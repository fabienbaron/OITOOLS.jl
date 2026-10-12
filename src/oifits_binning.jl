# Averaging spectral channels ("binning"): N adjacent channels, or a wavelength window, become
# ONE channel. Not grouping: `readoifits`' `spectralbin` gathers channels into polychromatic
# bins and keeps every one; this replaces them.
#
# Applied by the table readers while a table is still channels × rows, before uv is computed,
# so a row's channels are averaged together whatever a later filter removes.
#
# Two modes.
#
# :block (the default). Exact partitions: every channel is used once, and the new channel's
# spectral response is a top-hat whose width is the span averaged — what the smearing model
# assumes of a channel. Amplitudes (V², |V|, T3amp, flux) are averaged with weights 1/σ², phases
# (VISφ, T3φ) by the weighted circular mean arg Σ e^{iφ}/σ_φ². The error is σ_mean·√c with
# σ_mean = (Σ 1/σ²)^(-1/2): c = 1 for independent channels, c = n for no reduction at all, about
# 2-3 for an oversampled spectrograph whose neighbouring channels are correlated. The scatter of
# the channels about a straight line in λ fitted within each bin is returned as a χ² (a smooth
# chromatic trend is signal, not noise), so whether the errors (and c) are realistic can be judged
# rather than assumed.
#
# :pmoired. PMOIRED's `binning` (pmoired/oifits.py, `_binOI` and `_binVec`), for reproducing its
# results on the same data: a regular grid of nλ ÷ n channels, each a mean under a Gaussian kernel
# of FWHM half the new spacing with weights k/σ, phases by the circular mean of the weighted unit
# phasors, errors 1/⟨1/σ⟩ (not reduced), valid when one channel is. Its declared bandwidth, n
# times the old one, is not the kernel's.

"""
    BinSpec(binning; mode = :block, corr = 1)

How to average spectral channels. `binning` is the number of adjacent channels per new channel,
or a vector of `[λ_min, λ_max]` windows in metres (block mode only). `mode` is `:block` or
`:pmoired`; `corr` is c, the number of channels per independent element (block mode). `stats`
collects the scatter χ² per observable, as `(chi2, dof)` sums over every averaged channel.
"""
struct BinSpec
    binning::Union{Int,Vector{Vector{Float64}}}
    mode::Symbol
    corr::Float64
    stats::Dict{Symbol,Vector{Float64}}
end

function BinSpec(binning; mode::Symbol = :block, corr::Real = 1)
    mode in (:block, :pmoired) || throw(ArgumentError("binmode must be :block or :pmoired, got :$mode"))
    corr >= 1 || throw(ArgumentError("bincorr must be ≥ 1 (channels per independent element), got $corr"))
    b = binning isa Integer ? Int(binning) : [Float64.(collect(w)) for w in binning]
    b isa Int && b < 1 && throw(ArgumentError("binning must be ≥ 1 channels, got $b"))
    b isa Int || mode === :block ||
        throw(ArgumentError("binning by wavelength windows is block mode only"))
    return BinSpec(b, mode, Float64(corr), Dict{Symbol,Vector{Float64}}())
end

# The new channels of one wavelength table: their wavelengths, bandwidths and, in block mode, the
# old channels each one averages. Windows that hold no channel are dropped.
_bin_channels(spec::BinSpec, λ::AbstractVector, Δλ::AbstractVector) =
    _bin_channels(spec.binning, spec.mode, Float64.(λ), Float64.(Δλ))

function _bin_channels(nb::Int, mode::Symbol, λ::Vector{Float64}, Δλ::Vector{Float64})
    n = length(λ)
    o = sortperm(λ)
    if mode === :pmoired
        λn = collect(range(minimum(λ), maximum(λ); length = max(1, n ÷ nb)))
        return λn, nb .* _interp(λn, λ[o], Δλ[o]), Vector{Vector{Int}}()
    end
    return _block_channels([o[i:min(i + nb - 1, n)] for i in 1:nb:n], λ, Δλ)
end

_bin_channels(windows::Vector{Vector{Float64}}, ::Symbol, λ::Vector{Float64}, Δλ::Vector{Float64}) =
    _block_channels(filter(!isempty, [findall(x -> w[1] <= x <= w[2], λ) for w in windows]), λ, Δλ)

function _block_channels(groups::Vector{Vector{Int}}, λ::Vector{Float64}, Δλ::Vector{Float64})
    λn = Float64[sum(λ[g]) / length(g) for g in groups]
    Δn = Float64[maximum(λ[g] .+ Δλ[g] ./ 2) - minimum(λ[g] .- Δλ[g] ./ 2) for g in groups]
    return λn, Δn, groups
end

# Piecewise-linear interpolation, flat beyond the ends (numpy's `interp`).
function _interp(x::AbstractVector{Float64}, xp::AbstractVector{Float64},
                 fp::AbstractVector{Float64})
    map(x) do xi
        xi <= xp[1] && return fp[1]
        xi >= xp[end] && return fp[end]
        j = searchsortedlast(xp, xi)
        t = (xi - xp[j]) / (xp[j+1] - xp[j])
        (1 - t) * fp[j] + t * fp[j+1]
    end
end

_valid(v, e, f) = !f && isfinite(v) && isfinite(e) && e > 0

"""
    _bin_table(spec, λ, Δλ, val, err, flag; phase = false) -> (λn, Δλn, val, err, flag)

Average the channels (rows of the channels × records matrices `val`, `err`, `flag`) of one
table. `phase = true` for a phase in degrees. In block mode, with `stat` set, the scatter χ² and
its degrees of freedom are added to `spec.stats[stat]`.
"""
function _bin_table(spec::BinSpec, λ, Δλ, val::AbstractMatrix{T}, err::AbstractMatrix{T},
                    flag::AbstractMatrix{Bool}; phase::Bool = false,
                    stat::Union{Nothing,Symbol} = nothing) where {T}
    λn, Δn, groups = _bin_channels(spec, λ, Δλ)
    spec.mode === :pmoired &&
        return (T.(λn), T.(Δn), _pmoired_bin(λn, λ, val, err, flag; phase)...)
    v, e, f, χ2, dof = _bin_block(groups, Float64.(λ), val, err, flag, spec.corr, phase)
    if stat !== nothing
        st = get!(spec.stats, stat, [0.0, 0.0])
        st[1] += χ2; st[2] += dof
    end
    return T.(λn), T.(Δn), v, e, f
end

# The block average of every record, behind a function barrier so it compiles for the concrete
# matrix types; returns the averaged values, errors and flags, and the scatter χ² and its dof.
function _bin_block(groups::Vector{Vector{Int}}, λ::Vector{Float64}, val::AbstractMatrix{T},
                    err::AbstractMatrix{T}, flag::AbstractMatrix{Bool}, corr::Float64,
                    phase::Bool) where {T}
    nr = size(val, 2)
    v = fill(T(NaN), length(groups), nr); e = fill(T(NaN), length(groups), nr)
    f = trues(length(groups), nr)
    χ2, dof = 0.0, 0.0
    for r in 1:nr, (k, g) in enumerate(groups)
        ok = [i for i in g if _valid(val[i, r], err[i, r], flag[i, r])]
        isempty(ok) && continue
        w = [1 / Float64(err[i, r])^2 for i in ok]
        W = sum(w)
        x = Float64.(val[ok, r])
        m = if phase
            rad = deg2rad.(x)
            rad2deg(angle(sum(w .* cis.(rad))))
        else
            sum(w .* x) / W
        end
        if length(ok) >= 3
            χ, ν = _detrended_chi2(λ[ok], x, w, phase ? m : nothing)
            χ2 += χ; dof += ν
        end
        c = min(corr, length(ok))
        v[k, r] = T(m); e[k, r] = T(sqrt(c / W)); f[k, r] = false
    end
    return v, e, f, χ2, dof
end

# The scatter of a bin's channels about a straight line in λ (weighted), with n - 2 degrees of
# freedom. About a line rather than the mean, because a spectrum or a visibility that changes
# smoothly across the bin is signal, not noise; what is left is the channel-to-channel scatter
# the errors should describe. A phase is unwrapped about the bin's mean `m` first.
function _detrended_chi2(l::Vector{Float64}, x::Vector{Float64}, w::Vector{Float64}, m)
    y = m === nothing ? x : mod.(x .- m .+ 180, 360) .- 180
    W = sum(w); lm = sum(w .* l) / W; ym = sum(w .* y) / W
    sll = sum(w .* (l .- lm) .^ 2)
    b = sll > 0 ? sum(w .* (l .- lm) .* (y .- ym)) / sll : 0.0
    r = y .- ym .- b .* (l .- lm)
    return sum(w .* r .^ 2), length(x) - 2
end

# PMOIRED's `_binVec_flag` / `_binVec`, transcribed: per record, the valid channels only, a
# Gaussian kernel of FWHM half the new grid's local spacing, weights k/σ.
function _pmoired_bin(λn, λ, val::AbstractMatrix{T}, err::AbstractMatrix{T},
                      flag::AbstractMatrix{Bool}; phase::Bool = false) where {T}
    m, nr = length(λn), size(val, 2)
    v = fill(T(NaN), m, nr); e = fill(T(NaN), m, nr); f = trues(m, nr)
    G = _gradient(λn)
    for r in 1:nr
        ok = [i for i in axes(val, 1) if !flag[i, r]]
        any_ok = [any(!flag[i, r] for i in axes(val, 1) if
                      exp(-(λ[i] - λn[j])^2 / (2 * (0.5 * G[j] / 2.35482)^2)) > 0) for j in 1:m]
        isempty(ok) && continue
        x = Float64.(val[ok, r]); s = Float64.(err[ok, r]); l = Float64.(λ[ok])
        for j in 1:m
            k = exp.(-(l .- λn[j]) .^ 2 ./ (2 * (0.5 * G[j] / 2.35482)^2))
            no = sum(k ./ s)
            y = if phase
                (no != 0 && isfinite(no)) ? rad2deg(angle(sum(k ./ s .* cis.(deg2rad.(x))) / no)) :
                                            rad2deg(angle(sum(k .* cis.(deg2rad.(x))) / sum(k)))
            else
                (no != 0 && isfinite(no)) ? sum(k ./ s .* x) / no : sum(k .* x) / sum(k)
            end
            # the error: 1 / (the same weighted mean of 1/σ)
            inv = (no != 0 && isfinite(no)) ? sum(k ./ s .* (1 ./ s)) / no : sum(k .* (1 ./ s)) / sum(k)
            v[j, r] = T(y); e[j, r] = T(1 / inv)
            f[j, r] = !any_ok[j] || !isfinite(y)
        end
        if phase
            col = Float64.(v[:, r])
            good = isfinite.(col)
            any(good) && (v[good, r] = T.(_unwrap360(col[good] .+ 180) .- 180))
        end
    end
    return v, e, f
end

# numpy's `gradient` for a 1-D grid: one-sided at the ends, centred inside.
function _gradient(x::AbstractVector)
    n = length(x)
    n == 1 && return [one(eltype(x))]
    g = similar(x, Float64)
    g[1] = x[2] - x[1]; g[n] = x[n] - x[n-1]
    for i in 2:n-1
        g[i] = (x[i+1] - x[i-1]) / 2
    end
    return g
end

# numpy's `unwrap(…, period = 360)`.
function _unwrap360(p::AbstractVector)
    out = copy(p)
    for i in 2:length(out)
        d = out[i] - out[i-1]
        out[i] -= 360 * round(d / 360)
    end
    return out
end

"""
    binning_scatter(spec::BinSpec) -> Dict{Symbol,Float64}

The scatter check of a block-mode read: per observable, χ²/dof of the channels about a straight
line in λ fitted within each new channel (bins of three channels or more), summed over every new
channel. Near 1, the file's errors describe the
channel-to-channel scatter and c = 1 is right. Well below 1, neighbouring channels agree better
than their errors say: the channels are correlated (raise `bincorr`) or the errors inflated. Well
above 1, the errors are underestimated.
"""
binning_scatter(spec::BinSpec) =
    Dict(k => (s[2] > 0 ? s[1] / s[2] : NaN) for (k, s) in spec.stats)
