# The measurement model: how the model predicts what the instrument delivered, acting on the
# model only (`PLAN_C1.md`). Here, the reference of a differential VISPHI: per spectral vector,
# a polynomial in λ fitted to the phase residuals over the CONTINUUM channels — the vector's
# channels away from the model's lines, or within given continuum ranges — weighted by the data's
# errors, and subtracted at every channel. It is PMOIRED's DPHI residual,
# (d − P d) − (m − P m) = (d − m) − P (d − m), with the data's own errors, and each vector's
# coefficients counted out of the degrees of freedom.

"""
    MeasurementModel

How a fit references the model to the data (see [`measurement_model`](@ref)).
"""
struct MeasurementModel
    visphi    :: Symbol                          # :absolute, or :continuum
    order     :: Int                             # polynomial order; -1 = PMOIRED's "auto"
    continuum :: Vector{Tuple{Float64,Float64}}  # λ ranges (m) the continuum is taken from; empty: all
    kernel    :: Float64                         # PMOIRED's spectral kernel, FWHM in channels; 0 none
    pmoired   :: Bool                            # PMOIRED's DPHI as it computes it (see below)
end

"""
    measurement_model(; visphi = :absolute, order = :auto, continuum = nothing,
                      spectral_kernel = 0, pmoired = false) -> MeasurementModel

- `visphi = :continuum` — each spectral vector's phase residuals minus a polynomial in λ of
  `order` (`:auto`: PMOIRED's max(2, ⌊ptp(λ_continuum)/0.2 µm⌋)) fitted over its continuum
  channels: those outside the model's spectral lines (± 1.5 FWHM for a Gaussian, ± 3 widths for a
  Lorentzian, PMOIRED's widths), and inside `continuum` (λ ranges in metres) when given.
  `:absolute` (default) compares VISPHI as it is.
- `spectral_kernel = k` — PMOIRED's `wl kernel`: before the comparison (and before
  `visphi = :continuum`), the model is smoothed along each spectral vector with a Gaussian of
  FWHM `k` channels (2⌊2k⌋ + 3 taps, normalised, edges mirrored): the complex visibility and the
  complex triple product as phasors, V² separately. PMOIRED's setupFit uses k = 2 for GRAVITY.
  Only the channels present are smoothed; at the edges of a filtered window the values are
  mirrored, while PMOIRED smooths across them.
- `pmoired = true` — `visphi = :continuum` computed as PMOIRED's DPHI (`computeDiffPhiOI`): the
  continuum is every channel (Gaussian and Lorentzian lines included, `truncexp` lines
  excluded), and the data's polynomial is weighted by 1/σ² while the model's is unweighted:
  (d − P_w d) − (m − P_u m). The default fits one weighted polynomial to the residuals over the
  channels away from the lines.
"""
function measurement_model(; visphi::Symbol = :absolute, order = :auto, continuum = nothing,
                           spectral_kernel::Real = 0, pmoired::Bool = false)
    visphi in (:absolute, :continuum) || throw(ArgumentError("visphi is :absolute or :continuum, not :$visphi"))
    o = order === :auto ? -1 : Int(order)
    o >= -1 || throw(ArgumentError("order is :auto or ≥ 0"))
    c = continuum === nothing ? Tuple{Float64,Float64}[] :
        [(Float64(first(r)), Float64(last(r))) for r in continuum]
    spectral_kernel >= 0 || throw(ArgumentError("spectral_kernel is a FWHM in channels, ≥ 0"))
    return MeasurementModel(visphi, o, c, Float64(spectral_kernel), pmoired)
end

_references(mm) = mm !== nothing && (mm.visphi === :continuum || mm.kernel > 0)
_smooths(ref) = ref !== nothing && ref.mm.kernel > 0

# What the referencing needs of a dataset and a model, built once per fit: the VIS spectral
# vectors (positions in the VIS arrays, by wavelength), each point's λ (µm), and where the
# model's lines keep their centre and width among the resolver's values.
struct _RefContext
    mm     :: MeasurementModel
    groups :: Vector{Vector{Int}}
    lam    :: Vector{Float64}                   # µm, per VIS point
    lines  :: Vector{Tuple{Int,Int,Symbol}}     # (centre, width, profile) indices in all_names
    sg     :: SpectralGroups                    # every observable's vectors, for the kernel
end

function _ref_context(mm::MeasurementModel, fm, data)
    lines = Tuple{Int,Int,Symbol}[]
    n2i = Dict(n => i for (i, n) in enumerate(fm.all_names))
    for n in fm.all_names
        m = match(r"^(.+),line_(.+)_wl0$", n)
        m === nothing && continue
        for p in (:gaussian, :lorentzian, :truncexp, :trunclor)
            k = "$(m[1]),line_$(m[2])_$p"
            haskey(n2i, k) && push!(lines, (n2i[n], n2i[k], p))
        end
    end
    sg = spectral_groups(data)
    return _RefContext(mm, sg.vis, Float64.(data.vis_lam) .* 1e6, lines, sg)
end

# The continuum channels of a vector (positions within it), PMOIRED's exclusions around the lines.
function _continuum_mask(ctx::_RefContext, g, vals)
    λ = view(ctx.lam, g)
    keep = trues(length(g))
    if !isempty(ctx.mm.continuum)
        keep .= [any(lo * 1e6 <= l <= hi * 1e6 for (lo, hi) in ctx.mm.continuum) for l in λ]
    end
    for (ic, iw, p) in ctx.lines
        c, w = _plain(vals[ic]), abs(_plain(vals[iw]))
        if ctx.mm.pmoired && p !== :truncexp
            continue                       # PMOIRED's DPHI leaves these lines in the continuum
        elseif p === :gaussian
            keep .&= abs.(λ .- c) .>= 1.5w / 1000
        elseif p === :lorentzian
            keep .&= abs.(λ .- c) .>= 3w / 1000
        elseif p === :truncexp
            d = 2w / 1000
            keep .&= .!((λ .<= c + 1.3d) .& (λ .>= c - d / 8))
        end
    end
    return keep
end

# The vector's polynomial order: the given one, or PMOIRED's "auto" from the continuum's span.
_ref_order(ctx, λc) = ctx.mm.order >= 0 ? ctx.mm.order :
    max(2, floor(Int, (isempty(λc) ? 0.0 : maximum(λc) - minimum(λc)) / 0.2))

"""
    _reference_phases!(r, J, ctx, vals, σ, sw) -> number of coefficients removed

The VISPHI rows `r` (and their Jacobian `J`, or `nothing`) of one dataset, referenced per spectral
vector: the phase residuals unwrapped along λ (as PMOIRED unwraps the phase before its fit), then
rᵢ − Bᵢβ with B = diag(sw/σ)·[1, t, t², …] (t the vector's λ centred and scaled) and β the
least-squares fit of r over the continuum channels; J likewise, column by column (unwrapping
adds constants, which it does not see). `vals` are the
resolver's values (the lines' centres and widths), `σ` the VISPHI errors, `sw` √weight. A vector
with fewer continuum channels than coefficients takes all its channels; with fewer channels
still, its mean.
"""
function _reference_phases!(r, J, ctx::_RefContext, vals, σ, sw, dphi = nothing)
    ncoef = 0
    for g in ctx.groups
        n = length(g)
        n == 0 && continue
        keep = _continuum_mask(ctx, g, vals)
        λ = ctx.lam[g]
        k = _ref_order(ctx, λ[keep])
        count(keep) <= k && (keep .= true)
        k = min(k, n - 1)
        μ, s = sum(λ) / n, max(maximum(λ) - minimum(λ), eps())
        t = (λ .- μ) ./ s
        # unwrap the residuals (in degrees: r σ / sw) along the vector
        for i in 2:n
            a, b = g[i - 1], g[i]
            (isfinite(σ[a]) && isfinite(σ[b]) && σ[b] > 0) || continue
            Δ = r[b] * σ[b] / sw - r[a] * σ[a] / sw
            r[b] -= 360 * round(Δ / 360) * sw / σ[b]
        end
        if ctx.mm.pmoired && dphi !== nothing
            # PMOIRED: (d − P_w d) − (m − P_u m) = δ − P_u δ + (P_w − P_u) d, δ = m − d, in
            # degrees; the data's phase unwrapped along the vector, as PMOIRED unwraps it
            fin = [isfinite(σ[i]) && σ[i] > 0 for i in g]
            keep .&= fin
            count(keep) <= k && (keep .= fin)
            A = [t[i]^p for i in 1:n, p in 0:k]
            dd = Float64.(dphi[g])
            for i in 2:n; dd[i] = dd[i - 1] + mod(dd[i] - dd[i - 1] + 180, 360) - 180; end
            δ = [fin[i] ? r[g[i]] * σ[g[i]] / sw : 0.0 for i in 1:n]
            Ac = A[keep, :]
            Wc = [1 / σ[g[i]] for i in 1:n][keep]
            e = δ .- A * (Ac \ δ[keep]) .+ A * ((Ac .* Wc) \ (dd[keep] .* Wc) .- Ac \ dd[keep])
            for i in 1:n
                fin[i] && (r[g[i]] = sw * e[i] / σ[g[i]])
            end
            if J !== nothing
                Jδ = [fin[i] ? J[g[i], j] * σ[g[i]] / sw : 0.0 for i in 1:n, j in axes(J, 2)]
                Je = Jδ .- A * (Ac \ Jδ[keep, :])
                for i in 1:n
                    fin[i] && (J[g[i], :] .= sw .* Je[i, :] ./ σ[g[i]])
                end
            end
            ncoef += k + 1
            continue
        end
        B = [sw / σ[g[i]] * t[i]^p for i in 1:n, p in 0:k]
        Bc = B[keep, :]
        F = qr(Bc)
        β = F \ r[g[keep]]
        r[g] .-= B * β
        if J !== nothing
            Jc = J[g[keep], :]
            J[g, :] .-= B * (F \ Jc)
        end
        ncoef += k + 1
    end
    return ncoef
end

# ── PMOIRED's spectral kernel (`_applyWlKernel`, `_convolve`) ────────────────────────────

# The kernel: a Gaussian of FWHM `k` channels over 2⌊2k⌋ + 3 taps, normalised.
function _pm_kernel(k::Real)
    n = 2 * floor(Int, 2k) + 3
    x = (0:(n - 1)) .- (n - 1) / 2
    g = exp.(-x .^ 2 ./ (2 * (k / 2.35482)^2))
    return g ./ sum(g)
end

# PMOIRED's `_convolve`: `y` mirrored at both ends over as many samples as the kernel has (or
# `y` has), convolved ("same"), the padding cut off. Linear in `y`, generic in its element type.
function _pm_convolve(y::AbstractVector, ker::AbstractVector)
    n, m = length(y), length(ker)
    k = min(m, n)
    p = [reverse(y[1:k]); y; reverse(y[(n - k + 1):n])]
    h = (m - 1) ÷ 2
    out = similar(y, promote_type(eltype(y), eltype(ker)), n)
    for i in 1:n
        acc = zero(eltype(out))
        c = i + k                                    # centre in the padded array
        for j in 1:m
            q = c - h + j - 1
            1 <= q <= length(p) && (acc += ker[m - j + 1] * p[q])
        end
        out[i] = acc
    end
    return out
end

# A quantity per point (and its Jacobian rows, or nothing), smoothed within each vector `groups`.
function _pm_smooth!(a, da, groups, ker)
    for g in groups
        # (numpy's "same" changes length when the padded vector is shorter than the kernel;
        # such a short vector is left as it is)
        (length(g) < 2 || 3length(g) < length(ker)) && continue
        a[g] = _pm_convolve(a[g], ker)
        if da !== nothing
            for j in axes(da, 2)
                da[g, j] = _pm_convolve(da[g, j], ker)
            end
        end
    end
    return a
end

# The residual rows (and Jacobian rows, or `J = nothing`) of the observables formed from the
# model smoothed by PMOIRED's kernel: V, V² and the triple product per spectral vector, then the
# rows as `residuals_flat` forms them. The flux rows are the unsmoothed ones (the flux term
# compares one model flux).
function _smoothed_rows!(r, J, model, x, data, w, ref, jac_cache)
    ker = _pm_kernel(ref.mm.kernel)
    _wl  = hasproperty(data, :uv_lam) ? data.uv_lam : nothing
    _mjd = hasproperty(data, :uv_mjd) ? data.uv_mjd : nothing
    if J === nothing
        V = eval_model(model, x, data.uv; wl = _wl, mjd = _mjd); dV = nothing
    else
        V, dV = _model_jac(jac_cache, model, x, data.uv; wl = _wl, mjd = _mjd)
    end
    sub(idx) = V[idx]
    dsub(idx) = dV === nothing ? nothing : dV[idx, :]
    np = length(x)
    sg = ref.sg
    # V² per V² point, smoothed on its own
    v2 = abs2.(sub(data.indx_v2))
    dv2 = dV === nothing ? nothing : 2 .* real.(conj.(sub(data.indx_v2)) .* dsub(data.indx_v2))
    _pm_smooth!(v2, dv2, sg.v2, ker)
    # the complex triple product, smoothed as a phasor
    V1, V2, V3 = sub(data.indx_t3_1), sub(data.indx_t3_2), sub(data.indx_t3_3)
    t3 = V1 .* V2 .* V3
    dt3 = dV === nothing ? nothing :
          (V2 .* V3) .* dsub(data.indx_t3_1) .+ (V1 .* V3) .* dsub(data.indx_t3_2) .+ (V1 .* V2) .* dsub(data.indx_t3_3)
    _pm_smooth!(t3, dt3, sg.t3, ker)
    # the complex visibility, smoothed as a phasor
    vis = sub(data.indx_vis); dvis = dsub(data.indx_vis)
    _pm_smooth!(vis, dvis, sg.vis, ker)
    ε = eps(Float64)
    k = 0
    put!(val, dval, d, e, s) = for i in eachindex(val)
        k += 1
        r[k] = s * (val[i] - d[i]) / e[i]
        J === nothing || (J[k, :] .= s .* dval[i, :] ./ e[i])
    end
    putphase!(c, dc, d, e, s) = for i in eachindex(c)
        k += 1
        r[k] = s * _mod360_1(angle(c[i]) * (180.0 / π) - d[i]) / e[i]
        J === nothing || (J[k, :] .= s .* (180.0 / π) .* imag.(conj(c[i]) ./ max(abs2(c[i]), ε) .* dc[i, :]) ./ e[i])
    end
    putamp!(c, dc, d, e, s) = for i in eachindex(c)
        k += 1
        a = abs(c[i])
        r[k] = s * (a - d[i]) / e[i]
        J === nothing || (J[k, :] .= s .* real.(conj(c[i]) ./ max(a, ε) .* dc[i, :]) ./ e[i])
    end
    w[1] > 0 && data.nv2 > 0 && put!(v2, dv2, data.v2, data.v2_err, sqrt(w[1]))
    w[2] > 0 && data.nt3amp > 0 && putamp!(t3, dt3, data.t3amp, data.t3amp_err, sqrt(w[2]))
    w[3] > 0 && data.nt3phi > 0 && putphase!(t3, dt3, data.t3phi, data.t3phi_err, sqrt(w[3]))
    w[4] > 0 && data.nvisamp > 0 && putamp!(vis, dvis, data.visamp, data.visamp_err, sqrt(w[4]))
    w[5] > 0 && data.nvisphi > 0 && putphase!(vis, dvis, data.visphi, data.visphi_err, sqrt(w[5]))
    if w[6] > 0 && data.nflux > 0
        throw(ArgumentError("spectral_kernel: the flux term is not smoothed yet; leave OI_FLUX out"))
    end
    return r
end
