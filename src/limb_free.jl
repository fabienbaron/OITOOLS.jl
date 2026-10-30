# limb_free.jl — the free-form limb-darkened disc: I(μ) set by its values at N nodes.
#
#   "star,ldspline"   = 8.4    diameter (mas); a natural cubic spline in μ through the nodes
#   "star,ldpiecewise"= 8.4    the same with straight lines between the nodes
#   "star,ld1" … "star,ld<N-1>"  I(μ_k)/I(1) at μ_k = (k-1)/(N-1); the node at μ = 1 is 1
#
# N is the number of `ld<k>` keys plus one. The profile is linear in the node values,
#     I(μ) = Σ_k c_k b_k(μ),    c_N = 1,
# so V(ζ) = Σ c_k G_k(ζ) / Σ c_k G_k(0), with G_k(ζ) = ∫₀¹ b_k(μ) J₀(ζ√(1-μ²)) μ dμ and
# ζ = πθρ. G_k depends on the basis and N only: it is tabulated once per (basis, N) on a grid
# in ζ with its derivative and read by cubic Hermite interpolation, whose derivative is the
# one returned, so the value and the gradient a fitter sees agree.

using Printf

# ── basis on the nodes ───────────────────────────────────────────────────────────────────

"""
    ld_nodes(n) -> Vector{Float64}

The μ of a free-form disc's `n` nodes: evenly spaced from 0 (the limb) to 1 (the centre).
"""
ld_nodes(n::Integer) = collect(range(0.0, 1.0; length = n))

# Second derivatives of the natural cubic spline through `y` on `x`.
function _natural_m(x::AbstractVector{Float64}, y::AbstractVector)
    n = length(x); h = diff(x)
    m = zeros(promote_type(Float64, eltype(y)), n)
    n < 3 && return m
    A = zeros(n - 2, n - 2); r = zeros(eltype(m), n - 2)
    for i in 2:(n - 1)
        A[i-1, i-1] = 2 * (h[i-1] + h[i])
        i > 2     && (A[i-1, i-2] = h[i-1])
        i < n - 1 && (A[i-1, i] = h[i])
        r[i-1] = 6 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1])
    end
    m[2:(n-1)] = A \ r
    return m
end

# The profile through the node values `y` at μ (any number type in y).
function _ld_interp(basis::Symbol, x::AbstractVector{Float64}, y::AbstractVector, μ::Real)
    n = length(x)
    k = clamp(searchsortedlast(x, μ), 1, n - 1)
    h = x[k+1] - x[k]
    if basis === :piecewise
        a = (μ - x[k]) / h
        return (1 - a) * y[k] + a * y[k+1]
    end
    m = _natural_m(x, y)
    a = x[k+1] - μ; b = μ - x[k]
    return (m[k] * a^3 + m[k+1] * b^3) / (6h) + (y[k] / h - m[k] * h / 6) * a +
           (y[k+1] / h - m[k+1] * h / 6) * b
end

# b_k at each of `μs`, as a length(μs) × n matrix.
function _ld_basis(basis::Symbol, n::Int, μs::AbstractVector{Float64})
    x = ld_nodes(n)
    B = zeros(length(μs), n)
    for k in 1:n
        e = zeros(n); e[k] = 1.0
        if basis === :piecewise
            for (q, μ) in enumerate(μs)
                B[q, k] = _ld_interp(:piecewise, x, e, μ)
            end
        else
            m = _natural_m(x, e)
            for (q, μ) in enumerate(μs)
                j = clamp(searchsortedlast(x, μ), 1, n - 1); h = x[j+1] - x[j]
                a = x[j+1] - μ; b = μ - x[j]
                B[q, k] = (m[j] * a^3 + m[j+1] * b^3) / (6h) + (e[j] / h - m[j] * h / 6) * a +
                          (e[j+1] / h - m[j+1] * h / 6) * b
            end
        end
    end
    return B
end

# ── the tables ───────────────────────────────────────────────────────────────────────────

struct _LDFreeTable
    h    :: Float64              # step in ζ
    zmax :: Float64              # last tabulated ζ
    g    :: Matrix{Float64}      # G_k(ζ_j), (nζ, n)
    dg   :: Matrix{Float64}      # G_k'(ζ_j)
    g0   :: Vector{Float64}      # G_k(0)
end

const _LDFREE_H = 0.02
const _LDFREE_ZMAX = 512.0       # beyond it G is integrated at each ζ instead of read
const _LDFREE_ZCAP = 1.0e4       # beyond it V is taken as 0 (|V| < 1e-5 for a disc whose limb is not brighter than its centre)
const _LDFREE_TABLES = Dict{Tuple{Symbol,Int},Base.RefValue{_LDFreeTable}}()
const _LDFREE_LOCK = ReentrantLock()

# Quadrature in μ good to ζ ≤ zmax: each interval between nodes is cut into pieces even in
# s = √(1-μ²), each short enough in ζs to need few points, 16 Gauss-Legendre points a piece.
# J₀(ζ√(1-μ²)) is a function of μ², so it is smooth in μ up to μ = 1.
function _ld_quadrature(n::Int, zmax::Float64)
    xg, wg = gausslegendre(16)
    x = ld_nodes(n)
    μq = Float64[]; wq = Float64[]
    for k in 1:(n - 1)
        sa, sb = sqrt(1 - x[k]^2), sqrt(1 - x[k+1]^2)
        m = ceil(Int, zmax * (sa - sb) / 3) + 1
        for j in 1:m
            s1 = sa + (sb - sa) * (j - 1) / m; s2 = sa + (sb - sa) * j / m
            μ1, μ2 = sqrt(max(0.0, 1 - s1^2)), sqrt(max(0.0, 1 - s2^2))
            for (t, w) in zip(xg, wg)
                push!(μq, (μ1 + μ2) / 2 + (μ2 - μ1) / 2 * t)
                push!(wq, (μ2 - μ1) / 2 * w)
            end
        end
    end
    return μq, wq
end

function _build_ld_table(basis::Symbol, n::Int, zmax::Float64)
    μq, wq = _ld_quadrature(n, zmax)
    B = _ld_basis(basis, n, μq)
    sq = sqrt.(max.(0.0, 1 .- μq .^ 2))
    Bw = B .* (wq .* μq)                       # b_k μ dμ at each node
    nz = round(Int, zmax / _LDFREE_H) + 1
    g = zeros(nz, n); dg = zeros(nz, n)
    Threads.@threads for j in 1:nz
        ζ = (j - 1) * _LDFREE_H
        r0 = Vector{Float64}(undef, length(μq)); r1 = similar(r0)
        @inbounds for q in eachindex(μq)
            r0[q] = besselj0(ζ * sq[q]); r1[q] = -sq[q] * besselj1(ζ * sq[q])
        end
        g[j, :] = Bw' * r0
        dg[j, :] = Bw' * r1
    end
    return _LDFreeTable(_LDFREE_H, (nz - 1) * _LDFREE_H, g, dg, vec(sum(Bw; dims = 1)))
end

# The table for (basis, n) reaching ζ, built or extended (doubling) under the lock. Readers
# take the table once and keep it, so an extension replaces it without disturbing them.
function _ld_table(basis::Symbol, n::Int, ζ::Real)::_LDFreeTable
    r = get(_LDFREE_TABLES, (basis, n), nothing)
    r !== nothing && (r[].zmax >= ζ || r[].zmax >= _LDFREE_ZMAX) && return r[]
    lock(_LDFREE_LOCK)
    try
        r2 = get(_LDFREE_TABLES, (basis, n), nothing)
        r2 !== nothing && (r2[].zmax >= ζ || r2[].zmax >= _LDFREE_ZMAX) && return r2[]
        zmax = r2 === nothing ? 32.0 : r2[].zmax
        while zmax < min(ζ, _LDFREE_ZMAX); zmax *= 2; end
        t = _build_ld_table(basis, n, min(zmax, _LDFREE_ZMAX))
        r2 === nothing ? (_LDFREE_TABLES[(basis, n)] = Ref(t)) : (r2[] = t)
        return t
    finally
        unlock(_LDFREE_LOCK)
    end
end

_pv(x::ForwardDiff.Dual) = _pv(ForwardDiff.value(x))
_pv(x) = x

# G_k(ζ) and G_k'(ζ) for every k, into `G` and `dG`. Hermite on the table up to its end, the
# integral itself beyond.
function _ld_G!(G, dG, basis::Symbol, n::Int, t::_LDFreeTable, ζ)
    z = _pv(ζ)
    if !isfinite(z)
        fill!(G, NaN); fill!(dG, NaN)
    elseif z > _LDFREE_ZCAP
        fill!(G, zero(eltype(G))); fill!(dG, zero(eltype(dG)))
    elseif z <= t.zmax
        j = min(floor(Int, z / t.h), size(t.g, 1) - 2) + 1
        h = t.h
        s = (ζ - (j - 1) * h) / h
        h00 = 2s^3 - 3s^2 + 1; h10 = s^3 - 2s^2 + s; h01 = -2s^3 + 3s^2; h11 = s^3 - s^2
        d00 = (6s^2 - 6s) / h; d10 = 3s^2 - 4s + 1; d01 = (-6s^2 + 6s) / h; d11 = 3s^2 - 2s
        @inbounds for k in 1:n
            a, b, da, db = t.g[j, k], t.g[j+1, k], t.dg[j, k], t.dg[j+1, k]
            G[k]  = h00 * a + h10 * h * da + h01 * b + h11 * h * db
            dG[k] = d00 * a + d10 * da + d01 * b + d11 * db
        end
    else
        μq, wq = _ld_quadrature(n, Float64(z))
        Bw = _ld_basis(basis, n, μq) .* (wq .* μq)
        sq = sqrt.(max.(0.0, 1 .- μq .^ 2))
        for k in 1:n
            G[k]  = sum(Bw[q, k] * besselj0(ζ * sq[q]) for q in eachindex(μq))
            dG[k] = -sum(Bw[q, k] * sq[q] * besselj1(ζ * sq[q]) for q in eachindex(μq))
        end
    end
    return G, dG
end

# ── kernels ──────────────────────────────────────────────────────────────────────────────

_ld_basis_of(K::Symbol) = K === :ldpiecewise ? :piecewise : :spline

# The largest finite ζ of a call, which the table has to reach.
function _ld_reach(θ, ρ)
    z = 0.0
    for i in eachindex(ρ)
        v = _pv(π * _at(θ, i) * ρ[i] / MAS2RAD)
        isfinite(v) && (z = max(z, Float64(v)))
    end
    return min(z, _LDFREE_ZMAX)
end

# Values only, for `eval_model`: params = [θ, c_1 … c_{n-1}], each a scalar or one per ρ.
function _vis_ldfree(K::Symbol, params::AbstractVector, ρ::AbstractVector)
    n = length(params)
    basis = _ld_basis_of(K)
    T = promote_type(mapreduce(eltype, promote_type, params), eltype(ρ))
    V = Vector{T}(undef, length(ρ))
    isempty(ρ) && return V
    t = _ld_table(basis, n, _ld_reach(params[1], ρ))
    G = Vector{T}(undef, n); dG = Vector{T}(undef, n)
    @inbounds for i in eachindex(ρ)
        ζ = π * _at(params[1], i) * ρ[i] / MAS2RAD
        _ld_G!(G, dG, basis, n, t, ζ)
        num = G[n]; den = t.g0[n]
        for k in 1:(n - 1)
            c = _at(params[k+1], i)
            num += c * G[k]; den += c * t.g0[k]
        end
        V[i] = num / den
    end
    return V
end

_vis_ldspline(params::AbstractVector, ρ::AbstractVector)   = _vis_ldfree(:ldspline, params, ρ)
_vis_ldpiecewise(params::AbstractVector, ρ::AbstractVector) = _vis_ldfree(:ldpiecewise, params, ρ)

# The value, ∂/∂θ, ∂/∂c_k and ∂/∂ρ, as the other discs' `_vis_jac!`.
function _vis_jac_ldfree!(W, D, Dρ, K::Symbol, p, ρ)
    isempty(ρ) && return W
    n = length(p); basis = _ld_basis_of(K)
    t = _ld_table(basis, n, _ld_reach(p[1], ρ))
    T = eltype(W)
    G = Vector{T}(undef, n); dG = Vector{T}(undef, n)
    @inbounds for i in eachindex(ρ)
        θ = _at(p[1], i); r = ρ[i]
        _ld_G!(G, dG, basis, n, t, π * θ * r / MAS2RAD)
        num = G[n]; dnum = dG[n]; den = t.g0[n]
        for k in 1:(n - 1)
            c = _at(p[k+1], i)
            num += c * G[k]; dnum += c * dG[k]; den += c * t.g0[k]
        end
        V = num / den; dζ = dnum / den
        W[i] = V
        D[i, 1] = dζ * π * r / MAS2RAD
        Dρ[i] = dζ * π * θ / MAS2RAD
        for k in 1:(n - 1)
            D[i, k+1] = (G[k] - V * t.g0[k]) / den
        end
    end
    return W
end
_vis_jac!(W, D, Dρ, ::Val{:ldspline}, p, ρ)   = _vis_jac_ldfree!(W, D, Dρ, :ldspline, p, ρ)
_vis_jac!(W, D, Dρ, ::Val{:ldpiecewise}, p, ρ) = _vis_jac_ldfree!(W, D, Dρ, :ldpiecewise, p, ρ)

# ── keys ─────────────────────────────────────────────────────────────────────────────────

const _LD_NODE_KEY = r"^ld(\d+)$"
const _LD_RATIO_KEY = r"^ldq(\d+)$"

# The node keys a free-form disc `cn` has, ld1 … ld<k> while they run unbroken.
function _ld_node_suffixes(cn::AbstractString, pd::AbstractDict)
    out = String[]
    k = 1
    while haskey(pd, "$cn,ld$k"); push!(out, "ld$k"); k += 1; end
    return out
end

"""
    ld_profile(model_dict, component, μ) -> Vector{Float64}

The intensity across a disc, I(μ)/I(1), at each `μ`, for any disc law (`ud`, `ldlin`,
`ldquad`, `ldsqrt`, `ldpow`, `ldclaret4`) and for the free-form discs (`ldspline`,
`ldpiecewise`). Coefficients that are expressions are evaluated.
"""
function ld_profile(model_dict::AbstractDict, component::AbstractString, μ::AbstractVector{<:Real})
    md = Dict{String,Any}(String(k) => v for (k, v) in model_dict)
    cn = String(component)
    num(s) = _resolve_numeric("$cn,$s", md)
    has(s) = haskey(md, "$cn,$s")
    μv = Float64.(μ)
    if has("ldspline") || has("ldpiecewise")
        basis = has("ldspline") ? :spline : :piecewise
        y = [num(s) for s in _ld_node_suffixes(cn, md)]
        push!(y, 1.0)
        x = ld_nodes(length(y))
        return [_ld_interp(basis, x, y, m) for m in μv]
    elseif has("ldlin")
        u = num("u"); return @. 1 - u * (1 - μv)
    elseif has("ldquad")
        u, w = num("u"), num("w"); return @. 1 - u * (1 - μv) - w * (1 - μv)^2
    elseif has("ldsqrt")
        u, w = num("u"), num("w"); return @. 1 - u * (1 - μv) - w * (1 - sqrt(μv))
    elseif has("ldpow")
        α = num("alpha"); return @. μv^α
    elseif has("ldclaret4")
        c = [num("c$k") for k in 1:4]
        return [1 - sum(c[k] * (1 - m^(k / 2)) for k in 1:4) for m in μv]
    elseif has("ud")
        return ones(length(μv))
    end
    error("$cn is not a disc")
end

const _LD_LAW_KEYS = ("ud", "ldlin", "ldquad", "ldsqrt", "ldpow", "ldclaret4", "ldspline",
                      "ldpiecewise")
const _LD_LAW_COEFFS = Dict("ud" => String[], "ldlin" => ["u"], "ldquad" => ["u", "w"],
                            "ldsqrt" => ["u", "w"], "ldpow" => ["alpha"],
                            "ldclaret4" => ["c1", "c2", "c3", "c4"])

# Every key a disc's law owns, the diameter first: the law's coefficients, or a free-form
# disc's nodes and their ratios.
function _ld_law_keys(cn::AbstractString, pd::AbstractDict)
    i = findfirst(k -> haskey(pd, "$cn,$k"), _LD_LAW_KEYS)
    i === nothing && return String[]
    law = _LD_LAW_KEYS[i]
    rest = haskey(_LD_LAW_COEFFS, law) ? _LD_LAW_COEFFS[law] :
           [k for k in (s[(length(cn) + 2):end] for s in keys(pd) if startswith(s, cn * ","))
            if occursin(_LD_NODE_KEY, k) || occursin(_LD_RATIO_KEY, k)]
    return vcat([law], sort(rest))
end

"""
    ld_free!(model_dict, free, component; n = 5, basis = :spline, monotonic = false,
             lb = nothing, ub = nothing) -> model_dict

Make the disc `component` a free-form limb-darkened disc with `n` nodes, in place: its
diameter and whether it is free are kept, and the nodes take the values of the profile it had
(any law, or a free-form disc with another `n` or basis), so a fit continues from it.

`basis` is `:spline` (`c,ldspline`) or `:piecewise` (`c,ldpiecewise`). The nodes are
`c,ld1` … `c,ld<n-1>`, I(μ_k)/I(1) at μ_k = (k-1)/(n-1), all free.

With `monotonic = true` the free parameters are the ratios `c,ldq<k>` = I(μ_k)/I(μ_{k+1}) in
[0, 1], and each `c,ld<k>` is the expression multiplying them, so I(μ) never increases toward
the limb at the nodes. `lb` and `ub`, when given, receive the bounds of the new parameters and
lose those of the removed ones.
"""
function ld_free!(md::AbstractDict, free::AbstractVector, component::AbstractString;
                  n::Integer = 5, basis::Symbol = :spline, monotonic::Bool = false,
                  lb = nothing, ub = nothing)
    basis in (:spline, :piecewise) ||
        throw(ArgumentError("basis must be :spline or :piecewise, not :$basis"))
    n >= 2 || throw(ArgumentError("a free-form disc needs at least 2 nodes"))
    cn = String(component)
    old = _ld_law_keys(cn, md)
    isempty(old) && throw(ArgumentError("$cn is not a disc"))
    dkey = "$cn,$(old[1])"
    diam = md[dkey]
    wasfree = dkey in free
    dlb = lb === nothing ? nothing : get(lb, dkey, nothing)
    dub = ub === nothing ? nothing : get(ub, dkey, nothing)
    x = ld_nodes(n)
    y = clamp.(ld_profile(md, cn, x), 0.0, Inf)
    y ./= y[end] > 0 ? y[end] : 1.0
    for s in old
        k = "$cn,$s"
        delete!(md, k); filter!(!=(k), free)
        lb === nothing || delete!(lb, k); ub === nothing || delete!(ub, k)
    end
    nkey = "$cn," * (basis === :spline ? "ldspline" : "ldpiecewise")
    md[nkey] = diam
    wasfree && push!(free, nkey)
    dlb === nothing || (lb[nkey] = dlb)
    dub === nothing || (ub[nkey] = dub)
    for k in 1:(n - 1)
        key = "$cn,ld$k"
        if monotonic
            q = "$cn,ldq$k"
            md[q] = y[k+1] > 0 ? clamp(y[k] / y[k+1], 0.0, 1.0) : 1.0
            md[key] = k == n - 1 ? "\$$q" : "\$$q * \$$cn,ld$(k + 1)"
            push!(free, q)
            lb === nothing || (lb[q] = 0.0); ub === nothing || (ub[q] = 1.0)
        else
            md[key] = y[k]
            push!(free, key)
            lb === nothing || (lb[key] = 0.0); ub === nothing || (ub[key] = 2.0)
        end
    end
    return md
end

# ── error bars by variational inference (method in the VarInf extension) ────────────────

"""
    LDProfileVIResult

What [`fit_ld_profile_vi`](@ref) returns. `μ` is the grid; the profiles are `I(μ)/I(1)`
(`I_*`) and `I(μ)/⟨I⟩`, normalised by the mean intensity over the disc (`F_*`), which the
visibilities measure directly. `*_mean`/`*_std` are over the posterior samples (columns of
`I_samples`, `F_samples`); `*_lo`/`*_hi` their 16th and 84th percentiles. `diameter` and
`diameter_std` likewise. The samples are drawn about `centre`: `:mean`, the variational mean,
or `:map`, the posterior's maximum, when the variational mean diverged (χ²r more than twice the
maximum's). `model` is the model with the disc as a `ldpiecewise` disc on the grid at that
centre (a point of the posterior; the samples' mean need not be one), `free` its free
parameters, and `chi2r` the data's reduced χ² there; `chi2r_map` is the same at the
posterior's maximum, where the run starts. `error_scale` is the factor the data's error bars were multiplied by (1 unless
`scale_errors`). With bagging, `bootstrap` is the number of replicates pooled and `single` the
run on the full data alone; the samples and their statistics are the pool's, the centre,
`model` and `chi2r` the full data's. Without, `bootstrap` is 0 and `single` is `nothing`.
"""
struct LDProfileVIResult
    component    :: String
    method       :: Symbol
    prior        :: Symbol
    μ            :: Vector{Float64}
    I_mean       :: Vector{Float64}
    I_std        :: Vector{Float64}
    I_lo         :: Vector{Float64}
    I_hi         :: Vector{Float64}
    I_samples    :: Matrix{Float64}
    F_mean       :: Vector{Float64}
    F_std        :: Vector{Float64}
    F_lo         :: Vector{Float64}
    F_hi         :: Vector{Float64}
    F_samples    :: Matrix{Float64}
    diameter     :: Float64
    diameter_std :: Float64
    diameter_samples :: Vector{Float64}
    chi2r        :: Float64
    chi2r_map    :: Float64
    error_scale  :: Float64
    centre       :: Symbol
    ndata        :: Int
    model        :: Dict{String,Any}
    free         :: Vector{String}
    bootstrap    :: Int
    single       :: Union{Nothing,LDProfileVIResult}
    seconds      :: Float64
end

function Base.show(io::IO, ::MIME"text/plain", r::LDProfileVIResult)
    @printf(io, "LDProfileVIResult: %s, %s, %s prior, %d samples, χ²r = %.3f (%.1f s)%s\n",
            r.component, r.method, r.prior, size(r.I_samples, 2), r.chi2r, r.seconds,
            r.centre === :map ? "; drawn about the posterior's maximum" : "")
    r.bootstrap > 0 && @printf(io, "  bagged over %d bootstrap replicates\n", r.bootstrap)
    @printf(io, "  diameter = %.4f ± %.4f mas", r.diameter, r.diameter_std)
    r.error_scale == 1 ? println(io) : @printf(io, "   (error bars × %.3f)\n", r.error_scale)
    for i in round.(Int, range(1, length(r.μ); length = min(length(r.μ), 9)))
        @printf(io, "  μ = %.3f   I/I(1) = %.3f ± %.3f   I/⟨I⟩ = %.3f ± %.3f\n",
                r.μ[i], r.I_mean[i], r.I_std[i], r.F_mean[i], r.F_std[i])
    end
end

"""
    fit_ld_profile_vi(model_dict, component, data; kwargs...) -> LDProfileVIResult

The intensity profile I(μ) of the disc `component`, with error bars, by variational
inference (needs `using VarInf`). The disc's profile becomes a field on `ngrid` nodes even in
μ, joined by straight lines (`ldpiecewise`), log I = log I₀ + f, with I₀ the disc's current
profile (any law, or a free-form disc) and f a Gaussian process; the diameter is
log-normal about its current value. The run starts at the posterior's maximum, found from the
lowest minima of a scan of the diameter over ±22 %; posterior samples come from MGVI or geoVI.
The rest of the model is held at its values.

Keywords:
- `ngrid = 24`: nodes of the profile.
- `prior = :fixed`: f is a Gaussian process with a Matérn-3/2 kernel of length `length_scale`
  (in μ, default 0.3) and amplitude `amplitude` (in log I, default 0.3). `prior = :learned`:
  f is VarInf's correlated field, its power spectrum inferred, with the log-normal amplitude
  prior `fluct = (mean, std)` (default `(0.3, 0.2)`) and the slope prior `slope = (mean, std)`
  (default `(-4.0, 1.0)`), on a grid padded to twice its length.
- `diameter_std = 0.1`: the prior's relative width on the diameter.
- `method = :geovi` or `:mgvi`; `n_iterations = 6`, `n_samples = 4` per iteration;
  `n_draws = 50`: antithetic pairs drawn about the final mean for the statistics.
- `weights`: observables, as in [`fit_model`](@ref); by default V² and T3amp. Closure phases
  of a symmetric disc are 0 or 180° and add only the side of a null, which V² gives; their
  wrapped residuals are steps the samplers cannot follow.
- `scale_errors = false`: with `true`, error bars too small for the free-form profile's best
  fit (χ²r > 1 there) are multiplied by √χ²r before sampling. The band assumes the error bars;
  without this, underestimated ones make the profile follow the calibration's systematics.
- `bootstrap = 0`: with `n ≥ 2`, bagged VI. The run is repeated on `n` bootstrap replicates of
  the data (configurations drawn with replacement), each from the full data's maximum and
  with its error scale, and their samples pooled with equal weight, `max(4, cld(n_draws, n))`
  antithetic pairs each. The pool's spread is the posterior's plus the spread of its centre
  between replicates, which grows when the error bars or the model do not describe the data.
  The replicates run on all threads.
- `seed`, `verb`.
"""
function fit_ld_profile_vi end
