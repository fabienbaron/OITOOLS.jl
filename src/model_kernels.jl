# model_kernels.jl
#
# Component visibilities and their closed-form derivatives, point by point: the kernels the
# model Jacobian is assembled from (PLAN_ModelJacobian.md, J2), and the value-only versions
# `eval_model` uses.
#
# Every limb-darkening law in OITOOLS is a sum of the same building block,
#
#     B_ν(ζ) = ζ^{-ν} J_ν(ζ),      ζ = π θ ρ / MAS2RAD,
#
# with coefficients linear in the law's coefficients c:
#
#     V = Σ_k a_k(c) B_{ν_k}(ζ) / N(c),
#
# so its derivatives are
#
#     ∂V/∂ζ   = Σ_k a_k(c) (−ζ B_{ν_k+1}(ζ)) / N            (d/dζ[ζ^{-ν} J_ν] = −ζ^{-ν} J_{ν+1})
#     ∂V/∂c_j = (Σ_k ∂a_k/∂c_j B_{ν_k} − V ∂N/∂c_j) / N
#
# and ∂ζ/∂θ = π ρ / MAS2RAD, ∂ζ/∂ρ = π θ / MAS2RAD. The uniform disc is the one-term case.
#
# B_ν is where the time goes. A non-integer-order `besselj` costs ~350 ns a point against ~20 ns
# for `besselj1` (measured, 2026-10-08), so: a series below ζ = 2, where its alternating terms
# stay small; closed forms for the half-integer orders above it (spherical Bessel functions);
# `besselj0`/`besselj1` for orders 1 and 2. The quarter orders, and ldpow's moving order, use a
# longer series up to ζ = 12 and Hankel's asymptotic expansion beyond (`_hankel_asym`), both of
# which also give ∂/∂ν exactly; the library only for orders above `_ASYM_NU_MAX`.
#
# Every function here is generic in the number type, so ForwardDiff can differentiate the
# value kernels — which is how the closed-form derivatives are tested.

using SpecialFunctions: besselj, besselj0, besselj1, gamma, digamma
using ForwardDiff

# ── B_ν ──────────────────────────────────────────────────────────────────────────────────

# B_ν(0) = 1/(2^ν Γ(ν+1)).
_bnu0(ν) = 1 / (2^ν * gamma(ν + 1))

# The series in q = −x²/4: B_ν(x) = Σ_k c_k q^k, c_k = 1/(2^ν k! Γ(ν+k+1)), thirteen terms. Below
# x = 2 (|q| ≤ 1) the thirteenth is ~1e-19 of the first and nothing cancels. The coefficients
# are computed once per order, when the method is generated, and the sum is Horner's rule.
@generated function _bnu_series(x, ::Val{ν}) where {ν}
    c = Vector{Float64}(undef, 13)
    c[1] = _bnu0(ν)
    for k in 1:12
        c[k + 1] = c[k] / (k * (ν + k))
    end
    return :(evalpoly(-(x * x) / 4, $(Tuple(c))))
end

const _SQRT2OPI = sqrt(2 / π)

# Each order used, at x ≥ 2 (the series covers x < 2). Half-integer orders in closed form:
#   J_{3/2}(x) = √(2/(πx)) (sin x / x − cos x)
#   J_{5/2}(x) = √(2/(πx)) ((3/x² − 1) sin x − 3 cos x / x)
@inline _bnu_large(::Val{1.0}, x) = besselj1(x) / x
@inline _bnu_large(::Val{2.0}, x) = (2 * besselj1(x) / x - besselj0(x)) / (x * x)
# J₃ = 4J₂/x − J₁, one upward step from J₂ = 2J₁/x − J₀: above x = 2 the step is stable.
@inline function _bnu_large(::Val{3.0}, x)
    j1 = besselj1(x)
    j2 = 2 * j1 / x - besselj0(x)
    return (4 * j2 / x - j1) / (x * x * x)
end
@inline function _bnu_large(::Val{1.5}, x)
    s, c = sincos(x)
    return _SQRT2OPI * (s - x * c) / (x * x * x)
end
@inline function _bnu_large(::Val{2.5}, x)
    s, c = sincos(x)
    return _SQRT2OPI * ((3 - x * x) * s - 3 * x * c) / (x * x * x * x * x)
end
# Up to ζ = 12 the series still holds: its terms peak near (ζ/2)^ζ/Γ(ζ/2)² and the sum loses
# ~1e-13 of B_ν(0) at ζ = 12 (measured against 512-bit arithmetic), so the cancellation is
# harmless there. Coefficients as in `_bnu_series`, as many as the last term needs to fall
# below 1e-18 of the first at ζ = 12.
const _SERIES_XMAX = 12.0
@generated function _bnu_series_long(x, ::Val{ν}) where {ν}
    c = [_bnu0(ν)]
    k = 0
    while abs(c[end]) * 36.0^k > 1e-18 * c[1] && k < 80
        k += 1
        push!(c, c[end] / (k * (ν + k)))
    end
    return :(evalpoly(-(x * x) / 4, $(Tuple(c))))
end

# Hankel's expansion (DLMF 10.17.3), J_ν(x) = √(2/(πx)) (P cos ω − Q sin ω), ω = x − (ν/2 + 1/4)π,
# with P = Σ (−1)^m a_{2m} x^{−2m}, Q = Σ (−1)^m a_{2m+1} x^{−2m−1},
# a_k = a_{k−1} (4ν² − (2k−1)²)/(8k). The a_k are polynomials in ν, so ∂J_ν/∂ν comes term by
# term: ∂a_k = ∂a_{k−1} (4ν² − (2k−1)²)/(8k) + a_{k−1} ν/k, and ∂ω/∂ν = −π/2. The sum stops at
# its smallest term (it is asymptotic) or once both series are below 1e-17. For x ≥ 12 and
# ν ≤ `_ASYM_NU_MAX` that is within ~5e-11 of B_ν(0), relative, and ~1e-13 at ν ≤ 3 (measured
# against 512-bit arithmetic). A half-integer ν ends the a_k series; ∂a_k carries on, which
# is why the stopping rule watches both.
const _ASYM_NU_MAX = 7.0
@inline function _hankel_asym(ν, x)
    a = one(ν + x); da = zero(a); P = one(a); dP = zero(a); Q = zero(a); dQ = zero(a)
    xk = one(a); last = oftype(abs(a), Inf)
    for k in 1:40
        c  = (4 * ν * ν - (2k - 1)^2) / (8k)
        da = da * c + a * ν / k
        a *= c
        xk /= x
        t, dt = a * xk, da * xk
        mag = max(abs(t), abs(dt))
        mag > last && break
        last = mag
        sg = iseven(k ÷ 2) ? 1 : -1
        if isodd(k)
            Q += sg * t; dQ += sg * dt
        else
            P += sg * t; dP += sg * dt
        end
        mag < 1e-17 && break
    end
    sω, cω = sincos(x - (ν / 2 + oftype(ν, 0.25)) * π)
    A = sqrt(2 / (π * x))
    return (A * (P * cω - Q * sω), A * (dP * cω - dQ * sω + (P * sω + Q * cω) * π / 2))
end

@inline function _bnu_large(::Val{ν}, x) where {ν}
    ν > _ASYM_NU_MAX && return besselj(ν, x) / x^ν
    x < _SERIES_XMAX && return _bnu_series_long(x, Val(ν))
    return first(_hankel_asym(ν, x)) / x^ν
end

"""
    _bnu(Val(ν), x)

`x^{-ν} J_ν(x)` for a fixed order `ν`, accurate to Float64 resolution at every `x ≥ 0`.
"""
@inline _bnu(::Val{ν}, x) where {ν} = x < 2 ? _bnu_series(x, Val(ν)) : _bnu_large(Val(ν), x)

const _B1   = Val(1.0)
const _B2   = Val(2.0)
const _B3   = Val(3.0)
const _B32  = Val(1.5)
const _B52  = Val(2.5)
const _B54  = Val(1.25)
const _B94  = Val(2.25)
const _B74  = Val(1.75)
const _B114 = Val(2.75)

# The laws' constant basis factors (as in model_chainrules.jl): √(π/2), and wⱼΓ(νⱼ+1)2^νⱼ for
# b = 1/2 and 3/2.
const _K32 = sqrt(π / 2)
const _K54 = gamma(2.25) * 2^1.25 / 2.5
const _K74 = gamma(2.75) * 2^1.75 / 3.5

# At ζ below this the laws return exactly 1 and zero derivatives, as the value kernels always
# have (the limit is 1; the guard keeps it exact).
const _ZETA_MIN = 1e-8

# ── the laws, one point at a time ────────────────────────────────────────────────────────
#
# `_X_v(ζ, c...)` the value; `_X_d(ζ, c...)` the value, ∂/∂ζ, and ∂/∂c for each coefficient.

@inline _ud_v(ζ) = ζ < _ZETA_MIN ? one(ζ) : 2 * _bnu(_B1, ζ)
@inline _ud_d(ζ) = ζ < _ZETA_MIN ? (one(ζ), zero(ζ)) :
                                   (2 * _bnu(_B1, ζ), -2 * ζ * _bnu(_B2, ζ))

@inline function _ldlin_v(ζ, u)
    ζ < _ZETA_MIN && return one(ζ + u)
    return ((1 - u) * _bnu(_B1, ζ) + u * _K32 * _bnu(_B32, ζ)) / (0.5 - u / 6)
end
@inline function _ldlin_d(ζ, u)
    z = zero(ζ + u)
    ζ < _ZETA_MIN && return (one(z), z, z)
    N = 0.5 - u / 6
    b1, b32 = _bnu(_B1, ζ), _bnu(_B32, ζ)
    V = ((1 - u) * b1 + u * _K32 * b32) / N
    dζ = ((1 - u) * (-ζ * _bnu(_B2, ζ)) + u * _K32 * (-ζ * _bnu(_B52, ζ))) / N
    du = ((-b1 + _K32 * b32) + V / 6) / N
    return (V, dζ, du)
end

@inline function _ldquad_v(ζ, u, w)
    ζ < _ZETA_MIN && return one(ζ + u + w)
    num = (1 - u - w) * _bnu(_B1, ζ) + (u + 2w) * _K32 * _bnu(_B32, ζ) - 2w * _bnu(_B2, ζ)
    return num / (0.5 - u / 6 - w / 12)
end
@inline function _ldquad_d(ζ, u, w)
    z = zero(ζ + u + w)
    ζ < _ZETA_MIN && return (one(z), z, z, z)
    N = 0.5 - u / 6 - w / 12
    b1, b32, b2 = _bnu(_B1, ζ), _bnu(_B32, ζ), _bnu(_B2, ζ)
    V = ((1 - u - w) * b1 + (u + 2w) * _K32 * b32 - 2w * b2) / N
    dζ = ((1 - u - w) * (-ζ * b2) + (u + 2w) * _K32 * (-ζ * _bnu(_B52, ζ)) -
          2w * (-ζ * _bnu(_B3, ζ))) / N
    du = ((-b1 + _K32 * b32) + V / 6) / N
    dw = ((-b1 + 2 * _K32 * b32 - 2 * b2) + V / 12) / N
    return (V, dζ, du, dw)
end

@inline function _ldsqrt_v(ζ, c, d)
    ζ < _ZETA_MIN && return one(ζ + c + d)
    num = (1 - c - d) * _bnu(_B1, ζ) + c * _K32 * _bnu(_B32, ζ) + d * _K54 * _bnu(_B54, ζ)
    return num / (0.5 - c / 6 - d / 10)
end
@inline function _ldsqrt_d(ζ, c, d)
    z = zero(ζ + c + d)
    ζ < _ZETA_MIN && return (one(z), z, z, z)
    N = 0.5 - c / 6 - d / 10
    b1, b32, b54 = _bnu(_B1, ζ), _bnu(_B32, ζ), _bnu(_B54, ζ)
    V = ((1 - c - d) * b1 + c * _K32 * b32 + d * _K54 * b54) / N
    dζ = -ζ * ((1 - c - d) * _bnu(_B2, ζ) + c * _K32 * _bnu(_B52, ζ) +
               d * _K54 * _bnu(_B94, ζ)) / N
    dc = ((-b1 + _K32 * b32) + V / 6) / N
    dd = ((-b1 + _K54 * b54) + V / 10) / N
    return (V, dζ, dc, dd)
end

@inline function _ldclaret4_v(ζ, c1, c2, c3, c4)
    ζ < _ZETA_MIN && return one(ζ + c1 + c2 + c3 + c4)
    num = (1 - c1 - c2 - c3 - c4) * _bnu(_B1, ζ) + c1 * _K54 * _bnu(_B54, ζ) +
          c2 * _K32 * _bnu(_B32, ζ) + c3 * _K74 * _bnu(_B74, ζ) + c4 * 2 * _bnu(_B2, ζ)
    return num / (0.5 - c1 / 10 - c2 / 6 - 3c3 / 14 - c4 / 4)
end
@inline function _ldclaret4_d(ζ, c1, c2, c3, c4)
    z = zero(ζ + c1 + c2 + c3 + c4)
    ζ < _ZETA_MIN && return (one(z), z, z, z, z, z)
    N = 0.5 - c1 / 10 - c2 / 6 - 3c3 / 14 - c4 / 4
    b1, b54, b32, b74, b2 = _bnu(_B1, ζ), _bnu(_B54, ζ), _bnu(_B32, ζ), _bnu(_B74, ζ), _bnu(_B2, ζ)
    V = ((1 - c1 - c2 - c3 - c4) * b1 + c1 * _K54 * b54 + c2 * _K32 * b32 +
         c3 * _K74 * b74 + c4 * 2 * b2) / N
    dζ = -ζ * ((1 - c1 - c2 - c3 - c4) * b2 + c1 * _K54 * _bnu(_B94, ζ) +
               c2 * _K32 * _bnu(_B52, ζ) + c3 * _K74 * _bnu(_B114, ζ) +
               c4 * 2 * _bnu(_B3, ζ)) / N
    d1 = ((-b1 + _K54 * b54) + V / 10) / N
    d2 = ((-b1 + _K32 * b32) + V / 6) / N
    d3 = ((-b1 + _K74 * b74) + V * 3 / 14) / N
    d4 = ((-b1 + 2 * b2) + V / 4) / N
    return (V, dζ, d1, d2, d3, d4)
end

# ── vector kernels: W, ∂W/∂params (columns in the parser's order), ∂W/∂ρ ──────────────
#
# A parameter is a scalar, or a vector with one value per uv point (a chromatic parameter).

@inline _at(x::Real, i) = x
@inline _at(x::AbstractVector, i) = x[i]

"""
    _vis_jac!(W, D, Dρ, Val(kind), params, ρ)

The visibility `W` of an analytic component of `kind`, its derivatives with respect to the
kind's parameters (`D[:, j]`, in the order of `_ANALYTIC_PARAM_SUFFIXES[kind]`) and with
respect to ρ (`Dρ`, for the inclination/PA chain), at every ρ, in place.
"""
function _vis_jac!(W, D, Dρ, ::Val{:ud}, p, ρ)
    @inbounds for i in eachindex(ρ)
        θ = _at(p[1], i); r = ρ[i]
        V, dζ = _ud_d(π * θ * r / MAS2RAD)
        W[i] = V; D[i, 1] = dζ * π * r / MAS2RAD; Dρ[i] = dζ * π * θ / MAS2RAD
    end
    return W
end

for (K, vd, n) in ((:ldlin, :_ldlin_d, 1), (:ldquad, :_ldquad_d, 2), (:ldsqrt, :_ldsqrt_d, 2),
                   (:ldclaret4, :_ldclaret4_d, 4))
    cs = [:(_at(p[$(j + 1)], i)) for j in 1:n]
    @eval function _vis_jac!(W, D, Dρ, ::Val{$(QuoteNode(K))}, p, ρ)
        @inbounds for i in eachindex(ρ)
            θ = _at(p[1], i); r = ρ[i]
            out = $vd(π * θ * r / MAS2RAD, $(cs...))
            W[i] = out[1]
            D[i, 1] = out[2] * π * r / MAS2RAD
            Dρ[i] = out[2] * π * θ / MAS2RAD
            for j in 1:$n
                D[i, j + 1] = out[j + 2]
            end
        end
        return W
    end
end

# Gaussian: V = exp(−t²/ln2), t = π (FWHM/2) ρ / MAS2RAD (see `_vis_gaussian` for the half-width).
@inline _gauss_v(fwhm, r) = exp(-(π * (fwhm / 2) * r / MAS2RAD)^2 / log(2))
function _vis_jac!(W, D, Dρ, ::Val{:gaussian}, p, ρ)
    @inbounds for i in eachindex(ρ)
        f = _at(p[1], i); r = ρ[i]
        t = π * (f / 2) * r / MAS2RAD
        V = exp(-t^2 / log(2))
        g = -2t / log(2) * V                     # dV/dt
        W[i] = V; D[i, 1] = g * π * r / (2 * MAS2RAD); Dρ[i] = g * π * (f / 2) / MAS2RAD
    end
    return W
end

# Uniform ring: (dout² V_ud(dout) − din² V_ud(din)) / (dout² − din²).
function _vis_jac!(W, D, Dρ, ::Val{:ring}, p, ρ)
    @inbounds for i in eachindex(ρ)
        a = _at(p[1], i); b = _at(p[2], i); r = ρ[i]      # inner, outer diameter
        Va, ga = _ud_d(π * a * r / MAS2RAD); Vb, gb = _ud_d(π * b * r / MAS2RAD)
        Δ = b^2 - a^2
        V = (b^2 * Vb - a^2 * Va) / Δ
        W[i] = V
        D[i, 1] = (-2a * Va - a^2 * ga * π * r / MAS2RAD + 2a * V) / Δ
        D[i, 2] = (2b * Vb + b^2 * gb * π * r / MAS2RAD - 2b * V) / Δ
        Dρ[i] = (b^2 * gb * π * b / MAS2RAD - a^2 * ga * π * a / MAS2RAD) / Δ
    end
    return W
end

# Gaussian ring: (σo² Vo − σi² Vi) / (σo² − σi²), σ = FWHM/2.
function _vis_jac!(W, D, Dρ, ::Val{:gaussian_ring}, p, ρ)
    @inbounds for i in eachindex(ρ)
        fi = _at(p[1], i); fo = _at(p[2], i); r = ρ[i]
        si, so = fi / 2, fo / 2
        ti, to = π * si * r / MAS2RAD, π * so * r / MAS2RAD
        Vi, Vo = exp(-ti^2 / log(2)), exp(-to^2 / log(2))
        gi, go = -2ti / log(2) * Vi, -2to / log(2) * Vo       # dV/dt
        Δ = so^2 - si^2
        V = (so^2 * Vo - si^2 * Vi) / Δ
        W[i] = V
        # d/dσ, then d/dFWHM = d/dσ / 2
        dsi = (-2si * Vi - si^2 * gi * π * r / MAS2RAD + 2si * V) / Δ
        dso = (2so * Vo + so^2 * go * π * r / MAS2RAD - 2so * V) / Δ
        D[i, 1] = dsi / 2; D[i, 2] = dso / 2
        Dρ[i] = (so^2 * go * π * so / MAS2RAD - si^2 * gi * π * si / MAS2RAD) / Δ
    end
    return W
end

function _vis_jac!(W, D, Dρ, ::Val{:point}, p, ρ)
    fill!(W, one(eltype(W))); fill!(Dρ, zero(eltype(Dρ)))
    return W
end

function _vis_jac!(W, D, Dρ, ::Val{:resolved}, p, ρ)
    fill!(W, zero(eltype(W))); fill!(Dρ, zero(eltype(Dρ)))
    return W
end

# ── power-law limb darkening: the order moves with α ─────────────────────────────────────
#
# V = Γ(ν+1) 2^ν B_ν(ζ), ν = α/2 + 1. Below ζ = 12 as the series
#     V = Σ_k q^k t_k,  q = −ζ²/4,  t_k = 1/(k! (ν+1)(ν+2)…(ν+k)),
# whose ν-derivative is exact term by term: ∂t_k/∂ν = −t_k Σ_{j≤k} 1/(ν+j); it stops once a
# term is below 1e-17 (V is O(1)). Above it Hankel's expansion (`_hankel_asym`), with
#     ∂V/∂ν = V (ψ(ν+1) + log 2 − log ζ) + Γ(ν+1) 2^ν ∂J_ν/∂ν / ζ^ν,
# and the library with `_dbesselj_dnu` (model_chainrules.jl) for ν + 1 > `_ASYM_NU_MAX`.
# `G` = Γ(ν+1) 2^ν and `ψ` = ψ(ν+1) + log 2 are passed in so a loop at fixed α computes them
# once. Returns (V, ∂V/∂ζ, ∂V/∂α).
@inline _ldpow_consts(ν) = (gamma(ν + 1) * 2^ν, digamma(ν + 1) + log(oftype(ν, 2)))

@inline function _ldpow_series(ζ, ν)
    q  = -(ζ * ζ) / 4
    t  = one(ζ + ν); s = zero(t)          # t_k and Σ_{j≤k} 1/(ν+j)
    V  = one(t); dq = zero(t); dν = zero(t)
    qk = one(t)                           # q^k
    for k in 1:80
        t  /= k * (ν + k)
        s  += 1 / (ν + k)
        dq += k * t * qk                  # Σ k t_k q^{k−1}
        qk *= q
        tq  = t * qk
        V  += tq
        dν -= tq * s
        abs(tq) < 1e-17 && break
    end
    return (V, dq * (-ζ / 2), dν / 2)
end

@inline function _ldpow_d(ζ, α, G, ψ)
    ν = α / 2 + 1
    z = zero(ζ + α)
    ζ < _ZETA_MIN && return (one(z), z, z)
    ζ < _SERIES_XMAX && return _ldpow_series(ζ, ν)
    zν = ζ^ν
    if ν + 1 > _ASYM_NU_MAX
        Jν, dJ, Jν1 = besselj(ν, ζ), _dbesselj_dnu(ν, ζ), besselj(ν + 1, ζ)
    else
        Jν, dJ = _hankel_asym(ν, ζ)
        Jν1    = first(_hankel_asym(ν + 1, ζ))
    end
    V = G * Jν / zν
    return (V, -G * Jν1 / zν, (V * (ψ - log(ζ)) + G * dJ / zν) / 2)
end
@inline _ldpow_d(ζ, α) = _ldpow_d(ζ, α, _ldpow_consts(α / 2 + 1)...)

# The value alone, for `vis_ldpow`; `G` = Γ(ν+1) 2^ν as in `_ldpow_d`.
@inline _ldpow_v(ζ, α) = _ldpow_v(ζ, α, first(_ldpow_consts(α / 2 + 1)))
@inline function _ldpow_v(ζ, α, G)
    ν = α / 2 + 1
    ζ < _ZETA_MIN && return one(ζ + α)
    if ζ < _SERIES_XMAX                   # `_ldpow_series`'s recurrence, value only
        q = -(ζ * ζ) / 4
        t = one(ζ + ν); V = one(t); qk = one(t)
        for k in 1:80
            t  /= k * (ν + k)
            qk *= q
            tq  = t * qk
            V  += tq
            abs(tq) < 1e-17 && break
        end
        return V
    end
    J = ν > _ASYM_NU_MAX ? besselj(ν, ζ) : first(_hankel_asym(ν, ζ))
    return G * J / ζ^ν
end

function _vis_jac!(W, D, Dρ, ::Val{:ldpow}, p, ρ)
    isempty(ρ) && return W
    αc = _at(p[2], firstindex(ρ))
    G, ψ = _ldpow_consts(αc / 2 + 1)
    @inbounds for i in eachindex(ρ)
        θ = _at(p[1], i); r = ρ[i]; α = _at(p[2], i)
        if α != αc                       # a chromatic α
            αc = α; G, ψ = _ldpow_consts(α / 2 + 1)
        end
        V, dζ, dα = _ldpow_d(π * θ * r / MAS2RAD, α, G, ψ)
        W[i] = V; D[i, 1] = dζ * π * r / MAS2RAD; D[i, 2] = dα; Dρ[i] = dζ * π * θ / MAS2RAD
    end
    return W
end

# ── crescent: two offset uniform discs, complex, in (u, v) ───────────────────────────────
#
# W = (b² V_b(ρ) e^{iφ} − a² V_a(ρ) e^{−iφ}) / (b² − a²), a = crin, b = crout,
# φ = off·g, off = croff (b − a)/2, g = π (u sin PA + v cos PA)/MAS2RAD (PA = crprojang).
# Derivatives in (crin, crout, croff, crprojang) into the columns of D, and ∂W/∂u, ∂W/∂v.
function _crescent_jac!(W, D, Du, Dv, p, u, v)
    @inbounds for i in eachindex(u)
        a = _at(p[1], i); b = _at(p[2], i); cf = _at(p[3], i); pa = _at(p[4], i) * π / 180
        ui, vi = u[i], v[i]
        r = sqrt(ui^2 + vi^2)
        sp, cp = sincos(pa)
        g   = π * (ui * sp + vi * cp) / MAS2RAD
        off = cf * (b - a) / 2
        Eo, Ei = cis(off * g), cis(-off * g)
        Va, ga = _ud_d(π * a * r / MAS2RAD); Vb, gb = _ud_d(π * b * r / MAS2RAD)
        Δ  = b^2 - a^2
        Wi = (b^2 * Vb * Eo - a^2 * Va * Ei) / Δ
        W[i] = Wi
        dWdoff = im * g * (b^2 * Vb * Eo + a^2 * Va * Ei) / Δ
        dWdg   = im * off * (b^2 * Vb * Eo + a^2 * Va * Ei) / Δ
        D[i, 1] = (-2a * Va * Ei - a^2 * ga * π * r / MAS2RAD * Ei) / Δ + Wi * 2a / Δ -
                  dWdoff * cf / 2
        D[i, 2] = (2b * Vb * Eo + b^2 * gb * π * r / MAS2RAD * Eo) / Δ - Wi * 2b / Δ +
                  dWdoff * cf / 2
        D[i, 3] = dWdoff * (b - a) / 2
        D[i, 4] = dWdg * π * (ui * cp - vi * sp) / MAS2RAD * π / 180
        # through ρ (both discs) and through g
        dWdρ = r == 0 ? zero(Wi) :
               (b^2 * gb * π * b / MAS2RAD * Eo - a^2 * ga * π * a / MAS2RAD * Ei) / Δ
        Du[i] = (r == 0 ? zero(Wi) : dWdρ * ui / r) + dWdg * π * sp / MAS2RAD
        Dv[i] = (r == 0 ? zero(Wi) : dWdρ * vi / r) + dWdg * π * cp / MAS2RAD
    end
    return W
end

# ── a written visibility function: dual numbers through the user's expression ────────────
#
# The expression is arbitrary Julia, so it is differentiated rather than derived: one call of
# the compiled function on dual numbers, one partial for the baseline B (Mλ) and one per
# parameter, broadcast over all points at once. The one place AD remains in the Jacobian.
struct _VisFuncTag end

function _visfunc_jac!(W, D, DB, vis_fn, params, B, wl, mjd)
    _visfunc_jac!(W, D, DB, vis_fn, params, B, wl, mjd, Val(length(params) + 1))
end

function _visfunc_jac!(W, D, DB, vis_fn, params, B, wl, mjd, ::Val{NP}) where {NP}
    T  = promote_type(eltype(B), map(x -> eltype(x), params)...)
    Dl = ForwardDiff.Dual{_VisFuncTag,T,NP}
    seed(x, k) = Dl(T(x), ForwardDiff.Partials(ntuple(j -> j == k ? one(T) : zero(T), Val(NP))))
    Bd = [seed(b, 1) for b in B]
    pd = ntuple(j -> params[j] isa AbstractVector ? [seed(x, j + 1) for x in params[j]] :
                                                    seed(params[j], j + 1), length(params))
    Vd = vis_fn(Bd, wl, mjd, pd...)
    _visfunc_unpack!(W, D, DB, Vd, B, Val(NP))
    return W
end

function _visfunc_unpack!(W, D, DB, Vd, B, ::Val{NP}) where {NP}
    n = length(B)
    @inbounds for i in 1:n
        z = Vd isa AbstractVector ? Vd[i] : Vd      # a constant expression returns a scalar
        re, imv = real(z), imag(z)
        W[i]  = complex(ForwardDiff.value(re), ForwardDiff.value(imv))
        DB[i] = complex(ForwardDiff.partials(re, 1), ForwardDiff.partials(imv, 1))
        for j in 2:NP
            D[i, j - 1] = complex(ForwardDiff.partials(re, j), ForwardDiff.partials(imv, j))
        end
    end
    return W
end

# ── radial profiles: the Hankel transform as a cached matrix ─────────────────────────────
#
#     H_n(B_j) = Σ_k w_k I_k J_n(2π B_j r_k) r_k = (K_n I)_j,  K_n[j,k] = w_k r_k J_n(2π B_j r_k)
#
# K_n depends on the baselines and the nodes, not on the profile, so for a component with no
# inclination or position angle (whose B is the data's own) it is built once per uv set and
# the transform is a matrix-vector product. With the nodes r = r_max·s and weights g of
# `_hankel_nodes`, K_n[j,k] = g_k J_n(2π B_j r_k); its B-derivative, for the inclination/PA
# chain and the azimuthal modes, is K'_n[j,k] = g_k 2π r_k J_n'(2π B_j r_k), J_0' = −J_1,
# J_n' = (J_{n−1} − J_{n+1})/2.

"""
    _HankelLayout

What the probe (`_hankel_probe`, hankel_nodes.jl) found in a radial profile, in s = r/r_max:
its breakpoints and their kinds (`:jump`, `:kink`, `:transition`), whether the centre is
singular and with what exponent, and for each breakpoint its derivative with respect to each
profile parameter (zero for one that does not move).
"""
struct _HankelLayout
    breaks :: Vector{Float64}
    kinds  :: Vector{Symbol}
    sing   :: Bool
    pexp   :: Float64
    dbdp   :: Matrix{Float64}          # length(breaks) × number of profile parameters
end

"""
    HankelKernelCache{T}()

A radial-profile component's quadrature nodes, per node count, and per uv set (compared by
value) its node count and transform matrices, per order. The last few uv sets are kept. Also
what the probe found in the profile, and the multiple of the automatic count the convergence
check settled on, for baselines that move (inclination, position angle). Guarded by a lock,
since a `FlatModel` is shared between the threads of a batch likelihood.
"""
mutable struct HankelKernelCache{T<:AbstractFloat}
    entries :: Vector{Any}          # _HankelUV{TB,TK}: baselines of any element type
    nodes   :: Dict{Tuple{Int,Int},NTuple{3,Vector{T}}}  # (n, n₀) => (r, μ, g), r = r_max·s
    layout  :: Union{Nothing,_HankelLayout}  # nothing = not probed yet
    mult    :: Int                  # the check's multiple of `_hankel_n`
    warned  :: Bool                 # the not-converged warning was issued
    last_pp :: Vector{Float64}      # the parameters the breakpoints were last located for
    breaks  :: Vector{Float64}      # and where they were: the nodes and matrices are for these
    nbuilt  :: Int
    lock    :: ReentrantLock
end
HankelKernelCache{T}() where {T} =
    HankelKernelCache{T}(Any[], Dict{Tuple{Int,Int},NTuple{3,Vector{T}}}(), nothing, 1, false,
                         Float64[], Float64[], 0, ReentrantLock())

# One uv set: its baselines, the node count chosen for them, and the matrices built for them,
# keyed by (order, derivative).
struct _HankelUV{TB,TK}
    B  :: Vector{TB}
    n  :: Int
    Ks :: Dict{Tuple{Int,Bool},Matrix{TK}}
end

const _HANKEL_CACHE_SIZE = 4

@inline _besselj_int(n, x) = n == 0 ? besselj0(x) : n == 1 ? besselj1(x) : besselj(n, x)

function _build_hkernel(r::AbstractVector, g::AbstractVector, B::AbstractVector, n::Int,
                        deriv::Bool)
    Tk = promote_type(eltype(B), eltype(r), eltype(g))
    twoπ = Tk(2π)
    K = Matrix{Tk}(undef, length(B), length(r))
    @inbounds for k in eachindex(r), j in eachindex(B)
        x = twoπ * B[j] * r[k]
        K[j, k] = if !deriv
            g[k] * _besselj_int(n, x)
        else
            Jp = n == 0 ? -besselj1(x) : (_besselj_int(n - 1, x) - _besselj_int(n + 1, x)) / 2
            g[k] * twoπ * r[k] * Jp
        end
    end
    return K
end

# K_n (or K'_n) for these baselines: from the cache when `cache` is given, built otherwise.
# With a cache the nodes are a function of the baselines (`_hankel_setup`: a component whose
# nodes move with its parameters passes no cache), so the baselines alone key a uv set.
function _hkernel(cache, r, g, B::AbstractVector{TB}, n::Int, deriv::Bool) where {TB}
    cache === nothing && return _build_hkernel(r, g, B, n, deriv)
    TK = promote_type(TB, eltype(r), eltype(g))
    lock(cache.lock) do
        for e in cache.entries
            e isa _HankelUV{TB,TK} || continue
            if length(e.B) == length(B) && e.B == B
                return get!(() -> (cache.nbuilt += 1; _build_hkernel(r, g, B, n, deriv)),
                            e.Ks, (n, deriv))
            end
        end
        length(cache.entries) >= _HANKEL_CACHE_SIZE && popfirst!(cache.entries)
        K = _build_hkernel(r, g, B, n, deriv); cache.nbuilt += 1
        push!(cache.entries, _HankelUV{TB,TK}(Vector{TB}(B), length(r), Dict((n, deriv) => K)))
        return K
    end
end

"""
    _hankel_jac!(W, D, DB, Dψ, r, g, I, dI, B, ψ, modes, amps, phis, cache, jumps=())

A radial-profile component with optional azimuthal modes: its visibility `W` (complex), the
derivatives `D` — columns: the profile's parameters (through `dI = ∂I/∂p`, `n × np`), then
(amplitude, angle in degrees) for each mode — and `DB = ∂W/∂B`, `Dψ = ∂W/∂ψ` for the
inclination/PA chain. `r`, `g` are the nodes and weights of `_hankel_setup`. `modes` are the
orders; `amps`, `phis` their values (scalars or per point). `cache` is a `HankelKernelCache`,
or `nothing` when B changes every call.

`jumps` are the profile's jumps that move with its parameters, as `(bs, br, ΔI, dbdp)`: the
place in s and in mas, I(b⁺) − I(b⁻), and ∂b/∂p in s. ∂I/∂p does not see them; each adds
−ΔI·Jₙ(2πB b)·b·∂b/∂p to ∂Hₙ/∂p and −ΔI·b·∂b/∂p to ∂N/∂p (in s ds, as g is).
"""
function _hankel_jac!(W, D, DB, Dψ, r, g, I, dI, B, ψ, modes, amps, phis, cache, jumps = ())
    np = size(dI, 2)
    N  = dot(g, I)
    dN = dI' * g                                     # ∂N/∂p
    for b in jumps
        dN .-= b.ΔI * b.bs .* b.dbdp
    end
    K0 = _hkernel(cache, r, g, B, 0, false)
    V0 = (K0 * I) ./ N
    dV = (K0 * dI .+ _hankel_jumpterm(0, B, jumps, np, eltype(K0)) .- V0 .* dN') ./ N
    @views D[:, 1:np] .= dV
    W  .= V0
    DB .= (_hkernel(cache, r, g, B, 0, true) * I) ./ N
    fill!(Dψ, zero(eltype(Dψ)))
    for (m, n) in enumerate(modes)
        Hn  = (_hkernel(cache, r, g, B, n, false) * I) ./ N
        dHn = (_hkernel(cache, r, g, B, n, false) * dI .+
               _hankel_jumpterm(n, B, jumps, np, eltype(Hn)) .- Hn .* dN') ./ N
        HnB = (_hkernel(cache, r, g, B, n, true) * I) ./ N
        ph  = (-im)^n
        @inbounds for j in eachindex(B)
            a = _at(amps[m], j); φ = _at(phis[m], j) * π / 180
            s, c = sincos(n * (ψ[j] + φ - π / 2))
            W[j]  += a * ph * Hn[j] * c
            DB[j] += a * ph * HnB[j] * c
            Dψ[j] += -a * ph * Hn[j] * n * s
            for l in 1:np
                D[j, l] += a * ph * dHn[j, l] * c
            end
            D[j, np + 2m - 1] = ph * Hn[j] * c
            D[j, np + 2m]     = -a * ph * Hn[j] * n * s * π / 180
        end
    end
    return W
end

# ∂Hₙ/∂p from jumps that move with the parameters (see `_hankel_jac!`), `nB × np`.
function _hankel_jumpterm(n, B, jumps, np, ::Type{Tk}) where {Tk}
    M = zeros(Tk, length(B), np)
    for b in jumps
        twoπbr = Tk(2π) * Tk(b.br)
        @inbounds for j in eachindex(B)
            c = -Tk(b.ΔI) * _besselj_int(n, twoπbr * B[j]) * Tk(b.bs)
            for l in 1:np
                M[j, l] += c * Tk(b.dbdp[l])
            end
        end
    end
    return M
end
