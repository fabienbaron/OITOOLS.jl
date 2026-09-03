# model_chainrules.jl
#
# Differentiable (ForwardDiff / Zygote / AutoMALA-compatible) implementations of the
# limb-darkened disk visibility functions, plus explicit ChainRules rrules for each. These are
# what `parse_model.jl` builds a FlatModel from.
#
# Imported here rather than inherited. These names used to arrive from a module-level `using`
# in a sibling file that happened to be included first; depending on include order for a name
# is how a reordering becomes an UndefVarError at precompile time.
using SpecialFunctions
#
# This file provides:
#   - Scalar-parameter variants  vis_ldXXX(θ, params..., ρ) returning a Vector
#   - rrules for those scalar variants
#   - Safe (ifelse-based) drop-in replacements for the full (param, uv) API
#
# Conventions
#   θ        : angular diameter in mas
#   ρ        : baseline / wavelength  (cycles / rad), vector of length nuv
#   C        : mas-to-radian conversion constant  =  360*3600*1000 / (2π)
#   ζ        : dimensionless argument  π θ ρ / C   (one entry per baseline)
#
# Bessel-order recurrence used throughout:
#   d/dx [ x^{-ν} J_ν(x) ] = -x^{-ν} J_{ν+1}(x)       … (R1)
#   J_ν(0) / 0^ν → 1 / (2^ν Γ(ν+1))  so  x^{-ν} J_ν(x) → 1/(2^ν Γ(ν+1))  as x→0

using ChainRulesCore

# ---------------------------------------------------------------------------
# Shared constant
# ---------------------------------------------------------------------------

const MAS2RAD = 2.0626480624709636e8   # 1 radian in milli-arcseconds

# ---------------------------------------------------------------------------
# Per-point parameters
# ---------------------------------------------------------------------------
#
# A parameter may be a scalar, or a vector holding one value per uv point —
# the latter is what the resolver produces for an expression referencing the
# implicit variables `$WL`, `$MJD` or `$B` (a chromatic diameter, a
# time-variable limb darkening, ...).  Every formula below is written with
# broadcasting, so the same body serves both cases; only the reduction in the
# pullback differs, hence `_dparam`.

const VisParam = Union{Real,AbstractVector}

# Cotangent for a parameter: summed over uv points if the parameter is shared
# (scalar), left element-wise if each uv point has its own value.
_dparam(θ::Real, v) = sum(v)
_dparam(::AbstractVector, v) = v

# ---------------------------------------------------------------------------
# 1.  Uniform disk
# ---------------------------------------------------------------------------
#
# V_ud(θ, ρ) = 2 J₁(t) / t ,  t = π θ ρ / C
#
# ∂V/∂θ  = (∂V/∂t) * (∂t/∂θ)
#         = [ 2(t J₀(t) - 2 J₁(t)) / t² ]  *  (π ρ / C)
#   (the analytic derivative of the uniform-disc visibility)
#
# Note: jinc(y) = 2 J₁(πy) / (πy), i.e. V = jinc(θ ρ / C)

function vis_ud(θ::VisParam, ρ::AbstractVector)
    t = @. π * θ * ρ / MAS2RAD
    return @. ifelse(t < 1e-8, one(t), 2 * besselj1(t) / t)
end

function ChainRulesCore.rrule(::typeof(vis_ud), θ::VisParam, ρ::AbstractVector)
    t    = @. π * θ * ρ / MAS2RAD
    safe = t .> 1e-8
    V    = @. ifelse(safe, 2 * besselj1(t) / t, one(t))

    function vis_ud_pullback(ȳ)
        dV_dt = @. ifelse(safe, 2*(t*besselj0(t) - 2*besselj1(t)) / t^2, zero(t))
        dt_dθ = @. π * ρ / MAS2RAD
        ∂θ    = _dparam(θ, ȳ .* dV_dt .* dt_dθ)
        return NoTangent(), ∂θ, NoTangent()   # ρ is treated as data
    end
    return V, vis_ud_pullback
end

# ---------------------------------------------------------------------------
# 2.  Linear limb-darkening
# ---------------------------------------------------------------------------
#
# Profile: I(μ) ∝ 1 - u(1 - μ)
#
# V(θ, u, ρ) = [ (1-u) B₁(ζ)  +  u √(π/2) B₃₂(ζ) ] / N(u)
#
#   B₁(ζ)   = J₁(ζ) / ζ
#   B₃₂(ζ)  = J_{3/2}(ζ) / ζ^{3/2}
#   N(u)    = 0.5 - u/6
#   ζ       = π θ ρ / C
#
# Derivatives
# -----------
# Via recurrence (R1):
#   dB₁/dζ  = -J₂(ζ)  / ζ
#   dB₃₂/dζ = -J_{5/2}(ζ) / ζ^{3/2}
#
# ∂V/∂θ = [ (1-u)(-J₂/ζ) + u√(π/2)(-J_{5/2}/ζ^{3/2}) ] / N  *  (πρ/C)
#
# ∂V/∂u  (quotient rule):
#   ∂V_num/∂u = -B₁ + √(π/2) B₃₂
#   ∂N/∂u     = -1/6
#   ∂V/∂u     = (∂V_num/∂u · N + V_num/6) / N²
#             = (-B₁ + √(π/2) B₃₂) / N  +  V / (6N)

const _SQRTPIO2 = sqrt(π/2)

function vis_ldlin(θ::VisParam, u::VisParam, ρ::AbstractVector)
    N = @. 0.5 - u/6
    ζ = @. π * θ * ρ / MAS2RAD
    B1  = @. ifelse(ζ < 1e-8, one(ζ)/2, besselj1(ζ)/ζ)
    B32 = @. ifelse(ζ < 1e-8, one(ζ)/3, besselj(1.5,ζ)/ζ^1.5)
    # small-ζ limits: J1(ζ)/ζ → 1/2,  J_{3/2}(ζ)/ζ^{3/2} → 1/(3√(π/2)/2)… 
    # More precisely: x^{-ν}J_ν(x) → 1/(2^ν Γ(ν+1)) as x→0
    # J1(x)/x  → 1/2,   J_{3/2}(x)/x^{3/2} → 1/(2^{3/2}Γ(5/2)) = 1/(2√2·3√π/4)
    # At ζ=0: V_num = (1-u)/2 + u√(π/2)/(2^{3/2}Γ(5/2)) = (1-u)/2 + u/3
    #         N·V → (0.5-u/6)·V → finite, V → 1 (normalised)
    # The ifelse guards prevent NaN; the exact limit is 1.0.
    V   = @. ((1-u)*B1 + u*_SQRTPIO2*B32) / N
    return @. ifelse(ζ < 1e-8, one(V), V)
end

# Forward pass variant: returns (V, dV_dθ, dV_du) as element-wise vectors.
# Used by the analytic gradient in gc_polychromatic_fit.jl — avoids calling
# the pullback closure N_uv times.
function vis_ldlin_fwd(θ::Real, u::Real, ρ::AbstractVector)
    N    = 0.5 - u/6
    ζ    = @. π * θ * ρ / MAS2RAD
    safe = ζ .> 1e-8

    B1   = @. ifelse(safe, besselj1(ζ)/ζ,        0.5)
    B32  = @. ifelse(safe, besselj(1.5,ζ)/ζ^1.5,  1/(2^1.5*gamma(2.5)))
    J2   = @. ifelse(safe, besselj(2,ζ),          zero(ζ))
    J52  = @. ifelse(safe, besselj(2.5,ζ),        zero(ζ))

    dB1_dζ  = @. ifelse(safe, -J2  / ζ,           zero(ζ))
    dB32_dζ = @. ifelse(safe, -J52 / ζ^1.5,       zero(ζ))

    V_num   = @. (1-u)*B1 + u*_SQRTPIO2*B32
    V       = @. ifelse(safe, V_num/N, one(ζ))

    # dV/dθ element-wise: chain through ζ = π*θ*ρ/C, dζ/dθ = π*ρ/C
    dV_dζ   = @. ifelse(safe, ((1-u)*dB1_dζ + u*_SQRTPIO2*dB32_dζ) / N, zero(ζ))
    dV_dθ   = @. dV_dζ * π * ρ / MAS2RAD

    # dV/du element-wise (quotient rule)
    dV_du   = @. ifelse(safe, (-B1 + _SQRTPIO2*B32)/N + V/(6*N), zero(ζ))

    return V, dV_dθ, dV_du
end

# In-place variant: writes V, dV_dθ, dV_du directly into pre-allocated vectors.
# Used by channel_chi2_fg! via GCWorkspace to avoid allocation.
function vis_ldlin_fwd!(V::AbstractVector, dV_dθ::AbstractVector, dV_du::AbstractVector,
                         θ::Real, u::Real, ρ::AbstractVector)
    N    = 0.5 - u/6
    @inbounds for j in eachindex(ρ)
        ζ    = π * θ * ρ[j] / MAS2RAD
        safe = ζ > 1e-8
        if safe
            B1   = besselj1(ζ)/ζ
            B32  = besselj(1.5,ζ)/ζ^1.5
            J2   = besselj(2,ζ)
            J52  = besselj(2.5,ζ)
            dB1  = -J2  / ζ
            dB32 = -J52 / ζ^1.5
            V_num     = (1-u)*B1 + u*_SQRTPIO2*B32
            V[j]      = V_num / N
            dVdζ      = ((1-u)*dB1 + u*_SQRTPIO2*dB32) / N
            dV_dθ[j]  = dVdζ * π * ρ[j] / MAS2RAD
            dV_du[j]  = (-B1 + _SQRTPIO2*B32)/N + V[j]/(6*N)
        else
            V[j]     = 1.0
            dV_dθ[j] = 0.0
            dV_du[j] = 0.0
        end
    end
    return nothing
end

function ChainRulesCore.rrule(::typeof(vis_ldlin), θ::VisParam, u::VisParam, ρ::AbstractVector)
    N    = @. 0.5 - u/6
    ζ    = @. π * θ * ρ / MAS2RAD
    safe = ζ .> 1e-8

    B1   = @. ifelse(safe, besselj1(ζ)/ζ,         0.5)
    B32  = @. ifelse(safe, besselj(1.5,ζ)/ζ^1.5,   1/(2^1.5*gamma(2.5)))
    V_num= @. (1-u)*B1 + u*_SQRTPIO2*B32
    V    = @. ifelse(safe, V_num/N, one(ζ))

    # Precompute for gradient
    J2   = @. ifelse(safe, besselj(2,ζ),           zero(ζ))
    J52  = @. ifelse(safe, besselj(2.5,ζ),         zero(ζ))
    dB1_dζ  = @. -J2  / ζ    # from (R1): d[ζ^{-1}J₁]/dζ = -ζ^{-1}J₂
    dB32_dζ = @. -J52 / ζ^1.5

    function vis_ldlin_pullback(ȳ)
        ȳ_safe = @. ifelse(safe, ȳ, zero(ȳ))   # gradient is 0 in the ζ≈0 limit

        # ∂θ
        dV_dζ = @. ((1-u)*dB1_dζ + u*_SQRTPIO2*dB32_dζ) / N
        dζ_dθ = @. π * ρ / MAS2RAD
        ∂θ    = _dparam(θ, ȳ_safe .* dV_dζ .* dζ_dθ)

        # ∂u
        dV_du = @. (-B1 + _SQRTPIO2*B32) / N  +  V / (6*N)
        ∂u    = _dparam(u, ȳ_safe .* dV_du)

        return NoTangent(), ∂θ, ∂u, NoTangent()
    end
    return V, vis_ldlin_pullback
end

# ---------------------------------------------------------------------------
# 3.  Quadratic limb-darkening
# ---------------------------------------------------------------------------
#
# Profile: I(μ) ∝ 1 - u(1-μ) - w(1-μ)²
#
# V_num = (1-u-w) B₁ + (u+2w) √(π/2) B₃₂ - 2w B₂
#   B₂(ζ) = J₂(ζ) / ζ²
# N(u,w) = 0.5 - u/6 - w/12
#
# ∂V/∂θ: same structure via dB₁/dζ, dB₃₂/dζ as above, plus
#   dB₂/dζ = d[ζ^{-2}J₂]/dζ = -ζ^{-2}J₃  (from R1 with ν=2)
#
# ∂V/∂u  = (∂V_num/∂u · N + V_num/6)  / N²
#   ∂V_num/∂u = -B₁ + √(π/2) B₃₂
#
# ∂V/∂w  = (∂V_num/∂w · N + V_num/12) / N²
#   ∂V_num/∂w = -B₁ + 2√(π/2) B₃₂ - 2B₂

function vis_ldquad(θ::VisParam, u::VisParam, w::VisParam, ρ::AbstractVector)
    N   = @. 0.5 - u/6 - w/12
    ζ   = @. π * θ * ρ / MAS2RAD
    B1  = @. ifelse(ζ < 1e-8, 0.5,                  besselj1(ζ)/ζ)
    B32 = @. ifelse(ζ < 1e-8, 1/(2^1.5*gamma(2.5)),  besselj(1.5,ζ)/ζ^1.5)
    B2  = @. ifelse(ζ < 1e-8, 0.125,                 besselj(2,ζ)/ζ^2)
    # J₂(x)/x² → 1/8 as x→0  (x^{-ν}J_ν → 1/(2^ν Γ(ν+1)), ν=2: 1/(4·2)=1/8)
    V_num = @. (1-u-w)*B1 + (u+2w)*_SQRTPIO2*B32 - 2w*B2
    return @. ifelse(ζ < 1e-8, one(ζ), V_num/N)
end

function ChainRulesCore.rrule(::typeof(vis_ldquad), θ::VisParam, u::VisParam, w::VisParam, ρ::AbstractVector)
    N    = @. 0.5 - u/6 - w/12
    ζ    = @. π * θ * ρ / MAS2RAD
    safe = ζ .> 1e-8

    B1  = @. ifelse(safe, besselj1(ζ)/ζ,          0.5)
    B32 = @. ifelse(safe, besselj(1.5,ζ)/ζ^1.5,    1/(2^1.5*gamma(2.5)))
    B2  = @. ifelse(safe, besselj(2,ζ)/ζ^2,         0.125)
    J3  = @. ifelse(safe, besselj(3,ζ),             zero(ζ))
    J52 = @. ifelse(safe, besselj(2.5,ζ),           zero(ζ))
    J2  = @. ifelse(safe, besselj(2,ζ),             zero(ζ))

    dB1_dζ  = @. -J2  / ζ
    dB32_dζ = @. -J52 / ζ^1.5
    dB2_dζ  = @. -J3  / ζ^2

    V_num = @. (1-u-w)*B1 + (u+2w)*_SQRTPIO2*B32 - 2w*B2
    V     = @. ifelse(safe, V_num/N, one(ζ))

    function vis_ldquad_pullback(ȳ)
        ȳs = @. ifelse(safe, ȳ, zero(ȳ))

        dV_dζ = @. ((1-u-w)*dB1_dζ + (u+2w)*_SQRTPIO2*dB32_dζ - 2w*dB2_dζ) / N
        ∂θ = _dparam(θ, ȳs .* dV_dζ .* π .* ρ ./ MAS2RAD)

        dVnum_du = @. -B1 + _SQRTPIO2*B32
        dVnum_dw = @. -B1 + 2*_SQRTPIO2*B32 - 2*B2

        ∂u = _dparam(u, ȳs .* ((dVnum_du .* N .+ V_num ./ 6) ./ N^2))
        ∂w = _dparam(w, ȳs .* ((dVnum_dw .* N .+ V_num ./ 12) ./ N^2))

        return NoTangent(), ∂θ, ∂u, ∂w, NoTangent()
    end
    return V, vis_ldquad_pullback
end

# ---------------------------------------------------------------------------
# 3b. Square-root and Claret four-parameter limb darkening
# ---------------------------------------------------------------------------
#
# Both are sums of half-integer powers of mu, so both are the same shape of thing as
# `vis_ldquad` above and share its derivation. For a profile written as
#
#     I(mu) = sum_j a_j mu^{b_j}
#
# the normalised visibility is
#
#     V = [ sum_j a_j w_j F_j(zeta) ] / [ sum_j a_j w_j ],
#     w_j = 1/(b_j+2),   F_j = Gamma(nu_j+1) 2^{nu_j} J_{nu_j}(zeta)/zeta^{nu_j},  nu_j = b_j/2+1
#
# with F_j -> 1 as zeta -> 0. Substituting quadratic limb darkening reproduces `vis_ldquad`'s
# coefficients and its N exactly, and both laws below were checked against a direct numerical
# integration of their own I(mu) to 1.2e-10.
#
# The Bessel ORDERS here are fixed numbers, unlike `vis_ldpow` where nu = alpha/2+1 moves with a
# fitted parameter. That is why neither of these needs the ForwardDiff specialisation ldpow
# carries: `besselj` only sees a Dual in its argument, never in its order. Both are also linear
# in their coefficients, so those derivatives are just the basis functions.

# w_j * Gamma(nu_j+1) * 2^{nu_j} for b = 1/2 and b = 3/2 (b = 0, 1, 2 give 1, sqrt(pi/2), 2).
const _LDF12 = gamma(2.25) * 2^1.25 / 2.5      # 1.07790027…
const _LDF32 = gamma(2.75) * 2^1.75 / 3.5      # 1.54567277…

# x^{-nu} J_nu(x) -> 1/(2^nu Gamma(nu+1)) as x -> 0
_b0(ν) = 1 / (2^ν * gamma(ν + 1))

"""
    vis_ldsqrt(θ, c, d, ρ)

Square-root limb darkening, `I(μ) = 1 - c(1-μ) - d(1-√μ)`, for a disc of diameter `θ` in mas.

`c` is the linear coefficient and `d` the square-root one, in the order Claret tabulates them.
"""
function vis_ldsqrt(θ::VisParam, c::VisParam, d::VisParam, ρ::AbstractVector)
    N   = @. 0.5 - c/6 - d/10
    ζ   = @. π * θ * ρ / MAS2RAD
    B1  = @. ifelse(ζ < 1e-8, _b0(1.0),    besselj1(ζ)/ζ)
    B32 = @. ifelse(ζ < 1e-8, _b0(1.5),    besselj(1.5, ζ)/ζ^1.5)
    B54 = @. ifelse(ζ < 1e-8, _b0(1.25),   besselj(1.25, ζ)/ζ^1.25)
    V_num = @. (1-c-d)*B1 + c*_SQRTPIO2*B32 + d*_LDF12*B54
    return @. ifelse(ζ < 1e-8, one(ζ), V_num/N)
end

function ChainRulesCore.rrule(::typeof(vis_ldsqrt), θ::VisParam, c::VisParam, d::VisParam,
                              ρ::AbstractVector)
    N    = @. 0.5 - c/6 - d/10
    ζ    = @. π * θ * ρ / MAS2RAD
    safe = ζ .> 1e-8

    B1  = @. ifelse(safe, besselj1(ζ)/ζ,           _b0(1.0))
    B32 = @. ifelse(safe, besselj(1.5, ζ)/ζ^1.5,   _b0(1.5))
    B54 = @. ifelse(safe, besselj(1.25, ζ)/ζ^1.25, _b0(1.25))
    # d/dζ [ζ^{-ν} J_ν] = -ζ^{-ν} J_{ν+1}
    dB1  = @. ifelse(safe, -besselj(2.0,  ζ)/ζ,       zero(ζ))
    dB32 = @. ifelse(safe, -besselj(2.5,  ζ)/ζ^1.5,   zero(ζ))
    dB54 = @. ifelse(safe, -besselj(2.25, ζ)/ζ^1.25,  zero(ζ))

    V_num = @. (1-c-d)*B1 + c*_SQRTPIO2*B32 + d*_LDF12*B54
    V     = @. ifelse(safe, V_num/N, one(ζ))

    function vis_ldsqrt_pullback(ȳ)
        ȳs = @. ifelse(safe, ȳ, zero(ȳ))
        dV_dζ = @. ((1-c-d)*dB1 + c*_SQRTPIO2*dB32 + d*_LDF12*dB54) / N
        ∂θ = _dparam(θ, ȳs .* dV_dζ .* π .* ρ ./ MAS2RAD)
        dVnum_dc = @. -B1 + _SQRTPIO2*B32
        dVnum_dd = @. -B1 + _LDF12*B54
        ∂c = _dparam(c, ȳs .* ((dVnum_dc .* N .+ V_num ./ 6)  ./ N^2))
        ∂d = _dparam(d, ȳs .* ((dVnum_dd .* N .+ V_num ./ 10) ./ N^2))
        return NoTangent(), ∂θ, ∂c, ∂d, NoTangent()
    end
    return V, vis_ldsqrt_pullback
end

"""
    vis_ldclaret4(θ, c1, c2, c3, c4, ρ)

Claret's four-parameter limb darkening, `I(μ) = 1 - Σₖ cₖ (1 - μ^{k/2})`, for a disc of
diameter `θ` in mas.

The law most precise work uses, and the one Claret's tables are usually quoted in.
"""
function vis_ldclaret4(θ::VisParam, c1::VisParam, c2::VisParam, c3::VisParam, c4::VisParam,
                       ρ::AbstractVector)
    N    = @. 0.5 - c1/10 - c2/6 - 3*c3/14 - c4/4
    ζ    = @. π * θ * ρ / MAS2RAD
    B1   = @. ifelse(ζ < 1e-8, _b0(1.0),  besselj1(ζ)/ζ)
    B54  = @. ifelse(ζ < 1e-8, _b0(1.25), besselj(1.25, ζ)/ζ^1.25)
    B32  = @. ifelse(ζ < 1e-8, _b0(1.5),  besselj(1.5,  ζ)/ζ^1.5)
    B74  = @. ifelse(ζ < 1e-8, _b0(1.75), besselj(1.75, ζ)/ζ^1.75)
    B2   = @. ifelse(ζ < 1e-8, _b0(2.0),  besselj(2.0,  ζ)/ζ^2)
    V_num = @. (1-c1-c2-c3-c4)*B1 + c1*_LDF12*B54 + c2*_SQRTPIO2*B32 +
               c3*_LDF32*B74 + c4*2*B2
    return @. ifelse(ζ < 1e-8, one(ζ), V_num/N)
end

function ChainRulesCore.rrule(::typeof(vis_ldclaret4), θ::VisParam, c1::VisParam, c2::VisParam,
                              c3::VisParam, c4::VisParam, ρ::AbstractVector)
    N    = @. 0.5 - c1/10 - c2/6 - 3*c3/14 - c4/4
    ζ    = @. π * θ * ρ / MAS2RAD
    safe = ζ .> 1e-8

    B1  = @. ifelse(safe, besselj1(ζ)/ζ,           _b0(1.0))
    B54 = @. ifelse(safe, besselj(1.25, ζ)/ζ^1.25, _b0(1.25))
    B32 = @. ifelse(safe, besselj(1.5,  ζ)/ζ^1.5,  _b0(1.5))
    B74 = @. ifelse(safe, besselj(1.75, ζ)/ζ^1.75, _b0(1.75))
    B2  = @. ifelse(safe, besselj(2.0,  ζ)/ζ^2,    _b0(2.0))

    dB1  = @. ifelse(safe, -besselj(2.0,  ζ)/ζ,       zero(ζ))
    dB54 = @. ifelse(safe, -besselj(2.25, ζ)/ζ^1.25,  zero(ζ))
    dB32 = @. ifelse(safe, -besselj(2.5,  ζ)/ζ^1.5,   zero(ζ))
    dB74 = @. ifelse(safe, -besselj(2.75, ζ)/ζ^1.75,  zero(ζ))
    dB2  = @. ifelse(safe, -besselj(3.0,  ζ)/ζ^2,     zero(ζ))

    V_num = @. (1-c1-c2-c3-c4)*B1 + c1*_LDF12*B54 + c2*_SQRTPIO2*B32 +
               c3*_LDF32*B74 + c4*2*B2
    V     = @. ifelse(safe, V_num/N, one(ζ))

    function vis_ldclaret4_pullback(ȳ)
        ȳs = @. ifelse(safe, ȳ, zero(ȳ))
        dV_dζ = @. ((1-c1-c2-c3-c4)*dB1 + c1*_LDF12*dB54 + c2*_SQRTPIO2*dB32 +
                    c3*_LDF32*dB74 + c4*2*dB2) / N
        ∂θ = _dparam(θ, ȳs .* dV_dζ .* π .* ρ ./ MAS2RAD)
        # linear in every coefficient; dN/dc_k is the constant beside it in N
        d1 = @. -B1 + _LDF12*B54
        d2 = @. -B1 + _SQRTPIO2*B32
        d3 = @. -B1 + _LDF32*B74
        d4 = @. -B1 + 2*B2
        ∂c1 = _dparam(c1, ȳs .* ((d1 .* N .+ V_num ./ 10)      ./ N^2))
        ∂c2 = _dparam(c2, ȳs .* ((d2 .* N .+ V_num ./ 6)       ./ N^2))
        ∂c3 = _dparam(c3, ȳs .* ((d3 .* N .+ V_num .* 3 ./ 14) ./ N^2))
        ∂c4 = _dparam(c4, ȳs .* ((d4 .* N .+ V_num ./ 4)       ./ N^2))
        return NoTangent(), ∂θ, ∂c1, ∂c2, ∂c3, ∂c4, NoTangent()
    end
    return V, vis_ldclaret4_pullback
end

# ---------------------------------------------------------------------------
# 4.  Power-law limb-darkening
# ---------------------------------------------------------------------------
#
# Profile: I(μ) ∝ μ^α
#
# V(θ, α, ρ) = Γ(ν+1) · J_ν(ζ) · (ζ/2)^{-ν},   ν = α/2 + 1
#            = Γ(ν+1) · [ζ^{-ν} J_ν(ζ)] · 2^ν
#
# ∂V/∂θ  (through ζ, using R1):
#   d/dζ [ζ^{-ν}J_ν] = -ζ^{-ν}J_{ν+1}
#   ∂V/∂θ = Γ(ν+1)·2^ν · (-J_{ν+1}(ζ)/ζ^ν) · (πρ/C)
#
# ∂V/∂α  (through ν = α/2+1, dν/dα = 1/2):
#   ∂V/∂ν = ψ(ν+1)·V
#           + Γ(ν+1)·2^ν·[ ∂/∂ν(J_ν/ζ^ν) ]
#   where ∂/∂ν(J_ν/ζ^ν) = ∂J_ν/∂ν / ζ^ν - log(ζ)·J_ν/ζ^ν
#
# ∂J_ν(x)/∂ν is computed via the series (see _dbesselj_dnu below).

# Truncated series for ∂J_ν(x)/∂ν:
#
#   J_ν(x) = Σ_{k≥0} (-1)^k (x/2)^{ν+2k} / [k! Γ(ν+k+1)]
#
#   ∂J_ν/∂ν = log(x/2)·J_ν(x)
#             - Σ_{k≥0} (-1)^k (x/2)^{ν+2k} ψ(ν+k+1) / [k! Γ(ν+k+1)]
#
# The series is only usable for small x, and it fails *silently and catastrophically*
# rather than gracefully. Its terms alternate in sign and peak near k ≈ x/2 at a
# magnitude of order (x/2)^x, so the answer is the residue of a cancellation that costs
# more digits than Float64 has. Measured absolute error against the true derivative,
# with 30 terms:
#
#     x       1      10      15      20        30        60       120
#     err  1e-13   2e-12   3e-11   1.4e-6     1.6e+5    1e+24     4e+42
#
# and adding terms does not help — at x = 30 it is 6e-5 with 60 terms or with 500,
# because the limit is cancellation, not truncation. The original comment here assumed
# ζ = πθρ/C stays below 10 ("diameters ≲ 5 mas and baselines ≲ 100 m/λ"), which holds
# for the α Cen PIONIER data (ζ ≤ 11) but not in general: 8.5 mas on a 330 m CHARA
# baseline in H reaches ζ ≈ 26, and a 43 mas supergiant on the VLTI reaches ζ ≈ 39.
# There the series returns garbage, the power-law gradient follows it, and the fit
# fails without any warning.
#
# Above x = 15 we therefore differentiate the library Bessel in ν directly, with a
# 5-point central stencil. That is uniformly accurate (≤ 2e-11 for x from 1 to 300)
# but ~2.6x slower than the series, so the series keeps the common small-x case where
# it is both faster and slightly more accurate.
const _DBESSEL_SERIES_XMAX = 15.0

function _dbesselj_dnu(ν::Real, x::Real; nterms::Int=30)
    x == 0.0 && return zero(x)

    if x > _DBESSEL_SERIES_XMAX
        h = 1e-3
        return (besselj(ν-2h, x) - 8*besselj(ν-h, x) +
                8*besselj(ν+h, x) - besselj(ν+2h, x)) / (12h)
    end

    lhalf = log(x/2)
    half_x = x/2
    # accumulate the digamma series
    term0_coeff = half_x^ν / gamma(ν+1)
    s = digamma(ν+1) * term0_coeff    # k=0 term of the ψ-sum
    c = term0_coeff                    # running (x/2)^{ν+2k}/(k!Γ(ν+k+1))
    for k in 1:nterms
        c *= -(half_x^2) / (k * (ν+k))   # ratio between successive terms
        # Note: Γ(ν+k+1)/Γ(ν+k) = ν+k, so c_{k}/c_{k-1} = -(x/2)²/(k(ν+k))
        # (the 1/(ν+k) comes from ν+k+1-1 = ν+k in the ratio of Gamma functions)
        s += digamma(ν+k+1) * c
    end
    return lhalf * besselj(ν, x) - s
end

function vis_ldpow(θ::VisParam, α::VisParam, ρ::AbstractVector)
    ν   = @. α/2 + 1
    ζ   = @. π * θ * ρ / MAS2RAD
    Vz  = @. ifelse(ζ < 1e-8,
                    one(ζ),                              # limit: V → 1 as ζ→0
                    gamma(ν+1) * besselj(ν,ζ) * (ζ/2)^(-ν))
    return Vz
end

function ChainRulesCore.rrule(::typeof(vis_ldpow), θ::VisParam, α::VisParam, ρ::AbstractVector)
    ν    = @. α/2 + 1
    ζ    = @. π * θ * ρ / MAS2RAD
    safe = ζ .> 1e-8

    Jν   = @. ifelse(safe, besselj(ν, ζ),   zero(ζ))
    Jν1  = @. ifelse(safe, besselj(ν+1, ζ), zero(ζ))
    pow  = @. ifelse(safe, (ζ/2)^(-ν),      one(ζ))
    Gν1  = @. gamma(ν+1)

    V = @. ifelse(safe, Gν1 * Jν * pow, one(ζ))

    # ∂V/∂ζ = Γ(ν+1)·2^ν · d/dζ[ζ^{-ν}J_ν] = Γ(ν+1)·2^ν·(-J_{ν+1}/ζ^ν)
    #        = -Γ(ν+1) · J_{ν+1} · (ζ/2)^{-ν} / (ζ/2)   [absorbing 2^ν factor]
    # More directly: ∂V/∂ζ = Γ(ν+1)·(-J_{ν+1}(ζ)/ζ^ν)·2^ν
    #                       = -Gν1 · Jν1 · 2^ν / ζ^ν
    ζ_safe = @. ifelse(safe, ζ, one(ζ))   # avoid 0^ν below
    dV_dζ  = @. ifelse(safe, -Gν1 * Jν1 * 2^ν / ζ_safe^ν, zero(ζ))

    # ∂V/∂α: compute ∂J_ν/∂ν via series at each baseline (scalar loop is fine
    # since nuv is small compared to cost of Bessel evaluations)
    dJν_dν = map((z,s) -> s ? _dbesselj_dnu(ν, z) : 0.0, ζ, safe)

    # ∂/∂ν [Γ(ν+1)] = ψ(ν+1)·Γ(ν+1)
    # ∂V/∂ν = ψ(ν+1)·V + Γ(ν+1)·(ζ/2)^{-ν}·[dJν_dν - log(ζ/2)·Jν]
    #
    # The log is log(ζ/2), NOT log(ζ): it comes from ∂/∂ν (ζ/2)^{-ν} =
    # -log(ζ/2)·(ζ/2)^{-ν}. Using log(ζ) makes ∂V/∂α too small by exactly
    # log(2)·V/2, which is enough to stop a gradient-based fit ~0.015 short in α
    # while every derivative-free method converges correctly. See the gradient
    # test in test/test_model_gradients.jl.
    ψν1 = @. digamma(ν+1)
    dV_dν = @. ifelse(safe,
        ψν1*V + Gν1 * 2^ν * (dJν_dν - log(ζ_safe/2)*Jν) / ζ_safe^ν,
        zero(ζ))
    # dν/dα = 1/2
    dV_dα = dV_dν ./ 2

    function vis_ldpow_pullback(ȳ)
        ȳs = @. ifelse(safe, ȳ, zero(ȳ))
        ∂θ = _dparam(θ, ȳs .* dV_dζ .* π .* ρ ./ MAS2RAD)
        ∂α = _dparam(α, ȳs .* dV_dα)
        return NoTangent(), ∂θ, ∂α, NoTangent()
    end
    return V, vis_ldpow_pullback
end

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
#
# These match the original signatures exactly and are safe for ForwardDiff /
# Zygote because they use ifelse instead of mutating findall/setindex!.

