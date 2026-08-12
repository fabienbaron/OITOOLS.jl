# hankel.jl
#
# Numerical Hankel transform for radial profile → visibility conversion.
#
# All work in dimensionless / mas units:
#   r : radius in mas
#   B : spatial frequency in cycles/mas (= baseline/wavelength in SI,
#       converted by the caller)
#
# The Hankel transform of order n is defined as:
#   Hₙ(B) = ∫₀^∞ I(r) Jₙ(2π B r) r dr
#
# Normalized visibility:
#   V(B) = Hₙ(B) / N,   N = H₀(0) = ∫₀^∞ I(r) r dr

using SpecialFunctions      # besselj, besselj0, besselj1
using LinearAlgebra
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)


##############################################################################
# SECTION 1: Core quadrature
##############################################################################

"""
    trapz(y, x) -> Float64

Composite trapezoidal quadrature of `y` over `x`.
"""
function trapz(y::AbstractVector, x::AbstractVector)
    @assert length(y) == length(x) "trapz: y and x must have the same length"
    s = zero(promote_type(eltype(y), eltype(x)))
    @inbounds for i in 1:(length(x)-1)
        s += (y[i] + y[i+1]) * (x[i+1] - x[i])
    end
    return s / 2
end

"""
    trapz_weights(x) -> Vector{Float64}

Return the trapezoid quadrature weights for nodes `x`, so that
`sum(w .* y) ≈ trapz(y, x)` for any `y`.
"""
function trapz_weights(x::AbstractVector)
    n = length(x)
    return [(i == 1)   ? (x[2] - x[1]) / 2 :
            (i == n)   ? (x[n] - x[n-1]) / 2 :
                         (x[i+1] - x[i-1]) / 2
            for i in 1:n]
end


##############################################################################
# SECTION 2: Hankel kernel
##############################################################################

"""
    hankel_transform(I, r, B; n=0) -> Vector

Compute the nth-order Hankel transform:
    Hₙ(B[j]) = ∫ I(r) Jₙ(2π B[j] r) r dr   for each j

Uses composite trapezoid rule on the supplied `r` grid.

Arguments
---------
I  : radial profile values at grid points r  (length Nr)
r  : radius grid in mas, must be non-negative and increasing  (length Nr)
B  : spatial frequency grid in cycles/mas  (length nB)
n  : Bessel order (integer ≥ 0, default 0)

Returns
-------
H  : Vector of length nB
"""
function hankel_transform(I::AbstractVector, r::AbstractVector,
                           B::AbstractVector; n::Int=0)
    Nr = length(r)
    nB = length(B)
    @assert length(I) == Nr "hankel_transform: I and r must have the same length"

    w  = trapz_weights(r)          # trapezoid weights, length Nr

    # Functional style: map over B to avoid mutating an output array
    # (needed for Enzyme compatibility)
    function _ht_single(Bj)
        T   = promote_type(eltype(I), typeof(Bj))
        acc = zero(T)
        # T(2π), not 2π: the latter is Int*Irrational => Float64, which would promote the
        # whole transform (and every Hankel component downstream) out of Float32.
        twoπ = T(2π)
        @inbounds for k in 1:Nr
            arg = twoπ * Bj * r[k]
            Jn  = (n == 0) ? besselj0(arg) :
                  (n == 1) ? besselj1(arg) :
                             besselj(n, arg)
            acc += w[k] * I[k] * Jn * r[k]
        end
        return acc
    end

    return map(_ht_single, B)
end

"""
    hankel_norm(I, r) -> T

Compute the normalisation constant N = ∫ I(r) r dr = H₀(B=0).
"""
hankel_norm(I::AbstractVector, r::AbstractVector) = trapz(I .* r, r)

"""
    hankel_vis(I, r, B; n=0) -> Vector

Normalised Hankel visibility Hₙ(B) / N on the supplied baseline grid.
"""
function hankel_vis(I::AbstractVector, r::AbstractVector,
                    B::AbstractVector; n::Int=0)
    H = hankel_transform(I, r, B; n)
    N = hankel_norm(I, r)
    return H ./ N
end


# ChainRules rrules — must be included after hankel_norm/transform/vis are defined
include("hankel_chainrules.jl")

##############################################################################
# SECTION 3: Coarse-grid / interpolation pathway
##############################################################################

"""
    hankel_vis_interp(I, r, B_query; n=0, nB_coarse=40) -> Vector

Compute the normalised Hankel visibility at arbitrary query baselines
`B_query` by:
  1. Evaluating the transform on a coarse uniform grid of `nB_coarse` points
  2. Linearly interpolating to `B_query`

Faster than a direct call for large |B_query| with many repeats.
"""
function hankel_vis_interp(I::AbstractVector, r::AbstractVector,
                            B_query::AbstractVector;
                            n::Int=0, nB_coarse::Int=40)
    B_min = minimum(B_query)
    B_max = maximum(B_query)
    B_coarse = range(B_min, B_max; length=nB_coarse) |> collect
    H_coarse = hankel_vis(I, r, B_coarse; n)
    # Linear interpolation (extrapolates constantly beyond range)
    return [_interp1(B_coarse, H_coarse, b) for b in B_query]
end

function _interp1(x::AbstractVector, y::AbstractVector, xi::Real)
    xi <= x[1]   && return y[1]
    xi >= x[end] && return y[end]
    k = searchsortedfirst(x, xi) - 1
    t = (xi - x[k]) / (x[k+1] - x[k])
    return y[k] + t * (y[k+1] - y[k])
end


##############################################################################
# SECTION 4: Full complex visibility with azimuthal variations
##############################################################################

"""
    hankel_vis_full(I, r, B, PA; n_orders=Int[], amps=Float64[], phis=Float64[],
                    nB_coarse=40) -> Vector{Complex{T}}

Full complex visibility including azimuthal variation terms.

Arguments
---------
I        : radial profile on r grid
r        : radius grid (mas)
B        : baseline magnitudes (cycles/mas)
PA       : position angles (radians), same length as B
n_orders : azimuthal orders (integer list)
amps     : amplitudes of each azimuthal term (same length as n_orders)
phis     : phase offsets in radians for each term (same length as n_orders)

The element type follows `I`, `r` and `B`: Float32 inputs give a `Vector{ComplexF32}`.

`nB_coarse` is accepted for signature compatibility with `hankel_vis_interp` but is not
used here — this routine evaluates the transform directly on `B`.
"""
function hankel_vis_full(I::AbstractVector, r::AbstractVector,
                          B::AbstractVector, PA::AbstractVector;
                          n_orders::AbstractVector{<:Integer} = Int[],
                          amps::AbstractVector{<:Real}       = Float64[],
                          phis::AbstractVector{<:Real}       = Float64[],
                          nB_coarse::Int                     = 40)
    N    = hankel_norm(I, r)
    V    = complex.(hankel_transform(I, r, B; n=0) ./ N)   # H₀/N, real→complex

    for (ni, ai, φi) in zip(n_orders, amps, phis)
        Hn = hankel_transform(I, r, B; n=ni) ./ N
        V .+= ai .* (-1im)^ni .* Hn .* cos.(ni .* (PA .+ φi))
    end
    return V
end


##############################################################################
# SECTION 5: Profile expression evaluator
##############################################################################

"""
    eval_profile(expr_str, params, r) -> Vector{Float64}

Evaluate a radial profile string expression on radius grid `r`.

`\$R` in the expression is replaced by the radius array.
Other `\$name` references are replaced by scalar values from `params`.

Example:
    eval_profile("exp(-(\$R / \$scale)^2)", Dict("scale"=>2.0), r)
"""
function eval_profile(expr_str::AbstractString,
                       params::Dict{String},
                       r::AbstractVector)::Vector{Float64}
    # Replace $R with the Julia variable name _R_ (reserved)
    s = replace(expr_str, "\$R"  => "_R_")
    s = replace(s,        "\$MU" => "_MU_")   # cos(zenith angle), useful for LD

    # Replace $param_name with numeric literal
    for (k, v) in params
        s = replace(s, "\$$k" => string(Float64(v)))
    end

    # Build a function (_R_, _MU_) -> profile and call it
    # _MU_ = sqrt(1 - (r/(r_max))^2) for limb-darkening profiles
    r_max = maximum(r)
    _MU_  = @. sqrt(max(0.0, 1.0 - (r / r_max)^2))

    # Auto-broadcast: add dots before operators so expressions work on vectors
    s = replace(s, r"(?<![.])(\^)" => ".^")
    s = replace(s, r"(?<![.])(\*)" => ".*")
    s = replace(s, r"(?<![.])(\/)" => "./")

    fn_expr = Meta.parse("(_R_, _MU_) -> @. $s")
    fn      = eval(fn_expr)
    result  = Base.invokelatest(fn, r, _MU_)

    # result may be scalar if the expression doesn't reference _R_
    if result isa Number
        return fill(Float64(result), length(r))
    else
        return Float64.(result)
    end
end


##############################################################################
# SECTION 5b: Compiled profile closure (AD-transparent)
##############################################################################

"""
    compile_profile(expr_str, param_names) -> RuntimeGeneratedFunction

Compile a profile expression string into an AD-transparent closure.

The returned function has signature `(r, mu, param1, param2, ...) -> Vector`,
where each `param_i` corresponds to the names in `param_names` (in order).

`\$R` and `\$MU` are replaced by the first two arguments.
Other `\$name` references are replaced by the corresponding scalar argument.

Example:
    fn = compile_profile("exp(-(\$R / \$scale)^2 / 2)", ["scale"])
    I  = fn(r, mu, 2.0)   # equivalent to @. exp(-(r / 2.0)^2 / 2)
"""
function compile_profile(expr_str::AbstractString,
                          param_names::Vector{String})
    # Replace $R -> _R_, $MU -> _MU_, $param -> param_ident
    s = replace(expr_str, "\$R"  => "_R_")
    s = replace(s,        "\$MU" => "_MU_")
    for name in param_names
        ident = replace(name, "," => "__")   # mangle to valid Julia identifier
        s = replace(s, "\$$name" => ident)
    end

    # Build argument list: (_R_, _MU_, param1, param2, ...)
    arg_symbols = [Symbol(replace(n, "," => "__")) for n in param_names]
    args = Expr(:tuple, :_R_, :_MU_, arg_symbols...)

    fn_expr = Meta.parse("$args -> @. $s")
    return @RuntimeGeneratedFunction(fn_expr)
end


##############################################################################
# SECTION 6: Analytic ground truths (useful for testing and validation)
##############################################################################

"""Analytic normalised visibility for a uniform disk of diameter D (mas)."""
function analytic_ud_vis(D::Real, B::AbstractVector)
    x = @. π * D * B
    return @. ifelse(x < 1e-8, one(x), 2 * besselj1(x) / x)
end

"""Analytic ∂V/∂D for a uniform disk (via Bessel recurrence)."""
function analytic_ud_dvis_dD(D::Real, B::AbstractVector)
    x   = @. π * D * B
    # d/dD [2J₁(πDB)/(πDB)] = π B · d/dx[2J₁(x)/x]
    # d/dx[2J₁(x)/x] = 2[xJ₀(x) - 2J₁(x)] / x²   (standard recurrence)
    dVdx = @. ifelse(x < 1e-8, zero(x),
                     2*(x*besselj0(x) - 2*besselj1(x)) / x^2)
    return @. dVdx * π * B
end

"""Analytic normalised visibility for a Gaussian profile I(r)=exp(-r²/(2σ²))."""
function analytic_gauss_vis(σ::Real, B::AbstractVector)
    return @. exp(-2π^2 * σ^2 * B^2)
end

"""Analytic ∂V/∂σ for a Gaussian."""
function analytic_gauss_dvis_dσ(σ::Real, B::AbstractVector)
    return @. -4π^2 * σ * B^2 * exp(-2π^2 * σ^2 * B^2)
end

"""Analytic normalised visibility for a thin ring of radius R (mas)."""
analytic_ring_vis(R::Real, B::AbstractVector) = @. besselj0(2π * B * R)

"""Analytic ∂V/∂R for a thin ring."""
analytic_ring_dvis_dR(R::Real, B::AbstractVector) = @. -2π * B * besselj1(2π * B * R)


##############################################################################
# SECTION 7: Finite-difference gradient helper
##############################################################################

"""
    fd_gradient(f, p; h=1e-6) -> Vector

Central finite-difference gradient of scalar→scalar `f` at vector `p`.
"""
function fd_gradient(f, p::AbstractVector; h::Float64=1e-6)
    g = similar(p)
    for i in eachindex(p)
        pp = copy(p); pp[i] += h
        pm = copy(p); pm[i] -= h
        g[i] = (f(pp) - f(pm)) / (2h)
    end
    return g
end

"""
    fd_jacobian(f, p; h=1e-6) -> Matrix

Central finite-difference Jacobian of vector→vector `f` at vector `p`.
Returns matrix of shape (length(f(p)), length(p)).
"""
function fd_jacobian(f, p::AbstractVector; h::Float64=1e-6)
    f0 = f(p)
    J  = zeros(length(f0), length(p))
    for i in eachindex(p)
        pp = copy(p); pp[i] += h
        pm = copy(p); pm[i] -= h
        J[:, i] = (f(pp) - f(pm)) / (2h)
    end
    return J
end
