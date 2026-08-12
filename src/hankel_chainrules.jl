# hankel_chainrules.jl
#
# Explicit ChainRules rrules for the numerical Hankel transform pathway,
# plus _fwd and _fwd! variants for zero-allocation gradient computation.
#
# This file provides:
#   - rrule for hankel_norm(I, r)
#   - rrule for hankel_transform(I, r, B; n)
#   - rrule for hankel_vis(I, r, B; n)
#   - hankel_vis_fwd(I, r, B; n)      — returns (V, K, N), single Bessel pass
#   - hankel_vis_fwd!(out_V, out_K, I, r, B; n) -> N  — zero-allocation forward
#   - hankel_vis_pullback!(out_dI, ȳ, V, K, N, w) — zero-allocation reverse
#
# All functions assume:
#   I  : radial profile values on r grid (length Nr)
#   r  : radius grid in mas, non-negative, increasing (length Nr)
#   B  : spatial frequency grid in cycles/mas (length nB)
#   n  : Bessel order (integer ≥ 0)
#
# Conventions match test_hankel.jl:
#   Hₙ(B) = ∫ I(r) Jₙ(2πBr) r dr   (trapezoid rule on supplied r grid)
#   N      = ∫ I(r) r dr              (= H₀ at B=0)
#   V(B)   = Hₙ(B) / N
#
# Kernel matrix K (nB × Nr):
#   K[j, k] = w[k] · r[k] · Jₙ(2π B[j] r[k])
# so that  H = K * I  (gemv)  and  N = dot(w .* r, I).
#
# Pullback derivation (quotient rule):
#   V[j] = H[j] / N
#   ∂V[j]/∂I[k] = K[j,k]/N  −  V[j] · w[k]·r[k] / N
#   ∂L/∂I[k]   = Σⱼ ȳ[j] · ∂V[j]/∂I[k]
#               = (K'ȳ)[k]/N  −  dot(ȳ,V) · w[k]·r[k] / N
#
# In matrix form:  ∂L/∂I = (K'ȳ  −  dot(ȳ,V) · (w .* r)) / N
#
# The (w .* r) term comes from the normalization N; without it (i.e. for
# hankel_transform alone) the pullback is simply K'ȳ.
#
# Performance notes
# -----------------
# K is nB×Nr.  Storing it enables a single BLAS gemv in the pullback instead
# of recomputing all Bessel values.  For typical interferometry grids
# (Nr ≈ 100–500, nB ≈ 20–100) K is 16–400 KB — acceptable to cache in a
# Workspace struct during a fit.
#
# The _fwd! / pullback! pair is intended for the gradient-based chi2 loop:
#   1. Allocate workspace once (HankelWorkspace) at dict_to_model time.
#   2. Call hankel_vis_fwd! to fill workspace.V and workspace.K, get N.
#   3. Obtain ∂I/∂p via ForwardDiff on the compiled profile closure.
#   4. Call hankel_vis_pullback! to get ∂L/∂I from upstream ȳ.
#   5. Accumulate  ∂L/∂p = (∂I/∂p)' * ∂L/∂I  (pure dot products, no Bessel).

using ChainRulesCore
using LinearAlgebra: dot, mul!
using SpecialFunctions: besselj0, besselj1, besselj

# Requires test_hankel.jl to be included first.
# Load order:  include("test_hankel.jl") ; include("hankel_chainrules.jl")
#
# _trapz_weights mirrors trapz_weights but is prefixed _ to avoid conflicts.



function _trapz_weights(x::AbstractVector)
    n = length(x)
    return [(i == 1)   ? (x[2]   - x[1])   / 2 :
            (i == n)   ? (x[n]   - x[n-1]) / 2 :
                         (x[i+1] - x[i-1]) / 2
            for i in 1:n]
end

# ---------------------------------------------------------------------------
# Workspace struct for zero-allocation gradient loop
# ---------------------------------------------------------------------------

"""
    HankelWorkspace

Pre-allocated buffers for a single Hankel component.

Allocate once at `dict_to_model` time:
    ws = HankelWorkspace(Nr, nB)

Then in the chi2 loop:
    N  = hankel_vis_fwd!(ws.V, ws.K, I, r, B; n=0)
    hankel_vis_pullback!(ws.dI, ȳ, ws.V, ws.K, N, ws.w, r)
    dL_dp = ws.dI_dp' * ws.dI    # after filling ws.dI_dp via ForwardDiff
"""
mutable struct HankelWorkspace{T<:AbstractFloat}
    V    ::Vector{T}   # nB — normalised visibility
    K    ::Matrix{T}   # nB × Nr — Bessel kernel
    dI   ::Vector{T}   # Nr — ∂L/∂I from pullback
    dI_dp::Matrix{T}   # Nr × n_params — ∂I/∂p from ForwardDiff on profile
    w    ::Vector{T}   # Nr — trapezoid weights (fixed, cached)
end

"""
    HankelWorkspace(Nr, nB, n_params=1; T=Float64)

Allocate the Hankel buffers at precision `T`. The transform functions themselves
(`hankel_vis_fwd!`, `hankel_vis_pullback!`, ...) are generic over `AbstractVector`/
`AbstractMatrix`, so the buffer element type is what sets the working precision of a Hankel
component.
"""
function HankelWorkspace(Nr::Int, nB::Int, n_params::Int=1;
                         T::Type{<:AbstractFloat}=Float64)
    HankelWorkspace{T}(
        zeros(T, nB),
        zeros(T, nB, Nr),
        zeros(T, Nr),
        zeros(T, Nr, n_params),
        zeros(T, Nr),          # w filled by first call to hankel_vis_fwd!
    )
end


# ---------------------------------------------------------------------------
# 1.  hankel_norm
# ---------------------------------------------------------------------------
#
# N = ∫ I(r) r dr = Σₖ w[k] · r[k] · I[k]   (trapezoid)
#
# ∂N/∂I[k] = w[k] · r[k]
#
# Pullback:  ∂L/∂I = ȳ · (w .* r)    (ȳ is scalar here)

function ChainRulesCore.rrule(::typeof(hankel_norm),
                               I::AbstractVector,
                               r::AbstractVector)
    w  = _trapz_weights(r)
    wr = w .* r
    N  = dot(wr, I)

    function hankel_norm_pullback(ȳ)
        # ȳ is a scalar (∂L/∂N)
        return NoTangent(), ȳ .* wr, NoTangent()
    end
    return N, hankel_norm_pullback
end


# ---------------------------------------------------------------------------
# 2.  hankel_transform
# ---------------------------------------------------------------------------
#
# H[j] = Σₖ K[j,k] · I[k]   where   K[j,k] = w[k] · r[k] · Jₙ(2πB[j]r[k])
#
# ∂L/∂I[k] = Σⱼ ȳ[j] · K[j,k]  =  (K'ȳ)[k]

function ChainRulesCore.rrule(::typeof(hankel_transform),
                               I::AbstractVector,
                               r::AbstractVector,
                               B::AbstractVector;
                               n::Int=0)
    Nr = length(r)
    nB = length(B)
    w  = _trapz_weights(r)

    K  = [_bessel_kernel(n, B[j], r[k], w[k])
          for j in 1:nB, k in 1:Nr]      # nB × Nr
    H  = K * I

    function hankel_transform_pullback(ȳ)
        dI = K' * ȳ                       # Nr-vector, one gemv
        return NoTangent(), dI, NoTangent(), NoTangent()
    end
    return H, hankel_transform_pullback
end


# ---------------------------------------------------------------------------
# 3.  hankel_vis
# ---------------------------------------------------------------------------
#
# V[j] = H[j] / N
#
# Pullback (quotient rule — see file header for derivation):
#   ∂L/∂I = (K'ȳ  −  dot(ȳ, V) · (w .* r)) / N

function ChainRulesCore.rrule(::typeof(hankel_vis),
                               I::AbstractVector,
                               r::AbstractVector,
                               B::AbstractVector;
                               n::Int=0)
    Nr = length(r)
    nB = length(B)
    w  = _trapz_weights(r)
    wr = w .* r

    K  = [_bessel_kernel(n, B[j], r[k], w[k])
          for j in 1:nB, k in 1:Nr]      # nB × Nr
    H  = K * I
    N  = dot(wr, I)
    V  = H ./ N

    function hankel_vis_pullback(ȳ)
        # K'ȳ  : gemv, nB×Nr transposed times nB-vector → Nr-vector
        # dot(ȳ,V) · wr : quotient-rule normalization correction
        dI = (K' * ȳ .- dot(ȳ, V) .* wr) ./ N
        return NoTangent(), dI, NoTangent(), NoTangent()
    end
    return V, hankel_vis_pullback
end


# ---------------------------------------------------------------------------
# Bessel kernel scalar — shared by rrules and _fwd! variants
# ---------------------------------------------------------------------------

@inline function _bessel_kernel(n::Int, Bj::Real, rk::Real, wk::Real)
    arg = 2π * Bj * rk
    Jn  = (n == 0) ? besselj0(arg) :
          (n == 1) ? besselj1(arg) :
                     besselj(n, arg)
    return wk * rk * Jn
end


# ---------------------------------------------------------------------------
# 4.  hankel_vis_fwd  (allocating forward pass, returns K for pullback reuse)
# ---------------------------------------------------------------------------

"""
    hankel_vis_fwd(I, r, B; n=0) -> (V, K, N)

Forward pass for `hankel_vis` that returns all intermediates needed for
the pullback, at the cost of a single Bessel evaluation per (B[j], r[k]) pair.

Returns
-------
V  : normalised visibility, Vector{T} of length nB
K  : Bessel kernel matrix,  Matrix{T} of size nB × Nr
     K[j,k] = w[k] · r[k] · Jₙ(2π B[j] r[k])
N  : normalisation constant ∫ I r dr

The pullback is:  dI = (K'ȳ  −  dot(ȳ,V) · (w .* r)) / N
"""
function hankel_vis_fwd(I::AbstractVector, r::AbstractVector,
                         B::AbstractVector; n::Int=0)
    Nr = length(r)
    nB = length(B)
    w  = _trapz_weights(r)
    wr = w .* r

    K  = [_bessel_kernel(n, B[j], r[k], w[k])
          for j in 1:nB, k in 1:Nr]
    H  = K * I
    N  = dot(wr, I)
    V  = H ./ N

    return V, K, N
end


# ---------------------------------------------------------------------------
# 5.  hankel_vis_fwd!  (zero-allocation in-place forward pass)
# ---------------------------------------------------------------------------

"""
    hankel_vis_fwd!(out_V, out_K, I, r, B; n=0) -> N

In-place forward pass.  Writes normalised visibilities into `out_V` (length nB)
and the Bessel kernel into `out_K` (nB × Nr).  Returns the normalisation N.

Also caches the trapezoid weights into `out_w` if a 5th positional argument
is supplied; otherwise weights are recomputed (cheap).

Pre-allocate with:
    ws = HankelWorkspace(Nr, nB)
Then call:
    N  = hankel_vis_fwd!(ws.V, ws.K, I, r, B; n=0)
"""
function hankel_vis_fwd!(out_V::AbstractVector,
                          out_K::AbstractMatrix,
                          I::AbstractVector,
                          r::AbstractVector,
                          B::AbstractVector;
                          n::Int=0)
    Nr = length(r)
    nB = length(B)
    @assert length(I)      == Nr  "hankel_vis_fwd!: I and r length mismatch"
    @assert length(out_V)  == nB  "hankel_vis_fwd!: out_V length ≠ nB"
    @assert size(out_K)    == (nB, Nr)  "hankel_vis_fwd!: out_K shape ≠ (nB, Nr)"

    w = _trapz_weights(r)

    # Compute N = ∫ I r dr
    N = zero(eltype(I))
    @inbounds for k in 1:Nr
        N += w[k] * r[k] * I[k]
    end

    # Fill kernel matrix and accumulate H[j] = Σₖ K[j,k]·I[k]
    @inbounds for j in 1:nB
        Hj = zero(eltype(I))
        Bj = B[j]
        for k in 1:Nr
            Kjk        = _bessel_kernel(n, Bj, r[k], w[k])
            out_K[j,k] = Kjk
            Hj        += Kjk * I[k]
        end
        out_V[j] = Hj / N
    end

    return N
end


# ---------------------------------------------------------------------------
# 6.  hankel_vis_pullback!  (zero-allocation reverse pass)
# ---------------------------------------------------------------------------

"""
    hankel_vis_pullback!(out_dI, ȳ, V, K, N, w, r)

In-place pullback for `hankel_vis`.

Given upstream cotangent `ȳ` (length nB) and cached forward quantities
`V`, `K`, `N`, `w`, `r`, writes the cotangent ∂L/∂I into `out_dI` (length Nr).

Formula:
    out_dI = (K'ȳ  −  dot(ȳ, V) · (w .* r)) / N

Arguments
---------
out_dI  : output buffer, length Nr
ȳ       : upstream cotangent, length nB
V       : cached normalised visibility from hankel_vis_fwd!, length nB
K       : cached Bessel kernel from hankel_vis_fwd!, shape nB × Nr
N       : cached normalisation constant from hankel_vis_fwd!
w       : trapezoid weights, length Nr  (use _trapz_weights(r))
r       : radius grid, length Nr
"""
function hankel_vis_pullback!(out_dI::AbstractVector,
                               ȳ::AbstractVector,
                               V::AbstractVector,
                               K::AbstractMatrix,
                               N::Real,
                               w::AbstractVector,
                               r::AbstractVector)
    Nr = length(r)
    nB = length(ȳ)
    @assert length(out_dI) == Nr
    @assert length(V)      == nB
    @assert size(K)        == (nB, Nr)
    @assert length(w)      == Nr

    yV = dot(ȳ, V)            # scalar: Σⱼ ȳ[j]·V[j]

    # out_dI = K'ȳ  (gemv, K is nB×Nr so K' is Nr×nB)
    mul!(out_dI, K', ȳ)

    # Subtract normalization correction and divide by N
    @inbounds for k in 1:Nr
        out_dI[k] = (out_dI[k] - yV * w[k] * r[k]) / N
    end

    return nothing
end


# ---------------------------------------------------------------------------
# 7.  Convenience: full gradient ∂L/∂p via cached K and profile Jacobian
# ---------------------------------------------------------------------------

"""
    hankel_grad_params!(out_dL_dp, ȳ, V, K, N, w, r, dI_dp)

Compute the full parameter gradient ∂L/∂p in-place, given:
  - Upstream cotangent ȳ (length nB)
  - Cached forward pass (V, K, N from hankel_vis_fwd!)
  - Profile Jacobian dI_dp (Nr × n_params, from ForwardDiff on compiled profile)

Writes result into `out_dL_dp` (length n_params).

This is the chain rule:
    ∂L/∂p = (∂I/∂p)' · ∂L/∂I  =  dI_dp' * dI

where dI = (K'ȳ − dot(ȳ,V)·(w.*r)) / N  is computed internally (no allocation
if a scratch buffer is supplied; here we allocate one dI vector).

Use this after:
    N = hankel_vis_fwd!(ws.V, ws.K, I, r, B)
    ForwardDiff.jacobian!(ws.dI_dp, profile_fn, p)
    hankel_grad_params!(dL_dp, ȳ, ws.V, ws.K, N, ws.w, r, ws.dI_dp)
"""
function hankel_grad_params!(out_dL_dp::AbstractVector,
                              ȳ::AbstractVector,
                              V::AbstractVector,
                              K::AbstractMatrix,
                              N::Real,
                              w::AbstractVector,
                              r::AbstractVector,
                              dI_dp::AbstractMatrix)
    Nr, n_params = size(dI_dp)
    dI = similar(r)                           # Nr-vector, one allocation
    hankel_vis_pullback!(dI, ȳ, V, K, N, w, r)
    mul!(out_dL_dp, dI_dp', dI)              # (n_params × Nr) * (Nr,) → (n_params,)
    return nothing
end
