# Levenberg-Marquardt with bounds that are honoured, not clipped.
#
# The iteration is LsqFit's (`LsqFit/src/levenberg_marquardt.jl`, MIT): the same damped normal
# equations, the same step-quality test and the same λ schedule and constants, so a fit that
# touches no bound takes the steps LsqFit would. Two things differ.
#
# BOUNDS. LsqFit solves the step for every parameter and then clips it at the box. A parameter
# on its bound whose gradient points out of the box is clipped to no motion, but the other
# components of the step were solved assuming it moved, so the step is poor, λ never falls and
# the fit crawls to the iteration cap. Measured on PMOIRED's FS CMa ring (EX2 step 4, `diamin`
# on its bound 0): 200 iterations, χ² 8577.07, "not converged"; with that parameter held, 26
# iterations to χ² 8574.98. Here such a parameter is HELD for the iteration: left out of the
# solve, so the others take the step of the problem that is actually being solved. It rejoins
# as soon as its gradient turns back into the box.
#
# STOPPING. LsqFit's two tests are kept (‖g‖∞ below an absolute `g_tol`, the step below
# `x_tol` relative), and a third is added in units of the noise: the χ² decrease the undamped
# Gauss-Newton step would still make, `gᵀ(JᵀJ)⁻¹g` over the free parameters, below
# `chi2_tol · max(1, χ²/n)`. Residuals here are already divided by their errors, so a decrease
# of 1e-3 is a thousandth of what one data point contributes on average; σ is scaled by χ²/n
# the same way, so it is equally a step far below every σ.

"""
    _lm_bounded(f, jac, x0, lower, upper; maxIter, x_tol, g_tol, chi2_tol)
        -> (; x, resid, jacobian, converged, reason, iterations, trace, held)

Minimise `‖f(x)‖²` from `x0` within `lower ≤ x ≤ upper` (`±Inf` for none). `f(x)` returns the
residual vector, `jac(x)` its Jacobian.

`reason` says what stopped it: `:chi2` (the step could no longer move χ² by `chi2_tol` in units
of the noise; `chi2_tol = 0` turns the test off), `:gradient`, `:step`, `:bounds` (every
parameter is held on a bound), or `:maxiter`; `converged` is false only for the last. `held` is
the indices of the parameters that end on a bound with the gradient pushing them out of it, the
ones the minimum is constrained by. `trace` has one entry per iteration,
`(; iteration, chi2, gnorm, lambda)`.

The parameters keep `x0`'s element type, and the arithmetic is done in the wider of that and
the residuals' type.
"""
function _lm_bounded(f, jac, x0::AbstractVector{<:AbstractFloat}, lower::AbstractVector,
                     upper::AbstractVector; maxIter::Integer = 200, x_tol::Real = 1e-8,
                     g_tol::Real = 1e-12, chi2_tol::Real = 1e-3)
    n = length(x0)
    length(lower) == n && length(upper) == n ||
        throw(ArgumentError("Bounds must be of the same length as the number of parameters."))
    all(lower .<= x0 .<= upper) || throw(ArgumentError("Initial guess must be within bounds."))
    x = copy(x0)
    # The first evaluation outside the loop: `_lm_loop` is compiled for the residual's and the
    # Jacobian's concrete types, however well the caller's closures infer.
    return _lm_loop(f, jac, x, f(x), jac(x), lower, upper, Int(maxIter), x_tol, g_tol, chi2_tol)
end

function _lm_loop(f, jac, x::AbstractVector{P}, r0::AbstractVector, J0::AbstractMatrix,
                  lower, upper, maxIter::Int, x_tol, g_tol, chi2_tol) where {P}
    T = promote_type(P, eltype(r0), eltype(J0))
    n, m = length(x), length(r0)
    # The residuals and the Jacobian are copied into buffers of T as they arrive, so every
    # product below is BLAS whatever mix of precisions the model returns.
    r = Vector{T}(r0); J = Matrix{T}(J0); rt = similar(r)
    lo = convert(Vector{P}, lower); hi = convert(Vector{P}, upper)
    # LsqFit's constants.
    λ = T(10)
    λ_increase, λ_decrease = T(10), T(0.1)
    min_step_quality, good_step_quality = T(1e-3), T(0.75)
    MAX_LAMBDA, MIN_LAMBDA, MIN_DIAGONAL = T(1e16), T(1e-16), T(1e-6)
    # Below eps(T)·‖x‖ a step cannot change x; LsqFit's 1e-8 is below Float32's resolution.
    xtol = max(T(x_tol), eps(T))
    gtol, ctol = T(g_tol), T(chi2_tol)

    g = zeros(T, n); held = falses(n)
    JJ = zeros(T, n, n); A = zeros(T, n, n); DtD = zeros(T, n)
    δ = zeros(T, n); rhs = zeros(T, n); xt = similar(x); Jδ = zeros(T, m)
    S = sum(abs2, r)
    trace = NamedTuple{(:iteration, :chi2, :gnorm, :lambda),Tuple{Int,T,T,T}}[]
    reason = :maxiter
    iter = 0
    fresh = true                     # J and JᵀJ belong to the current x
    while iter < maxIter
        if fresh
            mul!(g, transpose(J), r)
            mul!(JJ, transpose(J), J)
            @inbounds for j in 1:n
                DtD[j] = max(JJ[j, j], MIN_DIAGONAL)
            end
            fresh = false
        end
        _lm_held!(held, x, g, lo, hi)
        if all(held)
            reason = :bounds
            break
        end
        gnorm = zero(T)
        @inbounds for i in 1:n
            held[i] || (gnorm = max(gnorm, abs(g[i])))
        end
        if gnorm < gtol
            reason = :gradient
            break
        end
        # The χ² the undamped Gauss-Newton step would still remove, over the free parameters.
        # A singular JᵀJ (a parameter the data does not constrain) leaves the test to the others.
        if ctol > 0
            _lm_system!(A, rhs, JJ, g, held, zero(T), DtD)
            C = cholesky!(Symmetric(A), check = false)
            if issuccess(C)
                ldiv!(C, rhs)
                dec = zero(T)
                @inbounds for i in 1:n
                    held[i] || (dec -= g[i] * rhs[i])       # rhs = −(JᵀJ)⁻¹g
                end
                if dec < ctol * max(one(T), S / m)
                    reason = :chi2
                    break
                end
            end
        end

        # (JᵀJ + λ·diag(JᵀJ)) δ = −g over the free parameters; a held one has δ = 0 exactly.
        _lm_system!(A, rhs, JJ, g, held, λ, DtD)
        C = cholesky!(Symmetric(A), check = false)
        if issuccess(C)
            ldiv!(δ, C, rhs)
        else
            _lm_system!(A, rhs, JJ, g, held, λ, DtD)
            δ .= A \ rhs
        end
        # A free parameter may still reach a bound within the step: it stops there, and is held
        # from the next iteration if the gradient keeps pushing it out.
        @inbounds for i in 1:n
            xt[i] = clamp(x[i] + δ[i], lo[i], hi[i])
            δ[i] = xt[i] - x[i]
        end

        mul!(Jδ, J, δ)
        predicted = zero(T)
        @inbounds for k in 1:m
            predicted += abs2(Jδ[k] + r[k])
        end
        copyto!(rt, f(xt))
        St = sum(abs2, rt)
        ρ = (St - S) / (predicted - S)
        if St < S && ρ > min_step_quality
            x, xt = xt, x
            r, rt = rt, r
            S = St
            copyto!(J, jac(x))
            fresh = true
            ρ > good_step_quality && (λ = max(λ_decrease * λ, MIN_LAMBDA))
        else
            λ = min(λ_increase * λ, MAX_LAMBDA)
        end
        iter += 1
        push!(trace, (; iteration = iter, chi2 = S, gnorm, lambda = λ))
        if norm(δ) < xtol * (xtol + norm(x))
            reason = :step
            break
        end
    end
    mul!(g, transpose(J), r)
    _lm_held!(held, x, g, lo, hi)
    return (; x, resid = r, jacobian = J, converged = reason !== :maxiter, reason,
              iterations = iter, trace, held = findall(held))
end

# A = JᵀJ + λ·diag(DtD) and rhs = −g, with a held parameter's row and column replaced by the
# identity and its right-hand side by 0, so that it does not move and the others solve the
# problem with it fixed — the same system as removing it, without allocating a smaller one.
function _lm_system!(A, rhs, JJ, g, held, λ, DtD)
    n = length(g)
    @inbounds for j in 1:n, i in 1:n
        A[i, j] = (held[i] || held[j]) ? (i == j ? one(eltype(A)) : zero(eltype(A))) :
                  (i == j ? JJ[i, i] + λ * DtD[i] : JJ[i, j])
    end
    @inbounds for i in 1:n
        rhs[i] = held[i] ? zero(eltype(rhs)) : -g[i]
    end
    return A
end

# Which parameters sit on a bound with the descent direction −g pointing out of the box.
function _lm_held!(held, x, g, lower, upper)
    @inbounds for i in eachindex(x)
        held[i] = (x[i] <= lower[i] && g[i] > 0) || (x[i] >= upper[i] && g[i] < 0)
    end
    return held
end

# How the stopping reason reads in a result.
const _LM_REASONS = Dict(
    :chi2     => "χ² could no longer change by a thousandth of a data point's share",
    :gradient => "gradient below tolerance",
    :step     => "step below tolerance",
    :bounds   => "every free parameter is held on a bound",
    :maxiter  => "iteration limit reached")
