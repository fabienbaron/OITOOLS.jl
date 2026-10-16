# Posterior samplers for model fitting, beside nested sampling: Ensemble Slice Sampling (here,
# pure Julia) and non-reversible parallel tempering (`fit_model_pigeons`, in the Pigeons
# extension). Both from ROTIR (`src/ensemble_slice.jl`, `src/mcmc_diagnostics.jl`,
# `src/fit_pigeons.jl`), on OITOOLS' models: the likelihood is −(χ² + constraint penalty)/2 and
# the bounds are the prior, uniform, as in `fit_model_nested`.

using Random, Statistics, Printf

"""
    SamplerResult

What [`fit_model_ensemble`](@ref) and [`fit_model_pigeons`](@ref) return.

| field | |
|---|---|
| `x_opt` | the highest-likelihood sample (as the nested samplers report) |
| `list_free_params` | parameter names, in the order of `x_opt` and the posterior's columns |
| `chi2`, `chi2r`, `ndof` | the data χ² at `x_opt` (no constraint penalty), and per degree of freedom |
| `posterior` | equal-weight samples, `(n_samples, n_params)` |
| `median`, `q16`, `q84` | per parameter, from `posterior` |
| `logz`, `logzerr` | the evidence on `fit_model_nested`'s scale; `NaN` when the sampler gives none |
| `diagnostics` | a NamedTuple: R̂ and τ (ensemble), round trips and chains (Pigeons) |
| `result` | the sampler's own output |
| `model` | the compiled `FlatModel` |
| `sampler` | `:ensemble` or `:pigeons` |
"""
struct SamplerResult
    x_opt            :: Vector{Float64}
    list_free_params :: Vector{String}
    chi2             :: Float64
    chi2r            :: Float64
    ndof             :: Int
    posterior        :: Matrix{Float64}
    median           :: Vector{Float64}
    q16              :: Vector{Float64}
    q84              :: Vector{Float64}
    logz             :: Float64
    logzerr          :: Float64
    diagnostics      :: NamedTuple
    result           :: Any
    model            :: FlatModel
    sampler          :: Symbol
end

function Base.show(io::IO, r::SamplerResult)
    @printf(io, "SamplerResult (%s): χ²r = %.4f  (chi2=%.2f, ndof=%d), %d samples\n",
            r.sampler, r.chi2r, r.chi2, r.ndof, size(r.posterior, 1))
    isfinite(r.logz) && @printf(io, "  log(Z) = %.2f ± %.2f\n", r.logz, r.logzerr)
    for (i, p) in enumerate(r.list_free_params)
        @printf(io, "  %-20s = %10.4f  (median %.4f, %.4f … %.4f)\n", p, r.x_opt[i],
                r.median[i], r.q16[i], r.q84[i])
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Ensemble Slice Sampling (Karamanis & Beutler 2020, the `zeus` algorithm)
# ─────────────────────────────────────────────────────────────────────────────
#
# Gradient-free and affine-invariant. Derivative-free affine-invariant ensembles need N = 2d
# walkers (Chen 2025, arXiv:2505.02987): the directions are walker differences, so with K
# walkers the reachable set is a (K-1)-dimensional affine hull. There is no accept/reject: the
# shrinking procedure always returns a point inside the slice.

"""
    ess_step(logp, x, lp_x, d, logy; max_expand, maxshrink, rng) -> (x, lp, ne, nc, capped)

One slice-sampling step along `d` from `x`, whose log-density `lp_x` is known. Steps the
bracket out until both ends leave the slice, then shrinks toward the current point until a
sample inside it is found.

The expansion budget is split RANDOMLY between the two sides, `J ~ U{0..B-1}`, `K = B-1-J`
(Neal 2003, Fig. 3), which keeps a limited stepping-out reversible.
"""
function ess_step(logp, x::AbstractVector{Float64}, lp_x::Float64, d::AbstractVector{Float64},
                  logy::Float64; max_expand::Int = 10_000, maxshrink::Int = 10_000,
                  rng = Random.default_rng())
    nd = length(x)
    xp = Vector{Float64}(undef, nd)
    J = rand(rng, 0:(max_expand - 1))
    K = (max_expand - 1) - J
    u = rand(rng)
    L = -u; R = L + 1.0
    ne = 0; capped = false
    while J > 0
        @. xp = x + L * d
        logp(xp) > logy || break
        L -= 1.0; ne += 1; J -= 1
    end
    J == 0 && (capped = true)
    while K > 0
        @. xp = x + R * d
        logp(xp) > logy || break
        R += 1.0; ne += 1; K -= 1
    end
    K == 0 && (capped = true)
    nc = 0
    while nc < maxshrink
        t = L + rand(rng) * (R - L)
        @. xp = x + t * d
        lp = logp(xp)
        lp > logy && return copy(xp), lp, ne, nc, capped
        t < 0 ? (L = t) : (R = t)
        nc += 1
    end
    return copy(x), lp_x, ne, nc, true         # degenerate: stay put rather than return garbage
end

"""
    ensemble_slice(logp, x0; nsweep = 5000, burn = 2000, thin = 10, mu0 = 1.0, tune = true,
                   max_expand, maxshrink, patience = 5, tol = 0.05, rng, verb = true)
        -> (chain, logps, diagnostics)

Ensemble Slice Sampling of `logp` from walkers `x0`, `(nwalkers, ndim)`, an even number of them.
`chain` is `(n_kept_sweeps, nwalkers, ndim)` and `logps` `(n_kept_sweeps, nwalkers)`, the
walkers' log-densities.

Each half of the walkers is updated with directions drawn from the other half; the halves are
reshuffled every sweep (zeus's `shuffle_ensemble=True`). `mu` is tuned during burn-in toward
equal numbers of expansions and contractions, and frozen once the balance has held `patience`
times, as in zeus.
"""
function ensemble_slice(logp, x0::Matrix{Float64}; nsweep::Int = 5000, burn::Int = 2000,
                        thin::Int = 10, mu0::Float64 = 1.0, tune::Bool = true,
                        max_expand::Int = 10_000, maxshrink::Int = 10_000,
                        patience::Int = 5, tol::Float64 = 0.05,
                        rng = Random.default_rng(), verb::Bool = true)
    nw, nd = size(x0)
    isodd(nw) && error("ensemble_slice needs an even number of walkers (red-blue split)")
    nw < 2nd && @warn "nw = $nw is below the N = 2d = $(2nd) that affine-invariant ensembles " *
                      "need (Chen 2025); directions span at most nw/2 - 1 dimensions per half."
    X = copy(x0); mu = mu0; half = nw ÷ 2
    lpX = [logp(@view X[i, :]) for i in 1:nw]
    nev = Ref(nw)
    wrapped = z -> (nev[] += 1; logp(z))
    kept = Matrix{Float64}[]; keptlp = Vector{Float64}[]
    nexp = 0; ncon = 0; ncap = 0; tuning = tune; nbal = 0; order = collect(1:nw)
    for sweep in 1:(burn + nsweep)
        shuffle!(rng, order)
        Ne = 0; Nc = 0
        for (a, b) in ((1:half, half+1:nw), (half+1:nw, 1:half))
            act = @view order[a]; comp = @view order[b]
            for i in act
                j = rand(rng, comp); k = rand(rng, comp)
                while k == j && length(comp) > 1; k = rand(rng, comp); end
                d = (2.0 * mu) .* (@view(X[j, :]) .- @view(X[k, :]))   # zeus's factor 2
                if !any(!iszero, d)
                    # a degenerate direction, rescued with the complementary half's spread only
                    d = [max(std(@view X[comp, q]), 1e-12) * randn(rng) for q in 1:nd]
                end
                logy = lpX[i] - randexp(rng)
                xn, lpn, ne, nc, cap = ess_step(wrapped, @view(X[i, :]), lpX[i], d, logy;
                                                max_expand, maxshrink, rng)
                X[i, :] = xn; lpX[i] = lpn
                Ne += ne; Nc += nc; cap && (ncap += 1)
            end
        end
        nexp += Ne; ncon += Nc
        if tuning && sweep <= burn && (Ne + Nc) > 0
            mu *= 2 * max(Ne, 1) / (max(Ne, 1) + Nc)   # zeus floors Ne, which bounds the factor
            abs(Ne / (Ne + Nc) - 0.5) < tol && (nbal += 1)
            nbal > patience && (tuning = false)
        end
        if sweep > burn && (sweep - burn) % thin == 0
            push!(kept, copy(X)); push!(keptlp, copy(lpX))
        end
        verb && sweep % 500 == 0 &&
            (@printf("  sweep %5d/%d  mu %.4f%s  Ne/Nc %.2f  evals %d\n", sweep, burn + nsweep,
                     mu, tuning ? "" : " (frozen)", Ne / max(Nc, 1), nev[]); flush(stdout))
    end
    chain = Array{Float64,3}(undef, length(kept), nw, nd)
    logps = Matrix{Float64}(undef, length(kept), nw)
    for (si, Xs) in enumerate(kept)
        chain[si, :, :] = Xs; logps[si, :] = keptlp[si]
    end
    return chain, logps,
           (mu = mu, tuned = !tuning, n_expand = nexp, n_contract = ncon, n_eval = nev[],
            n_capped = ncap, nkept = length(kept) * nw, nsweep_kept = length(kept), nwalkers = nw)
end

"""
    rhat(chain) -> Vector

Gelman-Rubin R̂ per parameter from a `(nsweep, nwalkers, ndim)` chain, each walker a chain.
Above about 1.01 the chain has not converged. The walkers of an ensemble are coupled, so this R̂
can be low while the ensemble as a whole drifts.
"""
function rhat(chain::AbstractArray{<:Real,3})
    m, nw, d = size(chain); out = zeros(float(eltype(chain)), d)
    for p in 1:d
        W = mean(var(@view chain[:, w, p]) for w in 1:nw)
        B = m * var([mean(@view chain[:, w, p]) for w in 1:nw])
        out[p] = sqrt(((m - 1) / m * W + B / m) / max(W, eps()))
    end
    return out
end

"""
    iat(chain; per_walker = false) -> Float64

Integrated autocorrelation time in sweeps, averaged over parameters, from a `(nsweep, nwalkers,
ndim)` chain with Sokal's `t > 5τ` window. By default the walkers are concatenated into one
series per parameter; `per_walker = true` is biased low (−24 % at τ = 99 on an AR(1) series).
A τ above about `nsweep/50` is a lower bound.
"""
function iat(chain::AbstractArray{<:Real,3}; per_walker::Bool = false)
    m, nw, d = size(chain); out = Float64[]
    for p in 1:d
        if per_walker
            acc = zeros(m)
            for w in 1:nw
                y = Float64.(chain[:, w, p]); y .-= mean(y)
                n0 = max(sum(abs2, y) / m, eps())
                acc .+= [sum(y[1:m-t] .* y[1+t:m]) / (m - t) for t in 0:m-1] ./ n0
            end
            push!(out, _sokal_window(acc ./ nw))
        else
            y = Float64.(vec(chain[:, :, p])); y .-= mean(y); n = length(y)
            f = [sum(y[1:n-t] .* y[1+t:n]) for t in 0:min(n - 1, 10m)] ./ n
            push!(out, _sokal_window(f ./ max(f[1], eps())))
        end
    end
    return mean(out)
end

function _sokal_window(acf::AbstractVector)
    τ = 1.0
    for t in 2:length(acf)
        τ += 2 * acf[t]
        t > 5τ && break
    end
    return max(τ, 1.0)
end

"Flatten a `(nsweep, nwalkers, ndim)` chain to rows `(nsweep·nwalkers, ndim)`."
flatten_chain(c::AbstractArray{<:Real,3}) = reshape(c, size(c, 1) * size(c, 2), size(c, 3))

# ─────────────────────────────────────────────────────────────────────────────
# Shared by both samplers
# ─────────────────────────────────────────────────────────────────────────────

"""
    _box_transform(lo, hi) -> (to_θ, to_z, logabsjac, dθdz, dljac)

An unconstrained reparameterisation of a box and its log-Jacobian (ROTIR's): a scaled logistic
for two finite bounds, `lo + exp(z)` for a finite lower bound only, the identity otherwise.
`dθdz(z)` and `dljac(z)` are the derivatives the gradient of the transformed log-posterior needs.
"""
function _box_transform(lo::AbstractVector, hi::AbstractVector)
    n = length(lo)
    kinds = [isfinite(lo[j]) && isfinite(hi[j]) ? :box :
             isfinite(lo[j])                    ? :pos : :free for j in 1:n]
    logistic(x) = 1 / (1 + exp(-x))
    to_θ(z) = [kinds[j] === :box ? lo[j] + (hi[j] - lo[j]) * logistic(z[j]) :
               kinds[j] === :pos ? lo[j] + exp(z[j]) : z[j] for j in 1:n]
    to_z(θ) = [kinds[j] === :box ? log((θ[j] - lo[j]) / max(hi[j] - θ[j], eps())) :
               kinds[j] === :pos ? log(max(θ[j] - lo[j], eps())) : θ[j] for j in 1:n]
    function logabsjac(z)
        s = 0.0
        @inbounds for j in 1:n
            if kinds[j] === :box
                u = logistic(z[j]); s += log(hi[j] - lo[j]) + log(u) + log1p(-u)
            elseif kinds[j] === :pos
                s += z[j]
            end
        end
        return s
    end
    dθdz(z) = [kinds[j] === :box ? (u = logistic(z[j]); (hi[j] - lo[j]) * u * (1 - u)) :
               kinds[j] === :pos ? exp(z[j]) : 1.0 for j in 1:n]
    dljac(z) = [kinds[j] === :box ? 1 - 2 * logistic(z[j]) : kinds[j] === :pos ? 1.0 : 0.0
                for j in 1:n]
    return to_θ, to_z, logabsjac, dθdz, dljac
end

# The compiled model, the constraints and the bounds of a sampled fit, as `fit_model_nested`
# sets them up.
function _sampler_setup(model_dict::Dict{String}, list_free_params::Vector{String}, data::OIdata;
                        lb, ub, constraints, nB_workspace = nothing)
    nB = something(nB_workspace, size(data.uv, 2))
    cons = parse_constraints(constraints)
    augmented = copy(model_dict)
    constraint_keys = augment_constraints!(augmented, cons)
    fm = parse_model(augmented, list_free_params; nB_workspace = nB)
    constraint_idx = [findfirst(==(k), fm.all_names) for k in constraint_keys]
    lower = Float64[get(lb, p, -Inf) for p in list_free_params]
    upper = Float64[get(ub, p, Inf) for p in list_free_params]
    x0 = Float64[Float64(model_dict[p]) for p in list_free_params]
    all(lower .<= x0 .<= upper) ||
        error("the starting values are outside the bounds: " *
              join([p for (i, p) in enumerate(list_free_params) if !(lower[i] <= x0[i] <= upper[i])], ", "))
    return fm, cons, constraint_idx, lower, upper, x0
end

# −(χ² + penalty)/2: the log-likelihood both samplers sample, the bounds being the prior.
_sampler_loglike(fm, data, weights, cons, idx) =
    θ -> -0.5 * (chi2_flat(fm, θ, data; weights) + constraint_penalty(θ, cons, fm, idx))

# The same with its gradient, from the closed-form χ² gradient.
function _sampler_loglike_grad(fm, data, weights, cons, idx)
    return function (θ)
        c, g = chi2_flat_fg(fm, θ, data; weights)
        g = Vector{Float64}(g)
        pen = constraint_penalty_and_grad!(g, θ, cons, fm, idx)
        return -0.5 * (Float64(c) + pen), -0.5 .* g
    end
end

# The result both samplers return, from equal-weight samples Θ (rows) and their
# log-likelihoods.
function _sampler_result(Θ::AbstractMatrix, logls::AbstractVector, fm, list_free_params, data,
                         weights; logz, logzerr, diagnostics, result, sampler)
    ib = argmax(logls)
    x_opt = Vector{Float64}(Θ[ib, :])
    chi2 = Float64(chi2_flat(fm, x_opt, data; weights))
    ndof = _ndof(data, weights)
    q(pr) = [quantile(view(Θ, :, j), pr) for j in axes(Θ, 2)]
    return SamplerResult(x_opt, list_free_params, chi2, chi2 / ndof, ndof, Matrix{Float64}(Θ),
                         q(0.5), q(0.16), q(0.84), Float64(logz), Float64(logzerr), diagnostics,
                         result, fm, sampler)
end

# ─────────────────────────────────────────────────────────────────────────────
# fit_model_ensemble
# ─────────────────────────────────────────────────────────────────────────────

"""
    fit_model_ensemble(model_dict, list_free_params, data; lb, ub, weights, constraints,
                       nwalkers = 0, nsweep = 2000, burn = nsweep ÷ 2, thin = 0, seed,
                       verb = true) -> SamplerResult

Sample the posterior of a model's free parameters with Ensemble Slice Sampling
([`ensemble_slice`](@ref); Karamanis & Beutler 2020), gradient-free and affine-invariant.

The log-likelihood is −(χ² + constraint penalty)/2; the prior is uniform within the bounds,
which may be infinite (an improper prior).

- `nwalkers`: `0` for 2d (d free parameters), at least 4, rounded up to an even number.
- `nsweep` sweeps kept after `burn` of burn-in, every `thin`-th (`0`: about 2000 kept sweeps at
  most). The walkers start in a small ball about the model's values, inside the bounds.
- `seed` makes the run reproducible.

`diagnostics` carries R̂ per parameter ([`rhat`](@ref); above 1.01, not converged) and τ, the
integrated autocorrelation time in sweeps ([`iat`](@ref); above about `nsweep/50`, a lower
bound), and the number of evaluations. `logz` is `NaN`.
"""
function fit_model_ensemble(model_dict::Dict{String}, list_free_params::Vector{String},
                            data::OIdata;
                            lb = Dict{String,Float64}(), ub = Dict{String,Float64}(),
                            weights = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0], constraints = [],
                            nB_workspace = nothing, nwalkers::Integer = 0,
                            nsweep::Integer = 2000, burn::Integer = nsweep ÷ 2,
                            thin::Integer = 0, seed = nothing, verb::Bool = true)
    fm, cons, idx, lower, upper, x0 = _sampler_setup(model_dict, list_free_params, data;
                                                     lb, ub, constraints, nB_workspace)
    nd = length(x0)
    nw = nwalkers > 0 ? Int(nwalkers) : max(4, 2nd)
    isodd(nw) && (nw += 1)
    rng = seed === nothing ? Random.default_rng() : Random.Xoshiro(seed)
    loglike = _sampler_loglike(fm, data, weights, cons, idx)
    logp = θ -> begin
        all(lower[i] <= θ[i] <= upper[i] for i in eachindex(θ)) || return -Inf
        v = loglike(θ)
        isfinite(v) ? Float64(v) : -Inf
    end
    # A small ball about the start, inside the box: a walker starting at −∞ could never move.
    span = [isfinite(upper[i] - lower[i]) ? min(upper[i] - lower[i], max(abs(x0[i]), 1.0)) :
            max(abs(x0[i]), 1.0) for i in 1:nd]
    X0 = Matrix{Float64}(undef, nw, nd)
    for w in 1:nw, i in 1:nd
        lo = isfinite(lower[i]) ? lower[i] + 1e-9 * span[i] : -Inf
        hi = isfinite(upper[i]) ? upper[i] - 1e-9 * span[i] : Inf
        X0[w, i] = clamp(x0[i] + 0.02 * span[i] * randn(rng), lo, hi)
    end
    th = thin > 0 ? Int(thin) : max(1, nsweep ÷ 2000)
    verb && @printf("ensemble slice sampling: %d walkers, %d + %d sweeps, every %d kept\n",
                    nw, burn, nsweep, th)
    chain, logps, dg = ensemble_slice(logp, X0; nsweep = Int(nsweep), burn = Int(burn),
                                      thin = th, rng, verb)
    Θ = flatten_chain(chain)
    R = rhat(chain); τ = iat(chain)
    if verb
        @printf("  %d samples, %d evaluations; R̂ %.4f%s, τ %.0f sweeps%s\n", size(Θ, 1),
                dg.n_eval, maximum(R), maximum(R) > 1.01 ? " (NOT converged)" : "", τ,
                τ > nsweep / 50 ? " (lower bound)" : "")
    end
    return _sampler_result(Θ, vec(logps), fm, list_free_params, data, weights;
                           logz = NaN, logzerr = NaN,
                           diagnostics = (rhat = R, tau = τ, n_eval = dg.n_eval,
                                          nwalkers = nw, mu = dg.mu, tuned = dg.tuned),
                           result = (chain = chain, logps = logps, info = dg),
                           sampler = :ensemble)
end

# ─────────────────────────────────────────────────────────────────────────────
# fit_model_pigeons: declared here, implemented in the Pigeons extension
# ─────────────────────────────────────────────────────────────────────────────

function _fit_model_pigeons end

"""
    fit_model_pigeons(model_dict, list_free_params, data; lb, ub, weights, constraints,
                      n_chains = 10, n_rounds = 10, explorer = :slice, variational = false,
                      n_chains_variational = 0, reference_sigma = 3.0, seed,
                      multithreaded = Threads.nthreads() > 1, verb = true) -> SamplerResult

Sample the posterior with non-reversible parallel tempering (Syed et al. 2021) through
Pigeons.jl, for multimodal posteriors. Needs `using Pigeons`.

The free parameters are mapped to ℝ (a logistic for two finite bounds, `lo + exp(z)` for one,
the identity for none) and the ladder anneals from N(0, `reference_sigma`²) on those
coordinates. The log-likelihood is −(χ² + constraint penalty)/2; the bounds are the prior.

- `n_rounds` is a power: the run makes 2^`n_rounds` scans, each round tuning the ladder for the
  next. `n_chains` ≈ 2Λ (the global communication barrier) is the sizing rule.
- `explorer`: `:slice` (no gradient) or `:mala` (AutoMALA, on the closed-form χ² gradient).
- `variational = true` fits a Gaussian reference during the run; `n_chains_variational > 0`
  adds a leg with the fixed reference (the stabilised variant, Surjanovic et al. 2022). On AX
  Cir's binary (4 parameters, bounds ±60 mas) the fixed reference gives 0–1 round trips at 8–10
  rounds, the variational or stabilised reference 76–352.

`logz` is the stepping-stone evidence on `fit_model_nested`'s scale (the log volume of the
bounds subtracted), `NaN` when a bound is infinite. `logzerr` is 1/√(round trips) and does not
include the bias of a short run: on α Cen A's power-law fit (Nautilus: −712.528), −711.71 at 8
rounds (1 round trip), −712.29 at 10 (18), −712.536 at 12 (103). `diagnostics` holds the number
of round trips; with fewer than 3, other modes may not have been visited.
"""
function fit_model_pigeons(args...; kwargs...)
    isempty(methods(_fit_model_pigeons)) &&
        error("fit_model_pigeons needs Pigeons: `using Pigeons`")
    return _fit_model_pigeons(args...; kwargs...)
end
