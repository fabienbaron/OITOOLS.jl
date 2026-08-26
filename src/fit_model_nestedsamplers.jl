# ─────────────────────────────────────────────────────────────────────────────
# The NestedSamplers.jl backend of `fit_model_nested`
# ─────────────────────────────────────────────────────────────────────────────
#
# Included by ext/OITOOLSNestedSamplersExt.jl. Nothing here is public: `fit_model_nested` is
# declared in the core package and documents the keywords.
#
# Pure Julia, which is the whole point of it. This is the nested sampler that survives into a
# PackageCompiler build, where PythonCall and a conda environment cannot go.
#
# Two differences from the UltraNest backend, both forced by the library and both absorbed here
# rather than pushed onto the caller:
#
#   * The likelihood is SCALAR. UltraNest is driven with a vectorised likelihood over a whole
#     batch of live points; `NestedModel` takes one point at a time. The Python round trip goes
#     away with it, so this is not simply the slower of the two — measure before assuming.
#   * The samples come out WEIGHTED. `_resample_equal` turns them into the equally weighted
#     posterior that `NestedResult.posterior` is documented to hold, matching what UltraNest
#     returns directly.

# ─────────────────────────────────────────────────────────────────────────────
# Rejection sampling, a batch at a time, evaluated across threads
# ─────────────────────────────────────────────────────────────────────────────
#
# Why this exists, measured on the α Cen A power-law fit at 16 threads, 7 seeds:
#
#                          wall: min   median   max
#   Proposals.Rejection          3.55     7.49  26.25 s
#   ThreadedRejection(8)         2.10     2.45   5.43 s     3.1x on the median
#
# One `chi2_flat` on that fit is 142.9 µs and the serial run spends 42125 of them at
# 157.6 µs each, so essentially all of the time is the likelihood.
#
# Read the spread, not just the median. Static nested sampling with a MultiEllipsoid bound is
# wildly variable run to run -- a factor of seven across those seven seeds -- and batching
# narrows that as much as it shortens it. UltraNest, for comparison, is metronomic on the same
# fit: 1.94 to 2.01 s.
#
# Note what is NOT the problem. The scalar likelihood carries only 10% overhead above a bare
# `chi2_flat`, so batching to amortise per-call cost can win at most 10%. What batching really
# buys is somewhere to put the other fifteen threads: the stock proposal evaluates one
# candidate at a time and leaves them idle. UltraNest gets this for free because it hands over
# a whole batch of live points and `_batch_loglike` threads it.
#
# EQUIVALENT, NOT IDENTICAL. Within one call this returns exactly what `Rejection` would: the
# candidates come off the same rng in the same order and the batch is scanned in that order,
# so the accepted point is the same one. But the batch consumes draws past the acceptance,
# so the stream diverges on the next call and a run is not reproducible against a serial one.
# The distribution is unchanged -- the extra evaluations are wasted work, not bias -- which is
# what `test_nested.jl` asserts by comparing logZ rather than bit patterns.
#
# BATCH SIZE TRACKS ACCEPTANCE, NOT THREAD COUNT. The mean number of trials per accepted point
# on that fit is about 7.5 (acceptance ≈ 13%), so a batch of 8 wastes 11% of its evaluations
# while a batch of 16 wastes 105% and is slower despite the same 16 threads.
#
# Batching is not what closes the gap to UltraNest, and it cannot be: UltraNest reaches the
# same evidence in 10985 evaluations against this sampler's 42125, because a reactive sampler
# adds live points where the likelihood demands them instead of holding 400 throughout. That
# is an algorithmic difference, and no amount of parallelism touches it.
struct ThreadedRejection <: Proposals.AbstractProposal
    batch   ::Int
    maxiter ::Int
end
ThreadedRejection(batch::Integer = 8; maxiter::Integer = 100_000) =
    ThreadedRejection(Int(batch), Int(maxiter))

# The cadence `Nested` refits its bounds at; the same one `Rejection` asks for.
NestedSamplers.default_update_interval(::ThreadedRejection, ndims) = 1.5

Base.show(io::IO, p::ThreadedRejection) = print(io, "ThreadedRejection(", p.batch, ")")

function (prop::ThreadedRejection)(rng::AbstractRNG, point::AbstractVector, logl_star,
                                   bounds::Bounds.AbstractBoundingSpace, model)
    ncall = 0
    B  = prop.batch
    us = Vector{Vector{Float64}}(undef, B)
    vs = Vector{Vector{Float64}}(undef, B)
    ls = Vector{Float64}(undef, B)

    for _ in 1:cld(prop.maxiter, B)
        # Drawn here, serially: an AbstractRNG is not thread-safe, and drawing in order is
        # what keeps the accepted point the one `Rejection` would have picked.
        for i in 1:B
            us[i] = rand(rng, bounds)
        end
        keep = findall(Proposals.unitcheck, us)
        isempty(keep) && continue

        # Dynamic scheduling, deliberately. `:static` would pin the tasks and make
        # `threadid()` stable for free, but it throws outright when it finds itself inside
        # another threaded region -- so a caller who fits several models in parallel would get
        # an error instead of a fit. The per-thread slots in `_nested_pieces` are safe without
        # it: a task can only migrate at a yield point, and there is none inside the
        # likelihood between reading `threadid()` and writing its slot.
        Threads.@threads for k in eachindex(keep)
            @inbounds vs[k], ls[k] = prior_transform_and_loglikelihood(model, us[keep[k]])
        end
        ncall += length(keep)

        for k in eachindex(keep)                     # in draw order
            @inbounds ls[k] ≥ logl_star && return us[keep[k]], vs[k], ls[k], ncall
        end
    end
    throw(ErrorException("ThreadedRejection could not find a point above logl_star = " *
                         "$logl_star after $ncall likelihood calls; bounds = $bounds"))
end

"""
    _resample_equal(samples, weights, rng) -> Matrix

Systematic resampling of weighted nested samples to equal weight, so `NestedResult.posterior`
means the same thing whichever backend produced it.

Systematic rather than multinomial: one uniform draw positions a comb of `n` equally spaced
points through the cumulative weights, which has strictly lower variance than `n` independent
draws and cannot, as multinomial resampling can, drop the highest-weight sample entirely.
"""
function _resample_equal(samples::AbstractMatrix{<:Real}, weights::AbstractVector{<:Real},
                         rng::AbstractRNG)
    n = size(samples, 1)
    n == length(weights) ||
        error("resampling needs one weight per sample; got $(length(weights)) for $n samples")
    total = sum(weights)
    total > 0 || error("nested sampling returned zero total weight; the run did not converge")

    cw = cumsum(weights ./ total)
    cw[end] = 1.0                        # guard the last comb tooth against rounding
    idx = Vector{Int}(undef, n)
    i = 1
    offset = rand(rng)
    for j in 1:n
        u = (j - 1 + offset) / n
        while i < n && u > cw[i]
            i += 1
        end
        idx[j] = i
    end
    return Matrix{Float64}(samples[idx, :])
end

# The prior transform, the log-likelihood and the constraint handling are shared by both
# methods below and are identical to the UltraNest backend's, so that the two samplers are
# demonstrably fitting the same posterior and any disagreement is the sampler's.
#
# The likelihood also keeps the best point it has seen. The sampler reports the evidence and
# does not report a maximum-likelihood point, and recovering one afterwards by re-evaluating
# the retained samples costs a measured 20% of the run for a worse answer: it can only find
# the best of the few thousand samples KEPT, while this sees all forty-odd thousand
# evaluations, which is the set UltraNest's `maximum_likelihood` is drawn from too.
function _nested_pieces(fm, data, weights, vonmises, cons, constraint_idx, lbounds, ubounds)
    Δx = ubounds .- lbounds
    prior_transform = u -> u .* Δx .+ lbounds

    # One slot per thread rather than one shared best, because `ThreadedRejection` calls this
    # from several threads at once and a shared `Ref` would be a data race. `threadid()` is a
    # safe index here for the usual reason: a task can only migrate at a yield point, and
    # there is none between reading the id and writing the slot. Two tasks sharing a thread
    # cannot be inside the closure at once for the same reason.
    # `maxthreadid()`, not `nthreads()`: the latter counts only the `:default` pool, while
    # thread ids are numbered across every pool, so a task running on an interactive thread
    # has an id above `nthreads()` and would index off the end.
    nt        = Threads.maxthreadid()
    best_logl = fill(-Inf, nt)
    best_x    = [Float64[] for _ in 1:nt]

    loglikelihood = function (param)
        ll = -0.5 * (chi2_flat(fm, param, data; weights, vonmises) +
                     constraint_penalty(param, cons, fm, constraint_idx))
        t = Threads.threadid()
        # `>` and not `>=`, so a NaN likelihood never becomes the best point.
        @inbounds if ll > best_logl[t]
            best_logl[t] = ll
            best_x[t]    = Vector{Float64}(param)  # the transform hands out a fresh vector,
        end                                        # but copying makes that not matter
        return ll
    end
    return prior_transform, loglikelihood, best_logl, best_x
end

"""
    _pick_proposal(proposal, batch, ndims) -> proposal

Substitute [`ThreadedRejection`](@ref) for the proposal `Nested` would have chosen by itself,
and only there.

`:auto` picks `Proposals.Rejection` below 10 parameters and a Markov-chain proposal above it;
only the first is a rejection sampler and only the first can be batched. So this stays out of
the way in every other case, and an explicitly chosen proposal is always honoured — a caller
who asked for `RWalk` gets `RWalk`.
"""
function _pick_proposal(proposal, batch::Integer, ndims::Integer)
    proposal === :auto || return proposal
    (ndims < 10 && batch > 1 && Threads.nthreads() > 1) || return :auto
    return ThreadedRejection(batch)
end

function _nestedsamplers_run(fm, list_free_params, data, lbounds, ubounds;
                             weights, vonmises, cons, constraint_idx,
                             nactive, bounds, proposal, batch, dlogz, maxiter, maxcall,
                             rng, verb, cornerplot)

    prior_transform, loglikelihood, best_logl, best_x =
        _nested_pieces(fm, data, weights, vonmises, cons, constraint_idx, lbounds, ubounds)

    ndims  = length(list_free_params)
    prop   = _pick_proposal(proposal, batch, ndims)
    model  = NestedModel(loglikelihood, prior_transform)
    smplr  = Nested(ndims, nactive; bounds = bounds, proposal = prop)

    # `chain_type = Array` rather than `Chains`: the result is a plain
    # (n_samples, n_params + 1) matrix whose last column is the normalised weight, which is
    # exactly what is needed and costs no MCMCChains conversion.
    vals, state = NestedSamplers.StatsBase.sample(rng, model, smplr;
                                                  dlogz = dlogz, maxiter = maxiter,
                                                  maxcall = maxcall, chain_type = Array,
                                                  progress = verb)
    verb && println()                    # the progress line ends without one

    params = @view vals[:, 1:end-1]
    wts    = @view vals[:, end]

    # The best over every thread's slot.
    minx = best_x[argmax(best_logl)]
    isempty(minx) && error("nested sampling made no likelihood evaluation; the run failed")

    # The data chi2, without the constraint penalty that shaped the sampling — the convention
    # every other fitter reports.
    minf = chi2_flat(fm, minx, data; weights, vonmises)
    ndof = _ndof(data, weights)

    if verb
        # `state.it`, not `state.ncall`: the state handed back by `bundle_samples` has already
        # absorbed the live points and drops the call counter.
        @printf("logZ = %.4f ± %.4f   (%d samples, %d iterations, %s)\n",
                state.logz, state.logzerr, size(params, 1), state.it, string(prop))
        @printf("Best-fit χ² = %.4f   (χ²r = %.4f, ndof = %d)\n", minf, minf / ndof, ndof)
    end

    posterior = _resample_equal(params, wts, rng)

    if cornerplot
        if isempty(methods(plot_corner_makie))
            @warn "cornerplot=true, but drawing one needs PairPlots: `using PairPlots` " *
                  "enables it. The fit itself is unaffected; the " *
                  "result carries `posterior`, so the corner plot can be drawn later." maxlog = 1
        else
            plot_corner_makie(posterior, list_free_params)
        end
    end

    return NestedResult(minx, list_free_params, minf, minf / ndof, ndof,
                        state.logz, state.logzerr, posterior, state, fm, :nestedsamplers)
end

function _nestedsamplers_fit(model_dict        ::Dict{String},
                             list_free_params  ::Vector{String},
                             data              ::OIdata;
    lb            = Dict{String,Float64}(),
    ub            = Dict{String,Float64}(),
    weights       = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    constraints   = [],
    vonmises      = false,
    nB_workspace  = nothing,
    nactive       = 400,
    bounds        = Bounds.MultiEllipsoid,
    proposal      = :auto,
    batch         = 8,
    dlogz         = 0.2,
    maxiter       = Inf,
    maxcall       = Inf,
    rng           = Random.default_rng(),
    verb          = true,
    cornerplot    = true,
)
    for p in list_free_params
        haskey(lb, p) || error("Lower bound required for '$p' (nested sampling needs finite bounds)")
        haskey(ub, p) || error("Upper bound required for '$p' (nested sampling needs finite bounds)")
    end
    lbounds = Float64[lb[p] for p in list_free_params]
    ubounds = Float64[ub[p] for p in list_free_params]

    nB              = something(nB_workspace, size(data.uv, 2))
    cons            = parse_constraints(constraints)
    augmented       = copy(model_dict)
    constraint_keys = augment_constraints!(augmented, cons)
    fm              = parse_model(augmented, list_free_params; nB_workspace = nB)
    constraint_idx  = [findfirst(==(k), fm.all_names) for k in constraint_keys]

    return _nestedsamplers_run(fm, list_free_params, data, lbounds, ubounds;
                               weights, vonmises, cons, constraint_idx,
                               nactive, bounds, proposal, batch, dlogz, maxiter, maxcall,
                               rng, verb, cornerplot)
end

# Same fit against a pre-compiled `FlatModel`. Bounds may be Dicts or Vectors here, since a
# FlatModel already fixes the parameter order.
function _nestedsamplers_fit(model             ::FlatModel,
                             data              ::OIdata;
    lb            = Dict{String,Float64}(),
    ub            = Dict{String,Float64}(),
    weights       = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    constraints   = [],
    vonmises      = false,
    nactive       = 400,
    bounds        = Bounds.MultiEllipsoid,
    proposal      = :auto,
    batch         = 8,
    dlogz         = 0.2,
    maxiter       = Inf,
    maxcall       = Inf,
    rng           = Random.default_rng(),
    verb          = true,
    cornerplot    = true,
)
    list_free_params = model.list_free_params
    !isempty(constraints) && error("Constraints are only supported with the Dict form of " *
                                   "fit_model_nested; pass model_dict instead of FlatModel")

    lbounds = _bounds_vec(lb, list_free_params, NaN)
    ubounds = _bounds_vec(ub, list_free_params, NaN)
    for (i, p) in enumerate(list_free_params)
        isnan(lbounds[i]) && error("Lower bound required for '$p' (nested sampling needs finite bounds)")
        isnan(ubounds[i]) && error("Upper bound required for '$p' (nested sampling needs finite bounds)")
    end

    return _nestedsamplers_run(model, list_free_params, data, lbounds, ubounds;
                               weights, vonmises, cons = ModelConstraint[], constraint_idx = Int[],
                               nactive, bounds, proposal, batch, dlogz, maxiter, maxcall,
                               rng, verb, cornerplot)
end
