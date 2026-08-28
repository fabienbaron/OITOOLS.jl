# ─────────────────────────────────────────────────────────────────────────────
# The Nautilus.jl backend of `fit_model_nested`
# ─────────────────────────────────────────────────────────────────────────────
#
# Included by ext/OITOOLSNautilusExt.jl. Nothing here is public: `fit_model_nested` is declared
# in the core package and documents the keywords.
#
# Nautilus.jl is a Julia port of the `nautilus` Python package (Lange 2023, MNRAS 525, 3181):
# importance nested sampling with neural-network boosted space exploration. It is the third
# backend, and it exists because it is strictly better placed than the other two for this job:
#
#   * Pure Julia, so unlike `:ultranest` it needs no PythonCall and no conda environment, and
#     it survives into a PackageCompiler build.
#   * Unlike the Python package it accepts a single free parameter, so a one-parameter diameter
#     fit -- the shape test/test_nested.jl uses -- works here. Quadrature is still the better
#     tool for a 1-D integral.
#   * Importance nested sampling reuses every likelihood evaluation it ever makes instead of
#     discarding dead points, so it buys far more evidence precision per χ² call -- and the χ²
#     calls are the entire cost of an OIFITS model fit.
#
# MEASURED, on the α Cen A power-law limb-darkening fit (2 free parameters, 10 runs each,
# 4 threads, one process). All three backends agree on the best fit to four decimals --
# θ_LD = 8.4566 mas, α = 0.1358, with run-to-run spread of 1e-5 on both. They differ in what a
# run costs, what precision it buys, and how reliably it finishes:
#
#                    calls    wall     logZ        scatter  quoted   samples  failed
#   :nautilus        12240    3.72 s   -712.6875   0.0092   0.0192     1264    0/10
#   :ultranest       10993    2.19 s   -712.5819   0.1145   0.2892     7902    0/10
#   :nestedsamplers   5875    9.24 s   -712.6803   0.1842   0.1659     5875    1/10
#
# At matched call counts nautilus is about 12x tighter on the evidence than UltraNest, for
# 1.7x the wall clock -- the gap is per-call overhead (four networks trained per bound), which
# is only visible because this χ² is cheap at 143 µs. Both quote conservative errors: their
# stated uncertainty is roughly twice their actual scatter.
#
# Three things worth knowing before reading too much into this:
#
#   * The posterior is 2-D and unimodal, which is exactly where a plain ellipsoidal nested
#     sampler is already near-optimal and the neural bounds have nothing to add. The advantage
#     claimed in the paper is for higher dimensions and harder geometries; this fit does not
#     test that.
#   * UltraNest's mean sits ~0.11 above the other two, which agree with each other to 0.007.
#     Against UltraNest's own scatter that is a ~3σ offset. Not conclusive from ten runs, but
#     it is the kind of thing keeping three independent backends exists to surface.
#   * :nautilus returns the FEWEST posterior samples here (1264 against UltraNest's 7902),
#     because equal-weight resampling at `equal_weight_boost = 1` cannot exceed the effective
#     sample size. Raise `n_eff` for a thicker corner plot; raising `equal_weight_boost` alone
#     only duplicates points it already has.
#
# A SECOND FIT, harder geometry: uniform disc + skewed ring (the `crescent` component) on the
# 2004 Beauty Contest data, 6 free parameters, 5 runs each. This is the regime the first fit
# could not test, and the ordering inverts:
#
#                    calls    wall     logZ         scatter  quoted   samples
#   :nautilus        36580   21.2 s   -19275.593    0.011    0.015      1936
#   :ultranest      235323   13.9 s   -19275.403    0.301    0.462     20320
#
# 6.4x FEWER likelihood calls, and 27x tighter on the evidence. Both recover all six
# parameters to three or four decimals and their evidences agree to 0.19, inside UltraNest's
# own scatter. This is the advantage the paper claims, and it does not appear at all in the
# 2-D fit above.
#
# Wall clock still favours UltraNest here, and the reason is worth understanding rather than
# generalising: χ² on this dataset costs 23.6 µs (247 uv points), so it accounts for under a
# second of nautilus's 21 s. The rest is fixed overhead -- four networks trained per bound,
# which grows with dimension, not with data size. Equating the two runtimes puts the crossover
# at a likelihood cost of roughly 60 µs, i.e. about 600 uv points on this model. Above that,
# the 6.4x call advantage starts showing up in wall clock too; a dataset of a few thousand
# visibilities, or any model with a Hankel profile, is comfortably past it.
#
# The practical rule: choose :nautilus when the likelihood is expensive or the evidence needs
# to be precise, :ultranest when the likelihood is cheap and a rough evidence will do.
#
# :nautilus also returns the fewest posterior samples in both fits. Equal-weight resampling
# cannot exceed the effective sample size, so `n_eff` is the knob for a denser corner plot.
#
# :nestedsamplers failed outright on 1 run in 10 -- `ThreadedRejection` exhausted its budget
# without finding a point above the likelihood threshold. That is a pre-existing property of
# that backend, not of this one, and its 9.24 s mean reflects the same instability its own
# source comments describe.
#
# Structured exactly like the other two backends, so the three stay comparable.
#
# TWO DIFFERENCES FROM THE OTHER BACKENDS, both absorbed here:
#
#   * Nautilus reports no evidence ERROR of its own. The package's own FAQ gives
#     σ(log Z) ≈ 1/sqrt(N_eff), and that is what `logzerr` carries -- a different quantity from
#     UltraNest's bootstrap error, so do not read a disagreement between them as a bug.
#   * There is no maximum-likelihood point in the result either, so the best point is tracked
#     during the run, per thread, exactly as `_nested_pieces` does for NestedSamplers.

"""
    _nautilus_pieces(fm, data, weights, vonmises, cons, constraint_idx, lbounds, ubounds)

The prior transform and log-likelihood closures, plus the per-thread slots that track the best
point seen.

One slot per thread rather than a shared best, because with `threaded = true` the likelihood is
called from several threads at once and a shared `Ref` would be a data race. `threadid()` is a
safe index here for the usual reason: a task can only migrate at a yield point, and there is
none between reading the id and writing the slot. `maxthreadid()` rather than `nthreads()`,
since thread ids are numbered across every pool.
"""
function _nautilus_pieces(fm, data, weights, vonmises, cons, constraint_idx, lbounds, ubounds)
    Δx = ubounds .- lbounds
    prior_transform = u -> u .* Δx .+ lbounds

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
            best_x[t]    = Vector{Float64}(param)
        end
        return ll
    end
    return prior_transform, loglikelihood, best_logl, best_x
end

function _nautilus_run(fm, list_free_params, data, lbounds, ubounds;
                       weights, vonmises, cons, constraint_idx,
                       n_live, n_eff, f_live, n_networks, discard_exploration,
                       n_like_max, timeout, equal_weight_boost, threaded, seed,
                       verb, cornerplot)

    ndims = length(list_free_params)

    prior_transform, loglikelihood, best_logl, best_x =
        _nautilus_pieces(fm, data, weights, vonmises, cons, constraint_idx, lbounds, ubounds)

    smplr = Nautilus.Sampler(prior_transform, loglikelihood;
                             n_dim = ndims, n_live = n_live, n_networks = n_networks,
                             threaded = threaded, seed = seed)

    finished = Nautilus.run!(smplr; f_live = f_live, n_eff = n_eff,
                             n_like_max = n_like_max, timeout = timeout,
                             discard_exploration = discard_exploration, verbose = verb)
    verb && !finished &&
        @warn "nautilus stopped on n_like_max or timeout before converging; the evidence " *
              "and posterior are usable but less accurate than requested."

    logz = Nautilus.log_z(smplr)
    neff = Nautilus.n_eff(smplr)
    # Nautilus has no built-in evidence error; see the header.
    logzerr = neff > 0 ? 1 / sqrt(neff) : NaN

    minx = best_x[argmax(best_logl)]
    isempty(minx) && error("nautilus made no likelihood evaluation; the run failed")

    # The data chi2, without the constraint penalty that shaped the sampling — the convention
    # every other fitter reports.
    minf = chi2_flat(fm, minx, data; weights, vonmises)
    ndof = _ndof(data, weights)

    # Equal-weight samples, so `NestedResult.posterior` means the same thing whichever backend
    # produced it — and so PairPlots, which takes no per-sample weights, can draw it directly.
    posterior, _, _ = Nautilus.posterior(smplr; equal_weight = true,
                                         equal_weight_boost = equal_weight_boost)

    if verb
        @printf("logZ = %.4f ± %.4f   (%d samples, %d calls, N_eff = %.0f, η = %.3f)\n",
                logz, logzerr, size(posterior, 1), smplr.n_like, neff, Nautilus.eta(smplr))
        @printf("Best-fit χ² = %.4f   (χ²r = %.4f, ndof = %d)\n", minf, minf / ndof, ndof)
    end

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
                        logz, logzerr, posterior, smplr, fm, :nautilus)
end

function _nautilus_fit(model_dict        ::Dict{String},
                       list_free_params  ::Vector{String},
                       data              ::OIdata;
    lb                  = Dict{String,Float64}(),
    ub                  = Dict{String,Float64}(),
    weights             = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    constraints         = [],
    vonmises            = false,
    nB_workspace        = nothing,
    n_live              = 500,
    n_eff               = 2000,
    f_live              = 0.01,
    n_networks          = 4,
    discard_exploration = false,
    n_like_max          = Inf,
    timeout             = Inf,
    equal_weight_boost  = 1.0,
    threaded            = Threads.nthreads() > 1,
    seed                = nothing,
    verb                = true,
    cornerplot          = true,
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

    return _nautilus_run(fm, list_free_params, data, lbounds, ubounds;
                         weights, vonmises, cons, constraint_idx,
                         n_live, n_eff, f_live, n_networks, discard_exploration,
                         n_like_max, timeout, equal_weight_boost, threaded, seed,
                         verb, cornerplot)
end

# Same fit against a pre-compiled `FlatModel`. Bounds may be Dicts or Vectors here, since a
# FlatModel already fixes the parameter order.
function _nautilus_fit(model             ::FlatModel,
                       data              ::OIdata;
    lb                  = Dict{String,Float64}(),
    ub                  = Dict{String,Float64}(),
    weights             = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    constraints         = [],
    vonmises            = false,
    n_live              = 500,
    n_eff               = 2000,
    f_live              = 0.01,
    n_networks          = 4,
    discard_exploration = false,
    n_like_max          = Inf,
    timeout             = Inf,
    equal_weight_boost  = 1.0,
    threaded            = Threads.nthreads() > 1,
    seed                = nothing,
    verb                = true,
    cornerplot          = true,
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

    return _nautilus_run(model, list_free_params, data, lbounds, ubounds;
                         weights, vonmises, cons = ModelConstraint[], constraint_idx = Int[],
                         n_live, n_eff, f_live, n_networks, discard_exploration,
                         n_like_max, timeout, equal_weight_boost, threaded, seed,
                         verb, cornerplot)
end
