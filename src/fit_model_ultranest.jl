# ─────────────────────────────────────────────────────────────────────────────
# The UltraNest backend of `fit_model_nested`
# ─────────────────────────────────────────────────────────────────────────────
#
# Included by ext/OITOOLSUltraNestExt.jl, which is what supplies PythonCall. Nothing here is
# public: `fit_model_nested` and `fit_model_ultranest` are declared in the core package and
# document the keywords; this file is only how UltraNest is driven.
#
# Needs the Python `ultranest` package, declared in CondaPkg.toml.

function _ultranest_fit(model_dict         ::Dict{String},
                        list_free_params  ::Vector{String},
                        data              ::OIdata;
    lb                       = Dict{String,Float64}(),
    ub                       = Dict{String,Float64}(),
    weights                  = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    constraints              = [],
    vonmises                 = false,
    nB_workspace             = nothing,
    min_num_live_points      = 400,
    cluster_num_live_points  = 100,
    num_bootstraps           = 30,
    use_stepsampler          = false,
    nsteps                   = 400,
    frac_remain              = 0.001,
    log_interval             = 100,
    verb                     = true,
    cornerplot               = true,
)
    # ── 1. Validate bounds ──────────────────────────────────────────────────
    for p in list_free_params
        haskey(lb, p) || error("Lower bound required for '$p' (UltraNest needs finite bounds)")
        haskey(ub, p) || error("Upper bound required for '$p' (UltraNest needs finite bounds)")
    end
    lbounds = Float64[lb[p] for p in list_free_params]
    ubounds = Float64[ub[p] for p in list_free_params]

    # ── 2. Compile the model (once), with the constraint expressions ────────
    nB = something(nB_workspace, size(data.uv, 2))
    cons      = parse_constraints(constraints)
    augmented = copy(model_dict)
    constraint_keys = augment_constraints!(augmented, cons)
    fm = parse_model(augmented, list_free_params; nB_workspace=nB)
    constraint_idx  = [findfirst(==(k), fm.all_names) for k in constraint_keys]

    # ── 3. Prior transform: uniform [0,1]^n → [lb, ub] ─────────────────────
    Δx = ubounds .- lbounds
    function prior_transform(u::AbstractVector{<:Real})
        u .* Δx .+ lbounds
    end

    # UltraNest calls this with a numpy array, which arrives as a zero-copy `PyArray` and so
    # already satisfies `AbstractMatrix` — no conversion, no copy.
    #
    # The RETURN value does need help: a bare Julia array reaches Python as a
    # `juliacall.VectorValue`, and UltraNest calls numpy methods (`.transpose`) on it.
    # `.to_numpy()` hands back a real ndarray.
    prior_transform_vectorized = let trafo = prior_transform
        (U::AbstractMatrix{<:Real}) ->
            Py(reduce(vcat, (u -> trafo(u)').(eachrow(U)))).to_numpy()
    end

    # ── 4. Log-likelihood: -0.5 * (chi2 + constraint penalty) ──────────────
    # A soft constraint states which parameter combinations are admissible, so it belongs in
    # the density being sampled; the posterior and the evidence are then both the constrained
    # ones.
    loglikelihood = let fm=fm, data=data, weights=weights, vonmises=vonmises,
                        cons=cons, constraint_idx=constraint_idx
        param::AbstractVector{<:Real} ->
            -0.5 * (chi2_flat(fm, param, data; weights, vonmises) +
                    constraint_penalty(param, cons, fm, constraint_idx))
    end

    loglikelihood_vectorized = let loglikelihood = loglikelihood
        (X::AbstractMatrix{<:Real}) -> Py(_batch_loglike(loglikelihood, X)).to_numpy()
    end

    # ── 5. Run UltraNest ───────────────────────────────────────────────────
    ultranest = pyimport("ultranest")

    # UltraNest concatenates this with a Python list (`names + [...]`), so it must be a real
    # list: a Julia Vector arrives as a juliacall.VectorValue, which does not support `+`.
    param_names = pylist(list_free_params)

    # `verb` has to reach UltraNest itself, not just our own printing. Raising
    # `log_interval` alone leaves `show_status` and `viz_callback` at their defaults, so
    # the sampler keeps drawing its live status bar and "Mono-modal Volume" blocks —
    # thousands of lines from a fit the caller asked to be quiet, which is unusable in a
    # loop over models and drowns any surrounding output in a script. Its `ultranest`
    # Python logger is a third, independent channel ("Sampling N live points...", the
    # logZ error budget), silenced by raising its level for the duration of the run.
    if !verb
        log_interval = 1_000_000
    end
    show_status  = verb
    viz_callback = verb ? "auto" : false

    _logging  = pyimport("logging")
    _un_log   = _logging.getLogger("ultranest")
    _un_level = _un_log.level                      # restore after: the logger is global

    smplr = ultranest.ReactiveNestedSampler(
        param_names, loglikelihood_vectorized;
        transform       = prior_transform_vectorized,
        num_bootstraps  = num_bootstraps,
        vectorized      = true,
    )

    if use_stepsampler
        stepsampler_mod = pyimport("ultranest.stepsampler")
        smplr.stepsampler = stepsampler_mod.RegionSliceSampler(
            nsteps         = nsteps,
            adaptive_nsteps = "move-distance",
        )
    end

    verb || _un_log.setLevel(_logging.WARNING)
    result = try
        smplr.run(;
            min_num_live_points     = min_num_live_points,
            cluster_num_live_points = cluster_num_live_points,
            log_interval            = log_interval,
            frac_remain             = frac_remain,
            show_status             = show_status,
            viz_callback            = viz_callback,
        )
    finally
        verb || _un_log.setLevel(_un_level)
    end

    # ── 6. Extract results ─────────────────────────────────────────────────
    minx = pyconvert(Vector{Float64}, result["maximum_likelihood"]["point"])
    # The data chi2, without the constraint penalty that shaped the sampling — the convention
    # the other fitters report.
    minf = chi2_flat(fm, minx, data; weights, vonmises)

    if verb
        @printf("Best-fit χ² = %.4f\n", minf)
        smplr.print_results()
    end

    # ── 7. Corner plot ─────────────────────────────────────────────────────
    if cornerplot
        if isempty(methods(plot_ultranest_corner))
            @warn "cornerplot=true, but drawing needs matplotlib: `using PythonPlot` " *
                  "enables it. The fit itself is unaffected; the result carries `posterior`, " *
                  "so the corner plot can be drawn later." maxlog = 1
        else
            plot_ultranest_corner(result)
        end
    end

    # ── 8. Build posterior matrix ──────────────────────────────────────────
    samples_py = result["samples"]
    posterior  = pyconvert(Matrix{Float64}, samples_py)

    logz    = pyconvert(Float64, result["logz"])
    logzerr = pyconvert(Float64, result["logzerr"])
    ndof    = _ndof(data, weights)

    return NestedResult(
        minx, list_free_params, minf, minf / ndof, ndof,
        logz, logzerr, posterior, result, fm, :ultranest,
    )
end

# Same fit against a pre-compiled `FlatModel`. Bounds may be Dicts or Vectors here, since a
# FlatModel already fixes the parameter order.
function _ultranest_fit(model             ::FlatModel,
                        data              ::OIdata;
    lb                       = Dict{String,Float64}(),
    ub                       = Dict{String,Float64}(),
    weights                  = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    constraints              = [],
    vonmises                 = false,
    min_num_live_points      = 400,
    cluster_num_live_points  = 100,
    num_bootstraps           = 30,
    use_stepsampler          = false,
    nsteps                   = 400,
    frac_remain              = 0.001,
    log_interval             = 100,
    verb                     = true,
    cornerplot               = true,
)
    list_free_params = model.list_free_params
    fm = model

    # A constraint is evaluated as an extra parameter of the compiled model, so it needs the
    # model dict to augment and recompile. Accepting it here and dropping it would sample a
    # posterior that ignores what the caller asked for.
    !isempty(constraints) && error("Constraints are only supported with the Dict form of fit_model_ultranest; pass model_dict instead of FlatModel")

    # ── Validate bounds ─────────────────────────────────────────────────────
    lbounds = _bounds_vec(lb, list_free_params, NaN)
    ubounds = _bounds_vec(ub, list_free_params, NaN)
    for (i, p) in enumerate(list_free_params)
        isnan(lbounds[i]) && error("Lower bound required for '$p' (UltraNest needs finite bounds)")
        isnan(ubounds[i]) && error("Upper bound required for '$p' (UltraNest needs finite bounds)")
    end

    # ── Prior transform: uniform [0,1]^n → [lb, ub] ────────────────────────
    Δx = ubounds .- lbounds
    function prior_transform(u::AbstractVector{<:Real})
        u .* Δx .+ lbounds
    end

    # UltraNest calls this with a numpy array, which arrives as a zero-copy `PyArray` and so
    # already satisfies `AbstractMatrix` — no conversion, no copy.
    #
    # The RETURN value does need help: a bare Julia array reaches Python as a
    # `juliacall.VectorValue`, and UltraNest calls numpy methods (`.transpose`) on it.
    # `.to_numpy()` hands back a real ndarray.
    prior_transform_vectorized = let trafo = prior_transform
        (U::AbstractMatrix{<:Real}) ->
            Py(reduce(vcat, (u -> trafo(u)').(eachrow(U)))).to_numpy()
    end

    # ── Log-likelihood: -0.5 * chi2 ─────────────────────────────────────────
    loglikelihood = let fm=fm, data=data, weights=weights, vonmises=vonmises
        param::AbstractVector{<:Real} -> -0.5 * chi2_flat(fm, param, data;
                                                           weights, vonmises)
    end

    loglikelihood_vectorized = let loglikelihood = loglikelihood
        (X::AbstractMatrix{<:Real}) -> Py(_batch_loglike(loglikelihood, X)).to_numpy()
    end

    # ── Run UltraNest ────────────────────────────────────────────────────────
    ultranest = pyimport("ultranest")

    # `verb` has to reach UltraNest itself, not just our own printing. Raising
    # `log_interval` alone leaves `show_status` and `viz_callback` at their defaults, so
    # the sampler keeps drawing its live status bar and "Mono-modal Volume" blocks —
    # thousands of lines from a fit the caller asked to be quiet, which is unusable in a
    # loop over models and drowns any surrounding output in a script. Its `ultranest`
    # Python logger is a third, independent channel ("Sampling N live points...", the
    # logZ error budget), silenced by raising its level for the duration of the run.
    if !verb
        log_interval = 1_000_000
    end
    show_status  = verb
    viz_callback = verb ? "auto" : false

    _logging  = pyimport("logging")
    _un_log   = _logging.getLogger("ultranest")
    _un_level = _un_log.level                      # restore after: the logger is global

    smplr = ultranest.ReactiveNestedSampler(
        list_free_params, loglikelihood_vectorized;
        transform       = prior_transform_vectorized,
        num_bootstraps  = num_bootstraps,
        vectorized      = true,
    )

    if use_stepsampler
        stepsampler_mod = pyimport("ultranest.stepsampler")
        smplr.stepsampler = stepsampler_mod.RegionSliceSampler(
            nsteps         = nsteps,
            adaptive_nsteps = "move-distance",
        )
    end

    verb || _un_log.setLevel(_logging.WARNING)
    result = try
        smplr.run(;
            min_num_live_points     = min_num_live_points,
            cluster_num_live_points = cluster_num_live_points,
            log_interval            = log_interval,
            frac_remain             = frac_remain,
            show_status             = show_status,
            viz_callback            = viz_callback,
        )
    finally
        verb || _un_log.setLevel(_un_level)
    end

    # ── Extract results ──────────────────────────────────────────────────────
    minx = pyconvert(Vector{Float64}, result["maximum_likelihood"]["point"])
    minf = chi2_flat(fm, minx, data; weights, vonmises)

    if verb
        @printf("Best-fit χ² = %.4f\n", minf)
        smplr.print_results()
    end

    if cornerplot
        if isempty(methods(plot_ultranest_corner))
            @warn "cornerplot=true, but drawing needs matplotlib: `using PythonPlot` " *
                  "enables it. The fit itself is unaffected; the result carries `posterior`, " *
                  "so the corner plot can be drawn later." maxlog = 1
        else
            plot_ultranest_corner(result)
        end
    end

    samples_py = result["samples"]
    posterior  = pyconvert(Matrix{Float64}, samples_py)
    logz    = pyconvert(Float64, result["logz"])
    logzerr = pyconvert(Float64, result["logzerr"])
    ndof    = _ndof(data, weights)

    return NestedResult(
        minx, list_free_params, minf, minf / ndof, ndof,
        logz, logzerr, posterior, result, fm, :ultranest,
    )
end
