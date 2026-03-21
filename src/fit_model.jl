# fit_model.jl
#
# Flat-dict parametric model fitting via NLopt.
# Replacement for the legacy modelfit.jl stack (OImodel, OIcomponent, etc.).
#
# Load order (within OITOOLS module):
#   include("resolvers.jl")          # SharedUtils, RGF
#   include("hankel.jl")             # Hankel + chainrules
#   include("model_chainrules.jl")   # vis_ud, vis_ldlin, ...
#   include("readoifits.jl")         # OIdata struct
#   include("parse_model.jl")        # FlatModel, eval_model, eval_model_grad
#   include("chi2_flat.jl")          # chi2_flat, chi2_flat_fg
#   include("fit_model.jl")          # this file
#
# Interface
# ---------
#   fit_model(model_dict, list_free_params, data; lb, ub, weights, priors, method, ...)
#       → FitResult

using NLopt, LsqFit, ForwardDiff, Printf, PyCall, FFTW


# ─────────────────────────────────────────────────────────────────────────────
# display_model — inspect model setup before fitting
# ─────────────────────────────────────────────────────────────────────────────

"""
    display_model(model_dict, list_free_params; lb, ub)

Print a human-readable summary of the model setup, grouped by component.
Shows each parameter's value (or expression), whether it is free or fixed,
and its bounds.  Warns about common mistakes (value out of bounds, missing
bounds, flux fractions that don't sum to 1).

Call this *before* `fit_model` / `fit_model_ultranest` to verify your setup.
"""
function display_model(model_dict  ::Dict{String},
                       list_free_params  ::Vector{String};
                       lb = Dict{String,Float64}(),
                       ub = Dict{String,Float64}())

    fit_set = Set(list_free_params)

    # ── Discover components ─────────────────────────────────────────────────
    comp_params = Dict{String, Vector{String}}()   # comp_name → [keys...]
    for k in sort(collect(keys(model_dict)))
        idx = findfirst(',', k)
        if isnothing(idx)
            cn = "__global__"
        else
            cn = k[1:idx-1]
        end
        push!(get!(comp_params, cn, String[]), k)
    end

    # ── Identify component kinds ────────────────────────────────────────────
    n_warnings = 0

    for cn in sort(collect(keys(comp_params)))
        params = comp_params[cn]

        # Detect kind
        kind = "unknown"
        for gk in ["profile", "diamin", "fwhmin", "crin",
                   "ud", "fwhm", "ldlin", "ldquad", "ldpow", "resolved"]
            if "$cn,$gk" in params
                kind = gk == "profile" ? "hankel" :
                       gk == "fwhm"    ? "gaussian" :
                       gk == "diamin"  ? "ring" :
                       gk == "fwhmin"  ? "gaussian_ring" :
                       gk == "crin"    ? "crescent" : gk
                break
            end
        end
        if kind == "unknown" && "$cn,f" in params
            kind = "point"
        end

        # Header
        printstyled("\nComponent: $cn  [$kind]\n", bold=true)
        println("─" ^ 78)
        @printf("  %-24s  %-8s  %-30s  %s\n", "Parameter", "Status", "Value", "Bounds")
        println("  " * "─"^74)

        for k in params
            suffix = occursin(',', k) ? k[findfirst(',', k)+1:end] : k
            val    = model_dict[k]
            free   = k in fit_set

            # Format value column
            if val isa Number
                val_str = @sprintf("%.6g", val)
            elseif val isa String
                val_str = length(val) > 28 ? val[1:25] * "..." : val
            else
                val_str = string(val)
                if length(val_str) > 28
                    val_str = val_str[1:25] * "..."
                end
            end

            # Format status and bounds
            if free
                lo = get(lb, k, nothing)
                hi = get(ub, k, nothing)
                lo_str = isnothing(lo) ? "-∞" : @sprintf("%.4g", lo)
                hi_str = isnothing(hi) ? "+∞" : @sprintf("%.4g", hi)
                bounds_str = "[$lo_str, $hi_str]"
                status_str = "FREE"

                # Sanity checks
                if val isa Number
                    if !isnothing(lo) && val < lo
                        bounds_str *= "  ⚠ value < lb!"
                        n_warnings += 1
                    end
                    if !isnothing(hi) && val > hi
                        bounds_str *= "  ⚠ value > ub!"
                        n_warnings += 1
                    end
                    if !isnothing(lo) && !isnothing(hi) && lo >= hi
                        bounds_str *= "  ⚠ lb ≥ ub!"
                        n_warnings += 1
                    end
                end

                printstyled(@sprintf("  %-24s  %-8s  %-30s  %s\n",
                    suffix, status_str, val_str, bounds_str), color=:red)
            else
                printstyled(@sprintf("  %-24s  %-8s  %-30s\n",
                    suffix, "fixed", val_str), color=:default)
            end
        end
    end

    # ── Summary ─────────────────────────────────────────────────────────────
    println("\n" * "─"^78)
    n_free = length(list_free_params)
    n_total = length(model_dict)
    println("Total: $n_total parameters, $n_free free")
    println("Free:  $(join(list_free_params, ", "))")

    # Check flux sum at reference values
    flux_keys = [k for k in keys(model_dict) if endswith(k, ",f")]
    numeric_fluxes = [(k, model_dict[k]) for k in flux_keys if model_dict[k] isa Number]
    expr_fluxes    = [(k, model_dict[k]) for k in flux_keys if model_dict[k] isa String]
    if !isempty(numeric_fluxes) && isempty(expr_fluxes)
        fsum = sum(v for (_, v) in numeric_fluxes)
        if abs(fsum - 1.0) > 0.01
            @printf("⚠ Numeric flux fractions sum to %.4f (expected ≈ 1.0)\n", fsum)
            n_warnings += 1
        end
    elseif !isempty(expr_fluxes)
        println("  Note: some flux fractions are expressions — sum check skipped")
    end

    if n_warnings > 0
        printstyled("  $n_warnings warning(s) found — review before fitting\n", color=:yellow)
    else
        printstyled("  Setup looks good ✓\n", color=:green)
    end
    println()
    return nothing
end


# ─────────────────────────────────────────────────────────────────────────────
# FitResult
# ─────────────────────────────────────────────────────────────────────────────

struct FitResult
    x_opt      ::Vector{Float64}   # best-fit parameter values
    list_free_params ::Vector{String}    # parameter names, same order as x_opt
    chi2       ::Float64           # raw weighted chi2 at x_opt
    chi2r      ::Float64           # chi2 / ndof
    ndof       ::Int               # number of data points (weighted)
    n_evals    ::Int               # number of chi2 function evaluations
    ret        ::Symbol            # NLopt return code
    model      ::FlatModel         # the compiled model
end

function Base.show(io::IO, r::FitResult)
    @printf(io, "FitResult: χ²r = %.4f  (chi2=%.2f, ndof=%d, evals=%d, ret=%s)\n",
            r.chi2r, r.chi2, r.ndof, r.n_evals, r.ret)
    for (i, p) in enumerate(r.list_free_params)
        @printf(io, "  %-20s = %10.4f\n", p, r.x_opt[i])
    end
end


# ─────────────────────────────────────────────────────────────────────────────
# _dollarify: inject $ before "comp,param" references in an expression string
# ─────────────────────────────────────────────────────────────────────────────

function _dollarify(expr::String)
    replace(expr, r"\b([A-Za-z_][A-Za-z0-9_]*,[A-Za-z_][A-Za-z0-9_]*)\b" => s"$\1")
end


# ─────────────────────────────────────────────────────────────────────────────
# Prior helpers
# ─────────────────────────────────────────────────────────────────────────────

function prior_penalty_and_grad!(g, x, priors, model, prior_idx)
    isempty(priors) && return 0.0
    pen = 0.0
    for (i, (_, target, sigma)) in enumerate(priors)
        idx = prior_idx[i]
        val = model.resolver.fn(x)[idx]
        r   = (val - target) / sigma^2
        pen += (val - target)^2 / sigma^2
        dpv = ForwardDiff.jacobian(x_ -> [model.resolver.fn(x_)[idx]], x)
        g .+= 2 .* r .* vec(dpv)
    end
    return pen
end

function prior_penalty(x, priors, model, prior_idx)
    isempty(priors) && return 0.0
    sum((model.resolver.fn(x)[prior_idx[i]] - target)^2 / sigma^2
        for (i, (_, target, sigma)) in enumerate(priors))
end


# ─────────────────────────────────────────────────────────────────────────────
# Helper: convert Dict or Vector bounds to a Vector aligned with list_free_params
# ─────────────────────────────────────────────────────────────────────────────

function _bounds_vec(b, list_free_params, default)
    if b isa AbstractVector
        length(b) == length(list_free_params) || error("Bounds vector length ($(length(b))) must match list_free_params length ($(length(list_free_params)))")
        return Float64.(b)
    else
        return Float64[get(b, p, default) for p in list_free_params]
    end
end


# ─────────────────────────────────────────────────────────────────────────────
# fit_model — main entry point
# ─────────────────────────────────────────────────────────────────────────────

"""
    fit_model(model_dict, list_free_params, data; kwargs...) -> FitResult

Fit a parametric model to interferometric data using NLopt.

# Arguments
- `model_dict::Dict{String}` — flat parameter dictionary (values are numbers or expression strings)
- `list_free_params::Vector{String}` — names of the free parameters to optimize
- `data::OIdata` — interferometric data

# Keywords
- `lb`, `ub` — `Dict{String,Float64}` of lower/upper bounds per parameter (default: ±Inf)
- `weights` — observable weights [V2, T3amp, T3phi, visamp, visphi, flux, diffphase]
- `priors` — Vector of `(expr_str, target, sigma)` tuples for Gaussian penalties
- `method` — NLopt algorithm (default `:LD_LBFGS`; use `:LN_NELDERMEAD` for gradient-free)
- `maxeval` — maximum function evaluations (default 2000)
- `ftol_rel`, `xtol_rel` — convergence tolerances
- `vonmises` — use von Mises statistic for T3phi
- `nB_workspace` — Hankel workspace size (default: nuv from data)
- `verb` — print per-evaluation chi2 breakdown
"""
function fit_model(model_dict   ::Dict{String},
                   list_free_params   ::Vector{String},
                   data         ::OIdata;
    lb           = Dict{String,Float64}(),
    ub           = Dict{String,Float64}(),
    weights      = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    priors       = [],
    method       = :LD_LBFGS,
    maxeval      = 2000,
    ftol_rel     = 1e-8,
    xtol_rel     = 1e-6,
    vonmises     = false,
    nB_workspace = nothing,
    verb         = false,
)
    # ── 1. Augment model_dict with prior expressions ─────────────────────────
    augmented = copy(model_dict)
    prior_keys = String[]
    for (i, (expr, _, _)) in enumerate(priors)
        key = "__prior_$(i)__"
        augmented[key] = _dollarify(expr)
        push!(prior_keys, key)
    end

    # ── 2. Compile the model (once) ──────────────────────────────────────────
    nB = something(nB_workspace, size(data.uv, 2))
    fm = parse_model(augmented, list_free_params; nB_workspace=nB)

    # ── 3. Resolve prior indices ─────────────────────────────────────────────
    prior_idx = [findfirst(==(k), fm.all_names) for k in prior_keys]

    # ── 4. Detect gradient need ──────────────────────────────────────────────
    use_grad = startswith(string(method), "LD_")

    # ── 5. Build objective ───────────────────────────────────────────────────
    n_evals = Ref(0)
    function objective(x, g)
        n_evals[] += 1
        if use_grad
            chi2, grad = chi2_flat_fg(fm, x, data; weights, vonmises, verb)
            chi2 += prior_penalty_and_grad!(grad, x, priors, fm, prior_idx)
            g .= grad
        else
            chi2 = chi2_flat(fm, x, data; weights, vonmises, verb)
            chi2 += prior_penalty(x, priors, fm, prior_idx)
        end
        return chi2
    end

    # ── 6. Configure NLopt ───────────────────────────────────────────────────
    opt = Opt(method, length(list_free_params))
    min_objective!(opt, objective)
    ftol_rel!(opt, ftol_rel)
    xtol_rel!(opt, xtol_rel)
    maxeval!(opt, maxeval)
    lower_bounds!(opt, [get(lb, p, -Inf) for p in list_free_params])
    upper_bounds!(opt, [get(ub, p,  Inf) for p in list_free_params])

    # ── 7. Run ───────────────────────────────────────────────────────────────
    x0 = Float64[model_dict[p] for p in list_free_params]
    local minf, minx, ret
    try
        minf, minx, ret = optimize(opt, x0)
    catch e
        if e isa NLopt.ForcedStop
            minf = objective(opt.last_optimum_value isa Vector ? opt.last_optimum_value : x0, Float64[])
            minx = opt.last_optimum_value isa Vector ? opt.last_optimum_value : x0
            ret  = :MAXEVAL_REACHED
        else
            rethrow(e)
        end
    end

    # ── 8. Return result ─────────────────────────────────────────────────────
    ndof = _ndof(data, weights)
    return FitResult(minx, list_free_params, minf, minf / ndof, ndof, n_evals[], ret, fm)
end

"""
    fit_model(model::FlatModel, x0, data; kwargs...) -> FitResult

Fit a pre-compiled `FlatModel` starting from parameter vector `x0`.

This method skips the `dict_to_model` step, which is useful when fitting the
same model multiple times (e.g. bootstrap, grid search).

Bounds `lb` and `ub` can be `Dict{String,Float64}` (keyed by parameter name)
or `Vector{Float64}` (ordered to match `model.list_free_params`).
"""
function fit_model(model        ::FlatModel,
                   x0           ::AbstractVector{<:Real},
                   data         ::OIdata;
    lb           = Dict{String,Float64}(),
    ub           = Dict{String,Float64}(),
    weights      = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    priors       = [],
    method       = :LD_LBFGS,
    maxeval      = 2000,
    ftol_rel     = 1e-8,
    xtol_rel     = 1e-6,
    vonmises     = false,
    verb         = false,
)
    list_free_params = model.list_free_params
    fm = model

    # ── Prior support ────────────────────────────────────────────────────────
    # When priors are used, we would need to augment and recompile the model.
    # For now, only the Dict method supports priors.
    !isempty(priors) && error("Priors are only supported with the Dict form of fit_model; pass model_dict instead of FlatModel")

    # ── Detect gradient need ─────────────────────────────────────────────────
    use_grad = startswith(string(method), "LD_")

    # ── Build objective ──────────────────────────────────────────────────────
    n_evals = Ref(0)
    function objective(x, g)
        n_evals[] += 1
        if use_grad
            chi2, grad = chi2_flat_fg(fm, x, data; weights, vonmises, verb)
            g .= grad
        else
            chi2 = chi2_flat(fm, x, data; weights, vonmises, verb)
        end
        return chi2
    end

    # ── Configure NLopt ──────────────────────────────────────────────────────
    opt = Opt(method, length(list_free_params))
    min_objective!(opt, objective)
    ftol_rel!(opt, ftol_rel)
    xtol_rel!(opt, xtol_rel)
    maxeval!(opt, maxeval)
    lower_bounds!(opt, _bounds_vec(lb, list_free_params, -Inf))
    upper_bounds!(opt, _bounds_vec(ub, list_free_params,  Inf))

    # ── Run ──────────────────────────────────────────────────────────────────
    local minf, minx, ret
    try
        minf, minx, ret = optimize(opt, Float64.(x0))
    catch e
        if e isa NLopt.ForcedStop
            minf = objective(opt.last_optimum_value isa Vector ? opt.last_optimum_value : Float64.(x0), Float64[])
            minx = opt.last_optimum_value isa Vector ? opt.last_optimum_value : Float64.(x0)
            ret  = :MAXEVAL_REACHED
        else
            rethrow(e)
        end
    end

    ndof = _ndof(data, weights)
    return FitResult(minx, list_free_params, minf, minf / ndof, ndof, n_evals[], ret, fm)
end


# ─────────────────────────────────────────────────────────────────────────────
# UltraNest result
# ─────────────────────────────────────────────────────────────────────────────

struct UltraNestResult
    x_opt       ::Vector{Float64}   # maximum-likelihood parameter values
    list_free_params  ::Vector{String}    # parameter names, same order as x_opt
    chi2        ::Float64           # chi2 at x_opt
    chi2r       ::Float64           # chi2 / ndof
    ndof        ::Int
    logz        ::Float64           # log-evidence
    logzerr     ::Float64           # log-evidence uncertainty
    posterior   ::Matrix{Float64}   # (n_samples, n_params) posterior samples
    result      ::Any               # raw UltraNest result dict (PyObject)
    model       ::FlatModel
end

function Base.show(io::IO, r::UltraNestResult)
    @printf(io, "UltraNestResult: χ²r = %.4f  (chi2=%.2f, ndof=%d)\n",
            r.chi2r, r.chi2, r.ndof)
    @printf(io, "  log(Z) = %.2f ± %.2f\n", r.logz, r.logzerr)
    for (i, p) in enumerate(r.list_free_params)
        samples = r.posterior[:, i]
        med  = sort(samples)[length(samples) ÷ 2]
        lo   = sort(samples)[max(1, round(Int, 0.16 * length(samples)))]
        hi   = sort(samples)[min(length(samples), round(Int, 0.84 * length(samples)))]
        @printf(io, "  %-20s = %10.4f  (%.4f … %.4f)\n", p, r.x_opt[i], lo, hi)
    end
end


# ─────────────────────────────────────────────────────────────────────────────
# fit_model_ultranest — Bayesian nested sampling via UltraNest
# ─────────────────────────────────────────────────────────────────────────────

"""
    fit_model_ultranest(model_dict, list_free_params, data; kwargs...) -> UltraNestResult

Fit a parametric model to interferometric data using UltraNest nested sampling.

Requires PyCall and the Python `ultranest` package.

# Arguments
- `model_dict::Dict{String}` — flat parameter dictionary
- `list_free_params::Vector{String}` — names of the free parameters
- `data::OIdata` — interferometric data

# Keywords
- `lb`, `ub` — `Dict{String,Float64}` of lower/upper bounds (REQUIRED for all list_free_params)
- `weights` — observable weights [V2, T3amp, T3phi, visamp, visphi, flux, diffphase]
- `vonmises` — use von Mises statistic for T3phi
- `nB_workspace` — Hankel workspace size
- `min_num_live_points` — minimum number of live points (default 400)
- `cluster_num_live_points` — live points for clustering (default 100)
- `num_bootstraps` — number of bootstraps for evidence (default 30)
- `use_stepsampler` — use RegionSliceSampler (default false)
- `nsteps` — steps per slice for stepsampler (default 400)
- `frac_remain` — termination fraction (default 0.001)
- `log_interval` — logging interval (default 100)
- `verb` — print progress (default true)
- `cornerplot` — show corner plot (default true)
"""
function fit_model_ultranest(model_dict   ::Dict{String},
                             list_free_params   ::Vector{String},
                             data         ::OIdata;
    lb                       = Dict{String,Float64}(),
    ub                       = Dict{String,Float64}(),
    weights                  = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
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

    # ── 2. Compile the model (once) ─────────────────────────────────────────
    nB = something(nB_workspace, size(data.uv, 2))
    fm = parse_model(model_dict, list_free_params; nB_workspace=nB)

    # ── 3. Prior transform: uniform [0,1]^n → [lb, ub] ─────────────────────
    Δx = ubounds .- lbounds
    function prior_transform(u::AbstractVector{<:Real})
        u .* Δx .+ lbounds
    end

    prior_transform_vectorized = let trafo = prior_transform
        (U::AbstractMatrix{<:Real}) -> reduce(vcat, (u -> trafo(u)').(eachrow(U)))
    end

    # ── 4. Log-likelihood: -0.5 * chi2 ─────────────────────────────────────
    loglikelihood = let fm=fm, data=data, weights=weights, vonmises=vonmises
        param::AbstractVector{<:Real} -> -0.5 * chi2_flat(fm, param, data;
                                                           weights, vonmises)
    end

    loglikelihood_vectorized = let loglikelihood = loglikelihood
        (X::AbstractMatrix{<:Real}) -> loglikelihood.(eachrow(X))
    end

    # ── 5. Run UltraNest ───────────────────────────────────────────────────
    ultranest = pyimport("ultranest")

    param_names = list_free_params  # UltraNest wants a list of strings

    if !verb
        log_interval = 1_000_000
    end

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

    result = smplr.run(;
        min_num_live_points     = min_num_live_points,
        cluster_num_live_points = cluster_num_live_points,
        log_interval            = log_interval,
        frac_remain             = frac_remain,
    )

    # ── 6. Extract results ─────────────────────────────────────────────────
    minx = Float64.(result["maximum_likelihood"]["point"])
    minf = chi2_flat(fm, minx, data; weights, vonmises)

    if verb
        @printf("Best-fit χ² = %.4f\n", minf)
        smplr.print_results()
    end

    # ── 7. Corner plot ─────────────────────────────────────────────────────
    if cornerplot
        PyDict(pyimport("matplotlib")."rcParams")["font.size"] = [8]
        histogram_color = "black"
        contour_colors = py"{'colors': ['#0072B2','#56B4E9','#009E73','#F0E442'], 'linestyles': ['-', '-', '-', '-']}"
        pyimport("ultranest.plot").cornerplot(result;
            contour_kwargs = contour_colors,
            color          = histogram_color,
        )
    end

    # ── 8. Build posterior matrix ──────────────────────────────────────────
    samples_py = result["samples"]
    posterior  = convert(Matrix{Float64}, samples_py)

    logz    = Float64(result["logz"])
    logzerr = Float64(result["logzerr"])
    ndof    = _ndof(data, weights)

    return UltraNestResult(
        minx, list_free_params, minf, minf / ndof, ndof,
        logz, logzerr, posterior, result, fm,
    )
end

"""
    fit_model_ultranest(model::FlatModel, data; kwargs...) -> UltraNestResult

Bayesian nested sampling on a pre-compiled `FlatModel`.
Starting values are not needed (UltraNest samples the full prior volume).
Bounds `lb` and `ub` are required and can be Dicts or Vectors.
"""
function fit_model_ultranest(model        ::FlatModel,
                             data         ::OIdata;
    lb                       = Dict{String,Float64}(),
    ub                       = Dict{String,Float64}(),
    weights                  = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
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

    prior_transform_vectorized = let trafo = prior_transform
        (U::AbstractMatrix{<:Real}) -> reduce(vcat, (u -> trafo(u)').(eachrow(U)))
    end

    # ── Log-likelihood: -0.5 * chi2 ─────────────────────────────────────────
    loglikelihood = let fm=fm, data=data, weights=weights, vonmises=vonmises
        param::AbstractVector{<:Real} -> -0.5 * chi2_flat(fm, param, data;
                                                           weights, vonmises)
    end

    loglikelihood_vectorized = let loglikelihood = loglikelihood
        (X::AbstractMatrix{<:Real}) -> loglikelihood.(eachrow(X))
    end

    # ── Run UltraNest ────────────────────────────────────────────────────────
    ultranest = pyimport("ultranest")

    if !verb
        log_interval = 1_000_000
    end

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

    result = smplr.run(;
        min_num_live_points     = min_num_live_points,
        cluster_num_live_points = cluster_num_live_points,
        log_interval            = log_interval,
        frac_remain             = frac_remain,
    )

    # ── Extract results ──────────────────────────────────────────────────────
    minx = Float64.(result["maximum_likelihood"]["point"])
    minf = chi2_flat(fm, minx, data; weights, vonmises)

    if verb
        @printf("Best-fit χ² = %.4f\n", minf)
        smplr.print_results()
    end

    if cornerplot
        PyDict(pyimport("matplotlib")."rcParams")["font.size"] = [8]
        histogram_color = "black"
        contour_colors = py"{'colors': ['#0072B2','#56B4E9','#009E73','#F0E442'], 'linestyles': ['-', '-', '-', '-']}"
        pyimport("ultranest.plot").cornerplot(result;
            contour_kwargs = contour_colors,
            color          = histogram_color,
        )
    end

    samples_py = result["samples"]
    posterior  = convert(Matrix{Float64}, samples_py)
    logz    = Float64(result["logz"])
    logzerr = Float64(result["logzerr"])
    ndof    = _ndof(data, weights)

    return UltraNestResult(
        minx, list_free_params, minf, minf / ndof, ndof,
        logz, logzerr, posterior, result, fm,
    )
end


# ─────────────────────────────────────────────────────────────────────────────
# LsqFit result
# ─────────────────────────────────────────────────────────────────────────────

struct LsqFitResult
    x_opt       ::Vector{Float64}       # best-fit parameter values
    list_free_params  ::Vector{String}        # parameter names, same order as x_opt
    chi2        ::Float64               # raw chi2 at x_opt
    chi2r       ::Float64               # chi2 / ndof
    ndof        ::Int
    covar       ::Matrix{Float64}       # parameter covariance matrix
    stderror    ::Vector{Float64}       # 1σ uncertainties (sqrt of diag(covar))
    converged   ::Bool                  # did LsqFit converge?
    model       ::FlatModel
    lsqfit_result ::LsqFit.LsqFitResult  # raw LsqFit result (for confidence_interval etc.)
end

function Base.show(io::IO, r::LsqFitResult)
    @printf(io, "LsqFitResult: χ²r = %.4f  (chi2=%.2f, ndof=%d, converged=%s)\n",
            r.chi2r, r.chi2, r.ndof, r.converged)
    for (i, p) in enumerate(r.list_free_params)
        @printf(io, "  %-20s = %10.4f ± %.4f\n", p, r.x_opt[i], r.stderror[i])
    end
end


# ─────────────────────────────────────────────────────────────────────────────
# fit_model_lsqfit — Levenberg-Marquardt via LsqFit.jl
# ─────────────────────────────────────────────────────────────────────────────

"""
    fit_model_lsqfit(model_dict, list_free_params, data; kwargs...) -> LsqFitResult

Fit a parametric model to interferometric data using Levenberg-Marquardt
(LsqFit.jl).  Returns parameter covariance and 1σ uncertainties.

Uses the analytic Jacobian from `residuals_flat_jac` (Wirtinger chain rule
through the complex visibility Jacobian), so no finite-difference overhead.

# Arguments
- `model_dict::Dict{String}` — flat parameter dictionary
- `list_free_params::Vector{String}` — names of the free parameters to optimize
- `data::OIdata` — interferometric data

# Keywords
- `lb`, `ub` — `Dict{String,Float64}` of lower/upper bounds per parameter
- `weights` — observable weights [V2, T3amp, T3phi, visamp, visphi, flux, diffphase]
- `vonmises` — use von Mises statistic for T3phi
- `nB_workspace` — Hankel workspace size
- `maxIter` — maximum iterations (default 200)
- `verb` — print per-evaluation chi2 breakdown
"""
function fit_model_lsqfit(model_dict   ::Dict{String},
                          list_free_params   ::Vector{String},
                          data         ::OIdata;
    lb           = Dict{String,Float64}(),
    ub           = Dict{String,Float64}(),
    weights      = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    vonmises     = false,
    nB_workspace = nothing,
    maxIter      = 200,
    verb         = false,
)
    # ── 1. Compile the model (once) ──────────────────────────────────────────
    nB = something(nB_workspace, size(data.uv, 2))
    fm = parse_model(model_dict, list_free_params; nB_workspace=nB)
    nparams = length(list_free_params)

    # ── 2. Build the model function and Jacobian for LsqFit ─────────────────
    # LsqFit convention: model_fn(xdata, p) -> predicted ydata
    # We use xdata = [1] as a dummy — the data is captured in the closure.
    # The "ydata" is zeros (we want residuals = predicted - 0 = residuals).

    function model_fn(_, p)
        residuals_flat(fm, p, data; weights, vonmises)
    end

    function jacobian_fn(_, p)
        _, J = residuals_flat_jac(fm, p, data; weights, vonmises)
        return J
    end

    # ── 3. Initial parameter vector and bounds ───────────────────────────────
    x0 = Float64[model_dict[p] for p in list_free_params]
    xdata = [1]  # dummy
    ydata = zeros(length(model_fn(xdata, x0)))  # target = 0 (residuals)

    # Build lower/upper bound vectors
    lb_vec = Float64[get(lb, p, -Inf) for p in list_free_params]
    ub_vec = Float64[get(ub, p,  Inf) for p in list_free_params]

    # ── 4. Run LsqFit ───────────────────────────────────────────────────────
    fit = curve_fit(model_fn, jacobian_fn, xdata, ydata, x0;
                    lower=lb_vec, upper=ub_vec,
                    maxIter=maxIter, inplace=false)

    # ── 5. Extract results ───────────────────────────────────────────────────
    minx = fit.param
    r_final = fit.resid
    chi2_val = dot(r_final, r_final)
    ndof = _ndof(data, weights)

    # Covariance matrix from LsqFit (may fail if Jacobian is singular,
    # e.g. when parameters hit bounds)
    cov, σ = try
        vcov(fit), stderror(fit)
    catch
        @warn "Covariance estimation failed (singular Jacobian — parameters may be at bounds)"
        fill(NaN, nparams, nparams), fill(NaN, nparams)
    end

    if verb
        chi2_flat(fm, minx, data; weights, verb=true, vonmises)
    end

    return LsqFitResult(minx, list_free_params, chi2_val, chi2_val / ndof, ndof,
                         cov, σ, fit.converged, fm, fit)
end

"""
    fit_model_lsqfit(model::FlatModel, x0, data; kwargs...) -> LsqFitResult

Levenberg-Marquardt fit on a pre-compiled `FlatModel` starting from `x0`.
Bounds `lb` and `ub` can be Dicts or Vectors.
"""
function fit_model_lsqfit(model        ::FlatModel,
                          x0           ::AbstractVector{<:Real},
                          data         ::OIdata;
    lb           = Dict{String,Float64}(),
    ub           = Dict{String,Float64}(),
    weights      = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    vonmises     = false,
    maxIter      = 200,
    verb         = false,
)
    list_free_params = model.list_free_params
    fm = model
    nparams = length(list_free_params)

    function model_fn(_, p)
        residuals_flat(fm, p, data; weights, vonmises)
    end

    function jacobian_fn(_, p)
        _, J = residuals_flat_jac(fm, p, data; weights, vonmises)
        return J
    end

    xdata = [1]
    ydata = zeros(length(model_fn(xdata, Float64.(x0))))

    lb_vec = _bounds_vec(lb, list_free_params, -Inf)
    ub_vec = _bounds_vec(ub, list_free_params,  Inf)

    fit = curve_fit(model_fn, jacobian_fn, xdata, ydata, Float64.(x0);
                    lower=lb_vec, upper=ub_vec,
                    maxIter=maxIter, inplace=false)

    minx = fit.param
    r_final = fit.resid
    chi2_val = dot(r_final, r_final)
    ndof = _ndof(data, weights)

    cov, σ = try
        vcov(fit), stderror(fit)
    catch
        @warn "Covariance estimation failed (singular Jacobian — parameters may be at bounds)"
        fill(NaN, nparams, nparams), fill(NaN, nparams)
    end

    if verb
        chi2_flat(fm, minx, data; weights, verb=true, vonmises)
    end

    return LsqFitResult(minx, list_free_params, chi2_val, chi2_val / ndof, ndof,
                         cov, σ, fit.converged, fm, fit)
end


# ─────────────────────────────────────────────────────────────────────────────
# model_to_obs — compute observables from a FlatModel
# ─────────────────────────────────────────────────────────────────────────────

"""
    model_to_obs(model, x, data) -> (v2, t3amp, t3phi, visamp, visphi)

Evaluate a `FlatModel` at parameter vector `x` and return the model
observables matching the data structure.

Returns a NamedTuple with fields `v2`, `t3amp`, `t3phi` (degrees),
`visamp`, `visphi` (degrees).  Empty vectors are returned for
observable types not present in `data`.
"""
function model_to_obs(model::FlatModel, x::AbstractVector, data::OIdata)
    cvis = eval_model(model, x, data.uv)

    # V²
    v2_model = data.nv2 > 0 ? abs2.(cvis[data.indx_v2]) : Float64[]

    # Triple products
    if data.nt3amp > 0 || data.nt3phi > 0
        t3 = cvis[data.indx_t3_1] .* cvis[data.indx_t3_2] .* cvis[data.indx_t3_3]
        t3amp_model = abs.(t3)
        t3phi_model = angle.(t3) .* (180.0 / π)
    else
        t3amp_model = Float64[]
        t3phi_model = Float64[]
    end

    # Visibility amplitude & phase
    if data.nvisamp > 0 || data.nvisphi > 0
        Vvis = cvis[data.indx_vis]
        visamp_model = abs.(Vvis)
        visphi_model = angle.(Vvis) .* (180.0 / π)
    else
        visamp_model = Float64[]
        visphi_model = Float64[]
    end

    return (v2=v2_model, t3amp=t3amp_model, t3phi=t3phi_model,
            visamp=visamp_model, visphi=visphi_model)
end

"""
    model_to_residuals(model, x, data)

Compute normalised residuals `(model - data) / error` for each observable type.
Returns a NamedTuple with fields `v2`, `t3amp`, `t3phi`, `visamp`, `visphi`.
Phase residuals (t3phi, visphi) are wrapped to [-180, 180] before dividing by error.
"""
function model_to_residuals(model::FlatModel, x::AbstractVector, data::OIdata)
    obs = model_to_obs(model, x, data)

    v2_res     = data.nv2 > 0 ? (obs.v2 .- data.v2) ./ data.v2_err : Float64[]
    t3amp_res  = data.nt3amp > 0 ? (obs.t3amp .- data.t3amp) ./ data.t3amp_err : Float64[]
    t3phi_res  = data.nt3phi > 0 ? _mod360(obs.t3phi .- data.t3phi) ./ data.t3phi_err : Float64[]
    visamp_res = data.nvisamp > 0 ? (obs.visamp .- data.visamp) ./ data.visamp_err : Float64[]
    visphi_res = data.nvisphi > 0 ? _mod360(obs.visphi .- data.visphi) ./ data.visphi_err : Float64[]

    return (v2=v2_res, t3amp=t3amp_res, t3phi=t3phi_res,
            visamp=visamp_res, visphi=visphi_res)
end

# ─────────────────────────────────────────────────────────────────────────────
# model_to_image — synthesise an image from a FlatModel via inverse FFT
# ─────────────────────────────────────────────────────────────────────────────

"""
    model_to_sed(model, x, wl_grid) -> (total, components)

Evaluate the spectral energy distribution of a `FlatModel` by computing
`V(u=0, v=0)` at each wavelength — the zero-baseline visibility equals
the total flux.

Per-component fluxes are obtained by resolving the `f` / `spectrum`
parameters at each wavelength.

# Arguments
- `model::FlatModel` — compiled model from `dict_to_model`
- `x::AbstractVector` — current free-parameter values
- `wl_grid` — wavelength array in microns

# Returns
- `total::Vector{Float64}` — total flux at each wavelength (= V(0,0))
- `components::Dict{String,Vector{Float64}}` — per-component flux
"""
function model_to_sed(model::FlatModel,
                      x::AbstractVector,
                      wl_grid::AbstractVector)
    nwl = length(wl_grid)

    # Total flux = V(0,0) at each wavelength
    uv_zero = zeros(2, 1)   # single (u,v) = (0,0) point
    total = Vector{Float64}(undef, nwl)
    for (i, wl) in enumerate(wl_grid)
        V = eval_model(model, x, uv_zero; wl)
        total[i] = real(V[1])
    end

    # Per-component fluxes from the resolver
    comp_fluxes = Dict{String, Vector{Float64}}()
    for comp in model.components
        comp_fluxes[comp.comp_name] = zeros(nwl)
    end
    for (i, wl) in enumerate(wl_grid)
        pv = model.resolver.fn(x, wl, NaN)
        for comp in model.components
            f = comp.f_idx == 0 ? 1.0 : pv[comp.f_idx]
            comp_fluxes[comp.comp_name][i] = f
        end
    end

    return total, comp_fluxes
end


const _MAS2RAD_IMG = 2.0626480624709636e8   # 1 rad in mas (same as _MAS2RAD_PM)

"""
    model_to_image(model, x; nx=256, pixsize=0.1, oversample=1, normalize=true, wl=nothing)
        -> Matrix{Float64}

Synthesise an image from a `FlatModel` at parameter vector `x`.

The image is computed by evaluating the model visibility on the 2-D FFT
frequency grid and applying an inverse real FFT (`irfft`), which exploits
the Hermitian symmetry V(-u,-v) = conj(V(u,v)) to halve the number of
visibility evaluations.

# Arguments
- `model::FlatModel` — compiled model from `dict_to_model`
- `x::AbstractVector` — current free-parameter values

# Keywords
- `nx` — image size in pixels (default 256)
- `pixsize` — pixel size in mas (default 0.1)
- `oversample` — oversampling factor (default 1)
- `normalize` — normalize image to unit sum (default true)
- `wl` — wavelength in microns (scalar); enables `\$WL` in spectrum expressions

# Returns
`img` — `nx × nx` Matrix{Float64}
"""
function model_to_image(model::FlatModel,
                        x::AbstractVector;
                        nx::Int       = 256,
                        pixsize::Real = 0.1,
                        oversample::Int = 1,
                        normalize::Bool = true,
                        wl = nothing)

    ns = nx * oversample
    ps = pixsize / oversample   # pixel size after oversampling (mas)

    # ── Spatial frequency grids (cycles/rad) ────────────────────────────────
    # Natural FFT order (NOT shifted) — matches what irfft expects.
    # u-axis (dim 1): rfft gives frequencies 0, 1/(ns·ps), ..., 1/(2·ps)
    # v-axis (dim 2): fft  gives frequencies 0, 1/(ns·ps), ..., -1/(ns·ps)
    nrfft = ns ÷ 2 + 1
    freq_u = rfftfreq(ns) ./ ps .* _MAS2RAD_IMG   # cycles/rad, length nrfft
    freq_v = fftfreq(ns)  ./ ps .* _MAS2RAD_IMG   # cycles/rad, length ns

    # ── Build UV grid for the half-plane (rfft convention) ──────────────────
    nuv_half = nrfft * ns
    uv = Matrix{Float64}(undef, 2, nuv_half)
    k = 0
    for jv in 1:ns         # v index (full, natural order)
        for iu in 1:nrfft  # u index (rfft half)
            k += 1
            uv[1, k] = freq_u[iu]
            uv[2, k] = freq_v[jv]
        end
    end

    # ── Evaluate model visibility on the half-plane ─────────────────────────
    V_half = eval_model(model, x, uv; wl)

    # Reshape to (nrfft, ns) — rfft output layout
    V_mat = reshape(V_half, nrfft, ns)

    # ── Inverse real FFT → real image ───────────────────────────────────────
    # irfft(F, ns) on a 2D array does irfft along dim 1 and ifft along dim 2,
    # recovering the full ns×ns real image.  fftshift centers the result.
    img = fftshift(irfft(V_mat, ns))

    # Clamp negatives (numerical noise from Gibbs ringing)
    @. img = max(img, 0.0)

    # ── Downsample oversampled image by block-averaging ───────────────────
    if oversample > 1
        img_ds = zeros(nx, nx)
        for j in 1:nx
            for i in 1:nx
                s = 0.0
                for dj in 0:(oversample-1)
                    for di in 0:(oversample-1)
                        s += img[(i-1)*oversample + di + 1,
                                 (j-1)*oversample + dj + 1]
                    end
                end
                img_ds[i, j] = s / oversample^2
            end
        end
        img = img_ds
    end

    if normalize && sum(img) > 0
        img ./= sum(img)
    end

    return img
end

"""
    resample_data(data::OIdata)

Create a bootstrapped copy of the data by adding Gaussian noise
(drawn from the error bars) to the observables. Useful for bootstrap
uncertainty estimation.
"""
function resample_data(data::OIdata)
    d = deepcopy(data)
    d.v2      .+= d.v2_err      .* randn(length(d.v2))
    d.t3amp   .+= d.t3amp_err   .* randn(length(d.t3amp))
    d.t3phi   .+= d.t3phi_err   .* randn(length(d.t3phi))
    d.visamp  .+= d.visamp_err  .* randn(length(d.visamp))
    d.visphi  .+= d.visphi_err  .* randn(length(d.visphi))
    if length(d.flux) > 0
        d.flux .+= d.flux_err   .* randn(length(d.flux))
    end
    return d
end
