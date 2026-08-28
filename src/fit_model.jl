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

using NLopt, LsqFit, ForwardDiff, Printf, FFTW


# ─────────────────────────────────────────────────────────────────────────────
# Vectorized likelihood evaluation for UltraNest
# ─────────────────────────────────────────────────────────────────────────────

"""
Batch size below which evaluating the rows serially beats threading them.

`Threads.@threads` costs a few μs to spawn and join, which a one-row batch cannot amortise.
Measured on a 4-parameter analytic model against `2004-data1.oifits` (≈85 μs per χ²
evaluation, 16 threads), speedup of the threaded batch over the serial one:

    batch      1     2     5    15    20    50   128   400
    speedup 0.77  1.33  2.99  6.15  7.18  8.40  8.88  9.46

so the break-even sits between 1 and 2. The threshold is set above that rather than at it,
because the curve shifts right for cheaper likelihoods and the points lost to it are few:
UltraNest's batches are small but rarely 2 or 3 (median 15 without a step sampler).
"""
const _MIN_THREADED_BATCH = 4

"""
    _batch_loglike(f, X) -> Vector

Apply `f` to each row of the batch `X`, threading over rows when that is worth it.

`vectorized = true` only means UltraNest hands over a whole batch per call; evaluating it is
still up to us, and a serial broadcast leaves every thread but one idle. Note the batch is
one row per call whenever a step sampler is active (`StepSampler.__next__` evaluates a single
reshaped point), so the guard below is what keeps `use_stepsampler = true` from paying for
threading it cannot use.

`X` arrives as a `PyArray`: a zero-copy view onto numpy's memory, not a Julia array. Only the
thread holding the GIL may touch Python, so the batch is copied here, on the calling thread,
before any worker looks at it. At these batch sizes that copy is a few hundred bytes against
`n` χ² evaluations. Nothing inside `f` reaches Python.

The rows share one `FlatModel`, which is safe because evaluation does not write to it:
`eval_model` allocates its accumulator per call and never touches `HankelSpec.workspace` —
those buffers are allocated but currently unused, since `eval_model_grad` goes through
ForwardDiff rather than the cached-K chain rule they were meant for. Wiring that up would
make the workspace shared mutable state and this loop a race; `test_python_boundary.jl`
asserts the buffers stay untouched so that change fails loudly here rather than silently
corrupting a posterior. `bootstrap.jl` gives each worker its own model for the same reason.
"""
function _batch_loglike(f, X::AbstractMatrix)
    n = size(X, 1)
    (Threads.nthreads() == 1 || n < _MIN_THREADED_BATCH) && return f.(eachrow(X))

    Xj = Matrix{eltype(X)}(X)          # off the PyArray before the threads start
    v1 = f(view(Xj, 1, :))             # fixes the element type without assuming Float64
    out = Vector{typeof(v1)}(undef, n)
    out[1] = v1
    Threads.@threads for k in 2:n
        @inbounds out[k] = f(view(Xj, k, :))
    end
    return out
end


# ─────────────────────────────────────────────────────────────────────────────
# display_model — inspect model setup before fitting
# ─────────────────────────────────────────────────────────────────────────────

"""
    model_warnings(model_dict, list_free_params; lb, ub) -> Vector{String}

Everything wrong with a model that can be seen without fitting it, one sentence each.

Four checks, and every one of them describes a setup that runs happily and answers the wrong
question: a starting value outside its own bounds, bounds that exclude everything, and flux
fractions that do not sum to one. The last is the quiet one -- a second component added with
`f = 1` doubles the model's flux, and the only symptom is a χ² in the millions.

Shared with [`display_model`](@ref), which prints these to a terminal, so a GUI and a script
cannot disagree about what is wrong with the same model.

Flux fractions given as expressions are not checked: their values depend on the resolver, so
summing the literals would be summing the wrong things.
"""
function model_warnings(model_dict::Dict{String}, list_free_params::AbstractVector{<:AbstractString};
                        lb = Dict{String,Float64}(), ub = Dict{String,Float64}())
    w = String[]
    for p in list_free_params
        v = get(model_dict, p, nothing)
        v isa Number || continue
        lo, hi = get(lb, p, nothing), get(ub, p, nothing)
        lo !== nothing && v < lo && push!(w, "$p starts at $v, below its lower bound $lo")
        hi !== nothing && v > hi && push!(w, "$p starts at $v, above its upper bound $hi")
        lo !== nothing && hi !== nothing && lo >= hi &&
            push!(w, "$p has an empty box: lower $lo is not below upper $hi")
    end

    fkeys = [k for k in keys(model_dict) if endswith(k, ",f")]
    numeric = [model_dict[k] for k in fkeys if model_dict[k] isa Number]
    if !isempty(numeric) && length(numeric) == length(fkeys)
        fsum = sum(numeric)
        abs(fsum - 1.0) > 0.01 && push!(w,
            "flux fractions sum to $(round(fsum; digits = 4)), not 1 — every component is " *
            "carrying its own full flux, which inflates χ² rather than sharing the source " *
            "between them")
    end
    return w
end

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

"""
    FitResult

What every point-estimate fitter returns: `fit_model`, `fit_model_lsqfit` and the grid search.

| field | meaning |
|---|---|
| `x_opt` | best-fit values, in `list_free_params` order |
| `list_free_params` | the names, so `x_opt[i]` is always identifiable |
| `chi2`, `chi2r`, `ndof` | raw weighted χ², χ²/ndof, and the number of DATA POINTS |
| `n_evals` | χ² evaluations spent |
| `ret` | the optimiser's return code |
| `model` | the compiled `FlatModel` the fit ran against |

`ndof` counts data points and is **not** reduced by the number of free parameters, so it is not
a degrees-of-freedom in the model-comparison sense; compute AIC or BIC from `chi2` and the
parameter count rather than reading them off `chi2r`.

See also `LsqFitResult`, [`UltraNestResult`](@ref) and [`Chi2Map`](@ref).
"""
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
- `constraints` — relations between parameters that bounds cannot express, as
  [`ModelConstraint`](@ref)s or `(lhs, op, rhs[, tol])` tuples. Handed to NLopt as real
  nonlinear constraints, so they hold at the optimum; a `method` that cannot take them is
  wrapped in `:AUGLAG` rather than replaced. The reported `chi2` is the data χ², unaffected
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
    constraints  = [],
    method       = :LD_LBFGS,
    maxeval      = 2000,
    ftol_rel     = 1e-8,
    xtol_rel     = 1e-6,
    vonmises     = false,
    nB_workspace = nothing,
    verb         = false,
)
    # ── 1. Augment model_dict with prior and constraint expressions ──────────
    # Both reach the optimiser as extra parameters of the compiled model: the resolver is the
    # only thing that can evaluate an expression over parameter names, so anything needing such
    # a value has to be compiled alongside the model rather than beside it.
    augmented = copy(model_dict)
    prior_keys = String[]
    for (i, (expr, _, _)) in enumerate(priors)
        key = "__prior_$(i)__"
        augmented[key] = _dollarify(expr)
        push!(prior_keys, key)
    end
    cons = parse_constraints(constraints)
    constraint_keys = augment_constraints!(augmented, cons)

    # ── 2. Compile the model (once) ──────────────────────────────────────────
    nB = something(nB_workspace, size(data.uv, 2))
    fm = parse_model(augmented, list_free_params; nB_workspace=nB)

    # ── 3. Resolve prior and constraint indices ──────────────────────────────
    prior_idx      = [findfirst(==(k), fm.all_names) for k in prior_keys]
    constraint_idx = [findfirst(==(k), fm.all_names) for k in constraint_keys]

    # ── 4. Detect gradient need ──────────────────────────────────────────────
    use_grad = startswith(string(method), "LD_")

    # ── 5. Build objective ───────────────────────────────────────────────────
    n_evals = Ref(0)
    function objective(x, g)
        n_evals[] += 1
        if use_grad
            chi2, grad = chi2_flat_fg(fm, x, data; weights, vonmises, verb)
            chi2 += prior_penalty_and_grad!(grad, x, priors, fm, prior_idx)
            # NLopt hands back a zero-length array when it wants the value without the
            # gradient, which AUGLAG does at some of its evaluations.
            length(g) > 0 && (g .= grad)
        else
            chi2 = chi2_flat(fm, x, data; weights, vonmises, verb)
            chi2 += prior_penalty(x, priors, fm, prior_idx)
        end
        return chi2
    end

    # ── 6. Configure NLopt ───────────────────────────────────────────────────
    # Constraints go to NLopt itself rather than into the objective as a penalty, so that
    # `diamout > diamin` holds at the optimum instead of merely being encouraged there. When
    # the chosen algorithm cannot take them, `constrained_opt` wraps it in AUGLAG rather than
    # substituting an algorithm the caller did not ask for.
    x0 = Float64[model_dict[p] for p in list_free_params]
    opt, wrapped = constrained_opt(method, length(list_free_params), cons;
                                   ftol_rel, xtol_rel, maxeval)
    min_objective!(opt, objective)
    lower_bounds!(opt, [get(lb, p, -Inf) for p in list_free_params])
    upper_bounds!(opt, [get(ub, p,  Inf) for p in list_free_params])
    add_nlopt_constraints!(opt, cons, fm, constraint_idx, x0)
    wrapped && @info "$(length(cons)) constraint(s) with $method, which takes none directly: " *
                     "wrapped in AUGLAG with $method as the local optimiser." maxlog = 1

    # ── 7. Run ───────────────────────────────────────────────────────────────
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
    # A constraint the optimiser could not satisfy is a fact about the fit. Reporting it beats
    # a silently infeasible answer, which looks exactly like a feasible one.
    for (c, v) in zip(cons, constraint_violations(minx, cons, fm, constraint_idx))
        v > 1 && @warn "Fit ended outside a constraint" constraint = c violation_in_tol_units = v
    end

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
    constraints  = [],
    method       = :LD_LBFGS,
    maxeval      = 2000,
    ftol_rel     = 1e-8,
    xtol_rel     = 1e-6,
    vonmises     = false,
    verb         = false,
)
    list_free_params = model.list_free_params
    fm = model

    # ── Priors and constraints ───────────────────────────────────────────────
    # Both are evaluated as extra parameters of the compiled model, so both need the model dict
    # to augment and recompile. A FlatModel is already compiled, and silently dropping them
    # would return a fit that ignores what the caller asked for.
    !isempty(priors) && error("Priors are only supported with the Dict form of fit_model; pass model_dict instead of FlatModel")
    !isempty(constraints) && error("Constraints are only supported with the Dict form of fit_model; pass model_dict instead of FlatModel")

    # ── Detect gradient need ─────────────────────────────────────────────────
    use_grad = startswith(string(method), "LD_")

    # ── Build objective ──────────────────────────────────────────────────────
    n_evals = Ref(0)
    function objective(x, g)
        n_evals[] += 1
        if use_grad
            chi2, grad = chi2_flat_fg(fm, x, data; weights, vonmises, verb)
            # NLopt hands back a zero-length array when it wants the value without the gradient.
            length(g) > 0 && (g .= grad)
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
# Nested sampling result
# ─────────────────────────────────────────────────────────────────────────────

"""
    NestedResult

Result of [`fit_model_nested`](@ref), whichever sampler produced it.

| field | type | meaning |
|---|---|---|
| `x_opt` | `Vector{Float64}` | maximum-likelihood parameter values |
| `list_free_params` | `Vector{String}` | parameter names, same order as `x_opt` |
| `chi2`, `chi2r`, `ndof` | `Float64`, `Float64`, `Int` | chi² at `x_opt`, reduced chi², degrees of freedom |
| `logz`, `logzerr` | `Float64` | Bayesian log-evidence and its uncertainty |
| `posterior` | `Matrix{Float64}` | `(n_samples, n_params)` equally weighted posterior samples |
| `result` | `Any` | the sampler's own result object, passed through unconverted |
| `model` | `FlatModel` | the compiled model that was fitted |
| `backend` | `Symbol` | which sampler produced this — `:nautilus` or `:ultranest` |

`backend` is not decoration. Two independent nested samplers agreeing on `logz` is evidence;
one of them quoted without saying which is not, because the estimator, the bounding strategy
and the stopping rule all differ between them.

!!! note "`result` is backend-specific"

    Under `:ultranest` it is a `PythonCall.Py` wrapping UltraNest's result dict — index it with
    `result["key"]` and convert explicitly with `pyconvert(T, …)`. Under `:nautilus` it
    is the sampler state Nautilus returned. Every other field is a concrete
    Julia type and means the same thing either way; prefer those where they suffice.
"""
struct NestedResult
    x_opt       ::Vector{Float64}   # maximum-likelihood parameter values
    list_free_params  ::Vector{String}    # parameter names, same order as x_opt
    chi2        ::Float64           # chi2 at x_opt
    chi2r       ::Float64           # chi2 / ndof
    ndof        ::Int
    logz        ::Float64           # log-evidence
    logzerr     ::Float64           # log-evidence uncertainty
    posterior   ::Matrix{Float64}   # (n_samples, n_params) posterior samples
    result      ::Any               # the sampler's own result object
    model       ::FlatModel
    backend     ::Symbol            # :nautilus or :ultranest
end

"""
    UltraNestResult

Alias for [`NestedResult`](@ref), so code written against the UltraNest-only API keeps working.
"""
const UltraNestResult = NestedResult

function Base.show(io::IO, r::NestedResult)
    @printf(io, "NestedResult (%s): χ²r = %.4f  (chi2=%.2f, ndof=%d)\n",
            r.backend, r.chi2r, r.chi2, r.ndof)
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
# Nested sampling backends
# ─────────────────────────────────────────────────────────────────────────────
#
# Two samplers implement the same fit, and neither is a dependency of the package:
#
#   :nautilus         Nautilus.jl       -- pure Julia, importance nested sampling; reuses every
#                                          likelihood call, so it buys far more evidence
#                                          precision per call (and, in higher dimensions,
#                                          needs far fewer). See src/fit_model_nautilus.jl.
#   :nautilus         Nautilus.jl       -- pure Julia, so it survives into a compiled build
#   :ultranest        UltraNest         -- Python, reached through PythonCall
#
# `using` either one activates its extension, which registers the backend and adds a method to
# `_fit_nested`. Keeping both is deliberate: an evidence estimate checked against a second,
# independent implementation is the only cheap way to find out that one of them is wrong.

"""
The nested sampler [`fit_model_nested`](@ref) uses when no `backend` keyword is given.

`:none` until an extension registers one, after which the first loaded wins. Change it with
[`set_nested_backend!`](@ref) rather than assigning to this directly, which checks that the
backend is actually available.
"""
const NESTED_BACKEND = Ref(:none)

"Backends whose extension has loaded. Written by the extensions' `__init__`."
const NESTED_BACKENDS_LOADED = Set{Symbol}()

function _register_nested_backend!(b::Symbol)
    push!(NESTED_BACKENDS_LOADED, b)
    NESTED_BACKEND[] === :none && (NESTED_BACKEND[] = b)
    return b
end

"""
    nested_backend() -> Symbol

Which nested sampler [`fit_model_nested`](@ref) will use: `:nautilus`,
`:ultranest`, or `:none` when no extension is loaded.
"""
nested_backend() = NESTED_BACKEND[]

"""
    set_nested_backend!(b::Symbol) -> Symbol

Choose the nested sampler: `:nautilus` or `:ultranest`.

The backend's package must already be loaded, and this says so immediately rather than letting
a long run end in a `MethodError`.
"""
function set_nested_backend!(b::Symbol)
    if !(b in NESTED_BACKENDS_LOADED)
        loaded = isempty(NESTED_BACKENDS_LOADED) ? "none" :
                 join(sort(collect(NESTED_BACKENDS_LOADED)), ", ")
        error("Nested sampling backend :$b is not loaded (loaded: $loaded). " *
              "`using Nautilus` enables :nautilus, " *
              "`using PythonCall` :ultranest.")
    end
    return NESTED_BACKEND[] = b
end

"Backend implementations. `OITOOLSNautilusExt` and `OITOOLSUltraNestExt` add the methods."
function _fit_nested end

"""
    fit_model_nested(model_dict, list_free_params, data; kwargs...) -> NestedResult
    fit_model_nested(model::FlatModel, data; kwargs...)             -> NestedResult

Fit a parametric model by Bayesian nested sampling, returning the posterior and the evidence
as well as a best-fit point.

Unlike the gradient fitters this needs no starting guess — it samples the whole prior volume —
but it does need **finite bounds on every free parameter**, since those bounds are the prior.

```julia
using OITOOLS, Nautilus
r = fit_model_nested(model_dict, ["star,ud"], data;
                     lb = Dict("star,ud" => 0.5), ub = Dict("star,ud" => 10.0))
r.logz, r.logzerr        # the evidence, for comparing models
r.posterior              # equally weighted samples, (n_samples, n_params)
```

`backend` selects the sampler and defaults to [`nested_backend()`](@ref nested_backend); see
[`set_nested_backend!`](@ref) to change it for a session.

# Keywords, both backends

| keyword | default | meaning |
|---|---|---|
| `lb`, `ub` | — | `Dict{String,Float64}` bounds, **required** for every free parameter: they are the prior |
| `weights` | `[1,1,1,0,0,0,0]` | observable weights `[V2, T3amp, T3phi, visamp, visphi, flux, diffphase]` |
| `vonmises` | `false` | use the von Mises statistic for T3phi |
| `constraints` | `[]` | [`ModelConstraint`](@ref)s, or `(lhs, op, rhs[, tol])` tuples. They enter the log-likelihood as `-penalty/2`, so the posterior and `logz` are both the constrained ones |
| `nB_workspace` | from the data | Hankel workspace size |
| `verb` | `true` | print progress |
| `cornerplot` | `true` | draw the corner plot when a plotting backend is loaded |

# Keywords, `:nautilus`

| keyword | default | meaning |
|---|---|---|
| `n_live` | `500` | live points in the exploration phase |
| `n_eff` | `2000` | effective sample size to reach before stopping |
| `f_live` | `0.01` | fraction of evidence left in the live set at which exploration ends |
| `n_networks` | `4` | networks in the bounding ensemble; more is more robust and slower |
| `threaded` | `nthreads() > 1` | evaluate the likelihood across threads |
| `seed` | `nothing` | fix for a reproducible run |

# Keywords, `:ultranest`

| keyword | default | meaning |
|---|---|---|
| `min_num_live_points` | `400` | minimum live points |
| `cluster_num_live_points` | `100` | live points per cluster |
| `num_bootstraps` | `30` | bootstraps for the evidence error |
| `use_stepsampler` | `false` | use `RegionSliceSampler` |
| `nsteps` | `400` | steps per slice, with `use_stepsampler` |
| `frac_remain` | `0.001` | termination fraction |
| `log_interval` | `100` | logging interval |

Evidence from two different samplers is worth more than evidence from one: `logz` should agree
within `logzerr`, and `NestedResult.backend` records which produced any given number.

!!! note "Nautilus is importance nested sampling"

    `:nautilus` draws from neural-network-boosted importance shells rather than by rejection
    inside a bounding ellipsoid, which buys it far more effective samples per likelihood call.
    Measured on the α Cen A and α Cen B power-law fits, its `logzerr` is **0.019** against
    UltraNest's 0.18–0.31 on the same data, while `x_opt` agrees to four decimals — so the two
    can be cross-checked against each other and the Julia one gives the tighter evidence.

    `threaded` defaults to `Threads.nthreads() > 1`; `seed` makes a run reproducible.

"""
function fit_model_nested(args...; backend::Symbol = nested_backend(), kwargs...)
    if backend === :none
        error("""No nested sampler is loaded. Add one:
                     using Nautilus         # pure Julia, importance nested sampling
                     using Nautilus         # pure Julia, no Python
                     using PythonCall       # UltraNest, needs the Python package""")
    end
    backend in NESTED_BACKENDS_LOADED ||
        error("Nested sampling backend :$backend is not loaded; loaded: " *
              (isempty(NESTED_BACKENDS_LOADED) ? "none" :
               join(sort(collect(NESTED_BACKENDS_LOADED)), ", ")))
    return _fit_nested(Val(backend), args...; kwargs...)
end

"""
    fit_model_ultranest(args...; kwargs...) -> NestedResult

[`fit_model_nested`](@ref) pinned to the UltraNest backend. Needs `using PythonCall` and the
Python `ultranest` package (declared in `CondaPkg.toml`).
"""
fit_model_ultranest(args...; kwargs...) =
    fit_model_nested(args...; backend = :ultranest, kwargs...)


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
- `constraints` — relations between parameters that bounds cannot express, as
  [`ModelConstraint`](@ref)s or `(lhs, op, rhs[, tol])` tuples. They enter as extra rows of the
  residual vector, which is exactly where Levenberg-Marquardt can act on them — and so also
  enter the Jacobian, meaning `covar` and `stderror` account for them. The reported `chi2`
  excludes them, so it stays comparable with what every other fitter reports
- `verb` — print per-evaluation chi2 breakdown
"""
function fit_model_lsqfit(model_dict   ::Dict{String},
                          list_free_params   ::Vector{String},
                          data         ::OIdata;
    lb           = Dict{String,Float64}(),
    ub           = Dict{String,Float64}(),
    weights      = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    constraints  = [],
    vonmises     = false,
    nB_workspace = nothing,
    maxIter      = 200,
    verb         = false,
)
    # ── 1. Compile the model (once), with the constraint expressions ─────────
    nB = something(nB_workspace, size(data.uv, 2))
    cons      = parse_constraints(constraints)
    augmented = copy(model_dict)
    constraint_keys = augment_constraints!(augmented, cons)
    fm = parse_model(augmented, list_free_params; nB_workspace=nB)
    constraint_idx  = [findfirst(==(k), fm.all_names) for k in constraint_keys]
    nparams = length(list_free_params)

    # ── 2. Build the model function and Jacobian for LsqFit ─────────────────
    # LsqFit convention: model_fn(xdata, p) -> predicted ydata
    # We use xdata = [1] as a dummy — the data is captured in the closure.
    # The "ydata" is zeros (we want residuals = predicted - 0 = residuals).

    # Constraints append rows. Their count is fixed for a given model — one per constraint, or
    # one per wavelength where the resolver broadcasts — which is what `curve_fit` requires of a
    # residual vector whose length it fixes from the first call.
    function model_fn(_, p)
        r = residuals_flat(fm, p, data; weights, vonmises)
        isempty(cons) ? r : vcat(r, constraint_residuals(p, cons, fm, constraint_idx))
    end

    function jacobian_fn(_, p)
        _, J = residuals_flat_jac(fm, p, data; weights, vonmises)
        isempty(cons) ? J : vcat(J, constraint_jacobian(p, cons, fm, constraint_idx))
    end

    # ── 3. Initial parameter vector and bounds ───────────────────────────────
    x0 = Float64[model_dict[p] for p in list_free_params]
    xdata = [1]  # dummy
    ndata = length(residuals_flat(fm, x0, data; weights, vonmises))
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
    # The data rows only: a chi2 carrying the constraint penalty would not be comparable with
    # what any other fitter reports, nor with the chi2r read off a published table.
    r_final = view(fit.resid, 1:ndata)
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
    constraints  = [],
    vonmises     = false,
    maxIter      = 200,
    verb         = false,
)
    list_free_params = model.list_free_params
    fm = model
    nparams = length(list_free_params)

    # A constraint is evaluated as an extra parameter of the compiled model, so it needs the
    # model dict to augment and recompile. Accepting it here and dropping it would return a fit
    # that ignores what the caller asked for.
    !isempty(constraints) && error("Constraints are only supported with the Dict form of fit_model_lsqfit; pass model_dict instead of FlatModel")

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


"""
    model_to_flux(model, x; wl=nothing)

Total flux of the model at zero baseline: `real(V(0,0))`.

For polychromatic models, pass `wl` as a scalar or 1-element vector.
See also `model_to_sed` for the full spectral energy distribution.
"""
function model_to_flux(model::FlatModel, x::AbstractVector; wl=nothing)
    wl_arg = isnothing(wl) ? nothing : (wl isa Number ? [wl] : wl)
    V = eval_model(model, x, zeros(2, 1); wl=wl_arg)
    return real(V[1])
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
                        wl = nothing,
                        mjd = nothing)

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
    # `mjd` as well as `wl`: a model may depend on either, and `eval_model` has always taken
    # both — only this entry point did not pass it, so a time-dependent model rendered at one
    # arbitrary epoch with nothing saying which.
    V_half = eval_model(model, x, uv; wl, mjd)

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

# `perturb_data` (formerly `resample_data`) and the block bootstrap now live in
# bootstrap.jl.
