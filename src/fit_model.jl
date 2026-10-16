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

# One fit evaluates one model against one dataset many times, so it keeps the model Jacobian's
# buffers between evaluations. A polychromatic fit evaluates each channel's uv set in turn and
# takes no cache.
_jac_cache_kw(::OIdata) = (; jac_cache = ModelJacCache())
_jac_cache_kw(_) = (;)

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
those buffers are allocated but unused: `eval_model_grad` builds a workspace per call, and
the Hankel transform matrices it shares across calls sit in a cache behind a lock. Writing to
`HankelSpec.workspace` would make it shared mutable state and this loop a race; `test_python_boundary.jl`
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
    implicit_vars(key) -> Tuple

The implicit variables an expression stored under `key` can use. Each expression kind is
compiled with its own arguments, so the set is fixed by the key's parameter name:

| key | variables | why |
|---|---|---|
| `c,profile` | `\$R`, `\$MU`, `\$D`, `\$RMIN`, `\$RMAX`, `\$DMIN`, `\$DMAX` | evaluated on the radial grid, before any uv point |
| `c,visfunc` | `\$BL`, `\$WL`, `\$MJD` | evaluated per uv point, on its baseline |
| anything else | `\$WL`, `\$MJD` | resolved once per uv point |

A variable used outside its set raises `UndefVarError` when the model is evaluated, naming the
variable but not the expression.
"""
function implicit_vars(key::AbstractString)
    i = findfirst(',', key)
    par = i === nothing ? key : key[i+1:end]
    par == "profile" && return ("R", "MU", "D", "RMIN", "RMAX", "DMIN", "DMAX")
    par == "visfunc" && return ("BL", "WL", "MJD")
    return ("WL", "MJD")
end

"""
    missing_refs(model_dict) -> Vector{Tuple{String,String}}

The `(key, reference)` pairs where an expression names something it cannot see, sorted by key:
a parameter the model does not define, or an implicit variable outside its
`implicit_vars` set.

The lookup rule is the parser's and it differs by key, so this mirrors it rather than guessing:
a `profile` or `visfunc` expression qualifies a bare `\$a` with its own component and falls back
to a global `a`, while every other expression resolves a reference literally — a derived
`"star,f" = "\$ud * 0.1"` does NOT find `star,ud`.
"""
function missing_refs(model_dict::AbstractDict)
    out = Tuple{String,String}[]
    for key in sort!(String.(collect(keys(model_dict))))
        v = model_dict[key]
        v isa AbstractString || continue
        i = findfirst(',', key)
        comp = i === nothing ? "" : key[1:i-1]
        par  = i === nothing ? key : key[i+1:end]
        qualifies = par == "profile" || par == "visfunc"
        inscope = implicit_vars(key)
        for ref in extract_refs(String(v))
            if ref in IMPLICIT_VARS || ref == "BL"
                ref in inscope || push!(out, (key, ref))
                continue
            end
            found = if qualifies && !occursin(',', ref)
                haskey(model_dict, comp * "," * ref) || haskey(model_dict, ref)
            else
                haskey(model_dict, ref)
            end
            found || push!(out, (key, ref))
        end
    end
    return out
end

# One sentence for a `missing_refs` pair: an implicit variable used out of its scope reads
# differently from a parameter that does not exist, and the fix differs too.
function _missing_ref_message(key, ref)
    ref in IMPLICIT_VARS || ref == "BL" || return "$key refers to \$$ref, which the model does not define"
    ok = implicit_vars(key)
    return "$key uses \$$ref, which this kind of expression cannot see (it has " *
           join(("\$" * v for v in ok), ", ") * ")"
end

"""
    model_warnings(model_dict, list_free_params; lb, ub) -> Vector{String}

Problems in a model that can be found without fitting it, one sentence each:

- a free parameter's starting value outside its bounds, or bounds with `lb ≥ ub`;
- numeric fluxes summing to zero or less (the visibilities are divided by the total flux);
- every flux free, so their common scale is undetermined (fix one, or write one as
  1 − the others);
- an expression referring to a parameter that is not in the model, or to an implicit variable
  its kind of expression does not have (`\$WL` in a profile); evaluating the model then raises
  an `UndefVarError` that does not name the key;
- a key the parser ignores: a non-shape key no expression refers to (unused, or misspelt such
  as `star,diamm`), or a shape key of another kind than the one that set the component's kind
  (`ud` beside `diamin`, which make a ring);
- `\$WL` (microns) compared with or divided by a constant below 1e-3, which looks like metres.

Fluxes given as expressions are not checked. [`display_model`](@ref) prints the bound, flux and
reference warnings.
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

    append!(w, _flux_warnings(model_dict, list_free_params))

    for (key, ref) in missing_refs(model_dict)
        push!(w, _missing_ref_message(key, ref) * " — evaluating the model raises " *
                 "UndefVarError naming $ref but not the key using it")
    end

    # The keys the parser ignores, as `dict_to_model` finds them, without its log message.
    ignored = try
        Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
            _warn_unrecognised_keys(Dict{String,Any}(model_dict))
        end
    catch
        String[]
    end
    # `$WL` is in microns. Divided by, or set against, a constant below 1e-3 it was almost
    # certainly written in metres, and is off by 10⁶ (raised to whatever power follows).
    for k in sort(collect(keys(model_dict)))
        v = model_dict[k]
        v isa AbstractString || continue
        for m in eachmatch(r"\$WL\s*[-/<>]=?\s*\(?\s*([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)", v)
            c = tryparse(Float64, m[1])
            c !== nothing && 0 < c < 1e-3 &&
                push!(w, "$k compares \$WL with $(m[1]): \$WL is in microns, so that constant " *
                         "looks like metres (write $(round(c * 1e6; sigdigits = 6)))")
        end
    end

    for k in sort(ignored)
        cn, suffix = split(k, ',', limit = 2)
        kind = _identify_kind(String(cn), model_dict)
        if suffix in _GEOMETRY_KEYS
            from = get(_KIND_GEOMETRY, kind, "")
            push!(w, "$cn is read as $kind" * (isempty(from) ? "" : " (from $cn,$from)") *
                     "; $k is ignored")
        else
            # Not a shape key, so the kind is unchanged: a parameter added and not referenced
            # yet, or a misspelling. Either way what is missing is a reference to it.
            push!(w, "$k is not used: no expression refers to \$$k")
        end
    end
    return w
end

# The fluxes, which V is divided by the sum of: numeric ones summing to zero or less leave V
# undefined, and every flux free leaves their common scale undetermined.
function _flux_warnings(model_dict::AbstractDict, free::AbstractVector{<:AbstractString})
    w = String[]
    comps = try
        String.(_component_names(Dict{String,Any}(model_dict)))
    catch
        String[]
    end
    isempty(comps) && return w
    fkeys = [c * ",f" for c in comps]
    vals = [get(model_dict, k, 1.0) for k in fkeys]          # no `f` is a flux of 1
    if all(v -> v isa Number, vals) && sum(vals) <= 0
        push!(w, "the fluxes sum to $(sum(vals)): the visibilities are divided by the total " *
                 "flux, which must be positive")
    end
    if length(fkeys) > 1 && all(k -> k in free, fkeys)
        push!(w, "every flux is free ($(join(fkeys, ", "))): the visibilities are divided by " *
                 "their sum, so only their ratios are measured — fix one, or write one as " *
                 "1 − the others")
    end
    return w
end

"""
    display_model(model_dict, list_free_params; lb, ub)

Print the model by component: each parameter's value or expression, free or fixed, and its
bounds, with warnings for values outside their bounds, empty bounds, fluxes, and references to
missing parameters (see [`model_warnings`](@ref)).
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
        for gk in ["image", "profile", "diamin", "fwhmin", "crin",
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

    # The fluxes: V is divided by their sum
    for msg in _flux_warnings(model_dict, list_free_params)
        println("⚠ ", msg)
        n_warnings += 1
    end

    # An expression naming a parameter the model does not define raises at evaluation time,
    # out of the resolver's own module, so the name of the key at fault is only available here.
    for (key, ref) in missing_refs(model_dict)
        println("⚠ ", _missing_ref_message(key, ref))
        n_warnings += 1
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
    # the suffixes of `_REF_PATTERN`, an azimuthal mode's key among them
    replace(expr, r"\b([A-Za-z_][A-Za-z0-9_]*,(?:az (?:amp|projang)[0-9]+|[A-Za-z_][A-Za-z0-9_]*))\b" => s"$\1")
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
  [`ModelConstraint`](@ref)s or `(lhs, op, rhs[, tol])` tuples, passed to NLopt as nonlinear
  constraints; a `method` that does not support them is wrapped in `:AUGLAG`. The reported
  `chi2` is the data χ²
- `rv` — radial velocities ([`read_rv`](@ref)) fitted with the interferometry, against the
  orbit of the model's component on an orbit; their χ² is added, and `ndof` counts them
- `method` — NLopt algorithm (default `:LD_LBFGS`; use `:LN_NELDERMEAD` for gradient-free)
- `maxeval` — maximum function evaluations (default 2000)
- `ftol_rel`, `xtol_rel` — convergence tolerances
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
    rv           ::Union{Nothing,RVData} = nothing,
    method       = :LD_LBFGS,
    maxeval      = 2000,
    ftol_rel     = 1e-8,
    xtol_rel     = 1e-6,
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
    rvix = rv === nothing ? nothing : _rv_index(fm, _orbit_component(model_dict))
    rvchi2(x) = sum(abs2, rv_residuals(fm, x, rv, rvix))

    # ── 4. Detect gradient need ──────────────────────────────────────────────
    use_grad = startswith(string(method), "LD_")

    # ── 5. Build objective ───────────────────────────────────────────────────
    n_evals = Ref(0)
    jkw = _jac_cache_kw(data)
    function objective(x, g)
        n_evals[] += 1
        if use_grad
            chi2, grad = chi2_flat_fg(fm, x, data; weights, verb, jkw...)
            chi2 += prior_penalty_and_grad!(grad, x, priors, fm, prior_idx)
            if rv !== nothing
                chi2 += rvchi2(x)
                grad .+= ForwardDiff.gradient(rvchi2, x)
            end
            # NLopt hands back a zero-length array when it wants the value without the
            # gradient, which AUGLAG does at some of its evaluations.
            length(g) > 0 && (g .= grad)
        else
            chi2 = chi2_flat(fm, x, data; weights, verb)
            chi2 += prior_penalty(x, priors, fm, prior_idx)
            rv === nothing || (chi2 += rvchi2(x))
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

    ndof = _ndof(data, weights) + (rv === nothing ? 0 : length(rv))
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
    jkw = _jac_cache_kw(data)
    function objective(x, g)
        n_evals[] += 1
        if use_grad
            chi2, grad = chi2_flat_fg(fm, x, data; weights, verb, jkw...)
            # NLopt hands back a zero-length array when it wants the value without the gradient.
            length(g) > 0 && (g .= grad)
        else
            chi2 = chi2_flat(fm, x, data; weights, verb)
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
| `backend` | `Symbol` | the sampler: `:nautilus` or `:ultranest` |

`result` is backend-specific: under `:ultranest` a `PythonCall.Py` wrapping UltraNest's result
dict (index with `result["key"]`, convert with `pyconvert(T, …)`), under `:nautilus` the sampler
state Nautilus returned.
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
#                                          likelihood call. See src/fit_model_nautilus.jl.
#   :ultranest        UltraNest         -- Python, reached through PythonCall
#
# `using` either one activates its extension, which registers the backend and adds a method to
# `_fit_nested`. Two independent implementations let an evidence estimate be cross-checked.

"""
The nested sampler [`fit_model_nested`](@ref) uses when no `backend` keyword is given.

`:none` until an extension registers one; the first loaded is the default. Change it with
[`set_nested_backend!`](@ref), which checks that the backend is loaded.
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

The backend's package must be loaded; otherwise this raises an error.
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

Fit a parametric model by nested sampling: posterior samples, log-evidence and the
maximum-likelihood point. The starting values are not used. Every free parameter needs finite
bounds, which define the uniform prior.

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
| `lb`, `ub` | | `Dict{String,Float64}` bounds, required for every free parameter |
| `weights` | `[1,1,1,0,0,0,0]` | observable weights `[V2, T3amp, T3phi, visamp, visphi, flux, diffphase]` |
| `constraints` | `[]` | [`ModelConstraint`](@ref)s or `(lhs, op, rhs[, tol])` tuples, added to the log-likelihood as `-penalty/2` |
| `nB_workspace` | from the data | Hankel workspace size |
| `verb` | `true` | print progress |
| `cornerplot` | `true` | draw the corner plot when a plotting backend is loaded |
| `seed` | `nothing` | reproducible run (UltraNest: seeds numpy's generator, restored afterwards, and uses fixed batches of 512 proposals) |

# Keywords, `:nautilus`

| keyword | default | meaning |
|---|---|---|
| `n_live` | `500` | live points in the exploration phase |
| `n_eff` | `2000` | effective sample size to reach before stopping |
| `f_live` | `0.01` | fraction of evidence left in the live set at which exploration ends |
| `n_networks` | `4` | networks in the bounding ensemble |
| `threaded` | `nthreads() > 1` | evaluate the likelihood across threads |

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

`:nautilus` is importance nested sampling with neural-network bounds. On the α Cen A and B
power-law fits its `logzerr` is 0.019, against 0.18–0.31 for UltraNest, with `x_opt` equal to
four decimals.
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
    converged   ::Bool                  # stopped by a convergence test, not the iteration cap
    model       ::FlatModel
    lsqfit_result ::LsqFit.LsqFitResult  # in LsqFit's form, so confidence_interval etc. apply;
                                         # `length(trace)` is the iteration count
    reason      ::Symbol                # what stopped it: :chi2, :gradient, :step, :bounds, :maxiter
    at_bound    ::Vector{String}        # parameters the minimum holds on a bound
    chi2_rv     ::Float64               # the radial velocities' part of chi2 (0 without them)
end

function Base.show(io::IO, r::LsqFitResult)
    @printf(io, "LsqFitResult: χ²r = %.4f  (chi2=%.2f, ndof=%d, converged=%s: %s)\n",
            r.chi2r, r.chi2, r.ndof, r.converged, _LM_REASONS[r.reason])
    isempty(r.at_bound) || println(io, "  on a bound: ", join(r.at_bound, ", "))
    r.chi2_rv > 0 && @printf(io, "  χ² interferometry %.2f, radial velocities %.2f\n",
                             r.chi2 - r.chi2_rv, r.chi2_rv)
    for (i, p) in enumerate(r.list_free_params)
        @printf(io, "  %-20s = %10.4f ± %.4f\n", p, r.x_opt[i], r.stderror[i])
    end
end


# ─────────────────────────────────────────────────────────────────────────────
# fit_model_lsqfit — Levenberg-Marquardt (`_lm_bounded`), reported in LsqFit's form
# ─────────────────────────────────────────────────────────────────────────────

# Run `_lm_bounded` and wrap it as LsqFit's result, with no weights: `vcov` and `stderror` then
# scale by rss/(n − p) exactly as they did when `curve_fit` produced it.
# The residual and Jacobian functions the Levenberg-Marquardt loop calls: the data's rows
# written in place (`residuals_flat!`, `residuals_flat_jac!`), then the extra rows of each
# `extra` function (radial velocities, constraints; `nothing` for none), into buffers allocated
# once per fit. The residuals and the Jacobian have buffers of their own, since the loop
# evaluates both before it copies either; it copies each into its own storage, so the next call
# may overwrite them.
function _lm_residual_fns(fm, data, weights, x0, ndata, extra, extra_jac, jkw; ref = nothing)
    nx = [e === nothing ? 0 : length(e(x0)) for e in extra]
    n, np = ndata + sum(nx; init = 0), length(x0)
    rf = Vector{Float64}(undef, n)
    rj = Vector{Float64}(undef, ndata); J = Matrix{Float64}(undef, n, np)
    function model_fn(p)
        residuals_flat!(view(rf, 1:ndata), fm, p, data; weights, ref)
        k = ndata
        for (e, m) in zip(extra, nx)
            e === nothing && continue
            copyto!(view(rf, (k + 1):(k + m)), e(p)); k += m
        end
        return rf
    end
    function jacobian_fn(p)
        residuals_flat_jac!(rj, view(J, 1:ndata, :), fm, p, data; weights, ref, jkw...)
        k = ndata
        for (e, m) in zip(extra_jac, nx)
            e === nothing && continue
            copyto!(view(J, (k + 1):(k + m), :), e(p)); k += m
        end
        return J
    end
    return model_fn, jacobian_fn
end

function _lm_fit(model_fn, jacobian_fn, x0, lb_vec, ub_vec, maxIter)
    lm = _lm_bounded(model_fn, jacobian_fn, Float64.(x0), lb_vec, ub_vec; maxIter)
    fit = LsqFit.LsqFitResult(lm.x, lm.resid, lm.jacobian, lm.converged, lm.trace, Float64[])
    return fit, lm
end

"""
    fit_model_lsqfit(model_dict, list_free_params, data; kwargs...) -> LsqFitResult

Fit a parametric model by Levenberg-Marquardt. Returns the covariance and 1σ uncertainties,
both scaled to χ²ᵣ = 1.

A parameter that reaches a bound is held there and left out of the step. Besides LsqFit's
gradient and step tests, the iteration stops when the Gauss-Newton step would lower χ² by less
than `0.001 · max(1, χ²/n)` (n data points). `reason` gives the test that stopped it, and `at_bound` the
parameters held on a bound. The Jacobian is analytic (`residuals_flat_jac`).

# Arguments
- `model_dict::Dict{String}` — flat parameter dictionary
- `list_free_params::Vector{String}` — names of the free parameters to optimize
- `data::OIdata` — interferometric data

# Keywords
- `lb`, `ub` — `Dict{String,Float64}` of lower/upper bounds per parameter
- `weights` — observable weights [V2, T3amp, T3phi, visamp, visphi, flux, diffphase]
- `nB_workspace` — Hankel workspace size
- `maxIter` — maximum iterations (default 200)
- `constraints` — relations between parameters that bounds cannot express, as
  [`ModelConstraint`](@ref)s or `(lhs, op, rhs[, tol])` tuples, added as extra rows of the
  residual vector and the Jacobian, so `covar` and `stderror` include them. The reported `chi2`
  excludes them
- `measurement` — the measurement model ([`measurement_model`](@ref)); with `visphi =
  :continuum`, each spectral vector's VISPHI residuals are referenced to a polynomial fitted over
  its continuum channels, as PMOIRED's DPHI (ndof reduced by the coefficients)
- `rv` — radial velocities ([`read_rv`](@ref)) fitted with the interferometry, against the
  orbit of the model's component on an orbit (`c,orb Ka`, `Kb`, `gamma`, or `M` and `q`). They
  enter the residual vector after the data; `chi2` and `ndof` count them, `chi2_rv` is their part
- `verb` — print per-evaluation chi2 breakdown
"""
function fit_model_lsqfit(model_dict   ::Dict{String},
                          list_free_params   ::Vector{String},
                          data         ::OIdata;
    lb           = Dict{String,Float64}(),
    ub           = Dict{String,Float64}(),
    weights      = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    constraints  = [],
    rv           ::Union{Nothing,RVData} = nothing,
    measurement  ::Union{Nothing,MeasurementModel} = nothing,
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
    rvix = rv === nothing ? nothing : _rv_index(fm, _orbit_component(model_dict))

    x0 = Float64[model_dict[p] for p in list_free_params]
    lb_vec = Float64[get(lb, p, -Inf) for p in list_free_params]
    ub_vec = Float64[get(ub, p,  Inf) for p in list_free_params]
    return _lsqfit_compiled(fm, list_free_params, x0, data, lb_vec, ub_vec; weights, cons,
                            constraint_idx, rv, rvix, maxIter, verb, measurement)
end

# The Dict form's fit on a compiled model (`fm`, with the constraint expressions
# `constraint_idx` already in it): what that form runs after compiling, and what a caller fitting
# many starts of one model runs per start. `warn_cov = false` keeps a failed covariance silent
# (its σ are NaN either way).
function _lsqfit_compiled(fm, list_free_params, x0::Vector{Float64}, data, lb_vec, ub_vec;
                          weights, cons, constraint_idx, rv, rvix, maxIter, verb,
                          warn_cov::Bool = true, measurement = nothing)
    nparams = length(list_free_params)
    # ── 2. The residuals and their Jacobian, the data captured in the closures ──
    # The data rows, the radial velocities', then the constraints'. Their count is fixed for a
    # given model — one per constraint, or one per wavelength where the resolver broadcasts —
    # so the residual vector keeps the length it has at the start.
    ndata = length(residuals_flat(fm, x0, data; weights))
    rv_fn  = rv === nothing ? nothing : p -> rv_residuals(fm, p, rv, rvix)
    rv_jfn = rv === nothing ? nothing : p -> rv_jacobian(fm, p, rv, rvix)
    c_fn   = isempty(cons) ? nothing : p -> constraint_residuals(p, cons, fm, constraint_idx)
    c_jfn  = isempty(cons) ? nothing : p -> constraint_jacobian(p, cons, fm, constraint_idx)
    ref = _references(measurement) ? _ref_context(measurement, fm, data) : nothing
    model_fn, jacobian_fn = _lm_residual_fns(fm, data, weights, x0, ndata,
                                             (rv_fn, c_fn), (rv_jfn, c_jfn), _jac_cache_kw(data); ref)

    # ── 3. Levenberg-Marquardt ───────────────────────────────────────────────
    fit, lm = _lm_fit(model_fn, jacobian_fn, x0, lb_vec, ub_vec, maxIter)

    # ── 4. Extract results ───────────────────────────────────────────────────
    minx = fit.param
    # The data rows (radial velocities included): a chi2 carrying the constraint penalty would
    # not be comparable with what any other fitter reports, nor with a published chi2r.
    nrv = rv === nothing ? 0 : length(rv)
    r_final = view(fit.resid, 1:(ndata + nrv))
    chi2_val = dot(r_final, r_final)
    chi2_rv = nrv == 0 ? 0.0 : sum(abs2, view(fit.resid, (ndata + 1):(ndata + nrv)))
    # each spectral vector's reference coefficients are fitted: out of the degrees of freedom
    ncoef = ref === nothing ? 0 :
            _apply_reference!(copy(fit.resid[1:ndata]), nothing, ref, fm, minx, data, _pad_weights(weights))
    ndof = _ndof(data, weights) + nrv - ncoef

    # Covariance matrix from LsqFit (may fail if Jacobian is singular,
    # e.g. when parameters hit bounds)
    cov, σ = try
        vcov(fit), stderror(fit)
    catch
        warn_cov && @warn "Covariance estimation failed (singular Jacobian — parameters may be at bounds)"
        fill(NaN, nparams, nparams), fill(NaN, nparams)
    end

    if verb
        chi2_flat(fm, minx, data; weights, verb=true)
    end

    return LsqFitResult(minx, list_free_params, chi2_val, chi2_val / ndof, ndof,
                         cov, σ, fit.converged, fm, fit, lm.reason,
                         list_free_params[lm.held], chi2_rv)
end

"""
    fit_model_lsqfit(model::FlatModel, x0, data; kwargs...) -> LsqFitResult

Levenberg-Marquardt fit on a pre-compiled `FlatModel` starting from `x0`.
Bounds `lb` and `ub` can be Dicts or Vectors. Radial velocities `rv` need `rv_component`, the
component on an orbit.
"""
function fit_model_lsqfit(model        ::FlatModel,
                          x0           ::AbstractVector{<:Real},
                          data         ::OIdata;
    lb           = Dict{String,Float64}(),
    ub           = Dict{String,Float64}(),
    weights      = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    constraints  = [],
    rv           ::Union{Nothing,RVData} = nothing,
    rv_component ::AbstractString = "",
    measurement  ::Union{Nothing,MeasurementModel} = nothing,
    maxIter      = 200,
    verb         = false,
)
    list_free_params = model.list_free_params
    fm = model
    nparams = length(list_free_params)
    # The radial velocities' orbit: the component is named, the compiled model has no dict.
    rvix = rv === nothing ? nothing : _rv_index(fm, rv_component)

    # A constraint is evaluated as an extra parameter of the compiled model, so it needs the
    # model dict to augment and recompile. Accepting it here and dropping it would return a fit
    # that ignores what the caller asked for.
    !isempty(constraints) && error("Constraints are only supported with the Dict form of fit_model_lsqfit; pass model_dict instead of FlatModel")

    x0f = Float64.(x0)
    ndata = length(residuals_flat(fm, x0f, data; weights))
    rv_fn  = rv === nothing ? nothing : p -> rv_residuals(fm, p, rv, rvix)
    rv_jfn = rv === nothing ? nothing : p -> rv_jacobian(fm, p, rv, rvix)
    ref = _references(measurement) ? _ref_context(measurement, fm, data) : nothing
    model_fn, jacobian_fn = _lm_residual_fns(fm, data, weights, x0f, ndata,
                                             (rv_fn,), (rv_jfn,), _jac_cache_kw(data); ref)

    lb_vec = _bounds_vec(lb, list_free_params, -Inf)
    ub_vec = _bounds_vec(ub, list_free_params,  Inf)

    fit, lm = _lm_fit(model_fn, jacobian_fn, x0, lb_vec, ub_vec, maxIter)

    minx = fit.param
    r_final = fit.resid
    chi2_val = dot(r_final, r_final)
    nrv = rv === nothing ? 0 : length(rv)
    chi2_rv = nrv == 0 ? 0.0 : sum(abs2, view(r_final, (length(r_final) - nrv + 1):length(r_final)))
    ncoef = ref === nothing ? 0 :
            _apply_reference!(copy(r_final[1:ndata]), nothing, ref, fm, minx, data, _pad_weights(weights))
    ndof = _ndof(data, weights) + nrv - ncoef

    cov, σ = try
        vcov(fit), stderror(fit)
    catch
        @warn "Covariance estimation failed (singular Jacobian — parameters may be at bounds)"
        fill(NaN, nparams, nparams), fill(NaN, nparams)
    end

    if verb
        chi2_flat(fm, minx, data; weights, verb=true)
    end

    return LsqFitResult(minx, list_free_params, chi2_val, chi2_val / ndof, ndof,
                         cov, σ, fit.converged, fm, fit, lm.reason,
                         list_free_params[lm.held], chi2_rv)
end


# ─────────────────────────────────────────────────────────────────────────────
# model_to_obs — compute observables from a FlatModel
# ─────────────────────────────────────────────────────────────────────────────

"""
    model_to_obs(model, x, data) -> (v2, t3amp, t3phi, visamp, visphi)

The model observables at parameter vector `x`, matching the data: a NamedTuple with fields `v2`, `t3amp`, `t3phi` (degrees),
`visamp`, `visphi` (degrees).  Empty vectors are returned for
observable types not present in `data`.
"""
function model_to_obs(model::FlatModel, x::AbstractVector, data::OIdata)
    # Per-uv wavelength and MJD, so chromatic expressions (`$WL`, `$MJD`) resolve. Without
    # them the resolver is handed NaN and every observable comes back NaN -- silently, since a
    # NaN model is not an error. `chi2_flat` passes the same pair for the same reason.
    # A smeared model: the observables averaged over each channel (`eval_model_sub`).
    _smeared(model, data) && return _model_to_obs_smeared(model, x, data)
    cvis = eval_model(model, x, data.uv;
                      wl  = hasproperty(data, :uv_lam) ? data.uv_lam : nothing,
                      mjd = hasproperty(data, :uv_mjd) ? data.uv_mjd : nothing)

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

# `model_to_obs` for a smeared model: the observables averaged over each channel.
function _model_to_obs_smeared(model::FlatModel, x::AbstractVector, data::OIdata)
    V = first(_smeared_cvis(model, x, data, false))
    v2, t3, vis = _smeared_observables(V isa AbstractMatrix ? V : reshape(V, :, 1), data)
    e(n, y) = n > 0 ? y : Float64[]
    return (v2 = e(data.nv2, v2), t3amp = e(data.nt3amp + data.nt3phi, abs.(t3)),
            t3phi = e(data.nt3amp + data.nt3phi, angle.(t3) .* (180.0 / π)),
            visamp = e(data.nvisamp + data.nvisphi, abs.(vis)),
            visphi = e(data.nvisamp + data.nvisphi, angle.(vis) .* (180.0 / π)))
end

"""
    model_to_residuals(model, x, data)

Normalised residuals `(model - data) / error` as a NamedTuple with fields `v2`, `t3amp`,
`t3phi`, `visamp`, `visphi`. Phase differences are wrapped to [-180, 180] before division.
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

The spectral energy distribution of a `FlatModel`: the total flux and each component's flux
(its `f` / `spectrum` resolved) at each wavelength.

# Arguments
- `model::FlatModel` — compiled model from `dict_to_model`
- `x::AbstractVector` — current free-parameter values
- `wl_grid` — wavelengths in metres, as `uv_lam` and every `wl` argument; expressions see
  `\$WL` in microns.

# Returns
- `total::Vector{T}` — total flux at each wavelength, Σ f ([`model_flux`](@ref)), in the
  parameters' type
- `components::Dict{String,Vector{T}}` — per-component flux, in the same type
"""
function model_to_sed(model::FlatModel,
                      x::AbstractVector,
                      wl_grid::AbstractVector)
    nwl = length(wl_grid)

    # Total flux, Σ f, at each wavelength, in the parameters' type
    total = model_flux(model, x; wl = wl_grid)
    T = eltype(total)

    # Per-component fluxes from the resolver
    comp_fluxes = Dict{String, Vector{T}}()
    for comp in model.components
        comp_fluxes[comp.comp_name] = zeros(T, nwl)
    end
    for (i, wl) in enumerate(wl_grid)
        pv = model.resolver.fn(x, wl, NaN)
        for comp in model.components
            f = comp.f_idx == 0 ? one(T) : T(pv[comp.f_idx])
            comp_fluxes[comp.comp_name][i] = f
        end
    end

    return total, comp_fluxes
end


"""
    model_to_flux(model, x; wl=nothing)

Total flux of the model at zero baseline: `real(V(0,0))`.

For a chromatic model, `wl` (metres) is a scalar or a 1-element vector. See also
[`model_to_sed`](@ref).
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

The model visibility is evaluated on the FFT frequency grid (half of it, by Hermitian symmetry)
and transformed by an inverse real FFT.

# Arguments
- `model::FlatModel` — compiled model from `dict_to_model`
- `x::AbstractVector` — current free-parameter values

# Keywords
- `nx` — image size in pixels (default 256)
- `pixsize` — pixel size in mas (default 0.1)
- `oversample` — oversampling factor (default 1)
- `normalize` — normalize image to unit sum (default true)
- `wl` — wavelength in metres (scalar), as every `wl` argument; enables `\$WL` (microns in
  the expressions) in spectrum expressions

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
    # Not normalised by the flux when the image is not: its sum is then the model's flux.
    V_half = eval_model(model, x, uv; wl, mjd, normalize)

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
