# resolvers.jl
#
# Flat-dict parameter expression resolver for model fitting: free parameters in, every named
# parameter (free, fixed, derived) out, compiled once per model with RuntimeGeneratedFunctions.
#
# Interface:
#   resolver = RGF.build_resolver(model_dict, list_free_params)
#   resolver.fn(x, WL, MJD)               -> Vector of every value, by `all_names`
#   resolver.fill!(s, cols, x, WL, MJD)   -> writes them into typed storage (see below)

using RuntimeGeneratedFunctions

# ─────────────────────────────────────────────────────────────────────────────
# Shared utilities (used by all implementations)
# ─────────────────────────────────────────────────────────────────────────────

# Reserved implicit variables available in expressions (not parameter names).
# $R, $MU              — radial grid variables (profile expressions)
# $D                   — diameter at each grid radius, i.e. 2$R (profile expressions)
# $RMIN, $RMAX         — the grid's own bounds, in mas (profile expressions)
# $DMIN, $DMAX         — the same as diameters (profile expressions)
# $WL                  — wavelength in metres (per UV point, from data.uv_lam)
# $MJD                 — Modified Julian Date  (per UV point, from data.uv_mjd, which is
#                        stored in T: ~5.6 min resolution at the default Float32)
#
# The five grid variables are PMOIRED's, spelled the same way, so a profile written for one
# package reads in the other. They are substituted where the profile is compiled, by which
# point the grid is known.
const IMPLICIT_VARS = Set(["R", "MU", "D", "RMIN", "RMAX", "DMIN", "DMAX", "WL", "MJD"])

# A `$`-reference: a name, then optionally `,` and a suffix. The suffix is an identifier, or
# one of the two keys written with a space, PMOIRED's spelling of an azimuthal mode
# (`az amp1`, `az projang1`). No other space belongs to a reference: `$a,b c` is `$a,b`, then c.
const _REF_PATTERN = r"\$([A-Za-z_][A-Za-z0-9_]*(?:,(?:az (?:amp|projang)[0-9]+|[A-Za-z_][A-Za-z0-9_]*))?)"

# `$WL` in an expression is in MICRONS, the unit the field writes a wavelength in; every
# wavelength outside the expression language (the `wl` arguments, `uv_lam`, `model_to_sed`'s
# grid) is in metres. A compiled expression reads `$WL` as `_wl_expr(WL)`, in the resolver and
# in a visibility function, and nowhere else: inside the generated code it is part of the
# fused broadcast, so it allocates nothing and keeps the wavelength's own precision.
_wl_expr(wl::Real) = wl * oftype(float(wl), 1e6)

# A reference as a Julia identifier: `,` → `__`, the space of an azimuthal mode's key → `_`.
_mangle(name::AbstractString) = replace(name, "," => "__", " " => "_")

"""
    extract_refs(expr_str) -> Vector{String}

Return all \$param,name references in an expression string, as bare
"param,name" strings (dollar stripped).  Implicit variables (WL, MJD,
R, MU) are included — filter with `IMPLICIT_VARS` if needed.
"""
function extract_refs(expr_str::String)::Vector{String}
    refs = String[]
    for m in eachmatch(_REF_PATTERN, expr_str)
        push!(refs, m.captures[1])
    end
    return unique(refs)
end

"""
    topo_sort(derived::Dict{String,String}) -> Vector{String}

Topological sort of derived parameter expressions.
Returns the names in evaluation order (dependencies first).
"""
function topo_sort(derived::Dict{String,String})::Vector{String}
    # Build dependency graph
    deps = Dict{String,Vector{String}}()
    for (name, expr) in derived
        refs = extract_refs(expr)
        deps[name] = filter(r -> haskey(derived, r), refs)
    end
    # Kahn's algorithm
    in_deg = Dict(n => length(d) for (n, d) in deps)
    queue  = [n for (n, d) in in_deg if d == 0]
    result = String[]
    while !isempty(queue)
        n = popfirst!(queue)
        push!(result, n)
        for (m, d) in deps
            if n in d
                in_deg[m] -= 1
                in_deg[m] == 0 && push!(queue, m)
            end
        end
    end
    length(result) == length(derived) ||
        error("Cycle detected in derived parameter expressions: " *
              "$(setdiff(keys(derived), result))")
    return result
end

"""
    partition_dict(model_dict, list_free_params) -> (free, fixed, derived)

Split a flat parameter dictionary into:
  free    — Dict of parameters that will be optimised
  fixed   — Dict of numeric parameters that are held constant
  derived — Dict of string expressions (to be evaluated in topo order)
"""
function partition_dict(model_dict::Dict{String}, list_free_params)
    fit_set = Set(list_free_params)
    free    = Dict{String,Float64}()
    fixed   = Dict{String,Float64}()
    derived = Dict{String,String}()
    for (k, v) in model_dict
        if k in fit_set
            v isa String && error("Parameter '$k' is in list_free_params but has a formula value (\"$v\"); derived parameters cannot be optimized — remove it from list_free_params or give it a numeric starting value")
            free[k] = Float64(v)
        elseif v isa Number
            fixed[k] = Float64(v)
        elseif v isa String
            derived[k] = v
        else
            # Non-string, non-numeric (e.g. Bool markers like "resolved" => true)
            # Skip — these are handled by dict_to_model, not the resolver
        end
    end
    return free, fixed, derived
end

"""
    broadcastify(expr) -> Expr

Transform a Julia AST so that all `:call` nodes become `:.` (dot-broadcast).
This allows scalar expressions to work element-wise on vector arguments
like \$WL and \$MJD.
"""
function broadcastify(expr)
    if expr isa Expr && expr.head == :call
        fn   = expr.args[1]
        args = [broadcastify(a) for a in expr.args[2:end]]
        return Expr(:., fn, Expr(:tuple, args...))
    elseif expr isa Expr
        return Expr(expr.head, [broadcastify(a) for a in expr.args]...)
    else
        return expr
    end
end


##############################################################################
# The resolver: RuntimeGeneratedFunctions.jl
##############################################################################
module RGF

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

using ..OITOOLS: extract_refs, topo_sort, partition_dict, broadcastify, IMPLICIT_VARS,
                _REF_PATTERN, _mangle, _wl_expr

const mangle = _mangle

function preprocess(s::String)::String
    replace(s, _REF_PATTERN => s -> (n = s[2:end]; n == "WL" ? "_wl_expr(WL)" : mangle(n)))
end

"""
The compiled resolver of one model.

- `fn(x, WL, MJD)` returns every parameter's value, ordered as `all_names`: a scalar, or one
  value per uv point for a parameter whose expression depends on `\$WL`/`\$MJD`.
- `fill!(s, cols, x, WL, MJD)` writes the same values into typed storage: the scalars into
  `s` (a `Vector{T}`), the per-point ones into `cols` (a `Vector` of preallocated
  `Vector{T}`, one per uv point each). `slot[i]` says where `all_names[i]` went: `k > 0` for
  `s[k]`, `-j` for `cols[j]`.

`fill!` is what `eval_model` calls, so its result is concretely typed whatever the model;
`fn` serves constraints and priors, which index single values.

A parameter is per-point ("chromatic") when its expression mentions `\$WL` or `\$MJD`, or
refers to a chromatic parameter; the rest stay scalar, so a chromatic flux does not make a
fixed diameter a vector.
"""
struct Resolver
    all_names ::Vector{String}
    fn        ::Function
    fill!     ::Function
    slot      ::Vector{Int}
    nscalar   ::Int
    nchromatic::Int
end

function build_resolver(model_dict::Dict{String}, list_free_params::Vector{String})
    free, fixed, derived = partition_dict(model_dict, list_free_params)
    sorted = topo_sort(derived)

    # 1. Define the names that will be returned in the output vector
    all_names = [list_free_params; collect(keys(fixed)); sorted]

    # 2. Which derived expressions depend on WL or MJD, directly or through another derived
    #    parameter (in topological order, so a dependency is classified first).
    direct(v) = any(r -> r in IMPLICIT_VARS && r ∉ ("R", "MU"), extract_refs(v))
    chromatic = Set{String}()
    for name in sorted
        (direct(derived[name]) || any(r -> r in chromatic, extract_refs(derived[name]))) &&
            push!(chromatic, name)
    end
    uses_wl_mjd = !isempty(chromatic)

    # 3. The body shared by both functions: free parameters from x, fixed ones as literals.
    head = [[:($(Symbol(mangle(list_free_params[i]))) = x[$i]) for i in eachindex(list_free_params)];
            [:($(Symbol(mangle(name))) = $val) for (name, val) in fixed]]

    # `fn`: every derived expression broadcast when any is chromatic, as vector WL/MJD then
    # propagate; plain scalar code otherwise.
    fn_stmts = if uses_wl_mjd
        [:($(Symbol(mangle(name))) = $(broadcastify(Meta.parse(preprocess(derived[name])))))
         for name in sorted]
    else
        [:($(Symbol(mangle(name))) = $(Meta.parse(preprocess(derived[name]))))
         for name in sorted]
    end
    fn_body = quote
        $(head...)
        $(fn_stmts...)
        return $(Expr(:vect, [Symbol(mangle(n)) for n in all_names]...))
    end
    compiled_fn = @RuntimeGeneratedFunction(:((x, WL, MJD) -> $fn_body))
    # Wrap with default NaN for backward compatibility: fn(x) still works. The wavelength
    # arrives in metres; the expressions read it in microns (`preprocess`).
    wrapped_fn(x, WL=NaN, MJD=NaN) = compiled_fn(x, WL, MJD)

    # `fill!`: chromatic expressions broadcast in place into their column; the others scalar.
    slot = zeros(Int, length(all_names))
    ns = 0; nc = 0
    for (i, n) in enumerate(all_names)
        if n in chromatic
            nc += 1; slot[i] = -nc
        else
            ns += 1; slot[i] = ns
        end
    end
    fill_stmts = Expr[]
    for name in sorted
        sym = Symbol(mangle(name))
        ex  = Meta.parse(preprocess(derived[name]))
        if name in chromatic
            j = -slot[findfirst(==(name), all_names)]
            push!(fill_stmts, :($sym = cols[$j]))
            push!(fill_stmts, :($sym .= $(broadcastify(ex))))
        else
            push!(fill_stmts, :($sym = $ex))
        end
    end
    stores = [:(s[$(slot[i])] = $(Symbol(mangle(n))))
              for (i, n) in enumerate(all_names) if slot[i] > 0]
    fill_body = quote
        $(head...)
        $(fill_stmts...)
        $(stores...)
        return nothing
    end
    compiled_fill = @RuntimeGeneratedFunction(:((s, cols, x, WL, MJD) -> $fill_body))

    Resolver(all_names, wrapped_fn, compiled_fill, slot, ns, nc)
end

function (r::Resolver)(x::Vector{Float64}, wl=NaN, mjd=NaN)
    vals_vec = r.fn(x, wl, mjd)
    return Dict{String,Any}(zip(r.all_names, vals_vec))
end

end # module RGF
