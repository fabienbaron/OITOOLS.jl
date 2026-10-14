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
# $MJD                 — Modified Julian Date  (per UV point, from data.uv_mjd, Float64
#                        whatever the data's T)
#
# The five grid variables are PMOIRED's, spelled the same way, so a profile written for one
# package reads in the other. They are substituted where the profile is compiled, by which
# point the grid is known.
const IMPLICIT_VARS = Set(["R", "MU", "D", "RMIN", "RMAX", "DMIN", "DMAX", "WL", "MJD"])

# A `$`-reference: a name (PMOIRED's components are often `1`, `2`), then optionally `,` and a
# suffix. The suffix is an identifier, or one of PMOIRED's keys written with a space: an
# azimuthal mode (`az amp1`, `az projang1`) or an orbit's element (`orb P`). No other space
# belongs to a reference: `$a,b c` is `$a,b`, then c.
const _REF_PATTERN = r"\$([A-Za-z0-9_][A-Za-z0-9_]*(?:,(?:az (?:amp|projang)[0-9]+|orb [A-Za-z_][A-Za-z0-9_]*|[A-Za-z_][A-Za-z0-9_]*))?)"

# `$WL` in an expression is in MICRONS, the unit the field writes a wavelength in; every
# wavelength outside the expression language (the `wl` arguments, `uv_lam`, `model_to_sed`'s
# grid) is in metres. A compiled expression reads `$WL` as `_wl_expr(WL)`, in the resolver and
# in a visibility function, and nowhere else: inside the generated code it is part of the
# fused broadcast, so it allocates nothing and keeps the wavelength's own precision.
_wl_expr(wl::Real) = wl * oftype(float(wl), 1e6)

# A reference as a Julia identifier: `,` → `__`, the space of a two-word key → `_`, and a `_`
# before a name that starts with a digit (`2,orb P` → `_2__orb_P`).
_mangle(name::AbstractString) =
    (m = replace(name, "," => "__", " " => "_"); isdigit(first(m)) ? "_" * m : m)

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
                _REF_PATTERN, _mangle, _wl_expr, orbit_x, orbit_y, orbit_a

const mangle = _mangle

function preprocess(s::String)::String
    replace(s, _REF_PATTERN => s -> (n = s[2:end]; n == "WL" ? "_wl_expr(WL)" : mangle(n)))
end

# A per-point expression of time alone (`$MJD`, no `$WL`, no other per-point parameter), written
# as a loop that evaluates it as scalar code once per run of equal MJDs and copies it along the
# run: `uv_mjd` comes in runs, one per epoch, so an orbit's Kepler equation is solved per epoch.
# Inline in the generated body (a closure there would be boxed). `epoch_fill` writes into the
# column `sym`; `epoch_vector` returns a new vector (a scalar for a scalar MJD). The loop's own
# names (`#y#<sym>`, …) are fixed by `sym` and cannot be a mangled parameter's, so the same
# model always generates the same body and reuses its compiled code.
function epoch_loop(sym, ex, y)
    i, t = Symbol("#i#", sym), Symbol("#t#", sym)
    quote
        $i = firstindex(MJD)
        @inbounds while $i <= lastindex(MJD)
            $t = MJD[$i]
            $y = let MJD = $t; $ex; end
            $sym[$i] = $y
            $i += 1
            while $i <= lastindex(MJD) && MJD[$i] == $t
                $sym[$i] = $y
                $i += 1
            end
        end
    end
end

function epoch_fill(sym, ex)
    y = Symbol("#y#", sym)
    quote
        if MJD isa Number
            $sym .= $ex
        elseif !isempty(MJD)
            local $y
            $(epoch_loop(sym, ex, y))
        end
    end
end

function epoch_vector(sym, ex)
    y, v = Symbol("#y#", sym), Symbol("#v#", sym)
    quote
        $sym = if MJD isa Number || isempty(MJD)
            $(broadcastify(ex))
        else
            local $y = let MJD = MJD[firstindex(MJD)]; $ex; end
            $v = Vector{typeof($y)}(undef, length(MJD))
            $(epoch_loop(v, ex, y))
            $v
        end
    end
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
fixed diameter a vector. One of time alone (`\$MJD`, no `\$WL`, no chromatic parameter) is
computed once per run of equal MJDs and copied along it.

The fixed parameters' values are data, not code: the generated functions read them from a
vector, which `fn` and `fill!` pass in. The generated code depends only on the model's
structure — its keys, which are free, its expressions — so a model compiled again with other
fixed values (an edit in the parameter table, a grid over a fixed parameter) reuses the
compiled code: compiling it again is parsing, not code generation.
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
    # Those of time alone, evaluated once per run of equal MJDs.
    by_epoch = Set(name for name in chromatic
                   if (refs = extract_refs(derived[name]); "MJD" in refs && "WL" ∉ refs &&
                       !any(r -> r in chromatic, refs)))

    # 3. The body shared by both functions: free parameters from x, fixed ones from c. The
    # generated functions' own arguments are named `#x`, `#c`, `#s`, `#cols`, which no
    # parameter's mangled name can be (a global named `s` or `c` would otherwise overwrite one).
    X, C, S, COLS = Symbol("#x"), Symbol("#c"), Symbol("#s"), Symbol("#cols")
    fixed_names = sort(collect(keys(fixed)))
    cvals = Float64[fixed[n] for n in fixed_names]
    head = [[:($(Symbol(mangle(list_free_params[i]))) = $X[$i]) for i in eachindex(list_free_params)];
            [:($(Symbol(mangle(name))) = $C[$j]) for (j, name) in enumerate(fixed_names)]]

    # `fn`: every derived expression broadcast when any is chromatic, as vector WL/MJD then
    # propagate; plain scalar code otherwise.
    fn_stmts = if uses_wl_mjd
        [(ex = Meta.parse(preprocess(derived[name])); sym = Symbol(mangle(name));
          name in by_epoch ? epoch_vector(sym, ex) : :($sym = $(broadcastify(ex))))
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
    # `fn(x)` works too (WL and MJD default to NaN). The wavelength arrives in metres; the
    # expressions read it in microns (`preprocess`).
    compiled_fn = @RuntimeGeneratedFunction(:(($X, $C, WL, MJD) -> $fn_body))
    wrapped_fn(x, WL = NaN, MJD = NaN) = compiled_fn(x, cvals, WL, MJD)

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
            push!(fill_stmts, :($sym = $COLS[$j]))
            push!(fill_stmts, name in by_epoch ? epoch_fill(sym, ex) : :($sym .= $(broadcastify(ex))))
        else
            push!(fill_stmts, :($sym = $ex))
        end
    end
    stores = [:($S[$(slot[i])] = $(Symbol(mangle(n))))
              for (i, n) in enumerate(all_names) if slot[i] > 0]
    fill_body = quote
        $(head...)
        $(fill_stmts...)
        $(stores...)
        return nothing
    end
    compiled_fill = @RuntimeGeneratedFunction(:(($S, $COLS, $X, $C, WL, MJD) -> $fill_body))

    fill_c(s, cols, x, WL, MJD) = compiled_fill(s, cols, x, cvals, WL, MJD)

    Resolver(all_names, wrapped_fn, fill_c, slot, ns, nc)
end

function (r::Resolver)(x::Vector{Float64}, wl=NaN, mjd=NaN)
    vals_vec = r.fn(x, wl, mjd)
    return Dict{String,Any}(zip(r.all_names, vals_vec))
end

end # module RGF
