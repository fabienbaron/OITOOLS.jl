# resolvers.jl
#
# Flat-dict parameter expression resolver for OITOOLS2 model fitting.
#
# Three implementations sharing a common utility module:
#   1. HandRolled  — custom AST interpreter, zero dependencies
#   2. RGF         — RuntimeGeneratedFunctions.jl, compiled at parse time
#   3. Symbolic    — Symbolics.jl, CAS-based, enables analytic gradients
#
# Interface:
#   resolver = build_resolver(param_dict, fit_params)
#   all_vals = resolver(x::Vector{Float64}) -> Dict{String,Float64}

using RuntimeGeneratedFunctions

# ─────────────────────────────────────────────────────────────────────────────
# Shared utilities (used by all implementations)
# ─────────────────────────────────────────────────────────────────────────────
module SharedUtils

export extract_refs, topo_sort, partition_dict, broadcastify, IMPLICIT_VARS

# Reserved implicit variables available in expressions (not parameter names).
# $R, $MU  — radial grid variables (profile expressions)
# $WL      — wavelength in metres (per UV point, from data.uv_lam)
# $MJD     — Modified Julian Date  (per UV point, from data.uv_mjd)
const IMPLICIT_VARS = Set(["R", "MU", "WL", "MJD"])

"""
    extract_refs(expr_str) -> Vector{String}

Return all \$param,name references in an expression string, as bare
"param,name" strings (dollar stripped).  Implicit variables (WL, MJD,
R, MU) are included — filter with `IMPLICIT_VARS` if needed.
"""
function extract_refs(expr_str::String)::Vector{String}
    refs = String[]
    for m in eachmatch(r"\$([A-Za-z_][A-Za-z0-9_,]*)", expr_str)
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
    partition_dict(param_dict, fit_params) -> (free, fixed, derived)

Split a flat parameter dictionary into:
  free    — Dict of parameters that will be optimised
  fixed   — Dict of numeric parameters that are held constant
  derived — Dict of string expressions (to be evaluated in topo order)
"""
function partition_dict(param_dict::Dict{String}, fit_params)
    fit_set = Set(fit_params)
    free    = Dict{String,Float64}()
    fixed   = Dict{String,Float64}()
    derived = Dict{String,String}()
    for (k, v) in param_dict
        if k in fit_set
            free[k] = Float64(v)
        elseif v isa Number
            fixed[k] = Float64(v)
        elseif v isa String
            derived[k] = v
        else
            # Non-string, non-numeric (e.g. Bool markers like "resolved" => true)
            # Skip — these are handled by parse_model, not the resolver
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

end # module SharedUtils


##############################################################################
# IMPLEMENTATION 1: Hand-rolled AST interpreter
##############################################################################
module HandRolled

using ..SharedUtils

# ── Upgraded AST: Using Indices instead of Names ──────────────────────────────
abstract type Expr_ end

struct Lit      <: Expr_; val::Float64; end
struct RefIdx   <: Expr_; idx::Int;     end  # Store index, not name
struct BinOp    <: Expr_; op::Symbol; lhs::Expr_; rhs::Expr_; end
struct UnaryOp  <: Expr_; op::Symbol; arg::Expr_; end
struct FuncCall <: Expr_; fn::Symbol; arg::Expr_; end

# ── Parser logic (Internal) ───────────────────────────────────────────────────
mangle(name::AbstractString) = "__" * replace(name, "," => "__")
demangle(id::AbstractString) = replace(lstrip(id, '_'), "__" => ",")

function preprocess(s::String)::String
    replace(s, r"\$([A-Za-z_][A-Za-z0-9_,]*)" => s -> mangle(s[2:end]))
end

function parse_to_idx(jexpr, name_to_idx::Dict{String, Int})::Expr_
    if jexpr isa Number
        return Lit(Float64(jexpr))
    elseif jexpr isa Symbol
        s = string(jexpr)
        if startswith(s, "__")
            name = demangle(s)
            return RefIdx(name_to_idx[name])
        elseif s == "pi"; return Lit(π)
        elseif s == "e";  return Lit(ℯ)
        else; error("Unknown symbol: $s")
        end
    elseif jexpr isa Expr
        h, args = jexpr.head, jexpr.args
        if h == :call
            fn = args[1]
            if length(args) == 3
                return BinOp(fn, parse_to_idx(args[2], name_to_idx), parse_to_idx(args[3], name_to_idx))
            elseif length(args) == 2
                if fn in (:+, :-); return UnaryOp(fn, parse_to_idx(args[2], name_to_idx))
                else;              return FuncCall(fn, parse_to_idx(args[2], name_to_idx))
                end
            end
        end
    end
    error("Cannot parse: $jexpr")
end

# ── Optimized Evaluator: Working on a Vector ──────────────────────────────────
eval_v(node::Lit,      v::Vector{Float64}) = node.val
eval_v(node::RefIdx,   v::Vector{Float64}) = v[node.idx] # Direct index access
eval_v(node::UnaryOp,  v::Vector{Float64}) = node.op == :- ? -eval_v(node.arg, v) : eval_v(node.arg, v)

function eval_v(node::BinOp, v::Vector{Float64})
    l, r = eval_v(node.lhs, v), eval_v(node.rhs, v)
    node.op == :+ ? l + r : node.op == :- ? l - r :
    node.op == :* ? l * r : node.op == :/ ? l / r :
    node.op == :^ ? l ^ r : error("Op error")
end

function eval_v(node::FuncCall, v::Vector{Float64})
    val = eval_v(node.arg, v)
    fn = node.fn
    fn == :sqrt ? sqrt(val) : fn == :sin ? sin(val) :
    fn == :cos  ? cos(val)  : fn == :exp ? exp(val) :
    fn == :log  ? log(val)  : abs(val) # Add more as needed
end

# ── Resolver ──────────────────────────────────────────────────────────────────
struct Resolver
    all_names   ::Vector{String}
    fixed_idxs  ::Vector{Int}
    fixed_vals  ::Vector{Float64}
    free_idxs   ::Vector{Int}
    derived_idxs::Vector{Int}
    compiled    ::Vector{Expr_}
end

function build_resolver(param_dict::Dict{String}, fit_params::Vector{String})
    free, fixed, derived = partition_dict(param_dict, fit_params)
    sorted = topo_sort(derived)

    # Create a stable global indexing for this model
    all_names = [fit_params; collect(keys(fixed)); sorted]
    name_to_idx = Dict(name => i for (i, name) in enumerate(all_names))

    fixed_names = collect(keys(fixed))
    f_idxs = [name_to_idx[n] for n in fixed_names]
    f_vals = [fixed[n] for n in fixed_names]
    free_idxs = [name_to_idx[n] for n in fit_params]
    d_idxs = [name_to_idx[n] for n in sorted]

    instr = [parse_to_idx(Meta.parse(preprocess(derived[n])), name_to_idx) for n in sorted]

    Resolver(all_names, f_idxs, f_vals, free_idxs, d_idxs, instr)
end

function (r::Resolver)(x::Vector{Float64})::Dict{String,Float64}
    # Pre-allocate the results vector
    vals = Vector{Float64}(undef, length(r.all_names))

    # 1. Load Fixed
    for i in 1:length(r.fixed_idxs); vals[r.fixed_idxs[i]] = r.fixed_vals[i]; end

    # 2. Load Free
    for i in 1:length(r.free_idxs);  vals[r.free_idxs[i]] = x[i]; end

    # 3. Compute Derived (no dict lookups)
    for i in 1:length(r.derived_idxs)
        vals[r.derived_idxs[i]] = eval_v(r.compiled[i], vals)
    end

    return Dict(zip(r.all_names, vals))
end

end


##############################################################################
# IMPLEMENTATION 2: RuntimeGeneratedFunctions.jl
##############################################################################
module RGF

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

using ..SharedUtils

mangle(name::AbstractString) = replace(name, "," => "__")

function preprocess(s::String)::String
    replace(s, r"\$([A-Za-z_][A-Za-z0-9_,]*)" => s -> mangle(s[2:end]))
end

struct Resolver
    all_names ::Vector{String}
    fn        ::Function
end

function build_resolver(param_dict::Dict{String}, fit_params::Vector{String})
    free, fixed, derived = partition_dict(param_dict, fit_params)
    sorted = topo_sort(derived)

    # 1. Define the names that will be returned in the output vector
    all_names = [fit_params; collect(keys(fixed)); sorted]

    # 2. Check whether any derived expression references WL or MJD.
    #    If so, the resolver is "chromatic" and all derived expressions
    #    are compiled with dot-broadcasting so that vector WL/MJD propagate.
    uses_wl_mjd = any(v -> any(r -> r in IMPLICIT_VARS && r ∉ ("R", "MU"),
                               extract_refs(v)),
                      values(derived))

    # 3. Build the function body as a single block
    # We embed fixed values as literals (constants) for the compiler.
    derived_stmts = if uses_wl_mjd
        # Broadcast all derived expressions so vector WL/MJD propagate naturally
        [:($(Symbol(mangle(name))) = $(broadcastify(Meta.parse(preprocess(derived[name])))))
         for name in sorted]
    else
        # Pure scalar path — no broadcasting overhead
        [:($(Symbol(mangle(name))) = $(Meta.parse(preprocess(derived[name]))))
         for name in sorted]
    end

    body = quote
        # Extract free parameters from the input vector x
        $([:($(Symbol(mangle(fit_params[i]))) = x[$i]) for i in eachindex(fit_params)]...)

        # Embed fixed parameters as local constants
        $([:($(Symbol(mangle(name))) = $val) for (name, val) in fixed]...)

        # Compute derived parameters in topological order
        $(derived_stmts...)

        # Return all values (may be a mix of scalars and vectors when chromatic)
        return $(Expr(:vect, [Symbol(mangle(n)) for n in all_names]...))
    end

    # 4. Compile: (x, WL, MJD) -> Vector
    fn_expr = :((x, WL, MJD) -> $body)
    compiled_fn = @RuntimeGeneratedFunction(fn_expr)

    # Wrap with default NaN for backward compatibility: fn(x) still works
    wrapped_fn(x, WL=NaN, MJD=NaN) = compiled_fn(x, WL, MJD)

    Resolver(all_names, wrapped_fn)
end

function (r::Resolver)(x::Vector{Float64}, wl=NaN, mjd=NaN)
    vals_vec = r.fn(x, wl, mjd)
    return Dict{String,Any}(zip(r.all_names, vals_vec))
end

end # module RGF


##############################################################################
# IMPLEMENTATION 3: Symbolics.jl
##############################################################################
module Symbolic

using ..SharedUtils

# Lazy-load Symbolics
const _symbolics_loaded = Ref(false)

function _load_symbolics()
    _symbolics_loaded[] && return true
    try
        @eval using Symbolics
        _symbolics_loaded[] = true
        return true
    catch e
        @warn "Symbolics.jl not available: $e"
        return false
    end
end

mangle(name::AbstractString) = replace(name, "," => "__")

function preprocess(s::String)::String
    replace(s, r"\$([A-Za-z_][A-Za-z0-9_,]*)" => s -> mangle(s[2:end]))
end

struct Resolver
    all_names ::Vector{String}
    fn        ::Function
end

function build_resolver(param_dict::Dict{String}, fit_params::Vector{String})
    _load_symbolics() || error("Symbolics.jl required for Symbolic resolver")

    free, fixed, derived = partition_dict(param_dict, fit_params)
    sorted = topo_sort(derived)

    # 1. Create symbolic variables for the free parameters only
    sym_fit = [Symbolics.variable(Symbol(mangle(n))) for n in fit_params]

    # 2. Build environment: name -> symbolic var OR numeric constant
    # By putting 'fixed' values here as numbers, they are folded into the expressions.
    env = Dict{String, Any}()
    for (i, name) in enumerate(fit_params)
        env[name] = sym_fit[i]
    end
    for (name, val) in fixed
        env[name] = val
    end

    # 3. Build symbolic expressions for derived parameters
    for name in sorted
        expr_str = preprocess(derived[name])
        # We build the graph by evaluating the Julia expression
        # using the symbols/constants currently in 'env'
        sym_expr = @eval begin
            $([:($( Symbol(mangle(k)) ) = $v) for (k,v) in env]...)
            $(Meta.parse(expr_str))
        end
        env[name] = sym_expr
    end

    # 4. Define the output vector (everything the resolver returns)
    all_names = [fit_params; collect(keys(fixed)); sorted]
    all_exprs = [env[n] for n in all_names]

    # 5. Compile to a native function: f(x) -> Vector{Float64}
    # build_function returns a tuple (f_out_of_place, f_in_place); we take the first.
    f_compiled = Symbolics.build_function(all_exprs, sym_fit, expression=Val{false})[1]

    Resolver(all_names, f_compiled)
end

function (r::Resolver)(x::Vector{Float64})::Dict{String,Float64}
    # Single native function call replaces the substitution loop
    vals_vec = r.fn(x)
    return Dict(zip(r.all_names, vals_vec))
end

end # module Symbolic
