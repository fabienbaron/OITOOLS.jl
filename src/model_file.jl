# model_file.jl
#
# TOML model files: everything needed to reproduce a fit, in one human-editable text file.
#
# Load order (within OITOOLS module):
#   include("parse_model.jl")        # default_bounds
#   include("constraints.jl")        # ModelConstraint
#   include("model_file.jl")         # this file
#
# What a model file holds, and why each part is in it
# ---------------------------------------------------
# A model dict alone does not describe a fit. `fit_model(model, free, data; lb, ub, ...)` takes
# four more arguments that change the answer, and none of them lived anywhere but in the
# script that happened to call it:
#
#     [model]        the parameter dict, values or expressions — the only part that was ever
#                    portable, via the PMOIRED transpiler
#     free           which parameters the optimiser moves; `x[i] ↔ free[i]` throughout OITOOLS
#     [bounds]       the box each free parameter lives in
#     [[constraints]] relations bounds cannot express — see constraints.jl
#     [[priors]]     Gaussian pulls toward a value
#
# So a file that carried only the dict would round-trip a model and lose the fit. All five are
# written, and `read_model_file` returns them in the shape the fitters take, so the round trip
# is a fit, not just a model.
#
# Format choice
# -------------
# TOML, because the package already depends on it for every facility, combiner and target
# config, because it diffs and merges as text, and because a model file is something a user
# edits by hand between runs. Keys are quoted since `"star,ud"` contains a comma; values are
# TOML numbers, or strings when they are expressions.
#
#     free = ["star,ud", "disk,pa"]
#
#     [model]
#     "star,ud"   = 6.5
#     "disk,fwhm" = "$star,ud * 3"
#
#     [bounds]
#     "star,ud" = [0.0, 20.0]
#
#     [[constraints]]
#     param = "disk,diamout"
#     op    = ">"
#     value = "disk,diamin"
#     tol   = 0.001
#
#     [[priors]]
#     expr   = "star,ud"
#     target = 6.0
#     sigma  = 0.5
#
# `free` is a top-level key rather than a member of `[model]` so that `[model]` mirrors the
# model dict exactly: a model may legitimately hold a bare global key, and a global named
# `free` would otherwise be silently eaten by the free-parameter list.

using TOML

"""
    read_model_file(path) -> (; model, free, lb, ub, constraints, priors, name)

Read a TOML model file into the arguments a fitter takes.

The result destructures straight into any of the fitters:

```julia
m = read_model_file("binary.toml")
res = fit_model(m.model, m.free, data; m.lb, m.ub, m.constraints, m.priors)
```

Every section is optional. A file with only `[model]` reads back as a model with no free
parameters, which is what `model_to_image` and `model_to_sed` want. Bounds absent from the file
are **not** filled in from [`default_bounds`](@ref) — an absent bound means ±Inf, and silently
substituting a guess would make two fits from the same file differ by which version of the
defaults was in the package.

See also [`write_model_file`](@ref), [`ModelConstraint`](@ref) and [`default_bounds`](@ref).
"""
function read_model_file(path::AbstractString)
    isfile(path) || error("No model file at '$path'")
    doc = TOML.parsefile(path)

    haskey(doc, "model") || error("'$path' has no [model] section, so there is no model to read")
    model = Dict{String,Any}(String(k) => v for (k, v) in doc["model"])

    free = String[String(p) for p in get(doc, "free", String[])]
    for p in free
        haskey(model, p) ||
            error("'$path' lists '$p' as free, but [model] has no such parameter")
        model[p] isa Number ||
            error("'$path' lists '$p' as free, but its value is the expression " *
                  "$(repr(model[p])); a free parameter must be numeric")
    end

    lb = Dict{String,Float64}()
    ub = Dict{String,Float64}()
    for (k, v) in get(doc, "bounds", Dict{String,Any}())
        (v isa AbstractVector && length(v) == 2) ||
            error("Bounds for '$k' in '$path' must be a two-element array [lower, upper], got $(repr(v))")
        lb[String(k)], ub[String(k)] = Float64(v[1]), Float64(v[2])
        lb[String(k)] < ub[String(k)] ||
            error("Bounds for '$k' in '$path' are empty: lower $(v[1]) is not below upper $(v[2])")
    end

    constraints = parse_constraints(get(doc, "constraints", []))

    priors = Tuple{String,Float64,Float64}[]
    for p in get(doc, "priors", [])
        for k in ("expr", "target", "sigma")
            haskey(p, k) || error("A [[priors]] entry in '$path' has no '$k'")
        end
        Float64(p["sigma"]) > 0 ||
            error("Prior on '$(p["expr"])' in '$path' has sigma $(p["sigma"]); it must be positive")
        push!(priors, (String(p["expr"]), Float64(p["target"]), Float64(p["sigma"])))
    end

    return (; model, free, lb, ub, constraints, priors,
              name = String(get(doc, "name", "")))
end

"""
    write_model_file(path, model_dict; free, lb, ub, constraints, priors, name)

Write a model and its fit settings to a TOML file readable by [`read_model_file`](@ref).

```julia
write_model_file("binary.toml", model; free, lb, ub,
                 constraints = [ModelConstraint("disk,diamout", ">", "disk,diamin")])
```

Keys are sorted, so the same model always produces the same bytes and two files can be
diffed. Bounds that are infinite on both sides are omitted rather than written as `[-inf, +inf]`:
they constrain nothing, and every fitter already treats a missing bound that way.

A `FitResult` is not what gets written — the fitted values are, if the caller passes them in
`model_dict`. Writing the result of a fit back out is therefore
`write_model_file(path, merge(model, Dict(zip(res.list_free_params, res.x_opt))); free = res.list_free_params, ...)`.
"""
function write_model_file(path::AbstractString, model_dict::AbstractDict;
                          free        = String[],
                          lb          = Dict{String,Float64}(),
                          ub          = Dict{String,Float64}(),
                          constraints = ModelConstraint[],
                          priors      = [],
                          name        ::AbstractString = "")
    freev = String[String(p) for p in free]
    for p in freev
        haskey(model_dict, p) ||
            error("'$p' is listed as free but is not in the model")
        model_dict[p] isa Number ||
            error("'$p' is listed as free but its value is the expression $(repr(model_dict[p]))")
    end
    cs = parse_constraints(constraints)
    # Bounds arrive as either of the two shapes the fitters take. A vector is positional, so it
    # only has meaning against `free`; normalise both to dicts before anything reads them.
    lbd = _bounds_dict(lb, freev, -Inf)
    ubd = _bounds_dict(ub, freev,  Inf)

    open(path, "w") do io
        println(io, "# OITOOLS model file")
        isempty(name) || println(io, "name = ", _toml_str(name))
        # `free` must precede every table header, since a bare TOML key belongs to whichever
        # table it follows.
        println(io, "free = [", join((_toml_str(p) for p in freev), ", "), "]")

        println(io, "\n[model]")
        for k in sort(collect(keys(model_dict)))
            v = model_dict[k]
            println(io, _toml_key(k), " = ", v isa AbstractString ? _toml_str(v) : _toml_num(v))
        end

        # Both key sets, since a parameter bounded from one side only is still bounded.
        bkeys = sort([k for k in union(keys(lbd), keys(ubd))
                      if isfinite(get(lbd, k, -Inf)) || isfinite(get(ubd, k, Inf))])
        if !isempty(bkeys)
            println(io, "\n[bounds]")
            for k in bkeys
                println(io, _toml_key(k), " = [",
                        _toml_num(get(lbd, k, -Inf)), ", ", _toml_num(get(ubd, k, Inf)), "]")
            end
        end

        for c in cs
            println(io, "\n[[constraints]]")
            println(io, "param = ", _toml_str(c.lhs))
            println(io, "op    = ", _toml_str(_op_string(c.op)))
            println(io, "value = ", c.rhs isa Float64 ? _toml_num(c.rhs) : _toml_str(c.rhs))
            println(io, "tol   = ", _toml_num(c.tol))
        end

        for (expr, target, sigma) in priors
            println(io, "\n[[priors]]")
            println(io, "expr   = ", _toml_str(String(expr)))
            println(io, "target = ", _toml_num(Float64(target)))
            println(io, "sigma  = ", _toml_num(Float64(sigma)))
        end
    end
    return path
end

# `lb`/`ub` in the fitters may be a Dict keyed by name or a Vector aligned with the free list.
# A file is keyed by name, so a vector is read against `free` here rather than being rejected.
function _bounds_dict(b, free::Vector{String}, default::Float64)
    b isa AbstractDict && return Dict{String,Float64}(String(k) => Float64(v) for (k, v) in b)
    length(b) == length(free) ||
        error("A bounds vector must have one entry per free parameter: got $(length(b)) for $(length(free))")
    return Dict{String,Float64}(free[i] => Float64(b[i]) for i in eachindex(free))
end

# TOML bare keys cannot contain a comma, and every component parameter does. Quoting only what
# needs it keeps a file of bare globals readable.
_toml_key(k::AbstractString) =
    occursin(r"^[A-Za-z0-9_-]+$", k) ? String(k) : _toml_str(k)

# TOML basic strings, not Julia ones. Julia's `repr` escapes `$` — every model expression starts
# with one — and TOML has no `\$` escape, so a file written with `repr` does not parse back.
function _toml_str(v::AbstractString)
    io = IOBuffer()
    print(io, '"')
    for ch in v
        ch == '"'       ? print(io, "\\\"") :
        ch == '\\'      ? print(io, "\\\\") :
        ch == '\n'      ? print(io, "\\n")  :
        ch == '\t'      ? print(io, "\\t")  :
        ch == '\r'      ? print(io, "\\r")  :
        ch == '\b'      ? print(io, "\\b")  :
        ch == '\f'      ? print(io, "\\f")  :
        iscntrl(ch)     ? print(io, "\\u", string(UInt16(ch); base = 16, pad = 4)) :
                          print(io, ch)
    end
    print(io, '"')
    return String(take!(io))
end

# TOML spells the non-finite floats in lower case with a mandatory sign on infinity, so Julia's
# own `Inf`/`NaN` are not valid TOML and would not read back.
_toml_num(v::AbstractFloat) =
    isnan(v) ? "nan" : isinf(v) ? (v > 0 ? "+inf" : "-inf") : repr(Float64(v))
_toml_num(v::Integer) = string(v)
_toml_num(v) = repr(v)
