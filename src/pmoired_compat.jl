"""
    pmoired_to_julia(s) -> String

Convert a PMOIRED Python dict literal string to an equivalent Julia Dict
literal string.  Handles nested dicts, expression strings with `\$` references,
single-quoted strings, Python boolean/None literals, and list literals.

Transformation rules applied:
  { ... }           ->  Dict( ... )
  'key': value      ->  "key" => value         (outside string literals)
  :                 ->  =>                      (outside string literals, single colon only)
  'plain string'    ->  "plain string"
  `'expr with \$ref'` ->  `raw"expr with \$ref"`   (raw string prevents interpolation)
  True / False      ->  true / false
  None              ->  nothing

Example
-------
```python
# Python / PMOIRED:
param = {
    'star,ud':    3.2,
    'star,f':     0.7,
    'ring,udout': '\$star,ud * 8',
    'ring,f':     '1 - \$star,f',
    'ring,incl':  30.0,
}
fitOnly = ['star,ud', 'star,f', 'ring,udout', 'ring,f']
```

```julia
# Generated Julia:
param = Dict(
    "star,ud"    => 3.2,
    "star,f"     => 0.7,
    "ring,udout" => raw"\$star,ud * 8",
    "ring,f"     => raw"1 - \$star,f",
    "ring,incl"  => 30.0,
)
fitOnly = ["star,ud", "star,f", "ring,udout", "ring,f"]
```
"""
function pmoired_to_julia(s::AbstractString)::String
    buf   = IOBuffer()
    i     = firstindex(s)
    last  = lastindex(s)

    # ------------------------------------------------------------------ helpers
    # Advance iterator safely
    nextidx(s, i) = i > last ? i : nextind(s, i)

    # Read a quoted string (single or double), return (content, next_i)
    function read_string(s, start, delim)
        content = IOBuffer()
        i = nextidx(s, start)          # skip opening quote
        while i <= last
            c = s[i]
            if c == '\\'               # escape sequence — preserve verbatim
                write(content, c)
                i = nextidx(s, i)
                if i <= last
                    write(content, s[i])
                    i = nextidx(s, i)
                end
            elseif c == delim          # closing quote
                i = nextidx(s, i)
                break
            else
                write(content, c)
                i = nextidx(s, i)
            end
        end
        return String(take!(content)), i
    end

    # ------------------------------------------------------------------ main loop
    while i <= last
        c = s[i]

        # ---- single-quoted string ----------------------------------------
        if c == '\''
            content, i = read_string(s, i, '\'')
            if occursin('$', content)
                write(buf, "raw\"", content, "\"")
            else
                write(buf, '"', content, '"')
            end

        # ---- double-quoted string ----------------------------------------
        elseif c == '"'
            content, i = read_string(s, i, '"')
            if occursin('$', content)
                write(buf, "raw\"", content, "\"")
            else
                write(buf, '"', content, '"')
            end

        # ---- Python hash comment — preserve as Julia comment -------------
        elseif c == '#'
            while i <= last && s[i] != '\n'
                write(buf, s[i])
                i = nextidx(s, i)
            end

        # ---- { -> Dict( --------------------------------------------------
        elseif c == '{'
            write(buf, "Dict(")
            i = nextidx(s, i)

        # ---- } -> ) ------------------------------------------------------
        elseif c == '}'
            write(buf, ')')
            i = nextidx(s, i)

        # ---- : -> =>  (key-value separator; avoid :: ) -------------------
        elseif c == ':'
            next_i = nextidx(s, i)
            if next_i <= last && s[next_i] == ':'   # :: — pass through
                write(buf, "::")
                i = nextidx(s, next_i)
            else
                write(buf, " =>")
                i = next_i
            end

        # ---- everything else passes through verbatim --------------------
        else
            write(buf, c)
            i = nextidx(s, i)
        end
    end

    result = String(take!(buf))

    # ---- Python keywords → Julia keywords --------------------------------
    # Use word-boundary anchors to avoid substring replacement
    result = replace(result, r"\bTrue\b"  => "true")
    result = replace(result, r"\bFalse\b" => "false")
    result = replace(result, r"\bNone\b"  => "nothing")

    return result
end


"""
    pmoired_to_dict(s) -> Dict{String,Any}

Convert a PMOIRED Python dict literal string directly to a Julia `Dict`.
This is a convenience wrapper around `pmoired_to_julia` that evaluates the
resulting Julia string, so users don't need `eval(Meta.parse(...))`.

```julia
model_dict = pmoired_to_dict("{'star,ud': 3.2, 'ring,f': '1 - \\\$star,f'}")
```
"""
function pmoired_to_dict(s::AbstractString)::Dict{String,Any}
    julia_str = pmoired_to_julia(s)
    return _adopt_pmoired_keys(eval(Meta.parse(julia_str)))
end

"""
    _adopt_pmoired_keys(d) -> Dict

Rename the PMOIRED keys OITOOLS spells differently, and warn about the one that cannot be
carried across faithfully.

This is a transpiler of syntax, so keys otherwise arrive verbatim -- which is right for almost
all of them, since the two packages share their vocabulary. The exception is the spatial
kernel. PMOIRED writes it `"comp,spatial kernel"`, with a space, and applies it PER COMPONENT;
OITOOLS spells it `spatial_kernel` and applies one kernel to the whole model. A single
component converts exactly. More than one cannot: adopting one component's kernel would smooth
the others too, and adopting none would drop a part of the model. The keys are left in place in
that case, so the inspector reports them as unrecognised rather than the import quietly
changing what the model means.
"""
function _adopt_pmoired_keys(d::AbstractDict)
    out = Dict{String,Any}(String(k) => v for (k, v) in d)
    kern = [k for k in keys(out) if endswith(k, ",spatial kernel")]
    if length(kern) == 1
        out["spatial_kernel"] = out[kern[1]]
        delete!(out, kern[1])
        @info "PMOIRED import: '$(kern[1])' adopted as the model-wide `spatial_kernel`. " *
              "OITOOLS applies one kernel to the whole model, which is the same thing here " *
              "because there is one component carrying one."
    elseif length(kern) > 1
        @warn "PMOIRED import: $(length(kern)) per-component spatial kernels, which OITOOLS " *
              "cannot represent -- it applies a single kernel to the whole model. The keys " *
              "are left as they are and will be reported as unrecognised; smooth the " *
              "components with a `profile` each, or keep one kernel." components = kern
    end
    # The Lorentzian fraction only means anything beside a kernel, so it travels with it.
    for k in [k for k in keys(out) if endswith(k, ",spatial kernel frac lorentzian")]
        @warn "PMOIRED import: '$k' has no OITOOLS equivalent -- its kernel is Gaussian only."
    end
    return out
end


"""
    pmoired_to_julia_file(infile, outfile)

Read a Python/PMOIRED notebook snippet from `infile`, convert it, write to `outfile`.
Lines that look like `param = {...}` or `fitOnly = [...]` are converted;
other lines are passed through.
"""
function pmoired_to_julia_file(infile::AbstractString, outfile::AbstractString)
    open(outfile, "w") do out
        for line in eachline(infile)
            println(out, pmoired_to_julia(line))
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Julia -> PMOIRED
# ─────────────────────────────────────────────────────────────────────────────

# Geometry keys OITOOLS understands that PMOIRED does not. A model using these is still
# written out, but it will not mean the same thing on the other side, so exporting warns.
const _PMOIRED_UNSUPPORTED_SUFFIXES = ["ldlin", "ldquad", "ldpow", "resolved"]

_pmoired_value(v::Bool)            = v ? "True" : "False"
_pmoired_value(::Nothing)          = "None"
_pmoired_value(v::Integer)         = string(v)
_pmoired_value(v::Real)            = string(float(v))
_pmoired_value(v::AbstractString)  = "'" * replace(v, "'" => "\\'") * "'"
_pmoired_value(v)                  = string(v)

"""
    dict_to_pmoired(model_dict; check=true) -> String

Render an OITOOLS model dictionary as PMOIRED source — the inverse of
[`pmoired_to_dict`](@ref).

Keys are emitted in sorted order so the output is reproducible and diffable. Julia `true`,
`false` and `nothing` become `True`, `False` and `None`; string values are quoted with
single quotes, so `\$`-expressions survive unchanged (both packages use the same syntax).

With `check=true` (the default) a warning is emitted for keys that have no PMOIRED
equivalent — the `ldlin`/`ldquad`/`ldpow`/`resolved` geometries are OITOOLS additions. The
model is still written; the warning exists because the failure is otherwise silent on the
PMOIRED side.

See also [`pmoired_to_dict`](@ref), [`dict_to_pmoired_file`](@ref).
"""
function dict_to_pmoired(model_dict::AbstractDict; check::Bool = true)
    if check
        bad = String[]
        for k in keys(model_dict)
            i = findfirst(==(','), k)
            suffix = i === nothing ? k : k[nextind(k, i):end]
            (suffix in _PMOIRED_UNSUPPORTED_SUFFIXES) && push!(bad, k)
        end
        isempty(bad) || @warn "These keys have no PMOIRED equivalent and will not be " *
                              "interpreted the same way after import" keys=sort(bad)
    end
    ks = sort(collect(keys(model_dict)))
    body = join(("    '$(k)': $(_pmoired_value(model_dict[k]))" for k in ks), ",\n")
    return "{\n" * body * ",\n}"
end

"""
    dict_to_pmoired_file(model_dict, outfile; check=true)

Write `model_dict` to `outfile` as PMOIRED source. Returns `outfile`.
Inverse of [`pmoired_to_julia_file`](@ref).
"""
function dict_to_pmoired_file(model_dict::AbstractDict, outfile::AbstractString; check::Bool = true)
    open(outfile, "w") do io
        println(io, dict_to_pmoired(model_dict; check = check))
    end
    return outfile
end
