# The Model perspective's data layer: a model dict, described.
#
# Two views are built from here, and they answer different questions.
#
# `model_rows` answers "what will the optimiser see". A parameter is in exactly one of three
# states, and they are not independent: anything named in `list_free_params` must be numeric,
# because the resolver raises on a string there. Free and derived are therefore mutually
# exclusive, and the table has to make that unrepresentable rather than report it afterwards.
#
# `model_inspection` answers "what did the parser actually understand". That is the view worth
# having: an unrecognised geometry key does not error, it silently changes the component's
# type — `"ring,gaussian_ring"` with a `fwhm` present yields a plain Gaussian and an inert
# parameter. Two demos shipped that way. The parser can already name the keys it ignored;
# this surfaces them.
#
# Nothing here draws. It is plain Julia over the OITOOLS parser, so it can be tested by
# comparing against the parser rather than against a screenshot.

"""How a parameter reaches the optimiser."""
@enum ParamMode PARAM_FIXED PARAM_FREE PARAM_EXPR

"""
One row of the parameter table.

`fitindex` is the position in `list_free_params`, and therefore the index into every `x` an
optimiser touches and every `x_opt` a result reports — `x[i] ↔ list_free_params[i]` holds
throughout OITOOLS. Showing it is what makes an exported script legible.

`value` for a `PARAM_EXPR` row is the resolved number when the expression can be evaluated
without data, and `NaN` when it cannot (it references `\$WL` or `\$MJD`, say). The distinction
matters: a blank cell and a computed cell mean different things.
"""
struct ParamRow
    component :: String
    param     :: String
    key       :: String
    mode      :: ParamMode
    value     :: Float64
    expr      :: String
    lb        :: Float64
    ub        :: Float64
    fitindex  :: Int
    atbound   :: Bool
end

"The component name the parser gives to bare keys, e.g. `\"PA\"` rather than `\"star,PA\"`."
const GLOBAL_COMPONENT = "__global__"

_split_key(k::AbstractString) = (i = findfirst(',', k);
                                 i === nothing ? (GLOBAL_COMPONENT, String(k)) :
                                                 (String(k[1:i-1]), String(k[i+1:end])))

"""
    model_rows(model_dict, list_free_params; lb, ub, wl, mjd) -> Vector{ParamRow}

Describe every parameter in the model: its state, value, bounds and fit-vector index.

Bounds come from `default_bounds` unless supplied. Rows are ordered globals first, then
components in the parser's own order, and within a component the geometry key leads — the key
that determines what the component *is* should not be buried among its modifiers.

`wl` and `mjd` are the wavelength and epoch a chromatic expression is DISPLAYED at. Such an
expression has no single value — the resolver broadcasts it to one entry per uv point — so
without them the row cannot be resolved at all and reports nothing. Pass the dataset's first
wavelength and the panel shows the model as it stands there.
"""
function model_rows(model_dict::AbstractDict, list_free_params::AbstractVector;
                    lb = nothing, ub = nothing, wl = nothing, mjd = nothing)
    md   = Dict{String,Any}(String(k) => v for (k, v) in model_dict)
    free = String.(list_free_params)
    # A chromatic expression has no single value, so the row shows it at one wavelength. The
    # substitution is a COPY: `_resolve_numeric` resolves `\$WL` by looking `"WL"` up in the
    # dict it is handed, so putting the number there is all it takes -- and putting it in the
    # real dict would make WL a global parameter of the model.
    disp = md
    if wl !== nothing || mjd !== nothing
        disp = copy(md)
        wl  === nothing || (disp["WL"]  = Float64(wl))
        mjd === nothing || (disp["MJD"] = Float64(mjd))
    end
    dlb, dub = try
        default_bounds(md, free)
    catch
        (Dict{String,Float64}(), Dict{String,Float64}())
    end
    lbd = lb === nothing ? dlb : Dict{String,Float64}(String(k) => Float64(v) for (k, v) in lb)
    ubd = ub === nothing ? dub : Dict{String,Float64}(String(k) => Float64(v) for (k, v) in ub)

    comps = try
        String.(OITOOLS._component_names(md))
    catch
        String[]
    end
    geomkey = Dict{String,String}()
    for c in comps
        k = _geometry_key(c, md)
        k === nothing || (geomkey[c] = k)
    end

    order = Dict(c => i for (i, c) in enumerate(comps))
    function rank(key)
        comp, par = _split_key(key)
        comp == GLOBAL_COMPONENT && return (0, 0, par)
        (get(order, comp, typemax(Int)) , get(geomkey, comp, "") == par ? 0 : 1, par)
    end

    rows = ParamRow[]
    for key in sort(collect(keys(md)); by = rank)
        comp, par = _split_key(key)
        v = md[key]
        idx = something(findfirst(==(key), free), 0)

        if idx > 0
            # Free parameters must be numeric; the resolver raises on a string. A model that
            # says otherwise is malformed, and the row should say so rather than pretend.
            val = v isa Number ? Float64(v) : NaN
            l = Float64(get(lbd, key, -Inf))
            u = Float64(get(ubd, key,  Inf))
            near = isfinite(val) && (_sits_on(val, l) || _sits_on(val, u))
            push!(rows, ParamRow(comp, par, key, PARAM_FREE, val, "", l, u, idx, near))
        elseif v isa AbstractString
            resolved = try
                Float64(OITOOLS._resolve_numeric(key, disp))
            catch
                NaN
            end
            push!(rows, ParamRow(comp, par, key, PARAM_EXPR, resolved, String(v),
                                 NaN, NaN, 0, false))
        else
            push!(rows, ParamRow(comp, par, key, PARAM_FIXED, Float64(v), "", NaN, NaN, 0, false))
        end
    end
    return rows
end

"""
Is this value sitting on that bound?

Measured against the bound, not against the span. `default_bounds` is deliberately generous —
a diameter gets 0 to 1000 mas — so a span-relative test would flag every small diameter as
pinned and the warning would mean nothing. What is worth reporting is the optimiser having
pushed a parameter to its limit, and that shows up as the value landing on the bound to within
solver tolerance.
"""
_sits_on(v::Float64, b::Float64) = isfinite(b) && abs(v - b) <= 1e-6 * max(1.0, abs(b))

"Which key made the parser decide this component's type, or `nothing`."
function _geometry_key(comp::AbstractString, md::AbstractDict)
    for gk in OITOOLS._GEOMETRY_KEYS
        haskey(md, "$comp,$gk") && return String(gk)
    end
    return haskey(md, "$comp,f") ? "f" : nothing
end

"""
What the parser made of one component.

`geometry_key` is the key that decided `kind`. `_GEOMETRY_KEYS` is searched in order and the
first match wins, so `profile` beats `diamin` beats `ud` — which is why naming the deciding key
is worth more than naming the kind alone.
"""
struct ComponentInfo
    name         :: String
    kind         :: Symbol
    geometry_key :: String
    params       :: Vector{String}
    unrecognised :: Vector{String}
end

"""
    model_inspection(model_dict) -> (; components, unrecognised, broadcasting, globals)

Report what the parser understood, including the keys it ignored.

`broadcasting` is true when any expression references `\$WL` or `\$MJD`. That is a behavioural
switch, not a detail: if *any* derived expression is chromatic the resolver broadcasts *all* of
them, so every parameter silently becomes a per-uv-point vector. The user should see it flip.
"""
function model_inspection(model_dict::AbstractDict)
    md = Dict{String,Any}(String(k) => v for (k, v) in model_dict)

    unknown = String[]
    try
        # The parser already knows which keys it ignored; it warns. Ask it directly, with the
        # warning suppressed — this is a panel, not a log.
        unknown = String.(Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
            OITOOLS._warn_unrecognised_keys(md)
        end)
    catch
    end

    comps = ComponentInfo[]
    names = try
        String.(OITOOLS._component_names(md))
    catch
        String[]
    end
    for c in names
        kind = try
            OITOOLS._identify_kind(c, md)
        catch
            :unknown
        end
        gk = something(_geometry_key(c, md), "")
        pre = c * ","
        ps  = sort([k[(length(pre)+1):end] for k in keys(md) if startswith(k, pre)])
        cu  = sort([k for k in unknown if startswith(k, pre)])
        push!(comps, ComponentInfo(c, kind, gk, ps, cu))
    end

    globals = sort([k for k in keys(md) if !occursin(',', k)])

    chromatic = any(v isa AbstractString && (occursin("\$WL", v) || occursin("\$MJD", v))
                    for v in values(md))

    return (; components = comps, unrecognised = sort(unknown),
              broadcasting = chromatic, globals)
end

"""
    free_parameter_vector(rows) -> Vector{ParamRow}

The free rows in fit-vector order, so row `i` is `x[i]`.
"""
free_parameter_vector(rows::AbstractVector{ParamRow}) =
    sort(filter(r -> r.mode === PARAM_FREE, rows); by = r -> r.fitindex)


# ─────────────────────────────────────────────────────────────────────────────
# Model → image
# ─────────────────────────────────────────────────────────────────────────────
#
# Rendering a model is how you find out what you actually typed. A parameter table says
# `"disk,fwhmin" => 1.5`; an image says whether that is the ring you meant. It is also the only
# way to see what a chromatic model does, because the parameter values do not change on the
# page — the picture does.

"""
    model_depends_on(model_dict) -> (; wl, mjd)

Whether any expression in the model references `\\\$WL` or `\\\$MJD`.

Separately, not as one `broadcasting` flag: the panel offers a wavelength control and a time
control, and each should be live only if it changes the answer. Offering both for a model that
uses neither invites the user to conclude the render is broken when nothing moves.
"""
function model_depends_on(model_dict::AbstractDict)
    vals = [v for v in values(model_dict) if v isa AbstractString]
    return (; wl  = any(v -> occursin("\$WL",  v), vals),
              mjd = any(v -> occursin("\$MJD", v), vals))
end

"""
    render_model_image(model_dict, free, x; nx, pixsize, wl, mjd) -> (; image, flux, …)

The model as a picture, at one wavelength and one epoch.

`x` is the free-parameter vector, in `free` order — the same `x` every optimiser sees, so what
is rendered is what would be fitted. Pass an empty `free` to render a model with nothing free,
which is the usual case here.

`wl` is in METRES, matching `\\\$WL`; `mjd` in days. Both may be `nothing`, and are ignored by a
model that does not reference them — which [`model_depends_on`](@ref) reports, so the panel can
say so rather than let the user wonder.
"""
function render_model_image(model_dict::AbstractDict, free::AbstractVector, x::AbstractVector;
                            nx::Integer = 128, pixsize::Real = 0.1,
                            wl = nothing, mjd = nothing, oversample::Integer = 1)
    md = Dict{String,Any}(String(k) => v for (k, v) in model_dict)
    fm = parse_model(md, String.(free); nB_workspace = 1)
    img = model_to_image(fm, Float64.(x); nx = Int(nx), pixsize = Float64(pixsize),
                         oversample = Int(oversample), wl = wl, mjd = mjd)
    dep = model_depends_on(md)
    return (; image = Float64.(img), flux = sum(img), nx = Int(nx),
              pixsize = Float64(pixsize), wl, mjd, depends = dep,
              fov = Int(nx) * Float64(pixsize))
end

"""
    model_image_cube(model_dict, free, x, wls; nx, pixsize, mjd) -> Array{Float64,3}

The model rendered across a wavelength grid, one plane per channel.

This is what `simulate` takes as a cube, and what a chromatic model is for: the planes differ
only if the model actually depends on `\\\$WL`, and identical planes are the answer to "is this
model chromatic" rather than a fault.
"""
function model_image_cube(model_dict::AbstractDict, free::AbstractVector, x::AbstractVector,
                          wls; nx::Integer = 128, pixsize::Real = 0.1, mjd = nothing)
    ws = collect(wls)
    cube = Array{Float64,3}(undef, Int(nx), Int(nx), length(ws))
    for (k, w) in enumerate(ws)
        cube[:, :, k] = render_model_image(model_dict, free, x;
                                           nx, pixsize, wl = w, mjd).image
    end
    return cube
end


# ─────────────────────────────────────────────────────────────────────────────
# Running a fit from the panel
# ─────────────────────────────────────────────────────────────────────────────
#
# The panel holds strings; the fitters take a dict, a free list, two bound dicts, a constraint
# list and a prior list. Parsing happens here, in one place, so the panel never has to know
# which fitter wants which shape — and so a malformed row is reported against the row rather
# than as a failure somewhere inside NLopt.

"`key=value` lines to a model dict. A value that parses as a number is one; anything else is an expression."
function parse_model_lines(text::AbstractString)
    md = Dict{String,Any}()
    for line in split(String(text), "\n")
        isempty(strip(line)) && continue
        i = findfirst('=', line)
        i === nothing && continue
        k = String(strip(line[1:i-1]))
        v = String(strip(line[nextind(line, i):end]))
        isempty(k) && continue
        md[k] = something(tryparse(Float64, v), v)
    end
    return md
end

"""
    parse_free_lines(text) -> (free, lb, ub)

`key\\tlb\\tub` lines, in fit-vector order.

The order is the contract: `x[i] ↔ free[i]` throughout OITOOLS, so the rows arrive in the order
the panel shows them and are not sorted here. An empty bound is ±Inf, which every fitter treats
as unbounded — except nested sampling, which says so.
"""
function parse_free_lines(text::AbstractString)
    free = String[]
    lb = Dict{String,Float64}(); ub = Dict{String,Float64}()
    for line in split(String(text), "\n")
        isempty(strip(line)) && continue
        f = split(line, "\t")
        k = String(strip(f[1])); isempty(k) && continue
        push!(free, k)
        length(f) >= 2 && (lb[k] = something(tryparse(Float64, strip(f[2])), -Inf))
        length(f) >= 3 && (ub[k] = something(tryparse(Float64, strip(f[3])),  Inf))
    end
    return free, lb, ub
end

"`lhs\\top\\trhs\\ttol` lines to `ModelConstraint`s."
function parse_constraint_lines(text::AbstractString)
    out = ModelConstraint[]
    for line in split(String(text), "\n")
        isempty(strip(line)) && continue
        f = split(line, "\t")
        length(f) >= 3 || continue
        rhs_s = String(strip(f[3]))
        rhs = something(tryparse(Float64, rhs_s), rhs_s)
        tol = length(f) >= 4 ? something(tryparse(Float64, strip(f[4])), DEFAULT_CONSTRAINT_TOL) :
                               DEFAULT_CONSTRAINT_TOL
        push!(out, ModelConstraint(String(strip(f[1])), String(strip(f[2])), rhs; tol))
    end
    return out
end

"`expr\\ttarget\\tsigma` lines to the prior tuples `fit_model` takes."
function parse_prior_lines(text::AbstractString)
    out = Tuple{String,Float64,Float64}[]
    for line in split(String(text), "\n")
        isempty(strip(line)) && continue
        f = split(line, "\t")
        length(f) >= 3 || continue
        t = tryparse(Float64, strip(f[2])); s = tryparse(Float64, strip(f[3]))
        (t === nothing || s === nothing || s <= 0) && continue
        push!(out, (String(strip(f[1])), t, s))
    end
    return out
end

"""
    fitting_weights(; v2, t3amp, t3phi, cvis, flux, diffvis) -> Vector{Float64}

The SEVEN-element weight vector model fitting takes: `[V², T3amp, T3φ, visamp, visφ, flux, diffφ]`.

Seven, not imaging's three, and `cvis` is one box driving TWO entries — amplitude and phase are
one observable to a user and two to the criterion. Building it here is what stops the panel
having to know that.
"""
fitting_weights(; v2::Bool = true, t3amp::Bool = true, t3phi::Bool = true,
                  cvis::Bool = false, flux::Bool = false, diffvis::Bool = false) =
    Float64[v2 ? 1 : 0, t3amp ? 1 : 0, t3phi ? 1 : 0,
            cvis ? 1 : 0, cvis ? 1 : 0, flux ? 1 : 0, diffvis ? 1 : 0]
