# The Observe perspective's data layer: what is observable, from where, and when.
#
# This perspective is the odd one out. Exploring, Modeling and Imaging all CONSUME a dataset;
# Observing PRODUCES one. Everything here is upstream of `simulate`, and the ordering of the
# panels follows from that: an instrument has to be chosen before anything can be computed
# about it, and a target before anything can be computed about the night.
#
# Nothing here draws. As with `model.jl` and `imaging.jl`, that is what lets it be tested
# against the library rather than against a screenshot — and this perspective needs it more
# than the others, because a planning tool that is confidently wrong about when a target is up
# is worse than no planning tool.

"""
    config_catalog() -> (; facilities, combiners, wavelengths, targets, wave_combiner, unknown)

The configuration files on disk, classified.

A thin pass over [`list_configs`](@ref) that keeps `wave_combiner`, which is the only
machine-readable link between a spectral setup and the combiner it belongs to — the `combiner`
key inside each wavelength TOML. Without it the setup list cannot be filtered by combiner, and
a user is offered spectral setups the chosen instrument cannot produce.
"""
config_catalog() = list_configs()

"""
    facility_telescopes(name) -> Vector{String}

Telescope names for a facility, in the order its configuration lists them.

The order is not cosmetic: `compute_delays` and `in_delay` take a `config::Vector{Int}` indexing
into exactly this list, so a panel that reordered it would silently plan a different array.
"""
function facility_telescopes(name::AbstractString)
    f = read_facility_file(String(name))
    return String.(f.tel_names)
end

"""
    telescope_config(selected, all) -> Vector{Int}

The `config` vector `compute_delays` wants, from a set of selected telescope names.

One entry per telescope in the facility's own order, 1 for in and 0 for out.
"""
telescope_config(selected, all) = Int[(String(t) in String.(selected)) ? 1 : 0 for t in all]

"""
Height of a target's own bar, in the chart's 0–10 row units.

Thinner than the constraint rows are tall. At 2.0 the bar was a fifth of the whole chart and
crowded the az/alt annotations sitting just above and below it; 1.2 leaves them clear while
keeping it the most prominent row, which it should be — it is the answer.
"""
const TARGET_BAR_HEIGHT = 1.2

"""
Gap between a bar edge and the az/alt numbers printed against it, in row units.

The label heights are DERIVED from the bar rather than written down: they were once fixed at
3.3 and 0.7 against a taller bar, and when the bar was slimmed to $TARGET_BAR_HEIGHT they
stayed put and drifted away from the thing they annotate. One constant, two positions.

Small on purpose. These numbers belong to the bar and should read as attached to it; the
vertical times sit right against the ends already, and a horizontal number floating a third of
a row away looks like it belongs to the row above.
"""
const LABEL_PAD = 0.05

"Row centre of a target bar, and the two label heights that follow from its height."
const TARGET_ROW = 2.0
az_label_y()  = TARGET_ROW + TARGET_BAR_HEIGHT / 2 + LABEL_PAD
alt_label_y() = TARGET_ROW - TARGET_BAR_HEIGHT / 2 - LABEL_PAD

"""
Shortest run that gets start/end annotations, in hours.

A sliver is a few minutes wide and its six labels land on top of each other and on those of the
next run -- the pile of overlapping digits that reads as broken glyph rendering rather than as
crowding. The bar itself is still drawn; only its numbers are dropped.
"""
const MIN_LABELLED_RUN = 0.4

"""
Lowest elevation worth observing at, in degrees.

25 rather than the library's 30: it is what ASPRO defaults to, and the last five degrees are a
real part of a short night.
"""
const DEFAULT_ALT_LIMIT = 25.0

"""
Highest elevation, in degrees.

Not 90. Near the zenith the delay lines and the azimuth axis both have to move fast to track,
and most arrays have a limit; 85 is ASPRO's. It genuinely bites — Vega transits at 85.4 from
CHARA, so the cap trims the top of its window rather than being a formality.
"""
const DEFAULT_ALT_MAX = 85.0

"""
One target's night: the arrays every Observe view is drawn from.

`good_alt`, `good_delay`, `good_twilight` and `good_moon` are INDEX vectors into the time grid,
not masks. An empty one means the constraint excludes the whole night; that is the opposite of
the constraint not applying, and conflating the two is how a night with no darkness at all
reports a target as observable.
"""
struct NightPlan
    name          :: String
    ra            :: Float64        # degrees
    dec           :: Float64        # degrees
    date          :: DateTime
    lst           :: Vector{Float32}
    lst_midnight  :: Float64
    ha            :: Vector{Float32}
    alt           :: Vector{Float32}
    az            :: Vector{Float32}
    good_alt      :: Vector{Int}
    good_delay    :: Vector{Int}
    good_twilight :: Vector{Int}
    good_moon     :: Vector{Int}
    moon_sep      :: Vector{Float32}
    moon_fli      :: Float64
    delay_applied :: Bool           # false ⇒ good_delay is "not checked", not "all fine"
    alt_limit     :: Float64
    alt_max       :: Float64
    # The array and its POPs. Carried on the plan because the chart has to say what it was made
    # with: the same target on the same night gives a different window under different POPs, so
    # a Gantt without them is not reproducible from what it shows.
    config        :: Vector{Int}
    pop           :: Vector{Int}
    telescopes    :: Vector{String}
    # One entry per baseline, `(; name, good)`. Empty unless the delay check ran: this is the
    # detailed view's whole content, and one baseline can close the night on its own.
    baselines     :: Vector{NamedTuple}
end

"Indices where every supplied constraint holds at once — the hours actually worth booking."
observable_indices(p::NightPlan) =
    intersect(p.good_alt, p.good_delay, p.good_twilight, p.good_moon)

"Hours per night that survive every constraint, given the grid step in minutes."
observable_hours(p::NightPlan, step_minutes::Integer = 1) =
    length(observable_indices(p)) * step_minutes / 60

"""
    night_plan(facility, name, ra, dec, date; kwargs...) -> NightPlan

Everything about one target on one night: altitude, azimuth, and each constraint separately.

`ra` and `dec` are in DEGREES. That is worth stating because the library took RA in hours until
0.11 while documenting degrees, and a plan that is 15× out in right ascension looks plausible
rather than broken.

The delay-line check is OPT-IN, and that is not timidity. It depends on a POP configuration
that has to be searched for, and applying it with an unsearched one reports far less time than
is really available — or none at all. Measured on Vega from CHARA, 2026-06-21:

| array | POPs | minutes inside the delay limits |
|---|---|---|
| all six | 1,1,1,1,1,1 | 0 |
| S1 S2 E1 E2 | 1,1,1,1,1,1 | 119 |
| S1 S2 E1 E2 | `best_pop` → 1,3,3,4,1,1 | 380 |

A panel that applied it by default would report every target as unobservable and read as
broken. So `use_delay = false` leaves `good_delay` covering the whole night and records
`delay_applied = false`, which the UI must say out loud rather than let pass for a constraint
that was met.
"""
function night_plan(facility, name::AbstractString, ra::Real, dec::Real, date::DateTime;
                    config        = nothing,
                    pop           = nothing,
                    use_delay     ::Bool = false,
                    alt_limit     ::Real = DEFAULT_ALT_LIMIT,
                    alt_max       ::Real = DEFAULT_ALT_MAX,
                    moon_min_sep  ::Real = 30.0,
                    dark_offset   ::Real = 0.0,
                    step_minutes  ::Integer = 1,
                    delay_length  = nothing)
    f = facility isa AbstractString ? read_facility_file(String(facility)) : facility

    obs = night_observability(f, Float64(ra), Float64(dec), date;
                              alt_limit = Float64(alt_limit), alt_max = Float64(alt_max),
                              moon_min_sep = Float64(moon_min_sep),
                              dark_offset = Float64(dark_offset),
                              step_minutes = Int(step_minutes))

    cfg = config === nothing ? ones(Int, f.ntel) : Int.(config)
    pp  = pop === nothing ? ones(Int, f.ntel) : Int.(pop)
    # The delay check is also CHARA-only unless the caller supplies its own POP array and
    # airpath, so a facility without them reports "not checked" rather than a wrong answer.
    all_idx = collect(1:length(obs.ha))
    good_delay, applied, blines = if use_delay
        try
            (Int.(in_delay(f, Float64(dec), obs.ha, cfg, pp; delay_length).good_delay), true,
             baseline_delay_windows(f, Float64(dec), obs.ha, cfg, pp; delay_length))
        catch
            (all_idx, false, NamedTuple[])
        end
    else
        (all_idx, false, NamedTuple[])
    end

    return NightPlan(String(name), Float64(ra), Float64(dec), date,
                     obs.lst, obs.lst_midnight, obs.ha, obs.alt, obs.az,
                     Int.(obs.good_alt), good_delay,
                     Int.(obs.good_twilight), Int.(obs.good_moon),
                     obs.moon_sep, obs.moon_fli, applied,
                     Float64(alt_limit), Float64(alt_max),
                     cfg, pp, String.(f.tel_names), blines)
end

function Base.show(io::IO, p::NightPlan)
    h = round(observable_hours(p); digits = 2)
    print(io, "NightPlan(", p.name, " on ", Dates.format(p.date, "yyyy-mm-dd"),
          ": ", h, " h observable, moon ", round(Int, 100 * p.moon_fli), "% lit)")
end

"""
    plan_rows(plans; step_minutes = 1) -> String

One `name\\tha\\thb\\thours\\tminsep` row per target, for the panel's summary table.

`ha`/`hb` are the first and last observable LST in hours, so a table can show the window
without the caller re-deriving it from the index vectors.
"""
function plan_rows(plans, step_minutes::Integer = 1)
    rows = String[]
    for p in plans
        idx = observable_indices(p)
        if isempty(idx)
            push!(rows, join((p.name, "", "", "0.0",
                              string(round(minimum(p.moon_sep); digits = 1))), "\t"))
        else
            push!(rows, join((p.name,
                              string(round(p.lst[first(idx)]; digits = 2)),
                              string(round(p.lst[last(idx)];  digits = 2)),
                              string(round(observable_hours(p, step_minutes); digits = 2)),
                              string(round(minimum(p.moon_sep[idx]); digits = 1))), "\t"))
        end
    end
    return join(rows, "\n")
end


"""
    baseline_delay_windows(facility, dec, ha, config, pop; delay_length) -> Vector{NamedTuple}

Where each BASELINE is inside its delay limits, as `(; name, good)`.

`in_delay` returns only the AND over every baseline — one answer for the array. That is the
right thing for "can I observe now", and the wrong thing for "why not": a single baseline out
of fifteen can close the whole night, and the summary cannot say which. ASPRO's detailed view
is exactly this, one row per baseline, which is what makes it worth having.

The limit is the smaller of the two telescopes' delay lengths, matching `in_delay`.
"""
function baseline_delay_windows(facility, dec::Real, ha::AbstractVector, config, pop;
                                delay_length = nothing)
    f = facility isa AbstractString ? read_facility_file(String(facility)) : facility
    r = in_delay(f, Float64(dec), Float32.(ha), Int.(config), Int.(pop); delay_length)

    dlens = if delay_length !== nothing
        fill(Float32(delay_length), f.ntel)
    elseif !isempty(f.delay_lengths)
        Float32.(f.delay_lengths)
    else
        fill(45.7f0, f.ntel)
    end

    out = NamedTuple[]
    for b in 1:r.nbaselines
        dmax = min(dlens[r.baseline_stations[1, b]], dlens[r.baseline_stations[2, b]])
        good = [t for t in axes(r.delay_carts, 2)
                if -dmax <= r.delay_carts[b, t] <= dmax]
        # The cart position itself, not just whether it is inside: the Gantt answers "in delay
        # or not", and the delay plot answers "by how much" — a baseline running at 44 m of a
        # 45.7 m limit is feasible and about to stop being so.
        push!(out, (; name = String(r.baseline_names[b]),
                      good, carts = Float64.(r.delay_carts[b, :]), limit = Float64(dmax)))
    end
    return out
end

"""
    best_pops(facility, dec, ha, config; n = 5) -> Vector{NamedTuple}

The POP configurations that keep a target inside the delay limits longest.

`score` is the count of grid steps, so it equals minutes only when the grid step is one minute
— which is the default, and the reason the panel should state its step rather than assume it.
"""
function best_pops(facility, dec::Real, ha::AbstractVector, config; n::Integer = 5)
    f = facility isa AbstractString ? read_facility_file(String(facility)) : facility
    res = best_pop(f, Float64(dec), Float32.(ha), Int.(config))
    return collect(Iterators.take(res, n))
end

"""
    pop_rows(pops) -> String

`pop\tscore` rows for the panel's POP table.
"""
pop_rows(pops) = join((join((join(string.(p.pop), " "), string(p.score)), "\t") for p in pops), "\n")


# ─────────────────────────────────────────────────────────────────────────────
# Gantt geometry
# ─────────────────────────────────────────────────────────────────────────────
#
# The chart described as data, so the Makie drawing and the matplotlib original are two
# renderings of ONE description rather than two implementations of one idea. That is what makes
# "faithful port" checkable: the test compares what matplotlib actually drew against this, and
# a divergence is a failing assertion instead of a chart nobody notices is different.
#
# Coordinates are LST hours, unwrapped so the night is monotonically increasing across the
# 24→0 crossing. matplotlib works in dates and Makie in numbers; both convert from here.

"One bar of the chart: a contiguous run of a constraint, in LST hours."
struct GanttBar
    x0     :: Float64
    x1     :: Float64
    y      :: Float64
    height :: Float64
    color  :: String
    label  :: String     # only the first run of a category carries one
end

"One annotation: az and alt at the ends of an observing block, and the times themselves."
struct GanttLabel
    x    :: Float64
    y    :: Float64
    text :: String
    kind :: Symbol       # :time, :az or :alt
    side :: Symbol       # :start or :end of the block, which decides which way it is anchored
end

"""
    unwrap_lst(lst) -> Vector{Float64}

LST made monotonically increasing across the 24→0 crossing.

A night that starts at 22 h and ends at 4 h is six hours long, not minus eighteen; every bar
position downstream depends on getting that right.
"""
function unwrap_lst(lst)
    v = Float64.(collect(lst))
    for i in 2:length(v)
        while v[i] < v[i-1]
            v[i] += 24.0
        end
    end
    return v
end

"""
    pop_label(plan) -> String

The array and its POPs as one line, e.g. `"S1:P1 S2:P3 E1:P3 E2:P4"`.

Only the telescopes actually in the array: a POP for a telescope that is not observing is not a
setting, it is noise in a label that has to be read at a glance.
"""
function pop_label(p::NightPlan)
    isempty(p.config) && return ""
    on = findall(>(0), p.config)
    isempty(on) && return ""
    parts = String[]
    for i in on
        name = i <= length(p.telescopes) ? p.telescopes[i] : string("T", i)
        push!(parts, string(name, ":P", i <= length(p.pop) ? p.pop[i] : 1))
    end
    return join(parts, " ")
end

"""
    gantt_geometry(plan; detailed = false) -> (; bars, labels, bands, rows, midnight, xlim, ymax)

The whole chart as data: background bands, one bar per contiguous run, and the annotations.

`detailed` changes what is drawn, not merely how, and summary is the default because it answers
the question that gets asked first: when can I observe this, at all.

Detailed is ASPRO's, and it is about DELAYS: one row per baseline, plus altitude and moon. That
is the view worth having when the answer is "you can't", because the summary cannot say why —
a single baseline out of fifteen closes the whole night, and only a per-baseline row names it.
"""
function gantt_geometry(p::NightPlan; detailed::Bool = false, show_alt = nothing)
    # `show_alt` was the old name for this switch, and it defaulted the other way.
    detailed = show_alt === nothing ? detailed : Bool(show_alt)
    lst = unwrap_lst(p.lst)
    mid = Float64(p.lst_midnight)
    !isempty(lst) && mid < lst[1] && (mid += 24.0)

    bars = GanttBar[]
    labels = GanttLabel[]

    # Background twilight bands: nested, darkening inward, spanning the full height.
    bands = GanttBar[]
    if !isempty(lst)
        for (pad, colour) in ((1.5, "lightgray"), (2.0, "lightgray"), (3.0, "gray"))
            push!(bands, GanttBar(lst[1] + pad, lst[end] - pad, 5.0, 10.0, colour, ""))
        end
    end

    # One bar per contiguous run, so a window that is interrupted reads as interrupted rather
    # than as one long block with a hole nobody sees.
    function add!(indices, y, height, colour, lbl)
        for (k, (i0, i1)) in enumerate(index_runs(indices))
            push!(bars, GanttBar(lst[i0], lst[i1], y, height, colour, k == 1 ? lbl : ""))
        end
    end

    rows = Tuple{Float64,String}[]        # (y, tick label), bottom-up
    ymax = 10.0

    if detailed
        # Bottom row is the answer; everything above it is why. One row per baseline, so the
        # baseline that closes the night is the one whose bar is short.
        y = 2.0
        show_idx = sort(collect(observable_indices(p)))
        lbl = "Observable"
        push!(rows, (y, p.name))

        # Clear of the az number printed above the target bar, rather than a fixed 1.4: in the
        # detailed view the rows stack tightly, and at 1.4 the altitude bar's lower edge landed
        # exactly on `az_label_y()` with the digits printed across it.
        y = az_label_y() + 1.0
        # The elevation LIMITS row: when the target sits between alt_limit and alt_max, which
        # is what the two fields in the Time panel set. Named to match them.
        isempty(p.good_alt) || add!(p.good_alt, y, 0.9, "orange", "Elevation")
        push!(rows, (y, "elevation")); y += 1.1

        isempty(p.good_moon) || add!(p.good_moon, y, 0.9, "mediumpurple", "Moon sep.")
        push!(rows, (y, "moon")); y += 1.1

        # One colour per baseline, and the SAME colour Explore gives it: a baseline that closes
        # the night here is the one to go looking for in the V² plot, and having to translate
        # between two colour schemes to do that is a step that should not exist. The bar carries
        # its own name as the colour key; `update_gantt!` resolves it through
        # `baseline_color_map`, which is what the observable plots use.
        for b in p.baselines
            isempty(b.good) || add!(b.good, y, 0.9, b.name, "")
            push!(rows, (y, b.name))
            y += 1.1
        end
        ymax = max(10.0, y + 0.6)
    else
        show_idx = sort(collect(observable_indices(p)))
        lbl = "Observable"
        push!(rows, (2.0, p.name))
    end

    if !isempty(show_idx)
        add!(show_idx, 2.0, TARGET_BAR_HEIGHT, detailed ? "green" : "blue", lbl)
        # Every run is annotated, not just the outermost pair: each is one observing block, and
        # its own start and end times with their az/alt are what goes on a schedule.
        for (i0, i1) in index_runs(show_idx)
            # Slivers are drawn but not annotated: six labels inside a few minutes overprint
            # each other and the next run's, which looks like broken text rather than crowding.
            if lst[i1] - lst[i0] >= MIN_LABELLED_RUN
                azy, alty = az_label_y(), alt_label_y()
                push!(labels, GanttLabel(lst[i0], TARGET_ROW, _hhmm(lst[i0]), :time, :start))
                push!(labels, GanttLabel(lst[i1], TARGET_ROW, _hhmm(lst[i1]), :time, :end))
                push!(labels, GanttLabel(lst[i0], azy,  string(round(Int, p.az[i0])),  :az,  :start))
                push!(labels, GanttLabel(lst[i0], alty, string(round(Int, p.alt[i0])), :alt, :start))
                push!(labels, GanttLabel(lst[i1], azy,  string(round(Int, p.az[i1])),  :az,  :end))
                push!(labels, GanttLabel(lst[i1], alty, string(round(Int, p.alt[i1])), :alt, :end))
            end
        end
    end

    # The window matplotlib uses: midnight ± 8 h, which centres the night regardless of season.
    offset = 4.0
    xlim = (mid - 12 + offset, mid + 12 - offset)

    # Bands span whatever height the rows ended up needing.
    bands = [GanttBar(b.x0, b.x1, ymax/2, ymax, b.color, b.label) for b in bands]

    # The POPs only when the delay check actually ran: quoting a setting that did not affect
    # the chart would be worse than quoting none.
    subtitle = p.delay_applied ? pop_label(p) : ""

    return (; bars, labels, bands, midnight = mid, xlim, target = p.name, rows, ymax, detailed,
              subtitle, delay_applied = p.delay_applied)
end

"LST hours as `H:M`, matching the annotation format of the matplotlib original."
function _hhmm(h::Real)
    t = mod(Float64(h), 24.0)
    hh = floor(Int, t)
    # Truncated, not rounded: `gantt_onenight` formats a DateTime, and that drops the seconds
    # rather than rounding them up. Rounding here puts the two renderers a minute apart on any
    # block whose end falls past the half-minute.
    mm = floor(Int, (t - hh) * 60)
    return string(mod(hh, 24), ":", mm)
end


# ─────────────────────────────────────────────────────────────────────────────
# Delay-cart geometry (the `chara_plan` view)
# ─────────────────────────────────────────────────────────────────────────────

"""
    delay_plot_geometry(plan) -> (; curves, limit, alt, xlim, ylim, altlim, target)

The delay-cart chart as data: one curve per baseline, the ±limit lines, and the altitude trace.

Where the Gantt reduces each baseline to in-or-out, this keeps the number. That is the
difference between "you can observe" and "you can observe, with 1.7 m of cart to spare", and
only the second tells you whether a small change of POP or array would help.

`alt` is carried on its own scale — the original draws it on a twin axis — so the renderer
maps it rather than mixing metres and degrees on one axis.
"""
function delay_plot_geometry(p::NightPlan)
    lst = unwrap_lst(p.lst)
    curves = [(; name = b.name, x = lst, y = b.carts) for b in p.baselines]
    limit  = isempty(p.baselines) ? 45.7 : maximum(b.limit for b in p.baselines)

    ymax = isempty(curves) ? limit :
           max(limit, maximum(maximum(abs, c.y) for c in curves)) * 1.05
    xlim = isempty(lst) ? (0.0, 24.0) : (minimum(lst), maximum(lst))

    return (; curves, limit, alt = (; x = lst, y = Float64.(p.alt)),
              xlim, ylim = (-ymax, ymax), altlim = (-5.0, 90.0),
              alt_limit = p.alt_limit, alt_max = p.alt_max, target = p.name)
end


# ─────────────────────────────────────────────────────────────────────────────
# What a simulation is OF
# ─────────────────────────────────────────────────────────────────────────────
#
# `simulate` needs a sky, and it takes one of two kinds:
#
#     image      = a FITS path or an array — 2-D is grey, 3-D is one plane per channel
#     flat_model = a compiled model, plus `flat_params`, its current values
#
# Those are genuinely different inputs, not one input with a flag, and the panel has to say
# which it holds before it can say whether Simulate is possible at all. Validating here rather
# than at the call means a wrong file is reported when it is chosen, not after an observation
# has been computed.

"""
    simulate_source_info(kind, path) -> (; ok, kind, summary, ndims, nx, nwav, detail)

Describe the sky a simulation would use, and whether it is usable.

`kind` is `"image"`, `"cube"` or `"model"`. For the two image kinds `path` is a FITS file; for
`"model"` it is a TOML model file (see [`read_model_file`](@ref)) — or empty, meaning the model
currently held by the Modeling perspective.

`ok = false` comes with a `summary` saying why, which is the whole point: "3-D file chosen as a
grey image" is a mistake worth catching before it produces an OIFITS nobody can explain.
"""
function simulate_source_info(kind::AbstractString, path::AbstractString = "")
    k = lowercase(strip(String(kind)))
    p = String(path)

    no(msg) = (; ok = false, kind = k, summary = msg, ndims = 0, nx = 0, nwav = 0, detail = "")

    if k == "model"
        isempty(p) && return (; ok = true, kind = k,
                                summary = "the model held by the Modeling perspective",
                                ndims = 0, nx = 0, nwav = 0, detail = "")
        isfile(p) || return no("no model file at '$p'")
        m = try
            read_model_file(p)
        catch err
            return no("could not read '$(basename(p))': " * sprint(showerror, err))
        end
        ncomp = length(OITOOLS._component_names(m.model))
        return (; ok = true, kind = k,
                  summary = "$(basename(p)): $ncomp component(s), $(length(m.free)) free",
                  ndims = 0, nx = 0, nwav = 0,
                  detail = join(sort(collect(keys(m.model))), ", "))
    end

    isempty(p) && return no("no file chosen")
    isfile(p)  || return no("no file at '$p'")
    img = try
        readfits(p)
    catch err
        return no("could not read '$(basename(p))': " * sprint(showerror, err))
    end

    nd = ndims(img)
    if k == "image"
        nd == 2 || return no("'$(basename(p))' is $(nd)-D; a grey image must be 2-D " *
                             (nd == 3 ? "— choose Image cube instead" : ""))
        return (; ok = true, kind = k,
                  summary = "$(basename(p)): $(size(img,1))×$(size(img,2)) grey",
                  ndims = 2, nx = size(img, 1), nwav = 1, detail = "")
    elseif k == "cube"
        nd == 3 || return no("'$(basename(p))' is $(nd)-D; a cube must be 3-D " *
                             (nd == 2 ? "— choose Image instead" : ""))
        return (; ok = true, kind = k,
                  summary = "$(basename(p)): $(size(img,1))×$(size(img,2)) × $(size(img,3)) channels",
                  ndims = 3, nx = size(img, 1), nwav = size(img, 3), detail = "")
    end
    return no("unknown source '$kind'; expected image, cube or model")
end
