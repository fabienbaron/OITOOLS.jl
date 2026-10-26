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

One entry per telescope in the facility's own order: 0 for out, 1 for in, 2 for the delay
reference, which `selected` names with a trailing `*` (`"W2*"`).
"""
function telescope_config(selected, all)
    sel = String.(selected)
    return Int[(String(t) * "*") in sel ? 2 : String(t) in sel ? 1 : 0 for t in all]
end

"""
    central_telescope(facility, config) -> Int

The telescope of `config` (entries > 0) whose largest distance to the others, projected on the
plane of the equator, is smallest: the delay reference that keeps the other carts' swing over
a night smallest. That swing is `cos δ` times this projected distance. `0` with none selected.
"""
function central_telescope(f::FacilityConfig, config::AbstractVector{<:Integer})
    used = findall(>(0), config)
    isempty(used) && return 0
    l = f.lat * π / 180
    n, e, u = f.sta_xyz[:, 1], f.sta_xyz[:, 2], f.sta_xyz[:, 3]     # north, east, up
    x = u .* cos(l) .- n .* sin(l)
    far(r) = maximum(hypot(x[i] - x[r], e[i] - e[r]) for i in used)
    return used[argmin([far(r) for r in used])]
end

"""
Height of a target's own bar, in the chart's 0–10 row units.

Thinner than the constraint rows are tall, which leaves the az/alt annotations above and below
it clear while keeping it the most prominent row: it is the answer.
"""
const TARGET_BAR_HEIGHT = 0.77

"""
Gap between a bar edge and the az/alt numbers printed against it, in row units.

The label heights are derived from the bar, so they follow `TARGET_BAR_HEIGHT`: one constant,
two positions.

Small on purpose. These numbers belong to the bar and should read as attached to it; the
vertical times sit right against the ends already, and a horizontal number floating a third of
a row away looks like it belongs to the row above.
"""
const LABEL_PAD = 0.02

"Row centre of a target bar, and the two label heights that follow from its height."
const TARGET_ROW = 2.0
az_label_y()  = TARGET_ROW + TARGET_BAR_HEIGHT / 2 + LABEL_PAD

"How far above `az_label_y()` the POPs over a bar are centred: clear of the az numbers."
const POP_LABEL_LIFT = 0.36
alt_label_y() = TARGET_ROW - TARGET_BAR_HEIGHT / 2 - LABEL_PAD

"""
Shortest run that gets start/end annotations, in hours.

A sliver is a few minutes wide and its six labels land on top of each other -- the pile of
overlapping digits that reads as broken glyph rendering rather than as crowding. The bar
itself is still drawn; only its numbers are dropped.

Defined in the core beside `gantt_onenight`, which obeys the same rule: `test_gantt.jl`
compares the two renderers' labels string for string.
"""
const MIN_LABELLED_RUN = OITOOLS.GANTT_MIN_LABELLED_RUN

"""
Width a run needs before its az and elev numbers are printed, as a fraction of the plotted span.

The numbers read inward from the bar's own edges, so what has to be wide enough is the BAR, not
the gap to its neighbour. A fraction rather than a number of hours because the numbers are a
fixed width in pixels: what they cost in hours depends on how many hours the axis is showing.
"""
const LABEL_ROOM_FRAC = OITOOLS.GANTT_LABEL_ROOM_FRAC

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

"One stretch of a multi-POP night plan, as indices into its time grid (see `NightPlan`)."
const PlanSegment = NamedTuple{(:i0, :i1, :config, :pop, :tier),
                               Tuple{Int,Int,Vector{Int},Vector{Int},Int}}

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
    # Minutes between samples on the time grid. Carried because every duration on this plan is
    # a COUNT of samples, so the count means nothing without it -- and a caller that had to
    # remember to supply it did not: the panel samples every 10 minutes and the summary line
    # divided by 1, reporting a tenth of the observable time and reading as "not observable".
    step_minutes  :: Float64
    # The facility's longitude and IANA zone, carried so a chart can be drawn in civil local
    # time without the renderer having to re-read the configuration. `local_utc_offset` takes
    # either this or a FacilityConfig, since it reads exactly these two names.
    lon           :: Float64
    timezone      :: String
    site          :: String
    lst           :: Vector{Float32}
    # The UTC grid the LST above was sampled on, in decimal hours from `date` at 00:00 UT, so
    # it runs past 24 for the morning half of a night. Carried because LST and solar time do
    # NOT differ by a constant -- sidereal runs fast by 1/365.25 -- so an axis drawn in UTC or
    # local time cannot be relabelled from LST by an offset; it has to be mapped through the
    # pair, which is what `lst_to_utc` does.
    utc           :: Vector{Float64}
    lst_midnight  :: Float64
    # Local sidereal time at the facility at the moment this plan was computed -- the wall
    # clock in the control room, in the units the chart's x axis is drawn in. Recorded here
    # rather than read at draw time so the marker stays where the night was worked out: the
    # chart is redrawn on every zoom and pan, and a line that crept rightwards each time would
    # be telling you something the rest of the figure does not know.
    lst_now       :: Float64
    # Whether the clock above belongs on THIS chart. LST repeats every day, so the current LST
    # falls inside almost any night's span -- including one planned for next month, where a
    # line saying "you are here" is simply false. Decided once, when the plan is made, against
    # the noon-to-noon bracket that defines the night beginning on `date`.
    is_tonight    :: Bool
    ha            :: Vector{Float32}
    alt           :: Vector{Float32}
    az            :: Vector{Float32}
    good_alt      :: Vector{Int}
    # The same window under the mount's SOFT elevation ceiling, where the facility declares one.
    # Time between the two is observable and worth flagging, not lost; `gantt_geometry` draws it
    # as a separate bar and the chart dashes it.
    good_alt_soft :: Vector{Int}
    # The two halves of `good_alt`, which the detailed view draws as ASPRO does: Rise/Set is the
    # elevation limits alone, Horizon the telescopes' horizon profile alone (`has_horizon`
    # false: the facility has none).
    good_riseset  :: Vector{Int}
    good_horizon  :: Vector{Int}
    has_horizon   :: Bool
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
    # A multi-POP plan (`pop_plan`): the stretches observed, `(; i0, i1, config, pop, tier)` as
    # indices into the time grid, and the POP changes between them, `(c0, c1)`: the change
    # starts at `c0` and the new setting at `c1`, nothing being observed in between (`c0 == c1`
    # for a change that costs no time). Both empty for a plan with one POP setting, which `pop`
    # and `config` then describe.
    segments      :: Vector{PlanSegment}
    changes       :: Vector{Tuple{Int,Int}}
    # The order the panel lists the telescopes in (the facility file's `display_order`).
    display_order :: Vector{Int}
end

"Whether a plan has more than one POP setting or drops a telescope (see `pop_plan`)."
is_multi_pop(p::NightPlan) = !isempty(p.segments)

"Indices where every supplied constraint holds at once — the hours actually worth booking."
observable_indices(p::NightPlan) =
    intersect(p.good_alt, p.good_delay, p.good_twilight, p.good_moon)

"Hours per night that survive every constraint."
observable_hours(p::NightPlan) = length(observable_indices(p)) * p.step_minutes / 60

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
| S1 S2 E1 E2 | 1,1,1,1,1,1 | 266 |
| S1 S2 E1 E2 | `best_pop` → 1,2,2,3,1,1 | 388 |
| all six | `best_pop` → 1,1,2,3,1,5 | 185 |

A panel that applied it by default would report every target as unobservable and read as
broken. So `use_delay = false` leaves `good_delay` covering the whole night and records
`delay_applied = false`, which the UI must say out loud rather than let pass for a constraint
that was met.
"""
function night_plan(facility, name::AbstractString, ra::Real, dec::Real, date::DateTime;
                    config        = nothing,
                    pop           = nothing,
                    popplan       = nothing,
                    change_minutes::Real = 10,
                    use_delay     ::Bool = false,
                    alt_limit     ::Real = DEFAULT_ALT_LIMIT,
                    alt_max       ::Real = DEFAULT_ALT_MAX,
                    moon_min_sep  ::Real = 30.0,
                    dark_offset   ::Real = 0.0,
                    step_minutes  ::Integer = 1,
                    delay_length  = nothing)
    f = facility isa AbstractString ? read_facility_file(String(facility)) : facility

    # The horizon is the max over the telescopes IN USE, so the selection has to reach the
    # observability rather than only the delay check.
    cfg0 = config === nothing ? ones(Int, f.ntel) : Int.(config)
    obs = night_observability(f, Float64(ra), Float64(dec), date;
                              alt_limit = Float64(alt_limit), alt_max = Float64(alt_max),
                              config = cfg0,
                              moon_min_sep = Float64(moon_min_sep),
                              dark_offset = Float64(dark_offset),
                              step_minutes = Int(step_minutes))

    cfg = config === nothing ? ones(Int, f.ntel) : Int.(config)
    pp  = pop === nothing ? ones(Int, f.ntel) : Int.(pop)
    # The delay check is also CHARA-only unless the caller supplies its own POP array and
    # airpath, so a facility without them reports "not checked" rather than a wrong answer.
    all_idx = collect(1:length(obs.ha))
    segments = PlanSegment[]; changes = Tuple{Int,Int}[]
    good_delay, applied, blines = if use_delay && popplan isa PopPlan && !isempty(popplan.segments)
        segments, changes, inforce, settings = _plan_segments(popplan, obs.ha, Float64(step_minutes),
                                                           Float64(change_minutes))
        pp = first(settings)
        (sort(unique(reduce(vcat, [collect(sg.i0:sg.i1) for sg in segments]))), true,
         _piecewise_baselines(f, Float64(dec), obs.ha, cfg, settings, inforce; delay_length))
    elseif use_delay
        try
            (Int.(in_delay(f, Float64(dec), obs.ha, cfg, pp; delay_length).good_delay), true,
             baseline_delay_windows(f, Float64(dec), obs.ha, cfg, pp; delay_length))
        catch
            (all_idx, false, NamedTuple[])
        end
    else
        (all_idx, false, NamedTuple[])
    end

    # `now(UTC)`, not `date`: this is the clock, not the night being planned. `ra` only
    # affects the hour angle the call also returns, which is discarded here.
    utc_now = Dates.now(Dates.UTC)
    lst_now, _ = OITOOLS.hour_angle_calc(utc_now, f.lon, Float64(ra) / 15)

    # Is the clock on this night? The night that BEGINS on `date` runs from that afternoon to
    # the next morning, so the noon-to-noon bracket is the test -- a plan made at 01:00 for
    # "yesterday" is still the night in progress. Local time is UTC shifted by the longitude,
    # four minutes per degree, which is as exact as this decision needs.
    local_now = utc_now + Dates.Minute(round(Int, 4 * f.lon))
    noon = Dates.DateTime(Dates.Date(date), Dates.Time(12))
    is_tonight = noon <= local_now < noon + Dates.Day(1)

    return NightPlan(String(name), Float64(ra), Float64(dec), date, Float64(step_minutes),
                     Float64(f.lon), String(f.timezone), String(f.site),
                     obs.lst, Float64.(obs.utc), obs.lst_midnight,
                     Float64(first(lst_now)), is_tonight,
                     obs.ha, obs.alt, obs.az,
                     Int.(obs.good_alt), Int.(obs.good_alt_soft),
                     Int.(obs.good_riseset), Int.(obs.good_horizon), Bool(obs.has_horizon),
                     good_delay,
                     Int.(obs.good_twilight), Int.(obs.good_moon),
                     obs.moon_sep, obs.moon_fli, applied,
                     Float64(alt_limit), Float64(alt_max),
                     cfg, pp, String.(f.tel_names), blines,
                     segments, changes, telescope_display_order(f))
end

# A `PopPlan`'s segments on a night plan's time grid. Returns the segments as index ranges, the
# POP changes between settings (`(c0, c1)`: starting `change_minutes` before the new setting's
# first segment, at `c1`), the setting in force at every sample (0 during a change) and the
# distinct settings in order.
function _plan_segments(plan::PopPlan, ha::AbstractVector, step_minutes::Float64,
                        change_minutes::Float64)
    tol = step_minutes / 120                 # half a step, in hours
    hc = unwrap_hour_angles(ha)              # the branch `pop_plan` searched on
    segs = PlanSegment[]
    for sg in plan.segments
        idx = findall(h -> sg.ha0 - tol <= h <= sg.ha1 + tol, hc)
        isempty(idx) && continue
        push!(segs, PlanSegment((first(idx), last(idx), sg.config, sg.pop, sg.tier)))
    end
    settings = Vector{Int}[]
    starts = Int[]
    for sg in segs
        if isempty(settings) || last(settings) != sg.pop
            push!(settings, sg.pop); push!(starts, sg.i0)
        end
    end
    n = length(ha)
    inforce = zeros(Int, n)
    changes = Tuple{Int,Int}[]
    nchg = round(Int, change_minutes / step_minutes)
    for j in eachindex(settings)
        lo = j == 1 ? 1 : starts[j]
        hi = j == length(settings) ? n : max(starts[j+1] - nchg - 1, lo)
        inforce[lo:hi] .= j
        j < length(settings) && push!(changes, (hi + 1, starts[j+1]))
    end
    return segs, changes, inforce, settings
end

# Each baseline under the POP setting in force at every sample: its carts where a setting is in
# force and NaN during a change, and in delay where it is -- during a change, under the setting
# it changes to, so the Gantt rows run on through it.
function _piecewise_baselines(f, dec, ha, config, settings, inforce; delay_length = nothing)
    per = [baseline_delay_windows(f, dec, ha, config, pp; delay_length) for pp in settings]
    shown = copy(inforce)
    for i in length(shown)-1:-1:1
        shown[i] == 0 && (shown[i] = shown[i+1])
    end
    out = NamedTuple[]
    for b in eachindex(first(per))
        carts = fill(NaN, length(ha))
        good = Int[]
        for (i, j) in enumerate(inforce)
            j == 0 || (carts[i] = per[j][b].carts[i])
            jg = shown[i]
            jg > 0 && insorted(i, per[jg][b].good) && push!(good, i)
        end
        push!(out, (; name = per[1][b].name, good, carts, limit = per[1][b].limit))
    end
    return out
end

# Read once per facility file: `night_plan` asks on every Compute, and parsing the TOML again
# cost more than the rest of the plan's bookkeeping.
const DISPLAY_ORDER_CACHE = Dict{Tuple{String,Vector{String}},Vector{Int}}()

"""
    telescope_display_order(facility) -> Vector{Int}

The facility's telescopes in the order the Observing panel lists them, as indices: the facility
file's `display_order` when it names each telescope once, else the file's own order.
"""
function telescope_display_order(f::FacilityConfig)
    names = String.(f.tel_names)
    return copy(get!(() -> _read_display_order(f, names), DISPLAY_ORDER_CACHE, (f.name, names)))
end

function _read_display_order(f::FacilityConfig, names::Vector{String})
    order::Vector{String} = try
        d = OITOOLS.TOML.parsefile(OITOOLS._resolve_config(f.name))
        String[String(x) for x in get(d, "display_order", String[])]
    catch
        String[]
    end
    idx = Union{Nothing,Int}[findfirst(==(n), names) for n in order]
    (length(idx) == length(names) && !any(isnothing, idx) && allunique(idx)) ||
        return collect(eachindex(names))
    return Int[something(i) for i in idx]
end
telescope_display_order(name::AbstractString) = telescope_display_order(read_facility_file(String(name)))

"""
    pop_plan_text(plan) -> String
    parse_pop_plan_text(text) -> PopPlan

A `PopPlan` as one line of text, so the panel can hold it, pass it back and show it: segments
separated by `;`, each `ha0 ha1 config tier pop…` with `config` one digit per telescope and the
POPs in the facility's telescope order, e.g. `-3.5800 -0.9800 111111 1 1 1 3 4 1 4`. `minutes`
and `changes` are recomputed on reading: the minutes from each segment's span on the 1-minute
grid, the changes from consecutive segments whose POPs differ.
"""
pop_plan_text(plan::PopPlan) =
    join((@sprintf("%.4f %.4f %s %d %s", sg.ha0, sg.ha1, join(string.(sg.config)), sg.tier,
                   join(string.(sg.pop), " ")) for sg in plan.segments), ";")

function parse_pop_plan_text(text::AbstractString)
    segs = PopSegment[]
    for part in split(String(text), ';')
        f = split(strip(part))
        length(f) >= 5 || continue
        cfg = [parse(Int, string(c)) for c in f[3]]
        pops = parse.(Int, f[5:end])
        length(pops) == length(cfg) || throw(ArgumentError("POP plan segment \"$(part)\": $(length(cfg)) telescopes, $(length(pops)) POPs"))
        push!(segs, PopSegment(parse(Float64, f[1]), parse(Float64, f[2]), cfg, pops,
                               parse(Int, f[4])))
    end
    mins = zeros(Int, 5)
    for sg in segs
        mins[sg.tier] += round(Int, (sg.ha1 - sg.ha0) * 60) + 1
    end
    nchg = count(k -> segs[k].pop != segs[k-1].pop, 2:length(segs); init = 0)
    return PopPlan(segs, ntuple(i -> mins[i], 5), nchg)
end

"""
    pop_plan_summary(plan, facility, ra, config; clock = nothing) -> String

One line for the panel's table, in its telescope order: the first POPs, then `→` and the new
POPs at each change, and `−E1 22:34` where a telescope of `config` is dropped (`+E1` where it
comes back), at the time it happens: LST, or `clock(lst)` hours when `clock` is given. A
telescope a segment does not use is a `×` in its POPs.
"""
function pop_plan_summary(plan::PopPlan, f::FacilityConfig, ra::Real, config::AbstractVector{<:Integer};
                          clock = nothing)
    isempty(plan.segments) && return ""
    on = [i for i in telescope_display_order(f) if i <= length(config) && config[i] > 0]
    lst(h) = clock === nothing ? _hhmm(mod(Float64(ra) / 15 + h, 24)) :
                                 _hhmm(clock(Float64(ra) / 15 + h))
    pops(sg) = join((_pop_text(sg.pop, sg.config, i) for i in on), " ")
    out = pops(plan.segments[1])
    for k in 2:length(plan.segments)
        a, b = plan.segments[k-1], plan.segments[k]
        a.pop != b.pop && (out *= " → " * pops(b))
        lost = [i for i in on if b.config[i] == 0 && a.config[i] > 0]
        back = [i for i in on if b.config[i] > 0 && a.config[i] == 0]
        isempty(lost) || (out *= " −" * join(f.tel_names[lost], " −") * " " * lst(b.ha0))
        isempty(back) || (out *= " +" * join(f.tel_names[back], " +") * " " * lst(b.ha0))
        r = _reference(b.config)
        r != 0 && r != _reference(a.config) && a.pop == b.pop &&
            (out *= " Ref → " * f.tel_names[r] * " " * lst(b.ha0))
    end
    return out
end

"""
    night_pop_plan(facility, ra, dec, date; config, reference = 0, max_changes = 2,
                   max_drop = 0, change_minutes = 10, alt_limit, alt_max, moon_min_sep = 30,
                   dark_offset = 0) -> PopPlan

`pop_plan` over the hours `night_plan` counts as observable for these settings: above the
elevation limits and the horizon, in astronomical night, clear of the moon. Pass `night_plan`
the same `change_minutes`.
"""
function night_pop_plan(facility, ra::Real, dec::Real, date::DateTime;
                        config = nothing, reference::Integer = 0, max_changes::Integer = 2,
                        max_drop::Integer = 0, change_minutes::Real = 10,
                        alt_limit::Real = DEFAULT_ALT_LIMIT, alt_max::Real = DEFAULT_ALT_MAX,
                        moon_min_sep::Real = 30.0, dark_offset::Real = 0.0,
                        delay_length = nothing)
    f = facility isa AbstractString ? read_facility_file(String(facility)) : facility
    cfg = config === nothing ? ones(Int, f.ntel) : Int.(config)
    ha, ok = _night_sky(f, ra, dec, date, cfg; alt_limit, alt_max, moon_min_sep, dark_offset)
    ref = reference > 0 ? Int(reference) : something(findfirst(==(2), cfg), 0)
    return pop_plan(f, Float64(dec), ha, cfg; observable = ok, reference = ref,
                    max_changes, max_drop, change_minutes, delay_length)
end

# The night's hour angles and where `night_plan` counts the target as observable.
function _night_sky(f, ra, dec, date, cfg; alt_limit = DEFAULT_ALT_LIMIT,
                    alt_max = DEFAULT_ALT_MAX, moon_min_sep = 30.0, dark_offset = 0.0)
    obs = night_observability(f, Float64(ra), Float64(dec), date;
                              alt_limit = Float64(alt_limit), alt_max = Float64(alt_max),
                              config = cfg, moon_min_sep = Float64(moon_min_sep),
                              dark_offset = Float64(dark_offset))
    ok = falses(length(obs.ha))
    ok[intersect(obs.good_alt, obs.good_twilight, obs.good_moon)] .= true
    return obs.ha, ok
end

function Base.show(io::IO, p::NightPlan)
    h = round(observable_hours(p); digits = 2)
    print(io, "NightPlan(", p.name, " on ", Dates.format(p.date, "yyyy-mm-dd"),
          ": ", h, " h observable, moon ", round(Int, 100 * p.moon_fli), "% lit)")
end

"""
    plan_rows(plans) -> String

One `name\\tha\\thb\\thours\\tminsep` row per target, for the panel's summary table.

`ha`/`hb` are the first and last observable LST in hours, so a table can show the window
without the caller re-deriving it from the index vectors.
"""
function plan_rows(plans)
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
                              string(round(observable_hours(p); digits = 2)),
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
    coverage_pops(facility, dec, ha, config) -> Vector{NamedTuple}

The single POP settings that cover the earliest, middle and latest third of the hour angles `ha`
(the observable ones) best, as `(; pop, score, part)`: for each third, the setting with the most
time in delay over it, `score` in its minutes there. Among settings tied in a third, the ones
moving the fewest POPs from one third to the next. A setting best in several thirds is listed
once, `part` naming them all (`"early, middle"`).
"""
function coverage_pops(facility, dec::Real, ha::AbstractVector, config)
    f = facility isa AbstractString ? read_facility_file(String(facility)) : facility
    h = sort(Float64.(collect(ha)))
    out = NamedTuple{(:pop, :score, :part),Tuple{Vector{Int},Int,String}}[]
    isempty(h) && return out
    lo, hi = first(h), last(h)
    w = (hi - lo) / 3
    # each third's best settings, every one tied for the most minutes
    thirds = Tuple{String,Vector{NamedTuple{(:pop, :score),Tuple{Vector{Int},Int}}}}[]
    for (k, part) in enumerate(("early", "middle", "late"))
        a, b = lo + (k - 1) * w, k == 3 ? hi : lo + k * w
        sub = Float32[x for x in h if a <= x <= b]
        r = best_pop(f, Float64(dec), sub, Int.(config); n_best = typemax(Int))
        isempty(r) && continue
        push!(thirds, (part, [q for q in r if q.score == r[1].score]))
    end
    isempty(thirds) && return out
    # the picks that move the fewest POPs from third to third
    active = findall(>(0), config)
    moved(a, b) = count(i -> a[i] != b[i], active)
    cost = [zeros(Int, length(thirds[1][2]))]; back = [Int[]]
    for t in 2:length(thirds)
        prev, cur = thirds[t-1][2], thirds[t][2]
        cb = Int[]; bb = Int[]
        for q in cur
            v, i = findmin(i -> cost[t-1][i] + moved(prev[i].pop, q.pop), eachindex(prev))
            push!(cb, v); push!(bb, i)
        end
        push!(cost, cb); push!(back, bb)
    end
    choice = zeros(Int, length(thirds)); choice[end] = argmin(cost[end])
    for t in length(thirds):-1:2
        choice[t-1] = back[t][choice[t]]
    end
    for (t, (part, qs)) in enumerate(thirds)
        q = qs[choice[t]]
        j = findfirst(o -> o.pop == q.pop, out)
        if j === nothing
            push!(out, (; pop = q.pop, score = q.score, part = part))
        else
            out[j] = (; pop = out[j].pop, score = out[j].score + q.score, part = out[j].part * ", " * part)
        end
    end
    return out
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

The array and its POPs as one line, in the panel's telescope order, e.g.
`"E1:P3 E2:P4 S1:P1 S2:P3"`.

Only the telescopes actually in the array: a POP for a telescope that is not observing is not a
setting, it is noise in a label that has to be read at a glance.
"""
function pop_label(p::NightPlan)
    isempty(p.config) && return ""
    on = [i for i in p.display_order if i <= length(p.config) && p.config[i] > 0]
    isempty(on) && return ""
    parts = String[]
    for i in on
        name = i <= length(p.telescopes) ? p.telescopes[i] : string("T", i)
        push!(parts, string(name, ":P", i <= length(p.pop) ? p.pop[i] : 1, p.config[i] == 2 ? "*" : ""))
    end
    return join(parts, " ")
end

# One telescope's POP as the labels write it: `5`, `5*` for the delay reference, `×` unused.
_pop_text(pop, config, i) = config[i] > 0 ? string(pop[i], config[i] == 2 ? "*" : "") : "×"

# The delay reference of a configuration (its entry 2), or 0.
_reference(config) = something(findfirst(==(2), config), 0)

"""
    reference_change(plan, k) -> Int

The telescope that becomes the delay reference at segment `k` when it has the same POPs as
segment `k-1` (a change of reference with no POP change); `0` otherwise.
"""
function reference_change(p::NightPlan, k::Integer)
    k == 1 && return 0
    sg, prev = p.segments[k], p.segments[k-1]
    r = _reference(sg.config)
    return r != 0 && r != _reference(prev.config) && sg.pop == prev.pop ? r : 0
end

"""
    segment_label(plan, segment) -> String

One segment of a multi-POP plan: the POPs of the selected telescopes in the panel's order, `×`
for one the segment does not use and `*` after the delay reference's, e.g. `"× 5* 3 1 5 4"`.
"""
function segment_label(p::NightPlan, sg)
    order = [i for i in p.display_order if i <= length(p.config) && p.config[i] > 0]
    return join((_pop_text(sg.pop, sg.config, i) for i in order), " ")
end

"""
    transition_label(plan, k) -> String

The telescopes segment `k` loses (`−W2`) and regains (`+E1`) from segment `k-1`, and
`Ref → E2` where the delay reference changes with no POP change (`reference_change`); `""` for
the first segment, which loses nothing: one the plan never had is a `×` in its POPs. A reference
change with a POP change shows in the new POPs instead, `*` marking the reference.
"""
function transition_label(p::NightPlan, k::Integer)
    k == 1 && return ""
    sg, prev = p.segments[k], p.segments[k-1]
    on = [i for i in p.display_order if i <= length(p.config) && p.config[i] > 0]
    r = reference_change(p, k)
    return join(vcat(["−" * p.telescopes[i] for i in on if sg.config[i] == 0 && prev.config[i] > 0],
                     ["+" * p.telescopes[i] for i in on if sg.config[i] > 0 && prev.config[i] == 0],
                     r == 0 ? String[] : ["Ref → " * p.telescopes[r]]), " ")
end

"""
    segment_changes(plan, k) -> Vector{Bool}

For the POPs `segment_label` prints for segment `k` (the selected telescopes, in the panel's
order), which differ from the setting before the last POP change. All false under the first
setting, and for a telescope the segment does not use.
"""
function segment_changes(p::NightPlan, k::Integer)
    sg = p.segments[k]
    order = [i for i in p.display_order if i <= length(p.config) && p.config[i] > 0]
    j = k - 1
    while j >= 1 && p.segments[j].pop == sg.pop
        j -= 1
    end
    j == 0 && return collect(falses(length(order)))
    prev = p.segments[j].pop
    return [sg.config[i] > 0 && sg.pop[i] != prev[i] for i in order]
end

"""
    multi_pop_summary(plan) -> String

A multi-POP plan in one line: segments, POP changes, hours observed and how many of them with
every selected telescope.
"""
function multi_pop_summary(p::NightPlan)
    nseg = length(p.segments)
    nchg = length(p.changes)
    n = count(>(0), p.config)
    full = sum((sg.i1 - sg.i0 + 1 for sg in p.segments if sg.tier == 1); init = 0)
    h(x) = @sprintf("%.1f", x * p.step_minutes / 60)
    return string(nseg, nseg == 1 ? " segment · " : " segments · ", nchg,
                  nchg == 1 ? " change · " : " changes · ", h(length(observable_indices(p))),
                  " h (", h(full), " h with ", n, " telescopes)")
end

"""
One row of a night's POP timeline (`pop_plan_timeline`).

`kind` is `:setting` (observed under `pop`/`config`), `:change` (the carts move; `moved` lists
`(telescope, from, to)`), `:drop` / `:readd` (an instant: `tels` stop or start being used),
`:ref` (an instant: `tels` becomes the delay reference with the POPs unchanged),
or `:gap` (dark, nothing observed; `why` is `:elevation`, `:moon` or `:delay`, the first of
them that excludes it, in that order). `i0:i1` are the samples of the night's grid it covers, and
`t0`/`t1` its start and end on the chosen clock, unwrapped (they may pass 24); `t1` is the start
of the next sample, so consecutive rows share an edge, and `NaN` for an instant.
"""
struct TimelineRow
    kind   :: Symbol
    i0     :: Int
    i1     :: Int
    t0     :: Float64
    t1     :: Float64
    ha0    :: Float64
    ha1    :: Float64
    pop    :: Vector{Int}
    config :: Vector{Int}
    moved  :: Vector{Tuple{Int,Int,Int}}
    tels   :: Vector{Int}
    hours  :: Float64
    why    :: Symbol
end

"""
    pop_plan_timeline(p::NightPlan; system = :lst) -> Vector{TimelineRow}

The dark hours of `p` as consecutive rows: each stretch observed under one POP setting and set
of telescopes, each POP change, each telescope dropped or added, and the stretches where nothing
is observed. Times on the clock `system` (`:lst`, `:utc`, `:local`), as the Gantt labels them.
Empty when the night has no astronomical darkness.
"""
function pop_plan_timeline(p::NightPlan; system::Symbol = :lst)
    rows = TimelineRow[]
    isempty(p.good_twilight) && return rows
    n = length(p.lst)
    lo, hi = extrema(p.good_twilight)
    obs = falses(n); obs[observable_indices(p)] .= true
    seg = zeros(Int, n)
    for (k, sg) in enumerate(p.segments)
        seg[sg.i0:sg.i1] .= k
    end
    chg = falses(n)
    for (c0, c1) in p.changes
        chg[c0:c1-1] .= true
    end
    alt = falses(n); alt[p.good_alt] .= true
    moon = falses(n); moon[p.good_moon] .= true
    function state(i)
        chg[i] && return (:change, Int[], Int[])
        obs[i] || return (!alt[i] ? :elevation : !moon[i] ? :moon : :delay, Int[], Int[])
        seg[i] > 0 && return (:setting, p.segments[seg[i]].config, p.segments[seg[i]].pop)
        return (:setting, p.config, p.pop)
    end

    x = unwrap_lst(p.lst)
    edge(i) = i < n ? x[i+1] : x[i] + (n > 1 ? x[i] - x[i-1] : p.step_minutes / 60)
    clock(v) = _gantt_clock(p, v, system)
    hours(i0, i1) = (i1 - i0 + 1) * p.step_minutes / 60
    row(kind, i0, i1; pop = Int[], config = Int[], moved = Tuple{Int,Int,Int}[], tels = Int[],
        why = :none) =
        TimelineRow(kind, i0, i1, clock(x[i0]), kind in (:drop, :readd, :ref) ? NaN : clock(edge(i1)),
                    kind === :setting ? Float64(p.ha[i0]) : NaN,
                    kind === :setting ? Float64(p.ha[i1]) : NaN,
                    pop, config, moved, tels, kind in (:drop, :readd, :ref) ? 0.0 : hours(i0, i1), why)

    # runs of one state over the dark hours
    runs = Tuple{Int,Int,Any}[]
    i0 = lo; s0 = state(lo)
    for i in lo+1:hi+1
        s = i <= hi ? state(i) : nothing
        if s != s0
            push!(runs, (i0, i - 1, s0)); i0 = i; s0 = s
        end
    end

    used = [i for i in p.display_order if i <= length(p.config) && p.config[i] > 0]
    settings = [k for (k, r) in enumerate(runs) if r[3][1] === :setting]
    gapwhy = (:elevation, :moon, :delay)
    prev = nothing                                  # the last setting's (config, pop)
    for (k, (a, b, s)) in enumerate(runs)
        kind, cfg, pop = s
        if kind === :change
            before = something(findlast(<(k), settings), 0)
            after  = something(findfirst(>(k), settings), 0)
            moved = Tuple{Int,Int,Int}[]
            if before > 0 && after > 0
                pa, pb = runs[settings[before]][3][3], runs[settings[after]][3][3]
                moved = [(t, pa[t], pb[t]) for t in used if pa[t] != pb[t]]
            end
            push!(rows, row(:change, a, b; moved))
        elseif kind in gapwhy
            push!(rows, row(:gap, a, b; why = kind))
        else
            if prev !== nothing && prev[1] != cfg
                lost = [t for t in used if prev[1][t] > 0 && cfg[t] == 0]
                back = [t for t in used if prev[1][t] == 0 && cfg[t] > 0]
                isempty(lost) || push!(rows, row(:drop, a, a; tels = lost))
                isempty(back) || push!(rows, row(:readd, a, a; tels = back))
                r = _reference(cfg)
                r != 0 && r != _reference(prev[1]) && prev[2] == pop &&
                    push!(rows, row(:ref, a, a; tels = [r]))
            end
            push!(rows, row(:setting, a, b; pop = copy(pop), config = copy(cfg)))
            prev = (cfg, pop)
        end
    end
    return rows
end

"What a timeline gap row says, by the constraint that excludes it."
const GAP_LABELS = Dict(:elevation => "elevation-restricted", :moon => "too close to the Moon",
                        :delay => "out of delay")

# A duration in hours as `h:mm`.
function _hmm(h::Real)
    m = round(Int, 60 * Float64(h))
    return string(div(m, 60), ":", lpad(mod(m, 60), 2, '0'))
end

# A timeline row as the fields the panel shows: from, to, what, POPs, HA, duration, note.
function timeline_fields(r::TimelineRow, p::NightPlan)
    order = [i for i in p.display_order if i <= length(p.config) && p.config[i] > 0]
    to = isnan(r.t1) ? "" : _hhmm(r.t1)
    if r.kind === :setting
        pops = join((_pop_text(r.pop, r.config, i) for i in order), " ")
        return (_hhmm(r.t0), to, "setting", pops, @sprintf("%+.1f … %+.1f", r.ha0, r.ha1),
                _hmm(r.hours), "")
    elseif r.kind === :change
        note = join((string(p.telescopes[t], " ", a, "→", b) for (t, a, b) in r.moved), " · ")
        return (_hhmm(r.t0), to, "change", "", "", _hmm(r.hours), note)
    elseif r.kind === :gap
        return (_hhmm(r.t0), to, "gap", "", "", _hmm(r.hours), GAP_LABELS[r.why])
    else
        return (_hhmm(r.t0), "", String(r.kind), "", "", "", join(p.telescopes[r.tels], " "))
    end
end

"""
    timeline_text(rows, p::NightPlan; system = :lst) -> String

A timeline as a plain-text table, one row per line, the telescopes named in the header in the
panel's order.
"""
function timeline_text(rows::AbstractVector{TimelineRow}, p::NightPlan; system::Symbol = :lst)
    order = [i for i in p.display_order if i <= length(p.config) && p.config[i] > 0]
    what(f) = f[3] == "setting" ? f[4] :
              f[3] == "change"  ? "POP change  " * f[7] :
              f[3] == "gap"     ? f[7] :
              f[3] == "drop"    ? "drop " * f[7] :
              f[3] == "ref"     ? "Ref → " * f[7] : "add " * f[7]
    w = max(30, 2 + maximum((length(what(timeline_fields(r, p))) for r in rows); init = 0))
    lines = [rpad(uppercase(String(system)), 6) * rpad("from", 7) * rpad("to", 7) *
             rpad("POPs (" * join(p.telescopes[order], " ") * ")", w) * rpad("HA", 15) * "time"]
    for r in rows
        f = timeline_fields(r, p)
        push!(lines, rpad("", 6) * rpad(f[1], 7) * rpad(f[2], 7) * rpad(what(f), w) *
                     rpad(f[5], 15) * f[6])
    end
    return join(rstrip.(lines), "\n")
end

"""
    lst_to_utc(p::NightPlan, x) -> Float64

Decimal UTC hours from `p.date` at 00:00 UT for a point `x` on the chart's LST axis.

LST and solar time run at different rates -- sidereal gains about four minutes a day -- so this
is an interpolation through the night's own (LST, UTC) pairs rather than an offset. Outside the
sampled window it extrapolates from the end segment, which is what keeps a tick just past
sunrise from folding back into the evening.
"""
function lst_to_utc(p::NightPlan, x::Real)
    lst = unwrap_lst(p.lst)
    n = length(lst)
    (n == 0 || length(p.utc) != n) && return NaN
    n == 1 && return p.utc[1]
    x <= lst[1]   && return p.utc[1]   + (x - lst[1])   * (p.utc[2] - p.utc[1])   / (lst[2] - lst[1])
    x >= lst[end] && return p.utc[end] + (x - lst[end]) * (p.utc[end] - p.utc[end-1]) /
                                                          (lst[end] - lst[end-1])
    i = searchsortedlast(lst, x)
    i = clamp(i, 1, n - 1)
    d = lst[i+1] - lst[i]
    return d == 0 ? p.utc[i] : p.utc[i] + (x - lst[i]) * (p.utc[i+1] - p.utc[i]) / d
end

"""
    utc_to_lst(p::NightPlan, u) -> Float64

The inverse of [`lst_to_utc`](@ref): where a UTC hour falls on the chart's LST axis.
"""
function utc_to_lst(p::NightPlan, u::Real)
    lst = unwrap_lst(p.lst)
    n = length(lst)
    (n == 0 || length(p.utc) != n) && return NaN
    n == 1 && return lst[1]
    u <= p.utc[1]   && return lst[1]   + (u - p.utc[1])   * (lst[2] - lst[1]) / (p.utc[2] - p.utc[1])
    u >= p.utc[end] && return lst[end] + (u - p.utc[end]) * (lst[end] - lst[end-1]) /
                                                            (p.utc[end] - p.utc[end-1])
    i = searchsortedlast(p.utc, u)
    i = clamp(i, 1, n - 1)
    d = p.utc[i+1] - p.utc[i]
    return d == 0 ? lst[i] : lst[i] + (u - p.utc[i]) * (lst[i+1] - lst[i]) / d
end

"""
    local_utc_offset(facility, when::DateTime) -> Float64

Hours to add to UTC to get civil local time at `facility`.

There is no time-zone database here -- TimeZones.jl is not a dependency, and pulling one in for
one observatory would be a poor trade -- so the rule is written out for the zone CHARA sits in
and everything else falls back to MEAN SOLAR time from the longitude. That fallback is honest
rather than wrong: it is local time in the original sense, and it is what a facility without a
declared zone can support. `FacilityConfig` declares one through `timezone`.

The implemented zone is America/Los_Angeles under the post-2007 US rule: UTC-8, and UTC-7 from
02:00 local on the second Sunday in March to 02:00 local on the first Sunday in November.
"""
function local_utc_offset(facility, when::DateTime)
    tz = hasproperty(facility, :timezone) ? String(getproperty(facility, :timezone)) : ""
    if tz == "America/Los_Angeles"
        y = Dates.year(when)
        # Dates.tonext is inclusive-from-the-next-day, so start the search a day early to let
        # the 1st itself be the answer when it is already the right weekday.
        dst_on  = Dates.tonext(Dates.Date(y, 3, 1)  - Dates.Day(1), Dates.Sunday) + Dates.Day(7)
        dst_off = Dates.tonext(Dates.Date(y, 11, 1) - Dates.Day(1), Dates.Sunday)
        # The switch times are quoted in LOCAL clock, so compare in the offset in force either
        # side of them: -8 going in, -7 coming out.
        on_utc  = DateTime(dst_on)  + Dates.Hour(2 + 8)
        off_utc = DateTime(dst_off) + Dates.Hour(2 + 7)
        return (on_utc <= when < off_utc) ? -7.0 : -8.0
    end
    return hasproperty(facility, :lon) ? Float64(getproperty(facility, :lon)) / 15 : 0.0
end

"""
    gantt_time_label(system, facility) -> String

The x-axis label for a Gantt drawn in `system`, one of `:lst`, `:utc` or `:local`.
"""
function gantt_time_label(system::Symbol, facility = nothing)
    system === :utc && return "UTC"
    if system === :local
        # The SITE, not the IANA zone: "Mount Wilson" is what an observer calls the place, and
        # "America/Los_Angeles" is a database key that happens to contain it.
        site = facility !== nothing && hasproperty(facility, :site) ?
               String(getproperty(facility, :site)) : ""
        tz   = facility !== nothing && hasproperty(facility, :timezone) ?
               String(getproperty(facility, :timezone)) : ""
        isempty(site) || return "Local (" * site * ")"
        return isempty(tz) ? "Local mean solar time" : "Local (" * tz * ")"
    end
    return "LST (h)"
end

"""
    gantt_geometry(plan; detailed = false) -> (; bars, labels, bands, softbars, rows, midnight, now, transit, xlim, ymax)

The whole chart as data: background bands, one bar per contiguous run, and the annotations.

`detailed` changes what is drawn, not merely how, and summary is the default because it answers
the question that gets asked first: when can I observe this, at all.

Detailed is ASPRO's, and it is about DELAYS: one row per baseline, plus altitude and moon. That
is the view worth having when the answer is "you can't", because the summary cannot say why —
a single baseline out of fifteen closes the whole night, and only a per-baseline row names it.
"""
function gantt_geometry(p::NightPlan; detailed::Bool = false, show_alt = nothing,
                        time_system::Symbol = :lst)
    # `show_alt` was the old name for this switch, and it defaulted the other way.
    detailed = show_alt === nothing ? detailed : Bool(show_alt)
    lst = unwrap_lst(p.lst)
    mid = Float64(p.lst_midnight)
    !isempty(lst) && mid < lst[1] && (mid += 24.0)

    # The clock, unwrapped into the same turn of LST the night is drawn in. `NaN` when it does
    # not land inside the night at all -- which is most of the time, since the chart spans the
    # dark hours and not a whole sidereal day. Clamping it to an edge instead would put a line
    # on the chart that says "now" while meaning "not now, and off to one side".
    nowlst = Float64(p.lst_now)
    !isempty(lst) && nowlst < lst[1] && (nowlst += 24.0)
    inside = p.is_tonight && !isempty(lst) && lst[1] <= nowlst <= lst[end]
    nowlst = inside ? nowlst : NaN

    bars = GanttBar[]
    labels = GanttLabel[]

    # Background twilight bands, darkening inward -- and they mean the Sun's depression, not a
    # decoration. The axis spans the NAUTICAL window (`night_window`'s 102 degree default, which
    # is what `night_observability` samples over), and `good_twilight` is the ASTRONOMICAL one
    # at 108, so the two together are the chart's whole dark-time story: light grey where the
    # Sun is between 12 and 18 degrees down, darker where it is past 18.
    #
    # They used to be rectangles inset 1.5, 2.0 and 3.0 HOURS from the ends of the night,
    # nested so they darkened inward. Nothing astronomical was consulted, and nested darkening
    # bands on an observing chart are read as twilight -- so the plot asserted a boundary that
    # was never computed, and moved it with the length of the night rather than with the Sun.
    bands = GanttBar[]
    if !isempty(lst)
        push!(bands, GanttBar(lst[1], lst[end], 5.0, 10.0, "lightgray", ""))
        if !isempty(p.good_twilight)
            i0, i1 = extrema(p.good_twilight)
            push!(bands, GanttBar(lst[i0], lst[i1], 5.0, 10.0, "gray", ""))
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
    # A multi-POP plan: each segment's label (x, y, text), the POP changes, and in the detailed
    # view a line through every row at each change.
    # (x, y, text, which of its POPs changed from the previous segment's)
    seglabels = Tuple{Float64,Float64,String,Vector{Bool}}[]
    # Written vertically in orange inside the target bar: "POP" just after where each change
    # starts (side :start), "DROP" just before a telescope is lost (side :end).
    marks = Tuple{Float64,String,Symbol,Float64}[]          # (x, text, side, y)
    # (x, y, POPs, which changed, whether it is a change): the first setting at its first
    # segment and each change where it starts, as "→" and the new POPs.
    poplabels = Tuple{Float64,Float64,String,Vector{Bool},Bool}[]
    switches = Float64[]
    span = isempty(lst) ? 1.0 : lst[end] - lst[1] + 0.5
    # Where a segment's bar ends: at its last sample, or at the next segment's first when the
    # two touch (a drop with no POP change) or only a POP change is between them. The change is
    # not drawn as a gap; "POP" marks where it starts.
    touches(k) = k < length(p.segments) &&
                 (p.segments[k+1].i0 == p.segments[k].i1 + 1 ||
                  any(c -> c[1] == p.segments[k].i1 + 1 && c[2] == p.segments[k+1].i0, p.changes))
    seg_x1(k) = touches(k) ? lst[p.segments[k+1].i0] : lst[p.segments[k].i1]
    # A segment's label in the POPs row: its POPs where they fit, on a short one only the
    # telescopes it loses or regains.
    function seglabel(k)
        sg = p.segments[k]
        w = lst[sg.i1] - lst[sg.i0]
        w >= 2 * LABEL_ROOM_FRAC * span && return segment_label(p, sg)
        w >= 0.5 * LABEL_ROOM_FRAC * span && return transition_label(p, k)
        return ""
    end

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
        # ASPRO's two rows: Rise/Set is the Time panel's elevation limits, Horizon the profile
        # of the telescopes in use.
        isempty(p.good_riseset) || add!(p.good_riseset, y, 0.9, "orange", "Rise/Set")
        push!(rows, (y, "Rise/Set")); y += 1.1
        if p.has_horizon
            isempty(p.good_horizon) || add!(p.good_horizon, y, 0.9, "sienna", "Horizon")
            push!(rows, (y, "Horizon")); y += 1.1
        end

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
        # The POPs row: one block per segment, labelled with its POPs and its drops.
        if is_multi_pop(p)
            for (k, sg) in enumerate(p.segments)
                push!(bars, GanttBar(lst[sg.i0], seg_x1(k), y, 0.9,
                                     sg.tier == 1 ? "lightsteelblue" : "lightblue", ""))
                t = seglabel(k)
                isempty(t) || push!(seglabels, ((lst[sg.i0] + lst[sg.i1]) / 2, y, t,
                                                segment_changes(p, k)))
            end
            for (c0, _) in p.changes
                push!(switches, lst[c0])
            end
            push!(rows, (y, "POPs")); y += 1.1
        end
        ymax = max(10.0, y + 0.6)
    else
        show_idx = sort(collect(observable_indices(p)))
        lbl = "Observable"
        push!(rows, (2.0, p.name))
    end

    if !isempty(show_idx) && is_multi_pop(p)
        # One bar per segment, a second shade where a telescope is dropped, and the POP changes
        # between them. The summary view labels each segment above its bar.
        for (k, sg) in enumerate(p.segments)
            idx = intersect(show_idx, sg.i0:sg.i1)
            col = sg.tier == 1 ? (detailed ? "green" : "blue") :
                                 (detailed ? "yellowgreen" : "dodgerblue")
            nb = length(bars)
            add!(idx, 2.0, TARGET_BAR_HEIGHT, col, k == 1 ? lbl : "")
            if touches(k) && length(bars) > nb && sg.i1 in idx
                b = bars[end]
                bars[end] = GanttBar(b.x0, seg_x1(k), b.y, b.height, b.color, b.label)
            end
            if !detailed
                # The first setting where it starts, and a telescope's loss ("−E1") at the point
                # it is lost. The changes are labelled below, where each starts.
                # The first POPs where they start, a telescope not used being a "×"; later, the
                # telescopes lost ("−W2") and regained ("+E1") where it happens.
                if k == 1
                    push!(poplabels, (lst[sg.i0], az_label_y() + POP_LABEL_LIFT, segment_label(p, sg),
                                      segment_changes(p, k), false))
                else
                    t = transition_label(p, k)
                    isempty(t) || push!(seglabels, (lst[sg.i0], az_label_y() + POP_LABEL_LIFT, t, Bool[]))
                end
            end
            # A telescope lost here: "DROP" at the end of the segment before, and the time it is
            # lost at the start of this one, as the other times are.
            if k > 1 && any(i -> sg.config[i] == 0 && p.segments[k-1].config[i] > 0, eachindex(sg.config))
                push!(marks, (seg_x1(k - 1), "DROP", :end, TARGET_ROW))
                push!(labels, GanttLabel(lst[sg.i0], TARGET_ROW,
                                         _hhmm(_gantt_clock(p, lst[sg.i0], time_system)), :time, :start))
            elseif reference_change(p, k) != 0
                # a new reference alone: its time, as a drop's
                push!(labels, GanttLabel(lst[sg.i0], TARGET_ROW,
                                         _hhmm(_gantt_clock(p, lst[sg.i0], time_system)), :time, :start))
            end
        end
        # "POP" where a change starts, and the time it starts just before it.
        for (c0, _) in p.changes
            push!(marks, (lst[c0], "POP", :start, TARGET_ROW))
            push!(marks, (lst[c0], _hhmm(_gantt_clock(p, lst[c0], time_system)), :end, TARGET_ROW))
        end
        if !detailed
            for (c0, c1) in p.changes
                k = findfirst(sg -> sg.i0 >= c1, p.segments)
                k === nothing && continue
                push!(poplabels, (lst[c0], az_label_y() + POP_LABEL_LIFT,
                                  segment_label(p, p.segments[k]),
                                  segment_changes(p, k), true))
            end
        end
    elseif !isempty(show_idx)
        add!(show_idx, 2.0, TARGET_BAR_HEIGHT, detailed ? "green" : "blue", lbl)
    end
    if !isempty(show_idx)
        # Every run is annotated, not just the outermost pair: each is one observing block, and
        # its own start and end times with their az/alt are what goes on a schedule.
        # The az and alt numbers read INWARD from the bar's own edges -- the left pair starts
        # at the left edge, the right pair ends at the right one -- so they stay over the bar
        # they describe. Growing outward, as they did, they met in the gap between two runs and
        # were written over each other there.
        # A POP change is not a gap on the chart, so the runs either side of one are labelled
        # as one, at their outer ends.
        label_idx = isempty(p.changes) ? show_idx :
                    sort(union(show_idx, Int[i for (c0, c1) in p.changes for i in c0:c1-1]))
        for (i0, i1) in index_runs(label_idx)
            # Slivers are drawn but not annotated: six labels inside a few minutes overprint
            # each other and the next run's, which looks like broken text rather than crowding.
            if lst[i1] - lst[i0] >= MIN_LABELLED_RUN
                azy, alty = az_label_y(), alt_label_y()
                # The clock the AXIS is in, so a bar that starts at 22:00 on the axis says
                # 22:00 at its end. Positions stay in LST either way -- only the text changes.
                push!(labels, GanttLabel(lst[i0], TARGET_ROW,
                                         _hhmm(_gantt_clock(p, lst[i0], time_system)), :time, :start))
                push!(labels, GanttLabel(lst[i1], TARGET_ROW,
                                         _hhmm(_gantt_clock(p, lst[i1], time_system)), :time, :end))
                # Both pairs are over the same bar now, so what has to be wide enough is the
                # bar. The times are drawn on their side and cost their height in x, so they
                # are printed whatever the width.
                if lst[i1] - lst[i0] >= LABEL_ROOM_FRAC * span
                    push!(labels, GanttLabel(lst[i0], azy,  string(round(Int, p.az[i0])),  :az,  :start))
                    push!(labels, GanttLabel(lst[i0], alty, string(round(Int, p.alt[i0])), :alt, :start))
                    push!(labels, GanttLabel(lst[i1], azy,  string(round(Int, p.az[i1])),  :az,  :end))
                    push!(labels, GanttLabel(lst[i1], alty, string(round(Int, p.alt[i1])), :alt, :end))
                end
            end
        end
    end

    # The axis IS the night -- sunset to sunrise -- which is what ASPRO draws in "Night only"
    # and what makes a short summer night read as short. The samples already span the dark
    # window, so the range is the grid's own, padded so a bar ending at sunrise is not drawn
    # on the frame. A fixed window centred on midnight would be the alternative and it hides
    # exactly the seasonal difference an observer is looking for.
    # A summary label sits at the az numbers' height, just above its bar, or clear above them
    # when one is under it. Widths are estimated from the character count, a character being
    # taken as 0.8% of the plotted span, on the generous side for the labels' 9 pt.
    if !detailed
        cw = 0.008 * span
        azs = [l for l in labels if l.kind === :az]
        for (k, (x, ylab, t, ch)) in enumerate(seglabels)
            hw = cw * length(t) / 2
            under(l) = l.side === :start ? (l.x - hw < x < l.x + 4cw + hw) :
                                           (l.x - 4cw - hw < x < l.x + hw)
            seglabels[k] = (x, az_label_y() + (any(under, azs) ? POP_LABEL_LIFT : 0.08), t, ch)
        end
    end

    pad = 0.25
    xlim = isempty(lst) ? (mid - 8.0, mid + 8.0) : (first(lst) - pad, last(lst) + pad)

    # Transit is HA = 0, i.e. LST = RA. This axis is unwrapped and can run past 24 h, so bring
    # the target's RA onto the branch the night occupies; NaN when transit is not in the night,
    # which is the honest drawing of "it does not culminate while it is dark".
    ra_h    = p.ra / 15
    tr      = ra_h + 24 * round(((xlim[1] + xlim[2]) / 2 - ra_h) / 24)
    transit = (xlim[1] <= tr <= xlim[2]) ? tr : NaN

    # Bands span whatever height the rows ended up needing.
    bands = [GanttBar(b.x0, b.x1, ymax/2, ymax, b.color, b.label) for b in bands]

    # The POPs only when the delay check actually ran: quoting a setting that did not affect
    # the chart would be worse than quoting none.
    subtitle = !p.delay_applied ? "" : is_multi_pop(p) ? multi_pop_summary(p) : pop_label(p)

    # Time the target spends above the mount's SOFT ceiling: observable, and flagged rather
    # than hidden, which is what ASPRO's dashed extension says. Empty when the facility
    # declares no soft limit, since then there is nothing to distinguish.
    soft_only = setdiff(observable_indices(p), intersect(p.good_alt_soft, p.good_delay,
                                                         p.good_twilight, p.good_moon))
    softbars = GanttBar[]
    if !isempty(soft_only)
        for (i0, i1) in index_runs(sort(collect(soft_only)))
            # The bar's own height, so the dashed edge lands ON the bar's top and bottom:
            # the band is a property of the hours it covers, and drawn shorter it read as a
            # separate, smaller thing floating inside the window rather than as part of it.
            push!(softbars, GanttBar(lst[i0], lst[i1], TARGET_ROW, TARGET_BAR_HEIGHT, "orange", ""))
        end
    end
    # "EL", once inside each stretch above the soft limit, at its centre; not in one under 18
    # minutes, where it would sit on the times at its ends.
    ellabels = [((b.x0 + b.x1) / 2, TARGET_ROW) for b in softbars if b.x1 - b.x0 >= 0.3]
    # The culmination on the target's row (`transit` is its x, NaN when not in the night).
    transits = isfinite(transit) ? [(transit, TARGET_ROW)] : Tuple{Float64,Float64}[]

    return (; bars, labels, bands, softbars, midnight = mid, now = nowlst, transit, transits,
              xlim, target = p.name, rows, ymax, detailed, subtitle,
              delay_applied = p.delay_applied, seglabels, poplabels, ellabels, marks, switches,
              multi = is_multi_pop(p), stack = false)
end

"Rows apart, in the chart's units, of the plans `gantt_stack_geometry` stacks: a bar with its
numbers below and its POPs above."
const STACK_ROW_STEP = 1.6

"""
    gantt_stack_geometry(plans, labels, title = ""; time_system = :lst) -> geometry

Several plans of one target on one night, each a summary row as `gantt_geometry` draws it, the
first at the top: the chart comparing single POP settings (`coverage_pops`), or the candidate
plans of the POP strategies tab. Each row's tick is its entry of `labels`; a single setting's
POPs are written above its bar, in the panel's telescope order, and a multi-POP plan's row
carries its own POPs, changes and drops.
"""
function gantt_stack_geometry(plans::AbstractVector{NightPlan}, labels::AbstractVector{<:AbstractString},
                              title::AbstractString = "";
                              time_system::Symbol = :lst)
    n = length(plans)
    geos = [gantt_geometry(p; time_system) for p in plans]
    g1 = first(geos)
    # The renderer fits its view round the rows, so the chart is only as tall as they are.
    ymax = (n - 1) * STACK_ROW_STEP + 4.0
    sh(y, k) = y + (n - k) * STACK_ROW_STEP
    bars = GanttBar[]; tlabels = GanttLabel[]; softbars = GanttBar[]
    rows = Tuple{Float64,String}[]
    poplabels = Tuple{Float64,Float64,String,Vector{Bool},Bool}[]
    seglabels = Tuple{Float64,Float64,String,Vector{Bool}}[]
    marks = Tuple{Float64,String,Symbol,Float64}[]
    ellabels = Tuple{Float64,Float64}[]; transits = Tuple{Float64,Float64}[]
    for (k, (p, g)) in enumerate(zip(plans, geos))
        append!(bars, [GanttBar(b.x0, b.x1, sh(b.y, k), b.height, b.color, b.label) for b in g.bars])
        append!(softbars, [GanttBar(b.x0, b.x1, sh(b.y, k), b.height, b.color, b.label) for b in g.softbars])
        append!(tlabels, [GanttLabel(l.x, sh(l.y, k), l.text, l.kind, l.side) for l in g.labels])
        append!(ellabels, [(x, sh(y, k)) for (x, y) in g.ellabels])
        append!(transits, [(x, sh(y, k)) for (x, y) in g.transits])
        push!(rows, (sh(TARGET_ROW, k), labels[k]))
        idx = sort(collect(observable_indices(p)))
        if is_multi_pop(p)
            # its own POPs, changes and drops, as its summary row draws them
            append!(poplabels, [(x, sh(y, k), t, c, a) for (x, y, t, c, a) in g.poplabels])
            append!(seglabels, [(x, sh(y, k), t, c) for (x, y, t, c) in g.seglabels])
            append!(marks, [(x, t, sd, sh(y, k)) for (x, t, sd, y) in g.marks])
        elseif !isempty(idx)
            on = [i for i in p.display_order if i <= length(p.config) && p.config[i] > 0]
            x = unwrap_lst(p.lst)[first(idx)]
            push!(poplabels, (x, sh(az_label_y() + POP_LABEL_LIFT, k),
                              join((_pop_text(p.pop, p.config, i) for i in on), " "),
                              falses(length(on)) |> collect, false))
        end
    end
    bands = [GanttBar(b.x0, b.x1, ymax / 2, ymax, b.color, b.label) for b in g1.bands]
    return (; bars, labels = tlabels, bands, softbars, midnight = g1.midnight, now = g1.now,
              transit = g1.transit, transits, xlim = g1.xlim, target = g1.target, rows, ymax,
              detailed = false, subtitle = String(title),
              delay_applied = g1.delay_applied, seglabels,
              poplabels, ellabels, marks,
              switches = Float64[], multi = false, stack = true)
end

"LST hours as `H:M`, matching the annotation format of the matplotlib original."
# An LST position expressed on whichever clock the chart is labelled in.
function _gantt_clock(p::NightPlan, x::Real, system::Symbol)
    system === :lst && return x
    u = lst_to_utc(p, x)
    isfinite(u) || return x
    system === :utc && return u
    return u + local_utc_offset(p, Dates.DateTime(Dates.Date(p.date)) +
                                   Dates.Millisecond(round(Int, u * 3.6e6)))
end

"""
    hour_angle_clock(p::NightPlan, ha; system = :lst) -> Float64

Where hour angle `ha` of `p`'s target falls on the clock `system`, in hours (not wrapped): the
LST `ra/15 + ha` taken on the branch of the night's own LST axis, then as `_gantt_clock` maps it.
"""
function hour_angle_clock(p::NightPlan, ha::Real; system::Symbol = :lst)
    x = p.ra / 15 + Float64(ha)
    lst = unwrap_lst(p.lst)
    isempty(lst) && return x
    mid = (first(lst) + last(lst)) / 2
    x += 24 * round((mid - x) / 24)
    return _gantt_clock(p, x, system)
end

function _hhmm(h::Real)
    # Truncated, not rounded: `gantt_onenight` formats a DateTime, and that drops the seconds
    # rather than rounding them up. Rounding here puts the two renderers a minute apart on any
    # block whose end falls past the half-minute. Reached through milliseconds for the same
    # reason -- that is the resolution `hours_to_date` hands the other renderer, and flooring
    # the raw hours instead loses a minute wherever the float lands just below one (19:05
    # printed as 19:4).
    ms = round(Int, mod(Float64(h), 24.0) * 3_600_000)
    mm = div(ms, 60_000)
    return string(div(mm, 60) % 24, ":", lpad(mod(mm, 60), 2, '0'))
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

    # NaN where a multi-POP plan is changing POPs, which the renderer draws as a break. A
    # branch-free loop: a filtered generator here cost five times the rest of this function.
    function peak(y)
        m = 0.0
        @inbounds @simd for v in y
            m = max(m, ifelse(isfinite(v), abs(v), 0.0))
        end
        return m
    end
    ymax = isempty(curves) ? limit :
           max(limit, maximum(peak(c.y) for c in curves)) * 1.05
    xlim = isempty(lst) ? (0.0, 24.0) : (minimum(lst), maximum(lst))
    switches = [lst[c0] for (c0, _) in p.changes]

    return (; curves, limit, alt = (; x = lst, y = Float64.(p.alt)),
              xlim, ylim = (-ymax, ymax), altlim = (-5.0, 90.0),
              alt_limit = p.alt_limit, alt_max = p.alt_max, target = p.name, switches)
end


# ─────────────────────────────────────────────────────────────────────────────
# What a simulation is OF
# ─────────────────────────────────────────────────────────────────────────────
#
# `simulate` needs a sky, and it takes one of two kinds:
#
#     image      = a FITS path, read by its keywords: 2-D, a wavelength cube, a time series or
#                  4-D (a SPARCO model, star plus image, is a model)
#     flat_model = a compiled model, plus `flat_params`, its current values
#
# Those are genuinely different inputs, not one input with a flag, and the panel has to say
# which it holds before it can say whether Simulate is possible at all. Validating here rather
# than at the call means a wrong file is reported when it is chosen, not after an observation
# has been computed.

"""
    facility_subset(f, config) -> FacilityConfig

The same array with only the selected telescopes in it.

`simulate` observes with every telescope a `FacilityConfig` carries — there is no config vector
in its signature — so simulating a 4T night on a 6T array means handing it a 4T array. Every
per-telescope field is subset together; missing one leaves `simulate` indexing a name vector
shorter than `sta_xyz` several hundred lines later.

`config` is the panel's own encoding, 0 for out and non-zero for in, so it is the vector the
telescope checkboxes already produce.
"""
function facility_subset(f, config::AbstractVector{<:Integer})
    keep = findall(!iszero, config)
    (isempty(keep) || length(keep) == f.ntel) && return f
    length(config) == f.ntel || throw(ArgumentError(
        "the telescope selection has $(length(config)) entries but $(f.name) has $(f.ntel)"))
    length(keep) >= 3 || throw(ArgumentError(
        "closure phases need three telescopes; $(length(keep)) selected"))
    g = deepcopy(f)
    g.ntel      = length(keep)
    g.tel_names = f.tel_names[keep]
    g.sta_names = isempty(f.sta_names) ? f.sta_names : f.sta_names[keep]
    g.tel_diams = isempty(f.tel_diams) ? f.tel_diams : f.tel_diams[keep]
    g.tel_gain  = isempty(f.tel_gain)  ? f.tel_gain  : f.tel_gain[keep]
    # Renumbered rather than carried over: OI_ARRAY's station indices are positions in THIS
    # array, and a subset that kept 1,3,5 would write rows pointing at stations it does not have.
    g.sta_index = collect(1:length(keep))
    g.sta_xyz   = f.sta_xyz[keep, :]
    g.delay_lengths = isempty(f.delay_lengths) ? f.delay_lengths : f.delay_lengths[keep]
    return g
end

"""
    epochs_for_hour_angles(facility, ra_deg, date, ha_min, ha_max, step_minutes) -> Vector{DateTime}

The UTC timestamps at which a target sits at each hour angle on a grid — what `simulate` takes,
from what the panel actually asks for.

The panel asks for hour angles because that is what a uv track is parameterised by; `simulate`
takes UTC because that is what the sky knows. The two are related by the sidereal rate, which
is a constant, so each epoch is solved for INDEPENDENTLY from its own hour angle rather than by
stepping from the last: an epoch grid built by accumulation drifts by the sidereal-vs-solar
difference, 4 minutes per day, which is a whole step over a long track.

Anchored to the nearest occurrence to 00:00 UTC on `date`, so an hour-angle range straddling
transit gives one continuous track around that night rather than two halves a day apart.
"""
function epochs_for_hour_angles(f, ra_deg::Real, date::Date,
                                ha_min::Real, ha_max::Real, step_minutes::Real)
    step_minutes > 0 || throw(ArgumentError("the epoch step must be positive"))
    ha_max >= ha_min || throw(ArgumentError("the hour-angle range runs backwards"))
    base = DateTime(date)                       # 00:00 UTC
    _, h0 = OITOOLS.hour_angle_calc(base, f.lon, Float64(ra_deg) / 15)
    href = first(h0)
    # Sidereal hours per solar hour: the sky runs fast, so an hour of hour angle takes slightly
    # less than an hour of clock.
    SIDEREAL = 1.0027379093
    at(ha) = base + Millisecond(round(Int, ((mod(ha - href + 12, 24) - 12) / SIDEREAL) * 3.6e6))
    step_h = Float64(step_minutes) / 60
    return [at(ha) for ha in Float64(ha_min):step_h:Float64(ha_max)]
end

"""
    simulate_source_info(kind, path) -> (; ok, kind, summary, ndims, nx, nwav, detail)

Describe the sky a simulation would use, and whether it is usable.

`kind` is `"image"` or `"model"` (`"cube"`, from older sessions, is read as `"image"`). For an
image `path` is a FITS file of any dimension, read as `simulate` reads it ([`read_model_image`](@ref)):
the summary gives its size, wavelengths and epochs, and its pixel scale. For `"model"` it is a
TOML model file (see [`read_model_file`](@ref)) — or empty, meaning the model currently held by
the Modeling perspective. A model may hold image components (a SPARCO star plus image).

`ok = false` comes with a `summary` saying why — for an image whose axes the keywords do not
give, the reader's reason and where to set them.
"""
function simulate_source_info(kind::AbstractString, path::AbstractString = "")
    k = lowercase(strip(String(kind)))
    p = String(path)

    no(msg) = (; ok = false, kind = k, summary = msg, ndims = 0, nx = 0, nwav = 0, detail = "",
                 pixsize = 0.0)

    if k == "model"
        isempty(p) && return (; ok = true, kind = k,
                                summary = "the model held by the Modeling perspective",
                                ndims = 0, nx = 0, nwav = 0, detail = "", pixsize = 0.0)
        isfile(p) || return no("no model file at '$p'")
        m = try
            read_model_file(p)
        catch err
            return no("could not read '$(basename(p))': " * sprint(showerror, err))
        end
        ncomp = length(OITOOLS._component_names(m.model))
        return (; ok = true, kind = k,
                  summary = "$(basename(p)): $ncomp component(s), $(length(m.free)) free",
                  ndims = 0, nx = 0, nwav = 0, pixsize = 0.0,
                  detail = join(sort(collect(keys(m.model))), ", "))
    end

    k in ("image", "cube") || return no("unknown source '$kind'; expected image or model")
    isempty(p) && return no("no file chosen")
    isfile(p)  || return no("no file at '$p'")
    # Read as `simulate` will read it: by its keywords. A file with no pixel scale is still
    # described (the panel's pixel size is then the one used); any other refusal is the reason.
    px = fits_pixsize(p)
    mi = try
        read_model_image(p; pixsize = something(px, 1.0))
    catch err
        # the reader's own message already names the file
        why = err isa ArgumentError ? err.msg : sprint(showerror, err)
        return no(why * " — fix it in Tools › FITS Image Handling (Set keywords, Save FITS)")
    end
    nx, ny, nλ, nt = size(mi.image)
    λs = filter(isfinite, mi.wavelengths); ts = filter(isfinite, mi.mjds)
    wave = nλ == 1 ? "1 λ" : "$nλ λ" * (isempty(λs) ? "" :
           @sprintf(" (%.4g–%.4g µm)", minimum(λs) * 1e6, maximum(λs) * 1e6))
    time = nt == 1 ? "1 epoch" : "$nt epochs" * (isempty(ts) ? "" :
           @sprintf(" (MJD %.1f–%.1f)", minimum(ts), maximum(ts)))
    scale = px === nothing ? "no pixel scale in the header" : @sprintf("%.4g mas", px)
    return (; ok = true, kind = "image",
              summary = "$(basename(p)): $(nx)×$(ny) · $wave · $time · $scale",
              ndims = nλ > 1 && nt > 1 ? 4 : nλ > 1 || nt > 1 ? 3 : 2, nx, nwav = nλ,
              detail = "", pixsize = px === nothing ? 0.0 : px)
end
