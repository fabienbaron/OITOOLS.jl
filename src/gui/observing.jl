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
                    alt_limit     ::Real = 30.0,
                    alt_max       ::Real = 90.0,
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
    good_delay, applied = if use_delay
        try
            (Int.(in_delay(f, Float64(dec), obs.ha, cfg, pp; delay_length).good_delay), true)
        catch
            (all_idx, false)
        end
    else
        (all_idx, false)
    end

    return NightPlan(String(name), Float64(ra), Float64(dec), date,
                     obs.lst, obs.lst_midnight, obs.ha, obs.alt, obs.az,
                     Int.(obs.good_alt), good_delay,
                     Int.(obs.good_twilight), Int.(obs.good_moon),
                     obs.moon_sep, obs.moon_fli, applied)
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
    gantt_geometry(plan; show_alt = true) -> (; bars, labels, bands, midnight, xlim, target)

The whole chart as data: background bands, one bar per contiguous run, and the annotations.

`show_alt` is the detailed/summary switch, and it changes what is drawn rather than merely how.
Detailed breaks the constraints onto their own rows, the way ASPRO's "detailed output" does;
summary draws only their intersection — one bar, the hours actually bookable.
"""
function gantt_geometry(p::NightPlan; show_alt::Bool = true)
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

    if show_alt
        isempty(p.good_alt)  || add!(p.good_alt,  5.0, 1.5, "orange",       "Altitude")
        isempty(p.good_moon) || add!(p.good_moon, 7.0, 1.0, "mediumpurple", "Moon sep.")
        show_idx = sort(collect(p.good_delay))
        lbl = "In Delay"
    else
        show_idx = sort(collect(observable_indices(p)))
        lbl = "Observable"
    end

    if !isempty(show_idx)
        add!(show_idx, 2.0, 2.0, "blue", lbl)
        # Every run is annotated, not just the outermost pair: each is one observing block, and
        # its own start and end times with their az/alt are what goes on a schedule.
        for (i0, i1) in index_runs(show_idx)
            push!(labels, GanttLabel(lst[i0], 2.0, _hhmm(lst[i0]), :time, :start))
            push!(labels, GanttLabel(lst[i1], 2.0, _hhmm(lst[i1]), :time, :end))
            push!(labels, GanttLabel(lst[i0], 3.3, string(round(Int, p.az[i0])),  :az,  :start))
            push!(labels, GanttLabel(lst[i0], 0.7, string(round(Int, p.alt[i0])), :alt, :start))
            push!(labels, GanttLabel(lst[i1], 3.3, string(round(Int, p.az[i1])),  :az,  :end))
            push!(labels, GanttLabel(lst[i1], 0.7, string(round(Int, p.alt[i1])), :alt, :end))
        end
    end

    # The window matplotlib uses: midnight ± 8 h, which centres the night regardless of season.
    offset = 4.0
    xlim = (mid - 12 + offset, mid + 12 - offset)

    return (; bars, labels, bands, midnight = mid, xlim, target = p.name)
end

"LST hours as `H:M`, matching the annotation format of the matplotlib original."
function _hhmm(h::Real)
    t = mod(Float64(h), 24.0)
    hh = floor(Int, t)
    mm = round(Int, (t - hh) * 60)
    mm == 60 && (mm = 0; hh += 1)
    return string(mod(hh, 24), ":", mm)
end
