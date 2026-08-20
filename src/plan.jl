# plan.jl — CHARA observation planning utilities
#
# Functions:
#   night_observability  — compute dark-time LST/HA/alt/az grid for a target
#   in_delay             — determine when a target is within delay-line limits
#   best_pop             — brute-force POP configuration optimizer
#   obs_plan             — Gantt-style observability plot (à la ASPRO)
#   chara_plan           — delay vs LST plot with altitude overlay (à la chara_plan)

using Dates, LinearAlgebra

# ─── CHARA POP data (from telescopes.chara, Jun 2019) ────────────────────────

const CHARA_POP_ARRAY = let
    p = zeros(Float32, 6, 5)  # 6 telescopes (S1,S2,E1,E2,W1,W2) × 5 POPs
    # S1
    p[1,:] = [0.0, 36.5639684, 73.1289412, 109.7249843, 143.0358086]
    # S2
    p[2,:] = [-36.5447103, 0.0, 36.5686106, 73.1577372, 106.4675474]
    # E1
    p[3,:] = [0.0, 36.5856149, 73.1206645, 109.7081816, 143.0198375]
    # E2
    p[4,:] = [-73.1122149, -36.5385324, 0.0, 36.578125261, 69.922248720]
    # W1
    p[5,:] = [-73.1098438, -36.5368346, 0.0, 36.5936311, 69.9061629]
    # W2
    p[6,:] = [-143.055000, -106.455324, -69.911098, -33.306356, 0.0]
    #
    # There is deliberately no seventh row. This table used to carry one for S3, filled with
    # a copy of S1's offsets and a "update when known" note -- but S3 has no POPs, so there
    # are no numbers to fill in and a POP search including it was silently optimising against
    # S1. CHARA.toml lists the six telescopes above and nothing else, and both this table and
    # CHARA_AIRPATH are indexed by telescope index into that list, so a seventh row was
    # unreachable as well as wrong.
    p
end

# One entry per telescope, in CHARA.toml order (S1, S2, E1, E2, W1, W2). See the note in
# CHARA_POP_ARRAY for why there is no seventh.
const CHARA_AIRPATH = Float32[0.0, 573633.833, 4250500.0, 3680300.0, 1842354.0, 2405131.0] .* Float32(1e-6)

# ─── Low-precision Moon ephemeris (Meeus, Astronomical Algorithms) ────────────

"""
    moon_radec(jd)

Low-precision Moon RA/Dec (degrees) for a given Julian Date.
Accuracy ~1° in position — sufficient for separation checks.
"""
function moon_radec(jd::Float64)
    T = (jd - 2451545.0) / 36525.0  # centuries from J2000.0

    # Fundamental arguments (degrees)
    Lp = mod(218.3165 + 481267.8813 * T, 360.0)  # mean longitude
    D  = mod(297.8502 + 445267.1115 * T, 360.0)  # mean elongation
    M  = mod(357.5291 +  35999.0503 * T, 360.0)  # Sun mean anomaly
    Mp = mod(134.9634 + 477198.8676 * T, 360.0)  # Moon mean anomaly
    F  = mod( 93.2721 + 483202.0175 * T, 360.0)  # argument of latitude

    d2r = π / 180.0

    # Ecliptic longitude (λ) and latitude (β) — main terms
    λ = Lp +
        6.289 * sin(Mp * d2r) +
        1.274 * sin((2D - Mp) * d2r) +
        0.658 * sin(2D * d2r) +
        0.214 * sin(2Mp * d2r) -
        0.186 * sin(M * d2r) -
        0.114 * sin(2F * d2r)

    β = 5.128 * sin(F * d2r) +
        0.281 * sin((Mp + F) * d2r) +
        0.278 * sin((Mp - F) * d2r) +
        0.173 * sin((2D - F) * d2r)

    # Obliquity of ecliptic
    ε = 23.4393 - 0.0130 * T

    # Ecliptic to equatorial
    λr = λ * d2r
    βr = β * d2r
    εr = ε * d2r

    ra  = atan(sin(λr) * cos(εr) - tan(βr) * sin(εr), cos(λr))
    dec = asin(sin(βr) * cos(εr) + cos(βr) * sin(εr) * sin(λr))

    return mod(ra * 180.0 / π, 360.0), dec * 180.0 / π
end

"""
    moon_illumination(jd)

Fractional lunar illumination (0–1) for a given Julian Date.
"""
function moon_illumination(jd::Float64)
    T = (jd - 2451545.0) / 36525.0
    D = mod(297.8502 + 445267.1115 * T, 360.0)
    M = mod(357.5291 +  35999.0503 * T, 360.0)
    Mp = mod(134.9634 + 477198.8676 * T, 360.0)
    d2r = π / 180.0
    # Phase angle (simplified)
    i = 180.0 - D - 6.289 * sin(Mp * d2r) +
        2.100 * sin(M * d2r) -
        1.274 * sin((2D - Mp) * d2r) -
        0.658 * sin(2D * d2r) -
        0.214 * sin(2Mp * d2r)
    return (1.0 + cos(i * d2r)) / 2.0
end

"""
    angular_separation(ra1, dec1, ra2, dec2)

Angular separation in degrees between two positions (all in degrees).
"""
function angular_separation(ra1::Float64, dec1::Float64, ra2::Float64, dec2::Float64)
    d2r = π / 180.0
    Δra = (ra2 - ra1) * d2r
    δ1 = dec1 * d2r
    δ2 = dec2 * d2r
    cos_sep = sin(δ1) * sin(δ2) + cos(δ1) * cos(δ2) * cos(Δra)
    return acos(clamp(cos_sep, -1.0, 1.0)) * 180.0 / π
end

# datetime_to_jd lives in astrometry.jl (included before this file) so that simulate.jl
# can use it too.

# ─── Night observability ──────────────────────────────────────────────────────

"""
    night_observability(facility, ra, dec, obsdate; alt_limit=30.0, alt_max=90.0,
                        moon_min_sep=30.0, dark_offset=0.0, step_minutes=1)

Compute the observability grid for a target on a given night.
Internally uses Float32 for LST, HA, alt, az arrays.

Returns a NamedTuple with fields:
- `utc`: UTC hours vector (Float64)
- `lst`: local sidereal time, hours (Float32)
- `ha`:  hour angle, hours (Float32)
- `alt`: altitude, degrees (Float32)
- `az`:  azimuth, degrees (Float32)
- `lst_midnight`: LST at local midnight (Float64)
- `good_alt`: indices where `alt_limit < alt < alt_max`
- `moon_sep`: Moon–target separation, degrees (Float32)
- `moon_fli`: fractional lunar illumination (Float64, single value for the night)
- `good_moon`: indices where `moon_sep > moon_min_sep`

Arguments:
- `facility`: FacilityConfig (provides lat, lon)
- `ra`: right ascension in degrees
- `dec`: declination in degrees
- `obsdate`: DateTime of the observing evening
- `alt_limit`: minimum elevation in degrees (default 30)
- `alt_max`: maximum safe elevation in degrees (default 90, e.g. 80 for CHARA)
- `moon_min_sep`: minimum Moon separation in degrees (default 30)
- `dark_offset`: hours after sunset / before sunrise for deeper twilight (default 0)
- `step_minutes`: time resolution in minutes (default 1)
"""
function night_observability(facility::FacilityConfig, ra::Float64, dec::Float64,
                             obsdate::DateTime;
                             alt_limit::Float64=30.0, alt_max::Float64=90.0,
                             moon_min_sep::Float64=30.0,
                             dark_offset::Float64=0.0,
                             step_minutes::Int=1)
    lat, lon = facility.lat, facility.lon

    # LST at local midnight (next day, 7h UT ≈ midnight at CHARA longitude)
    # `ra` is in degrees (see the docstring and TargetConfig.raep0); hour_angle_calc wants
    # hours. Before 0.11 `ra` was passed through unconverted, so callers had to supply hours
    # -- which then made the angular_separation() call below, which needs degrees, wrong.
    lst_midnight, _ = hour_angle_calc(obsdate + Dates.Day(1) + Dates.Hour(7), lon, ra/15)
    lst_midnight = lst_midnight[1]

    # Dark window
    _, UTC_set  = sunrise_sunset(obsdate, lat, lon)
    UTC_rise, _ = sunrise_sunset(obsdate + Dates.Day(1), lat, lon)
    utc = collect(range(UTC_set + dark_offset, UTC_rise - dark_offset, step=1.0/60*step_minutes))

    # LST and hour angle — compute in Float64 then convert
    dates = hours_to_date(obsdate, utc)
    lst64, ha64 = hour_angle_calc(dates, lon, ra/15)
    lst = Float32.(lst64)
    ha  = Float32.(ha64)

    # Altitude and azimuth
    alt64, az64 = alt_az(dec, lat, ha64)
    alt_vec = Float32.(alt64)
    az_vec  = Float32.(az64)

    good_alt = findall((alt_vec .> alt_limit) .& (alt_vec .< alt_max))

    # Moon separation and illumination
    jd_mid = datetime_to_jd(dates[length(dates) ÷ 2 + 1])
    moon_fli = moon_illumination(jd_mid)

    moon_sep = Vector{Float32}(undef, length(utc))
    for i in eachindex(utc)
        jd = datetime_to_jd(dates[i])
        mra, mdec = moon_radec(jd)
        moon_sep[i] = Float32(angular_separation(ra, dec, mra, mdec))
    end

    good_moon = findall(moon_sep .> moon_min_sep)

    # Twilight filter: indices within 1.5h of sunset/sunrise (outermost light-gray band)
    # Use unwrapped LST to handle 24→0 crossing
    lst_unwrap = Float64.(lst)
    for i in 2:length(lst_unwrap)
        while lst_unwrap[i] < lst_unwrap[i-1]
            lst_unwrap[i] += 24.0
        end
    end
    tw_start = lst_unwrap[1] + 1.5
    tw_end   = lst_unwrap[end] - 1.5
    good_twilight = findall((lst_unwrap .>= tw_start) .& (lst_unwrap .<= tw_end))

    return (utc=utc, lst=lst, ha=ha, alt=alt_vec, az=az_vec,
            lst_midnight=lst_midnight, good_alt=good_alt,
            moon_sep=moon_sep, moon_fli=moon_fli, good_moon=good_moon,
            good_twilight=good_twilight)
end

# ─── Delay computation ────────────────────────────────────────────────────────

"""
    compute_delays(facility, dec, ha, config, pop; pop_array=CHARA_POP_ARRAY, airpath=CHARA_AIRPATH)

Compute delay cart positions for all baselines over a time grid.

Returns `(delay_carts, nbaselines, baseline_names, baseline_stations)` where
`delay_carts` is `nbaselines × length(ha)`.

Arguments:
- `facility`: FacilityConfig
- `dec`: declination in degrees
- `ha`: hour angle vector (hours)
- `config`: telescope configuration vector (0=unused, 1=use, 2=reference)
- `pop`: POP assignment vector (one per telescope)
"""
function compute_delays(facility::FacilityConfig, dec::Float64, ha::Vector{Float32},
                        config::Vector{Int}, pop::Vector{Int};
                        pop_array::Matrix{Float32}=CHARA_POP_ARRAY,
                        airpath::Vector{Float32}=CHARA_AIRPATH)
    nbaselines, baseline_xyz, baseline_stations, baseline_names = get_baselines(facility; config=config)
    delay_geo = Float32.(geometric_delay(facility.lat, Float64.(ha), dec, baseline_xyz))

    delay_airpath = Float32[airpath[baseline_stations[2, i]] - airpath[baseline_stations[1, i]]
                            for i in 1:nbaselines]
    delay_pop = Float32[pop_array[baseline_stations[2, i], pop[baseline_stations[2, i]]] -
                        pop_array[baseline_stations[1, i], pop[baseline_stations[1, i]]]
                        for i in 1:nbaselines]

    delay_carts = 0.5f0 .* (delay_geo .- delay_airpath .- delay_pop)
    return delay_carts, nbaselines, baseline_names, baseline_stations
end

"""
    in_delay(facility, dec, ha, config, pop; delay_length=nothing, kwargs...)

Determine when the target is within delay-line limits for all baselines.

Returns a NamedTuple with:
- `delay_carts`: delay cart positions (nbaselines × ntimes)
- `has_delay`: BitVector, true where all baselines are within limits
- `good_delay`: indices into the time grid where observing is feasible
- `nbaselines`, `baseline_names`, `baseline_stations`

The per-telescope delay limits come from `facility.delay_lengths`.
Pass `delay_length=43.0` to override with a uniform (conservative) value.
"""
function in_delay(facility::FacilityConfig, dec::Float64, ha::Vector{Float32},
                  config::Vector{Int}, pop::Vector{Int};
                  delay_length::Union{Nothing,Float64}=nothing,
                  pop_array::Matrix{Float32}=CHARA_POP_ARRAY,
                  airpath::Vector{Float32}=CHARA_AIRPATH)

    delay_carts, nbaselines, baseline_names, baseline_stations =
        compute_delays(facility, dec, ha, config, pop; pop_array=pop_array, airpath=airpath)

    # Per-baseline delay limit = min of the two telescopes' delay lengths
    if !isnothing(delay_length)
        dlens = fill(Float32(delay_length), facility.ntel)
    elseif !isempty(facility.delay_lengths)
        dlens = Float32.(facility.delay_lengths)
    else
        dlens = fill(45.7f0, facility.ntel)
    end

    has_delay = trues(size(delay_carts, 2))
    for b in 1:nbaselines
        dmax = min(dlens[baseline_stations[1, b]], dlens[baseline_stations[2, b]])
        for t in eachindex(has_delay)
            if delay_carts[b, t] < -dmax || delay_carts[b, t] > dmax
                has_delay[t] = false
            end
        end
    end

    good_delay = findall(has_delay)
    return (delay_carts=delay_carts, has_delay=has_delay, good_delay=good_delay,
            nbaselines=nbaselines, baseline_names=baseline_names,
            baseline_stations=baseline_stations)
end

# ─── Observability filter for simulate() ─────────────────────────────────────

"""
    observable_epochs(facility, target, dates; min_elevation=nothing, max_elevation=nothing,
                      pops=nothing, config=Int[], delay_length=nothing)

Select the epochs in `dates` at which `target` is actually observable from `facility`.

This is **opt-in and composable**: `simulate()` does not call it unless you ask. A simulation
made to test image reconstruction usually wants the full uv coverage regardless of whether a
real night could deliver it, so every constraint here is off by default and only switches on
when you pass the corresponding keyword.

- `min_elevation`, `max_elevation`: degrees. `nothing` (default) applies no elevation cut.
- `pops`: a POP configuration, one entry per telescope, values in `1:5`. `nothing` (default)
  applies **no delay-line check at all**. POPs are never chosen for you — run [`best_pop`](@ref)
  first if you want a recommendation, then pass the result here.
- `config`: telescope-use flags for the delay check (`1` use, `0` skip, `2` reference cart);
  defaults to using every telescope, i.e. all `ntel*(ntel-1)/2` baselines must fit.
- `delay_length`: override the per-telescope delay-line lengths with one uniform value (m).

`target.raep0` / `target.decep0` are in degrees, as elsewhere in the OIFITS-facing API.

Returns `(dates, mask, report)`:
- `dates`: the surviving epochs, in input order
- `mask::BitVector`: true where the epoch survived
- `report`: `(n_in, n_out, n_dropped_elevation, n_dropped_delay, elevation)`, where
  `elevation` is the altitude in degrees at every *input* epoch.

```julia
dates_ok, mask, rep = observable_epochs(facility, target, dates;
                                        min_elevation = 30.0,
                                        pops = [1,3,5,2,4,1])
simulate(facility, target, combiner, wavelength, dates_ok, "sim.oifits"; flat_model=m)
```
"""
function observable_epochs(facility::FacilityConfig, target::TargetConfig,
                           dates::AbstractVector{DateTime};
                           min_elevation::Union{Nothing,Real}=nothing,
                           max_elevation::Union{Nothing,Real}=nothing,
                           pops::Union{Nothing,AbstractVector{<:Integer}}=nothing,
                           config::AbstractVector{<:Integer}=Int[],
                           delay_length::Union{Nothing,Float64}=nothing)

    isempty(dates) && return (dates=dates, mask=BitVector(), report=(n_in=0, n_out=0,
                              n_dropped_elevation=0, n_dropped_delay=0, elevation=Float64[]))

    _, ha_hours = hour_angle_calc(collect(dates), facility.lon, target.raep0/15)
    alt, _ = alt_az(target.decep0, facility.lat, ha_hours)

    n = length(dates)
    keep = trues(n)

    n_drop_elev = 0
    if !isnothing(min_elevation) || !isnothing(max_elevation)
        lo = isnothing(min_elevation) ? -Inf : Float64(min_elevation)
        hi = isnothing(max_elevation) ?  Inf : Float64(max_elevation)
        for i in 1:n
            if alt[i] < lo || alt[i] > hi
                keep[i] = false
                n_drop_elev += 1
            end
        end
    end

    n_drop_delay = 0
    if !isnothing(pops)
        length(pops) == facility.ntel || throw(ArgumentError(
            "pops has $(length(pops)) entries but facility has $(facility.ntel) telescopes"))
        cfg = isempty(config) ? ones(Int, facility.ntel) : collect(Int, config)
        dl = in_delay(facility, target.decep0, Float32.(ha_hours), cfg, collect(Int, pops);
                      delay_length=delay_length)
        for i in 1:n
            if !dl.has_delay[i]
                keep[i] && (n_drop_delay += 1)
                keep[i] = false
            end
        end
    end

    report = (n_in=n, n_out=count(keep), n_dropped_elevation=n_drop_elev,
              n_dropped_delay=n_drop_delay, elevation=collect(alt))
    return (dates=dates[keep], mask=keep, report=report)
end

# ─── Best POP search ─────────────────────────────────────────────────────────

"""
    best_pop(facility, dec, ha, config; n_best=5, min_minutes=10, delay_length=nothing)

Brute-force search over all POP combinations.

Returns a vector of NamedTuples `(pop, score)` sorted by score (minutes observable),
keeping up to `n_best` results with score ≥ `min_minutes`.

Arguments:
- `config`: telescope configuration (0/1/2); only telescopes with config>0 are used
- `n_best`: number of top solutions to return
- `min_minutes`: discard solutions below this threshold
- `delay_length`: override per-telescope delay limits with a uniform value (m)
"""
function best_pop(facility::FacilityConfig, dec::Float64, ha::Vector{Float32},
                  config::Vector{Int};
                  n_best::Int=5, min_minutes::Int=10,
                  delay_length::Union{Nothing,Float64}=nothing,
                  pop_array::Matrix{Float32}=CHARA_POP_ARRAY,
                  airpath::Vector{Float32}=CHARA_AIRPATH)

    nbaselines, baseline_xyz, baseline_stations, baseline_names = get_baselines(facility; config=config)
    delay_geo = Float32.(geometric_delay(facility.lat, Float64.(ha), dec, baseline_xyz))

    delay_airpath = Float32[airpath[baseline_stations[2, i]] - airpath[baseline_stations[1, i]]
                            for i in 1:nbaselines]

    if !isnothing(delay_length)
        dlens = fill(Float32(delay_length), facility.ntel)
    elseif !isempty(facility.delay_lengths)
        dlens = Float32.(facility.delay_lengths)
    else
        dlens = fill(45.7f0, facility.ntel)
    end

    # Active telescope indices
    active = findall(config .> 0)
    npops = 5
    nactive = length(active)

    # Build all POP combos for active telescopes
    results = Tuple{Vector{Int}, Int}[]

    # Iterate over all POP combinations
    _pop_search!(results, facility, active, nbaselines, baseline_stations,
                 delay_geo, delay_airpath, dlens, pop_array, npops, nactive, min_minutes)

    sort!(results, by=x -> x[2], rev=true)

    # Build full pop vectors and return
    out = NamedTuple{(:pop, :score), Tuple{Vector{Int}, Int}}[]
    for (i, (pop_active, score)) in enumerate(results)
        i > n_best && break
        pop_full = ones(Int, facility.ntel)
        for (k, tel_idx) in enumerate(active)
            pop_full[tel_idx] = pop_active[k]
        end
        push!(out, (pop=pop_full, score=score))
    end
    return out
end

function _pop_search!(results, facility, active, nbaselines, baseline_stations,
                      delay_geo, delay_airpath, dlens, pop_array, npops, nactive, min_minutes)
    ntimes = size(delay_geo, 2)
    pop_active = ones(Int, nactive)

    # Full pop vector (all telescopes)
    pop_full = ones(Int, facility.ntel)

    function recurse(depth)
        if depth > nactive
            # Evaluate this POP combination
            for (k, tel_idx) in enumerate(active)
                pop_full[tel_idx] = pop_active[k]
            end

            score = 0
            for t in 1:ntimes
                all_ok = true
                for b in 1:nbaselines
                    i1 = baseline_stations[1, b]
                    i2 = baseline_stations[2, b]
                    dp = pop_array[i2, pop_full[i2]] - pop_array[i1, pop_full[i1]]
                    da = delay_airpath[b]
                    dc = 0.5f0 * (delay_geo[b, t] - da - dp)
                    dmax = min(dlens[i1], dlens[i2])
                    if dc < -dmax || dc > dmax
                        all_ok = false
                        break
                    end
                end
                score += all_ok
            end

            if score >= min_minutes
                push!(results, (copy(pop_active), score))
            end
            return
        end

        for p in 1:npops
            pop_active[depth] = p
            recurse(depth + 1)
        end
    end

    recurse(1)
end

# ─── Gantt-style observing plan plot ──────────────────────────────────────────

"""
    obs_plan(targetname, facility, ra, dec, obsdate, pop, config;
             alt_limit=30.0, alt_max=90.0, delay_length=nothing,
             dark_offset=0.0, step_minutes=1, figsize=(10,5), savefile="")

Produce a Gantt-style observability plot for one target on one night,
showing dark time, altitude window, and delay feasibility.

Pass `delay_length=43.0` for a conservative delay estimate.
Pass `savefile="path.png"` to save to file and close the figure.
"""

"""
    empty_night(facility, obsdate; figsize=(10,5), savefile="")

Render a Gantt plot showing only the twilight bands for a given night,
with no target. Useful as a blank canvas before a target is selected.
"""

# ─── CHARA-plan style delay plot ──────────────────────────────────────────────

"""
    chara_plan(targetname, facility, ra, dec, obsdate, pop, config;
               alt_limit=30.0, alt_max=90.0, delay_length=nothing,
               dark_offset=0.0, step_minutes=1)

Produce a delay-vs-LST plot with altitude overlay, similar to the
classic chara_plan software from GSU.

Each baseline's delay cart position is plotted vs LST, with the altitude
curve and elevation limit overlaid.

Pass `delay_length=43.0` for a conservative delay estimate.
"""

# ─── Display helpers ──────────────────────────────────────────────────────────

"""
    print_pop_results(facility, config, results)

Pretty-print the output of `best_pop`.
"""
function print_pop_results(facility::FacilityConfig, config::AbstractVector{Int}, results)
    active = findall(config .> 0)
    println("─── Best POP configurations ───")
    for (rank, r) in enumerate(results)
        pop_str = join(["$(facility.sta_names[i])-POP$(r.pop[i])" for i in active], "  ")
        println("  #$rank  ($( r.score) min):  $pop_str")
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Contiguous runs of an index vector
# ─────────────────────────────────────────────────────────────────────────────
#
# In the core rather than beside one of the renderers, because BOTH draw from it: a Gantt bar
# is one contiguous run, and a matplotlib chart and a Makie one that each found their own runs
# would be free to disagree about where an observing block starts.

"""
    index_runs(idx) -> Vector{Tuple{Int,Int}}

Split a sorted index vector into maximal runs of consecutive integers.

An observability window is not necessarily one block: a target can dip below the elevation
limit and come back, and a baseline routinely leaves and re-enters the delay range, which is
the normal shape of a POP configuration. Drawing such a set as a single bar from `idx[1]` to
`idx[end]` -- which is what this file used to do -- paints the gaps as observable.
"""
function index_runs(idx::AbstractVector{<:Integer})
    runs = Tuple{Int,Int}[]
    isempty(idx) && return runs
    v = sort(collect(idx))
    first_i = v[1]; prev = v[1]
    for k in @view v[2:end]
        if k == prev + 1
            prev = k
        else
            push!(runs, (first_i, prev)); first_i = k; prev = k
        end
    end
    push!(runs, (first_i, prev))
    return runs
end
