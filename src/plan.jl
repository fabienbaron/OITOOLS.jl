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
    p = zeros(Float32, 7, 5)  # 7 telescopes (S1,S2,E1,E2,W1,W2,S3) × 5 POPs
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
    # S3 (same POP offsets as S1 for now — update when known)
    p[7,:] = [0.0, 36.5639684, 73.1289412, 109.7249843, 143.0358086]
    p
end

const CHARA_AIRPATH = Float32[0.0, 573633.833, 4250500.0, 3680300.0, 1842354.0, 2405131.0, 0.0] .* Float32(1e-6)

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

"""
    datetime_to_jd(dt::DateTime)

Convert a Julia DateTime to Julian Date.
"""
function datetime_to_jd(dt::DateTime)
    Y = Dates.year(dt)
    M = Dates.month(dt)
    D = Dates.day(dt) + (Dates.hour(dt) + Dates.minute(dt)/60.0 + Dates.second(dt)/3600.0) / 24.0
    if M <= 2
        Y -= 1
        M += 12
    end
    A = floor(Int, Y / 100)
    B = 2 - A + floor(Int, A / 4)
    return floor(Int, 365.25 * (Y + 4716)) + floor(Int, 30.6001 * (M + 1)) + D + B - 1524.5
end

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
    lst_midnight, _ = hour_angle_calc(obsdate + Dates.Day(1) + Dates.Hour(7), lon, ra)
    lst_midnight = lst_midnight[1]

    # Dark window
    _, UTC_set  = sunrise_sunset(obsdate, lat, lon)
    UTC_rise, _ = sunrise_sunset(obsdate + Dates.Day(1), lat, lon)
    utc = collect(range(UTC_set + dark_offset, UTC_rise - dark_offset, step=1.0/60*step_minutes))

    # LST and hour angle — compute in Float64 then convert
    dates = hours_to_date(obsdate, utc)
    lst64, ha64 = hour_angle_calc(dates, lon, ra)
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
function obs_plan(targetname::AbstractString, facility::FacilityConfig,
                  ra::Float64, dec::Float64, obsdate::DateTime,
                  pop::Vector{Int}, config::Vector{Int};
                  alt_limit::Float64=30.0, alt_max::Float64=90.0,
                  delay_length::Union{Nothing,Float64}=nothing,
                  dark_offset::Float64=0.0, step_minutes::Int=1,
                  figsize=(10,5), savefile::AbstractString="",
                  show_alt::Bool=true)

    obs = night_observability(facility, ra, dec, obsdate;
                              alt_limit=alt_limit, alt_max=alt_max,
                              dark_offset=dark_offset, step_minutes=step_minutes)

    dr = in_delay(facility, dec, obs.ha, config, pop; delay_length=delay_length)

    gantt_onenight(targetname, obsdate, obs.lst, obs.lst_midnight,
                   obs.az, obs.alt, obs.good_alt, dr.good_delay;
                   good_twilight=obs.good_twilight,
                   figsize=figsize, savefile=savefile, show_alt=show_alt)
end

"""
    empty_night(facility, obsdate; figsize=(10,5), savefile="")

Render a Gantt plot showing only the twilight bands for a given night,
with no target. Useful as a blank canvas before a target is selected.
"""
function empty_night(facility::FacilityConfig, obsdate::DateTime;
                     dark_offset::Float64=0.0, step_minutes::Int=1,
                     figsize=(10,5), savefile::AbstractString="")
    lat, lon = facility.lat, facility.lon
    # Use a dummy RA to get LST grid (RA only shifts HA, not LST)
    ra_dummy = 0.0
    lst_midnight, _ = hour_angle_calc(obsdate + Dates.Day(1) + Dates.Hour(7), lon, ra_dummy)
    lst_midnight = lst_midnight[1]
    _, UTC_set  = sunrise_sunset(obsdate, lat, lon)
    UTC_rise, _ = sunrise_sunset(obsdate + Dates.Day(1), lat, lon)
    utc = collect(range(UTC_set + dark_offset, UTC_rise - dark_offset, step=1.0/60*step_minutes))
    lst, _ = hour_angle_calc(hours_to_date(obsdate, utc), lon, ra_dummy)
    lst = Float32.(lst)

    gantt_onenight("", obsdate, lst, lst_midnight,
                   Float32[], Float32[], Int[], Int[];
                   figsize=figsize, savefile=savefile)
end

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
function chara_plan(targetname::AbstractString, facility::FacilityConfig,
                    ra::Float64, dec::Float64, obsdate::DateTime,
                    pop::Vector{Int}, config::Vector{Int};
                    alt_limit::Float64=30.0, alt_max::Float64=90.0,
                    delay_length::Union{Nothing,Float64}=nothing,
                    dark_offset::Float64=0.0, step_minutes::Int=1)

    obs = night_observability(facility, ra, dec, obsdate;
                              alt_limit=alt_limit, alt_max=alt_max,
                              dark_offset=dark_offset, step_minutes=step_minutes)

    dr = in_delay(facility, dec, obs.ha, config, pop; delay_length=delay_length)

    # Effective delay limit for the plot
    if !isnothing(delay_length)
        dmax = Float64(delay_length)
    elseif !isempty(facility.delay_lengths)
        dmax = maximum(facility.delay_lengths)
    else
        dmax = 45.7
    end

    fig = figure(figsize=(12, 6))
    ax1 = fig.add_subplot(111)

    # Plot delay carts for each baseline
    for b in 1:dr.nbaselines
        ax1.plot(obs.lst, dr.delay_carts[b, :], label=dr.baseline_names[b])
    end

    ax1.axhline(y=dmax, color="gray", linestyle="--", alpha=0.5, label="Delay limit")
    ax1.axhline(y=-dmax, color="gray", linestyle="--", alpha=0.5)

    ax1.set_ylabel("Delay cart position (m)")
    ax1.set_xlabel("LST (hours)")

    # Altitude on secondary axis
    ax2 = ax1.twinx()
    ax2.plot(obs.lst, obs.alt, color="black", linestyle=":", linewidth=1.5, label="Altitude")
    ax2.axhline(y=alt_limit, color="red", linestyle="-", alpha=0.5, label="El. limit $(alt_limit)°")
    ax2.set_ylabel("Altitude (°)")
    ax2.set_ylim(-5, 90)

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(vcat(lines1, lines2), vcat(labels1, labels2),
               loc="upper left", fontsize=8, ncol=2)

    pop_str = join(["$(facility.sta_names[i]):P$(pop[i])" for i in findall(config .> 0)], " ")
    ax1.set_title("$(targetname) — $(Date(obsdate)) — POPs: $(pop_str)")

    ax1.set_xlim(minimum(obs.lst), maximum(obs.lst))
    tight_layout()
    return fig
end

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
