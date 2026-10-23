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
- `horizon`: the lowest usable elevation at each sample, degrees: the larger of `alt_limit`
  and the horizon profile of the telescopes `config` selects
- `good_alt`: indices where `horizon < alt < alt_max`
- `good_riseset`: indices where `alt_limit < alt < alt_max` (ASPRO's Rise/Set)
- `good_horizon`: indices where the target is above the horizon profile (ASPRO's Horizon); every
  index when the facility has no profile
- `has_horizon`: whether the facility has horizon profiles
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
                             config::AbstractVector{<:Integer}=Int[],
                             moon_min_sep::Float64=30.0,
                             dark_offset::Float64=0.0,
                             step_minutes::Int=1)
    lat, lon = facility.lat, facility.lon

    # LST at the local midnight closing `obsdate`'s evening -- the same night the grid below
    # covers. From the longitude rather than a fixed hour, which was CHARA's offset and wrong
    # anywhere else.
    #
    # `ra` is in degrees (see the docstring and TargetConfig.raep0); hour_angle_calc wants
    # hours. Before 0.11 `ra` was passed through unconverted, so callers had to supply hours
    # -- which then made the angular_separation() call below, which needs degrees, wrong.
    midnight_utc = Dates.DateTime(Dates.Date(obsdate)) +
                   Dates.Millisecond(round(Int, (24.0 - lon / 15) * 3.6e6))
    lst_midnight, _ = hour_angle_calc(midnight_utc, lon, ra/15)
    lst_midnight = lst_midnight[1]

    # Dark window of the night BEGINNING on `obsdate`, which is the night an observer means
    # when they name a date. Hours are measured from `obsdate` at 00:00 UT and run past 24 for
    # the morning half, which `hours_to_date` rolls over.
    t_set, t_rise = night_window(obsdate, lat, lon)
    base  = Dates.DateTime(Dates.Date(obsdate))
    hrs(t) = Dates.value(t - base) / 3.6e6
    utc = collect(range(hrs(t_set) + dark_offset, hrs(t_rise) - dark_offset,
                        step = 1.0/60*step_minutes))

    # LST and hour angle — compute in Float64 then convert
    dates = hours_to_date(obsdate, utc)
    lst64, ha64 = hour_angle_calc(dates, lon, ra/15)
    lst = Float32.(lst64)
    ha  = Float32.(ha64)

    # Altitude and azimuth
    alt64, az64 = alt_az(dec, lat, ha64)
    alt_vec = Float32.(alt64)
    az_vec  = Float32.(az64)

    # The horizon the telescopes actually look over, where the facility gives one, floored by
    # whatever the caller asked for: `alt_limit` is a preference and the terrain is not.
    profile = horizon_limit(facility, az_vec; config = config)
    horizon = max.(Float64(alt_limit), profile)
    good_alt = findall(i -> alt_vec[i] > horizon[i] && alt_vec[i] < alt_max, eachindex(alt_vec))
    good_riseset = findall(i -> alt_limit < alt_vec[i] < alt_max, eachindex(alt_vec))
    good_horizon = findall(i -> alt_vec[i] > profile[i], eachindex(alt_vec))
    # The same window under the SOFT ceiling, when the facility declares one. Time between the
    # two is observable and flagged, not lost, which is how ASPRO draws it.
    soft = isfinite(facility.alt_soft_max) ? min(alt_max, facility.alt_soft_max) : alt_max
    good_alt_soft = findall(i -> alt_vec[i] > horizon[i] && alt_vec[i] < soft, eachindex(alt_vec))

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

    # Astronomical night: the Sun 18° below the horizon, which is what "dark" means for an
    # observing window. A fixed pad off the ends of the grid is not that — it is unrelated to
    # the Sun, varies with season, latitude and the grid's own extent, and at CHARA in April
    # discarded nearly two hours of genuinely dark sky every night.
    a_set, a_rise = night_window(obsdate, lat, lon; zenith = 108.0)
    good_twilight = findall(d -> a_set <= d <= a_rise, dates)

    return (utc=utc, lst=lst, ha=ha, alt=alt_vec, az=az_vec, horizon=horizon,
            lst_midnight=lst_midnight, good_alt=good_alt, good_alt_soft=good_alt_soft,
            good_riseset=good_riseset, good_horizon=good_horizon,
            has_horizon=!isempty(facility.horizon_az),
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
                        pop_array::Union{Nothing,Matrix{Float32}}=nothing,
                        airpath::Union{Nothing,Vector{Float32}}=nothing,
                        beam_order::Union{Nothing,AbstractVector{<:Integer}}=nothing)
    nbaselines, baseline_xyz, baseline_stations, baseline_names = get_baselines(facility; config=config)
    delay_geo = Float32.(geometric_delay(facility.lat, Float64.(ha), dec, baseline_xyz))
    off = _delay_offsets(facility, config, pop; pop_array, airpath, beam_order)

    delay_carts = 0.5f0 .* (delay_geo .- Float32.(off))
    return delay_carts, nbaselines, baseline_names, baseline_stations
end

# The fixed optical path difference of each baseline: the air path ahead of the delay line,
# the beam sampling table for the channel that telescope feeds, and the POPs. Factored out
# because the sampled path and the analytic one must not be able to disagree about it.
#
# The facility's own tables win -- these belong in the configuration, beside the station
# coordinates they are measured against. The CHARA_* constants remain for a facility file that
# predates them, and an explicit argument still overrides both.
function _delay_offsets(facility::FacilityConfig, config::Vector{Int}, pop::Vector{Int};
                        pop_array=nothing, airpath=nothing, beam_order=nothing)
    pops = _pop_table(facility; pop_array)
    path = _fixed_path(facility, config; airpath, beam_order)
    nbl, _, sts, _ = get_baselines(facility; config = config)
    return Float64[(path[sts[2, i]] - path[sts[1, i]]) +
                   (pops[sts[2, i], pop[sts[2, i]]] - pops[sts[1, i], pop[sts[1, i]]])
                   for i in 1:nbl]
end

# The per-telescope halves of the offset above, which `best_pop` needs separately: it varies
# the POPs over a fixed path rather than evaluating one assignment.
_pop_table(facility::FacilityConfig; pop_array=nothing) =
    pop_array !== nothing ? pop_array :
    isempty(facility.pop_offsets) ? CHARA_POP_ARRAY : Float32.(facility.pop_offsets)

# A telescope's channel is its POSITION among the telescopes in use, so the same six in a
# different order feed different channels and get different paths -- which is why the BST is
# added per configuration rather than folded into `fixed_offset`. `beam_order` names that
# order; without one it is the facility's own, and no BST means no change.
function _fixed_path(facility::FacilityConfig, config::Vector{Int};
                     airpath=nothing, beam_order=nothing)
    path = airpath !== nothing ? airpath :
           (isempty(facility.fixed_offsets) || any(isnan, facility.fixed_offsets)) ?
           CHARA_AIRPATH : Float32.(facility.fixed_offsets)
    isempty(facility.bst) && return path
    used = beam_order === nothing ? findall(!iszero, config) : collect(beam_order)
    path = copy(path)
    for (chan, t) in enumerate(used)
        chan <= size(facility.bst, 2) || break
        path[t] += Float32(facility.bst[t, chan])
    end
    return path
end

# Per-telescope cart span, in metres of CART travel. The limit switches are the measurement;
# `delay_lengths` is their difference, for a facility file that gives no switches.
function _delay_spans(facility::FacilityConfig, delay_length::Union{Nothing,Float64})
    delay_length !== nothing && return fill(Float64(delay_length), facility.ntel)
    if length(facility.delay_front) == facility.ntel &&
       !any(isnan, facility.delay_front) && !any(isnan, facility.delay_back)
        return Float64.(facility.delay_back .- facility.delay_front)
    end
    isempty(facility.delay_lengths) || return Float64.(facility.delay_lengths)
    return fill(45.7, facility.ntel)
end

# w(h) = A cos h + B sin h + C, the delay of one baseline as a function of hour angle.
# `l` and `dec` in radians; the baseline vector is local (East, North, Up).
@inline function _w_coeffs(l::Float64, δ::Float64, bx::Float64, by::Float64, bz::Float64)
    cd, sd = cos(δ), sin(δ)
    return (cd * (cos(l) * bz - sin(l) * bx), -cd * by, sd * (cos(l) * bx + sin(l) * bz))
end

@inline _w_at(A, B, C, h) = A * cos(h * π / 12) + B * sin(h * π / 12) + C

# Hour angles, in hours, where `wlo <= w(h) <= whi`, restricted to `ha_range`. Exact: the
# crossings are the roots of R cos(h - φ) = W - C, and between consecutive roots w cannot
# cross a limit, so each sub-interval is settled by its midpoint.
function _w_intervals(A::Float64, B::Float64, C::Float64, wlo::Float64, whi::Float64,
                      ha_range::Tuple{Float64,Float64})
    R = hypot(A, B)
    (whi < C - R || wlo > C + R) && return Tuple{Float64,Float64}[]   # never within limits
    (wlo <= C - R && whi >= C + R) && return [ha_range]               # never outside them
    φ = atan(B, A)
    cuts = Float64[ha_range[1], ha_range[2]]
    for W in (wlo, whi)
        R == 0 && break
        k = (W - C) / R
        abs(k) > 1 && continue
        a = acos(clamp(k, -1.0, 1.0))
        for h in (φ + a, φ - a)
            x = (mod(h + π, 2π) - π) * 12 / π
            ha_range[1] < x < ha_range[2] && push!(cuts, x)
        end
    end
    sort!(cuts)
    out = Tuple{Float64,Float64}[]
    for i in 1:length(cuts)-1
        lo, hi = cuts[i], cuts[i+1]
        hi - lo < 1e-12 && continue
        if wlo <= _w_at(A, B, C, 0.5 * (lo + hi)) <= whi
            if !isempty(out) && lo - out[end][2] < 1e-12
                out[end] = (out[end][1], hi)          # the cut was not a crossing
            else
                push!(out, (lo, hi))
            end
        end
    end
    return out
end

# Intersection of two sorted, disjoint interval lists, and the total length of one.
function _isect(a::Vector{Tuple{Float64,Float64}}, b::Vector{Tuple{Float64,Float64}})
    out = Tuple{Float64,Float64}[]
    i = j = 1
    while i <= length(a) && j <= length(b)
        lo = max(a[i][1], b[j][1]); hi = min(a[i][2], b[j][2])
        hi > lo && push!(out, (lo, hi))
        a[i][2] < b[j][2] ? (i += 1) : (j += 1)
    end
    return out
end
_ilen(v::Vector{Tuple{Float64,Float64}}) = isempty(v) ? 0.0 : sum(x -> x[2] - x[1], v)

"""
    delay_ha_intervals(facility, dec, config, pop; kwargs...) -> Vector{Tuple{Float64,Float64}}

Hour-angle intervals, in hours, where every baseline is within its delay limits, with the edges
solved analytically. The delay of a baseline is a sinusoid in hour angle,

    w(h) = A cos h + B sin h + C,    A = cosδ (cos l b₃ − sin l b₁)
                                     B = −cosδ b₂
                                     C = sinδ (cos l b₁ + sin l b₃)

so `w(h) = W` inverts in closed form: with `R = √(A²+B²)` and `φ = atan(B, A)` it is
`R cos(h − φ) = W − C`, with zero, one or two roots. Each baseline gives at most four critical
hour angles, and the midpoint of each sub-interval is tested, as in ASPRO 2's
`DelayLineService.findHAIntervalsForBaseLine`.
"""
function delay_ha_intervals(facility::FacilityConfig, dec::Float64,
                            config::Vector{Int}, pop::Vector{Int};
                            delay_length::Union{Nothing,Float64}=nothing,
                            pop_array::Union{Nothing,Matrix{Float32}}=nothing,
                            airpath::Union{Nothing,Vector{Float32}}=nothing,
                            beam_order::Union{Nothing,AbstractVector{<:Integer}}=nothing,
                            ha_range::Tuple{Float64,Float64}=(-12.0, 12.0))
    nbl, bxyz, bst_idx, _ = get_baselines(facility; config = config)
    nbl == 0 && return Tuple{Float64,Float64}[]

    # The same fixed offsets and limits `in_delay` uses, so the two cannot disagree.
    off  = _delay_offsets(facility, config, pop; pop_array, airpath, beam_order)
    dlen = _delay_spans(facility, delay_length)

    l, δ = facility.lat * π / 180, dec * π / 180
    acc = [ha_range]
    for b in 1:nbl
        A, B, C = _w_coeffs(l, δ, bxyz[1, b], bxyz[2, b], bxyz[3, b])
        dmax = min(dlen[bst_idx[1, b]], dlen[bst_idx[2, b]])
        acc = _isect(acc, _w_intervals(A, B, C, off[b] - 2dmax, off[b] + 2dmax, ha_range))
        isempty(acc) && return Tuple{Float64,Float64}[]
    end
    return acc
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
                  pop_array::Union{Nothing,Matrix{Float32}}=nothing,
                  airpath::Union{Nothing,Vector{Float32}}=nothing,
                  beam_order::Union{Nothing,AbstractVector{<:Integer}}=nothing)

    delay_carts, nbaselines, baseline_names, baseline_stations =
        compute_delays(facility, dec, ha, config, pop; pop_array=pop_array, airpath=airpath,
                       beam_order=beam_order)

    # Per-baseline delay limit = min of the two telescopes' cart spans. The limit switches are
    # the measurement; `delay_lengths` is their difference and is what a facility file that
    # gives no switches carries instead.
    dlens = Float32.(_delay_spans(facility, delay_length))

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
    # The EXACT hour-angle edges, beside the sampled indices. `good_delay` can only resolve a
    # boundary to the sampling step; `intervals` is where it actually falls, so a caller that
    # wants to draw or quote an edge need not inherit the grid's resolution.
    intervals = delay_ha_intervals(facility, dec, config, pop;
                                   delay_length, pop_array, airpath, beam_order)
    return (delay_carts=delay_carts, has_delay=has_delay, good_delay=good_delay,
            intervals=intervals,
            nbaselines=nbaselines, baseline_names=baseline_names,
            baseline_stations=baseline_stations)
end

# ─── Observability filter for simulate() ─────────────────────────────────────

"""
    night_window(obsdate, lat, lon; zenith = 102.0) -> (t_set, t_rise)

The dark window of the night that begins on `obsdate`, as UTC `DateTime`s: the sunset and
sunrise bracketing the local midnight after `obsdate`'s evening, at any longitude.

`zenith` is the Sun's zenith angle defining the boundary: 90°50′ for geometric sunrise/sunset,
96 civil, 102 nautical, 108 astronomical.
"""
function night_window(obsdate, lat, lon; zenith::Float64 = 102.0)
    # Local midnight closing the evening of `obsdate`, in UTC. Mean solar time from the
    # longitude is ample for picking WHICH night; nothing here needs civil time or a zone.
    mid = Dates.DateTime(Dates.Date(obsdate)) +
          Dates.Millisecond(round(Int, (24.0 - lon / 15) * 3.6e6))
    at(d, h) = Dates.DateTime(d) + Dates.Millisecond(round(Int, h * 3.6e6))
    sets = Dates.DateTime[]; rises = Dates.DateTime[]
    for k in -1:1
        d = Dates.Date(mid) + Dates.Day(k)
        r, st = sunrise_sunset(Dates.DateTime(d), lat, lon; zenith = zenith)
        push!(sets, at(d, st)); push!(rises, at(d, r))
    end
    before = filter(t -> t <= mid, sets)
    after  = filter(t -> t >= mid, rises)
    # Polar day or night: no boundary brackets the midnight, so fall back to the nearest
    # events rather than throwing. The caller sees a window; whether it is dark is `good_alt`'s
    # business and the Sun's, not this function's.
    isempty(before) && (before = sets)
    isempty(after)  && (after  = rises)
    return (maximum(before), minimum(after))
end

"""
    horizon_limit(facility, az; config = Int[]) -> Float64

Minimum elevation, degrees, at azimuth `az`, for the telescopes `config` selects.

The maximum over the horizon profiles of the telescopes in use. `-Inf` when the facility has no
horizon profiles. (At Mount Wilson the profiles range from 18° in the south to 48° in the
north-west.)
"""
function horizon_limit(facility::FacilityConfig, az::Real; config::AbstractVector{<:Integer} = Int[])
    isempty(facility.horizon_az) && return -Inf
    return _horizon_at(facility, _horizon_used(facility, config), Float64(az))
end

# The same limit over a whole azimuth grid. The telescopes in use are derived ONCE: deriving
# them per sample allocated a fresh index vector for every minute of the night, which was more
# than half the cost of this check and 44% of `night_observability`.
function horizon_limit(facility::FacilityConfig, az::AbstractVector;
                       config::AbstractVector{<:Integer} = Int[])
    isempty(facility.horizon_az) && return fill(-Inf, length(az))
    used = _horizon_used(facility, config)
    return [_horizon_at(facility, used, Float64(a)) for a in az]
end

_horizon_used(facility::FacilityConfig, config::AbstractVector{<:Integer}) =
    isempty(config) ? collect(eachindex(facility.horizon_az)) : findall(!iszero, config)

function _horizon_at(facility::FacilityConfig, used::Vector{Int}, az::Float64)
    a = mod(az, 360.0)
    lim = -Inf
    for t in used
        t <= length(facility.horizon_az) || continue
        xs, ys = facility.horizon_az[t], facility.horizon_el[t]
        (isempty(xs) || length(xs) != length(ys)) && continue
        # The profiles are given in increasing azimuth and close at 360, so a plain scan is
        # enough; interpolating keeps a 2° sampling from quantising the limit.
        if a <= xs[1]
            lim = max(lim, ys[1])
        elseif a >= xs[end]
            lim = max(lim, ys[end])
        else
            i = searchsortedlast(xs, a)
            i = clamp(i, 1, length(xs) - 1)
            d = xs[i+1] - xs[i]
            lim = max(lim, d == 0 ? ys[i] : ys[i] + (a - xs[i]) * (ys[i+1] - ys[i]) / d)
        end
    end
    return lim
end

"""
    observable_epochs(facility, target, dates; min_elevation=nothing, max_elevation=nothing,
                      pops=nothing, config=Int[], delay_length=nothing)

The epochs in `dates` at which `target` is observable from `facility`. Each constraint is off
unless its keyword is given; `simulate` calls this only with `observability`.

- `min_elevation`, `max_elevation`: degrees. `nothing` (default) applies no elevation cut.
- `pops`: a POP configuration, one entry per telescope, values in `1:5`. `nothing` (default):
  no delay-line check. [`best_pop`](@ref) proposes configurations.
- `config`: telescope-use flags for the delay check (`1` use, `0` skip, `2` reference cart);
  defaults to using every telescope, i.e. all `ntel*(ntel-1)/2` baselines must fit.
- `delay_length`: override the per-telescope delay-line lengths with one uniform value (m).

`target.raep0` / `target.decep0` are in degrees.

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

# How a Gantt annotates one observing run. Here, in the core, because BOTH renderers obey it --
# `gantt_onenight` in the matplotlib extension and `gantt_geometry` in the GUI one -- and
# sibling extensions cannot import from each other, so neither could own it. `test_gantt.jl`
# compares the labels they produce string for string, and a rule applied on one side only shows
# up there as a chart that disagrees with itself.
#
# The shortest run, in hours, that gets annotated at all:
const GANTT_MIN_LABELLED_RUN = 0.4
# and the width a run needs to hold its az and elev numbers, as a fraction of the plotted span.
# A fraction rather than a number of hours because the numbers are a fixed width in PIXELS, so
# what they cost in hours depends on how many hours the axis is showing.
const GANTT_LABEL_ROOM_FRAC = 0.07

# ─── Best POP search ─────────────────────────────────────────────────────────

"""
    best_pop(facility, dec, ha, config; n_best=5, min_minutes=10, delay_length=nothing)

Search over all POP combinations.

Returns a vector of NamedTuples `(pop, score)` sorted by score (minutes observable),
keeping up to `n_best` results with score ≥ `min_minutes`. The score is the length of the
hour-angle intervals where every baseline is within its delay limits
([`delay_ha_intervals`](@ref)), within the span of `ha`.

Arguments:
- `config`: telescope configuration (0/1/2); only telescopes with config>0 are used
- `n_best`: number of top solutions to return
- `min_minutes`: discard solutions below this threshold
- `delay_length`: override per-telescope delay limits with a uniform value (m)
- `beam_order`: the order the telescopes feed the combiner's channels, which selects each
  one's beam sampling table entry; without one the facility's own order is used
"""
function best_pop(facility::FacilityConfig, dec::Float64, ha::Vector{Float32},
                  config::Vector{Int};
                  n_best::Int=5, min_minutes::Int=10,
                  delay_length::Union{Nothing,Float64}=nothing,
                  pop_array::Union{Nothing,Matrix{Float32}}=nothing,
                  airpath::Union{Nothing,Vector{Float32}}=nothing,
                  beam_order::Union{Nothing,AbstractVector{<:Integer}}=nothing)

    nbaselines, baseline_xyz, baseline_stations, baseline_names = get_baselines(facility; config=config)

    # The same tables `in_delay` resolves, so a POP this search recommends is scored against
    # the delays the plan will then compute from it.
    pops = _pop_table(facility; pop_array)
    path = _fixed_path(facility, config; airpath, beam_order)
    delay_fixed = Float64[path[baseline_stations[2, i]] - path[baseline_stations[1, i]]
                          for i in 1:nbaselines]
    dlens = _delay_spans(facility, delay_length)

    # Active telescope indices
    active = findall(config .> 0)
    npops = size(pops, 2)
    nactive = length(active)

    # The score is observable time within THIS night, so the search is restricted to the hour
    # angles the caller sampled rather than the whole ±12 h: a combination that is wonderful at
    # noon is worth nothing.
    ha_range = isempty(ha) ? (-12.0, 12.0) :
               (Float64(minimum(ha)), Float64(maximum(ha)))
    l, δ = facility.lat * π / 180, dec * π / 180

    results = Tuple{Vector{Int}, Int}[]
    _pop_search!(results, active, nbaselines, baseline_stations,
                 delay_fixed, dlens, pops, npops, nactive, min_minutes,
                 ha_range, l, δ, baseline_xyz)

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

# The POP search, as interval arithmetic rather than a scan.
#
# The POP choice enters a baseline's condition ONLY as a shift of its delay window:
#
#     off_ij + (pop_j - pop_i) - 2d  <=  w(h)  <=  off_ij + (pop_j - pop_i) + 2d
#
# so a baseline's feasible hour angles depend on the whole combination through the single pair
# (pop_i, pop_j). There are 25 such pairs, not 5^6, which is what makes the precompute below
# worth far more than making the inner loop faster: 15 baselines x 25 pairs of CLOSED-FORM
# solves replaces a scan of every combination against every time sample.
#
# What remains is then a depth-first assignment with pruning. Once telescopes 1..k are fixed
# every baseline among them is determined, so the running intersection can be tested at once
# and a subtree abandoned the moment it cannot reach `min_minutes` -- the same early exit
# ASPRO 2 makes when a baseline is incompatible with a W range.
function _pop_search!(results, active, nbaselines, baseline_stations,
                      delay_fixed, dlens, pop_array, npops, nactive, min_minutes,
                      ha_range, l, δ, bxyz)
    # tbl[b][pi, pj] -- feasible hour angles for baseline b with that POP pair.
    tbl = Vector{Matrix{Vector{Tuple{Float64,Float64}}}}(undef, nbaselines)
    for b in 1:nbaselines
        i1, i2 = baseline_stations[1, b], baseline_stations[2, b]
        A, B, C = _w_coeffs(l, δ, bxyz[1, b], bxyz[2, b], bxyz[3, b])
        dmax = min(dlens[i1], dlens[i2])
        m = Matrix{Vector{Tuple{Float64,Float64}}}(undef, npops, npops)
        for pi in 1:npops, pj in 1:npops
            off = delay_fixed[b] +
                  Float64(pop_array[i2, pj]) - Float64(pop_array[i1, pi])
            m[pi, pj] = _w_intervals(A, B, C, off - 2dmax, off + 2dmax, ha_range)
        end
        tbl[b] = m
    end

    # Baselines whose two telescopes are both among the first k assigned, grouped by the k at
    # which that becomes true -- so each baseline is intersected exactly once, as late as
    # possible and as early as it can be.
    pos = Dict(t => k for (k, t) in enumerate(active))
    at_depth = [Int[] for _ in 1:nactive]
    for b in 1:nbaselines
        i1, i2 = baseline_stations[1, b], baseline_stations[2, b]
        (haskey(pos, i1) && haskey(pos, i2)) || continue
        push!(at_depth[max(pos[i1], pos[i2])], b)
    end

    min_hours = min_minutes / 60
    pop_active = ones(Int, nactive)

    function recurse(depth, acc)
        if depth > nactive
            score = round(Int, _ilen(acc) * 60)
            score >= min_minutes && push!(results, (copy(pop_active), score))
            return
        end
        for p in 1:npops
            pop_active[depth] = p
            cur = acc
            ok = true
            for b in at_depth[depth]
                i1, i2 = baseline_stations[1, b], baseline_stations[2, b]
                cur = _isect(cur, tbl[b][pop_active[pos[i1]], pop_active[pos[i2]]])
                # Pruning: an intersection only shrinks, so a subtree that is already too
                # short cannot be rescued by the telescopes still to be assigned.
                if _ilen(cur) < min_hours
                    ok = false
                    break
                end
            end
            ok && recurse(depth + 1, cur)
        end
    end
    recurse(1, [ha_range])
    return results
end

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

# ─── Multi-POP night plans ────────────────────────────────────────────────────

"""
    PopSegment

One stretch of a [`PopPlan`](@ref): hour angles `ha0`–`ha1` (hours, on the plan's 1-minute
grid), the telescopes used (`config`, 0/1/2 as everywhere, full length) and the POPs (`pop`,
full length; a dropped telescope keeps the POP it is parked on). `tier` is 1 with every selected
telescope, 2 with one dropped that is not the reference, 3 with the reference dropped, 4 and 5
the same with two dropped.
"""
struct PopSegment
    ha0    :: Float64
    ha1    :: Float64
    config :: Vector{Int}
    pop    :: Vector{Int}
    tier   :: Int
end

"""
    PopPlan

The result of [`pop_plan`](@ref): `segments` in time order, `minutes` observed in each of the five
tiers (see [`PopSegment`](@ref)), and the number of POP `changes`.
"""
struct PopPlan
    segments :: Vector{PopSegment}
    minutes  :: NTuple{5,Int}
    changes  :: Int
end

"""
    pop_plan(facility, dec, ha, config; observable = trues(length(ha)), reference = 0,
             max_changes = 2, max_drop = 0, change_minutes = 10, min_minutes = 10,
             delay_length = nothing, pop_array = nothing, airpath = nothing) -> PopPlan

The sequence of POP settings, and of telescopes dropped from `config`, that observes a target
longest over the hour angles `ha` where `observable` is true.

Plans are ranked by the minutes observed with every selected telescope, then with one dropped
that is not the reference, then with the reference dropped, then the same with two dropped, then
by the fewest POP changes, and last by the fewest POPs moved over all the changes. Each
comparison only breaks ties in the previous one, so time with the full array is never given up
for time with fewer telescopes.

- `reference`: the telescope protected from being dropped while another drop gains the same
  time; by default the telescope `config` marks 2, if any. A telescope marked 2 is also the
  delay reference, as in [`in_delay`](@ref): a minute is then constrained by the baselines to
  it, and a segment that drops it has another kept telescope as its delay reference (marked 2
  in the segment's `config`).
- `max_changes`: POP changes allowed; `0` is the single setting [`best_pop`](@ref) finds.
- `max_drop`: telescopes that may be dropped at once (0, 1 or 2).
- `change_minutes`: the time a POP change takes, during which nothing is observed. Dropping or
  re-adding a telescope moves no POP and costs nothing.
- `min_minutes`: POP settings observable for less than this are not considered.

The search runs on a 1-minute grid over the observable span. The delay offsets are those of the
whole selection: a dropped telescope's beam keeps its channel. [`pop_plans`](@ref) gives the
best plan for each number of changes from the same search.
"""
pop_plan(facility::FacilityConfig, dec::Real, ha::AbstractVector{<:Real},
         config::AbstractVector{<:Integer}; kwargs...) = last(pop_plans(facility, dec, ha, config; kwargs...))

"""
    pop_plans(facility, dec, ha, config; kwargs...) -> Vector{PopPlan}

The best plan with at most 0, 1, …, `max_changes` POP changes, in that order, from one search;
the keywords are [`pop_plan`](@ref)'s, whose result is the last of these.
"""
function pop_plans(facility::FacilityConfig, dec::Real, ha::AbstractVector{<:Real},
                   config::AbstractVector{<:Integer};
                  observable::AbstractVector{Bool} = trues(length(ha)),
                  reference::Integer = something(findfirst(==(2), config), 0),
                  max_changes::Integer = 2, max_drop::Integer = 0,
                  change_minutes::Real = 10, min_minutes::Real = 10,
                  delay_length::Union{Nothing,Float64} = nothing,
                  pop_array::Union{Nothing,Matrix{Float32}} = nothing,
                  airpath::Union{Nothing,Vector{Float32}} = nothing)
    cfg = Int.(config)
    active = findall(>(0), cfg)
    nactive = length(active)
    length(observable) == length(ha) ||
        throw(ArgumentError("observable has $(length(observable)) entries for $(length(ha)) hour angles"))
    empty = [PopPlan(PopSegment[], (0, 0, 0, 0, 0), 0) for _ in 0:max(0, Int(max_changes))]
    (nactive < 2 || !any(observable)) && return empty
    max_drop = clamp(Int(max_drop), 0, min(2, nactive - 2))

    # The 1-minute grid over the observable span, and which of its minutes the sky allows: each
    # takes the observability of the nearest sample of `ha`.
    hs = Float64.(ha)
    obs_i = findall(observable)
    lo, hi = minimum(hs[obs_i]), maximum(hs[obs_i])
    grid = collect(lo:1/60:hi)
    T = length(grid)
    order = sortperm(hs)
    sky = falses(T)
    for (t, h) in enumerate(grid)
        j = searchsortedfirst(view(hs, order), h)
        cands = (order[clamp(j - 1, 1, length(hs))], order[clamp(j, 1, length(hs))])
        k = abs(hs[cands[1]] - h) <= abs(hs[cands[2]] - h) ? cands[1] : cands[2]
        sky[t] = observable[k]
    end

    # Each pair's in-delay minutes for every POP pair, from the same closed-form intervals
    # `best_pop` intersects. Every pair, whatever the reference: which of them constrain a
    # minute depends on the telescopes kept and on which of them is the reference.
    cfg1 = Int.(cfg .> 0)
    nbl, bxyz, bst, _ = get_baselines(facility; config = cfg1)
    pairidx = Dict{Tuple{Int,Int},Int}()
    for b in 1:nbl
        pairidx[(bst[1, b], bst[2, b])] = b; pairidx[(bst[2, b], bst[1, b])] = b
    end
    pops = _pop_table(facility; pop_array)
    path = _fixed_path(facility, cfg1; airpath)
    dlens = _delay_spans(facility, delay_length)
    npops = size(pops, 2)
    l, δ = facility.lat * π / 180, Float64(dec) * π / 180
    wrap(h) = mod(h + 12, 24) - 12
    inD = Array{BitVector}(undef, nbl, npops, npops)
    for b in 1:nbl
        i1, i2 = bst[1, b], bst[2, b]
        A, B, C = _w_coeffs(l, δ, bxyz[1, b], bxyz[2, b], bxyz[3, b])
        dmax = min(dlens[i1], dlens[i2])
        for pi in 1:npops, pj in 1:npops
            off = (path[i2] - path[i1]) + Float64(pops[i2, pj]) - Float64(pops[i1, pi])
            iv = _w_intervals(A, B, C, off - 2dmax, off + 2dmax, (-12.0, 12.0))
            inD[b, pi, pj] = BitVector(any(x -> x[1] <= wrap(h) <= x[2], iv) for h in grid)
        end
    end

    # The drops allowed, in order of preference, with their tiers: none; one telescope other than
    # the reference; the reference; two without the reference; two with it. With a delay
    # reference (`config` 2) a minute is constrained by the baselines to it, as in `in_delay`;
    # when the reference is dropped, another kept telescope takes its place, each one tried.
    # Without one, by every pair kept.
    dref = something(findfirst(==(2), cfg), 0)
    dsets = Vector{Int}[Int[]]
    max_drop >= 1 && append!(dsets, [[t] for t in active])
    max_drop >= 2 && append!(dsets, [[active[a], active[b]] for a in 1:nactive for b in a+1:nactive])
    drops = Tuple{Vector{Int},Int,Vector{Int},Int}[]      # (dropped, tier, baselines, reference)
    for d in dsets
        kept = setdiff(active, d)
        tier = isempty(d) ? 1 : (length(d) == 1 ? (reference in d ? 3 : 2) : (reference in d ? 5 : 4))
        if dref == 0
            push!(drops, (d, tier, [pairidx[(a, b)] for a in kept for b in kept if a < b], 0))
        else
            for r in (dref in kept ? [dref] : kept)
                push!(drops, (d, tier, [pairidx[(r, x)] for x in kept if x != r], r))
            end
        end
    end
    sort!(drops; by = x -> x[2])     # stable: preference within a tier is the facility's order
    keep_bl = [x[3] for x in drops]

    # States: every POP setting of the selected telescopes, with the best drop at each minute.
    # Settings that observe the same minutes in the same TIERS are one state, whichever drops
    # get them there: the score cannot tell them apart, and the first one found is reported.
    nst = npops^nactive
    nst > 5^8 && throw(ArgumentError("pop_plan: $(nactive) telescopes is more than the search handles"))
    seen = Dict{Vector{UInt8},Int}()
    state_pop = Vector{Vector{Int}}()
    state_gain = Vector{Vector{UInt8}}()      # per minute: 0 = nothing observed, else the drop
    state_tier = Vector{Vector{UInt8}}()      # per minute: 0, else the drop's tier
    state_members = Vector{Vector{Tuple{Vector{Int},Vector{UInt8}}}}()   # every (POPs, drops)
    p = ones(Int, nactive)
    full = ones(Int, facility.ntel)
    gain = zeros(UInt8, T); tier = zeros(UInt8, T)
    m = trues(T); mc = m.chunks
    left = falses(T); lc = left.chunks
    for code in 0:nst-1
        c = code
        for k in 1:nactive
            p[k] = c % npops + 1; c ÷= npops
        end
        for (k, t) in enumerate(active); full[t] = p[k]; end
        fill!(gain, 0x00); fill!(tier, 0x00)
        copyto!(lc, sky.chunks)
        nobs = 0
        for (di, d) in enumerate(drops)
            any(left) || break
            copyto!(mc, lc)
            for b in keep_bl[di]
                xc = inD[b, full[bst[1, b]], full[bst[2, b]]].chunks
                @inbounds for i in eachindex(mc)
                    mc[i] &= xc[i]
                end
                any(m) || break
            end
            # the minutes this drop gains, read off the set bits
            @inbounds for (w, chunk) in enumerate(mc)
                while chunk != 0
                    t = 64 * (w - 1) + trailing_zeros(chunk) + 1
                    gain[t] = UInt8(di); tier[t] = UInt8(d[2]); nobs += 1
                    chunk &= chunk - 1
                end
                lc[w] &= ~mc[w]
            end
        end
        nobs >= min_minutes || continue
        # Every setting of a state is kept: the times cannot tell them apart, but the POPs moved
        # at a change can, and the choice is made once the plan is known (`_pop_plan_trace`).
        if haskey(seen, tier)
            push!(state_members[seen[tier]], (copy(full), copy(gain)))
            continue
        end
        seen[copy(tier)] = length(state_pop) + 1
        push!(state_members, [(copy(full), copy(gain))])
        push!(state_pop, copy(full))
        push!(state_gain, copy(gain))
        push!(state_tier, copy(tier))
    end
    K = length(state_pop)
    K == 0 && return empty
    # the tier of each state at each minute, states contiguous for the inner loop below
    tiermat = Matrix{UInt8}(undef, K, T)
    for s in 1:K, t in 1:T
        tiermat[s, t] = state_tier[s][t]
    end

    # Dynamic programming over the minutes. A score is (minutes per tier..., -changes), compared
    # lexicographically, which adding a minute or a change preserves. Being in state s at minute
    # t comes either from s at t-1, or from another state at t-Δ-1 through a change occupying
    # the Δ minutes between.
    Δ = max(0, round(Int, change_minutes))
    kmax = Int(max_changes)
    NEG = (typemin(Int) ÷ 4, 0, 0, 0, 0, 0)
    tiers = [d[2] for d in drops]
    # one minute in tier `tr` (0: none)
    addmin(s, tr) = ntuple(i -> s[i] + (i == tr), Val(6))
    # Only minutes t-1 and t-Δ-1 are read, so the scores live in a ring of Δ+2 rows; what is
    # kept for every minute is how each state was entered and the two best states.
    nring = Δ + 2
    ring = [fill(NEG, K, kmax + 1) for _ in 1:nring]
    row(t) = ring[mod1(t, nring)]
    came = falses(T, kmax + 1, K)            # true: entered s at t by a change
    top = zeros(Int, T, kmax + 1, 2)         # the two best states at (t, k)
    topval = fill(NEG, T, kmax + 1, 2)
    for t in 1:T
        cur = row(t); prev = t > 1 ? row(t - 1) : cur
        r = t - Δ - 1
        for k in 1:kmax+1
            a = b = 0
            for s in 1:K
                stay = t == 1 ? (k == 1 ? (0, 0, 0, 0, 0, 0) : NEG) : prev[s, k]
                sw = NEG
                if k > 1 && r >= 1
                    j = top[r, k-1, 1] == s ? 2 : 1
                    if top[r, k-1, j] > 0
                        v = topval[r, k-1, j]
                        sw = (v[1], v[2], v[3], v[4], v[5], v[6] - 1)
                    end
                end
                tr = @inbounds tiermat[s, t]
                v = if sw > stay
                    came[t, k, s] = true; addmin(sw, tr)
                else
                    addmin(stay, tr)
                end
                cur[s, k] = v
                if a == 0 || v > topval[t, k, 1]
                    b = a; topval[t, k, 2] = topval[t, k, 1]
                    a = s; topval[t, k, 1] = v
                elseif b == 0 || v > topval[t, k, 2]
                    b = s; topval[t, k, 2] = v
                end
            end
            top[t, k, 1] = a; top[t, k, 2] = b
        end
    end

    # For each number of changes allowed, the best end and back to the start.
    last_row = row(T)
    plans = PopPlan[]
    for c in 0:kmax
        bs, bk = 1, 1
        for s in 1:K, k in 1:c+1
            last_row[s, k] > last_row[bs, bk] && ((bs, bk) = (s, k))
        end
        push!(plans, _pop_plan_trace(bs, bk, T, Δ, came, top, state_members, tiermat, drops,
                                     tiers, grid, cfg1, active))
    end
    return plans
end

# Back from state `bs` with `bk - 1` changes at the last minute to the plan that reached it, then
# the settings that move the fewest POPs over its changes. For each run of one state the
# candidates are every setting observing that run's own minutes in the same tiers -- of any
# state, as what a setting does outside the run does not matter to it.
function _pop_plan_trace(bs, bk, T, Δ, came, top, state_members, tiermat, drops, tiers, grid,
                         cfg1, active)
    path_s = zeros(Int, T)                   # the state at each minute; 0 during a change
    t, s, k = T, bs, bk
    while t >= 1
        path_s[t] = s
        if came[t, k, s]
            r = t - Δ - 1
            src = top[r, k-1, 1] == s ? top[r, k-1, 2] : top[r, k-1, 1]
            t, s, k = r, src, k - 1
        else
            t -= 1
        end
    end

    # The plan's states in order, and for each the setting that moves the fewest POPs over all
    # its changes: a shortest path through the states' settings, the cost of a step being the
    # POPs that differ.
    blocks = Tuple{Int,Int,Int}[]            # (state, first minute, last minute)
    for t in 1:T
        st = path_s[t]
        st == 0 && continue
        if isempty(blocks) || last(blocks)[1] != st || path_s[t-1] == 0
            push!(blocks, (st, t, t))
        else
            blocks[end] = (st, blocks[end][2], t)
        end
    end
    K = size(tiermat, 1)
    cands = map(blocks) do (st, t0, t1)
        obs = [t for t in t0:t1 if tiermat[st, t] != 0]
        same = [s for s in 1:K if all(t -> tiermat[s, t] == tiermat[st, t], obs)]
        reduce(vcat, (state_members[s] for s in same))
    end
    moved(a, b) = count(i -> a[i] != b[i], active)
    pick = Dict{Int,Tuple{Vector{Int},Vector{UInt8}}}()
    if !isempty(blocks)
        nb = length(blocks)
        cost = [zeros(Int, length(cands[1]))]
        back = [Int[]]
        for b in 2:nb
            prevm, curm = cands[b-1], cands[b]
            cb = Vector{Int}(undef, length(curm)); bb = Vector{Int}(undef, length(curm))
            for (j, (pj, _)) in enumerate(curm)
                best, arg = typemax(Int), 1
                for (i, (pi, _)) in enumerate(prevm)
                    v = cost[b-1][i] + moved(pi, pj)
                    v < best && ((best, arg) = (v, i))
                end
                cb[j] = best; bb[j] = arg
            end
            push!(cost, cb); push!(back, bb)
        end
        choice = Vector{Int}(undef, nb)
        choice[nb] = argmin(cost[nb])
        for b in nb:-1:2
            choice[b-1] = back[b][choice[b]]
        end
        for (b, (_, t0, t1)) in enumerate(blocks), t in t0:t1
            pick[t] = cands[b][choice[b]]
        end
    end

    # Segments: runs of one state and one drop, over the minutes something is observed.
    segs = PopSegment[]
    mins = zeros(Int, 5)
    t = 1
    while t <= T
        s = path_s[t]
        mem = s == 0 ? nothing : pick[t]
        g = s == 0 ? 0 : Int(mem[2][t])
        if g == 0
            t += 1; continue
        end
        u = t
        while u < T && path_s[u+1] == s && pick[u+1] === mem && mem[2][u+1] == g
            u += 1
        end
        c = copy(cfg1)
        for d in drops[g][1]; c[d] = 0; end
        drops[g][4] > 0 && (c[drops[g][4]] = 2)
        push!(segs, PopSegment(grid[t], grid[u], c, copy(mem[1]), tiers[g]))
        mins[tiers[g]] += u - t + 1
        t = u + 1
    end
    return PopPlan(segs, ntuple(i -> mins[i], 5), bk - 1)
end

"""
    print_pop_plan(facility, plan)

Pretty-print a [`PopPlan`](@ref): one line per segment, then the minutes per tier.
"""
function print_pop_plan(facility::FacilityConfig, plan::PopPlan)
    println("─── POP plan: $(length(plan.segments)) segments, $(plan.changes) POP changes ───")
    for sg in plan.segments
        used = findall(>(0), sg.config)
        pop_str = join(["$(facility.sta_names[i])-POP$(sg.pop[i])" for i in used], "  ")
        @printf("  HA %+6.2f … %+6.2f h  (%3d min):  %s\n", sg.ha0, sg.ha1,
                round(Int, (sg.ha1 - sg.ha0) * 60) + 1, pop_str)
    end
    labels = ("all telescopes", "one dropped", "reference dropped", "two dropped",
              "two dropped, reference among them")
    for (i, n) in enumerate(plan.minutes)
        n > 0 && println("  $(n) min with $(labels[i])")
    end
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

# ─────────────────────────────────────────────────────────────────────────────
# Contiguous runs of an index vector
# ─────────────────────────────────────────────────────────────────────────────
#
# In the core rather than beside one of the renderers, because BOTH draw from it: a Gantt bar
# is one contiguous run, and a matplotlib chart and a Makie one that each found their own runs
# would be free to disagree about where an observing block starts.

"""
    index_runs(idx) -> Vector{Tuple{Int,Int}}

Split an index vector into maximal runs of consecutive integers, as `(first, last)` pairs; used
to draw observability windows with gaps.
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
