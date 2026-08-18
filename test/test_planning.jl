# test_planning.jl — astrometry and observation-planning regressions.
#
# Why this exists: before 0.11 there was no test of any function in plan.jl or astrometry.jl,
# and no test of the SIMBAD coordinate path. Three bugs lived there undisturbed, all of the
# same kind — a unit convention that was wrong but plausible, so the output looked reasonable
# and nothing crashed:
#
#   1. `sunrise_sunset` fed degrees to `cos()` and a cosine to `acos()` after scaling it by
#      pi/180. It returned a ~12 h night at every latitude on every date of the year.
#   2. Callers reduced `ra_dec_from_simbad`'s sexagesimal components with a dot product
#      against [1, 1/60, 1/3600], which negates only the degrees field. Every southern
#      declination was wrong by 2x(arcmin + arcsec).
#   3. `night_observability` used its `ra` argument as hours for the hour angle and as
#      degrees for the Moon separation, in the same function body.
#
# The assertions below are therefore about *units and signs*, not about high accuracy: the
# underlying almanac algorithm is only good to a few minutes, and that is fine.

using OITOOLS, Test, Dates

@testset "planning and astrometry" begin

    CHARA_LAT, CHARA_LON = 34.2243836944, -118.0570313111

    @testset "sexagesimal_to_degrees" begin
        @test sexagesimal_to_degrees("46 28 00.5731825") ≈ 46.46682588402778 rtol = 1e-12
        # The sign must apply to the whole quantity, not just the leading field.
        @test sexagesimal_to_degrees("-46 28 00.5731825") ≈ -46.46682588402778 rtol = 1e-12
        # ...including when the degrees field is zero, where a naive dot product loses the
        # sign entirely and flips hemisphere.
        @test sexagesimal_to_degrees("-00 30 00.0") ≈ -0.5 rtol = 1e-12
        @test sexagesimal_to_degrees("+00 30 00.0") ≈  0.5 rtol = 1e-12
        # Separators and field counts
        @test sexagesimal_to_degrees("12:34:56.7") ≈ sexagesimal_to_degrees("12 34 56.7")
        @test sexagesimal_to_degrees("-12:30")     ≈ -12.5
        @test sexagesimal_to_degrees("15")         ≈ 15.0
        @test_throws ArgumentError sexagesimal_to_degrees("   ")
    end

    @testset "sunrise_sunset: seasonal variation" begin
        darkness(d, lat, lon) = begin
            _, set  = sunrise_sunset(d, lat, lon)
            rise, _ = sunrise_sunset(d + Day(1), lat, lon)
            rise - set
        end

        jun = darkness(DateTime(2026, 6, 21), CHARA_LAT, CHARA_LON)
        dec = darkness(DateTime(2026, 12, 21), CHARA_LAT, CHARA_LON)
        mar = darkness(DateTime(2026, 3, 21), CHARA_LAT, CHARA_LON)

        # The bug: all three were ~12.0 h. Nights at 34 N must be much shorter in June than
        # in December, with the equinox in between.
        @test dec - jun > 3.5
        @test jun < mar < dec
        # Sanity envelope for nautical twilight at this latitude.
        @test 6.5 < jun < 8.5
        @test 11.0 < dec < 13.5

        # Southern hemisphere must invert the seasons.
        jun_s = darkness(DateTime(2026, 6, 21), -24.627, -70.405)   # VLTI / Paranal
        dec_s = darkness(DateTime(2026, 12, 21), -24.627, -70.405)
        @test jun_s > dec_s

        # Local clock times must be plausible: dusk in the evening, dawn before morning.
        _, set = sunrise_sunset(DateTime(2026, 12, 21), CHARA_LAT, CHARA_LON)
        rise, _ = sunrise_sunset(DateTime(2026, 12, 22), CHARA_LAT, CHARA_LON)
        local_set  = mod(set  - 8, 24)      # Mt Wilson is UTC-8 (no DST applied here)
        local_rise = mod(rise - 8, 24)
        @test 16.5 < local_set  < 19.0
        @test  4.5 < local_rise <  7.0
    end

    @testset "sunrise_sunset: polar latitudes do not throw" begin
        # |cosH| > 1 there — acos would raise DomainError without the clamp. The returned
        # window is degenerate rather than meaningful; the contract is only that it is
        # finite and does not blow up, since both supported arrays are mid-latitude.
        for d in (DateTime(2026, 6, 21), DateTime(2026, 12, 21)), lat in (78.0, -78.0)
            r, s = sunrise_sunset(d, lat, 15.0)
            @test isfinite(r) && isfinite(s)
        end
    end

    @testset "night_observability takes RA in degrees" begin
        facility = read_facility_file("CHARA")
        obsdate  = DateTime(2026, 6, 3)
        # Vega: 18h36m56.3s = 279.2347 deg, +38d47m01s
        ra_deg, dec_deg = 279.23473479, 38.78368896

        obs = night_observability(facility, ra_deg, dec_deg, obsdate; alt_limit = 30.0)

        @test length(obs.ha) == length(obs.alt) == length(obs.az) == length(obs.moon_sep)
        @test !isempty(obs.good_alt)                    # Vega is circumpolar-ish from CHARA
        @test all(-12 .<= obs.ha .<= 12)
        @test all(-90 .<= obs.alt .<= 90)
        @test all(0 .<= obs.az .< 360)

        # Moon separation is an angle on the sky: it must be in [0,180]. When `ra` was
        # consumed as hours this silently compared an hours value against the Moon's degrees.
        @test all(0 .<= obs.moon_sep .<= 180)
        @test 0 <= obs.moon_fli <= 1

        # The altitudes must be the ones alt_az gives for the same hour angles — i.e. the
        # RA -> hour-angle conversion inside is self-consistent.
        alt_direct, _ = alt_az(dec_deg, facility.lat, Float64.(obs.ha))
        @test Float32.(alt_direct) ≈ obs.alt rtol = 1e-4

        # Feeding the same target as hours must NOT reproduce the degrees answer -- this is
        # what makes the convention observable rather than a matter of taste.
        obs_h = night_observability(facility, ra_deg/15, dec_deg, obsdate; alt_limit = 30.0)
        @test !(obs_h.ha ≈ obs.ha)
    end

    @testset "supporting astrometry" begin
        @test airmass(90.0) ≈ 1.0 rtol = 1e-6
        @test airmass(30.0) ≈ 2.0 rtol = 1e-2
        @test airmass(0.0) > 50                       # clamped at 0.5 deg, matches ASPRO
        @test datetime_to_mjd(DateTime(2000,1,1,12,0,0)) ≈ 51544.5 atol = 1e-6
        @test datetime_to_jd(DateTime(2000,1,1,12,0,0)) ≈ 2451545.0 atol = 1e-6
        @test angular_separation(0.0, 0.0, 0.0, 90.0) ≈ 90.0 rtol = 1e-9
        @test angular_separation(10.0, -5.0, 10.0, -5.0) ≈ 0.0 atol = 1e-9
        # alt_az at transit (ha=0) must put the target at 90 - |lat - dec|
        alt0, _ = alt_az(34.0, 34.0, 0.0)
        @test alt0 ≈ 90.0 rtol = 1e-6
    end
end
