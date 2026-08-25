# The Observe perspective's data layer.
#
# These assert ASTRONOMY, not self-consistency. A planning tool that agrees with itself and
# disagrees with the sky is worse than none, and the three bugs this perspective was blocked on
# (twilight in the wrong units, southern declinations with the wrong sign, RA in hours where
# degrees were documented) were all of that kind: plausible numbers, wrong answers. So the
# checks below are things independently known — Vega is a northern summer target, Canopus never
# rises from California, a June night is short — rather than comparisons against whatever the
# code happens to return.

@testset "Observe data layer" begin

    # Vega: high overhead at CHARA (dec +38.8 against latitude +34.2), a summer target.
    VEGA_RA, VEGA_DEC = 279.2347, 38.7837
    # Canopus: far south, and so never usefully above the horizon from California.
    CANOPUS_RA, CANOPUS_DEC = 95.988, -52.696
    JUNE, DECEMBER = DateTime(2026, 6, 21), DateTime(2026, 12, 21)

    @testset "the catalogue links a spectral setup to its combiner" begin
        c = config_catalog()
        @test "CHARA" in c.facilities
        @test !isempty(c.combiners)
        @test isempty(c.unknown)          # every shipped config classifies
        # wave_combiner is the ONLY machine-readable link between the two, and without it a
        # panel offers spectral setups the chosen instrument cannot produce.
        for (wave, comb) in c.wave_combiner
            @test wave in c.wavelengths
            @test comb in c.combiners
        end
        @test !isempty(c.wave_combiner)
    end

    @testset "telescope order is the facility's own" begin
        tels = facility_telescopes("CHARA")
        @test tels == ["S1", "S2", "E1", "E2", "W1", "W2"]
        # compute_delays indexes into exactly this order, so a config vector built from names
        # has to line up with it positionally.
        @test telescope_config(["S1", "E1"], tels) == [1, 0, 1, 0, 0, 0]
        @test telescope_config(String[], tels)     == zeros(Int, 6)
        @test telescope_config(tels, tels)         == ones(Int, 6)
    end

    @testset "Vega is a summer target from CHARA, and not a winter one" begin
        june = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE)
        dec  = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, DECEMBER)
        # dec +38.8 against latitude +34.2 transits within 5 degrees of the zenith
        @test maximum(june.alt) > 80
        @test observable_hours(june) > 3
        # In December it barely scrapes the 25-degree floor -- 28.5 at transit -- so it is
        # nominally "up" for a few minutes and useless. The check is the ratio, not zero:
        # asserting zero would have been an accident of the old 30-degree default.
        @test maximum(dec.alt) < 30
        @test length(dec.good_alt) < length(june.good_alt) / 10
        @test observable_hours(dec) < 0.5
    end

    @testset "a far-southern target never rises from a northern site" begin
        p = night_plan("CHARA", "Canopus", CANOPUS_RA, CANOPUS_DEC, DateTime(2026, 2, 1))
        # dec -52.7 from latitude +34.2 cannot exceed 90 - 34.2 - 52.7 ≈ 3 degrees
        @test maximum(p.alt) < 5
        @test observable_hours(p) == 0
    end

    @testset "the June night is short, and it is what limits the target" begin
        # alt_max = 90 on purpose: with the default 85 the zenith cap ALSO bites (Vega transits
        # at 85.4), and this test is about twilight being computed rather than assumed. The cap
        # gets its own test below.
        p = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE; alt_max = 90.0)
        @test length(p.good_alt) > length(p.good_twilight)
        @test observable_hours(p) ≈ length(p.good_twilight) / 60 rtol = 1e-6
        # and a June night at mid-latitude is short
        @test 3 < length(p.good_twilight) / 60 < 8
    end

    @testset "the defaults are ASPRO's, and the zenith cap is not decorative" begin
        @test DEFAULT_ALT_LIMIT == 25.0
        @test DEFAULT_ALT_MAX   == 85.0
        d = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE)
        @test d.alt_limit == 25.0 && d.alt_max == 85.0
        # Vega transits at 85.4 from CHARA, so the cap removes real time rather than nothing
        open_ = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE; alt_max = 90.0)
        @test observable_hours(d) < observable_hours(open_)
    end

    @testset "the delay check is opt-in, and says when it did not run" begin
        # Applying it with an unsearched POP reports far less time than is available, so it
        # must not be the default — and "not checked" must not look like "passed".
        p = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE)
        @test p.delay_applied == false
        @test p.good_delay == collect(1:length(p.ha))   # covers the night: not a constraint

        cfg = telescope_config(["S1", "S2", "E1", "E2"], facility_telescopes("CHARA"))
        q = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE; config = cfg, use_delay = true)
        @test q.delay_applied == true
        @test observable_hours(q) <= observable_hours(p)   # a constraint can only remove time
    end

    @testset "searching the POPs buys real time" begin
        cfg = telescope_config(["S1", "S2", "E1", "E2"], facility_telescopes("CHARA"))
        base = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE;
                          config = cfg, use_delay = true)               # POP 1 throughout
        pops = best_pops("CHARA", VEGA_DEC, base.ha, cfg; n = 3)
        @test length(pops) == 3
        @test issorted([p.score for p in pops]; rev = true)             # best first
        best = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE;
                          config = cfg, pop = pops[1].pop, use_delay = true)
        # This is the whole point of the POP search: an arbitrary configuration wastes the night
        @test observable_hours(best) > observable_hours(base)
        @test occursin("\t", pop_rows(pops))
    end

    @testset "the summary table reports a window, not just a total" begin
        rows = split(plan_rows([night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE),
                                night_plan("CHARA", "Canopus", CANOPUS_RA, CANOPUS_DEC, JUNE)]), "\n")
        @test length(rows) == 2
        vega = split(rows[1], "\t")
        @test vega[1] == "Vega"
        @test parse(Float64, vega[4]) > 3          # hours
        @test parse(Float64, vega[2]) != parse(Float64, vega[3])   # a real window
        # a target that never rises still gets a row, with zero hours rather than no row
        @test split(rows[2], "\t")[1] == "Canopus"
        @test parse(Float64, split(rows[2], "\t")[4]) == 0
    end

    @testset "SIMBAD resolution answers in the panel's own format" begin
        # The reply is parsed by ObserveTab.qml's applySimbad, which splits on tabs and
        # branches on field 1. Both branches are asserted here because a malformed reply
        # would be read by QML as a successful resolve and write NaNs into the row.
        @test startswith(GUI.shell_simbad(""), "bad\t")
        @test startswith(GUI.shell_simbad("   "), "bad\t")

        # The network is not a test dependency. This runs only when asked for, because a
        # suite that fails on a SIMBAD outage is a suite people learn to ignore.
        if get(ENV, "OITOOLS_SIMBAD_TESTS", "0") == "1"
            f = split(GUI.shell_simbad("Vega"), "\t")
            @test f[1] == "ok"
            @test length(f) == 9
            @test parse(Float64, f[2]) ≈ VEGA_RA  atol = 1e-3   # degrees, not hours
            @test parse(Float64, f[3]) ≈ VEGA_DEC atol = 1e-3
            @test f[8] == "* alf Lyr"
            @test startswith(f[9], "A0")

            # Southern declinations keep their sign: this is P0-2, and a sexagesimal parser
            # that applies the minus to the degrees field alone gets it wrong by 2×(′+″).
            g = split(GUI.shell_simbad("Alpha Cen A"), "\t")
            @test g[1] == "ok"
            @test parse(Float64, g[3]) < -60
            # α Cen A has no published R or I magnitude. That must arrive as an empty field,
            # not as 0.0 — magnitude zero is Vega-bright, and writing it for "unknown" would
            # silently make an unmeasured band the brightest in the row.
            t = simbad_target("Alpha Cen A")
            @test isnan(t.mags["R"])
            @test t.mags["K"] < 0
            @test t.plx > 700                       # the nearest star system

            @test_throws Exception simbad_target("NoSuchStarXYZ")
        end

        # Rejected before the network, so this costs nothing and runs everywhere
        @test_throws ArgumentError simbad_target("  ")
        @test issorted(indexin(["V", "J", "H", "K"], SIMBAD_BANDS))   # panel order
    end

    @testset "the moon is reported even when it does not exclude anything" begin
        p = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE)
        @test 0 <= p.moon_fli <= 1
        @test length(p.moon_sep) == length(p.alt)
        @test all(0 .<= p.moon_sep .<= 180)
        @test occursin("moon", sprint(show, p))
    end
end
