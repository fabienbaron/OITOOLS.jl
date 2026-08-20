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
        # in December it never clears the default 30-degree limit
        @test maximum(dec.alt) < 30
        @test isempty(dec.good_alt)
        @test observable_hours(dec) == 0
    end

    @testset "a far-southern target never rises from a northern site" begin
        p = night_plan("CHARA", "Canopus", CANOPUS_RA, CANOPUS_DEC, DateTime(2026, 2, 1))
        # dec -52.7 from latitude +34.2 cannot exceed 90 - 34.2 - 52.7 ≈ 3 degrees
        @test maximum(p.alt) < 5
        @test observable_hours(p) == 0
    end

    @testset "the June night is short, and it is what limits the target" begin
        p = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE)
        # Vega is above the limit all night, so the observable time IS the dark window — which
        # is the check that twilight is being computed at all rather than assumed.
        @test length(p.good_alt) > length(p.good_twilight)
        @test observable_hours(p) ≈ length(p.good_twilight) / 60 rtol = 1e-6
        # and a June night at mid-latitude is short
        @test 3 < length(p.good_twilight) / 60 < 8
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

    @testset "the moon is reported even when it does not exclude anything" begin
        p = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE)
        @test 0 <= p.moon_fli <= 1
        @test length(p.moon_sep) == length(p.alt)
        @test all(0 .<= p.moon_sep .<= 180)
        @test occursin("moon", sprint(show, p))
    end
end
