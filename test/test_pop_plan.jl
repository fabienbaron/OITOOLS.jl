# test_pop_plan.jl — multi-POP night plans (`pop_plan`).

using OITOOLS, Test, Dates

@testset "pop_plan" begin
    f = read_facility_file("CHARA")

    # The hours a target is above the elevation limit and the horizon, in astronomical night.
    function sky(ra, dec, date; cfg = ones(Int, 6))
        o = night_observability(f, ra, dec, date; alt_limit = 25.0, config = cfg)
        ok = falses(length(o.ha))
        ok[intersect(o.good_alt, o.good_twilight)] .= true
        return o, ok
    end

    # Every minute of every segment is in delay, by `in_delay` on the segment's own telescopes
    # and delay reference, with the whole selection's beam channels (a dropped telescope's beam
    # keeps its channel).
    function segments_in_delay(plan, dec, cfg)
        for sg in plan.segments
            hs = Float32.(collect(sg.ha0:1/60:sg.ha1))
            r = in_delay(f, Float64(dec), hs, sg.config, sg.pop; beam_order = findall(>(0), cfg))
            length(r.good_delay) == length(hs) || return false
        end
        return true
    end

    vega = (279.2347, 38.7837, DateTime(2026, 6, 21))
    hd217357 = (345.067176, -22.524347, DateTime(2026, 10, 26))    # SIMBAD

    @testset "one setting is best_pop's best" begin
        o, ok = sky(vega...)
        bp = best_pop(f, vega[2], o.ha[ok], ones(Int, 6); n_best = 1)
        p = pop_plan(f, vega[2], o.ha, ones(Int, 6); observable = ok, max_changes = 0)
        @test length(p.segments) == 1 && p.changes == 0
        @test p.segments[1].pop == bp[1].pop
        @test abs(p.minutes[1] - bp[1].score) <= 1
        @test segments_in_delay(p, vega[2], ones(Int, 6))
    end

    @testset "a change is made only when it gains more than it costs" begin
        o, ok = sky(vega...)
        cfg = ones(Int, 6)
        p0 = pop_plan(f, vega[2], o.ha, cfg; observable = ok, max_changes = 0)
        p1 = pop_plan(f, vega[2], o.ha, cfg; observable = ok, max_changes = 2)
        @test p1.changes >= 1
        @test p1.minutes[1] > p0.minutes[1] + 10
        @test segments_in_delay(p1, vega[2], cfg)
        # the 10 minutes of the change are not observed
        for (a, b) in zip(p1.segments[1:end-1], p1.segments[2:end])
            a.pop != b.pop && @test (b.ha0 - a.ha1) * 60 >= 10
        end
        # a change dearer than anything it could gain is not made
        pnever = pop_plan(f, vega[2], o.ha, cfg; observable = ok, max_changes = 2,
                          change_minutes = 600)
        @test pnever.changes == 0 && pnever.minutes == p0.minutes
    end

    @testset "time with every telescope is never given up" begin
        for t in (vega, hd217357)
            o, ok = sky(t...)
            a = pop_plan(f, t[2], o.ha, ones(Int, 6); observable = ok, max_changes = 2)
            b = pop_plan(f, t[2], o.ha, ones(Int, 6); observable = ok, max_changes = 2,
                         max_drop = 2)
            @test b.minutes[1] == a.minutes[1]
            @test sum(b.minutes) >= sum(a.minutes)
            @test segments_in_delay(b, t[2], ones(Int, 6))
        end
        # Vega: one change gives the full array most of the night, then E1 is dropped for the
        # last half hour, with no further change
        o, ok = sky(vega...)
        p = pop_plan(f, vega[2], o.ha, ones(Int, 6); observable = ok, max_changes = 2,
                     max_drop = 1)
        @test p.changes == 1
        @test p.minutes[1] > 300 && p.minutes[2] > 0
        @test p.segments[end].config[3] == 0 && p.segments[end].pop == p.segments[end-1].pop
    end

    @testset "HD 217357 from CHARA on 2026-10-26" begin
        # Nothing with all six telescopes while it is above the horizon; without E1 the whole
        # 72 minutes.
        o, ok = sky(hd217357...)
        dec = hd217357[2]
        p6 = pop_plan(f, dec, o.ha, ones(Int, 6); observable = ok, max_changes = 2)
        @test isempty(p6.segments)
        cfg = [1, 1, 1, 1, 1, 2]                       # W2 the reference
        p = pop_plan(f, dec, o.ha, cfg; observable = ok, max_changes = 2, max_drop = 1)
        @test p.minutes == (0, 72, 0, 0, 0)
        @test length(p.segments) == 1 && p.changes == 0
        @test p.segments[1].config == [1, 1, 0, 1, 1, 2]
        @test segments_in_delay(p, dec, cfg)
        # E1 as the reference: it is still the one dropped, as no other drop gains anything,
        # and the minutes count as the reference's
        pr = pop_plan(f, dec, o.ha, [1, 1, 2, 1, 1, 1]; observable = ok, max_changes = 2,
                      max_drop = 1)
        @test pr.minutes == (0, 0, 72, 0, 0)
        # and another kept telescope is the delay reference for that segment
        @test pr.segments[1].config[3] == 0 && count(==(2), pr.segments[1].config) == 1
        @test segments_in_delay(pr, dec, [1, 1, 2, 1, 1, 1])
    end

    @testset "unwrap_hour_angles" begin
        @test unwrap_hour_angles([11.9, -11.9, -11.8]) ≈ [11.9, 12.1, 12.2]
        @test unwrap_hour_angles([18.3, 18.7, -5.3]) ≈ [-5.7, -5.3, -5.3]
        @test unwrap_hour_angles([-3.0, -2.0, -1.0]) == [-3.0, -2.0, -1.0]
        @test eltype(unwrap_hour_angles(Float32[1, 2])) == Float32
        @test eltype(unwrap_hour_angles([1, 2])) == Float64
        @test isempty(unwrap_hour_angles(Float64[]))
    end

    # A night through LST 0 h for a target at RA 5.3 h, observable from HA -5.9 h. Hour angles
    # used to run [-12, 24), jumping from +18.7 to -5.3 at LST 0 h; the plan then held a segment
    # at HA +18.1 after one at -1.3 and counted two POP changes for one.
    capella = (79.172, 45.998, DateTime(2026, 9, 15))
    @testset "Capella from CHARA on 2026-09-15, through LST 0 h" begin
        ra, dec, d = capella
        o, ok = sky(capella...)
        @test all(h -> -12 <= h < 12, o.ha)
        @test unwrap_hour_angles(o.ha) == Float32.(o.ha)          # no jump in the night
        p = pop_plan(f, dec, o.ha, ones(Int, 6); observable = ok, max_changes = 2, max_drop = 1)
        @test p.changes == 1 && p.minutes == (230, 37, 0, 0, 0)
        @test issorted([sg.ha0 for sg in p.segments]) &&
              all(p.segments[k].ha0 > p.segments[k-1].ha1 for k in 2:length(p.segments))
        @test count(k -> p.segments[k].pop != p.segments[k-1].pop, 2:length(p.segments)) == p.changes
        @test segments_in_delay(p, dec, ones(Int, 6))
        # S1 dropped for the first 37 min, then every telescope; E1 5 → 2 for the last 17
        @test p.segments[1].config[1] == 0 && p.segments[1].pop == p.segments[2].pop
        @test p.segments[end].pop[3] == 2 && p.segments[2].pop[3] == 5
        # with W2 as the delay reference, one setting all night
        cfg = ones(Int, 6); cfg[6] = 2
        pw = pop_plan(f, dec, o.ha, cfg; observable = ok, max_changes = 2, max_drop = 1)
        @test pw.changes == 0 && pw.minutes == (277, 0, 0, 0, 0)
        @test only(pw.segments).pop == [1, 1, 4, 4, 1, 5]
    end

    # The every-cart-free delay model (no reference) against ASPRO 2, which checks every pair
    # with the same window. ASPRO's windows, held PoPs, converted to hour angle with its transit:
    #   xvfb-run -a java -cp .:Aspro2-26.09.jar AsproNight --config "CHARA 2026A" \
    #       --stations "E1 W2 W1 S2 S1 E2" --target Capella --ra 05:16:41.28 \
    #       --dec +45:59:52.8 --date 2026-09-16 --minElev 25 --pops 551215
    #   → 08:12:48–11:50:42 UTC, transit 13:28:51; --pops 451215 → from 08:52:52 (to twilight)
    @testset "Capella's delay limits agree with ASPRO 2" begin
        dec = capella[2]
        ha(utc) = (utc - (13 + 28 / 60 + 51 / 3600)) * 1.0027379
        iv = only(delay_ha_intervals(f, dec, ones(Int, 6), [1, 2, 5, 5, 1, 5]))
        @test iv[1] ≈ ha(8 + 12 / 60 + 48 / 3600) atol = 2 / 60     # within 2 minutes
        @test iv[2] ≈ ha(11 + 50 / 60 + 42 / 3600) atol = 2 / 60
        iv4 = only(delay_ha_intervals(f, dec, ones(Int, 6), [1, 2, 4, 5, 1, 5]))
        @test iv4[1] ≈ ha(8 + 52 / 60 + 52 / 3600) atol = 3 / 60
    end

    @testset "the reference is kept when another drop gains the same" begin
        o, ok = sky(vega...)
        # Vega's last half hour needs E1 out; with each other telescope as the reference the
        # drop is E1's, and with E1 as the reference it is still E1 only if no other drop works
        for r in (1, 2, 4, 5, 6)
            cfg = ones(Int, 6); cfg[r] = 2
            p = pop_plan(f, vega[2], o.ha, cfg; observable = ok, max_changes = 2, max_drop = 1)
            @test p.minutes[3] == 0
            @test all(sg -> sg.config[r] > 0, p.segments)
        end
    end

    @testset "pop_plans: the best for each number of changes, from one search" begin
        o, ok = sky(vega...)
        ps = pop_plans(f, vega[2], o.ha, ones(Int, 6); observable = ok, max_changes = 2, max_drop = 1)
        @test length(ps) == 3
        for (c, q) in enumerate(ps)
            @test q.changes <= c - 1
            @test q.minutes == pop_plan(f, vega[2], o.ha, ones(Int, 6); observable = ok,
                                        max_changes = c - 1, max_drop = 1).minutes
        end
        @test ps[end].minutes == pop_plan(f, vega[2], o.ha, ones(Int, 6); observable = ok,
                                          max_changes = 2, max_drop = 1).minutes
    end

    @testset "no POP is moved that did not have to be" begin
        # At each change, putting any one moved telescope back on its previous POP takes some
        # baseline out of delay while the new setting is in force -- over every segment that
        # keeps it, not only the first: otherwise the move was unnecessary.
        cases = (vega, (279.2347, 38.7837, DateTime(2027, 6, 23)),
                 (213.9153, 19.1824, DateTime(2026, 4, 15)), (310.358, 45.2803, DateTime(2026, 8, 20)))
        for t in cases, cfg in (ones(Int, 6), [1, 1, 1, 1, 1, 2]), kd in (0, 1)
            o, ok = sky(t...)
            p = pop_plan(f, t[2], o.ha, cfg; observable = ok, max_changes = 3, max_drop = kd)
            for k in 2:length(p.segments)
                a, b = p.segments[k-1], p.segments[k]
                a.pop == b.pop && continue
                blk = [sg for sg in p.segments[k:end]]
                j = findfirst(sg -> sg.pop != b.pop, blk)
                j === nothing || (blk = blk[1:j-1])
                for i in findall(>(0), cfg)
                    a.pop[i] == b.pop[i] && continue
                    q = copy(b.pop); q[i] = a.pop[i]
                    fits = all(blk) do sg
                        hs = Float32.(collect(sg.ha0:1/60:sg.ha1))
                        r = in_delay(f, Float64(t[2]), hs, sg.config, q; beam_order = findall(>(0), cfg))
                        length(r.good_delay) == length(hs)
                    end
                    @test !fits
                end
            end
        end
    end

    @testset "min_stretch" begin
        # RA 120° Dec +50° on 2026-09-15, W2 the reference, drop ≤ 1: S1 comes back after one
        # minute without S1; with a 15-minute stretch that minute is not observed.
        cfg = [1, 1, 1, 1, 1, 2]
        o, ok = sky(120.0, 50.0, DateTime(2026, 9, 15); cfg)
        a = pop_plan(f, 50.0, o.ha, cfg; observable = ok, max_drop = 1)
        b = pop_plan(f, 50.0, o.ha, cfg; observable = ok, max_drop = 1, min_stretch = 15)
        @test a.minutes == (120, 1, 0, 0, 0) && b.minutes == (120, 0, 0, 0, 0)
        @test length(b.segments) == 1 && segments_in_delay(b, 50.0, cfg)
        # every segment next to another under the same POPs lasts 15 minutes
        o, ok = sky(345.0, 60.0, DateTime(2026, 9, 15); cfg)
        p = pop_plan(f, 60.0, o.ha, cfg; observable = ok, max_drop = 1, min_stretch = 15)
        @test all(sg -> round(Int, (sg.ha1 - sg.ha0) * 60) + 1 >= 15, p.segments)
        # a setting held for longer than half the night leaves no room for a change
        o, ok = sky(vega...)
        p0 = pop_plan(f, vega[2], o.ha, ones(Int, 6); observable = ok, max_changes = 0)
        pl = pop_plan(f, vega[2], o.ha, ones(Int, 6); observable = ok, max_changes = 2,
                      min_stretch = 400)
        @test pl.changes == 0 && pl.minutes == p0.minutes
        # a change whose settings are held 60 minutes each
        p60 = pop_plan(f, vega[2], o.ha, ones(Int, 6); observable = ok, max_changes = 2,
                       min_stretch = 60)
        p1 = pop_plan(f, vega[2], o.ha, ones(Int, 6); observable = ok, max_changes = 2)
        @test p60.changes >= 1 && p0.minutes[1] < p60.minutes[1] <= p1.minutes[1]
        @test segments_in_delay(p60, vega[2], ones(Int, 6))
    end

    @testset "degenerate input" begin
        o, ok = sky(vega...)
        @test isempty(pop_plan(f, vega[2], o.ha, [1, 0, 0, 0, 0, 0]; observable = ok).segments)
        @test isempty(pop_plan(f, vega[2], o.ha, ones(Int, 6); observable = falses(length(ok))).segments)
        @test_throws ArgumentError pop_plan(f, vega[2], o.ha, ones(Int, 6); observable = ok[1:end-1])
    end
end

@testset "pop_plan_targets" begin
    f = read_facility_file("CHARA")
    date = DateTime(2026, 10, 26)
    function sky(ra, dec; cfg = ones(Int, 6), d = date)
        o = night_observability(f, ra, dec, d; alt_limit = 25.0, config = cfg)
        ok = falses(length(o.ha))
        ok[intersect(o.good_alt, o.good_twilight)] .= true
        return o, ok
    end
    target(name, ra, dec; cal = false) = (; name, ra, dec, cal)
    skies(ts; kw...) = reduce(hcat, [sky(t.ra, t.dec; kw...)[2] for t in ts])
    lst = sky(0.0, 0.0)[1].lst
    cfg = ones(Int, 6)
    wrap(h) = mod(h + 12, 24) - 12

    @testset "one target is pop_plans' plan" begin
        vega = target("Vega", 279.2347, 38.7837)
        o, ok = sky(vega.ra, vega.dec; d = DateTime(2026, 6, 21))
        ps = pop_plans_targets(f, [vega], o.lst, cfg; observable = reshape(ok, :, 1))
        pp = pop_plans(f, vega.dec, o.ha, cfg; observable = ok)
        @test [p.minutes[1] for p in ps] == [p.minutes[1] for p in pp]
        @test [p.changes for p in ps] == [p.changes for p in pp]
        # each target's part is in delay at every minute, with its own telescopes and POPs
        tp = target_pop_plan(ps[end], 1, vega.ra, cfg; ha_start = unwrap_hour_angles(o.ha)[1])
        @test sum(sg -> round(Int, (sg.ha1 - sg.ha0) * 60) + 1, tp.segments) == ps[end].minutes[1]
        for sg in tp.segments
            hs = Float32.(collect(sg.ha0:1/60:sg.ha1))
            @test length(in_delay(f, vega.dec, hs, cfg, sg.pop).good_delay) == length(hs)
        end
    end

    @testset "no change: the best setting by brute force" begin
        # S1 S2 E2 W2, W2 the delay reference; up all night, so each target's whole night counts
        c4 = [1, 1, 0, 1, 0, 2]
        ts = [target("A", 330.0, 20.0), target("B", 15.0, 45.0), target("C", 60.0, 10.0)]
        obs = trues(length(lst), 3)
        p = pop_plan_targets(f, ts, lst, c4; observable = obs, max_changes = 0)
        T = length(p.lst)
        best = 0.0
        for code in 0:624
            pop = ones(Int, 6)
            for (k, t) in enumerate((1, 2, 4, 6))
                pop[t] = (code ÷ 5^(k - 1)) % 5 + 1
            end
            sc = 0.0
            for t in ts
                iv = delay_ha_intervals(f, t.dec, c4, pop)
                sc += count(l -> any(x -> x[1] <= wrap(l - t.ra / 15) <= x[2], iv), p.lst) / T
            end
            best = max(best, sc)
        end
        @test p.score ≈ best atol = 1e-9
        @test p.score ≈ sum(p.minutes) / T atol = 1e-9
        @test length(p.segments) == 1 && p.changes == 0
    end

    @testset "two targets 6 h apart: one change serves both" begin
        a, b = target("A", 330.0, 20.0), target("B", 60.0, 20.0)
        obs = skies([a, b])
        p0, p1 = pop_plans_targets(f, [a, b], lst, cfg; observable = obs, max_changes = 1)
        alone(t, j) = pop_plans(f, t.dec, sky(t.ra, t.dec)[1].ha, cfg; observable = obs[:, j],
                                max_changes = 0)[1].minutes[1]
        @test p0.minutes == [103, 149]
        @test p1.changes == 1 && p1.minutes == [alone(a, 1), alone(b, 2)] == [106, 163]
        @test p1.segments[1].pop != p1.segments[2].pop
        @test (p1.segments[2].lst0 - p1.segments[1].lst1) * 60 >= 10
        @test p1.score > p0.score
    end

    @testset "a calibrator out of reach takes the science minutes with it" begin
        sci = target("S", 330.0, 30.0)
        cal = target("C", 15.0, 30.0; cal = true)          # 3 h later
        obs = skies([sci, cal])
        free = pop_plan_targets(f, [sci, cal], lst, cfg; observable = obs, max_changes = 0)
        tied = pop_plan_targets(f, [sci, cal], lst, cfg; observable = obs, max_changes = 0,
                                calibrators = [[2], Int[]], bracket_minutes = 30)
        @test 0 < tied.minutes[1] < free.minutes[1]
        # every science minute has the calibrator within 30 minutes
        cm = findall(tied.observed[:, 2])
        @test all(t -> any(u -> abs(u - t) <= 30, cm), findall(tied.observed[:, 1]))
        @test isempty(tied.ignored)
        wide = pop_plan_targets(f, [sci, cal], lst, cfg; observable = obs, max_changes = 0,
                                calibrators = [[2], Int[]], bracket_minutes = 600)
        @test wide.minutes[1] == free.minutes[1]
    end

    @testset "calibrators that cannot be used are ignored" begin
        sci = target("S", 330.0, 30.0)
        plain = pop_plans_targets(f, [sci], lst, cfg; observable = skies([sci]))
        # below the horizon all night
        low = target("L", 330.0, -80.0; cal = true)
        # up most of the night, and in delay under no setting of the six telescopes
        pole = target("P", 0.0, 85.0; cal = true)
        @test count(skies([low])) == 0 && count(skies([pole])) > 500
        for c in (low, pole)
            ps = pop_plans_targets(f, [sci, c], lst, cfg; observable = skies([sci, c]),
                                   calibrators = [[2], Int[]])
            @test [p.minutes[1] for p in ps] == [p.minutes[1] for p in plain]
            @test all(p -> p.ignored == [2] && p.minutes[2] == 0, ps)
        end
    end

    @testset "a list of calibrators only is scored as science" begin
        ts = [target("A", 330.0, 20.0; cal = true), target("B", 60.0, 20.0; cal = true)]
        sc = [(; t..., cal = false) for t in ts]
        obs = skies(ts)
        a = pop_plans_targets(f, ts, lst, cfg; observable = obs, max_changes = 1)
        b = pop_plans_targets(f, sc, lst, cfg; observable = obs, max_changes = 1)
        @test [p.minutes for p in a] == [p.minutes for p in b]
        @test [p.score for p in a] == [p.score for p in b]
    end

    @testset "min_stretch" begin
        a, b = target("A", 330.0, 20.0), target("B", 60.0, 20.0)
        obs = skies([a, b])
        p0 = pop_plan_targets(f, [a, b], lst, cfg; observable = obs, max_changes = 0)
        p1 = pop_plan_targets(f, [a, b], lst, cfg; observable = obs, max_changes = 1, min_stretch = 60)
        pl = pop_plan_targets(f, [a, b], lst, cfg; observable = obs, max_changes = 1, min_stretch = 400)
        @test p1.changes == 1 && p1.score > p0.score
        @test pl.changes == 0 && pl.minutes == p0.minutes
    end

    @testset "degenerate input" begin
        t = target("A", 330.0, 20.0)
        @test isempty(pop_plan_targets(f, [t], lst, [1, 0, 0, 0, 0, 0]; observable = skies([t])).segments)
        @test isempty(pop_plan_targets(f, [t], lst, cfg; observable = falses(length(lst), 1)).segments)
        @test_throws ArgumentError pop_plan_targets(f, [t], lst, cfg; observable = falses(3, 1))
        # a calibrator listed for a target must be marked cal
        @test_throws ArgumentError pop_plan_targets(f, [t, t], lst, cfg; observable = skies([t, t]),
                                                    calibrators = [[2], Int[]])
    end
end
