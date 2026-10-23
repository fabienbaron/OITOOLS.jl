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

    @testset "degenerate input" begin
        o, ok = sky(vega...)
        @test isempty(pop_plan(f, vega[2], o.ha, [1, 0, 0, 0, 0, 0]; observable = ok).segments)
        @test isempty(pop_plan(f, vega[2], o.ha, ones(Int, 6); observable = falses(length(ok))).segments)
        @test_throws ArgumentError pop_plan(f, vega[2], o.ha, ones(Int, 6); observable = ok[1:end-1])
    end
end
