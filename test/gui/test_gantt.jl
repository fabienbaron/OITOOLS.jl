# The Gantt port: Makie against the matplotlib original.
#
# The two renderers do not compute the chart. `gantt_geometry` does, and both draw what it
# says — so what is checked here is that matplotlib really draws that geometry. If it does, the
# Makie rendering of the same geometry is faithful by construction, and a future change to
# either drawing shows up as a failing assertion rather than as a chart nobody notices moved.
#
# Positions in x cannot be compared: matplotlib works in date units and the geometry in LST
# hours. Widths, heights, rows, colours, ordering and the annotation strings all can, and
# together they are the chart.

@testset "Gantt port" begin

    VEGA_RA, VEGA_DEC = 279.2347, 38.7837
    NIGHT = DateTime(2026, 6, 21)
    CFG   = [1, 1, 1, 1, 0, 0]        # S1 S2 E1 E2
    POP   = [1, 3, 3, 4, 1, 1]        # what best_pop finds for this night

    facility = read_facility_file("CHARA")
    obs = night_observability(facility, VEGA_RA, VEGA_DEC, NIGHT)
    dr  = in_delay(facility, VEGA_DEC, obs.ha, CFG, POP)

    plan = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, NIGHT;
                      config = CFG, pop = POP, use_delay = true,
                      alt_limit = 30.0, alt_max = 90.0)   # the library's, to match the call below
    geo  = gantt_geometry(plan; detailed = false)

    # SUMMARY mode is where the two renderers draw the same chart: the intersection of every
    # constraint, on one row. Detailed is not comparable, and deliberately so — ASPRO's detailed
    # view is one row per BASELINE, which `gantt_onenight` does not draw at all.
    mpl() = gantt_onenight("Vega", NIGHT, obs.lst, obs.lst_midnight,
                           obs.az, obs.alt, obs.good_alt, dr.good_delay;
                           good_twilight = obs.good_twilight, good_moon = obs.good_moon,
                           show_alt = false)

    @testset "the delay window agrees between the two paths" begin
        # night_plan calls in_delay itself; if it passed different arguments the whole
        # comparison below would be comparing two different nights.
        @test sort(collect(plan.good_delay)) == sort(collect(dr.good_delay))
        @test plan.delay_applied
    end

    @testset "matplotlib draws the geometry, bar for bar" begin
        got = mpl_gantt_bars(mpl)
        want = geometry_bars(geo)
        @test length(got) == length(want)
        for (g, w) in zip(got, want)
            @test g.width  ≈ w.width  rtol = 1e-6
            @test g.height ≈ w.height rtol = 1e-9
            @test g.y      ≈ w.y      rtol = 1e-9
            @test lowercase(g.color) == lowercase(w.color)
        end
    end

    @testset "summary draws one row, and it is the answer" begin
        bars = Dict(b.label => b for b in geo.bars if !isempty(b.label))
        @test haskey(bars, "Observable")
        @test bars["Observable"].y == 2.0
        @test bars["Observable"].height == 2.0
        # no per-constraint rows at all in this mode
        @test !any(b -> b.y != 2.0, geo.bars)
        @test only(geo.rows)[2] == "Vega"
    end

    @testset "annotations match, string for string" begin
        got = mpl_gantt_texts(mpl)
        want = [l.text for l in geo.labels]
        @test length(got) == length(want)
        @test got == want
        # six per observing block: start and end time, and az/alt at each end
        @test length(want) == 6 * length(index_runs(sort(collect(plan.good_delay))))
    end

    @testset "one bar per contiguous run, so a gap reads as a gap" begin
        # A window broken in two must draw as two bars. Drawing it as one is the failure that
        # makes an interrupted night look bookable.
        runs = index_runs(sort(collect(plan.good_delay)))
        delaybars = [b for b in geo.bars if b.y == 2.0]
        @test length(delaybars) == length(runs)
        # and only the first of them carries the legend label
        @test count(!isempty, [b.label for b in delaybars]) == 1
    end

    @testset "detailed is the delay view: one row per baseline" begin
        # ASPRO's detailed output. The point is diagnostic: the summary says "you cannot
        # observe", and only a per-baseline row says which baseline is the reason.
        det = gantt_geometry(plan; detailed = true)
        @test length(plan.baselines) == 6            # four telescopes -> six baselines
        # a row for the answer, one for altitude, one for the moon, and one per baseline
        @test length(det.rows) == 3 + length(plan.baselines)
        ticks = [r[2] for r in det.rows]
        @test ticks[1] == "Vega"
        @test "altitude" in ticks && "moon" in ticks
        for b in plan.baselines
            @test b.name in ticks
        end
        # the chart grew to fit them, rather than overprinting inside a fixed 0..10
        @test det.ymax > 10.0
        # and the answer row is still the intersection, no wider than any single baseline
        span(g, y) = sum(b.x1 - b.x0 for b in g.bars if b.y == y; init = 0.0)
        @test span(det, 2.0) <= span(geo, 2.0) + 1e-9
    end

    @testset "a baseline can close the night on its own" begin
        # The whole reason the detailed view exists: the array's window is the intersection of
        # the baselines', so it is only ever as good as the worst one.
        @test !isempty(plan.baselines)
        worst = minimum(length(b.good) for b in plan.baselines)
        @test length(plan.good_delay) <= worst
        @test sort(collect(reduce(intersect, [b.good for b in plan.baselines]))) ==
              sort(collect(plan.good_delay))
    end

    @testset "the elevation cap is applied, and is not a formality" begin
        # Vega transits at 85.4 degrees from CHARA, so an 85-degree cap removes the top of its
        # window. A default of 90 would make the control look decorative.
        open_ = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, NIGHT; alt_max = 90.0)
        capped = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, NIGHT; alt_max = 85.0)
        @test maximum(open_.alt) > 85
        @test length(capped.good_alt) < length(open_.good_alt)
        @test capped.alt_max == 85.0
        # and the floor too
        low  = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, NIGHT; alt_limit = 25.0)
        high = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, NIGHT; alt_limit = 60.0)
        @test length(high.good_alt) < length(low.good_alt)
    end

    @testset "LST is unwrapped across midnight" begin
        # A night from 22 h to 4 h is six hours long, not minus eighteen, and every bar
        # position downstream depends on it.
        @test unwrap_lst([22.0, 23.0, 0.5, 1.5]) ≈ [22.0, 23.0, 24.5, 25.5]
        @test issorted(unwrap_lst(plan.lst))
        @test geo.midnight >= unwrap_lst(plan.lst)[1]
    end

    @testset "the Makie renderer draws what the geometry says" begin
        fig = Makie.Figure()
        ax  = Makie.Axis(fig[1, 1])
        g   = build_gantt(fig, ax)
        update_gantt!(g, plan; detailed = false)

        @test length(g.bandrects[]) == length(geo.bands)
        @test length(g.barrects[])  == length(geo.bars)
        for (r, b) in zip(g.barrects[], geo.bars)
            @test r.widths[1] ≈ Float32(b.x1 - b.x0) rtol = 1e-5
            @test r.widths[2] ≈ Float32(b.height)    rtol = 1e-6
            @test r.origin[2] ≈ Float32(b.y - b.height/2) rtol = 1e-5
        end
        # the annotations are split by kind, and together they are all of them
        @test length(g.tstxt[]) + length(g.tetxt[]) + length(g.aztxt[]) + length(g.alttxt[]) ==
              length(geo.labels)
        # start and end times anchor opposite ways, so they are separate plots
        @test length(g.tstxt[]) == length(g.tetxt[])
        # the target names its own row, as a y tick rather than a title
        @test g.axis.yticks[][2] == ["Vega"]
        # one legend entry per named bar
        @test length(g.legtxt[]) == count(b -> !isempty(b.label), geo.bars)
    end
end
