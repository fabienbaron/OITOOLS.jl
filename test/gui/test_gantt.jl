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
                      config = CFG, pop = POP, use_delay = true)
    geo  = gantt_geometry(plan; show_alt = true)

    mpl() = gantt_onenight("Vega", NIGHT, obs.lst, obs.lst_midnight,
                           obs.az, obs.alt, obs.good_alt, dr.good_delay;
                           good_twilight = obs.good_twilight, good_moon = obs.good_moon)

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

    @testset "the rows are the ones the chart means" begin
        # These are not arbitrary: y = 2 is the target's own row (the only y-tick), 5 is the
        # altitude window drawn across the band, 7 the moon. A port that shuffled them would
        # still "look like a Gantt" and be wrong.
        bars = Dict(b.label => b for b in geo.bars if !isempty(b.label))
        @test bars["Altitude"].y  == 5.0
        @test bars["Moon sep."].y == 7.0
        @test bars["In Delay"].y  == 2.0
        @test bars["In Delay"].height == 2.0
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

    @testset "summary mode draws the intersection, not the delay window" begin
        det = gantt_geometry(plan; show_alt = true)
        sum_ = gantt_geometry(plan; show_alt = false)
        @test any(b -> b.label == "In Delay",   det.bars)
        @test any(b -> b.label == "Observable", sum_.bars)
        # summary has no per-constraint rows at all
        @test !any(b -> b.y == 7.0, sum_.bars)
        # and its bar covers no more than the delay one, since it is an intersection
        span(bars, y) = sum(b.x1 - b.x0 for b in bars if b.y == y; init = 0.0)
        @test span(sum_.bars, 2.0) <= span(det.bars, 2.0) + 1e-9
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
        update_gantt!(g, plan; show_alt = true)

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
