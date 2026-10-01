# The Gantt chart, in Makie.
#
# One rendering of `gantt_geometry`; the matplotlib `gantt_onenight` is the other. Neither
# computes the chart — both are handed it — which is what lets `test/gui/plotport.jl` assert
# that they agree rather than hope so.
#
# Static on purpose: the rectangles are not pickable. A Gantt is read, not manipulated, and
# adding hit-testing would buy interaction nobody asked for at the cost of a second input path
# through MakieArea.

"""
Point size for the numbers printed around each bar, before the screen scale is applied.

Smaller than the legend's, and deliberately so: four numbers meet at every bar end -- the time
above, the azimuth and the altitude below, and the opposite bar's three coming the other way --
so they are set in the gaps between bars rather than in open space. At the legend's size they
crowd the bar they are annotating.
"""
const GANTT_ANNOTATION_PT = 7.0

"Colours by name, so the geometry can stay renderer-independent."
const GANTT_COLORS = Dict(
    "lightgray"    => Makie.RGBAf(0.827, 0.827, 0.827, 0.75),
    "gray"         => Makie.RGBAf(0.502, 0.502, 0.502, 0.75),
    "orange"       => Makie.RGBAf(1.000, 0.647, 0.000, 1.0),
    "mediumpurple" => Makie.RGBAf(0.576, 0.439, 0.859, 1.0),
    "blue"         => Makie.RGBAf(0.000, 0.000, 1.000, 1.0),
    "green"        => Makie.RGBAf(0.000, 0.502, 0.000, 1.0),
)

_gantt_color(name) = get(GANTT_COLORS, name, Makie.RGBAf(0.5, 0.5, 0.5, 1.0))

"""
    build_gantt(figure, axis) -> NamedTuple

Create every plot the Gantt will ever need, before the window exists.

Same rule as `build_canvas`: once Qt owns the GL context, inserting a plot allocates buffers
with none bound. So the chart is a fixed set of plots whose Observables are rewritten, not a
figure that is rebuilt per target.
"""
function build_gantt(fig, ax)
    # Not zoomable. This is a night at a fixed scale, not a plot to explore: the x axis
    # IS the night and the y axis is a list of rows, so panning either only loses the
    # thing being read. Makie gives every Axis these by default.
    for it in (:scrollzoom, :dragpan, :rectanglezoom, :limitreset)
        Makie.deregister_interaction!(ax, it)
    end

    bandrects = Makie.Observable(Makie.Rect2f[])
    bandcols  = Makie.Observable(Makie.RGBAf[])
    barrects  = Makie.Observable(Makie.Rect2f[])
    barcols   = Makie.Observable(Makie.RGBAf[])
    midline   = Makie.Observable(Makie.Point2f[])

    # Start and end labels are separate plots because they anchor in opposite directions: a
    # time centred on the bar end overprints the az/alt numbers, which is what the original
    # avoids with ha="right" / ha="left".
    tspos     = Makie.Observable(Makie.Point2f[]);  tstxt   = Makie.Observable(String[])
    tepos     = Makie.Observable(Makie.Point2f[]);  tetxt   = Makie.Observable(String[])
    # Split by side, like the times: an az number centred ON the bar edge lands under the
    # vertical time already anchored there, and on a short run the two ends collide as well.
    # Anchored outward, each label leans away from the bar and away from its opposite number.
    azspos    = Makie.Observable(Makie.Point2f[]);  azstxt  = Makie.Observable(String[])
    azepos    = Makie.Observable(Makie.Point2f[]);  azetxt  = Makie.Observable(String[])
    altspos   = Makie.Observable(Makie.Point2f[]);  altstxt = Makie.Observable(String[])
    altepos   = Makie.Observable(Makie.Point2f[]);  altetxt = Makie.Observable(String[])

    legpos    = Makie.Observable(Makie.Point2f[]);  legtxt  = Makie.Observable(String[])
    legcols   = Makie.Observable(Makie.RGBAf[])

    # Local midnight, the date it turns over, and the target's culmination. Separate from the
    # red "now" line above: that one says where the clock is, these say where the NIGHT is.
    softrects    = Makie.Observable(Makie.Rect2f[])
    gridsegs     = Makie.Observable(Makie.Point2f[])
    midnightline = Makie.Observable(Makie.Point2f[])
    datepos      = Makie.Observable(Makie.Point2f[]); datetxt = Makie.Observable(String[])
    transitpos   = Makie.Observable(Makie.Point2f[])

    # Bands first, then the time rules, then bars, then text.
    #
    # The rules are drawn HERE rather than left to the axis grid, because Makie puts axis
    # decorations under plot content and the bands are opaque: the axis grid is invisible
    # exactly where the chart is darkest. Drawn between bands and bars they read over the
    # twilight and stay behind the bar, which is the order ASPRO uses.
    Makie.poly!(ax, bandrects; color = bandcols, strokewidth = 0)
    Makie.linesegments!(ax, gridsegs; color = Makie.RGBAf(1, 1, 1, 0.85),
                        linestyle = :dot, linewidth = 0.9)
    Makie.poly!(ax, barrects;  color = barcols,  strokewidth = 0)
    # Above the mount's soft ceiling: observable, but the drives complain. A WHITE wash over
    # whatever the bar's own colour is, plus a dashed edge -- so it stays the same bar, paler,
    # rather than becoming a differently coloured one. Washing rather than tinting is what lets
    # it work for every baseline colour in the detailed view.
    Makie.poly!(ax, softrects; color = (:white, 0.55), strokewidth = 1.1,
                strokecolor = :black, linestyle = :dash)
    Makie.lines!(ax, midline; color = :red, linewidth = 1.5)
    # Local midnight, under the bars rather than over them: it is a reference, not a datum.
    Makie.lines!(ax, midnightline; color = :black, linewidth = 1.2)
    # The culmination, which is the one instant on the row an observer aims for.
    Makie.scatter!(ax, transitpos; color = :yellow, marker = :diamond, markersize = 11,
                   strokewidth = 0.8, strokecolor = :black)

    sc = live_plot_scale()
    fs = Float32(9 * sc)                            # the legend
    afs = Float32(GANTT_ANNOTATION_PT * sc)         # the numbers around each bar
    # Times are rotated to vertical, as in the original: at one-minute resolution the starts
    # and ends of adjacent blocks would otherwise overprint each other.
    # Nudged outward in SCREEN pixels, not in data: the anchor stays the true bar end, which is
    # what the geometry says, while the glyphs sit clear of the bar and of the az/alt numbers —
    # the placement `ha="right"`/`ha="left"` gives in the original. A data-space offset would
    # move with the zoom and stop meaning the same thing.
    Makie.text!(ax, tspos; text = tstxt, rotation = Float32(π/2),
                align = (:center, :center), offset = (-0.75afs, 0), fontsize = afs)
    Makie.text!(ax, tepos; text = tetxt, rotation = Float32(π/2),
                align = (:center, :center), offset = (0.75afs, 0), fontsize = afs)
    # The numbers read INWARD: the start pair begins at the bar's left edge and the end pair
    # finishes at its right one, so each stays over the bar it describes. Outward, as they were,
    # two runs a few minutes apart wrote their numbers over each other in the gap between them.
    Makie.text!(ax, azspos;  text = azstxt,  align = (:left,  :bottom),
                offset = (0.2afs, 0),  fontsize = afs)
    Makie.text!(ax, azepos;  text = azetxt,  align = (:right, :bottom),
                offset = (-0.2afs, 0), fontsize = afs)
    Makie.text!(ax, altspos; text = altstxt, align = (:left,  :top),
                offset = (0.2afs, 0),  fontsize = afs)
    Makie.text!(ax, altepos; text = altetxt, align = (:right, :top),
                offset = (-0.2afs, 0), fontsize = afs)

    # The date the night turns over, set just above the axis at midnight: a Gantt that crosses
    # into the next day otherwise gives no clue which day a morning hour belongs to.
    Makie.text!(ax, datepos; text = datetxt, align = (:center, :bottom), fontsize = afs,
                color = :black)

    # A hand-drawn legend, for the reason build_canvas has one: a Makie Legend fixes its entry
    # count at construction, and this chart's varies with which constraints apply.
    Makie.scatter!(ax, legpos; color = legcols, marker = :rect, markersize = fs)
    Makie.text!(ax, legpos; text = legtxt, align = (:left, :center), fontsize = fs,
                offset = (fs, 0))

    # The target names the row it occupies, as a y tick — which is where the original puts it,
    # and it leaves the title free. y = 2 is the row the observable bar is drawn on.
    ax.ylabel = ""
    ax.xlabel = "LST (h)"
    # Dark enough that the white time rules read against it. The twilight bands sit on top and
    # step darker from here, so this is the chart's "outside the night" tone rather than paper.
    ax.backgroundcolor = Makie.RGBAf(0.88, 0.88, 0.88, 1.0)
    # The axis text follows `live_plot_scale`, as the annotations do, so the settings panel's
    # plot scale moves all of it together. Left at Makie's default the tick labels are more
    # than twice the size of the numbers they sit under.
    ax.xticklabelsize = Float32(8 * sc)
    ax.yticklabelsize = Float32(9 * sc)
    ax.xlabelsize     = Float32(9 * sc)
    ax.titlesize      = Float32(9 * sc)
    ax.yticks = ([2.0], [""])
    # No horizontal gridline. The y axis is categorical -- one row per target -- so a line
    # through the row divides nothing, and it runs straight through the rotated start/end times
    # drawn at the bar ends. The vertical grid stays: that one marks the hours, which is the
    # axis a reader actually measures against.
    ax.ygridvisible = false
    ax.yticksvisible = false
    Makie.ylims!(ax, 0, 10)

    return (; figure = fig, axis = ax, bandrects, bandcols, barrects, barcols, midline,
              softrects, gridsegs, midnightline, datepos, datetxt, transitpos,
              tspos, tstxt, tepos, tetxt, azspos, azstxt, azepos, azetxt, altspos, altstxt, altepos, altetxt,
              legpos, legtxt, legcols)
end

"""
    gantt_time_axis!(g, p::NightPlan, geo; system = :lst)

Label the x axis in `system` -- `:lst`, `:utc` or `:local` -- and rule it every 15 minutes.

The DATA stays in LST: moving the bars would mean re-deriving the whole geometry per time
system, and both renderers consume that geometry. Only the ticks move. LST and solar time run
at different rates, so a tick at a round UTC minute does NOT sit at a round LST one; each
label's position is mapped through the night's own (LST, UTC) pairs.

Ticks are hourly and labelled `HH:MM`, with unlabelled minor ticks every 15 minutes and a
dotted grid on both — the reading an observer does off this chart is "how long have I got",
which is a measurement against the time axis rather than against the bars.
"""
function gantt_time_axis!(g, p::NightPlan, geo; system::Symbol = :lst)
    ax = g.axis
    x0, x1 = geo.xlim

    # White dotted rules on a light background, which is how ASPRO draws them. A grey grid is
    # invisible against the twilight bands -- they are grey too -- and the bands are the whole
    # point of the chart, so it is the GRID that has to contrast with them rather than the
    # other way round. The axis background below is what makes white read outside the bands.
    # The rules themselves are `g.gridsegs`, drawn with the plots; the axis only carries the
    # ticks. See `build_gantt` for why they cannot be the axis's own grid.
    ax.xgridvisible       = false
    ax.xminorgridvisible  = false
    ax.xminorticksvisible = true
    ax.xminorticks        = Makie.IntervalsBetween(4)
    quarters = Float64[]

    if system === :lst
        h0 = ceil(x0); h1 = floor(x1)
        majors = collect(h0:1.0:h1)
        ax.xticks = (majors, [_gantt_hhmm(h) for h in majors])
        ax.xlabel = gantt_time_label(:lst, p)
        append!(quarters, (ceil(x0 * 4) / 4):0.25:x1)
        _gantt_rules!(g, quarters, geo)
        return g
    end

    u0 = lst_to_utc(p, x0); u1 = lst_to_utc(p, x1)
    (isfinite(u0) && isfinite(u1) && u1 > u0) || return g

    # Civil offset in hours at a given UTC hour of this night. Zero for UTC; for local it is
    # evaluated AT THE TICK, so a night spanning a daylight-saving change is right either side.
    off(u) = system === :utc ? 0.0 :
             local_utc_offset(p, Dates.DateTime(Dates.Date(p.date)) +
                                 Dates.Millisecond(round(Int, u * 3.6e6)))
    # Displayed clock back to UTC. The offset is piecewise constant, so one refinement settles
    # it everywhere except within an hour of the switch itself.
    to_utc(q) = (u = q - off(q); q - off(u))

    majors = Float64[]; labels = String[]
    # Quarter hours of the DISPLAYED clock, each mapped back onto the LST axis. Stepping in LST
    # instead would put the rules at ragged clock times, which is the opposite of the point.
    q = ceil((u0 + off(u0)) * 4) / 4
    while q <= u1 + off(u1)
        x = utc_to_lst(p, to_utc(q))
        if isfinite(x) && x0 <= x <= x1
            push!(quarters, x)
            if abs(q - round(q)) < 1e-9
                push!(majors, x); push!(labels, _gantt_hhmm(q))
            end
        end
        q += 0.25
    end

    isempty(majors) || (ax.xticks = (majors, labels))
    ax.xlabel = gantt_time_label(system, p)
    _gantt_rules!(g, quarters, geo)
    return g
end

"Vertical rules at `xs`, spanning the chart."
function _gantt_rules!(g, xs, geo)
    segs = Makie.Point2f[]
    for x in xs
        push!(segs, Makie.Point2f(x, 0), Makie.Point2f(x, geo.ymax))
    end
    g.gridsegs[] = segs
    return g
end

"Decimal hours as HH:MM, for a time axis."
function _gantt_hhmm(h)
    m = round(Int, mod(h, 24) * 60)
    return string(lpad(div(m, 60) % 24, 2, '0'), ":", lpad(mod(m, 60), 2, '0'))
end

"""
    update_gantt!(g, plan; detailed = false, time_system = :lst)

Draw one night. Allocates only the vectors handed to the Observables.
"""
function update_gantt!(g, p::NightPlan; detailed::Bool = false,
                       time_system::Symbol = :lst)
    geo = gantt_geometry(p; detailed, time_system)

    _rect(b) = Makie.Rect2f(b.x0, b.y - b.height/2, max(b.x1 - b.x0, 1e-6), b.height)

    g.bandrects[] = [_rect(b) for b in geo.bands]
    # Baselines colour by name through Explore's map, everything else by the named palette.
    # Built from the ROWS, not from the bars, so a baseline that is never in delay -- and so
    # draws no bar -- still consumes its colour and the others keep theirs.
    blmap = baseline_color_map([r[2] for r in geo.rows if occursin('-', r[2])])
    # `baseline_color_map` yields colour NAMES, as the observable plots take them; convert the
    # same way `canvas_data` does so the two end up byte-identical rather than merely similar.
    barcol(b) = haskey(GANTT_COLORS, b.color) ? GANTT_COLORS[b.color] :
                haskey(blmap, b.color)       ? Makie.RGBAf(Makie.to_color(blmap[b.color])) :
                                               _gantt_color(b.color)

    g.bandcols[]  = [_gantt_color(b.color) for b in geo.bands]
    g.barrects[]  = [_rect(b) for b in geo.bars]
    g.softrects[] = [_rect(b) for b in geo.softbars]
    g.barcols[]   = [barcol(b) for b in geo.bars]
    # The clock at the facility as the night was worked out, not local midnight. Empty when
    # `gantt_geometry` reports it outside the plotted hours: no line is the honest drawing of
    # "now is not on this chart", where a line at the edge would read as "now is right here".
    g.midline[]   = isfinite(geo.now) ?
                    [Makie.Point2f(geo.now, 0), Makie.Point2f(geo.now, geo.ymax)] :
                    Makie.Point2f[]

    for (pick, pos, txt) in ((l -> l.kind === :time && l.side === :start, g.tspos, g.tstxt),
                             (l -> l.kind === :time && l.side === :end,   g.tepos, g.tetxt),
                             (l -> l.kind === :az   && l.side === :start, g.azspos,  g.azstxt),
                             (l -> l.kind === :az   && l.side === :end,   g.azepos,  g.azetxt),
                             (l -> l.kind === :alt  && l.side === :start, g.altspos, g.altstxt),
                             (l -> l.kind === :alt  && l.side === :end,   g.altepos, g.altetxt))
        sel = filter(pick, geo.labels)
        pos[] = [Makie.Point2f(l.x, l.y) for l in sel]
        txt[] = [l.text for l in sel]
    end

    # No legend, in either view. Every row is named by its own y tick, so the legend restated
    # what the axis already said while covering the right-hand end of the chart -- exactly
    # where a target that stays up late puts its bar. The plots stay (removing them would mean
    # rebuilding the figure) and are simply fed nothing.
    named = eltype(geo.bars)[]
    x0 = geo.xlim[1] + 0.80 * (geo.xlim[2] - geo.xlim[1])   # top right, as in the original
    ytop = geo.ymax - 0.6
    g.legpos[]  = [Makie.Point2f(x0, ytop - 0.06 * geo.ymax * (k - 1)) for k in eachindex(named)]
    g.legtxt[]  = [b.label for b in named]
    g.legcols[] = [_gantt_color(b.color) for b in named]

    # One tick per row: in detailed mode those are the baselines, which is what makes the view
    # readable at all -- fifteen unlabelled bars say nothing.
    # The POPs the chart was made with. Same target, same night, different POPs gives a
    # different window, so a Gantt that does not say which is not reproducible from what it
    # shows. "delay lines not checked" is the honest label when they were not applied.
    g.axis.title = isempty(geo.subtitle) ?
                   (geo.delay_applied ? "" : "delay lines not checked") :
                   "POPs   " * geo.subtitle

    g.axis.yticks = ([r[1] for r in geo.rows], [r[2] for r in geo.rows])

    # Local midnight, and the date it turns over written just above the axis. `geo.midnight` is
    # on the same LST axis as everything else, so nothing needs converting to place it.
    g.midnightline[] = geo.xlim[1] <= geo.midnight <= geo.xlim[2] ?
                       [Makie.Point2f(geo.midnight, 0), Makie.Point2f(geo.midnight, geo.ymax)] :
                       Makie.Point2f[]
    if geo.xlim[1] <= geo.midnight <= geo.xlim[2]
        d0 = Dates.Date(p.date)
        g.datepos[] = [Makie.Point2f(geo.midnight, 0.012 * geo.ymax)]
        g.datetxt[] = [Dates.format(d0, "mm/dd") * " - " *
                       Dates.format(d0 + Dates.Day(1), "mm/dd")]
    else
        g.datepos[] = Makie.Point2f[]; g.datetxt[] = String[]
    end

    # The culmination, on the row the target occupies. y = 2 is that row, as `build_gantt` says.
    g.transitpos[] = isfinite(geo.transit) ?
                     [Makie.Point2f(geo.transit, 2.0)] : Makie.Point2f[]

    gantt_time_axis!(g, p, geo; system = time_system)
    Makie.xlims!(g.axis, geo.xlim[1], geo.xlim[2])
    Makie.ylims!(g.axis, 0, geo.ymax)
    return g
end


# ─────────────────────────────────────────────────────────────────────────────
# The delay-cart plot (`chara_plan`), in Makie
# ─────────────────────────────────────────────────────────────────────────────
#
# Same arrangement as the Gantt: `delay_plot_geometry` describes the chart and this draws it.
#
# Altitude gets a real second axis on the right, as the original's `twinx` does. Mapping it onto
# the metre axis was tried first and is worse for one reason: a curve with no scale beside it
# cannot be read. Degrees and metres share a plot here, so one of them has to say which is
# which, and a labelled axis is how.

"Distinct colours for the baselines. Beyond this many they repeat, which is honest: fifteen
baselines cannot each have a memorable colour, and the legend is what names them."
const DELAY_COLORS = [
    Makie.RGBAf(0.12, 0.47, 0.71, 1), Makie.RGBAf(1.00, 0.50, 0.05, 1),
    Makie.RGBAf(0.17, 0.63, 0.17, 1), Makie.RGBAf(0.84, 0.15, 0.16, 1),
    Makie.RGBAf(0.58, 0.40, 0.74, 1), Makie.RGBAf(0.55, 0.34, 0.29, 1),
    Makie.RGBAf(0.89, 0.47, 0.76, 1), Makie.RGBAf(0.50, 0.50, 0.50, 1),
    Makie.RGBAf(0.74, 0.74, 0.13, 1), Makie.RGBAf(0.09, 0.75, 0.81, 1),
]

"How many baselines the plot can draw. Six telescopes give fifteen, which is the practical cap."
const MAX_BASELINES = 15

"""
    build_delay_plot(figure, axis) -> NamedTuple

Every plot the delay chart can need, created before the window exists.

A fixed pool of `MAX_BASELINES` curves, hidden until used: the count changes with the array,
and creating a line per baseline on demand is what allocates GL buffers with no context bound.
"""
function build_delay_plot(fig, ax)
    # Not zoomable. This is a night at a fixed scale, not a plot to explore: the x axis
    # IS the night and the y axis is a list of rows, so panning either only loses the
    # thing being read. Makie gives every Axis these by default.
    for it in (:scrollzoom, :dragpan, :rectanglezoom, :limitreset)
        Makie.deregister_interaction!(ax, it)
    end

    curves = [Makie.Observable(Makie.Point2f[]) for _ in 1:MAX_BASELINES]
    lines  = [Makie.lines!(ax, curves[i];
                           color = DELAY_COLORS[mod1(i, length(DELAY_COLORS))],
                           visible = false, linewidth = 1.5) for i in 1:MAX_BASELINES]

    limitpts = Makie.Observable(Makie.Point2f[])
    Makie.lines!(ax, limitpts; color = (:gray, 0.6), linestyle = :dash, linewidth = 1.5)

    # The altitude axis: same cell, ticks on the right, x linked so the two never drift apart.
    ax2 = Makie.Axis(fig[1, 1]; yaxisposition = :right, ylabel = "altitude (°)",
                     ygridvisible = false, xgridvisible = false)
    Makie.hidespines!(ax2)
    Makie.hidexdecorations!(ax2)
    Makie.linkxaxes!(ax, ax2)
    Makie.ylims!(ax2, -5, 90)

    altpts = Makie.Observable(Makie.Point2f[])
    Makie.lines!(ax2, altpts; color = :black, linestyle = :dot, linewidth = 1.5)
    ellim = Makie.Observable(Makie.Point2f[])
    Makie.lines!(ax2, ellim; color = (:red, 0.6), linewidth = 1.5)

    fs = Float32(9 * live_plot_scale())
    legpos = Makie.Observable(Makie.Point2f[]); legtxt = Makie.Observable(String[])
    legcol = Makie.Observable(Makie.RGBAf[])
    Makie.scatter!(ax, legpos; color = legcol, marker = :rect, markersize = fs)
    Makie.text!(ax, legpos; text = legtxt, align = (:left, :center), fontsize = fs,
                offset = (fs, 0))

    ax.xlabel = "LST (h)"
    ax.ylabel = "delay cart position (m)"
    return (; figure = fig, axis = ax, altaxis = ax2, curves, lines, limitpts, altpts, ellim,
              legpos, legtxt, legcol)
end

"""
    update_delay_plot!(d, plan)

Draw one night's delay carts. Curves beyond this array's baseline count are hidden, not removed.
"""
function update_delay_plot!(d, p::NightPlan)
    geo = delay_plot_geometry(p)
    n = min(length(geo.curves), MAX_BASELINES)

    for i in 1:MAX_BASELINES
        if i <= n
            c = geo.curves[i]
            d.curves[i][] = [Makie.Point2f(c.x[k], c.y[k]) for k in eachindex(c.x)]
            d.lines[i].visible[] = true
        else
            d.lines[i].visible[] = false
        end
    end

    x0, x1 = geo.xlim
    # Both limit lines as one polyline with a NaN break, so the pair costs one plot not two.
    d.limitpts[] = [Makie.Point2f(x0, geo.limit), Makie.Point2f(x1, geo.limit),
                    Makie.Point2f(NaN, NaN),
                    Makie.Point2f(x0, -geo.limit), Makie.Point2f(x1, -geo.limit)]

    # Altitude in DEGREES on its own axis, so the dotted curve can actually be read off.
    lo, hi = geo.ylim
    d.altpts[] = [Makie.Point2f(geo.alt.x[k], geo.alt.y[k]) for k in eachindex(geo.alt.x)]
    d.ellim[]  = [Makie.Point2f(x0, geo.alt_limit), Makie.Point2f(x1, geo.alt_limit)]
    Makie.ylims!(d.altaxis, geo.altlim[1], geo.altlim[2])

    ytop = hi - 0.06 * (hi - lo)
    d.legpos[] = [Makie.Point2f(x0 + 0.02 * (x1 - x0), ytop - 0.055 * (hi - lo) * (k - 1))
                  for k in 1:n]
    d.legtxt[] = [geo.curves[k].name for k in 1:n]
    d.legcol[] = [DELAY_COLORS[mod1(k, length(DELAY_COLORS))] for k in 1:n]

    d.axis.title = geo.target * " — delay carts (limit ±" * string(round(geo.limit; digits = 1)) *
                   " m, elevation limit " * string(round(Int, geo.alt_limit)) * "°)"
    Makie.xlims!(d.axis, x0, x1)
    Makie.ylims!(d.axis, lo, hi)
    return d
end

# ── the radial-profile preview ───────────────────────────────────────────────
#
# Two panels of one figure: I(r) as it is written, and its Hankel transform beside it. Editing
# a profile and seeing only the brightness distribution tells you half of what you need -- the
# visibility signature is what the data constrains, and two profiles that look similar in I(r)
# can be obvious apart in V(B).
#
# Built before the window, like every other plot here: after it, allocating GL buffers with no
# context bound is the failure this whole pre-create design exists to avoid.

"Baselines the preview transform is evaluated on, in Mλ."
const PROFILE_B_MAX = 400.0
const PROFILE_NB    = 300

"""
    build_profile_plot(fig, iax, vax) -> profile panel

`I(r)` on the left axis and `V(B)` on the right, each a single line whose data is an Observable.
"""
function build_profile_plot(fig, iax, vax)
    rvals = Makie.Observable(Float32[0, 1])
    ivals = Makie.Observable(Float32[0, 0])
    bvals = Makie.Observable(Float32[0, 1])
    vvals = Makie.Observable(Float32[1, 1])

    Makie.lines!(iax, rvals, ivals; color = :black, linewidth = 2)
    Makie.lines!(vax, bvals, vvals; color = :black, linewidth = 2)
    # V = 0 is where a resolved component crosses over, and reading a null off the curve is
    # most of what this panel is for.
    Makie.hlines!(vax, [0.0]; color = (:red, 0.5), linewidth = 1)

    iax.xlabel = "r (mas)";      iax.ylabel = "I(r)"
    vax.xlabel = "B (Mλ)";       vax.ylabel = "V(B)"

    # Four ticks, not `style_axis!`'s ten. Ten is matplotlib's density and right for a figure
    # that fills a window; this pair lives in a ~200 px preview column beside the editor, where
    # ten labels overlap into an unreadable smear. Same reason the type is smaller: the panel
    # is here to show the SHAPE of the profile and where its nulls fall, and a number that
    # cannot be read is worse than one that is not drawn.
    for a in (iax, vax)
        a.xticks[] = Makie.LinearTicks(4)
        a.yticks[] = Makie.LinearTicks(4)
        a.xticklabelsize[] = 9;  a.yticklabelsize[] = 9
        a.xlabelsize[]     = 10; a.ylabelsize[]     = 10
        a.xminorticksvisible[] = false
        a.yminorticksvisible[] = false
    end
    return (; figure = fig, iaxis = iax, vaxis = vax, rvals, ivals, bvals, vvals)
end

"""
    update_profile_plot!(p, r, I, B, V)

Assignment only, as everywhere else on a live canvas: `r` and `B` in mas and Mλ.
"""
function update_profile_plot!(p, r, I, B, V)
    p.rvals[] = Float32.(r); p.ivals[] = Float32.(I)
    p.bvals[] = Float32.(B); p.vvals[] = Float32.(V)
    Makie.autolimits!(p.iaxis)
    Makie.autolimits!(p.vaxis)
    return p
end
