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
    "sienna"       => Makie.RGBAf(0.627, 0.322, 0.176, 1.0),
    # a multi-POP plan: segments with a telescope dropped, and the POPs row
    "dodgerblue"     => Makie.RGBAf(0.118, 0.565, 1.000, 1.0),
    "yellowgreen"    => Makie.RGBAf(0.604, 0.804, 0.196, 1.0),
    "lightsteelblue" => Makie.RGBAf(0.690, 0.769, 0.871, 1.0),
    "lightblue"      => Makie.RGBAf(0.678, 0.847, 0.902, 1.0),
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
    # A multi-POP plan: a line through every row at each POP change, the segment labels, and the
    # POPs where each setting starts. "EL" marks the stretches above the soft elevation limit.
    switchsegs   = Makie.Observable(Makie.Point2f[])
    segpos       = Makie.Observable(Makie.Point2f[]); segtxt = Makie.Observable(Any[])
    segoff       = Makie.Observable(Makie.Vec2f[])
    poppos       = Makie.Observable(Makie.Point2f[]); poptxt = Makie.Observable(Any[])
    elpos        = Makie.Observable(Makie.Point2f[]); eltxt  = Makie.Observable(String[])
    markpos      = Makie.Observable(Makie.Point2f[]); marktxt = Makie.Observable(String[])
    markepos     = Makie.Observable(Makie.Point2f[]); marketxt = Makie.Observable(String[])

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
    Makie.linesegments!(ax, switchsegs; color = (:black, 0.6), linewidth = 1.0,
                        linestyle = :dash)
    Makie.lines!(ax, midline; color = :red, linewidth = 1.5)
    # Local midnight, under the bars rather than over them: it is a reference, not a datum.
    Makie.lines!(ax, midnightline; color = Makie.RGBAf(0.0, 0.39, 0.0, 1.0), linewidth = 1.2,
                 linestyle = :dash)
    # The culmination, which is the one instant on the row an observer aims for.
    Makie.scatter!(ax, transitpos; color = :yellow, marker = :diamond, markersize = 11,
                   strokewidth = 0.8, strokecolor = :black)

    # Text sizes follow `scale`, which every draw sets to the live plot scale (`_gantt_restyle!`):
    # the chart is built before the window reports its scale.
    scale = Makie.Observable(Float32(live_plot_scale()))
    fs  = Makie.lift(s -> 9s, scale)                                    # the legend
    afs = Makie.lift(s -> Float32(GANTT_ANNOTATION_PT) * s, scale)      # the numbers around each bar
    nfs = Makie.lift(s -> Float32(GANTT_ANNOTATION_PT - 1) * s, scale)  # their az and elevation
    off(k) = Makie.lift(a -> (Float32(k) * a, 0.0f0), afs)
    # Times are rotated to vertical, as in the original: at one-minute resolution the starts
    # and ends of adjacent blocks would otherwise overprint each other.
    # Nudged INTO the bar in screen pixels, not in data: the anchor stays the true bar end, which
    # is what the geometry says, and a data-space offset would move with the zoom. Inside rather
    # than outside, so the end of one block and the start of the next do not overprint across a
    # 10-minute POP change.
    Makie.text!(ax, tspos; text = tstxt, rotation = Float32(π/2), color = :white,
                align = (:center, :center), offset = off(0.75), fontsize = afs)
    Makie.text!(ax, tepos; text = tetxt, rotation = Float32(π/2), color = :white,
                align = (:center, :center), offset = off(-0.75), fontsize = afs)
    # The numbers read INWARD: the start pair begins at the bar's left edge and the end pair
    # finishes at its right one, so each stays over the bar it describes. Outward, as they were,
    # two runs a few minutes apart wrote their numbers over each other in the gap between them.
    Makie.text!(ax, azspos;  text = azstxt,  align = (:left,  :bottom),
                offset = off(0.2),  fontsize = nfs)
    Makie.text!(ax, azepos;  text = azetxt,  align = (:right, :bottom),
                offset = off(-0.2), fontsize = nfs)
    Makie.text!(ax, altspos; text = altstxt, align = (:left,  :top),
                offset = off(0.2),  fontsize = nfs)
    Makie.text!(ax, altepos; text = altetxt, align = (:right, :top),
                offset = off(-0.2), fontsize = nfs)

    # The date the night turns over, set just above the axis at midnight: a Gantt that crosses
    # into the next day otherwise gives no clue which day a morning hour belongs to.
    Makie.text!(ax, datepos; text = datetxt, align = (:center, :bottom), fontsize = afs,
                color = :black)
    # A transition set after a POPs label is anchored at that label's start and offset in
    # pixels (`_draw_gantt!`), so it stays clear of it whatever the panel's width.
    segplot = Makie.text!(ax, segpos; text = segtxt, align = (:center, :center), fontsize = fs,
                          offset = segoff, color = :black)
    Makie.text!(ax, poppos; text = poptxt, align = (:left, :center), fontsize = fs,
                color = :black)
    Makie.text!(ax, elpos; text = eltxt, rotation = Float32(π/2), align = (:center, :center),
                fontsize = fs, color = GANTT_MARK_COLOR)
    # "POP" just after where a change starts, the time it starts and "DROP" just before.
    Makie.text!(ax, markpos; text = marktxt, rotation = Float32(π/2), align = (:center, :center),
                offset = off(0.75), fontsize = afs, color = GANTT_MARK_COLOR)
    Makie.text!(ax, markepos; text = marketxt, rotation = Float32(π/2), align = (:center, :center),
                offset = off(-0.75), fontsize = afs, color = GANTT_MARK_COLOR)

    # A hand-drawn legend, for the reason build_canvas has one: a Makie Legend fixes its entry
    # count at construction, and this chart's varies with which constraints apply.
    Makie.scatter!(ax, legpos; color = legcols, marker = :rect, markersize = fs)
    Makie.text!(ax, legpos; text = legtxt, align = (:left, :center), fontsize = fs,
                offset = Makie.lift(f -> (f, 0.0f0), fs))

    # The target names the row it occupies, as a y tick — which is where the original puts it,
    # and it leaves the title free. y = 2 is the row the observable bar is drawn on.
    ax.ylabel = ""
    ax.xlabel = "LST (h)"
    # Dark enough that the white time rules read against it. The twilight bands sit on top and
    # step darker from here, so this is the chart's "outside the night" tone rather than paper.
    ax.backgroundcolor = Makie.RGBAf(0.88, 0.88, 0.88, 1.0)
    ax.yticks = ([2.0], [""])
    # No horizontal gridline. The y axis is categorical -- one row per target -- so a line
    # through the row divides nothing, and it runs straight through the rotated start/end times
    # drawn at the bar ends. The vertical grid stays: that one marks the hours, which is the
    # axis a reader actually measures against.
    ax.ygridvisible = false
    ax.yticksvisible = false
    Makie.ylims!(ax, 0, 10)

    # The rows' vertical extent, in chart units; the view is fitted round it (`_gantt_view!`).
    yspan = Ref((0.0, 10.0))
    g = (; figure = fig, axis = ax, bandrects, bandcols, barrects, barcols, midline,
           softrects, gridsegs, midnightline, datepos, datetxt, transitpos,
           switchsegs, segpos, segtxt, segoff, segplot, poppos, poptxt, elpos, eltxt, markpos, marktxt, markepos, marketxt,
           tspos, tstxt, tepos, tetxt, azspos, azstxt, azepos, azetxt, altspos, altstxt, altepos, altetxt,
           legpos, legtxt, legcols, yspan, scale)
    _gantt_restyle!(g)
    # A resized panel shows more or less around the same rows, at the same scale.
    Makie.on(_ -> _gantt_view!(g), ax.scene.viewport)
    return g
end

"""
Pixels per chart unit, at plot scale 1. The chart keeps this scale whatever the panel's height,
so a bar and the gaps between rows keep their size: a taller panel shows more space round the
rows, not taller bars.
"""
const GANTT_PX_PER_UNIT = 60.0

"The height things that span the chart -- bands, rules, lines -- are drawn to, either way of 0."
const GANTT_FAR = 1.0e3

# Fit the view round the rows at `GANTT_PX_PER_UNIT`, centred, or to the rows when they need more.
# The text at the live plot scale. The axis text follows it as the annotations do, so the
# settings panel's plot scale moves all of it together; left at Makie's default the tick labels
# are more than twice the size of the numbers they sit under.
function _gantt_restyle!(g)
    sc = Float32(live_plot_scale())
    g.scale[] == sc || (g.scale[] = sc)
    ax = g.axis
    ax.xticklabelsize = 8sc
    ax.yticklabelsize = 9sc
    ax.xlabelsize     = 9sc
    ax.titlesize      = 9sc
    return g
end

function _gantt_view!(g)
    lo, hi = g.yspan[]
    vh = Float64(g.axis.scene.viewport[].widths[2])
    h = vh > 0 ? max(hi - lo, vh / (GANTT_PX_PER_UNIT * live_plot_scale())) : hi - lo
    c = (lo + hi) / 2
    y0, y1 = c - h / 2, c + h / 2
    Makie.ylims!(g.axis, y0, y1)
    # The date the night turns over sits at the foot of the view.
    isempty(g.datepos[]) || (g.datepos[] = [Makie.Point2f(first(g.datepos[])[1], y0 + 0.012 * (y1 - y0))])
    return g
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

    t = gantt_time_ticks(p, geo.xlim...; system)
    t === nothing && return g
    (isempty(t.majors) && system !== :lst) || (ax.xticks = (t.majors, t.labels))
    ax.xlabel = t.xlabel
    _gantt_rules!(g, t.quarters, geo)
    return g
end

"""
    gantt_time_ticks(p::NightPlan, x0, x1; system = :lst) -> (; majors, labels, quarters, xlabel)

The time axis of `gantt_time_axis!` between LST `x0` and `x1`: the hourly ticks' positions in
LST and their `HH:MM` labels in `system`, the quarter hours, and the axis label. `nothing` when
the night has no UTC mapping for `:utc` or `:local`.
"""
function gantt_time_ticks(p::NightPlan, x0, x1; system::Symbol = :lst)
    if system === :lst
        majors = collect(ceil(x0):1.0:floor(x1))
        return (; majors, labels = [_gantt_hhmm(h) for h in majors],
                  quarters = collect((ceil(x0 * 4) / 4):0.25:x1), xlabel = gantt_time_label(:lst, p))
    end

    u0 = lst_to_utc(p, x0); u1 = lst_to_utc(p, x1)
    (isfinite(u0) && isfinite(u1) && u1 > u0) || return nothing

    # Civil offset in hours at a given UTC hour of this night. Zero for UTC; for local it is
    # evaluated AT THE TICK, so a night spanning a daylight-saving change is right either side.
    off(u) = system === :utc ? 0.0 :
             local_utc_offset(p, Dates.DateTime(Dates.Date(p.date)) +
                                 Dates.Millisecond(round(Int, u * 3.6e6)))
    # Displayed clock back to UTC. The offset is piecewise constant, so one refinement settles
    # it everywhere except within an hour of the switch itself.
    to_utc(q) = (u = q - off(q); q - off(u))

    majors = Float64[]; labels = String[]; quarters = Float64[]
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
    return (; majors, labels, quarters, xlabel = gantt_time_label(system, p))
end

"Vertical rules at `xs`, spanning the chart."
function _gantt_rules!(g, xs, geo)
    segs = Makie.Point2f[]
    for x in xs
        push!(segs, Makie.Point2f(x, -GANTT_FAR), Makie.Point2f(x, GANTT_FAR))
    end
    g.gridsegs[] = segs
    return g
end

"Decimal hours as HH:MM, for a time axis."
function _gantt_hhmm(h)
    m = round(Int, mod(h, 24) * 60)
    return string(lpad(div(m, 60) % 24, 2, '0'), ":", lpad(mod(m, 60), 2, '0'))
end

# A segment's POPs with those that changed from the previous setting in dark red (`changed` is
# `segment_changes`, one entry per POP printed), or a transition label with the telescopes lost
# ("−W2") in dark red and regained ("+E1") in dark green.
const GANTT_CHANGE_COLOR = Makie.RGBAf(0.55, 0.0, 0.0, 1.0)
const GANTT_REGAIN_COLOR = Makie.RGBAf(0.0, 0.39, 0.0, 1.0)
# "EL", "POP" and the time a telescope is lost, written vertically inside the bar.
const GANTT_MARK_COLOR = Makie.RGBAf(1.0, 0.55, 0.0, 1.0)
function _segment_rich(t::AbstractString, changed::AbstractVector{Bool})
    drop, pops = startswith(t, "−") || startswith(t, "+") ? (t, "") : ("", t)
    parts = Any[]
    for (j, w) in enumerate(split(drop))
        j > 1 && push!(parts, " ")
        push!(parts, Makie.rich(String(w); color = startswith(w, "+") ? GANTT_REGAIN_COLOR :
                                                                         GANTT_CHANGE_COLOR))
    end
    for (j, v) in enumerate(split(pops))
        j > 1 && push!(parts, " ")
        push!(parts, j <= length(changed) && changed[j] ?
                     Makie.rich(String(v); color = GANTT_CHANGE_COLOR) : String(v))
    end
    return Makie.rich(parts...)
end

"""
    update_gantt!(g, plan; detailed = false, time_system = :lst)

Draw one night. Allocates only the vectors handed to the Observables.
"""
function update_gantt!(g, p::NightPlan; detailed::Bool = false,
                       time_system::Symbol = :lst)
    return _draw_gantt!(g, p, gantt_geometry(p; detailed, time_system); time_system)
end

"""
    update_gantt_stack!(g, plans, labels[, title]; time_system = :lst)

Draw several plans of one target as stacked rows, the first at the top
(`gantt_stack_geometry`).
"""
update_gantt_stack!(g, plans::AbstractVector{NightPlan}, labels, title...; time_system::Symbol = :lst) =
    _draw_gantt!(g, first(plans), gantt_stack_geometry(plans, labels, title...; time_system); time_system)

# Draw a geometry, `p` giving the night (date, clock) the axis is labelled from.
function _draw_gantt!(g, p::NightPlan, geo; time_system::Symbol = :lst)
    _gantt_restyle!(g)
    _rect(b) = Makie.Rect2f(b.x0, b.y - b.height/2, max(b.x1 - b.x0, 1e-6), b.height)

    # The twilight bands run the full height of whatever the view shows.
    g.bandrects[] = [Makie.Rect2f(b.x0, -GANTT_FAR, max(b.x1 - b.x0, 1e-6), 2GANTT_FAR) for b in geo.bands]
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
    g.switchsegs[] = reduce(vcat, ([Makie.Point2f(x, -GANTT_FAR), Makie.Point2f(x, GANTT_FAR)]
                                   for x in geo.switches); init = Makie.Point2f[])
    # A transition lifted onto the POPs labels' line goes one space after the POPs label
    # before it when it would cover it, the widths measured in the labels' font.
    font = Makie.to_font(g.segplot.fonts[], g.segplot.font[])
    pxfs = Float32(9 * g.scale[])
    width(t) = Float64(Makie.widths(Makie.text_bb(t, font, pxfs))[1])
    pxu = g.axis.scene.viewport[].widths[1] / (geo.xlim[2] - geo.xlim[1])
    segpos = Makie.Point2f[]; segoff = Makie.Vec2f[]
    for ((x, y, t, _), (px, pt)) in zip(geo.seglabels, geo.segafter)
        dx = isempty(pt) ? 0.0 : width(pt) + width(" ") + width(t) / 2
        if dx > 0 && (x - px) * pxu < dx
            push!(segpos, Makie.Point2f(px, y)); push!(segoff, Makie.Vec2f(dx, 0))
        else
            push!(segpos, Makie.Point2f(x, y)); push!(segoff, Makie.Vec2f(0, 0))
        end
    end
    g.segoff[] = segoff
    g.segpos[] = segpos
    g.segtxt[] = Any[_segment_rich(l[3], l[4]) for l in geo.seglabels]
    g.poppos[] = [Makie.Point2f(l[1], l[2]) for l in geo.poplabels]
    g.poptxt[] = Any[l[5] ? Makie.rich(Makie.rich("→ "; color = GANTT_CHANGE_COLOR),
                                       _segment_rich(l[3], l[4])) :
                            _segment_rich(l[3], l[4]) for l in geo.poplabels]
    g.markpos[]  = [Makie.Point2f(m[1], m[4]) for m in geo.marks if m[3] === :start]
    g.marktxt[]  = [m[2] for m in geo.marks if m[3] === :start]
    g.markepos[] = [Makie.Point2f(m[1], m[4]) for m in geo.marks if m[3] === :end]
    g.marketxt[] = [m[2] for m in geo.marks if m[3] === :end]
    g.elpos[]  = [Makie.Point2f(x, y) for (x, y) in geo.ellabels]
    g.eltxt[]  = fill("EL", length(geo.ellabels))
    g.barcols[]   = [barcol(b) for b in geo.bars]
    # The clock at the facility as the night was worked out, not local midnight. Empty when
    # `gantt_geometry` reports it outside the plotted hours: no line is the honest drawing of
    # "now is not on this chart", where a line at the edge would read as "now is right here".
    g.midline[]   = isfinite(geo.now) ?
                    [Makie.Point2f(geo.now, -GANTT_FAR), Makie.Point2f(geo.now, GANTT_FAR)] :
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
                   (geo.multi || geo.stack) ? geo.subtitle : "POPs   " * geo.subtitle

    g.axis.yticks = ([r[1] for r in geo.rows], [r[2] for r in geo.rows])

    # Local midnight, and the date it turns over written just above the axis. `geo.midnight` is
    # on the same LST axis as everything else, so nothing needs converting to place it.
    g.midnightline[] = geo.xlim[1] <= geo.midnight <= geo.xlim[2] ?
                       [Makie.Point2f(geo.midnight, -GANTT_FAR), Makie.Point2f(geo.midnight, GANTT_FAR)] :
                       Makie.Point2f[]
    if geo.xlim[1] <= geo.midnight <= geo.xlim[2]
        d0 = Dates.Date(p.date)
        g.datepos[] = [Makie.Point2f(geo.midnight, 0.012 * geo.ymax)]
        g.datetxt[] = [Dates.format(d0, "mm/dd") * " - " *
                       Dates.format(d0 + Dates.Day(1), "mm/dd")]
    else
        g.datepos[] = Makie.Point2f[]; g.datetxt[] = String[]
    end

    # The culmination, on each row the target occupies.
    g.transitpos[] = [Makie.Point2f(x, y) for (x, y) in geo.transits]

    gantt_time_axis!(g, p, geo; system = time_system)
    Makie.xlims!(g.axis, geo.xlim[1], geo.xlim[2])
    # The rows' extent: the detailed view's whole height, a summary's rows with their labels.
    g.yspan[] = geo.detailed ? (0.0, Float64(geo.ymax)) :
                (minimum(r[1] for r in geo.rows) - TARGET_BAR_HEIGHT / 2 - 0.35,
                 maximum(r[1] for r in geo.rows) + TARGET_BAR_HEIGHT / 2 + 0.55)
    _gantt_view!(g)
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
    # A multi-POP plan's changes, as on the Gantt.
    switchpts = Makie.Observable(Makie.Point2f[])
    Makie.linesegments!(ax, switchpts; color = (:black, 0.6), linewidth = 1.0, linestyle = :dash)

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

    # The legend, in rows under the chart: an axis of its own in the row below, unit square
    # in x and one unit per row in y, with no decorations.
    # Text sizes follow `scale`, set to the live plot scale at every draw (`_delay_restyle!`).
    scale = Makie.Observable(Float32(live_plot_scale()))
    fs = Makie.lift(s -> 9s, scale)
    legax = Makie.Axis(fig[2, 1]; height = 2fs[])
    Makie.hidedecorations!(legax)
    Makie.hidespines!(legax)
    for it in (:scrollzoom, :dragpan, :rectanglezoom, :limitreset)
        Makie.deregister_interaction!(legax, it)
    end
    legpos = Makie.Observable(Makie.Point2f[]); legtxt = Makie.Observable(String[])
    legcol = Makie.Observable(Makie.RGBAf[])
    Makie.scatter!(legax, legpos; color = legcol, marker = :rect, markersize = fs)
    Makie.text!(legax, legpos; text = legtxt, align = (:left, :center), fontsize = fs,
                offset = Makie.lift(f -> (f, 0.0f0), fs))

    ax.xlabel = "LST (h)"
    ax.ylabel = "delay cart position (m)"
    d = (; figure = fig, axis = ax, altaxis = ax2, legaxis = legax, legfs = fs, curves, lines,
           limitpts, altpts, ellim, legpos, legtxt, legcol, switchpts, scale)
    _delay_restyle!(d)
    return d
end

# The Gantt's tick size, and axis labels a step above the legend's, at the live plot scale.
function _delay_restyle!(d)
    sc = Float32(live_plot_scale())
    d.scale[] == sc || (d.scale[] = sc)
    for a in (d.axis, d.altaxis)
        a.xlabelsize = a.ylabelsize = 10.5f0 * sc
        a.xticklabelsize = a.yticklabelsize = 8sc
    end
    return d
end

"Baselines per legend row under the delay chart."
const DELAY_LEGEND_COLUMNS = 8

"""
    update_delay_plot!(d, plan; time_system = :lst)

Draw one night's delay carts. Curves beyond this array's baseline count are hidden, not removed.
The x axis is labelled in `time_system` as the Gantt's is (`gantt_time_ticks`).
"""
function update_delay_plot!(d, p::NightPlan; time_system::Symbol = :lst)
    _delay_restyle!(d)
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

    lo, hi = geo.ylim
    d.switchpts[] = reduce(vcat, ([Makie.Point2f(x, lo), Makie.Point2f(x, hi)]
                                  for x in geo.switches); init = Makie.Point2f[])
    # Altitude in DEGREES on its own axis, so the dotted curve can actually be read off.
    d.altpts[] = [Makie.Point2f(geo.alt.x[k], geo.alt.y[k]) for k in eachindex(geo.alt.x)]
    d.ellim[]  = [Makie.Point2f(x0, geo.alt_limit), Makie.Point2f(x1, geo.alt_limit)]
    Makie.ylims!(d.altaxis, geo.altlim[1], geo.altlim[2])

    ncol = max(1, min(n, DELAY_LEGEND_COLUMNS))
    nrow = max(1, cld(n, ncol))
    d.legpos[] = [Makie.Point2f((mod(k - 1, ncol) + 0.1) / ncol, div(k - 1, ncol) + 0.5)
                  for k in 1:n]
    d.legaxis.height = 1.6f0 * d.legfs[] * nrow
    Makie.limits!(d.legaxis, 0, 1, nrow, 0)
    d.legtxt[] = [geo.curves[k].name for k in 1:n]
    d.legcol[] = [DELAY_COLORS[mod1(k, length(DELAY_COLORS))] for k in 1:n]

    d.axis.title = geo.target * " — delay carts (limit ±" * string(round(geo.limit; digits = 1)) *
                   " m, elevation limit " * string(round(Int, geo.alt_limit)) * "°)"
    t = gantt_time_ticks(p, x0, x1; system = time_system)
    if t !== nothing
        isempty(t.majors) || (d.axis.xticks = (t.majors, t.labels))
        d.axis.xlabel = t.xlabel
    end
    d.axis.xminorticksvisible = true
    d.axis.xminorticks = Makie.IntervalsBetween(4)
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
