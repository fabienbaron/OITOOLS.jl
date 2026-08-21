# The Gantt chart, in Makie.
#
# One rendering of `gantt_geometry`; the matplotlib `gantt_onenight` is the other. Neither
# computes the chart — both are handed it — which is what lets `test/gui/plotport.jl` assert
# that they agree rather than hope so.
#
# Static on purpose: the rectangles are not pickable. A Gantt is read, not manipulated, and
# adding hit-testing would buy interaction nobody asked for at the cost of a second input path
# through MakieArea.

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

    # Bands first, then bars on top, then text: the z-order of the original.
    Makie.poly!(ax, bandrects; color = bandcols, strokewidth = 0)
    Makie.poly!(ax, barrects;  color = barcols,  strokewidth = 0)
    Makie.lines!(ax, midline; color = :red, linewidth = 1.5)

    fs = Float32(9 * live_plot_scale())
    # Times are rotated to vertical, as in the original: at one-minute resolution the starts
    # and ends of adjacent blocks would otherwise overprint each other.
    # Nudged outward in SCREEN pixels, not in data: the anchor stays the true bar end, which is
    # what the geometry says, while the glyphs sit clear of the bar and of the az/alt numbers —
    # the placement `ha="right"`/`ha="left"` gives in the original. A data-space offset would
    # move with the zoom and stop meaning the same thing.
    Makie.text!(ax, tspos; text = tstxt, rotation = Float32(π/2),
                align = (:center, :center), offset = (-0.75fs, 0), fontsize = fs)
    Makie.text!(ax, tepos; text = tetxt, rotation = Float32(π/2),
                align = (:center, :center), offset = (0.75fs, 0), fontsize = fs)
    Makie.text!(ax, azspos;  text = azstxt,  align = (:right,  :bottom),
                offset = (-0.4fs, 0), fontsize = fs)
    Makie.text!(ax, azepos;  text = azetxt,  align = (:left,   :bottom),
                offset = (0.4fs, 0),  fontsize = fs)
    Makie.text!(ax, altspos; text = altstxt, align = (:right,  :top),
                offset = (-0.4fs, 0), fontsize = fs)
    Makie.text!(ax, altepos; text = altetxt, align = (:left,   :top),
                offset = (0.4fs, 0),  fontsize = fs)

    # A hand-drawn legend, for the reason build_canvas has one: a Makie Legend fixes its entry
    # count at construction, and this chart's varies with which constraints apply.
    Makie.scatter!(ax, legpos; color = legcols, marker = :rect, markersize = fs)
    Makie.text!(ax, legpos; text = legtxt, align = (:left, :center), fontsize = fs,
                offset = (fs, 0))

    # The target names the row it occupies, as a y tick — which is where the original puts it,
    # and it leaves the title free. y = 2 is the row the observable bar is drawn on.
    ax.ylabel = ""
    ax.xlabel = "LST (h)"
    ax.yticks = ([2.0], [""])
    # No horizontal gridline. The y axis is categorical -- one row per target -- so a line
    # through the row divides nothing, and it runs straight through the rotated start/end times
    # drawn at the bar ends. The vertical grid stays: that one marks the hours, which is the
    # axis a reader actually measures against.
    ax.ygridvisible = false
    ax.yticksvisible = false
    Makie.ylims!(ax, 0, 10)

    return (; figure = fig, axis = ax, bandrects, bandcols, barrects, barcols, midline,
              tspos, tstxt, tepos, tetxt, azspos, azstxt, azepos, azetxt, altspos, altstxt, altepos, altetxt,
              legpos, legtxt, legcols)
end

"""
    update_gantt!(g, plan; detailed = false)

Draw one night. Allocates only the vectors handed to the Observables.
"""
function update_gantt!(g, p::NightPlan; detailed::Bool = false)
    geo = gantt_geometry(p; detailed)

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
    g.barcols[]   = [barcol(b) for b in geo.bars]
    g.midline[]   = [Makie.Point2f(geo.midnight, 0), Makie.Point2f(geo.midnight, geo.ymax)]

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
