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
    azpos     = Makie.Observable(Makie.Point2f[]);  aztxt   = Makie.Observable(String[])
    altpos    = Makie.Observable(Makie.Point2f[]);  alttxt  = Makie.Observable(String[])

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
    Makie.text!(ax, azpos;   text = aztxt,   align = (:center, :top),    fontsize = fs)
    Makie.text!(ax, altpos;  text = alttxt,  align = (:center, :bottom), fontsize = fs)

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
              tspos, tstxt, tepos, tetxt, azpos, aztxt, altpos, alttxt,
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
    g.bandcols[]  = [_gantt_color(b.color) for b in geo.bands]
    g.barrects[]  = [_rect(b) for b in geo.bars]
    g.barcols[]   = [_gantt_color(b.color) for b in geo.bars]
    g.midline[]   = [Makie.Point2f(geo.midnight, 0), Makie.Point2f(geo.midnight, geo.ymax)]

    for (pick, pos, txt) in ((l -> l.kind === :time && l.side === :start, g.tspos, g.tstxt),
                             (l -> l.kind === :time && l.side === :end,   g.tepos, g.tetxt),
                             (l -> l.kind === :az,                        g.azpos, g.aztxt),
                             (l -> l.kind === :alt,                       g.altpos, g.alttxt))
        sel = filter(pick, geo.labels)
        pos[] = [Makie.Point2f(l.x, l.y) for l in sel]
        txt[] = [l.text for l in sel]
    end

    # One legend entry per NAMED bar, which is the first run of each category — the same rule
    # matplotlib applies with "_nolegend_".
    named = [b for b in geo.bars if !isempty(b.label)]
    x0 = geo.xlim[1] + 0.80 * (geo.xlim[2] - geo.xlim[1])   # top right, as in the original
    ytop = geo.ymax - 0.6
    g.legpos[]  = [Makie.Point2f(x0, ytop - 0.06 * geo.ymax * (k - 1)) for k in eachindex(named)]
    g.legtxt[]  = [b.label for b in named]
    g.legcols[] = [_gantt_color(b.color) for b in named]

    # One tick per row: in detailed mode those are the baselines, which is what makes the view
    # readable at all -- fifteen unlabelled bars say nothing.
    g.axis.yticks = ([r[1] for r in geo.rows], [r[2] for r in geo.rows])
    Makie.xlims!(g.axis, geo.xlim[1], geo.xlim[2])
    Makie.ylims!(g.axis, 0, geo.ymax)
    return g
end
