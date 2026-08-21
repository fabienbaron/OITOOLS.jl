# The live plot surface: draw once, then only ever push data.
#
# WHY THIS IS NOT JUST `draw!`.
#
# GLMakie allocates a plot's GPU buffers the moment the plot is inserted into a scene that is
# already on screen. Under QMLMakie the GL context belongs to Qt's render thread and is current
# only inside Qt's render callback, so an insertion from a QML callback -- which is where every
# user action arrives -- allocates with no context bound:
#
#     gl_renderobject: glGenBuffers returned invalid id. OpenGL Context active?
#
# It does not always bite: with nothing else competing for the context the redraw often gets a
# usable one. Opening a ComboBox popup or a FileDialog gives that surface its own context, and
# from then on every insertion fails. That is why "the first plot works and nothing after it
# does", and why a scripted test that never opens a dialog passes on a machine where the GUI is
# plainly broken. test/gui_click.sh reproduces it deliberately.
#
# There is no fix from the Julia side of that boundary: QML.jl exposes no QOpenGLContext, and
# QMLMakie's `native_switch_context!` for a QMLWindow is a deliberate no-op because GL work is
# only ever meant to happen inside the render callback.
#
# So: create every plot ONCE, before the window exists, and afterwards only assign to
# Observables. Updating an Observable uploads at the next render, on the render thread, with a
# valid context -- which is what QMLMakie's own Lorenz example does. Nothing here inserts a
# plot, a Legend or a Colorbar after `loadqml`.
#
# The offline builders (`uv_figure`, `observable_figure`, and `draw!` underneath them) are
# unchanged and still create plots freely: they build a fresh Figure that is not on screen, so
# the constraint does not apply, and they remain what the port harness checks against. The
# tests in test/runtests.jl assert that this file and `draw!` produce the same points and the
# same colours, because two implementations that drift is the failure this project already had
# once.

"""
Pre-built plot surface for the shell.

Every field is either a plot object created before display or an Observable feeding one. The
only thing that ever happens after the window opens is assignment.
"""
struct LiveCanvas
    figure      :: Any
    axis        :: Any
    # scatter
    points      :: Makie.Observable{Vector{Makie.Point2f}}
    colors      :: Makie.Observable{Vector{Makie.RGBAf}}
    markersize  :: Makie.Observable{Float32}
    scatterplot :: Any
    # error bars, empty for uv coverage
    errpoints   :: Makie.Observable{Vector{Makie.Point2f}}
    errlow      :: Makie.Observable{Vector{Float32}}
    errhigh     :: Makie.Observable{Vector{Float32}}
    errplot     :: Any
    # legend, drawn by hand into its own decoration-free axis
    legaxis     :: Any
    legmarks    :: Makie.Observable{Vector{Makie.Point2f}}
    legcolors   :: Makie.Observable{Vector{Makie.RGBAf}}
    legtextpos  :: Makie.Observable{Vector{Makie.Point2f}}
    legtext     :: Makie.Observable{Vector{String}}
    legfontsize :: Makie.Observable{Float32}
    # colorbar for the continuous colour modes
    cbarlimits  :: Makie.Observable{Tuple{Float32,Float32}}
    cbarlabel   :: Makie.Observable{String}
    colorbar    :: Any
    # reconstructed image, drawn by the Image perspective
    imagedata   :: Makie.Observable{Matrix{Float32}}
    imagex      :: Makie.Observable{Vector{Float32}}
    imagey      :: Makie.Observable{Vector{Float32}}
    imageplot   :: Any
    # Colormap for the image and its colorbar, so both move together.
    colormap    :: Makie.Observable{Any}
    # Colormap the colorbar draws. The scatter's continuous modes carry oiplot's ramps while
    # the image carries whatever the colormap buttons last chose, so the two cannot share one.
    cbarmap     :: Makie.Observable{Any}
    # The per-group panel grid: one axis per baseline / triplet / station.
    panels      :: Any
    # Width and height of the view `autolimits!` chose for the current data. The zoom clamp
    # measures against this, so "10x out" means ten times the data extent rather than ten
    # times whatever the user last left on screen.
    homespan    :: Base.RefValue{Tuple{Float64,Float64}}
end

"""
    build_canvas(figure, axis) -> LiveCanvas

Create every plot the shell will ever need. **Call before the window is shown.**
"""
function build_canvas(fig, ax)
    points   = Makie.Observable(Makie.Point2f[])
    colors   = Makie.Observable(Makie.RGBAf[])
    msize    = Makie.Observable(6.0f0)
    errpts   = Makie.Observable(Makie.Point2f[])
    errlo    = Makie.Observable(Float32[])
    errhi    = Makie.Observable(Float32[])

    # Error bars first so the markers sit on top, matching oiplot's z-order.
    errplot = Makie.errorbars!(ax, errpts, errlo, errhi; color = OIPLOT_ECOLOR)
    sc      = Makie.scatter!(ax, points; color = colors, marker = OIPLOT_MARKER,
                             markersize = msize)

    # The legend is drawn by hand rather than with Makie.Legend.
    #
    # A Legend block cannot be resized in place: its entry count is fixed at construction, and
    # changing it means assigning `entrygroups`, which rebuilds the block's internal plots --
    # an allocation, which is the one thing forbidden here. A fixed pool of 28 entries was
    # tried and is worse: blank entries still reserve their slot, so the legend overflowed the
    # figure and clipped the real labels.
    #
    # Two plots in a decoration-free axis have neither problem. The number of entries is just
    # the length of an Observable, so it can be anything from zero upwards, and zero draws
    # nothing at all.
    legax = Makie.Axis(fig[2, 1]; limits = (0, 1, 0, 1))
    Makie.hidedecorations!(legax)
    Makie.hidespines!(legax)
    legmarks = Makie.Observable(Makie.Point2f[])
    legcols  = Makie.Observable(Makie.RGBAf[])
    legtpos  = Makie.Observable(Makie.Point2f[])
    legtxt   = Makie.Observable(String[])
    legfs    = Makie.Observable(Float32(UVPLOT_LEGEND_SIZE * live_plot_scale()))
    legms    = Makie.lift(f -> f * LEGEND_MARKER_RATIO, legfs)
    Makie.scatter!(legax, legmarks; color = legcols, marker = OIPLOT_MARKER, markersize = legms)
    Makie.text!(legax, legtpos; text = legtxt, align = (:left, :center), fontsize = legfs)

    cblim = Makie.Observable((0.0f0, 1.0f0))
    # Seeded with the label it will eventually carry, then blanked by `_show_colorbar!` at the
    # end of this function. This is glyph-atlas pre-warming, and it is the same argument as
    # everything else built here: Makie keeps one texture atlas, and a glyph first met AFTER
    # the window exists has to extend it inside Qt's GL context. Measured, `μ` came out with
    # strokes of neighbouring atlas cells through it while `λ` -- already drawn by the axis
    # labels -- was clean. Meeting the character up front is what avoids the extension.
    cblab = Makie.Observable(WAVELENGTH_LABEL)
    cmap  = Makie.Observable{Any}(LIVE_COLORMAP)
    cbmap = Makie.Observable{Any}(WAV_COLORMAP)
    cbar  = Makie.Colorbar(fig[1, 2]; colormap = cbmap, limits = cblim, label = cblab)

    # The reconstructed image, created here and hidden, like everything else on this canvas.
    # Building a plot after the window exists allocates GL buffers with no context bound, which
    # is the failure this whole pre-create design is here to avoid.
    imdata = Makie.Observable(zeros(Float32, 2, 2))
    imx    = Makie.Observable(Float32[0, 1])
    imy    = Makie.Observable(Float32[0, 1])
    implot = Makie.heatmap!(ax, imx, imy, imdata; colormap = cmap, visible = false)

    panels = _build_panels(fig)

    canvas = LiveCanvas(fig, ax, points, colors, msize, sc,
                        errpts, errlo, errhi, errplot,
                        legax, legmarks, legcols, legtpos, legtxt, legfs,
                        cblim, cblab, cbar,
                        imdata, imx, imy, implot, cmap, cbmap, panels,
                        Base.RefValue((1.0, 1.0)))
    _install_zoom!(canvas)

    # The legend axis is a drawing surface, not a plot. Left interactive it responds to the
    # wheel and to drags like any other axis, and since it holds no data there is nothing on
    # screen to show that it moved -- the labels simply drift out of the visible unit square.
    for i in (:scrollzoom, :dragpan, :rectanglezoom, :limitreset)
        Makie.deregister_interaction!(legax, i)
    end
    set_legend!(canvas, Pair{String,Makie.RGBAf}[])
    _show_colorbar!(canvas, false)
    show_panels!(canvas, false)
    return canvas
end

"""
    show_image!(canvas, img, pixsize; label = "")

Put a reconstructed image on the canvas, in place of the scatter view.

Axes are in milliarcseconds and East is to the LEFT, which is the convention every figure in
this package uses and the opposite of a bare `heatmap` — an image drawn without the flip is
mirrored, and a mirrored reconstruction is not obviously wrong to the eye.

Nothing is created here. The heatmap already exists; this sets its data and swaps which plots
are visible.
"""
function show_image!(c::LiveCanvas, img::AbstractMatrix, pixsize::Real; label::AbstractString = "")
    nx = size(img, 1)
    half = nx * pixsize / 2
    # Descending x: East left. Makie takes the coordinate vectors as cell centres.
    c.imagex[] = Float32.(range(half, -half; length = nx))
    c.imagey[] = Float32.(range(-half, half; length = nx))
    c.imagedata[] = Float32.(img)

    c.points[] = Makie.Point2f[]
    c.errpoints[] = Makie.Point2f[]
    c.errlow[] = Float32[]; c.errhigh[] = Float32[]
    set_legend!(c, Pair{String,Makie.RGBAf}[])

    c.scatterplot.visible[] = false
    c.errplot.visible[] = false
    c.imageplot.visible[] = true

    mx = maximum(img)
    c.cbarlimits[] = (0.0f0, Float32(mx > 0 ? mx : 1))
    c.cbarlabel[] = isempty(label) ? "flux / pixel" : String(label)
    c.cbarmap[] = c.colormap[]
    _show_colorbar!(c, true)

    c.axis.xlabel = "α (mas)"
    c.axis.ylabel = "δ (mas)"
    Makie.limits!(c.axis, half, -half, -half, half)
    c.axis.aspect = 1
    # The zoom bounds are multiples of the home view, and for an image the home view is this
    # one -- not whatever scatter plot the canvas last held. Without this, zooming a
    # reconstruction is bounded against the uv coverage that happened to precede it.
    c.homespan[] = (2 * Float64(half), 2 * Float64(half))
    return c
end

"Return the canvas to the scatter view after an image has been shown."
function hide_image!(c::LiveCanvas)
    c.imageplot.visible[] = false
    c.scatterplot.visible[] = true
    c.errplot.visible[] = true
    c.axis.aspect = nothing
    return c
end

"""
Wavelength ramp for the `wav` colour option: matplotlib's `gist_rainbow`, which runs red at
short λ through green and blue to magenta at long λ. A spectral ramp for a spectral axis, so
the colour of a point reads as its wavelength rather than as a position on an arbitrary scale.

Resampled to 256 steps because `ramp_colors` bins onto whatever ramp it is given, and Makie's
`gist_rainbow` carries only 101 colours.
"""
const WAV_COLORMAP = Makie.resample_cmap(:gist_rainbow, 256)

"MJD ramp, matching `oiplot.jl`'s `cmap=\"plasma\"`."
const MJD_COLORMAP = :plasma

"Colormap the image starts with, and the fallback ramp for `ramp_colors`."
const LIVE_COLORMAP = :viridis

"""
    live_plot_scale() -> Float64

How much larger the on-screen plot is than the publication figure.

oiplot's 12 pt labels and 6 pt markers are right for a paper column and too small for a
window being read from a normal sitting distance -- and the port harness pins those numbers,
so they cannot simply be raised. This multiplies them for the live canvas only; `draw!` and
the offline builders are untouched and still match oiplot exactly.

Override with `OITOOLSGUI_PLOT_SCALE=2.0` (or 1.0 to get the publication sizes on screen).
"""
function live_plot_scale()
    # A value dialled in the settings panel wins over everything: it is the most recent thing
    # the user said, and the other two are a startup variable and a guess from the screen.
    PLOT_SCALE_USER[] > 0 && return PLOT_SCALE_USER[]
    v = tryparse(Float64, get(ENV, "OITOOLSGUI_PLOT_SCALE", ""))
    (v === nothing || !isfinite(v) || v <= 0) || return clamp(v, 0.5, 5.0)
    # From the screen's PHYSICAL density, on its own anchors -- not from the UI scale.
    #
    # Judged by eye: 1.19 on a 92.6 dpi desktop and 1.7 on a ~189 dpi laptop. Keeping this
    # separate from the chrome factor is what lets either be retuned without disturbing the
    # other, which tying them together did.
    dpi = SCREEN_DPI_LIVE[]
    dpi > 0 || return clamp(PLOT_SCALE_AT_REFERENCE * UI_SCALE_LIVE[] / UI_SCALE_REFERENCE,
                            0.5, 5.0)
    return clamp(PLOT_SCALE_AT_REF_DPI * sqrt(dpi / PLOT_REF_DPI), 0.5, 5.0)
end

"Colorbar label for wavelength colouring. Also the glyph-atlas seed -- see `build_canvas`."
const WAVELENGTH_LABEL = "λ (μm)"

"""
Legend swatch size, as a fraction of the legend font size.

A swatch stands for a data point, so it should not tower over one. At 1.1 it did: measured at
ui scale 0.875 the swatch came out 13.1 px against a 7.14 px marker on the plot, nearly twice
the thing it labels. 0.8 keeps it a little larger than the data point, which it needs to be to
read as a colour at this size, without dominating the row.
"""
const LEGEND_MARKER_RATIO = 0.8f0

"""
Width of one legend character as a fraction of the axis, at `UVPLOT_LEGEND_SIZE`.

Measured off a rendered legend: seven characters of a baseline name span 0.047 of the axis.
Only `set_legend!`'s centring uses it.
"""
const CHAR_W = 0.0067

# Neither Legend nor Colorbar has a `visible` attribute, so they are hidden by collapsing the
# row or column they occupy. Deleting and re-adding them is not an option: recreating a block
# allocates GPU buffers.
#
# No try/catch around these assignments. A failure to hide leaves an empty legend frame and a
# meaningless 0..1 colorbar on screen, which a swallowed exception would make invisible.
"""
Zoom per wheel detent, and how far the view may travel from the data.

Both spans are multiples of the one `autolimits!` chose, so they are absolute stops measured
against the DATA rather than against wherever the user last left the view: zooming in stops at
4x magnification (a quarter of the data extent on screen) and zooming out at 3x the extent.
Scrolling past either does nothing, which is the intended feel -- the view cannot be lost.

The ceiling is not only a convenience. Without one a few wheel clicks reach 1e49, where tick
spacing falls below what a Float64 can represent, every tick lands on the same value, and
Makie throws `AssertionError: vmin != vmax` out of its tick locator.
"""
const ZOOM_PER_DETENT = 1.1
const ZOOM_MIN_SPAN   = 1 / 4          # zoom in  ⇒ at most 4x magnification
const ZOOM_MAX_SPAN   = 3.0            # zoom out ⇒ at most 3x the data extent

"One wheel detent, in the scroll units Qt delivers for a discrete mouse wheel."
const WHEEL_DETENT = 120.0

"""
    _install_zoom!(canvas)

Replace the axis's stock scroll zoom with one that normalises the event and bounds the result.

Makie's `ScrollZoom` takes the scroll value as a count of steps and raises `(1 - speed)` to it,
which assumes something delivers small numbers. QMLMakie does not, and not by its own choice:
it scales `angleDelta` by `wheelfactor` ONLY when `pixelDelta` is zero, on the reading that a
non-zero `pixelDelta` means fine-grained pixel scrolling that needs no scaling. Under libinput
a discrete wheel reports BOTH as ±120, so the guard is false, 120 arrives unscaled, and
`0.9^120` zooms out by 270000x per click. Measured on this hardware: eight clicks reach 1e43.

A touchpad sends many small deltas instead, which is why the same code behaves on a laptop and
runs away on a desktop mouse -- so the fix has to normalise the EVENT, not lower the speed. A
speed small enough to tame 120 units per detent would make a touchpad barely zoom at all.

Dividing by one detent gives both devices the same meaning: a wheel click is one step, and a
touchpad's small deltas are fractions of a step that accumulate as the finger moves.
"""
function _install_zoom!(c::LiveCanvas)
    ax = c.axis
    Makie.deregister_interaction!(ax, :scrollzoom)
    Makie.register_interaction!(ax, :scrollzoom) do event::Makie.ScrollEvent, axis
        # Anchor on the pointer, so whatever is under the cursor stays under it.
        mp = Makie.mouseposition(axis.scene)
        zoom_step!(c, Float64(event.y) / WHEEL_DETENT; at = (Float64(mp[1]), Float64(mp[2])))
        return Makie.Consume(true)
    end
    return nothing
end

"""
    zoom_step!(canvas, steps; at = nothing) -> Bool

Zoom by `steps` wheel detents about the data point `at`, honouring the bounds. Positive steps
zoom in. Returns whether the view actually moved, which is `false` once a stop is reached.

Separate from the interaction that calls it so the bounds can be tested without a mouse: a
synthetic scroll event is not a faithful stand-in for a wheel, and the arithmetic is the part
worth pinning down.
"""
function zoom_step!(c::LiveCanvas, steps::Real; at = nothing)
    steps == 0 && return false
    ax = c.axis
    fl = ax.finallimits[]
    ox, oy = Float64(fl.origin[1]), Float64(fl.origin[2])
    wx, wy = Float64(fl.widths[1]), Float64(fl.widths[2])
    (isfinite(wx) && isfinite(wy) && wx > 0 && wy > 0) || return false

    hx, hy = c.homespan[]
    (hx > 0 && hy > 0) || return false
    # One factor for both axes. Scaling them differently would quietly change the aspect
    # ratio, and uv coverage is drawn with DataAspect, where that is a lie about the sky.
    factor = ZOOM_PER_DETENT^(-steps)
    nwx, nwy = wx * factor, wy * factor

    # Refuse the whole step rather than clamping it to the bound: overzooming does NOTHING.
    #
    # Clamping was tried and is worse. The bound cannot be landed on exactly -- `finallimits`
    # is a Rect2f, so the span round-trips through Float32, and under DataAspect the x and y
    # spans finish straddling it by about 1e-7 in opposite directions (measured) -- so a clamp
    # goes on requesting invisible corrections for ever and the view never comes to rest.
    # Refusing stops within one detent of the bound, which is 10% and invisible.
    if nwx < ZOOM_MIN_SPAN * hx || nwy < ZOOM_MIN_SPAN * hy ||
       nwx > ZOOM_MAX_SPAN * hx || nwy > ZOOM_MAX_SPAN * hy
        return false
    end

    fx, fy = if at === nothing
        (0.5, 0.5)
    else
        (clamp((Float64(at[1]) - ox) / wx, 0.0, 1.0),
         clamp((Float64(at[2]) - oy) / wy, 0.0, 1.0))
    end
    # `limits!`, not an assignment to `finallimits`. While `ax.limits` is still automatic,
    # Makie recomputes the view from the plots and a directly-assigned rectangle is discarded:
    # measured, zooming IN crept a little and then snapped back to the data extent on the next
    # step, over and over, while zooming out happened to survive because autolimits only ever
    # expand. `limits!` sets an explicit view, and `update_canvas!` calls `autolimits!` on
    # every redraw, so nothing stays pinned once the data changes.
    x0 = ox + fx * (wx - nwx)
    y0 = oy + fy * (wy - nwy)
    # Emit each pair in the axis's own direction. `finallimits` always reports positive
    # widths, so a reversed axis is a FLAG rather than an ordering, and handing `limits!` an
    # ascending pair silently clears it -- which mirrored the reconstruction east-for-west the
    # first time anyone zoomed one, since `show_image!` reverses x to put East on the left.
    xlo, xhi = ax.xreversed[] ? (x0 + nwx, x0) : (x0, x0 + nwx)
    ylo, yhi = ax.yreversed[] ? (y0 + nwy, y0) : (y0, y0 + nwy)
    Makie.limits!(ax, xlo, xhi, ylo, yhi)
    return true
end

"""
Colormaps offered for the reconstructed image, in the order the buttons show them.

Named as matplotlib names them, because that is what `imdisp` takes: whatever reads best here
can be copied straight into `imdisp(image; colormap = "gist_heat")` and give the same picture.
`gist_heat` is `imdisp`'s own default.

Each has its `_r` reverse. Which way round works is not a matter of taste — it depends on
whether the figure ends up on a screen or on a white page, and faint structure that is obvious
one way can vanish the other.
"""
const IMAGE_COLORMAPS = ["viridis"      => :viridis,
                         "gist_earth"   => :gist_earth,
                         "gist_earth_r" => Makie.Reverse(:gist_earth),
                         "gist_gray"    => :gist_gray,
                         "gist_gray_r"  => Makie.Reverse(:gist_gray),
                         "gist_heat"    => :gist_heat,
                         "gist_heat_r"  => Makie.Reverse(:gist_heat)]

"""
    set_colormap!(canvas, name) -> Bool

Switch the image colormap by its button label. Unknown names are ignored rather than thrown:
this arrives from QML, and a typo there should not take the window down.
"""
function set_colormap!(c::LiveCanvas, name::AbstractString)
    i = findfirst(p -> first(p) == name, IMAGE_COLORMAPS)
    i === nothing && return false
    c.colormap[] = last(IMAGE_COLORMAPS[i])
    c.cbarmap[]  = c.colormap[]
    return true
end

"Button labels for the image colormaps, newline-separated, for QML to build a row from."
image_colormap_names() = join(first.(IMAGE_COLORMAPS), "\n")

"""
Most panels the per-group view will draw.

21 covers a six-telescope closure set (20 triplets) and its 15 baselines whole. Beyond that the
panels are too small to read anyway, and `update_panels!` says how many it dropped rather than
quietly showing a subset.

Chosen for legibility, not cost: measured, a grid of 25 axes costs the same 130 ms at
construction as a grid of 15, because the time is first-figure overhead rather than per-axis
work -- against a launch already several seconds long.
"""
const MAX_PANELS = 21

"""
    _build_panels(figure) -> NamedTuple

Every panel the view can ever need, created before the window exists.

Same rule as the rest of this file: once Qt owns the GL context, adding a plot allocates buffers
with none bound. So the pool is fixed and redraws only assign to it -- axes are REPOSITIONED
into the shape the data wants, which Makie allows, rather than created to fit it.
"""
function _build_panels(fig)
    gl = Makie.GridLayout(fig[3, 1])
    axes = Any[]; pts = Makie.Observable{Vector{Makie.Point2f}}[]
    cols = Makie.Observable{Makie.RGBAf}[]
    epts = Makie.Observable{Vector{Makie.Point2f}}[]
    elo  = Makie.Observable{Vector{Float32}}[]
    ehi  = Makie.Observable{Vector{Float32}}[]
    for k in 1:MAX_PANELS
        a = Makie.Axis(gl[fld(k - 1, 5) + 1, mod(k - 1, 5) + 1])
        p = Makie.Observable(Makie.Point2f[])
        c = Makie.Observable(Makie.RGBAf(0, 0, 0, 1))
        ep = Makie.Observable(Makie.Point2f[])
        el = Makie.Observable(Float32[]); eh = Makie.Observable(Float32[])
        Makie.errorbars!(a, ep, el, eh; color = OIPLOT_ECOLOR)
        Makie.scatter!(a, p; color = c, marker = OIPLOT_MARKER, markersize = 5)
        push!(axes, a); push!(pts, p); push!(cols, c)
        push!(epts, ep); push!(elo, el); push!(ehi, eh)
    end
    return (; grid = gl, axes, points = pts, colors = cols,
              errpoints = epts, errlow = elo, errhigh = ehi)
end

"An Axis has no `visible`; like every Block it draws into a Scene, and that Scene has one."
_set_axis_visible!(ax, on) = (getfield(ax, :blockscene).visible[] = on; nothing)

"""
    show_panels!(canvas, on)

Swap between the single plot and the panel grid. Rows are collapsed rather than the contents
destroyed, the same way the legend and colorbar are hidden.
"""
function show_panels!(c::LiveCanvas, on::Bool)
    _set_axis_visible!(c.axis, !on)
    Makie.rowsize!(c.figure.layout, 1, on ? Makie.Fixed(0) : Makie.Auto())
    Makie.rowsize!(c.figure.layout, 3, on ? Makie.Auto() : Makie.Fixed(0))
    on && set_legend!(c, Pair{String,Makie.RGBAf}[])
    on && _show_colorbar!(c, false)
    for a in c.panels.axes
        _set_axis_visible!(a, false)
    end
    return nothing
end

"""
    update_panels!(canvas, d, kind) -> Int

Fill the grid with one panel per group and return how many were drawn.

Independent y ranges: each panel autoscales to its own group. A shared range would make the
panels comparable, but it also flattens the one thing this view is for -- the structure within
a single baseline's spectrum, which sits well inside the spread across all of them.
"""
function update_panels!(c::LiveCanvas, d, kind::Symbol)
    pd = panel_data(d, kind)
    n  = min(length(pd.panels), MAX_PANELS)
    sc = live_plot_scale()

    # A near-square grid sized to what is actually there: six baselines in a 3x2 read far
    # better than six scattered along a row of five.
    ncol = max(1, ceil(Int, sqrt(n)))
    nrow = max(1, cld(n, ncol))

    for k in 1:MAX_PANELS
        a = c.panels.axes[k]
        if k > n
            _set_axis_visible!(a, false)
            c.panels.points[k][] = Makie.Point2f[]
            c.panels.errpoints[k][] = Makie.Point2f[]
            continue
        end
        g = pd.panels[k]
        c.panels.grid[fld(k - 1, ncol) + 1, mod(k - 1, ncol) + 1] = a
        _set_axis_visible!(a, true)
        style_axis!(a; scale = sc * 0.75)      # smaller type: many panels, little room each
        # Three ticks, not the ten `style_axis!` sets to match matplotlib's density on a full
        # figure. At panel size those ten overprint each other into a grey smear and the axis
        # stops being readable at all -- which is worse than a coarse scale.
        a.xticks[] = Makie.LinearTicks(3)
        a.yticks[] = Makie.LinearTicks(3)
        a.title = g.name
        a.titlesize = OIPLOT_LABELSIZE * sc * 0.7
        # Labels only on the edges, or every panel spends its width on repeated axis text.
        a.xlabel = (fld(k - 1, ncol) + 1 == nrow) ? pd.xlabel : ""
        a.ylabel = (mod(k - 1, ncol) == 0)        ? pd.ylabel : ""
        c.panels.colors[k][]    = g.color
        c.panels.points[k][]    = Makie.Point2f.(g.x, g.y)
        c.panels.errpoints[k][] = Makie.Point2f.(g.x, g.y)
        c.panels.errlow[k][]    = Float32.(g.err)
        c.panels.errhigh[k][]   = Float32.(g.err)
        Makie.autolimits!(a)
    end
    for r in 1:5, cc in 1:5
        # Collapse the rows and columns the layout no longer uses.
        r <= nrow || Makie.rowsize!(c.panels.grid, r, Makie.Fixed(0))
        cc <= ncol || Makie.colsize!(c.panels.grid, cc, Makie.Fixed(0))
    end
    for r in 1:nrow;  Makie.rowsize!(c.panels.grid, r, Makie.Auto()); end
    for cc in 1:ncol; Makie.colsize!(c.panels.grid, cc, Makie.Auto()); end
    return n
end

"""
    set_legend!(canvas, entries; ncol = UVPLOT_LEGEND_NCOL)

Lay `label => colour` pairs out in columns and push them to the legend axis. An empty vector
clears it and collapses the row, so there is no empty frame when a colour mode has no groups.
"""
function set_legend!(c::LiveCanvas, entries::AbstractVector; ncol::Integer = UVPLOT_LEGEND_NCOL)
    n = length(entries)
    if n == 0
        c.legmarks[]   = Makie.Point2f[]
        c.legcolors[]  = Makie.RGBAf[]
        c.legtextpos[] = Makie.Point2f[]
        c.legtext[]    = String[]
        Makie.rowsize!(c.figure.layout, 2, Makie.Fixed(0))
        return nothing
    end
    rows  = cld(n, max(1, min(Int(ncol), n)))
    # Columns ACTUALLY filled, which is not always the column budget. Six entries in a
    # five-column budget need two rows, and two rows hold six entries in three columns, not
    # five. Centring against the budget would then reserve two empty columns on the right and
    # push the whole legend left of the plot it belongs to.
    cols  = cld(n, rows)
    # Column pitch from the longest label, not 1/cols.
    #
    # Spreading the columns evenly across the whole axis made the legend as wide as the window
    # whatever it contained, so five short baseline names sat in five widely separated columns
    # with the plot squeezed above them. Sizing the pitch to the text and centring the block
    # keeps it compact. The constants are in axis fractions and approximate a character's
    # width; they only need to be close, because the block is centred either way.
    maxlen = maximum(length(first(e)) for e in entries)
    pitch  = min(1.0 / cols, 0.045 + 0.0135 * maxlen)
    # The block ends at the last column's TEXT, not at the last column's pitch: a marker plus
    # its label is narrower than the pitch, which carries the gap to the next column too.
    # CHAR_W is a measured glyph width at the legend font size, and is deliberately not the
    # 0.0135 in `pitch` — that one sets a comfortable gap between columns, this one has to be
    # the real width or the block is centred against space the text does not occupy.
    width  = (cols - 1) * pitch + 0.022 + CHAR_W * maxlen
    x0     = max(0.0, (1.0 - width) / 2)
    marks = Vector{Makie.Point2f}(undef, n)
    cs    = Vector{Makie.RGBAf}(undef, n)
    tpos  = Vector{Makie.Point2f}(undef, n)
    txt   = Vector{String}(undef, n)
    for (i, e) in enumerate(entries)
        col = (i - 1) ÷ rows          # fill down each column, as oiplot's ncol legend does
        row = (i - 1) % rows
        x = x0 + col * pitch
        y = 1 - (row + 0.5) / rows
        marks[i] = Makie.Point2f(x, y)
        tpos[i]  = Makie.Point2f(x + 0.022, y)
        cs[i]    = last(e)
        txt[i]   = first(e)
    end
    c.legmarks[]   = marks
    c.legcolors[]  = cs
    c.legtextpos[] = tpos
    c.legtext[]    = txt
    # One row of legend per ~14 points of height, so many baselines get room without a
    # dataset with three of them wasting half the window.
    Makie.rowsize!(c.figure.layout, 2, Makie.Fixed((14.0 * rows + 8) * live_plot_scale()))
    return nothing
end

function _show_colorbar!(c::LiveCanvas, on::Bool)
    # `Colorbar` has no `visible` attribute (checked, not assumed). Width 0 plus a collapsed
    # column is NOT enough on its own: the ticks, their labels and the axis line still draw,
    # leaving a ghost 0..1 scale floating to the right of every baseline-coloured plot. Each
    # decoration has to be switched off by name.
    sc = live_plot_scale()
    c.colorbar.width[]            = on ? 14 * sc : 0
    c.colorbar.ticksvisible[]     = on
    c.colorbar.ticklabelsvisible[] = on
    c.colorbar.labelvisible[]     = on
    c.colorbar.ticklabelsize[]    = OIPLOT_TICKLABELSIZE * sc
    c.colorbar.labelsize[]        = OIPLOT_LABELSIZE * sc
    # Hide the whole Block. Switching decorations off one at a time does not do it -- ticks,
    # tick labels, label and all four spines were named individually and a bare vertical line
    # still ran down the side of every baseline-coloured plot, because at width 0 what remains
    # is the gradient's own outline. A Block draws into its own Scene, and that Scene has the
    # `visible` a Colorbar lacks.
    getfield(c.colorbar, :blockscene).visible[] = on
    on || (c.cbarlabel[] = "")
    Makie.colsize!(c.figure.layout, 2, on ? Makie.Auto() : Makie.Fixed(0))
    # Sit beside the axis, not adrift from it. The default column gap is generous for a
    # figure with several panels and looks detached with one.
    Makie.colgap!(c.figure.layout, 1, on ? 8 : 0)
    # Hug the colorbar rather than sit centred in the cell. An isotropic plot -- uv coverage --
    # is square, so under DataAspect it fills only part of a wide cell and a colorbar pinned to
    # the cell's right edge ends up hundreds of pixels adrift of the data it describes.
    c.axis.halign[] = on ? :right : :center
    return nothing
end

"""
    ramp_colors(values; colormap = WAV_COLORMAP) -> Vector{RGBAf}

Map numbers onto colours here rather than handing Makie a numeric vector.

The live scatter's `color` Observable has to keep ONE element type for the lifetime of the
plot -- swapping a `Vector{RGBAf}` for a `Vector{Float64}` when the user picks "wav" would
change the attribute's type and force the plot to be rebuilt, which is the very thing this
file exists to avoid. So every colour mode resolves to RGBAf, and the colorbar is a standalone
block carrying its own limits instead of being attached to the plot.
"""
function ramp_colors(values::AbstractVector{<:Real}; colormap = WAV_COLORMAP)
    cm = Makie.to_colormap(colormap)
    isempty(values) && return Makie.RGBAf[]
    lo, hi = extrema(values)
    n = length(cm)
    hi == lo && return fill(Makie.RGBAf(cm[1]), length(values))
    idx = clamp.(round.(Int, (values .- lo) ./ (hi - lo) .* (n - 1)) .+ 1, 1, n)
    return [Makie.RGBAf(cm[i]) for i in idx]
end

"""
    canvas_data(d, kind; color, conjugate, logscale) -> NamedTuple

Everything needed to draw, computed and nothing drawn.

Shared with `draw!` through the same helpers (`uv_point_labels`, `group_names`,
`baseline_color_map`, `OBS_SPECS`), and asserted equal to it in the tests.
"""
function canvas_data(d, kind::Symbol; color::Union{Nothing,Symbol} = nothing,
                     conjugate::Bool = true, logscale::Bool = false)
    if kind === :uv
        color = color === nothing ? :baseline : color
        names, info = uv_point_labels(d)
        u = Float64.(d.uv[1, :]) ./ 1e6
        v = Float64.(d.uv[2, :]) ./ 1e6
        x, y, inf = conjugate ? (vcat(u, -u), vcat(v, -v), vcat(info, info)) : (u, v, info)

        legend = Pair{String,Makie.RGBAf}[]
        cvals  = Float64[]
        clabel = ""
        cmap   = WAV_COLORMAP
        cols = if color === :baseline
            cmap = baseline_color_map(names)
            legend = [n => Makie.RGBAf(Makie.to_color(cmap[n])) for n in sort(unique(names))]
            cs = [Makie.RGBAf(Makie.to_color(cmap[n])) for n in names]
            conjugate ? vcat(cs, cs) : cs
        elseif color === :wav
            w = Float64.(d.uv_lam) .* 1e6
            cvals = conjugate ? vcat(w, w) : w
            clabel = WAVELENGTH_LABEL
            ramp_colors(cvals; colormap = cmap)
        elseif color === :mjd
            m = Float64.(d.uv_mjd)
            cvals = conjugate ? vcat(m, m) : m
            clabel = "MJD"
            cmap   = MJD_COLORMAP
            ramp_colors(cvals; colormap = cmap)
        else
            fill(Makie.RGBAf(Makie.to_color("#1f77b4")), length(x))
        end

        return (; x, y, err = Float64[], colors = cols, info = inf, legend, cvals, clabel, cmap,
                xlabel = "U (Mλ)", ylabel = "V (Mλ)",
                title = "uv coverage — " * basename(d.filename),
                isotropic = true, logscale = false, ylims = nothing, markersize = 6.0f0)
    end

    spec = get(OBS_SPECS, kind, nothing)
    spec === nothing && throw(ArgumentError(
        "unknown plot $(repr(kind)); use :uv or one of $(sort(collect(keys(OBS_SPECS))))"))
    color = color === nothing ? spec.default_color : color

    if kind === :visphi && occursin("differential", lowercase(String(d.phityp)))
        throw(ArgumentError(
            "this dataset has phityp=\"differential\": its phase is measured against the " *
            "other channels, so it belongs against WAVELENGTH — use the `diffphi` view."))
    end

    names = group_names(d, spec)
    x = Float64.(getfield(d, spec.x)) .* spec.xscale
    y = Float64.(getfield(d, spec.y))
    e = Float64.(getfield(d, spec.yerr))
    isempty(x) && throw(ArgumentError("no $(kind) data in $(basename(d.filename))"))
    info = [obs_info(d, spec, i, names) for i in eachindex(x)]

    legend = Pair{String,Makie.RGBAf}[]
    cvals  = Float64[]
    clabel = ""
    cmap   = WAV_COLORMAP
    cols = if color === :wav || color === :mjd
        v = Float64.(getfield(d, color === :wav ? spec.lam : spec.mjd))
        color === :wav && (v = v .* 1e6)
        cvals  = v
        clabel = color === :wav ? WAVELENGTH_LABEL : "MJD"
        cmap   = color === :wav ? WAV_COLORMAP : MJD_COLORMAP
        ramp_colors(v; colormap = cmap)
    elseif color === :none
        fill(Makie.RGBAf(Makie.to_color("#1f77b4")), length(x))
    else
        cmap = baseline_color_map(names)
        legend = [n => Makie.RGBAf(Makie.to_color(cmap[n])) for n in sort(unique(names))]
        [Makie.RGBAf(Makie.to_color(cmap[n])) for n in names]
    end

    ylims = nothing
    lo_err = e
    yv = y
    if logscale
        # Non-positive points are DROPPED, not clamped to the axis floor. This is what
        # matplotlib does under `set_yscale("log")`, so it is what `plot_v2(...; logplot=true)`
        # already does and what the port tests compare against. Clamping instead piles every
        # noise-dominated V² onto the bottom of the axis, where it reads as a measurement at
        # that value rather than as one the log axis cannot show.
        keep = [isfinite(v) && v > 0 for v in y]
        any(keep) || throw(ArgumentError(
            "cannot use a log scale: no positive $(kind) values in $(basename(d.filename))"))
        x    = x[keep]
        yv   = y[keep]
        e    = e[keep]
        info = info[keep]
        cols = cols[keep]
        isempty(cvals) || (cvals = cvals[keep])
        l, h = extrema(yv)
        # The lower whisker still reaches below zero on a noisy point. Makie transforms the
        # bar's bounding box, so log10 would be called on the negative; truncate at the floor.
        floorv = l / 10
        lo_err = yv .- max.(yv .- e, floorv)
        ylims  = (l / 2, h * 2)
    end

    return (; x, y = yv, err = e, errlow = lo_err, colors = cols, info, legend, cvals, clabel,
            cmap, xlabel = spec.xlabel, ylabel = spec.ylabel,
            title = string(kind) * " — " * basename(d.filename),
            isotropic = false, logscale, ylims, markersize = 7.0f0)
end

"""
    update_canvas!(canvas, d, kind; color, conjugate, logscale, markersize) -> Vector{String}

Push a dataset onto the live surface. Assignment only — no plot is created, so this is safe to
call from a QML callback with no GL context bound. Returns the per-point info strings.
"""
function update_canvas!(c::LiveCanvas, d, kind::Symbol;
                        color::Union{Nothing,Symbol} = nothing, conjugate::Bool = true,
                        logscale::Bool = false, markersize::Union{Nothing,Real} = nothing)
    pd = canvas_data(d, kind; color = color, conjugate = conjugate, logscale = logscale)

    ax = c.axis
    ax.xlabel[] = pd.xlabel
    ax.ylabel[] = pd.ylabel
    ax.title[]  = pd.title
    ax.aspect[] = pd.isotropic ? Makie.DataAspect() : nothing
    sc = live_plot_scale()
    style_axis!(ax; scale = sc)
    c.legfontsize[] = Float32(UVPLOT_LEGEND_SIZE * sc)

    # `identity` first, always. Changing the scale transforms whatever points are loaded AT
    # THAT MOMENT, so setting log10 before the new data arrives calls log10 on the previous
    # plot's values -- a negative V², a closure phase -- and throws DomainError. Going the
    # other way has the same problem: pushing a closure phase while the axis is still log10.
    # `identity` accepts anything, so it is the only safe state in which to swap the data.
    ax.yscale[] = identity

    # The panel's size is absolute, in points, and does NOT take the plot scale on top: it is
    # set by eye against what is on screen, so scaling it again would move it as soon as the
    # window changed monitors.
    c.markersize[] = MARKER_SIZE_USER[] > 0 ? Float32(MARKER_SIZE_USER[]) :
                     Float32(something(markersize, pd.markersize) * sc)
    c.points[]     = Makie.Point2f.(pd.x, pd.y)
    c.colors[]     = pd.colors

    if isempty(pd.err)
        c.errpoints[] = Makie.Point2f[]
        c.errlow[]    = Float32[]
        c.errhigh[]   = Float32[]
    else
        c.errpoints[] = Makie.Point2f.(pd.x, pd.y)
        c.errlow[]    = Float32.(hasproperty(pd, :errlow) ? pd.errlow : pd.err)
        c.errhigh[]   = Float32.(pd.err)
    end

    set_legend!(c, pd.legend)

    if isempty(pd.cvals)
        _show_colorbar!(c, false)
    else
        lo, hi = extrema(pd.cvals)
        if hi == lo
            # A single-valued colour axis: one wavelength, or one epoch. Makie needs lo < hi,
            # and `lo + 1` -- the obvious filler -- draws a colorbar spanning a whole invented
            # micron, which reads as a range the data does not have. A narrow symmetric band
            # around the value says "everything is at this value" instead.
            pad = max(abs(lo) * 1.0f-3, 1.0f-6)
            c.cbarlimits[] = (Float32(lo - pad), Float32(lo + pad))
        else
            c.cbarlimits[] = (Float32(lo), Float32(hi))
        end
        c.cbarlabel[]  = pd.clabel
        c.cbarmap[]    = hasproperty(pd, :cmap) ? pd.cmap : WAV_COLORMAP
        _show_colorbar!(c, true)
    end

    # Limits from the data while the axis is still linear, then log10 last of all -- by which
    # point every loaded value is positive, because canvas_data clamped them to the floor.
    Makie.autolimits!(ax)
    let fl = ax.finallimits[]
        c.homespan[] = (Float64(fl.widths[1]), Float64(fl.widths[2]))
    end
    if pd.logscale
        pd.ylims !== nothing && Makie.ylims!(ax, pd.ylims...)
        ax.yscale[] = log10
        # style_axis! has just forced LinearTicks, which on a log axis places ten evenly
        # SPACED labels at evenly VALUED positions -- so they bunch into the top decade and
        # leave the rest of the axis unlabelled. LogTicks puts them one per decade.
        ax.yticks[] = Makie.LogTicks(Makie.LinearTicks(OIPLOT_XTICKS))
    end
    return pd.info
end
