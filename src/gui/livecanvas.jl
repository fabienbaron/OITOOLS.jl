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
    legms    = Makie.lift(f -> f * 1.1f0, legfs)
    Makie.scatter!(legax, legmarks; color = legcols, marker = OIPLOT_MARKER, markersize = legms)
    Makie.text!(legax, legtpos; text = legtxt, align = (:left, :center), fontsize = legfs)

    cblim = Makie.Observable((0.0f0, 1.0f0))
    cblab = Makie.Observable("")
    cbar  = Makie.Colorbar(fig[1, 2]; colormap = LIVE_COLORMAP, limits = cblim, label = cblab)

    # The reconstructed image, created here and hidden, like everything else on this canvas.
    # Building a plot after the window exists allocates GL buffers with no context bound, which
    # is the failure this whole pre-create design is here to avoid.
    imdata = Makie.Observable(zeros(Float32, 2, 2))
    imx    = Makie.Observable(Float32[0, 1])
    imy    = Makie.Observable(Float32[0, 1])
    implot = Makie.heatmap!(ax, imx, imy, imdata; colormap = LIVE_COLORMAP, visible = false)

    canvas = LiveCanvas(fig, ax, points, colors, msize, sc,
                        errpts, errlo, errhi, errplot,
                        legax, legmarks, legcols, legtpos, legtxt, legfs,
                        cblim, cblab, cbar,
                        imdata, imx, imy, implot)
    set_legend!(canvas, Pair{String,Makie.RGBAf}[])
    _show_colorbar!(canvas, false)
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
    _show_colorbar!(c, true)

    c.axis.xlabel = "α (mas)"
    c.axis.ylabel = "δ (mas)"
    Makie.limits!(c.axis, half, -half, -half, half)
    c.axis.aspect = 1
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

"Colormap used for the continuous colour modes; also what the colorbar shows."
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
    v = tryparse(Float64, get(ENV, "OITOOLSGUI_PLOT_SCALE", ""))
    (v === nothing || !isfinite(v) || v <= 0) && return 1.7
    return clamp(v, 0.5, 5.0)
end

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
    on || (c.cbarlabel[] = "")
    Makie.colsize!(c.figure.layout, 2, on ? Makie.Auto() : Makie.Fixed(0))
    # Sit beside the axis, not adrift from it. The default column gap is generous for a
    # figure with several panels and looks detached with one.
    Makie.colgap!(c.figure.layout, 1, on ? 8 : 0)
    return nothing
end

"""
    ramp_colors(values; colormap = LIVE_COLORMAP) -> Vector{RGBAf}

Map numbers onto colours here rather than handing Makie a numeric vector.

The live scatter's `color` Observable has to keep ONE element type for the lifetime of the
plot -- swapping a `Vector{RGBAf}` for a `Vector{Float64}` when the user picks "wav" would
change the attribute's type and force the plot to be rebuilt, which is the very thing this
file exists to avoid. So every colour mode resolves to RGBAf, and the colorbar is a standalone
block carrying its own limits instead of being attached to the plot.
"""
function ramp_colors(values::AbstractVector{<:Real}; colormap = LIVE_COLORMAP)
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
        cols = if color === :baseline
            cmap = baseline_color_map(names)
            legend = [n => Makie.RGBAf(Makie.to_color(cmap[n])) for n in sort(unique(names))]
            cs = [Makie.RGBAf(Makie.to_color(cmap[n])) for n in names]
            conjugate ? vcat(cs, cs) : cs
        elseif color === :wav
            w = Float64.(d.uv_lam) .* 1e6
            cvals = conjugate ? vcat(w, w) : w
            clabel = "λ (μm)"
            ramp_colors(cvals)
        elseif color === :mjd
            m = Float64.(d.uv_mjd)
            cvals = conjugate ? vcat(m, m) : m
            clabel = "MJD"
            ramp_colors(cvals)
        else
            fill(Makie.RGBAf(Makie.to_color("#1f77b4")), length(x))
        end

        return (; x, y, err = Float64[], colors = cols, info = inf, legend, cvals, clabel,
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
    cols = if color === :wav || color === :mjd
        v = Float64.(getfield(d, color === :wav ? spec.lam : spec.mjd))
        color === :wav && (v = v .* 1e6)
        cvals  = v
        clabel = color === :wav ? "λ (μm)" : "MJD"
        ramp_colors(v)
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
            xlabel = spec.xlabel, ylabel = spec.ylabel,
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

    c.markersize[] = Float32(something(markersize, pd.markersize) * sc)
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
        _show_colorbar!(c, true)
    end

    # Limits from the data while the axis is still linear, then log10 last of all -- by which
    # point every loaded value is positive, because canvas_data clamped them to the floor.
    Makie.autolimits!(ax)
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
