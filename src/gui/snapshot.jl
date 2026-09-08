# ── "Save PNG": what every perspective's plot area writes to a file ──────────
#
# Nothing here reads the on-screen framebuffer, and it cannot: QMLMakie swaps the screen's last
# postprocessor for one that renders into Qt's FBO, so `Makie.colorbuffer` on the live screen
# hands back the raw buffer instead of the composited frame -- measured, it comes out as noise
# -- and `GLFW.SwapBuffers` has no method for a `QMLWindow` at all, so `Makie.save` on the live
# figure throws before it gets that far.
#
# A figure that was never displayed has neither problem. `Makie.save` opens its own hidden
# GLMakie screen for it, which works while Qt holds the context -- measured in the same
# session. So a snapshot REBUILDS the panel offscreen, through the very builders the window
# uses, and saves that. What lands in the file is then what the equivalent script would draw,
# which is the property the command log is for.

"""
The image a canvas is currently showing, as `(image, pixsize)`, or `nothing` when it is
showing a scatter plot instead.

`pixsize` is read back off the coordinate vector rather than tracked separately, so the two
cannot fall out of step. `show_image!` spans `nx * pixsize` end to end with `nx` samples, so
the pixel size is that span over `nx` -- not the sample spacing, which is the span over
`nx - 1` and would hand the snapshot a field of view `nx / (nx - 1)` too wide.
"""
function _canvas_image(c)
    c === nothing && return nothing
    c.imageplot.visible[] || return nothing
    x = c.imagex[]
    (length(x) < 2 || isempty(c.imagedata[])) && return nothing
    return (image = c.imagedata[],
            pixsize = abs(Float64(x[end]) - Float64(x[1])) / length(x))
end

"""The most recent fit that produced a χ² map, or `nothing` if none did."""
function _last_chi2_map(sh::ShellState)
    for f in Iterators.reverse(sh.fits)
        m = try f.map catch; nothing end
        m === nothing || return m
    end
    return nothing
end

"""
    _snapshot_image_figure(img, pixsize, colormap, size; label, title) -> Figure

`show_image!`'s view of an image, rebuilt in a figure of its own: milliarcseconds, East to the
LEFT, the colormap the buttons last chose, and the colorbar label the canvas carries -- which
is how a posterior mean or spread says which of the two it is.

`title` comes from the live axis rather than being written here: `show_image!` never touches
it, so the panel still carries whatever the window titled it with.
"""
function _snapshot_image_figure(img, pixsize, colormap, size; label = "flux / pixel",
                                title = "")
    nx = Base.size(img, 1)
    half = nx * pixsize / 2
    fig = Makie.Figure(; size, fonts = PLOT_FONTS)
    ax = Makie.Axis(fig[1, 1]; aspect = 1, xlabel = "α (mas)", ylabel = "δ (mas)", title)
    style_axis!(ax; scale = live_plot_scale())
    image_minorticks!(ax, nx * pixsize, nx * pixsize; scale = live_plot_scale())
    # The canvas anchors the colour scale at zero rather than at the image minimum, so a
    # snapshot left on Makie's automatic range would shade the same reconstruction differently
    # from the panel it was taken from.
    mx = maximum(img)
    hm = Makie.heatmap!(ax, range(-half, half; length = nx), range(-half, half; length = nx),
                        img; colormap, colorrange = (0.0, mx > 0 ? Float64(mx) : 1.0))
    Makie.Colorbar(fig[1, 2], hm; label)
    Makie.limits!(ax, half, -half, -half, half)
    return fig
end

"""
    _live_panel(sh, which) -> NamedTuple or nothing

The on-screen figure and axis a snapshot is meant to reproduce, or `nothing` when that panel
has none -- headless, or a plot the window never built.

Every panel is either a `LiveCanvas` or a NamedTuple carrying `figure` and `axis`, so one
lookup serves all of them.
"""
function _live_panel(sh::ShellState, which::AbstractString)
    p = which == "explore"   ? sh.canvas      :
        which == "image"     ? sh.imcanvas    :
        which == "model"     ? sh.modelcanvas :
        which == "chi2map"   ? sh.chi2map     :
        which == "residuals" ? sh.residplot   :
        which == "sed"       ? sh.sedplot     :
        which == "gantt"     ? sh.gantt       :
        which == "delay"     ? sh.delayplot   : nothing
    p === nothing && return nothing
    (hasproperty(p, :figure) && hasproperty(p, :axis)) || return nothing
    return (; figure = p.figure, axis = p.axis)
end

"""
    _live_size(live, fallback) -> (w, h)

The size to rebuild at: the on-screen scene's own, when there is one.

Taking it from the panel instead of from QML is what makes the file the picture on screen
rather than a differently-shaped one. Makie sizes text, ticks and the legend in POINTS, so a
figure built at twice the scene's width does not scale up -- it gets twice as much room for
the same 12 pt labels, which relocates the ticks, rewraps the legend and thins every marker.
Resolution is bought with `px_per_unit` at save time, where it costs nothing but pixels.
"""
function _live_size(live, fallback)
    live === nothing && return fallback
    sz = try Makie.size(live.figure.scene) catch; nothing end
    sz === nothing && return fallback
    w, h = round(Int, Float64(sz[1])), round(Int, Float64(sz[2]))
    return (w >= 64 && h >= 64) ? (w, h) : fallback
end

"""
    _match_view!(dst, src) -> dst

Give the rebuilt axis the view the live one is showing, zoom and pan included.

`finallimits` is in DATA coordinates on a log axis as much as on a linear one (measured), so
one path serves both. It always reports positive widths, so a reversed axis is a flag rather
than an ordering -- handing `limits!` an ascending pair on one silently clears the reversal,
which is how a saved reconstruction would come out mirrored.
"""
function _match_view!(dst, src)
    (dst === nothing || src === nothing) && return dst
    fl = src.finallimits[]
    x0, y0 = Float64(fl.origin[1]), Float64(fl.origin[2])
    wx, wy = Float64(fl.widths[1]), Float64(fl.widths[2])
    all(isfinite, (x0, y0, wx, wy)) && wx > 0 && wy > 0 || return dst
    xlo, xhi = dst.xreversed[] ? (x0 + wx, x0) : (x0, x0 + wx)
    ylo, yhi = dst.yreversed[] ? (y0 + wy, y0) : (y0, y0 + wy)
    Makie.limits!(dst, xlo, xhi, ylo, yhi)
    return dst
end

"""
    _snapshot_figure(sh, which, size, detailed) -> Figure or nothing

Rebuild one perspective's plot in a fresh figure. `nothing` means the panel has nothing on it
yet, which is a message rather than an error.
"""
function _snapshot_figure(sh::ShellState, which::AbstractString, size, detailed::Bool)
    live = _live_panel(sh, which)
    size = _live_size(live, size)
    liveax = live === nothing ? nothing : live.axis

    if which == "explore"
        e = current_dataset(sh)
        e === nothing && return nothing
        d = e.data[1, 1]
        fig = Makie.Figure(; size, fonts = PLOT_FONTS)
        ax = Makie.Axis(fig[1, 1])
        style_axis!(ax)
        # Rebuilt through the canvas builders the window itself uses, not through `plot_into!`:
        # residuals and the model/image overplot are `update_canvas!` arguments, so a figure
        # assembled any other way saves the plain data plot whatever the ticks say. Creating
        # plots is free here — the figure has no GL context to allocate into.
        c = build_canvas(fig, ax)
        if sh.panels
            show_panels!(c, true)
            update_panels!(c, d, sh.kind)
        else
            rset = residual_set(sh)
            rset isa String && return rset
            oset = overlay_set(sh)
            oset isa String && return oset
            show_panels!(c, false)
            update_canvas!(c, d, sh.kind; color = sh.color, logscale = sh.logy,
                           residual = residual_for(rset, sh.kind),
                           overlay  = residual_for(oset, sh.kind))
            # The zoom, once the data is loaded: `update_canvas!` autolimits, so matching the
            # view has to come after it. The panel grid is left alone -- each of its axes is
            # framed on its own group and none of them is zoomable.
            _match_view!(ax, liveax)
        end
        return fig

    elseif which == "image" || which == "model"
        c = which == "image" ? sh.imcanvas : sh.modelcanvas
        im = _canvas_image(c)
        im === nothing && return nothing
        fig = _snapshot_image_figure(im.image, im.pixsize, c.colormap[], size;
                                     label = c.cbarlabel[],
                                     title = liveax === nothing ? "" : String(liveax.title[]))
        _match_view!(fig.content[findfirst(b -> b isa Makie.Axis, fig.content)], liveax)
        return fig

    elseif which == "chi2map"
        m = _last_chi2_map(sh)
        m === nothing && return nothing
        fig = Makie.Figure(; size, fonts = PLOT_FONTS)
        ax = Makie.Axis(fig[1, 1])
        style_axis!(ax)
        update_chi2_map!(build_chi2_map(fig, ax), m)
        _match_view!(ax, liveax)
        return fig

    elseif which == "residuals"
        r = _current_residuals(sh)
        r isa String && return nothing
        fig = Makie.Figure(; size, fonts = PLOT_FONTS)
        ax = Makie.Axis(fig[1, 1])
        style_axis!(ax)
        update_residuals!(build_residuals(fig, ax), r.data, r.res)
        _match_view!(ax, liveax)
        return fig

    elseif which == "sed"
        sh.sedplot === nothing && return nothing
        isempty(sh.sedplot.total[]) && return nothing
        fig = Makie.Figure(; size, fonts = PLOT_FONTS)
        sax = Makie.Axis(fig[1, 1])
        style_axis!(sax)
        panel = build_sed(fig, sax)
        # Copied from the live panel rather than recomputed: the file must be the picture on
        # screen, and a second `model_to_sed` could be of a model edited in between.
        panel.total[] = sh.sedplot.total[]
        for i in eachindex(panel.comps)
            panel.comps[i][]  = sh.sedplot.comps[i][]
            panel.labels[i][] = sh.sedplot.labels[i][]
            panel.colors[i][] = sh.sedplot.colors[i][]
        end
        Makie.autolimits!(panel.axis)
        _match_view!(panel.axis, liveax)
        return fig

    elseif which == "gantt" || which == "delay"
        sh.plan === nothing && return nothing
        fig = Makie.Figure(; size, fonts = PLOT_FONTS)
        ax = Makie.Axis(fig[1, 1])
        style_axis!(ax)
        if which == "gantt"
            update_gantt!(build_gantt(fig, ax), sh.plan; detailed)
        else
            isempty(sh.plan.baselines) && return nothing
            update_delay_plot!(build_delay_plot(fig, ax), sh.plan)
        end
        _match_view!(ax, liveax)
        return fig
    end
    return "! no such plot: " * String(which)
end

"""
Pixels per figure unit in a saved PNG.

The file is the panel's own layout at twice the resolution, which is what a figure meant for a
paper or a slide needs and what scaling the figure itself cannot give.
"""
const SNAPSHOT_PX_PER_UNIT = 2

"""
    shell_save_figure(which, path, width, height, detailed) -> String

Write one perspective's plot to a PNG. Returns `""` on success or a message beginning with `!`.

`which` is `explore`, `image`, `model`, `chi2map`, `residuals`, `sed`, `gantt` or `delay` --
one per plot area, not one per tab, because a tab showing two plots would otherwise have to
guess which one was meant.

`width`/`height` are the panel's, and are the FALLBACK: with a window up, the size comes from
the live scene itself (`_live_size`), which is the only measurement that makes the file the
same shape as the picture. `detailed` is the Gantt's own row-detail switch and is ignored by
every other plot.
"""
function shell_save_figure(which, path, width, height, detailed)
    sh = _shell()
    # `String`, not the `SubString` `strip` returns, and this is load-bearing: Makie's
    # `save(::String, ::FigureLike)` does not accept a SubString, so a name already ending in
    # `.png` (which skips the concatenation below) missed it and fell through to FileIO's
    # backend loop -- which tried ImageMagick on a Figure and finally wrote the file through
    # `show(io, MIME"image/png", fig)`. That path takes no screen options and renders the
    # figure by itself, so what landed in the file was not the panel that was rebuilt.
    p = String(strip(String(path)))
    isempty(p) && return "! no file name"
    endswith(lowercase(p), ".png") || (p = p * ".png")
    sz = (max(320, round(Int, Float64(width))), max(240, round(Int, Float64(height))))

    fig = try
        _snapshot_figure(sh, String(which), sz, Bool(detailed))
    catch err
        msg = "! could not build the figure: " * _cause(err)
        console!(sh, msg; kind = :err); return msg
    end
    fig isa AbstractString && (console!(sh, fig; kind = :err); return fig)
    fig === nothing && return "! this plot has nothing on it yet"

    # Resolution is bought here rather than by enlarging the figure: `px_per_unit` renders the
    # SAME layout at twice the pixels, where a bigger figure would have re-laid it out with
    # the same point sizes and produced a different picture.
    out = Makie.size(fig.scene)
    ppu = SNAPSHOT_PX_PER_UNIT
    # Named rather than left to `current_backend()`: this figure was never displayed, and
    # GLMakie is the backend that opens a hidden screen for one. Which backend is "current"
    # depends on what the session loaded last, and the file must not.
    try
        Makie.save(p, fig; backend = GLMakie, px_per_unit = ppu)
    catch
        # One retry at the figure's own resolution. `px_per_unit` is a SCREEN option, and a
        # stale environment can route the write through FileIO's MIME path instead, which
        # takes none -- better a correctly laid-out file at 1:1 than no file at all.
        ppu = 1
        try
            Makie.save(p, fig)
        catch err2
            msg = "! could not write " * p * ": " * _cause(err2)
            console!(sh, msg; kind = :err); return msg
        end
    end
    console!(sh, "save(\"" * p * "\", figure; px_per_unit = $(ppu))"; kind = :cmd)
    console!(sh, "  wrote $(p), $(round(Int, out[1] * ppu))×$(round(Int, out[2] * ppu))")
    return ""
end
