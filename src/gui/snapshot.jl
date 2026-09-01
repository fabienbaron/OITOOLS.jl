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

`pixsize` is read back off the coordinate vector rather than tracked separately: `show_image!`
builds it as an evenly spaced range of cell centres, so the spacing IS the pixel size and the
two cannot fall out of step.
"""
function _canvas_image(c)
    c === nothing && return nothing
    c.imageplot.visible[] || return nothing
    x = c.imagex[]
    (length(x) < 2 || isempty(c.imagedata[])) && return nothing
    return (image = c.imagedata[], pixsize = abs(Float64(x[2]) - Float64(x[1])))
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
    _snapshot_image_figure(img, pixsize, colormap, size) -> Figure

`show_image!`'s view of an image, rebuilt in a figure of its own: milliarcseconds, East to the
LEFT, and the colormap the buttons last chose.
"""
function _snapshot_image_figure(img, pixsize, colormap, size)
    nx = Base.size(img, 1)
    half = nx * pixsize / 2
    fig = Makie.Figure(; size)
    ax = Makie.Axis(fig[1, 1]; aspect = 1, xlabel = "α (mas)", ylabel = "δ (mas)")
    style_axis!(ax; scale = live_plot_scale())
    image_minorticks!(ax, nx * pixsize, nx * pixsize; scale = live_plot_scale())
    hm = Makie.heatmap!(ax, range(half, -half; length = nx), range(-half, half; length = nx),
                        img; colormap)
    Makie.Colorbar(fig[1, 2], hm; label = "flux / pixel")
    Makie.limits!(ax, half, -half, -half, half)
    return fig
end

"""
    _snapshot_figure(sh, which, size, detailed) -> Figure or nothing

Rebuild one perspective's plot in a fresh figure. `nothing` means the panel has nothing on it
yet, which is a message rather than an error.
"""
function _snapshot_figure(sh::ShellState, which::AbstractString, size, detailed::Bool)
    if which == "explore"
        e = current_dataset(sh)
        e === nothing && return nothing
        d = e.data[1, 1]
        fig = Makie.Figure(; size)
        ax = Makie.Axis(fig[1, 1])
        if sh.panels
            c = build_canvas(fig, ax)
            show_panels!(c, true)
            update_panels!(c, d, sh.kind)
        else
            # The same call `refresh_plot!` makes when there is no window: on this path
            # creating plots is free, because the figure has no GL context to allocate into.
            plot_into!(fig, ax, d, sh.kind; color = sh.color, logscale = sh.logy,
                       extras = Any[])
            Makie.autolimits!(ax)
        end
        return fig

    elseif which == "image" || which == "model"
        c = which == "image" ? sh.imcanvas : sh.modelcanvas
        im = _canvas_image(c)
        im === nothing && return nothing
        return _snapshot_image_figure(im.image, im.pixsize, c.colormap[], size)

    elseif which == "chi2map"
        m = _last_chi2_map(sh)
        m === nothing && return nothing
        fig = Makie.Figure(; size)
        update_chi2_map!(build_chi2_map(fig, Makie.Axis(fig[1, 1])), m)
        return fig

    elseif which == "residuals"
        r = _current_residuals(sh)
        r isa String && return nothing
        fig = Makie.Figure(; size)
        update_residuals!(build_residuals(fig, Makie.Axis(fig[1, 1])), r.data, r.res)
        return fig

    elseif which == "sed"
        sh.sedplot === nothing && return nothing
        isempty(sh.sedplot.total[]) && return nothing
        fig = Makie.Figure(; size)
        panel = build_sed(fig, Makie.Axis(fig[1, 1]))
        # Copied from the live panel rather than recomputed: the file must be the picture on
        # screen, and a second `model_to_sed` could be of a model edited in between.
        panel.total[] = sh.sedplot.total[]
        for i in eachindex(panel.comps)
            panel.comps[i][]  = sh.sedplot.comps[i][]
            panel.labels[i][] = sh.sedplot.labels[i][]
            panel.colors[i][] = sh.sedplot.colors[i][]
        end
        Makie.autolimits!(panel.axis)
        return fig

    elseif which == "gantt" || which == "delay"
        sh.plan === nothing && return nothing
        fig = Makie.Figure(; size)
        ax = Makie.Axis(fig[1, 1])
        if which == "gantt"
            update_gantt!(build_gantt(fig, ax), sh.plan; detailed)
        else
            isempty(sh.plan.baselines) && return nothing
            update_delay_plot!(build_delay_plot(fig, ax), sh.plan)
        end
        return fig
    end
    return "! no such plot: " * String(which)
end

"""
    shell_save_figure(which, path, width, height, detailed) -> String

Write one perspective's plot to a PNG. Returns `""` on success or a message beginning with `!`.

`which` is `explore`, `image`, `model`, `chi2map`, `residuals`, `sed`, `gantt` or `delay` --
one per plot area, not one per tab, because a tab showing two plots would otherwise have to
guess which one was meant.

`width`/`height` come from the panel on screen, so the file is framed as the window frames it
rather than at some size chosen here. `detailed` is the Gantt's own row-detail switch and is
ignored by every other plot.
"""
function shell_save_figure(which, path, width, height, detailed)
    sh = _shell()
    p = strip(String(path))
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

    try
        Makie.save(p, fig)
    catch err
        msg = "! could not write " * p * ": " * _cause(err)
        console!(sh, msg; kind = :err); return msg
    end
    console!(sh, "> save(\"" * p * "\", figure)"; kind = :cmd)
    console!(sh, "  wrote $(p), $(sz[1])×$(sz[2])")
    return ""
end
