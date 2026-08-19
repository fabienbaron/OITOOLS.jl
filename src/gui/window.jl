# The window: Qt callback registration and `gui()`.
#
# Split from the rest of src/gui/ because everything else describes what to draw, while this
# builds the thing that draws it and hands control to Qt's event loop.

function __init__()
    # Callbacks must be registered before any QML file that calls them is loaded.
    QML.@qmlfunction(shell_open, shell_set_view, shell_select_dataset,
                     shell_dataset_names, shell_script, shell_export, shell_pick,
                     shell_console, shell_ready, shell_reset_zoom, shell_pick_text)
    return nothing
end

"""
Folder the file dialog opens in: next to the last dataset loaded if there is one, otherwise
the bundled `data/` directory, otherwise the working directory. Opening in the filesystem
root or the home directory is a small tax paid on every single open, and it makes automated
clicking guesswork.
"""
function _initial_folder(session::Session)
    if !isempty(session.datasets)
        d = dirname(abspath(session.datasets[end].path))
        isdir(d) && return "file://" * d
    end
    bundled = joinpath(pkgdir(OITOOLS), "test", "gui", "data")
    isdir(bundled) && return "file://" * bundled
    return "file://" * pwd()
end

"Tab names, in the order Main.qml lists them."
const TAB_NAMES = ("explore", "observe", "model", "image")

"""
Which perspective the window opens on, from `\$OITOOLSGUI_TAB` (default Explore).

Its reason for existing is testing: a tab that is never current is constructed but never laid
out, so any layout warning it would emit never appears. Naming a tab is what makes an
automated run render it. An unknown name is reported rather than silently ignored, since a
typo would otherwise look like the setting having no effect.
"""
function _initial_tab()
    want = lowercase(strip(get(ENV, "OITOOLSGUI_TAB", "")))
    isempty(want) && return 0
    i = findfirst(==(want), TAB_NAMES)
    if i === nothing
        @warn "OITOOLSGUI_TAB is not a tab name; opening on Explore" got = want options = TAB_NAMES
        return 0
    end
    return i - 1        # QML indexes from zero
end

function OITOOLS.gui(session::Session = Session();
                        qmlfile::AbstractString = joinpath(pkgdir(OITOOLS), "src", "gui", "qml", "Main.qml"),
                        autoquit_ms::Integer = 0,
                        on_ready = nothing)
    check_qt_conflict()
    # Backstop for sessions started by hand rather than through bin/oitoolsgui.jl. GLMakie is
    # imported by now but has not built a context yet, and Mesa reads the variables at context
    # creation, so this is still early enough. A no-op when the launcher already ran it.
    READY_HOOK[] = on_ready
    gfx = configure_graphics!()

    # QML works the scale out from Screen.logicalPixelDensity; this is only the override.
    # 0.0 means "decide for yourself", which is the default -- see src/scaling.jl.
    uiscale = ui_scale_override()
    @info "UI scale: $(uiscale.reason)"

    fig = Makie.Figure()
    ax  = Makie.Axis(fig[1, 1]; xlabel = "", ylabel = "", title = "OITOOLS")
    style_axis!(ax)
    # Build every plot BEFORE the window exists. Once QML is up, GLMakie allocates GPU buffers
    # on insertion and the context belongs to Qt's render thread -- see src/livecanvas.jl.
    canvas = build_canvas(fig, ax)
    sh  = ShellState(session, fig, ax, nothing, String[], Any[], 0, :uv, :baseline,
                     "no dataset loaded", String[], canvas, "")
    SHELL[] = sh
    install_interactions!(sh)

    if !isempty(session.datasets)
        sh.current = 1
        sh.status  = refresh_plot!(sh)
    end
    # Seed the console before the window exists, so it opens with a transcript rather than an
    # empty box -- including the graphics and scale decisions, which are exactly what you want
    # to see when the window looks wrong.
    console!(sh, "OITOOLSGUI ready")
    console!(sh, "graphics: $(gfx.reason)")
    console!(sh, "ui scale: $(uiscale.reason)")
    for (i, d) in enumerate(session.datasets)
        console!(sh, "load_dataset!(session, \"$(d.path)\")"; kind = :cmd)
        i == sh.current && console!(sh, sh.status)
    end
    isempty(session.datasets) && console!(sh, "no dataset loaded — use Open OIFITS…")

    # `fig`, not `fig.scene`: QMLMakie's own examples pass the Figure, and MakieArea resolves
    # the scene itself (Makie.get_scene). Handing it the root Scene directly skips whatever
    # the figure-level path sets up.
    QML.loadqml(qmlfile;
                plot            = fig,
                autoQuitMs      = Int(autoquit_ms),
                initialTab      = _initial_tab(),
                uiScaleOverride = uiscale.scale,
                initialFolder   = _initial_folder(session),
                fullscreenOnStart = get(ENV, "OITOOLSGUI_FULLSCREEN", "") == "1",
                initialStatus   = sh.status)
    QML.exec()

    # An automated run has no way to read the console pane once the window is gone, so hand it
    # over on the way out. Only active when asked for; production runs never write a file.
    dump = get(ENV, "OITOOLSGUI_CONSOLE_DUMP", "")
    if !isempty(dump)
        try
            open(dump, "w") do io
                for l in sh.console; println(io, l); end
                # Not console output: the last identified point, so an automated run can
                # check that picking works without reading a screenshot.
                println(io, "PICKTEXT: ", sh.picktext)
            end
        catch err
            @warn "could not write console dump" path = dump exception = err
        end
    end
    return session
end
