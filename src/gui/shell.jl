# The application shell.
#
# Four tabs — Observe, Explore, Model, Image — over one Session, plus three regions that live
# OUTSIDE the tabs and must never be per-tab:
#
#   * the context bar, because the current dataset is what the four perspectives share; if
#     each tab picked its own, this would be four programs in one window;
#   * the task tray, because a reconstruction started in Image has to survive switching to
#     Explore to look at residuals;
#   * the command log, because it records across modes and the exported script is one story.
#
# Julia owns all state. QML holds no model of its own; it reads through `state` and calls back
# for everything else, so the GUI cannot drift from what an exported script would do.

"""
    ShellState

Live application state: the session, the single Makie figure the window shows, and the
per-point descriptions backing click-to-identify.

One figure for the lifetime of the window, re-plotted in place (`plot_into!`). Swapping
figures would tear down and rebuild the GL screen on every view change.
"""
mutable struct ShellState
    session :: Session
    figure  :: Any
    axis    :: Any
    plotobj :: Any
    info    :: Vector{String}
    extras  :: Vector{Any}  # legend / colorbar blocks of the current draw, removed on redraw
    current :: Int          # index into session.datasets, 0 = none
    kind    :: Symbol
    color   :: Symbol
    status  :: String
    console :: Vector{String}   # what the bottom pane shows: commands and outcomes, in order
    canvas  :: Any              # LiveCanvas once the window exists; nothing when headless
    picktext:: String           # last identified point; QML polls this
end

# Cap the console so a long session cannot grow the buffer without bound. The oldest lines go
# first; the command log itself (session.log) is never trimmed, so script export is unaffected.
const CONSOLE_MAX = 500

"""
    console!(sh, text; kind = :info)

Append one line to the console pane.

`kind` picks the marker, and the markers are the point: `:cmd` lines are the OITOOLS call that
just ran and are replayable, everything else is commentary. That keeps the pane readable as a
transcript rather than a log dump, and it matches what `export_script` will write -- the `>`
lines are the script, in order.

    > load_dataset!(session, "iota_peg4t.oifits")
      1520 V2 points, 1 wav-bin x 1 time-bin
    ! cannot plot t3phi/wav: BoundsError ...
"""
function console!(sh::ShellState, text::AbstractString; kind::Symbol = :info)
    mark = kind === :cmd ? "> " : kind === :err ? "! " : "  "
    stamp = Dates.format(Dates.now(), "HH:MM:SS")
    push!(sh.console, string(stamp, " ", mark, text))
    length(sh.console) > CONSOLE_MAX && deleteat!(sh.console, 1:(length(sh.console) - CONSOLE_MAX))
    return text
end

"""Everything the console pane shows, oldest first."""
shell_console() = join(_shell().console, "\n")

const SHELL = Ref{Union{Nothing,ShellState}}(nothing)

_shell() = SHELL[] === nothing ? error("no GUI running") : SHELL[]

"Current dataset entry, or `nothing`."
function current_dataset(sh::ShellState = _shell())
    (sh.current < 1 || sh.current > length(sh.session.datasets)) && return nothing
    return sh.session.datasets[sh.current]
end

"""
    refresh_plot!(shell) -> String

Re-draw the current view into the existing axis. Returns a status line for the context bar.
"""
function refresh_plot!(sh::ShellState = _shell())
    e = current_dataset(sh)
    e === nothing && return "no dataset loaded"
    d = e.data[1, 1]
    if sh.canvas === nothing
        # Headless/offline: no window, so creating plots is free. This is the path the unit
        # tests take.
        sh.plotobj, sh.info = plot_into!(sh.figure, sh.axis, d, sh.kind;
                                         color = sh.color, extras = sh.extras)
        Makie.autolimits!(sh.axis)
    else
        # Live: assignment only. See src/livecanvas.jl for why nothing may be created here.
        sh.info    = update_canvas!(sh.canvas, d, sh.kind; color = sh.color)
        sh.plotobj = sh.canvas.scatterplot
    end
    return string(e.name, " — ", basename(e.path), " — ", sh.kind,
                  " (", length(sh.info), " points)")
end

# ── callbacks reachable from QML ──────────────────────────────────────────────

"""
    _cause(err) -> String

One-line description of the *innermost* cause of a possibly-wrapped exception.

Wrappers hide the useful part. ComputePipeline (which drives Makie's attributes) raises a
`ResolveException` whose first line only names the node that failed --

    Failed to resolve gl_renderobject:

-- and the actual error is nested in its `.error` field. Reporting `first(split(...))` of the
wrapper therefore produced a status message that named a GPU internal and said nothing about
what went wrong, for every possible underlying failure. `LoadError` and `CapturedException`
wrap the same way.
"""
function _cause(err)
    e = err
    outer = nothing
    for _ in 1:8                                   # bounded: wrappers can nest, cycles cannot
        inner = if e isa LoadError
            e.error
        elseif hasproperty(e, :error) && getproperty(e, :error) isa Exception
            e.error
        elseif hasproperty(e, :ex) && getproperty(e, :ex) isa Exception
            e.ex
        elseif e isa TaskFailedException
            e.task.exception
        else
            nothing
        end
        inner === nothing && break
        outer === nothing && (outer = e)
        e = inner
    end
    line = first(split(sprint(showerror, e), "\n"))
    # Keep the wrapper's node name when there was one: "gl_renderobject" is not the cause,
    # but it does say where the failure surfaced, which is worth one word.
    if outer isa Exception && hasproperty(outer, :start)
        return string(getproperty(outer, :start).name, ": ", line)
    end
    return String(line)
end

# ── automation hook ──────────────────────────────────────────────────────────
#
# A function run once from a QML timer just after the window appears. It exists so the GUI can
# be driven without a human: actions fired from here run on Qt's GUI thread, after the GL
# screen is live, which is precisely the context that plain unit tests cannot reproduce and
# where the interesting failures happen (a plot drawn before `loadqml` succeeds even when
# every later one fails).
#
# Production runs leave it `nothing`, so the timer calls a no-op.
const READY_HOOK = Ref{Any}(nothing)

"""Run the automation hook, if one was installed. Never throws: a failing hook must not take
the window down with it, or the transcript explaining the failure dies too."""
function shell_ready()
    f = READY_HOOK[]
    f === nothing && return ""
    sh = _shell()
    try
        f()
        console!(sh, "automation hook finished")
    catch err
        console!(sh, "automation hook failed: " * _cause(err); kind = :err)
    end
    return sh.status
end

"""Load an OIFITS chosen in the file dialog. Returns a status string for the UI."""
function shell_open(path::AbstractString)
    sh = _shell()
    p = String(path)
    startswith(p, "file://") && (p = p[8:end])
    console!(sh, "load_dataset!(session, \"$p\")"; kind = :cmd)
    try
        load_dataset!(sh.session, p; warn = false, verbose = false)
        sh.current = length(sh.session.datasets)
        sh.status = refresh_plot!(sh)
        console!(sh, sh.status)
    catch err
        sh.status = "could not load $(basename(p)): " * _cause(err)
        console!(sh, sh.status; kind = :err)
    end
    return sh.status
end

"""Switch the view (`"uv"`, `"v2"`, `"t3phi"`, `"t3amp"`) and/or the colour encoding."""
function shell_set_view(kind::AbstractString, color::AbstractString)
    sh = _shell()
    old_kind, old_color = sh.kind, sh.color
    sh.kind  = Symbol(kind)
    sh.color = Symbol(color)
    # An exception thrown out of a QML callback terminates the application. Report it in the
    # status bar and roll back instead, so a dataset without T3 or a colour mode that does
    # not apply is an inconvenience rather than a crash.
    console!(sh, "plot $(kind), coloured by $(color)"; kind = :cmd)
    try
        sh.status = refresh_plot!(sh)
        console!(sh, sh.status)
    catch err
        sh.kind, sh.color = old_kind, old_color
        sh.status = "cannot plot $(kind)/$(color): " * _cause(err)
        console!(sh, sh.status; kind = :err)
        try; refresh_plot!(sh); catch; end
    end
    return sh.status
end

"""Select dataset `i` (1-based) from the context bar."""
function shell_select_dataset(i::Integer)
    sh = _shell()
    (1 <= i <= length(sh.session.datasets)) || return sh.status
    sh.current = Int(i)
    try
        sh.status = refresh_plot!(sh)
    catch err
        sh.status = "cannot plot dataset $i: " * _cause(err)
    end
    return sh.status
end

"""Dataset names for the context-bar dropdown, newline-separated (QML reads a plain string)."""
shell_dataset_names() = join([e.name * " — " * basename(e.path)
                              for e in _shell().session.datasets], "\n")

"""The command log as a runnable script, for the log panel and 'Export script'."""
shell_script() = export_script(_shell().session)

"""
Reset the plot back to the whole dataset, undoing any zoom or pan.

`autolimits!`, not `reset_limits!`: the latter honours an explicitly-set limits attribute,
which the log-scale path sets, so it would restore the zoom rather than clear it. This only
recomputes limits -- no plot is created -- so it is safe from a QML callback.

Deliberately NOT written to the console: zoom is view state, and §1.2 of the design keeps
view state out of the command log so that the exported script stays a description of what was
computed rather than of where the user happened to be looking.
"""
function shell_reset_zoom()
    sh = _shell()
    try
        Makie.autolimits!(sh.axis)
    catch err
        sh.status = "could not reset the view: " * _cause(err)
        console!(sh, sh.status; kind = :err)
    end
    return sh.status
end

"""Write the exported script to `path`."""
function shell_export(path::AbstractString)
    p = String(path); startswith(p, "file://") && (p = p[8:end])
    export_script(_shell().session, p)
    return "exported to " * basename(p)
end

"""
    shell_pick(; tol = 0.02) -> String

Identify the plotted point nearest the cursor, or `""` if the cursor is not near one.

`Makie.pick` is unavailable here: GPU picking reads a framebuffer GLMakie fills during its own
render loop, and QMLMakie builds its screen with `start_renderloop = false` because Qt drives
rendering, so that buffer is never populated and `pick` returns `nothing` everywhere.

The search therefore runs over the points already held. Distances are in axis fractions, not
data units: U in Mλ and V² dimensionless are not comparable, and a plain Euclidean distance
would be dominated by whichever axis carries the larger numbers. `tol` is that fraction of the
axis diagonal, so 2% is roughly a finger's width on screen. `oiplot`'s `onclickidentify`
normalises the same way.

`sh.info` is index-aligned with the drawn points by construction, including uv coverage's
conjugate half, where the points and the info vector are duplicated together.
"""
function shell_pick(; tol::Real = 0.02)
    sh = _shell()
    isempty(sh.info) && return ""
    sh.canvas === nothing && return ""
    pts = sh.canvas.points[]
    isempty(pts) && return ""
    try
        scene = sh.axis.scene
        # `Makie.mouseposition(scene)` gives data coordinates for THIS scene. `events.mouseposition`
        # is window-relative, so feeding it to `to_world` on the axis scene ignores the axis's
        # viewport offset and lands somewhere else entirely.
        world = Makie.mouseposition(scene)
        lims  = sh.axis.finallimits[]
        wx, wy = Float64(lims.widths[1]), Float64(lims.widths[2])
        (wx > 0 && wy > 0) || return ""
        cx, cy = Float64(world[1]), Float64(world[2])

        best, bestd = 0, Inf
        for (i, q) in enumerate(pts)
            dx = (Float64(q[1]) - cx) / wx
            dy = (Float64(q[2]) - cy) / wy
            d  = dx * dx + dy * dy
            if d < bestd
                bestd, best = d, i
            end
        end
        (best < 1 || best > length(sh.info)) && return ""
        return sqrt(bestd) <= tol ? sh.info[best] : ""
    catch
        return ""
    end
end


"""
    install_interactions!(sh)

Wire left-click identification and right-click view reset, reading Makie's own event stream.

**Why not QML.** `MakieArea` fills itself with a `MouseArea` that has
`acceptedButtons: Qt.AllButtons` and becomes active as soon as the area takes focus. A
`MouseArea` takes an *exclusive* pointer grab on press, which cancels any competing
`TapHandler` -- so a QML handler of ours fires only until the plot is first clicked, and never
again. That is exactly the "clicking does nothing" symptom, and no arrangement of QML handlers
fixes it, because the grab is not ours to lose.

Makie already receives those events (that is how zoom and pan work), so the reliable place to
listen is the scene. Nothing here creates a plot, so it is safe under the live GL constraint.

Right click is a click, not a drag: Makie binds right-drag to panning, so resetting on press
would fight it. The reset happens on RELEASE, and only if the pointer has barely moved since
the press -- which is the same distinction a TapHandler's drag threshold makes.
"""
function install_interactions!(sh::ShellState)
    scene  = sh.axis.scene
    events = Makie.events(scene)
    downat = Ref((0.0, 0.0))

    Makie.on(events.mousebutton, priority = 100) do ev
        pos = Tuple(Float64.(events.mouseposition[]))
        if ev.button == Makie.Mouse.left && ev.action == Makie.Mouse.press
            sh.picktext = shell_pick()
            if get(ENV, "OITOOLSGUI_DEBUG_PICK", "") == "1"
                console!(sh, "pick at $(round.(pos; digits=1)): ninfo=$(length(sh.info)) " *
                             "-> $(isempty(sh.picktext) ? "<nothing>" : sh.picktext)")
            end
        elseif ev.button == Makie.Mouse.right
            if ev.action == Makie.Mouse.press
                downat[] = pos
            elseif ev.action == Makie.Mouse.release &&
                   hypot(pos[1] - downat[][1], pos[2] - downat[][2]) < 5
                shell_reset_zoom()
            end
        end
        return Makie.Consume(false)      # never swallow: zoom and pan still need these
    end
    return nothing
end

"""Last identified point. QML polls this rather than asking for a pick itself."""
shell_pick_text() = _shell().picktext

"""
    check_qt_conflict()

Fail early, and legibly, if a second Qt is already in the process.

Two Qt builds in one process do not coexist. QML.jl's is **Qt 6.10.2** from its JLL, and a
foreign one makes its libraries fail with

    libQt6DBus.so: undefined symbol: _ZN14QObjectPrivateC2E16QtPrivate_6_10_2

which names the symbol rather than the cause, so the check below turns it into a sentence.

Historically the foreign Qt came from OITOOLS itself: `using OITOOLS` imported matplotlib
unconditionally, matplotlib probed for an interactive backend, found **PySide6** in OITOOLS'
CondaPkg environment (`ultranest → getdist → pyside6 → qt6-main`) and mapped its **Qt 6.11.1**
into the process. The workaround was to set `MPLBACKEND=Agg` before OITOOLS loaded.

OITOOLS keeps all matplotlib in an extension (`OITOOLSPythonPlotExt`), so importing OITOOLS
starts no probe and maps no Qt, and `MPLBACKEND` needs no forcing.

The check remains because the failure mode is still reachable: it fires if the session loads
PythonPlot itself, or any other package carrying a system or conda Qt. The remedy is to keep
that package out of the GUI process, or to force it to a non-interactive backend before
loading it.
"""
function check_qt_conflict()
    Sys.islinux() || return nothing
    isfile("/proc/self/maps") || return nothing
    foreign = String[]
    for line in eachline("/proc/self/maps")
        occursin("libQt6Core", line) || continue
        path = String(split(line)[end])
        occursin("/artifacts/", path) || push!(foreign, path)
    end
    isempty(foreign) && return nothing
    error("""
          A second Qt is already loaded, so QML cannot start:
              $(first(unique(foreign)))

          QML.jl needs its own Qt (6.10.2, from its JLL) to be the only one in the process.
          The usual source is matplotlib picking an interactive backend and mapping the Qt
          that ships with PySide6 in a conda environment.

          OITOOLS does not do this — its matplotlib lives in an extension — so something
          else in this session loaded one. Start the GUI in a process that does not import
          PythonPlot:
              julia --project=. -e 'using OITOOLS, GLMakie, QMLMakie, QML; gui()'
          or use bin/oitoolsgui.jl. If you do need matplotlib in the same process, force it
          to a non-interactive backend before loading it:
              ENV["MPLBACKEND"] = "Agg"
          """)
end

"""
    gui(session = Session(); qmlfile, autoquit_ms = 0) -> Session

Open the application window. Returns when it closes.

Requires the Qt stack, which is a weak dependency:

    using QML, QMLMakie, GLMakie      # enables this method

Calls [`check_qt_conflict`](@ref) and [`configure_graphics!`](@ref) first: the one refuses to
start beside a foreign Qt, the other gives Mesa the driver hint it needs on WSL. Both are
no-ops on a healthy Linux desktop, and `bin/oitoolsgui.jl` already runs the second one early,
which is where it belongs.
"""
function gui end
