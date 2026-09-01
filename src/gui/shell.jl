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
    logy    :: Bool             # log y axis, where the observable allows one
    panels  :: Bool             # one panel per baseline/triplet/station, against wavelength
    status  :: String
    console :: Vector{String}   # what the bottom pane shows: commands and outcomes, in order
    canvas  :: Any              # LiveCanvas once the window exists; nothing when headless
    picktext:: String           # last identified point; QML polls this
    imaging :: Any              # the last ImagingResult, or nothing
    imcanvas:: Any              # the Image perspective's own LiveCanvas
    ftcache :: Any              # (; key, ft, summary) for the Image perspective's plan
    gantt   :: Any              # the Observe perspective's Gantt chart
    delayplot :: Any            # ...and its delay-cart chart
    modelcanvas :: Any          # the Model perspective's rendering canvas
    plan    :: Any              # the NightPlan behind the current Gantt, for the hover readout
    fits    :: Vector{Any}      # completed fits, newest last
    chi2map :: Any              # the Model perspective's χ² map chart
    profileplot :: Any       # the radial-profile preview: I(r) beside its transform
    residplot :: Any            # normalised model residuals, shown in place of the χ² map
    sedplot :: Any              # the model SED, sharing that same rectangle
    enginelog :: String         # everything the last reconstruction printed
    fitlog    :: String          # ...and everything the last model fit printed
    job       :: Any             # the running GuiJob, or nothing
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
                                         color = sh.color, logscale = sh.logy,
                                         extras = sh.extras)
        Makie.autolimits!(sh.axis)
    elseif sh.panels
        # One panel per group, against wavelength. The pool of axes already exists; this only
        # fills it, for the same reason nothing else here creates a plot.
        show_panels!(sh.canvas, true)
        n = update_panels!(sh.canvas, d, sh.kind)
        sh.info = String[]                      # picking is per-point on the single plot only
        total = length(panel_data(d, sh.kind).panels)
        return string(e.name, " — ", basename(e.path), " — ", sh.kind,
                      " (", n, " of ", total, " ", grouping_noun(sh.kind), "s)")
    else
        # Live: assignment only. See src/livecanvas.jl for why nothing may be created here.
        show_panels!(sh.canvas, false)
        sh.info    = update_canvas!(sh.canvas, d, sh.kind; color = sh.color,
                                    logscale = sh.logy)
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
function shell_set_view(kind::AbstractString, color::AbstractString, logy::Bool = false,
                        panels::Bool = false)
    sh = _shell()
    old_kind, old_color, old_log, old_pan = sh.kind, sh.color, sh.logy, sh.panels
    sh.kind  = Symbol(kind)
    sh.color = Symbol(color)
    # Silently ignored where it would be wrong rather than refused: a phase takes both signs,
    # and a log axis on one drops half the points without saying so.
    sh.logy  = logy && Symbol(kind) in LOG_Y_KINDS
    # uv coverage has no groups to split by -- it is geometry, not an observable.
    sh.panels = panels && haskey(OBS_SPECS, Symbol(kind))
    # An exception thrown out of a QML callback terminates the application. Report it in the
    # status bar and roll back instead, so a dataset without T3 or a colour mode that does
    # not apply is an inconvenience rather than a crash.
    console!(sh, "plot $(kind), coloured by $(color)" * (sh.logy ? ", log y" : "") *
                 (sh.panels ? ", per " * grouping_noun(sh.kind) : ""); kind = :cmd)
    try
        sh.status = refresh_plot!(sh)
        console!(sh, sh.status)
    catch err
        sh.kind, sh.color, sh.logy, sh.panels = old_kind, old_color, old_log, old_pan
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
    shell_shift_date(iso, days, months) -> String

Move a date by whole days and/or months, returning ISO `yyyy-mm-dd`.

In Julia because this is calendar arithmetic, not string editing: adding a month to 31 January
has to land on 28 February, and `Dates.Month` clamps to the end of the month while a naive
`+30 days` would not. An empty or unparseable input comes back as today, which is the only
useful answer to "step from nothing".
"""
function shell_shift_date(iso::AbstractString, days::Integer, months::Integer)
    d = try
        Dates.Date(String(iso))
    catch
        Dates.today()
    end
    (days == 0 && months == 0) && return string(Dates.today())
    return string(d + Dates.Month(months) + Dates.Day(days))
end

"""
    shell_fit_model(model_lines, free_lines, constraint_lines, prior_lines,
                    v2, t3amp, t3phi, cvis, flux, diffvis, optimiser, maxeval) -> String

Fit the panel's model to the current dataset. Returns a one-line summary.

The optimiser choice decides more than speed. NLopt takes the constraints as real nonlinear
constraints, so they hold at the optimum; Levenberg-Marquardt and nested sampling have no such
machinery and apply a soft penalty that a steep χ² can overrule. That difference is stated in
the console rather than left to be discovered.
"""
function shell_fit_model(model_lines::AbstractString, free_lines::AbstractString,
                         constraint_lines::AbstractString, prior_lines::AbstractString,
                         v2::Bool, t3amp::Bool, t3phi::Bool,
                         cvis::Bool, flux::Bool, diffvis::Bool,
                         optimiser::AbstractString, maxeval::Integer,
                         gridp1::AbstractString = "", gridp2::AbstractString = "",
                         gridn::Integer = 60, gridn2::Integer = 0)
    sh = _shell()
    e = current_dataset(sh)
    e === nothing && return "! no dataset loaded"

    md = parse_model_lines(model_lines)
    isempty(md) && return "! no model"
    free, lb, ub = parse_free_lines(free_lines)
    isempty(free) && return "! no free parameters"

    cons = try
        parse_constraint_lines(constraint_lines)
    catch err
        msg = "! constraint: " * _cause(err); console!(sh, msg); return msg
    end
    priors  = parse_prior_lines(prior_lines)
    weights = fitting_weights(; v2, t3amp, t3phi, cvis, flux, diffvis)
    all(iszero, weights) && return "! no observables selected"

    data = e.data[1, 1]
    opt  = String(optimiser)

    isempty(cons) || console!(sh, "  constraints: " *
        (startswith(opt, "L") && opt != "lsqfit" ?
         "enforced by NLopt (they hold at the optimum)" :
         "SOFT penalty with this fitter — a steep chi2 can overrule them"))

    # The grid is not a call to `fit_model`, so it does not claim to be one: the console line is
    # the call that was actually made, which is what makes the exported script reproduce it.
    gp1, gp2 = String(gridp1), String(gridp2)
    if opt == "grid"
        length(free) >= 2 || return "! grid search needs two free parameters"
        isempty(gp1) && (gp1 = free[1])
        isempty(gp2) && (gp2 = free[findfirst(!=(gp1), free)])
        # A step count per axis: the two parameters rarely deserve the same resolution -- a
        # diameter may need 200 points where its darkening coefficient needs 40, and forcing
        # them equal spends the whole budget on the axis that did not need it.
        n2 = gridn2 > 1 ? Int(gridn2) : Int(gridn)
        console!(sh, "> chi2_map(model, $(free), data, $(repr(gp1)), $(repr(gp2)); " *
                     "n1 = $(Int(gridn)), n2 = $(n2), weights = $(weights))")
    else
        console!(sh, "> fit_model(model, $(free), data; " *
                     "weights = $(weights), method = :$(opt))")
    end

    # `verb = true` throughout: the console under the results table shows what the optimiser
    # printed, so it has to print. The fit runs on a worker — see `GuiJob` — so the window
    # stays alive and its output can be read while it is produced.
    sh.fitlog = ""
    optsym = Symbol(opt)
    start_job!(sh, :fit, function (_stop)
        t0 = time()
        themap = nothing
        # The chi2 of the model as it stands, measured BEFORE the optimiser touches it and
        # with the fit's own weights, so "start" and "final" are the same quantity and their
        # ratio means something. Divided by the fit's own ndof below, not by one computed here.
        chi2_start = try
            # Qualified: `chi2_flat` is not exported, and an unqualified call here resolves to
            # nothing and is swallowed by this very catch.
            OITOOLS.chi2_flat(parse_model(md, free),
                              Float64[Float64(md[k]) for k in free], data; weights)
        catch err
            # Say so rather than only showing a dash. A starting model that cannot be
            # evaluated is worth knowing about before reading the fit that came out of it.
            console!(sh, "  starting chi2 not available: " * _cause(err))
            NaN
        end
        r = if opt == "grid"
            themap = chi2_map(md, free, data, gp1, gp2;
                              lb, ub, weights, n1 = Int(gridn),
                              n2 = gridn2 > 1 ? Int(gridn2) : Int(gridn))
            FitResult(themap)
        elseif opt == "lsqfit"
            fit_model_lsqfit(md, free, data; lb, ub, weights, constraints = cons,
                             maxIter = Int(maxeval), verb = true)
        elseif opt == "nested"
            # Whichever sampler is loaded. `nested_backend()` picks it; the panel says which,
            # so a logz in the fit table can be attributed to the sampler that produced it.
            fit_model_nested(md, free, data; lb, ub, weights, constraints = cons,
                             verb = true, cornerplot = false)
        else
            fit_model(md, free, data; lb, ub, weights, priors, constraints = cons,
                      method = optsym, maxeval = Int(maxeval), verb = true)
        end
        # The map travels with the result rather than being drawn here: `update_chi2_map!`
        # touches Makie, which belongs to the GUI thread.
        return (; result = r, map = themap, seconds = time() - t0, optimiser = opt, free,
                  chi2_start)
    end)
    return "fitting…"
end

"""
    finish_fit!(sh, res) -> String

Record a finished fit and draw whatever it produced. **GUI thread only**: `update_chi2_map!` is
a Makie call.
"""
function finish_fit!(sh::ShellState, res)
    if !res.ok
        msg = "! fit failed: " * _cause(res.err)
        sh.fitlog = isempty(res.output) ? msg : res.output * "\n" * msg
        console!(sh, msg)
        return msg
    end
    f = res.result
    r = f.result
    # The grid search prints nothing -- it is a comprehension over `chi2_flat`, not an optimiser
    # with a trace -- so the console gets the one thing worth saying about it.
    sh.fitlog = isempty(strip(res.output)) && f.optimiser == "grid" ?
        "grid search evaluates χ² on a fixed grid and has no iteration trace to report.\n" *
        "the surface itself is the output: see the χ² map below." :
        res.output

    # The map belongs to the fit that produced it, not to the shell: two grid fits over
    # different parameter pairs both stay reachable, and selecting an older fit can redraw it.
    f.map === nothing || sh.chi2map === nothing || update_chi2_map!(sh.chi2map, f.map)
    push!(sh.fits, (; result = r, optimiser = f.optimiser, seconds = f.seconds,
                      free = f.free, map = f.map, chi2_start = f.chi2_start))
    for (k, v) in zip(r.list_free_params, r.x_opt)
        console!(sh, @sprintf("    %-22s = %.6g", k, v))
    end
    it, ev = _fit_work(r)
    start  = isfinite(f.chi2_start) ? f.chi2_start / max(r.ndof, 1) : NaN
    line = @sprintf("chi2r %s -> %.4f   ndof %d   %d free   %s   %.2f s",
                    isfinite(start) ? @sprintf("%.4f", start) : "?",
                    r.chi2r, r.ndof, length(f.free),
                    it == 0 && ev == 0 ? "" :
                    it > 0 && ev > 0   ? string(it, " iters / ", ev, " evals") :
                    it > 0             ? string(it, " iters") : string(ev, " evals"),
                    f.seconds)
    console!(sh, "  " * line)
    return line
end

"""
    shell_chi2_map_info() -> String

`p1\tp2\tbest1\tbest2\tchi2r` for the most recent grid fit, or `""` when the last fit was not
one. QML uses it to decide whether the map panel has anything to show.
"""
function shell_chi2_map_info()
    sh = _shell()
    isempty(sh.fits) && return ""
    f = last(sh.fits)
    m = hasproperty(f, :map) ? f.map : nothing
    m === nothing && return ""
    b = argmin(m.chi2)
    return join((m.p1, m.p2, string(m.v1[b[1]]), string(m.v2[b[2]]),
                 string(minimum(m.chi2) / m.ndof)), "\t")
end

"""
    shell_free_names() -> String

The free parameters, newline-separated, in fit-vector order — the two axis pickers for the
grid choose from these, and only these, because `chi2_map` maps free parameters.
"""
function shell_free_names()
    m = _model()
    isempty(m.free) && return ""
    return join(m.free, "\n")
end

"The fitted values as `key=value` lines, so the table can adopt them."
function shell_fit_values()
    sh = _shell()
    isempty(sh.fits) && return ""
    r = last(sh.fits).result
    return join((string(k, "=", v) for (k, v) in zip(r.list_free_params, r.x_opt)), "\n")
end

"""
    _fit_params(r) -> String

A fit's optimised parameters as `name = value` pairs, with an uncertainty where the fitter
produced one.

The point of a row per fit is comparing fits, and two fits that differ only in χ² are two
numbers, not two models. What separates them is where the parameters landed.

`± σ` appears only when it was actually computed: Levenberg-Marquardt carries `stderror` from
the Jacobian, nested sampling has a posterior to take a standard deviation of, and the plain
NLopt methods have neither. Writing `±` on a number that came from nowhere would be worse than
leaving it off.
"""
function _fit_params(r)
    names = r.list_free_params
    vals  = r.x_opt
    # `stderror` on the LM result; a posterior column-wise for nested sampling.
    errs = if hasproperty(r, :stderror) && r.stderror !== nothing &&
              length(r.stderror) == length(vals)
        Float64.(r.stderror)
    elseif hasproperty(r, :posterior) && r.posterior isa AbstractMatrix &&
           size(r.posterior, 2) == length(vals)
        [std(@view r.posterior[:, j]) for j in eachindex(vals)]
    else
        nothing
    end
    parts = String[]
    for (j, n) in enumerate(names)
        v = @sprintf("%.5g", vals[j])
        push!(parts, errs === nothing || !isfinite(errs[j]) ? "$n = $v" :
                     "$n = $v ± " * @sprintf("%.3g", errs[j]))
    end
    return join(parts, "   ")
end

"""
    shell_fit_output() -> String

Everything the last model fit printed, for the console under the results table.

Same mechanism as the reconstructors' output; `_capture_stdout` documents why the obvious
approaches do not work inside a Qt event loop.
"""
shell_fit_output() = (sh = SHELL[]; sh === nothing ? "" : sh.fitlog)

"""
    _fit_work(result) -> (iterations, evaluations)

How much work a fit took, as whichever of the two counts its optimiser actually reports. `0`
means "this optimiser does not report that", not "zero".

They are different quantities and are not interchangeable: NLopt counts chi2 EVALUATIONS and
never reports iterations; Levenberg--Marquardt counts ITERATIONS, each of which costs a
residual and a Jacobian. Presenting either under the other's name would misstate the cost of
the fit, so both are carried and the panel labels whichever it has.
"""
function _fit_work(r)
    if r isa FitResult                      # NLopt, and the grid search, which reports its size
        return (0, r.n_evals)
    elseif r isa LsqFitResult
        # LsqFit reports no count of its own; the trace is stored so that its length is one.
        tr = try length(r.lsqfit_result.trace) catch; 0 end
        return (tr, 0)
    end
    return (0, 0)
end

"""
    shell_fit_rows() -> String

Completed fits as
`label\toptimiser\tchi2r_start\tchi2r\twork\tndof\tnfree\taic\tbic\tparams` rows, newest
first.

`chi2r_start` is the model as it stood when the fit was launched, measured with the fit's own
weights and divided by the fit's own `ndof`, so the two chi2r columns are the same quantity and
the fit's progress is the ratio between them. It is `NaN` when the starting model could not be
evaluated, which the panel renders as an em dash.

`work` is `_fit_work` rendered with its unit, because iterations and evaluations are different
quantities and only one of them exists for any given optimiser.

AIC and BIC are computed here rather than read off `chi2r`, because `ndof` counts DATA POINTS
and is not reduced by the free-parameter count — so a comparison that used it as a degrees-of-
freedom would favour the larger model for the wrong reason.
"""
function shell_fit_rows()
    sh = _shell()
    rows = String[]
    for (i, f) in enumerate(reverse(sh.fits))
        r = f.result
        k = length(f.free)
        aic = r.chi2 + 2k
        bic = r.chi2 + k * log(max(r.ndof, 1))
        # The starting chi2 divided by the FIT'S ndof, so the two chi2r columns are the same
        # quantity and the fit's progress is the ratio between them.
        start = isfinite(f.chi2_start) ? f.chi2_start / max(r.ndof, 1) : NaN
        it, ev = _fit_work(r)
        work = it > 0 && ev > 0 ? string(it, " / ", ev) :
               it > 0           ? string(it, " it")     :
               ev > 0           ? string(ev, " ev")     : "—"
        push!(rows, join((string("fit ", length(sh.fits) - i + 1), f.optimiser,
                          _numv(start), string(round(r.chi2r; digits = 4)),
                          work, string(r.ndof), string(k),
                          string(round(aic; digits = 1)), string(round(bic; digits = 1)),
                          _fit_params(r)), "\t"))
    end
    return join(rows, "\n")
end

# ── the live model ───────────────────────────────────────────────────────────
#
# ONE model, owned by the session, which the panel reads and edits through these calls. It
# replaces a round trip through text: the table used to serialise itself to `key=value` lines
# that the fitter re-parsed, so the dict the inspector described and the dict the fitter saw
# were two different objects that merely tended to agree. They cannot disagree now.

"The model the panel is editing, created empty on first use so there is always one."
function _model(sh::ShellState = _shell())
    isempty(sh.session.models) &&
        push!(sh.session.models, ModelEntry("untitled", Dict{String,Any}(), String[],
                                            Dict{String,Float64}(), Dict{String,Float64}(), Any[]))
    return sh.session.models[end]
end

"""
    _display_wl_mjd() -> (wl, mjd)

The wavelength and epoch a chromatic parameter is shown at: the loaded dataset's first, or
`nothing` when no dataset is loaded and there is nothing to be chromatic against.
"""
function _display_wl_mjd()
    e = current_dataset(_shell())
    e === nothing && return (nothing, nothing)
    d = e.data[1, 1]
    w = hasproperty(d, :uv_lam) && !isempty(d.uv_lam) ? Float64(minimum(d.uv_lam)) : nothing
    t = hasproperty(d, :uv_mjd) && !isempty(d.uv_mjd) ? Float64(minimum(d.uv_mjd)) : nothing
    return (w, t)
end

"""
    shell_model_rows() -> String

The parameter table, one row per line, as
`comp \t param \t key \t mode \t value \t expr \t lb \t ub \t fitindex \t atbound`.

Ten fields, no `kind`: that belongs to a COMPONENT, not to one of its parameters, and
`shell_model_components` reports it there along with the key that decided it.

Straight from `model_rows`, which resolves derived values and flags a parameter sitting on its
bound — so the panel reports what the PARSER made of the model rather than what was typed.
"""
function shell_model_rows()
    m = _model()
    isempty(m.dict) && return ""
    # `nothing`, not an empty Dict: `model_rows` reads a supplied dict as authoritative, so
    # handing it an empty one suppresses `default_bounds` and every bound comes back 0/0.
    # The dataset's first wavelength and epoch, so a chromatic expression resolves to the
    # number it actually takes there rather than to nothing.
    w, t = _display_wl_mjd()
    rows = try
        model_rows(m.dict, m.free;
                   lb = isempty(m.lb) ? nothing : m.lb,
                   ub = isempty(m.ub) ? nothing : m.ub,
                   wl = w, mjd = t)
    catch err
        console!(_shell(), "! model: " * _cause(err); kind = :err)
        return ""
    end
    return join((join((r.component, r.param, r.key, _mode_name(r.mode),
                       _numv(r.value), r.expr, _num(r.lb), _num(r.ub),
                       string(r.fitindex), string(r.atbound)), "\t")
                 for r in rows), "\n")
end

_mode_name(m) = m === PARAM_FREE ? "free" : m === PARAM_EXPR ? "expr" : "fixed"

# Bounds: a non-finite bound is no bound, and the panel leaves those fields blank.
_num(x) = isfinite(x) ? string(x) : "0.0"

# Values: a non-finite value is a value the resolver could not produce -- a chromatic
# expression with no wavelength to show it at, or a reference to a key that is not there. The
# panel renders NaN as an em dash; reporting it as 0.0 would put a number in the table that
# the model does not hold.
_numv(x) = isfinite(x) ? string(x) : "NaN"

"""
    shell_model_components() -> String

`name \t kind \t geometry_key \t nparams \t nunrecognised` per component, from
`model_inspection` — including the key that DECIDED the kind, which is the part that catches a
`gaussian_ring` silently becoming a plain Gaussian.
"""
function shell_model_components()
    m = _model()
    isempty(m.dict) && return ""
    insp = try
        model_inspection(m.dict)
    catch err
        console!(_shell(), "! inspection: " * _cause(err); kind = :err)
        return ""
    end
    return join((join((c.name, string(c.kind), c.geometry_key,
                       string(length(c.params)), string(length(c.unrecognised))), "\t")
                 for c in insp.components), "\n")
end

"""
    shell_model_warnings() -> String

What is wrong with the model as it stands, one per line, or `""` when nothing is.

The same checks `display_model` prints to a terminal — [`model_warnings`](@ref) is shared, so
the panel and a script cannot disagree about the same model.
"""
function shell_model_warnings()
    m = _model()
    isempty(m.dict) && return ""
    w = try
        model_warnings(m.dict, m.free; lb = m.lb, ub = m.ub)
    catch err
        return "could not check the model: " * _cause(err)
    end
    return join(w, "\n")
end

"""
    shell_model_inspection() -> String

Three lines: unrecognised keys, whether the resolver broadcasts, and the globals. Each is a
fact about the model the user cannot see any other way.
"""
function shell_model_inspection()
    m = _model()
    isempty(m.dict) && return "\n\n"
    insp = try
        model_inspection(m.dict)
    catch
        return "\n\n"
    end
    return join((join(insp.unrecognised, " "),
                 string(insp.broadcasting),
                 join(insp.globals, " ")), "\n")
end

"""
    shell_model_depends() -> String

Whether the model references `\$WL` and `\$MJD`, then the wavelength and epoch the panel
displays chromatic values at, as `wl<TAB>mjd<TAB>wl_value<TAB>mjd_value`. The first two are
`1` or `0`; the last two are empty when no dataset is loaded.

Two flags rather than one: the render panel offers a wavelength control and a time control,
and each should be live only if it changes the answer. A control that moves nothing reads as a
broken render rather than as a model that does not vary.
"""
function shell_model_depends()
    m = _model()
    isempty(m.dict) && return "0\t0\t\t"
    d = model_depends_on(m.dict)
    w, t = _display_wl_mjd()
    return string(d.wl ? 1 : 0, '\t', d.mjd ? 1 : 0, '\t',
                  w === nothing ? "" : string(w), '\t',
                  t === nothing ? "" : string(t))
end

"""
    shell_model_sed() -> String

Draw the model's SED over the loaded dataset's wavelengths, and return the summary.

The data's own wavelengths, not a grid invented here: the band that was observed is the band
the model was constrained over, and a curve drawn outside it is extrapolation dressed as a
measurement.

Refused for a model that references neither `\$WL` nor `\$MJD`. Such a model is the same at
every wavelength, and a panel of flat lines is not a spectrum.
"""
function shell_model_sed()
    sh = _shell()
    m  = _model(sh)
    e  = current_dataset(sh)
    e === nothing && return "! no dataset loaded"
    isempty(m.dict) && return "! no model"
    model_depends_on(m.dict).wl ||
        return "! this model references no \$WL, so its SED is flat — nothing to plot"
    return try
        d  = e.data[1, 1]
        # Metres, the unit `\$WL` carries in the resolver and `uv_lam` carries in the data.
        # The panel's x axis is in microns, so the conversion happens at the drawing, not here.
        wl = sort(unique(Float64.(d.uv_lam)))
        length(wl) < 2 && return "! this dataset has one wavelength, so there is no spectrum"
        fm = parse_model(Dict{String,Any}(String(k) => v for (k, v) in m.dict),
                         String.(m.free); nB_workspace = 1)
        x = Float64[_row_value(m, k) for k in m.free]
        total, comps = model_to_sed(fm, x, wl)
        line = update_sed!(sh.sedplot, wl, total, comps)
        console!(sh, "model_to_sed(model, x, wl_grid)"; kind = :cmd)
        console!(sh, line)
        line
    catch err
        msg = "! SED not available: " * _cause(err)
        console!(sh, msg; kind = :err)
        msg
    end
end

"""
    shell_model_residuals() -> String

Draw the current model's normalised residuals and return the one-line summary.

The model as it stands, not the last fit: after "Adopt" the two are the same, and before a fit
the residuals of a hand-built model are worth seeing -- that is how you find out whether the
starting point is anywhere near the data before spending an optimiser on it.

Built with `parse_model` and `model_to_residuals`, the same pair `shell_model_chi2` uses, so
the picture and the chi2 in the corner are the same model evaluated the same way.

`quiet` suppresses the console lines. The panel redraws itself after every edit while it is on
screen -- residuals of the model as it was two keystrokes ago would be the wrong picture -- and
a transcript line per keystroke is not a transcript.
"""
function shell_model_residuals(quiet::Bool = false)
    sh = _shell()
    r = _current_residuals(sh; quiet)
    r isa String && return r
    line = update_residuals!(sh.residplot, r.data, r.res)
    if !quiet
        console!(sh, "model_to_residuals(model, x, data)"; kind = :cmd)
        console!(sh, line)
    end
    return line
end

"""
    _current_residuals(sh) -> (; data, res) or an error String

Evaluate the model in the table against the loaded data.

Shared by the panel and by "Save PNG", so the file and the screen show the same numbers rather
than two evaluations that could have drifted apart between them.
"""
function _current_residuals(sh::ShellState; quiet::Bool = false)
    m = _model(sh)
    e = current_dataset(sh)
    e === nothing && return "! no dataset loaded"
    isempty(m.dict) && return "! no model"
    return try
        d  = e.data[1, 1]
        fm = parse_model(Dict{String,Any}(String(k) => v for (k, v) in m.dict),
                         String.(m.free); nB_workspace = 1)
        x = Float64[_row_value(m, k) for k in m.free]
        (; data = d, res = model_to_residuals(fm, x, d))
    catch err
        msg = "! residuals not available: " * _cause(err)
        quiet || console!(sh, msg; kind = :err)
        msg
    end
end

"""
    shell_model_chi2() -> String

Reduced chi2 of the model as it stands, against the current dataset, or `""` when either is
missing. Cheap for analytic components, which is what makes it worth recomputing on every edit.
"""
function shell_model_chi2()
    sh = _shell()
    m = _model(sh)
    e = current_dataset(sh)
    (e === nothing || isempty(m.dict)) && return ""
    return try
        d = e.data[1, 1]
        # `parse_model`, the same call `render_model_image` makes: it takes the free list too,
        # because which parameters are free changes what the flat model IS.
        fm = parse_model(Dict{String,Any}(String(k) => v for (k, v) in m.dict),
                         String.(m.free); nB_workspace = 1)
        x = Float64[_row_value(m, k) for k in m.free]
        chi2 = model_to_chi2(fm, x, d)
        @sprintf("%.4f", chi2 / max(1, _chi2_ndof(d)))
    catch err
        # A model mid-edit is not an error, it is simply not evaluable yet — but say so once,
        # or a permanently blank chi2 looks like the number rather than the model being absent.
        console!(_shell(), "  chi2 not available: " * _cause(err))
        ""
    end
end

_row_value(m, key) = get(m.dict, key, 0.0) isa Real ? Float64(m.dict[key]) : 0.0

"Data points the chi2 is taken over. Counts VALID points, which is what `ndof` means here."
_chi2_ndof(d) = d.nv2 + d.nt3amp + d.nt3phi + d.nvisamp + d.nvisphi

"""
    shell_set_param(key, field, value) -> String

Edit one cell of the model: `field` is `value`, `mode`, `lb`, `ub` or `expr`.

`mode` is where the exclusions live. A parameter listed in `list_free_params` MUST be numeric —
the resolver raises on a string there — so going free drops any expression, and going derived
drops it from the fit vector. Doing that here rather than in QML keeps the rule in one place.
"""
function shell_set_param(key::AbstractString, field::AbstractString, value)
    sh = _shell(); m = _model(sh); k = String(key)
    haskey(m.dict, k) || return "! no such parameter: " * k
    f = String(field)
    if f == "value"
        m.dict[k] = Float64(value)
    elseif f == "expr"
        m.dict[k] = String(value)
        filter!(!=(k), m.free)                 # derived and free are exclusive
    elseif f == "lb"
        m.lb[k] = Float64(value)
    elseif f == "ub"
        m.ub[k] = Float64(value)
    elseif f == "mode"
        mode = String(value)
        if mode == "free"
            m.dict[k] isa Real || (m.dict[k] = 0.0)   # a free parameter must be numeric
            k in m.free || push!(m.free, k)
        else
            filter!(!=(k), m.free)
        end
    else
        return "! unknown field: " * f
    end
    return ""
end

"""
    shell_open_model(path)   -> String
    shell_save_model(path)   -> String
    shell_import_pmoired(path) -> String
    shell_export_pmoired(path) -> String

Model files. TOML carries the model, the free list, the bounds, the constraints and the priors
together; PMOIRED carries the model alone, and export warns about what it cannot represent
rather than writing something that means a different thing on the other side.
"""
function shell_open_model(path::AbstractString)
    sh = _shell()
    r = try
        read_model_file(String(path))
    catch err
        msg = "! could not read model: " * _cause(err); console!(sh, msg; kind = :err); return msg
    end
    m = _model(sh)
    m.dict = r.model; m.free = collect(r.free)
    m.lb = Dict{String,Float64}(r.lb); m.ub = Dict{String,Float64}(r.ub)
    m.name = isempty(r.name) ? basename(String(path)) : r.name
    console!(sh, "> read_model_file(\"" * String(path) * "\")")
    console!(sh, "  $(m.name): $(length(m.dict)) keys, $(length(m.free)) free")
    return m.name
end

function shell_save_model(path::AbstractString)
    sh = _shell(); m = _model(sh)
    isempty(m.dict) && return "! no model to save"
    try
        write_model_file(String(path), m.dict; free = m.free, lb = m.lb, ub = m.ub, name = m.name)
    catch err
        msg = "! could not write model: " * _cause(err); console!(sh, msg; kind = :err); return msg
    end
    console!(sh, "> write_model_file(\"" * String(path) * "\", model; free = $(m.free))")
    return ""
end

function shell_import_pmoired(path::AbstractString)
    sh = _shell()
    md = try
        pmoired_to_dict(read(String(path), String))
    catch err
        msg = "! PMOIRED import: " * _cause(err); console!(sh, msg; kind = :err); return msg
    end
    m = _model(sh)
    m.dict = md; m.free = String[]; empty!(m.lb); empty!(m.ub)
    m.name = basename(String(path))
    console!(sh, "> pmoired_to_dict(read(\"" * String(path) * "\", String))")
    # Import is a transpiler, not a validator, so say at once what the parser made of it.
    insp = try model_inspection(md) catch; nothing end
    insp === nothing || isempty(insp.unrecognised) ||
        console!(sh, "  ! keys the parser does not recognise: " *
                     join(insp.unrecognised, ", "); kind = :err)
    return m.name
end

function shell_export_pmoired(path::AbstractString)
    sh = _shell(); m = _model(sh)
    isempty(m.dict) && return "! no model to export"
    try
        dict_to_pmoired_file(String(path), m.dict)
    catch err
        msg = "! PMOIRED export: " * _cause(err); console!(sh, msg; kind = :err); return msg
    end
    console!(sh, "> dict_to_pmoired_file(\"" * String(path) * "\", model)")
    return ""
end

"""
    shell_model_image(model_json, nx, pixsize, wl, mjd) -> String

Render the current model and put it on the Model perspective's canvas.

`model_json` is `key=value` pairs, one per line — the model dict as the panel holds it. A value
that parses as a number is one; anything else is an expression, which is exactly the
distinction the parser itself makes.

`wl` in metres and `mjd` in days; either ≤ 0 means "not specified", which a model that does not
reference them ignores anyway.
"""
function shell_model_image(model_json::AbstractString, nx::Integer, pixsize::Real,
                           wl::Real, mjd::Real)
    sh = _shell()
    md = Dict{String,Any}()
    for line in split(String(model_json), "\n")
        isempty(strip(line)) && continue
        i = findfirst('=', line)
        i === nothing && continue
        k = strip(line[1:i-1]); v = strip(line[nextind(line, i):end])
        md[String(k)] = something(tryparse(Float64, v), String(v))
    end
    isempty(md) && return "! no model"

    r = try
        render_model_image(md, String[], Float64[];
                           nx = Int(nx), pixsize = Float64(pixsize),
                           wl  = wl  > 0 ? Float64(wl)  : nothing,
                           mjd = mjd > 0 ? Float64(mjd) : nothing)
    catch err
        msg = "! render failed: " * _cause(err); console!(sh, msg); return msg
    end

    sh.modelcanvas === nothing || show_image!(sh.modelcanvas, r.image, r.pixsize)
    console!(sh, "> model_to_image(model, x; nx = $(r.nx), pixsize = $(r.pixsize)" *
                 (wl > 0 ? ", wl = $(wl)" : "") * ")")
    dep = r.depends.wl || r.depends.mjd ? "" :
          "   (this model uses neither \$WL nor \$MJD, so λ and MJD change nothing)"
    line = "rendered $(r.nx)×$(r.nx), FOV $(round(r.fov; digits=2)) mas, flux $(round(r.flux; digits=4))" * dep
    console!(sh, "  " * line)
    return line
end

"""
    shell_sim_source(kind, path) -> String

Validate the sky a simulation would use: `"ok\t<summary>"` or `"bad\t<why>"`.

Checked when the file is chosen rather than when Simulate is pressed, so a 3-D cube picked as a
grey image is caught before an observation is computed from it.
"""
function shell_sim_source(kind::AbstractString, path::AbstractString)
    i = simulate_source_info(String(kind), String(path))
    return (i.ok ? "ok\t" : "bad\t") * i.summary
end

"""
    shell_sim_band(setup) -> String

`band\tλ_µm` for a spectral setup — the photometric band its channels fall in, and their mean
wavelength. `"\t"` when the setup is unknown.

Which band an observation is made in is a property of the instrument, not of the target, so the
panel should not make the user look it up: MIRC-X is H, MYSTIC is K, and a magnitude typed for
the wrong one feeds a noise model that is wrong by whole magnitudes.

The band comes from the MEAN of the setup's channels through `band_for_wavelength`, not from
the setup's name. Names are not a reliable guide -- `MIRCX_LOWJ` and `MIRCX_LOWH` differ by one
letter and two bands -- while the wavelengths in the file are the thing itself.
"""
function shell_sim_band(setup::AbstractString)
    isempty(strip(setup)) && return "\t"
    # Only the FILE READ is guarded: an unknown setup is an ordinary answer, and the panel
    # shows it as "no band". Everything after it is arithmetic on a WaveConfig, and wrapping
    # that in a catch too is how a plain typo -- `w.lambda` for `w.λ` -- came back as "unknown
    # setup" for every input instead of failing where it was written.
    w = try
        read_wave_file(String(strip(setup)))
    catch
        return "\t"
    end
    λ = sum(w.λ) / length(w.λ)
    (isfinite(λ) && λ > 0) || return "\t"
    return string(band_for_wavelength(λ).name, '\t', round(λ * 1e6; digits = 3))
end

"""
    shell_simbad(name) -> String

Resolve one target: `"ok\tra\tdec\tV\tJ\tH\tK\tmain_id\tsptype"`, or `"bad\t<why>"`.

Coordinates in DEGREES. A magnitude SIMBAD does not have comes back empty rather than as a
number, because a missing magnitude and a magnitude of zero are very different things to an
SNR estimate — Vega is roughly zero in every band.

Blocking: this is a network call on Qt's thread, and the caller is expected to have said so
before making it.
"""
function shell_simbad(name::AbstractString)
    sh = SHELL[]
    say(t) = sh === nothing ? nothing : console!(sh, t)
    n = strip(String(name))
    isempty(n) && return "bad\tno name"
    say("> simbad_target(\"$n\")")
    t = try
        simbad_target(n)
    catch err
        msg = "bad\t" * _cause(err)
        say("! " * _cause(err))
        return msg
    end
    m(b) = (v = get(t.mags, b, NaN); isnan(v) ? "" : string(round(v; digits = 3)))
    say(@sprintf("  %s: ra %.5f  dec %+.5f  %s  plx %.2f mas",
                 t.main_id, t.ra, t.dec, t.sptype, t.plx))
    return join(("ok", string(round(t.ra; digits = 6)), string(round(t.dec; digits = 6)),
                 m("V"), m("J"), m("H"), m("K"), t.main_id, t.sptype), "\t")
end

"""
    shell_nested_backend() -> String

`backend\tlabel`: the nested sampler in force and how to name it in the panel, or `"\t"` when
neither extension is loaded.

The Modeling tab needs both halves. The key decides whether the entry is offered at all, and
the label has to name the sampler rather than saying "nested sampling", because two samplers
that disagree on an evidence are only distinguishable if the interface admits which one ran.
"""
function shell_nested_backend()
    b = nested_backend()
    b === :none && return "\t"
    label = b === :nautilus  ? "Nautilus.jl" :
            b === :ultranest ? "UltraNest"   : String(b)
    return string(b, '\t', label)
end

"""
    shell_plot_scale() -> String

`effective\tuser\tmarker\tzoom`: the plot scale actually in force, the user's override behind
it (zero when it is computed rather than chosen), the marker size (zero for the per-view
default) and the zoom factor per wheel detent.

Two scale numbers rather than one because they answer different questions -- the panel opens
showing what is being drawn, while "Save defaults" stores the override, so a scale that was
computed from this screen's DPI is not pinned onto the next screen.
"""
shell_plot_scale() = join((live_plot_scale(), PLOT_SCALE_USER[], MARKER_SIZE_USER[],
                           zoom_per_detent()), '\t')

"""
    shell_set_plot_scale(x)  -> String
    shell_set_marker_size(x) -> String

Apply a value from the settings panel and redraw. Zero restores the computed default.

Both redraw immediately rather than waiting for the next action: the panel is being used to
judge a size by eye, so the plot has to change while the user is looking at it.
"""
function shell_set_plot_scale(x)
    v = set_plot_scale!(Float64(x))
    sh = SHELL[]
    sh === nothing && return ""
    console!(sh, @sprintf("plot scale: %.2f%s", live_plot_scale(), v == 0 ? " (auto)" : ""))
    current_dataset(sh) === nothing || try
        sh.status = refresh_plot!(sh)
    catch err
        console!(sh, "could not restyle: " * _cause(err); kind = :err)
    end
    return ""
end

"""
    shell_set_zoom_step(x) -> String

How far one wheel detent zooms. Zero restores the default.

No redraw: this changes what the NEXT scroll does, not what is on screen, and redrawing to
show a setting that has not been used yet would be a lie about the current view.
"""
function shell_set_zoom_step(x)
    v = set_zoom_step!(Float64(x))
    sh = SHELL[]
    sh === nothing && return ""
    console!(sh, @sprintf("zoom step: %.2fx per detent%s", zoom_per_detent(),
                          v == 0 ? " (default)" : ""))
    return ""
end

function shell_set_marker_size(x)
    v = set_marker_size!(Float64(x))
    sh = SHELL[]
    sh === nothing && return ""
    console!(sh, v == 0 ? "marker size: auto" : @sprintf("marker size: %.0f pt", v))
    current_dataset(sh) === nothing || try
        sh.status = refresh_plot!(sh)
    catch err
        console!(sh, "could not restyle: " * _cause(err); kind = :err)
    end
    return ""
end

"""
    gui_settings_file() -> String

Where the appearance defaults live: `oitools/gui.toml` under the platform's per-user config
directory.

Deliberately not `LocalPreferences.toml`: these are the user's window, not the project's, and
they should follow the user across every project that launches the GUI.
"""
function gui_settings_file()
    base = Sys.iswindows() ? get(ENV, "APPDATA", homedir()) :
           get(ENV, "XDG_CONFIG_HOME", joinpath(homedir(), ".config"))
    return joinpath(base, "oitools", "gui.toml")
end

"""
    shell_save_settings(payload) -> String

Write the settings panel's values to `gui_settings_file()`, and say where they went.

`payload` is `key\tvalue` lines, the same tab-separated shape the rest of this bridge uses.
Values are stored as they arrive: QML owns what the keys mean, and a key this version does not
recognise is kept rather than dropped, so an older build cannot silently erase a newer one's
settings.
"""
function shell_save_settings(payload)
    sh = SHELL[]
    path = gui_settings_file()
    try
        d = isfile(path) ? TOML.parsefile(path) : Dict{String,Any}()
        for line in split(String(payload), '\n')
            isempty(line) && continue
            f = split(line, '\t')
            length(f) == 2 || continue
            v = tryparse(Float64, f[2])
            d[f[1]] = v === nothing ? f[2] : v
        end
        mkpath(dirname(path))
        open(io -> TOML.print(io, d), path, "w")
        sh === nothing || console!(sh, "saved appearance defaults to " * path)
        return path
    catch err
        sh === nothing || console!(sh, "could not save settings: " * _cause(err); kind = :err)
        return ""
    end
end

"""
    shell_load_settings() -> String

Read the saved defaults back as `key\tvalue` lines, empty if none have been saved.

Called from QML before the first plot is drawn. A missing or unreadable file is not an error:
the built-in defaults are a perfectly good answer, and a corrupt settings file should not stop
the window from opening.
"""
function shell_load_settings()
    path = gui_settings_file()
    isfile(path) || return ""
    try
        d = TOML.parsefile(path)
        return join((string(k, '\t', v) for (k, v) in sort(collect(d), by = first)), '\n')
    catch
        return ""
    end
end

"""
One selectable component, as the panel offers it.

`id` is what `shell_add_component` takes and is NOT always a parser kind: "ring, custom I(r)"
and "custom image function" are both `:hankel`, differing only in which keys they seed, and the
five limb-darkening laws are five kinds under one category. The parser knows kinds; the panel
offers shapes, and the two are not the same list.
"""
struct ComponentChoice
    category :: String   # grouping key
    catlabel :: String   # how the group is named
    id       :: String   # what shell_add_component takes
    label    :: String   # how the entry is named within its group
    name     :: String   # seed for the component's name
end

"""
The component menu, in order.

Categories are alphabetical, and a category with one entry is offered as itself rather than as
a menu of one. Grouping follows what the shapes ARE rather than how the parser spells them: a
uniform disc and a limb-darkened disc are both discs and differ only in the law across the
face, which is the sub-choice; a ring is a ring whether its cross-section is flat, Gaussian or
written out.
"""
const COMPONENT_CHOICES = ComponentChoice[
    # Simplest first, then by how much structure each adds, with the two written-out kinds
    # last: those are the ones you reach for when none of the named shapes fits, so they are
    # the ones you look for deliberately rather than land on by scrolling.
    ComponentChoice("point",    "point source", "point", "point source", "pt"),

    ComponentChoice("disc", "disc", "ud",        "uniform",                    "star"),
    ComponentChoice("disc", "disc", "ldlin",     "linear limb darkening",      "star"),
    ComponentChoice("disc", "disc", "ldquad",    "quadratic limb darkening",   "star"),
    ComponentChoice("disc", "disc", "ldsqrt",    "square-root limb darkening", "star"),
    ComponentChoice("disc", "disc", "ldpow",     "power-law limb darkening",   "star"),
    ComponentChoice("disc", "disc", "ldclaret4", "Claret four-parameter",      "star"),

    # The three ring cross-sections. `ring_profile` is the same `:hankel` machinery as the
    # custom image function, seeded with `diamin`/`diamout` so the profile is evaluated on the
    # annulus between them -- which is a thin ring smeared by I(r), since a Hankel transform is
    # a superposition of thin rings. PMOIRED builds its `profile` rings on exactly that grid.
    ComponentChoice("ring", "ring", "ring",          "uniform",     "ring"),
    ComponentChoice("ring", "ring", "gaussian_ring", "Gaussian",    "ring"),
    ComponentChoice("ring", "ring", "ring_profile",  "custom I(r)", "ring"),

    ComponentChoice("gaussian", "Gaussian", "gaussian", "Gaussian", "gauss"),
    ComponentChoice("resolved", "resolved flux", "resolved", "resolved flux", "bg"),
    ComponentChoice("crescent", "crescent", "crescent", "crescent", "crescent"),

    # You write I(r); the Hankel transform turns it into V(B). The image is what is specified,
    # which is what the name says.
    ComponentChoice("image", "custom image function", "image_func",
                    "custom image function", "disk"),

    # You write V(B) directly, in $B. The other direction, and the only kind whose expression
    # is not an image at all.
    ComponentChoice("visfunc", "custom visibility function", "vis_func",
                    "custom visibility function", "comp"),
]

"""
    shell_component_kinds() -> String

What "+ component" offers, as `category\tcatlabel\tid\tlabel\tname` lines, in menu order.

`COMPONENT_CHOICES` is the table; the panel groups consecutive lines sharing a `category` and
offers a group of one as a plain entry rather than as a menu of one.

The last field seeds the name. A component name is the prefix of every key it owns and of
every `\$name,suffix` an expression refers to, so it is read far more often than it is typed;
`star` and `ring` say what the model is where `c1` and `c2` say only how many there are.
"""
function shell_component_kinds()
    return join((join((c.category, c.catlabel, c.id, c.label, c.name), '\t')
                 for c in COMPONENT_CHOICES), "\n")
end

# Starting values for a freshly added component. Not physics -- somewhere sane for a fit to
# begin, in the spirit of `default_bounds`. Sizes in mas; the outer of a pair is put clear of
# the inner so a new ring is a ring rather than an error.
const _COMPONENT_SEEDS = Dict{String,Float64}(
    "diamout" => 4.0, "fwhmout" => 4.0, "crout" => 4.0,
    # A crescent, not an off-centre hole. `croff` is the fraction of the way to internal
    # tangency, so 0.5 in a 2-to-4 mas annulus moves the hole halfway and still reads as a
    # ring; 0.8 in a 3-to-4 one opens the limb and looks like the thing it is named after.
    "crin" => 3.0, "croff" => 0.8, "crprojang" => 0.0,
    "u" => 0.3, "w" => 0.1, "alpha" => 0.15, "resolved" => 1.0,
    "c1" => 0.5, "c2" => -0.2, "c3" => 0.3, "c4" => -0.1)

_component_seed(suffix) = get(_COMPONENT_SEEDS, suffix, 2.0)

"""
    _share_flux!(m, newkey) -> Float64

The flux fraction a component being added should start at, having rescaled the others so the
fractions still sum to 1.

`display_model` warns when numeric flux fractions do not sum to 1, and every component starting
at 1.0 is the commonest way to earn that warning: two components then each carry a whole
source, which inflates chi2 rather than splitting the light between them. Adding a second
component means "share the source", not "double it".

Existing fractions keep their ratio to each other -- a deliberate 4:1 stays 4:1 -- and the new
one takes an equal share. Nothing is touched when any existing fraction is an expression:
renormalising around a value the resolver has yet to compute would be a guess, and the warning
is the better thing to see.
"""
function _share_flux!(m::ModelEntry, newkey::AbstractString)
    others = [k for k in keys(m.dict) if endswith(k, ",f") && k != newkey]
    isempty(others) && return 1.0
    all(m.dict[k] isa Real for k in others) || return 1.0
    share = 1.0 / (length(others) + 1)
    tot = sum(Float64(m.dict[k]) for k in others)
    for k in others
        m.dict[k] = tot > 0 ? (1.0 - share) * Float64(m.dict[k]) / tot : share
    end
    return share
end

"""
    shell_add_component(name, kind) -> String

Add a component of `kind`, with the keys the parser identifies that kind by, and a flux
fraction. Returns `""` on success or a message beginning with `!`.

The geometry key is written even when the kind consumes no geometry parameters: `_identify_kind`
looks the key up by name, so a `resolved` component without a `resolved` key is not a resolved
component -- it is a parse error.
"""
function shell_add_component(name, kind)
    sh = _shell(); m = _model(sh)
    nm = strip(String(name))
    isempty(nm) && return "! a component needs a name"
    occursin(',', nm) && return "! a component name cannot contain a comma"
    nm == GLOBAL_COMPONENT && return "! that name is reserved for globals"
    k = Symbol(String(kind))
    k === :ld && return "! choose a limb-darkening law: " *
                        join((l.key for l in LD_LAWS), ", ")
    # Panel ids that are not parser kinds. `image_func` and `ring_profile` are both `:hankel`
    # and differ only in the keys seeded below: the first is a profile about the centre, the
    # second is one evaluated on an annulus, which is the same machinery pointed at a ring.
    k === :hankel || k === :image_func || k === :ring_profile || k === :vis_func ||
        haskey(_ANALYTIC_PARAM_SUFFIXES, k) ||
        return "! unknown component kind: " * String(kind)
    pre = nm * ","
    any(startswith(key, pre) for key in keys(m.dict)) &&
        return "! there is already a component called " * nm

    if k === :vis_func
        # V(B) written out. The seed is a uniform disc's own visibility, which is a real model
        # rather than a placeholder: it evaluates, it fits, and it is the thing a user is most
        # likely to be about to modify.
        # `ifelse`, because V is evaluated at B = 0: the zero-baseline point IS the total flux,
        # and `chi2_flat`'s flux term and `model_to_sed` both ask for it. 2J₁(x)/x is 0/0 there
        # and a NaN at one uv point poisons the whole chi2, so the seed shows the guard a
        # written-out visibility generally needs.
        m.dict[pre * "visfunc"] =
            "ifelse(\$B < 1e-8, 1.0, 2*besselj1(pi*\$d*\$B/206.264806) / " *
            "(pi*\$d*\$B/206.264806))"
        m.dict[pre * "d"]       = _component_seed("ud")
        m.dict[pre * "f"]       = _share_flux!(m, pre * "f")
        push!(m.free, pre * "d")
        console!(sh, "  added $(nm): V(B) written in \$B (Mλ); $(pre)d is free")
        console!(sh, "  a component's flux fraction assumes V(0) = 1, and V IS evaluated at " *
                     "B = 0; guard any 0/0 there as the seed does")
        return ""
    end

    if k === :ring_profile
        # A ring whose cross-section you write. The grid runs over the annulus rather than from
        # the centre, so `$R` is a radius within the ring and a flat profile is a uniform ring
        # -- exactly the grid PMOIRED builds for its own `profile` rings.
        din, dout = _component_seed("diamin"), _component_seed("diamout")
        m.dict[pre * "profile"] = "1.0"
        m.dict[pre * "diamin"]  = din
        m.dict[pre * "diamout"] = dout
        m.dict[pre * "nr"]      = 100.0
        m.dict[pre * "f"]       = _share_flux!(m, pre * "f")
        push!(m.free, pre * "diamin"); push!(m.free, pre * "diamout")
        console!(sh, "  added $(nm): ring cross-section on r ∈ [" * string(din/2) * ", " *
                     string(dout/2) * "] mas, 100 points; both diameters are free")
        return ""
    end

    if k === :hankel || k === :image_func
        # A profile has no fixed parameter list, so nothing here writes one: what goes in is the
        # expression and the r grid, and `compile_profile` discovers the rest from the string as
        # it is edited.
        #
        # The seed is a flat profile, which is a uniform disk of diameter `udout`. It parses,
        # evaluates and fits, where an empty string or a half-written formula does none of the
        # three -- and a component that cannot be evaluated cannot be previewed either, so
        # there would be nothing on screen to edit against.
        m.dict[pre * "profile"] = "1.0"
        m.dict[pre * "udout"]   = _component_seed("diamout")
        m.dict[pre * "nr"]      = 100.0
        m.dict[pre * "f"]       = _share_flux!(m, pre * "f")
        # `udout` and not `profile`: a free parameter must be numeric, and the profile is a
        # string. `udout` is the one the grid is built from and the one worth fitting.
        push!(m.free, pre * "udout")
        console!(sh, "  added $(nm): radial profile on r ∈ [0, " *
                     string(_component_seed("diamout") / 2) * "] mas, 100 points; " *
                     "$(pre)udout is free")
        return ""
    end

    gk = get(_KIND_GEOMETRY, k, "")
    suffixes = copy(_ANALYTIC_PARAM_SUFFIXES[k])
    isempty(gk) || gk in suffixes || pushfirst!(suffixes, gk)
    for s in suffixes
        m.dict[pre * s] = _component_seed(s)
    end
    # Every component carries a flux fraction: it is what `:point` is identified BY, and what
    # makes a second component mean anything relative to the first.
    m.dict[pre * "f"] = _share_flux!(m, pre * "f")
    # The geometry parameter starts free, because it is the one the component was added to
    # measure. A component whose every parameter is fixed contributes nothing to the fit
    # vector, so adding one and pressing Run would fit exactly what was already there.
    # Bounds need no seeding here: `model_rows` falls back to `default_bounds` for any key
    # without one.
    isempty(gk) || push!(m.free, pre * gk)
    console!(sh, "  added $(nm): $(k), keys " * join((pre * s for s in suffixes), ", ") *
                 (isempty(gk) ? "" : "; $(pre * gk) is free"))
    return ""
end

"""
    shell_expression_keywords(key) -> String

The implicit variables an expression for `key` may use, as `name\tmeaning\tenabled\treason`.

They are SCOPED, and nothing else says so. `\$R` and `\$MU` are substituted by
`compile_profile` and exist only inside a profile string; `\$WL` and `\$MJD` come from the uv
points and exist only outside one, because a profile is evaluated on the radial grid before any
uv point is in hand. Crossing the boundary is not caught when it is typed -- it surfaces later
as `UndefVarError: WL not defined in OITOOLS`, which names neither the rule nor the expression
that broke it, and which a plain typo produces too.

So the four are always listed and the out-of-scope pair is disabled with the reason, rather
than being absent: a keyword that is simply missing reads as one that does not exist.
"""
function shell_expression_keywords(key::AbstractString)
    inprofile = endswith(String(key), ",profile")
    why_out = inprofile ?
        "not in a profile: the profile is evaluated on the radial grid, before any uv point" :
        "only inside a profile string, where compile_profile substitutes it"
    rows = (("R",    "radius over the profile grid, mas",            inprofile),
            ("MU",   "sqrt(1 - (r/r_max)^2), for limb darkening",    inprofile),
            ("D",    "diameter at each grid radius, i.e. 2\$R",      inprofile),
            ("RMIN", "the grid's inner radius, mas",                 inprofile),
            ("RMAX", "the grid's outer radius, mas",                 inprofile),
            ("DMIN", "the grid's inner diameter, mas",               inprofile),
            ("DMAX", "the grid's outer diameter, mas",               inprofile),
            ("WL",   "wavelength in METRES, one per uv point",       !inprofile),
            ("MJD",  "Modified Julian Date, one per uv point (~5.6 min resolution " *
                     "unless the file was read with T = Float64)",   !inprofile))
    return join((string("\$", n, '\t', m, '\t', ok ? "1" : "0", '\t', ok ? "" : why_out)
                 for (n, m, ok) in rows), "\n")
end

"""
    shell_check_expression(key, expr) -> String

Every `\$` reference in `expr`, classified, as `ref\tclass\tmessage`.

`class` is `ok`, `scope` or `unknown`. This is the check the parser does not do when the
expression is written: `dict_to_model` accepts all three happily and the mistake only appears
when something evaluates, as an `UndefVarError` naming an internal module. An out-of-scope
implicit variable and a misspelt parameter are indistinguishable there, and distinguishing them
is the whole point of doing it here.
"""
function shell_check_expression(key::AbstractString, expr::AbstractString)
    m = _model()
    k = String(key)
    inprofile = endswith(k, ",profile")
    comp = (i = findfirst(==(','), k)) === nothing ? "" : k[1:i-1]
    out = String[]
    for ref in OITOOLS.extract_refs(String(expr))
        if ref in OITOOLS.IMPLICIT_VARS
            # The grid variables belong to a profile; the uv-point ones cannot appear in one,
            # because a profile is evaluated on the radial grid before any uv point is in hand.
            gridvar = ref in ("R", "MU", "D", "RMIN", "RMAX", "DMIN", "DMAX")
            ok = gridvar ? inprofile : !inprofile
            push!(out, join((ref, ok ? "ok" : "scope",
                             ok ? "" : (gridvar ?
                                        "\$$(ref) only means anything inside a profile string" :
                                        "\$$(ref) is a uv-point quantity and is not available " *
                                        "inside a profile")), '\t'))
        elseif haskey(m.dict, ref)
            push!(out, join((ref, "ok", ""), '\t'))
        elseif !isempty(comp) && haskey(m.dict, comp * "," * ref)
            # A profile resolves a bare name against its own component first, then the globals
            # -- so inside `disk,profile`, `\$scale` means `disk,scale` if that exists.
            push!(out, join((ref, "ok", "resolves to " * comp * "," * ref), '\t'))
        else
            push!(out, join((ref, "unknown",
                             inprofile ?
                             "no parameter of this name: a profile treats every name that is " *
                             "not \$R or \$MU as one, so this would have to be created" :
                             "no parameter or global of this name"), '\t'))
        end
    end
    return join(out, "\n")
end

"""
    shell_profile_params(key, expr) -> String

The parameters a profile expression needs but the model does not have, as `name\tseed` lines.

`compile_profile` treats every name that is not `\$R` or `\$MU` as a parameter, resolved against
the component first and then the globals, so a profile is the one place where writing a name
CREATES the need for a parameter rather than referring to one. Listing them is what turns a
typo into a visible extra entry instead of an `UndefVarError` at the first evaluation.
"""
function shell_profile_params(key::AbstractString, expr::AbstractString)
    m = _model()
    k = String(key)
    endswith(k, ",profile") || return ""
    comp = k[1:findfirst(==(','), k)-1]
    out = String[]
    for ref in OITOOLS.extract_refs(String(expr))
        ref in OITOOLS.IMPLICIT_VARS && continue
        (haskey(m.dict, ref) || haskey(m.dict, comp * "," * ref)) && continue
        push!(out, string(ref, '\t', _component_seed(ref)))
    end
    return join(out, "\n")
end

"""
    shell_add_profile_params(key, expr) -> String

Create those parameters on the component, seeded. Returns `""` or a message beginning with `!`.

Component-qualified rather than global: `compile_profile` looks there first, and a global would
be silently shared with every other profile that happens to use the same name.
"""
function shell_add_profile_params(key::AbstractString, expr::AbstractString,
                                  seeds::AbstractString = "")
    sh = _shell(); m = _model(sh)
    k = String(key)
    endswith(k, ",profile") || return "! not a profile: " * k
    comp = k[1:findfirst(==(','), k)-1]
    # A template's own seeds win over `_component_seed`'s generic 2.0: `Rin = 1.26` describes
    # the shape the template is for, where 2.0 describes nothing and leaves the preview showing
    # something the user did not ask for.
    want = Dict{String,Float64}()
    for kv in split(seeds, ',')
        isempty(kv) && continue
        nv = split(kv, '=')
        length(nv) == 2 || continue
        v = tryparse(Float64, nv[2])
        v === nothing || (want[String(nv[1])] = v)
    end

    made = String[]
    for line in split(shell_profile_params(k, expr), "\n")
        isempty(line) && continue
        name, seed = split(line, '\t')
        nk = comp * "," * name
        m.dict[nk] = get(want, String(name), parse(Float64, seed))
        push!(made, nk)
    end
    isempty(made) && return ""
    console!(sh, "  created " * join(made, ", "))
    return ""
end

"""
    shell_profile_grid(component) -> String

The radial grid a profile is evaluated on, as `key\tr_min\tr_max\tnr`, or `!` and a reason.

Which key set `r_max` is worth saying out loud: four different ones can, they are tried in the
order `udout`, `diamout`, `diam`, `r_max`, and the first three are DIAMETERS that get halved
while the fourth is already a radius. A grid that came out twice the expected size is otherwise
a silent puzzle.
"""
function shell_profile_grid(component::AbstractString)
    m = _model(); cn = String(component)
    haskey(m.dict, cn * ",profile") || return "! " * cn * " is not a profile component"
    src, rmax = "", 0.0
    for (suffix, halve) in (("udout", true), ("diamout", true), ("diam", true), ("r_max", false))
        kk = cn * "," * suffix
        haskey(m.dict, kk) || continue
        v = try
            Float64(OITOOLS._resolve_numeric(kk, m.dict))
        catch
            return "! " * kk * " does not resolve to a number"
        end
        src, rmax = suffix, halve ? v / 2 : v
        break
    end
    isempty(src) && return "! no grid key: one of udout, diamout, diam or r_max is required"
    rmin = try
        Float64(OITOOLS._hankel_r_min(cn, m.dict))
    catch
        0.0
    end
    nr = try
        OITOOLS._hankel_Nr(cn, m.dict)
    catch
        100
    end
    # Whether the grid is big enough for the profile written on it, as I(r_max)/max(I).
    #
    # A profile is evaluated ONLY on [r_min, r_max]: everything beyond is not attenuated, it is
    # absent. So a grid that stops while the profile is still bright models a sharply cut ring
    # rather than the one that was typed, and nothing else on screen says so — the double
    # sigmoid at its published values sits at 54% of peak on the grid a new component starts
    # with. Reported rather than judged: the number says how bad it is, and 1% is only the
    # threshold at which it is worth mentioning.
    edge = 0.0
    expr = get(m.dict, cn * ",profile", nothing)
    if expr isa AbstractString && rmax > rmin
        try
            names = [r for r in OITOOLS.extract_refs(expr) if r ∉ OITOOLS.IMPLICIT_VARS]
            vals = Float64[]
            for n in names
                kk = haskey(m.dict, cn * "," * n) ? cn * "," * n : n
                push!(vals, Float64(OITOOLS._resolve_numeric(kk, m.dict)))
            end
            rr = collect(range(rmin, rmax; length = nr))
            mm = sqrt.(max.(0.0, 1 .- (rr ./ rmax) .^ 2))
            II = OITOOLS.compile_profile(expr, String.(names))(rr, mm, vals...)
            II isa AbstractVector || (II = fill(Float64(II), length(rr)))
            # Only when the profile is still FALLING at the edge. A flat profile is bright
            # all the way out and that is not a mistake -- it is a uniform disc, truncated by
            # design, and it is what a new component starts as. What matters is a shape cut off
            # mid-decline, which is the ring being modelled as something narrower than written.
            pk = maximum(abs, II)
            k = max(1, length(II) - 5)
            falling = length(II) > 5 && abs(II[end]) < abs(II[k])
            pk > 0 && falling && (edge = abs(II[end]) / pk)
        catch
            edge = 0.0            # will not evaluate: `shell_profile_curves` reports why
        end
    end
    return join((src, string(rmin), string(rmax), string(nr), string(edge)), '\t')
end

"""
    shell_profile_curves(component) -> String

Evaluate a profile and its Hankel transform onto the preview, returning `""` or a `!` message.

`hankel_transform` wants r in mas and B in cycles/mas, which is the same reciprocal pair, so
the only conversion here is from the baselines an observer reads -- Mλ -- into cycles/mas:
one cycle per radian is 4.8481e-9 cycles per mas.

Everything is read from the dict rather than from a parsed `FlatModel`, because the point is to
preview an expression WHILE it is being written, when the model as a whole may not yet parse.
"""
function shell_profile_curves(component::AbstractString)
    sh = _shell(); m = _model(sh)
    sh.profileplot === nothing && return ""          # headless: nothing to draw on
    cn = String(component)
    expr = get(m.dict, cn * ",profile", nothing)
    expr isa AbstractString || return "! " * cn * " is not a profile component"

    g = shell_profile_grid(cn)
    startswith(g, "!") && return g
    _, rmin_s, rmax_s, nr_s = split(g, '\t')
    rmin = parse(Float64, rmin_s); rmax = parse(Float64, rmax_s); nr = parse(Int, nr_s)
    rmax > rmin || return "! the grid is empty: r_max is not greater than r_min"

    names = [r for r in OITOOLS.extract_refs(expr) if r ∉ OITOOLS.IMPLICIT_VARS]
    vals = Float64[]
    for n in names
        k = haskey(m.dict, cn * "," * n) ? cn * "," * n : (haskey(m.dict, n) ? n : "")
        isempty(k) && return "! \$" * n * " has no parameter behind it yet"
        v = try
            Float64(OITOOLS._resolve_numeric(k, m.dict))
        catch err
            return "! " * k * ": " * _cause(err)
        end
        push!(vals, v)
    end

    r  = collect(range(rmin, rmax; length = nr))
    mu = sqrt.(max.(0.0, 1 .- (r ./ rmax) .^ 2))
    I = try
        OITOOLS.compile_profile(expr, String.(names))(r, mu, vals...)
    catch err
        return "! the profile will not evaluate: " * _cause(err)
    end
    # A profile that does not mention `$R` or `$MU` compiles to `@. <constant>`, which has no
    # array operand and so returns a SCALAR. That is a legitimate profile -- a flat one is a
    # uniform disk, and it is what a new component starts as -- so it is spread over the grid
    # rather than rejected. `hankel_vis` asserts on the length and would otherwise refuse it.
    I isa AbstractVector || (I = fill(Float64(I), length(r)))
    length(I) == length(r) ||
        return "! the profile returned $(length(I)) values for $(length(r)) grid points"
    all(isfinite, I) || return "! the profile is not finite on this grid"

    bml = collect(range(0.0, PROFILE_B_MAX; length = PROFILE_NB))
    V = try
        OITOOLS.hankel_vis(I, r, bml .* 1e6 .* 4.84813681109536e-9)
    catch err
        return "! the transform failed: " * _cause(err)
    end
    update_profile_plot!(sh.profileplot, r, I, bml, V)
    return ""
end

"""
    shell_flux_constraint() -> String

The constraint that makes the component fluxes behave as flux ratios, as `lhs\top\trhs\ttol`,
or `!` and a reason.

Built from the components actually present rather than typed, because the expression is one
term per component and gets long and easy to mistype the moment there are more than two -- and
a normalisation that silently omits a component is worse than none.

Only the components whose `f` is a plain number: one already derived from the others is the
better way to do this and needs no constraint on top, and adding one would pull against an
identity that already holds exactly.
"""
function shell_flux_constraint()
    m = _model()
    fkeys = sort!([k for k in keys(m.dict) if endswith(k, ",f") && m.dict[k] isa Real])
    length(fkeys) < 2 && return "! a flux constraint needs at least two components with " *
                                "numeric flux fractions"
    lhs = join(("\$" * k for k in fkeys), " + ")
    return join((lhs, "==", "1", "0.001"), '\t')
end

"""
    shell_default_bounds() -> String

Fill the lower and upper bounds of every free parameter from `default_bounds`, returning a
one-line summary or a message beginning with `!`.

The current dataset is passed when there is one, which is what makes this worth a button: the
angular ceiling then comes from the coverage itself -- twice the largest scale the shortest
baseline senses -- rather than from `DEFAULT_MAX_SIZE_MAS`. A source bigger than that is
resolved out and its size is not constrained by this data at all, so a bound drawn from the uv
coverage is a statement about the observation rather than a guess.

Only the free parameters: a bound on a fixed or derived one constrains nothing, and writing it
would put a number in the table that no fitter reads.
"""
function shell_default_bounds()
    sh = _shell(); m = _model(sh)
    isempty(m.free) && return "! no free parameters to bound"
    e = current_dataset(sh)
    data = e === nothing ? nothing : e.data
    lb, ub = try
        default_bounds(m.dict, m.free; data)
    catch err
        msg = "! default_bounds: " * _cause(err); console!(sh, msg; kind = :err); return msg
    end
    n = 0
    for k in m.free
        haskey(lb, k) || continue
        m.lb[k] = Float64(lb[k]); m.ub[k] = Float64(ub[k]); n += 1
    end
    src = data === nothing ? "no dataset loaded, so sizes use the default ceiling" :
                             "sizes scaled from the uv coverage"
    console!(sh, "> default_bounds(model, free" * (data === nothing ? "" : "; data") * ")"; kind = :cmd)
    console!(sh, "  bounded $(n) free parameter(s); " * src)
    return ""
end

"""
    shell_check_constraints(lines) -> String

Evaluate the constraints, as `1`/`0` per line in the order given, or `!` and a reason.

`lines` is the panel's own `lhs\top\trhs\ttol` text rather than anything stored here, so what
is checked is what is on screen -- including a constraint being edited and not yet used in a
fit.
"""
function shell_check_constraints(lines::AbstractString)
    m = _model()
    cons = parse_constraint_lines(lines)
    isempty(cons) && return "! no constraints to check"
    ok = try
        check_constraints(cons, m.dict; verb = false)
    catch err
        return "! " * _cause(err)
    end
    return join((b ? "1" : "0" for b in ok), "\n")
end

"""
Radial profiles worth starting from, as `key`, `label`, `expression`, `seeds` and a suggested
grid radius in mas.

A template carries three things because pasting the expression alone fixes only the typing.
Its parameters are seeded at values that mean something for the shape rather than at
`_component_seed`'s generic 2.0, and it names an `r_max` wide enough that the profile has
actually died away by the edge of the grid -- the double sigmoid at the paper's own numbers
sits at 54% of its peak on the grid a new component is created with, which truncates the ring
being modelled without saying so.

A shape with an exact analytic kind is deliberately NOT here. A uniform disc through a
numerical Hankel transform of a flat profile is `ud` computed worse -- measured at 8e-5 against
the closed form, plus a grid to get wrong -- and `\$MU^\$alpha` is `ldpow` the same way. A
template earns its place only where there is no analytic component to reach for instead.

`seeds` is `name=>value` in the order they should be created. `r_max` is a RADIUS; `udout` is
twice it.
"""
# Parameter names are lowercase throughout, and none of them is a function.
#
# Two rules, both learned the hard way. A name starting with a capital R reads as the implicit
# `$R` and was the case that broke `compile_profile` before it substituted whole tokens; and a
# name that is also a Julia function -- `sin` was the one here -- becomes an argument that
# SHADOWS it, so an expression using both `$sin` and `sin(...)` fails with "objects of type
# Float64 are not callable", which names nothing the user wrote. `rin`/`rout`/`win`/`wout`
# avoid both, and mean the same thing in every template.
const PROFILE_TEMPLATES = [
    (key = "doughnut", label = "Doughnut",
     expr = "1 - ((\$R - \$rbar)/(2*(\$rout - \$rin)))^2",
     seeds = ["rbar" => 1.5, "rin" => 1.0, "rout" => 2.0], rmax = 2.5,
     note = "parabolic ring: one width, symmetric rims"),

    (key = "double_sigmoid", label = "Double sigmoid",
     expr = "1/(1+exp(-(\$R-\$rin)/\$win)) * 1/(1+exp((\$R-\$rout)/\$wout))",
     seeds = ["rin" => 1.26, "rout" => 1.25, "win" => 0.26, "wout" => 0.43], rmax = 4.0,
     note = "inner and outer rims with independent sharpness; dark cavity"),

    (key = "alpha_double_sigmoid", label = "α double sigmoid",
     expr = "(\$a + (1-\$a)/(1+exp(-(\$R-\$rin)/\$win))) * 1/(1+exp((\$R-\$rout)/\$wout))",
     seeds = ["a" => 0.17, "rin" => 1.25, "rout" => 1.24, "win" => 0.10, "wout" => 0.50],
     rmax = 4.0,
     note = "as above, with α letting the cavity carry flux"),

    (key = "ring_gauss", label = "Ring ⊛ Gaussian",
     expr = "exp(-(\$R^2 + \$r0^2)/(2*\$s^2)) * besseli(0, \$R*\$r0/\$s^2)",
     seeds = ["r0" => 2.0, "s" => 0.3], rmax = 4.0,
     note = "the EXACT convolution of a thin ring with a Gaussian, not an approximation"),

    (key = "power_law", label = "Power law",
     expr = "(\$R/\$r0)^(-\$p)",
     seeds = ["r0" => 1.0, "p" => 1.5], rmax = 4.0,
     note = "falling envelope; r0 keeps the origin finite"),

]

"""
    shell_profile_templates() -> String

The templates, as `key\tlabel\tnote` lines, for the menus that offer them.
"""
shell_profile_templates() =
    join((string(t.key, '\t', t.label, '\t', t.note) for t in PROFILE_TEMPLATES), "\n")

"""
    shell_profile_template(key) -> String

One template as `expression\tr_max\tname=value,name=value…`, or `!` and a reason.

The caller puts the expression in the editor rather than committing it: choosing a shape is a
draft, and Apply/Revert is already there to arbitrate.
"""
function shell_profile_template(key::AbstractString)
    k = String(key)
    i = findfirst(t -> t.key == k, PROFILE_TEMPLATES)
    i === nothing && return "! no such profile template: " * k
    t = PROFILE_TEMPLATES[i]
    seeds = join((string(n, "=", v) for (n, v) in t.seeds), ",")
    return join((t.expr, string(t.rmax), seeds), '\t')
end

"""
    shell_component_geometry(name) -> String

Which optional geometry a component carries, as `position\torientation` (`1`/`0` each).

`f` is not offered: every component has one, and one without it is not a component the parser
would weight. `x`/`y` and `incl`/`pa` are the ones a component may or may not have, and their
absence is the reason an inclined ring could not be built here at all.
"""
function shell_component_geometry(name::AbstractString)
    m = _model(); cn = String(name)
    pos = haskey(m.dict, cn * ",x")    || haskey(m.dict, cn * ",y")
    ori = haskey(m.dict, cn * ",incl") || haskey(m.dict, cn * ",pa")
    return string(pos ? "1" : "0", '\t', ori ? "1" : "0")
end

"""
    shell_set_component_geometry(name, which, on) -> String

Add or remove a pair of optional geometry parameters. `which` is `position` or `orientation`.

Both members of a pair together, like the azimuthal modes: an inclination with no position
angle only means something by accident -- it inclines about whichever axis the sky happens to
give -- and a position angle with no inclination rotates a circle. Seeded at zero, which is the
value that leaves the component exactly as it was, so switching this on cannot move a fit.
"""
function shell_set_component_geometry(name::AbstractString, which::AbstractString, on)
    sh = _shell(); m = _model(sh)
    cn = String(name)
    any(startswith(k, cn * ",") for k in keys(m.dict)) ||
        return "! no component called " * cn
    keys_ = String(which) == "position"    ? ("x", "y") :
            String(which) == "orientation" ? ("incl", "pa") :
            return "! unknown geometry: " * String(which)
    if Bool(on)
        for suffix in keys_
            k = cn * "," * suffix
            haskey(m.dict, k) || (m.dict[k] = 0.0)
        end
        console!(sh, "  " * cn * ": added " * join((cn * "," * s for s in keys_), " and "))
    else
        for suffix in keys_
            k = cn * "," * suffix
            delete!(m.dict, k); delete!(m.lb, k); delete!(m.ub, k)
            filter!(!=(k), m.free)
        end
        console!(sh, "  " * cn * ": removed " * join((cn * "," * s for s in keys_), " and "))
    end
    return ""
end

"""
The limb-darkening laws, as one shape wearing different coefficients.

They are a sub-choice, not five separate components: every one of them is a disc of some
diameter, and which law you use is a question about the atmosphere rather than about the
geometry. Offering them as five entries beside "uniform ring" and "crescent" put two different
questions in one list and made limb darkening most of it.

`coeffs` is what the law needs BESIDES the diameter, in the parser's order.
"""
const LD_LAWS = [
    (key = "ldlin",     label = "linear",         coeffs = ["u"],
     note = "I(μ) = 1 − u(1−μ)"),
    (key = "ldquad",    label = "quadratic",      coeffs = ["u", "w"],
     note = "I(μ) = 1 − u(1−μ) − w(1−μ)²"),
    (key = "ldsqrt",    label = "square root",    coeffs = ["u", "w"],
     note = "I(μ) = 1 − u(1−μ) − w(1−√μ)"),
    (key = "ldpow",     label = "power law",      coeffs = ["alpha"],
     note = "I(μ) = μ^α  (Hestroffer)"),
    (key = "ldclaret4", label = "four-parameter", coeffs = ["c1", "c2", "c3", "c4"],
     note = "I(μ) = 1 − Σₖ cₖ(1 − μ^{k/2})  (Claret)"),
]

"""
    shell_ld_laws() -> String

The laws, as `key\tlabel\tnote` lines.
"""
shell_ld_laws() = join((string(l.key, '\t', l.label, '\t', l.note) for l in LD_LAWS), "\n")

"""
    shell_ld_law(component) -> String

The law a component currently uses, as `key\tlabel`, or `""` when it is not a limb-darkened
disc at all.
"""
function shell_ld_law(component::AbstractString)
    m = _model(); cn = String(component)
    i = findfirst(l -> haskey(m.dict, cn * "," * l.key), LD_LAWS)
    i === nothing && return ""
    return string(LD_LAWS[i].key, '\t', LD_LAWS[i].label)
end

"""
    shell_set_ld_law(component, law) -> String

Change which limb-darkening law a component uses, keeping its diameter.

The diameter is the same physical quantity in every law and carries over, along with whether it
was being fitted. The coefficients do NOT: `u` is the linear coefficient in the quadratic and
the square-root laws alike, but the two laws are different shapes and a value fitted under one
is not a value under the other. They are reseeded, and the console says so, rather than being
carried across to look like a measurement that survived the change.
"""
function shell_set_ld_law(component::AbstractString, law::AbstractString)
    sh = _shell(); m = _model(sh)
    cn = String(component); want = String(law)
    j = findfirst(l -> l.key == want, LD_LAWS)
    j === nothing && return "! unknown limb-darkening law: " * want
    i = findfirst(l -> haskey(m.dict, cn * "," * l.key), LD_LAWS)
    i === nothing && return "! " * cn * " is not a limb-darkened disc"
    old, new = LD_LAWS[i], LD_LAWS[j]
    old.key == want && return ""

    dkey = cn * "," * old.key
    diam = m.dict[dkey]
    wasfree = dkey in m.free
    lb = get(m.lb, dkey, nothing); ub = get(m.ub, dkey, nothing)

    for k in vcat([old.key], old.coeffs)
        kk = cn * "," * k
        delete!(m.dict, kk); delete!(m.lb, kk); delete!(m.ub, kk)
        filter!(!=(kk), m.free)
    end

    nkey = cn * "," * new.key
    m.dict[nkey] = diam
    wasfree && push!(m.free, nkey)
    lb === nothing || (m.lb[nkey] = lb)
    ub === nothing || (m.ub[nkey] = ub)
    for c in new.coeffs
        m.dict[cn * "," * c] = _component_seed(c)
    end
    console!(sh, "  " * cn * ": " * old.label * " → " * new.label *
                 "; diameter kept at " * string(diam) * ", coefficients reseeded")
    return ""
end

"""
    shell_az_modes(component) -> String

The azimuthal modes a component carries, as `order\tamp\tprojang` lines, lowest order first.

Read from the model rather than from anything the panel keeps, so the list cannot drift from
what the parser will build.
"""
function shell_az_modes(component::AbstractString)
    m = _model(); cn = String(component)
    out = Tuple{Int,Float64,Float64}[]
    for k in keys(m.dict)
        mm = match(Regex("^" * cn * ",az amp(\\d+)\$"), k)
        mm === nothing && continue
        n = parse(Int, mm.captures[1])
        amp = m.dict[k]
        phi = get(m.dict, cn * ",az projang" * string(n), 0.0)
        push!(out, (n, amp isa Real ? Float64(amp) : NaN,
                       phi isa Real ? Float64(phi) : NaN))
    end
    sort!(out; by = first)
    return join((join((string(n), string(a), string(p)), '\t') for (n, a, p) in out), "\n")
end

"""
    shell_add_az_mode(component) -> String

Add the next azimuthal mode, as the PAIR the parser requires.

`dict_to_model` errors when an `az amp<N>` has no matching `az projang<N>`, so the two are
written together and removed together -- exposing one without the other builds a model that
cannot be evaluated. Seeded at zero amplitude, which is the component exactly as it was: adding
the structure should not move a fit, only make the asymmetry available to one.
"""
function shell_add_az_mode(component::AbstractString)
    sh = _shell(); m = _model(sh); cn = String(component)
    any(startswith(k, cn * ",") for k in keys(m.dict)) ||
        return "! no component called " * cn
    haskey(m.dict, cn * ",profile") ||
        return "! azimuthal modes belong to a radial-profile component; " * cn * " is not one"
    n = 1
    while haskey(m.dict, cn * ",az amp" * string(n)); n += 1; end
    m.dict[cn * ",az amp" * string(n)]     = 0.0
    m.dict[cn * ",az projang" * string(n)] = 0.0
    console!(sh, "  " * cn * ": added az mode " * string(n) * " (amp and projang, both 0)")
    return ""
end

"""
    shell_remove_az_mode(component, order) -> String

Remove one mode, both keys together, for the reason they are added together.
"""
function shell_remove_az_mode(component::AbstractString, order)
    sh = _shell(); m = _model(sh); cn = String(component)
    n = order isa Integer ? Int(order) : Int(round(Float64(order)))
    ka, kp = cn * ",az amp" * string(n), cn * ",az projang" * string(n)
    haskey(m.dict, ka) || return "! " * cn * " has no az mode " * string(n)
    for k in (ka, kp)
        delete!(m.dict, k); delete!(m.lb, k); delete!(m.ub, k)
        filter!(!=(k), m.free)
    end
    console!(sh, "  " * cn * ": removed az mode " * string(n))
    return ""
end

"""
    shell_remove_component(name) -> String

Delete a component and every key belonging to it, including its bounds and its place in the
free list.

Expressions elsewhere that referenced it are deliberately NOT rewritten: a dangling `\$name,ud`
shows up in the inspector's unrecognised list, which is a visible problem, where a silent
substitution would be an invisible one.
"""
function shell_remove_component(name)
    sh = _shell(); m = _model(sh)
    nm = String(name); pre = nm * ","
    ks = [k for k in keys(m.dict) if startswith(k, pre)]
    isempty(ks) && return "! no component called " * nm
    for k in ks
        delete!(m.dict, k); delete!(m.lb, k); delete!(m.ub, k)
    end
    filter!(p -> !startswith(p, pre), m.free)
    console!(sh, "  removed $(nm): $(length(ks)) keys")
    return ""
end

"""
    shell_rename_component(old, new) -> String

Rename a component, carrying its parameters, its bounds and its place in the free list with it.

Unlike removal, expressions ARE rewritten. `\$old,ud` and `\$new,ud` name the same parameter of
the same component, so leaving the old spelling behind would break a model the user never
edited -- where after a removal a dangling reference is the honest report that something the
model depended on is gone. The reference form is `\$name,suffix` and a name cannot contain a
comma, so replacing `\$old,` matches the whole name and never a longer one that merely starts
with it.
"""
function shell_rename_component(old, new)
    sh = _shell(); m = _model(sh)
    o = String(old); nm = strip(String(new))
    isempty(nm) && return "! a component needs a name"
    occursin(',', nm) && return "! a component name cannot contain a comma"
    nm == GLOBAL_COMPONENT && return "! that name is reserved for globals"
    nm == o && return ""
    opre = o * ","; npre = nm * ","
    ks = [k for k in keys(m.dict) if startswith(k, opre)]
    isempty(ks) && return "! no component called " * o
    any(startswith(k, npre) for k in keys(m.dict)) &&
        return "! there is already a component called " * nm

    for k in ks
        nk = npre * k[nextind(k, lastindex(opre)):end]
        m.dict[nk] = m.dict[k];                      delete!(m.dict, k)
        haskey(m.lb, k) && (m.lb[nk] = m.lb[k];      delete!(m.lb, k))
        haskey(m.ub, k) && (m.ub[nk] = m.ub[k];      delete!(m.ub, k))
        replace!(m.free, k => nk)
    end
    nref = 0
    for k in collect(keys(m.dict))
        v = m.dict[k]
        v isa AbstractString && occursin("\$" * opre, v) || continue
        m.dict[k] = replace(v, "\$" * opre => "\$" * npre)
        nref += 1
    end
    console!(sh, "  renamed $(o) → $(nm): $(length(ks)) keys" *
                 (nref > 0 ? ", $(nref) expression(s) rewritten" : ""))
    return ""
end

"""
    _parse_options(text) -> Dict{String,String}

The engine settings panel's `key\tvalue` lines.

Values stay strings here and are converted by whichever engine branch reads them: the panel
sends what its controls hold, and only the engine knows whether a given key is an integer, a
weight or a file path.
"""
function _parse_options(text::AbstractString)
    o = Dict{String,String}()
    for line in split(String(text), '\n')
        isempty(line) && continue
        f = split(line, '\t')
        length(f) == 2 || continue
        o[String(f[1])] = String(f[2])
    end
    return o
end

"""
    shell_job_poll() -> String

`state\toutput` — what the running job has printed, and whether it is still going.

`state` is `running`, `done`, or `idle`. On `done` the job has already been wound up and its
result applied: the canvas is drawn, the fit recorded, the console written. QML calls this from
a Timer, which is what makes the output appear while it is produced rather than at the end.

The capture file is written unbuffered, so `output` is everything printed up to this instant
rather than up to the last buffer flush.
"""
function shell_job_poll()
    sh = SHELL[]
    (sh === nothing || sh.job === nothing) && return "idle\t"
    j = sh.job
    istaskdone(j.task) || return "running\t" * job_output(sh)

    stopped = j.stop[]
    kind = j.kind
    res  = job_finish!(sh)
    if stopped
        # Abandoned rather than aborted; see `shell_job_stop`. Whatever came back is discarded,
        # because the user asked to stop looking at it.
        line = "stopped"
        kind === :image ? (sh.enginelog = res.output) : (sh.fitlog = res.output)
        console!(sh, "  " * line)
    else
        line = kind === :image ? finish_reconstruct!(sh, res) : finish_fit!(sh, res)
    end
    return "done\t" * line
end

"""
    shell_job_stop() -> String

Ask the running job to stop.

None of the engines takes a cancellation token -- there is no callback API anywhere in the
reconstruction or fitting code -- so this cannot abort the computation. What it does is
ABANDON it: the flag is set, an interrupt is thrown at the task in case it is at a yield point,
the GUI returns to idle immediately, and whatever the worker eventually produces is discarded.

The CPU may stay busy until the current run ends. Saying so is better than a Stop button that
appears to work and silently leaves a thread running, and better than no button at all.
"""
function shell_job_stop()
    sh = SHELL[]
    (sh === nothing || sh.job === nothing) && return "nothing is running"
    sh.job.stop[] = true
    # Only lands if the task is at a yield point; a blocking ccall will not see it. Harmless
    # when it does not.
    try
        Base.schedule(sh.job.task, InterruptException(); error = true)
    catch
    end
    console!(sh, "> stop requested")
    return "stopping…"
end

"True while a job is running, so the panel can disable what must not be touched."
shell_job_running() = (sh = SHELL[]; sh !== nothing && sh.job !== nothing)

"""
    GuiJob

A long-running Julia call, off the GUI thread, with its output readable while it runs.

Reconstructions and fits used to run ON the GUI thread. The window froze for the duration, the
engine's output could only appear once it was over, and Stop could not do anything because
nothing was listening. Moving the call to a worker fixes all three at once.

The split matters: the worker computes and NOTHING else. Every canvas update stays on the GUI
thread, in `job_finish!`, because GL calls from a worker are the crash §5.4b of the design
documents. `stop` is a flag the GUI sets and the finish path reads.
"""
mutable struct GuiJob
    task      :: Task
    kind      :: Symbol           # :image or :fit — decides what `job_finish!` does with it
    path      :: String           # capture file, read while it grows
    io        :: Base.Filesystem.File
    saved     :: Any              # the `Base.stdout` to put back
    savedfd   :: Cint             # ...and the descriptor
    stop      :: Base.RefValue{Bool}
    started   :: Float64
end

"""
    start_job!(sh, kind, f) -> String

Redirect output to a file, run `f` on a worker thread, and remember the job.

The redirect is set up HERE, on the GUI thread, and taken down in `job_finish!`, also on the
GUI thread. `Base.stdout` is a global: assigning it from the worker would race with whatever
the GUI thread is doing.
"""
function start_job!(sh::ShellState, kind::Symbol, f)
    path, tmp = mktemp(); close(tmp)
    # A `Filesystem.File`, NOT the `IOStream` mktemp hands back. An IOStream buffers, and the
    # whole point here is that the console fills while the run is going: measured, a 54 kB
    # VMLMB trace arrived in a single lump at the end through an IOStream, and line by line
    # through this. Each write is a syscall, which is affordable for a progress log.
    io = Base.Filesystem.open(path, Base.Filesystem.JL_O_WRONLY |
                                    Base.Filesystem.JL_O_CREAT |
                                    Base.Filesystem.JL_O_TRUNC, 0o644)
    saved   = Base.stdout
    savedfd = Sys.isunix() ? ccall(:dup, Cint, (Cint,), 1) : Cint(-1)
    flush(stdout)
    Base.stdout = IOContext(io, :color => true)
    savedfd >= 0 && ccall(:dup2, Cint, (Cint, Cint),
                          Base.cconvert(Cint, Base.fd(io)), Cint(1))
    stop = Ref(false)
    task = Threads.@spawn f(stop)
    sh.job = GuiJob(task, kind, path, io, saved, savedfd, stop, time())
    return ""
end

"Whatever the running job has printed so far. Buffered, so it arrives in chunks."
function job_output(sh::ShellState)
    j = sh.job
    j === nothing && return ""
    return try; read(j.path, String); catch; ""; end
end

"""
    job_finish!(sh) -> (; ok, result, err, output)

Put stdout back, collect what the worker produced, and delete the job. **GUI thread only.**
"""
function job_finish!(sh::ShellState)
    j = sh.job
    j === nothing && return (; ok = false, result = nothing, err = nothing, output = "")
    Base.stdout = j.saved
    if j.savedfd >= 0
        ccall(:dup2, Cint, (Cint, Cint), j.savedfd, Cint(1))
        ccall(:close, Cint, (Cint,), j.savedfd)
    end
    try; close(j.io); catch; end
    out = try; read(j.path, String); catch; ""; end
    rm(j.path; force = true)
    res, err = nothing, nothing
    try
        res = fetch(j.task)
    catch e
        err = e isa TaskFailedException ? e.task.exception : e
    end
    sh.job = nothing
    return (; ok = err === nothing, result = res, err, output = out)
end

"""
    _capture_stdout(f) -> (result, output)

Run `f` with `stdout` on a temp file, and return everything it printed.

The engines have no callback API -- they report by printing -- so this is what puts their own
account in front of the user. Getting it took four attempts, and the three that failed are
properties of this environment rather than carelessness:

  * `redirect_stdout()` onto a PIPE hangs. Not the cooperative-scheduling deadlock it looks
    like: it hangs with the reader on a spawned thread too. Inside `QML.exec()` Qt owns the
    thread and libuv's event loop is never pumped, so a libuv pipe read cannot progress at all.
    Measured both ways in a standalone QML app; both variants work in a plain Julia process,
    which is what makes it Qt's doing.
  * `redirect_stdout(io)` onto a file throws `ArgumentError: invalid stdio type` here, because
    QML.jl has already put the Core streams in place.
  * `dup2` on descriptor 1 alone captures nothing from Julia: `println` goes through the libuv
    handle bound at startup, not through the descriptor.

Assigning the binding is what QML.jl itself does (`Base.stdout = Core.stdout`), and a plain
`IOStream` needs no event loop. `dup2` is layered on top so that C-level printing -- OptimPack
and the MaxEnt kernel write with `printf` -- lands in the same file rather than escaping to the
terminal.

The captured text keeps its ANSI escapes: the console renders them, so the colours OITOOLS
chose survive the trip instead of being flattened to black.
"""
function _capture_stdout(f)
    path, io = mktemp()
    old = Base.stdout
    saved = Sys.isunix() ? ccall(:dup, Cint, (Cint,), 1) : Cint(-1)
    result = nothing
    try
        # `IOContext(:color => true)`, not the bare stream. `printstyled` emits ANSI escapes
        # only when `get(stdout, :color, false)` is true, and an `IOStream` says false -- so
        # capturing to a plain file silently threw away every colour OITOOLS chose. chi2_flat
        # prints V2 red, T3amp blue and T3phi green, and that is worth keeping.
        Base.stdout = IOContext(io, :color => true)
        # Same file for the C side. `Base.fd` yields a `RawFD`, which has no `Int32`
        # constructor; `cconvert` is the way down to the Cint the syscall wants.
        saved >= 0 && ccall(:dup2, Cint, (Cint, Cint),
                            Base.cconvert(Cint, Base.fd(io)), Cint(1))
        result = f()
    finally
        Base.stdout = old
        if saved >= 0
            ccall(:dup2, Cint, (Cint, Cint), saved, Cint(1))
            ccall(:close, Cint, (Cint,), saved)
        end
        try; close(io); catch; end
    end
    out = try; read(path, String); catch; ""; end
    rm(path; force = true)
    return (result, out)
end


"""
    shell_reset_image_zoom() -> String

Return the imaging display to the whole image.

Separate from [`shell_reset_zoom`](@ref), which resets the Exploring canvas: the two are
different axes on different figures, and a right click on one must not move the other.
"""
function shell_reset_image_zoom()
    sh = _shell()
    cv = sh.imcanvas
    cv === nothing && return ""
    try
        Makie.autolimits!(cv.axis)
    catch err
        console!(sh, "could not reset the view: " * _cause(err); kind = :err)
    end
    return ""
end

"""
    shell_show_start_image(nx, pixsize, mode, startkind, startfwhm, startseed) -> String

Draw the starting image on the imaging canvas, without running anything.

The start is not a formality: a Dirac, a Gaussian and a random field lead a reconstruction to
different places, and the width of a Gaussian start is a real prior. Being able to look at it
before spending a run is the point.
"""
function shell_show_start_image(nx::Integer, pixsize::Real, mode::AbstractString,
                                startkind::AbstractString, startfwhm::Real = AUTO_FWHM,
                                startseed::Integer = 1, startpath::AbstractString = "")
    sh = _shell()
    e = current_dataset(sh)
    e === nothing && return "! no dataset loaded"
    setup = ImagingSetup(; nx = Int(nx), pixsize = Float64(pixsize), mode = Symbol(mode),
                           startkind = Symbol(startkind), startfwhm = Float64(startfwhm),
                           startseed = Int(startseed), startpath = String(startpath))
    img = try
        ft, _, _ = ensure_ft!(sh.ftcache, e.data, setup)
        Float64.(start_image(setup, ft))
    catch err
        msg = "! could not build the starting image: " * _cause(err)
        console!(sh, msg; kind = :err); return msg
    end
    sh.imcanvas === nothing || show_image!(sh.imcanvas, img, setup.pixsize;
                                           label = "starting image")
    console!(sh, "> start_image(setup, ft)   # $(startkind), $(Int(nx))×$(Int(nx))")
    return "showing the starting image"
end

"""
    shell_save_image(path) -> String

Write the reconstructed image to FITS, with the pixel size in the header.

`pixsize` is passed so the file carries its own scale: an image without one is a picture, not
a measurement, and nothing downstream can put it back.
"""
function shell_save_image(path::AbstractString)
    sh = _shell()
    r = sh.imaging
    r === nothing && return "! nothing reconstructed yet"
    p = String(path); startswith(p, "file://") && (p = p[8:end])
    try
        writefits(r.image, p; pixsize = r.setup.pixsize)
    catch err
        msg = "! could not write the image: " * _cause(err)
        console!(sh, msg; kind = :err); return msg
    end
    console!(sh, "> writefits(image, \"" * p * "\"; pixsize = $(r.setup.pixsize))")
    return "saved " * basename(p)
end

"""
    shell_engine_output() -> String

Everything the last reconstruction printed.

The engines report their progress by printing -- there is no callback API -- so this is
captured from stdout for the duration of the run and handed back whole. It is the engine's own
account of what it did: MaxEnt's α and entropy per iteration, ADMM's primal and dual residuals,
the sampler's temperature and acceptance.
"""
shell_engine_output() = (sh = SHELL[]; sh === nothing ? "" : sh.enginelog)

"""
    shell_recenter_image() -> String

Recentre the reconstructed image on its centroid and redraw.

Nothing in the reconstruction holds the image still: V² and closure phase are both invariant
under translation, so the image is free to drift anywhere in the field and routinely does.
`recenter` circularly shifts the centroid back to the middle. That leaves the fit untouched for
the same invariance -- as long as no flux wraps around the edge, which is what a centred image
is for in the first place.
"""
function shell_recenter_image()
    sh = SHELL[]
    sh === nothing && return ""
    r = sh.imaging
    r === nothing && return "nothing reconstructed yet"
    img = recenter(r.image)
    # ImagingResult is immutable, so rebuild it. Every other field survives the shift: the
    # χ², the point count and the timings all describe the run that produced the pixels, and
    # moving them does not change what that run did.
    sh.imaging = ImagingResult(img, r.chi2, r.chi2_start, r.ndof, r.flux, r.maxiter,
                               r.seconds, r.setup, r.weights, r.breakdown)
    sh.imcanvas === nothing || show_image!(sh.imcanvas, img, r.setup.pixsize)
    console!(sh, "> image = recenter(image)")
    return "recentred on the centroid"
end

"""
    shell_image_colormaps() -> String
    shell_image_colormap(name) -> String

The colormaps offered for the reconstructed image, and the setter for one of them.

Applies to the imaging canvas — the one the reconstruction is drawn on — not to the explore
view, whose colours mean baselines rather than flux.
"""
shell_image_colormaps() = image_colormap_names()

function shell_image_colormap(name)
    sh = SHELL[]
    sh === nothing && return ""
    cv = sh.imcanvas
    cv === nothing && return "no image canvas"
    if set_colormap!(cv, String(name))
        console!(sh, "colormap: " * String(name))
        return String(name)
    end
    return "unknown colormap"
end

"""
    shell_ui_scale(x) -> String

Record the UI scale QML worked out from the screen, so anything drawn in Julia matches it.

Called once from `Component.onCompleted`, which runs before the first plot is drawn — the
canvas is built earlier still, but every redraw restyles from the live value.
"""
function shell_ui_scale(x, dpi = 0.0)
    before = live_plot_scale()
    v = set_ui_scale!(Float64(x), Float64(dpi))
    sh = SHELL[]
    sh === nothing && return ""
    console!(sh, @sprintf("plot scale: %.2f (ui %.2f, %.0f dpi)",
                         live_plot_scale(), v, SCREEN_DPI_LIVE[]))

    # Restyle. `gui()` draws the first plot before `exec()`, so by the time QML reports the
    # scale the figure already exists at whatever the default was -- and `makieArea.update()`
    # only re-renders it, it does not re-run `update_canvas!`. Without this the very first
    # view keeps the wrong font and marker sizes until something else forces a redraw.
    #
    # Safe here despite the window being up: `update_canvas!` only assigns to Observables and
    # creates no plots, which is what the whole pre-built canvas design is for.
    live_plot_scale() ≈ before && return ""
    current_dataset(sh) === nothing && return ""
    try
        sh.status = refresh_plot!(sh)
    catch err
        console!(sh, "could not restyle for the new scale: " * _cause(err); kind = :err)
    end
    return ""
end

"Facilities on disk, newline-separated. CHARA first: it is the only one the delay-line and POP checks cover."
function shell_facilities()
    c = config_catalog()
    fac = sort(String.(c.facilities); by = f -> (f != "CHARA", f))
    return join(fac, "\n")
end

"Telescope names for a facility, in the order its configuration lists them."
shell_telescopes(name::AbstractString) = join(facility_telescopes(String(name)), "\n")

"""
    shell_gantt(facility, name, ra, dec, dateiso, telescopes, pops, use_delay, detailed) -> String

Compute one night and draw it on the Observe canvas. Returns a one-line summary.

`ra` and `dec` are DEGREES, `telescopes` and `pops` are space-separated. `use_delay` is opt-in
because an unsearched POP configuration reports far less time than is really available — see
[`night_plan`](@ref).
"""
function shell_gantt(facility::AbstractString, name::AbstractString, ra::Real, dec::Real,
                     dateiso::AbstractString, telescopes::AbstractString,
                     pops::AbstractString, use_delay::Bool, detailed::Bool,
                     alt_limit::Real = DEFAULT_ALT_LIMIT, alt_max::Real = DEFAULT_ALT_MAX)
    sh = _shell()
    date = try
        DateTime(String(dateiso))
    catch
        return "! not a date: $(dateiso)"
    end
    tels = facility_telescopes(String(facility))
    sel  = filter(!isempty, String.(split(String(telescopes))))
    cfg  = isempty(sel) ? ones(Int, length(tels)) : telescope_config(sel, tels)
    pp   = isempty(strip(String(pops))) ? nothing :
           [parse(Int, x) for x in split(String(pops))]

    p = try
        night_plan(String(facility), String(name), Float64(ra), Float64(dec), date;
                   config = cfg, pop = pp, use_delay = use_delay,
                   alt_limit = Float64(alt_limit), alt_max = Float64(alt_max))
    catch err
        msg = "! " * _cause(err); console!(sh, msg); return msg
    end

    sh.gantt === nothing || update_gantt!(sh.gantt, p; detailed)
    # The delay chart is drawn from the same plan, so both views are always of one night.
    sh.delayplot === nothing || isempty(p.baselines) || update_delay_plot!(sh.delayplot, p)

    sh.plan = p          # the hover readout reads this; see shell_gantt_hover
    h = observable_hours(p)
    console!(sh, "> night_plan(\"$(facility)\", \"$(name)\", $(round(Float64(ra); digits=4)), " *
                 "$(round(Float64(dec); digits=4)), DateTime(\"$(dateiso)\"))")

    # Zero hours with the delay check on and no POPs searched is almost always the POPs, not
    # the sky: an arbitrary configuration can close the night outright, and an empty chart with
    # no explanation reads as a broken tool rather than an unsearched one.
    why = if p.delay_applied && h == 0 && isempty(strip(String(pops)))
        "  — POPs not searched; try Search POPs, or fewer telescopes"
    elseif !p.delay_applied
        "  (delay not checked)"
    else
        ""
    end
    line = @sprintf("%s: %.2f h observable%s   moon %d%% lit",
                    name, h, why, round(Int, 100 * p.moon_fli))
    console!(sh, "  " * line)

    # `summary \t dark window`. The summary is for the console, which keeps a transcript; the
    # dark window is a readout the panel shows on its own. Returning both from the one call
    # avoids a second pass over the night to recompute what this one already knows.
    dark = isempty(p.good_twilight) ? "no astronomical night" :
           string(_hhmm(p.lst[first(p.good_twilight)]), " – ",
                  _hhmm(p.lst[last(p.good_twilight)]), " LST")
    return line * "\t" * dark
end

"The best POP configurations for this array and declination, as `pop\tscore` rows."
function shell_best_pops(facility::AbstractString, dec::Real, dateiso::AbstractString,
                         ra::Real, telescopes::AbstractString, n::Integer = 5)
    sh = SHELL[]
    say(t) = sh === nothing ? nothing : console!(sh, t)
    date = try
        Dates.DateTime(String(dateiso))
    catch
        say("! POP search: '$(dateiso)' is not a date"); return ""
    end
    try
        f = read_facility_file(String(facility))
        obs = night_observability(f, Float64(ra), Float64(dec), date)
        tels = facility_telescopes(String(facility))
        sel  = filter(!isempty, String.(split(String(telescopes))))
        cfg  = isempty(sel) ? ones(Int, length(tels)) : telescope_config(sel, tels)
        say("> best_pop(facility, $(round(Float64(dec); digits=4)), ha, $(cfg))")
        res = best_pops(f, Float64(dec), obs.ha, cfg; n = Int(n))
        if isempty(res)
            say("  no POP configuration reaches this target with these telescopes")
            return ""
        end
        say("  best: POPs $(join(res[1].pop, " ")) — $(res[1].score) steps")

        # `pops \t score \t haRange` — the table promises an HA range, and it is worth the
        # extra in_delay per candidate (microseconds): a configuration is chosen for WHEN it
        # works, not only for how long.
        rows = String[]
        for r in res
            span = try
                g = in_delay(f, Float64(dec), obs.ha, cfg, Int.(r.pop)).good_delay
                isempty(g) ? "—" :
                    string(round(obs.ha[first(g)]; digits = 2), " … ",
                           round(obs.ha[last(g)];  digits = 2), " h")
            catch
                "—"
            end
            push!(rows, join((join(string.(r.pop), " "), string(r.score), span), "\t"))
        end
        return join(rows, "\n")
    catch err
        say("! POP search failed: " * _cause(err))
        return ""
    end
end

"""
    shell_ft_setup(nx, pixsize, mode) -> String

Build the Fourier plan for this geometry and describe it, reusing the last one when nothing
changed.

Called whenever nx, the pixel size or the transform changes, so the readout beside the
selector always describes the plan a Run would actually use. Building an NFFT plan is the
expensive part of setting up a reconstruction, which is why it is keyed and not repeated.
"""
function shell_ft_setup(nx::Integer, pixsize::Real, mode::AbstractString)
    sh = _shell()
    e = current_dataset(sh)
    e === nothing && (sh.ftcache = nothing; return "no dataset")
    setup = ImagingSetup(; nx = Int(nx), pixsize = Float64(pixsize), mode = Symbol(mode))
    try
        ft, summary, rebuilt = ensure_ft!(sh.ftcache, e.data, setup)
        sh.ftcache = (; key = (objectid(e.data), setup.nx, setup.pixsize, setup.mode),
                        ft, summary)
        rebuilt && console!(sh, "> setup_ft(data, $(setup.nx), $(setup.pixsize); " *
                                "mode = \"$(mode)\")")
        return summary
    catch err
        msg = "! plan failed: " * _cause(err)
        console!(sh, msg)
        sh.ftcache = nothing
        return msg
    end
end

"""
    shell_reconstruct(nx, pixsize, mode, startkind, v2, t3amp, t3phi, maxiter) -> String

Run one image reconstruction on the current dataset and put the result on the canvas.

Positivity is the only thing shaping the image: it is VMLMB's `lower = 0` rather than a
regulariser, so it is in force with the regulariser list empty. That is a real reconstruction —
an unregularised maximum-likelihood image — and it is the configuration in which a broken plan,
criterion or weighting shows up as a chi2 that does not move.

Runs on the calling thread, which is Qt's. A reconstruction takes seconds and the window will
not repaint during it; the caller is expected to have said so before calling.
"""
function shell_reconstruct(nx::Integer, pixsize::Real, mode::AbstractString,
                           startkind::AbstractString,
                           v2::Bool, t3amp::Bool, t3phi::Bool, maxiter::Integer,
                           regularizers::AbstractString = "",
                           startfwhm::Real = AUTO_FWHM, startseed::Integer = 1,
                           from_previous::Bool = false,
                           engine::AbstractString = "vmlmb",
                           options::AbstractString = "")
    sh = _shell()
    e = current_dataset(sh)
    e === nothing && return "no dataset loaded"

    weights = imaging_weights(; v2, t3amp, t3phi)
    all(iszero, weights) && return "no observables selected"

    eng = Symbol(engine)
    haskey(IMAGING_ENGINES, eng) ||
        return "! $(engine) is not wired; pick one of " *
               join(sort(string.(keys(IMAGING_ENGINES))), ", ")
    opts = _parse_options(options)

    setup = ImagingSetup(; engine = eng,
                           nx = Int(nx), pixsize = Float64(pixsize),
                           mode = Symbol(mode), startkind = Symbol(startkind),
                           startfwhm = Float64(startfwhm), startseed = Int(startseed))
    regs = try
        parse_regularizers(regularizers)
    catch err
        msg = "! " * _cause(err); console!(sh, msg); return msg
    end
    regtext = isempty(regs) ? "positivity only" :
              "regularizers = [" * join((string(r) for r in regs), ", ") * "]"
    startdesc = from_previous ? "x_start = previous image" : "x_start = $(startkind)"
    # The call that was actually made, named by the engine that made it. A console that said
    # `reconstruct` for every engine would be a log of something that did not happen.
    console!(sh, "> $(IMAGING_ENGINES[eng])(x_start, data, ft; maxiter = $(Int(maxiter)))  " *
                 "# $(Int(nx))×$(Int(nx)) at $(round(Float64(pixsize); digits=4)) mas, " *
                 "$startdesc, " * regtext)
    isempty(opts) || console!(sh, "  " * join(("$k = $v" for (k, v) in sort(collect(opts), by = first)), ", "))
    # Continue from the last image rather than from a fresh start.
    xprev = nothing
    if from_previous
        sh.imaging === nothing && return "no previous image to continue from"
        xprev = sh.imaging.image
        size(xprev, 1) == Int(nx) ||
            return "the previous image is $(size(xprev,1))×$(size(xprev,2)), not $(Int(nx))×$(Int(nx))"
    end

    cached = sh.ftcache
    ftkey  = (objectid(e.data), setup.nx, setup.pixsize, setup.mode)
    ft     = (cached !== nothing && cached.key == ftkey) ? cached.ft : nothing
    # `verb = true`: the console below the panel shows what the engine printed, so it has to
    # print. The run itself goes to a worker thread — see `GuiJob` — so the window stays alive,
    # the output can be read while it is produced, and Stop has something to act on.
    sh.enginelog = ""
    start_job!(sh, :image, function (_stop)
        reconstruct_image(e.data, setup; weights, maxiter = Int(maxiter),
                          regularizers = regs, ft, x_start = xprev, options = opts,
                          verb = true)
    end)
    return "reconstructing…"
end

"""
    finish_reconstruct!(sh, res) -> String

Take a finished reconstruction and put it on screen. **GUI thread only** — `show_image!` is a
GL call, and §5.4b of the design has the measurement for what happens when those run on a
worker.
"""
function finish_reconstruct!(sh::ShellState, res)
    r = res.result
    sh.enginelog = res.output
    if !res.ok
        msg = "! reconstruction failed: " * _cause(res.err)
        sh.enginelog = isempty(res.output) ? msg : res.output * "\n" * msg
        console!(sh, msg)
        return msg
    end
    sh.imcanvas === nothing || show_image!(sh.imcanvas, r.image, r.setup.pixsize)
    sh.imaging = r

    # Total flux is reported because nothing constrains it: V² and closure phase are both
    # invariant under a global scaling, so an image fitted from them alone drifts away from
    # unit flux, and a user who does not know that reads the number as a bug.
    line = @sprintf("chi2r %.4g -> %.4g   flux %.4f   %d iter, %.1f s",
                    chi2r_start(r), chi2r(r), r.flux, r.maxiter, r.seconds)
    console!(sh, "  " * line)
    isempty(r.breakdown) ||
        console!(sh, "  " * join((@sprintf("%s %.3f%s", b.name, b.chi2r, b.used ? "" : " (unused)")
                                  for b in r.breakdown), "   "))
    return line
end

"""
    shell_chi2_breakdown() -> String

The last run's reduced chi2 per observable, as `name\tchi2r\tn\tused` rows.

Separate from the run so the panel can re-read it without re-running anything.
"""
function shell_chi2_breakdown()
    r = _shell().imaging
    r === nothing && return ""
    return join((join((b.name, string(round(b.chi2r; digits = 4)), string(b.n),
                       b.used ? "1" : "0"), "\t") for b in r.breakdown), "\n")
end

"""
    shell_image_defaults() -> String

Image geometry the current dataset suggests, as `"pixsize=…,nx=…"`.

`auto_pixsize` samples the longest baseline about three times over. A constant in its place is
not a neutral default: too coarse and the model cannot represent what the data resolves, too
fine and the reconstruction is unconstrained at the pixel scale.
"""
function shell_image_defaults()
    e = current_dataset()
    e === nothing && return ""
    s = imaging_defaults(e.data)
    return "pixsize=$(round(s.pixsize; digits = 4)),nx=$(s.nx)"
end

"The last reconstruction, or an empty string when there is none."
shell_imaging_summary() =
    _shell().imaging === nothing ? "" : sprint(show, _shell().imaging)

"""
    shell_gantt_hover() -> String

What the Gantt says at the cursor, as lines, or `""` when the pointer is off a row.

Read from Makie's event stream rather than from a QML mouse handler, the same route picking
takes: QMLMakie forwards mouse moves into the scene, so the position is already in the axis's
own coordinates and no conversion between the two coordinate systems is needed.

Reports what was COMPUTED, not what was typed. A row shows the elevation window in force and
the night's totals; the cursor's x gives the instant — LST, hour angle, altitude, azimuth, and
whether the delay lines reach at that moment, which is the question a Gantt is read to answer.
"""
function shell_gantt_hover()
    sh = SHELL[]
    (sh === nothing || sh.gantt === nothing || sh.plan === nothing) && return ""
    p = sh.plan
    ax = sh.gantt.axis
    mp = try
        Makie.mouseposition(ax.scene)
    catch
        return ""
    end
    x, y = Float64(mp[1]), Float64(mp[2])
    (isfinite(x) && isfinite(y)) || return ""

    # Summary geometry: the rows are what the hover needs to name, and in the detailed view
    # the extra rows are baselines whose own labels are already on the y axis.
    geo = gantt_geometry(p)
    # Which row: the nearest, and only if the pointer is actually within half a row of it.
    isempty(geo.rows) && return ""
    _, k = findmin(r -> abs(r[1] - y), geo.rows)
    row = geo.rows[k]
    abs(row[1] - y) > 0.6 && return ""

    # Nearest sample in time. The grid is regular, so this is the instant under the cursor.
    isempty(p.lst) && return ""
    _, i = findmin(t -> abs(t - x), p.lst)

    out = String[]
    push!(out, row[2])
    push!(out, @sprintf("LST %s   HA %+.2f h", _hhmm(p.lst[i]), p.ha[i]))
    push!(out, @sprintf("alt %.1f°   az %.1f°", p.alt[i], p.az[i]))
    push!(out, @sprintf("elevation window %.0f–%.0f°", p.alt_limit, p.alt_max))
    inside(v, idx) = i in idx ? v : "no " * v
    push!(out, join((inside("dark", p.good_twilight),
                     inside("above limits", p.good_alt),
                     inside("moon ok", p.good_moon),
                     p.delay_applied ? inside("in delay", p.good_delay) : "delay not checked"),
                    " · "))
    push!(out, @sprintf("night: %.2f h observable   moon %d%% lit",
                        observable_hours(p), round(Int, 100 * p.moon_fli)))
    isempty(p.pop) || push!(out, "POPs " * pop_label(p))
    return join(out, "\n")
end

"""
    shell_grouping_noun(kind) -> String

What one panel would contain for this observable: `baseline`, `triplet` or `station`. The tick
box is labelled from it, so it never reads "per baseline" over a grid of closure triangles.
"""
shell_grouping_noun(kind::AbstractString) =
    haskey(OBS_SPECS, Symbol(kind)) ? grouping_noun(Symbol(kind)) : ""

"""
    shell_plot_kinds() -> String

`kind=1` or `kind=0` for every entry of the Explore plot menu, comma-separated.

Computed HERE rather than in QML so the menu cannot disagree with `canvas_data`: the two rules
that are easy to get wrong both live in the same file as the plotting. A kind greyed out for a
file that would in fact plot is as unhelpful as one offered that throws.

A "differential" plot is one against WAVELENGTH, so it needs more than one wavelength in the
data — which is a property of the file, not of how it was binned on load.

`visphi` and `diffphi` read the SAME field. Which of the two is meaningful is decided by
`PHITYP`: a differential phase is measured against the other channels, so it belongs against
wavelength and `canvas_data` refuses it against baseline. Exactly one of the pair is available
for any given file.
"""
function shell_plot_kinds()
    e = current_dataset()
    e === nothing && return ""
    d = e.data[1, 1]
    a = try
        observable_availability(e.data)
    catch
        return ""
    end
    diff = occursin("differential", lowercase(String(d.phityp)))
    # Distinct wavelengths in the data, NOT `size(e.data, 1)` -- that is the spectral BINNING
    # dimension, and `readoifits` leaves it at 1 unless asked to bin, so a 119-channel file
    # reported as monochromatic and `diffvisamp` was greyed out on exactly the files it is for.
    poly = length(unique(d.uv_lam)) > 1
    flags = ("uv"         => true,                      # geometry, present whenever data is
             "v2"         => a.v2,
             "t3phi"      => a.t3phi,
             "t3amp"      => a.t3amp,
             "visamp"     => d.nvisamp > 0,
             "visphi"     => d.nvisphi > 0 && !diff,
             "diffphi"    => d.nvisphi > 0 && diff,
             "diffvisamp" => d.nvisamp > 0 && poly,
             "flux"       => a.flux)
    return join((string(k, "=", v ? 1 : 0) for (k, v) in flags), ",")
end

"""
    shell_observables() -> String

Which observables the current dataset holds, as `"v2=1,t3amp=1,…"`.

The Image and Model panels both gate their tick boxes on this. Empty when no dataset is
loaded, which is the state where every box should be off rather than optimistically on.
"""
function shell_observables()
    e = current_dataset()
    e === nothing && return ""
    try
        return observable_flags_string(e.data)
    catch err
        console!(_shell(), "! could not read observables: " * _cause(err))
        return ""
    end
end

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
