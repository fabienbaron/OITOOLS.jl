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
                         gridn::Integer = 60)
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
        console!(sh, "> chi2_map(model, $(free), data, $(repr(gp1)), $(repr(gp2)); " *
                     "n1 = $(Int(gridn)), n2 = $(Int(gridn)), weights = $(weights))")
    else
        console!(sh, "> fit_model(model, $(free), data; " *
                     "weights = $(weights), method = :$(opt))")
    end

    themap = nothing
    t = @elapsed r = try
        if opt == "grid"
            themap = chi2_map(md, free, data, gp1, gp2;
                              lb, ub, weights, n1 = Int(gridn), n2 = Int(gridn))
            FitResult(themap)
        elseif opt == "lsqfit"
            fit_model_lsqfit(md, free, data; lb, ub, weights, constraints = cons,
                             maxIter = Int(maxeval))
        elseif opt == "ultranest"
            fit_model_ultranest(md, free, data; lb, ub, weights, constraints = cons,
                                verb = false, cornerplot = false)
        else
            fit_model(md, free, data; lb, ub, weights, priors, constraints = cons,
                      method = Symbol(opt), maxeval = Int(maxeval))
        end
    catch err
        msg = "! fit failed: " * _cause(err); console!(sh, msg); return msg
    end

    # The map belongs to the fit that produced it, not to the shell: two grid fits over
    # different parameter pairs both stay reachable, and selecting an older fit can redraw it.
    themap === nothing || sh.chi2map === nothing || update_chi2_map!(sh.chi2map, themap)
    push!(sh.fits, (; result = r, optimiser = opt, seconds = t, free, map = themap))
    for (k, v) in zip(r.list_free_params, r.x_opt)
        console!(sh, @sprintf("    %-22s = %.6g", k, v))
    end
    line = @sprintf("chi2r %.4f   ndof %d   %d free   %.2f s",
                    r.chi2r, r.ndof, length(free), t)
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
    shell_fit_rows() -> String

Completed fits as `label\toptimiser\tchi2r\tndof\tnfree\taic\tbic` rows, newest first.

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
        push!(rows, join((string("fit ", length(sh.fits) - i + 1), f.optimiser,
                          string(round(r.chi2r; digits = 4)), string(r.ndof), string(k),
                          string(round(aic; digits = 1)), string(round(bic; digits = 1))), "\t"))
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
    rows = try
        model_rows(m.dict, m.free;
                   lb = isempty(m.lb) ? nothing : m.lb,
                   ub = isempty(m.ub) ? nothing : m.ub)
    catch err
        console!(_shell(), "! model: " * _cause(err); kind = :err)
        return ""
    end
    return join((join((r.component, r.param, r.key, _mode_name(r.mode),
                       _num(r.value), r.expr, _num(r.lb), _num(r.ub),
                       string(r.fitindex), string(r.atbound)), "\t")
                 for r in rows), "\n")
end

_mode_name(m) = m === PARAM_FREE ? "free" : m === PARAM_EXPR ? "expr" : "fixed"
_num(x) = isfinite(x) ? string(x) : "0.0"

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
    shell_component_kinds() -> String

The component kinds "+ component" offers, as `kind\tlabel` lines.

Taken from the parser's own table rather than a list written here, so a kind OITOOLS gains --
or loses -- cannot leave the two disagreeing about what can be built.
"""
function shell_component_kinds()
    labels = Dict(:point => "point source", :ud => "uniform disk", :gaussian => "Gaussian",
                  :ldlin => "linear limb darkening", :ldquad => "quadratic limb darkening",
                  :ldpow => "power-law limb darkening", :ring => "uniform ring",
                  :gaussian_ring => "Gaussian ring", :crescent => "crescent",
                  :resolved => "fully resolved")
    order = [:point, :ud, :gaussian, :ldlin, :ldquad, :ldpow,
             :ring, :gaussian_ring, :crescent, :resolved]
    return join((string(k, '\t', get(labels, k, string(k)))
                 for k in order if haskey(_ANALYTIC_PARAM_SUFFIXES, k)), "\n")
end

# Starting values for a freshly added component. Not physics -- somewhere sane for a fit to
# begin, in the spirit of `default_bounds`. Sizes in mas; the outer of a pair is put clear of
# the inner so a new ring is a ring rather than an error.
const _COMPONENT_SEEDS = Dict{String,Float64}(
    "diamout" => 4.0, "fwhmout" => 4.0, "crout" => 4.0,
    "croff" => 0.5, "crprojang" => 0.0,
    "u" => 0.3, "w" => 0.1, "alpha" => 0.15, "resolved" => 1.0)

_component_seed(suffix) = get(_COMPONENT_SEEDS, suffix, 2.0)

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
    haskey(_ANALYTIC_PARAM_SUFFIXES, k) || return "! unknown component kind: " * String(kind)
    pre = nm * ","
    any(startswith(key, pre) for key in keys(m.dict)) &&
        return "! there is already a component called " * nm

    gk = get(_KIND_GEOMETRY, k, "")
    suffixes = copy(_ANALYTIC_PARAM_SUFFIXES[k])
    isempty(gk) || gk in suffixes || pushfirst!(suffixes, gk)
    for s in suffixes
        m.dict[pre * s] = _component_seed(s)
    end
    # Every component carries a flux fraction: it is what `:point` is identified BY, and what
    # makes a second component mean anything relative to the first.
    m.dict[pre * "f"] = 1.0
    console!(sh, "  added $(nm): $(k), keys " * join((pre * s for s in suffixes), ", "))
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
    r = try
        reconstruct_image(e.data, setup; weights, maxiter = Int(maxiter),
                          regularizers = regs, ft, x_start = xprev, options = opts)
    catch err
        msg = "! reconstruction failed: " * _cause(err)
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
