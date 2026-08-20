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
    fits    :: Vector{Any}      # completed fits, newest last
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
    else
        # Live: assignment only. See src/livecanvas.jl for why nothing may be created here.
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
function shell_set_view(kind::AbstractString, color::AbstractString, logy::Bool = false)
    sh = _shell()
    old_kind, old_color, old_log = sh.kind, sh.color, sh.logy
    sh.kind  = Symbol(kind)
    sh.color = Symbol(color)
    # Silently ignored where it would be wrong rather than refused: a phase takes both signs,
    # and a log axis on one drops half the points without saying so.
    sh.logy  = logy && Symbol(kind) in LOG_Y_KINDS
    # An exception thrown out of a QML callback terminates the application. Report it in the
    # status bar and roll back instead, so a dataset without T3 or a colour mode that does
    # not apply is an inconvenience rather than a crash.
    console!(sh, "plot $(kind), coloured by $(color)" * (sh.logy ? ", log y" : ""); kind = :cmd)
    try
        sh.status = refresh_plot!(sh)
        console!(sh, sh.status)
    catch err
        sh.kind, sh.color, sh.logy = old_kind, old_color, old_log
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
                         optimiser::AbstractString, maxeval::Integer)
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

    console!(sh, "> fit_model(model, $(free), data; " *
                 "weights = $(weights), method = :$(opt))")
    t = @elapsed r = try
        if opt == "lsqfit"
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

    push!(sh.fits, (; result = r, optimiser = opt, seconds = t, free))
    for (k, v) in zip(r.list_free_params, r.x_opt)
        console!(sh, @sprintf("    %-22s = %.6g", k, v))
    end
    line = @sprintf("chi2r %.4f   ndof %d   %d free   %.2f s",
                    r.chi2r, r.ndof, length(free), t)
    console!(sh, "  " * line)
    return line
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
    return line
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
                           from_previous::Bool = false)
    sh = _shell()
    e = current_dataset(sh)
    e === nothing && return "no dataset loaded"

    weights = imaging_weights(; v2, t3amp, t3phi)
    all(iszero, weights) && return "no observables selected"

    setup = ImagingSetup(; nx = Int(nx), pixsize = Float64(pixsize),
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
    console!(sh, "> reconstruct(x_start, data, ft; maxiter = $(Int(maxiter)))  " *
                 "# $(Int(nx))×$(Int(nx)) at $(round(Float64(pixsize); digits=4)) mas, " *
                 "$startdesc, " * regtext)
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
                          regularizers = regs, ft, x_start = xprev)
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
