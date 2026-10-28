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
    # :none, :model or :image -- plot (model - data)/sigma for the SELECTED observable rather
    # than the observable itself. A modifier on the current plot, not a panel of its own: a
    # residual belongs on the same axes as the quantity it is a residual OF.
    residual :: Symbol
    # :none, :model or :image -- draw the model's or the reconstruction's OBSERVABLES on top of
    # the data, rather than replacing them. Independent of `residual`: seeing the curve through
    # the points and seeing what is left over are different questions.
    overlay  :: Symbol
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
    # Which channel of `imaging` the panel is showing. A cube reconstruction has one image per
    # wavelength bin and the canvas draws one at a time, so this is the answer to "one of
    # what?" for everything that reads the result: the canvas, the residual overlay, the model
    # observables, the saved PNG. Always 1 for a grey run.
    #
    # Last in the struct, and with a default below, because the field was added after three
    # places were already building a state positionally.
    imchannel :: Int
end

"""
    ShellState(fields...)   # every field but `imchannel`

Build a shell without naming the view state, which defaults to the first channel.

Written as a `Vararg` on the field count rather than by listing the fields: three places build
a state positionally, and a constructor that spelled out thirty-odd names would be one rename
away from being silently wrong about which argument is which.
"""
ShellState(args::Vararg{Any,N}) where {N} =
    N == fieldcount(ShellState) - 1 ? ShellState(args..., 1) :
                                      throw(MethodError(ShellState, args))

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

    > load_dataset!(session, "Iota_Peg4T.oifits")
      1520 V2 points, 1 wav-bin x 1 time-bin
    ! cannot plot t3phi/wav: BoundsError ...
"""
function console!(sh::ShellState, text::AbstractString; kind::Symbol = :info)
    mark = kind === :cmd ? "> " : kind === :err ? "! " : "  "
    stamp = Dates.format(Dates.now(), "HH:MM:SS")
    # A marker already in the text is not repeated. Callers build an error as `"! ..."` because
    # that leading `!` is what the panels test to decide something failed, and they then pass
    # the same string here — which read as `! ! could not write ...` in the pane. Only `>` and
    # `!` are de-duplicated: an info line's marker is two spaces, and a message that genuinely
    # starts with an indent should keep it.
    dup = (kind === :cmd || kind === :err) && startswith(text, mark)
    push!(sh.console, string(stamp, " ", dup ? "" : mark, text))
    length(sh.console) > CONSOLE_MAX && deleteat!(sh.console, 1:(length(sh.console) - CONSOLE_MAX))
    return text
end

"""Everything the console pane shows, oldest first."""
shell_console() = join(_shell().console, "\n")

const SHELL = Ref{Union{Nothing,ShellState}}(nothing)

# Which clock the Gantt's x axis is labelled in: `:lst`, `:utc` or `:local`. Module state rather
# than a ShellState field because the OFFSCREEN renderer needs it too -- "Save PNG" rebuilds the
# chart in a second figure, and a saved image labelled in a different clock from the one on
# screen is a bug report waiting to happen. The panel sends it on every draw regardless.
const GANTT_TIME_SYSTEM = Ref{Symbol}(:lst)
# Whether the Gantt's target is a calibrator, for the same reason.
const GANTT_CAL = Ref(false)

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
        return string(e.name, " — ", source_label(e), " — ", sh.kind,
                      " (", n, " of ", total, " ", grouping_noun(sh.kind), "s)")
    else
        # Live: assignment only. See src/livecanvas.jl for why nothing may be created here.
        show_panels!(sh.canvas, false)
        # Residuals and the model/image overlay, both derived from the same two producers.
        # A failure is RETURNED, not just recorded: writing it to `sh.status` here and then
        # falling through to the usual "N points" line at the end overwrote it, so a residual
        # that could not be computed showed the plain data plot and said nothing at all.
        rset = residual_set(sh)
        rset isa String && return rset
        oset = overlay_set(sh)
        oset isa String && return oset
        sh.info    = update_canvas!(sh.canvas, d, sh.kind; color = sh.color,
                                    logscale = sh.logy,
                                    residual = residual_for(rset, sh.kind),
                                    overlay  = residual_for(oset, sh.kind),
                                    source   = source_label(e))
        sh.plotobj = sh.canvas.scatterplot
    end
    suffix = sh.residual === :none ? "" : " — " * String(sh.residual) * " residuals"
    suffix *= sh.overlay === :none ? "" : " + " * String(sh.overlay) * " overplotted"
    return string(e.name, " — ", source_label(e), " — ", sh.kind,
                  " (", length(sh.info), " points)", suffix)
end

"""
    residual_set(sh) -> NamedTuple or an error String

`(model - data)/sigma` per observable, from the model or from the last reconstruction.

Both producers key their result identically -- `v2`, `t3amp`, `t3phi`, `visamp`, `visphi` --
so the caller picks a field by the plot kind and never learns which made it.
"""
function residual_set(sh::ShellState)
    sh.residual === :none && return nothing
    if sh.residual === :model
        r = _current_residuals(sh; quiet = true)
        r isa String && return r
        return r.res
    end
    im = sh.imaging
    im === nothing && return "! no reconstruction yet — run one in Imaging first"
    e = current_dataset(sh)
    e === nothing && return "! no dataset loaded"
    return try
        ft = setup_ft(e.data, im.setup.nx, im.setup.pixsize; mode = String(im.setup.mode))
        # `ImagingResult.image` is always Float64, and the panel's plans follow `readoifits`,
        # which is Float32 by default. Without the cast this is a MethodError inside
        # `image_to_vis` — which is exactly why the image overlay drew nothing.
        #
        # The channel the panel is SHOWING, and only that one. `image_to_obs(::Matrix, ::Matrix
        # {OIdata}, …)` reduces to `ft[1], data[1]`, so handing it a plane against a binned
        # dataset silently answers for bin 1 whatever the user is looking at.
        image_to_residuals(to_ft_precision(result_plane(im, sh.imchannel), ft), ft, e.data)
    catch err
        "! image residuals not available: " * _cause(err)
    end
end

"""
    overlay_set(sh) -> NamedTuple or an error String

The model's or the reconstruction's OBSERVABLES, keyed like the residual set so `residual_for`
picks a kind out of either.

`model_to_obs` and `image_to_obs` are what the two engines already use to compare against the
data, so an overlay drawn from them is the same prediction the chi2 is computed from — not a
second evaluation that could disagree with it.
"""
function overlay_set(sh::ShellState)
    sh.overlay === :none && return nothing
    e = current_dataset(sh)
    e === nothing && return "! no dataset loaded"
    if sh.overlay === :model
        m = _model()
        return try
            fm = parse_model(Dict{String,Any}(String(k) => v for (k, v) in m.dict),
                             String.(m.free); nB_workspace = 1)
            x = Float64[_row_value(m, k) for k in m.free]
            model_to_obs(fm, x, e.data[1, 1])
        catch err
            "! model observables not available: " * _cause(err)
        end
    end
    im = sh.imaging
    im === nothing && return "! no reconstruction yet — run one in Imaging first"
    return try
        ft = setup_ft(e.data, im.setup.nx, im.setup.pixsize; mode = String(im.setup.mode))
        # `ImagingResult.image` is always Float64, and the panel's plans follow `readoifits`,
        # which is Float32 by default. Without the cast this is a MethodError inside
        # `image_to_vis` — which is exactly why the image overlay drew nothing.
        image_to_obs(to_ft_precision(result_plane(im, sh.imchannel), ft), ft, e.data)
    catch err
        "! image observables not available: " * _cause(err)
    end
end

"""
    residual_for(res, kind) -> Vector or nothing

The residual vector matching a plot kind, or `nothing` where the kind has none.

`t3amp_max` and `t3phi_max` differ from their plain forms only in which baseline they are drawn
against, so they share a residual. `uv`, `flux` and the differential views have none: the first
is geometry, and the others are not what `model_to_residuals` reports.
"""
function residual_for(res, kind::Symbol)
    res === nothing && return nothing
    k = kind === :t3amp_max ? :t3amp : kind === :t3phi_max ? :t3phi : kind
    return hasproperty(res, k) ? getproperty(res, k) : nothing
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

"""
    shell_about() -> String

What the About panel shows, as `key\tvalue` lines: `version` (OITOOLS), `julia`, and `license`,
the licence named at the top of the package's `LICENSE.md` with its version line when it has one
(blank when the file is not shipped).
"""
function shell_about()
    lic = ""
    try
        f = joinpath(pkgdir(OITOOLS), "LICENSE.md")
        if isfile(f)
            l = filter(!isempty, strip.(readlines(f)))
            lic = isempty(l) ? "" : length(l) > 1 && startswith(l[2], "Version") ?
                  "$(l[1]), $(l[2])" : String(l[1])
        end
    catch
    end
    return join(["version\t$(pkgversion(OITOOLS))", "julia\t$(VERSION)", "license\t$(lic)"], "\n")
end

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

# The window-walk locator: with `OITOOLSGUI_WALK` naming a folder, a QML timer asks for the
# next request written there and writes back where the named control is, so a walk driven from
# outside (xdotool) clicks controls by name rather than by pixel. Inert otherwise.

"The pending walk request (`<OITOOLSGUI_WALK>/request`, removed once read), or \"\"."
function shell_walk_request()
    dir = get(ENV, "OITOOLSGUI_WALK", "")
    f = joinpath(dir, "request")
    (isempty(dir) || !isfile(f)) && return ""
    s = read(f, String)
    rm(f; force = true)
    return String(strip(s))
end

"Write the answer to a walk request to `<OITOOLSGUI_WALK>/reply`, whole or not at all."
function shell_walk_reply(s::AbstractString)
    dir = get(ENV, "OITOOLSGUI_WALK", "")
    isempty(dir) && return ""
    tmp = joinpath(dir, "reply.tmp")
    write(tmp, s)
    mv(tmp, joinpath(dir, "reply"); force = true)
    return ""
end

"""
    shell_open(path) -> String

Load what the file dialog chose: one OIFITS file, or several -- one path per line -- read as
one dataset. Returns a status string for the UI.
"""
function shell_open(path::AbstractString)
    sh = _shell()
    ps = [strip_file_url(String(l)) for l in split(path, '\n') if !isempty(strip(l))]
    isempty(ps) && return "! no file chosen"
    p = length(ps) == 1 ? ps[1] : ps
    console!(sh, "load_dataset!(session, $(_literal(p)))"; kind = :cmd)
    try
        load_dataset!(sh.session, p; warn = false, verbose = false)
        sh.current = length(sh.session.datasets)
        sh.status = refresh_plot!(sh)
        console!(sh, sh.status)
    catch err
        sh.status = "could not load $(basename(ps[1]))" *
                    (length(ps) > 1 ? " + $(length(ps) - 1) more: " : ": ") * _cause(err)
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
shell_dataset_names() = join([e.name * " — " *
                              (e.derived_from === nothing ? source_label(e) :
                                                            e.derived_from * " › filtered")
                              for e in _shell().session.datasets], "\n")

"""Index (1-based) of the dataset every perspective is looking at; 0 with none loaded."""
shell_current_dataset() = _shell().current

# ── the Filter section (Explore) ──────────────────────────────────────────────────────

# The dataset a filter applies to: the current one, or the unfiltered dataset it came from.
# A filter's settings describe the whole cut, so applying them again replaces the cut rather
# than cutting what an earlier one left.
function _filter_root(sh::ShellState)
    e = current_dataset(sh)
    for _ in 1:length(sh.session.datasets)
        (e === nothing || e.derived_from === nothing) && break
        i = findfirst(x -> x.name == e.derived_from, sh.session.datasets)
        i === nothing && break
        e = sh.session.datasets[i]
    end
    return e
end

"""
    shell_filter_choices() -> String

What the Filter section offers, from the unfiltered dataset and its file, as `kind\t…` lines:
`extent_wav\tlo\thi` (µm), `extent_bl` (Mλ), `extent_blm` (m), `night\tlo:hi\tlabel` (the MJD
window that keeps one night, from `_nights`),
`insname\tname`, `target\tname`, `station\tname`, `baseline\tA-B`, `triangle\tA-B-C`, and
`has\tobservable` for each of v2, t3amp, t3phi, visamp, visphi, flux the data hold.
"""
function shell_filter_choices()
    sh = _shell()
    e = _filter_root(sh)
    e === nothing && return ""
    data = vec(e.data)
    lines = String[]
    allof(f) = reduce(vcat, (Float64.(getfield(d, f)) for d in data))
    lam, bl, mjd = allof(:uv_lam), allof(:uv_baseline), allof(:uv_mjd)
    g(v) = @sprintf("%.4g", v)
    if !isempty(lam)
        push!(lines, join(("extent_wav", g(minimum(lam) * 1e6), g(maximum(lam) * 1e6)), "\t"))
        push!(lines, join(("extent_bl", g(minimum(bl) / 1e6), g(maximum(bl) / 1e6)), "\t"))
        blm = bl .* lam
        push!(lines, join(("extent_blm", g(minimum(blm)), g(maximum(blm))), "\t"))
    end
    ns = _nights(mjd)
    for (i, n) in enumerate(ns)
        a, b = _night_window(ns, i)
        push!(lines, string("night\t", @sprintf("%.5f:%.5f", a, b), "\t",
                            Date(1858, 11, 17) + Day(floor(Int, n.lo)), " UT  ·  ", n.n, " uv points"))
    end
    hdrs = [h for h in (try OITOOLS.load_oifits(f) catch; nothing end
                        for f in (isempty(e.files) ? [e.path] : e.files)) if h !== nothing]
    for n in unique(String[i.insname for h in hdrs for i in h.instr]); push!(lines, "insname\t" * n); end
    for n in unique(String[t.target for h in hdrs for t in h.target]); push!(lines, "target\t" * n);  end
    d1 = first(data)
    name(k) = d1.sta_name[findfirst(==(k), d1.sta_index)]
    used = Set{Int}(); pairs = Set{Vector{Int}}(); tris = Set{Vector{Int}}()
    for d in data
        for c in eachcol(d.v2_sta_index);  push!(pairs, sort(collect(c))); union!(used, c); end
        for c in eachcol(d.vis_sta_index); push!(pairs, sort(collect(c))); union!(used, c); end
        for c in eachcol(d.t3_sta_index);  push!(tris, sort(collect(c)));  union!(used, c); end
    end
    for (o, n) in (("v2", :nv2), ("t3amp", :nt3amp), ("t3phi", :nt3phi), ("visamp", :nvisamp),
                   ("visphi", :nvisphi), ("flux", :nflux))
        any(d -> getfield(d, n) > 0, data) && push!(lines, "has\t" * o)
    end
    for k in sort(collect(used)); push!(lines, "station\t" * name(k)); end
    for p in sort(collect(pairs)); push!(lines, "baseline\t" * join(name.(p), "-")); end
    for t in sort(collect(tris));  push!(lines, "triangle\t" * join(name.(t), "-")); end
    return join(lines, "\n")
end

"""
    shell_filter_preview(spec) -> String

What the Filter section's settings would remove, without applying them: "Filtered 37/494 V²
· 12/260 T3", over every bin (a T3 point holds its amplitude and phase together, so one count). When the plot shows the dataset being filtered (its first bin,
as Explore plots it), those points are drawn hollow and grey on it.

The instrument, target, table and flag choices are read options, decided when the file is
read: they are not counted here, only on Apply.
"""
function shell_filter_preview(spec::AbstractString)
    sh = _shell()
    root = _filter_root(sh)
    root === nothing && return ""
    k = filter_spec_kwargs(spec)
    if k isa String
        sh.canvas === nothing || (sh.canvas.prevpoints[] = Makie.Point2f[])
        return k
    end
    data = vec(root.data)
    idxs = try
        [OITOOLS.set_data_filter(d; k.filter...) for d in data]
    catch err
        sh.canvas === nothing || (sh.canvas.prevpoints[] = Makie.Point2f[])
        return "! " * _cause(err)
    end
    parts = String[]
    for (label, field, slot) in (("V²", :v2, 3), ("T3", :t3phi, 4), ("VIS", :visamp, 2), ("flux", :flux, 5))
        n = sum(length(getfield(d, field)) for d in data)
        n > 0 || continue
        push!(parts, string(sum(length(ix[slot]) for ix in idxs), "/", n, " ", label))
    end
    msg = "Filtered " * join(parts, " · ")
    haskey(k.read, :binning) && (msg *= " · read again averaging every $(k.read[:binning]) channels")
    isempty(k.drop) || (msg *= " · without " * join(string.(k.drop), ", "))
    _set_filter_preview!(sh, root, first(idxs))
    return msg
end

# The removed points of the plot on screen, when it shows `root`'s first bin as data (not
# panels, not residuals, whose y is not the data's).
function _set_filter_preview!(sh::ShellState, root, idx)
    c = sh.canvas
    c === nothing && return nothing
    if current_dataset(sh) !== root || sh.panels || sh.residual !== :none
        c.prevpoints[] = Makie.Point2f[]
        return nothing
    end
    d = root.data[1, 1]
    kind = sh.kind
    pts = Makie.Point2f[]
    if kind === :uv
        keep(n, bad) = setdiff(1:n, bad)
        kept = Set{Int}()
        union!(kept, d.indx_v2[keep(length(d.indx_v2), idx[3])])
        t3k = keep(length(d.indx_t3_1), idx[4])
        union!(kept, d.indx_t3_1[t3k], d.indx_t3_2[t3k], d.indx_t3_3[t3k])
        union!(kept, d.indx_vis[keep(length(d.indx_vis), idx[2])])
        gone = union(setdiff(1:d.nuv, kept), idx[1])
        u = Float64.(d.uv[1, gone]) ./ 1e6; v = Float64.(d.uv[2, gone]) ./ 1e6
        pts = vcat(Makie.Point2f.(u, v), Makie.Point2f.(-u, -v))
    elseif haskey(OBS_SPECS, kind)
        spec = OBS_SPECS[kind]
        slot = spec.y === :v2 ? 3 : spec.y in (:t3phi, :t3amp) ? 4 :
               spec.y === :flux ? 5 : 2
        x = Float64.(getfield(d, spec.x)) .* spec.xscale
        y = Float64.(getfield(d, spec.y))
        r = sort(idx[slot])
        pts = Makie.Point2f.(x[r], y[r])
        sh.logy && filter!(p -> isfinite(p[2]) && p[2] > 0, pts)
    end
    c.prevpoints[] = pts
    return nothing
end

"""
    shell_apply_filter(spec) -> String

Filter the dataset on screen -- or, when that is itself a filter's result, the dataset it came
from -- with the Filter section's settings, into a new dataset shown as "data1 › filtered".
Logged, so the exported script makes the same cut.
"""
function shell_apply_filter(spec::AbstractString)
    sh = _shell()
    root = _filter_root(sh)
    root === nothing && return "! no dataset loaded"
    k = filter_spec_kwargs(spec)
    k isa String && return k
    (isempty(k.read) && isempty(k.filter) && isempty(k.prep) && isempty(k.drop)) &&
        return "! nothing to filter: every setting is at its default"
    nlog = length(sh.session.log)
    entry = try
        filter_dataset!(sh.session, root.name; read = k.read, prep = k.prep, drop = k.drop,
                        spec = String(spec), k.filter...)
    catch err
        msg = "! could not filter $(root.name): " * _cause(err)
        console!(sh, msg; kind = :err)
        return msg
    end
    for l in sh.session.log[nlog+1:end]
        console!(sh, l.code; kind = :cmd)
    end
    sh.current = length(sh.session.datasets)
    sh.canvas === nothing || (sh.canvas.prevpoints[] = Makie.Point2f[])
    sh.status = try
        refresh_plot!(sh)
    catch err
        "filtered into $(entry.name), but cannot plot it: " * _cause(err)
    end
    console!(sh, sh.status)
    return sh.status
end

"""The Filter section's settings that made the current dataset; empty for a file as read."""
function shell_filter_spec()
    e = current_dataset(_shell())
    return e === nothing ? "" : e.filterspec
end

"""Remove the filter preview from the plot (the section was reset or closed)."""
function shell_filter_preview_clear()
    sh = _shell()
    sh.canvas === nothing || (sh.canvas.prevpoints[] = Makie.Point2f[])
    return ""
end

"""
    shell_dataset_summary() -> String

The current dataset's [`dataset_summary`](@ref) as `label\tvalue` lines, for Explore's Dataset
section; `""` with none loaded.
"""
function shell_dataset_summary()
    e = current_dataset(_shell())
    e === nothing && return ""
    rows = try
        dataset_summary(e)
    catch err
        ["summary" => "! " * _cause(err)]
    end
    return join((string(first(r), "\t", last(r)) for r in rows), "\n")
end

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
    p = strip_file_url(path)
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
    drawing = Ref(false)

    # The selection tools: while one is chosen, a left drag draws a box or a lasso instead of
    # zooming, and its release selects what is inside. With none chosen nothing here runs.
    Makie.on(events.mouseposition, priority = 100) do _
        (drawing[] && sh.canvas !== nothing) || return Makie.Consume(false)
        p = Makie.Point2f(Makie.mouseposition(scene))
        o = sh.canvas.seloutline[]
        if SEL_TOOL[] === :box && !isempty(o)
            a = o[1]
            sh.canvas.seloutline[] = [a, Makie.Point2f(p[1], a[2]), p, Makie.Point2f(a[1], p[2]), a]
        else
            # `[p]`, not `p`: a Point2f is itself a vector, and vcat would splice in its two
            # coordinates rather than append the point.
            sh.canvas.seloutline[] = vcat(o, [p])
        end
        return Makie.Consume(true)
    end

    Makie.on(events.mousebutton, priority = 100) do ev
        pos = Tuple(Float64.(events.mouseposition[]))
        if SEL_TOOL[] !== :none && ev.button == Makie.Mouse.left && sh.canvas !== nothing &&
           (ev.action == Makie.Mouse.press ? Makie.is_mouseinside(scene) : drawing[])
            if ev.action == Makie.Mouse.press
                drawing[] = true
                sh.canvas.seloutline[] = [Makie.Point2f(Makie.mouseposition(scene))]
            else
                drawing[] = false
                poly = sh.canvas.seloutline[]
                sh.canvas.seloutline[] = Makie.Point2f[]
                select_polygon!(sh, poly)
            end
            return Makie.Consume(true)
        end
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

# ── the selection tools (Explore) ──────────────────────────────────────────────────

# `:none`, `:box` or `:lasso`; and the selection: the dataset and plot kind it was made on,
# which observable (`set_data_filter`'s discard-list slot: 1 uv, 2 VIS, 3 V², 4 T3, 5 flux)
# and the indices into it.
const SEL_TOOL  = Ref(:none)
const SELECTION = Ref{Any}(nothing)

"""
    shell_set_select_tool(mode) -> String

`"none"`, `"box"` or `"lasso"`. With a tool chosen, a left drag on the plot selects instead of
zooming; with none, the plot's zoom, pan, pick and reset are untouched.
"""
function shell_set_select_tool(mode::AbstractString)
    sh = _shell()
    t = Symbol(mode)
    t in (:none, :box, :lasso) || return "! the selection tools are none, box and lasso"
    SEL_TOOL[] = t
    if sh.canvas !== nothing
        t === :none ? Makie.activate_interaction!(sh.canvas.axis, :rectanglezoom) :
                      Makie.deactivate_interaction!(sh.canvas.axis, :rectanglezoom)
    end
    return ""
end

"""
    select_polygon!(sh, poly) -> Int

Select the points of the plot on screen inside `poly` (data coordinates; a box is its four
corners), ring them, and keep them as the selection. What the mouse handler calls on release,
and what a test calls without a mouse. Returns how many were selected.
"""
function select_polygon!(sh::ShellState, poly::AbstractVector)
    e = current_dataset(sh)
    (e === nothing || sh.panels || sh.residual !== :none || length(poly) < 3) && return 0
    slot, idx, pts = selection_indices(e.data[1, 1], sh.kind, poly)
    SELECTION[] = isempty(idx) ? nothing : (; dataset = e.name, kind = sh.kind, slot, idx)
    sh.canvas === nothing || (sh.canvas.selpoints[] = pts)
    return length(idx)
end

const _SLOT_NAMES = ("uv point", "VIS point", "V² point", "T3 point", "flux point")
const _SLOT_KEYS  = ("remove_uv", "remove_vis", "remove_v2", "remove_t3", "remove_flux")

"""What is selected, for the Select section: "12 V² points on data1", or empty."""
function shell_selection_text()
    s = SELECTION[]
    s === nothing && return ""
    return string(length(s.idx), " ", _SLOT_NAMES[s.slot], length(s.idx) == 1 ? "" : "s",
                  " selected on ", s.dataset)
end

"""Forget the selection and its rings."""
function shell_clear_selection()
    sh = _shell()
    SELECTION[] = nothing
    sh.canvas === nothing || (sh.canvas.selpoints[] = Makie.Point2f[]; sh.canvas.seloutline[] = Makie.Point2f[])
    return ""
end

# A point's identity across filtered copies of a dataset, so a selection made on a filtered view
# names the same points in the unfiltered dataset the Filter section applies to: its uv
# coordinates (one pair, or a triangle's three), wavelength and stations. Not its epoch alone:
# a Float32 MJD can be the same for a whole night. Flux has no uv, so its value stands in.
function _point_keys(d, slot)
    uvat(k) = (d.uv[1, k], d.uv[2, k])
    slot == 1 && return [(uvat(i), d.uv_lam[i]) for i in 1:size(d.uv, 2)]
    slot == 2 && return [(uvat(d.indx_vis[i]), d.vis_lam[i], sort(d.vis_sta_index[:, i]))
                         for i in eachindex(d.visamp)]
    slot == 3 && return [(uvat(d.indx_v2[i]), d.v2_lam[i], sort(d.v2_sta_index[:, i]))
                         for i in eachindex(d.v2)]
    slot == 4 && return [(uvat(d.indx_t3_1[i]), uvat(d.indx_t3_2[i]), uvat(d.indx_t3_3[i]),
                          d.t3_lam[i]) for i in eachindex(d.t3phi)]
    return [(d.flux_mjd[i], d.flux_lam[i], d.flux_sta_index[i], d.flux[i]) for i in eachindex(d.flux)]
end

"""
    shell_selection_spec() -> String

The selection as Filter-section lines -- `remove_v2\t3,17,42` -- indexed into the unfiltered
dataset the section filters, so "Add to filter" can add them to its hand-picked points. A
selection made on a filtered copy is matched to the source point by point. `""` with nothing
selected; a line starting with "!" when it cannot be (a binned source, whose bins a single
plotted bin does not index).
"""
function shell_selection_spec()
    sh = _shell()
    s = SELECTION[]
    s === nothing && return ""
    i = findfirst(e -> e.name == s.dataset, sh.session.datasets)
    i === nothing && return "! the dataset the selection was made on is gone"
    e = sh.session.datasets[i]
    root = e
    for _ in 1:length(sh.session.datasets)
        root.derived_from === nothing && break
        j = findfirst(x -> x.name == root.derived_from, sh.session.datasets)
        j === nothing && break
        root = sh.session.datasets[j]
    end
    size(root.data) == (1, 1) ||
        return "! hand-picked points need a single-bin dataset; $(root.name) has $(length(root.data)) bins"
    idx = if root === e
        s.idx
    else
        rk = _point_keys(root.data[1, 1], s.slot)
        where = Dict(k => n for (n, k) in enumerate(rk))
        sk = _point_keys(e.data[1, 1], s.slot)
        sort!([where[sk[i]] for i in s.idx if haskey(where, sk[i])])
    end
    return string(_SLOT_KEYS[s.slot], "\t", join(idx, ","))
end

"""
    shell_export_oifits(path) -> String

Write the current dataset -- filtered or binned, as it is -- to an OIFITS file through
[`write_oifits`](@ref). The original file is never the target: writing over the file a dataset
was read from would replace the data with a copy of itself less whatever was filtered out.
"""
function shell_export_oifits(path::AbstractString)
    sh = _shell()
    e = current_dataset(sh)
    e === nothing && return "! no dataset loaded"
    p = strip_file_url(path)
    any(f -> OITOOLS._same_file(p, f), [e.path; e.files]) &&
        return "! that is a file the dataset was read from; choose another name"
    try
        write_oifits(p, e.data; overwrite = true)
    catch err
        msg = "! could not write OIFITS: " * _cause(err); console!(sh, msg; kind = :err); return msg
    end
    console!(sh, "> write_oifits(\"" * p * "\", " * e.name * ".data; overwrite = true)")
    return "saved " * basename(p)
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
    _collapse_cr(s) -> String

Show a carriage-return progress line the way a terminal does: keep only what follows the last
`\r` on each line.

Samplers redraw a progress bar in place with `\r` rather than emitting a line per update.
Captured verbatim into a `TextArea` that string is ONE line, thousands of characters wide — and
the console does not wrap (`OutputConsole` sets `NoWrap` on purpose, so columns of numbers stay
in columns). The result is a horizontal rule where the output should be. Nautilus is where this
was first noticed; anything with a progress bar does it.
"""
_collapse_cr(s::AbstractString) =
    join((last(split(line, '\r')) for line in split(String(s), '\n')), '\n')

"""
    shell_fit_model(free_lines, constraint_lines, prior_lines,
                    v2, t3amp, t3phi, visamp, visphi, flux, diffvis, optimiser, maxeval,
                    gridp1, gridp2, gridn, gridn2,
                    uncertainty, nboot, bootmode, granularity, seed, threaded, fit_rv,
                    sampler_opts) -> String

Fit the model the panel is editing (the session's live model) to the current dataset. Returns
a one-line summary.

The optimiser choice decides more than speed. NLopt takes the constraints as real nonlinear
constraints, so they hold at the optimum; Levenberg-Marquardt and nested sampling have no such
machinery and apply a soft penalty that a steep χ² can overrule. That difference is stated in
the console rather than left to be discovered.

`uncertainty` is the panel's Uncertainty selector, and decides the `±` the fits table shows:

| choice | σ |
|---|---|
| `"none"` | none |
| `"jacobian"` | Levenberg-Marquardt's `stderror` (none from the other optimisers) |
| `"bootstrap"` | [`bootstrap_fit`](@ref) after the fit, started from its optimum, with `nboot`, `bootmode`, `granularity`, `seed` and `threaded`; the inner fitter is LM when the optimiser is, NLopt otherwise |
| `"posterior"` | the standard deviation of the posterior (nested sampling, the ensemble sampler, Pigeons) |

The two samplers beside nested sampling: `"ensemble"` ([`fit_model_ensemble`](@ref)) and
`"pigeons"` ([`fit_model_pigeons`](@ref), when Pigeons is loaded). `sampler_opts` gives their
settings as `key\tvalue` lines: `sweeps` (the ensemble), `rounds`, `chains`, `explorer`
(`slice` or `mala`) and `reference` (`fixed`, `variational`, or `stabilised`: a variational leg
and a fixed one, as many chains each), Pigeons'; `seed` is the panel's.

Bootstrap replicates are fitted without constraints or priors, so a bootstrap is refused when
the panel has either, as it is after a grid search or nested sampling.
"""
function shell_fit_model(free_lines::AbstractString,
                         constraint_lines::AbstractString, prior_lines::AbstractString,
                         v2::Bool, t3amp::Bool, t3phi::Bool,
                         visamp::Bool, visphi::Bool, flux::Bool, diffvis::Bool,
                         optimiser::AbstractString, maxeval::Integer,
                         gridp1::AbstractString = "", gridp2::AbstractString = "",
                         gridn::Integer = 60, gridn2::Integer = 0,
                         uncertainty::AbstractString = "jacobian", nboot::Integer = 200,
                         bootmode::AbstractString = "replacement",
                         granularity::AbstractString = "config", seed::Integer = 1,
                         threaded::Bool = true, fit_rv::Bool = false,
                         sampler_opts::AbstractString = "")
    sh = _shell()
    e = current_dataset(sh)
    e === nothing && return "! no dataset loaded"

    # A copy: the fit runs on a worker, and the panel may edit the live model meanwhile.
    md = copy(_model(sh).dict)
    isempty(md) && return "! no model"
    free, lb, ub = parse_free_lines(free_lines)
    isempty(free) && return "! no free parameters"

    cons = try
        parse_constraint_lines(constraint_lines)
    catch err
        msg = "! constraint: " * _cause(err); console!(sh, msg); return msg
    end
    priors  = parse_prior_lines(prior_lines)
    weights = fitting_weights(; v2, t3amp, t3phi, visamp, visphi, flux, diffvis)
    all(iszero, weights) && return "! no observables selected"

    data = e.data[1, 1]
    opt  = String(optimiser)
    so = _parse_options(sampler_opts)
    samplers = ("nested", "ensemble", "pigeons")
    opt == "pigeons" && Base.get_extension(OITOOLS, :OITOOLSPigeonsExt) === nothing &&
        return "! parallel tempering needs Pigeons"
    rv = fit_rv ? RVSET[] : nothing
    if fit_rv
        rv === nothing && return "! no radial-velocity table open"
        (opt == "grid" || opt in samplers) &&
            return "! radial velocities are fitted with Levenberg-Marquardt or NLopt"
        isempty(_orbit_name(_model(sh))) &&
            return "! radial velocities need a component on an orbit (Geometry › orbit)"
        String(uncertainty) == "bootstrap" && opt != "lsqfit" &&
            return "! a bootstrap with radial velocities refits with Levenberg-Marquardt"
    end

    unc = String(uncertainty)
    unc in ("none", "jacobian", "bootstrap", "posterior") ||
        return "! unknown uncertainty: " * unc
    bmode, bgran = Symbol(bootmode), Symbol(granularity)
    bfitter = opt == "lsqfit" ? :lsqfit : :nlopt
    if unc == "bootstrap"
        (opt == "grid" || opt in samplers) &&
            return "! bootstrap refits with Levenberg-Marquardt or NLopt; choose one of those"
        (isempty(cons) && isempty(priors)) ||
            return "! bootstrap replicates are fitted without constraints or priors; remove them or choose another uncertainty"
        bmode in (:replacement, :halfsample, :weights, :pmoired) ||
            return "! unknown bootstrap mode: " * String(bootmode)
        bgran in (:config, :epoch, :point) ||
            return "! unknown bootstrap granularity: " * String(granularity)
        nboot >= 2 || return "! nboot must be at least 2"
    end

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
    elseif opt == "ensemble"
        console!(sh, "> fit_model_ensemble(model, $(free), data; weights = $(weights), " *
                     "nsweep = $(_optint(so, "sweeps", 2000)), seed = $(Int(seed)))")
    elseif opt == "pigeons"
        console!(sh, "> fit_model_pigeons(model, $(free), data; weights = $(weights), " *
                     "n_rounds = $(_optint(so, "rounds", 10)), n_chains = $(_optint(so, "chains", 10)), " *
                     "explorer = :$(_optstr(so, "explorer", "slice"))" *
                     (_optstr(so, "reference", "stabilised") == "variational" ? ", variational = true" :
                      _optstr(so, "reference", "stabilised") == "stabilised" ?
                      ", n_chains_variational = $(_optint(so, "chains", 10))" : "") *
                     ", seed = $(Int(seed)))")
    else
        console!(sh, "> fit_model(model, $(free), data; " *
                     "weights = $(weights), method = :$(opt)" * (rv === nothing ? "" : ", rv") * ")")
    end
    unc == "bootstrap" &&
        console!(sh, "> bootstrap_fit(model, $(free), data; nboot = $(Int(nboot)), " *
                     "mode = :$(bmode), granularity = :$(bgran), fitter = :$(bfitter), " *
                     "weights = $(weights), seed = $(Int(seed)), threaded = $(threaded))")

    # `verb = true` throughout: the console under the results table shows what the optimiser
    # printed, so it has to print. The fit runs on a worker — see `GuiJob` — so the window
    # stays alive and its output can be read while it is produced.
    sh.fitlog = ""
    optsym = Symbol(opt)
    GRID_LIVE[] = nothing
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
                              n2 = gridn2 > 1 ? Int(gridn2) : Int(gridn),
                              on_grid = g -> (GRID_LIVE[] = g))
            FitResult(themap)
        elseif opt == "lsqfit"
            fit_model_lsqfit(md, free, data; lb, ub, weights, constraints = cons, rv,
                             maxIter = Int(maxeval), verb = true)
        elseif opt == "nested"
            # Whichever sampler is loaded. `nested_backend()` picks it; the panel says which,
            # so a logz in the fit table can be attributed to the sampler that produced it.
            fit_model_nested(md, free, data; lb, ub, weights, constraints = cons,
                             verb = true, cornerplot = false)
        elseif opt == "ensemble"
            fit_model_ensemble(md, free, data; lb, ub, weights, constraints = cons,
                               nsweep = _optint(so, "sweeps", 2000), seed = Int(seed),
                               verb = true)
        elseif opt == "pigeons"
            fit_model_pigeons(md, free, data; lb, ub, weights, constraints = cons,
                              n_rounds = _optint(so, "rounds", 10),
                              n_chains = _optint(so, "chains", 10),
                              explorer = Symbol(_optstr(so, "explorer", "slice")),
                              variational = _optstr(so, "reference", "stabilised") == "variational",
                              n_chains_variational = _optstr(so, "reference", "stabilised") == "stabilised" ?
                                                     _optint(so, "chains", 10) : 0,
                              seed = Int(seed), multithreaded = threaded, verb = true)
        else
            fit_model(md, free, data; lb, ub, weights, priors, constraints = cons, rv,
                      method = optsym, maxeval = Int(maxeval), verb = true)
        end
        boot = nothing
        sigma = if unc == "bootstrap"
            # From the optimum, so the full-data fit `bootstrap_fit` starts with is a few
            # iterations rather than the whole fit again.
            mb = copy(md)
            for (k, v) in zip(r.list_free_params, r.x_opt); mb[k] = v; end
            boot = bootstrap_fit(mb, free, data; nboot = Int(nboot), mode = bmode,
                                 granularity = bgran, fitter = bfitter, lb, ub, weights, rv,
                                 method = bfitter == :nlopt ? optsym : :LD_LBFGS,
                                 maxeval = Int(maxeval), seed = Int(seed), threaded,
                                 verb = true)
            boot.sigma
        elseif unc == "jacobian" && opt != "grid"
            _jacobian_sigma(r)
        elseif unc == "posterior"
            _posterior_sigma(r)
        else
            nothing
        end
        label = unc == "bootstrap" ? string(opt, " + bootstrap ", Int(nboot)) : opt
        # The map travels with the result rather than being drawn here: `update_chi2_map!`
        # touches Makie, which belongs to the GUI thread.
        return (; result = r, map = themap, seconds = time() - t0, optimiser = label, free,
                  chi2_start, sigma, boot)
    end)
    return "fitting…"
end

# ── radial velocities (Modeling) ──────────────────────────────────────────────────────

# The radial-velocity table the Modeling perspective fits with the interferometry, when one is
# open (`read_rv`).
const RVSET = Ref{Any}(nothing)

# The component on an orbit in the model, or "".
function _orbit_name(m)
    for cn in OITOOLS._component_names(m.dict)
        get(m.dict, "$cn,x", nothing) == "orbit" && return cn
    end
    return ""
end

"""
    shell_rv_open(path) -> String

Read a radial-velocity table (`read_rv`: CSV `MJD,RV,sigma,star`) for the fits. When the model
has a component on an orbit without velocity elements, adds them (`c,orb Ka`, `Kb`, `gamma`,
fixed, at 10, 10 and 0 km/s). Returns the table's summary, or a message beginning with `!`.
"""
function shell_rv_open(path::AbstractString)
    sh = _shell()
    p = strip_file_url(String(path))
    rv = try
        OITOOLS.read_rv(p)
    catch err
        return console!(sh, "! " * _cause(err); kind = :err)
    end
    RVSET[] = rv
    console!(sh, "rv = read_rv($(_literal(p)))"; kind = :cmd)
    m = _model(sh); cn = _orbit_name(m)
    if !isempty(cn) && !(haskey(m.dict, "$cn,orb M") && haskey(m.dict, "$cn,orb q"))
        added = String[]
        for (k, v) in ("Ka" => 10.0, "Kb" => 10.0, "gamma" => 0.0)
            key = "$cn,orb $k"
            haskey(m.dict, key) || (m.dict[key] = v; push!(added, key))
        end
        isempty(added) || console!(sh, "  added " * join(added, ", ") * " for the velocities")
    end
    return shell_rv_summary()
end

"The open radial-velocity table: `name: N velocities (a, b), MJD range`, or \"\"."
function shell_rv_summary()
    rv = RVSET[]
    rv === nothing && return ""
    na = count(==(:a), rv.star)
    return @sprintf("%s: %d velocities (%d a, %d b), MJD %.1f–%.1f", basename(rv.path),
                    length(rv), na, length(rv) - na, minimum(rv.mjd), maximum(rv.mjd))
end

"The model's component on an orbit, or \"\"."
shell_orbit_component() = _orbit_name(_model())

"Close the radial-velocity table."
shell_rv_close() = (RVSET[] = nothing; "")

"""
    finish_fit!(sh, res) -> String

Record a finished fit and draw whatever it produced. **GUI thread only**: `update_chi2_map!` is
a Makie call.
"""
function finish_fit!(sh::ShellState, res)
    if !res.ok
        msg = "! fit failed: " * _cause(res.err)
        sh.fitlog = isempty(res.output) ? msg : _collapse_cr(res.output) * "\n" * msg
        console!(sh, msg)
        return msg
    end
    f = res.result
    r = f.result
    r isa LsqFitResult && r.chi2_rv > 0 &&
        console!(sh, @sprintf("  χ² interferometry %.1f, radial velocities %.1f (%d)",
                              r.chi2 - r.chi2_rv, r.chi2_rv, RVSET[] === nothing ? 0 : length(RVSET[])))
    # The grid search prints nothing -- it is a comprehension over `chi2_flat`, not an optimiser
    # with a trace -- so the console gets the one thing worth saying about it.
    sh.fitlog = isempty(strip(res.output)) && f.optimiser == "grid" ?
        "grid search evaluates χ² on a fixed grid and has no iteration trace to report.\n" *
        "the surface itself is the output: see the χ² map below." :
        _collapse_cr(res.output)

    # The map belongs to the fit that produced it, not to the shell: two grid fits over
    # different parameter pairs both stay reachable, and selecting an older fit can redraw it.
    f.map === nothing || sh.chi2map === nothing || update_chi2_map!(sh.chi2map, f.map)
    # A companion search chooses no uncertainty; its best node keeps whatever its fit carries.
    sigma = haskey(f, :sigma) ? f.sigma : _result_sigma(r)
    push!(sh.fits, (; result = r, optimiser = f.optimiser, seconds = f.seconds,
                      free = f.free, map = f.map, chi2_start = f.chi2_start, sigma,
                      boot = get(f, :boot, nothing)))
    for (j, (k, v)) in enumerate(zip(r.list_free_params, r.x_opt))
        console!(sh, sigma === nothing || !isfinite(sigma[j]) ?
                     @sprintf("    %-22s = %.6g", k, v) :
                     @sprintf("    %-22s = %.6g ± %.3g", k, v, sigma[j]))
    end
    for line in _reversed_rings(_model(sh).dict, r.list_free_params, r.x_opt)
        console!(sh, "  " * line)
    end
    if f.map isa CompanionSearch
        COMPANION_LIVE[] = nothing
        s = f.map; b = s.best
        line = @sprintf("%s %+.3f, %+.3f mas (sep %.3f mas, PA %.2f°)   %s,f %.5f (Δm %.2f)   χ²r %.4f, %.4f without   %.1fσ   %d nodes, %.1f s",
                        isempty(s.refined) ? "best node" : "best refined minimum", b.x, b.y, b.sep, b.pa, s.companion, b.f, b.dmag, b.chi2r, s.chi2r_null,
                        s.nsigma, s.nnodes, f.seconds)
        console!(sh, "  " * line)
        for (k, q) in enumerate(s.refined)
            console!(sh, @sprintf("    minimum %d: %+.3f ± %.3f, %+.3f ± %.3f mas   f %.5f   χ²r %.4f   (from node %+.2f, %+.2f)",
                                  k, q.x, q.stderror[end - 1], q.y, q.stderror[end], q.f, q.chi2r,
                                  q.start...))
        end
        return line
    end
    if f.map isa MultiStartResult
        m = f.map
        nmin = count(==(:minimum), m.status)
        line = @sprintf("%d starts, %d kept, %d minima; best χ²r %.4f from %d start(s)   %.1f s",
                        size(m.starts, 2), nmin, length(m.minima), r.chi2r,
                        length(m.minima[1].starts), f.seconds)
        console!(sh, "  " * line)
        for (k, q) in enumerate(m.minima[1:min(10, end)])
            console!(sh, @sprintf("    minimum %d: χ²r %.4f   %d start(s)", k, q.fit.chi2r,
                                  length(q.starts)))
        end
        return line
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

`p1\tp2\tbest1\tbest2\tchi2r\tkind` for the most recent fit that drew a map, or `""` when the
last fit did not. `kind` is `grid` for a χ² map over two parameters and `companion` for a
companion search, whose axes are the companion's position. QML uses it to decide whether the
map panel has anything to show, and what its caption and toggle are.
"""
function shell_chi2_map_info()
    sh = _shell()
    isempty(sh.fits) && return ""
    f = last(sh.fits)
    m = hasproperty(f, :map) ? f.map : nothing
    m === nothing && return ""
    if m isa CompanionSearch
        return join(("x", "y", string(m.best.x), string(m.best.y), string(m.best.chi2r),
                     "companion"), "\t")
    end
    if m isa MultiStartResult
        isempty(m.minima) && return ""
        names = [m.explored; filter(n -> n ∉ m.explored, [m.list_free_params; m.held])]
        v = MULTISTART_VIEW[]
        a = v.p1 in names ? v.p1 : names[1]
        b = v.p2 in names ? v.p2 : names[min(2, length(names))]
        q = m.minima[clamp(v.selected, 1, length(m.minima))]
        val = OITOOLS._ms_values(m, q)
        return join((a, b, string(val[a]), string(val[b]), string(q.fit.chi2r), "minima"), "\t")
    end
    b = argmin(m.chi2)
    return join((m.p1, m.p2, string(m.v1[b[1]]), string(m.v2[b[2]]),
                 string(minimum(m.chi2) / m.ndof), "grid"), "\t")
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

"""
The fitted values as `key=value` lines, so the table can adopt them. After a companion search,
the best node: its fitted parameters and the companion's position, which the grid set.
"""
function shell_fit_values()
    sh = _shell()
    isempty(sh.fits) && return ""
    f = last(sh.fits)
    r = f.result
    kv = [string(k, "=", v) for (k, v) in zip(r.list_free_params, r.x_opt)]
    if f.map isa MultiStartResult && !isempty(f.map.held)
        # the held parameters of the selected minimum, which its start set
        m = f.map.minima[clamp(MULTISTART_VIEW[].selected, 1, length(f.map.minima))]
        v = OITOOLS._ms_values(f.map, m)
        append!(kv, [string(h, "=", v[h]) for h in f.map.held])
    end
    if f.map isa CompanionSearch
        # A refined minimum's fit has the position among its parameters already.
        for (a, v) in ((",x", f.map.best.x), (",y", f.map.best.y))
            k = f.map.companion * a
            k in r.list_free_params || push!(kv, string(k, "=", v))
        end
    end
    return join(kv, "\n")
end

"""
    _fit_correlation(f) -> (names, matrix) or nothing

The correlation matrix of a fits entry, from the same source as its `±`: the bootstrap's after
a bootstrap; LM's covariance after a Jacobian fit (and for a refined companion minimum); the
nested posterior's. `nothing` when the entry has no uncertainty -- "skip", NLopt, the grid, an
unrefined companion search -- since a correlation without a covariance has nothing behind it.
"""
function _fit_correlation(f)
    r = f.result
    get(f, :boot, nothing) isa BootstrapResult && return (f.boot.list_free_params, f.boot.correlation)
    f.sigma === nothing && return nothing
    if r isa LsqFitResult
        s = r.stderror
        return (r.list_free_params, r.covar ./ (s * s'))
    elseif hasproperty(r, :posterior) && r.posterior isa AbstractMatrix && size(r.posterior, 1) > 1
        return (r.list_free_params, cor(r.posterior))
    end
    return nothing
end

"""
    shell_fit_correlation() -> String

The latest fit's correlation matrix for the panel under the fits table: the parameter names on
the first line, then one row of the matrix per line, tab-separated, to two decimals. `""` when
the fit has no uncertainty to correlate (see `_fit_correlation`).
"""
function shell_fit_correlation()
    sh = _shell()
    isempty(sh.fits) && return ""
    c = _fit_correlation(last(sh.fits))
    c === nothing && return ""
    names, R = c
    rows = [join((@sprintf("%.2f", R[i, j]) for j in axes(R, 2)), "\t") for i in axes(R, 1)]
    return join(vcat(join(names, "\t"), rows), "\n")
end

"""
    _fit_samples(f) -> (samples, names, best) or nothing

The samples a corner plot of a fits entry is drawn from: a bootstrap's kept replicates, or a
nested sampler's posterior, with the fit's optimum to mark. `nothing` for a fit with neither.
"""
function _fit_samples(f)
    b = get(f, :boot, nothing)
    b isa BootstrapResult && return (b.samples[b.mask, :], b.list_free_params, b.x_opt)
    r = f.result
    hasproperty(r, :posterior) && r.posterior isa AbstractMatrix && size(r.posterior, 1) > 1 &&
        return (r.posterior, r.list_free_params, r.x_opt)
    return nothing
end

# The corner PNG on screen, so "Save corner…" writes the picture that was looked at.
const CORNER_PNG = Ref("")

"""
    shell_corner_plot(width, height, px_per_unit = 1) -> String

Draw the latest fit's corner plot (PairPlots, through `plot_corner_makie`) to a PNG and return
its path, or a message beginning with `!`. The figure is `width` × `height` units, written at
`px_per_unit` pixels per unit.

The figure is never put on screen as a Makie scene: no plot may be created once the QML window
holds the GL context, and a corner plot is a grid of new axes. It is built offscreen and saved,
as "Save PNG" does, and the panel shows the file. Each call writes a new file name, so the
panel's `Image` cannot show a cached older one.
"""
function shell_corner_plot(width::Real = 800, height::Real = 800, px_per_unit::Real = 1)
    sh = _shell()
    isempty(sh.fits) && return "! no fit yet"
    Base.get_extension(OITOOLS, :OITOOLSPairPlotsExt) === nothing &&
        return "! corner plots need PairPlots"
    src = _fit_samples(last(sh.fits))
    src === nothing &&
        return "! the latest fit has no samples: fit with Uncertainty › bootstrap, or nested sampling"
    S, names, best = src
    fig = try
        pd = plot_corner_makie(S, names; truth = best)
        Makie.resize!(pd.figure, max(320, round(Int, width)), max(320, round(Int, height)))
        pd.figure
    catch err
        msg = "! could not draw the corner plot: " * _cause(err)
        console!(sh, msg; kind = :err); return msg
    end
    dir = joinpath(tempdir(), "oitoolsgui-corner")
    mkpath(dir)
    isempty(CORNER_PNG[]) || rm(CORNER_PNG[]; force = true)
    p = joinpath(dir, string("corner-", getpid(), "-", time_ns(), ".png"))
    try
        Makie.save(p, fig; backend = GLMakie, px_per_unit = Float64(px_per_unit))
    catch err
        msg = "! could not write the corner plot: " * _cause(err)
        console!(sh, msg; kind = :err); return msg
    end
    CORNER_PNG[] = p
    console!(sh, "> plot_corner_makie(" *
                 (get(last(sh.fits), :boot, nothing) isa BootstrapResult ? "boot" : "result") *
                 "; truth = fit.x_opt)")
    return p
end

"""
    shell_save_corner(path) -> String

Copy the corner plot on screen to `path` (`.png` added if missing). `""` on success.
"""
function shell_save_corner(path::AbstractString)
    sh = _shell()
    isfile(CORNER_PNG[]) || return "! no corner plot drawn yet"
    p = String(strip(strip_file_url(String(path))))
    isempty(p) && return "! no file name"
    endswith(lowercase(p), ".png") || (p = p * ".png")
    try
        cp(CORNER_PNG[], p; force = true)
    catch err
        return "! could not write " * p * ": " * _cause(err)
    end
    console!(sh, "> save(\"" * p * "\", plot_corner_makie(...).figure)")
    return ""
end

# ── the companion search (Model › Optimiser) ─────────────────────────────────────────

"""
    shell_companion_candidates() -> String

The model's components that can be the searched companion, newline-separated: those with
`f`, `x` and `y` keys. The search moves one over a grid of positions.
"""
function shell_companion_candidates()
    m = _model()
    names = String[]
    for k in keys(m.dict)
        i = findfirst(',', k)
        i === nothing && continue
        c = k[1:prevind(k, i)]
        (c in names || !all(haskey(m.dict, c * s) for s in (",f", ",x", ",y"))) && continue
        push!(names, c)
    end
    return join(sort(names), "\n")
end

# NaN from an empty field, which means "the search's own default".
_or_default(v) = (v isa Real && isfinite(v)) ? Float64(v) : nothing

"""
    shell_companion_grid(companion, step, rmin, rmax, mask) -> String

How large a search will be, before it runs: "1024 nodes · step 1.81 mas · r 1.81–36.6 mas",
from [`companion_grid`](@ref) with the same defaults the search takes. An empty field (NaN)
takes the default, and the readout shows what that default is.
"""
function shell_companion_grid(companion, step, rmin, rmax, mask)
    sh = _shell()
    e = current_dataset(sh)
    e === nothing && return "! no dataset loaded"
    cn = String(companion)
    isempty(cn) && return "! the model has no component with f, x and y to search for"
    g = try
        companion_grid(_model(sh).dict, e.data[1, 1]; companion = cn, step = _or_default(step),
                       rmin = _or_default(rmin), rmax = _or_default(rmax),
                       mask = something(_or_default(mask), 0.01))
    catch err
        return "! " * _cause(err)
    end
    return @sprintf("%d nodes · step %.3g mas · r %.3g–%.3g mas", length(g.nodes), g.step,
                    g.rmin, g.rmax)
end

# `x` to two significant figures, rounded by `r` (floor or ceil).
function _round2(x::Real, r)
    (isfinite(x) && x > 0) || return Float64(x)
    p = 10.0^(floor(log10(x)) - 1)
    return r(x / p) * p
end

# A value from `_round2`, written with its two significant figures: 2.0, 1.2, 0.85, 43.
function _fmt2(x::Real)
    (isfinite(x) && x > 0) || return @sprintf("%g", x)
    return @sprintf("%.*f", max(0, 1 - floor(Int, log10(x))), x)
end

"""
    shell_companion_defaults(companion) -> String

The grid the search would take by default, as the dialog prefills it: `"step\trmin\trmax"`,
in mas, to two significant figures: the step rounded down, so the grid is no coarser than the
default λ/(2 B_max), and the radii outwards (the inner one down, the outer one up), which only
widens the annulus. `""` when there is no grid to compute (no dataset, or a
model the grid cannot be laid around).
"""
function shell_companion_defaults(companion)
    sh = _shell()
    e = current_dataset(sh)
    (e === nothing || isempty(String(companion))) && return ""
    g = try
        companion_grid(_model(sh).dict, e.data[1, 1]; companion = String(companion))
    catch
        return ""
    end
    # The step rounded down, so that it is never coarser than λ/(2 B_max).
    return join((_fmt2(_round2(g.step, floor)), _fmt2(_round2(g.rmin, floor)),
                 _fmt2(_round2(g.rmax, ceil))), "\t")
end

"""
    shell_companion_search(free_lines, v2, t3amp, t3phi, visamp, visphi, flux, diffvis,
                           companion, step, rmin, rmax, mask, fmax, fitter,
                           refit_primary, refine) -> String

Run [`companion_search`](@ref) on the panel's model and the current dataset, on a worker, as a
fit is run. The panel's free parameters are fitted at every node, less the companion's
position, which the grid sets. `fitter` is the per-node fitter: `"bradi"` or `"halley"` when
only the companion's flux is free, `"vmlmb"` for any free set (and the only one when more than
the flux is). With `refit_primary` false the primary's free parameters are held at their
current values: only the companion's are fitted at a node, and the significance is against this
primary as it stands. `refine` is `companion_search`'s: the best `refine` local minima fitted
again with the position free (0 for none). The map lands in the χ² map panel, the best node (or
refined minimum) in the fits table, and "Adopt" writes it into the model.
"""
function shell_companion_search(free_lines::AbstractString,
                                v2::Bool, t3amp::Bool, t3phi::Bool,
                                visamp::Bool, visphi::Bool, flux::Bool, diffvis::Bool,
                                companion, step, rmin, rmax, mask, fmax,
                                fitter::AbstractString = "bradi", refit_primary::Bool = true,
                                refine::Integer = 0)
    sh = _shell()
    e = current_dataset(sh)
    e === nothing && return "! no dataset loaded"
    md = copy(_model(sh).dict)
    isempty(md) && return "! no model"
    cn = String(companion)
    isempty(cn) && return "! the model has no component with f, x and y to search for"
    free, lb, ub = parse_free_lines(free_lines)
    # An empty bound is ±Inf here, and passed on it would override the search's own default.
    filter!(p -> isfinite(p.second), lb); filter!(p -> isfinite(p.second), ub)
    kx, ky, kf = cn * ",x", cn * ",y", cn * ",f"
    dropped = filter(k -> k in (kx, ky), free)
    free = filter(k -> !(k in (kx, ky)), free)
    if !refit_primary
        held = filter(k -> !startswith(k, cn * ","), free)
        isempty(held) || console!(sh, "  primary held at its current values: $(join(held, ", "))")
        free = filter(k -> startswith(k, cn * ","), free)
    end
    kf in free || return "! free $kf: the search fits the companion's flux at every node"
    (flux || diffvis) && return "! the search fits V², T3 and VIS; untick OI_FLUX and differential"
    weights = fitting_weights(; v2, t3amp, t3phi, visamp, visphi, flux, diffvis)
    all(iszero, weights) && return "! no observables selected"
    data = e.data[1, 1]
    fz = String(fitter)
    fz in ("bradi", "halley", "vmlmb") || return "! unknown per-node fitter: " * fz
    kw = (; companion = cn, free, lb, ub, fmax = something(_or_default(fmax), 0.5),
            step = _or_default(step), rmin = _or_default(rmin), rmax = _or_default(rmax),
            mask = something(_or_default(mask), 0.01), weights,
            solver = fz == "halley" ? :halley : :bradi,
            path = fz == "vmlmb" ? :general : :auto, refine = max(0, Int(refine)))
    isempty(dropped) || console!(sh, "  $(join(dropped, " and ")) set by the grid, not fitted")
    console!(sh, "> companion_search(model, data; " *
                 join((string(k, " = ", repr(v)) for (k, v) in pairs(kw)), ", ") * ")")
    sh.fitlog = ""
    COMPANION_LIVE[] = nothing
    start_job!(sh, :fit, function (_stop)
        t0 = time()
        chi2_start = try
            OITOOLS.chi2_flat(parse_model(md, free), Float64[Float64(md[k]) for k in free],
                              data; weights)
        catch
            NaN
        end
        s = companion_search(md, data; kw...,
                             on_grid = g -> (COMPANION_LIVE[] = (; g..., companion = cn)))
        return (; result = s.fit, map = s, seconds = time() - t0, optimiser = "companion search",
                  free, chi2_start)
    end)
    return "searching…"
end

# The grid search in progress: its axes and the χ² array `chi2_map` fills, set by its `on_grid`
# hook on the worker and drawn from the GUI thread by `shell_grid_progress`.
const GRID_LIVE = Ref{Any}(nothing)

"""
    shell_grid_progress() -> String

Draw the running grid search into the map panel, and say how far it is: `"done\ttotal"`, or
`""` before the grid exists. Called by the panel's job poll, on the GUI thread.
"""
function shell_grid_progress()
    sh = _shell()
    g = GRID_LIVE[]
    (g === nothing || sh.chi2map === nothing) && return ""
    done, total = update_grid_progress!(sh.chi2map, g)
    return string(done, "\t", total)
end

# The search in progress: its grid and the arrays its threads fill, set by `companion_search`'s
# `on_grid` hook on the worker and drawn from the GUI thread by `shell_companion_progress`.
const COMPANION_LIVE = Ref{Any}(nothing)

"""
    shell_companion_progress() -> String

Draw the running companion search into the map panel, and say how far it is: `"done\ttotal"`,
or `""` before the grid exists. Called by the panel's job poll, on the GUI thread, which is where
Makie has to be touched; the first call, with nothing fitted yet, puts up the axes and the
grid's extent.
"""
function shell_companion_progress()
    sh = _shell()
    g = COMPANION_LIVE[]
    (g === nothing || sh.chi2map === nothing) && return ""
    done, total = update_companion_progress!(sh.chi2map, g, g.companion)
    return string(done, "\t", total)
end

"""
    shell_companion_map_view(which) -> String

Show the last companion search's χ²ᵣ map (`"chi2r"`) or its fitted-flux map (`"flux"`).
"""
function shell_companion_map_view(which::AbstractString)
    sh = _shell()
    isempty(sh.fits) && return ""
    m = last(sh.fits).map
    (m isa CompanionSearch && sh.chi2map !== nothing) || return ""
    update_chi2_map!(sh.chi2map, m; show = which == "flux" ? :flux : :chi2r)
    return ""
end

# ── multi-start search (Model Optimization) ───────────────────────────────────────────

# One explored parameter per line, `name<TAB>kind<TAB>a<TAB>b<TAB>c`: kind `grid` (min, max,
# step), `uniform` (min, max) or `normal` (mean, σ). Lines of another kind are not explored.
function _parse_explore_lines(lines::AbstractString)
    grid, urand, nrand = Pair{String,Any}[], Pair{String,Any}[], Pair{String,Any}[]
    for l in split(lines, '\n')
        f = split(l, '\t')
        length(f) >= 2 || continue
        name, kind = String(f[1]), String(f[2])
        num(i) = (i <= length(f) ? tryparse(Float64, strip(f[i])) : nothing)
        if kind == "grid"
            v = (num(3), num(4), num(5))
            any(isnothing, v) && throw(ArgumentError("$name: grid needs min, max and step"))
            push!(grid, name => v)
        elseif kind == "uniform"
            v = (num(3), num(4))
            any(isnothing, v) && throw(ArgumentError("$name: uniform needs min and max"))
            push!(urand, name => v)
        elseif kind == "normal"
            v = (num(3), num(4))
            any(isnothing, v) && throw(ArgumentError("$name: normal needs mean and σ"))
            push!(nrand, name => v)
        end
    end
    return grid, urand, nrand
end

# One condition per line, `lhs op rhs` (`<=`, `>=`, `<`, `>`, `=`), as `constrain` takes them.
function _parse_condition_text(text::AbstractString)
    out = ModelConstraint[]
    for l in split(String(text), '\n')
        isempty(strip(l)) && continue
        m = match(r"^(.*?)\s*(<=|>=|<|>|=)\s*(.*)$", strip(l))
        m === nothing && throw(ArgumentError("condition \"$(strip(l))\": lhs op rhs, op one of <= >= < > ="))
        lhs, rhs = String(strip(m.captures[1])), String(strip(m.captures[3]))
        (isempty(lhs) || isempty(rhs)) &&
            throw(ArgumentError("condition \"$(strip(l))\": both sides are needed"))
        push!(out, ModelConstraint(lhs, m.captures[2],
                                   something(tryparse(Float64, rhs), rhs)))
    end
    return out
end

# The latest multi-start view: its two parameters ("" for the defaults) and the selected minimum.
const MULTISTART_VIEW = Ref((p1 = "", p2 = "", selected = 1))

"""
    shell_multistart(free_lines, constraint_lines, explore_lines, constrain_lines,
                     v2, t3amp, t3phi, visamp, visphi, flux, diffvis,
                     nstarts, seed, deltachi2, maxeval, fit_rv) -> String

Run [`multistart_fit`](@ref) on the current dataset and model as a job: the free rows' bounds,
the constraints, the observables and the radial velocities as the Fit button uses them;
`explore_lines` the explored parameters (see `_parse_explore_lines`), `constrain_lines` the
conditions on the starts in the constraint panel's form. The global minimum becomes the fit's
result, the search its map (the Minima view).
"""
function shell_multistart(free_lines::AbstractString, constraint_lines::AbstractString,
                          explore_lines::AbstractString, constrain_lines::AbstractString,
                          v2::Bool, t3amp::Bool, t3phi::Bool, visamp::Bool, visphi::Bool,
                          flux::Bool, diffvis::Bool, nstarts::Integer, seed::Integer,
                          deltachi2::Real, maxeval::Integer, fit_rv::Bool)
    sh = _shell()
    e = current_dataset(sh)
    e === nothing && return "! no dataset loaded"
    md = copy(_model(sh).dict)
    isempty(md) && return "! no model"
    free, lb, ub = parse_free_lines(free_lines)
    isempty(free) && return "! no free parameters"
    cons, conds, (grid, urand, nrand) = try
        parse_constraint_lines(constraint_lines), _parse_condition_text(constrain_lines),
        _parse_explore_lines(explore_lines)
    catch err
        msg = "! " * _cause(err); console!(sh, msg); return msg
    end
    isempty(grid) && isempty(urand) && isempty(nrand) && return "! explore at least one parameter"
    for (n, _) in [grid; urand; nrand]
        n in free || get(md, n, nothing) isa Real ||
            return "! $n is explored but is neither free nor a fixed number"
    end
    isempty(grid) && nstarts < 1 && return "! give the number of starts"
    weights = fitting_weights(; v2, t3amp, t3phi, visamp, visphi, flux, diffvis)
    all(iszero, weights) && return "! no observables selected"
    rv = fit_rv ? RVSET[] : nothing
    fit_rv && rv === nothing && return "! no radial-velocity table open"
    dc = isfinite(deltachi2) && deltachi2 > 0 ? Float64(deltachi2) : nothing
    data = e.data[1, 1]
    kw = (; grid, rand = urand, randn = nrand, nstarts = isempty(grid) ? Int(nstarts) : nothing,
            constrain = conds, lb, ub, weights, constraints = cons, rv, deltachi2 = dc,
            seed = Int(seed), maxIter = Int(maxeval))
    shown = [string(k, " = ", k === :rv ? "rv" : repr(v)) for (k, v) in pairs(kw)
             if !(v === nothing || (v isa AbstractVector && isempty(v)))]
    console!(sh, "> multistart_fit(model, $(free), data; " * join(shown, ", ") * ")")
    sh.fitlog = ""
    MULTISTART_VIEW[] = (p1 = "", p2 = "", selected = 1)
    start_job!(sh, :fit, function (_stop)
        t0 = time()
        chi2_start = try
            OITOOLS.chi2_flat(parse_model(md, free), Float64[Float64(md[k]) for k in free],
                              data; weights)
        catch
            NaN
        end
        r = multistart_fit(md, free, data; kw..., verbose = true)
        isempty(r.minima) && error("no fit converged with uncertainties")
        best = r.minima[1].fit
        return (; result = best, map = r, seconds = time() - t0, optimiser = "multi-start",
                  free, chi2_start, sigma = _jacobian_sigma(best))
    end)
    return "searching…"
end

"""
    shell_multistart_minima() -> String

The latest multi-start search's minima, one per line: `χ²ᵣ<TAB>starts<TAB>name=value …`, lowest
χ²ᵣ first; `""` when the latest fit was not a multi-start search.
"""
function shell_multistart_minima()
    sh = _shell()
    isempty(sh.fits) && return ""
    r = last(sh.fits).map
    r isa MultiStartResult || return ""
    names = [r.list_free_params; r.held]
    return join((string(@sprintf("%.4f", m.fit.chi2r), "\t", length(m.starts), "\t",
                        (v = OITOOLS._ms_values(r, m);
                         join((@sprintf("%s=%.6g", k, v[k]) for k in names), " ")))
                 for m in r.minima), "\n")
end

"""
    shell_multistart_select(k) -> String

Make minimum `k` (1 = the global one) of the latest multi-start search the result that Adopt
writes into the table and the χ² readout reports, and ring it in the Minima view.
"""
function shell_multistart_select(k::Integer)
    sh = _shell()
    isempty(sh.fits) && return "! no fit"
    f = last(sh.fits)
    r = f.map
    r isa MultiStartResult || return "! the latest fit is not a multi-start search"
    1 <= k <= length(r.minima) || return "! no minimum $k"
    m = r.minima[k].fit
    sh.fits[end] = merge(f, (; result = m, sigma = _jacobian_sigma(m)))
    v = MULTISTART_VIEW[]
    MULTISTART_VIEW[] = (; v.p1, v.p2, selected = Int(k))
    sh.chi2map === nothing || update_chi2_map!(sh.chi2map, r; v.p1, v.p2, selected = k)
    console!(sh, @sprintf("  minimum %d selected: χ²r %.4f, from %d start(s)", k, m.chi2r,
                          length(r.minima[k].starts)))
    return ""
end

"""
    shell_multistart_view(p1, p2) -> String

Draw the latest multi-start search over the parameters `p1`, `p2`: each start joined to the
minimum it reached, the minima coloured by χ²ᵣ.
"""
function shell_multistart_view(p1::AbstractString, p2::AbstractString)
    sh = _shell()
    isempty(sh.fits) && return ""
    r = last(sh.fits).map
    (r isa MultiStartResult && sh.chi2map !== nothing) || return ""
    MULTISTART_VIEW[] = (p1 = String(p1), p2 = String(p2), selected = MULTISTART_VIEW[].selected)
    update_chi2_map!(sh.chi2map, r; p1 = String(p1), p2 = String(p2),
                     selected = MULTISTART_VIEW[].selected)
    return ""
end

"""
    shell_multistart_names() -> String

The parameters the latest multi-start search can be drawn over, newline-separated: the explored
ones first, then the other free ones.
"""
function shell_multistart_names()
    sh = _shell()
    isempty(sh.fits) && return ""
    r = last(sh.fits).map
    r isa MultiStartResult || return ""
    return join([r.explored; filter(n -> n ∉ r.explored, [r.list_free_params; r.held])], "\n")
end

"""
    shell_multistart_pmoired(text) -> String

A PMOIRED `expl` dict (`{'grid': {'B,x': (-12, 12, 1.5)}, 'rand': {...}, 'randn': {...}}`) as
the dialog's explore lines (`name<TAB>kind<TAB>a<TAB>b<TAB>c`); a global whose name is not an
identifier renamed as the model import renames it (`MJD0-59800` → `MJD0_59800`). `"! …"` when
the text is not such a dict.
"""
function shell_multistart_pmoired(text::AbstractString)
    d = try
        Core.eval(OITOOLS, Meta.parse(OITOOLS.pmoired_to_julia(String(text))))
    catch err
        return "! not a PMOIRED expl dict: " * _cause(err)
    end
    d isa AbstractDict || return "! not a PMOIRED expl dict"
    lines = String[]
    for (pk, kind) in ("grid" => "grid", "rand" => "uniform", "randn" => "normal")
        haskey(d, pk) || continue
        d[pk] isa AbstractDict || return "! expl['$pk'] is not a dict"
        for (n, v) in sort(collect(d[pk]); by = first)
            name = String(n)
            (occursin(',', name) || occursin(r"^[A-Za-z0-9_]+$", name)) ||
                (name = replace(name, r"[^A-Za-z0-9_]" => "_"))
            push!(lines, join([name; kind; string.(collect(v))], "\t"))
        end
    end
    isempty(lines) && return "! expl has no grid, rand or randn"
    return join(lines, "\n")
end

"""
    _jacobian_sigma(r), _posterior_sigma(r), _result_sigma(r) -> Vector{Float64} or nothing

The uncertainty a fit result carries by itself: Levenberg-Marquardt's `stderror` from the
Jacobian, or the standard deviation of a nested-sampling posterior; `_result_sigma` takes
whichever there is. The plain NLopt methods carry neither, and get `nothing` rather than a `±`
on a number that came from nowhere.
"""
function _jacobian_sigma(r)
    hasproperty(r, :stderror) && r.stderror !== nothing &&
        length(r.stderror) == length(r.x_opt) || return nothing
    return Float64.(r.stderror)
end

function _posterior_sigma(r)
    hasproperty(r, :posterior) && r.posterior isa AbstractMatrix &&
        size(r.posterior, 2) == length(r.x_opt) || return nothing
    return [std(@view r.posterior[:, j]) for j in axes(r.posterior, 2)]
end

function _result_sigma(r)
    j = _jacobian_sigma(r)
    return j === nothing ? _posterior_sigma(r) : j
end

"""
    _fit_params(r, errs) -> String

A fit's optimised parameters as `name = value` pairs, with `± σ` from `errs` where it holds
a finite one (`errs === nothing` for none).

The point of a row per fit is comparing fits, and two fits that differ only in χ² are two
numbers, not two models. What separates them is where the parameters landed.
"""
function _fit_params(r, errs)
    names = r.list_free_params
    vals  = r.x_opt
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
    elseif r isa SamplerResult              # the ensemble counts its evaluations; Pigeons does not
        return (0, get(r.diagnostics, :n_eval, 0))
    end
    return (0, 0)
end

"""
    _reversed_rings(dict, free, x) -> Vector{String}

One line per profile component whose inner diameter ended above its outer one at `x`, saying
the model reads it as the ring between them. The fit is not wrong for it: the values are
reported as fitted, and this says what they mean.
"""
function _reversed_rings(dict, free, x)
    d = copy(dict)
    for (k, v) in zip(free, x)
        d[k] = v
    end
    out = String[]
    for k in keys(d)
        endswith(k, ",profile") || continue
        cn = k[1:end-length(",profile")]
        haskey(d, cn * ",diamin") || continue
        a, b = try
            OITOOLS._hankel_r_min(cn, d), OITOOLS._hankel_r_max(cn, d)
        catch
            continue
        end
        a > b || continue
        push!(out, @sprintf("%s: diamin %.4g > outer %.4g mas, read as the ring between them",
                             cn, 2a, 2b))
    end
    return sort!(out)
end

"""
    shell_fit_rows() -> String

Completed fits as
`label\toptimiser\tchi2r_start\tchi2r\twork\tndof\tnfree\tparams` rows, newest first.

`chi2r_start` is the model as it stood when the fit was launched, measured with the fit's own
weights and divided by the fit's own `ndof`, so the two chi2r columns are the same quantity and
the fit's progress is the ratio between them. It is `NaN` when the starting model could not be
evaluated, which the panel renders as an em dash.

`work` is `_fit_work` rendered with its unit, because iterations and evaluations are different
quantities and only one of them exists for any given optimiser. A fit that reports it did not
converge -- Levenberg-Marquardt stopped by its iteration limit -- has ` · not converged` after
it: its values are where it stopped, not a minimum.
A Levenberg-Marquardt minimum that holds parameters on their bounds has ` · on a bound: ` and
their names, `; `-separated: the minimum there is set by the bound, not by the data.

No AIC or BIC. Both are information criteria for choosing BETWEEN models, and this package has
a better answer for that question: `fit_model_nested` returns the log-evidence, which compares
models without assuming the likelihood is Gaussian in the parameters or that the sample is
large. A column that is a worse version of something the package already does well invites the
comparison to be made with the worse one.
"""
function shell_fit_rows()
    sh = _shell()
    rows = String[]
    for (i, f) in enumerate(reverse(sh.fits))
        r = f.result
        k = length(f.free)
        # The starting chi2 divided by the FIT'S ndof, so the two chi2r columns are the same
        # quantity and the fit's progress is the ratio between them.
        start = isfinite(f.chi2_start) ? f.chi2_start / max(r.ndof, 1) : NaN
        it, ev = _fit_work(r)
        work = it > 0 && ev > 0 ? string(it, " / ", ev) :
               it > 0           ? string(it, " it")     :
               ev > 0           ? string(ev, " ev")     : "—"
        hasproperty(r, :converged) && r.converged === false && (work *= " · not converged")
        hasproperty(r, :at_bound) && !isempty(r.at_bound) &&
            (work *= " · on a bound: " * join(r.at_bound, "; "))
        push!(rows, join((string("fit ", length(sh.fits) - i + 1), f.optimiser,
                          _numv(start), string(round(r.chi2r; digits = 4)),
                          work, string(r.ndof), string(k),
                          _fit_params(r, f.sigma)), "\t"))
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
    # The bounds set in the session, `default_bounds` for the rest (`model_rows`). The
    # dataset's first wavelength and epoch, so a chromatic expression resolves to the number
    # it actually takes there rather than to nothing.
    w, t = _display_wl_mjd()
    rows = try
        model_rows(m.dict, m.free; lb = m.lb, ub = m.ub, wl = w, mjd = t)
    catch err
        console!(_shell(), "! model: " * _cause(err); kind = :err)
        return ""
    end
    # Smearing is a setting, edited through its own menu (`shell_set_smearing`), not a row.
    return join((join((r.component, r.param, r.key, _mode_name(r.mode),
                       _numv(r.value), r.expr, _num(r.lb), _num(r.ub),
                       string(r.fitindex), string(r.atbound)), "\t")
                 for r in rows if !(r.key in _SMEAR_KEYS)), "\n")
end

const _SMEAR_KEYS = ("smear", "smear_mode")

"""
    shell_smearing() -> String

The model's bandwidth smearing as the Smearing menu shows it, `"samples\tmode"`: samples is
`off`, `auto` or a whole number of sub-wavelengths per channel; mode `incoherent` or
`coherent`.
"""
function shell_smearing()
    d = _model().dict
    sm = get(d, "smear", 0)
    n = sm == "auto" ? "auto" :
        (sm isa Real && isinteger(sm) && sm > 1) ? string(Int(sm)) : "off"
    return n * "\t" * string(get(d, "smear_mode", "incoherent"))
end

"""
    shell_set_smearing(samples, mode) -> String

Set the model's `"smear"` and `"smear_mode"` from the Smearing menu: `samples` is `off`,
`auto` or a whole number ≥ 2. Off removes both keys, so a model without smearing is saved
without them. They are model keys rather than a fit option because every evaluation against
data reads them -- a fit, χ², the residuals, the companion search -- and a saved model or an
exported script has to evaluate the same way.
"""
function shell_set_smearing(samples::AbstractString, mode::AbstractString)
    sh = _shell(); m = _model(sh)
    n, md = String(samples), String(mode)
    md in ("incoherent", "coherent") || return "! smearing mode is incoherent or coherent"
    if n == "off"
        delete!(m.dict, "smear"); delete!(m.dict, "smear_mode")
        console!(sh, "  bandwidth smearing off"; kind = :cmd)
        return ""
    end
    v = n == "auto" ? "auto" : tryparse(Int, n)
    (v == "auto" || (v isa Int && v >= 2)) ||
        return "! smearing takes auto or a whole number of sub-wavelengths ≥ 2"
    m.dict["smear"] = v
    md == "incoherent" ? delete!(m.dict, "smear_mode") : (m.dict["smear_mode"] = md)
    console!(sh, "  bandwidth smearing: $(v) sub-wavelengths per channel, $(md)"; kind = :cmd)
    return ""
end

_mode_name(m) = m === PARAM_FREE ? "free" : m === PARAM_EXPR ? "expr" : "fixed"

# Bounds: an infinite bound is no bound. It travels as JavaScript writes it, so the panel's
# `parseFloat` reads ±Infinity, leaves the field blank and sends it to the fit as no bound.
# NaN is the bound of a row that is not free, which shows none.
_num(x) = isfinite(x) || isnan(x) ? string(x) : (x > 0 ? "Infinity" : "-Infinity")

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
    shell_orbit_view() -> String

Draw the Orbit view for the model's component on an orbit: the sky track over one period, the
positions at the current dataset's nights, and the open radial velocities against phase with the
model's curves. Returns a one-line summary, or a message beginning with `!`.
"""
function shell_orbit_view()
    sh = _shell(); m = _model(sh)
    cn = _orbit_name(m)
    isempty(cn) && return "! no component on an orbit (Geometry › orbit)"
    return try
        fm = parse_model(Dict{String,Any}(String(k) => v for (k, v) in m.dict), String.(m.free);
                         nB_workspace = 1)
        vals = fm.resolver.fn(Float64[_row_value(m, k) for k in m.free])
        val(k) = (i = findfirst(==("$cn,orb $k"), fm.all_names); i === nothing ? nothing : Float64(vals[i]))
        a = something(val("a"), (val("M") !== nothing && val("plx") !== nothing) ?
                      orbit_a(val("P"), val("M"), val("plx")) : nothing, NaN)
        el = (; P = val("P"), MJD0 = val("MJD0"), e = val("e"), omega = val("omega"),
                OMEGA = val("OMEGA"), incl = val("incl"), a)
        e = current_dataset(sh)
        mjds = e === nothing ? Float64[] :
               sort(unique(round.(reduce(vcat, (Float64.(d.uv_mjd) for d in vec(e.data))); digits = 1)))
        rv = RVSET[]
        K = if val("Ka") !== nothing && val("Kb") !== nothing
            (; Ka = val("Ka"), Kb = val("Kb"), gamma = something(val("gamma"), 0.0))
        elseif val("M") !== nothing && val("q") !== nothing
            ka, kb = orbit_semiamplitudes(el.P, el.e, el.incl, val("M"), val("q"))
            (; Ka = ka, Kb = kb, gamma = something(val("gamma"), 0.0))
        else
            nothing
        end
        update_orbit_view!(ORBIT_CANVAS[], el, mjds, rv, K)
        @sprintf("%s: P %.4g d, e %.3f, a %.4g mas, i %.1f°; %d nights%s", cn, el.P, el.e, el.a,
                 el.incl, length(mjds),
                 rv === nothing ? "" : K === nothing ? "; RVs open, no Ka/Kb" : "; $(length(rv)) RVs")
    catch err
        msg = "! orbit not available: " * _cause(err)
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
    shell_model_chi2(v2, t3amp, t3phi, visamp, visphi, flux, diffvis) -> String

Reduced chi2 of the model as it stands, against the current dataset, or `""` when either is
missing. Cheap for analytic components, which is what makes it worth recomputing on every edit.

The seven flags are the panel's Observables ticks, and the number is `chi2_flat` over
`_ndof(data, weights)` with the weights `shell_fit_model` builds from the same ticks -- so after
"Adopt" it is the χ²r the fit reported, not a figure over observables the fit never saw.
"""
function shell_model_chi2(v2::Bool = true, t3amp::Bool = true, t3phi::Bool = true,
                          visamp::Bool = false, visphi::Bool = false, flux::Bool = false,
                          diffvis::Bool = false)
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
        weights = fitting_weights(; v2, t3amp, t3phi, visamp, visphi, flux, diffvis)
        all(iszero, weights) && return ""
        chi2 = model_to_chi2(fm, x, d; weights)
        @sprintf("%.4f", chi2 / OITOOLS._ndof(d, weights))
    catch err
        # A model mid-edit is not an error, it is simply not evaluable yet — but say so once,
        # or a permanently blank chi2 looks like the number rather than the model being absent.
        console!(_shell(), "  chi2 not available: " * _cause(err))
        ""
    end
end

_row_value(m, key) = get(m.dict, key, 0.0) isa Real ? Float64(m.dict[key]) : 0.0

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
rather than writing something that means a different thing on the other side. Opening or
importing a model that does not set bandwidth smearing keeps the smearing in force.
"""
function shell_open_model(path::AbstractString)
    sh = _shell()
    r = try
        read_model_file(String(path))
    catch err
        msg = "! could not read model: " * _cause(err); console!(sh, msg; kind = :err); return msg
    end
    m = _model(sh)
    old = m.dict
    m.dict = r.model; m.free = collect(r.free)
    m.lb = Dict{String,Float64}(r.lb); m.ub = Dict{String,Float64}(r.ub)
    m.name = basename(String(path))
    console!(sh, "> read_model_file(\"" * String(path) * "\")")
    console!(sh, "  $(m.name): $(length(m.dict)) keys, $(length(m.free)) free" *
                 (isempty(r.name) ? "" : ", named in the file \"" * r.name * "\""))
    _keep_smearing!(sh, old, m.dict)
    return m.name
end

# Smearing is chosen for the data as much as for the model, and a model file or a PMOIRED
# import usually does not set it: one that does not keeps the smearing in force, rather than
# switching it off unseen. One that sets it has it. The keys go into the model, so a saved
# model or an exported script evaluates as the panel did.
function _keep_smearing!(sh, old::AbstractDict, new::AbstractDict)
    (haskey(new, "smear") || !haskey(old, "smear")) && return
    for k in _SMEAR_KEYS
        haskey(old, k) && (new[k] = old[k])
    end
    console!(sh, "  bandwidth smearing kept: $(old["smear"]) sub-wavelengths per channel, " *
                 string(get(old, "smear_mode", "incoherent")); kind = :cmd)
    return
end

function shell_save_model(path::AbstractString)
    sh = _shell(); m = _model(sh)
    isempty(m.dict) && return "! no model to save"
    try
        write_model_file(String(path), m.dict; free = m.free, lb = m.lb, ub = m.ub,
                         name = shell_model_description())
    catch err
        msg = "! could not write model: " * _cause(err); console!(sh, msg; kind = :err); return msg
    end
    m.name = basename(String(path))
    console!(sh, "> write_model_file(\"" * String(path) * "\", model; free = $(m.free))")
    return ""
end

function shell_import_pmoired(path::AbstractString)
    src = try
        read(String(path), String)
    catch err
        msg = "! PMOIRED import: " * _cause(err); console!(_shell(), msg; kind = :err); return msg
    end
    return _import_pmoired!(src, basename(String(path)),
                            "pmoired_to_dict(read(\"" * String(path) * "\", String))")
end

"""
    shell_import_pmoired_text(text) -> String

Import a PMOIRED `param` dict pasted from a notebook, as [`shell_import_pmoired`](@ref) imports a
`.py` file: the model is replaced, and the name is returned, or a message beginning `!`. The
console line holds the text as a raw string, so a `\$WL` in it stays a `\$WL` when the script
runs.
"""
function shell_import_pmoired_text(text::AbstractString)
    s = strip(String(text))
    isempty(s) && return "! nothing to import: paste a PMOIRED param dict"
    return _import_pmoired!(String(s), "pasted PMOIRED model",
                            "pmoired_to_dict(raw\"\"\"" * s * "\"\"\")")
end

function _import_pmoired!(src::String, name::String, call::String)
    sh = _shell()
    md = try
        pmoired_to_dict(src)
    catch err
        msg = "! PMOIRED import: " * _cause(err); console!(sh, msg; kind = :err); return msg
    end
    m = _model(sh)
    old = m.dict
    m.dict = md; m.free = String[]; empty!(m.lb); empty!(m.ub)
    m.name = name
    console!(sh, "> " * call)
    _keep_smearing!(sh, old, m.dict)
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
    shell_model_image(nx, pixsize, wl, mjd) -> String

Render the current model and put it on the Model perspective's canvas.

The model is the one Julia holds, not a copy assembled from the panel's rows: the table does
not show every key (an image component's file is in its card), and a model rebuilt from the
rows alone would be a different model.

`wl` in metres and `mjd` in days; either ≤ 0 means "not specified", which a model that does not
reference them ignores anyway.
"""
function shell_model_image(nx::Integer, pixsize::Real, wl::Real, mjd::Real)
    sh = _shell()
    md = _model(sh).dict
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

Validate the sky a simulation would use: `"ok\t<summary>\t<pixsize>"` or `"bad\t<why>\t0"`.

Checked when the file is chosen rather than when Simulate is pressed, so a 3-D cube picked as a
grey image is caught before an observation is computed from it.

The third field is the pixel size the FITS header declares, in mas, or `0` when it declares
none — the panel fills its box with it rather than making the user retype a number the file
already carries. A model has no pixel size and reports `0` too.
"""
function shell_sim_source(kind::AbstractString, path::AbstractString)
    i = simulate_source_info(String(kind), String(path))
    return (i.ok ? "ok\t" : "bad\t") * i.summary * "\t" * string(i.pixsize)
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

# The image arguments of a simulation, as the console line writes them: `, image = "…",
# pixsize = …` for an image, nothing for a model.
_sky_args(sky) = haskey(sky, :image) ?
    ", image = $(_literal(String(sky.image))), pixsize = $(sky.pixsize)" : ""

"""
    _simulate_sky(sh, o) -> NamedTuple or String

The sky a simulation is of, as the keyword arguments `simulate` and `simulate_from_oifits`
both take: `(; image, pixsize)` for a FITS file, `(; flat_model, flat_params)` for a model.

One resolver for both entry points, so "Use existing OIFITS as template" changes the COVERAGE
and nothing else — the sky the two produce cannot differ.

An empty model path means the model the Modeling perspective is holding, which is the whole
point of the "From Modeling tab" button: no file, no round trip, and what is simulated is what
is on screen next door.
"""
function _simulate_sky(sh::ShellState, o::AbstractDict)
    kind = String(get(o, "source", "image"))
    path = String(get(o, "path", ""))

    if kind == "model"
        dict, free, value = if isempty(path)
            m = _model(sh)
            (Dict{String,Any}(String(k) => v for (k, v) in m.dict), String.(m.free),
             k -> _row_value(m, k))
        else
            isfile(path) || return "! no model file at '$(path)'"
            r = try
                read_model_file(path)
            catch err
                return "! could not read the model: " * _cause(err)
            end
            (Dict{String,Any}(String(k) => v for (k, v) in r.model), String.(r.free),
             k -> Float64(r.model[k]))
        end
        isempty(dict) && return "! no model to simulate"
        return try
            (; flat_model = parse_model(dict, free; nB_workspace = 1),
               flat_params = Float64[value(k) for k in free])
        catch err
            "! the model does not parse: " * _cause(err)
        end
    end

    # an image of any dimension, by path: `simulate` reads its keywords; the panel's pixel size,
    # pre-filled from the header, is passed as the override
    i = simulate_source_info("image", path)
    i.ok || return "! " * i.summary
    return (; image = path, pixsize = _optreal(o, "pixsize", 0.1))
end

"""
    shell_simulate(outfile, facility, telescopes, combiner, wavelength,
                   name, ra, dec, dateiso, ha_min, ha_max, step_minutes, options) -> String

Write an OIFITS of a simulated observation. `"ok\t<path>\t<summary>"`, or a message
beginning with `!`.

TWO simulators behind one button, chosen by the `copy` option:

  * **the observing calculation** (default) — `simulate` builds the uv coverage from the array,
    the hour angles and the wavelength table, and the error bars from the photon budget of the
    combiner at the target's magnitude. Everything about the observation is computed.
  * **copying a file's structure** (`copy = 1`) — `simulate_from_oifits` re-observes the sky
    through the LOADED dataset's uv points and copies its signal-to-noise. Nothing is computed
    about the array or the night; what comes back is that dataset, of a different source. It is
    the right mode for "what would this dataset look like if the star were a binary", and the
    only one whose result is comparable to the original point by point.

`options` is `key\tvalue` lines, as the Image panel's engine options are: `source`
(image | model; `cube` is read as image), `path`, `pixsize`, `image_channels` (1: the channels
of the image, or of the model's image cube, instead of the spectral setup — `wave_from_image`),
`mag`, `mag_ao`, `noise`, `debias`, `n_samples`, `seed`, `observability`, `alt_limit`,
`alt_max`, `copy`, `copy_mode`.
"""
function shell_simulate(outfile::AbstractString, facility::AbstractString,
                        telescopes::AbstractString, combiner::AbstractString,
                        wavelength::AbstractString, name::AbstractString,
                        ra::Real, dec::Real, dateiso::AbstractString,
                        ha_min::Real, ha_max::Real, step_minutes::Real,
                        options::AbstractString = "")
    sh = _shell()
    out = String(strip(String(outfile)))
    isempty(out) && return "! no output file"

    # Never write over a file the session has open — in EITHER mode. The copy path would
    # replace the very data it is simulating from, and the observing path would replace an
    # unrelated dataset the user is still looking at. The file picker asks before replacing a
    # file, but "replace?" is a question about a name; this is a question about someone's
    # observations, and the answer is no.
    for e in sh.session.datasets
        if any(f -> OITOOLS._same_file(out, f), [e.path; e.files])
            msg = "! that is $(e.name)'s own file ($(source_label(e))) — simulating into it " *
                  "would overwrite the loaded data; choose another output file"
            console!(sh, msg); return msg
        end
    end

    o = _parse_options(options)

    sky = _simulate_sky(sh, o)
    sky isa String && (console!(sh, sky; kind = :err); return sky)

    frac = get(o, "noise", "") == "fractional" ?
           FractionalNoise(v2 = _optreal(o, "frac_v2", 1) / 100, visamp = _optreal(o, "frac_visamp", 1) / 100,
                           t3amp = _optreal(o, "frac_t3amp", 1) / 100, visphi = _optreal(o, "frac_visphi", 0.5),
                           t3phi = _optreal(o, "frac_t3phi", 1), flux = _optreal(o, "frac_flux", 1) / 100) :
           nothing
    noise = frac === nothing && _optbool(o, "noise", true)
    # `systematics` draws the calibration error correlated per baseline per night, and
    # `fringe_tracker` engages one on an instrument that declares it. Both default ON, matching
    # `simulate`; both reach only the observing path, since the copy path takes its error bars
    # from the file it copies rather than from the noise model.
    systematics = _optbool(o, "systematics", true)
    ftrack      = _optbool(o, "fringe_tracker", true)
    seed  = _optint(o, "seed", 1)
    t0 = time()

    if _optbool(o, "copy", false)
        e = current_dataset(sh)
        e === nothing && return "! Use existing OIFITS as template needs a dataset loaded in Exploring"
        length(e.files) > 1 &&
            return "! Use existing OIFITS as template takes one file; $(e.name) is $(length(e.files)) files read as one"
        mode = get(o, "copy_mode", "copy_snr")
        frac === nothing ||
            return "! fractional noise needs a simulated observation; a copied structure keeps its own noise"
        try
            simulate_from_oifits(e.path, out; mode = mode, noise = noise, seed = seed, sky...)
        catch err
            why = _cause(err)
            # One failure is worth naming, because it is a property of the FILE rather than of
            # the simulation: this path re-writes the input through OIFITS.jl, and several
            # perfectly readable files cannot be written back by it (a missing CATEGORY column
            # in OI_TARGET, a missing FOV in OI_ARRAY, a CFITSIO error). Everything
            # else — a cube with the wrong number of planes, say — is the caller's to fix and
            # must not be dressed up as the file's fault.
            hint = (occursin("CFITSIO", why) || occursin("FitsError", why) ||
                    occursin("MissingColumn", why)) ?
                   " — $(basename(e.path)) is readable but cannot be written back by OIFITS.jl" : ""
            msg = "! simulating from $(basename(e.path)): " * why * hint
            console!(sh, msg); return msg
        end
        console!(sh, "simulate_from_oifits(\"$(e.path)\", \"$(out)\"; mode = \"$(mode)\"" *
                     _sky_args(sky) * (noise ? "" : ", noise = false") * ")"; kind = :cmd)
    else
        f = try
            facility_subset(read_facility_file(String(facility)),
                            [parse(Int, x) for x in split(String(telescopes)) if !isempty(x)])
        catch err
            msg = "! " * _cause(err); console!(sh, msg; kind = :err); return msg
        end
        date = try
            Date(DateTime(String(dateiso)))
        catch
            return "! not a date: $(dateiso)"
        end
        tgt = read_obs_file("default_obs")
        tgt.target = String(name); tgt.raep0 = Float64(ra); tgt.decep0 = Float64(dec)
        dates = try
            epochs_for_hour_angles(f, Float64(ra), date, Float64(ha_min), Float64(ha_max),
                                   Float64(step_minutes))
        catch err
            msg = "! " * _cause(err); console!(sh, msg; kind = :err); return msg
        end
        obs = _optbool(o, "observability", false) ?
              (; min_elevation = _optreal(o, "alt_limit", DEFAULT_ALT_LIMIT),
                 max_elevation = _optreal(o, "alt_max", DEFAULT_ALT_MAX)) : nothing
        # the image's own channels (an image, or a model holding a cube), else the setup
        imgch = _optbool(o, "image_channels", false)
        wave = try
            !imgch ? read_wave_file(String(wavelength)) :
            haskey(sky, :image) ? wave_from_image(sky.image) : wave_from_image(sky.flat_model)
        catch err
            msg = "! " * _cause(err); console!(sh, msg; kind = :err); return msg
        end
        try
            simulate(f, tgt, read_comb_file(String(combiner)), wave,
                     dates, out;
                     mag = _optreal(o, "mag", 2.0), mag_ao = _optreal(o, "mag_ao", 2.0),
                     noise = something(frac, noise), debias = _optbool(o, "debias", true),
                     systematics = systematics, fringe_tracker = ftrack,
                     n_samples = _optint(o, "n_samples", 100), seed = seed,
                     observability = obs, sky...)
        catch err
            msg = "! simulation failed: " * _cause(err)
            console!(sh, msg; kind = :err); return msg
        end
        # The command log has to be runnable, so every NON-DEFAULT switch has to appear in it.
        kw = String[]
        isempty(_sky_args(sky)) || push!(kw, _sky_args(sky)[3:end])
        frac === nothing || push!(kw, "noise = " * _frac_literal(frac))
        frac === nothing && !noise && push!(kw, "noise = false")
        systematics || push!(kw, "systematics = false")
        ftrack      || push!(kw, "fringe_tracker = false")
        obs === nothing || push!(kw, "observability = $(obs)")
        console!(sh, "simulate(facility, target, combiner, " * (imgch ? ":image" : "wavelength") *
                     ", dates, \"$(out)\"" *
                     (isempty(kw) ? "" : "; " * join(kw, ", ")) * ")   # $(length(dates)) epochs, " *
                     "$(f.ntel) telescopes"; kind = :cmd)
    end

    summary = try
        d = readoifits(out; filter_bad_data = false, verbose = false, warn = false)[1, 1]
        @sprintf("%d V², %d T3, %d channels", d.nv2, d.nt3amp, length(unique(d.uv_lam)))
    catch err
        return "! wrote $(out) but it does not read back: " * _cause(err)
    end
    console!(sh, @sprintf("  wrote %s — %s, %.1f s", out, summary, time() - t0))
    return "ok\t" * out * "\t" * summary
end

# A FractionalNoise as the code that rebuilds it: the fields away from their defaults.
function _frac_literal(f::FractionalNoise)
    d = FractionalNoise()
    kw = ["$(k) = $(_g(getfield(f, k)))" for k in fieldnames(FractionalNoise) if getfield(f, k) != getfield(d, k)]
    return "FractionalNoise(" * join(kw, ", ") * ")"
end

# The target list's own columns, in the panel's order.
const TARGET_COLUMNS = ["name", "ra", "dec", "V", "J", "H", "K", "simbad", "cal", "for"]

"""
    shell_save_targets(path, rows) -> String
    shell_load_targets(path) -> String

The Observing panel's target list as a CSV file, one target per row, with a header naming the
columns: `name, ra, dec, V, J, H, K, simbad, cal, for`, coordinates in decimal degrees, an empty
magnitude for one not known, `simbad` `1` for a row whose values are SIMBAD's, `cal` `1` for a
calibrator, `for` the names of the science targets a calibrator is for, `;`-separated. `rows` and
the result are the panel's own form, one target per line with those ten fields tab-separated,
followed by one `column=value` field for each other column the file has, which is written back.
A row of eight fields is not a calibrator, one of nine is for no target.

A file with a header is read by column name, in any order and case; one without is read by
position in that order, and a file of names alone, one per line, loads as unresolved targets.
Each returns `""` or a message beginning with `!`, the loader its rows.
"""
function shell_save_targets(path::AbstractString, rows::AbstractString)
    sh = SHELL[]
    p = String(strip(strip_file_url(String(path))))
    isempty(p) && return "! no file name"
    endswith(lowercase(p), ".csv") || (p *= ".csv")
    q(x) = occursin(r"[\",]", x) ? "\"" * replace(x, "\"" => "\"\"") * "\"" : x
    lines = [split(l, '\t') for l in split(String(rows), '\n'; keepempty = false)]
    filter!(f -> length(f) >= 8, lines)
    cal(f) = length(f) >= 9 && f[9] == "1" ? "1" : "0"
    calfor(f) = length(f) >= 10 ? String(f[10]) : ""
    # the other columns, in the order first met
    extra(f) = [String(first(split(x, '=', limit = 2))) => String(last(split(x, '=', limit = 2)))
                for x in f[11:end] if occursin('=', x)]
    keys_ = unique(k for f in lines for (k, _) in extra(f))
    try
        open(p, "w") do io
            println(io, join([TARGET_COLUMNS; q.(keys_)], ","))
            for f in lines
                ex = Dict(extra(f))
                println(io, join([q(String(f[1])); String.(f[2:8]); cal(f); q(calfor(f));
                                  [q(get(ex, k, "")) for k in keys_]], ","))
            end
        end
    catch err
        return "! could not write " * p * ": " * _cause(err)
    end
    sh === nothing || console!(sh, "wrote the target list to " * p)
    return ""
end

function shell_load_targets(path::AbstractString)
    p = String(strip(strip_file_url(String(path))))
    lines = try
        readlines(p)
    catch err
        return "! could not read " * p * ": " * _cause(err)
    end
    # One CSV field list, honouring double quotes.
    function fields(l)
        out = String[]; buf = IOBuffer(); inq = false; i = firstindex(l)
        while i <= lastindex(l)
            c = l[i]
            if inq
                if c == '"' && nextind(l, i) <= lastindex(l) && l[nextind(l, i)] == '"'
                    write(buf, '"'); i = nextind(l, i)
                elseif c == '"'
                    inq = false
                else
                    write(buf, c)
                end
            elseif c == '"'
                inq = true
            elseif c == ','
                push!(out, String(take!(buf)))
            else
                write(buf, c)
            end
            i = nextind(l, i)
        end
        push!(out, String(take!(buf)))
        return String.(strip.(out))
    end
    data = [fields(strip(l)) for l in lines if !(isempty(strip(l)) || startswith(strip(l), "#"))]
    isempty(data) && return "! no targets in " * p
    # A header is a first line with a `name` column. Each own column's position in the file
    # (0 when absent); the others are carried as `column=value`.
    head = lowercase.(first(data))
    if "name" in head
        col = [something(findfirst(==(lowercase(c)), head), 0) for c in TARGET_COLUMNS]
        others = [k for k in eachindex(head) if !(head[k] in lowercase.(TARGET_COLUMNS))]
        names = first(data)
        data = data[2:end]
    else
        col = collect(eachindex(TARGET_COLUMNS)); others = Int[]; names = String[]
    end
    rows = String[]
    for f in data
        val(c) = 0 < col[c] <= length(f) ? f[col[c]] : ""
        num(c) = tryparse(Float64, val(c)) !== nothing ? val(c) : ""
        name = val(1)
        isempty(name) && continue
        ra, dec = num(2), num(3)
        simbad = val(8) == "1" && !isempty(ra) && !isempty(dec) ? "1" : "0"
        ex = [names[k] * "=" * (k <= length(f) ? f[k] : "") for k in others]
        push!(rows, join([name, isempty(ra) ? "0" : ra, isempty(dec) ? "0" : dec,
                          num(4), num(5), num(6), num(7), simbad, val(9) == "1" ? "1" : "0",
                          val(10), ex...], "\t"))
    end
    isempty(rows) && return "! no targets in " * p
    sh = SHELL[]
    sh === nothing || console!(sh, "read $(length(rows)) targets from " * p)
    return join(rows, "\n")
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

`effective\tuser\tmarker\tzoom\tzoomuser`: the plot scale actually in force, the user's
override behind it (zero when it is computed rather than chosen), the marker size (zero for the
per-view default), the zoom factor per wheel detent in force, and the user's zoom override
(zero for the default).

Two scale numbers rather than one because they answer different questions -- the panel opens
showing what is being drawn, while "Save config" stores the override, so a scale that was
computed from this screen's DPI is not pinned onto the next screen.
"""
shell_plot_scale() = join((live_plot_scale(), PLOT_SCALE_USER[], MARKER_SIZE_USER[],
                           zoom_per_detent(), ZOOM_STEP_USER[]), '\t')

"The plots' font (`PLOT_FONT`), which Makie fixes when a figure is built."
shell_plot_font() = PLOT_FONT

"""
    shell_diagnostics() -> String

Julia's thread count, and the Mesa driver variables when any is set (`configure_graphics!`
sets them under WSL), for the settings panel's version line.
"""
function shell_diagnostics()
    parts = [string(Threads.nthreads(), Threads.nthreads() == 1 ? " thread" : " threads")]
    for k in ("GALLIUM_DRIVER", "LIBGL_ALWAYS_SOFTWARE", "MESA_LOADER_DRIVER_OVERRIDE")
        haskey(ENV, k) && push!(parts, k * "=" * ENV[k])
    end
    return join(parts, "  ·  ")
end

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
    shell_version() -> String

The OITOOLS version, and the Julia it is running on, as one line for the settings panel.

`pkgversion` reads the loaded package's own `Project.toml`, so a `dev`ed checkout reports what
that checkout declares rather than whatever a registry last published — which is the useful
answer when a bug report has to say which code was running.
"""
function shell_version()
    v = try
        string(pkgversion(@__MODULE__))
    catch
        "unknown"
    end
    return "OITOOLS v" * v * "  ·  Julia " * string(VERSION)
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
    shell_settings_path() -> String

Where the appearance config is kept, whether or not anything has been saved there.

The panel shows it on a line of its own. That line is always present, and says what it is:
empty until a button was pressed, it grew the panel when a message appeared in it, and a bare
truncated path beside the buttons reads as a stray fragment rather than as an answer.
"""
shell_settings_path() = gui_settings_file()

"""
    shell_reset_settings() -> String

Delete the saved appearance defaults, and say where they were.

The panel's reset has to remove the FILE and not merely the values on screen: settings are
applied at startup, so a window restored to its built-in look while the file still says
otherwise would come back tweaked at the next launch — which is precisely the confusion the
button exists to end.

Returns the path removed, or `""` when there was nothing saved, which is also the answer when
the settings have never been touched.
"""
function shell_reset_settings()
    sh = SHELL[]
    path = gui_settings_file()
    isfile(path) || return ""
    try
        rm(path)
        sh === nothing || console!(sh, "removed saved appearance defaults: " * path)
        return path
    catch err
        sh === nothing || console!(sh, "could not remove settings: " * _cause(err); kind = :err)
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
The laws across a disc's face: uniform, and the limb-darkening laws, as one shape wearing
different coefficients. The law menu on a disc's card and the add-component dialog's disc
entries both read this list.

They are a sub-choice, not six separate components: every one of them is a disc of some
diameter, and which law you use is a question about the atmosphere rather than about the
geometry. Offering them as separate entries beside "uniform ring" and "crescent" put two different
questions in one list and made limb darkening most of it.

`coeffs` is what the law needs BESIDES the diameter, in the parser's order.
"""
const DISC_LAWS = [
    (key = "ud",        label = "uniform",        coeffs = String[],
     note = "I(μ) = 1"),
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
What each kind the parser reads is called on screen, by `model_inspection`'s kind: the
section headers (`disk · Gaussian`) and the model's description both read it. A disc is not
here: its law names it, from `DISC_LAWS`. Written so it cannot be taken for a key --
`\$disk,gaussian` would not resolve.
"""
const KIND_LABELS = Dict{String,String}(
    "point" => "point", "gaussian" => "Gaussian", "ring" => "ring",
    "gaussian_ring" => "Gaussian ring", "crescent" => "crescent", "hankel" => "written profile",
    "image" => "image", "visfunc" => "visibility function", "resolved" => "resolved flux")

"""
    shell_kind_labels() -> String

`KIND_LABELS` as `kind\tlabel` lines, for the section headers.
"""
shell_kind_labels() = join((string(k, '\t', v) for (k, v) in sort!(collect(KIND_LABELS))), "\n")

# A kind as the description names it: a disc by its law ("power law disc"), anything else from
# `KIND_LABELS`, and an unknown kind as itself.
function _kind_label(kind::AbstractString)
    haskey(KIND_LABELS, kind) && return KIND_LABELS[kind]
    i = findfirst(l -> l.key == kind, DISC_LAWS)
    return i === nothing ? String(kind) : DISC_LAWS[i].label * " disc"
end

"""
    shell_model_description() -> String

The model as `name: kind` per component, in the model's order and joined by ", ": "star: power
law disc", or "star: uniform disc, rim: written profile". Built from what the parser reads, so it
follows a law switch or a rename the moment either happens and cannot describe a model that is
no longer there. `""` for an empty model, or one that does not parse yet.

The toolbar shows it after the model's source, and Save writes it as the file's `name`.
"""
function shell_model_description()
    m = _model()
    isempty(m.dict) && return ""
    insp = try model_inspection(m.dict) catch; return "" end
    return join((string(c.name, ": ", _kind_label(string(c.kind))) for c in insp.components), ", ")
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

    # the disc laws, from the list the law menu reads, so the two cannot disagree
    (ComponentChoice("disc", "disc", l.key, l.label, "star") for l in DISC_LAWS)...,

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

    # An image or cube from a FITS file: a reconstruction saved with Save FITS…, a ROTIR cube.
    # Not written out like the two around it, but read; added through
    # `shell_add_image_component`, which takes the file.
    ComponentChoice("fitsimage", "image (FITS file)", "fits_image", "image (FITS file)", "img"),

    # You write V(B) directly, in $BL. The other direction, and the only kind whose expression
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

The visibilities are divided by the total flux, so the sum does not change V; sharing it out
keeps every `f` readable as the component's fraction of the light, which is what the table, the
SED and a fitted value then say.

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
    problem = _new_component_problem(m, nm)
    isempty(problem) || return problem
    k = Symbol(String(kind))
    k === :fits_image && return "! an image component is added with its file: " *
                                "shell_add_image_component(name, path)"
    # Panel ids that are not parser kinds. `image_func` and `ring_profile` are both `:hankel`
    # and differ only in the keys seeded below: the first is a profile about the centre, the
    # second is one evaluated on an annulus, which is the same machinery pointed at a ring.
    k === :hankel || k === :image_func || k === :ring_profile || k === :vis_func ||
        haskey(_ANALYTIC_PARAM_SUFFIXES, k) ||
        return "! unknown component kind: " * String(kind)
    pre = nm * ","

    if k === :vis_func
        # V(B) written out. The seed is a uniform disc's own visibility, which is a real model
        # rather than a placeholder: it evaluates, it fits, and it is the thing a user is most
        # likely to be about to modify.
        # `ifelse`, because V is evaluated at B = 0: the zero-baseline point IS the total flux,
        # and `chi2_flat`'s flux term and `model_to_sed` both ask for it. 2J₁(x)/x is 0/0 there
        # and a NaN at one uv point poisons the whole chi2, so the seed shows the guard a
        # written-out visibility generally needs.
        m.dict[pre * "visfunc"] =
            "ifelse(\$BL < 1e-8, 1.0, 2*besselj1(pi*\$d*\$BL/206.264806) / " *
            "(pi*\$d*\$BL/206.264806))"
        m.dict[pre * "d"]       = _component_seed("ud")
        m.dict[pre * "f"]       = _share_flux!(m, pre * "f")
        push!(m.free, pre * "d")
        console!(sh, "  added $(nm): V(B) written in \$BL (Mλ); $(pre)d is free")
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

# Why `nm` cannot name a new component, or "".
function _new_component_problem(m, nm)
    isempty(nm) && return "! a component needs a name"
    occursin(',', nm) && return "! a component name cannot contain a comma"
    nm == GLOBAL_COMPONENT && return "! that name is reserved for globals"
    any(startswith(key, nm * ",") for key in keys(m.dict)) &&
        return "! there is already a component called " * nm
    return ""
end

"""
    shell_add_image_component(name, path) -> String

Add a FITS image or cube as the component `name`: `name,image` holds the file's absolute path
and `name,f` a share of the flux. Returns `""` or a message beginning with `!`.

The file is read here, so one the reader refuses (a 3-D file that does not say what its third
axis is, a plane with no flux) is refused before it enters the model. A file with no pixel
scale is added with `name,pixsize = 1.0` and a console line saying the number is a
placeholder: refusing it would leave no way to give the scale, which is the card's job.
Nothing is free: an image has no shape parameter, and `x`/`y` come from the position toggle.
"""
function shell_add_image_component(name, path)
    sh = _shell(); m = _model(sh)
    nm = strip(String(name))
    problem = _new_component_problem(m, nm)
    isempty(problem) || return problem
    p = abspath(expanduser(String(path)))
    hdr_ps = try fits_pixsize(p) catch; nothing end
    mi = try
        read_model_image(p; pixsize = hdr_ps === nothing ? 1.0 : nothing)
    catch err
        return "! " * _cause(err)
    end
    pre = nm * ","
    m.dict[pre * "image"] = p
    hdr_ps === nothing && (m.dict[pre * "pixsize"] = 1.0)
    m.dict[pre * "f"] = _share_flux!(m, pre * "f")
    console!(sh, "  added $(nm): image $(basename(p)), " * _image_summary(mi))
    hdr_ps === nothing &&
        console!(sh, "! $(basename(p)) has no pixel scale; $(pre)pixsize = 1.0 mas is a " *
                     "placeholder, set it in the image card"; kind = :err)
    return ""
end

"""
    shell_image_info(name) -> String

What the image card shows, as `path\tsummary\tpixsize\tpixsource\tnorm\tfilenorm`, or `!`
and a reason when the file cannot be read.

`pixsource` is `model` when `name,pixsize` overrides the header, `header` otherwise, and
`none` when the header has no scale (the value shown is then the model's placeholder).
`norm` is the model's `name,norm` or `""`; `filenorm` is what the file itself says
(`plane` or `relative`), which applies when `norm` is empty.
"""
function shell_image_info(name::AbstractString)
    m = _model(); cn = String(name)
    p = get(m.dict, cn * ",image", nothing)
    p isa AbstractString || return "! " * cn * " is not an image component"
    ov = get(m.dict, cn * ",pixsize", nothing)
    facts = image_facts(p)
    if facts === nothing
        return "! " * (isfile(p) ? (try read_model_image(p; pixsize = 1.0); "" catch err;
                                    _cause(err) end)
                                 : "no such file: " * p)
    end
    hdr_ps = try fits_pixsize(p) catch; nothing end
    ps, src = ov isa Real ? (Float64(ov), hdr_ps === nothing ? "none" : "model") :
                            (something(hdr_ps, NaN), hdr_ps === nothing ? "none" : "header")
    nm = get(m.dict, cn * ",norm", "")
    # Ten significant digits: a scale stored in degrees comes back as 0.05000000000000004 mas.
    return join((p, facts.summary, string(round(ps; sigdigits = 10)), src, string(nm),
                 string(facts.normalization)), '\t')
end

"""
    shell_set_image(name, field, value) -> String

Change one thing about an image component: `file` (a path), `pixsize` (mas, `""` to use the
header's) or `norm` (`plane`, `relative`, or `""` for what the file says). Returns `""` or `!`
and a reason; a file the reader refuses is not put in the model.
"""
function shell_set_image(name::AbstractString, field::AbstractString, value::AbstractString)
    sh = _shell(); m = _model(sh)
    cn = String(name); k = cn * ","
    haskey(m.dict, k * "image") || return "! " * cn * " is not an image component"
    v = strip(String(value))
    if field == "file"
        p = abspath(expanduser(v))
        ov = get(m.dict, k * "pixsize", nothing)
        mi = try
            read_model_image(p; pixsize = ov isa Real ? ov : nothing)
        catch err
            return "! " * _cause(err)
        end
        m.dict[k * "image"] = p
        console!(sh, "  $(k)image = $(basename(p)), " * _image_summary(mi))
    elseif field == "pixsize"
        if isempty(v)
            fits_pixsize(m.dict[k * "image"]) === nothing &&
                return "! the file has no pixel scale, so " * k * "pixsize is needed"
            delete!(m.dict, k * "pixsize")
            console!(sh, "  $(k)pixsize removed: the header's scale applies")
        else
            x = tryparse(Float64, v)
            (x === nothing || !(x > 0)) && return "! a pixel scale is a positive number of mas"
            m.dict[k * "pixsize"] = x
            console!(sh, "  $(k)pixsize = $(x) mas")
        end
    elseif field == "norm"
        if isempty(v)
            delete!(m.dict, k * "norm")
            console!(sh, "  $(k)norm removed: the file's normalisation applies")
        elseif v in ("plane", "relative")
            m.dict[k * "norm"] = String(v)
            console!(sh, "  $(k)norm = $(v)")
        else
            return "! norm is plane or relative"
        end
    else
        return "! unknown image field: " * String(field)
    end
    return ""
end

"""
    shell_expression_keywords(key) -> String

Every implicit variable, as `name\tshort\tenabled\treason\tmeaning`: a few words for beside
its button, whether an expression for `key` can use it, why not when it cannot, and the full
description.

They are SCOPED, and nothing else says so. Each kind of expression is compiled with its own
arguments -- a profile with the radial grid, a visibility function with each uv point's
baseline, wavelength and epoch, every other expression with the wavelength and epoch alone --
and `OITOOLS.implicit_vars` is that rule. Crossing it is not caught when it is typed: it surfaces later as
`UndefVarError: WL not defined in OITOOLS`, which names neither the rule nor the expression that
broke it, and which a plain typo produces too.

So every variable is always listed and the out-of-scope ones are disabled with the reason,
rather than being absent: a keyword that is simply missing reads as one that does not exist.
"""
function shell_expression_keywords(key::AbstractString)
    inscope = OITOOLS.implicit_vars(String(key))
    rows = (("R",    "radius, mas",            "radius over the profile grid, mas"),
            ("MU",   "μ, for limb darkening",  "sqrt(1 - (r/r_max)^2), for limb darkening"),
            ("D",    "diameter, 2R",           "diameter at each grid radius, i.e. 2\$R"),
            ("RMIN", "grid inner radius",      "the grid's inner radius, mas"),
            ("RMAX", "grid outer radius",      "the grid's outer radius, mas"),
            ("DMIN", "grid inner diameter",    "the grid's inner diameter, mas"),
            ("DMAX", "grid outer diameter",    "the grid's outer diameter, mas"),
            ("BL",   "baseline, Mλ",           "baseline length in Mλ, one per uv point"),
            ("WL",   "wavelength, µm",         "wavelength in MICRONS, one per uv point"),
            ("MJD",  "date, MJD",              "Modified Julian Date, one per uv point"))
    return join((string("\$", n, '\t', short, '\t', n in inscope ? "1" : "0", '\t',
                        n in inscope ? "" : _implicit_scope_reason(n), '\t', long)
                 for (n, short, long) in rows), "\n")
end

"""
    shell_expression_params(key) -> String

The model's parameters an expression for `key` can refer to, one `\$comp,name` per line: the
key's own component first, then the globals, then the other components. Written in full, which
resolves in every kind of expression. The names are the parser's (`all_names`), so a profile
string is not offered and a `spectrum` is offered as the `f` it is; `key` itself is left out.
Empty when the model does not compile, since the checker below the field says why.
"""
function shell_expression_params(key::AbstractString)
    m = _model(); k = String(key)
    isempty(m.dict) && return ""
    names = try
        dict_to_model(m.dict, String[]).all_names
    catch
        return ""
    end
    i = findfirst(',', k)
    comp = i === nothing ? "" : k[1:i-1]
    self = i !== nothing && k[i+1:end] == "spectrum" ? comp * ",f" : k
    rank(n) = (j = findfirst(',', n); j === nothing ? 1 : n[1:j-1] == comp ? 0 : 2)
    keep = sort([String(n) for n in names if n != self && n != k]; by = n -> (rank(n), n))
    return join(("\$" * n for n in keep), "\n")
end

# Where an implicit variable CAN be used, for the line explaining why it cannot be used here.
function _implicit_scope_reason(name)
    name == "BL" && return "only in a visibility function, which is the one place a baseline exists"
    name in ("WL", "MJD") &&
        return "not in a profile: the profile is evaluated on the radial grid, before any uv point"
    return "only inside a profile string, where compile_profile substitutes it"
end

"""
    shell_check_expression(key, expr) -> String

Every `\$` reference in `expr`, classified, as `ref\tclass\tmessage`.

`class` is `ok`, `scope` or `unknown`. This is the check the parser does not do when the
expression is written: `dict_to_model` accepts all three happily and the mistake only appears
when something evaluates, as an `UndefVarError` naming an internal module. An out-of-scope
implicit variable and a misspelt parameter are indistinguishable there, and distinguishing them
is the whole point of doing it here.

The lookup follows the parser, which differs by key: a `profile` or `visfunc` expression finds
a bare `\$a` on its own component first and then among the globals, while any other expression
looks a reference up literally -- `\$ud` in `star,f` is NOT `star,ud`, and says so.
"""
function shell_check_expression(key::AbstractString, expr::AbstractString)
    m = _model()
    k = String(key)
    i = findfirst(==(','), k)
    comp = i === nothing ? "" : k[1:i-1]
    par  = i === nothing ? k : k[i+1:end]
    qualifies = par == "profile" || par == "visfunc"
    inscope = OITOOLS.implicit_vars(k)
    out = String[]
    for ref in OITOOLS.extract_refs(String(expr))
        if ref in OITOOLS.IMPLICIT_VARS || ref == "BL"
            ok = ref in inscope
            push!(out, join((ref, ok ? "ok" : "scope",
                             ok ? "" : "\$$(ref): " * _implicit_scope_reason(ref)), '\t'))
        elseif haskey(m.dict, ref)
            push!(out, join((ref, "ok", ""), '\t'))
        elseif qualifies && !occursin(',', ref) && haskey(m.dict, comp * "," * ref)
            push!(out, join((ref, "ok", "resolves to " * comp * "," * ref), '\t'))
        elseif !qualifies && !isempty(comp) && !occursin(',', ref) &&
               haskey(m.dict, comp * "," * ref)
            push!(out, join((ref, "unknown",
                             "write \$" * comp * "," * ref * ": only a profile or a visibility " *
                             "function looks inside its own component"), '\t'))
        else
            push!(out, join((ref, "unknown",
                             par == "profile" ?
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

The radial grid a profile is evaluated on, as `key\tr_min\tr_max\tnr\tedge\tnodes\tfound`,
or `!` and a reason: `nr` is `c,nr` (0 when the count is chosen from the data and the
profile), `nodes` the count the model's quadrature uses on the preview's baselines, and
`found` what the probe saw in the profile (a jump, a kink, a singular centre), or empty.

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
    # Reversed diameters are the ring between them, as the model reads them.
    rmin > rmax && ((rmin, rmax) = (rmax, rmin))
    nr = try
        OITOOLS._hankel_Nr(cn, m.dict)
    catch
        0
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
    nodes, found = 0, ""
    expr = get(m.dict, cn * ",profile", nothing)
    if expr isa AbstractString && rmax > rmin
        try
            names = [r for r in OITOOLS.extract_refs(expr) if r ∉ OITOOLS.IMPLICIT_VARS]
            vals = Float64[]
            for n in names
                kk = haskey(m.dict, cn * "," * n) ? cn * "," * n : n
                push!(vals, Float64(OITOOLS._resolve_numeric(kk, m.dict)))
            end
            fn = OITOOLS.compile_profile(expr, String.(names); rmin = rmin, rmax = rmax)
            # The quadrature the model uses, on the preview's baselines: how many nodes, and
            # what the profile shows (jumps, kinks, a singular centre) that it adapted to.
            _, nodes, lay, _, _ = OITOOLS.hankel_profile_vis(fn, vals, rmin, rmax,
                                                             _profile_baselines(); nr)
            found = OITOOLS._hankel_describe(lay, rmax)
            found == "nothing unusual" && (found = "")
            rr = collect(range(rmin, rmax; length = PROFILE_NR_EDGE))
            mm = sqrt.(max.(0.0, 1 .- (rr ./ rmax) .^ 2))
            II = fn(rr, mm, vals...)
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
    return join((src, string(rmin), string(rmax), string(nr), string(edge), string(nodes),
                 found), '\t')
end

# The even grid the edge check and the drawn I(r) use: for the eye, not for the transform.
const PROFILE_NR_EDGE = 400

# The preview's baselines, Mλ, and the same in cycles/mas (one cycle per radian is
# 4.8481e-9 cycles per mas).
_profile_bml() = collect(range(0.0, PROFILE_B_MAX; length = PROFILE_NB))
_profile_baselines() = _profile_bml() .* 1e6 .* 4.84813681109536e-9

"""
    shell_profile_curves(component) -> String

Evaluate a profile and its Hankel transform onto the preview, returning `""` or a `!` message.

V(B) comes from `hankel_profile_vis`, the quadrature a model component uses, so the preview is
what a fit sees. It wants r in mas and B in cycles/mas, which is the same reciprocal pair, so
the only conversion here is from the baselines an observer reads -- Mλ -- into cycles/mas:
one cycle per radian is 4.8481e-9 cycles per mas. I(r) is drawn on an even grid, for the eye.

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
    fields = split(g, '\t')
    rmin = parse(Float64, fields[2]); rmax = parse(Float64, fields[3]); nr = parse(Int, fields[4])
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

    r  = collect(range(rmin, rmax; length = PROFILE_NR_EDGE))
    mu = sqrt.(max.(0.0, 1 .- (r ./ rmax) .^ 2))
    fn = try
        OITOOLS.compile_profile(expr, String.(names); rmin = rmin, rmax = rmax)
    catch err
        return "! the profile will not compile: " * _cause(err)
    end
    I = try
        fn(r, mu, vals...)
    catch err
        return "! the profile will not evaluate: " * _cause(err)
    end
    # A profile that does not mention `$R` or `$MU` compiles to `@. <constant>`, which has no
    # array operand and so returns a SCALAR. That is a legitimate profile -- a flat one is a
    # uniform disk, and it is what a new component starts as -- so it is spread over the grid
    # rather than rejected.
    I isa AbstractVector || (I = fill(Float64(I), length(r)))
    length(I) == length(r) ||
        return "! the profile returned $(length(I)) values for $(length(r)) grid points"
    # A centre where the profile is infinite (a power law) is integrable, and the quadrature
    # handles it; only the drawing leaves those points out.
    fin = isfinite.(I)
    any(fin) || return "! the profile is not finite anywhere on its grid"

    bml = _profile_bml()
    V = try
        # the model's own quadrature, so that the preview is what a fit sees
        first(OITOOLS.hankel_profile_vis(fn, vals, rmin, rmax, _profile_baselines(); nr))
    catch err
        return "! the transform failed: " * _cause(err)
    end
    update_profile_plot!(sh.profileplot, r[fin], I[fin], bml, V)
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
    # "=", not "==": `_op_string` canonicalises equality to a single "=", and that is what the
    # operator menu offers. Handing QML "==" left it with no matching entry, so the combo fell
    # back to its first item and the button silently produced `sum < 1`.
    return join((lhs, "=", "1", "0.001"), '\t')
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
twice it. An `r_max` of 0 suggests none: a profile written on the grid's own radii, as the
doughnut is with `\$RMIN`/`\$RMAX`, scales with whatever grid it is given.
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
     expr = "1 - (2*(\$R - (\$RMIN + \$RMAX)/2)/(\$RMAX - \$RMIN))^2",
     seeds = Pair{String,Float64}[], rmax = 0.0,
     note = "PMOIRED's doughnut: brightest mid-grid, zero at r_min and r_max"),

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
    orb = get(m.dict, cn * ",x", nothing) == "orbit"
    return string(pos ? "1" : "0", '\t', ori ? "1" : "0", '\t', orb ? "1" : "0")
end

# The elements a component put on an orbit starts with (PMOIRED's keys); MJD0 is the current
# dataset's mean date when there is one.
const _ORBIT_START = ["P" => 10.0, "MJD0" => 60000.0, "e" => 0.0, "omega" => 0.0,
                      "OMEGA" => 0.0, "incl" => 45.0, "a" => 5.0]

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
    if String(which) == "orbit"
        if Bool(on)
            mjd = try
                e = current_dataset(sh)
                e === nothing ? NaN : sum(Float64(d.mean_mjd) for d in e.data) / length(e.data)
            catch
                NaN
            end
            for (k, v) in _ORBIT_START
                key = "$cn,orb $k"
                haskey(m.dict, key) || (m.dict[key] = k == "MJD0" && isfinite(mjd) ? round(mjd; digits = 2) : v)
            end
            for c in ("x", "y")
                k = "$cn,$c"; m.dict[k] = "orbit"; filter!(!=(k), m.free)
                delete!(m.lb, k); delete!(m.ub, k)
            end
            console!(sh, "  $cn: on an orbit ($cn,x = $cn,y = \"orbit\"; elements $cn,orb P, MJD0, " *
                         "e, omega, OMEGA, incl, a)")
        else
            for k in [k for k in keys(m.dict) if startswith(k, "$cn,orb ")]
                delete!(m.dict, k); delete!(m.lb, k); delete!(m.ub, k); filter!(!=(k), m.free)
            end
            m.dict["$cn,x"] = 0.0; m.dict["$cn,y"] = 0.0
            console!(sh, "  $cn: off its orbit (x, y = 0)")
        end
        return ""
    end
    keys_ = String(which) == "position"    ? ("x", "y") :
            String(which) == "orientation" ? ("incl", "pa") :
            return "! unknown geometry: " * String(which)
    Bool(on) && keys_[1] == "incl" && haskey(m.dict, cn * ",image") &&
        return "! " * cn * " is an image, used at the orientation it was made with: it takes " *
               "no incl or pa"
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
    shell_ld_laws() -> String

The laws, as `key\tlabel\tnote` lines.
"""
shell_ld_laws() = join((string(l.key, '\t', l.label, '\t', l.note) for l in DISC_LAWS), "\n")

"""
    shell_ld_law(component) -> String

The law a component currently uses, as `key\tlabel`, or `""` when it is not a disc at all.
"""
function shell_ld_law(component::AbstractString)
    m = _model(); cn = String(component)
    i = findfirst(l -> haskey(m.dict, cn * "," * l.key), DISC_LAWS)
    i === nothing && return ""
    return string(DISC_LAWS[i].key, '\t', DISC_LAWS[i].label)
end

"""
    shell_set_ld_law(component, law) -> String

Change which law a disc uses (uniform or a limb-darkening law), keeping its diameter.

The diameter is the same physical quantity in every law and carries over, along with whether it
was being fitted. The coefficients do NOT: `u` is the linear coefficient in the quadratic and
the square-root laws alike, but the two laws are different shapes and a value fitted under one
is not a value under the other. They are reseeded, and the console says so, rather than being
carried across to look like a measurement that survived the change.
"""
function shell_set_ld_law(component::AbstractString, law::AbstractString)
    sh = _shell(); m = _model(sh)
    cn = String(component); want = String(law)
    j = findfirst(l -> l.key == want, DISC_LAWS)
    j === nothing && return "! unknown disc law: " * want
    i = findfirst(l -> haskey(m.dict, cn * "," * l.key), DISC_LAWS)
    i === nothing && return "! " * cn * " is not a disc"
    old, new = DISC_LAWS[i], DISC_LAWS[j]
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

# ── spectral lines (component card) ──────────────────────────────────────────────────

const _LINE_PROFILES = ("lorentzian", "gaussian", "truncexp", "trunclor")
const _C_KMS = 2.998e5          # PMOIRED's speed of light in km/s, as EX5 writes it

"""
    shell_lines(component) -> String

The component's spectral lines, one per line: `name<TAB>profile<TAB>velocity<TAB>rest`,
`velocity` the global its centre is written with (`""` when the centre is a value of its own),
`rest` the rest wavelength (µm) of a centre written by velocity, else the centre itself.
"""
function shell_lines(component::AbstractString)
    m = _model(); cn = String(component)
    out = String[]
    for n in OITOOLS._line_names(m.dict, cn)
        prof = something(findfirst(p -> haskey(m.dict, "$cn,line_$(n)_$p"), _LINE_PROFILES), 0)
        w0 = get(m.dict, "$cn,line_$(n)_wl0", nothing)
        v = w0 isa AbstractString ? (r = match(r"\$(v_[A-Za-z0-9_]+)", w0); r === nothing ? "" : r[1]) : ""
        rest = w0 isa Real ? string(w0) :
               w0 isa AbstractString && (r = match(r"^\s*([0-9.eE+-]+)\s*\*", w0)) !== nothing ? String(r[1]) : ""
        push!(out, join((n, prof == 0 ? "" : _LINE_PROFILES[prof], v, rest), '\t'))
    end
    return join(out, "\n")
end

"""
    shell_add_line(component, name, profile, wl0, f, width) -> String

Add a spectral line to a component: `c,line_<name>_wl0` (µm), `_f` (peak, < 0 absorption) and
`_<profile>` (width in nm), PMOIRED's keys, as fixed rows of the parameter table.
"""
function shell_add_line(component::AbstractString, name::AbstractString, profile::AbstractString,
                        wl0::Real, f::Real, width::Real)
    sh = _shell(); m = _model(sh); cn = String(component); n = strip(String(name))
    any(startswith(k, cn * ",") for k in keys(m.dict)) || return "! no component called " * cn
    occursin(r"^[A-Za-z0-9]+$", n) || return "! a line's name is letters and digits"
    n in OITOOLS._line_names(m.dict, cn) && return "! $cn already has a line $n"
    p = String(profile)
    p in _LINE_PROFILES || return "! unknown profile " * p
    (isfinite(wl0) && wl0 > 0) || return "! the centre is a wavelength in µm"
    (isfinite(width) && width != 0) || return "! the width is in nm, not 0"
    m.dict["$cn,line_$(n)_wl0"] = Float64(wl0)
    m.dict["$cn,line_$(n)_f"]   = Float64(f)
    m.dict["$cn,line_$(n)_$p"]  = Float64(width)
    console!(sh, "  $cn: line $n ($p) at $(wl0) µm, peak $(f), width $(width) nm")
    return ""
end

"""
    shell_remove_line(component, name) -> String

Remove a spectral line, every key of it, with the velocity global it was written with.
"""
function shell_remove_line(component::AbstractString, name::AbstractString)
    sh = _shell(); m = _model(sh); cn = String(component); n = String(name)
    pre = "$cn,line_$(n)_"
    ks = [k for k in keys(m.dict) if startswith(k, pre) &&
          occursin(OITOOLS._LINE_KEY, k[(length(cn) + 2):end])]
    isempty(ks) && return "! $cn has no line $n"
    v = _line_velocity_global(m, cn, n)
    for k in [ks; isempty(v) ? String[] : [v]]
        delete!(m.dict, k); delete!(m.lb, k); delete!(m.ub, k); filter!(!=(k), m.free)
    end
    console!(sh, "  $cn: removed line $n")
    return ""
end

_line_velocity_global(m, cn, n) =
    (w0 = get(m.dict, "$cn,line_$(n)_wl0", nothing);
     w0 isa AbstractString && (r = match(r"\$(v_[A-Za-z0-9_]+)", w0)) !== nothing ? String(r[1]) : "")

"""
    shell_line_by_velocity(component, name, on, rest) -> String

Write a line's centre as a velocity (`rest·(1 + \$v_<c>_<name>/c)`, `rest` in µm, a new global
`v_<c>_<name>` in km/s set to the centre's current shift — EX5's reparametrisation), or back as
a value of its own (`on = false`, the centre the velocity gives).
"""
function shell_line_by_velocity(component::AbstractString, name::AbstractString, on::Bool, rest::Real)
    sh = _shell(); m = _model(sh); cn = String(component); n = String(name)
    k = "$cn,line_$(n)_wl0"
    haskey(m.dict, k) || return "! $cn has no line $n"
    v = "v_$(cn)_$(n)"
    if on
        w0 = m.dict[k]
        w0 isa Real || return "! the centre is already an expression"
        (isfinite(rest) && rest > 0) || return "! the rest wavelength is in µm"
        haskey(m.dict, v) && return "! a parameter $v exists already"
        m.dict[v] = (w0 / rest - 1) * _C_KMS
        m.dict[k] = "$(rest) * (1 + \$$(v) / $(_C_KMS))"
        console!(sh, "  $cn: line $n's centre written as $(rest) µm at velocity $v")
    else
        g = _line_velocity_global(m, cn, n)
        isempty(g) && return "! the centre is not written with a velocity"
        r = match(r"^\s*([0-9.eE+-]+)\s*\*", m.dict[k])
        r === nothing && return "! the centre is not of the form rest * (1 + \$v/c)"
        rest0 = parse(Float64, r[1])
        m.dict[k] = rest0 * (1 + Float64(m.dict[g]) / _C_KMS)
        for kk in (g,); delete!(m.dict, kk); delete!(m.lb, kk); delete!(m.ub, kk); filter!(!=(kk), m.free); end
        console!(sh, "  $cn: line $n's centre back to a value, $(m.dict[k]) µm")
    end
    return ""
end

"""
    shell_remove_component(name) -> String

Delete a component and every key belonging to it, including its bounds and its place in the
free list.

Expressions elsewhere that referenced it are deliberately NOT rewritten: a dangling `\$name,ud`
is listed by `model_warnings`, which the panel shows, where a silent substitution would be an
invisible one.
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
    shell_remove_param(key) -> String

Delete one parameter that is the model's own rather than the parser's: a global, or a key on a
component whose name the parser gives no meaning to (`disk,slope`, added with
[`shell_add_param`](@ref)). Its bounds and its place in the free list go with it.

A component's OTHER keys are what the parser identifies it by, so deleting `disk,fwhm` on its
own does not leave a smaller disc, it leaves a `disk,f` the parser can no longer place;
removing the component is the operation that exists for that. A global belongs to no
component, so [`shell_remove_component`](@ref) can never reach one, and without this a global
outlives every component that referenced it and cannot be cleared at all.

Expressions that still refer to it are not rewritten, for the same reason a removed component's
references are not. Unlike that case the result IS reported: `model_warnings` lists every
expression naming a parameter the model does not define, which is what the panel shows.
"""
function shell_remove_param(key)
    sh = _shell(); m = _model(sh); k = String(key)
    haskey(m.dict, k) || return "! no such parameter: " * k
    i = findfirst(',', k)
    i !== nothing && OITOOLS.reserved_suffix(k[i+1:end]) &&
        return "! " * k * " is part of what makes the component what it is; remove the " *
               "component instead"
    delete!(m.dict, k); delete!(m.lb, k); delete!(m.ub, k)
    filter!(!=(k), m.free)
    nref = count(v -> v isa AbstractString && k in OITOOLS.extract_refs(v), values(m.dict))
    console!(sh, "  removed " * (i === nothing ? "the global " : "") * k *
                 (nref > 0 ? ": $(nref) expression(s) still refer to it" : ""))
    return ""
end

# Why `name` cannot be a new parameter on component `cn` ("" for a global), or "" when it can.
function _param_name_problem(m, cn::AbstractString, nm::AbstractString)
    isempty(nm) && return "a parameter needs a name"
    occursin(r"^[A-Za-z_][A-Za-z0-9_]*$", nm) ||
        return "a name is letters, digits and underscores, not starting with a digit: " *
               "it has to survive being written as \$" * nm
    (nm in OITOOLS.IMPLICIT_VARS || nm == "BL") &&
        return "\$" * nm * " is a built-in variable"
    if isempty(cn)
        nm in ("spatial_kernel", "smear", "smear_mode") &&
            return nm * " is a setting the parser reads by name"
        haskey(m.dict, nm) && return "there is already a global called " * nm
    else
        any(k -> startswith(k, cn * ","), keys(m.dict)) || return "no component called " * cn
        OITOOLS.reserved_suffix(nm) &&
            return "on a component, " * nm * " already means something to the parser"
        haskey(m.dict, cn * "," * nm) && return "there is already a parameter " * cn * "," * nm
    end
    return ""
end

"""
    shell_remove_unused() -> String

Delete every key the parser ignores -- the ones Validation reports as not used or ignored --
with their bounds and their places in the free list. Returns `""` or a message beginning `!`.

These are keys `model_inspection` reports as unrecognised: on a component, a name no kind
reads and no expression refers to; it is either a parameter added and never used, or a
misspelt shape key. Either way the model is the same without it. The console names each key
removed, so a misspelling deleted here can be typed back correctly.
"""
function shell_remove_unused()
    sh = _shell(); m = _model(sh)
    insp = try
        model_inspection(m.dict)
    catch err
        return "! " * _cause(err)
    end
    ks = [String(k) for k in insp.unrecognised if haskey(m.dict, k)]
    isempty(ks) && return "! nothing to remove: the parser uses every key"
    for k in ks
        delete!(m.dict, k); delete!(m.lb, k); delete!(m.ub, k)
        filter!(!=(k), m.free)
    end
    console!(sh, "  removed what nothing used: " * join(ks, ", "))
    return ""
end

"""
    shell_param_removable(key) -> String

`"1"` when [`shell_remove_param`](@ref) would delete `key`, otherwise `"0"`: the panel enables
its Remove button from this, so the button and the refusal cannot disagree.
"""
function shell_param_removable(key)
    k = String(key)
    haskey(_model().dict, k) || return "0"
    i = findfirst(',', k)
    return i === nothing || !OITOOLS.reserved_suffix(k[i+1:end]) ? "1" : "0"
end

"""
    shell_param_name_problem(component, name) -> String

Why `name` cannot be added to `component` (`""` or `"__global__"` for a global), or `""`. The
dialog asks as the name is typed, so OK is never pressed on a name that would be refused.
"""
function shell_param_name_problem(component, name)
    cn = String(component); cn == GLOBAL_COMPONENT && (cn = "")
    return _param_name_problem(_model(), cn, strip(String(name)))
end

"""
    shell_add_param(component, name, value) -> String

Create a fixed parameter at `value`: a global when `component` is `""` or `"__global__"`,
otherwise `component,name`. Returns `""` or a message beginning with `!`.

This is how a model gets a quantity that no component kind seeds: SPARCO's flux fraction `fs`,
a spectral slope, the `d0` and `p` of a diameter written as `\$d0 * (\$WL/1.6)^\$p`. It
starts fixed and referenced by nothing, because what it means is decided by the expression
that uses it; until one does, the parser lists it among the keys it ignored, which is true.

A name the parser already reads on a component is refused rather than added: `ud` on a
Gaussian is not a new parameter, it turns the Gaussian into a uniform disc.
"""
function shell_add_param(component, name, value)
    sh = _shell(); m = _model(sh)
    cn = String(component); cn == GLOBAL_COMPONENT && (cn = "")
    nm = strip(String(name))
    msg = _param_name_problem(m, cn, nm)
    isempty(msg) || return "! " * msg
    key = isempty(cn) ? String(nm) : cn * "," * nm
    m.dict[key] = Float64(value)
    console!(sh, "  added " * (isempty(cn) ? "the global " : "") * key * " = " *
                 string(Float64(value)) * "; nothing refers to it yet -- write \$" * key *
                 " in an expression")
    return ""
end

"""
    shell_sparco_fluxes(star, env, index, lambda0) -> String

SPARCO's chromatic fluxes on two components of the model: `star` follows a power law of fixed
`index` (−4, Rayleigh–Jeans, for a photosphere in the near infrared) and `env` (a disc, a ring,
an image…) its own slope, the global `fs` being the star's share of the flux at `lambda0`
(µm). Writes

    star,f = \$fs·(\$WL/λ0)^index / N        env,f = (1 − \$fs)·(\$WL/λ0)^\$env,slope / N
    N      = \$fs·(\$WL/λ0)^index + (1 − \$fs)·(\$WL/λ0)^\$env,slope

so the two sum to 1 at every wavelength, and frees `fs` in [0, 1] and `env,slope` in [−4, 4]
(created at 0.5 and 0 when absent; kept when present). Other components keep their own fluxes.
Returns `""` or a message beginning with `!`.
"""
function shell_sparco_fluxes(star, env, index::Real, lambda0::Real)
    sh = _shell(); m = _model(sh)
    st, en = String(star), String(env)
    comps = OITOOLS._component_names(m.dict)
    (st in comps && en in comps) || return "! choose two components of the model"
    st == en && return "! the star and the environment must be two components"
    isfinite(index) || return "! the star index must be a number"
    lambda0 > 0 || return "! λ0 must be a positive wavelength in µm"
    slope = en * ",slope"
    if !haskey(m.dict, "fs")
        msg = _param_name_problem(m, "", "fs")
        isempty(msg) || return "! fs: " * msg
        m.dict["fs"] = 0.5
    end
    if !haskey(m.dict, slope)
        msg = _param_name_problem(m, en, "slope")
        isempty(msg) || return "! $(slope): " * msg
        m.dict[slope] = 0.0
    end
    l0, k = @sprintf("%g", lambda0), @sprintf("%g", index)
    sp = "\$fs * (\$WL/$l0)^($k)"
    ep = "(1 - \$fs) * (\$WL/$l0)^\$$slope"
    norm = "($sp + $ep)"
    for (key, num) in ((st * ",f", sp), (en * ",f", ep))
        m.dict[key] = "$num / $norm"
        filter!(!=(key), m.free)
    end
    for (key, lo, hi) in (("fs", 0.0, 1.0), (slope, -4.0, 4.0))
        m.dict[key] isa Real || (m.dict[key] = key == "fs" ? 0.5 : 0.0)   # free means numeric
        m.lb[key] = lo; m.ub[key] = hi
        key in m.free || push!(m.free, key)
    end
    others = setdiff(comps, [st, en])
    console!(sh, "  SPARCO fluxes: $(st),f and $(en),f sum to 1 at every wavelength; fs and $(slope) free" *
                 (isempty(others) ? "" : " ($(join(others, ", ")) keep their own fluxes)"))
    return ""
end

"""
    shell_link_params(name, keys) -> String

Tie parameters together: create the global `name` and make each of `keys` (tab-separated) the
expression `\$name`. Returns `""` or a message beginning with `!`.

The global takes the first parameter's value and bounds, and is free when any of them was, so
linking two fitted position angles leaves one fitted position angle. The linked rows lose their
bounds and their place in the fit vector, as any derived row does.

Through a global rather than by making one row follow another, because a link should survive
either side being removed: `disk,pa = \$ring,pa` dangles when the ring goes, `\$pa` does not.
"""
function shell_link_params(name, keys_)
    sh = _shell(); m = _model(sh)
    nm = String(strip(String(name)))
    ks = [String(k) for k in split(String(keys_), '\t') if !isempty(k)]
    length(unique(ks)) < 2 && return "! linking needs at least two parameters"
    for k in ks
        haskey(m.dict, k) || return "! no such parameter: " * k
        m.dict[k] isa Real || return "! " * k * " is already an expression"
    end
    msg = _param_name_problem(m, "", nm)
    isempty(msg) || return "! " * msg

    k1 = ks[1]
    lb, ub = get(m.lb, k1, nothing), get(m.ub, k1, nothing)
    if lb === nothing || ub === nothing
        dlb, dub = default_bounds(m.dict, [k1])
        lb === nothing && (lb = get(dlb, k1, nothing))
        ub === nothing && (ub = get(dub, k1, nothing))
    end
    vals = [Float64(m.dict[k]) for k in ks]
    wasfree = any(k -> k in m.free, ks)

    m.dict[nm] = vals[1]
    lb === nothing || (m.lb[nm] = Float64(lb))
    ub === nothing || (m.ub[nm] = Float64(ub))
    for k in ks
        m.dict[k] = "\$" * nm
        delete!(m.lb, k); delete!(m.ub, k)
        filter!(!=(k), m.free)
    end
    wasfree && push!(m.free, nm)

    console!(sh, "  linked " * join(ks, ", ") * " through the global " * nm *
                 (wasfree ? " (free)" : " (fixed)"))
    allequal(vals) ||
        console!(sh, "  they held " * join(string.(vals), ", ") * "; " * nm * " starts at " *
                     string(vals[1]) * ", the first")
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
    # Collapsed on the live path too, or a progress bar is a horizontal rule for the whole run
    # and only tidies itself up at the end.
    istaskdone(j.task) || return "running\t" * _collapse_cr(job_output(sh))

    stopped = j.stop[]
    kind = j.kind
    res  = job_finish!(sh)
    if stopped
        # Abandoned rather than aborted; see `shell_job_stop`. Whatever came back is discarded,
        # because the user asked to stop looking at it.
        line = "stopped"
        kind === :image ? (sh.enginelog = _collapse_cr(res.output)) :
                          (sh.fitlog   = _collapse_cr(res.output))
        console!(sh, "  " * line)
    else
        line = kind === :image  ? finish_reconstruct!(sh, res) :
               kind === :epochs ? finish_epochs!(sh, res)      : finish_fit!(sh, res)
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
    kind      :: Symbol           # :image, :fit or :epochs — decides what `job_finish!` does with it
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

Return the imaging display to the whole image, framed exactly as `show_image!` framed it.

Separate from [`shell_reset_zoom`](@ref), which resets the Exploring canvas: the two are
different axes on different figures, and a right click on one must not move the other.

`reset_view!`, not `autolimits!`. The right click on this canvas arrives through the QML
`MouseArea` laid over it, so this function decides the view the click returns to, and
`autolimits!` pads the heatmap's cell edges by the axis's 5% autolimit margin: 1.117x the field
of view at 64 pixels, measured, which reads as the reset not quite working.
"""
function shell_reset_image_zoom()
    sh = _shell()
    cv = sh.imcanvas
    cv === nothing && return ""
    try
        reset_view!(cv)
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
    cube = try
        ft, _, _ = ensure_ft!(sh.ftcache, e.data, setup)
        Float64.(start_image(setup, ft))
    catch err
        msg = "! could not build the starting image: " * _cause(err)
        console!(sh, msg; kind = :err); return msg
    end
    # A start is a cube, one plane per channel, and the canvas draws a plane: the channel the
    # slider is on, clamped, since the slider may still be set for a result with more channels.
    nw = _nchannels(cube)
    w  = clamp(sh.imchannel, 1, nw)
    sh.imcanvas === nothing || show_image!(sh.imcanvas, _plane(cube, w), setup.pixsize;
                                           label = nw > 1 ? "starting image, channel $w of $nw" :
                                                            "starting image",
                                           title = image_title(setup; kind = "start image",
                                                               channel = w, nchannels = nw))
    console!(sh, "> start_image(setup, ft)   # $(startkind), $(Int(nx))×$(Int(nx))")
    return "showing the starting image"
end

"""
    shell_vi_prior_sample(nx, pixsize, options) -> String

Draw one image from the VI sky prior and put it on the imaging canvas, without running
anything. The Show-start-image of variational inference.

A prior is the whole of what VI adds to the data, and until it can be looked at it is a set of
numbers in boxes. This draws `xi ~ N(0, I)` through the same `SkyModelParams` a run would
build — `vi_sky_params`, one definition for both — so what appears is the prior the next Run
will use, support map included.

Each call is a new draw. What matters about a prior is the range of images it allows, and a
button that returned the same picture every time would hide exactly that.
"""
function shell_vi_prior_sample(nx::Integer, pixsize::Real, options::AbstractString = "")
    sh = _shell()
    e = current_dataset(sh)
    e === nothing && return "! no dataset loaded"
    ext = Base.get_extension(OITOOLS, :OITOOLSVarInfExt)
    ext === nothing && return "! variational inference needs VarInf; install it and restart"
    o = _parse_options(options)
    img = try
        p = vi_sky_params(ext, o, Int(nx), Float64(pixsize), _mean_lam(e.data[1, 1]))
        vi_prior_draw(ext, p)
    catch err
        msg = "! could not draw from the prior: " * _cause(err)
        console!(sh, msg; kind = :err); return msg
    end
    sh.imcanvas === nothing || show_image!(sh.imcanvas, img, Float64(pixsize);
                                           label = "prior sample",
                                           title = image_title(
                                               ImagingSetup(; nx = Int(nx),
                                                              pixsize = Float64(pixsize));
                                               kind = "prior sample"))
    support = haskey(o, "vi_prior") && !isempty(o["vi_prior"]) ?
              "weight = " * basename(o["vi_prior"]) :
              "disc R = " * string(round(_optreal(o, "vi_radius", nx * pixsize / 5);
                                         digits = 3)) * " mas"
    console!(sh, "sky_forward(randn(_latent_size(p)), p)   # $(support)"; kind = :cmd)
    console!(sh, @sprintf("  prior sample: flux %.4g, peak %.4g", sum(img), maximum(img)))
    return "showing a sample from the sky prior"
end

"""
    shell_set_overlay_mode(mode) -> String

Draw the model's or the reconstruction's observables ON TOP of the data: `""`, `"model"` or
`"image"`.

The companion of [`shell_set_residual_mode`](@ref) and the other half of the same question: an
overplot shows whether the prediction goes through the points, a residual shows what is left
when it does not. Mutually exclusive with the residual view, because a residual already has the
model subtracted — overplotting it there would draw the prediction against its own difference.
"""
function shell_set_overlay_mode(mode::AbstractString)
    sh = _shell()
    m = mode == "model" ? :model : mode == "image" ? :image : :none
    sh.overlay = m
    m === :none || (sh.residual = :none)
    sh.status = refresh_plot!(sh)
    m === :none || console!(sh, m === :model ? "model_to_obs(model, x, data)" :
                                               "image_to_obs(x, ft, data)"; kind = :cmd)
    return sh.status
end

"""
    shell_set_residual_mode(mode) -> String

Switch the plot between the observable and its residual: `""`, `"model"` or `"image"`.

A modifier on whatever the kind menu has selected, so `v2` becomes `v2` residuals and `t3phi`
becomes `t3phi` residuals, each against its own baseline and in its own groups. That is the
point of doing it here rather than in a panel of its own: the residual is read against the same
axes as the quantity it came from, and every view that works for the data — colour-by, per
group, picking — works for it unchanged.
"""
function shell_set_residual_mode(mode::AbstractString)
    sh = _shell()
    m = mode == "model" ? :model : mode == "image" ? :image : :none
    sh.residual = m
    m === :none || (sh.overlay = :none)   # a residual has the model subtracted already
    sh.status = refresh_plot!(sh)
    m === :none || console!(sh, m === :model ? "model_to_residuals(model, x, data)" :
                                               "image_to_residuals(x, ft, data)"; kind = :cmd)
    return sh.status
end

"""
    shell_sparco_converged() -> String

The SPARCO parameters the last run converged to, as `key<TAB>value` lines keyed by the panel's
own names, or `""` when the last run was not a SPARCO one.

The two SPARCO engines report differently and neither reports in the panel's vocabulary, so the
translation lives here rather than in QML:

  * the hybrid returns `param_names`, a Dict over its own free-parameter spellings
    (`star,f0`, `bg,f0`, `d_env`, `star,di`);
  * the annealing sampler returns `params`, one vector PER CHAIN, against the fixed tuple
    `SPARCO_PARAM_NAMES` -- so the values worth reporting are the best chain's, not chain 1's.

`star,di` is deliberately NOT called `ud`: it is the star's spectral index in
`(lambda/lambda0)^(-star,di)`, whereas the annealing engine's `ud` is a uniform-disc diameter
in mas. Two quantities, two keys.
"""
function shell_sparco_converged()
    sh = _shell()
    r = sh.imaging
    r === nothing && return ""
    e = r.extra
    e === nothing && return ""
    out = String[]

    # The hybrid: a Dict of its own free-parameter names.
    if e isa NamedTuple && haskey(e, :param_names) && e.param_names isa AbstractDict
        m = Dict("star,f0" => "f_star", "bg,f0" => "f_bg",
                 "d_env" => "env_indx", "star,di" => "star_di")
        for (k, v) in e.param_names
            gk = get(m, String(k), nothing)
            gk === nothing && continue
            push!(out, string(gk, '\t', Float64(v)))
        end

    # The annealing sampler: per-chain vectors against a fixed name tuple.
    elseif e isa NamedTuple && haskey(e, :params) && haskey(e, :param_names) &&
           e.params isa AbstractVector && !isempty(e.params)
        best = haskey(e, :best_chain) ? Int(e.best_chain) : 1
        v = e.params[clamp(best, 1, length(e.params))]
        v isa AbstractVector || return ""
        names = String.(collect(e.param_names))
        for (i, nm) in enumerate(names)
            i <= length(v) || break
            # lambda0 is never sampled; echoing it back as "converged" would suggest it was.
            nm == "lambda0" && continue
            push!(out, string(nm, '\t', Float64(v[i])))
        end
    end
    return join(sort(out), "\n")
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
    p = strip_file_url(path)
    try
        # The cube and its wavelengths together: a saved cube whose third axis carries no WCS
        # is a set of channels nothing can identify, and the result is the only thing that
        # knows which bin each plane came from.
        writefits(r.image, p; pixsize = r.setup.pixsize, wavelengths = r.wavelengths)
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
    # Channel by channel, with the SAME shift for all of them: `recenter` finds its own
    # centroid, so recentring each plane independently would align the channels to each other
    # and destroy exactly the differential astrometry a cube is reconstructed to show.
    img = _recenter_cube(r.image, sh.imchannel)
    # ImagingResult is immutable, so rebuild it. Every other field survives the shift: the
    # χ², the point count and the timings all describe the run that produced the pixels, and
    # moving them does not change what that run did.
    #
    # Rebuilt with EVERY field, not through the short constructor: that one fills `extra` and
    # `ensemble` with `nothing`, so recentring a SQUEEZE or VI result silently threw away its
    # diagnostics and its posterior — the run's most expensive output, discarded by a shift.
    sh.imaging = ImagingResult(img, r.chi2, r.chi2_start, r.ndof, r.flux, r.maxiter,
                               r.seconds, r.setup, r.weights, r.breakdown, r.extra,
                               r.ensemble, r.wavelengths)
    sh.imcanvas === nothing ||
        show_image!(sh.imcanvas, result_plane(sh.imaging, sh.imchannel), r.setup.pixsize;
                    title = _result_title(sh, sh.imaging))
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
    shell_central_telescope(facility, telescopes) -> String

The name of the delay reference to use among `telescopes` (space-separated names):
[`central_telescope`](@ref)'s choice. `""` with none.
"""
function shell_central_telescope(facility::AbstractString, telescopes::AbstractString)
    f = read_facility_file(String(facility))
    sel = filter(!isempty, String.(split(String(telescopes))))
    r = central_telescope(f, telescope_config(sel, String.(f.tel_names)))
    return r == 0 ? "" : String(f.tel_names[r])
end

"""
    shell_telescope_order(name) -> String

The facility's telescopes in the order the Observing panel lists them, one per line: the
facility file's `display_order` when it names each telescope once, else the file's own order.
"""
function shell_telescope_order(name::AbstractString)
    f = read_facility_file(String(name))
    return join(String.(f.tel_names)[telescope_display_order(f)], "\n")
end

# A clock name from the panel: `:lst`, `:utc` or `:local`, LST for anything else.
function _time_system(name::AbstractString)
    t = Symbol(lowercase(String(name)))
    return t in (:lst, :utc, :local) ? t : :lst
end

# The part of the night the panel uses, from its `"from-to"` UT text (`"3:30-8"`, `"3:30-"`):
# an empty side is the start or the end of the night (NaN), an empty text the whole night.
function _ut_window(night::AbstractString)
    t = strip(String(night))
    (isempty(t) || t == "-") && return nothing
    hours(x) = isempty(strip(x)) ? NaN :
               (f = split(strip(x), ':'); parse(Int, f[1]) + (length(f) > 1 ? parse(Int, f[2]) / 60 : 0.0))
    a, b = split(t, '-')
    return (hours(a), hours(b))
end

# The night the Observe panel draws, from its controls: `(; plan, date, config)`, or a message
# beginning with `!`. `plan_text` is a multi-POP plan as `pop_plan_text`; empty, the POPs `pops`.
function _night_plan_for(facility, name, ra, dec, dateiso, telescopes, pops, use_delay,
                         alt_limit, alt_max, plan_text, change_minutes; ut_window = nothing)
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
    return try
        plan = isempty(strip(String(plan_text))) ? nothing : parse_pop_plan_text(plan_text)
        (; plan = night_plan(String(facility), String(name), Float64(ra), Float64(dec), date;
                             config = cfg, pop = pp, use_delay = use_delay, popplan = plan,
                             change_minutes = Float64(change_minutes),
                             alt_limit = Float64(alt_limit), alt_max = Float64(alt_max), ut_window),
           date, config = cfg)
    catch err
        "! " * _cause(err)
    end
end

"""
    shell_gantt(facility, name, ra, dec, dateiso, telescopes, pops, use_delay, detailed,
                alt_limit, alt_max, time_system, plan_text, change_minutes, stack_mode, cal,
                night = "") -> String

Compute one night and draw it on the Observe canvas. Returns a one-line summary.

`ra` and `dec` are DEGREES, `telescopes` and `pops` are space-separated. `use_delay` is opt-in
because an unsearched POP configuration reports far less time than is really available — see
[`night_plan`](@ref). `time_system` is `"lst"`, `"utc"` or `"local"` and labels the x axis;
the bars do not move, only the ticks. `cal` draws a calibrator's bars in the calibrator colours.
`night` is the part of the night used, `"from-to"` in UT hours (`"3:30-8"`); empty, all of it.
"""
function shell_gantt(facility::AbstractString, name::AbstractString, ra::Real, dec::Real,
                     dateiso::AbstractString, telescopes::AbstractString,
                     pops::AbstractString, use_delay::Bool, detailed::Bool,
                     alt_limit::Real = DEFAULT_ALT_LIMIT, alt_max::Real = DEFAULT_ALT_MAX,
                     time_system::AbstractString = "lst", plan_text::AbstractString = "",
                     change_minutes::Real = 10, stack_mode::AbstractString = "",
                     cal::Bool = false, night::AbstractString = "")
    sh = _shell()
    uw = try
        _ut_window(night)
    catch
        return "! night: not from-to UT hours: $(night)"
    end
    np = _night_plan_for(facility, name, ra, dec, dateiso, telescopes, pops, use_delay,
                         alt_limit, alt_max, plan_text, change_minutes; ut_window = uw)
    np isa String && (startswith(np, "! not a date") || console!(sh, np); return np)
    p, date, cfg = np.plan, np.date, np.config

    tsys = Symbol(lowercase(String(time_system)))
    tsys in (:lst, :utc, :local) || (tsys = :lst)
    GANTT_TIME_SYSTEM[] = tsys
    GANTT_CAL[] = cal
    # With `stack_mode` "coverage", the single POP settings best over the earliest, middle and
    # latest third of the observable hours, one row each, in the summary view. `p` becomes the
    # first, which the delay chart shows; the detailed view shows `p` as given.
    GANTT_STACK[] = nothing
    mode = String(stack_mode)
    if mode == "coverage" && use_delay && !detailed && isempty(strip(String(plan_text)))
        stack, labels = try
            f = read_facility_file(String(facility))
            ha, ok = _night_sky(f, ra, dec, date, cfg; alt_limit, alt_max, ut_window = uw)
            res = any(ok) ? coverage_pops(f, Float64(dec), ha[ok], cfg) : []
            console!(sh, "> coverage_pops(facility, $(round(Float64(dec); digits=4)), ha, $(cfg))")
            ([night_plan(String(facility), String(name), Float64(ra), Float64(dec), date;
                         config = cfg, pop = Int.(r.pop), use_delay = true,
                         alt_limit = Float64(alt_limit), alt_max = Float64(alt_max),
                         ut_window = uw) for r in res],
             [string(name, " ", r.part) for r in res])
        catch err
            console!(sh, "! POP rows: " * _cause(err)); (NightPlan[], String[])
        end
        if !isempty(stack)
            p = first(stack)
            GANTT_STACK[] = (; kind = :coverage, plans = stack, labels,
                               title = "best time coverage POPs")
        end
    end
    if sh.gantt !== nothing
        GANTT_STACK[] !== nothing ?
            draw_gantt_stack!(sh.gantt, GANTT_STACK[]; time_system = tsys, cal) :
            update_gantt!(sh.gantt, p; detailed, time_system = tsys, cal)
    end
    # The delay chart is drawn from the same plan, so both views are always of one night.
    sh.delayplot === nothing || isempty(p.baselines) || update_delay_plot!(sh.delayplot, p; time_system = tsys)

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

"""
    shell_pop_timeline(facility, name, ra, dec, dateiso, telescopes, pops, use_delay,
                       alt_limit, alt_max, time_system, plan_text, change_minutes,
                       as_text = false, night = "") -> String

The night `shell_gantt` draws for the same arguments, as its POP timeline
([`pop_plan_timeline`](@ref)): one line per row, `from \t to \t kind \t pops \t ha \t time \t
note`, times `H:MM` on the clock `time_system`, `kind` one of `setting change drop readd gap`.
With `as_text`, the table as [`timeline_text`](@ref) writes it. `""` for a night with no
darkness; a message beginning with `!` on failure.
"""
function shell_pop_timeline(facility::AbstractString, name::AbstractString, ra::Real, dec::Real,
                            dateiso::AbstractString, telescopes::AbstractString,
                            pops::AbstractString, use_delay::Bool,
                            alt_limit::Real = DEFAULT_ALT_LIMIT, alt_max::Real = DEFAULT_ALT_MAX,
                            time_system::AbstractString = "lst", plan_text::AbstractString = "",
                            change_minutes::Real = 10, as_text::Bool = false,
                            night::AbstractString = "")
    uw = try
        _ut_window(night)
    catch
        return "! night: not from-to UT hours: $(night)"
    end
    np = _night_plan_for(facility, name, ra, dec, dateiso, telescopes, pops, use_delay,
                         alt_limit, alt_max, plan_text, change_minutes; ut_window = uw)
    np isa String && return np
    tsys = _time_system(time_system)
    rows = pop_plan_timeline(np.plan; system = tsys)
    as_text && return timeline_text(rows, np.plan; system = tsys)
    return join((join(timeline_fields(r, np.plan), "\t") for r in rows), "\n")
end

# A list of targets from the panel, `name \t ra \t dec \t cal \t for` per line, on one night:
# `(; f, date, cfg, uw, targets, cals, obs, lst, sky)` with each target's
# `night_observability` and where it is observable (`sky`, samples × targets), or a message
# beginning with `!`. `for` names the science targets a calibrator is for, `;`-separated; with
# `needs_cal` false no calibrator is tied to a target.
function _targets_night(facility, rows, dateiso, telescopes, alt_limit, alt_max, night, needs_cal)
    date = try
        DateTime(String(dateiso))
    catch
        return "! not a date: $(dateiso)"
    end
    uw = try
        _ut_window(night)
    catch
        return "! night: not from-to UT hours: $(night)"
    end
    f = read_facility_file(String(facility))
    tels = facility_telescopes(String(facility))
    sel = filter(!isempty, String.(split(String(telescopes))))
    cfg = isempty(sel) ? ones(Int, length(tels)) : telescope_config(sel, tels)
    targets = NamedTuple{(:name, :ra, :dec, :cal),Tuple{String,Float64,Float64,Bool}}[]
    fors = Vector{String}[]
    for line in split(String(rows), '\n'; keepempty = false)
        x = split(line, '\t')
        length(x) >= 3 || continue
        ra, dec = tryparse(Float64, x[2]), tryparse(Float64, x[3])
        (isempty(strip(x[1])) || ra === nothing || dec === nothing) && continue
        push!(targets, (; name = String(strip(x[1])), ra, dec, cal = length(x) >= 4 && x[4] in ("1", "true")))
        push!(fors, length(x) >= 5 ? filter(!isempty, String.(strip.(split(x[5], ';')))) : String[])
    end
    isempty(targets) && return "! no targets"
    cals = [Int[] for _ in targets]
    if needs_cal
        for (j, t) in enumerate(targets), nm in fors[j]
            t.cal || continue
            i = findfirst(s -> s.name == nm && !s.cal, targets)
            i === nothing || push!(cals[i], j)
        end
    end
    obs = [night_observability(f, t.ra, t.dec, date; alt_limit = Float64(alt_limit),
                               alt_max = Float64(alt_max), config = cfg, ut_window = uw)
           for t in targets]
    sky = falses(length(first(obs).lst), length(targets))
    for (j, o) in enumerate(obs)
        sky[intersect(o.good_alt, o.good_twilight, o.good_moon), j] .= true
    end
    return (; f, date, cfg, uw, targets, cals, obs, lst = first(obs).lst, sky)
end

# A shared plan as text, `lst0 lst1 config 1 POPs` per segment, `;`-separated: the form
# `pop_plan_text` gives a target's plan, so the panel reads the first POPs the same way.
_shared_plan_text(plan::TargetsPopPlan, cfg) =
    join((@sprintf("%.4f %.4f %s 1 %s", sg.lst0, sg.lst1, join(string.(Int.(cfg .> 0))),
                   join(string.(sg.pop), " ")) for sg in plan.segments), ";")

# The POP sequence of a shared plan in the panel's telescope order, `→` before each change and
# the time it starts on the clock `clock` (LST hours → hours).
function _shared_plan_summary(plan::TargetsPopPlan, f, cfg, clock)
    on = [i for i in telescope_display_order(f) if i <= length(cfg) && cfg[i] > 0]
    pops(sg) = join((_pop_text(sg.pop, cfg, i) for i in on), " ")
    isempty(plan.segments) && return ""
    out = pops(plan.segments[1])
    for k in 2:length(plan.segments)
        a, b = plan.segments[k-1], plan.segments[k]
        a.pop == b.pop && continue
        out *= " → " * pops(b) * " " * _hhmm(clock(a.lst1 + 1 / 60))
    end
    return out
end

"""
    shell_gantt_targets(facility, rows, dateiso, telescopes, autopops, pops, max_changes,
                        change_minutes, min_stretch, bracket_minutes, needs_cal, alt_limit,
                        alt_max, time_system, selected, detailed, night = "",
                        plan_text = "") -> String

A list of targets on one night under one shared POP plan, drawn on the Observe canvas: a row per
target, calibrators in the calibrator colours, and the POPs row on top
(`gantt_targets_geometry`). `rows` is one target per line, `name \\t ra \\t dec \\t cal \\t for`
(degrees; `cal` `1` for a calibrator; `for` the science targets it is for, `;`-separated).
With `autopops` the plan is `pop_plan_targets`' with at most `max_changes` changes (each
`change_minutes`), settings held `min_stretch` minutes, and calibrators within
`bracket_minutes` of their targets when `needs_cal`; `plan_text`, one of
`shell_pop_plans_targets`' plans, is drawn instead when the same search still gives it. Without
`autopops`, every target under the POPs `pops`. `selected` (1-based in `rows`) is the target
the detailed view and the delay chart show. `night` as in [`shell_gantt`](@ref).

Returns `summary \\t dark window \\t the selected target's plan (pop_plan_text) \\t its first
POPs`, or a message beginning with `!`.
"""
function shell_gantt_targets(facility::AbstractString, rows::AbstractString,
                             dateiso::AbstractString, telescopes::AbstractString, autopops::Bool,
                             pops::AbstractString, max_changes::Integer, change_minutes::Real,
                             min_stretch::Real, bracket_minutes::Real, needs_cal::Bool,
                             alt_limit::Real, alt_max::Real, time_system::AbstractString,
                             selected::Integer, detailed::Bool, night::AbstractString = "",
                             plan_text::AbstractString = "")
    sh = _shell()
    tn = try
        _targets_night(facility, rows, dateiso, telescopes, alt_limit, alt_max, night, needs_cal)
    catch err
        "! " * _cause(err)
    end
    tn isa String && (console!(sh, tn); return tn)
    f, cfg, targets = tn.f, tn.cfg, tn.targets
    n = length(targets)
    tsys = _time_system(time_system)
    shared = try
        if autopops
            console!(sh, "> pop_plans_targets(facility, targets, lst, $(cfg); observable, " *
                         "calibrators = $(tn.cals), bracket_minutes = $(bracket_minutes), " *
                         "max_changes = $(max_changes), change_minutes = $(change_minutes), " *
                         "min_stretch = $(min_stretch))")
            ps = pop_plans_targets(f, targets, tn.lst, cfg; observable = tn.sky,
                                   calibrators = tn.cals, bracket_minutes,
                                   max_changes = max(0, Int(max_changes)), change_minutes,
                                   min_stretch)
            k = findfirst(p -> _shared_plan_text(p, cfg) == String(plan_text), ps)
            ps[something(k, length(ps))]
        else
            nothing
        end
    catch err
        msg = "! POPs for the list: " * _cause(err)
        console!(sh, msg); return msg
    end
    pp = isempty(strip(String(pops))) ? ones(Int, length(cfg)) : [parse(Int, x) for x in split(String(pops))]
    plans = NightPlan[]
    for (k, t) in enumerate(targets)
        popplan = shared === nothing ? nothing :
                  target_pop_plan(shared, k, t.ra, cfg; ha_start = unwrap_hour_angles(tn.obs[k].ha)[1])
        push!(plans, night_plan(f, t.name, t.ra, t.dec, tn.date; config = cfg,
                                pop = shared === nothing || isempty(shared.segments) ? pp :
                                      first(shared.segments).pop,
                                popplan, use_delay = true, change_minutes = Float64(change_minutes),
                                alt_limit = Float64(alt_limit), alt_max = Float64(alt_max),
                                ut_window = tn.uw))
    end
    if shared === nothing
        # one setting for the whole night, as a plan of one segment
        x = unwrap_lst(first(plans).lst)
        shared = TargetsPopPlan([SharedSegment(first(x), last(x), pp)], x, falses(length(x), n),
                                [round(Int, observable_hours(p) * 60) for p in plans], 0, 0.0, Int[])
    end
    sel = clamp(Int(selected), 1, n)
    p = plans[sel]
    labels = [t.cal ? t.name * " Cal" : t.name for t in targets]
    cal = [t.cal for t in targets]
    title = string(n, " targets · ", shared.changes, shared.changes == 1 ? " change" : " changes",
                   autopops ? @sprintf(" · score %.2f", shared.score) : "")
    GANTT_TIME_SYSTEM[] = tsys
    GANTT_CAL[] = cal[sel]
    GANTT_STACK[] = detailed ? nothing :
                    (; kind = :targets, plans, labels, cal, shared, config = cfg, title)
    if sh.gantt !== nothing
        detailed ? update_gantt!(sh.gantt, p; detailed = true, time_system = tsys, cal = cal[sel]) :
                   draw_gantt_stack!(sh.gantt, GANTT_STACK[]; time_system = tsys)
    end
    sh.delayplot === nothing || isempty(p.baselines) || update_delay_plot!(sh.delayplot, p; time_system = tsys)
    sh.plan = p
    hours = join((@sprintf("%s %.1f", t.name, observable_hours(q)) for (t, q) in zip(targets, plans)), " · ")
    line = string(title, ": ", hours)
    console!(sh, "  " * line)
    isempty(shared.ignored) ||
        console!(sh, "  calibrators not used: " * join((targets[j].name for j in shared.ignored), ", "))
    dark = isempty(p.good_twilight) ? "no astronomical night" :
           string(_hhmm(p.lst[first(p.good_twilight)]), " – ", _hhmm(p.lst[last(p.good_twilight)]), " LST")
    ptext = autopops ? pop_plan_text(target_pop_plan(shared, sel, targets[sel].ra, cfg;
                                                     ha_start = unwrap_hour_angles(tn.obs[sel].ha)[1])) : ""
    first_pops = isempty(shared.segments) ? pp : first(shared.segments).pop
    return join((line, dark, ptext, join(string.(first_pops), " ")), "\t")
end

"""
    shell_pop_plans_targets(facility, rows, dateiso, telescopes, max_changes, change_minutes,
                            min_stretch, bracket_minutes, needs_cal, alt_limit, alt_max,
                            time_system, night = "") -> String

The POP strategies tab's shared plans for a list of targets (`rows` as in
[`shell_gantt_targets`](@ref)): the best with at most 0 … `max_changes` changes, each distinct
plan once, best first. One row per plan, `summary \\t changes \\t score \\t hours of the science
targets \\t hours of all \\t plan`, the summary being the POP sequence with each change's time on
the clock `time_system`, `plan` the text `shell_gantt_targets` takes. `""` when no plan observes
anything; a message beginning with `!` on failure.
"""
function shell_pop_plans_targets(facility::AbstractString, rows::AbstractString,
                                 dateiso::AbstractString, telescopes::AbstractString,
                                 max_changes::Integer, change_minutes::Real, min_stretch::Real,
                                 bracket_minutes::Real, needs_cal::Bool, alt_limit::Real,
                                 alt_max::Real, time_system::AbstractString,
                                 night::AbstractString = "")
    sh = SHELL[]
    say(t) = sh === nothing ? nothing : console!(sh, t)
    try
        tn = _targets_night(facility, rows, dateiso, telescopes, alt_limit, alt_max, night, needs_cal)
        tn isa String && (say(tn); return tn)
        f, cfg, targets = tn.f, tn.cfg, tn.targets
        say("> pop_plans_targets(facility, targets, lst, $(cfg); observable, calibrators = $(tn.cals), " *
            "bracket_minutes = $(bracket_minutes), max_changes = $(max_changes), " *
            "change_minutes = $(change_minutes), min_stretch = $(min_stretch))")
        ps = pop_plans_targets(f, targets, tn.lst, cfg; observable = tn.sky, calibrators = tn.cals,
                               bracket_minutes, max_changes = max(0, Int(max_changes)),
                               change_minutes, min_stretch)
        filter!(p -> !isempty(p.segments), ps)
        unique!(p -> _shared_plan_text(p, cfg), ps)
        isempty(ps) && (say("  no POP plan observes these targets"); return "")
        sort!(ps; by = p -> (p.score, -p.changes), rev = true)
        tsys = _time_system(time_system)
        p0 = night_plan(f, "", 0.0, 0.0, tn.date; config = cfg,
                        alt_limit = Float64(alt_limit), alt_max = Float64(alt_max))
        clock(x) = tsys === :lst ? mod(x, 24) : hour_angle_clock(p0, x; system = tsys)
        sci = findall(t -> !t.cal, targets)
        isempty(sci) && (sci = collect(eachindex(targets)))
        out = String[]
        for p in ps
            push!(out, join((_shared_plan_summary(p, f, cfg, clock), string(p.changes),
                             @sprintf("%.2f", p.score), @sprintf("%.1f", sum(p.minutes[sci]) / 60),
                             @sprintf("%.1f", sum(p.minutes) / 60), _shared_plan_text(p, cfg)), "\t"))
        end
        say("  best: " * _shared_plan_summary(first(ps), f, cfg, clock))
        return join(out, "\n")
    catch err
        msg = "! POP plans for the list: " * _cause(err)
        say(msg)
        return msg
    end
end

"The POP strategies tab's candidate plans, drawn as stacked Gantt rows (`shell_plans_gantt`)."
const PLANS_GANTT = Ref{Any}(nothing)

"""
    shell_plans_gantt(facility, name, ra, dec, dateiso, telescopes, rows, use_delay,
                      alt_limit, alt_max, time_system, change_minutes, night = "") -> String

Draw candidate plans of one target as stacked Gantt rows on the POP strategies tab's chart, the
first at the top. `rows` is one candidate per line, `label \t pops \t plan_text`: a single
setting by its POPs (`plan_text` empty), a multi-POP plan by its `pop_plan_text`. Returns the
chart's height in pixels at which every row keeps the Gantt's own scale (`GANTT_PX_PER_UNIT`),
or a message beginning with `!`.
"""
function shell_plans_gantt(facility::AbstractString, name::AbstractString, ra::Real, dec::Real,
                           dateiso::AbstractString, telescopes::AbstractString,
                           rows::AbstractString, use_delay::Bool,
                           alt_limit::Real = DEFAULT_ALT_LIMIT, alt_max::Real = DEFAULT_ALT_MAX,
                           time_system::AbstractString = "lst", change_minutes::Real = 10,
                           night::AbstractString = "")
    g = PLANS_GANTT[]
    uw = try
        _ut_window(night)
    catch
        return "! night: not from-to UT hours: $(night)"
    end
    plans = NightPlan[]; labels = String[]
    for line in split(String(rows), '\n'; keepempty = false)
        f = split(line, '\t')
        length(f) >= 2 || continue
        np = _night_plan_for(facility, name, ra, dec, dateiso, telescopes, f[2], use_delay,
                             alt_limit, alt_max, length(f) >= 3 ? f[3] : "", change_minutes;
                             ut_window = uw)
        np isa String && return np
        push!(plans, np.plan); push!(labels, String(f[1]))
    end
    isempty(plans) && return "! no candidates"
    rows_px = ((length(plans) - 1) * STACK_ROW_STEP + 4.0) * GANTT_PX_PER_UNIT * live_plot_scale()
    px = string(ceil(Int, rows_px + 70 * live_plot_scale()))     # and the time axis and title
    g === nothing && return px
    update_gantt_stack!(g, plans, labels, string(name, " · ", length(plans), " candidates");
                        time_system = _time_system(time_system))
    return px
end

"""
    shell_pop_plans(facility, dec, dateiso, ra, telescopes, reference, max_changes, max_drop,
                    change_minutes, alt_limit, alt_max, time_system = "lst", min_stretch = 0,
                    night = "") -> String

The POP strategies tab's plans for one target, best first: the best with at most 0 … `max_changes` POP
changes and 0 … `max_drop` telescopes dropped, each distinct plan once, ranked as `pop_plan`
ranks them. One row per plan, `summary \t changes \t dropped \t hours with every telescope
\t hours in all \t plan`, `plan` being the `pop_plan_text` that `shell_gantt` takes.
`reference` is the telescope protected from being dropped (a name, or empty). The summary's
times are on the clock `time_system`. `min_stretch` as in `pop_plan`; `night` as in
[`shell_gantt`](@ref). `""` when no plan observes anything.
"""
function shell_pop_plans(facility::AbstractString, dec::Real, dateiso::AbstractString,
                         ra::Real, telescopes::AbstractString, reference::AbstractString,
                         max_changes::Integer, max_drop::Integer, change_minutes::Real,
                         alt_limit::Real = DEFAULT_ALT_LIMIT, alt_max::Real = DEFAULT_ALT_MAX,
                         time_system::AbstractString = "lst", min_stretch::Real = 0,
                         night::AbstractString = "")
    sh = SHELL[]
    say(t) = sh === nothing ? nothing : console!(sh, t)
    date = try
        Dates.DateTime(String(dateiso))
    catch
        say("! POP plans: '$(dateiso)' is not a date"); return ""
    end
    tsys = _time_system(time_system)
    try
        f = read_facility_file(String(facility))
        tels = facility_telescopes(String(facility))
        sel  = filter(!isempty, String.(split(String(telescopes))))
        cfg  = isempty(sel) ? ones(Int, length(tels)) : telescope_config(sel, tels)
        ref  = something(findfirst(==(String(reference)), tels), 0)
        ha, ok = _night_sky(f, ra, dec, date, cfg; alt_limit, alt_max, ut_window = _ut_window(night))
        kc, kd = max(0, Int(max_changes)), max(0, Int(max_drop))
        say("> pop_plans(facility, $(round(Float64(dec); digits=4)), ha, $(cfg); observable, " *
            "reference = $(ref), max_changes = $(kc), max_drop = 0:$(kd), change_minutes = $(change_minutes), " *
            "min_stretch = $(min_stretch))")
        plans = PopPlan[]
        for d in 0:kd
            for pl in pop_plans(f, Float64(dec), ha, cfg; observable = ok, reference = ref,
                                max_changes = kc, max_drop = d, change_minutes, min_stretch)
                isempty(pl.segments) && continue
                t = pop_plan_text(pl)
                any(q -> pop_plan_text(q) == t, plans) || push!(plans, pl)
            end
        end
        isempty(plans) && (say("  no POP plan observes this target"); return "")
        sort!(plans; by = q -> (q.minutes..., -q.changes), rev = true)
        clock = if tsys === :lst
            nothing
        else
            p0 = night_plan(f, "", Float64(ra), Float64(dec), date; config = cfg,
                            alt_limit = Float64(alt_limit), alt_max = Float64(alt_max))
            x -> hour_angle_clock(p0, x - Float64(ra) / 15; system = tsys)
        end
        rows = String[]
        for pl in plans
            dropped = unique(i for sg in pl.segments for i in eachindex(cfg) if cfg[i] > 0 && sg.config[i] == 0)
            push!(rows, join((pop_plan_summary(pl, f, ra, cfg; clock), string(pl.changes),
                              isempty(dropped) ? "—" : join(tels[dropped], " "),
                              @sprintf("%.1f", pl.minutes[1] / 60), @sprintf("%.1f", sum(pl.minutes) / 60),
                              pop_plan_text(pl)), "\t"))
        end
        say("  best: " * pop_plan_summary(first(plans), f, ra, cfg))
        return join(rows, "\n")
    catch err
        say("! POP plans failed: " * _cause(err))
        return ""
    end
end

"""
    shell_best_pops(facility, dec, dateiso, ra, telescopes, n = 5,
                    alt_limit = DEFAULT_ALT_LIMIT, alt_max = DEFAULT_ALT_MAX,
                    time_system = "lst", night = "") -> String

The best POP configurations for this target and night, as `pop\tscore\thaRange\ttimeRange`
rows, the time range on the clock `time_system`. They are
scored over the hours the target is observable: above `alt_limit` and the horizon, below
`alt_max`, in astronomical night (within `night`, as in [`shell_gantt`](@ref)), clear of the
moon.
"""
function shell_best_pops(facility::AbstractString, dec::Real, dateiso::AbstractString,
                         ra::Real, telescopes::AbstractString, n::Integer = 5,
                         alt_limit::Real = DEFAULT_ALT_LIMIT, alt_max::Real = DEFAULT_ALT_MAX,
                         time_system::AbstractString = "lst", night::AbstractString = "")
    sh = SHELL[]
    say(t) = sh === nothing ? nothing : console!(sh, t)
    date = try
        Dates.DateTime(String(dateiso))
    catch
        say("! POP search: '$(dateiso)' is not a date"); return ""
    end
    try
        f = read_facility_file(String(facility))
        tels = facility_telescopes(String(facility))
        sel  = filter(!isempty, String.(split(String(telescopes))))
        cfg  = isempty(sel) ? ones(Int, length(tels)) : telescope_config(sel, tels)
        obs = night_observability(f, Float64(ra), Float64(dec), date;
                                  alt_limit = Float64(alt_limit), alt_max = Float64(alt_max),
                                  config = cfg, ut_window = _ut_window(night))
        ha = obs.ha[intersect(obs.good_alt, obs.good_twilight, obs.good_moon)]
        if isempty(ha)
            say("  POP search: the target is not observable this night")
            return ""
        end
        say("> best_pop(facility, $(round(Float64(dec); digits=4)), ha, $(cfg))")
        res = best_pops(f, Float64(dec), ha, cfg; n = Int(n))
        if isempty(res)
            say("  no POP configuration reaches this target with these telescopes")
            return ""
        end
        say("  best: POPs $(join(res[1].pop, " ")) — $(res[1].score) steps")

        # `pops \t score \t haRange \t timeRange` — the table promises an HA range, and it is
        # worth the extra in_delay per candidate (microseconds): a configuration is chosen for
        # WHEN it works, not only for how long. The time range is on the clock `time_system`.
        tsys = _time_system(time_system)
        p0 = night_plan(f, "", Float64(ra), Float64(dec), date; config = cfg,
                        alt_limit = Float64(alt_limit), alt_max = Float64(alt_max))
        rows = String[]
        for r in res
            span, clk = try
                g = in_delay(f, Float64(dec), ha, cfg, Int.(r.pop)).good_delay
                isempty(g) ? ("—", "—") :
                    (string(round(ha[first(g)]; digits = 2), " … ",
                            round(ha[last(g)];  digits = 2), " h"),
                     string(_hhmm(hour_angle_clock(p0, ha[first(g)]; system = tsys)), " – ",
                            _hhmm(hour_angle_clock(p0, ha[last(g)]; system = tsys))))
            catch
                ("—", "—")
            end
            push!(rows, join((join(string.(r.pop), " "), string(r.score), span, clk), "\t"))
        end
        return join(rows, "\n")
    catch err
        say("! POP search failed: " * _cause(err))
        return ""
    end
end

"""
    shell_bin_info() -> String

How the current dataset is binned and what it could be, as
`nbins <TAB> nchannels <TAB> has_diffvis <TAB> mergeable`.

`nchannels` is the number of DISTINCT wavelengths in the file, which is the most bins that
mean anything -- past that a bin holds no data. `mergeable` says whether the file carries
several OI_WAVELENGTH tables sharing one grid, the case where reading per channel without
merging silently multiplies the data by the number of tables.
"""
function shell_bin_info()
    e = current_dataset(_shell())
    e === nothing && return ""
    nb = size(e.data, 1)
    lams = Float64[]
    for d in e.data
        append!(lams, Float64.(d.uv_lam))
    end
    nch = length(unique(round.(lams; sigdigits = 10)))
    hasdv = any(d -> d.nvisphi > 0, e.data)
    return string(nb, '\t', nch, '\t', hasdv ? 1 : 0, '\t',
                  get(e.readopts, :merge_oi_wavelength, false) ? 1 : 0)
end

"""
    shell_rebin(mode, nbins) -> String

Re-read the current dataset with a different spectral binning, and log the `readoifits` call.

`mode` is `single`, `channels` or `groups`:

  * `single` -- one bin over the whole band, which is what every monochromatic engine wants.
  * `channels` -- one image plane per spectral channel, ALWAYS with `merge_oi_wavelength`.
    That is not a preference. A file carrying several OI_WAVELENGTH tables on one grid --
    v1295 Aql has seven -- otherwise yields one bin per channel PER TABLE, and since the grids
    coincide every measurement lands in every matching bin and is counted once per table.
    `readoifits` warns about the tables but not about that consequence, which inflates chi2 and
    over-weights the data against the regularisers.
  * `groups` -- `nbins` equal windows across the band. A window then holds several channels,
    and that is what forces `use_vis = false`: the cross-channel differential term pairs one
    baseline per bin, so a bin holding many channels makes the pairing ambiguous and
    `image_to_chi2` raises a DimensionMismatch. Dropping OI_VIS costs nothing in the fit --
    the imaging criterion is V2, T3amp and T3phi, and a reconstruction with and without it is
    bit-identical -- but it is a real observable to lose from the reported chi2, so it is said
    out loud rather than done quietly.

Re-reading, not re-slicing: binning happens inside `readoifits` and there is no way to rebin an
`OIdata` after the fact. The entry keeps its name and its place, so everything pointing at it
still does.
"""
function shell_rebin(mode::AbstractString, nbins::Integer = 0)
    sh = _shell()
    e = current_dataset(sh)
    e === nothing && return "! no dataset loaded"
    # Re-reading the file would silently drop the cut; rebinning comes before filtering.
    e.derived_from === nothing ||
        return "! $(e.name) is a filtered copy of $(e.derived_from): rebin $(e.derived_from), " *
               "then filter it again"
    m = String(mode)

    kw = Dict{Symbol,Any}()
    for (k, v) in e.readopts
        k in (:polychromatic, :merge_oi_wavelength, :spectralbin, :use_vis) || (kw[k] = v)
    end

    if m == "single"
        # nothing to add: one bin is readoifits' own default
    elseif m == "channels"
        kw[:polychromatic] = true
        kw[:merge_oi_wavelength] = true
    elseif m == "groups"
        n = Int(nbins)
        n > 1 || return "! grouping needs at least 2 bins"
        lams = Float64[]
        for d in e.data
            append!(lams, Float64.(d.uv_lam))
        end
        isempty(lams) && return "! this dataset carries no wavelengths"
        lo, hi = extrema(lams)
        hi > lo || return "! this dataset is at a single wavelength; there is nothing to bin"
        nch = length(unique(round.(lams; sigdigits = 10)))
        n <= nch || return "! only $(nch) distinct wavelengths; $(n) bins would leave some empty"
        edges = range(lo - 1e-12, hi + 1e-12; length = n + 1)
        kw[:spectralbin] = [[edges[i], edges[i+1]] for i in 1:n]
        if any(d -> d.nvisphi > 0, e.data)
            kw[:use_vis] = false
            console!(sh, "  grouped bins hold several channels, so OI_VIS is dropped: the " *
                         "differential term pairs one baseline per bin. The fit is unchanged " *
                         "— it uses V2, T3amp and T3phi — only the reported chi2 loses a term.")
        end
    else
        return "! unknown binning mode: " * m
    end

    data = try
        OITOOLS.readoifits(sources(e); warn = false, verbose = false, kw...)
    catch err
        msg = "! could not re-read " * source_label(e) * ": " * _cause(err)
        console!(sh, msg; kind = :err); return msg
    end
    data isa AbstractArray || return "! readoifits returned nothing for " * source_label(e)

    e.data = data
    e.readopts = kw
    e.spectralbin = get(kw, :spectralbin, nothing)

    kwstr = _kwargs(Pair{Symbol,Any}[k => v for (k, v) in sort(collect(kw), by = first)];
                    defaults = _READOIFITS_DEFAULTS)
    call = isempty(kwstr) ? "readoifits($(_literal(sources(e))))" :
                            "readoifits($(_literal(sources(e))); $kwstr)"
    log!(sh.session, "$(e.name) = $call"; note = "rebin $(source_label(e))", binding = e.name)
    console!(sh, "$(e.name) = $call"; kind = :cmd)

    sh.ftcache = nothing          # the plan is per bin; a rebin invalidates every cell
    line = "$(size(data, 1)) wavelength bin" * (size(data, 1) == 1 ? "" : "s") *
           " × $(size(data, 2)) epoch" * (size(data, 2) == 1 ? "" : "s")
    console!(sh, line)
    sh.status = refresh_plot!(sh)
    return line
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
                           options::AbstractString = "",
                           startpath::AbstractString = "")
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

    # `startpath` is carried through: without it a FITS start worked from "Show start image",
    # which passes it, and silently fell back to whatever `start_image` does with an empty path
    # from Run, which did not -- the two buttons disagreeing about the starting image.
    setup = ImagingSetup(; engine = eng,
                           nx = Int(nx), pixsize = Float64(pixsize),
                           mode = Symbol(mode), startkind = Symbol(startkind),
                           startfwhm = Float64(startfwhm), startseed = Int(startseed),
                           startpath = String(startpath))
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
    console!(sh, "> $(engine_call_name(eng, opts))(x_start, data, ft; maxiter = $(Int(maxiter)))  " *
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
# The title for a finished reconstruction, with whichever channel the panel is showing.
# Here rather than in `image_title` because the channel lives on the shell, not on the result.
function _result_title(sh::ShellState, r)
    n = result_channels(r)
    λ = (n > 1 && !isempty(r.wavelengths) && 1 <= sh.imchannel <= length(r.wavelengths)) ?
        r.wavelengths[sh.imchannel] : NaN
    return image_title(r.setup; channel = sh.imchannel, nchannels = n, wavelength = λ)
end

function finish_reconstruct!(sh::ShellState, res)
    r = res.result
    sh.enginelog = res.output
    if !res.ok
        msg = "! reconstruction failed: " * _cause(res.err)
        sh.enginelog = isempty(res.output) ? msg : res.output * "\n" * msg
        console!(sh, msg)
        return msg
    end
    sh.imcanvas === nothing || show_image!(sh.imcanvas, result_plane(r, sh.imchannel),
                                          r.setup.pixsize; title = _result_title(sh, r))
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

"""
    shell_tempering_diagnostics() -> String

The tempered sampler's health, as `name=value` pairs, or `""` when the last run was not one.

Plain numbers, not a verdict: the panel owns the thresholds, because they are what a user reads
and argues with, and `demos/pigeons_diagnostics.md` is where they are written down. `NaN` for a
quantity Pigeons could not supply, which the panel shows as "—" rather than as a zero — a swap
acceptance of 0 means a dead rung, and printing that for "not measured" would invent a fault.

`restarts` is the field that decides whether the run means anything at all: zero round trips
means the chains never communicated, and no other number rescues that.
"""
function shell_tempering_diagnostics()
    r = _shell().imaging
    (r === nothing || !(r.extra isa NamedTuple) || !haskey(r.extra, :pt_health)) && return ""
    h = r.extra.pt_health
    fields = (:restarts, :barrier, :swap_min, :swap_mean, :rho_max, :rho_mean,
              :explorer_min, :explorer_mean, :logz, :logz_delta)
    parts = [string(f, "=", getfield(h, f)) for f in fields]
    # The ladder, one acceptance per adjacent pair, so a single bad rung can be located rather
    # than merely reported as "min(α) = 0".
    push!(parts, "rungs=" * join(h.swap_per_rung, "|"))
    return join(parts, ",")
end

"""
    shell_channels() -> String

The cube's channels, as `count<TAB>current<TAB>label1|label2|...`, or `""` for a grey run.

`""` is the answer the panel needs for "there is nothing to scrub": a grey reconstruction has
one image and a slider over it would be a control that cannot do anything.

The labels are wavelengths in µm because that is what an observer identifies a channel by; the
INDEX is what every other call takes, since two bins can round to the same µm string.
"""
function shell_channels()
    r = _shell().imaging
    r === nothing && return ""
    n = result_channels(r)
    n > 1 || return ""
    labels = isempty(r.wavelengths) ? ["channel $w" for w in 1:n] :
             [string(round(l * 1e6; digits = 4), " µm") for l in r.wavelengths]
    return string(n, '\t', _shell().imchannel, '\t', join(labels, '|'))
end

"""
    shell_show_channel(w) -> String

Draw channel `w` of the cube, and remember it.

Remembered rather than passed through, because the channel is not only a drawing choice: the
residual overlay, the model observables and the saved PNG all answer "one of what?" with it,
and a panel where the picture showed one channel while the residuals showed another would be
worse than no channel control at all.
"""
function shell_show_channel(w::Integer)
    sh = _shell()
    r  = sh.imaging
    r === nothing && return "! nothing reconstructed yet"
    n = result_channels(r)
    sh.imchannel = clamp(Int(w), 1, n)
    sh.imcanvas === nothing ||
        show_image!(sh.imcanvas, result_plane(r, sh.imchannel), r.setup.pixsize;
                    title = _result_title(sh, r))
    lab = isempty(r.wavelengths) ? "channel $(sh.imchannel)" :
          string(round(r.wavelengths[sh.imchannel] * 1e6; digits = 4), " µm")
    return "channel $(sh.imchannel) of $n — $lab"
end

"""
    shell_result_ensemble() -> String

What the last reconstruction offers beyond one image, as `nsamples<TAB>source`, or `""`.

`""` means the engine returned a point estimate and the panel's mean/spread/sample controls
have nothing to show — which is every engine but the samplers.
"""
function shell_result_ensemble()
    r = _shell().imaging
    (r === nothing || r.ensemble === nothing) && return ""
    return string(length(r.ensemble.samples), '\t', r.ensemble.source)
end

"""
    shell_show_result(mode, index) -> String

Put one view of the last reconstruction on the canvas and return a one-line description.

`mode` is `result`, `mean`, `sigma` or `sample`:

  * `result` is the image the engine actually returned. For annealing that is the BEST chain,
    not the mean of them, so it and `mean` are different images and both are offered.
  * `mean` averages the ensemble.
  * `sigma` is the spread across it — for annealing, chain-to-chain disagreement rather than a
    within-chain posterior width. The colourbar is relabelled, because a spread map shares no
    units with a flux map and an unchanged legend would say it did.
  * `sample` is one member, `index` counted from 1.

GUI thread only: `show_image!` is a Makie call.
"""
function shell_show_result(mode::AbstractString, index::Integer = 1)
    sh = _shell()
    r  = sh.imaging
    r === nothing && return "! nothing reconstructed yet"
    sh.imcanvas === nothing && return "! no canvas"
    m = String(mode)
    px = r.setup.pixsize

    if m == "result"
        show_image!(sh.imcanvas, result_plane(r, sh.imchannel), px; title = _result_title(sh, r))
        return "showing the reconstruction"
    end

    e = r.ensemble
    e === nothing && return "! this engine returned one image, not an ensemble"

    if m == "mean"
        show_image!(sh.imcanvas, e.mean, px; label = "mean flux / pixel")
        return "showing the mean of " * e.source
    elseif m == "sigma"
        e.sigma === nothing && return "! only one member; there is no spread to show"
        show_image!(sh.imcanvas, e.sigma, px; label = "spread / pixel")
        return "showing the spread across " * e.source
    elseif m == "sample"
        n = length(e.samples)
        i = clamp(Int(index), 1, n)
        show_image!(sh.imcanvas, e.samples[i], px)
        return "showing member $(i) of $(n)"
    end
    return "! unknown result view: " * m
end

"The last reconstruction, or an empty string when there is none."
shell_imaging_summary() =
    _shell().imaging === nothing ? "" : sprint(show, _shell().imaging)

# The plans the Gantt shows as stacked rows, `(; kind, plans, labels, title, …)`: `kind`
# `:coverage` for one target's best settings over the night's thirds (`stack_mode`), `:targets`
# for a list of targets under one shared POP plan, with `cal`, `shared` (the
# `TargetsPopPlan`) and `config`; `nothing` when it shows one plan.
const GANTT_STACK = Ref{Any}(nothing)

# Draw a `GANTT_STACK` value on a Gantt.
function draw_gantt_stack!(g, st; time_system::Symbol = :lst, cal::Bool = false)
    st.kind === :coverage && return update_gantt_stack!(g, st.plans, st.labels, st.title; time_system, cal)
    st.kind === :targets && return update_gantt_targets!(g, st.plans, st.labels, st.cal, st.shared,
                                                         st.config, st.title; time_system)
    error("unknown Gantt stack: $(st.kind)")
end

# The geometry of a `GANTT_STACK` value, for the hover.
function _stack_geometry(st)
    st.kind === :coverage && return gantt_stack_geometry(st.plans, st.labels, st.title)
    st.kind === :targets && return gantt_targets_geometry(st.plans, st.labels, st.cal, st.shared,
                                                          st.config, st.title)
    error("unknown Gantt stack: $(st.kind)")
end

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
    stack = GANTT_STACK[]
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
    geo = stack === nothing ? gantt_geometry(p) : _stack_geometry(stack)
    # Which row: the nearest, and only if the pointer is actually within half a row of it.
    isempty(geo.rows) && return ""
    _, k = findmin(r -> abs(r[1] - y), geo.rows)
    row = geo.rows[k]
    abs(row[1] - y) > 0.6 && return ""
    # the shared POPs row of a list: the segment under the cursor
    if stack !== nothing && k > length(stack.plans)
        sp = stack.shared
        x0 = unwrap_lst(first(stack.plans).lst)[1]
        onx(l) = l + 24 * ceil((x0 - 2 - l) / 24)
        j = findfirst(sg -> onx(sg.lst0) <= x <= onx(sg.lst0) + (sg.lst1 - sg.lst0), sp.segments)
        out = ["POPs", @sprintf("%d changes · score %.2f", sp.changes, sp.score)]
        if j !== nothing
            sg = sp.segments[j]
            on = [i for i in first(stack.plans).display_order if stack.config[i] > 0]
            push!(out, @sprintf("segment %d of %d: %s", j, length(sp.segments),
                                join((_pop_text(sg.pop, stack.config, i) for i in on), " ")))
            push!(out, string("LST ", _hhmm(mod(sg.lst0, 24)), " – ", _hhmm(mod(sg.lst1, 24))))
        end
        return join(out, "\n")
    end
    stack === nothing || (p = stack.plans[k])

    # Nearest sample in time. The grid is regular, so this is the instant under the cursor.
    isempty(p.lst) && return ""
    _, i = findmin(t -> abs(t - x), p.lst)

    out = String[]
    push!(out, row[2])
    push!(out, @sprintf("LST %s   HA %+.2f h", _hhmm(p.lst[i]), p.ha[i]))
    push!(out, @sprintf("alt %.1f°   az %.1f°", p.alt[i], p.az[i]))
    push!(out, @sprintf("elevation window %.0f–%.0f°", p.alt_limit, p.alt_max))
    inside(v, idx) = i in idx ? v : "no " * v
    push!(out, join(filter(!isempty, [inside("dark", p.good_twilight),
                     inside("within limits", p.good_riseset),
                     p.has_horizon ? inside("above horizon", p.good_horizon) : "",
                     inside("moon ok", p.good_moon),
                     p.delay_applied ? inside("in delay", p.good_delay) : "delay not checked"]),
                    " · "))
    push!(out, @sprintf("night: %.2f h observable   moon %d%% lit",
                        observable_hours(p), round(Int, 100 * p.moon_fli)))
    if is_multi_pop(p)
        k = findfirst(sg -> sg.i0 <= i <= sg.i1, p.segments)
        g = findfirst(r -> r[1] <= i < r[2], p.changes)
        g !== nothing ? push!(out, "POP change") :
        k !== nothing ? push!(out, @sprintf("segment %d of %d: %s", k, length(p.segments),
                                            segment_label(p, p.segments[k]))) : nothing
        push!(out, multi_pop_summary(p))
    else
        isempty(p.pop) || push!(out, "POPs " * pop_label(p))
    end
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
    shell_compare_available() -> String

Whether a model and a reconstruction exist to compare against, as `"model=1,image=0"`.

The residual and overplot ticks are disabled on this rather than being allowed to fail: with no
model there is nothing to subtract or draw, and a tick that reports an error is a worse answer
than one that is visibly not available yet.
"""
function shell_compare_available()
    sh = SHELL[]
    sh === nothing && return "model=0,image=0"
    m = try
        mm = _model()
        !isempty(mm.dict)
    catch
        false
    end
    return string("model=", m ? 1 : 0, ",image=", sh.imaging === nothing ? 0 : 1)
end

"""
    shell_optional_engines() -> String

Which optional engines this build can actually run, as `"vi=1,tempering=0,corner=1"`: the two
stochastic imaging engines, and the Model tab's corner plot (PairPlots).

Both are gated on a package extension rather than a dependency, so whether they work is a
property of the session and not of the data. The panel greys the entry and says which package
is missing; before this the two flags were hardcoded `false` in QML, so an engine stayed greyed
even with its package installed — the reason the VI entry has never been reachable.
"""
function shell_optional_engines()
    vi   = Base.get_extension(OITOOLS, :OITOOLSVarInfExt)  !== nothing
    temp = Base.get_extension(OITOOLS, :OITOOLSPigeonsExt) !== nothing
    corner = Base.get_extension(OITOOLS, :OITOOLSPairPlotsExt) !== nothing
    return string("vi=", vi ? 1 : 0, ",tempering=", temp ? 1 : 0, ",corner=", corner ? 1 : 0)
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

Calls [`check_qt_conflict`](@ref) (refuses to start beside another Qt) and
[`configure_graphics!`](@ref) first.
"""
function gui end
