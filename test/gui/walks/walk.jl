# Drive the OITOOLS window from outside, by control name: real clicks and keys (xdotool) on a
# virtual screen (Xvfb), at positions the window itself reports (`walkFind` in Main.qml,
# enabled by `OITOOLSGUI_WALK`). Run by hand, not by the suites: it needs Xvfb, xdotool,
# ImageMagick's `import` and the PMOIRED example data.
#
#     julia test/gui/walks/ex1.jl
#
# Queries: an objectName, `text:Fit` (a visible control whose text is "Fit"), `popup:auto` (the
# same inside an open popup), `a/b` (b inside a).

module Walk

export start, finish, locate, click, shiftclick, typein, choose, tick, value, waitfor, shot,
       check, stage, press, logtext, fitresult, pickfile, spin, fitnow, pastepmoired

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))

mutable struct Window
    name::String
    out::String          # screenshots, log, report
    io::String           # request / reply
    display::String
    xvfb::Base.Process
    gui::Base.Process
    wid::String
    log::String
    results::Vector{Tuple{String,Bool,String}}
end

_x(w::Window, args...) = run(addenv(`xdotool $args`, "DISPLAY" => w.display))

"""
    start(name; data_dir, display = ":117", screen = "1920x1440x24", env = []) -> Window

A fresh Xvfb and a fresh window (`bin/oitoolsgui.jl`, an empty session) under `script`, so Qt's
warnings reach the log line by line; `cmd` replaces the launch (an AppImage, say); `env` adds variables (`"OITOOLSGUI_SCALE" => "1.6"`). Output in `~/.cache/oitools-claude/walks/<name>/`.
"""
function start(name; data_dir = "", display = ":117", screen = "1920x1440x24", env = Pair{String,String}[],
               cmd = "cd $(REPO) && julia --project=bin bin/oitoolsgui.jl")
    out = joinpath(homedir(), ".cache", "oitools-claude", "walks", name)
    io = joinpath(out, "io")
    rm(out; recursive = true, force = true); mkpath(io)
    isfile("/tmp/.X$(display[2:end])-lock") && error("display $display is taken")
    xvfb = run(pipeline(`Xvfb $display -screen 0 $screen`; stdout = devnull, stderr = devnull);
               wait = false)
    sleep(2)
    log = joinpath(out, "run.log")
    # A config directory of its own, so the user's saved settings (start tab, full screen,
    # scales) do not reach the walk; `env` can name another.
    cfg = joinpath(out, "config"); mkpath(cfg)
    env = ["DISPLAY" => display, "OITOOLSGUI_WALK" => io, "OITOOLSGUI_DATA_DIR" => data_dir,
           "XDG_CONFIG_HOME" => cfg, env...]
    gui = run(pipeline(addenv(`script -qfec $cmd $log`, env...); stdout = devnull, stderr = devnull);
              wait = false)
    wid = ""
    for _ in 1:150
        wid = try
            readchomp(addenv(`xdotool search --name OITOOLS`, "DISPLAY" => display))
        catch
            ""
        end
        isempty(wid) || break
        sleep(2)
    end
    isempty(wid) && error("no window after 300 s; see $log")
    w = Window(name, out, io, display, xvfb, gui, first(split(wid, "\n")), log, [])
    # the window answers once its QML timers run
    for _ in 1:60
        try
            locate(w, "text:Exploring"; timeout = 2); break
        catch
            sleep(1)
        end
    end
    sleep(2)
    return w
end

"""
    request(w, s; timeout = 10) -> String

One locator request, answered by the window's timer.
"""
function request(w::Window, s; timeout = 10)
    reply = joinpath(w.io, "reply")
    rm(reply; force = true)
    tmp = joinpath(w.io, "request.tmp")
    write(tmp, s); mv(tmp, joinpath(w.io, "request"); force = true)
    t0 = time()
    while !isfile(reply)
        time() - t0 > timeout && error("no reply to \"$s\" in $(timeout) s")
        sleep(0.05)
    end
    r = read(reply, String); rm(reply; force = true)
    startswith(r, "!\t") && error(r[3:end])
    return r
end

"Where a control is and what it shows: `(name, x, y, w, h, enabled, checked, current, text)`."
function locate(w::Window, q; timeout = 10, reveal = false)
    f = split(request(w, (reveal ? "reveal " : "locate ") * q; timeout), "\t"; limit = 9)
    return (name = f[1], x = parse(Int, f[2]), y = parse(Int, f[3]), w = parse(Int, f[4]),
            h = parse(Int, f[5]), enabled = f[6] == "true", checked = f[7] == "true",
            current = f[8], text = f[9])
end

dump(w::Window) = request(w, "dump")

"Wait until `q` is visible and enabled (the window may still be building it)."
function ready(w::Window, q; timeout = 30)
    t0 = time()
    err = nothing
    while time() - t0 < timeout
        try
            l = locate(w, q; reveal = true)
            l.enabled && l.w > 0 && return l
        catch e
            err = e
        end
        sleep(0.25)
    end
    error("$q not ready after $(timeout) s" * (err === nothing ? "" : ": " * sprint(showerror, err)))
end

function click(w::Window, q; pause = 0.6, mods = String[], at = :centre)
    l = ready(w, q)
    x = l.x + (at === :left ? min(10, l.w ÷ 2) : l.w ÷ 2)
    y = l.y + l.h ÷ 2
    for m in mods; _x(w, "keydown", m); end
    _x(w, "mousemove", "--window", w.wid, string(x), string(y), "click", "1")
    for m in mods; _x(w, "keyup", m); end
    sleep(pause)
    return l
end

shiftclick(w::Window, q; kw...) = click(w, q; mods = ["shift"], kw...)

press(w::Window, k...; pause = 0.3) = (_x(w, "key", "--window", w.wid, k...); sleep(pause))

"Click a field, replace its text with `s`, press Return (or Tab, where Return would accept a dialog)."
function typein(w::Window, q, s; pause = 0.6, enter = "Return")
    click(w, q; pause = 0.2)
    press(w, "ctrl+a")
    _x(w, "type", "--window", w.wid, "--delay", "20", "--", string(s))
    press(w, enter; pause)
end

"Open a combo box and pick the entry whose text is `item`."
function choose(w::Window, q, item; pause = 0.8)
    locate(w, q).current == item && return
    click(w, q; pause = 0.5)
    click(w, "popup:" * item; pause)
    got = locate(w, q).current
    got == item || error("$q shows \"$got\" after choosing \"$item\"")
end

"Set a tick box (by its objectName, else its text) to `on`."
function tick(w::Window, label, on::Bool)
    q = occursin(":", label) ? label :
        (try locate(w, label; timeout = 2); label catch; "text:" * label end)
    locate(w, q).checked == on || click(w, q)
    locate(w, q).checked == on || error("$label did not change to $on")
end

value(w::Window, q) = locate(w, q).text

"Type `path` into the open file picker and confirm; a pending \"Discard\" question is answered first."
function pickfile(w::Window, path; pause = 3)
    try
        locate(w, "text:Discard"; timeout = 2)
        click(w, "text:Discard"; pause = 1)
    catch
    end
    typein(w, "pickerName", path; pause)
end

"Set an editable SpinBox to `n`."
function spin(w::Window, q, n)
    click(w, q; pause = 0.3)
    press(w, "ctrl+a")
    _x(w, "type", "--window", w.wid, "--", string(n))
    press(w, "Tab")
end

"Press Fit and wait until the button reads \"Fit\" again."
function fitnow(w::Window; timeout = 300)
    click(w, "fitButton"; pause = 1)
    waitfor(() -> value(w, "fitButton") == "Fit", "the fit"; timeout)
    sleep(1)
end

"Paste `text` into the Import PMOIRED box and import it."
function pastepmoired(w::Window, text)
    click(w, "text:Import PMOIRED…"; pause = 1)
    try
        locate(w, "text:Discard"; timeout = 2)
        click(w, "text:Discard"; pause = 1)
    catch
    end
    click(w, "pmoiredPaste"; pause = 0.2)
    press(w, "ctrl+a")
    _x(w, "type", "--window", w.wid, "--delay", "2", "--", replace(text, "\n" => " "))
    click(w, "text:Import"; pause = 2)
end

"""
    fitresult(w, n) -> (chi2r, Dict(name => (value, σ)))

Row `fit n` of the Fits table, as shown.
"""
function fitresult(w::Window, n)
    chi2r = parse(Float64, value(w, "fit:fit $n:chi2r"))
    p = Dict{String,Tuple{Float64,Float64}}()
    # parameters apart by three spaces; inside one, non-breaking spaces ("rim,az amp1 = …")
    for item in split(value(w, "fit:fit $n:params"), "   "; keepempty = false)
        m = match(r"^(.+) = (\S+) ± (\S+)$", strip(replace(item, '\u00a0' => ' ')))
        m === nothing || (p[m[1]] = (parse(Float64, m[2]), parse(Float64, m[3])))
    end
    return chi2r, p
end

function waitfor(f, what; timeout = 300, every = 0.5)
    t0 = time()
    while time() - t0 < timeout
        try
            f() && return true
        catch
        end
        sleep(every)
    end
    error("timed out after $(timeout) s waiting for $what")
end

"The terminal output so far (the console's transcript and Qt's messages)."
logtext(w::Window) = read(w.log, String)

function shot(w::Window, name)
    run(addenv(`import -window $(w.wid) $(joinpath(w.out, name * ".png"))`, "DISPLAY" => w.display))
    return name
end

"Record a comparison: `got` read from the window, `want` the script's, within `tol`."
function check(w::Window, label, got::Real, want::Real; tol)
    ok = abs(got - want) <= tol
    push!(w.results, (label, ok, "$(got) (script $(want) ± $(tol))"))
    println(ok ? "PASS  " : "FAIL  ", rpad(label, 44), got, "   script ", want)
    return ok
end

"Run one named step; an error is recorded and stops the walk with a screenshot."
function stage(f, w::Window, label)
    println("── ", label)
    try
        f()
    catch e
        msg = first(split(sprint(showerror, e), "\n"))
        push!(w.results, (label, false, msg))
        println("FAIL  ", label, ": ", msg)
        try shot(w, "fail") catch end
        rethrow()
    end
end

const QT_WARN = r"Qt Warning(?!: failed to create compose table)|qml:[0-9]+|TypeError|ReferenceError|Binding loop|Cannot assign|QQml"

"Stop the window, then Xvfb; scan the log for Qt warnings; print the results."
function finish(w::Window)
    try
        run(pipeline(`pkill -TERM -P $(getpid(w.gui))`; stderr = devnull); wait = true)
    catch
    end
    sleep(3)
    kill(w.gui); kill(w.xvfb)
    warns = filter(l -> occursin(QT_WARN, l), readlines(w.log))
    println("\n── $(w.name): ", count(r -> r[2], w.results), " of ", length(w.results), " passed; ",
            isempty(warns) ? "no Qt warnings" : "$(length(warns)) Qt warning lines:")
    foreach(l -> println("   ", l), first(warns, 20))
    for r in w.results
        r[2] || println("   FAIL ", r[1], ": ", r[3])
    end
    open(joinpath(w.out, "report.txt"), "w") do io
        for r in w.results; println(io, r[2] ? "PASS" : "FAIL", "\t", r[1], "\t", r[3]); end
        foreach(l -> println(io, "WARN\t", l), warns)
    end
    return all(r -> r[2], w.results) && isempty(warns)
end

end
