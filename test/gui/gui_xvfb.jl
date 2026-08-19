# Automated GUI smoke test. Needs a display -- any display.
#
#     # on a machine with a real screen (this is the one that matters):
#     julia --project=bin test/gui_xvfb.jl
#
#     # headless, e.g. the dev box or CI:
#     xvfb-run -a --server-args="-screen 0 1400x1000x24" julia --project=bin test/gui_xvfb.jl
#
# Run it BOTH ways. Under Xvfb the GL driver is llvmpipe, which is forgiving; a real driver is
# not necessarily. The delete-and-recreate redraw below passes on llvmpipe and fails on WSL's
# D3D12 with "Failed to resolve gl_renderobject", so a green run here does NOT mean the GUI
# works on your machine -- which is exactly why each step is named and reported separately.
#
# WHY. Everything else in test/ runs without a display, which means it cannot touch the code
# path where the GUI actually breaks. A plot drawn in `gui()` BEFORE `loadqml` succeeds even
# when every subsequent one fails, so "the first plot works and nothing after it does" is
# invisible to a headless test by construction. This drives the window after it is open, from
# Qt's own GUI thread, via the `on_ready` hook.
#
# It asserts on Julia-side state rather than pixels: the console transcript, the session, and
# the number of points the axis holds. A screenshot diff would fail for font reasons and would
# not say which action broke.

# OITOOLSGUI first, and prefer_native_wayland! before GLMakie: GLFW.jl hard-codes X11, so
# without this the test would run Qt on Wayland and GLMakie on XWayland under a Wayland
# compositor -- two windowing systems, which is not what it is meant to be checking.
using OITOOLS, GLMakie, QMLMakie, QML, Test

# The GUI is an extension, so its types are not `OITOOLS.Session`; reach them the way
# test/gui/runtests.jl does.
const GUI = Base.get_extension(OITOOLS, :OITOOLSGUIExt)
using .GUI

const DATA = joinpath(@__DIR__, "data")
const MONO = joinpath(DATA, "2004-data1.oifits")   # nwav=1: the degenerate colour case
const POLY = joinpath(DATA, "v1295Aql.oifits")     # nwav=5, 4932 uv points

# Each step records (name, ok, detail). Steps run in order and never throw out of the hook,
# so one failure still leaves a full transcript rather than a dead window.
const STEPS = Tuple{String,Bool,String}[]

# Function first: `step("name") do ... end` passes the closure as argument ONE.
function step(f, name)
    try
        f()
        push!(STEPS, (name, true, ""))
    catch err
        push!(STEPS, (name, false, first(split(sprint(showerror, err), "\n"))))
    end
end

"Number of points currently on the axis — the check that a redraw actually redrew."
function npoints(sh)
    n = 0
    for p in sh.axis.scene.plots
        try
            v = p[1][]
            v isa AbstractVector && (n += length(v))
        catch
        end
    end
    return n
end

function scripted()
    sh = GUI.SHELL[]

    # 1. The first plot after the window exists. This is the one that fails today.
    step("open mono after window is up") do
        s = GUI.shell_open(MONO)
        occursin("could not load", s) && error(s)
        npoints(sh) > 0 || error("axis is empty after load")
    end

    # 2. Redraw with the same kind: pure replot, no new observable type.
    step("replot same view") do
        s = GUI.shell_set_view("uv", "baseline")
        occursin("cannot plot", s) && error(s)
    end

    # 3. Colour mode change — degenerate here, since mono has one wavelength.
    step("mono + colour by wav (degenerate)") do
        s = GUI.shell_set_view("uv", "wav")
        occursin("cannot plot", s) && error(s)
    end

    # 4. A second file, which is the other half of the bug report.
    step("open poly as a second dataset") do
        s = GUI.shell_open(POLY)
        occursin("could not load", s) && error(s)
        length(sh.session.datasets) == 2 || error("expected 2 datasets")
    end

    # 5. Colour by wavelength for real, and by MJD.
    step("poly + colour by wav") do
        s = GUI.shell_set_view("uv", "wav")
        occursin("cannot plot", s) && error(s)
    end
    step("poly + colour by mjd") do
        s = GUI.shell_set_view("uv", "mjd")
        occursin("cannot plot", s) && error(s)
    end

    # 6. Kind changes: different plot types, error bars, log scale.
    for kind in ("v2", "t3phi", "t3amp", "uv")
        step("switch to $kind") do
            s = GUI.shell_set_view(kind, "baseline")
            occursin("cannot plot", s) && error(s)
            npoints(sh) > 0 || error("axis empty after switching to $kind")
        end
    end

    # 7. Back to the first dataset.
    step("select dataset 1") do
        s = GUI.shell_select_dataset(1)
        occursin("cannot plot", s) && error(s)
    end
end

session = Session()
gui(session; autoquit_ms = 9000, on_ready = scripted)

# ── report ───────────────────────────────────────────────────────────────────
# The console is the window's own record of what happened; print it whenever a run produces
# no steps, because that means the hook never ran and the steps cannot say why.
let sh = GUI.SHELL[]
    if sh !== nothing && (isempty(STEPS) || any(s -> !s[2], STEPS))
        println("\n── console transcript ──")
        for l in sh.console; println("  ", l); end
        println("  READY_HOOK installed: ", GUI.READY_HOOK[] !== nothing)
    end
end
println("\n── GUI smoke test ──")
for (name, ok, detail) in STEPS
    println("  ", ok ? "PASS  " : "FAIL  ", rpad(name, 38), detail)
end
nfail = count(s -> !s[2], STEPS)
println(nfail == 0 ? "\nall $(length(STEPS)) steps passed" : "\n$nfail of $(length(STEPS)) steps FAILED")

@testset "GUI on a virtual screen" begin
    @test !isempty(STEPS)                       # the hook ran at all
    for (name, ok, detail) in STEPS
        @testset "$name" begin
            @test ok || error(detail)
        end
    end
end
exit(nfail == 0 ? 0 : 1)
