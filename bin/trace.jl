# What `bin/compile.jl` traces to build the sysimage.
#
# PackageCompiler records every method this script compiles and bakes those specialisations
# into the image, so what belongs here is exactly the work a user waits for on a cold start:
# loading the stack, building the figures `gui()` builds, and drawing into them once.
#
# NOT `gui()` itself. That opens a Qt window and hands control to its event loop, which would
# never return and would hang the build. Everything below the window is reachable without it.
#
# NO PYTHON. PythonPlot and PythonCall are dependencies of `bin/` — the port harness needs the
# first, the UltraNest backend the second — but neither is loaded here and neither is passed to
# `create_sysimage`. A compiled build that dragged in a conda environment would defeat the
# point: nested sampling comes from Nautilus.jl and plotting from Makie.
#
# Needs a display for the GL context. `bin/compile.jl` supplies one.

using OITOOLS

configure_graphics!()

using GLFW_jll
# Same order as bin/oitoolsgui.jl: GLMakie's platform first, then Qt follows it.
configure_qt_platform!(; match_x11 = !prefer_native_wayland!().applied)

using GLMakie, QMLMakie, QML
using Nautilus, PairPlots
using Random: Xoshiro

const GUI = Base.get_extension(OITOOLS, :OITOOLSGUIExt)
const MK  = Base.get_extension(OITOOLS, :OITOOLSMakieExt)
(GUI === nothing || MK === nothing) && error("extensions did not load; nothing to trace")
using .GUI

const FILE = joinpath(pkgdir(OITOOLS), "test", "oifits_for_tests", "2004-data1.oifits")
const OUT  = joinpath(mktempdir(), "trace.png")

"""
    traced(f, what)

Run one section of the trace, reporting rather than throwing if it fails.

A sysimage build is twelve minutes and everything before the failure is discarded, so a single
unlucky step must not destroy it. Used for the STOCHASTIC section only — a deterministic step
that fails is a real fault and should stop the build. The cost of catching is a thinner image, so the failure is
reported loudly: a warning here means whatever `what` names will be compiled in the user's
first session instead.
"""
function traced(f, what)
    try
        f()
    catch err
        @warn "trace section failed; its methods will NOT be in the sysimage" section = what exception = err
    end
    return nothing
end

# Both shapes: `setup_ft` and `reconstruct` take the whole `(nwav, ntime)` array `readoifits`
# returns, while every plotting function takes one cell of it.
cube = readoifits(FILE; warn = false, verbose = false)
data = cube[1, 1]

# ── The plots, saved rather than merely built ────────────────────────────────
#
# `Makie.save` is what forces the render, and rendering is where most of the compilation is.
# Building a Figure and stopping leaves the expensive half untraced.
for p in (uvplot_makie(data),
          plot_v2_makie(data),
          plot_v2_makie(data; logscale = true),
          plot_t3phi_makie(data),
          plot_obs_makie(data),
          imdisp_makie(rand(32, 32); pixsize = 0.2, colorbar = true))
    Makie.save(OUT, p.figure)
end

# ── The GUI's own panels ─────────────────────────────────────────────────────
#
# Each is built before the window exists, so its compilation is time the user spends looking at
# nothing. `build_chi2_map` alone was measured at 1.72 s on a cold session — the single most
# expensive step in `gui()`, because `contour!` is a recipe compiled nowhere else.
fig = Makie.Figure(fonts = OITOOLS.PLOT_FONTS)
ax  = Makie.Axis(fig[1, 1])
OITOOLS.style_axis!(ax)
canvas = build_canvas(fig, ax)
update_canvas!(canvas, data, :uv; color = :baseline)
update_canvas!(canvas, data, :v2; color = :wav)
update_canvas!(canvas, data, :t3phi; color = :mjd)
GUI.show_image!(canvas, zeros(Float32, 16, 16), 0.3)
GUI.hide_image!(canvas)

let f = Makie.Figure(), a = Makie.Axis(f[1, 1])
    build_delay_plot(f, a)
end
let f = Makie.Figure(), a = Makie.Axis(f[1, 1])
    build_gantt(f, a)
end

# ── Model fitting, and the χ² map the Modeling tab draws ─────────────────────
md   = Dict{String,Any}("s,ud" => 6.0, "s,f" => 1.0)
free = ["s,ud"]
fit_model(md, free, data; lb = Dict("s,ud" => 1.0), ub = Dict("s,ud" => 12.0), verb = false)
fit_model_lsqfit(md, free, data; lb = Dict("s,ud" => 1.0), ub = Dict("s,ud" => 12.0), verb = false)

let f = Makie.Figure(), a = Makie.Axis(f[1, 1])
    c = build_chi2_map(f, a)
    # `ldpow` + `alpha`, not `ud` + `alpha`: a uniform disc has no limb-darkening exponent,
    # and the parser warns that the key is unused rather than erroring — exactly the silent
    # misparse a trace file should not be teaching.
    m = chi2_map(Dict{String,Any}("s,ldpow" => 6.0, "s,alpha" => 0.2, "s,f" => 1.0),
                 ["s,ldpow", "s,alpha"], data, "s,ldpow", "s,alpha";
                 lb = Dict("s,ldpow" => 5.0, "s,alpha" => 0.0),
                 ub = Dict("s,ldpow" => 8.0, "s,alpha" => 0.5), n1 = 6, n2 = 6)
    update_chi2_map!(c, m)
end

# Residuals and a corner plot: the Modeling tab's two output figures.
let fm = dict_to_model(md, free)
    Makie.save(OUT, plot_residuals_makie(fm, [6.8], data).figure)
end

# Nested sampling, the pure-Julia backend — the one a compiled build actually has.
#
# Seeded, and with enough live points to be stable. At `nactive = 60` on a one-parameter fit
# the bounding ellipsoid collapses onto the peak and rejection sampling exhausts its 100 000
# draws without landing above the likelihood threshold — which killed a whole build. That is
# a property of rejection sampling at tight settings, not of the batching: stock
# `Proposals.Rejection` throws the same way. A trace is not the place to run a sampler at its
# limit.
traced("nested sampling + corner plot") do
    r = fit_model_nested(md, free, data; backend = :nestedsamplers,
                         lb = Dict("s,ud" => 5.0), ub = Dict("s,ud" => 9.0),
                         nactive = 200, rng = Xoshiro(1234),
                         verb = false, cornerplot = false)
    Makie.save(OUT, plot_corner_makie(r.posterior, r.list_free_params).figure)
end

# ── Imaging ──────────────────────────────────────────────────────────────────
#
# `reconstruct` is the single most expensive thing to compile in the package: measured from the
# GUI, the first Run spent one to two MINUTES compiling before ten seconds of actual work. Two
# iterations are enough, because the cost is compilation and not iteration. The regularised
# criterion is a separate specialisation, so both are traced.
let nx = 32, pixsize = 0.5
    ft = setup_ft(cube, nx, pixsize)
    x0 = gaussian2d(nx, nx, nx ÷ 6)
    reconstruct(x0, cube, ft; maxiter = 2, verb = false)
    reconstruct(x0, cube, ft; maxiter = 2, verb = false,
                regularizers = [["centering", 1e4], ["tv", 1e2]])
end

@info "trace complete"
