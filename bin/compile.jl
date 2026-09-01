# Build a sysimage for the GUI.
#
#     xvfb-run -a julia --project=bin bin/compile.jl        # headless
#     julia --project=bin bin/compile.jl                    # on a machine with a screen
#
# then launch with
#
#     julia --project=bin --sysimage bin/oitoolsgui.so bin/oitoolsgui.jl [file.oifits]
#
# WHAT THIS BUYS. A sysimage holds the packages already loaded and the methods already
# compiled, so the time between the command and the window is the time to map a file rather
# than to compile Makie, Qt and the reconstruction driver. It does NOT make the package
# faster once running.
#
# WHAT IT COSTS. The image is roughly a gigabyte and is tied to this machine's Julia version
# and this project's Manifest. Rebuild it after changing either, or after changing OITOOLS
# itself — a stale sysimage silently runs the OLD code, which is the one failure mode of this
# whole approach worth being afraid of. `bin/oitoolsgui.jl` cannot detect that for you.
#
# NO PYTHON IN THE IMAGE. The package list below names neither PythonPlot nor PythonCall, and
# `bin/trace.jl` loads neither, so the image carries no conda environment. Nested sampling
# comes from Nautilus.jl and plotting from Makie; `plot_ultranest_corner` and the
# matplotlib figures are simply absent from a session started on this image, and calling one
# raises a MethodError naming the function. Load PythonPlot in that session to get them back —
# the sysimage does not prevent it, it just does not carry it.
#
# A DISPLAY IS NEEDED to build, because the trace opens a GL context. Xvfb is enough.
#
# ONE ORDERING CAVEAT, handled for you. `bin/oitoolsgui.jl` decides GLMakie's windowing
# platform before Qt's, because only GLMakie's can fail: GLFW.jl initialises from its own
# `__init__`, which in a sysimage runs before any user code, so `JULIA_GLFW_PLATFORM` arrives
# too late and GLFW stays on X11. The launcher notices and pins Qt to xcb to match, so the two
# halves agree either way. On a Wayland session that means a sysimage run uses XWayland where
# the plain launcher uses Wayland.
#
# To get native Wayland out of a sysimage, export the variable before Julia starts:
#
#     JULIA_GLFW_PLATFORM=wayland julia --project=bin --sysimage bin/oitoolsgui.so bin/oitoolsgui.jl
#
# A DISPLAY IS NEEDED to build, because the trace opens a GL context. Xvfb is enough.

using Pkg
using PackageCompiler

const HERE    = @__DIR__
const SYSIMG  = get(ENV, "OITOOLS_SYSIMAGE", joinpath(HERE, "oitoolsgui.so"))
const TRACE   = joinpath(HERE, "trace.jl")

# Named explicitly rather than taken from the project, which also holds PythonPlot, PythonCall
# and PackageCompiler itself — none of which belong in the image.
# Makie is NOT listed: PackageCompiler requires every named package to be a direct dependency
# of the project, and `bin/` reaches Makie through GLMakie. `include_transitive_dependencies`
# defaults to true, so it lands in the image regardless.
const PACKAGES = [:OITOOLS, :GLMakie, :QMLMakie, :QML, :GLFW_jll,
                  :Nautilus, :PairPlots]

isfile(TRACE) || error("missing $TRACE")

@info "Building sysimage" sysimage = SYSIMG packages = PACKAGES
@info "This takes tens of minutes and several GB of scratch space."

t = @elapsed create_sysimage(PACKAGES;
                             sysimage_path = SYSIMG,
                             precompile_execution_file = TRACE,
                             incremental = true)

@info "Done" seconds = round(t, digits = 1) size_MB = round(filesize(SYSIMG) / 1e6, digits = 1)
println("""

Launch with:

    julia --project=$(relpath(HERE)) --sysimage $(relpath(SYSIMG)) $(relpath(joinpath(HERE, "oitoolsgui.jl"))) [file.oifits]

Rebuild it whenever OITOOLS, the Manifest or Julia itself changes: a stale sysimage runs the
code it was built from, not the code on disk.
""")
