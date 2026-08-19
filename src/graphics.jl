# Graphics-environment setup.
#
# This lives in the core package, not in the GUI extension, because of a hard ordering
# constraint: Mesa reads its driver variables when it dlopens the driver, at first GL context
# creation, so `configure_graphics!` has to run BEFORE `using GLMakie`. Anything defined in the
# GUI extension only exists AFTER that, since loading GLMakie is what triggers the extension.
#
# `prefer_native_wayland!` has the same ordering constraint but needs `GLFW_jll` for its ccall,
# so it lives in `ext/OITOOLSGLFWExt.jl` — a separate, one-function extension triggered by
# GLFW_jll alone, which a launcher can therefore load on its own.
#
# Exactly one platform needs help: WSL. There is no /dev/dri render node there — the GPU is
# reached through /dev/dxg and D3D12 — so Mesa's usual probe fails, falls through to zink,
# finds no Vulkan driver installed, and quietly lands on llvmpipe. Nothing errors out. The
# only symptom is that Makie feels slow, which reads as Makie's fault rather than the
# driver's. Pointing Mesa at its d3d12 gallium driver gets the real GPU back.
#
# Mesa reads these variables with getenv() when it dlopens the driver, and that happens at
# first GL context creation. Julia's `ENV[...] =` goes through setenv(), so C libraries do
# see the change — provided we get there first. Hence a function the launcher calls
# explicitly, above `using GLMakie`, rather than anything in `__init__`: a module's
# `__init__` runs whenever the package is loaded, which may be long after the caller opened a
# window, and a library that rewrites the process environment behind your back is no fun to
# debug either.
#
# Every check below is a guard that switches the workaround *off*: not WSL, a real /dev/dri,
# no /dev/dxg, no d3d12 driver on disk, or any of these variables set by hand — each of those
# means do nothing. If WSL ever grows a render node, this file stops acting on its own.

"Environment variable that disables everything in this file."
const NO_GPU_SETUP = "OITOOLSGUI_NO_GPU_SETUP"

"An environment variable is 'set' only if it says something other than off."
_isset(name) = !(lowercase(strip(get(ENV, name, ""))) in ("", "0", "false", "no"))

"True if this process is running under WSL."
function is_wsl()
    Sys.islinux() || return false
    haskey(ENV, "WSL_DISTRO_NAME") && return true
    isfile("/proc/version") || return false
    return occursin("microsoft", lowercase(read("/proc/version", String)))
end

"Directories Mesa searches for gallium DRI drivers. The layout differs per distribution."
function dri_driver_dirs()
    dirs = String[]
    haskey(ENV, "LIBGL_DRIVERS_PATH") &&
        append!(dirs, split(ENV["LIBGL_DRIVERS_PATH"], ':'; keepempty = false))
    append!(dirs, ["/usr/lib/dri", "/usr/lib64/dri", "/usr/lib/x86_64-linux-gnu/dri"])
    return dirs
end

"Path to Mesa's D3D12 gallium driver, or `nothing` if this Mesa was not built with it."
function d3d12_driver()
    for d in dri_driver_dirs()
        p = joinpath(d, "d3d12_dri.so")
        isfile(p) && return p
    end
    return nothing
end

"""
True if WSL is exposing an NVIDIA GPU.

File existence rather than `nvidia-smi`: this runs on every start, and the subprocess would
cost more than everything else in this file put together — and can block on a busy GPU.
"""
has_wsl_nvidia() = isfile("/usr/lib/wsl/lib/libnvidia-ml.so.1")

"""
    configure_graphics!(; verbose = true) -> NamedTuple

Give Mesa the one hint it cannot work out for itself on WSL, and do nothing anywhere else.

Returns `(applied, reason, vars)`: whether anything was set, the one-line explanation of why
or why not, and the variables set. Call it **before** the first OpenGL context exists —
`bin/oitoolsgui.jl` calls it above `using GLMakie`, and [`gui`](@ref) calls it again as a
backstop for sessions started by hand. The second call is a no-op, since by then the
variables are set and set variables are left alone.

Nothing here is forced on you. Any of

    OITOOLSGUI_NO_GPU_SETUP=1     # skip all of it
    GALLIUM_DRIVER=llvmpipe       # or any other explicit choice
    LIBGL_ALWAYS_SOFTWARE=1       # compare against software rendering

wins, which is what you want when the question is whether a bug is the driver's fault.

To see which renderer you actually ended up on, start the GUI with `QSG_INFO=1`; Qt reports
the GL vendor and renderer as it builds the scene graph.
"""
function configure_graphics!(; verbose::Bool = true)
    vars = Dict{String,String}()
    no(reason) = (applied = false, reason = reason, vars = vars)

    _isset(NO_GPU_SETUP)         && return no("disabled by $NO_GPU_SETUP")
    is_wsl()                     || return no("not WSL; the platform picks its own driver")
    isdir("/dev/dri")            && return no("/dev/dri exists; Mesa can find the GPU itself")
    ispath("/dev/dxg")           || return no("no /dev/dxg; this WSL exposes no GPU")
    _isset("LIBGL_ALWAYS_SOFTWARE") &&
        return no("LIBGL_ALWAYS_SOFTWARE is set; leaving software rendering alone")
    haskey(ENV, "GALLIUM_DRIVER") &&
        return no("GALLIUM_DRIVER already set to $(ENV["GALLIUM_DRIVER"])")

    driver = d3d12_driver()
    if driver === nothing
        verbose && @warn """
            No d3d12_dri.so in $(join(dri_driver_dirs(), ", ")).

            On WSL that leaves OpenGL on llvmpipe — software rendering on the CPU. The GUI
            will work and will be slow. Installing a Mesa built with the d3d12 gallium
            driver (Arch: `mesa`) is what makes the GPU reachable.
            """
        return no("no d3d12_dri.so found; staying on software rendering")
    end

    vars["GALLIUM_DRIVER"] = "d3d12"
    # The default D3D12 adapter is whichever Windows hands over first, which on a hybrid
    # laptop is the integrated GPU. The name is matched as a substring of the adapter
    # description, and is only worth setting once we know a discrete NVIDIA card is there.
    has_wsl_nvidia() && !haskey(ENV, "MESA_D3D12_DEFAULT_ADAPTER_NAME") &&
        (vars["MESA_D3D12_DEFAULT_ADAPTER_NAME"] = "NVIDIA")
    merge!(ENV, vars)

    verbose && @info "WSL detected: routing OpenGL through D3D12 instead of llvmpipe" driver vars
    return (applied = true, reason = "WSL with a D3D12 GPU", vars = vars)
end
