# Graphics-environment setup.
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

# ── native Wayland for GLMakie ───────────────────────────────────────────────

"Environment variable that keeps GLMakie on XWayland (the GLFW.jl default)."
const FORCE_X11 = "OITOOLSGUI_GLFW_X11"

"""
    prefer_native_wayland!(; verbose = true) -> NamedTuple

On a Wayland session, make GLMakie use Wayland rather than XWayland.

GLFW.jl hard-codes X11 on Linux -- `Init(; platform = Sys.islinux() ? PLATFORM_X11 : ...)` --
so under a Wayland compositor `using GLMakie` either fails with "X11: The DISPLAY environment
variable is missing", or, if DISPLAY happens to be set, silently runs on **XWayland**.

The second case is the one that matters here. Qt picks native Wayland while GLMakie is on
XWayland, so the application is running two windowing systems at once, with two EGL display
connections and two surface lifetimes, and Qt Quick popups become `xdg_popup` surfaces whose
create/destroy sequence differs from the X11 child-window equivalent. That is a bad place to
be for a GUI whose plot surface already had one GL-context bug.

`glfwInit` is idempotent, so initialising it here with the Wayland hint makes GLFW.jl's later
`Init()` a no-op and its X11 hint irrelevant. **Must run before `using GLMakie`.**

Set `\$FORCE_X11=1` to skip it and get the XWayland behaviour back -- the A/B switch for
deciding whether a bug is the split's fault.
"""
function prefer_native_wayland!(; verbose::Bool = true)
    no(reason) = (applied = false, reason = reason)
    Sys.islinux()             || return no("not Linux")
    haskey(ENV, "WAYLAND_DISPLAY") || return no("not a Wayland session; GLFW's X11 default is right")
    _isset(FORCE_X11)         && return no("$FORCE_X11 is set; leaving GLMakie on XWayland")

    # GLFW_jll is a direct dependency: loading it only dlopens libglfw, which needs no display
    # (glfwInit is the part that does), so this stays safe on a headless machine. The library
    # is named directly in the ccall because a ccall cannot take it from a local variable.
    # Silence xkbcommon's Compose-file complaint, which this function causes.
    #
    #   xkbcommon: ERROR: [XKB-679] No Compose file for locale "en_US.UTF-8"
    #
    # GLFW's Wayland path uses xkbcommon for keyboard handling; its X11 path does not. So the
    # error appears the moment we switch to Wayland and never before, which makes it look like
    # a new fault rather than what it is: xkbcommon_jll's artifact ships no X11 locale data, so
    # it cannot find a Compose table. Keyboard input works regardless -- Compose sequences are
    # for typing accented characters, which nothing here does. Point it at the system table
    # when there is one, and leave any existing setting alone.
    if !haskey(ENV, "XCOMPOSEFILE")
        for f in ("/usr/share/X11/locale/en_US.UTF-8/Compose", "/usr/share/X11/locale/C/Compose")
            if isfile(f)
                ENV["XCOMPOSEFILE"] = f
                break
            end
        end
    end

    PLATFORM = Cint(0x00050003)   # GLFW_PLATFORM
    WAYLAND  = Cint(0x00060003)   # GLFW_PLATFORM_WAYLAND
    ccall((:glfwInitHint, GLFW_jll.libglfw), Cvoid, (Cint, Cint), PLATFORM, WAYLAND)
    if ccall((:glfwInit, GLFW_jll.libglfw), Cint, ()) == 1
        verbose && @info "GLFW pre-initialised on native Wayland; GLMakie will not use XWayland"
        return (applied = true, reason = "GLFW forced to native Wayland")
    end
    verbose && @warn "forced-Wayland glfwInit failed; GLFW.jl will fall back to X11"
    return no("forced-Wayland glfwInit failed")
end
