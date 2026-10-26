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
    configure_qt_platform!(; match_x11 = false, verbose = true) -> (; applied, reason, vars)

Keep Qt on the same windowing system as GLMakie.

Qt follows the session (Wayland or X11) on its own. With `match_x11 = true`, set when GLMakie
could not use Wayland, this sets `QT_QPA_PLATFORM=xcb` so that Qt uses X11 too.
[`oitoolsgui`](@ref) passes `match_x11 = !prefer_native_wayland!().applied`, which is true when
`\$FORCE_X11` is set or GLFW was already initialised (as in a PackageCompiler sysimage). An
explicit `QT_QPA_PLATFORM` is left alone.
"""
function configure_qt_platform!(; match_x11::Bool = false, verbose::Bool = true)
    vars = Dict{String,String}()
    no(reason) = (applied = false, reason = reason, vars = vars)

    haskey(ENV, "QT_QPA_PLATFORM") &&
        return no("QT_QPA_PLATFORM already set to $(ENV["QT_QPA_PLATFORM"])")
    haskey(ENV, "WAYLAND_DISPLAY") ||
        return no("not a Wayland session; Qt picks its own platform")

    match_x11 || return no("Wayland session; Qt and GLMakie both on Wayland")

    vars["QT_QPA_PLATFORM"] = "xcb"
    merge!(ENV, vars)
    verbose && @info "GLMakie is on X11, so Qt is pinned to xcb to match it"
    return (applied = true, reason = "Qt pinned to xcb to match GLMakie on X11", vars = vars)
end

"""
    configure_graphics!(; verbose = true) -> NamedTuple

Set the Mesa driver variables needed on WSL; elsewhere, do nothing.

Returns `(applied, reason, vars)`: whether anything was set, a one-line reason, and the
variables set. Call it before the first OpenGL context exists: [`oitoolsgui`](@ref) calls
it before loading GLMakie, and [`gui`](@ref) calls it again. Variables already set are left alone,
so any of these takes precedence:

    OITOOLSGUI_NO_GPU_SETUP=1     # skip all of it
    GALLIUM_DRIVER=llvmpipe       # or any other explicit choice
    LIBGL_ALWAYS_SOFTWARE=1       # compare against software rendering


With `QSG_INFO=1` Qt prints the GL vendor and renderer.
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

import TOML

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
    edid_native_mode(edid) -> (; width, height, width_mm, height_mm) or nothing

The native mode of a display and its image size, from the first detailed timing descriptor
of its EDID (bytes 54-71). `nothing` when the block is too short or the descriptor is not a
timing (pixel clock 0) or gives no size.
"""
function edid_native_mode(edid::AbstractVector{UInt8})
    length(edid) >= 72 || return nothing
    d = @view edid[55:72]
    (d[1] == 0 && d[2] == 0) && return nothing
    w  = Int(d[3])  | (Int(d[5] >> 4) << 8)
    h  = Int(d[6])  | (Int(d[8] >> 4) << 8)
    wm = Int(d[13]) | (Int(d[15] >> 4) << 8)
    hm = Int(d[14]) | (Int(d[15] & 0x0f) << 8)
    (w > 0 && h > 0 && wm > 0 && hm > 0) || return nothing
    return (; width = w, height = h, width_mm = wm, height_mm = hm)
end

"""
    panel_dpis(; root = "/sys/class/drm") -> Vector{(; width, height, dpi)}

Every connected display with a readable EDID under `root`: its native mode in device pixels
and its density in device pixels per inch. Empty off Linux.
"""
function panel_dpis(; root::AbstractString = "/sys/class/drm")
    out = @NamedTuple{width::Int, height::Int, dpi::Float64}[]
    Sys.islinux() && isdir(root) || return out
    for c in readdir(root; join = true)
        st, ed = joinpath(c, "status"), joinpath(c, "edid")
        isfile(st) && isfile(ed) || continue
        try
            strip(read(st, String)) == "connected" || continue
            m = edid_native_mode(read(ed))
            m === nothing && continue
            push!(out, (; width = m.width, height = m.height, dpi = m.width / (m.width_mm / 25.4)))
        catch err
            @debug "unreadable display connector" c err
        end
    end
    return out
end

# One density for a set of displays: theirs when they agree to within 1 dpi, else 0.
_one_dpi(ps) = isempty(ps) ? 0.0 :
               maximum(p -> p.dpi, ps) - minimum(p -> p.dpi, ps) < 1 ? first(ps).dpi : 0.0

"""
    panel_dpi(width, height; root = "/sys/class/drm") -> Float64

Density, in device pixels per inch, of the connected display whose native mode is `width` x
`height` device pixels (to 2 pixels, for a size rounded through a fractional scale), from its
EDID. Zero when no connected display has that mode, when the
displays that do disagree, or off Linux.

Under X11 Qt takes the screen's size in millimetres from the X server, which can report a
nominal figure (a configured DPI, XWayland's default) rather than the panel's.
"""
panel_dpi(width::Real, height::Real; root::AbstractString = "/sys/class/drm") =
    _one_dpi(filter(p -> abs(p.width - width) <= 2 && abs(p.height - height) <= 2,
                    panel_dpis(; root)))

"Variables through which Qt is told its scale; any of them set means it is not ours to choose."
const QT_SCALE_VARS = ("QT_SCALE_FACTOR", "QT_SCREEN_SCALE_FACTORS",
                       "QT_ENABLE_HIGHDPI_SCALING", "QT_AUTO_SCREEN_SCALE_FACTOR")

"""
    window_scale_for(dpi) -> Float64

The whole-window scale a panel of `dpi` device pixels per inch gets: `dpi / 150` to the nearest
0.25, and never below 1. 1.25 for a 189 dpi 16" 2560x1600 panel, 1 for a 141 dpi 15.6"
1920x1080 one and for a desktop monitor.
"""
window_scale_for(dpi::Real) = max(1.0, round(4 * dpi / 150) / 4)

# Xft.dpi / 96: the scale Qt on X11 applies by itself. 1 when unset or unreadable.
function _x11_desktop_scale()
    xrdb = Sys.which("xrdb")
    xrdb === nothing && return 1.0
    try
        for l in eachline(pipeline(`$xrdb -query`; stderr = devnull))
            m = match(r"^Xft\.dpi:\s*([0-9.]+)", l)
            m === nothing && continue
            v = tryparse(Float64, m[1])
            return v === nothing || v <= 0 ? 1.0 : v / 96
        end
    catch err
        @debug "xrdb -query failed" err
    end
    return 1.0
end

# Windows' display scale (its dpi / 96) as read under WSL from the registry through `reg.exe`;
# 0 when it cannot be read. WSLg gives Qt a scale of 1 whatever Windows uses.
function _windows_desktop_scale()
    reg = Sys.which("reg.exe")
    reg === nothing && isfile("/mnt/c/Windows/System32/reg.exe") &&
        (reg = "/mnt/c/Windows/System32/reg.exe")
    reg === nothing && return 0.0
    for (key, name) in (("HKCU\\Control Panel\\Desktop\\WindowMetrics", "AppliedDPI"),
                        ("HKCU\\Control Panel\\Desktop", "LogPixels"))
        try
            out = read(pipeline(`$reg query $key /v $name`; stderr = devnull), String)
            m = match(r"REG_DWORD\s+0x([0-9a-fA-F]+)", out)
            m === nothing && continue
            v = parse(Int, m[1]; base = 16)
            v > 0 && return v / 96
        catch err
            @debug "reg.exe query failed" key name err
        end
    end
    return 0.0
end

"""
    wsl_panel_dpi() -> Float64

Under WSL, the panel density Windows' display scale stands for: 120 dpi per 100%, the step
Windows recommends, so 180 at 150%. 0 outside WSL or when the registry cannot be read. WSLg
reports a nominal 96-100 dpi for any display and exposes no EDID.
"""
function wsl_panel_dpi()
    is_wsl() || return 0.0
    return 120 * _windows_desktop_scale()
end

"""
    configure_qt_scale!(; verbose = true) -> (; applied, reason, vars)

Scale all of Qt on Linux, through `QT_SCALE_FACTOR`, when the desktop does not.

Qt draws the controls' padding, borders and indicators at the desktop's scale. On X11 at 100%
a dense panel gets them at one device pixel per point while the window's own scale enlarges
the text, so text crowds its fields; Windows at 150% draws the same panel's controls at 1.5.
This sets the factor `window_scale_for` gives the panel's EDID density, divided by what X11
already applies (`Xft.dpi`). Under WSL, whose WSLg gives Qt a scale of 1, the density is
`wsl_panel_dpi`, from Windows' display scale (`AppliedDPI` in the registry).

`window_scale` in `gui_settings_file()` wins when positive (1 turns it off). Left alone: any
of `QT_SCALE_VARS` already set; a Wayland session other than WSLg, whose compositor gives Qt
its scale; WSL when the registry cannot be read; displays that disagree in density or have no
EDID. Must run before `QML` is loaded, which creates the `QGuiApplication`; `app/launcher.sh`
does the same for the bundle.
"""
function configure_qt_scale!(; verbose::Bool = true)
    vars = Dict{String,String}()
    no(reason) = (applied = false, reason = reason, vars = vars)
    Sys.islinux() || return no("not Linux")
    for v in QT_SCALE_VARS
        haskey(ENV, v) && return no("$v already set")
    end
    saved = 0.0
    try
        path = gui_settings_file()
        if isfile(path)
            v = get(TOML.parsefile(path), "window_scale", 0.0)
            v isa Real && isfinite(v) && (saved = Float64(v))
        end
    catch err
        @debug "could not read the saved window scale" err
    end
    f = if saved > 0
        clamp(saved, 0.5, 4.0)
    elseif is_wsl()
        dpi = wsl_panel_dpi()
        dpi > 0 || return no("WSL; Windows' display scale unreadable")
        window_scale_for(dpi)
    else
        x11 = lowercase(get(ENV, "QT_QPA_PLATFORM", "")) == "xcb" || !haskey(ENV, "WAYLAND_DISPLAY")
        x11 || return no("Wayland session; the compositor scales Qt")
        dpi = _one_dpi(panel_dpis())
        dpi > 0 || return no("no single panel density from the EDID")
        window_scale_for(dpi) / _x11_desktop_scale()
    end
    abs(f - 1) < 0.01 && return no("scale 1")
    vars["QT_SCALE_FACTOR"] = string(round(f; digits = 3))
    merge!(ENV, vars)
    verbose && @info "Qt scaled to match the panel" QT_SCALE_FACTOR = vars["QT_SCALE_FACTOR"]
    return (applied = true, reason = (saved > 0 ? "saved window_scale" :
                                  is_wsl() ? "Windows' display scale" : "from the panel's EDID"),
            vars = vars)
end
