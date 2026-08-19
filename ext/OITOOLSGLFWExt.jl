# `prefer_native_wayland!`, in its own extension.
#
# It is one function, and it is here rather than in `OITOOLSGUIExt` for an ordering reason:
# it must run BEFORE `using GLMakie`, so it cannot live in an extension that GLMakie triggers.
# GLFW_jll alone is enough for it, and loading GLFW_jll only dlopens libglfw — no display
# needed — so a launcher can pull this in, call the function, and only then load GLMakie:
#
#     using OITOOLS
#     configure_graphics!()        # core: no dependency at all
#     using GLFW_jll               # activates this extension
#     prefer_native_wayland!()
#     using GLMakie, QMLMakie, QML # activates OITOOLSGUIExt
#
# `configure_graphics!` and the rest of the graphics-environment logic are in
# `src/graphics.jl`, in the core package.
module OITOOLSGLFWExt

using OITOOLS
using OITOOLS: _isset
using GLFW_jll

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
function OITOOLS.prefer_native_wayland!(; verbose::Bool = true)
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

end # module
