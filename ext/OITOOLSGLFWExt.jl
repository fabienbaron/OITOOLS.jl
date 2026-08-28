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
    # An explicit request for XWayland. `configure_qt_platform!` no longer pins xcb by itself,
    # so this only fires when the caller set the variable — and then GLFW's X11 default is
    # already what is wanted.
    lowercase(get(ENV, "QT_QPA_PLATFORM", "")) == "xcb" &&
        return no("Qt is on xcb; GLFW's X11 default already matches it")

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

    # GLFW.jl's own override, rather than a race against it.
    #
    # `GLFW.Init` hard-codes `PLATFORM_X11` on Linux (glfw3.jl:499) -- which is why GLMakie
    # lands on XWayland under a Wayland compositor -- but two lines later it honours
    # `JULIA_GLFW_PLATFORM`, and `HandlePlatformSelection` maps "wayland" to PLATFORM_WAYLAND.
    # Setting that is the supported route and needs no ccall.
    #
    # This function used to call `glfwInitHint` + `glfwInit` here instead, getting in first so
    # that GLFW.jl's later `Init` became a no-op. That works only if nothing has initialised
    # GLFW yet, and in a PackageCompiler sysimage carrying GLMakie something has: GLFW.jl's
    # `__init__` runs at image load, before any user code. The variable has no such ordering
    # problem — exported from a shell it is set before Julia starts at all.
    already = ccall((:glfwGetPlatform, GLFW_jll.libglfw), Cint, ())
    if already != Cint(0)
        name = already == Cint(0x00060004) ? "X11" :
               already == Cint(0x00060003) ? "Wayland" :
               already == Cint(0x00060002) ? "Cocoa" :
               already == Cint(0x00060001) ? "Win32" :
               already == Cint(0x00060005) ? "null" : "0x$(string(already, base = 16))"
        already == Cint(0x00060003) &&
            return (applied = false, reason = "GLFW is already on Wayland; nothing to do")
        verbose && @warn """
            GLFW is already initialised on $name, so it is too late to move it to Wayland.

            GLMakie will run on $name while Qt may not, which is the mismatch this function
            exists to avoid. The usual cause is a sysimage built with GLMakie in it: its
            `__init__` runs before any user code, and `glfwInit` cannot be undone.

            Export the variable before starting Julia instead:

                JULIA_GLFW_PLATFORM=wayland julia --sysimage ... bin/oitoolsgui.jl

            The launcher handles this on its own: it passes `match_x11 = true` to
            `configure_qt_platform!` when this returns not-applied, so Qt follows GLMakie to
            xcb rather than the two ending up split.
            """
        return no("GLFW already initialised on $name; too late to choose a platform")
    end

    ENV["JULIA_GLFW_PLATFORM"] = "wayland"
    verbose && @info "GLFW will initialise on native Wayland; GLMakie will not use XWayland"
    return (applied = true, reason = "JULIA_GLFW_PLATFORM=wayland set before GLFW.jl loads")
end

end # module
