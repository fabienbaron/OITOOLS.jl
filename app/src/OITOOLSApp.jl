# The application entry point.
#
# WHAT IS DIFFERENT ABOUT AN APPLICATION, and it is not obvious: a bundle is a sysimage, and
# every package `__init__` in a sysimage runs BEFORE any user code. Measured against
# `bin/oitoolsgui.so`, before a single line of the script had executed:
#
#     GLFW.INITIALIZED[]        : true
#     QML QGuiApplication built : true
#
# So two of the three ordering hooks that `bin/oitoolsgui.jl` relies on cannot work here:
# `prefer_native_wayland!` needs to run before `glfwInit`, and `configure_qt_platform!` before
# Qt starts, and both have already happened. Only `configure_graphics!` survives, because Mesa
# reads its variables when the GL CONTEXT is created and that happens later, inside `gui`.
#
# The consequence is the reason `launcher.sh` exists. Left alone on a Wayland session, GLFW
# takes its hardcoded X11 and Qt follows the session to Wayland -- the SPLIT configuration,
# two windowing systems and two EGL connections in one process, which is the state
# `prefer_native_wayland!` was written to prevent. The launcher sets the variables before the
# process starts, which is the only moment that still works.
#
# `oitoolsgui()` is deliberately not used. It offers to install missing optional packages, and
# an application cannot: its depot is read-only and its user never installed Julia. Every
# optional engine is an ordinary dependency in Project.toml instead, so nothing is greyed out
# and nothing is ever offered.

module OITOOLSApp

using OITOOLS
using GLFW_jll
using GLMakie, QMLMakie, QML
using Nautilus                       # "Nested sampling" in the Model panel
using Pigeons                        # "Tempering" in the Image panel
using VarInf                         # "Variational inference" in the Image panel

const GUI = Base.get_extension(OITOOLS, :OITOOLSGUIExt)

"""
    julia_main() -> Cint

The bundle's entry point.

Returns a status code and lets nothing escape: there is no REPL behind an application, so an
uncaught exception would take the process down with its message going nowhere a user can read.
Failures are written to [`crash_log_path`](@ref) as well as to stderr.
"""
Base.@ccallable function julia_main()::Cint
    try
        run_app(copy(ARGS))
        return 0
    catch err
        report_crash(err, catch_backtrace())
        return 1
    end
end

"""
    run_app(files; autoquit_ms = 0) -> Session

Open the window, with any files given on the command line already loaded.

Separated from [`julia_main`](@ref) so that it can be driven from a test without building a
bundle: `autoquit_ms` closes the window by itself, and `\$OITOOLS_AUTOQUIT_MS` does the same to
a bundle that is already built, which is how the bundle is checked in a container where there
is nobody to close a window.
"""
function run_app(files::AbstractVector{<:AbstractString} = String[]; autoquit_ms::Integer = 0)
    autoquit_ms = autoquit_ms > 0 ? autoquit_ms :
                  something(tryparse(Int, get(ENV, "OITOOLS_AUTOQUIT_MS", "")), 0)
    GUI === nothing && error("OITOOLSGUIExt did not load; GLMakie, QMLMakie and QML are all needed")

    # Still effective in an application: Mesa reads these when the GL context is created, which
    # has not happened yet. The two platform hooks are not, hence launcher.sh -- see the notes
    # at the top of this file.
    configure_graphics!()
    warn_if_split_windowing()
    use_writable_cache_dir!()
    register_qml_modules!()

    session = GUI.Session()
    for f in files
        if isfile(f)
            GUI.load_dataset!(session, f; warn = false, verbose = false)
        else
            @warn "not a file, skipping" file = f
        end
    end
    OITOOLS.gui(session; autoquit_ms)
    return session
end

# ── the two things a bundle needs and a checkout does not ────────────────────

"""
    warn_if_split_windowing()

Say so when GLFW and Qt have landed on different windowing systems.

It cannot be fixed from here -- both are already initialised by the time this runs -- so this
reports rather than repairs, and names the variable that would have prevented it. The state is
reachable only when the launcher was bypassed, which is exactly when nobody is expecting it.
"""
function warn_if_split_windowing()
    Sys.islinux() || return nothing
    haskey(ENV, "WAYLAND_DISPLAY") || return nothing        # an X11 session: both take X11
    glfw_wayland = lowercase(get(ENV, "JULIA_GLFW_PLATFORM", "")) == "wayland"
    qt_x11       = lowercase(get(ENV, "QT_QPA_PLATFORM", "")) == "xcb"
    (glfw_wayland || qt_x11) && return nothing              # the launcher settled it
    @warn """
        Qt is on Wayland and GLMakie is on XWayland, so this process is running two windowing
        systems at once. Start through the launcher script, or set one of these before the
        application starts -- neither can be changed from inside it:

            JULIA_GLFW_PLATFORM=wayland     both on Wayland (preferred)
            QT_QPA_PLATFORM=xcb             both on X11
        """
    return nothing
end

"""
    use_writable_cache_dir!()

Point Makie's font-atlas cache at a per-user directory.

Makie caches the atlas in a scratch space inside the DEPOT, which in an installed bundle is
read-only: it would then try to download the atlas from GitHub on every launch and, failing
that, re-render every glyph with a warning, having nowhere to store the result. A per-user
cache directory makes the first launch the only slow one -- and the bundle ships a prebuilt
atlas, so ideally not even that.
"""
function use_writable_cache_dir!()
    dir = joinpath(cache_home(), "oitools")
    try
        mkpath(dir)
        seed_font_atlas!(dir)
        GLMakie.Makie.makie_cache_dir[] = dir
    catch err
        @debug "could not set a writable Makie cache directory" dir err
    end
    return dir
end

"""
    seed_font_atlas!(dir)

Copy the bundle's prebuilt font atlas into the per-user cache, for files not already there.

`build.jl` stages the atlas the precompile trace produced, so a bundle has one; a checkout does
not, and `resource` returns `nothing` there, which is why this is a no-op rather than an error.
Files already in `dir` are left alone: the user's own cache is newer than the bundle's by
definition once Makie has written to it.
"""
function seed_font_atlas!(dir)
    src = OITOOLS.resource("makie-cache")
    src === nothing && return 0
    n = 0
    for f in readdir(src)
        dst = joinpath(dir, f)
        isfile(dst) && continue
        cp(joinpath(src, f), dst)
        n += 1
    end
    return n
end

"""
    register_qml_modules!() -> String or nothing

Tell Qt where the bundled `Makie` QML module is.

QMLMakie registers its own module from `__init__` with
`QML.add_import_path(joinpath(@__DIR__, "qml"))`, and `@__DIR__` is expanded at PRECOMPILE
time -- so the path in the image is the package directory of the machine that built the
bundle. It does not exist on the machine running it, Qt finds no module, and every `.qml` that
says `import Makie` fails with `module "Makie" is not installed` before the window appears.

`build.jl` stages those 16 kB beside the other resources; this points Qt at them. A checkout
needs none of it -- `resource` returns `nothing` there and QMLMakie's own path is live.
"""
function register_qml_modules!()
    dir = OITOOLS.resource("qml-modules")
    dir === nothing && return nothing
    QML.add_import_path(dir)
    return dir
end

"Base directory for per-user caches, following the platform's own convention."
cache_home() = Sys.iswindows() ? get(ENV, "LOCALAPPDATA", homedir()) :
               Sys.isapple()   ? joinpath(homedir(), "Library", "Caches") :
               get(ENV, "XDG_CACHE_HOME", joinpath(homedir(), ".cache"))

"Base directory for per-user state, which is where a crash log belongs."
state_home() = Sys.iswindows() ? get(ENV, "LOCALAPPDATA", homedir()) :
               Sys.isapple()   ? joinpath(homedir(), "Library", "Logs") :
               get(ENV, "XDG_STATE_HOME", joinpath(homedir(), ".local", "state"))

"""
    crash_log_path() -> String

Where an application failure is recorded.

An application has no console to print a stacktrace to, and the one a user can find is worth
more than the one they cannot. Appended to, not truncated, so a failure that only happens
every so often still has its predecessors beside it.
"""
crash_log_path() = joinpath(state_home(), "oitools", "crash.log")

function report_crash(err, bt)
    msg = sprint() do io
        println(io, "─"^72)
        println(io, "OITOOLS ", pkgversion(OITOOLS), " crashed at ", Libc.strftime(time()))
        println(io, "─"^72)
        showerror(io, err, bt)
        println(io)
    end
    print(stderr, msg)
    path = crash_log_path()
    try
        mkpath(dirname(path))
        open(io -> print(io, msg), path, "a")
        println(stderr, "\nThis was also written to $path")
    catch
        # Nothing useful is left to do if even the log cannot be written; stderr already has it.
    end
    return path
end

end # module
