#!/usr/bin/env julia
#
# Launcher.
#
#     julia --project=bin bin/oitoolsgui.jl [file.oifits]
#
# Use --project=bin, not --project=. — bin/Project.toml is where GLMakie, QMLMakie and QML are
# real dependencies, so bin/Manifest.toml pins them. Under --project=. those three are weak
# dependencies and are not loadable at all; if they resolve anyway it is from your default
# environment, at whatever versions happen to be there.
#
# The ordering below is the whole point of this file, and it is forced by one fact: Mesa and
# GLFW both read their configuration when the first OpenGL context is created, so both hints
# have to be set above `using GLMakie`. That is why neither function lives in the GUI
# extension — loading GLMakie is what creates it, so anything inside it is already too late:
#
#   configure_graphics!    core package,          src/graphics.jl      (needs nothing)
#   prefer_native_wayland! GLFW_jll extension,    ext/OITOOLSGLFWExt.jl (needs libglfw only)
#
# The GUI's own types (Session, load_dataset!) do live in OITOOLSGUIExt, so they are reached
# after the fact through `Base.get_extension`.
#
# MPLBACKEND does not need forcing here: OITOOLS keeps matplotlib in an extension
# (OITOOLSPythonPlotExt), so importing it starts no backend probe and maps no second Qt.
# `check_qt_conflict` runs inside `gui` in case something else in the session loads one.

# A CondaPkg resolve that gets interrupted leaves its lock file behind, and every later start
# then blocks on it FOREVER -- the process stays alive, no window ever appears, and the only
# clue is one CondaPkg info line. That is indistinguishable from a graphics hang, and it cost
# a long debugging session. So: notice it, say what it is, and say how to clear it.
let lock = joinpath(@__DIR__, "..", ".CondaPkg", "lock")
    if isfile(lock)
        age = round(Int, time() - mtime(lock))
        @warn """
            A CondaPkg lock file exists and startup will block on it until it is released.

                $lock   (age $(age)s)

            If no other Julia process is resolving a Conda environment right now, it is stale
            -- delete it and start again:

                rm -f $lock
            """
    end
end

using OITOOLS

configure_graphics!()          # before the first GL context — see src/graphics.jl

using GLFW_jll                 # activates OITOOLSGLFWExt; dlopens libglfw, needs no display
prefer_native_wayland!()       # before `using GLMakie` — GLFW.jl would otherwise pick XWayland

using GLMakie, QMLMakie, QML   # these three activate OITOOLSGUIExt, which defines gui()

# An extension's exports do not reach the caller on their own, so name the module and pull them
# in. This is the same route test/gui/runtests.jl uses.
const GUI = Base.get_extension(OITOOLS, :OITOOLSGUIExt)
GUI === nothing && error("OITOOLSGUIExt did not load; GLMakie, QMLMakie and QML are all needed")
using .GUI

session = Session()
for f in ARGS
    isfile(f) ? load_dataset!(session, f; warn = false, verbose = false) :
                @warn "not a file, skipping" f
end
gui(session)
