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
# The one ordering constraint here is `configure_graphics!()`: Mesa reads its driver
# variables when the first OpenGL context is created, so the call has to happen above
# `using GLMakie`. See src/graphics.jl.
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

configure_graphics!()        # before the first GL context — see src/graphics.jl
prefer_native_wayland!()     # before `using GLMakie` — GLFW.jl would otherwise pick XWayland

using GLMakie, QMLMakie, QML   # these three activate OITOOLSGUIExt, which defines gui()

session = Session()
for f in ARGS
    isfile(f) ? load_dataset!(session, f; warn = false, verbose = false) :
                @warn "not a file, skipping" f
end
gui(session)
