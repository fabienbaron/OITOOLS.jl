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
# The window is opened by `oitoolsgui` (src/gui_launcher.jl), which sets the Mesa, GLFW and Qt
# variables in the order they must come and then loads GLMakie, QMLMakie, QML and the optional
# packages. UltraNest needs PythonCall, which is loaded by hand:
#
#     julia --project=bin -e 'using PythonCall; include("bin/oitoolsgui.jl")'

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
oitoolsgui(ARGS...)
