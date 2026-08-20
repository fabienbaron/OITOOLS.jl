# Drive `filepicker_min.qml`: the file picker with nothing else around it.
#
#     PICKER_MODE=close julia --project=bin test/gui/filepicker_min.jl
#
# `PICKER_MODE` is one of `none`, `close` or `timer` (default `close`, which is what the real
# GUI does). The transcript goes to stdout and to $PICKER_LOG if set, so an automated run can
# read it after the window is gone.
#
# This deliberately loads NOTHING from the GUI: no Makie, no GL, no session. QML.jl and
# QtQuick.Dialogs are the entire surface under test.

using QML

const LOGLINES = String[]
const LOGFILE  = get(ENV, "PICKER_LOG", "")

function picker_log(s)
    line = String(s)
    push!(LOGLINES, line)
    println(line); flush(stdout)
    isempty(LOGFILE) || open(LOGFILE, "a") do io; println(io, line); end
    return nothing
end

"""
Directory listing for the `inline` mode's picker, as `kind\tname` lines.

The inline picker is the candidate that cannot fail the way QtQuick.Dialogs does: it is an
in-window `Popup`, so there is no second window to be left behind. That means it needs Julia
for what the toolkit would otherwise supply.
"""
function picker_listdir(path)
    dir = abspath(String(path))
    isdir(dir) || return ""
    rows = String[]
    parent = dirname(dir)
    parent == dir || push!(rows, "up\t..")
    names = try; sort(readdir(dir)); catch; String[]; end
    for n in names
        startswith(n, ".") && continue
        full = joinpath(dir, n)
        isdir(full) ? push!(rows, "dir\t" * n) : push!(rows, "file\t" * n)
    end
    return join(rows, "\n")
end

"Absolute, normalised path, so QML never has to do path arithmetic."
picker_join(a, b) = abspath(joinpath(String(a), String(b)))

picker_home() = abspath(String(get(ENV, "PWD", pwd())))

"Hold the GUI thread, the way a real dataset load does."
function picker_block(seconds)
    t0 = time()
    while time() - t0 < Float64(seconds)
    end
    return nothing
end

QML.@qmlfunction picker_log picker_block picker_listdir picker_join picker_home

mode   = get(ENV, "PICKER_MODE", "close")
folder = "file://" * abspath(joinpath(@__DIR__, "data"))
picker_log("=== filepicker_min: mode=$mode folder=$folder ===")

qmlfile = joinpath(@__DIR__, "filepicker_min.qml")
QML.loadqml(qmlfile; pickerMode = mode, initialFolder = folder)

# No autoquit here: the driving script kills the process. Reaching for QML.jl's timer API
# would add a second thing under test to a file whose whole point is having only one.
QML.exec()
picker_log("=== event loop returned ===")
