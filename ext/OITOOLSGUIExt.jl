# The OITOOLS GUI: four perspectives over one session, in a QML window drawing with Makie.
#
#     using OITOOLS, GLMakie, QMLMakie, QML
#     gui()
#
# The whole GUI lives here rather than in the core package because Makie and Qt together cost
# about 6.4 s of load time, and `using OITOOLS` is 1.69 s. A reduction script, a CI run, or any
# session that never opens a window should not pay that. Gating the GUI behind its own stack
# means they do not.
#
# The trigger lists five packages, but a caller only ever names three: GLMakie pulls in Makie
# and GLFW_jll, so `using GLMakie, QMLMakie, QML` loads all five and fires this extension. Makie
# and GLFW_jll are named because the sources below use them directly, and an extension may only
# `using` the parent's dependencies and its own triggers.
#
# Sources live in `src/gui/` and are included from here. `Session`, `ShellState` and
# `LiveCanvas` are therefore extension types: functions can be forward-declared in the parent
# (`function gui end`), types cannot. `gui()` builds its own `Session`, so callers never need to
# name one; tests reach the rest through `Base.get_extension(OITOOLS, :OITOOLSGUIExt)`.
module OITOOLSGUIExt

using OITOOLS
using OITOOLS: OIdata, OBS_PLOT_SPECS, oiplot_colors, canonical_color, as_datavec
using Printf, Dates
using Makie
using GLFW_jll                    # prefer_native_wayland! only; needs no display to load
using QML, QMLMakie, GLMakie

# Exported so tests and scripts can `using .GUI` after fetching the module with
# `Base.get_extension`. `gui` itself is exported by OITOOLS, since it is declared there.
export Session, DatasetEntry, ModelEntry, ImageEntry, Selection, LogEntry
export log!, export_script, load_dataset!, filter_dataset!, dataset
export uv_figure, observable_figure, point_info, uv_point_labels, plot_into!
export OBS_SPECS, group_names, baseline_names, triplet_names, station_names
export baseline_color_map, style_axis!, add_baseline_legend!
export LiveCanvas, build_canvas, update_canvas!, canvas_data
export ShellState, check_qt_conflict, configure_graphics!, ui_scale_override
export prefer_native_wayland!

const GUIDIR = joinpath(pkgdir(OITOOLS), "src", "gui")

include(joinpath(GUIDIR, "graphics.jl"))      # WSL OpenGL setup; before any GL context
include(joinpath(GUIDIR, "scaling.jl"))       # UI scale policy; QML does the detecting
include(joinpath(GUIDIR, "commandlog.jl"))
include(joinpath(GUIDIR, "session.jl"))
include(joinpath(GUIDIR, "actions_data.jl"))
include(joinpath(GUIDIR, "plots.jl"))
include(joinpath(GUIDIR, "livecanvas.jl"))    # the live, allocation-free drawing surface
include(joinpath(GUIDIR, "shell.jl"))
include(joinpath(GUIDIR, "window.jl"))        # gui(): builds the window and runs Qt

end # module
