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
# The trigger lists four packages, but a caller only ever names three: GLMakie pulls in Makie,
# so `using GLMakie, QMLMakie, QML` loads all four and fires this extension. Makie is named
# because the sources below use it directly, and an extension may only `using` the parent's
# dependencies and its own triggers.
#
# `configure_graphics!` and `prefer_native_wayland!` are NOT here. Both must run before the
# first GL context exists, which is before `using GLMakie` — and therefore before this
# extension can exist at all. The first lives in the core package (`src/graphics.jl`), the
# second in `ext/OITOOLSGLFWExt.jl`, triggered by GLFW_jll alone.
#
# Sources live in `src/gui/` and are included from here. `Session`, `ShellState` and
# `LiveCanvas` are therefore extension types: functions can be forward-declared in the parent
# (`function gui end`), types cannot. `gui()` builds its own `Session`, so callers never need to
# name one; tests reach the rest through `Base.get_extension(OITOOLS, :OITOOLSGUIExt)`.
module OITOOLSGUIExt

using OITOOLS
using OITOOLS: OIdata, OBS_PLOT_SPECS, oiplot_colors, canonical_color, as_datavec
# Not exported by the parent, and both are needed by src/gui/imaging.jl: the precision
# cast that keeps a starting image loadable by the Fourier plan, and the criterion the
# panel reports before and after a run.
using OITOOLS: to_ft_precision, image_to_chi2
using Printf, Dates
using Makie
using QML, QMLMakie, GLMakie

# Exported so tests and scripts can `using .GUI` after fetching the module with
# `Base.get_extension`. `gui` itself is exported by OITOOLS, since it is declared there.
export Session, DatasetEntry, ModelEntry, ImageEntry, Selection, LogEntry
export log!, export_script, load_dataset!, filter_dataset!, dataset
export uv_figure, observable_figure, point_info, uv_point_labels, plot_into!
export OBS_SPECS, group_names, baseline_names, triplet_names, station_names
export baseline_color_map, style_axis!, add_baseline_legend!
export LiveCanvas, build_canvas, update_canvas!, canvas_data
export ParamRow, ParamMode, PARAM_FIXED, PARAM_FREE, PARAM_EXPR
export model_rows, model_inspection, ComponentInfo, free_parameter_vector
export ImagingSetup, ImagingResult, imaging_defaults, imaging_weights, fov,
       start_image, reconstruct_image, observable_availability, observable_flags_string,
       parse_regularizers
export picker_list, picker_places, picker_join, picker_parent, picker_start,
       picker_kind, picker_would_overwrite, picker_matches
export GLOBAL_COMPONENT
export ShellState, check_qt_conflict, ui_scale_override

const GUIDIR = joinpath(pkgdir(OITOOLS), "src", "gui")

include(joinpath(GUIDIR, "scaling.jl"))       # UI scale policy; QML does the detecting
include(joinpath(GUIDIR, "commandlog.jl"))
include(joinpath(GUIDIR, "session.jl"))
include(joinpath(GUIDIR, "actions_data.jl"))
include(joinpath(GUIDIR, "model.jl"))       # Model perspective: parameter table, parser inspection
include(joinpath(GUIDIR, "imaging.jl"))     # Image perspective: geometry, start image, reconstruct
include(joinpath(GUIDIR, "filepicker.jl"))  # the in-window file picker's directory listing
include(joinpath(GUIDIR, "plots.jl"))
include(joinpath(GUIDIR, "livecanvas.jl"))    # the live, allocation-free drawing surface
include(joinpath(GUIDIR, "shell.jl"))
include(joinpath(GUIDIR, "window.jl"))        # gui(): builds the window and runs Qt

end # module
