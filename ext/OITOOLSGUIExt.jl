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
# The toolkit-free half of the drawing layer, which moved into the core package so that the
# matplotlib front-end, the Makie one and this live canvas cannot drift apart.
using OITOOLS: PlotData, ObsSpec, OBS_SPECS, LOG_Y_KINDS,
               OIPLOT_COLORS, OIPLOT_ECOLOR, OIPLOT_MARKER,
               OIPLOT_LABELSIZE, OIPLOT_TICKLABELSIZE,
               OIPLOT_LEGEND_SIZE, OIPLOT_LEGEND_NCOL,
               UVPLOT_LEGEND_SIZE, UVPLOT_LEGEND_NCOL, OIPLOT_XTICKS,
               OIPLOT_MINOR_INTERVALS, imdisp_tickinterval,
               PLOT_FONT, PLOT_FONTS, PLOT_GLYPHS,
               baseline_color_map, baseline_names, triplet_names, station_names,
               group_names, grouping_noun, panel_data, point_info, obs_info,
               uv_point_labels, _colors_for
# The Makie drawing layer, which belongs to OITOOLSMakieExt now. `using GLMakie` activates it,
# so it is always present by the time this extension exists — but one extension cannot import
# from another, so these arrive through the functions the parent package declares.
using OITOOLS: prewarm_glyphs!, style_axis!, add_baseline_legend!,
               draw!, plot_into!, image_minorticks!
# Not exported by the parent, and both are needed by src/gui/imaging.jl: the precision
# cast that keeps a starting image loadable by the Fourier plan, and the criterion the
# panel reports before and after a run.
using OITOOLS: to_ft_precision, image_to_chi2
# The component-kind tables, so "+ component" builds the exact key set the parser identifies a
# kind by. Duplicating them here would mean the GUI could offer a component OITOOLS cannot read.
using OITOOLS: _KIND_GEOMETRY, _ANALYTIC_PARAM_SUFFIXES
using Printf, Dates, TOML
using Statistics: std          # posterior spread, for the fit table's ± column
using Makie
using QML, QMLMakie, GLMakie

# Exported so tests and scripts can `using .GUI` after fetching the module with
# `Base.get_extension`. `gui` itself is exported by OITOOLS, since it is declared there.
export Session, DatasetEntry, ModelEntry, ImageEntry, Selection, LogEntry
export log!, export_script, load_dataset!, filter_dataset!, dataset
export point_info, uv_point_labels, plot_into!
export OBS_SPECS, group_names, baseline_names, triplet_names, station_names
export panel_data, grouping_noun
export baseline_color_map, style_axis!, add_baseline_legend!
export LiveCanvas, build_canvas, update_canvas!, canvas_data, zoom_step!,
       set_colormap!, image_colormap_names, IMAGE_COLORMAPS,
       show_panels!, update_panels!, MAX_PANELS,
       ZOOM_MIN_SPAN, ZOOM_MAX_SPAN, ZOOM_PER_DETENT
export ParamRow, ParamMode, PARAM_FIXED, PARAM_FREE, PARAM_EXPR
export model_rows, model_inspection, ComponentInfo, free_parameter_vector
export model_depends_on, render_model_image, model_image_cube
export parse_model_lines, parse_free_lines, parse_constraint_lines, parse_prior_lines,
       fitting_weights
export ImagingSetup, ImagingResult, imaging_defaults, imaging_weights, fov,
       start_image, reconstruct_image, observable_availability, observable_flags_string,
       parse_regularizers, ft_summary, ensure_ft!, AUTO_FWHM, chi2_breakdown, chi2r, chi2r_start,
       IMAGING_ENGINES, POLYCHROMATIC_ENGINES, run_engine, prior_image, result_ensemble
export NightPlan, config_catalog, facility_telescopes, telescope_config,
       night_plan, observable_indices, observable_hours, plan_rows, best_pops, pop_rows,
       GanttBar, GanttLabel, gantt_geometry, unwrap_lst, build_gantt, update_gantt!,
       pop_label,
       baseline_delay_windows, DEFAULT_ALT_LIMIT, DEFAULT_ALT_MAX, TARGET_BAR_HEIGHT,
       delay_plot_geometry, build_delay_plot, update_delay_plot!, simulate_source_info
export build_chi2_map, update_chi2_map!
export build_residuals, update_residuals!, residual_series,
       RESIDUAL_KINDS, RESIDUAL_SPECS, RESIDUAL_COLORS
export build_sed, update_sed!, SED_MAX_COMPONENTS
export ZOOM_STEP_USER, zoom_per_detent, set_zoom_step!
export picker_list, picker_places, picker_join, picker_parent, picker_start,
       picker_kind, picker_would_overwrite, picker_matches, picker_examples
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
include(joinpath(GUIDIR, "observing.jl"))   # Observe perspective: configs, targets, nights
include(joinpath(GUIDIR, "gantt.jl"))       # the Gantt chart, in Makie
include(joinpath(GUIDIR, "livecanvas.jl"))    # the live, allocation-free drawing surface
include(joinpath(GUIDIR, "chi2map.jl"))       # the grid-search chi2 surface, in Makie
include(joinpath(GUIDIR, "residuals.jl"))     # normalised model residuals, in Makie
include(joinpath(GUIDIR, "sed.jl"))           # the model SED, in Makie
include(joinpath(GUIDIR, "shell.jl"))
include(joinpath(GUIDIR, "snapshot.jl"))      # "Save PNG": rebuilds a panel offscreen
include(joinpath(GUIDIR, "window.jl"))        # gui(): builds the window and runs Qt

# ── precompilation ───────────────────────────────────────────────────────────
#
# The core package's workload cannot reach any of this: these methods only exist once Makie and
# Qt are loaded, so they are compiled in the user's session on first use — which is exactly
# when a window is already open and every pause is visible.
#
# Only what runs headless is here. Building a Figure, an Axis and every plot on it is
# scene-graph construction and needs no GL; rendering does, and nothing below renders. `gui()`
# itself and anything touching QML are deliberately absent — they need a display.
using PrecompileTools

@setup_workload begin
    _pcfile = joinpath(pkgdir(OITOOLS), "test", "gui", "data", "2004-data1.oifits")
    @compile_workload begin
        if isfile(_pcfile)
            data = readoifits(_pcfile; verbose = false, warn = false)
            d = data[1, 1]

            # The redraw path: turning an OIdata into points and colours is what every plot
            # change runs, and it is pure Julia.
            for kind in (:uv, :v2, :t3phi)
                canvas_data(d, kind; color = :baseline)
            end
            uv_point_labels(d)
            observable_flags_string(data)

            # The panels' own data layers.
            imaging_defaults(data)

            # The whole imaging path as the Run button takes it. Precompiling core's
            # `reconstruct` is NOT enough: the call goes through this extension's
            # `reconstruct_image`, and inference through that wrapper is its own
            # specialisation — measured at 7.7 s on first use with only the core workload in
            # place, against 0.2 s of actual work. Two iterations compile all of it.
            _s = ImagingSetup(; nx = 16, pixsize = 0.5, startkind = :gaussian)
            reconstruct_image(data, _s; maxiter = 2)
            reconstruct_image(data, _s; maxiter = 2, regularizers = "centering,1e4")
            start_image(ImagingSetup(; nx = 16, pixsize = 0.5, startkind = :dirac),
                        setup_ft(data, 16, 0.5))
            parse_regularizers("centering,1e5;l1l2,1e7,1e-3")
            md = Dict{String,Any}("s,ud" => 3.0, "s,f" => 1.0)
            model_rows(md, ["s,ud"])
            model_inspection(md)
            picker_list(dirname(_pcfile), "*.oifits", false)
            picker_places()

            # The Observe perspective. `build_gantt` alone is 1.02 s of compilation on first
            # use, and it runs while the user is waiting for a window — the same argument as
            # for `build_canvas` below.
            plan = night_plan("CHARA", "Vega", 279.2347, 38.7837, DateTime(2026, 6, 21);
                              config = [1, 1, 1, 1, 0, 0], use_delay = true, step_minutes = 10)
            gantt_geometry(plan; show_alt = true)
            gantt_geometry(plan; show_alt = false)
            plan_rows([plan])

            # The POP search, which AutoPOPs now runs on every Compute rather than only when
            # asked. 11.6 ms of work at six telescopes and far more than that to compile the
            # first time, and it happens while the user is waiting for a Gantt.
            let f = read_facility_file("CHARA"),
                obs = night_observability(f, 279.2347, 38.7837, DateTime(2026, 6, 21))
                best_pops(f, 38.7837, obs.ha, [1, 1, 1, 1, 0, 0]; n = 3)
            end
            gfig = Makie.Figure()
            gax  = Makie.Axis(gfig[1, 1])
            gg   = build_gantt(gfig, gax)
            update_gantt!(gg, plan)

            # The canvas itself. `build_canvas` creates one of every plot the shell can show,
            # and that construction is a large part of the delay before the first window
            # appears.
            fig = Makie.Figure()
            ax  = Makie.Axis(fig[1, 1])
            # `style_axis!` belongs to OITOOLSMakieExt and is compiled by its workload. It is
            # not called here because sibling extensions have no guaranteed load order: at
            # precompile time this one can exist before that one does, and the call would then
            # find a function with no methods. At run time both are loaded long before `gui()`.
            c = build_canvas(fig, ax)
            # `update_canvas!` calls `style_axis!`, which belongs to OITOOLSMakieExt. Sibling
            # extensions are not activated while this one is being precompiled, so the guard
            # is what decides whether these get compiled here or in the user's first session.
            # At RUN time the ordering is not in doubt: OITOOLSMakieExt triggers on Makie,
            # which GLMakie pulls in before QML and QMLMakie are anywhere near loaded.
            if !isempty(methods(style_axis!))
                update_canvas!(c, d, :uv; color = :baseline)
                update_canvas!(c, d, :v2; color = :wav)
            end
            show_image!(c, zeros(Float32, 8, 8), 0.3)
            hide_image!(c)

            # The delay-cart chart and the χ² map, for the same reason as the two above: both
            # are built before the window exists, so their compilation is time the user spends
            # looking at nothing. Measured on a cold session, `build_chi2_map` alone was 1.72 s
            # -- the most expensive single step in `gui()` -- because `contour!` is a recipe
            # that had never been compiled anywhere else in the package.
            dfig = Makie.Figure()
            dax  = Makie.Axis(dfig[1, 1])
            build_delay_plot(dfig, dax)

            ofig = Makie.Figure()
            oax  = Makie.Axis(ofig[1, 1])
            om   = build_chi2_map(ofig, oax)
            # Drawing one too: `update_chi2_map!` is where the heatmap, the contour levels and
            # the marker actually get their data, and that is a separate specialisation.
            let mdl = Dict{String,Any}("s,ud" => 3.0, "s,f" => 1.0),
                mp = chi2_map(mdl, ["s,ud", "s,f"], d, "s,ud", "s,f";
                              range1 = (1.0, 5.0), range2 = (0.5, 1.5), n1 = 4, n2 = 4)
                update_chi2_map!(om, mp)
            end
        end
    end
end

end # module
