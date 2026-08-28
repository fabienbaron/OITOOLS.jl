# Plotting in Makie, without Qt and without Python.
#
# Loaded when the caller has Makie — which any backend brings: `using CairoMakie` gets vector
# figures with no GPU and no display, `using GLMakie` gets a window, and the GUI gets it
# through QMLMakie. Choosing the backend is the caller's job and this file never mentions one.
#
# Why it is not part of OITOOLSGUIExt, where this code started: everything here draws a figure
# and hands it back, so none of it needs Qt. Splitting it out is what lets a script — or a
# PackageCompiler build, which cannot contain PythonPlot — plot at all. The GUI extension
# keeps only what is genuinely interactive: the live canvas, the Gantt and χ² panels whose
# Observables are rewritten in place, and the window itself.
#
# `import`, not `using`: the definitions inside src/oiplot_makie.jl are written as plain
# `function draw!(...)`, and importing the names makes those extend OITOOLS' functions rather
# than defining new ones local to this module. That is the same arrangement
# OITOOLSPythonPlotExt uses, and it is what lets the file be read as the drawing layer and
# nothing else.
module OITOOLSMakieExt

using OITOOLS
using Makie
using Printf

import OITOOLS: prewarm_glyphs!, style_axis!, add_baseline_legend!,
                draw!, plot_into!, uvplot_makie, plot_observable_makie,
                imdisp_makie, image_into!, image_minorticks!, imdisp_multi_makie, plot_residuals_makie,
                plot_obs_makie, plot_v2_multifile_makie, plot_facility_makie,
                plot_v2_makie, plot_t3phi_makie, plot_t3amp_makie,
                plot_visamp_makie, plot_visphi_makie, plot_flux_makie, plot_diffphi_makie

# The toolkit-free descriptions this layer draws from. They live in the core package so that
# the matplotlib front-end, this one, and the GUI's live canvas cannot disagree about which
# array feeds an observable or what colour a baseline is.
using OITOOLS: OIdata, PlotData, ObsSpec, OBS_SPECS, LOG_Y_KINDS,
               OIPLOT_COLORS, OIPLOT_ECOLOR, OIPLOT_MARKER,
               OIPLOT_LABELSIZE, OIPLOT_TICKLABELSIZE,
               OIPLOT_LEGEND_SIZE, OIPLOT_LEGEND_NCOL,
               UVPLOT_LEGEND_SIZE, UVPLOT_LEGEND_NCOL, OIPLOT_XTICKS,
               OIPLOT_MINOR_INTERVALS, imdisp_tickinterval,
               PLOT_FONT, PLOT_FONTS, PLOT_GLYPHS,
               baseline_color_map, baseline_names, triplet_names, station_names,
               group_names, grouping_noun, panel_data, point_info, obs_info,
               uv_point_labels, _colors_for, mod360, FlatModel

include(joinpath(pkgdir(OITOOLS), "src", "oiplot_makie.jl"))

# ── precompilation ───────────────────────────────────────────────────────────
#
# The core package's workload cannot reach any of this: these methods exist only once Makie is
# loaded. Everything here runs headless — building a Figure and drawing into it needs no
# backend, and no backend is named in this file — so it is safe in any environment.
#
# OITOOLSGUIExt has its own workload for the interactive pieces and deliberately does NOT call
# into this one: sibling extensions have no guaranteed load order.
using PrecompileTools

@setup_workload begin
    file = joinpath(pkgdir(OITOOLS), "test", "oifits_for_tests", "2004-data1.oifits")
    @compile_workload begin
        if isfile(file)
            d = OITOOLS.readoifits(file; warn = false, verbose = false)[1, 1]
            fig = Makie.Figure()
            ax  = Makie.Axis(fig[1, 1])
            style_axis!(ax)
            # `draw!` is the single drawing implementation and by far the most expensive thing
            # to compile here: one uv view and one observable view cover both of its branches.
            draw!(fig, ax, d, :uv; color = :baseline)
            draw!(fig, ax, d, :v2; color = :wav)
            add_baseline_legend!(fig, ax, baseline_names(d))
        end
    end
end

end # module
