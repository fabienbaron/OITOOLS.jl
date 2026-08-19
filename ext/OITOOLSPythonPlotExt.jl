# matplotlib plotting for OITOOLS.
#
# Loaded automatically when the caller has PythonPlot. Everything here used to be in the core
# load path, which meant `using OITOOLS` imported matplotlib whether or not anything was ever
# plotted — see the note beside the stub declarations in src/OITOOLS.jl.
#
# `import`, not `using`: the definitions inside oiplot.jl are written as plain
# `function uvplot(...)`, and importing the names makes those extend OITOOLS' functions rather
# than defining new ones local to this module. That is what lets the file move unchanged.
module OITOOLSPythonPlotExt

using OITOOLS
using PythonPlot, PythonCall, LaTeXStrings, Statistics, LinearAlgebra, Dates, Printf

import OITOOLS: chara_plan, gantt_onenight, imdisp, imdisp_multi, obs_plan, plot_diffphi,
                plot_facility, plot_flux, plot_multi, plot_obs, plot_residuals, plot_t3amp,
                plot_t3amp_residuals, plot_t3phi, plot_t3phi_residuals, plot_v2,
                plot_v2_multifile, plot_v2_residuals, plot_visamp, plot_visamp_residuals,
                plot_visphi, plot_visphi_residuals, set_oiplot_defaults, uvplot,
                plot_ultranest_corner,
                empty_night, monitor_update!

# Internals oiplot.jl reaches for.
using OITOOLS: OIdata, FlatModel, SqueezeMonitor, FacilityConfig,
               image_to_obs, image_to_residuals, model_to_obs, model_to_residuals,
               hours_to_date, night_observability, in_delay, sunrise_sunset,
               alt_az, datetime_to_jd, moon_radec, moon_illumination, angular_separation,
               best_pop, cdg, recenter, hour_angle_calc, mod360, _is_obstype,
               SPARCO_PARAM_NAMES,
               # plot metadata, which stayed in core so non-matplotlib front-ends can use it
               ObsPlotSpec, OBS_PLOT_SPECS, canonical_color, as_datavec,
               oiplot_colors, oiplot_markers

include(joinpath(pkgdir(OITOOLS), "src", "oiplot.jl"))

end # module
