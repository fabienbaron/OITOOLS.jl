# OITOOLS.jl
#
# Dependencies:
#   emmt/OIFITS.jl (new API: OIDataSet, read(OIDataSet, …), OIFITS.select(ds, OI_VIS2), …)
#   FITSIO        (for readfits / writefits)
#   All other deps unchanged
#
# Float precision:
#   OIdata is now OIdata{T<:AbstractFloat}.
#   Default T=Float32; use readoifits(file; T=Float64) for double precision.
#   NFFT plans, DFT matrices, and chi2 functions propagate T automatically.

module OITOOLS

using PrecompileTools

version() = println("OITOOLS v$(pkgversion(OITOOLS))")

include("graphics.jl")           # GL driver environment; must be callable before GLMakie
include("readoifits.jl")
include("oiplot_specs.jl")   # observable/plot metadata: no toolkit, used by every front-end
include("vis_functions.jl")
include("model_chainrules.jl")

# ── Flat-dict model fitting stack ────────────────────────────────────────────
include("resolvers.jl")          # SharedUtils, RGF, HandRolled, Symbolic
include("hankel.jl")             # Hankel transform + chainrules + profile compiler
include("parse_model.jl")        # FlatModel, eval_model, eval_model_grad
include("chi2_flat.jl")          # chi2_flat, chi2_flat_fg
include("constraints.jl")        # ModelConstraint: NLopt constraints, and penalties elsewhere
include("model_file.jl")         # TOML model files: model + free + bounds + constraints + priors
include("fit_model.jl")          # fit_model (NLopt), fit_model_ultranest, display_model
include("bootstrap.jl")          # block bootstrap + parametric Monte Carlo


include("utils.jl")
include("oichi2.jl")
include("sparco_flat.jl")
include("astrometry.jl")
include("atmosphere.jl")     # bands, seeing, AO Strehl — used by simulate.jl
include("vonmises.jl")
include("simulate.jl")
include("oifitslib.jl")
include("maximent.jl")
include("oimem.jl")
include("proximal_operators.jl")
include("bsdmm.jl")
include("squeeze.jl")
include("plan.jl")
include("pmoired_compat.jl")

# ── Plotting: declared here, implemented in ext/OITOOLSPythonPlotExt.jl ──────
#
# Every matplotlib call in the package lives in src/oiplot.jl, which is loaded only when the
# caller has PythonPlot. These are the function objects the extension adds methods to.
#
# Why: importing OITOOLS used to import matplotlib, which probes for an interactive backend,
# finds PySide6 in the CondaPkg environment and maps a Qt into the process. That cost a
# ~3 s load, broke `Pkg.test()` and the docs build (each project building its own 1.5 GB
# Python environment), and made OITOOLS unusable alongside a Qt-based Julia GUI. None of that
# is needed to read an OIFITS file.
#
# Calling one of these without PythonPlot raises a MethodError naming the function; see the
# note in the package README.
function chara_plan end
function gantt_onenight end
function imdisp end
function imdisp_multi end
function obs_plan end
function plot_diffphi end
function plot_facility end
function plot_flux end
function plot_multi end
function plot_obs end
function plot_residuals end
function plot_t3amp end
function plot_t3amp_residuals end
function plot_t3phi end
function plot_t3phi_residuals end
function plot_v2 end
function plot_v2_multifile end
function plot_v2_residuals end
function plot_visamp end
function plot_visamp_residuals end
function plot_visphi end
function plot_visphi_residuals end
function set_oiplot_defaults end
function plot_ultranest_corner end

# The GUI. Defined by OITOOLSGUIExt, which loads with GLMakie + QMLMakie + QML:
#     using OITOOLS, GLMakie, QMLMakie, QML
#     gui()
# A function can be declared here and given methods by an extension; the GUI's types
# (Session, ShellState, LiveCanvas) cannot, and stay inside the extension.
function gui end
function uvplot end
function empty_night end

# Native Wayland for GLMakie. Defined by OITOOLSGLFWExt, which loads with GLFW_jll alone —
# deliberately, because this has to be callable BEFORE `using GLMakie`. See src/graphics.jl.
function prefer_native_wayland! end

export configure_graphics!, prefer_native_wayland!, is_wsl

# ── Reading OIFITS data ─────────────────────────────────────────────────────
export OIdata
export readoifits, readoifits_multiepochs, list_oifits_targets
export readfits, writefits
export oifits_prep, updatefits_aspro
export filter_data, set_data_filter

# ── Plotting ─────────────────────────────────────────────────────────────────
export set_oiplot_defaults, uvplot, plot_v2, plot_t3phi, plot_t3amp,
       plot_ultranest_corner, gui,
       plot_visamp, plot_visphi, plot_diffphi, plot_flux, plot_obs, plot_multi,
       plot_v2_residuals, plot_t3phi_residuals, plot_t3amp_residuals,
       plot_visamp_residuals, plot_visphi_residuals, plot_residuals,
       plot_v2_multifile,
       imdisp, imdisp_multi, plot_facility

# ── Model fitting (flat-dict interface) ──────────────────────────────────────
export FlatModel, dict_to_model, parse_model, eval_model, eval_model_grad, display_model, default_bounds, max_angular_scale, DEFAULT_MAX_SIZE_MAS
export fit_model, fit_model_lsqfit, fit_model_ultranest
export ModelConstraint, parse_constraints, check_constraints, DEFAULT_CONSTRAINT_TOL
export read_model_file, write_model_file
export FitResult, LsqFitResult, UltraNestResult
export model_to_vis, model_to_obs, model_to_residuals, model_to_chi2, model_to_chi2_fg,
       model_to_image, model_to_sed, model_to_flux
export cvis_to_chi2_f, cvis_to_chi2_fg, cvis_to_chi2_noalloc, model_and_image_to_chi2_fg

# ── Uncertainty estimation by resampling ─────────────────────────────────────
export bootstrap_fit, bootstrap_driver, BootstrapResult
export data_blocks, DataBlocks, resample_blocks
export block_counts, apply_block_counts, block_weights, apply_block_weights
export perturb_data
export resample_data          # deprecated alias for perturb_data

# ── Visibility functions ─────────────────────────────────────────────────────
export visibility_ud, visibility_ldlin, visibility_ldquad, visibility_ldquad_tri,
       visibility_ldpow, visibility_ldsquareroot,
       visibility_annulus, visibility_ellipse_uniform, visibility_ellipse_quad,
       visibility_thin_ring, visibility_Gaussian_ring, visibility_Gaussian_ring_az,
       visibility_Lorentzian_ring, visibility_GaussianLorentzian_ring_az

# ── Image reconstruction ─────────────────────────────────────────────────────
export setup_ft, setup_dft, setup_nfft, ft_info, OIft, NFFTCell, DFTCell
export image_to_vis,
       image_to_v2, image_to_t3phi, image_to_t3amp, image_to_obs,
       image_to_residuals, image_to_chi2, image_to_chi2_fg
export vis_to_v2, vis_to_t3, observables
export crit_f, crit_fg, image_to_crit
export reconstruct
export gaussian2d, recenter
export scale_regularizers, scale_tv_epsilon
export lcurve, lcurve_elbow, lcurve_normalize

# ── SPARCO ──────────────────────────────────────────────────────────────────
export reconstruct_sparco, reconstruct_hybrid
export reconstruct_sparco_flat, chi2_sparco_flat_f, optimize_sparco_flat_parameters
export model_and_image_to_vis
# Deprecated (error on call):
export reconstruct_sparco_gray, reconstruct_sparco_multi
export chi2_sparco_multi_f, optimize_sparco_multi_parameters

# ── BSDMM reconstruction ──────────────────────────────────────────────────
export reconstruct_bsdmm
export reconstruct_squeeze, reconstruct_squeeze_tempered, SqueezeSparco, default_nelements
export prox_tv, prox_l2smooth, prox_centering!
export tv_norm, l2smooth_norm
export prox_group_sparsity, group_sparsity_norm
export prox_grouptv, grouptv_norm

# ── Deprecated (still exported for backward compatibility) ──────────────────
export setup_dft_polychromatic, setup_nfft_multiepochs, setup_nfft_polychromatic
export chi2_f, chi2_fg, chi2_polychromatic_f
export crit_polychromatic_fg, crit_multitemporal_fg
export reconstruct_multitemporal, reconstruct_polychromatic

# ── Maximum entropy (BSMEM) ──────────────────────────────────────────────────
export reconstruct_bsmem
export auto_pixsize, gaussian_prior, gaussian_prior_cube

# ── Simulation ───────────────────────────────────────────────────────────────
export simulate, simulate_from_oifits, get_uv
export read_facility_file, read_obs_file, read_comb_file, read_wave_file
export list_configs
export predict_errors
export FacilityConfig, TargetConfig, CombinerConfig, WaveConfig, AOConfig
export strehl_ratio, coupling_efficiency, fried_parameter, seeing_from_r0, atm_transmission
export PhotometricBand, PHOTOMETRIC_BANDS, band_for_wavelength, band_by_name, zero_point_flux
export gantt_onenight
export sunrise_sunset, alt_az, opd_limits, airmass, datetime_to_jd, datetime_to_mjd
export query_target_from_simbad, ra_dec_from_simbad, magnitudes_from_simbad, sexagesimal_to_degrees
export simbad_target, SIMBAD_BANDS

# ── Observation planning (CHARA) ──────────────────────────────────────────────
export night_observability, compute_delays, in_delay, best_pop, observable_epochs
export index_runs
export obs_plan, chara_plan, print_pop_results
export moon_radec, moon_illumination, angular_separation

# ── Writing OIFITS ───────────────────────────────────────────────────────────
export oifits_check, oifits_merge, oifits_filter, oifits_fix_tdim

# ── PMOIRED compatibility ────────────────────────────────────────────────────
export pmoired_to_julia, pmoired_to_julia_file, pmoired_to_dict, dict_to_pmoired, dict_to_pmoired_file

# ── Precompilation workload ──────────────────────────────────────────────────
include("precompile.jl")

end
