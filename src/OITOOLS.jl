# OITOOLS.jl
#
# Dependencies:
#   emmt/OIFITS.jl (new API: OIDataSet, read(OIDataSet, …), OIFITS.select(ds, OI_VIS2), …)
#   FITSIO        (for readfits / writefits)
#   All other deps unchanged
#
# Float precision:
#   OIdata is now OIdata{T<:AbstractFloat}.
#   Default T=Float64 preserves full backward compatibility.
#   Use readoifits(file; T=Float32) for half-memory datasets.
#   NFFT plans, DFT matrices, and chi2 functions propagate T automatically.

module OITOOLS

include("readoifits.jl")
include("vis_functions.jl")
include("model_chainrules.jl")

# ── Flat-dict model fitting stack ────────────────────────────────────────────
include("resolvers.jl")          # SharedUtils, RGF, HandRolled, Symbolic
include("hankel.jl")             # Hankel transform + chainrules + profile compiler
include("parse_model.jl")        # FlatModel, eval_model, eval_model_grad
include("chi2_flat.jl")          # chi2_flat, chi2_flat_fg
include("fit_model.jl")          # fit_model (NLopt), fit_model_ultranest, display_model

include("write_oifits_ha.jl")
include("write_oifits_obs.jl")
include("utils.jl")
include("oichi2.jl")
include("oiplot.jl")
include("astrometry.jl")
include("vonmises.jl")
include("simulate.jl")
include("oifitslib.jl")
include("maximent.jl")
include("oimem.jl")
include("pmoired_compat.jl")

# readoifits
export OIdata
export readoifits, readoifits_multiepochs, readfits, writefits
export oifits_prep, updatefits_aspro, readoifits_multicolors, list_oifits_targets
export remove_redundant_uv!, filter_data, set_data_filter


# oiplot
export set_oiplot_defaults, uvplot, plot_v2,
       plot_diffphi, plot_visphi, plot_t3phi, imdisp,
       imdisp_temporal, plot_v2_residuals, plot_t3phi_residuals, plot_t3amp_residuals,
       plot_t3amp, plot_v2_multifile, imdisp_polychromatic,
       plot_flux, plot_facility

# oichi2
export setup_dft, setup_dft_polychromatic, setup_nfft, setup_nfft_multiepochs,
       mod360, vis_to_v2, vis_to_t3, image_to_v2, image_to_t3phi, image_to_t3amp,
       observables, image_to_obs, image_to_vis_dft, image_to_vis, chi2_f, chi2_vis_nfft_f,
       chi2_vis_dft_fg, chi2_vis_nfft_fg, gaussian2d, cdg
export regularization, reg_centering, tvsq, tv, l1l2, l1l2w, l1hyp, l2sq,
       compactness, radial_variance, entropy, reg_support
export crit_fg, crit_f, crit_polychromatic_fg, chi2_fg, chi2_polychromatic_f,
       crit_multitemporal_fg, reconstruct, reconstruct_multitemporal,
       setup_radial_reg

# vis_functions
export bb, visibility_ud, visibility_ldpow, visibility_ldquad, visibility_ldquad_alt,
       visibility_ldlin, visibility_annulus, visibility_ellipse_quad,
       visibility_ellipse_uniform, visibility_thin_ring, visibility_Gaussian_ring,
       visibility_Gaussian_ring_az, visibility_ldsquareroot, visibility_Lorentzian_ring,
       visibility_GaussianLorentzian_ring_az
# model_chainrules — AD-compatible visibility functions with ChainRules rrules
export vis_ud, vis_ldlin, vis_ldquad, vis_ldpow,
       visibility_ud_d, visibility_ldlin_d, visibility_ldquad_d, visibility_ldpow_d

# resolvers — flat-dict parameter expression resolvers (internals not exported)

# hankel — numerical Hankel transform for radial profiles
export trapz, trapz_weights, hankel_transform, hankel_norm, hankel_vis,
       hankel_vis_interp, hankel_vis_full,
       eval_profile, compile_profile

# parse_model — compile flat param dicts into evaluable models
export FlatModel, AnalyticSpec, HankelSpec, AbstractComponentSpec,
       parse_model, eval_model, eval_model_grad

# chi2_flat — chi-squared for flat parametric models
export chi2_flat, chi2_flat_fg, residuals_flat, residuals_flat_jac

# fit_model — new flat-dict model fitting (NLopt + LsqFit + UltraNest)
export FitResult, LsqFitResult, UltraNestResult, display_model,
       fit_model, fit_model_lsqfit, fit_model_ultranest,
       model_to_obs, model_to_image, model_to_sed,
       resample_data

export get_uv, get_uv_indxes, prep_arrays, read_obs_file,
       read_comb_file, read_wave_file, simulate, simulate_from_oifits,
       get_v2_baselines, v2mapt3, get_t3_baselines, hour_angle_calc
export write_oi_header, write_oi_array, write_oi_target, write_oi_wavelength,
       write_oi_vis2, write_oi_t3
export setup_nfft_polychromatic, reconstruct_polychromatic
export optimize_sparco_parameters
export FacilityConfig, TargetConfig, CombinerConfig, WaveConfig,
       facility_info, obsv_info, combiner_info, wave_info,
       read_facility_file
export disk
export chi2_sparco_f, chi2_sparco_fg, reconstruct_sparco_gray

# simulate
export hours_to_date, sunrise_sunset, mjd_to_utdate, dates_to_jd, jd_to_hour_angle,
       opd_limits, alt_az, geometric_delay, cart_delay
export query_target_from_simbad, ra_dec_from_simbad, get_baselines, recenter
export gantt_onenight

# von Mises
export gaussianwrapped_to_vonmises_fast, logbesselI0

# oifitslib (wraps oifitslib C tools)
export oifits_check, oifits_merge, oifits_filter

# maximent / oimem — MaxEnt image reconstruction
export MaximENTParams, maxent_step!, reconstruct!,
       maxent_setup, maxent_reconstruct!, reconstruct_bsmem,
       auto_pixsize, gaussian_prior

# pmoired_compat — convert PMOIRED Python dicts to Julia
export pmoired_to_julia, pmoired_to_julia_file

end
