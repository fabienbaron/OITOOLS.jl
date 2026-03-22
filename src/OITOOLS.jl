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

version() = println("OITOOLS v$(pkgversion(OITOOLS))")

include("readoifits.jl")
include("vis_functions.jl")
include("model_chainrules.jl")

# ── Flat-dict model fitting stack ────────────────────────────────────────────
include("resolvers.jl")          # SharedUtils, RGF, HandRolled, Symbolic
include("hankel.jl")             # Hankel transform + chainrules + profile compiler
include("parse_model.jl")        # FlatModel, eval_model, eval_model_grad
include("chi2_flat.jl")          # chi2_flat, chi2_flat_fg
include("fit_model.jl")          # fit_model (NLopt), fit_model_ultranest, display_model


include("utils.jl")
include("oichi2.jl")
include("sparco_flat.jl")
include("oiplot.jl")
include("astrometry.jl")
include("vonmises.jl")
include("simulate.jl")
include("oifitslib.jl")
include("maximent.jl")
include("oimem.jl")
include("plan.jl")
include("pmoired_compat.jl")

# ── Reading OIFITS data ─────────────────────────────────────────────────────
export OIdata
export readoifits, readoifits_multiepochs, list_oifits_targets
export readfits, writefits
export oifits_prep, updatefits_aspro
export filter_data, set_data_filter

# ── Plotting ─────────────────────────────────────────────────────────────────
export set_oiplot_defaults, uvplot, plot_v2, plot_t3phi, plot_t3amp,
       plot_visamp, plot_visphi, plot_diffphi, plot_flux, plot_obs, plot_multi,
       plot_v2_residuals, plot_t3phi_residuals, plot_t3amp_residuals,
       plot_visamp_residuals, plot_visphi_residuals, plot_residuals,
       plot_v2_multifile,
       imdisp, imdisp_multi, plot_facility

# ── Model fitting (flat-dict interface) ──────────────────────────────────────
export FlatModel, dict_to_model, parse_model, eval_model, eval_model_grad, display_model
export fit_model, fit_model_lsqfit, fit_model_ultranest
export FitResult, LsqFitResult, UltraNestResult
export model_to_vis, model_to_obs, model_to_residuals, model_to_chi2, model_to_chi2_fg,
       model_to_image, model_to_sed, resample_data

# ── Visibility functions ─────────────────────────────────────────────────────
export visibility_ud, visibility_ldlin, visibility_ldquad, visibility_ldquad_alt,
       visibility_ldpow, visibility_ldsquareroot,
       visibility_annulus, visibility_ellipse_uniform, visibility_ellipse_quad,
       visibility_thin_ring, visibility_Gaussian_ring, visibility_Gaussian_ring_az,
       visibility_Lorentzian_ring, visibility_GaussianLorentzian_ring_az

# ── Image reconstruction ─────────────────────────────────────────────────────
export setup_ft, setup_dft, setup_nfft
export image_to_vis,
       image_to_v2, image_to_t3phi, image_to_t3amp, image_to_obs,
       image_to_residuals, image_to_chi2, image_to_chi2_fg
export vis_to_v2, vis_to_t3, observables
export crit_f, crit_fg, image_to_crit
export reconstruct
export gaussian2d, recenter

# ── SPARCO ──────────────────────────────────────────────────────────────────
export reconstruct_sparco_gray, reconstruct_sparco_multi
export chi2_sparco_multi_f, optimize_sparco_multi_parameters
export reconstruct_sparco_flat, chi2_sparco_flat_f, optimize_sparco_flat_parameters

# ── Deprecated (still exported for backward compatibility) ──────────────────
export setup_dft_polychromatic, setup_nfft_multiepochs, setup_nfft_polychromatic
export chi2_f, chi2_fg, chi2_polychromatic_f
export crit_polychromatic_fg, crit_multitemporal_fg
export reconstruct_multitemporal, reconstruct_polychromatic

# ── Maximum entropy (BSMEM) ──────────────────────────────────────────────────
export reconstruct_bsmem
export auto_pixsize, gaussian_prior

# ── Simulation ───────────────────────────────────────────────────────────────
export simulate, simulate_from_oifits, get_uv
export read_facility_file, read_obs_file, read_comb_file, read_wave_file
export FacilityConfig, TargetConfig, CombinerConfig, WaveConfig
export gantt_onenight
export sunrise_sunset, alt_az, opd_limits
export query_target_from_simbad, ra_dec_from_simbad, magnitudes_from_simbad

# ── Observation planning (CHARA) ──────────────────────────────────────────────
export night_observability, compute_delays, in_delay, best_pop
export obs_plan, chara_plan, print_pop_results
export moon_radec, moon_illumination, angular_separation

# ── Writing OIFITS ───────────────────────────────────────────────────────────
export oifits_check, oifits_merge, oifits_filter

# ── PMOIRED compatibility ────────────────────────────────────────────────────
export pmoired_to_julia, pmoired_to_julia_file, pmoired_to_dict

end
