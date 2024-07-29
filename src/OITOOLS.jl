module OITOOLS
include("readoifits.jl");
include("vis_functions.jl");
include("modelfit.jl");
include("write_oifits_ha.jl");
include("write_oifits_obs.jl");
include("utils.jl");
include("oichi2.jl");
include("oiplot.jl");
include("astrometry.jl");
include("vonmises.jl")
include("simulate.jl");
include("oifitslib.jl");

#readoifits
export OIdata
export readoifits, readoifits_multiepochs, readfits, writefits
export oifits_prep, updatefits_aspro,readoifits_multicolors, list_oifits_targets
export remove_redundant_uv!,filter_data,set_data_filter

#modelfit
export pos_fixed,spectrum_powerlaw,spectrum_gray,model_to_image, limbdarkened_disk
export OIparam, OIcomponent,OImodel
export create_component,create_model,update_model, model_to_vis, dispatch_params, model_to_obs, model_to_chi2,visfunc_to_chi2, get_model_bounds,get_model_params,get_model_pnames,fit_model_ultranest, fit_model_levenberg, fit_model_nlopt, fit_visfunc_nlopt, resample_data, bootstrap_fit

#oiplot
export set_oiplot_defaults, uvplot, onclickidentify, plot_v2,plot_v2_timelapse, plot_diffphi, plot_visphi, plot_t3phi, plot_v2_and_t3phi_wav, imdisp, imdisp_temporal, plot_v2_residuals, plot_t3phi_residuals, plot_t3amp_residuals, plot_v2_model_vs_func, plot_v2, plot_t3amp,plot_v2_multifile,imdisp_polychromatic
#oichi2
export setup_dft, setup_nfft,setup_nfft_multiepochs, mod360, vis_to_v2,vis_to_t3, image_to_v2, image_to_t3phi, image_to_t3amp, image_to_obs, image_to_vis_dft,image_to_vis,chi2_f,chi2_f,chi2_vis_nfft_f,chi2_vis_dft_fg,chi2_vis_nfft_fg,gaussian2d,cdg
export regularization, reg_centering,tvsq,tv,l1l2,l1l2w,l1hyp,l2sq,compactness,radial_variance,entropy, reg_support
export crit_fg,crit_f, crit_polychromatic_fg, chi2_fg,chi2_f,chi2_polychromatic_f,crit_multitemporal_fg,reconstruct,reconstruct_multitemporal, setup_radial_reg
#vis_functions
export bb, visibility_ud, visibility_ldpow, visibility_ldquad, visibility_ldquad_alt, visibility_ldlin, visibility_annulus, visibility_ellipse_quad, visibility_ellipse_uniform,visibility_thin_ring,visibility_Gaussian_ring, visibility_Gaussian_ring_az,visibility_ldsquareroot, visibility_Lorentzian_ring, visibility_GaussianLorentzian_ring_az
export get_uv,get_uv_indxes,prep_arrays,read_array_file,read_obs_file,read_comb_file,read_wave_file,simulate,simulate_from_oifits,vis_to_t3_conj,get_v2_baselines,v2mapt3,get_t3_baselines,hour_angle_calc
export write_oi_header,write_oi_array,write_oi_target,write_oi_wavelength,write_oi_vis2,write_oi_t3
export setup_nfft_polychromatic,imdisp_polychromatic,reconstruct_polychromatic, image_to_vis_polychromatic_nfft,chi2_polychromatic_nfft_f
export optimize_sparco_parameters
export facility_info,obsv_info,combiner_info,wave_info,error_struct,read_facility_file,read_obs_file,read_wave_file,read_comb_file,define_errors
export disk, setup_nfft_t4, vis_to_t4
export chi2_sparco_f, chi2_sparco_f_alt, chi2_sparco_fg,reconstruct_sparco_gray
# simulate
export facility_info, obsv_info, combiner_info, wave_info, error_struct
export hours_to_date, sunrise_sunset, hour_angle_calc, mjd_to_utdate,dates_to_jd,jd_to_hour_angle,opd_limits,alt_az,geometric_delay,cart_delay
export query_target_from_simbad, ra_dec_from_simbad, get_baselines, recenter
export gantt_onenight
#von mises
export gaussianwrapped_to_vonmises_fast,logbesselI0

#oifitslib
export oifits_check, oifits_merge, oifits_filter
end
