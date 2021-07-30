__precompile__(true)

module OITOOLS
include("readoifits.jl");
include("write_oifits_ha.jl");
include("write_oifits_obs.jl");
include("utils.jl");
include("oichi2.jl");
include("oiplot.jl");
include("vis_functions.jl");
include("modelfit.jl");
include("astrometry.jl");
include("vonmises.jl")
include("simulate.jl");

#readoifits
export OIdata
export readoifits, readoifits_multiepochs, readfits, writefits
export oifits_prep, updatefits_aspro,readoifits_multicolors
export remove_redundant_uv!,filter_data,set_data_filter
precompile(readoifits, (String,))


#modelfit
export pos_fixed,spectrum_powerlaw,spectrum_gray
export OIparam, OIcomponent,OImodel
export create_component,create_model,update_model, model_to_cvis,dispatch_params,model_to_chi2,visfunc_to_chi2, get_model_bounds,get_model_params,get_model_pnames,fit_model_ultranest, fit_model_levenberg, fit_model_nlopt, resample_data, bootstrap_fit

#oiplot
export uvplot, onclickidentify, v2plot,v2plot_timelapse, diffphiplot, visphiplot, t3phiplot, imdisp, imdisp_temporal, v2plot_modelvsdata, v2plot_modelvsfunc,v2plot,t3phiplot,t3ampplot,v2plot_multifile,imdisp_multiwave,imdisp_polychromatic
#oichi2
export setup_dft, setup_nfft,setup_nfft_multiepochs, mod360, cvis_to_v2,cvis_to_t3, image_to_cvis_dft,image_to_cvis_nfft,chi2_dft_f,chi2_nfft_f,chi2_vis_nfft_f,chi2_vis_dft_fg,chi2_vis_nfft_fg,gaussian2d,cdg,reg_centering,tvsq,tv,regularization,crit_dft_fg,crit_nfft_fg, chi2_dft_fg,chi2_nfft_fg,crit_nfft_fg,crit_multitemporal_nfft_fg,reconstruct,reconstruct_multitemporal
#vis_functions
export visibility_ud, visibility_ldpow, visibility_ldquad, visibility_ldlin, visibility_annulus, visibility_ellipse_quad, visibility_ellipse_uniform,visibility_thin_ring,visibility_Gaussian_ring, visibility_Gaussian_ring_az,visibility_ldsquareroot, visibility_Lorentzian_ring, visibility_GaussianLorentzian_ring_az
export get_uv,get_uv_indxes,prep_arrays,read_array_file,read_obs_file,read_comb_file,read_wave_file,simulate,simulate_from_oifits,cvis_to_t3_conj,get_v2_baselines,v2mapt3,get_t3_baselines,hour_angle_calc
export write_oi_header,write_oi_array,write_oi_target,write_oi_wavelength,write_oi_vis2,write_oi_t3
export setup_nfft_polychromatic,imdisp_polychromatic,reconstruct_polychromatic, image_to_cvis_polychromatic_nfft,chi2_polychromatic_nfft_f
export facility_info,obsv_info,combiner_info,wave_info,error_struct,read_facility_file,read_obs_file,read_wave_file,read_comb_file,define_errors
export setup_nfft_t4, cvis_to_t4
export chi2_sparco_nfft_f,chi2_sparco_nfft_f_alt, chi2_sparco_nfft_fg,reconstruct_sparco_gray
# simulate
export facility_info, obsv_info, combiner_info, wave_info, error_struct
export hours_to_date, sunrise_sunset, hour_angle_calc, mjd_to_utdate,dates_to_jd,jd_to_hour_angle,opd_limits,alt_az,geometric_delay,cart_delay
export query_target_from_simbad, ra_dec_from_simbad, get_baselines, recenter
export gantt_onenight
#von mises
export gaussianwrapped_to_vonmises_fast,logbesselI0
end
