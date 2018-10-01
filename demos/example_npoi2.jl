# Example to load one of NPOI files with multiple targets
using OITOOLS
filename = "./data/HD140573.oifits";
targetname =  "FKV0582";
data = (readoifits(filename, targetname=targetname, filter_bad_data = true, filter_v2_snr_threshold=1.0))[1,1];
v2plot(data.v2_baseline,data.v2,data.v2_err)
f_chi2, params, cvis_model = fit_model_v2(data, visibility_ud, [1.0]);# diameter is the parameter
v2_model = cvis_to_v2(cvis_model, data.indx_v2)
v2plot_modelvsdata(data.v2_baseline,data.v2,data.v2_err, v2_model,logplot=true);
v2plot_modelvsfunc(data.v2_baseline,data.v2,data.v2_err, visibility_ud,params , yrange=[1e-4, 1.0], logplot=true);
