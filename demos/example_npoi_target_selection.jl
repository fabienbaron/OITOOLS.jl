# Example to load one of NPOI files with multiple targets
using OITOOLS
filename = "./data/HD140573.oifits";
targetname =  "FKV0582";
data = (readoifits(filename, targetname=targetname, filter_bad_data = true, filter_v2_snr_threshold=1.0))[1,1];
v2plot(data)

model = create_model(create_component(type="ldpow", name="Model"));

minf, minx, cvis_model, result = fit_model_nlopt(data, model, weights=[1.0,0,0]);

v2_model = vis_to_v2(cvis_model, data.indx_v2)
v2plot_modelvsdata(data, v2_model,logplot=true);
v2plot_modelvsfunc(data, visibility_ud,params , yrange=[1e-4, .1]);
