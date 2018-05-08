# Example to load one of NPOI files with multiple targets
include("oitools.jl");

filename = "HD48329_oidata.fits";
targetname =  "FKV0254";
data = (readoifits_filter(filename, targetname=targetname))[1,1];

# Then we need to filter bad data
good = find(  (data.v2_err.>0) .& (data.v2_err.<1.0) .& (data.v2.>-0.2) .& (data.v2.<1.2) )
data.v2 = data.v2[good]
data.v2_err = data.v2_err[good]
data.v2_baseline = data.v2_baseline[good]
data.nv2 = length(data.v2)
v2plot(data,logplot=true);

f_chi2, params, cvis_model = fit_model_v2(data, visibility_ud, [1.0]);# diameter is the parameter
