#
# BOOTSTRAP EXAMPLE
#

include("oitools.jl");
oifitsfile = "AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1];

f_chi2, params_mode, cvis_model = fit_model_v2(data, visibility_ud, [8.0]);# diameter is the parameter
println("Nominal value: $(params_mode)")
nbootstraps = 10000
params = zeros(Float64,nbootstraps)
println("Boostraping...")
for k=1:nbootstraps
 f_chi2, paropt, ~ = fit_model_v2(bootstrap_v2_data(data), visibility_ud, [8.0], verbose = false, calculate_vis = false);# diameter is the parameter
 params[k]=paropt[1];
end
h = plt[:hist](params,50);
params_err = std(params, corrected=false)
println("Error bar: $(params_err)")
