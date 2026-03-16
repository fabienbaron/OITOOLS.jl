# Example to load one of NPOI files with multiple targets
using OITOOLS
filename = "./data/HD140573.oifits";
targetname =  "FKV0582";
data = (readoifits(filename, targetname=targetname, filter_bad_data = true, filter_v2_snr_threshold=1.0))[1,1];
v2plot(data)

# Fit a power-law limb-darkened disk
param = Dict{String,Any}(
    "star,ldpow" => 1.0,
    "star,alpha" => 0.2,
    "star,f"     => 1.0,
)
fit_params = ["star,ldpow", "star,alpha"]
lb = Dict("star,ldpow" => 0.1, "star,alpha" => 0.0)
ub = Dict("star,ldpow" => 40.0, "star,alpha" => 1.0)

result = fit_model(param, fit_params, data;
    lb=lb, ub=ub, weights=[1.0, 0, 0])
println(result)

obs = model_to_obs(result.model, result.x_opt, data)
v2plot_modelvsdata(data, obs.v2, logplot=true)
