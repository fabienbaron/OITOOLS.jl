#
# BOOTSTRAP EXAMPLE
#
using OITOOLS
using Statistics

# VLTI example
oifitsfile = "./data/AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1];

# Model: linear limb-darkened disk
param = Dict{String,Any}(
    "star,ldlin" => 8.0,
    "star,u"     => 0.3,
    "star,f"     => 1.0,
)
fit_params = ["star,ldlin", "star,u"]
lb = Dict("star,ldlin" => 5.0, "star,u" => 0.0)
ub = Dict("star,ldlin" => 12.0, "star,u" => 1.0)

# First fit to get best-fit values
result = fit_model(param, fit_params, data;
    lb=lb, ub=ub, weights=[1.0, 0.0, 0.0])
println("Best fit: ", result.x_opt, "  χ²r = ", result.chi2r)

# Bootstrap resampling
nbootstraps = 1000
npar = length(fit_params)
boot_params = zeros(nbootstraps, npar)

for i in 1:nbootstraps
    # Resample data with replacement
    boot_data = resample_data(data)
    # Re-fit on bootstrapped data
    r = fit_model(param, fit_params, boot_data;
        lb=lb, ub=ub, weights=[1.0, 0.0, 0.0],
        method=:LN_NELDERMEAD, maxeval=3000, verb=false)
    boot_params[i, :] = r.x_opt
end

# Summarise
for j in 1:npar
    m = mean(boot_params[:, j])
    s = std(boot_params[:, j])
    println("$(fit_params[j]):  mean = $(round(m, digits=4)) ± $(round(s, digits=4))")
end
