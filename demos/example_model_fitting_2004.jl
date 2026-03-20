# Model fitting example: multiple components
# Trying to model-fit the 2004 Beauty Contest data set
using OITOOLS

data = readoifits("./data/2004-data1.oifits")[1,1];

# Setup model: unresolved star + Gaussian ring
model_dict = Dict{String,Any}(
    "star,ud"       => 0.0,     # unresolved star (fixed at 0 diameter)
    "star,f"        => 0.5,     # star flux fraction (free)
    "ring,gaussian_ring" => 1.0, # ring radius (free)
    "ring,fwhm"     => 0.5,     # ring FWHM (free)
    "ring,incl"     => 0.0,     # inclination (fixed for now)
    "ring,pa"       => 0.0,     # position angle (fixed for now)
    "ring,f"        => "1.0 - \$star,f", # ring flux = complement of star flux
)

list_free_params = ["star,f", "ring,gaussian_ring", "ring,fwhm"]

lb = Dict("star,f" => 0.0, "ring,gaussian_ring" => 0.1, "ring,fwhm" => 0.01)
ub = Dict("star,f" => 1.0, "ring,gaussian_ring" => 10.0, "ring,fwhm" => 5.0)

display_model(model_dict, list_free_params; lb=lb, ub=ub)

# UltraNest nested sampling for robust posterior exploration
result = fit_model_ultranest(model_dict, list_free_params, data;
    lb=lb, ub=ub,
    min_num_live_points=400, cluster_num_live_points=200)

println(result)
