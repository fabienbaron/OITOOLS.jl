# Model fitting example: multiple components
# Trying to model-fit the 2004 Beauty Contest data set
using OITOOLS
# Activates the pure-Julia nested sampler that `fit_model_nested` runs. `using
# PythonCall` instead would select the UltraNest backend.
using NestedSamplers

data = readoifits("./data/2004-data1.oifits")[1,1];

# Setup model: unresolved star + Gaussian ring.
#
# NOTE: this demo used to write "ring,gaussian_ring" => radius together with
# "ring,fwhm" => width. `gaussian_ring` has never been a geometry key in parse_model.jl, so
# the parser matched on `fwhm` instead and silently fitted a *plain Gaussian* — a different
# model from the one the comments described. The parser now warns about keys it ignores.
#
# The Gaussian ring OITOOLS implements is a difference of two Gaussians, parameterised by
# the inner and outer FWHM rather than by radius and width, so the numbers below are the
# equivalent of the old radius/width pair (fwhmin = 2R - W, fwhmout = 2R + W).
model_dict = Dict{String,Any}(
    "star,ud"       => 0.0,     # unresolved star (fixed at 0 diameter)
    "star,f"        => 0.5,     # star flux fraction (free)
    "ring,fwhmin"   => 1.5,     # inner FWHM, mas (free)
    "ring,fwhmout"  => 2.5,     # outer FWHM, mas (free)
    "ring,incl"     => 0.0,     # inclination (fixed for now)
    "ring,pa"       => 0.0,     # position angle (fixed for now)
    "ring,f"        => "1.0 - \$star,f", # ring flux = complement of star flux
)

list_free_params = ["star,f", "ring,fwhmin", "ring,fwhmout"]

lb = Dict("star,f" => 0.0, "ring,fwhmin" => 0.1, "ring,fwhmout" => 0.2)
ub = Dict("star,f" => 1.0, "ring,fwhmin" => 10.0, "ring,fwhmout" => 20.0)

display_model(model_dict, list_free_params; lb=lb, ub=ub)

# UltraNest nested sampling for robust posterior exploration
result = fit_model_nested(model_dict, list_free_params, data;
    lb=lb, ub=ub,
    min_num_live_points=400, cluster_num_live_points=200)

println(result)
