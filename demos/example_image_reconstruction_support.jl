#
# Image reconstruction with soft support constraint
# Example with Polaris MIRCX data
#
using OITOOLS
using PythonPlot
set_oiplot_defaults();
oifitsfile = "./data/polaris.oifits"
data = readoifits(oifitsfile);

# Let's find an adequate mask by fitting a UD to the data
weights=[1.0,1.0,0.0]
model = Dict{String,Any}(
    "disc,ud" => 2.0,
    "disc,f"  => 1.0,
)
free_params = ["disc,ud"]
lb = Dict("disc,ud" => 1.0)
ub = Dict("disc,ud" => 5.0)

# Global minimization via UltraNest
result_un = fit_model_ultranest(model, free_params, data[1];
    lb=lb, ub=ub, weights=weights)
# Local refinement
result = fit_model(model, free_params, data[1];
    lb=lb, ub=ub, weights=weights)
println(result)

# Now imaging
pixsize = 0.05 # size of a pixel in milliarcseconds
nx = 128 # width of image (number of pixels)
weights=[1.0,1.0,1.0]
mask  = disk(npix=nx, diameter=result.x_opt[1]/pixsize+1) # binary mask
prior = model_to_image(result.model, result.x_opt, nx=nx, pixsize=pixsize, oversample=1).*mask
regularizers = [["centering", 1e2], ["l1l2", 2e8, 1e-4], ["support", 1.0, mask]];
ft = setup_ft(data, nx, pixsize);
x = reconstruct(prior, data, ft, regularizers = regularizers, verb = true, maxiter=500, weights=weights);
x = reconstruct(x.*mask, data, ft, regularizers = regularizers, verb = true, maxiter=500, weights=weights);
chi2 = image_to_chi2(x, ft, data, weights=weights, verb=true);
imdisp(x, pixsize=pixsize, beamsize=0.5*1/maximum(data[1].uv_baseline)*180/pi*3600*1000);
savefig("polaris_image.png")
