#
# Pi Gruis, monochromatic reconstruction
#
using OITOOLS
# Activates the pure-Julia nested sampler that `fit_model_nested` runs. `using
# PythonCall` instead would select the UltraNest backend.
using Nautilus
oifitsfile = "./data/pigru.oifits"
pixsize = .5 # size of a pixel in milliarcseconds
nx = 64 # width of image (number of pixels)
data = readoifits(oifitsfile);
# Fourier transform setup
ft = setup_ft(data, nx, pixsize);

#initial image from model fitting
weights=[1.0,0.0,0.0]
model = Dict{String,Any}(
    "disc,ldquad" => 10.0,
    "disc,u1"     => 0.2,
    "disc,u2"     => 0.2,
    "disc,f"      => 1.0,
)
free_params = ["disc,ldquad", "disc,u1", "disc,u2"]
lb = Dict("disc,ldquad" => 5.0, "disc,u1" => -1.0, "disc,u2" => -1.0)
ub = Dict("disc,ldquad" => 30.0, "disc,u1" => 2.0, "disc,u2" => 1.0)

# Global minimization
result = fit_model_nested(model, free_params, data[1];
    lb=lb, ub=ub, weights=weights)
println(result)

x_start = model_to_image(result.model, result.x_opt, nx=nx, pixsize=pixsize)

image_to_chi2(x_start, ft, data, verb=true); # Evaluate chi2
regularizers = [["centering", 1e4], ["l1l2", 7e7, 1e-3]];
weights = [1.0, 0.0, 1.0]; # disable T3amp for VLTI/PIONIER
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = true, maxiter=500, weights=weights);
imdisp(x,pixsize=pixsize)
# Uncomment if you want to write the result
#writefits(x,"reconstruction.fits")
