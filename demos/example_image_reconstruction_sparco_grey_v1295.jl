# Noura Paper "SPARCO"-type reconstruction
#

using OITOOLS
oifitsfile = "./data/2019_v1295Aql.WL_SMOOTH.A.oifits"
pixsize = 0.125 #mas/pixel
nx = 64 # pixels
data = readoifits(oifitsfile, filter_bad_data=true);
ft = setup_ft(data, nx, pixsize);
params_start = [0.444, 0, 0, 4, 1.6e-6] #V2:12.36 T3A: 4.14 T3P: 1.61
x_start = rand(nx,nx)
x_start /= sum(x_start)
weights=[1.0, 0.0, 1.0] # For this dataset, the T3amp are a bit fishy
chi2_sparco_f(x_start,  params_start, ft[1], data[1], weights=weights)

# Initial image, no regularization except positivity
regularizers = []
params_start = [0.5, 0, 0, 4, 1.6e-6] 
params, x = reconstruct_sparco_gray(x_start, params_start, data[1], ft[1], regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=200); #grey environment
minchi2, params,ret = optimize_sparco_parameters(params_start, x, ft[1], data[1]; weights = weights )
params, x = reconstruct_sparco_gray(x, params, data[1], ft[1], regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=200); #grey environment
imdisp(x, pixsize = pixsize, figtitle="v1295 Aql - A")
# The idea is to slowly decrease the regularization level
regularizers = [["tvsq", 1e7]]
params, x = reconstruct_sparco_gray(x, params, data[1], ft[1], regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=200); #grey environment
minchi2, params,ret = optimize_sparco_parameters(params_start, x, ft[1], data[1]; weights = weights )
regularizers = [["tvsq", 1e5]]
params, x = reconstruct_sparco_gray(x, params, data[1], ft[1], regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=200); #grey environment
minchi2, params,ret = optimize_sparco_parameters(params_start, x, ft[1], data[1]; weights = weights )
regularizers = [["tvsq", 1e3]]
params, x = reconstruct_sparco_gray(x, params, data[1], ft[1], regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=200); #grey environment
minchi2, params,ret = optimize_sparco_parameters(params_start, x, ft[1], data[1]; weights = weights )
regularizers = [["tvsq", 1e1]]
params, x = reconstruct_sparco_gray(x, params, data[1], ft[1], regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=200); #grey environment
minchi2, params,ret = optimize_sparco_parameters(params_start, x, ft[1], data[1]; weights = weights )
imdisp(x, pixsize = pixsize, figtitle="v1295 Aql - A")
chi2_sparco_f(x,  params, ft[1], data[1], weights=weights)
