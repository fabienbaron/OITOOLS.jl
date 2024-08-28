# Noura Paper "SPARCO"-type reconstruction
#

using OITOOLS
oifitsfile = "./data/2019_v1295Aql.WL_SMOOTH.A.oifits"
pixsize = 0.25 #mas/pixel
nx = 32 # pixels
data = readoifits(oifitsfile, filter_bad_data=true)[1,1];
ft = setup_nfft(data, nx, pixsize);
params_start = [0.444, 0, 0, 4, 1.6e-6] #V2:12.36 T3A: 4.14 T3P: 1.61
x_start = rand(nx,nx)
x_start /= sum(x_start)
chi2_sparco_f(x_start,  params_start, ft, data, weights=[1.0,0.0,1.0])

# Initial image, no regularization except positivity
regularizers = []
params_start = [0.5, 0, 0, 4, 1.6e-6] 
params, x = reconstruct_sparco_gray(x_start, params_start, data, ft, regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=200); #grey environment
minchi2, params,ret = optimize_sparco_parameters(params_start, x, ft, data; weights = [1.0,1.0,1.0] )
params, x = reconstruct_sparco_gray(x, params, data, ft, regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=200); #grey environment
imdisp(x, pixsize = pixsize, figtitle="v1295 Aql - A")
regularizers = [["tvsq", 1e7]]
params, x = reconstruct_sparco_gray(x, params, data, ft, regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=200); #grey environment
minchi2, params,ret = optimize_sparco_parameters(params_start, x, ft, data; weights = [1.0,1.0,1.0] )
regularizers = [["tvsq", 1e5]]
params, x = reconstruct_sparco_gray(x, params, data, ft, regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=200); #grey environment
minchi2, params,ret = optimize_sparco_parameters(params_start, x, ft, data; weights = [1.0,1.0,1.0] )
regularizers = [["tvsq", 1e3]]
params, x = reconstruct_sparco_gray(x, params, data, ft, regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=200); #grey environment
minchi2, params,ret = optimize_sparco_parameters(params_start, x, ft, data; weights = [1.0,1.0,1.0] )
regularizers = [["tvsq", 1e1]]
params, x = reconstruct_sparco_gray(x, params, data, ft, regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=200); #grey environment
minchi2, params,ret = optimize_sparco_parameters(params_start, x, ft, data; weights = [1.0,1.0,1.0] )
imdisp(x, pixsize = pixsize, figtitle="v1295 Aql - A")
g=similar([params_start;x_start])
chi2_f = chi2_sparco_fg([params_start;x_start],  g, ft, data, 5, weights=[1.0, 0.0, 1.0])