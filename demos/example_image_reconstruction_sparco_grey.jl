#
# Basic use of "SPARCO"-type reconstruction
#
using PyPlot
include("../src/OITOOLS.jl");using Main.OITOOLS
# using OITOOLS
oifitsfile = "/home/baron/Downloads/2019_v1295Aql.WL_SMOOTH.A.oifits"
pixsize = 0.1
nx = 128
data = readoifits(oifitsfile, filter_bad_data=true)[1,1];
ft = setup_nfft(data, nx, pixsize);
# Parameters
# 1: flux fraction of star at λ0
# 2: flux fraction of background at λ0
# 3: stellar angular diameter
# 4: spectral index of the environment+4
# 5: λ0
params_start=[0.46, 0, 0, 0, 1.6e-6] #V2:12.36 T3A: 4.14 T3P: 1.61
x_start = gaussian2d(nx,nx,nx/6);
x_start = vec(x_start)/sum(x_start);
#regularizers = [["l1l2", 7e8, 1e-3]];
regularizers = []
params, x = reconstruct_sparco_gray(x_start, params_start, data, ft, regularizers=regularizers, weights=[1.0, 0.0, 1.0], verb=true, maxiter=200); #grey environment
imdisp(x, pixscale = pixsize)

savefig("v1295_A.png")
writefits(x, "v1295_A.fits")
