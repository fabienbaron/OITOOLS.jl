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

# Frank
using OITOOLS, DelimitedFiles
oifitsfile = "/home/baron/Downloads/2019_v1295Aql.WL_SMOOTH.A.oifits"
data = readoifits(oifitsfile, filter_bad_data=true, use_vis=false)[1,1];
# Remove polychromatic point source contribution from data (with SPARCO squeeze, OITOOLS modeling, etc.)
...
model = create_model(create_component(type="ud-polychromatic", name="Star"));
dispatch_params([8.0,0.1], model);
cvis_star = model_to_cvis(model, data)
v2_star = cvis_to_v2(cvis_star, data.indx_v2)

# Convert V2 to V with "correct" error prop -- based off Sivia's Bayesian Tutorial  3.98a/3.98b
VV = 0.5*( data.v2 +  sqrt.(data.v2.^2+2*data.v2_err.^2))
Vamp = sqrt.(VV)
#Vamp_err = sqrt.(1.0./(1.0./VV+2*(3*VV-data.v2)./(data.v2_err.^2)))  # = σ
Re = Vamp
Im = 0*Re
weights = 1.0./VV+2*(3*VV-data.v2)./(data.v2_err.^2)
u = data.uv[1,data.indx_v2]
v = data.uv[2,data.indx_v2]
open("frank.txt"; write=true) do f
    #write(f, "u v Re Im weights\n")
    writedlm(f, hcat([0;u], [0;v], [1.0;Re], [0.0;Im], [10000.;weights]))
end
