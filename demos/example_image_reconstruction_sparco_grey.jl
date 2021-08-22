#
# Very Basic Image reconstruction code
#
include("../src/OITOOLS.jl");using Main.OITOOLS
# using OITOOLS
oifitsfile = "./data/MWC480.oifits"
pixsize = 0.15
nx = 64
data = readoifits(oifitsfile, filter_bad_data=true)[1,1];
ft = setup_nfft(data, nx, pixsize);
# Parameters
# 1: flux fraction of star at λ0
# 2: flux fraction of background at λ0
# 3: stellar angular diameter
# 4: spectral index of the enrironment
# 5: λ0
params=[0.5, 0.1, 0.4, 1.0, 1.6e-6] #V2:12.36 T3A: 4.14 T3P: 1.61
x_start = gaussian2d(nx,nx,nx/6);
x_start = vec(x_start)/sum(x_start);
δ=1e-6
numerical_g = zeros(4)
for i=1:4
    orig = params[i]
    params[i] = orig + 2*δ
    f2r = chi2_sparco_nfft_f(x_start, ft, data, params,verb=false)
    params[i] = orig + δ
    f1r = chi2_sparco_nfft_f(x_start, ft, data, params,verb=false)
    params[i] = orig - δ
    f1l = chi2_sparco_nfft_f(x_start, ft, data, params,verb=false)
    params[i] = orig - 2*δ
    f2l = chi2_sparco_nfft_f(x_start, ft, data, params,verb=false)
    numerical_g[i] = (f2l-f2r+8*(f1r-f1l))
    params[i] = orig
end
numerical_g ./= (12*δ)

start=[params;x_start];
g=zeros(Float64, length(params)+nx*nx);
chi2_sparco_nfft_fg(start, g, ft, data, length(params),verb=false);
g[1:4]

sol = reconstruct_sparco_gray(start, params, data, ft, verb=true); #grey environment
