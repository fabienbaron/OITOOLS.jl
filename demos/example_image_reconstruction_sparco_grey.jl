#
# Very Basic Image reconstruction code
#
#include("../src/OITOOLS.jl");using Main.OITOOLS
using OITOOLS
oifitsfile = "./data/MWC480.oifits"
pixsize = 0.15
nx = 64
data = readoifits(oifitsfile, filter_bad_data=true)[1,1];
ft = setup_nfft(data, nx, pixsize);

params=[0.5, 0.4, 1.6e-6, 0.0] #V2:12.36 T3A: 4.14 T3P: 1.61
#x_start=vec(readfits("../../squeeze/mwc480_grey.fits")); 0.1mas/pix

x_start = gaussian2d(nx,nx,nx/6);
x_start = vec(x_start)/sum(x_start);

chi2_sparco_nfft_f_alt(x_start, ft, data, params)
start=[x_start;params[1]; params[2]];
chi2_sparco_nfft_f(start, ft, data, params)
#g=zeros(nx*nx+2);
#chi2_sparco_nfft_fg(start, g, ft, data, params)
sol = reconstruct_sparco_gray(start, params, data, ft, verb=true); #grey environment
x_sol=sol[1:end-2];
params[1] = sol[end-1]; params[2] = sol[end]
imdisp(x_sol)
