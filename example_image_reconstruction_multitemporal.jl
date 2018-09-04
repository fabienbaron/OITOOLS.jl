#
# Very Basic Image reconstruction code
#
include("oitools.jl");
oifitsfiles = ["AZCYG_JUN01_2018.oifits","AZCYG_JUL_2018.oifits","AZCYG_AUG25_2018.oifits"]
epochs_weights = [1.0, 1.0, 1.0];
nepochs, tepochs, data = readoifits_multiepochs(oifitsfiles, filter_bad_data=true, force_full_t3 = true);
pixsize = 0.2;
nx = 64;
fftplan_multi = setup_nfft_multiepochs(data, nx, pixsize);

#initial image is a simple Gaussian
x_start_single = [ exp(-((i-(nx+1)/2)^2+(j-(nx+1)/2)^2)/(2*(nx/6)^2)) for i=1:nx for j=1:nx]
x_start_single /= sum(x_start_single);
x_start = repmat(x_start_single, 3)


mu = 1e4
g = zeros(Float64, length(x_start))
chi2_multitemporal_tv_fg(x_start, g, mu, epochs_weights, fftplan_multi, data, verb = true);

imdisp_temporal(x, nepochs, pixscale=pixsize)
writefits(reshape(x,nx,nx),"reconstruction.fits")
