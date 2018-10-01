#
# Image reconstruction code with temporal regularization
#
using OITOOLS
oifitsfiles = ["./data/AZCYG_JUN01_2018.oifits","./data/AZCYG_JUL_2018.oifits","./data/AZCYG_AUG25_2018.oifits"]
printcolors = [:red, :green, :blue];
nepochs, tepochs, data = readoifits_multiepochs(oifitsfiles, filter_bad_data=true, force_full_t3 = true);
pixsize = 0.2;
nx = 64;
fftplan = setup_nfft_multiepochs(data, nx, pixsize);

#initial image is a simple Gaussian
x_start_single = [ exp(-((i-(nx+1)/2)^2+(j-(nx+1)/2)^2)/(2*(nx/6)^2)) for i=1:nx for j=1:nx]
x_start_single /= sum(x_start_single);
x_start = repmat(x_start_single, nepochs)

regularizers = [   [ ["centering", 1e4], ["tv", 7e3] ],  #epoch 1
                   [ ["tv", 7e3] ],  #epoch 2
                   [ ["tv", 7e3] ], #epoch 3
                   [ ["temporal_tvsq", 1e2] ] ]; #transtemporal

x = reconstruct_multitemporal(x_start, data, fftplan, printcolor= printcolors, regularizers = regularizers);

imdisp_temporal(x, nepochs, pixscale=pixsize)
writefits(reshape(x,nx,nx),"reconstruction.fits")
