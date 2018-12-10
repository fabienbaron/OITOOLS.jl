#
# Image reconstruction code with temporal regularization
#
using OITOOLS
using FITSIO
oifitsfiles=["./data/AZCYGJUN01.oifits","./data/AZCYGJUL.oifits","./data/AZCYGAUG.oifits"]
printcolors = [:red, :green, :blue]
nepochs, tepochs, data = readoifits_multiepochs(oifitsfiles, filter_bad_data=true, force_full_t3 = true);
pixsize = 0.1;
nx = 64;
fftplan = setup_nfft_multiepochs(data, nx, pixsize);;

#initial image is a simple Gaussian
prior_im=read(FITS("./data/4masdisk.fits")[1])
prior_im=convert(Array{Float64,2},prior_im)
prior_im=vec(prior_im)/sum(prior_im)
x_start=prior_im=repeat(prior_im,nepochs)


regularizers = [   [ ["centering", 1e4], ["tv", 2e3] ],  #epoch 1
                   [ ["tv", 5e2] ],  #epoch 2
                   [ ["tv", 5e2] ], #epoch 3
                   [ ["temporal_tvsq", 3e4] ] ]; #transtemporal

x = reconstruct_multitemporal(x_start, data, fftplan, printcolor= printcolors, regularizers = regularizers);

imdisp_temporal(x, nepochs, pixscale=pixsize, cmap="gist_heat")
writefits(reshape(x,nx,nx,3),"reconstruction.fits")
