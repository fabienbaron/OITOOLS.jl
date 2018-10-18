#
# Image reconstruction code with temporal regularization
#
using OITOOLS
using FITSIO
oifitsfile=["./data/spotchange1.oifits","./data/spotchange2.oifits","./data/spotchange3.oifits"]
printcolors = [:red, :green, :blue]
nepochs, tepochs, data = readoifits_multiepochs(oifitsfile, filter_bad_data=true, force_full_t3 = true);
pixsize = 0.1;
nx = 64;
fftplan = setup_nfft_multiepochs(data, nx, pixsize);;


#initial image is a simple Gaussian
#x_start_single = [ exp(-((i-(nx+1)/2)^2+(j-(nx+1)/2)^2)/(2*(nx/6)^2)) for i=1:nx for j=1:nx]
#x_start_single /= sum(x_start_single);
#x_start = repeat(x_start_single, nepochs);

start_im=read(FITS("./data/4masdisk.fits")[1])
start_im=convert(Array{Float64,2},start_im)
start_im=vec(start_im)/sum(start_im)
x_start=prior_im=repeat(start_im,nepochs)


regularizers = [   [ ["centering", 1e4], ["tv", 2e4] ],  #epoch 1
                   [ ["tv", 2e4] ],  #epoch 2
                   [ ["tv", 2e4] ], #epoch 3
                   [ ["temporal_tvsq", 3e8] ] ]; #transtemporal

x = reconstruct_multitemporal(x_start, data, fftplan, printcolor= printcolors, regularizers = regularizers);

imdisp_temporal(x, nepochs, pixscale=pixsize)

source_file1 = "./data/4masdisk_spotchange1.fits";
source_image1=readfits(source_file1);  source_image1=vec(source_image1);
source_file2 = "./data/4masdisk_spotchange2.fits";
source_image2=readfits(source_file2); source_image2=vec(source_image2);
source_file3 = "./data/4masdisk_spotchange3.fits";
source_image3=readfits(source_file3); source_image3=vec(source_image3);
source_images=vcat(convert(Array{Float64,1},source_image1),convert(Array{Float64,1},source_image2),convert(Array{Float64,1},source_image3))

writefits(reshape(x,nx,nx,3),"reconstruction.fits")
