#
# Image reconstruction code with spectral regularization
# This requires a private example file at the moment
#
using OITOOLS
using FITSIO

oifitsfile = "./data/MWC480.oifits"
data = vec(readoifits(oifitsfile, filter_bad_data = true, polychromatic = true)) # vec is to get rid of degenerate (temporal) dimension

nx = 64 #number of pixels (side)
pixsize = 0.2 # mas/pixel

fftplan = setup_nfft_polychromatic(data, nx, pixsize);
nwav = length(fftplan)

# Setup regularization
regularizers = [   [ ["centering", 1e4], ["tv", 7e2] ]]  # Frame 1 is centered
for i=1:nwav-1
    push!(regularizers,[["centering", 1e4], ["tv",7e2]]) # Total variation for all
end
push!(regularizers,[ ["transspectral_tvsq", 1e6] ] ); #transspectral regularization ties the frames together


#x_start = vec(repeat(ones(nx*nx)/(nx*nx), nwav));
#x_start = vec(rand(nx*nx, nwav));
pointsource = zeros(nx,nx); pointsource[div(nx+1,2), div(nx+1,2)] = 1.0;
x_start = zeros(nx, nx, nwav);
for i=1:nwav
    x_start[:,:,i]=pointsource
end
x_start = vec(x_start)
x = reconstruct_polychromatic(x_start, data, fftplan, regularizers = regularizers, maxiter = 500);

imdisp_polychromatic(reshape(x,nx*nx,nwav).^.25, pixscale=pixsize)
