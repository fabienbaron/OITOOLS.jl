#
# Image reconstruction code with spectral regularization
# This requires a private example file at the moment
#
using OITOOLS
oifitsfile = "./data/MWC480.oifits"

# Step 1: reconstruct a "gray" image
data = readoifits(oifitsfile, filter_bad_data = true, polychromatic = false)[1,1]
nx = 64 #number of pixels (side)
pixsize = 0.1 # mas/pixel
ft = setup_nfft(data, nx, pixsize);
regularizers = [["centering", 1e7]] 
pointsource = zeros(nx,nx); 
pointsource[div(nx+1,2), div(nx+1,2)] = 1.0;
x = copy(pointsource);
x = reconstruct(x, data, ft, regularizers = regularizers, maxiter = 500, verb=true);
imdisp(x.^.2, pixsize=pixsize)

# Step 2: polychromatic reconstruction
data = vec(readoifits(oifitsfile, filter_bad_data = true, polychromatic = true)) # vec is to get rid of degenerate (temporal) dimension
nx = 64 #number of pixels (side)
pixsize = 0.2 # mas/pixel
ft = setup_nfft_polychromatic(data, nx, pixsize);
nwavs = length(ft)

# Setup regularization
regularizers = [   [ ["centering", 1e4], ["tv", 1e2] ]]  # Frame 1 is centered
for i=1:nwavs-1
    push!(regularizers,[["tv",1e2]]) # Total variation for all
end

regularizers = [   [ ["centering", 1e4], ["l1l2", 1e5, 0.4] ]]  # Frame 1 is centered
for i=1:nwavs-1
    push!(regularizers,[["l1l2",1e5,0.4]]) # Total variation for all
end

# Uncomment the desired transspectral regularization
# push!(regularizers,[ ["transspectral_tvsq", 1e5] ] );
push!(regularizers,[["transspectral_structnorm", 1e2], ["transspectral_tv", 1e2] ] );

pointsource = zeros(nx,nx); 
pointsource[div(nx+1,2), div(nx+1,2)] = 1.0;
x_start = repeat(pointsource, 1,1,nwavs);

x = copy(x_start);
for i=1:10
    x = reconstruct_polychromatic(x, data, ft, regularizers = regularizers, maxiter = 200, verb=false);
end
imdisp_polychromatic(x.^.2, pixsize=pixsize)

# Check chi2
chi2 = chi2_polychromatic_f(x, ft, data, verb = true)