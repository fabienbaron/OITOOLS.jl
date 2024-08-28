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
pointsource[nx÷2+1, nx÷2+1] = 1.0;
x = copy(pointsource);
for i=1:30
    x = reconstruct(x, data, ft, regularizers = regularizers, maxiter = 1000, verb=false);
end
imdisp(x.^.2, pixsize=pixsize)
plot_v2_residuals(x, data, ft)

x_mono = copy(x)
# Step 2: polychromatic reconstruction
data = vec(readoifits(oifitsfile, filter_bad_data = true, polychromatic = true)) # vec is to get rid of degenerate (temporal) dimension
nx = 64 #number of pixels (side)
pixsize = 0.1 # mas/pixel
ft = setup_nfft_polychromatic(data, nx, pixsize);
nwavs = length(ft)

# # Setup regularization
#  regularizers = [   [ ["centering", 1e4], ["tv", 1e2] ]]  # Frame 1 is centered
#  for i=1:nwavs-1
#      push!(regularizers,[["tv",1e2]]) # Total variation for all
#  end

regularizers = [[ ["centering", 1e7]]]  # Frame 1 is centered

for i=1:nwavs-1
    push!(regularizers,[]) # Total variation for all
end


# for i=1:nwavs-1
#      push!(regularizers,[["l1l2",1e5,0.4]]) # Total variation for all
#  end

# Uncomment the desired transspectral regularization
push!(regularizers,[ ["transspectral_structnorm", 1e2], ["transspectral_tvsq", 1e3] ] );
#push!(regularizers,[["transspectral_structnorm", 1e7], ["transspectral_tv", 1e4] ] );

pointsource = zeros(nx,nx); 
pointsource[nx÷2+1, nx÷2+1] = 1.0;
x_start = repeat(pointsource, 1,1,nwavs);
x = copy(x_start);

x = repeat(x_mono,1,1,nwavs)

for i=1:3
    x = reconstruct_polychromatic(x, data, ft, regularizers = regularizers, maxiter = 400, verb=true);
end
imdisp_polychromatic(x.^.2, pixsize=pixsize)

# Check chi2
chi2 = chi2_polychromatic_f(x, ft, data, verb = true)