#
# Image reconstruction of the BC2026 data
#
using OITOOLS
oifitsfile = "./data/BC2026/OBJECT1_N.oifits"

# In this example we will not try a polychromatic reconstruction
data = readoifits(oifitsfile)[1,1];

# Fourier transform setup
pixsize = 1.0 # size of a pixel in milliarcseconds
nx = 128 # width of image (number of pixels)
ft = setup_nfft(data, nx, pixsize);
#initial image is a simple Gaussian
x_start = gaussian2d(nx,nx,nx/6);
image_to_chi2(x_start, ft, data, verb=true); # Evaluate chi2
regularizers = [["centering", 1e4], ["l1l2", 1e9, 1e-3]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = true, maxiter=500);
imdisp(x,pixsize=pixsize)
writefits(x, "object1N_gray.fits")