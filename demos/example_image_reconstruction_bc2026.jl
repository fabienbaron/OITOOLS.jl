#
# Very Basic Image reconstruction code
#
using OITOOLS
oifitsfile = "./data/BC2026/OBJECT1_N.oifits"
data = readoifits(oifitsfile)[1,1];
# Fourier transform setup
pixsize = 1.0 # size of a pixel in milliarcseconds
nx = 128 # width of image (number of pixels)
ft = setup_nfft(data, nx, pixsize);
# Once could also use the DFT instead of NFFT, but NFFT is much faster for larger images 
#dft = setup_dft(data, nx, pixsize);
#initial image is a simple Gaussian
x_start = gaussian2d(nx,nx,nx/6);
chi2_f(x_start, ft, data, verb=true); # Evaluate chi2
regularizers = [["centering", 1e4], ["l1l2", 7e8, 1e-3]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = true, maxiter=500);
imdisp(x,pixsize=pixsize)
# Uncomment if you want to write the result
#writefits(x,"reconstruction.fits")
