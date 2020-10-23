#
# Very Basic Image reconstruction code
#
using OITOOLS
oifitsfile = "./data/2004-data1.oifits"
pixsize = 0.2
nx = 64
data = readoifits(oifitsfile)[1,1];
dft = setup_dft(data, nx, pixsize);
#initial image is a simple Gaussian
x_start = gaussian2d(nx,nx,nx/6);
x_start = vec(x_start)/sum(x_start);
regularizers = [["centering", 1e4], ["l1l2", 7e6, 1e-3]];
x = reconstruct(x_start, data, dft, regularizers = regularizers, verb = true, maxiter=500);
imdisp(x,pixscale=pixsize)
writefits(reshape(x,nx,nx),"reconstruction.fits")
