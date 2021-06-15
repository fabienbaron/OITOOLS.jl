#
# Very Basic Image reconstruction code
#
#using OITOOLS
include("../src/OITOOLS.jl"); using Main.OITOOLS
oifitsfile = "./data/2004-data1.oifits"
pixsize = 0.202
nx = 64
data = readoifits(oifitsfile)[1,1];
data.t3phi_vonmises_err = gaussianwrapped_to_vonmises_fast.(data.t3phi_err/180*pi)
data.t3phi_vonmises_chi2_offset = log(2*pi) .+ 2*(logbesselI0.(data.t3phi_vonmises_err)-log.(data.t3phi_err/180*pi))
dft = setup_dft(data, nx, pixsize);

x = vec(readfits("./data/2004-64.fits"));
chi2_dft_f(x, dft, data)
chi2_dft_f(x, dft, data, vonmises=true)

# #initial image is a simple Gaussian
# x_start = gaussian2d(nx,nx,nx/6);
# x_start = vec(x_start)/sum(x_start);
# regularizers = [["centering", 1e4], ["l1l2", 7e6, 1e-3]];
# x = reconstruct(x_start, data, dft, regularizers = regularizers, verb = true, maxiter=500);
# imdisp(x,pixscale=pixsize)
# writefits(reshape(x,nx,nx),"reconstruction.fits")
