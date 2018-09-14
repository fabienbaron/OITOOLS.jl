#
# Very Basic Image reconstruction code
#
include("oitools.jl");
include("../MLIR/EPLL.jl")
using OptimPack;
oifitsfile = "2004-data1.oifits"
pixsize = 0.1
nx = 137
data = readoifits(oifitsfile)[1,1];
fftplan = setup_nfft(data, nx, pixsize);
#initial image is a simple Gaussian

x_start = Array{Float64}(nx, nx);
     for i=1:nx
       for j=1:nx
         x_start[i,j] = exp(-((i-(nx+1)/2)^2+(j-(nx+1)/2)^2)/(2*(nx/6)^2));
       end
     end
 x_start = vec(x_start)/sum(x_start);


regularizers = [["centering", 1e4], ["tv", 1e5]];
x = reconstruct(x_start, data, fftplan, regularizers = regularizers); # good start, no central star



regularizers = [["centering", 1e4], ["EPLL", 1.0, load("GMM_8x8_100_5000.jld", "GMM")]];
x = reconstruct(x, data, fftplan, regularizers = regularizers, maxiter = 10);
writefits(reshape(x,nx,nx),"reconstruction.fits")
imdisp(x,pixscale=pixsize)
