#
# Very Basic Image reconstruction code
#
include("oitools.jl");
oifitsfile = "2004-data1.oifits"
pixsize = 0.2
nx = 64
data = readoifits(oifitsfile)[1,1];
ft = setup_nfft(data, nx, pixsize);
#initial image is a simple Gaussian

x_start = Array{Float64}(nx, nx);
     for i=1:nx
       for j=1:nx
         x_start[i,j] = exp(-((i-(nx+1)/2)^2+(j-(nx+1)/2)^2)/(2*(nx/6)^2));
       end
     end
 x_start = vec(x_start)/sum(x_start);

regularizers = [["centering", 1e4], ["tv", 7e3]];
x = reconstruct(x_start, data, ft, regularizers = regularizers);
imdisp(x,pixscale=pixsize)
writefits(reshape(x,nx,nx),"reconstruction.fits")
