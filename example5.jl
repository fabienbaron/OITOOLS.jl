#
# Very Basic Image reconstruction code
#

include("oitools.jl")
using OptimPack
nw = 1;# monochromatic mode
oifitsfile = "2004-data1.oifits"
pixsize = 0.1
nx = 128
data = (readoifits(oifitsfile))[1,1];
dft = setup_dft(data, nx, pixsize);
#nfft = setup_nfft(data, nx, pixsize);
#initial image is a simple Gaussian
 x_start = Array{Float64}(nx, nx);
     for i=1:nx
       for j=1:nx
         x_start[i,j] = exp(-((i-(nx+1)/2)^2+(j-(nx+1)/2)^2)/(2*(nx/6)^2));
       end
     end
 x_start = vec(x_start)/sum(x_start);
#x_start = rand(nx*nx);
#x_start /=sum(x_start);
mu = 1e3
crit = (x,g)->chi2_centered_fg(x, g, mu, dft, data);
x_sol = OptimPack.vmlmb(crit, x_start, verb=true, lower=0, maxiter=80, blmvm=false);
f = FITS("reconst.fits", "w"); write(f, reshape(x_sol,(nx,nx))); close(f);
