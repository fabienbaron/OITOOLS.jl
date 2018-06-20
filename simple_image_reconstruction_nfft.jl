#
# Very Basic Image reconstruction code
#
include("oitools.jl");
using OptimPack;
oifitsfile = "2004-data1.oifits"
pixsize = 0.2
nx = 64
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
crit = (x,g)->chi2_centered_tv_nfft_fg(x, g, fftplan, data, mu = 7e3);
@time x_sol = OptimPack.vmlmb(crit, x_start, verb=true, lower=0, maxiter=2000, blmvm=false, gtol=(1e-8,1e-8));
imdisp(x_sol,pixscale=pixsize)
#f = FITS("reconst.fits", "w"); write(f, reshape(x_sol,(nx,nx))); close(f);
