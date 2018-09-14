#
# Very Basic Image reconstruction code
#
include("oitools.jl");
using OptimPack;
#oifitsfile = "betlyr6t.oifits"
oifitsfile = "2004-data1.oifits"

pixsize = 0.1
nx = 137
#data = readoifits(oifitsfile, filter_bad_data = true)[1,1];
data = readoifits(oifitsfile)[1,1];
fftplan = setup_nfft(data, nx, pixsize);
#initial image is a simple Gaussian
 x_start = Array{Float64}(nx, nx);
     for i=1:nx
       for j=1:nx
         x_start[i,j] = exp(-((i-(nx+1)/2)^2+(j-(nx+1)/2)^2)/(2*(nx/8)^2));
       end
     end
 x_start = vec(x_start)/sum(x_start);

# TOTAL VARIATION
crit = (x,g)->chi2_centered_tv_nfft_fg(x, g, fftplan, data, mu = 1e4);
x_sol = OptimPack.vmlmb(crit, x_start, verb=true, lower=0, maxiter=200, blmvm=false, gtol=(1e-8,1e-8));
for i=1:4
    x_sol = OptimPack.vmlmb(crit, x_sol, verb=true, lower=0, maxiter=200, blmvm=false, gtol=(1e-8,1e-8));
end
imdisp(x_sol,pixscale=pixsize)
f = FITS("reconst137_TV.fits", "w"); write(f, reshape(x_sol,(nx,nx))); close(f);

#EPLL
include("../MLIR/EPLL.jl")
using JLD;
dict = load("GMM_YSO.jld", "GMM");
crit = (x,g)->chi2_centered_EPLL_nfft_fg(x, g, fftplan, data, dict, mu=1e1);
@time x_sol = OptimPack.vmlmb(crit, x_start, verb=true, lower=0, maxiter=200, blmvm=false, gtol=(1e-8,1e-8));
imdisp(x_sol,pixscale=pixsize)
f = FITS("reconst137_EPLL.fits", "w"); write(f, reshape(x_sol,(nx,nx))); close(f);
