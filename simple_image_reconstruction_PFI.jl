#
# Very Basic Image reconstruction code
#
include("readoifits_vis.jl");
include("oichi2.jl");
include("oiplot.jl");

oifitsfile = "PFI_Y12_1200m_3m_10_11mu.10hrx10nights.HIGHSNR.oifits"
data = readoifits(oifitsfile)[1,1];
#Debug
# nuv=100000;
# indx=Int.(ceil.(data.nuv*rand(nuv)))
# data.nuv= nuv;
# data.nvisamp=nuv;
# data.nvisphi=nuv;
# data.visamp=data.visamp[indx];
# data.visamp_err=data.visamp_err[indx];
# data.visphi=data.visphi[indx];
# data.visphi_err=data.visphi_err[indx];
# data.uv=data.uv[:,indx];

pixsize = 0.7
nx = 801
#initial image is a simple Gaussian
x_start = Array{Float64}(nx, nx);
for i=1:nx
   for j=1:nx
     x_start[i,j] = exp(-((i-(nx+1)/2)^2+(j-(nx+1)/2)^2)/(2*(nx/6)^2));
   end
 end
x_start = vec(x_start)/sum(x_start);

using OptimPack

# NFFT reconstruction
fftplan = setup_nfft(data.uv, nx, pixsize);
crit = (x,g)->chi2_vis_nfft_fg(x, g, fftplan, data);
x_sol = OptimPack.vmlmb(crit, x_start, verb=true, lower=0, maxiter=20, blmvm=false, gtol=(1e-6,1e-6));
x_sol = OptimPack.vmlmb(crit, x_sol, verb=true, lower=0, maxiter=20, blmvm=false, gtol=(1e-6,1e-6));
imdisp(x_sol)
f = FITS("reconst.fits", "w"); write(f, reshape(x_sol,(nx,nx))); close(f);
