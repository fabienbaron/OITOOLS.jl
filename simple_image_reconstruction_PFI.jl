#
# Very Basic Image reconstruction code
#
using FITSIO
using NFFT
FFTW.set_num_threads(8);
using OptimPack
include("readoifits_vis.jl")
include("oichi2.jl")
include("oiplot.jl")
oifitsfile = "PFI_Y12_1200m_3m_10_11mu.10hrx10nights.HIGHSNR.oifits"
pixsize = 4.0
nx = 128
data = readoifits(oifitsfile)[1,1];
#initial image is a simple Gaussian
x_start = Array{Float64}(nx, nx);
 for i=1:nx
       for j=1:nx
         x_start[i,j] = exp(-((i-(nx+1)/2)^2+(j-(nx+1)/2)^2)/(2*(nx/6)^2));
       end
     end
 x_start = vec(x_start)/sum(x_start);

fftplan = setup_nfft(data.uv, nx, pixsize);
#tic();
crit = (x,g)->chi2_vis_nfft_fg(x, g, fftplan, data);
x_sol = OptimPack.vmlmb(crit, x_start, verb=true, lower=0, maxiter=80, blmvm=false);
#toc();

#f = FITS("reconst.fits", "w"); write(f, reshape(x_sol,(nx,nx))); close(f);
