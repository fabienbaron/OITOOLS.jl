# Comparing NFFT vs DFT on complex visibility reconstruction
include("oitools.jl")
oifitsfile = "2004-data1.oifits";
data = readoifits(oifitsfile)[1,1];
xtrue = vec(read(FITS("2004-64.fits")[1]));
xtrue /=sum(xtrue);

# Generate fake data
cvis_true = setup_dft(data.uv, 64, 0.1)*xtrue;
data.visamp_err = 0.001*ones(data.nuv)
data.visamp = abs.(cvis_true)+ data.visamp_err.*randn(data.nuv)
data.nvisamp = data.nuv
data.nvisphi = data.nuv
data.visphi_err= 0.5*ones(size(cvis_true));
data.visphi = angle.(cvis_true)*180./pi + data.visphi_err.*randn(size(cvis_true));

# Setup both DFT and NFFT for reconstruction
nx = 256
pixsize = 0.025

#initial image is a simple Gaussian
x_start = Array{Float64}(nx, nx);
     for i=1:nx
       for j=1:nx
         x_start[i,j] = exp(-((i-(nx+1)/2)^2+(j-(nx+1)/2)^2)/(2*(nx/6)^2));
       end
     end
 x_start = vec(x_start)/sum(x_start);

using OptimPack

#DFT reconstruction
dft = setup_dft(data.uv, nx, pixsize);
tic();
crit = (x,g)->chi2_vis_dft_fg(x, g, dft, data);
x_sol = OptimPack.vmlmb(crit, x_start, verb=true, lower=0, maxiter=80, blmvm=false);
toc();
imdisp(x_sol)

# NFFT reconstruction
fftplan = setup_nfft(data.uv, nx, pixsize);
tic();
crit = (x,g)->chi2_vis_nfft_fg(x, g, fftplan, data);
x_sol = OptimPack.vmlmb(crit, x_start, verb=true, lower=0, maxiter=80, blmvm=false);
toc();
imdisp(x_sol)
