using LinearAlgebra, Statistics, OITOOLS

oifitsfile = "./data/2004-data1.oifits"
pixsize = 0.2# size of a pixel in milliarcseconds
nx = 64 # width of image (number of pixels)
data = readoifits(oifitsfile)[1,1];
dft = setup_dft(data, nx, pixsize);
ft = setup_nfft(data, nx, pixsize);
x = 10*vec(gaussian2d(nx,nx,nx/6));

# V2
dftg = (xx,gg)->chi2_dft_fg(xx, gg, dft, data, weights = [1,0,0], verb=false, vonmises=false);
nfftg= (xx,gg)->chi2_nfft_fg(xx, gg, ft, data, weights = [1,0,0], verb=false, vonmises=false);
g_dft=similar(x); dftg(x, g_dft); g_dft[:] = (g_dft .- sum(vec(x).*g_dft) / sum(x) ) / sum(x);
g_nfft=similar(x); nfftg(x, g_nfft); g_nfft[:] = (g_nfft .- sum(vec(x).*g_nfft) / sum(x) ) / sum(x);
print("V2 ", norm(g_dft-g_nfft), mean(g_dft./g_nfft))

# T3AMP
dftg = (xx,gg)->chi2_dft_fg(xx, gg, dft, data, weights = [0,1,0], verb=false, vonmises=false);
nfftg= (xx,gg)->chi2_nfft_fg(xx, gg, ft, data, weights = [0,1,0], verb=false, vonmises=false);
g_dft=similar(x); dftg(x, g_dft); g_dft[:] = (g_dft .- sum(vec(x).*g_dft) / sum(x) ) / sum(x);
g_nfft=similar(x); nfftg(x, g_nfft); g_nfft[:] = (g_nfft .- sum(vec(x).*g_nfft) / sum(x) ) / sum(x);
norm(g_dft-g_nfft)
mean(g_dft./g_nfft)

# T3PHI
dftg = (xx,gg)->chi2_dft_fg(xx, gg, dft, data, weights = [0,0,1], verb=false, vonmises=false);
nfftg= (xx,gg)->chi2_nfft_fg(xx, gg, ft, data, weights = [0,0,1], verb=false, vonmises=false);
g_dft=similar(x); dftg(x, g_dft); g_dft[:] = (g_dft .- sum(vec(x).*g_dft) / sum(x) ) / sum(x);
g_nfft=similar(x); nfftg(x, g_nfft); g_nfft[:] = (g_nfft .- sum(vec(x).*g_nfft) / sum(x) ) / sum(x);
norm(g_dft-g_nfft)
mean(g_dft./g_nfft)

# T3PHI - von Mises
data.t3phi_vonmises_err = gaussianwrapped_to_vonmises_fast.(data.t3phi_err/180*pi);
data.t3phi_vonmises_chi2_offset = log(2*pi) .+ 2*(logbesselI0.(data.t3phi_vonmises_err)-log.(data.t3phi_err/180*pi));
dftg = (xx,gg)->chi2_dft_fg(xx, gg, dft, data, weights = [0,0,1], verb=false, vonmises=true);
nfftg= (xx,gg)->chi2_nfft_fg(xx, gg, ft, data, weights = [0,0,1], verb=false, vonmises=true);
g_dft=similar(x); dftg(x, g_dft); g_dft[:] = (g_dft .- sum(vec(x).*g_dft) / sum(x) ) / sum(x);
g_nfft=similar(x); nfftg(x, g_nfft); g_nfft[:] = (g_nfft .- sum(vec(x).*g_nfft) / sum(x) ) / sum(x);
norm(g_dft-g_nfft)
mean(g_dft./g_nfft)

# TIMING
x0 = 10*vec(gaussian2d(nx,nx,nx/6));
@time x_dft = reconstruct(x0, data, dft, verb = false, maxiter=500);
@time x_nfft = reconstruct(x0, data, ft, verb = false, maxiter=500);
