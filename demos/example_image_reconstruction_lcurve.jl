#
# Image reconstruction using total variation and l-curve
#
using OITOOLS
oifitsfile = "./data/2004-data1.oifits"
pixsize = 0.1
nx = 137
data = readoifits(oifitsfile)[1,1];
fftplan = setup_nfft(data, nx, pixsize);
#initial image is a simple Gaussian
x_start = gaussian2d(nx,nx,nx/6);
x_start = vec(x_start)/sum(x_start);


# L-CURVE
# in this example we're looking for the best total variation weight value
#
tv_weights = [1e1, 1e2, 1e3, 2e3, 5e3, 1e4]
lcurve_chi2 = zeros(length(tv_weights))
lcurve_reg = zeros(length(tv_weights))
for i=1:length(tv_weights)
   regularizers = [["centering", 1e4], ["tv", tv_weights[i]]];
   x = reconstruct(x_start, data, fftplan, regularizers = regularizers, verb = false);
   g = Array{Float64}(size(x));
   for t=1:3 # make sure we converged
       x = reconstruct(x, data, fftplan, regularizers = regularizers, verb = false);
   end
   lcurve_chi2[i] = chi2_nfft_fg(x, g, fftplan, data);
   lcurve_reg[i] = regularization(x,g, regularizers=regularizers);
   imdisp(x,pixscale=pixsize)
end
clf(); loglog(lcurve_reg, lcurve_chi2); scatter(lcurve_reg, lcurve_chi2); xlabel("Regularization"); ylabel("Chi2")
