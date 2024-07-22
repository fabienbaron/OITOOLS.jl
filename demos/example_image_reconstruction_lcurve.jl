#
# Image reconstruction using total variation and l-curve
#
using OITOOLS, PyPlot
oifitsfile = "./data/2004-data1.oifits"
pixsize = 0.2
nx = 64
data = readoifits(oifitsfile)[1,1];
fftplan = setup_nfft(data, nx, pixsize);
#initial image is a simple Gaussian
x_start = gaussian2d(nx,nx,nx/6);

# L-CURVE
# in this example we're looking for the best total variation weight value
#
tv_weights = [1e1, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5]
lcurve_chi2 = zeros(length(tv_weights))
lcurve_reg = zeros(length(tv_weights))
for i=1:length(tv_weights)
   regularizers = [["centering", 1e3], ["tv", tv_weights[i]]];
   x = reconstruct(x_start, data, fftplan, regularizers = regularizers, verb = false, maxiter=500);
   g = similar(x);
      for t=1:3 # uncomment to make sure we converged
          x = reconstruct(x, data, fftplan, regularizers = regularizers, verb = false, maxiter=100);
     end
   lcurve_chi2[i] = chi2_fg(x, g, fftplan, data);
   lcurve_reg[i] = regularization(x,g, regularizers=regularizers)/tv_weights[i];
   #imdisp(x,pixscale=pixsize)
end


# True L-shape
clf();
loglog(lcurve_chi2, lcurve_reg);
scatter(lcurve_chi2, lcurve_reg);
for i=1:length(tv_weights)
text(lcurve_chi2[i]*1.1, lcurve_reg[i], "μ=$(tv_weights[i])", size="small")
end
ylabel("Regularization");
xlabel("Chi2")


clf();
lcurve_reg .*= tv_weights
loglog(lcurve_reg, lcurve_chi2);
scatter(lcurve_reg, lcurve_chi2);
for i=1:length(tv_weights)
text(lcurve_reg[i], lcurve_chi2[i]*1.1, "μ=$(tv_weights[i])", size="xx-small")
end
xlabel("Regularization");
ylabel("Chi2")
