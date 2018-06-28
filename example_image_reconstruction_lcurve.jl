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


# L-CURVE
# in this example we're looking for the best total variation weight value
#
tv_weights = [1e2, 1e3, 2e3, 3e3, 5e3, 6e3, 7e3, 8e3, 1e4, 1e5]
lcurve_chi2 = zeros(length(tv_weights))
lcurve_reg = zeros(length(tv_weights))
for i=1:length(tv_weights)
   regularizers = [["centering", 1e4], ["tv", tv_weights[i]]];
   x = reconstruct(x_start, fftplan, regularizers = regularizers, verb = false);
   g = Array{Float64}(size(x));
   for t=1:10 # make sure we converged
       x = reconstruct(x, fftplan, regularizers = regularizers, verb = false);
   end
   lcurve_chi2[i] = chi2_nfft_fg(x, g, fftplan, data);
   lcurve_reg[i] = regularization(x,g, regularizers=regularizers);
   imdisp(x,pixscale=pixsize)
end
clf(); loglog(lcurve_reg, lcurve_chi2); scatter(lcurve_reg, lcurve_chi2); xlabel("Regularization"); ylabel("Chi2")
