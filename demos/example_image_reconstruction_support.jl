#
# Image reconstruction with soft support constrain 
# Example with Polaris MIRCX data 
#
using OITOOLS
using PyPlot
set_oiplot_defaults();
oifitsfile = "./data/polaris.oifits"
data = readoifits(oifitsfile)[1,1];

# Let's find an adequate mask by fitting the data
weights=[1.0,1.0,0.0]
disc = create_component(type="ud", name="disc");
disc.vis_params[1].val=2.0
disc.vis_params[1].minval=1.0
disc.vis_params[1].maxval=5.0
model=create_model(disc)
# Global minimization 
minf, minx, cvis_model, result = fit_model_ultranest(data, model, weights=weights); #interesting local minimum with v2 only!
# Local minimization 
minf, minx, cvis_model, result = fit_model_nlopt(data, model, weights=weights);

# Now imaging
pixsize = 0.05 # size of a pixel in milliarcseconds
nx = 128 # width of image (number of pixels)
weights=[1.0,1.0,1.0]
mask  = disk(npix=nx, diameter=minx[1]/pixsize+1) # binary mask
prior = model_to_image(model, nx=nx, pixsize=pixsize, oversample=1).*mask 
regularizers = [["centering", 1e2], ["l1l2", 2e8, 1e-4], ["support", 1.0, mask]];
ft = setup_nfft(data, nx, pixsize);
x = reconstruct(prior, data, ft, regularizers = regularizers, verb = true, maxiter=500, weights=weights);
x = reconstruct(x.*mask, data, ft, regularizers = regularizers, verb = true, maxiter=500, weights=weights);
chi2 = chi2_f(x, ft, data, weights=weights, verb=true); 
imdisp(x, pixsize=pixsize, beamsize=0.5*1/maximum(data.uv_baseline)*180/pi*3600*1000);
savefig("polaris_image.png")