#
# Pi Gruis, monochromatic reconstruction
#
using OITOOLS
oifitsfile = "./data/pigru.oifits"
pixsize = .5 # size of a pixel in milliarcseconds
nx = 64 # width of image (number of pixels)
data = readoifits(oifitsfile)[1,1];
# Fourier transform setup
ft = setup_nfft(data, nx, pixsize);
# Once could also use the DFT instead of NFFT, but NFFT is much faster for larger images 
#dft = setup_dft(data, nx, pixsize);

#initial image from model fitting
weights=[1.0,0.0,0.0]
disc = create_component(type="ldquad", name="disc");
disc.vis_params[1].val=10.0
disc.vis_params[1].minval=5.0
disc.vis_params[1].maxval=30.0
model=create_model(disc)
minf, minx, cvis_model, result = fit_model_ultranest(data, model, weights=weights); #interesting local minimum with v2 only!
x_start = model_to_image(model, nx=nx, pixsize=pixsize)

#x_start = vec(x_start)/sum(x_start)
chi2_f(x_start, ft, data, verb=true); # Evaluate chi2
regularizers = [["centering", 1e4], ["l1l2", 7e7, 1e-3]];
weights = [1.0, 0.0, 1.0]; # disable T3amp for VLTI/PIONIER
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = true, maxiter=500, weights=weights);
imdisp(x,pixsize=pixsize)
# Uncomment if you want to write the result
#writefits(x,"reconstruction.fits")
