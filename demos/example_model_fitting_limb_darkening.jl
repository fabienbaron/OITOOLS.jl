using OITOOLS
#
# EXAMPLE 2: fit uniform disc and limb-darkening law to data
#

# Check https://arxiv.org/abs/1610.06185 for official results

oifitsfile = "./data/AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1]; # data can be split by wavelength, time, etc.
#
# LIMB DARKENED DISC FITTING
#
model = create_model(create_component(type="ldpow", name="Model"));

# You can check which parameters are free just by displaying the models

display(model);

# We can choose between three main packages for optimization
# NLopt: several local and global optimizers
# Default is Nelder Mead
minf, minx, cvis_model, result = fit_model_nlopt(data, model, weights=[1.0,0,0]);
# One could use ISRES as a global optimization strategy
minf, minx, cvis_model, result = fit_model_nlopt(data, model, fitter=:GN_ISRES, weights=[1.0,0,0]);

# LsqFit: Levenberg-Marquardt (fast but inaccurate statistical uncertainties)
minf, minx, cvis_model, result = fit_model_levenberg(data, model, weights=[1.0,0,0]);

# UltraNest: Nested Sampling (best for statistical uncertainties)
minf, minx, cvis_model, result = fit_model_ultranest(data, model, weights=[1.0,0,0]);

# Note: the result structures are different and depend on the packages

# How to compute complex visibilities
# Here for Hestroffer power law with diameter = 8.0 and limb-darkening parameter 0.1
model = create_model(create_component(type="ldlin", name="Model"));
dispatch_params([8.0,0.1], model);
cvis_model = model_to_vis(model, data)


# How to compute chi2 for given model (note: you don't need to compute cvis first)
chi2v2 = model_to_chi2(data, model, weights=[1.0,0,0], verb=true)
chi2v2 = model_to_chi2(data, model, [8.3,0.2], weights=[1.0,0,0], verb=true)

# Plot model vs data
v2_model, t3amp_model, t3phi_model = model_to_obs(model,data);
plot_v2_residuals(data, v2_model,logplot=true);

# Visualize the model as an image
img = model_to_image(model, pixsize=0.1)
imdisp(img, pixsize=0.1)

diam = 3.4
ldd = 0.3

model = create_model(create_component(type="ldpow", name="Model"));
dispatch_params([diam,ldd], model);
disk_fft = model_to_image(model, pixsize=0.1, nx=64, oversample=1)
disk_img = limbdarkened_disk(diam/2, [3.0, ldd], nx=64, pixsize=0.1)




