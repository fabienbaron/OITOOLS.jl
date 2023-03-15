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
model = create_model(create_component(type="ldpow", name="Model 1"));

# You can check which parameters are free just by displaying the models

display(model);

# We can choose between three main packages for optimization
# NLopt: several local and global optimizers
# Default is Nelder Mead
minf, minx, cvis_model, result = fit_model_nlopt(model, data, weights=[1.0,0,0]);
# One could use ISRES as a global optimization strategy
minf, minx, cvis_model, result = fit_model_nlopt(model, data, fitter=:GN_ISRES, weights=[1.0,0,0]);

# LsqFit: Levenberg-Marquardt (fast but inaccurate statistical uncertainties)
minf, minx, cvis_model, result = fit_model_levenberg(model, data, weights=[1.0,0,0]);

# UltraNest: Nested Sampling (best for statistical uncertainties)
minf, minx, cvis_model, result = fit_model_ultranest(model, data, weights=[1.0,0,0]);

# Note: the result structures are different and depend on the packages

# Let's look at the ultranest results




# How to compute complex visibilities
# Here for Hestroffer power law with diameter = 8.0 and limb-darkening parameter 0.1
model = create_model(create_component(type="ldlin", name="Model 2"));
dispatch_params([8.0,0.1], model);

# How to compute chi2 for given model (note: you don't need to compute cvis first)
chi2v2 = model_to_chi2(model, data, [8.3,0.2], weights=[1.0,0,0], verb=true)

v2_model, t3amp_model, t3phi_model = model_to_obs(model, data)
# Could have also called v2_model = vis_to_v2(cvis_model, data.indx_v2);

# Evaluate model at data points
cvis_model = model_to_cvis(model, data);

# Evaluate at any point
model_to_cvis(model, uv, λ)


# Plot model vs data
v2plot_model_vs_data(data, v2_model,logplot=true);

# Visualize the model as an image
img = model_to_image(model)
imdisp(img)

# TODO: evaluate model at any (u,v) point?
