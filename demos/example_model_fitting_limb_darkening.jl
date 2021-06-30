using OITOOLS
#
# EXAMPLE 2: fit uniform disc and limb-darkening law to data
#

# Check https://arxiv.org/abs/1610.06185 for official results

oifitsfile = "./data/AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1]; # data can be split by wavelength, time, etc.
# If you want to check the data:
# uvplot(data)
# uvplot(data, bywavelength=true)
# v2plot(data,logplot=true);
#t3phiplot(data);

#
# UNIFORM DISC FITTING
#
model = create_model(create_component(type="ud", name="Model"));

# You can check which parameters are free just by displaying the models

display(model);

# We can choose three main packages for optimization
# NLopt: several local and global optimizers
minf, minx, cvis_model, result = fit_model_nlopt(data, model, chi2_weights=[1.0,0,0]);

# LsqFit: Levenberg-Marquardt (fast but inaccurate statistical uncertainties)
#minf, minx, cvis_model, result = fit_model_levenberg(data, model, chi2_weights=[1.0,0,0]);

# UltraNest: Nested Sampling (best for statistical uncertainties)
minf, minx, cvis_model, result = fit_model_ultranest(data, model, chi2_weights=[1.0,0,0]);

# Note: the result structures are different and depend on the packages

# How to compute complex visibilities
# Here for Hestroffer power law with diameter = 8.0 and limb-darkening parameter 0.1
model = create_model(create_component(type="ldlin", name="Model"));
dispatch_params([8.0,0.1], model);
cvis = model_to_cvis(model, data)

# How to compute chi2 for given model (note: you don't need to compute cvis first)
chi2v2 = model_to_chi2(data, model, [8.3,0.2], chi2_weights=[1.0,0,0])

# Plot model vs data
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data, v2_model,logplot=true);


#chi2v2 = model_to_chi2(data, visibility_ud, [8.306655883789062], weights=[1.0,0,0])
#chi2v2 = model_to_chi2(data, visibility_ldquad, [8.410184002460234, 0.1917095059503015, -0.03404992047596201], weights=[1.0,0,0])

# # Example of fitting with bound constraints
# f_chi2, params, cvis_model = fit_model(data, visibility_ud, [8.0], lbounds=[7.5], hbounds=[8.2], weights=[1.0,0,0]);# will stop at upper bound
