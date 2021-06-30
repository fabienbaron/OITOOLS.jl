using OITOOLS
#
# EXAMPLE 2: fit uniform disc and limb-darkening law to data
#

# Check https://arxiv.org/abs/1610.06185 for official results

# New interface -- coming soon


printstyled("WARNING --- DEPRECATED --- ", color=:red)

oifitsfile = "./data/AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1]; # data can be split by wavelength, time, etc.
uvplot(data)
uvplot(data, bywavelength=true)
v2plot(data,logplot=true);
#t3phiplot(data);


model=create_model()

minf, minx, cvis_model, result = fit_model_ultranest(data, model);
minf, minx, cvis_model, result = fit_model_levenberg(data, model);
minf, minx, cvis_model, result = fit_model_nlopt(data, model);



# Example of visibilities, here for Hestroffer with limb-darkening parameter 0.1
cvis = visibility_ldpow([8.0,0.1], data.uv);

# Fit uniform disc and plot
f_chi2, params, cvis_model = fit_model(data, visibility_ud, [8.0],weights=[1.0,0,0]);# diameter is the parameter, chi2 ~ 15.23

# Plot model vs data
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data, v2_model,logplot=true);

# Plot model vs function
v2plot_modelvsfunc(data, visibility_ud,params,logplot=true);

# Fit limb-darkened disc (Hestroffer) and plot
f_chi2, params, cvis_model = fit_model(data, visibility_ldpow, [8.0,0.], weights=[1.0,0,0]);#diameter, ld1, ld2 coeffs
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data, v2_model,logplot=true);

# Fit limb-darkened disc (linear) and plot
f_chi2, params, cvis_model = fit_model(data, visibility_ldlin, [8.0,0.1], weights=[1.0,0,0]);#diameter, ld1, ld2 coeffs
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data, v2_model,logplot=true);

# Fit limb-darkened disc (quadratic) and plot
f_chi2, params, cvis_model = fit_model(data, visibility_ldquad, [8.0,0.1,0.1], weights=[1.0,0,0]);#diameter, ld1, ld2 coeffs
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data, v2_model,logplot=true);

# Fit limb-darkened disc (square root) and plot
f_chi2, params, cvis_model = fit_model(data, visibility_ldsquareroot, [8.0,0.1,0.1], weights=[1.0,0,0]);#diameter, ld1, ld2 coeffs
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data, v2_model,logplot=true);

# Directly compare chi2 for given law
chi2v2 = model_to_chi2(data, visibility_ud, [8.306655883789062], weights=[1.0,0,0])
chi2v2 = model_to_chi2(data, visibility_ldquad, [8.410184002460234, 0.1917095059503015, -0.03404992047596201], weights=[1.0,0,0])

# Example of fitting with bound constraints
f_chi2, params, cvis_model = fit_model(data, visibility_ud, [8.0], lbounds=[7.5], hbounds=[8.2], weights=[1.0,0,0]);# will stop at upper bound
