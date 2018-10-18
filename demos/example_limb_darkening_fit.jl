using OITOOLS
#
# EXAMPLE 2: fit uniform disc and limb-darkening law to data
#
oifitsfile = "./data/AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1]; # data can be split by wavelength, time, etc.
uvplot(data)
v2plot(data,logplot=true);
t3phiplot(data);

# Example of visibilities, here for Hestroffer with limb-darkening parameter 0.1
cvis = visibility_ldpow([8.0,0.1], data.v2_baseline);

# Fit uniform disc and plot
f_chi2, params, cvis_model = fit_model_v2(data, visibility_ud, [8.0]);# diameter is the parameter

# Plot model vs data
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data.v2_baseline,data.v2,data.v2_err, v2_model,logplot=true);

# Fit limb-darkened disc (Hestroffer) and plot
f_chi2, params, cvis_model = fit_model_v2(data, visibility_ldpow, [8.0,0.]);#diameter, ld1, ld2 coeffs
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data.v2_baseline,data.v2,data.v2_err, v2_model,logplot=true);

# Fit limb-darkened disc (quadratic) and plot
f_chi2, params, cvis_model = fit_model_v2(data, visibility_ldquad, [8.0,0.1,0.1]);#diameter, ld1, ld2 coeffs
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data.v2_baseline,data.v2,data.v2_err, v2_model,logplot=true);
