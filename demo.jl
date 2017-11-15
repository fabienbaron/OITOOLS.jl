include("oitools.jl");

#
# EXAMPLE 1: read image and data, and compare the observables
#

# read the image file
fitsfile = "2004true.fits";
pixsize = 0.1; # in mas/pixel
x_true = (read((FITS(fitsfile))[1])); nx = (size(x_true))[1]; x_true=vec(x_true);

# display the image
imdisp(x_true)

# read the data file
oifitsfile = "2004-data1.oifits";
data = (readoifits(oifitsfile))[1,1]; # data can be split by wavelength, time, etc.

# Display the data
uvplot(data.uv);
v2plot(data.v2_baseline,data.v2,data.v2_err,logplot=true);
t3phiplot(data.t3_baseline,data.t3phi,data.t3phi_err);


# Compare data to image
#

# Setup Fourier transform via DFT (NFFT also possible)
dft = setup_dft(data, nx, pixsize);
# This computes the complete chi2
f_chi2 = chi2(x_true, dft, data);

# Compute |V|^2 observables and plot
cvis_model = image_to_cvis_dft(x_true, dft);
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data.v2_baseline,data.v2,data.v2_err, v2_model);


#
# EXAMPLE 2: fit uniform disc and limb-darkening law to data
#
oifitsfile = "AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1]; # data can be split by wavelength, time, etc.
uvplot(data.uv)
v2plot(data.v2_baseline,data.v2,data.v2_err,logplot=true);
t3phiplot(data.t3_baseline,data.t3phi,data.t3phi_err);

# Fit uniform disc and plot
f_chi2, params, cvis_model = fit_model_v2(data, visibility_ud, [8.0]);# diameter is the parameter
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data.v2_baseline,data.v2,data.v2_err, v2_model,logplot=true);

# Fit limb-darkened disc (quadratic) and plot
f_chi2, params, cvis_model = fit_model_v2(data, visibility_ldquad, [8.0,0.1,0.1]);#diameter, ld1, ld2 coeffs
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data.v2_baseline,data.v2,data.v2_err, v2_model,logplot=true);
