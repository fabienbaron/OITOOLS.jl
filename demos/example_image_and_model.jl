using OITOOLS
#
# EXAMPLE 1: read image and data, and compare the observables
#

# read the image file
fitsfile = "./data/2004true.fits";
pixsize = 0.101; # in mas/pixel
x_true = readfits(fitsfile); nx = (size(x_true))[1]; x_true=vec(x_true);

# display the image
imdisp(x_true, pixscale = pixsize, tickinterval = 1.0, beamsize = 1.0, beamlocation = [0.85, 0.85])

# read the data file
oifitsfile = "./data/2004-data1.oifits";
data = (readoifits(oifitsfile))[1,1]; # data can be split by wavelength, time, etc.

# Display the data
uvplot(data);
v2plot(data,logplot=true);# Alternatively, one can do v2plot(data.v2_baseline,data.v2,data.v2_err,logplot=true);
t3phiplot(data);


# Compare data to image
#


# DFT method

# Setup Fourier transform via DFT
dft = setup_dft(data, nx, pixsize);
# This computes the complete chi2
f_chi2 = chi2_dft_f(x_true, dft, data);
# Compute |V|^2 observables and plot
cvis_model = image_to_cvis_dft(x_true, dft);
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data.v2_baseline,data.v2,data.v2_err, v2_model);

# NFFT method
ft = setup_nfft(data, nx, pixsize);
# This computes the complete chi2
f_chi2 = chi2_nfft_f(x_true, ft, data);
# Compute |V|^2 observables and plot
cvis_model = image_to_cvis_nfft(x_true, ft);
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data.v2_baseline,data.v2,data.v2_err, v2_model);
