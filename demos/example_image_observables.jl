using OITOOLS
#
# EXAMPLE 1: read image and data, and compare the observables
#

# read the image file
fitsfile = "./data/2004true.fits";
pixsize = 0.101; # in mas/pixel
x_true = readfits(fitsfile); 
nx = (size(x_true))[1];

# display the image
imdisp(x_true, pixsize = pixsize, tickinterval = 1.0, beamsize = 1.0, beamlocation = [0.85, 0.85], use_colorbar = true)

# read the data file
oifitsfile = "./data/2004-data1.oifits";
data = readoifits(oifitsfile); # data can be split by wavelength, time, etc.

# Display the data
uvplot(data);
plot_v2(data,logplot=true);
plot_t3phi(data);
plot_t3amp(data);

# Compare data to image
#

# DFT method
data = data[1,1]
# Setup Fourier transform via DFT
dft = setup_dft(data, nx, pixsize);
# This computes the complete chi2
f_chi2 = chi2_f(x_true, dft, data);
# Compute |V|^2 observables and plot
@time cvis_model = image_to_vis(x_true, dft);
v2_model = vis_to_v2(cvis_model, data.indx_v2);
plot_v2_residuals(data, v2_model);

# NFFT method
ft = setup_nfft(data, nx, pixsize);
# This computes the complete chi2
f_chi2 = chi2_f(x_true, ft, data);
# Compute |V|^2 observables and plot
@time cvis_model = image_to_vis(x_true, ft);
v2_model = vis_to_v2(cvis_model, data.indx_v2);
plot_v2_residuals(data, v2_model);
