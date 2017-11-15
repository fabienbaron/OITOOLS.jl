using FITSIO
#include("checkpackages.jl")
include("readoifits.jl")
include("oichi2.jl")
include("oiplot.jl")

#
# Main
#

#read model fits file
fitsfile = "2004true.fits";

pixsize = 0.1; # in mas/pixel
x_true = (read((FITS(fitsfile))[1])); nx = (size(x_true))[1]; x_true=vec(x_true);

oifitsfile = "2004-data1.oifits";
data = (readoifits(oifitsfile))[1,1]; # data can be split by wavelength, time, etc.

# setup DFT (nfft also possible)
dft = setup_dft(data, nx, pixsize);

# this computes the complete chi2
f_true = chi2(x_true, dft, data);


#uv plot
uvplot(data.uv)

#v2 plot
v2plot(data.v2_baseline,data.v2,data.v2_err);

#t3phi plot
t3phiplot(data.t3_baseline,data.t3phi,data.t3phi_err);

#v2 model vs data
# compute observables cvis then v2
cvis_model = image_to_cvis_dft(x_true, dft);
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data.v2_baseline,data.v2,data.v2_err, v2_model);

# the image
imdisp(x_true)
