using FITSIO
include("readoifits.jl")
include("setupft.jl")
include("oichi2.jl")
include("oiplot.jl")

#
# Main
#

#read model fits file
pixellation = 0.1; # in mas
fitsfile = "azcyg2014mean_mean.fits";
#fitsfile = "2004true.fits";
scale_rad = pixellation * (pi / 180.0) / 3600000.0;
x_true = (read((FITS(fitsfile))[1])); nx = (size(x_true))[1]; x_true=vec(x_true);

oifitsfile = "azcyg2014.oifits";
#oifitsfile="2004-data1.oifits";
data = read_oifits(oifitsfile);

# setup DFT (nfft also possible soon)
dft = setup_ft(data, nx);

# this computes the complete chi2
f_true = chi2(x_true, dft, data);


# now plot some stuff
PyPlot.show()
#uv plot
uvplot(data.uv[1,:],data.uv[2,:])

#t3phi plot
t3phiplot(data.baseline_t3,data.t3phi_data,data.t3phi_data_err);

#v2 model vs data
# compute observables cvis then v2
cvis_model = image_to_cvis(x_true, dft);
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data.baseline_v2,data.v2_data,data.v2_data_err, v2_model);

# the image
imdisp(x_true)
readline()
