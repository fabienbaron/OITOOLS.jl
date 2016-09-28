using FITSIO
using PyPlot,PyCall
PyPlot.show()
#@pyimport mpl_toolkits.axes_grid1 as axgrid

include("readoifits.jl")
include("setupft.jl")
include("oichi2.jl")
include("oiplot.jl")
##########################################
##########################################
#
# Code actually starts
#
pixellation = 0.1; # in mas
fitsfile = "2004true.fits";
oifitsfile = "2004-data1.oifits";
nw = 1;# monochromatic mode
#read model fits file
scale_rad = pixellation * (pi / 180.0) / 3600000.0;
x_true = (read((FITS(fitsfile))[1])); nx = (size(x_true))[1]; x_true=vec(x_true);

data = read_oifits(oifitsfile);
dft = setup_ft(data, nx);
#init required because of OptimPack way
f_true = chi2(x_true, dft, data);

# now plot some stuff

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
