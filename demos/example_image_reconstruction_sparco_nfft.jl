#
# Very Basic Image reconstruction code
#
include("../src/OITOOLS.jl");using Main.OITOOLS
oifitsfile = "./data/MWC480.oifits"
pixsize = 0.1
nx = 64
data = readoifits(oifitsfile, filter_bad_data=true)[1,1];
ft = setup_nfft(data, nx, pixsize);
params=[4.358e-01, 2.073e-01, 1.6e-6, 0.0] #V2:12.36 T3A: 4.14 T3P: 1.61
x=vec(readfits("../../squeeze/mwc480_grey.fits"));
chi2_sparco_nfft_f(x, ft, data, params)
