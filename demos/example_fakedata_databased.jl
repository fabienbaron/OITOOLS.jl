using OITOOLS
#
#EXAMPLE 9
#Simulate an observation using an input image and input oifits file

image_file="./data/2004true.fits"
pixsize=0.101
oifitsin="./data/2004-data1.oifits"
out_file="!./data/2004testsimulation.oifits"

simulate_obs(oifitsin,out_file,image_file,pixsize,dft=true,nfft=false)

#Compare simulated data to impate

#display the image
fitsfile = "./data/2004true.fits";
pixsize = 0.101; # in mas/pixel
x_true = readfits(fitsfile); nx = (size(x_true))[1]; x_true=vec(x_true);
imdisp(x_true, pixscale = pixsize, tickinterval = 1.0, beamsize = 1.0, beamlocation = [0.85, 0.85])

#read the data file

oifitsfile = "./data/2004testsimulation.oifits";
data = (readoifits(oifitsfile))[1,1]; # data can be split by wavelength, time, etc.
uvplot(data);
v2plot(data,logplot=true);# Alternatively, one can do v2plot(data.v2_baseline,data.v2,data.v2_err,logplot=true);
t3phiplot(data);

# NFFT method
ft = setup_nfft(data, nx, pixsize);
# This computes the complete chi2
f_chi2 = chi2_nfft_f(x_true, ft, data);
# Compute |V|^2 observables and plot
cvis_model = image_to_cvis_nfft(x_true, ft);
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data.v2_baseline,data.v2,data.v2_err, v2_model);
