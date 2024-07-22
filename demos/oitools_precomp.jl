# # How to use:
# using PackageCompiler;  create_sysimage(:OITOOLS, sysimage_path="oitools.so", precompile_execution_file="oitools_precomp.jl"); exit()
# # Then you can launch Julia with
# julia  -t auto -Joitools.so

using OITOOLS
using AstroTime
oifitsfile = "./data/AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1]; 

model = create_model(create_component(type="ud", name="Model"));
minf, minx, cvis_model, result = fit_model_nlopt(data, model, weights=[1.0,0,0]);
minf, minx, cvis_model, result = fit_model_levenberg(data, model, weights=[1.0,0,0]);
minf, minx, cvis_model, result = fit_model_ultranest(data, model, weights=[1.0,1.0,1.0]);
#Simulate an observation using an analytic model of your target
#In this example, observations start at UT 2018-08-13 at 3:00:00AM and last until 2018-08-13 8:30:00AM, with a period of 15 minutes
dates = collect(from_utc(2018,8,13,3,0,0.0):15minutes:from_utc(2018,8,13,8,30,0.0))

# Model info - here a simple limb-darkened disk
model = create_model(create_component(type="ldlin", name="Model"));
dispatch_params([3.0,0.15], model);

# Target info
target = read_obs_file("./data/default_obs.txt"); # read defaults (for OIFITS header)
target.target[1] = "EZ Psc"
target.raep0[1] =  [20, 57, 59.4437981]'*[1.0, 1/60., 1/3600] # set ra
target.decep0[1] = [46, 28, 00.5731825]'*[1.0, 1/60., 1/3600] # set dec
# Simulation info
facility    = read_facility_file("./data/CHARA_new.txt");
combiner    = read_comb_file("./data/MIRC.txt");
wavelength        = read_wave_file("./data/MIRC_LOWH.txt");
v2m=1.0/100; v2a=1e-5; t3ampm=1.0/100; t3ampa=1e-6; t3phim=0.0; t3phia=0.5;
errors      = define_errors(v2m,v2a,t3ampm,t3ampa,t3phim,t3phia);

# Output file
out_file="./data/simulation-limb_darkened_disk.oifits";

simulate(facility, target, combiner, wavelength, dates, errors, out_file, model=model);

#Check simulated data
data = (readoifits(out_file))[1,1]; # data can be split by wavelength, time, etc.
uvplot(data)
uvplot(data, color="wav")
uvplot(data, color="mjd")
v2plot(data)
t3phiplot(data)
# Image reconstruction

oifitsfile = "./data/2004-data1.oifits"
pixsize = 0.2 # size of a pixel in milliarcseconds
nx = 64 # width of image (number of pixels)
data = readoifits(oifitsfile)[1,1];
ft = setup_nfft(data, nx, pixsize);
#initial image is a simple Gaussian
x_start = gaussian2d(nx,nx,nx/6);
x_start = vec(x_start)/sum(x_start);
chi2_f(x_start, ft, data, verb=true); # Starting chi2
regularizers = [["centering", 1e3], ["l1l2", 7e6, 1e-3]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = true, maxiter=5);
