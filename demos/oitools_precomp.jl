# # How to use:
# using PackageCompiler;  create_sysimage(:OITOOLS, sysimage_path="oitools.so", precompile_execution_file="oitools_precomp.jl"); exit()
# # Then you can launch Julia with
# julia  -t auto -Joitools.so

using OITOOLS
using AstroTime
oifitsfile = "./data/AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1];

# Model fitting with new flat-dict interface
param_ud = Dict{String,Any}("star,ud" => 8.0, "star,f" => 1.0)
result = fit_model(param_ud, ["star,ud"], data, weights=[1.0,0,0])

param_ld = Dict{String,Any}("star,ldpow" => 8.0, "star,alpha" => 0.15, "star,f" => 1.0)
result = fit_model(param_ld, ["star,ldpow", "star,alpha"], data, weights=[1.0,0,0])

#Simulate an observation using an analytic model of your target
#In this example, observations start at UT 2018-08-13 at 3:00:00AM and last until 2018-08-13 8:30:00AM, with a period of 15 minutes
dates = collect(from_utc(2018,8,13,3,0,0.0):15minutes:from_utc(2018,8,13,8,30,0.0))

# Model info - here a simple limb-darkened disk
param = Dict{String,Any}(
    "star,ldlin" => 3.0,
    "star,u"     => 0.15,
    "star,f"     => 1.0,
)
model = parse_model(param, String[])

# Target info
target = read_obs_file("default_obs"); # read defaults (for OIFITS header)
target.target = "EZ Psc"
target.raep0 =  [20, 57, 59.4437981]'*[1.0, 1/60., 1/3600] # set ra
target.decep0 = [46, 28, 00.5731825]'*[1.0, 1/60., 1/3600] # set dec
# Simulation info
facility    = read_facility_file("CHARA_new");
combiner    = read_comb_file("MIRC");
wavelength        = read_wave_file("MIRC_LOWH");

# Output file
out_file="./data/simulation-limb_darkened_disk.oifits";

simulate(facility, target, combiner, wavelength, dates, out_file, flat_model=model, flat_params=Float64[]);

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
image_to_chi2(x_start, ft, data, verb=true); # Starting chi2
regularizers = [["centering", 1e3], ["l1l2", 7e6, 1e-3]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = true, maxiter=5);
