include("../src/OITOOLS.jl"); using Main.OITOOLS

using AstroTime
#Simulate an observation using an input image, given telescope parameters, and input image or image filename and observation times
#In this example, observations start at UT 2018-08-13 at 3:00:00AM and last until 2018-08-13 8:30:00AM, with a period of 15 minutes
dates = collect(TAIEpoch(2018,8,13,3,0,0.0):15minutes:TAIEpoch(2018,8,13,8,30,0.0))

# Input image defining the target brightness distribution
image_file="./data/2004true.fits";
pixsize=0.101;

# Target info
target = read_obs_file("./data/default_obs.txt"); # read defaults (for OIFITS header)
target.target[1] = "AZ Cyg"
target.raep0[1] =  [20, 57, 59.4437981]'*[1.0, 1/60., 1/3600] # set ra
target.decep0[1] = [46, 28, 00.5731825]'*[1.0, 1/60., 1/3600] # set dec


# Simulation info
facility    = read_facility_file("./data/CHARA_new.txt");
combiner    = read_comb_file("./data/MIRC.txt");
wavelength        = read_wave_file("./data/MIRC_LOWH.txt");
v2m=1.0/100; v2a=1e-5; t3ampm=1.0/100; t3ampa=1e-6; t3phim=0.0; t3phia=0.5;
errors      = define_errors(v2m,v2a,t3ampm,t3ampa,t3phim,t3phia);

errors      = define_errors(0.0,1e-3,0.0,1e-6,0.0,0.1);
# Output file
out_file="./data/simulation-2004image.oifits";

simulate(facility, target, combiner, wavelength, dates, errors, out_file, image=image_file, pixsize=pixsize);

#Check simulated data
data = (readoifits(out_file, filter_bad_data=false))[1,1]; # data can be split by wavelength, time, etc.
#uvplot(data)
uvplot(data, color="wav")
#uvplot(data, color="mjd")

#image_file="./data/2004true.fits";
x_true = readfits(image_file)
x_true /= sum(x_true)
ft = setup_nfft(data, size(x_true,1), pixsize);
f_chi2 = chi2_f(x_true, ft, data, verb=true);
image_to_vis(x, ft)