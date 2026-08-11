using OITOOLS
using Dates
#Simulate an observation using an input image, given telescope parameters, and input image or image filename and observation times
#In this example, observations start at UT 2018-08-13 at 3:00:00AM and last until 2018-08-13 8:30:00AM, with a period of 15 minutes
dates = collect(DateTime(2018,8,13,3,0,0):Minute(15):DateTime(2018,8,13,8,30,0))

# Input image defining the target brightness distribution
image_file="./data/2004true.fits";
pixsize=0.101;

# Target info
target = read_obs_file("default_obs"); # read defaults (for OIFITS header)
target.target = "AZ Cyg"
target.raep0 =  15*[20, 57, 59.4437981]'*[1.0, 1/60., 1/3600] # set ra  (hms -> degrees)
target.decep0 =   [46, 28, 00.5731825]'*[1.0, 1/60., 1/3600] # set dec (dms -> degrees)


# Simulation info
facility    = read_facility_file("CHARA");
combiner    = read_comb_file("MIRCX");
wavelength        = read_wave_file("MIRCX_LOWH");

# Output file
out_file="./data/simulation-2004image.oifits";

simulate(facility, target, combiner, wavelength, dates, out_file, image=image_file, pixsize=pixsize);

#Check simulated data
data = (readoifits(out_file, filter_bad_data=false, T=Float64))[1,1]; # data can be split by wavelength, time, etc.
#uvplot(data)
uvplot(data, color="wav")
#uvplot(data, color="mjd")

#image_file="./data/2004true.fits";
x_true = readfits(image_file)
x_true /= sum(x_true)
ft = setup_nfft(data, size(x_true,1), pixsize);
f_chi2 = image_to_chi2(x_true, ft, data, verb=true);
# With noise off the chi2 above is 0 to machine precision; with noise on it should sit near 1
# for each observable, which is the end-to-end check that the error bars match the noise.