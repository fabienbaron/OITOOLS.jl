using OITOOLS
using AstroTime
#Simulate an observation using an analytic model of your target
#In this example, observations start at UT 2018-08-13 at 3:00:00AM and last until 2018-08-13 8:30:00AM, with a period of 15 minutes
dates = collect(from_utc(2018,8,13,3,0,0.0):15minutes:from_utc(2018,8,13,8,30,0.0))

# Model info - here a simple limb-darkened disk
model = create_model(create_component(type="ldlin", name="Model"));
dispatch_params([3.0,0.15], model);

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

# Output file
out_file="./data/simulation-limb_darkened_disk.oifits";

simulate(facility, target, combiner, wavelength, dates, errors, out_file, model=model);

#Check simulated data
data = (readoifits(out_file))[1,1]; # data can be split by wavelength, time, etc.
uvplot(data)
uvplot(data, color="wav")
uvplot(data, color="mjd")
