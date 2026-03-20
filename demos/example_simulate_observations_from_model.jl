using OITOOLS
using AstroTime
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
params = Float64[]

# Target info
target = read_obs_file("default_obs"); # read defaults (for OIFITS header)
target.target = "AZ Cyg"
target.raep0 =  [20, 57, 59.4437981]'*[1.0, 1/60., 1/3600] # set ra
target.decep0 = [46, 28, 00.5731825]'*[1.0, 1/60., 1/3600] # set dec


# Simulation info
facility    = read_facility_file("CHARA_new");
combiner    = read_comb_file("MIRC");
wavelength        = read_wave_file("MIRC_LOWH");

# Output file
out_file="./data/simulation-limb_darkened_disk.oifits";

simulate(facility, target, combiner, wavelength, dates, out_file, flat_model=model, flat_params=params);

#Check simulated data
data = (readoifits(out_file))[1,1]; # data can be split by wavelength, time, etc.
# uvplot(data)
# uvplot(data, color="wav")
# uvplot(data, color="mjd")
