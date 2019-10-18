using OITOOLS

#EXAMPLE 8
#Simulate an observation using an input image, given telescope parameters, and input hour angles

dates=[2018 8 13 5 13 56.7; 2018 8 13 5 15 56.7; 2018 8 13 5 20 56.7; 2018 8 13 5 25 56.7; 2018 8 13 5 30 56.7;
       2018 8 13 6 15 56.7; 2018 8 13 6 20 56.7; 2018 8 13 6 25 56.7; 2018 8 13 6 30 56.7; 2018 8 13 7 15 56.7;
       2018 8 13 7 20 56.7; 2018 8 13 7 25 56.7; 2018 8 13 7 30 56.7;]

# Object info
image_file="./data/2004true.fits";
pixsize=0.101;
out_file="./data/2004testsimulation.oifits";
obs = read_obs_file("./data/default_obs.txt"); # read defaults
obs.raep0[1] =  [20, 57, 59.4437981]'*[1.0, 1/60., 1/3600] # UPDATE ra
obs.decep0[1] = [46, 28, 00.5731825]'*[1.0, 1/60., 1/3600] # UPDATE DEC

facility    = read_facility_file("./data/CHARA.txt");
combiner    = read_comb_file("./data/MIRC.txt");
wave        = read_wave_file("./data/MIRC_LOWH.txt");
v2m=1.0/100; v2a=1e-5; t3ampm=1.0/100; t3ampa=1e-6; t3phim=0.0; t3phia=0.5;
errors      = define_errors(v2m,v2a,t3ampm,t3ampa,t3phim,t3phia);

lst, hour_angles = hour_angle_calc(dates,facility.lon[1],obs.raep0[1]);
simulate_ha(facility, obs, combiner, wave, hour_angles, image_file, pixsize, errors, out_file);

#Compare simulated data to original data
data = (readoifits(out_file))[1,1]; # data can be split by wavelength, time, etc.
uvplot(data, fancy = true)
