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
obs = read_obs_file("./data/example_obs_config.txt"); # read defaults
obs.raep0[1] =  [20, 57, 59.4437981]'*[1.0, 1/60., 1/3600] # UPDATE ra
obs.decep0[1] = [46, 28, 00.5731825]'*[1.0, 1/60., 1/3600] # UPDATE DEC

facility    = read_facility_file("./data/example_facility_config.txt");
combiner    = read_comb_file("./data/example_combiner_config.txt");
wave        = read_wave_file("./data/example_wave_config.txt");
v2m=1.0/100; v2a=1e-5; t3ampm=1.0/100; t3ampa=1e-6; t3phim=0.0; t3phia=0.5;
errors      = define_errors(v2m,v2a,t3ampm,t3ampa,t3phim,t3phia);

lst, hour_angles = hour_angle_calc(dates,facility.lon[1],obs.raep0[1]);
simulate_ha(facility, obs, combiner, wave, hour_angles, image_file, pixsize, errors, out_file);

#Compare simulated data to original data
data = (readoifits(out_file))[1,1]; # data can be split by wavelength, time, etc.
uvplot(data, fancy = true)

#
# DEBUGGING - WORK IN PROGRESS AFTER THIS POINT
#
ha = collect(range(-12, 12, step=1.0/60))
alt_limit = 35; # observe above 25 degrees elevation
alt, az = alt_az(obs.decep0[1], facility.lat[1], ha); # result will be in degrees
good_alt = findall(alt.>alt_limit);
ha = ha[good_alt[1]:good_alt[end]];
print("HA range based on CHARA telescope elevation limit: from ", ha[1], " to ",  ha[end] )

h = ha' * pi / 12; #✓
δ = obs.decep0[1]/180*pi
l = facility.lat[1]/180*pi; #✓
λ = wave.lam#

ntel = facility.ntel[1]
station_xyz= hcat([facility.sta_xyz[(i*3-2):i*3] for i=1:ntel]...)'
nv2,v2_baselines,v2_stations,v2_indx,baseline_name         = get_v2_baselines(ntel,station_xyz,facility.tel_names);
nuv,uv,u_M,v_M,w_M = get_uv(l, h, λ, δ, v2_baselines)
delay_geo   = geometric_delay(l,h,δ,v2_baselines)

# POP DELAY
include("../src/popranges.jl");
pop_tel1 = 1
pop_tel2 = 1
delay_pop = zeros(nv2)
delay_airpath = zeros(nv2)
for i=1:nv2
    delay_pop[i]     = pop_array[ v2_stations[2, i] , pop_tel2] -  pop_array[ v2_stations[1, i] , pop_tel1]
    delay_airpath[i] = airpath[v2_stations[2, i]] - airpath[v2_stations[1, i]]
end

# AIRMASS
delay_carts = -0.5*( delay_geo .+ delay_airpath .+ delay_pop )
has_delay = (delay_carts.>-43).&(delay_carts.<43)
has_delay_all_baselines = vec(prod(has_delay, dims=1));
findall(has_delay_all_baselines)
