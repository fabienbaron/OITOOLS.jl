using OITOOLS

# Night limits for a given day
day = 17; month = 10; year = 2019;
lat = 34.2243836944;
lon = -118.0570313111;
ra = [20, 57, 59.4437981]'*[1.0, 1/60., 1/3600];
dec = [46, 28, 00.5731825]'*[1.0, 1/60., 1/3600];
~, UTC_set = sunrise_sunset(day, month, year, lat, lon;zenith=90.833333333);
UTC_rise, ~ = sunrise_sunset(day+1, month, year, lat, lon;zenith=90.833333333);
# Nights
#UT_set+1 to UT_rise-1
dark_UTCs = collect(range(UTC_set, UTC_rise, step=1.0/60));
#dark_hours = collect(range(UT_set+1, UT_rise-1, step=1.0/60));

#dark_hms = hours_to_hms(dark_hours)
dark_dates = Array(hcat([[year;month;day;(hours_to_hms(dark_UTCs))[i, :]] for i=1:size(dark_UTCs,1)]...)');
lst, ha = hour_angle_calc(dark_dates,lon,ra, ldir="E");

# Elevation limits
alt_limit = 45; # observe above 45 degrees elevation
alt, az = alt_az(dec, lat, ha); # result will be in degrees
plot(dark_UTCs, alt)
# good_alt = findall(alt.>alt_limit);
# ha = ha[good_alt];
# alt= alt[good_alt];
# az = az[good_alt]
# print("HA range based on CHARA telescope elevation limit: from ", ha[1], " to ",  ha[end] )

obs = read_obs_file("./data/default_obs"); # read defaults
obs.raep0[1] =  ra
obs.decep0[1] = dec
facility    = read_facility_file("./data/CHARA.txt");
#facility    = read_facility_file("./data/example_S2E1.txt");
combiner    = read_comb_file("./data/MIRC.txt");
wave        = read_wave_file("./data/MIRC_LOWH.txt");

h = ha' * pi / 12; #✓
δ = obs.decep0[1]/180*pi
l = facility.lat[1]/180*pi; #✓
λ = wave.lam#
ntel = facility.ntel[1]
station_xyz= hcat([facility.sta_xyz[(i*3-2):i*3] for i=1:ntel]...)'
nv2,v2_baselines,v2_stations,v2_indx,baseline_name         = get_v2_baselines(ntel,station_xyz,facility.tel_names);
#nuv,uv,u_M,v_M,w_M = get_uv(l, h, λ, δ, v2_baselines); uvplot(uv);
delay_geo   = geometric_delay(l,h,δ,v2_baselines)

# POP DELAY
include("../src/popranges.jl");
pop = [5,5,5,5,5,5]
delay_pop = zeros(nv2)
delay_airpath = zeros(nv2)
for i=1:nv2
    itel2 = v2_stations[2, i]
    itel1 = v2_stations[1, i]
    delay_pop[i]     = pop_array[ itel2 , pop[itel2]] -  pop_array[ itel1 , pop[itel1]]
    delay_airpath[i] = airpath[itel2] - airpath[itel1]
end
delay_carts = 0.5*( delay_geo .- delay_airpath .- delay_pop )

clf()
plot(dark_UTCs, delay_carts[5,:], label="Delay")
plot(dark_UTCs, alt, label = "Altitude")
legend()
xlabel("UT")
ylabel("EL/DELAY")
tight_layout()

has_delay = (delay_carts.>-43).&(delay_carts.<43)
#has_delay_all_baselines = vec(prod(has_delay, dims=1));
#findall(has_delay_all_baselines.>0)
