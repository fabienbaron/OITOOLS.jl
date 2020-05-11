using OITOOLS
using PyPlot

#Note: to install astroquery
# using Conda
# Conda.add("astropy", channel="astropy")
# Conda.add("astroquery", channel="astropy")

# To query simbad  uncomment the next line
#radec = ra_dec_from_simbad("Az Cyg")
radec = ([20.0, 57.0, 59.4437], [46.0, 28.0, 0.573])
ra = radec[1]'*[1.0, 1/60., 1/3600];
dec = radec[2]'*[1.0, 1/60., 1/3600];


# Night limits for a given day
day = 21; month = 8; year = 2018;
#obsdate = Date(2019, 10, 17);

# Observatory location
lat = 34.2243836944;
lon = -118.0570313111;

#
#lst_midnight, ~ = hour_angle_calc(Array([year, month, day+1, 0.0, 0.0, 0.0 ]'),lon, 0.0, timezone="local")
lst_midnight, ~ = hour_angle_calc(Array([year, month, day+1, 7.0, 0.0, 0.0 ]'),lon, 0.0)

# Dark Observability
~, UTC_set = sunrise_sunset(day, month, year, lat, lon);
UTC_rise, ~ = sunrise_sunset(day+1, month, year, lat, lon);
dark_offset = 0; #in hours, set offset to 1 or 2 for deeper dark
utc = collect(range(UTC_set+dark_offset, UTC_rise-dark_offset, step=1.0/60));
dates = Array(hcat([[year;month;day;(hours_to_hms(utc))[i, :]] for i=1:size(utc,1)]...)');
lst, ha = hour_angle_calc(dates,lon,ra);

# Elevation limits
alt_limit = 45; # observe above 45 degrees elevation
alt, az = alt_az(dec, lat, ha); # result will be in degrees
good_alt = findall(alt.>alt_limit);
# ha = ha[good_alt];
# alt= alt[good_alt];
# az = az[good_alt]
# print("HA range based on CHARA telescope elevation limit: from ", ha[1], " to ",  ha[end] )

plot(utc, alt, label="Object altitude", color=:blue)
plt.axhline(y=alt_limit, label="Limit altitude", color=:red)
legend()


# Gantt plot
fig = figure()
ax = fig.add_subplot(111)
xlim(lst_midnight[1]-12.0, lst_midnight[1]+12.0)
ylim(0, 10)
grid()
ax.broken_barh([(lst[1]-1, 24+lst[end] - lst[1]+2)],(0, 10), color=:lightgray)
ax.broken_barh([(lst[1], 24+lst[end] - lst[1])],(0, 10), color=:gray)
ax.broken_barh([(lst[good_alt[1]], 24+lst[good_alt[end]] - lst[good_alt[1]])],(9.5, 0.3), color=:orange, label="Horizon")
legend()



facility    = read_facility_file("./data/CHARA.txt");
combiner    = read_comb_file("./data/MIRC.txt");
wave        = read_wave_file("./data/MIRC_LOWH.txt");

ntel = facility.ntel[1]
station_xyz= hcat([facility.sta_xyz[(i*3-2):i*3] for i=1:ntel]...)'
nv2,v2_baselines,v2_stations,v2_indx,baseline_name         = get_v2_baselines(ntel,station_xyz,facility.tel_names);
#nuv,uv,u_M,v_M,w_M = get_uv(l, h, λ, δ, v2_baselines); uvplot(uv);
delay_geo   = geometric_delay(lat/180*pi,ha'*pi/12,dec/180*pi,v2_baselines)

# POP DELAY
include("../src/popranges.jl");
pop = [1,2,2,3,2,5]
delay_pop = zeros(nv2)
delay_airpath = zeros(nv2)
for i=1:nv2
    itel2 = v2_stations[2, i]
    itel1 = v2_stations[1, i]
    delay_pop[i]     = pop_array[ itel2 , pop[itel2]] -  pop_array[ itel1 , pop[itel1]]
    delay_airpath[i] = airpath[itel2] - airpath[itel1]
end
delay_carts = 0.5*( delay_geo .- delay_airpath .- delay_pop )
has_delay = (delay_carts.>-43).&(delay_carts.<43);
has_delay_all_baselines = vec(prod(has_delay, dims=1));
findall(has_delay_all_baselines.>0)
# on first order length(findall(has_delay_all_baselines.>0)) is good


clf()
plot(dark_UTCs, delay_carts[1,:], label="Delay")
plot(dark_UTCs, alt, label = "Altitude")
legend()
xlabel("UT")
ylabel("EL/DELAY")
tight_layout()
