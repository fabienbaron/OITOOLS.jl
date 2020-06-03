using PyPlot, PyCall, Dates
include("../src/OITOOLS.jl")
using Main.OITOOLS
# To query simbad  uncomment the next line
targetname = "AZ Cyg"
radec = ra_dec_from_simbad(targetname)
ra, dec = (radec[1]'*[1.0, 1/60., 1/3600], radec[2]'*[1.0, 1/60., 1/3600])

# Night limits for a given day
obsdate = DateTime(2020, 6, 03);

# Observatory location
facility = read_facility_file("./data/CHARA.txt");
lat, lon = facility.lat[1], facility.lon[1]; #CHARA

lst_midnight, ~ = hour_angle_calc(obsdate+Dates.Day(1)+Dates.Hour(7),lon, ra)

# Dark Observability
~, UTC_set = sunrise_sunset(obsdate, lat, lon);
UTC_rise, ~ = sunrise_sunset(obsdate+Dates.Day(1), lat, lon);
dark_offset = 0; #in hours, set offset to 1 or 2 for deeper dark
utc = collect(range(UTC_set+dark_offset, UTC_rise-dark_offset, step=1.0/60));
lst, ha = hour_angle_calc(hours_to_date(obsdate, utc),lon,ra);

# Elevation limits
alt_limit = 45; # observe above 45 degrees elevation
alt, az = alt_az(dec, lat, ha); # result will be in degrees

plot(utc, alt, label="Object altitude", color=:blue)
plt.axhline(y=alt_limit, label="Limit altitude", color=:red)
xlabel("UTC"); title("Altitude vs time")
legend()
good_alt = findall(alt.>alt_limit);
printstyled("HA range based on CHARA telescope elevation limit: from ", ha[good_alt][1], " to ",  ha[good_alt][end], color=:red )

# Gantt plot
fig = figure(figsize=(20,10))
ax = fig.add_subplot(111)
ax.xaxis.set_major_locator(pyimport("matplotlib.dates").HourLocator(interval = 2))
ax.xaxis.set_major_formatter(pyimport("matplotlib.dates").DateFormatter("%H:%M"))
xlim(hours_to_date(obsdate, lst_midnight[1]-12), hours_to_date(obsdate+Dates.Day(1), lst_midnight[1]-12) ); xlabel("LST")
ylim(0, 10)
title(string("Observing night: ", Date(obsdate), " - ", Date(obsdate+Dates.Day(1)), " Target: ", targetname))
fig.autofmt_xdate(bottom=0.2,rotation=30,ha="right"); ax.xaxis_date(); grid()
plt.axvline(x=hours_to_date(obsdate, lst_midnight[1]), color=:red) # nextday transition
tight_layout()
start_date = hours_to_date(obsdate,lst[1]+1.5)
end_date   = hours_to_date(obsdate, lst[end]-1.5 )
ax.barh(5, end_date - start_date, left=start_date, height=10, align="center", color=:lightgray, alpha = 0.75)
start_date = hours_to_date(obsdate,lst[1]+2)
end_date   = hours_to_date(obsdate, lst[end]-2 )
ax.barh(5, end_date - start_date, left=start_date, height=10, align="center", color=:lightgray, alpha = 0.75)
start_date = hours_to_date(obsdate,lst[1]+3)
end_date   = hours_to_date(obsdate, lst[end]-3 )
ax.barh(5, end_date - start_date, left=start_date, height=10, align="center", color=:gray, alpha = 0.75)

start_date = hours_to_date(obsdate,lst[good_alt[1]])
end_date   = hours_to_date(obsdate, lst[good_alt[end]])
ax.barh(5, end_date - start_date, left=start_date, height=1.5, align="center", color=:orange, label="Altitude")
legend()

ntel = facility.ntel[1]
station_xyz= hcat([facility.sta_xyz[(i*3-2):i*3] for i=1:ntel]...)'
nv2,v2_baselines,v2_stations,v2_indx,baseline_name         = get_v2_baselines(ntel,station_xyz,facility.tel_names);
delay_geo   = geometric_delay(lat,ha,dec,v2_baselines);

# POP DELAY
include("../src/popranges.jl");
pop = [1,2,4,5,2,5]
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
has_delay_times = findall(has_delay_all_baselines.>0)
# on first order length(findall(has_delay_all_baselines.>0)) is good
start_date = hours_to_date(obsdate,lst[has_delay_times[1]])
end_date   = hours_to_date(obsdate, lst[has_delay_times[end]])
ax.barh(2, end_date - start_date, left=start_date, height=1.5, align="center", color=:blue, label="In Delay")
legend()

# LST vs delay
figure()
plot(lst, delay_carts[1,:], label="Delay")
plt.axhline(y=alt_limit, label="Limit altitude", color=:red)
plot(lst, alt, label = "Altitude")
legend()
xlabel("UT")
ylabel("EL/DELAY")
tight_layout()
