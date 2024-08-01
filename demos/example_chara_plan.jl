#
# This demonstrates how to calculate delays in a way that is 
# analogous with the chara_plan software used at GSU
#

using OITOOLS, Dates
using PyPlot

# To query simbad  uncomment the next line
targetname = "AZ Cyg"
radec = ra_dec_from_simbad(targetname)
ra, dec = (radec[1]'*[1.0, 1/60., 1/3600], radec[2]'*[1.0, 1/60., 1/3600])

# Night limits for a given date
obsdate = DateTime(2020, 6, 03);

# Observatory location
facility = read_facility_file("./data/CHARA.txt");
lat, lon = facility.lat[1], facility.lon[1]; #CHARA

lst_midnight, _ = hour_angle_calc(obsdate+Dates.Day(1)+Dates.Hour(7),lon, ra)
lst_midnight = lst_midnight[1];

# Dark Observability
_, UTC_set = sunrise_sunset(obsdate, lat, lon);
UTC_rise, _ = sunrise_sunset(obsdate+Dates.Day(1), lat, lon);
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

# In this example we're finding the best pops

tel_names= ["S1"  "S2"  "E1"  "E2"  "W1"  "W2"]
config = [1, 1, 1, 1, 1, 2]# telescopes to use = 1, to not use = 0, reference = 2
nbaselines, baseline_xyz , baseline_stations, baseline_name = get_baselines(facility, config=config);
delay_geo   = geometric_delay(lat,ha,dec,baseline_xyz);
include("../src/popranges.jl");
delay_airpath = [airpath[baseline_stations[2, i]] - airpath[baseline_stations[1, i]] for i=1:nbaselines]
delay_pop = zeros(nbaselines, 1, 5,5,5,5,5,5);
for i=1:5
    for j=1:5
        for k=1:5
            for l=1:5
                for m=1:5
                    for n=1:5
                        pop = [i,j,k,l,m,n]
                        delay_pop[:,1, i,j,k,l,m,n] = [pop_array[ baseline_stations[2, i] , pop[baseline_stations[2, i]]] -  pop_array[ baseline_stations[1, i] , pop[baseline_stations[1, i]]] for i=1:nbaselines]
                    end
                end
            end
        end
    end
end
delay_carts = 0.5*( delay_geo .- delay_airpath .- delay_pop ); # baselines x lst x pop1 x pop2 x pop3 xpop4 x pop5x pop6
has_delay = prod((delay_carts.>-43).&(delay_carts.<43), dims=1); # 1 = delay on all baselines
pop_score = reshape(sum(has_delay,dims=2),5,5,5,5,5,5); # assumes non-discontinous delay block, may be wrong

# Best 5 POP configs
minimum_minutes_on_object = 10; # we'll assume we need at least 10 minute on an object
indx_good = findall(pop_score.>minimum_minutes_on_object)
good_pops_loc = sortperm(pop_score[indx_good],rev=true)
good_pops = indx_good[good_pops_loc]
good_pops = good_pops[1:min(5, length(good_pops))] # only keep the 5 best ones
println("Good POP solutions:");
for j=1:length(good_pops)
println("\nSolution ranked #$j")
[print(string(tel_names[i], " - POP ", good_pops[j][i])," | ") for i=1:length(tel_names)];
end

best_pops = findmax(pop_score)
println("Best POP solution:");
[println(string(tel_names[i], " - POP ", best_pops[2][i])) for i=1:length(tel_names)];
#good_delay = findall(vec(has_delay[:,:,best_pops[2]]).>0)


# ALTERNATIVE EXAMPLE, IF WE REQUEST SPECIFIC POPS
pop = [1,2,4,5,2,5]
delay_airpath = [airpath[baseline_stations[2, i]] - airpath[baseline_stations[1, i]] for i=1:nbaselines]
delay_pop = [pop_array[ baseline_stations[2, i] , pop[baseline_stations[2, i]]] -  pop_array[ baseline_stations[1, i] , pop[baseline_stations[1, i]]] for i=1:nbaselines]
delay_carts = 0.5*( delay_geo .- delay_airpath .- delay_pop )
has_delay = (delay_carts.>-43).&(delay_carts.<43);
good_delay = findall(vec(prod(has_delay, dims=1)).>0);
gantt_onenight(targetname,obsdate, lst, lst_midnight, az, alt, good_alt, good_delay);tight_layout();

# LST vs delay
figure()
plot(lst, delay_carts[4,:], label="E2")
plot(lst, delay_carts[3,:], label="E1")
plt.axhline(y=alt_limit, label="Limit altitude = $(alt_limit)°", color=:red)
plot(lst, alt, label = "Altitude")
legend()
ylim(-5,90)
xlim(minimum(lst), maximum(lst))
xlabel("UT")
ylabel("EL/DELAY")
tight_layout()
