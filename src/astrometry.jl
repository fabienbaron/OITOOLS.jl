using AstroTime

function hours_to_date(obsdate, hours)
    h= floor.(hours)
    min = div.((hours-h)*3600.0, 60.0)
    sec = rem.((hours-h)*3600.0, 60.0)
    isec =  floor.(sec)
    msec = round.((sec-floor.(sec))*1000)
    return Dates.DateTime.(Dates.year(obsdate),Dates.month(obsdate),Dates.day(obsdate),h,min,isec,msec)
end



function hour_angle_calc(dates::Union{DateTime, Array{DateTime, 1}}, longitude::Float64, ra::Float64;ldir="E")
return hour_angle_calc(from_utc.(dates), longitude, ra, ldir=ldir)
end

function hour_angle_calc(dates, longitude::Float64, ra::Float64;ldir="E")

"""
This function calculates and returns the hour angle for the desired object given a RA, time, and longitude
of observations. This function assumes UTC.
"""

alpha= 1
if ldir == "W"
    alpha = -1.
elseif ldir == "E"
    alpha = 1.
end

if typeof(dates) == Epoch{InternationalAtomicTime, Float64} # transform single epochs into Array
    dates = [dates]
end

Y= year.(dates)
M= month.(dates)
D = day.(dates)
H= hour.(dates) # CHARA UTC offset
MIN= minute.(dates)
SEC= second.(dates)+millisecond.(dates)/1000

h_ad = alpha*longitude/15 #longitude in degrees, Measures the hours offset due to longitude

#Below: the code calculates first the Julian Date given the input time and then determines the GMST based
#on this JD.  This is then converted to LST and finally to Local Hour Angle (LHA or HA).  The final result
#is in terms of hours for both LST and HA.
jdn = floor.((1461*(Y .+4800 .+(M.-14)/12))/4+(367*(M .-2-12*((M.-14)/12)))/12-(3*((Y .+4900 +(M.-14)/12)/100))/4 .+D.-32075)
jdn = jdn + ((H.-12)/24)+(MIN/1440)+(SEC/86400)
jd0 = jdn .-2451545.0
t = jd0/36525.0
gmst = 24110.54841 .+ 8640184.812866*t + 0.093104*t.^2 - 6.2E-6*t.^3 + (1.00273790935 .+ 5.9e-11*t).*(H*3600 + MIN*60 + SEC) #seconds
gmst = gmst/3600 #hours

gmst_over = findall(gmst.>24)
gmst[gmst_over] -= (24*floor.(gmst[gmst_over]/24))
lst = gmst .+ h_ad
lst_under = findall(lst.<0)
lst_over = findall(lst.>24)
lst[lst_under] .+= 24
lst[lst_over] -= (24*floor.(lst[lst_over]/24))
hour_angle = lst .-ra
hour_angle[findall(hour_angle.<-12)] .+= 24; # had to add this, normal ???? FB
return lst,hour_angle
end

function mjd_to_utdate(mjd::Float64) # used in interactive plots in oiplot.jl to return utc of points TODO: replace this
    """
    This function calculates and returns the UTC given a specific mjd.  It is a julia implementation of the codes  by Peter J. Acklam found
        on http://www.radiativetransfer.org/misc/atmlabdoc/atmlab/time/mjd2date.html and linked pages, except for the error handling (for now)
    """
    jd = mjd + 2400000.5;
    #get year month and days
    intjd = floor(Int64,jd+0.5);
    a = intjd+32044;
    b = floor(Int64, (4 * a + 3) / 146097);
    c = a - floor(Int64, (b * 146097) / 4);
    d = floor(Int64,(4 * c + 3) / 1461);
    e = c - floor(Int64,(1461 * d) / 4);
    m = floor(Int64, (5 * e + 2) / 153);
    day   = e - floor(Int64, (153 * m + 2) / 5) + 1;
    month = m + 3 - 12 * floor(Int64, m / 10);
    year  = b * 100 + d - 4800 + floor(Int64, m / 10);
    #get hour min seconds
    fracmjd = mjd-floor(Int64, mjd); #days;
    secs = 86400*fracmjd;
    hour = trunc(Int64,secs/3600);
    secs= secs - 3600*hour;  #remove hour
    mins = Int(trunc(secs/60)); #minutes
    secs = secs - 60*mins; #seconds
    msecs = round(Int64,1000*(secs - floor(Int64, secs)));
    secs = floor(Int64, secs);
    return Dates.DateTime(year, month, day, hour, mins, secs, msecs)
end


function opd_limits(base, alt, az)
    s = hcat(cos.(alt).*cos.(az),cos.(alt).*sin.(az),sin.(alt))'
    opd = dot(base,s) #in meters (?)
    return opd
end

function alt_az(dec_deg,lat_deg, ha_hours) #returns alt, az in degrees
    dec = dec_deg*pi/180;
    ha = ha_hours*pi/12
    lat = lat_deg*pi/180
    # Simple version
    alt = asin.(sin(dec)*sin(lat).+cos(dec)*cos(lat)*cos.(ha))
    az = atan.((-cos(dec)*sin.(ha))./(sin(dec)*cos(lat).-cos(dec)*cos.(ha)*sin(lat)))
    return alt*180/pi, mod.(az*180/pi.+360, 360)
end

function geometric_delay(lat_deg,h_deg,δ_deg,baselines)
    l = lat_deg/180*pi
    h = h_deg'*pi/12
    δ = δ_deg/180*pi
    Δgeo =  -(sin(l)*cos(δ)*cos.(h).-cos(l)*sin(δ) ).*baselines[1,:]-(cos(δ)*sin.(h)) .* baselines[2,:]+ (cos(l)*cos(δ)*cos.(h) .+ sin(l)*sin(δ)).*baselines[3,:]
    return Δgeo
end

function dayofyear(day, month, year)
N = floor(275 * month / 9) - floor((month + 9) / 12) * (1 + floor((year - 4 * floor(year / 4) + 2) / 3)) + day - 30
return N
end

function sunrise_sunset(obsdate::DateTime, latitude, longitude;zenith=102.0)
    return sunrise_sunset(from_utc(obsdate), latitude, longitude,zenith=zenith)
end

function sunrise_sunset(obsdate, latitude, longitude;zenith=102.0)
#Source:
#	Almanac for Computers, 1990
#	published by Nautical Almanac Office
#	United States Naval Observatory
#	Washington, DC 20392
#
#Inputs:
#	day, month, year:      date of sunrise/sunset
#	latitude, longitude:   location for sunrise/sunset
#	zenith:                Sun's zenith for sunrise/sunset
#	  offical sunrise/sunset     = 90 degrees 50'
#	  civil        twilight = 96 degrees
#	  nautical     twilight  = 102 degrees
#	  astronomical twilight = 108 degrees

#	NOTE: longitude is positive for East and negative for West

Y   = year.(obsdate)
M  = month.(obsdate)
D    = day.(obsdate)

dtr=pi/180.
rtd = 180.0/pi
#1. first calculate the day of the year
N = dayofyear(D, M, Y)

#2. convert the longitude to hour value and calculate an approximate time
lngHour = longitude / 15
t_rise = N .+ ((6 - lngHour) / 24)
t_set = N .+ ((18 - lngHour) / 24)

#3. calculate the Sun's mean anomaly
M_rise = 0.9856 * t_rise .- 3.289
M_set = 0.9856 * t_set .- 3.289

#4. calculate the Sun's true longitude
L_rise = mod.(M_rise +  1.916 * sin(M_rise*dtr) +   0.020 * sin(2 * M_rise*dtr) .+ 282.634, 360)
L_set  = mod.(M_set  +  1.916 * sin(M_set*dtr)  +   0.020 * sin(2 * M_set*dtr)  .+ 282.634, 360)
#	NOTE: L potentially needs to be adjusted into the range [0,360) by adding/subtracting 360

#5a. calculate the Sun's right ascension
RA_rise = mod.(rtd*atan(0.91764 * tan(L_rise*dtr)),360)
RA_set = mod.(rtd*atan(0.91764 * tan(L_set*dtr)),360)

#	NOTE: RA potentially needs to be adjusted into the range [0,360) by adding/subtracting 360
#5b. right ascension value needs to be in the same quadrant as L
RA_rise += 90*(floor( L_rise/90) -  floor(RA_rise/90))
RA_set  += 90*(floor( L_set/90) -  floor(RA_set/90))

#5c. right ascension value needs to be converted into hours
RA_rise = RA_rise / 15
RA_set = RA_set / 15

#6. calculate the Sun's declination
sinDec_rise = 0.39782 * sin(L_rise*dtr)
cosDec_rise = cos(asin(sinDec_rise))

sinDec_set = 0.39782 * sin(L_set*dtr)
cosDec_set = cos(asin(sinDec_set))


#7a. calculate the Sun's local hour angle
cosH_rise = (cos(zenith) - (sinDec_rise * sin(dtr*latitude))) / (cosDec_rise * cos(dtr*latitude))
#if (cosH_rise >  1)
#  print("the sun never rises on this location");
#elseif (cosH_rise < -1)
#  print("the sun never sets on this location");
#end
cosH_set = (cos(zenith) - (sinDec_set * sin(dtr*latitude))) / (cosDec_set * cos(dtr*latitude))

#7b. finish calculating H and convert into hours
H_rise = (360 - rtd*acos(dtr*cosH_rise))/15
H_set = rtd*acos(dtr*cosH_set)/15

#8. calculate local mean time of rising/setting
T_rise = H_rise + RA_rise - 0.06571 * t_rise - 6.622
T_set = H_set + RA_set - 0.06571 * t_set - 6.622

#9. adjust back to UTC
UT_rise = T_rise - lngHour
UT_set = T_set - lngHour

#	NOTE: UT potentially needs to be adjusted into the range [0,24) by adding/subtracting 24
#10. convert UT value to local time zone of latitude/longitude
#localT_rise = UT_rise + localOffse
#localT_rise = UT_rise + localOffse
return UT_rise, UT_set
end


# Aspro-type plot
function gantt_onenight(targetname,obsdate, lst, lst_midnight, az, alt, good_alt, good_delay)
# Gantt plot
fig = figure(figsize=(20,10))
ax = fig.add_subplot(111)
ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=1))
ax.xaxis.set_minor_locator(matplotlib.dates.MinuteLocator(interval=15))
ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
offset = 4# 0 if we want to see delay during the day, 4 to keep mostly dark times
xlim(hours_to_date(obsdate, lst_midnight[1]-12+offset), hours_to_date(obsdate+Dates.Day(1), lst_midnight[1]-12-offset) ); xlabel("LST")
ylim(0, 10);
df = DateFormat("h:m");
title(string("Observing night: ", Date(obsdate), " - ", Date(obsdate+Dates.Day(1)), " -- Target: ", targetname))
fig.autofmt_xdate(bottom=0.2,rotation=30,ha="right"); ax.xaxis_date(); grid()
plt.axvline(x=hours_to_date(obsdate, lst_midnight[1]), color=:red) # nextday transition # add label?
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
ax.barh(5, end_date - start_date, left=start_date, height=1.5, align="center", color=:orange, label="Altitude",zorder=3)
start_date = hours_to_date(obsdate,lst[good_delay[1]])
end_date   = hours_to_date(obsdate, lst[good_delay[end]])
ax.barh(2, end_date - start_date, left=start_date, height=2, align="center", color=:blue, label="In Delay", zorder=3)
text(start_date, 2, Dates.format(start_date, dateformat"H:M") , rotation=90,va="center", ha="right", color=:black)
text(end_date, 2, Dates.format(end_date, dateformat"H:M"), rotation=90,va="center",ha="left",color=:black)
text(start_date, 3.3, round(Int64,az[good_delay[1]]),va="top", ha="center", color=:black) # add label to show it's az ?
text(start_date, 0.7, round(Int64,alt[good_delay[1]]),va="bottom",ha="center",color=:black)
text(end_date, 3.3, round(Int64,az[good_delay[end]]),va="top", ha="right", color=:black)
text(end_date, 0.7, round(Int64,alt[good_delay[end]]),va="bottom",ha="right",color=:black)
yticks([2],[targetname])
legend()
tight_layout()
end

function get_baselines(facility; config = []) # similar to get_v2_baselines in simulate.jl, but here for planning
    # determine baselines and make necessary arrays
    N = facility.ntel[1]
    station_xyz= hcat([facility.sta_xyz[(i*3-2):i*3] for i=1:N]...)'

    # Use all scopes by default, without a specific reference cart
    if config == []
        config = ones(N)
    end

    # Example for CHARA, the facility.sta_names ["S1"  "S2"  "E1"  "E2"  "W1"  "W2"]
    # config = [1, 1, 0, 1, 1, 2]# telescopes to use = 1, to not use = 0, reference = 2
    iref  = findall(config .==2 )
    use_ref = ( length(iref) == 1 ) # we discard several refs as a mistake
    itels =  findall(config .== 1 )

    list_one = []
    if use_ref
        list_one = iref
    else
        list_one = itels
    end
    list_two = itels

    # Note: this should be done with push!
    nbaselines = Int64(N*(N-1)/2);
    baseline_xyz = Array{Float64}(undef,3,nbaselines);
    baseline_stations  = Array{Int64}(undef,2,nbaselines);
    ind = 1
    for i in list_one
      for j in list_two
          if i!=j
            baseline_xyz[:,ind] .= station_xyz[j,:]-station_xyz[i,:];
            baseline_stations[:,ind] = [i,j];
            ind += 1
            end
        end
    end

    nbaselines = ind-1;
    baseline_xyz = baseline_xyz[:,1:nbaselines]
    baseline_stations = baseline_stations[:,1:nbaselines];

    baseline_names = Array{String}(undef,nbaselines);
    for i=1:nbaselines
        baseline_names[i]=string(facility.sta_names[baseline_stations[1,i]],"-",facility.sta_names[baseline_stations[2,i]])
    end

    return nbaselines, baseline_xyz, baseline_stations, baseline_names
end
