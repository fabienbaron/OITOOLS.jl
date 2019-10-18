function hours_to_hms(hours)
    h= floor.(hours)
    min = div.((hours-h)*3600.0, 60.0)
    sec = rem.((hours-h)*3600.0, 60.0)
    return hcat(h, min, sec)
end



function hour_angle_calc(dates::Union{Array{Any,2},Array{Float64,2}},longitude::Float64, ra::Float64;dst="no",ldir="E",timezone="UTC")

"""
This function calculates and returns the hour angle for the desired object given a RA, time, and longitude
of observations.  The program assumes UTC but changing the timezone argument to "local" will adjust the
input times accordingly.

Arguments:
Manual Inputs:
* date [2018  3 5 21 13 56.7; 2018 3 5 21 14 12.7] day,month,year,hour,min,sec: Correction to UTC is done within code if necessary/dictated.
* longitude: degree (in decimal form) value of the longitudinal location of the telescope.
* ra: [h,m,s] array of the right ascension of the desired object.

Optional Inputs:
* dst: ["yes","no"] whether daylight savings time needs to be accounted for.
* ldir: ["W","E"] defines the positive direction for change in longitude (should be kept to W for best applicability).
* timezone: ["UTC","local"] states whether the time inputs are UTC or not.  If not the code will adjust the times as needed with the given longitude.

Accuracy:
* With preliminary testing the LST returned is accurate to within a few minutes when compared with other calculators (MORE TESTING AND QUANTITATIVE ERROR NEEDS TO BE ESTABLISHED).
"""


if ldir == "W"
    alpha = -1.
elseif ldir == "E"
    alpha = 1.
end

years=Int.(dates[:,1])
months=Int.(dates[:,2])
days = Int.(dates[:,3])
hours=Int.(dates[:,4])
minutes=Int.(dates[:,5])
seconds=Float64.(dates[:,6])

h_ad = alpha*longitude/15 #Measures the hours offset due to longitude


if timezone != "UTC"
    if dst=="no"
        hours .-= h_ad
    elseif dst=="yes"
        h_ad += 1
        hours .-= h_ad
    end
    wrap_hours_over = findall(hours.>24)
    hours[wrap_hours_over] .-= 24
    days[wrap_hours_over] .+= 1

    wrap_hours_under = findall(hours.<0)
    hours[wrap_hours_under] +=24
    days[wrap_hours_under] -= 1
end

#Below: the code calculates first the Julian Date given the input time and then determines the GMST based
#on this JD.  This is then converted to LST and finally to Local Hour Angle (LHA or HA).  The final result
#is in terms of hours for both LST and HA.
jdn = floor.((1461*(years .+4800 .+(months.-14)/12))/4+(367*(months .-2-12*((months.-14)/12)))/12-(3*((years .+4900 +(months.-14)/12)/100))/4 .+days.-32075)
jdn = jdn + ((hours.-12)/24)+(minutes/1440)+(seconds/86400)
jd0 = jdn .-2451545.0
#gmst = 6.697374558 + (0.06570982441908*(jd0-(hour/24.))) + (1.00273790935*hour) + 0.000026*(jd0/36525)^2
t = jd0/36525.0
H = 24110.54841 .+ 8640184.812866*t + 0.093104*t.^2 - 6.2E-6*t.^3
w = 1.00273790935 .+ 5.9e-11*t
tt = hours*3600 + minutes*60 + seconds
gmst = H + w.*tt #seconds
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


function mjd_to_utc(mjd)
    """
    This function calculates and returns the UTC given a specific mjd.  It is a julia implementation of the codes  by Peter J. Acklam found
        on http://www.radiativetransfer.org/misc/atmlabdoc/atmlab/time/mjd2date.html and linked pages, except for the error handling (for now)
    """
    jd = mjd + 2400000.5;
    #get year month and days
    intjd = floor(jd+0.5);
    a = intjd+32044;
    b = floor((4 * a + 3) / 146097);
    c = a - floor((b * 146097) / 4);
    d = floor((4 * c + 3) / 1461);
    e = c - floor((1461 * d) / 4);
    m = floor((5 * e + 2) / 153);
    day   = e - floor((153 * m + 2) / 5) + 1;
    month = m + 3 - 12 * floor(m / 10);
    year  = b * 100 + d - 4800 + floor(m / 10);
    #get hour min seconds
    fracmjd = mjd-floor(mjd); #days;
    secs = 86400*fracmjd;
    hour = trunc(secs/3600);
    secs= secs-3600*hour;  #remove hour
    mins = trunc(secs/60); #minutes
    secs = secs -60*mins; #seconds
    return secs,mins,hour,day,month,year
end

# function dates_to_jd(dates::Union{Array{Any,2},Array{Float64,2}})
#     """
#     This function calculates and returns the hour angle for the desired object given a RA, time, and longitude
#     of observations.  The program assumes UTC but changing the timezone argument to "local" will adjust the
#     input times accordingly.
#
#     Arguments:
#     Manual Inputs:
#     * date [2018  3 5 21 13 56.7; 2018 3 5 21 14 12.7] day,month,year,hour,min,sec: Correction to UTC is done within code if necessary/dictated.
#     * longitude: degree (in decimal form) value of the longitudinal location of the telescope.
#     * ra: [h,m,s] array of the right ascension of the desired object.
#
#     Optional Inputs:
#     * dst: ["yes","no"] whether daylight savings time needs to be accounted for.
#     * ldir: ["W","E"] defines the positive direction for change in longitude (should be kept to W for best applicability).
#     * timezone: ["UTC","local"] states whether the time inputs are UTC or not.  If not the code will adjust the times as needed with the given longitude.
#
#     Accuracy:
#     * With preliminary testing the LST returned is accurate to within a few minutes when compared with other calculators (MORE TESTING AND QUANTITATIVE ERROR NEEDS TO BE ESTABLISHED).
#     """
#     if ldir == "W"
#         alpha = -1.
#     elseif ldir == "E"
#         alpha = 1.
#     end
#     years=Int.(dates[:,1])
#     months=Int.(dates[:,2])
#     days = Int.(dates[:,3])
#     hours=Int.(dates[:,4])
#     minutes=Int.(dates[:,5])
#     seconds=Float64.(dates[:,6])
#     h_ad = alpha*longitude/15 #Measures the hours offset due to longitude
#     if timezone != "UTC"
#         if dst=="no"
#             hours .-= h_ad
#         elseif dst=="yes"
#             h_ad += 1
#             hours .-= h_ad
#         end
#         wrap_hours_over = findall(hours.>24)
#         hours[wrap_hours_over] .-= 24
#         days[wrap_hours_over] .+= 1
#
#         wrap_hours_under = findall(hours.<0)
#         hours[wrap_hours_under] +=24
#         days[wrap_hours_under] -= 1
#     end
#
#     #Below: the code calculates first the Julian Date given the input time and then determines the GMST based
#     #on this JD.  This is then converted to LST and finally to Local Hour Angle (LHA or HA).  The final result
#     #is in terms of hours for both LST and HA.
#     jdn = floor.((1461*(years .+4800 .+(months.-14)/12))/4+(367*(months .-2-12*((months.-14)/12)))/12-(3*((years .+4900 +(months.-14)/12)/100))/4 .+days.-32075)
#     jdn = jdn + ((hours.-12)/24)+(minutes/1440)+(seconds/86400)
#     jd0 = jdn .-2451545.0
#     return jd0
# end


# function jd_to_hour_angle(jd::Array{Float64,1},tel_longitude::Float64, obj_ra::Float64;dst="no",ldir="W",timezone="UTC")
#
# """
# This function calculates and returns the hour angle for the desired object given the object RA, the JD time, and longitude
# of observations.
# Arguments:
# Manual Inputs:
# * longitude: degree (in decimal form) value of the longitudinal location of the telescope.
# * ra: [h,m,s] array of the right ascension of the desired object.
# Accuracy:
# * With preliminary testing the LST returned is accurate to within a few minutes when compared with other calculators (MORE TESTING AND QUANTITATIVE ERROR NEEDS TO BE ESTABLISHED).
# """
#
# t = jd/36525.0 # TODO check that correct jd here (= rjd - 51545.0) ???
# # Greenwich Mean Sidereal Time at Oh UT (in seconds)
# H = 24110.54841 .+ 8640184.812866*t + 0.093104*t.^2 - 6.2e-6*t.^3
# w = 1.00273790935 .+ 5.9e-11*t  # what's this ???
# tt = hours*3600 + minutes*60 + seconds # what's this ???
# gmst = H + w.*tt #seconds
# # Greenwich Mean Sidereal Time in hours
# gmst = gmst/3600 #hours
# gmst_over = findall(gmst.>24)
# gmst[gmst_over] -= (24*floor.(gmst[gmst_over]/24))
# # Local Mean Sidereal Time in hours at 0h UT (longitude correction)
# lst = gmst .+ longitude/15
# lst_under = findall(lst.<0)
# lst_over = findall(lst.>24)
# lst[lst_under] .+= 24
# lst[lst_over] -= (24*floor.(lst[lst_over]/24))
# # HA of star at 0h UT on RJD
# hour_angle = lst .-obj_ra
# return lst, hour_angle
# end


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

function geometric_delay(l,h,δ,baselines)
    Δgeo =  -(sin(l)*cos(δ)*cos.(h).-cos(l)*sin(δ) ).*baselines[1,:]-(cos(δ)*sin.(h)) .* baselines[2,:]+ (cos(l)*cos(δ)*cos.(h) .+ sin(l)*sin(δ)).*baselines[3,:]
    return Δgeo
end


#function cart_delay(baselines)
#    Δcarts = -0.5*(geometric_delay(l,h,δ,baselines) - airpath_delay(baselines) + pop_delay(baselines))
#    return Δcarts
#end


function dayofyear(day, month, year)
N = floor(275 * month / 9) - floor((month + 9) / 12) * (1 + floor((year - 4 * floor(year / 4) + 2) / 3)) + day - 30
return N
end

function sunrise_sunset(day, month, year, latitude, longitude;zenith=102.0)
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

dtr=pi/180.
rtd = 180.0/pi
#1. first calculate the day of the year
N = dayofyear(day, month, year)

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
