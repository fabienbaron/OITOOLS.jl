#
# This file includes general purpose functions
#

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
