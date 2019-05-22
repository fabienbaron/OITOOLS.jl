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

using SparseArrays
using LinearAlgebra
function x_start_from_V2_dft(data, dft)
# estimate V and 1/sigma_V^2 from V2 and V2_err using equation 3.98a in Data Analysis (Sivia/Skilling)
V2 = 0.5*( data.v2 +  sqrt.(data.v2.^2+2*data.v2_err.^2))
V = sqrt.(V2)
W = spdiagm(0=>(1.0./V2+2*(3*V2-data.v2)./(data.v2_err.^2))) # 1/sigma^2
H = dft[1:data.nv2, :];
y = real(H'*(W*V)); y=y.*(y.>=0);imshow(reshape(y,64,64))
y = real(H'*(W*H)+1e9*sparse(1.0I, 4096,4096))\(real(H'*(W*V))); y=y.*(y.>=0);imshow(reshape(y,64,64))
y=y.*(y.>=0)

nx = 64; o = ones(nx); D_1D = spdiagm(-1=>-o[1:nx-1],0=>o); D = [kron(spdiagm(0=>ones(nx)), D_1D) ;  kron(D_1D, spdiagm(0=>ones(nx)))];
DtD = D'*D;
y = real(H'*(W*H)+1e9*DtD+1e9*sparse(1.0I, 4096,4096))\(real(H'*(W*V))); y=y.*(y.>=0);imshow(reshape(y,64,64));chi2_dft_f(y, dft, data)
end
