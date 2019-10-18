#
# This file includes general purpose functions
#

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

using PyCall

function query_simbad(targetname)
    return pyimport("astroquery.simbad").Simbad.query_object(targetname)
end

function ra_dec_from_simbad(targetname)
res=query_simbad(targetname)
ra = [parse(Float64, i) for i in split(get(get(res, "RA"),0))]
dec = [parse(Float64, i) for i in split(get(get(res, "DEC"),0))]
return ra, dec
end
