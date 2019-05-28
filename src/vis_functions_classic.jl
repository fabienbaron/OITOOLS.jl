#
# Classic visibility functions (uniform disc, etc.)
#

using SpecialFunctions

#param  : model parameter
#param[1] = diameter in mas
#v_r = radius in uv space
function visibility_ud(param, v_r;tol=1e-14)
t = param[1]/2.0626480624709636e8*pi*v_r;
V = 2.0*besselj1.(t)./t
indx= findall(abs.(t).<tol)
if indx !=[]
    V[indx].=1.0;
end
return V
end

function dvisibility_ud(param, v_r;tol=1e-14)
dt_dp = pi*v_r/2.0626480624709636e8
t= param[1]*dt_dp
dV_dt = (t.*besselj0.(t)-2*besselj1.(t))./t.^2
return dV_dt.*dt_dp
end


#power law
function visibility_ldpow(param, v_r)
maxk=200;
theta = param[1]/2.0626480624709636e8;
#f=k->gamma(param[2]/2+2)/(gamma(param[2]/2+k+2)*gamma(k+1))*(-0.25*(pi*theta*v_r).^2).^k; # note: can lead to NaN due to typical gamma blow up
f=k->(-1)^k*exp.(lgamma(param[2]/2+2)-lgamma(param[2]/2+k+2)-lgamma(k+1).+2k*log.(0.5*pi*theta*v_r));
V = sum(map(f,collect(0:maxk)));
return V
end

#quadratic law
function visibility_ldquad(param,v_r)
theta = param[1]/2.0626480624709636e8;
zeta = (pi*v_r*theta);
V = ((1.0-param[2]-param[3])*(besselj1.(zeta)./zeta)+((theta+2*param[3])/sqrt(2/pi))*(
(sqrt.(2 ./(pi*zeta)).*((sin.(zeta)./zeta)-cos.(zeta)))./zeta.^(3/2))-2*param[3]*
(besselj.(2, zeta)./zeta.^2))./(0.5-param[2]/6-param[3]/12)
return V
end

#linear law
function visibility_ldlin(param,v_r)
theta = param[1]/2.0626480624709636e8;
zeta = (pi*v_r*theta);
V = ((1.0-param[2])*(besselj1.(zeta)./zeta)+theta/sqrt(2/pi)*(
(sqrt.(2 ./(pi*zeta)).*((sin.(zeta)./zeta)-cos.(zeta)))./zeta.^(3/2)))./(0.5-param[2]/6)
return V
end

# visibility of an annulus of unit flux
function visibility_annulus(r_in, r_out, v_r;tol=1e-10)
if abs(r_out - r_in)<tol #
  return zeros(Complex{Float64},length(v_r)); #shouldn't we use thin ring here ?
else
  return (1.0+0*im).*(visibility_ud([2*r_out],v_r)*r_out^2-visibility_ud([2*r_in],v_r)*r_in^2)/(r_out^2-r_in^2);
end

function visibility_thinring(param,v_r)
return besselj0.(pi*param[1]/2.0626480624709636e8*v_r)
end

# function visibility_ldsqrt(param,v_r)
# theta = param[1]/2.0626480624709636e8;
# bs1 = bessel(0, pi*v_r*theta, 2)
# bs2 = bessel(0.5, pi*v_r*theta, 2)
# bs3 = bessel(0.25, pi*v_r*theta, 2)
# x1 = ((1D0-alpha-beta)*(bs1(2)))/(pi*a*rho)
# x2 = (alpha*(sqrt(pi/2D0))*bs2(2))/((pi*a*rho)**(1.5D0))
# x3 = (beta*((2D0)**(0.25D0))*(gamma(1.25D0))*bs3(2))/((pi*a*rho)**(1.25D0))
# x4 = (0.5D0)-((1D0/6D0)*alpha)-((1D0/10D0)*beta)
# return (x1+x2+x3)/x4
# end


end

function init_bounds(visfunc)
#return default lower and upper bounds on parameters
lbounds = []
hbounds = []
if visfunc==visibility_ud
 lbounds = [0.0]
 hbounds = [1e9]
end

if visfunc==visibility_ldpow
 lbounds = [0.0, 0.0]
 hbounds = [1e9, 3.0]
end

if visfunc == visibility_ldquad
 lbounds = [0.0, -1.0, -1.0]
 hbounds = [1e9, 1.0, 1.0]
end

if visfunc == visibility_ldlin
 lbounds = [0.0, -1.0]
 hbounds = [1e9, 1.0]
end

return lbounds, hbounds
end
