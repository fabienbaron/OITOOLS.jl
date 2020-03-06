#
# Classic visibility functions (uniform disc, etc.)
#

using SpecialFunctions

#param  : model parameter
#param[1] = diameter in mas
#ρ = radius in uv space
function visibility_ud(param, uv;tol=1e-14)
ρ=sqrt.(uv[1,:].^2+uv[2,:].^2)
t = param[1]/2.0626480624709636e8*pi*ρ;
V = 2.0*besselj1.(t)./t
indx= findall(abs.(t).<tol)
if indx !=[]
    V[indx].=1.0;
end
return V
end

function dvisibility_ud(param, uv;tol=1e-14)
ρ=sqrt.(uv[1,:].^2+uv[2,:].^2)
dt_dp = pi*ρ/2.0626480624709636e8
t= param[1]*dt_dp
dV_dt = (t.*besselj0.(t)-2*besselj1.(t))./t.^2
return dV_dt.*dt_dp
end

function visibility_ellipse_uniform(param, uv) # i: inclination (deg), ϕ: semi-major axis orientation (deg)
ϵ = param[2];
ϕ = param[3]/180*pi;
ρ = sqrt.(ϵ^2*(uv[1,:].*cos.(ϕ)-uv[2,:].*sin.(ϕ) ).^2+(uv[2,:].*cos.(ϕ)+uv[1,:].*sin.(ϕ)).^2)
t = param[1]/2.0626480624709636e8*pi*ρ;
V = 2.0*besselj1.(t)./t
indx= findall(abs.(t).<tol)
if indx !=[]
    V[indx].=1.0;
end
return V
end

function visibility_ellipse_quad(param, uv) # i: inclination (deg), ϕ: semi-major axis orientation (deg)
ϵ = param[4];
ϕ = param[5]/180*pi;
ρ = sqrt.(ϵ^2*(uv[1,:].*cos.(ϕ)-uv[2,:].*sin.(ϕ) ).^2+(uv[2,:].*cos.(ϕ)+uv[1,:].*sin.(ϕ)).^2)
theta = param[1]/2.0626480624709636e8;
zeta = (pi*ρ*theta);
V = ((1.0-param[2]-param[3])*(besselj1.(zeta)./zeta)+((theta+2*param[3])/sqrt(2/pi))*(
(sqrt.(2 ./(pi*zeta)).*((sin.(zeta)./zeta)-cos.(zeta)))./zeta.^(3/2))-2*param[3]*
(besselj.(2, zeta)./zeta.^2))./(0.5-param[2]/6-param[3]/12)
return V
end

#power law
function visibility_ldpow(param, uv; maxk = 200)
ρ=sqrt.(uv[1,:].^2+uv[2,:].^2)
theta = param[1]/2.0626480624709636e8;
#f=k->gamma(param[2]/2+2)/(gamma(param[2]/2+k+2)*gamma(k+1))*(-0.25*(pi*theta*ρ).^2).^k; # note: can lead to NaN due to typical gamma blow up
f=k->(-1)^k*exp.(2k*log.(0.5*pi*theta*ρ).+(logabsgamma(param[2]/2+2)[1]-logabsgamma(param[2]/2+k+2)[1]-logabsgamma(k+1)[1]));
V = sum(map(f,collect(0:maxk)));
return V
end

#quadratic law
function visibility_ldquad(param,uv)
ρ=sqrt.(uv[1,:].^2+uv[2,:].^2)
theta = param[1]/2.0626480624709636e8;
zeta = (pi*ρ*theta);
V = ((1.0-param[2]-param[3])*(besselj1.(zeta)./zeta)+((theta+2*param[3])/sqrt(2/pi))*(
(sqrt.(2 ./(pi*zeta)).*((sin.(zeta)./zeta)-cos.(zeta)))./zeta.^(3/2))-2*param[3]*
(besselj.(2, zeta)./zeta.^2))./(0.5-param[2]/6-param[3]/12)
return V
end

#linear law
function visibility_ldlin(param,uv)
ρ=sqrt.(uv[1,:].^2+uv[2,:].^2)
theta = param[1]/2.0626480624709636e8;
zeta = (pi*ρ*theta);
V = ((1.0-param[2])*(besselj1.(zeta)./zeta)+theta/sqrt(2/pi)*(
(sqrt.(2 ./(pi*zeta)).*((sin.(zeta)./zeta)-cos.(zeta)))./zeta.^(3/2)))./(0.5-param[2]/6)
return V
end

# visibility of an annulus of unit flux
function visibility_annulus(r_in, r_out, uv;tol=1e-10)
ρ=sqrt.(uv[1,:].^2+uv[2,:].^2)
if abs(r_out - r_in)<tol #
  return zeros(Complex{Float64},length(ρ)); #shouldn't we use thin ring here ?
else
  return (1.0+0*im).*(visibility_ud([2*r_out],ρ)*r_out^2-visibility_ud([2*r_in],ρ)*r_in^2)/(r_out^2-r_in^2);
end

function visibility_thinring(param,uv)
ρ=sqrt.(uv[1,:].^2+uv[2,:].^2)
return besselj0.(pi*param[1]/2.0626480624709636e8*ρ)
end

# function visibility_ldsqrt(param,ρ)
# theta = param[1]/2.0626480624709636e8;
# bs1 = bessel(0, pi*ρ*theta, 2)
# bs2 = bessel(0.5, pi*ρ*theta, 2)
# bs3 = bessel(0.25, pi*ρ*theta, 2)
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

if visfunc == visibility_ellipse_uniform
 lbounds = [0.0,0.0,-180.0]
 hbounds = [1e9,1.0,180.0]
end


if visfunc == visibility_ellipse_quad
 lbounds = [0.0, -1.0, -1.0, 0.0, -180.0]
 hbounds = [1e9, 1.0, 1.0, 0.0, 180.0]
end

return lbounds, hbounds
end
