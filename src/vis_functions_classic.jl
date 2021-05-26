#
# Classic visibility functions (uniform disc, etc.)
#
using SpecialFunctions

#param  : model parameter
#param[1] = diameter in mas
#ρ = radius in uv space
function visibility_ud(param, uv::Array{Float64,2})
ρ=sqrt.(uv[1,:].^2+uv[2,:].^2)
V = jinc.(param[1]/2.0626480624709636e8*ρ)
return V
end

# function dvisibility_ud(param, uv::Array{Float64,2})
# ρ=sqrt.(uv[1,:].^2+uv[2,:].^2)
# dt_dp = pi*ρ/2.0626480624709636e8
# t= param[1]*dt_dp
# dV_dt = (t.*besselj0.(t)-2*besselj1.(t))./t.^2   # ouch singularity!
# return dV_dt.*dt_dp
# end

function visibility_ellipse_uniform(param, uv::Array{Float64,2}) # ϕ: semi-major axis orientation (deg)
ϵ = param[2];
ϕ = param[3]/180*pi;
ρ = sqrt.(ϵ^2*(uv[1,:].*cos.(ϕ)-uv[2,:].*sin.(ϕ) ).^2+(uv[2,:].*cos.(ϕ)+uv[1,:].*sin.(ϕ)).^2)
t = param[1]/2.0626480624709636e8*ρ;
V = jinc.(t)
return V
end


function visibility_ellipse_quad(param, uv::Array{Float64,2};tol=1e-6) #ϵ: ellipticity  i: inclination (deg), ϕ: semi-major axis orientation (deg)
ϵ = param[4];
ϕ = param[5]/180*pi;
ρ = sqrt.(ϵ^2*(uv[1,:].*cos.(ϕ)-uv[2,:].*sin.(ϕ) ).^2+(uv[2,:].*cos.(ϕ)+uv[1,:].*sin.(ϕ)).^2)
θ = param[1]/2.0626480624709636e8;
ζ = pi*ρ*θ;
V = ( (1.0-param[2]-param[3])*besselj1.(ζ)./ζ +(param[2]+2*param[3])/sqrt(2/pi)*besselj.(1.5,ζ)./ζ.^(3/2)-2*param[3]*besselj.(2, ζ)./ζ.^2)/(.5-param[2]/6-param[3]/12)
indx= findall(abs.(ζ).<tol) #get rid of possible nan
V[findall(.!(isfinite.(V)))].=1.0;
return V
end

#power law
function visibility_ldpow(param, uv::Array{Float64,2}; maxk = 200)
ρ=sqrt.(uv[1,:].^2+uv[2,:].^2)
θ = param[1]/2.0626480624709636e8;
#f=k->gamma(param[2]/2+2)/(gamma(param[2]/2+k+2)*gamma(k+1))*(-0.25*(pi*θ*ρ).^2).^k; # note: can lead to NaN due to typical gamma blow up
#f=k->(-1)^k*exp.(2k*log.(0.5*pi*θ*ρ).+(loggamma(param[2]/2+2)-loggamma(param[2]/2+k+2)-loggamma(k+1)));
#V = sum(map(f,collect(0:maxk)));
V = gamma(param[2]/2+2)*besselj.(param[2]/2+1, pi*θ*ρ).*(0.5*pi*θ*ρ).^-(param[2]/2+1)
V[findall(.!(isfinite.(V)))].=1.0;
return V
end

#quadratic law
function visibility_ldquad(param,uv::Array{Float64,2};tol=1e-6)
ρ=sqrt.(uv[1,:].^2+uv[2,:].^2)
θ = param[1]/2.0626480624709636e8;
ζ = (pi*ρ*θ);
V = ( (1.0-param[2]-param[3])*besselj1.(ζ)./ζ +(param[2]+2*param[3])/sqrt(2/pi)*besselj.(1.5,ζ)./ζ.^1.5-2*param[3]*besselj.(2, ζ)./ζ.^2)/(.5-param[2]/6-param[3]/12)
indx = findall(abs.(ζ).<tol) #get rid of possible nan
V[findall(.!(isfinite.(V)))].=1.0;
return V
end

#square root law
function visibility_ldsquareroot(param,uv::Array{Float64,2};tol=1e-6)
ρ=sqrt.(uv[1,:].^2+uv[2,:].^2)
θ = param[1]/2.0626480624709636e8;
ζ = pi*ρ*θ;
V = ((1.0-param[2]-param[3])*besselj.(1, ζ)./ζ + param[2]*sqrt(pi/2)* besselj.(1.5, ζ)./ζ.^1.5 + param[3]*2^0.25*gamma(1.25)*besselj.(1.25, ζ)./ζ.^1.25)/(0.5-param[2]/6-param[3]/10)
V[findall(.!(isfinite.(V)))].=1.0;
return V
end

#linear law
function visibility_ldlin(param,uv::Array{Float64,2};tol=1e-6)
ρ=sqrt.(uv[1,:].^2+uv[2,:].^2)
θ = param[1]/2.0626480624709636e8;
ζ = (pi*ρ*θ);
V = ( (1.0-param[2])*besselj1.(ζ)./ζ +param[2]/sqrt(2/pi)*besselj.(1.5,ζ)./ζ.^(3/2))/(.5-param[2]/6)
indx = findall(abs.(ζ).<tol) #get rid of possible nan
V[findall(.!(isfinite.(V)))].=1.0;
return V
end

# visibility of an annulus of unit flux
function visibility_annulus(r_in, r_out, uv::Array{Float64,2};tol=1e-10)
ρ=sqrt.(uv[1,:].^2+uv[2,:].^2)
if abs(r_out - r_in)<tol #
  return zeros(Complex{Float64},length(ρ)); #shouldn't we use thin ring here ?
else
  return (1.0+0*im).*(visibility_ud([2*r_out],ρ)*r_out^2-visibility_ud([2*r_in],ρ)*r_in^2)/(r_out^2-r_in^2);
end
end

function visibility_thin_ring(param,uv::Array{Float64,2}) # i: inclination (deg), ϕ: semi-major axis orientation (deg)
# Parameters
#  1: radius of the ring (radius and not diameter)
#  2: ϕ = position angle
#  3: i = inclination
ϕ = param[2]/180*pi;
i = param[3]/180*pi;
θ = param[1]/2.0626480624709636e8;
uu = uv[1,:].*cos(ϕ) + uv[2,:].*sin(ϕ)
vv = (-uv[1,:].*sin(ϕ) + uv[2,:].*cos(ϕ))*cos(i)
ρ = sqrt.( uu.^2 + vv.^2)
return besselj0.(2*pi*θ*ρ)
end

function visibility_Gaussian_ring(param,uv::Array{Float64,2}) # i: inclination (deg), ϕ: semi-major axis orientation (deg)
# Parameters
#  1: radius of the ring
#  2: ϕ = position angle
#  3: i = inclination
#  4: ratio = half flux radius / ring radius
θ = param[1]/2.0626480624709636e8;
ϕ = param[2]/180*pi;
i = param[3]/180*pi;
uu = uv[1,:].*cos(ϕ) + uv[2,:].*sin(ϕ)
vv = (-uv[1,:].*sin(ϕ) + uv[2,:].*cos(ϕ))*cos(i)
ρ = sqrt.( uu.^2 + vv.^2)
return besselj0.(2*pi*θ*ρ).*exp.(-pi^2/log(2)*(θ*param[4])^2*(uv[1,:].^2+uv[2,:].^2))
end


function visibility_Lorentzian_ring(param,uv::Array{Float64,2}) # i: inclination (deg), ϕ: semi-major axis orientation (deg)
# Parameters
#  1: radius of the ring
#  2: ϕ = position angle
#  3: i = inclination
#  4: ratio = half flux radius / ring radius
# Note:
θ = param[1]/2.0626480624709636e8;
ϕ = param[2]/180*pi;
i = param[3]/180*pi;
uu = uv[1,:].*cos(ϕ) + uv[2,:].*sin(ϕ)
vv = (-uv[1,:].*sin(ϕ) + uv[2,:].*cos(ϕ))*cos(i)
ρ = sqrt.( uu.^2 + vv.^2)
return besselj0.(2*pi*θ*ρ).*exp.(-2*pi/sqrt(3)*θ*param[4]*sqrt.(uv[1,:].^2+uv[2,:].^2))
end

function visibility_Gaussian_ring_az(param,uv::Array{Float64,2}) # i: inclination (deg), ϕ: semi-major axis orientation (deg)
#  1: θ = radius of the ring
#  2: ϕ = position angle
#  3: i = inclination
#  4: ratio = half flux radius / ring radius
#  5,6,7,8: azimuthal coefficients
θ = param[1]/2.0626480624709636e8;
ϕ = param[2]/180*pi;
i = param[3]/180*pi;
uu = uv[1,:].*cos(ϕ) + uv[2,:].*sin(ϕ)
vv = (-uv[1,:].*sin(ϕ) + uv[2,:].*cos(ϕ))*cos(i)
ρ = sqrt.( uu.^2 + vv.^2)
#cα = uu./ρ
#sα = vv./ρ
#s2α = 2*cα.*sα = 2*uu.*vv./ρ.^2
#c2α = (uu.*uu - vv.*vv)/ρ^2
#az_modulation = - im*(param[5]*cos(α) + param[6]*sin(α))*besselj1.(2*pi*θ*ρ) + (param[7]*cos(2α)+param[8]*sin(2α))*besselj.(2,2*pi*θ*ρ)
az_modulation = - im*(param[5]*uu + param[6]*vv).*(besselj1.(2*pi*θ*ρ)./ρ) - (2*param[7]*uu.*vv+param[8]*(uu.*uu-vv.*vv)).*(besselj.(2,2*pi*θ*ρ)./ρ.^2)
az_modulation[findall(.!(isfinite.(az_modulation)))].=0.0; # prevents failure at ρ=0
return (besselj0.(2*pi*θ*ρ) + az_modulation).*exp.(-pi^2/log(2)*(θ*param[4])^2*(uv[1,:].^2+uv[2,:].^2))
end

function visibility_GaussianLorentzian_ring_az(param,uv::Array{Float64,2}) # i: inclination (deg), ϕ: semi-major axis orientation (deg)
#  1: θ = radius of the ring
#  2: ϕ = position angle
#  3: i = inclination
#  4: ratio = half flux radius / ring radius
#  5,6,7,8: azimuthal coefficients
#  9 : ratio Lorentzian/Gaussian
θ = param[1]/2.0626480624709636e8;
ϕ = param[2]/180*pi;
i = param[3]/180*pi;
uu = uv[1,:].*cos(ϕ) + uv[2,:].*sin(ϕ)
vv = (-uv[1,:].*sin(ϕ) + uv[2,:].*cos(ϕ))*cos(i)
ρ = sqrt.( uu.^2 + vv.^2)
ρρ = sqrt.(uv[1,:].^2+uv[2,:].^2) # original ρ, no incl
#cα = uu./ρ
#sα = vv./ρ
#s2α = 2*cα.*sα = 2*uu.*vv./ρ.^2
#c2α = (uu.*uu - vv.*vv)/ρ^2
#az_modulation = - im*(param[5]*cos(α) + param[6]*sin(α))*besselj1.(2*pi*θ*ρ) + (param[7]*cos(2α)+param[8]*sin(2α))*besselj.(2,2*pi*θ*ρ)
az_modulation = - im*(param[5]*uu + param[6]*vv).*(besselj1.(2*pi*θ*ρ)./ρ) - (2*param[7]*uu.*vv+param[8]*(uu.*uu-vv.*vv)).*(besselj.(2,2*pi*θ*ρ)./ρ.^2)
az_modulation[findall(.!(isfinite.(az_modulation)))].=0.0; # prevents failure at ρ=0
return (besselj0.(2*pi*θ*ρ) + az_modulation).*(param[9]*exp.(-pi^2/log(2)*(θ*param[4])^2*ρρ.^2)+(1-param[9])*exp.(-2*pi/sqrt(3)*θ*param[4]*ρρ))
end

#
# Overloaded functions of ρ   - TO REMOVE
#
# function visibility_ud(param, ρ::Array{Float64,1})
# t = param[1]/2.0626480624709636e8*pi*ρ;
# V = 2.0*besselj1.(t)./t
# V[findall(.!(isfinite.(V)))].=1.0;
# return V
# end
#
# function dvisibility_ud(param, ρ::Array{Float64,1})
# dt_dp = pi*ρ/2.0626480624709636e8
# t= param[1]*dt_dp
# dV_dt = (t.*besselj0.(t)-2*besselj1.(t))./t.^2
# return dV_dt.*dt_dp
# end
#
#
#
# function visibility_ldpow(param, ρ::Array{Float64,1}; maxk = 200)
# θ = param[1]/2.0626480624709636e8;
# #f=k->gamma(param[2]/2+2)/(gamma(param[2]/2+k+2)*gamma(k+1))*(-0.25*(pi*θ*ρ).^2).^k; # note: can lead to NaN due to typical gamma blow up
# f=k->(-1)^k*exp.(2k*log.(0.5*pi*θ*ρ).+(logabsgamma(param[2]/2+2)[1]-logabsgamma(param[2]/2+k+2)[1]-logabsgamma(k+1)[1]));
# V = sum(map(f,collect(0:maxk)));
# return V
# end
#
# #quadratic law
# function visibility_ldquad(param,ρ::Array{Float64,1})
# θ = param[1]/2.0626480624709636e8;
# ζ = (pi*ρ*θ);
# V = ((1.0-param[2]-param[3])*(besselj1.(ζ)./ζ)+((θ+2*param[3])/sqrt(2/pi))*(
# (sqrt.(2 ./(pi*ζ)).*((sin.(ζ)./ζ)-cos.(ζ)))./ζ.^(3/2))-2*param[3]*
# (besselj.(2, ζ)./ζ.^2))./(0.5-param[2]/6-param[3]/12)
# return V
# end
#
# #linear law
# function visibility_ldlin(param,ρ::Array{Float64,1})
# θ = param[1]/2.0626480624709636e8;
# ζ = (pi*ρ*θ);
# V = ((1.0-param[2])*(besselj1.(ζ)./ζ)+θ/sqrt(2/pi)*(
# (sqrt.(2 ./(pi*ζ)).*((sin.(ζ)./ζ)-cos.(ζ)))./ζ.^(3/2)))./(0.5-param[2]/6)
# return V
# end
#
# # visibility of an annulus of unit flux
# function visibility_annulus(r_in, r_out, ρ::Array{Float64,1};tol=1e-10)
#     if abs(r_out - r_in)<tol #
#         return zeros(Complex{Float64},length(ρ)); #shouldn't we use thin ring here ?
#     else
#         return (1.0+0*im).*(visibility_ud([2*r_out],ρ)*r_out^2-visibility_ud([2*r_in],ρ)*r_in^2)/(r_out^2-r_in^2);
#     end
# end
#
# function visibility_thinring(param,ρ::Array{Float64,1})
# return besselj0.(pi*param[1]/2.0626480624709636e8*ρ)
# end

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


if visfunc == visibility_ldsquareroot
 lbounds = [0.0, -1.0, -1.0]
 hbounds = [1e9, 1.0, 1.0]
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
