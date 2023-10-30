#
# Classic visibility functions (uniform disc, etc.)
#
using SpecialFunctions

#param  : model parameter
#param[1] = diameter in mas
#ρ = radius in uv space
function visibility_ud(param, uv::Array{Float64,2})
ρ = sqrt.(uv[1,:].^2+uv[2,:].^2)
V = jinc.(param[1]/2.0626480624709636e8*ρ)
return V
end

function dvisibility_ud(param, uv::Array{Float64,2})
ρ=sqrt.(uv[1,:].^2+uv[2,:].^2)
dt_dp = pi*ρ/2.0626480624709636e8
t = param[1]*dt_dp .+ eps(Float64) # need eps() to remove singularity in next line
dV_dt = 2.0*(t.*besselj0.(t)-2*besselj1.(t))./t.^2   
return dV_dt.*dt_dp
end

function visibility_ellipse_uniform(param, uv::Array{Float64,2}) # ϕ: semi-major axis orientation (deg)
ϵ = param[2];
ϕ = param[3]/180*pi;
ρ = sqrt.(ϵ^2*(uv[1,:].*cos.(ϕ)-uv[2,:].*sin.(ϕ) ).^2+(uv[2,:].*cos.(ϕ)+uv[1,:].*sin.(ϕ)).^2)
t = param[1]/2.0626480624709636e8*ρ;
V = jinc.(t)
return V
end

function visibility_sigmoid(param, uv::Array{Float64,2})
ρ=sqrt.(uv[1,:].^2+uv[2,:].^2)
V = -im*pi*csch.(pi*param[1]/2.0626480624709636e8*ρ)
return V
end



function visibility_ellipse_quad(param, uv::Array{Float64,2};tol=1e-6) #ϵ: ellipticity  i: inclination (deg), ϕ: position angle (deg)
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
function visibility_ldpow(param, uv::Array{Float64,2})
ρ=sqrt.(uv[1,:].^2+uv[2,:].^2)
θ = param[1]/2.0626480624709636e8;
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

#quadratic law, alternative formulation with "triangular sampling"
#  Kipping et al. 2013, https://arxiv.org/abs/1308.0009# eq 15 and 16
function visibility_ldquad_tri(param,uv::Array{Float64,2};tol=1e-6)
ρ=sqrt.(uv[1,:].^2+uv[2,:].^2)
θ = param[1]/2.0626480624709636e8;
u1 = 2*sqrt(param[2])*param[3]
u2 = sqrt(param[2])*(1-2*param[3])
ζ = (pi*ρ*θ);
V = ( (1.0-u1-u2)*besselj1.(ζ)./ζ +(u1+2*u2)/sqrt(2/pi)*besselj.(1.5,ζ)./ζ.^1.5-2*u2*besselj.(2, ζ)./ζ.^2)/(.5-u1/6-u2/12)
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
function visibility_annulus(param, uv::Array{Float64,2};tol=1e-10) # bugged
if param[2] - param[1]<tol #
  return zeros(Complex{Float64},size(uv,2)); #shouldn't we use thin ring here ?
else
  return (visibility_ud([2*param[2]],uv)*param[2]^2-visibility_ud([2*param[1]],uv)*param[1]^2)/(param[2]^2-param[1]^2);
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



function visibility_thin_ring_az(param,uv::Array{Float64,2}) # i: inclination (deg), ϕ: semi-major axis orientation (deg)
#  1: θ = radius of the ring
#  2: ϕ = position angle
#  3: i = inclination
#  4,5,6,7: azimuthal coefficients
θ = param[1]/2.0626480624709636e8;
ϕ = param[2]/180*pi;
i = param[3]/180*pi;
u = @view uv[1,:]
v = @view uv[2,:]
uu = u.*cos(ϕ) + v.*sin(ϕ)
vv = (-u.*sin(ϕ) + v.*cos(ϕ))*cos(i)
ρ = sqrt.( uu.^2 + vv.^2)
az_modulation = - im*(param[4]*uu + param[5]*vv).*(besselj1.(2*pi*θ*ρ)./ρ) - (2*param[6]*uu.*vv+param[7]*(uu.*uu-vv.*vv)).*(besselj.(2,2*pi*θ*ρ)./ρ.^2)
az_modulation[findall(.!(isfinite.(az_modulation)))].=0.0; # prevents failure at ρ=0
return (besselj0.(2*pi*θ*ρ) + az_modulation)
end


function visibility_Gaussian_ring(param,uv::Array{Float64,2}) # i: inclination (deg), ϕ: semi-major axis orientation (deg)
# Parameters
#  1: radius of the ring
#  2: ϕ = position angle
#  3: i = inclination
#  4: w = ring FWHM
θ = param[1]/2.0626480624709636e8;
ϕ = param[2]/180*pi;
i = param[3]/180*pi;
w = param[4]/2.0626480624709636e8;
uu = uv[1,:].*cos(ϕ) + uv[2,:].*sin(ϕ)
vv = (-uv[1,:].*sin(ϕ) + uv[2,:].*cos(ϕ))*cos(i)
ρ = sqrt.( uu.^2 + vv.^2)
return besselj0.(2*pi*θ*ρ).*exp.(-pi^2/log(2)*w^2*(uv[1,:].^2+uv[2,:].^2))
end



function visibility_Gaussian(param,uv::Array{Float64,2}) # i: inclination (deg), ϕ: semi-major axis orientation (deg)
# Parameters
#  1: θ=FWHM
#  2: i: inclination (used for aspect ratio too)
#  3: ϕ = position angle
θ = param[1]/2.0626480624709636e8;
i = param[2]/180*pi;
ϕ = param[3]/180*pi;
uu = uv[1,:].*cos(ϕ) + uv[2,:].*sin(ϕ)
vv = (-uv[1,:].*sin(ϕ) + uv[2,:].*cos(ϕ))*cos(i)
return exp.(-pi^2/log(2)*θ^2*(uu.^2+vv.^2))
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
w = param[4]/2.0626480624709636e8;
uu = uv[1,:].*cos(ϕ) + uv[2,:].*sin(ϕ)
vv = (-uv[1,:].*sin(ϕ) + uv[2,:].*cos(ϕ))*cos(i)
ρ = sqrt.( uu.^2 + vv.^2)
return besselj0.(2*pi*θ*ρ).*exp.(-2*pi/sqrt(3)*w*sqrt.(uv[1,:].^2+uv[2,:].^2))
end

function visibility_Gaussian_ring_az(param,uv::Array{Float64,2}) # i: inclination (deg), ϕ: semi-major axis orientation (deg)
#  1: θ = radius of the ring
#  2: ϕ = position angle
#  3: i = inclination
#  4: w = ring FWHM
#  5,6,7,8: azimuthal coefficients
θ = param[1]/2.0626480624709636e8;
ϕ = param[2]/180*pi;
i = param[3]/180*pi;
w = param[4]/2.0626480624709636e8;
uu = uv[1,:].*cos(ϕ) + uv[2,:].*sin(ϕ)
vv = (-uv[1,:].*sin(ϕ) + uv[2,:].*cos(ϕ))*cos(i)
ρ = sqrt.( uu.^2 + vv.^2)
az_modulation = - im*(param[5]*uu + param[6]*vv).*(besselj1.(2*pi*θ*ρ)./ρ) - (2*param[7]*uu.*vv+param[8]*(uu.*uu-vv.*vv)).*(besselj.(2,2*pi*θ*ρ)./ρ.^2)
az_modulation[findall(.!(isfinite.(az_modulation)))].=0.0; # prevents failure at ρ=0
return (besselj0.(2*pi*θ*ρ) + az_modulation).*exp.(-pi^2/log(2)*w^2*(uv[1,:].^2+uv[2,:].^2))
end

function visibility_GaussianLorentzian_ring_az(param,uv::Array{Float64,2}) # i: inclination (deg), ϕ: semi-major axis orientation (deg)
#  1: θ = radius of the ring
#  2: ϕ = position angle
#  3: i = inclination
#  4: w = ring FWHM
#  5,6,7,8: azimuthal coefficients
#  9 : ratio Lorentzian/Gaussian
θ = param[1]/2.0626480624709636e8;
ϕ = param[2]/180*pi;
i = param[3]/180*pi;
w = param[4]/2.0626480624709636e8;
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
return (besselj0.(2*pi*θ*ρ) + az_modulation).*(param[9]*exp.(-pi^2/log(2)*w^2*ρρ.^2)+(1-param[9])*exp.(-2*pi/sqrt(3)*w*ρρ))
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
 lbounds = [0.0, 0.0, -1.0]
 hbounds = [1e9, 2.0, 1.0]
end

if visfunc == visibility_ldquad_alt
 lbounds = [0.0, 0.0, 0.0]
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
