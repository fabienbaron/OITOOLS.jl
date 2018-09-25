#
# This is a helper file to fit SATLAS profiles to interferometric data
#
#
using DelimitedFiles
using SpecialFunctions

struct SATLAS_MODEL
  flux::Array{Float64,1}
  mu::Array{Float64,1}
end

function read_satlas(file)
sat1 = readdlm(file)
f_satlas = sat1[:,6];
push!(f_satlas,1.0);
mu_satlas = sat1[:,1];
push!(mu_satlas,1.0);
return SATLAS_MODEL(f_satlas,mu_satlas)
end

using PyCall
sp=pyimport("scipy.special")
#hankel = sp[:struve](order, x)


function IA(x)
    return (x.^2).*besselj1.(x)-(pi*x./2).*(besselj1.(x).*sp[:struve](0, x).-besselj0.(x).*sp[:struve](1, x))
end

function visibility_satlas_any(params::Array{Float64,1}, v_r::Array{Float64,1}, satmodel::SATLAS_MODEL)
# v_r is the baseline radius
# diameter of star in mas
# SATLAS model from read_satlas()
r = reverse(0.5*params[1]*(1.0 .- satmodel.mu.^2));
flux = reverse(satmodel.flux);
#a = (satmodel.flux[2:end]-satmodel.flux[1:end-1])./(r[2:end]- r[1:end-1]);
#b = satmodel.flux[1:end-1]-a.*r[1:end-1];
b = 0.5*(flux[1:end-1] + flux[2:end]).*(r[2:end].^2-r[1:end-1].^2)*pi
nuv = length(v_r);
vis = zeros(Complex{Float64},nuv);
totflux=sum(b);
for i=1:length(r)-1
#  vis += a[i]*(  IA(2*r[i+1]*v_r) - IA(2*r[i]*v_r) );
  vis += b[i]*visibility_annulus(r[i],r[i+1],v_r);
end
vis /= totflux;
return vis
end

function visibility_satlas_img(params::Array{Float64,1}, v_r::Array{Float64,1},satmodel::SATLAS_MODEL; npix=256, pixsize=0.05)
  # Note: NFFT could be used here, of course
x2=(repeat(collect(1:npix),1,npix).-(npix+1)/2).^2;
r=sqrt.(x2+x2');
rstar = 0.5*params[1]/pixsize # rstar in pixels
mask=findall(r.>rstar);
r[mask].=rstar;
#new mu values to interpolate:
mu_star = sqrt.(1.0 .- (r/rstar).^2);
indx_lo = vec(Int.(floor.(mu_star*1000 .+1)))
indx_hi = vec(Int.(floor.(mu_star*1000).+2)); indx_hi[end]=1001;
mu_lo = satmodel.mu[indx_lo]
mu_hi = satmodel.mu[indx_hi]
f_lo = satmodel.flux[indx_lo]
f_hi = satmodel.flux[indx_hi]
f = f_lo+(vec(mu_star)-mu_lo).*((f_hi-f_lo)./(mu_hi-mu_lo))
flux = sum(f);
nuv = length(v_r);
dft = zeros(Complex{Float64}, nuv, npix*npix);
xvals = -2 * pi * im * pixsize * (pi / 180.0) / 3600000.0*[(mod(i-1,npix)+1) for i=1:npix*npix];
for uu=1:nuv
    dft[uu,:] = exp.(v_r[uu] * xvals);
end
return dft * f / flux;
end

using NLopt
function fit_satlas_v2(data::OIdata, visfunc, init_param::Array{Float64,1}, satmodel::SATLAS_MODEL)
  indx= data.indx_v2;
  nv2 = length(data.v2[indx]);
  chisq=(param,g)->norm((abs2.(visfunc(param,data.v2_baseline,satmodel))-data.v2[indx])./data.v2_err[indx])^2/nv2;
  opt = Opt(:LN_NELDERMEAD, 1);
  min_objective!(opt, chisq);
  xtol_rel!(opt,1e-5);
  bounds = [0e0];
  lower_bounds!(opt, bounds);
  (minf,minx,ret) = optimize(opt, init_param);
  println("Chi2: $minf \t parameters:$minx \t \t $ret")
  cvis_model = visfunc(minx,sqrt.(data.uv[1,:].^2+data.uv[2,:].^2),satmodel)
  return (minf,minx,cvis_model)
end
