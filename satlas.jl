using DataFrames
using SpecialFunctions

type SATLAS_MODEL
  flux::Array{Float64,1}
  mu::Array{Float64,1}
end

function read_satlas(file)
sat1 = DataFrames.readtable(file, header = false, separator = ' ', nastrings=["NA", "na", "n/a", "missing"]);
names!(sat1, [Symbol("col$i") for i in 1:9]);
rename!(sat1, [:col1, :col6], [:mu, :flux]);
f_satlas = convert(Array,sat1[:6]);
push!(f_satlas,1.0);
mu_satlas = convert(Array,sat1[:1]);
push!(mu_satlas,1.0);
return SATLAS_MODEL(f_satlas,mu_satlas)
end

using PyCall
sp=pyimport("scipy.special")
#hankel = sp[:struve](order, x)


function IA(x)
    return (x.^2).*besselj1.(x)-(pi*x./2).*(besselj1.(x).*sp[:struve](0, x).-besselj0.(x).*sp[:struve](1, x))
end

function complex_vis_satlas_any(v_r::Array{Float64,1}, diameter::Real, satmodel::SATLAS_MODEL)
# v_r is the baseline radius
# diameter of star in mas
# SATLAS model from read_satlas()
r = reverse(0.5*diameter*(1-satmodel.mu.^2));
satmodel.flux = reverse(satmodel.flux);
a = (satmodel.flux[2:end]-satmodel.flux[1:end-1])./(r[2:end]- r[1:end-1]);
#b = satmodel.flux[1:end-1]-a.*r[1:end-1];
b = 0.5*(satmodel.flux[1:end-1] + satmodel.flux[2:end]).*(r[2:end].^2-r[1:end-1].^2)*pi
nuv = length(v_r);
vis = zeros(Complex64,nuv);

totflux=sum(b);
for i=1:length(r)-1
#  vis += a[i]*(  IA(2*r[i+1]*v_r) - IA(2*r[i]*v_r) );
  vis += b[i]*visibility_annulus(r[i],r[i+1],v_r);
end
vis /= totflux;
return vis
end


include("oichi2.jl")
function complex_vis_satlas_img(v_r::Array{Float64,1}, diameter::Real, satmodel::SATLAS_MODEL; npix=256, pixsize=0.05)
x2=(repmat(collect(1:npix),1,npix)-(npix+1)/2).^2;
r=sqrt.(x2+x2');
rstar = 0.5*diameter/pixsize # rstar in pixels
mask=find(r.>rstar);
r[mask]=rstar;
#new mu values to interpolate:
mu_star = sqrt.(1-(r/rstar).^2);
indx_lo = vec(Int.(floor.(mu_star*1000+1)))
indx_hi = vec(Int.(floor.(mu_star*1000)+2)); indx_hi[end]=1001;
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
