include("../OITOOLS.jl/oitools.jl")
using NFFT
function cvis_to_t3(cvis, indx1, indx2, indx3)
  t3 = cvis[indx1].*cvis[indx2].*conj(cvis[indx3]);
  t3amp = abs.(t3);
  t3phi = angle.(t3)*180./pi;
  return t3, t3amp, t3phi
end

δ=+18.594234/180*pi #ce tau 10 mas 32w 0.37 pix

#include("chara_config.jl")
include("npoi_config.jl")
hour_angles = linspace(-3.75,3.75,13);

staxyz=Array{Float64}(3,N)
for i=1:N
        staxyz[:,i]=station_xyz[i][:];
end

# V2 Baselines
nv2 = Int64(N*(N-1)/2);
v2_baselines = Array{Float64}(3,nv2);
v2_stations  = Array{Int64}(2,nv2);
v2_stations_nonredun=Array{Int64}(2,nv2);
v2_indx      = Array{Int64}(nv2);
baseline_list = Array{String}(nv2);
indx = 1
for i=1:N
  for j=i+1:N
    #if(j!=i)
        v2_baselines[:,indx] = station_xyz[j]-station_xyz[i];
        v2_stations[:,indx] = [i,j];
        v2_indx[indx] = indx;
        indx += 1
#    end
  end
end

# T3 Baselines
nt3 = binomial(N,3);# V2 Baselines
function v2mapt3(i,j)
    # really dumb way of doing this
    # there is probably an analytic formula for this
    listi = find(v2_stations[1,:].==i);
    return listi[findfirst(v2_stations[2,listi].==j)]
end

t3_baselines=Array{Float64}(3,3,nt3)
t3_stations  = Array{Int64}(3, nt3)
t3_indx_1  = Array{Int64}(nt3)
t3_indx_2  = Array{Int64}(nt3)
t3_indx_3  = Array{Int64}(nt3)
indx = 1
for i=1:N
  for j=i+1:N-1
    for k=j+1:N
        t3_baselines[:, 1, indx] = station_xyz[j]-station_xyz[i];
        t3_baselines[:, 2, indx] = station_xyz[k]-station_xyz[j];
        t3_baselines[:, 3, indx] = station_xyz[i]-station_xyz[k];
        t3_stations[:,indx]=[i,j,k];
        t3_indx_1[indx] = v2mapt3(i,j);
        t3_indx_2[indx] = v2mapt3(j,k);
        t3_indx_3[indx] = v2mapt3(i,k);
        indx += 1
    end
  end
end

# Merge all uv

#Use following expression only if there are missing baselines somewhere
#vis_baselines = hcat(v2_baselines, t3_baselines[:, 1, :], t3_baselines[:, 2, :], t3_baselines[:, 3, :])

#Expression to use for pure simulation where the full complement of v2 and t3 are created
vis_baselines = deepcopy(v2_baselines)
nuv = size(vis_baselines, 2);

# Baselines to proj Baselines (uv m)
# Now compute the UV coordinates, again according to the APIS++ standards.
#δ=+46.4668383/180*pi
nhours = length(hour_angles);
l=lat/180*pi
h = hour_angles' * pi / 12;
u_M = -sin(l)*sin.(h) .* vis_baselines[1,:] .+ cos.(h) .* vis_baselines[2,:]+cos(l)*sin.(h).*vis_baselines[3,:];
v_M = (sin(l)*sin(δ)*cos.(h).+cos(l)*cos(δ)) .* vis_baselines[1,:] + sin(δ)*sin.(h) .* vis_baselines[2,:]+(-cos(l)*cos.(h)*sin(δ).+sin(l)*cos(δ)) .* vis_baselines[3,:];
w_M =  (-sin(l)*cos(δ)* cos.(h).+cos(l)*sin(δ)).* vis_baselines[1,:] - cos(δ)* sin.(h) .* vis_baselines[2,:]  .+ (cos(l)*cos.(h)*cos(δ).+sin(l)sin(δ)) .* vis_baselines[3,:];

u_M = -u_M
v_M = -v_M

v2_indx_M = repmat(v2_indx,1,nhours)+repmat(Int64.(collect(linspace(0,nuv*(nhours-1),nhours)))',nv2)
t3_indx_1_M = repmat(t3_indx_1,1,nhours)+repmat(Int64.(collect(linspace(0,nuv*(nhours-1),nhours)))',nt3)
t3_indx_2_M = repmat(t3_indx_2,1,nhours)+repmat(Int64.(collect(linspace(0,nuv*(nhours-1),nhours)))',nt3)
t3_indx_3_M = repmat(t3_indx_3,1,nhours)+repmat(Int64.(collect(linspace(0,nuv*(nhours-1),nhours)))',nt3)

# proj baselines to (uv wav)
nw=length(λ)
u = reshape((1./λ)'.*vec(u_M), (nuv,nhours,length(λ)));
v = reshape((1./λ)'.*vec(v_M), (nuv,nhours,length(λ)));
w = reshape((1./λ)'.*vec(w_M), (nuv,nhours,length(λ)));
uv = vcat(vec(u)',vec(v)')

v2_indx_w = vec(repmat(vec(v2_indx_M),1,nw)+repmat(Int64.(collect(linspace(0,nuv*nhours*(nw-1),nw)))',nv2*nhours))
t3_indx_1_w = vec(repmat(vec(t3_indx_1_M),1,nw)+repmat(Int64.(collect(linspace(0,nuv*nhours*(nw-1),nw)))',nt3*nhours))
t3_indx_2_w = vec(repmat(vec(t3_indx_2_M),1,nw)+repmat(Int64.(collect(linspace(0,nuv*nhours*(nw-1),nw)))',nt3*nhours))
t3_indx_3_w = vec(repmat(vec(t3_indx_3_M),1,nw)+repmat(Int64.(collect(linspace(0,nuv*nhours*(nw-1),nw)))',nt3*nhours))
t3_baseline = (sqrt.(uv[1,t3_indx_1_w].^2 + uv[2,t3_indx_1_w].^2).*
    sqrt.(uv[1,t3_indx_2_w].^2 + uv[2,t3_indx_2_w].^2).*
    sqrt.(uv[1,t3_indx_3_w].^2 + uv[2,t3_indx_3_w].^2)).^(1./3.);


fitsfile = "2004true137.fits";
pixsize=0.1
x = (read((FITS(fitsfile))[1])); x=x[:,end:-1:1]; nx = (size(x))[1]; x=vec(x)/sum(x);


#nfft_plan = setup_nfft(-uv, nx, pixsize)
#cvis_model = image_to_cvis_nfft(x, nfft_plan)

dft = setup_dft(-uv, nx, pixsize);
cvis_model = image_to_cvis_dft(x, dft);
v2_model = cvis_to_v2(cvis_model, v2_indx_w);
t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, t3_indx_1_w, t3_indx_2_w, t3_indx_3_w);

#Add noise

v2_model_err = 0.00333*v2_model+1e-5
v2_model += v2_model_err.*randn(length(v2_model))

t3amp_model_err =0.007*t3amp_model+1e-6
t3amp_model += t3amp_model_err.*randn(length(t3amp_model))

t3phi_model_err = zeros(length(t3phi_model))+2. # degree  -- there is another way of setting this with Haniff formula
t3phi_model += t3phi_model_err.*randn(length(t3phi_model))

telnames=tel_names
sta_names=telnames
sta_index=Int64.(collect(linspace(1,N,N)))

target_id=1;
target="SuperStar!"; #get from user
raep0=314.49769042;#get from user
decep0=46.46683833;#get from user
equinox=2000;
ra_err=-1;
dec_err=-1;
sysvel=-1;
veltyp="LSR";
veldef="OPTICAL";
pmra=-1;
pmdec=-1;
pmra_err=-1;
pmdec_err=-1;
parallax=-1;
para_err=-1;
spectyp="Unknown";

eff_wave= λ
eff_band = δλ

target_id_vis2=ones(nv2*nhours)
time_vis2=ones(nv2*nhours) #change
mjd_vis2=ones(nv2*nhours)*58297.  #change
int_time_vis2=ones(nv2*nhours) #exchange
flag_vis2=fill(false,nv2*nhours,1)
#need to get vis2,vis2err,u,v,sta_index from DATA
target_id_t3=ones(nt3*nhours)
time_t3=ones(nt3*nhours) #change
mjd_t3=ones(nt3*nhours)*58297.  #change
int_time_t3=ones(nt3*nhours) #exchange;
flag_t3=fill(false,(nt3*nhours),1);


ucoord_vis2 = u_M[v2_indx_M]
vcoord_vis2 = v_M[v2_indx_M]
u1coord = u_M[t3_indx_1_M]
v1coord = v_M[t3_indx_1_M]
u2coord = u_M[t3_indx_2_M]
v2coord = v_M[t3_indx_2_M]

v2_model=reshape(reshape(v2_model,(nhours,nv2,nw)),(nhours*nv2,nw));
v2_model_err=reshape(reshape(v2_model_err,(nhours,nv2,nw)),(nhours*nv2,nw));

t3amp_model=reshape(reshape(t3amp_model,(nhours,nt3,nw)),(nhours*nt3,nw));
t3amp_model_err=reshape(reshape(t3amp_model_err,(nhours,nt3,nw)),(nhours*nt3,nw));

t3phi_model=reshape(reshape(t3phi_model,(nhours,nt3,nw)),(nhours*nt3,nw));
t3phi_model_err=reshape(reshape(t3phi_model_err,(nhours,nt3,nw)),(nhours*nt3,nw));

v2_model_stations=repmat(v2_stations,1,nhours);
t3_model_stations=repmat(t3_stations,1,nhours);
v2_model=transpose(v2_model);
v2_model_err=transpose(v2_model_err);
t3amp_model=transpose(t3amp_model);
t3amp_model_err=transpose(t3amp_model_err);
t3phi_model=transpose(t3phi_model);
t3phi_model_err=transpose(t3phi_model_err);

oi_array=[telnames,sta_names,sta_index,tel_diams,staxyz];
target_array=[target_id,target,raep0,decep0,equinox,ra_err,dec_err,sysvel,veltyp,veldef,pmra,pmdec,pmra_err,pmdec_err,parallax,para_err,spectyp];
wave_array=[eff_wave,eff_band];
vis2_array=[target_id_vis2,time_vis2,mjd_vis2,int_time_vis2,v2_model,v2_model_err,vec(ucoord_vis2),vec(vcoord_vis2),v2_model_stations,flag_vis2];
t3_array=[target_id_t3,time_t3,mjd_t3,int_time_t3,t3amp_model,t3amp_model_err,t3phi_model,t3phi_model_err,vec(u1coord),vec(v1coord),vec(u2coord),vec(v2coord),t3_model_stations,flag_t3];


outfilename ="!2004fake137.oifits"
f = fits_create_file(outfilename);
write_oi_header(f,1);
write_oi_array(f,oi_array);
write_oi_target(f,target_array);
write_oi_wavelength(f,wave_array);
write_oi_vis2(f,vis2_array);
write_oi_t3(f,t3_array);
fits_close_file(f);

#
# Double check
#
function cvis_to_t3(cvis, indx1, indx2, indx3)
  t3 = cvis[indx1].*cvis[indx2].*cvis[indx3];
  t3amp = abs.(t3);
  t3phi = angle.(t3)*180./pi;
  return t3, t3amp, t3phi
end
oifitsfile = "2004fake137.oifits"
data = readoifits(oifitsfile)[1,1];

fitsfile = "2004true137.fits";
pixsize = 0.1
x_true = (read((FITS(fitsfile))[1])); nx = (size(x_true))[1]; x_true=vec(x_true);
dft = setup_dft(data.uv, nx, pixsize);
f_chi2 = chi2(x_true, dft, data);


#
# Imaging
#
using OptimPack
pixsize = 0.2
nx = 64
dft = setup_dft(data.uv, nx, pixsize);

x_start = Array{Float64}(nx, nx);
     for i=1:nx
       for j=1:nx
         x_start[i,j] = exp(-((i-(nx+1)/2)^2+(j-(nx+1)/2)^2)/(2*(nx/6)^2));
       end
     end
 x_start = vec(x_start)/sum(x_start);
#crit = (x,g)->chi2_centered_fg(x, g, dft, data);
crit = (x,g)->chi2_centered_l2_fg(x, g, 1e0, dft, data );
crit = (x,g)->chi2_TVSQ_fg(x, g, 1e0, dft, data );

x_sol = OptimPack.vmlmb(crit, x_start, verb=true, lower=0, maxiter=80, blmvm=false);
