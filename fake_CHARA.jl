include("oitools.jl")
include("write_oifits_ha.jl")
PyDict(pyimport("matplotlib")["rcParams"])["font.family"]=["serif"]
#PyDict(pyimport("matplotlib")["rcParams"])["mathtext.fontset"]=["custom"]
PyDict(pyimport("matplotlib")["rcParams"])["xtick.major.size"]=[6]
PyDict(pyimport("matplotlib")["rcParams"])["ytick.major.size"]=[6]
PyDict(pyimport("matplotlib")["rcParams"])["xtick.minor.size"]=[6]
PyDict(pyimport("matplotlib")["rcParams"])["ytick.minor.size"]=[6]
PyDict(pyimport("matplotlib")["rcParams"])["xtick.major.width"]=[1]
PyDict(pyimport("matplotlib")["rcParams"])["ytick.major.width"]=[1]
PyDict(pyimport("matplotlib")["rcParams"])["xtick.minor.width"]=[1]
PyDict(pyimport("matplotlib")["rcParams"])["ytick.minor.width"]=[1]
PyDict(pyimport("matplotlib")["rcParams"])["lines.markeredgewidth"]=[1]
PyDict(pyimport("matplotlib")["rcParams"])["legend.numpoints"]=[1]
#PyDict(pyimport("matplotlib")["rcParams"])["legend.frameon"]=["False"]
PyDict(pyimport("matplotlib")["rcParams"])["legend.handletextpad"]=[0.3]
using NFFT

function cvis_to_t3(cvis, indx1, indx2, indx3)
  t3 = cvis[indx1].*cvis[indx2].*conj(cvis[indx3]);
  t3amp = abs.(t3);
  t3phi = angle.(t3)*180./pi;
  return t3, t3amp, t3phi
end

δ=+18.594234/180*pi #ce tau 10 mas 32w 0.37 pix
#δ=δ/180*pi # mu cep 20 mas 0.37 pixsize
hour_angles=collect(linspace(-6,6,10)) #[-]'
#hour_angles=collect(linspace(-3,3,6)) #[-]'
name = "CHARA"
lat = 34.2243836944
lon = -118.0570313111
alt = 1731.264

N = 6
tel_diams = ones(Float64, N)
#tel_names =["S1","S3","S2","E1","E2","W3","W1","W2"]
tel_names =["S1","S2","E1","E2","W1","W2"]
#extra w
station_xyz = [[0  ,  0 ,  0 ],
[-5.744414576 , 33.585327515 ,  0.634532856],
[ 125.3331333, 305.9284973  ,  -5.9190997],
[70.3891451 , 269.7146871  ,  -2.8025644],
[-175.0684101 , 216.3272464  ,  -10.7975261],
[-69.0845925, 199.3424346, 0.4706086]]

staxyz=Array{Float64}(3,N)
for i=1:N
        staxyz[:,i]=station_xyz[i][:]
end


# V2 Baselines
nv2 = Int64(N*(N-1)/2);
v2_baselines = Array{Float64}(3,nv2)
v2_stations  = Array{Int64}(2,nv2)
v2_stations_nonredun=Array{Int64}(2,nv2)
v2_indx      = Array{Int64}(nv2)
baseline_list = Array{String}(nv2)
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
    listi = find(v2_stations[1,:].==i)
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


#hour_angles = linspace(-5,-1.8,)'

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
# add spectral info
#λ = collect(linspace(6.4E-7,9.4E-7,16))
λ = [1.745909E-06, 1.717131E-06, 1.684602E-06, 1.650831E-06 ,1.615819E-06,  1.579565E-06,1.542070E-06 , 1.508262E-06];
#λ = [1.745909E-06, 1.717131E-06]
nw=length(λ)
u = reshape((1./λ)'.*vec(u_M), (nuv,nhours,length(λ)));
v = reshape((1./λ)'.*vec(v_M), (nuv,nhours,length(λ)));
w = reshape((1./λ)'.*vec(w_M), (nuv,nhours,length(λ)));
uv = vcat(vec(u)',vec(v)')

v2_indx_w = vec(repmat(vec(v2_indx_M),1,nw)+repmat(Int64.(collect(linspace(0,nuv*nhours*(nw-1),nw)))',nv2*nhours))
t3_indx_1_w = vec(repmat(vec(t3_indx_1_M),1,nw)+repmat(Int64.(collect(linspace(0,nuv*nhours*(nw-1),nw)))',nt3*nhours))
t3_indx_2_w = vec(repmat(vec(t3_indx_2_M),1,nw)+repmat(Int64.(collect(linspace(0,nuv*nhours*(nw-1),nw)))',nt3*nhours))
t3_indx_3_w = vec(repmat(vec(t3_indx_3_M),1,nw)+repmat(Int64.(collect(linspace(0,nuv*nhours*(nw-1),nw)))',nt3*nhours))



# Load file
# Monochromatic monotemporal image

fitsfile = "2004true137.fits";
pixsize=0.1
#fitsfile = "oifits-sim-test_64.fits"

x = (read((FITS(fitsfile))[1])); x=x[:,end:-1:1]; nx = (size(x))[1]; x=vec(x)/sum(x);

nfft_plan = setup_nfft(-uv, nx, pixsize)
cvis_model = image_to_cvis_nfft(x, nfft_plan)
v2_model = cvis_to_v2(cvis_model, v2_indx_w);
t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, t3_indx_1_w, t3_indx_2_w, t3_indx_3_w);


t3_baseline = (sqrt.(uv[1,t3_indx_1_w].^2 + uv[2,t3_indx_1_w].^2).*
    sqrt.(uv[1,t3_indx_2_w].^2 + uv[2,t3_indx_2_w].^2).*
    sqrt.(uv[1,t3_indx_3_w].^2 + uv[2,t3_indx_3_w].^2)).^(1./3.);

#Add noise

v2_model_err = 0.01*v2_model+1e-5
v2_model += v2_model_err.*randn(length(v2_model))

t3amp_model_err =0.05*t3amp_model+1e-6
t3amp_model += t3amp_model_err.*randn(length(t3amp_model))

t3phi_model_err = zeros(length(t3phi_model))+0.5 # degree  -- there is another way of setting this with Haniff formula
t3phi_model += t3phi_model_err.*randn(length(t3phi_model))
#=
v2_model_err = 1e-8

t3amp_model_err = abs(0.1*t3amp_model)


t3phi_model_err = 1. # degree  -- there is another way of setting this with Haniff formula
=#

#=
 v2_baselines=sqrt.(uv[1,:].^2+uv[2,:].^2)
 v2_baselines_new=reshape(v2_baselines,(nv2,nhours,nw))

v2_model=reshape(v2_model,(nv2,nhours,nw))
=#

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

eff_wave= λ #[1.1e-6,1.2e-6,1.3e-6]#, 1.684602E-06, 1.650831E-06, 1.615819E-06, 1.579565E-06, 1.542070E-06, 1.508262E-06]
eff_band=[2.877857E-08, 3.252876E-08, 3.377056E-08, 3.501202E-08, 3.625382E-08, 3.749551E-08, 3.749551E-08, 3.380819E-08]

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

v2_model=reshape(reshape(v2_model,(nhours,nv2,nw)),(nhours*nv2,nw))
v2_model_err=reshape(reshape(v2_model_err,(nhours,nv2,nw)),(nhours*nv2,nw))

t3amp_model=reshape(reshape(t3amp_model,(nhours,nt3,nw)),(nhours*nt3,nw))
t3amp_model_err=reshape(reshape(t3amp_model_err,(nhours,nt3,nw)),(nhours*nt3,nw))

t3phi_model=reshape(reshape(t3phi_model,(nhours,nt3,nw)),(nhours*nt3,nw))
t3phi_model_err=reshape(reshape(t3phi_model_err,(nhours,nt3,nw)),(nhours*nt3,nw))

v2_model_stations=repmat(v2_stations,1,nhours)
t3_model_stations=repmat(t3_stations,1,nhours)

#v2_model=collect(linspace(1,0,length(v2_model)))
v2_model=transpose(v2_model)
v2_model_err=transpose(v2_model_err)
t3amp_model=transpose(t3amp_model)
t3amp_model_err=transpose(t3amp_model_err)
t3phi_model=transpose(t3phi_model)
t3phi_model_err=transpose(t3phi_model_err)

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



getchar();


staxyzplot=Array{Float64}(3,N)
for i=1:N
        staxyzplot[:,i]=station_xyzplot[i][:]
end


#tel_names =["S1","S3","S2","E1","W1","W2","E2"]
#tel_names =["S1","S2","E1","W1","W2","E2"]

fig = figure("Telescope Configuration",figsize=(32,20),facecolor="White")
clf();
ax=gca()
minorticks_on
ax[:locator_params](axis ="y", nbins=20)
ax[:locator_params](axis ="x", nbins=20)
axis("equal")
plot(staxyzplot[1,:],staxyzplot[2,:],markersize=20,marker="o",linestyle="None",markeredgecolor="green",markerfacecolor="green")
ax[:set_ylim](-30)
title("Telescope Configuration")
xlabel("Meters")
ylabel("Meters")
grid("on")
for i=1:N
    annotate(telnames[i],(staxyzplot[1,i]+10,staxyzplot[2,i]+3),size=17,color="Black")
end
annotate("",xy=[200,100],xytext=[200,25],arrowprops=Dict("facecolor"=>"black","arrowstyle"=>"simple"))
annotate("N",(195,110),size=20,color="Black")
savefig(telconfigfile)
#

#----------------------------------------------------------------------

for i=1:nv2
    baseline_list[i]=string(tel_names[v2_stations[1,i]],"-",tel_names[v2_stations[2,i]])
end


colors=["black","darkgreen","firebrick","forestgreen","peru","navy","gray","darkorange","lime","orange",
"fuchsia","saddlebrown","red","darkslateblue","blueviolet","indigo","blue","dodgerblue",
"sienna","olive","purple","darkorchid","tomato","darkturquoise","steelblue","seagreen","darkgoldenrod","darkseagreen"]


fig = figure("UV plot",figsize=(32,20),facecolor="White")
clf();
ax = axes()
minorticks_on
markeredgewidth=0.1
ax[:locator_params](axis ="y", nbins=20)
ax[:locator_params](axis ="x", nbins=20)
axis("equal")
for i=1:nv2
    #for j=1:length(λ)
        scatter(u[i,:,:]/1e6,v[i,:,:]/1e6,color=colors[i],label=baseline_list[i])
        scatter(-u[i,:,:]/1e6,-v[i,:,:]/1e6,color=colors[i])
        println,i

end
ax[:legend](bbox_to_anchor=[1.01,1.0],loc=2,borderaxespad=0)
title("UV coverage")
xlabel(L"U (M$\lambda$)")
ylabel(L"V (M$\lambda$)")
grid("on")
savefig(uvfigfile)

#-----------------------------------------
v2baseline_plot=sqrt.(uv[1,:].^2+uv[2,:].^2)
v2baseline_plot=reshape(v2baseline_plot,(nv2, nhours, nw))
v2_plot=reshape(v2_model',(nv2, nhours, nw))
v2_err_plot=reshape(v2_model_err',(nv2, nhours, nw))

for i=1:nv2
    baseline_list[i]=string(tel_names[v2_stations[1,i]],"-",tel_names[v2_stations[2,i]])
end


fig = figure("V2 data",figsize=(30,15),facecolor="White");
clf();
ax = gca();
ax[:set_yscale]("log")
ax[:set_ylim](1e-8,3.0)
minorticks_on
markeredgewidth=0.1
    for i=1:nv2
     errorbar(vec(v2baseline_plot[i,:,:])/1E6,vec(v2_plot[i,:,:]),yerr=vec(v2_err_plot[i,:,:]),fmt="o", fillstyle="full",markersize=4,ecolor="gainsboro",color=colors[i],markeredgecolor=colors[i],label=baseline_list[i])
    end
    ax[:legend](loc=0,borderaxespad=0,markerscale=3)

title("Squared Visibility Amplitude Data")
xlabel(L"Baseline (M$\lambda$)")
ylabel("Squared Visibility Amplitudes")
grid("on")
savefig(v2figfile)

end
