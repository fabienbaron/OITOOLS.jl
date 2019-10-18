#
# Image reconstruction code with temporal regularization
#
using OITOOLS
using FITSIO
include("spot_motion.jl")
#parameters for simulation of moving spots
#star params
star_params=[64,0.1,4,1,90] #params=[npix,scale,d_mean_star,majorminor,pa_star]
spot_params=[0.75 0.75 1 0 90 0.45;0.75 0.5 2 145 35 0.45] # params[x_d_spot,y_d_spot,dist_from_edge_mas,pa_loc,pa_spot,ratio]
ra=20+57/60+59.44/3600
longitude=118.0570313111
error_array=[1.0/100,1e-5,1.0/100,1e-6,0.0,0.25] # v2m,v2a,t3ampm,t3ampa,t3phim,t3phia
config_files=["./data/CHARA.txt","./data/default_obs.txt","./data/MIRC.txt","./data/MIRC_LOWH.txt"]

year=2018
hours=[8 9 10; 7 9 10; 5 6 7]
minutes=[5,10,15,30,35,40]
simulation_params=[0.1,90,4,6,15,year,hours,minutes] #speed (mas/day) span,days/month,start_mont,start_day

#make simulations
files,models=make_spot_move(star_params,spot_params,ra,longitude,error_array,config_files,simulation_params,hours,minutes)

oifitsfiles=files[:,1]
num_files=size(oifitsfiles)[1]
#oifitsfiles=["./data/spotchange1.oifits","./data/spotchange2.oifits","./data/spotchange3.oifits"]
printcolors =repeat([:red, :green, :blue],Int(num_files/3))
nepochs, tepochs, data = readoifits_multiepochs(oifitsfiles, filter_bad_data=true, force_full_t3 = true);
pixsize = Float64(star_params[2])
nx = Int(star_params[1]);
fftplan = setup_nfft_multiepochs(data, nx, pixsize);;

#initial image is a simple Gaussian
#prior_im=read(FITS("./data/4masdisk.fits")[1])
prior_im=make_disk(nx,pixsize,star_params[3],star_params[4],star_params[5])
prior_im=convert(Array{Float64,2},prior_im)
prior_im=vec(prior_im)/sum(prior_im)
x_start=prior_im=repeat(prior_im,nepochs)
regularizers = [   [ ["centering", 1e4], ["tv", 2e6] ]]  #epoch 1
for i=1:nepochs-1
    push!(regularizers,[["tv",2e6]])
end
push!(regularizers,[ ["temporal_tvsq", 3e8] ] ); #transtemporal

x = reconstruct_multitemporal(x_start, data, fftplan, printcolor= printcolors, regularizers = regularizers);

imdisp_temporal(x, nepochs, pixscale=pixsize, cmap="gist_heat",name="Reconstructions")

imdisp_temporal(vec(models),nepochs,pixscale=pixsize,cmap="gist_heat",name="Models")
