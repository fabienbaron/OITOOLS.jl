# EXAMPLE: classic grid search for a resolved binary
#           Note: - this is the classic version with a loop on ra,dec of the secondary
#                 - the primary visibilities are loaded from a file 
#                 - NLopt is used with the Nelder-Mead simplex to fit the flux and diameter of secondary for (ra,dec) points

using OITOOLS, SpecialFunctions, NLopt

function binary_ud_primary_centered(vis_primary::Array{Complex{Float64},1}, params::Array{Float64,1}, ra::Float64, dec::Float64, uv::Array{Float64,2}, uv_baseline::Array{Float64,1})  # flux ratio is primary/secondary
    t = params[1]/2.062648062e8*pi*uv_baseline .+1e-8;
    vis_secondary_centered = 2.0*besselj1.(t)./t
    V = (vis_primary .+ params[2] * vis_secondary_centered .* cis.(2*pi/206264806.2*(uv[1,:]*ra - uv[2,:]*dec)))/(1.0+params[2]);
return V
end

# Load Data
#
oifitsfile = "./data/2011Sep02.lam_And_prepped.oifits";
data = (readoifits(oifitsfile))[1,1]; # data can be split by wavelength, time, etc.

# Load visibilities for the primary (or compute them here from scratch)
using HDF5
vis_primary = h5read("./data/cvis_primary.h5","cvis_primary1")

#
# GRID SEARCH
#
gridside = 100; # this will be rounded to respect steps
gridstep = 0.2 #mas/pixel
hwidth = gridside*gridstep/2;
ra = collect(range(-hwidth, hwidth , step = gridstep)); # in mas
dec = collect(range(-hwidth, hwidth, step = gridstep));  ; # in mas

init_flux_ratio = .01 ; #initial guess for the flux ratio
init_diameter_secondary = 0.001; # initial guess for diameter of the secondary

visfunc=(params, uv)->binary_ud_primary_centered(vis_primary, params, 0.0, 0.0, uv, data.uv_baseline)  # flux ratio is primary/secondary
minchi2 = model_to_chi2(data,visfunc,[1.0,0.0])

chi2_map = minchi2*ones(length(ra), length(dec));
best_par = []

#
# Grid search: single thread version, only optimizes flux ratio and diameter of secondary
#
for i=1:length(ra)
    print("New row: ra = $(ra[i]) mas\n")
    for j=1:length(dec)
        visfunc=(params,uv)->binary_ud_primary_centered(vis_primary, params, ra[i], dec[j], uv, data.uv_baseline)  # flux ratio is primary/secondary 
        #chi2_map[i,j] = model_to_chi2(data,visfunc,[init_diameter_secondary,init_flux_ratio]) # if we don't want to fit
        chi2_map[i,j], opt_params, _ =  fit_model(data, visfunc, [init_diameter_secondary, init_flux_ratio], lbounds=[0, 0], hbounds=[1.0, .2], calculate_vis = false, verbose=false);
        if chi2_map[i,j] < minchi2
            minchi2 = chi2_map[i,j];
            best_par = opt_params;
            print("New best chi2: $(chi2_map[i,j]) at ra=$(ra[i]) mas and dec=$(dec[j]) mas, diameter = $(opt_params[1]), flux = $(opt_params[2])\n");
        end
    end
    # Update chi2_map view # comment out next line to speed up search
    imdisp(chi2_map, pixscale = gridstep, colormap = "gist_heat_r", figtitle="Companion search: lighter is more probable")
end
elapsed = time() - start
print(elapsed)

#
# Grid search: multi-threaded version, only optimizes flux ratio and diameter of secondary
#
chi2_map = minchi2*ones(length(ra), length(dec));
best_par = []
start = time()
Threads.@threads for i = 1:length(ra)
    for j=1:length(dec)
        visfunc=(params,uv)->binary_ud_primary_centered(vis_primary, params, ra[i], dec[j], uv, data.uv_baseline)  # flux ratio is primary/secondary 
        chi2_map[i,j], opt_params, _ =  fit_model(data, visfunc, [init_diameter_secondary, init_flux_ratio], lbounds=[0, 0], hbounds=[.4, .2], calculate_vis = false, verbose=false);
    end
end
imdisp(chi2_map, pixscale = gridstep, colormap = "gist_heat_r", figtitle="Companion search: lighter is more probable")
minchi2, radec = findmin(chi2_map)
i=radec[1]; j = radec[2]
visfunc=(params,uv)->binary_ud_primary_centered(vis_primary, params, ra[i], dec[j], uv, data.uv_baseline)  # flux ratio is primary/secondary 
chi2_map[i,j], opt_params, _ =  fit_model(data, visfunc, [init_diameter_secondary, init_flux_ratio], lbounds=[0, 0], hbounds=[.4, .2], calculate_vis = false, verbose=false);
print("Best chi2: $(chi2_map[i,j]) at ra=$(ra[i]) mas and dec=$(dec[j]) mas, diameter = $(opt_params[1]), flux = $(opt_params[2])\n")
elapsed = time() - start
print(elapsed)


function binary_ud_primary_centered_radec(vis_primary::Array{Complex{Float64},1}, params::Array{Float64,1}, uv::Array{Float64,2}, uv_baseline::Array{Float64,1})  # flux ratio is primary/secondary
    t = params[1]/2.0626480624709636e8*pi*uv_baseline .+1e-8;
    vis_secondary_centered = 2.0*besselj1.(t)./t
    V = (vis_primary .+ params[2] * vis_secondary_centered .* cis.(2*pi/206264806.2*(uv[1,:]*params[3] - uv[2,:]*params[4])))/(1.0+params[2]);
return V
end

#
# Grid search: single thread version, optimizes flux ratio, diameter of secondary and location
#
start = time()
for i=1:length(ra)
    print("New row: ra = $(ra[i]) mas\n")
    for j=1:length(dec)
        visfunc=(params,uv)->binary_ud_primary_centered_radec(vis_primary, params, uv, data.uv_baseline)  # flux ratio is primary/secondary 
        chi2_map[i,j], opt_params, _ =  fit_model(data, visfunc, [init_diameter_secondary, init_flux_ratio, ra[i], dec[j]], lbounds=[0.0, 0.0, ra[i]-2*gridstep, dec[j]-2*gridstep], hbounds=[4.0, .2, ra[i]+2*gridstep, dec[j]+2*gridstep], calculate_vis = false, verbose=false);
        if chi2_map[i,j] < minchi2
            minchi2 = chi2_map[i,j];
            best_par = opt_params;
            print("New best chi2: $(chi2_map[i,j]) at Initial: ra=$(ra[i]) mas and dec=$(dec[j]) mas / Final: ra=$(opt_params[3]) mas and dec=$(opt_params[4]) mas, diameter = $(opt_params[1]), flux = $(opt_params[2])\n");
        end
    end
    # Update chi2_map view # comment out next line to speed up search
    imdisp(chi2_map, pixscale = gridstep, colormap = "gist_heat_r", figtitle="Companion search: lighter is more probable")
end
elapsed = time() - start
print(elapsed)
