#
# EXAMPLE: classic grid search for a resolved binary
#           Output: diameters of the two stars + flux ratio + location of secondary 
#           Note: - this is the classic version with a loop on ra,dec of the secondary
#                 - NLopt is used with the Nelder-Mead simplex to fit the flux and diameter of secondary for (ra,dec) points

using OITOOLS, SpecialFunctions, Statistics, ProgressMeter

function binary_ud_fixedradec(params::Array{Float64,1}, ra::Float64, dec::Float64, uv::Array{Float64,2}, uv_baseline::Array{Float64,1})  # flux ratio is primary/secondary
     # params[1]: ud primary
     # params[2]: ud_secondary
     # params[3]: flux ratio
     t1 = params[1]/2.0626480624709636e8*pi*uv_baseline .+1e-8;
     vis_primary_centered = 2.0*besselj1.(t1)./t1
     t2 = params[2]/2.0626480624709636e8*pi*uv_baseline .+1e-8;
     vis_secondary_centered = 2.0*besselj1.(t2)./t2
     V = (vis_primary_centered .+ params[3] * vis_secondary_centered .* cis.(-2*pi/206264806.2*(uv[1,:]*ra + uv[2,:]*dec)))/(1.0+params[3]);
 return V
 end

function binary_ud(params::Array{Float64,1}, uv::Array{Float64,2}, uv_baseline::Array{Float64,1})  # flux ratio is primary/secondary
    # params[1]: ud primary
    # params[2]: ud_secondary
    # params[3]: flux ratio
    # params[4]: ra secondary
    # params[5]: dec secondary
    t1 = params[1]/2.0626480624709636e8*pi*uv_baseline .+1e-8;
    vis_primary_centered = 2.0*besselj1.(t1)./t1
    t2 = params[2]/2.0626480624709636e8*pi*uv_baseline .+1e-8;
    vis_secondary_centered = 2.0*besselj1.(t2)./t2
    V = (vis_primary_centered .+ params[3] * vis_secondary_centered .* cis.(-2*pi/206264806.2*(uv[1,:]*params[4] + uv[2,:]*params[5])))/(1.0+params[3]);
return V
end

function binary_ud_bws(params::Array{Float64,1}, uv::Array{Float64,2}, uv_baseline::Array{Float64,1})  # flux ratio is primary/secondary
    # params[1]: ud primary
    # params[2]: ud_secondary
    # params[3]: flux ratio
    # params[4]: ra secondary
    # params[5]: dec secondary
    # params[6]: bandwidth smearing
    t1 = params[1]/2.0626480624709636e8*pi*uv_baseline .+1e-8;
    vis_primary_centered = 2.0*besselj1.(t1)./t1
    t2 = params[2]/2.0626480624709636e8*pi*uv_baseline .+1e-8;
    vis_secondary_centered = 2.0*besselj1.(t2)./t2
    t3 = pi * params[6] * abs.(uv[1,:] * params[4] + uv[2,:] * params[5])/206264806.2 .+ 1e-13;
    vis_bw = sin.(t3) ./ t3; # faster than sinc() function at the moment
    V = (vis_primary_centered .+ params[3] * vis_bw .* vis_secondary_centered .* cis.(-2*pi/206264806.2*(uv[1,:]*params[4] + uv[2,:]*params[5])))/(1.0+params[3]);
return V
end



# Load Data
#
oifitsfile = "./data/AXCir.oifits";
data = readoifits(oifitsfile)[1,1]; # data can be split by wavelength, time, etc.

#
# GRID SEARCH
#
gridside = 200; # this will be rounded to respect steps
gridstep = 0.5 #mas/pixel
hwidth = gridside*gridstep/2;
ra = collect(range(-hwidth, hwidth , step = gridstep)); # in mas
dec = collect(range(-hwidth, hwidth, step = gridstep));  ; # in mas

init_diameter_primary = 1.0; # initial guess for diameter of the secondary
init_diameter_secondary = 0.0; # initial guess for diameter of the secondary
init_flux_ratio = .01 ; #initial guess for the flux ratio

minchi2 = 1e4; # placeholder value
chi2_map = ones(length(ra), length(dec));
best_par = []


# Example on how to establish a Chi2 map for a given/fixed flux ratio

model = create_model(create_component(type="ud", name="Primary"));
minf, minx, cvis_model, result = fit_model_nlopt(data, model, weights=[1.0,0,0]);
init_diameter_primary = minx[1]
bws = mean(data.v2_dlam./data.v2_lam)*2  # *2 seems needed for good results! look into this
# Estimate rmin, rmax, bws   # for bws we will need to look into OIFITS
p = Progress(length(ra), desc="Exploring RA/DEC range");
Threads.@threads for i=1:length(ra)
    #print("New row: ra = $(ra[i]) mas\n")
    for j=1:length(dec)
    visfunc=(params,uv)->binary_ud_bws(params, uv, data.uv_baseline)  # flux ratio is primary/secondary
    chi2_map[i,j] = visfunc_to_chi2(data, visfunc,[init_diameter_primary, init_diameter_secondary, init_flux_ratio, ra[i], dec[j], bws], weights=[1.0,1.0,1.0])    
    end
    next!(p)
end
finish!(p)

imdisp(chi2_map, pixsize = gridstep, colormap = "gist_heat_r", figtitle="Binary search: lighter is more probable")
best = findmin(chi2_map)
min_chi2 = best[1]
best_ra = ra[best[2][1]]
best_dec = dec[best[2][2]]
scatter(-best_ra, best_dec, alpha=0.5)

#
# Grid search: flux ratio, diameters of primary and secondary and location of secondary
#
minchi2 = 1e4; # placeholder value
chi2_map = ones(length(ra), length(dec));
visfunc=(params,uv)->binary_ud(params, uv, data.uv_baseline)  # flux ratio is primary/secondary
# To use bandwidth smearing:
# init_fracband = 0.02
# visfunc=(params,uv)->binary_ud_bws(params, uv, data.uv_baseline)  # flux ratio is primary/secondary
# Then change a line to read:
#         chi2_map[i,j], opt_params, _ =  fit_visfunc_nlopt(data, visfunc, [init_diameter_primary, init_diameter_secondary, init_flux_ratio, ra[i], dec[j], init_fracband], lbounds=[0, 0, 0, ra[i]-2*gridstep, dec[j]-2*gridstep], hbounds=[5.0, 5.0, 20.0, ra[i]+2*gridstep, dec[j]+2*gridstep], calculate_vis = false, verbose=false);
p = Progress(length(ra), desc="Exploring RA/DEC range");
Threads.@threads for i=1:length(ra)
    for j=1:length(dec)
        chi2_map[i,j], opt_params, _, _ =  fit_visfunc_nlopt(data, visfunc, [init_diameter_primary, init_diameter_secondary, init_flux_ratio, ra[i], dec[j]], lbounds=[0, 0, 0, ra[i]-2*gridstep, dec[j]-2*gridstep], hbounds=[5.0, 5.0, 20.0, ra[i]+2*gridstep, dec[j]+2*gridstep], calculate_vis = false, verbose=false);
    end
    next!(p)
end
finish!(p)

imdisp(chi2_map, pixsize = gridstep, colormap = "gist_heat_r", figtitle="Binary search: lighter is more probable")
minchi2, radec = findmin(chi2_map)
i=radec[1]; j = radec[2]
chi2_map[i,j], opt_params, _ =  fit_model_old(data, visfunc, [init_diameter_primary, init_diameter_secondary, init_flux_ratio, ra[i], dec[j]], lbounds=[0, 0, 0, ra[i]-2*gridstep, dec[j]-2*gridstep], hbounds=[5.0, 5.0, 20.0, ra[i]+2*gridstep, dec[j]+2*gridstep], calculate_vis = false, verbose=false);
print("Best chi2: $(chi2_map[i,j]) at Initial: ra=$(ra[i]) mas and dec=$(dec[j]) mas / Final: ra=$(opt_params[4]) mas and dec=$(opt_params[5]) mas, diameter_primary = $(opt_params[1]), diameter_secondary = $(opt_params[2]), flux_ratio = $(opt_params[3])\n");





    




