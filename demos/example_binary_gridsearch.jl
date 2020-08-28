#
# EXAMPLE: clasic grid search for a resolved binary
#           Output: diameters of the two stars + flux ratio + location of secondary 
#           Note: - this is the classic version with a loop on ra,dec of the secondary
#                 - NLopt is used with the Nelder-Mead simplex to fit the flux and diameter of secondary for (ra,dec) points

using OITOOLS, SpecialFunctions, NLopt

function binary_ud_fixedradec(params::Array{Float64,1}, ra::Float64, dec::Float64, uv::Array{Float64,2}, uv_baseline::Array{Float64,1})  # flux ratio is primary/secondary
     # params[1]: ud primary
     # params[2]: ud_secondary
     # params[3]: flux ratio
     t1 = params[1]/2.0626480624709636e8*pi*uv_baseline .+1e-8;
     vis_primary_centered = 2.0*besselj1.(t1)./t1
     t2 = params[2]/2.0626480624709636e8*pi*uv_baseline .+1e-8;
     vis_secondary_centered = 2.0*besselj1.(t2)./t2
     V = (vis_primary_centered .+ params[3] * vis_secondary_centered .* cis.(2*pi/206264806.2*(uv[1,:]*ra - uv[2,:]*dec)))/(1.0+params[3]);
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
    V = (vis_primary_centered .+ params[3] * vis_secondary_centered .* cis.(2*pi/206264806.2*(uv[1,:]*params[4] - uv[2,:]*params[5])))/(1.0+params[3]);
return V
end

# Load Data
#
oifitsfile = "./data/iota_peg4t.oifits";
data = (readoifits(oifitsfile))[1,1]; # data can be split by wavelength, time, etc.

#
# GRID SEARCH
#
gridside = 100; # this will be rounded to respect steps
gridstep = 0.4 #mas/pixel
hwidth = gridside*gridstep/2;
ra = collect(range(-hwidth, hwidth , step = gridstep)); # in mas
dec = collect(range(-hwidth, hwidth, step = gridstep));  ; # in mas

init_diameter_primary = 1.0; # initial guess for diameter of the secondary
init_diameter_secondary = 1.0; # initial guess for diameter of the secondary
init_flux_ratio = .5 ; #initial guess for the flux ratio

minchi2 = 1e4; # placeholder value
chi2_map = ones(length(ra), length(dec));
best_par = []

#
# Grid search: single thread version, only optimizes flux ratio and diameters of and primary secondary
#
start = time()
for i=1:length(ra)
    print("New row: ra = $(ra[i]) mas\n")
    for j=1:length(dec)
        visfunc=(params,uv)->binary_ud_fixedradec(params, ra[i], dec[j], uv, data.uv_baseline)  # flux ratio is primary/secondary 
        chi2_map[i,j], opt_params, ~ =  fit_model(data, visfunc, [init_diameter_primary, init_diameter_secondary, init_flux_ratio], lbounds=[0, 0, 0.0001], hbounds=[10, 10, 1.0], calculate_vis = false, verbose=false);
        if chi2_map[i,j] < minchi2
            minchi2 = chi2_map[i,j];
            best_par = opt_params;
            print("New best chi2: $(chi2_map[i,j]) at ra=$(ra[i]) mas and dec=$(dec[j]) mas, diameter_primary = $(opt_params[1]), diameter_secondary = $(opt_params[2]), flux_ratio = $(opt_params[3])\n");
        end
    end
    # Update chi2_map view # comment out next line to speed up search
    if i%2==1
        imdisp(chi2_map, pixscale = gridstep)
    end
end
elapsed = time() - start


#
# Grid search: single thread version, optimizes flux ratio, diameters of primary and secondary and location of secondary
#
minchi2 = 1e4; # placeholder value
chi2_map = ones(length(ra), length(dec));

start = time()
best_par = []
visfunc=(params,uv)->binary_ud(params, uv, data.uv_baseline)  # flux ratio is primary/secondary
for i=1:length(ra)
    print("New row: ra = $(ra[i]) mas\n")
    for j=1:length(dec)
        #chi2_map[i,j]= model_to_chi2(data,visfunc,[init_diameter_primary, init_diameter_secondary, init_flux_ratio, ra[i], dec[j]]) # if we don't want to fit
         chi2_map[i,j], opt_params, ~ =  fit_model(data, visfunc, [init_diameter_primary, init_diameter_secondary, init_flux_ratio, ra[i], dec[j]], lbounds=[0, 0, 0, ra[i]-2*gridstep, dec[j]-2*gridstep], hbounds=[5.0, 5.0, 20.0, ra[i]+2*gridstep, dec[j]+2*gridstep], calculate_vis = false, verbose=false);
         if chi2_map[i,j] < minchi2
             minchi2 = chi2_map[i,j];
             best_par = opt_params;
             print("New best chi2: $(chi2_map[i,j]) at Initial: ra=$(ra[i]) mas and dec=$(dec[j]) mas / Final: ra=$(opt_params[4]) mas and dec=$(opt_params[5]) mas, diameter_primary = $(opt_params[1]), diameter_secondary = $(opt_params[2]), flux_ratio = $(opt_params[3])\n");
         end
    end
    # Update chi2_map view # comment out next line to speed up search
    imdisp(chi2_map, pixscale = gridstep, colormap = "gist_heat_r")
end
elapsed = time() - start
print(elapsed)



#
# Grid search: multiple thread version, optimizes flux ratio, diameters of primary and secondary and location of secondary
#
minchi2 = 1e4; # placeholder value
chi2_map = ones(length(ra), length(dec));

start = time()
best_par = []
visfunc=(params,uv)->binary_ud(params, uv, data.uv_baseline)  # flux ratio is primary/secondary
Threads.@threads for i=1:length(ra)
    for j=1:length(dec)
         chi2_map[i,j], opt_params, ~ =  fit_model(data, visfunc, [init_diameter_primary, init_diameter_secondary, init_flux_ratio, ra[i], dec[j]], lbounds=[0, 0, 0, ra[i]-2*gridstep, dec[j]-2*gridstep], hbounds=[5.0, 5.0, 20.0, ra[i]+2*gridstep, dec[j]+2*gridstep], calculate_vis = false, verbose=false);
    end
end
elapsed = time() - start
print(elapsed)
imdisp(chi2_map, pixscale = gridstep, colormap = "gist_heat_r")
minchi2, radec = findmin(chi2_map)
i=radec[1]; j = radec[2]
chi2_map[i,j], opt_params, ~ =  fit_model(data, visfunc, [init_diameter_primary, init_diameter_secondary, init_flux_ratio, ra[i], dec[j]], lbounds=[0, 0, 0, ra[i]-2*gridstep, dec[j]-2*gridstep], hbounds=[5.0, 5.0, 20.0, ra[i]+2*gridstep, dec[j]+2*gridstep], calculate_vis = false, verbose=false);
print("Best chi2: $(chi2_map[i,j]) at Initial: ra=$(ra[i]) mas and dec=$(dec[j]) mas / Final: ra=$(opt_params[4]) mas and dec=$(opt_params[5]) mas, diameter_primary = $(opt_params[1]), diameter_secondary = $(opt_params[2]), flux_ratio = $(opt_params[3])\n");
