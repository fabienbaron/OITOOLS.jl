using OITOOLS
using LinearAlgebra
using PyPlot

#
# EXAMPLE: fast grid search for an unresolved binary, using closure phases only
#           Note: - this is not meant to be precise
#                 - a parallel implementation (splitting the grid into smaller chunks) would be an easy extension for better performance

function binary_dirac_primary_centered(flux_ratio::Float64, side::Array{Float64,1}, uv::Array{Float64,2} )  # flux ratio is primary/secondary
# Note: we use the fact the grid is a circulant matrix
nside = length(side)
nuv = size(uv,2)
# V = (1.0 .+ flux_ratio * cis.( (2*pi)* (uv[1,:] * grid_ra + uv[2,:] * grid_dec)))/(1.0+flux_ratio);  # this would be slower
V = (1.0 .+ flux_ratio * repeat( cis.(uv[1,:]*side'), 1, nside).*reshape(repeat(cis.(uv[2,:]*side'), nside), (nuv,nside*nside)) )/(1.0+flux_ratio);
return V
end

function chi2_map(cvis::Array{Complex{Float64},2}, data::OIdata)
    @inbounds t3 = cvis[data.indx_t3_1,:].*cvis[data.indx_t3_2,:].*cvis[data.indx_t3_3,:];
    #NOTE: closure phases are often the only thing needed here
    #@inbounds chi2_v2 = sum( ((abs2.(cvis[data.indx_v2,:]) .- data.v2)./data.v2_err).^2, dims=1)
    #@inbounds chi2_t3amp = sum(((abs.(t3) .- data.t3amp)./data.t3amp_err).^2, dims=1);
    @inbounds chi2_t3phi = sum( (mod360.(angle.(t3)*(180.0/pi) .- data.t3phi)./data.t3phi_err).^2, dims=1);
    return chi2_t3phi/data.nt3phi #(chi2_v2 + chi2_t3amp + chi2_t3phi)/(data.nv2+data.nt3amp+data.nt3phi)
end

#oifitsfile = "./data/HD196867_oifits_merged.fits";
oifitsfile = "./data/HD_189037_2019Aug07.fits";
data = (readoifits(oifitsfile))[1,1]; # data can be split by wavelength, time, etc.


#
# SQUARE GRID SEARCH
#
# ng = size of grid in pixels
nside = 100;
# pixsize = grid mesh size in mas per grid point
pixsize = 4.0; # 0.1 mas/grid pixel
hwidth = nside*pixsize/2;
side = collect(range(-hwidth,hwidth, length = nside)/ 206264806.2*2*pi);
#grid_ra = vec(repeat(side,1,nside))'; # -10 to 10 mas, 0.1 mas precision
#grid_dec = vec(repeat(side,1,nside)')';

# Square grid search  - binary

#flux_ratio = range(0.05,1.0, length=20); # Every 5%
flux_ratio = [0.01, 0.02, 0.05, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8]; # Semilog strategy
minchi2 = ones(length(flux_ratio));
radec = Array{Array{Float64,1}}(undef,length(flux_ratio));
for i=1:length(flux_ratio)
    cvis = binary_dirac_primary_centered(flux_ratio[i], side, data.uv);
    chi2 = reshape(chi2_map(cvis, data),nside,nside);
    minchi2[i], loc = findmin(chi2)
    radec[i] = range(-hwidth,hwidth, length = nside)[[loc[1], loc[2]]]
    print("Flux ratio: $(flux_ratio[i]) \t Chi2_min: $(minchi2[i]) \t RA,DEC: $(radec[i]) mas \t r = $(norm(radec[i])) mas\t θ = $(180/pi*atan(-radec[i][1],radec[i][2]))°\n")
end

best_chi2, ibest = findmin(minchi2)
cvis = binary_dirac_primary_centered(flux_ratio[ibest], side, data.uv);
best_chi2_map = reshape(chi2_map(cvis, data),nside,nside);
imshow(best_chi2_map, extent=[hwidth,-hwidth,-hwidth,hwidth],ColorMap(:gist_heat)) # Display map, slow
scatter(radec[ibest][2], radec[ibest][1])

print("Best solution: $best_chi2 \t Flux_ratio = $(flux_ratio[ibest]) \t RA, DEC: $(radec[ibest]) mas  \t r = $(norm(radec[ibest])) mas\t θ = $(180/pi*atan(-radec[ibest][1],radec[ibest][2]))°\n")

# Once secondary found
sec_ra = radec[ibest][1]/ 206264806.2
sec_dec = radec[ibest][2]/ 206264806.2
vis_binary = (1.0 .+ flux_ratio[ibest] *cis.(2*pi*(data.uv[1,:]*sec_ra-data.uv[2,:]*sec_dec)))/(1.0+flux_ratio[ibest]);
chi2_map(reshape(vis_binary,length(vis_binary),1), data)

#
# Continue to next step to look for triple
#

# Square grid search  - triple
