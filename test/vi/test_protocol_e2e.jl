# Step-1b validation: a very small reconstruct_hybrid run on 2004 data, then
# a very small reconstruct_pointsource run. Both must succeed end-to-end; they
# exercise the entire chain of refactored helpers through the protocol.

using OITOOLS, VarInf
const OIVI = Base.get_extension(OITOOLS, :OITOOLSVarInfExt)
using .OIVI
using OITOOLS
using Statistics: mean
using Random

Random.seed!(42)

println("=== Step 1b end-to-end test ===")

oifitsfile = "/home/baron/SOFTWARE/OITOOLS.jl/demos/data/BC2004/2004-data1.oifits"
npix = 24
pixsize = 0.4

println("\n-- Loading data --")
data = readoifits(oifitsfile; filter_bad_data=true, verbose=false, warn=false)
ft = setup_ft(data, npix, pixsize)
nf = size(data, 1)
freq = [3e8 / mean(data[1, 1].uv_lam)]

p = SkyModelParams(npix, pixsize, freq;
                   R_mas=npix * pixsize / 2, u=0.2,
                   spatial_slope_prior=(-4.0, 1.0),
                   spatial_fluct_prior=(1.6487, 2.1612),
                   spectral_slope_prior=(-4.0, 0.0),
                   spectral_fluct_prior=(1e-30, 0.0))

frozen = nf == 1 ? frozen_spectral_range(p) : nothing

println("\n-- reconstruct_hybrid (2 MGVI + 1 GeoVI iters, 2 samples) --")
center, image_mean, image_std, samples = OIVI.reconstruct_hybrid(
    p, ft, data;
    n_mgvi=2, n_geovi=1, n_samples=2,
    map_maxiter=30, kl_maxiter=8, cg_maxiter=20,
    geo_newton_maxiter=3, geo_cg_maxiter=10,
    frozen_ranges=frozen, verb=true)

println("\nimage_mean range: ", extrema(image_mean))
println("image_std  range: ", extrema(image_std))
println("samples drawn:    ", length(samples))
@assert all(isfinite, image_mean) "image_mean has non-finite values"
@assert all(isfinite, image_std)  "image_std has non-finite values"
@assert maximum(image_mean) > 0   "image_mean is non-positive"
println("\nreconstruct_hybrid PASSED")

println("\n-- reconstruct_pointsource (2 MGVI + 1 GeoVI iters) --")
N = 2
pos_init = [0.5 -0.5; 0.0 0.0]
ps = PointSourceParams(N, data, ft;
                       pos_init=pos_init, pos_std=0.5, logf_std=1.0)
z_ps, posterior, ps_samples = reconstruct_pointsource(
    ps, data;
    n_mgvi=2, n_geovi=1, n_samples=2,
    map_maxiter=30, kl_maxiter=8, cg_maxiter=20,
    geo_newton_maxiter=3, geo_cg_maxiter=10, verb=true)

@assert all(isfinite, posterior.x_mas) "posterior x_mas has non-finite values"
@assert all(isfinite, posterior.flux)  "posterior flux has non-finite values"
println("\nreconstruct_pointsource PASSED")

println("\n=== Step 1b end-to-end test PASSED ===")
