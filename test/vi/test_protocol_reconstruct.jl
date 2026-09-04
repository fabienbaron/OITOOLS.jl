# Step-3a test: exercise the new generic reconstruct_*(prob; ...) overloads on
# both InterferometricProblem and PointSourceProblem, and verify they produce
# valid output (finite, expected sizes) — algorithm correctness is already
# covered by the SkyModelParams-typed versions in test_protocol_e2e.jl.

using OITOOLS, VarInf
const OIVI = Base.get_extension(OITOOLS, :OITOOLSVarInfExt)
using .OIVI
using OITOOLS
using Statistics: mean
using LinearAlgebra
using Random

Random.seed!(7)

println("=== Step 3a generic-reconstruct test ===")

oifitsfile = "/home/baron/SOFTWARE/OITOOLS.jl/demos/data/BC2004/2004-data1.oifits"
npix = 24
pixsize = 0.4
data = readoifits(oifitsfile; filter_bad_data=true, verbose=false, warn=false)
ft = setup_ft(data, npix, pixsize)
nf = size(data, 1)
freq = [3e8 / mean(data[1, 1].uv_lam)]

# ---------------------------------------------------------------------------
# InterferometricProblem
# ---------------------------------------------------------------------------

p = SkyModelParams(npix, pixsize, freq;
                   R_mas=npix * pixsize / 2, u=0.2,
                   spatial_slope_prior=(-4.0, 1.0),
                   spatial_fluct_prior=(1.6487, 2.1612),
                   spectral_slope_prior=(-4.0, 0.0),
                   spectral_fluct_prior=(1e-30, 0.0))

obs_vec = [ObservationConfig(ft[c, 1], data[c, 1]) for c in 1:nf]
sigma_vec = vcat([vcat(o.v2_err, o.t3amp_err, o.t3phi_err) for o in obs_vec]...)
ctx = ObsContext(obs_vec, sigma_vec, nothing)
prob = InterferometricProblem(p, ft, data, ctx)
n = latent_size(prob)

println("\n-- reconstruct_map(prob; ...) --")
z_map = reconstruct_map(prob; maxiter=15, verb=true)
@assert length(z_map) == n
@assert all(isfinite, z_map)
println("  z_map: |z| = $(round(norm(z_map), digits=4))")

println("\n-- reconstruct_mgvi(prob; ...) (1 iter, 2 samples) --")
cb_count = Ref(0)
cb = (pr, z, s, i) -> (cb_count[] += 1; nothing)
z, samples = reconstruct_mgvi(prob;
                              z0=z_map,
                              n_iterations=1, n_samples=2,
                              kl_maxiter=5, cg_maxiter=10,
                              iter_callback=cb, verb=true)
@assert length(z) == n
@assert all(isfinite, z)
@assert length(samples) == 2
@assert cb_count[] == 1
println("  callback fired: $(cb_count[]) times")

println("\n-- reconstruct_geovi(prob; ...) (1 iter, 1 sample) --")
cb_count[] = 0
z, samples = reconstruct_geovi(prob;
                                z0=z_map,
                                n_iterations=1, n_samples=1,
                                kl_maxiter=3, cg_maxiter=10,
                                geo_newton_maxiter=2, geo_cg_maxiter=5,
                                iter_callback=cb, verb=true)
@assert length(z) == n
@assert all(isfinite, z)
@assert length(samples) == 2  # antithetic pair
@assert cb_count[] == 1

println("\n-- OIVI.reconstruct_hybrid(prob; ...) (1 MGVI + 1 GeoVI, 1 sample) --")
cb_count[] = 0
frozen = frozen_spectral_range(p)
z, samples = OIVI.reconstruct_hybrid(prob;
                                 z0=z_map,
                                 n_mgvi=1, n_geovi=1, n_samples=1,
                                 kl_maxiter=3, cg_maxiter=10,
                                 geo_newton_maxiter=2, geo_cg_maxiter=5,
                                 frozen_ranges=frozen,
                                 iter_callback=cb, verb=true)
@assert length(z) == n
@assert all(isfinite, z)
@assert length(samples) == 2
@assert cb_count[] == 2

# ---------------------------------------------------------------------------
# PointSourceProblem
# ---------------------------------------------------------------------------

println("\n-- PointSourceProblem: OIVI.reconstruct_hybrid(prob; ...) --")
N = 2
pos_init = [0.5 -0.5; 0.0 0.0]
ps = PointSourceParams(N, data, ft;
                        pos_init=pos_init, pos_std=0.5, logf_std=1.0)
prob_ps = PointSourceProblem(ps, data)
n_ps = latent_size(prob_ps)

z_ps_map = reconstruct_map(prob_ps; maxiter=15, verb=false)
@assert length(z_ps_map) == n_ps

cb_count[] = 0
z_ps, ps_samples = OIVI.reconstruct_hybrid(prob_ps;
                                       z0=z_ps_map,
                                       n_mgvi=1, n_geovi=1, n_samples=1,
                                       kl_maxiter=3, cg_maxiter=10,
                                       geo_newton_maxiter=2, geo_cg_maxiter=5,
                                       iter_callback=cb, verb=true)
@assert length(z_ps) == n_ps
@assert all(isfinite, z_ps)
@assert length(ps_samples) == 2
@assert cb_count[] == 2

println("\n=== Step 3a test PASSED ===")
