# Step-1a smoke test: confirm InterferometricProblem and PointSourceProblem
# implementations of the AbstractInferenceProblem protocol return bit-identical
# results compared to the original free-function calls they delegate to.
#
# Run with:  julia --project=. test/test_protocol_smoke.jl

using OITOOLS, VarInf
const OIVI = Base.get_extension(OITOOLS, :OITOOLSVarInfExt)
using .OIVI
using .OIVI: geovi_transformation, geovi_right_sqrt_metric,
                   geovi_left_sqrt_metric, ps_geovi_transformation,
                   ps_geovi_right_sqrt_metric, ps_geovi_left_sqrt_metric,
                   energy_fg, ps_energy_fg, total_latent_size,
                   ps_latent_size, _latent_size
using OITOOLS
using LinearAlgebra
using Statistics: mean
using Random

Random.seed!(0)

println("=== Step 1a smoke test ===")
println()

# ---------------------------------------------------------------------------
# InterferometricProblem round-trip
# ---------------------------------------------------------------------------

oifitsfile = "/home/baron/SOFTWARE/OITOOLS.jl/demos/data/BC2004/2004-data1.oifits"
npix = 32
pixsize = 0.2

println("Loading interferometric data ($oifitsfile)...")
data = readoifits(oifitsfile; filter_bad_data=true, verbose=false, warn=false)
ft = setup_ft(data, npix, pixsize)
nf = size(data, 1)
d1 = data[1, 1]
freq = [3e8 / mean(d1.uv_lam)]

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
println("  latent_size(prob)  = $(latent_size(prob))")
println("  data_size(prob)    = $(data_size(prob))")

n = latent_size(prob)
nd = data_size(prob)
z = 0.01 .* randn(n)
v_lat = randn(n)
v_dat = randn(nd)

# energy_and_gradient
e_proto, g_proto = energy_and_gradient(prob, z)
e_direct, g_direct = energy_fg(z, p, ft, data, ctx)
@assert e_proto === e_direct "energy mismatch: $e_proto vs $e_direct"
@assert g_proto == g_direct  "gradient mismatch"
println("  energy_and_gradient: ✓  (E=$(round(e_proto, digits=4)))")

# transformation
t_proto  = transformation(prob, z)
t_direct = geovi_transformation(z, p, ctx)
@assert t_proto == t_direct "transformation mismatch"
println("  transformation:      ✓  (||T||=$(round(norm(t_proto), digits=4)))")

# right_sqrt_metric
r_proto  = right_sqrt_metric(prob, z, v_lat)
r_direct = geovi_right_sqrt_metric(z, v_lat, p, ctx)
@assert r_proto == r_direct "right_sqrt_metric mismatch"
println("  right_sqrt_metric:   ✓")

# left_sqrt_metric
l_proto  = left_sqrt_metric(prob, z, v_dat)
l_direct = geovi_left_sqrt_metric(z, v_dat, p, ctx)
@assert l_proto == l_direct "left_sqrt_metric mismatch"
println("  left_sqrt_metric:    ✓")

# Adjoint sanity: <R·v_lat, v_dat> == <v_lat, L·v_dat>
ip_R = dot(r_proto, v_dat)
ip_L = dot(v_lat, l_proto)
println("  ⟨R·v, w⟩ = $ip_R")
println("  ⟨v, L·w⟩ = $ip_L")
@assert abs(ip_R - ip_L) / max(abs(ip_R), 1.0) < 1e-9 "adjoint identity broken"
println("  adjoint identity:    ✓")

println()

# ---------------------------------------------------------------------------
# PointSourceProblem round-trip
# ---------------------------------------------------------------------------

println("Building PointSourceProblem (3 sources)...")
N = 3
pos_init = [0.5 -0.3 0.0;
            0.0  0.4 -0.2]
ps = PointSourceParams(N, data, ft;
                       pos_init=pos_init, pos_std=0.5,
                       logf_std=1.0)

prob_ps = PointSourceProblem(ps, data)
println("  latent_size(prob_ps) = $(latent_size(prob_ps))  (expected $(ps_latent_size(ps)))")
println("  data_size(prob_ps)   = $(data_size(prob_ps))")

n_ps = latent_size(prob_ps)
nd_ps = data_size(prob_ps)
z_ps = 0.1 .* randn(n_ps)
v_ps_lat = randn(n_ps)
v_ps_dat = randn(nd_ps)

e_p_proto, g_p_proto = energy_and_gradient(prob_ps, z_ps)
e_p_direct, g_p_direct = ps_energy_fg(z_ps, ps, data)
@assert e_p_proto === e_p_direct
@assert g_p_proto == g_p_direct
println("  energy_and_gradient: ✓  (E=$(round(e_p_proto, digits=4)))")

t_ps_proto = transformation(prob_ps, z_ps)
t_ps_direct = ps_geovi_transformation(z_ps, ps)
@assert t_ps_proto == t_ps_direct
println("  transformation:      ✓")

r_ps_proto = right_sqrt_metric(prob_ps, z_ps, v_ps_lat)
r_ps_direct = ps_geovi_right_sqrt_metric(z_ps, v_ps_lat, ps)
@assert r_ps_proto == r_ps_direct
println("  right_sqrt_metric:   ✓")

l_ps_proto = left_sqrt_metric(prob_ps, z_ps, v_ps_dat)
l_ps_direct = ps_geovi_left_sqrt_metric(z_ps, v_ps_dat, ps)
@assert l_ps_proto == l_ps_direct
println("  left_sqrt_metric:    ✓")

ip_R_ps = dot(r_ps_proto, v_ps_dat)
ip_L_ps = dot(v_ps_lat, l_ps_proto)
println("  ⟨R·v, w⟩ = $ip_R_ps")
println("  ⟨v, L·w⟩ = $ip_L_ps")
@assert abs(ip_R_ps - ip_L_ps) / max(abs(ip_R_ps), 1.0) < 1e-9
println("  adjoint identity:    ✓")

println()
println("=== Step 1a smoke test PASSED ===")
