using OITOOLS, VarInf
const OIVI = Base.get_extension(OITOOLS, :OITOOLSVarInfExt)
using .OIVI
using .OIVI: amplitude_spectrum_adjoint, amplitude_spectrum_jvp,
                   sky_forward_jvp, _unpack, _spectral_latent_size
using LinearAlgebra
using FFTW
using OITOOLS
using Statistics: mean
using FiniteDifferences
using Test

# Test pixel size: 0.5e-9 rad = 0.5e-9 / MAS2RAD mas
const test_pixsize = 0.5e-9 / (π / (180.0 * 3600.0 * 1000.0))

@testset "Harmonic Smoothing" begin
    nx, ny = 64, 64
    dx = 1.0
    sigma = 3.0
    kernel = make_smoothing_kernel(nx, ny, dx, sigma)

    delta = zeros(nx, ny)
    delta[nx÷2+1, ny÷2+1] = 1.0
    smoothed = harmonic_smooth(delta, kernel)

    @test argmax(smoothed) == CartesianIndex(nx÷2+1, ny÷2+1)
    @test minimum(smoothed) > -1e-10
    @test abs(sum(smoothed) - 1.0) < 1e-10
    a = randn(nx, ny)
    b = randn(nx, ny)
    Sa = harmonic_smooth(a, kernel)
    Sb = harmonic_smooth(b, kernel)
    @test abs(dot(Sa, b) - dot(a, Sb)) / abs(dot(Sa, b)) < 1e-10

    println("  Smoothing tests passed.")
end

@testset "Fourier Mode Distributor" begin
    npix = 32
    dx = 1.0e-9
    bin_index, mult, mode_lengths, rel_log, log_vol, total_volume, n_bins =
        fourier_mode_distributor(npix, dx)

    @test sum(mult) ≈ npix^2
    @test minimum(bin_index) == 1
    @test maximum(bin_index) == n_bins
    @test issorted(mode_lengths)
    @test mode_lengths[1] ≈ 0.0 atol=1e-20
    @test mult[1] ≈ 1.0
    @test rel_log[1] ≈ 0.0
    @test length(log_vol) == n_bins - 2
    @test total_volume ≈ (npix * dx)^2

    println("  Fourier mode distributor tests passed.")
end

@testset "CorrFieldConfig" begin
    cfg = CorrFieldConfig(slope_prior=(-4.0, 1.0), fluct_prior=(1.0, 0.5))
    @test !cfg.use_iwp
    @test cfg.slope_mean == -4.0          # slope is a normal prior, stored as-is
    # fluct is value-space lognormal (mean, std), stored internally as (μ, σ);
    # check it round-trips to the input lognormal (mean=1.0, std=0.5).
    @test exp(cfg.fluct_mean + cfg.fluct_std^2 / 2) ≈ 1.0
    @test sqrt(exp(cfg.fluct_std^2) - 1) * exp(cfg.fluct_mean + cfg.fluct_std^2 / 2) ≈ 0.5

    cfg_iwp = CorrFieldConfig(slope_prior=(-4.0, 1.0), fluct_prior=(1.0, 0.5),
                               flex_prior=(0.37, 0.2), asp_prior=(0.14, 0.07))
    @test cfg_iwp.use_iwp

    println("  CorrFieldConfig tests passed.")
end

@testset "Limb Weight" begin
    D = limb_weight(64, 64, 10.0; u=0.5)
    @test size(D) == (64, 64)
    @test maximum(D) ≈ 1.0
    @test minimum(D) >= 0.0
    @test D[33, 33] ≈ 1.0 atol=0.05

    println("  Limb weight tests passed.")
end

@testset "Amplitude Spectrum (no IWP)" begin
    npix = 16
    freq = [1.5e14]

    p = SkyModelParams(npix, test_pixsize, freq;
                       R_mas=2.0, u=0.0,
                       slope_prior=(-4.0, 0.5),
                       fluct_prior=(1.0, 0.5))

    @test !p.spatial.use_iwp
    cfg = p.spatial
    amp = amplitude_spectrum(0.0, 0.0, 0.0, 0.0, zeros(0, 2), cfg, p)
    @test length(amp) == p.n_bins
    @test amp[1] ≈ Float64(p.npix^2)
    @test all(amp .> 0)

    # FD test for adjoint
    g_amp = randn(p.n_bins)
    v_slope, v_fluct = randn(), randn()

    g_xi_slope, g_xi_fluct, _, _, _, _ =
        amplitude_spectrum_adjoint(g_amp, 0.1, 0.2, 0.0, 0.0, zeros(0, 2), cfg, p)

    ε = 1e-7
    amp_p = amplitude_spectrum(0.1 + ε * v_slope, 0.2 + ε * v_fluct,
                               0.0, 0.0, zeros(0, 2), cfg, p)
    amp_0 = amplitude_spectrum(0.1, 0.2, 0.0, 0.0, zeros(0, 2), cfg, p)
    fd_dot = dot(g_amp, (amp_p .- amp_0) ./ ε)
    adj_dot = g_xi_slope * v_slope + g_xi_fluct * v_fluct
    rel_err = abs(fd_dot - adj_dot) / (abs(fd_dot) + abs(adj_dot) + 1e-20)
    @test rel_err < 1e-4
    println("  Amplitude spectrum adjoint: rel_err = $rel_err")
end

@testset "Amplitude Spectrum (with IWP)" begin
    npix = 16
    freq = [1.5e14]

    p = SkyModelParams(npix, test_pixsize, freq;
                       R_mas=2.0, u=0.0,
                       slope_prior=(-4.0, 0.5),
                       fluct_prior=(1.0, 0.5),
                       flex_prior=(0.37, 0.2),
                       asp_prior=(0.14, 0.07))

    @test p.spatial.use_iwp
    cfg = p.spatial
    n_iwp_steps = p.n_bins - 2
    xi_spec = 0.1 * randn(n_iwp_steps, 2)

    amp = amplitude_spectrum(0.1, 0.2, -0.1, 0.3, xi_spec, cfg, p)
    @test length(amp) == p.n_bins
    @test all(amp .> 0)

    # FD test for adjoint
    g_amp = randn(p.n_bins)
    v_slope, v_fluct, v_flex, v_asp = randn(4)
    v_spec = randn(n_iwp_steps, 2)
    xi_s, xi_f, xi_fx, xi_a = 0.1, 0.2, -0.1, 0.3

    g_xi_slope, g_xi_fluct, g_xi_flex, g_xi_asp, g_xi_spectrum, _ =
        amplitude_spectrum_adjoint(g_amp, xi_s, xi_f, xi_fx, xi_a, xi_spec, cfg, p)

    ε = 1e-7
    amp_p = amplitude_spectrum(xi_s + ε*v_slope, xi_f + ε*v_fluct,
                               xi_fx + ε*v_flex, xi_a + ε*v_asp,
                               xi_spec .+ ε.*v_spec, cfg, p)
    amp_0 = amplitude_spectrum(xi_s, xi_f, xi_fx, xi_a, xi_spec, cfg, p)
    fd_dot = dot(g_amp, (amp_p .- amp_0) ./ ε)
    adj_dot = g_xi_slope * v_slope + g_xi_fluct * v_fluct +
              g_xi_flex * v_flex + g_xi_asp * v_asp + dot(g_xi_spectrum, v_spec)
    rel_err = abs(fd_dot - adj_dot) / (abs(fd_dot) + abs(adj_dot) + 1e-20)
    @test rel_err < 1e-4
    println("  Amplitude spectrum (IWP) adjoint: rel_err = $rel_err")
end

@testset "Sky Model Forward (no IWP)" begin
    npix = 32
    freq = [1.5e14, 1.6e14, 1.7e14]

    p = SkyModelParams(npix, test_pixsize, freq;
                       R_mas=2.0, u=0.5,
                       slope_prior=(-4.0, 0.5),
                       fluct_prior=(1.0, 0.5))

    xi = randn(_latent_size(p))
    image = sky_forward(xi, p)

    @test size(image) == (npix, npix, 3)
    @test minimum(image) >= 0.0
    @test image[1, 1, 1] ≈ 0.0 atol=1e-15

    println("  Sky forward (no IWP) tests passed.")
end

@testset "Sky Model Forward (with IWP)" begin
    npix = 32
    freq = [1.5e14, 1.6e14]

    p = SkyModelParams(npix, test_pixsize, freq;
                       R_mas=2.0, u=0.5,
                       slope_prior=(-4.0, 0.5),
                       fluct_prior=(1.0, 0.5),
                       flex_prior=(0.37, 0.2),
                       asp_prior=(0.14, 0.07))

    xi = randn(_latent_size(p))
    image = sky_forward(xi, p)

    @test size(image) == (npix, npix, 2)
    @test minimum(image) >= 0.0
    @test image[1, 1, 1] ≈ 0.0 atol=1e-15

    println("  Sky forward (with IWP) tests passed.")
end

@testset "Sky Model Adjoint - no IWP (finite difference)" begin
    npix = 16
    freq = [1.5e14, 1.6e14]

    p = SkyModelParams(npix, test_pixsize, freq;
                       R_mas=2.0, u=0.3,
                       slope_prior=(-4.0, 0.5),
                       fluct_prior=(1.0, 0.5),
                       ch_slope_prior=(-3.0, 0.5),
                       ch_fluct_prior=(0.03, 0.015))

    n = _latent_size(p)
    xi = 0.1 * randn(n)
    g_image = randn(npix, npix, length(freq))

    g_xi = sky_adjoint(g_image, xi, p)

    v = randn(n)
    ε = 1e-6
    image_plus = sky_forward(xi .+ ε .* v, p)
    image_0 = sky_forward(xi, p)
    Jv = (image_plus .- image_0) ./ ε

    lhs = dot(g_xi, v)
    rhs = dot(g_image, Jv)
    rel_err = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-20)

    @test rel_err < 1e-4
    println("  Adjoint (no IWP) test: rel_err = $rel_err")
end

@testset "Sky Model Adjoint - with IWP (finite difference)" begin
    npix = 16
    freq = [1.5e14, 1.6e14]

    p = SkyModelParams(npix, test_pixsize, freq;
                       R_mas=2.0, u=0.3,
                       slope_prior=(-4.0, 0.5),
                       fluct_prior=(1.0, 0.5),
                       flex_prior=(0.37, 0.2),
                       asp_prior=(0.14, 0.07),
                       ch_slope_prior=(-3.0, 0.5),
                       ch_fluct_prior=(0.03, 0.015),
                       ch_flex_prior=(0.37, 0.2),
                       ch_asp_prior=(0.14, 0.07))

    n = _latent_size(p)
    xi = 0.1 * randn(n)
    g_image = randn(npix, npix, length(freq))

    g_xi = sky_adjoint(g_image, xi, p)

    v = randn(n)
    ε = 1e-6
    image_plus = sky_forward(xi .+ ε .* v, p)
    image_0 = sky_forward(xi, p)
    Jv = (image_plus .- image_0) ./ ε

    lhs = dot(g_xi, v)
    rhs = dot(g_image, Jv)
    rel_err = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-20)

    @test rel_err < 1e-4
    println("  Adjoint (with IWP) test: rel_err = $rel_err")
end

@testset "Sky Forward JVP (finite difference)" begin
    npix = 16
    freq = [1.5e14, 1.6e14]

    p = SkyModelParams(npix, test_pixsize, freq;
                       R_mas=2.0, u=0.3,
                       slope_prior=(-4.0, 0.5),
                       fluct_prior=(1.0, 0.5),
                       ch_slope_prior=(-3.0, 0.5),
                       ch_fluct_prior=(0.03, 0.015))

    n = _latent_size(p)
    xi = 0.1 * randn(n)
    v = randn(n)

    image, d_image = sky_forward_jvp(xi, v, p)

    ε = 1e-7
    image_p = sky_forward(xi .+ ε .* v, p)
    Jv_fd = (image_p .- image) ./ ε

    rel_err = norm(d_image .- Jv_fd) / (norm(Jv_fd) + 1e-20)
    @test rel_err < 1e-4
    println("  JVP test: rel_err = $rel_err")
end

@testset "Adjoint-JVP consistency" begin
    npix = 16
    freq = [1.5e14, 1.6e14]

    p = SkyModelParams(npix, test_pixsize, freq;
                       R_mas=2.0, u=0.3,
                       slope_prior=(-4.0, 0.5),
                       fluct_prior=(1.0, 0.5),
                       flex_prior=(0.37, 0.2),
                       asp_prior=(0.14, 0.07),
                       ch_slope_prior=(-3.0, 0.5),
                       ch_fluct_prior=(0.03, 0.015),
                       ch_flex_prior=(0.37, 0.2),
                       ch_asp_prior=(0.14, 0.07))

    n = _latent_size(p)
    xi = 0.1 * randn(n)
    v = randn(n)
    g_image = randn(npix, npix, length(freq))

    g_xi = sky_adjoint(g_image, xi, p)
    _, d_image = sky_forward_jvp(xi, v, p)

    lhs = dot(g_xi, v)
    rhs = dot(g_image, d_image)
    rel_err = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-20)
    @test rel_err < 1e-8
    println("  Adjoint-JVP consistency: rel_err = $rel_err")
end

# ============================================================================
# Observe tests — require real NFFT plans from OIFITS data
# ============================================================================

const oifitsfile = "/home/baron/SOFTWARE/OITOOLS.jl/demos/data/BC2004/2004-data1.oifits"

function _make_obs_and_image()
    data = readoifits(oifitsfile; filter_bad_data=true, verbose=false, warn=false)
    ft = setup_ft(data, 32, 0.3)
    obs = ObservationConfig(ft, data)
    npix = 32; pixsize = 0.3
    freq = [3e8 / mean(data[1,1].uv_lam)]
    p = SkyModelParams(npix, pixsize, freq; R_mas=3.0, u=0.2,
                       spatial_slope_prior=(-4.0, 1.0),
                       spatial_fluct_prior=(1.6487, 2.1612),
                       spectral_slope_prior=(-1000.0, 0.0),
                       spectral_fluct_prior=(1e-30, 0.0))
    xi = 0.1 * randn(_latent_size(p))
    image = sky_forward(xi, p)[:, :, 1]
    return obs, image
end

@testset "observe_adjoint (finite difference)" begin
    obs, image = _make_obs_and_image()
    g_obs = randn(obs.nv2 + obs.nt3amp + obs.nt3phi)

    g_image = observe_adjoint(g_obs, image, obs)

    # FD Jacobian-vector product: <J'·g_obs, v> should equal <g_obs, J·v>
    fdm = central_fdm(5, 1)
    v = randn(size(image))
    Jv = jvp(fdm, img -> observe(img, obs)[1], (image, v))

    lhs = dot(g_image, v)
    rhs = dot(g_obs, Jv)
    rel_err = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-20)
    @test rel_err < 1e-5
    println("  observe_adjoint FD test: rel_err = $rel_err")
end

@testset "observe_jvp (finite difference)" begin
    obs, image = _make_obs_and_image()
    v = randn(size(image))

    _, _, d_obs_jvp = observe_jvp(image, v, obs)

    fdm = central_fdm(5, 1)
    d_obs_fd = jvp(fdm, img -> observe(img, obs)[1], (image, v))

    rel_err = norm(d_obs_jvp .- d_obs_fd) / (norm(d_obs_fd) + 1e-20)
    @test rel_err < 1e-5
    println("  observe_jvp FD test: rel_err = $rel_err")
end

@testset "observe_adjoint / observe_jvp consistency" begin
    obs, image = _make_obs_and_image()
    v_image = randn(size(image))
    v_obs = randn(obs.nv2 + obs.nt3amp + obs.nt3phi)

    g_image = observe_adjoint(v_obs, image, obs)
    _, _, d_obs = observe_jvp(image, v_image, obs)

    lhs = dot(g_image, v_image)
    rhs = dot(v_obs, d_obs)
    rel_err = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-20)
    @test rel_err < 1e-10
    println("  observe adjoint-JVP consistency: rel_err = $rel_err")
end

# ============================================================================
# Hybrid model tests — image + parametric component
# ============================================================================

using .OIVI: _split_latent, _combined_cvis

function _make_hybrid_obs_and_image()
    data = readoifits(oifitsfile; filter_bad_data=true, verbose=false, warn=false)
    ft = setup_ft(data, 32, 0.3)

    # Simple uniform disk model with free flux and diameter
    model_dict = Dict{String,Any}(
        "star,ud" => 1.0,
        "star,f"  => 0.3
    )
    list_free = ["star,ud", "star,f"]
    model = dict_to_model(model_dict, list_free)

    obs = ObservationConfig(ft, data; model=model)
    npix = 32; pixsize = 0.3
    freq = [3e8 / mean(data[1,1].uv_lam)]
    p = SkyModelParams(npix, pixsize, freq; R_mas=3.0, u=0.2,
                       spatial_slope_prior=(-4.0, 1.0),
                       spatial_fluct_prior=(1.6487, 2.1612),
                       spectral_slope_prior=(-1000.0, 0.0),
                       spectral_fluct_prior=(1e-30, 0.0))
    xi = 0.1 * randn(_latent_size(p))
    image = sky_forward(xi, p)[:, :, 1]
    params = [1.0, 0.3]  # ud=1mas, f=0.3
    return obs, image, p, xi, params
end

@testset "Hybrid observe_adjoint (finite difference)" begin
    obs, image, p, xi, params = _make_hybrid_obs_and_image()
    g_obs = randn(obs.nv2 + obs.nt3amp + obs.nt3phi)

    g_image, g_params = observe_adjoint(g_obs, image, obs; params=params)

    fdm = central_fdm(5, 1)

    # Check image gradient
    v_img = randn(size(image))
    Jv_img = jvp(fdm, img -> observe(img, obs; params=params)[1], (image, v_img))
    lhs_img = dot(g_image, v_img)
    rhs_img = dot(g_obs, Jv_img)
    rel_err_img = abs(lhs_img - rhs_img) / (abs(lhs_img) + abs(rhs_img) + 1e-20)
    @test rel_err_img < 1e-5
    println("  Hybrid observe_adjoint (image) FD test: rel_err = $rel_err_img")

    # Check params gradient
    v_par = randn(length(params))
    Jv_par = jvp(fdm, par -> observe(image, obs; params=par)[1], (params, v_par))
    lhs_par = dot(g_params, v_par)
    rhs_par = dot(g_obs, Jv_par)
    rel_err_par = abs(lhs_par - rhs_par) / (abs(lhs_par) + abs(rhs_par) + 1e-20)
    @test rel_err_par < 1e-5
    println("  Hybrid observe_adjoint (params) FD test: rel_err = $rel_err_par")
end

@testset "Hybrid observe_jvp (finite difference)" begin
    obs, image, p, xi, params = _make_hybrid_obs_and_image()
    v_img = randn(size(image))
    v_par = randn(length(params))

    _, _, d_obs_jvp = observe_jvp(image, v_img, obs;
                                   params=params, d_params=v_par)

    fdm = central_fdm(5, 1)

    # Image tangent
    d_obs_img_fd = jvp(fdm, img -> observe(img, obs; params=params)[1], (image, v_img))
    # Params tangent
    d_obs_par_fd = jvp(fdm, par -> observe(image, obs; params=par)[1], (params, v_par))
    d_obs_fd = d_obs_img_fd .+ d_obs_par_fd

    rel_err = norm(d_obs_jvp .- d_obs_fd) / (norm(d_obs_fd) + 1e-20)
    @test rel_err < 1e-5
    println("  Hybrid observe_jvp FD test: rel_err = $rel_err")
end

@testset "Hybrid adjoint-JVP consistency" begin
    obs, image, p, xi, params = _make_hybrid_obs_and_image()
    v_img = randn(size(image))
    v_par = randn(length(params))
    v_obs = randn(obs.nv2 + obs.nt3amp + obs.nt3phi)

    g_image, g_params = observe_adjoint(v_obs, image, obs; params=params)
    _, _, d_obs = observe_jvp(image, v_img, obs;
                               params=params, d_params=v_par)

    lhs = dot(g_image, v_img) + dot(g_params, v_par)
    rhs = dot(v_obs, d_obs)
    rel_err = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-20)
    @test rel_err < 1e-10
    println("  Hybrid adjoint-JVP consistency: rel_err = $rel_err")
end

@testset "Extended latent energy_fg (finite difference)" begin
    obs, image, p, xi, params = _make_hybrid_obs_and_image()
    data = readoifits(oifitsfile; filter_bad_data=true, verbose=false, warn=false)
    ft = setup_ft(data, 32, 0.3)

    z = vcat(xi, params)
    @test total_latent_size(p, obs) == length(z)

    e, g = energy_fg(z, p, ft, data, obs)

    fdm = central_fdm(5, 1)
    g_fd = grad(fdm, z0 -> energy_fg(z0, p, ft, data, obs)[1], z)[1]

    rel_err = norm(g .- g_fd) / (norm(g_fd) + 1e-20)
    @test rel_err < 1e-4
    println("  Extended latent energy_fg FD test: rel_err = $rel_err")
end

@testset "Hybrid backward compatibility (model=nothing)" begin
    # When model=nothing, hybrid code path should give same results as image-only
    data = readoifits(oifitsfile; filter_bad_data=true, verbose=false, warn=false)
    ft = setup_ft(data, 32, 0.3)
    obs_plain = ObservationConfig(ft, data)
    obs_hybrid = ObservationConfig(ft, data; model=nothing)

    npix = 32; pixsize = 0.3
    freq = [3e8 / mean(data[1,1].uv_lam)]
    p = SkyModelParams(npix, pixsize, freq; R_mas=3.0, u=0.2,
                       spatial_slope_prior=(-4.0, 1.0),
                       spatial_fluct_prior=(1.6487, 2.1612),
                       spectral_slope_prior=(-1000.0, 0.0),
                       spectral_fluct_prior=(1e-30, 0.0))
    xi = 0.1 * randn(_latent_size(p))
    image = sky_forward(xi, p)[:, :, 1]

    obs1, err1 = observe(image, obs_plain)
    obs2, err2 = observe(image, obs_hybrid)
    @test obs1 ≈ obs2
    @test err1 ≈ err2

    g_obs = randn(length(obs1))
    g1 = observe_adjoint(g_obs, image, obs_plain)
    g2 = observe_adjoint(g_obs, image, obs_hybrid)
    @test g1 ≈ g2

    println("  Backward compatibility test passed.")
end

println("\nAll tests passed!")
