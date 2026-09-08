using OITOOLS, VarInf
const OIVI = Base.get_extension(OITOOLS, :OITOOLSVarInfExt)
using .OIVI
using .OIVI: ps_cvis, ps_cvis_jvp, ps_cvis_adjoint, ps_observables,
                   ps_unpack, ps_unpack_tangent, ps_adjoint_chain!,
                   ps_geovi_transformation, ps_geovi_right_sqrt_metric,
                   ps_geovi_left_sqrt_metric,
                   _obs_to_g_cvis, _cvis_to_obs_jvp, MAS2RAD
using VarInf: _posterior_metric_mul
using LinearAlgebra
using Printf
using OITOOLS
using Statistics: mean
using FiniteDifferences
using Test

# ============================================================================
# Test setup: create a small synthetic dataset
# ============================================================================

function make_test_ps()
    # Use the OITOOLS demo data for testing
    oifitsfile = joinpath(dirname(pathof(OITOOLS)), "..", "demos", "data",
                          "2019_v1295Aql.WL_SMOOTH.A.oifits")
    if !isfile(oifitsfile)
        @warn "Test OIFITS file not found, skipping point source tests"
        return nothing, nothing, nothing
    end

    data = readoifits(oifitsfile; filter_bad_data=true, verbose=false,
                      warn=false, polychromatic=true, merge_oi_wavelength=true, T=Float64)
    ft = setup_ft(data, 64, 0.2)
    nwav = size(data, 1)

    N = 3
    pos_init = [0.0 0.5 -0.3;
                0.0 0.3 -0.4]
    ps = PointSourceParams(N, data, ft;
                           pos_init=pos_init, pos_std=0.5, logf_std=0.5)
    return ps, data, ft
end


# ============================================================================
# Tests
# ============================================================================

@testset "Point Source Model" begin
    ps, data, ft = make_test_ps()
    if ps === nothing
        @warn "Skipping point source tests (no test data)"
        return
    end

    N, nf = ps.N, ps.nf
    n = ps_latent_size(ps)
    @test n == 2 * N + N * nf

    z = randn(n) * 0.3
    x_rad, y_rad, flux = ps_unpack(z, ps)
    @test length(x_rad) == N
    @test length(y_rad) == N
    @test size(flux) == (N, nf)
    @test all(flux .> 0)  # exp ensures positivity

    println("  Latent layout tests passed.")
end


@testset "Point Source DFT" begin
    ps, data, ft = make_test_ps()
    ps === nothing && return

    N, nf = ps.N, ps.nf
    n = ps_latent_size(ps)
    z = randn(n) * 0.3
    x_rad, y_rad, flux = ps_unpack(z, ps)

    for c in 1:nf
        obs = ps.obs_vec[c]
        cvis = ps_cvis(x_rad, y_rad, flux[:, c], obs)

        # V(0) should be 1 (flux normalization)
        nuv = size(obs.uv, 2)
        @test length(cvis) == nuv

        # All |V| <= 1 for normalized point sources
        @test all(abs.(cvis) .<= 1.0 + 1e-10)

        # Observable extraction
        obs_model, obs_err = ps_observables(cvis, obs)
        @test length(obs_model) == obs.nv2 + obs.nt3amp + obs.nt3phi
        @test all(obs_err .> 0)
    end

    println("  DFT forward tests passed.")
end


@testset "Point Source JVP vs Finite Differences" begin
    ps, data, ft = make_test_ps()
    ps === nothing && return

    N, nf = ps.N, ps.nf
    n = ps_latent_size(ps)
    z = randn(n) * 0.3
    v = randn(n) * 0.01
    x_rad, y_rad, flux = ps_unpack(z, ps)
    dx_rad, dy_rad, dflux = ps_unpack_tangent(v, ps, flux)

    eps_fd = 1e-7
    for c in 1:min(nf, 2)  # test first 2 channels for speed
        obs = ps.obs_vec[c]
        cvis, d_cvis = ps_cvis_jvp(x_rad, y_rad, flux[:, c],
                                     dx_rad, dy_rad, dflux[:, c], obs)

        # Finite-difference check on full z -> cvis(channel c) map
        function cvis_of_z(zz)
            xr, yr, fl = ps_unpack(zz, ps)
            return ps_cvis(xr, yr, fl[:, c], obs)
        end

        cvis_plus  = cvis_of_z(z .+ eps_fd .* v)
        cvis_minus = cvis_of_z(z .- eps_fd .* v)
        d_cvis_fd = (cvis_plus .- cvis_minus) ./ (2 * eps_fd)

        rel_err = norm(d_cvis .- d_cvis_fd) / (norm(d_cvis_fd) + 1e-30)
        @test rel_err < 1e-4
        @printf("  ch%d JVP rel error: %.2e\n", c, rel_err)
    end

    println("  JVP tests passed.")
end


@testset "Point Source VJP vs Finite Differences" begin
    ps, data, ft = make_test_ps()
    ps === nothing && return

    N, nf = ps.N, ps.nf
    n = ps_latent_size(ps)
    z = randn(n) * 0.3
    x_rad, y_rad, flux = ps_unpack(z, ps)

    eps_fd = 1e-7
    for c in 1:min(nf, 2)
        obs = ps.obs_vec[c]
        cvis = ps_cvis(x_rad, y_rad, flux[:, c], obs)
        nuv = length(cvis)
        g_cvis = randn(ComplexF64, nuv)

        g_x, g_y, g_flux = ps_cvis_adjoint(g_cvis, x_rad, y_rad, flux[:, c], obs)

        # Build full gradient w.r.t. z via chain rule
        g_z_adj = zeros(n)
        ps_adjoint_chain!(g_z_adj, g_x, g_y, g_flux, flux, ps, c)

        # Finite-difference gradient of Re(g_cvis^H * cvis(z))
        function loss_of_z(zz)
            xr, yr, fl = ps_unpack(zz, ps)
            cv = ps_cvis(xr, yr, fl[:, c], obs)
            return real(dot(g_cvis, cv))
        end

        g_z_fd = zeros(n)
        for j in 1:n
            z_p = copy(z); z_p[j] += eps_fd
            z_m = copy(z); z_m[j] -= eps_fd
            g_z_fd[j] = (loss_of_z(z_p) - loss_of_z(z_m)) / (2 * eps_fd)
        end

        rel_err = norm(g_z_adj .- g_z_fd) / (norm(g_z_fd) + 1e-30)
        @test rel_err < 1e-4
        @printf("  ch%d VJP rel error: %.2e\n", c, rel_err)
    end

    println("  VJP tests passed.")
end


@testset "JVP-VJP Consistency" begin
    ps, data, ft = make_test_ps()
    ps === nothing && return

    N, nf = ps.N, ps.nf
    n = ps_latent_size(ps)
    z = randn(n) * 0.3
    v = randn(n) * 0.1

    x_rad, y_rad, flux = ps_unpack(z, ps)
    dx_rad, dy_rad, dflux = ps_unpack_tangent(v, ps, flux)

    for c in 1:min(nf, 2)
        obs = ps.obs_vec[c]

        # JVP: z-space → data-space
        cvis, d_cvis = ps_cvis_jvp(x_rad, y_rad, flux[:, c],
                                     dx_rad, dy_rad, dflux[:, c], obs)
        _, _, d_obs = _cvis_to_obs_jvp(cvis, d_cvis, obs)

        # VJP: data-space → z-space
        w = randn(length(d_obs))
        g_cvis = _obs_to_g_cvis(w, cvis, obs)
        g_x, g_y, g_flux = ps_cvis_adjoint(g_cvis, x_rad, y_rad, flux[:, c], obs)
        g_z = zeros(n)
        ps_adjoint_chain!(g_z, g_x, g_y, g_flux, flux, ps, c)

        # Consistency: <w, J*v> == <J'*w, v>
        lhs = dot(w, d_obs)
        rhs = dot(g_z, v)
        rel_err = abs(lhs - rhs) / (abs(lhs) + 1e-30)
        @test rel_err < 1e-6
        @printf("  ch%d JVP-VJP consistency: %.2e\n", c, rel_err)
    end

    println("  JVP-VJP consistency tests passed.")
end


@testset "Energy Gradient vs Finite Differences" begin
    ps, data, ft = make_test_ps()
    ps === nothing && return

    n = ps_latent_size(ps)
    z = randn(n) * 0.3

    e0, g0 = ps_energy_fg(z, ps, data; weights=[1.0, 0.0, 1.0])
    @test isfinite(e0)
    @test length(g0) == n

    eps_fd = 1e-6
    g_fd = zeros(n)
    for j in 1:n
        z_p = copy(z); z_p[j] += eps_fd
        z_m = copy(z); z_m[j] -= eps_fd
        e_p, _ = ps_energy_fg(z_p, ps, data; weights=[1.0, 0.0, 1.0])
        e_m, _ = ps_energy_fg(z_m, ps, data; weights=[1.0, 0.0, 1.0])
        g_fd[j] = (e_p - e_m) / (2 * eps_fd)
    end

    rel_err = norm(g0 .- g_fd) / (norm(g_fd) + 1e-30)
    @test rel_err < 1e-4
    @printf("  Energy gradient rel error: %.2e\n", rel_err)

    println("  Energy gradient tests passed.")
end


@testset "GeoVI Metric Symmetry" begin
    ps, data, ft = make_test_ps()
    ps === nothing && return

    n = ps_latent_size(ps)
    z = randn(n) * 0.3
    v = randn(n)
    w = randn(n)

    prob = PointSourceProblem(ps, data)
    Mv = _posterior_metric_mul(prob, z, v)
    Mw = _posterior_metric_mul(prob, z, w)

    # M should be symmetric: <Mv, w> == <v, Mw>
    lhs = dot(Mv, w)
    rhs = dot(v, Mw)
    rel_err = abs(lhs - rhs) / (abs(lhs) + 1e-30)
    @test rel_err < 1e-6
    @printf("  Metric symmetry rel error: %.2e\n", rel_err)

    # M should be positive definite: <v, Mv> > 0
    vMv = dot(v, Mv)
    @test vMv > 0

    println("  GeoVI metric symmetry tests passed.")
end


@testset "GeoVI Transformation JVP-VJP Consistency" begin
    ps, data, ft = make_test_ps()
    ps === nothing && return

    n = ps_latent_size(ps)
    n_data = length(ps.sigma_vec)
    z = randn(n) * 0.3
    v = randn(n) * 0.1
    w = randn(n_data)

    Jv = ps_geovi_right_sqrt_metric(z, v, ps)
    JTw = ps_geovi_left_sqrt_metric(z, w, ps)

    lhs = dot(w, Jv)
    rhs = dot(JTw, v)
    rel_err = abs(lhs - rhs) / (abs(lhs) + 1e-30)
    @test rel_err < 1e-6
    @printf("  Transformation JVP-VJP consistency: %.2e\n", rel_err)

    println("  Transformation consistency tests passed.")
end
