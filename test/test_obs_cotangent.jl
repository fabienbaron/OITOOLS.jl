# The adjoint of the cvis -> observables map, exercised on ARBITRARY cotangents.
#
# `cvis_to_chi2_fg` only ever asks for one direction: the chi2 residual. Variational inference
# asks for the same operator applied to arbitrary directions, because the Fisher metric is the
# forward map linearised about a point rather than the gradient at the data. That is why the
# scatter is factored out of `_accumulate_g_cvis!` and why it is tested here on cotangents that
# have nothing to do with chi2 -- a bug reachable only from the general path would otherwise
# sit behind a chi2 suite that passes.

using OITOOLS, Test, LinearAlgebra, Random, FiniteDifferences
using OITOOLS: scatter_obs_cotangent!, _pad_weights, _chi2_terms, _accumulate_g_cvis!

@testset "observable cotangent adjoint" begin
    Random.seed!(4)
    dm = readoifits(joinpath(@__DIR__, "..", "demos", "data",
                             "2019_v1295Aql.WL_SMOOTH.A.oifits");
                    merge_oi_wavelength = true, filter_bad_data = true, T = Float64)
    d = dm[1, 1]; nx, ps = 20, 0.35
    cell = setup_ft(dm, nx, ps)[1, 1]

    # Observables as the scatter defines them: phases in RADIANS, which is the convention the
    # cotangents are per. The stored t3phi/visphi are degrees; mixing the two is a factor
    # 180/pi that no chi2 test would catch, because chi2's own residual carries the conversion.
    function obs_of(x)
        V = OITOOLS.image_to_vis(x, cell.uv)
        (v2     = abs2.(V[d.indx_v2]),
         t3amp  = abs.(V[d.indx_t3_1] .* V[d.indx_t3_2] .* V[d.indx_t3_3]),
         t3phi  = angle.(V[d.indx_t3_1] .* V[d.indx_t3_2] .* V[d.indx_t3_3]),
         visamp = abs.(V[d.indx_vis]),
         visphi = angle.(V[d.indx_vis]))
    end

    # A linear functional of the observables with NOTHING to do with chi2.
    c = (v2     = randn(d.nv2),     t3amp  = randn(d.nt3amp), t3phi  = randn(d.nt3phi),
         visamp = randn(d.nvisamp), visphi = randn(d.nvisphi))
    L(x) = (o = obs_of(x);
            dot(c.v2, o.v2) + dot(c.t3amp, o.t3amp) + dot(c.t3phi, o.t3phi) +
            dot(c.visamp, o.visamp) + dot(c.visphi, o.visphi))

    @testset "gradient of an arbitrary functional matches finite differences" begin
        x = abs.(randn(nx, nx)) .+ 0.1
        V = OITOOLS.image_to_vis(x, cell.uv)
        g_cvis = zeros(ComplexF64, length(V))
        scatter_obs_cotangent!(g_cvis, V, d;
            g_v2 = c.v2, g_t3amp = c.t3amp, g_t3phi = c.t3phi,
            g_visamp = c.visamp, g_visphi = c.visphi)
        g = real.(adjoint(cell.uv) * conj(g_cvis))       # the NFFT convention, see chi2_flat
        flux = sum(x); g = (g .- sum(x .* g) / flux) ./ flux   # image is normalised by its flux

        for _ in 1:3
            dir = randn(nx, nx)
            num = central_fdm(5, 1)(t -> L(x .+ t .* dir), 0.0)
            @test isapprox(dot(g, dir), num; rtol = 1e-6)
        end
    end

    @testset "chi2 is the special case, not a second implementation" begin
        # Feed the scatter the chi2 cotangents by hand and require the chi2 path to agree.
        x = abs.(randn(nx, nx)) .+ 0.1
        V = OITOOLS.image_to_vis(x, cell.uv)
        w = _pad_weights([1.0, 1, 1, 1, 1, 0, 0], Float64)
        t = _chi2_terms(V, d, w, false)

        viachi2 = zeros(ComplexF64, length(V))
        _accumulate_g_cvis!(viachi2, V, t, d, w)

        general = zeros(ComplexF64, length(V))
        scatter_obs_cotangent!(general, V, d;
            g_v2 = t.r_v2, scale_v2 = 2 * w[1],
            g_t3amp = t.r_t3amp, scale_t3amp = w[2],
            g_t3phi = t.r_t3phi, scale_t3phi = w[3],
            g_visamp = t.r_visamp, scale_visamp = w[4],
            g_visphi = t.r_visphi, scale_visphi = w[5],
            V1 = t.V1, V2 = t.V2, V3 = t.V3, Vvis = t.Vvis)
        @test general == viachi2
    end

    @testset "an absent cotangent contributes nothing" begin
        x = abs.(randn(nx, nx)) .+ 0.1
        V = OITOOLS.image_to_vis(x, cell.uv)
        none = zeros(ComplexF64, length(V))
        scatter_obs_cotangent!(none, V, d)          # every cotangent defaulted to nothing
        @test all(iszero, none)

        only_v2 = zeros(ComplexF64, length(V))
        scatter_obs_cotangent!(only_v2, V, d; g_v2 = c.v2)
        @test count(!iszero, only_v2) > 0
        @test all(iszero, only_v2[setdiff(eachindex(only_v2), d.indx_v2)])
    end
end
