# The precision split: Float32 plans, Float64 accumulation.
#
# This guards a failure that is silent by construction. Variational inference accepts Float32
# Fourier plans, but the ADJOINT SOURCE must still be accumulated at Float64: a uv point
# collects contributions from several observables — it appears in V² and in all three legs of a
# T3 — which are large and of opposite sign, so summing them at the plan's precision loses the
# cancellation. Measured against a Float64 reference in the image kernel, gradient error at
# nx = 512 is 0.57 with a Float32 accumulator and 8.0e-04 with a Float64 one. Nothing throws
# either way; the reconstruction merely stops converging properly.
#
# So the assertion that matters is `eltype(_obs_to_g_cvis(...)) === ComplexF64` even when the
# visibilities are ComplexF32.

using OITOOLS, VarInf, Test, Random, LinearAlgebra

const E = Base.get_extension(OITOOLS, :OITOOLSVarInfExt)
const PFILE = joinpath(pkgdir(OITOOLS), "demos", "data", "BC2004", "2004-data1.oifits")
const PNX, PPS = 32, 0.5

_setup(T) = let d = readoifits(PFILE; T = T, filter_bad_data = true, verbose = false, warn = false)
    (d, setup_ft(d, PNX, PPS))
end

_params(d) = let l = unique(d[1, 1].uv_lam)
    E.SkyModelParams(PNX, Float64(PPS), [Float64(sum(l) / length(l))]; R_mas = PNX * PPS / 5)
end

@testset "both plan precisions are accepted" begin
    for T in (Float32, Float64)
        _, ft = _setup(T)
        @test E._check_ft_precision(ft) === T
        obs = E.ObservationConfig(ft[1, 1], _setup(T)[1][1, 1])
        @test E.plan_precision(obs) === T
    end
end

@testset "the adjoint source is accumulated at Float64 whatever the plans are" begin
    for T in (Float32, Float64)
        data, ft = _setup(T)
        d   = data[1, 1]
        obs = E.ObservationConfig(ft[1, 1], d)
        img = abs.(randn(Random.Xoshiro(4), PNX, PNX)) .+ 1e-3

        # The visibilities follow the plans ...
        cvis = E._combined_cvis(img, obs)
        @test eltype(cvis) === Complex{T}

        # ... and the accumulator does not. This is the whole point.
        nobs = d.nv2 + d.nt3amp + d.nt3phi
        g_cvis = E._obs_to_g_cvis(ones(Float64, nobs), cvis, obs)
        @test eltype(g_cvis) === ComplexF64

        # The gradient comes back at full precision regardless of the plans, because the
        # rounding to the plan's type happens once, inside the adjoint.
        g_img = E._g_cvis_to_g_image(g_cvis, img, obs)
        @test eltype(g_img) === Float64
        @test all(isfinite, g_img)
    end
end

@testset "a Float32 gradient tracks the Float64 one" begin
    # Not bitwise: the transform itself is Float32, whose own error is ~7.5e-06. What is being
    # checked is that nothing worse than that leaks in — a Float32 accumulator would show up
    # here as a relative difference of order 1e-1, four orders above this bound.
    gs = map((Float32, Float64)) do T
        data, ft = _setup(T)
        p    = _params(data)
        prob = E._build_interferometric_problem(p, ft, data; weights = [1.0, 1.0, 1.0], verb = false)
        z    = randn(Random.Xoshiro(7), VarInf.latent_size(prob))
        last(E.energy_and_gradient(prob, z))
    end
    @test norm(gs[1] .- gs[2]) / norm(gs[2]) < 1e-4
end
