using OITOOLS, Test, Random, Statistics, LinearAlgebra

# The posterior samplers beside nested sampling: Ensemble Slice Sampling (core) and parallel
# tempering (`fit_model_pigeons`, the Pigeons extension), both from ROTIR. Checked on targets
# with known answers, then against Nautilus on α Cen A's power-law fit. Pigeons is a weak
# dependency CI does not install; its tests skip without it.

const _SAMPLERS_HAVE_PIGEONS = try
    @eval using Pigeons
    true
catch
    false
end
const _SAMPLERS_HAVE_NAUTILUS = try
    @eval using Nautilus
    true
catch
    false
end

@testset "posterior samplers" begin
    @testset "ensemble slice sampling: a correlated Gaussian" begin
        C = [1.0 0.8; 0.8 2.0]; P = inv(C)
        logp = x -> -0.5 * dot(x, P * x)
        rng = Xoshiro(3)
        chain, logps, dg = OITOOLS.ensemble_slice(logp, 0.1 .* randn(rng, 8, 2); nsweep = 6000,
                                                  burn = 1000, thin = 1, rng, verb = false)
        F = OITOOLS.flatten_chain(chain)
        @test size(chain) == (6000, 8, 2) && size(logps) == (6000, 8)
        @test logps[end, 3] ≈ logp(chain[end, 3, :])
        @test cov(F) ≈ C rtol = 0.1
        @test maximum(OITOOLS.rhat(chain)) < 1.01
        @test dg.tuned
        @test_throws ErrorException OITOOLS.ensemble_slice(logp, zeros(3, 2); verb = false)
    end

    @testset "autocorrelation time of AR(1)" begin
        ρ = 0.9; τ = (1 + ρ) / (1 - ρ)                     # 19 sweeps
        rng = Xoshiro(5); m, nw = 20_000, 4
        c = zeros(m, nw, 1)
        for w in 1:nw, t in 2:m
            c[t, w, 1] = ρ * c[t-1, w, 1] + sqrt(1 - ρ^2) * randn(rng)
        end
        @test OITOOLS.iat(c) ≈ τ rtol = 0.15
    end

    @testset "box transform and the gradient MALA uses" begin
        lo = [0.0, 1.0, -Inf]; hi = [2.0, Inf, Inf]
        to_θ, to_z, ljac, dθdz, dljac = OITOOLS._box_transform(lo, hi)
        z = [0.3, -0.7, 1.2]
        @test to_z(to_θ(z)) ≈ z
        θ = to_θ(z); @test 0 < θ[1] < 2 && θ[2] > 1 && θ[3] == z[3]
        h = 1e-6; e(j) = [i == j ? h : 0.0 for i in 1:3]
        @test dθdz(z) ≈ [(to_θ(z + e(j))[j] - to_θ(z - e(j))[j]) / 2h for j in 1:3] rtol = 1e-6
        @test dljac(z) ≈ [(ljac(z + e(j)) - ljac(z - e(j))) / 2h for j in 1:3] atol = 1e-6
        # the log-likelihood's closed-form gradient, with a constraint active
        d = readoifits(joinpath(@__DIR__, "..", "demos", "data", "AlphaCenA.oifits");
                       T = Float64, warn = false, verbose = false)[1, 1]
        md = Dict{String,Any}("s,ldpow" => 8.4, "s,alpha" => 0.15, "s,f" => 1.0)
        free = ["s,ldpow", "s,alpha"]
        fm, cons, idx, lower, upper, x0 = OITOOLS._sampler_setup(md, free, d;
            lb = Dict("s,ldpow" => 7.0, "s,alpha" => 0.0), ub = Dict("s,ldpow" => 10.0, "s,alpha" => 1.0),
            constraints = [("\$s,alpha", "<=", 0.1, 0.01)])
        w = [1.0, 0, 0, 0, 0, 0, 0]
        ll = OITOOLS._sampler_loglike(fm, d, w, cons, idx)
        llg = OITOOLS._sampler_loglike_grad(fm, d, w, cons, idx)
        v, g = llg(x0)
        @test v ≈ ll(x0)
        fd = [(ll(x0 .+ [i == j ? 1e-6 : 0.0 for i in 1:2]) - ll(x0 .- [i == j ? 1e-6 : 0.0 for i in 1:2])) / 2e-6
              for j in 1:2]
        @test g ≈ fd rtol = 1e-4
    end

    d = readoifits(joinpath(@__DIR__, "..", "demos", "data", "AlphaCenA.oifits");
                   T = Float64, warn = false, verbose = false)[1, 1]
    md = Dict{String,Any}("s,ldpow" => 8.4, "s,alpha" => 0.15, "s,f" => 1.0)
    free = ["s,ldpow", "s,alpha"]
    lb = Dict("s,ldpow" => 7.0, "s,alpha" => 0.0); ub = Dict("s,ldpow" => 10.0, "s,alpha" => 1.0)
    w = [1.0, 0, 0, 0, 0, 0, 0]
    # α Cen A's power-law fit, Levenberg-Marquardt: θ 8.4566, α 0.13581
    lm = fit_model_lsqfit(md, free, d; lb, ub, weights = w)

    re = fit_model_ensemble(md, free, d; lb, ub, weights = w, nsweep = 1500, seed = 11, verb = false)
    @testset "fit_model_ensemble on α Cen A" begin
        @test re isa SamplerResult && re.sampler === :ensemble
        @test size(re.posterior, 2) == 2 && re.list_free_params == free
        @test re.x_opt ≈ lm.x_opt rtol = 1e-4
        @test re.chi2 ≈ lm.chi2 rtol = 1e-6
        @test re.median ≈ lm.x_opt rtol = 1e-3
        # LM's σ is scaled by √χ²r (χ²r 4.3 here); the posterior's spread is not
        @test std(re.posterior; dims = 1)[:] ≈ lm.stderror ./ sqrt(lm.chi2r) rtol = 0.15
        @test maximum(re.diagnostics.rhat) < 1.05 && isnan(re.logz)
        r2 = fit_model_ensemble(md, free, d; lb, ub, weights = w, nsweep = 1500, seed = 11, verb = false)
        @test r2.posterior == re.posterior                          # seeded: the same run
        @test_throws ErrorException fit_model_ensemble(merge(md, Dict{String,Any}("s,alpha" => 2.0)),
                                                       free, d; lb, ub, weights = w, verb = false)
    end

    if _SAMPLERS_HAVE_PIGEONS
        @testset "fit_model_pigeons on α Cen A" begin
            # 12 rounds: the stepping-stone evidence settles with tens of round trips
            rp = fit_model_pigeons(md, free, d; lb, ub, weights = w, n_rounds = 12, seed = 2,
                                   multithreaded = false, verb = false)
            @test rp isa SamplerResult && rp.sampler === :pigeons
            @test rp.median ≈ lm.x_opt rtol = 1e-3
            @test std(rp.posterior; dims = 1)[:] ≈ lm.stderror ./ sqrt(lm.chi2r) rtol = 0.15
            @test rp.diagnostics.round_trips >= 30
            # the evidence on the nested samplers' scale
            if _SAMPLERS_HAVE_NAUTILUS
                rn = fit_model_nested(md, free, d; backend = :nautilus, lb, ub, weights = w,
                                      verb = false, cornerplot = false, seed = 4)
                @test abs(rp.logz - rn.logz) <= 3 * hypot(rp.logzerr, rn.logzerr)
                @test mean(rp.posterior; dims = 1)[:] ≈ mean(rn.posterior; dims = 1)[:] rtol = 1e-3
            end
            # MALA, on the closed-form gradient
            rm = fit_model_pigeons(md, free, d; lb, ub, weights = w, n_rounds = 8, seed = 2,
                                   explorer = :mala, multithreaded = false, verb = false)
            @test rm.median ≈ lm.x_opt rtol = 1e-3
            # an infinite bound: no evidence
            ri = fit_model_pigeons(md, free, d; lb, ub = Dict("s,ldpow" => Inf, "s,alpha" => 1.0),
                                   weights = w, n_rounds = 7, seed = 2, multithreaded = false,
                                   verb = false)
            @test isnan(ri.logz)
            @test ri.median ≈ lm.x_opt rtol = 1e-3
        end
    else
        @info "Pigeons not available — skipping fit_model_pigeons (weak dependency)"
        @test_throws ErrorException fit_model_pigeons(md, free, d; lb, ub, weights = w)
    end
end
