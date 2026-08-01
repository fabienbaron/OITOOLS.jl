# Tests for the resampling-based uncertainty estimation (src/bootstrap.jl)

using OITOOLS, Test, Random, Statistics

const _BOOT_FILE = joinpath(@__DIR__, "oifits_for_tests", "AlphaCenA.oifits")

@testset "bootstrap" begin

data = readoifits(_BOOT_FILE)[1, 1]

@testset "data_blocks partitions the data" begin
    for g in (:config, :epoch, :point)
        b = data_blocks(data; granularity=g)
        @test length(b) > 0
        # every row of every table belongs to exactly one block
        @test sort(vcat(b.idx_v2...))   == collect(1:length(data.v2))
        @test sort(vcat(b.idx_t3...))   == collect(1:length(data.t3phi))
        @test sort(vcat(b.idx_vis...))  == collect(1:length(data.visamp))
        @test sort(vcat(b.idx_flux...)) == collect(1:length(data.flux))
    end
    # granularity ordering: point ≥ config ≥ epoch
    @test length(data_blocks(data; granularity=:point)) >=
          length(data_blocks(data; granularity=:config)) >=
          length(data_blocks(data; granularity=:epoch))
    # a :config block never mixes configurations, and holds whole spectral vectors
    b = data_blocks(data; granularity=:config)
    for idx in b.idx_v2
        isempty(idx) && continue
        @test all(data.v2_sta_index[:, i] == data.v2_sta_index[:, idx[1]] for i in idx)
        @test all(data.v2_mjd[i] ≈ data.v2_mjd[idx[1]] for i in idx)
        # all wavelengths of that (MJD, baseline) are in the same block
        @test length(unique(data.v2_lam[idx])) == length(idx)
    end
end

@testset "apply_block_counts with unit counts is the identity" begin
    b = data_blocks(data; granularity=:config)
    d = apply_block_counts(data, b, ones(Int, length(b)))
    @test d.v2      == data.v2
    @test d.v2_err  == data.v2_err
    @test d.indx_v2 == data.indx_v2
    @test d.t3phi   == data.t3phi
    @test d.t3amp   == data.t3amp
    @test d.indx_t3_1 == data.indx_t3_1
    @test d.nv2     == data.nv2
    @test d.nt3phi  == data.nt3phi
    @test d.uv      == data.uv          # uv table is carried over untouched
end

@testset "block_counts statistics" begin
    n, R = 60, 4000
    rng = Xoshiro(7)
    for (mode, evar) in ((:replacement, 1.0), (:pmoired, 0.5), (:halfsample, 0.25))
        C = hcat([block_counts(n, mode; rng=rng) for _ in 1:R]...)
        # every scheme preserves the expected amount of data...
        @test isapprox(mean(C), mode === :halfsample ? 0.5 : 1.0, atol=0.02)
        # ...but they differ in resampling variance, which sets the width of the
        # bootstrap distribution: :pmoired is half that of :replacement.
        @test isapprox(var(vec(C)), evar, rtol=0.05)
    end
    @test_throws ErrorException block_counts(10, :nonsense)
end

@testset "block_weights / apply_block_weights" begin
    rng = Xoshiro(5)
    n = 200
    W = hcat([block_weights(n; rng=rng) for _ in 1:2000]...)
    # mean 1 and variance 1: the same first two moments as the multiplicity of a
    # draw with replacement, which is what makes the multiplier bootstrap match
    # the nonparametric one
    @test isapprox(mean(W), 1.0, atol=0.02)
    @test isapprox(var(vec(W)), 1.0, rtol=0.05)
    @test all(W .> 0)

    b = data_blocks(data; granularity=:config)
    w = block_weights(length(b); rng=Xoshiro(6))
    d = apply_block_weights(data, b, w)
    @test d.v2 == data.v2                       # values untouched
    @test length(d.v2) == length(data.v2)       # nothing added or dropped
    @test d.v2_err != data.v2_err               # only the weighting changed
    # error bars scale as 1/sqrt(w) within each block
    for (i, idx) in enumerate(b.idx_v2)
        isempty(idx) && continue
        @test all(isapprox.(d.v2_err[idx], data.v2_err[idx] ./ sqrt(w[i])))
    end
    @test_throws ErrorException apply_block_weights(data, b, ones(length(b) + 1))
    @test_throws ErrorException apply_block_weights(data, b, zeros(length(b)))

    dr = resample_blocks(data, b; mode=:weights, rng=Xoshiro(6))
    @test dr.v2_err == d.v2_err
end

@testset "resample_blocks" begin
    b = data_blocks(data; granularity=:config)
    d1 = resample_blocks(data, b; mode=:replacement, rng=Xoshiro(1))
    d2 = resample_blocks(data, b; mode=:replacement, rng=Xoshiro(1))
    d3 = resample_blocks(data, b; mode=:replacement, rng=Xoshiro(2))
    @test d1.v2 == d2.v2                       # reproducible given the RNG
    @test d1.v2 != d3.v2                       # and actually random
    @test d1.nv2 == length(d1.v2)              # counts stay consistent
    @test length(d1.indx_v2) == length(d1.v2)
    @test maximum(d1.indx_v2) <= size(d1.uv, 2)
    @test d1.uv == data.uv
    # values are never modified, only selected
    @test issubset(Set(d1.v2), Set(data.v2))
    # half-sampling drops about half the data
    dh = resample_blocks(data, b; mode=:halfsample, rng=Xoshiro(3))
    @test 0.35 < (length(dh.v2) + length(dh.t3phi)) /
                 (length(data.v2) + length(data.t3phi)) < 0.65
end

@testset "perturb_data" begin
    d = perturb_data(data; rng=Xoshiro(11))
    @test length(d.v2) == length(data.v2)
    @test d.v2 != data.v2
    @test d.v2_err == data.v2_err
    @test d.uv == data.uv
    @test d.nv2 == data.nv2
    # residuals should scale like the error bars
    @test isapprox(std((d.v2 .- data.v2) ./ data.v2_err), 1.0, atol=0.15)
    @test perturb_data(data; rng=Xoshiro(11)).v2 == d.v2
end

@testset "Float32 vs Float64 data" begin
    f32 = readoifits(_BOOT_FILE)[1, 1]                  # readoifits default
    f64 = readoifits(_BOOT_FILE; T=Float64)[1, 1]
    @test eltype(f32.v2) === Float32
    @test eltype(f64.v2) === Float64

    # the whole pipeline runs in either precision
    for d in (f32, f64)
        b = data_blocks(d)
        @test length(b) > 0
        r = resample_blocks(d, b; mode=:replacement, rng=Xoshiro(1))
        @test eltype(r.v2) === eltype(d.v2)             # precision is preserved
        @test r.nv2 == length(r.v2)
        w = resample_blocks(d, b; mode=:weights, rng=Xoshiro(1))
        @test eltype(w.v2_err) === eltype(d.v2_err)
        @test all(isfinite, w.v2_err)
        @test eltype(perturb_data(d; rng=Xoshiro(1)).v2) === eltype(d.v2)
    end

    # MJDs are stored in Float64 whatever T is, because Float32 cannot resolve
    # them: near MJD 55000 its representable values are 3.9e-3 d = 5.6 min apart,
    # which would merge the blocks of exposures taken minutes apart.
    for d in (f32, f64)
        @test eltype(d.v2_mjd)   === Float64
        @test eltype(d.t3_mjd)   === Float64
        @test eltype(d.vis_mjd)  === Float64
        @test eltype(d.flux_mjd) === Float64
        @test eltype(d.uv_mjd)   === eltype(d.v2)   # uv_mjd follows T: it feeds $MJD
    end
    @test f32.v2_mjd == f64.v2_mjd
    @test f32.mean_mjd == f64.mean_mjd

    # so the block structure is identical in both precisions, at every granularity
    for g in (:config, :epoch, :point)
        @test length(data_blocks(f32; granularity=g)) ==
              length(data_blocks(f64; granularity=g))
        @test data_blocks(f32; granularity=g).keys == data_blocks(f64; granularity=g).keys
    end
end

@testset "bootstrap_fit" begin
    model = Dict{String,Any}("star,ldlin" => 8.0, "star,u" => 0.3, "star,f" => 1.0)
    free  = ["star,ldlin", "star,u"]
    lb    = Dict("star,ldlin" => 5.0, "star,u" => 0.0)
    ub    = Dict("star,ldlin" => 12.0, "star,u" => 1.0)
    w     = [1.0, 0.0, 0.0]
    kw    = (lb=lb, ub=ub, weights=w, nboot=40, seed=99, verb=false)

    b = bootstrap_fit(model, free, data; kw...)
    @test size(b.samples) == (40, 2)
    @test b.nfailed == 0
    @test all(isfinite, b.median)
    @test all(b.sigma .> 0)
    @test size(b.covar) == (2, 2)
    @test isapprox(b.correlation[1, 1], 1.0)
    # bootstrap centre agrees with the full-data fit
    @test isapprox(b.median[1], b.x_opt[1], atol=5 * b.sigma[1])
    # reproducible
    @test bootstrap_fit(model, free, data; kw...).samples == b.samples

    # these data are correlated: block bootstrap must exceed a per-point bootstrap
    p = bootstrap_fit(model, free, data; granularity=:point, kw...)
    @test b.sigma[1] > p.sigma[1]

    # every scheme should agree on well-behaved data to within its own noise
    for mode in (:halfsample, :weights, :pmoired)
        r = bootstrap_fit(model, free, data; mode=mode, kw...)
        @test all(r.sigma .> 0)
        @test isapprox(r.median[1], b.x_opt[1], atol=5 * r.sigma[1])
    end
end

end # testset bootstrap
