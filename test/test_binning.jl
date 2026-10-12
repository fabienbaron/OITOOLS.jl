# Averaging spectral channels on reading (`readoifits(…; binning)`, src/oifits_binning.jl).
# The tolerances follow from the definitions: a weighted mean, its formal error, a χ² that has
# its expected value. `binmode = :pmoired` was checked against PMOIRED's own `_binVec_flag` on a
# BC2026 table (values to 1e-8, flags identical; PLAN_Tools.md), which the suite cannot run.
using OITOOLS, Test, Random, Statistics

@testset "averaging channels" begin
    O = OITOOLS
    λ = collect(range(8e-6, 13e-6; length = 120))
    Δλ = fill(λ[2] - λ[1], 120)

    @testset "block: weights, errors, the scatter check" begin
        rng = Xoshiro(7)
        nrec = 400
        σ = 0.02
        val = 0.5 .+ σ .* randn(rng, 120, nrec)
        err = fill(σ, 120, nrec); flag = falses(120, nrec)
        spec = O.BinSpec(10)
        λn, Δn, v, e, f = O._bin_table(spec, λ, Δλ, val, err, flag; stat = :v2)
        @test length(λn) == 12 && !any(f)
        @test λn ≈ [sum(λ[i:i+9]) / 10 for i in 1:10:120]
        @test Δn ≈ fill(10Δλ[1], 12)                       # the span averaged
        @test all(e .≈ σ / sqrt(10))                        # c = 1: σ/√n
        # the averaged values scatter by their stated error
        @test 0.9 < std(vec(v .- 0.5)) / (σ / sqrt(10)) < 1.1
        @test 0.9 < O.binning_scatter(spec)[:v2] < 1.1      # independent noise: χ²/dof ≈ 1

        # c: no reduction at c = n, part of it in between; capped at n
        @test all(O._bin_table(O.BinSpec(10; corr = 10), λ, Δλ, val, err, flag)[4] .≈ σ)
        @test all(O._bin_table(O.BinSpec(10; corr = 2.5), λ, Δλ, val, err, flag)[4] .≈ σ * sqrt(2.5 / 10))
        @test all(O._bin_table(O.BinSpec(10; corr = 50), λ, Δλ, val, err, flag)[4] .≈ σ)

        # a smooth chromatic trend is signal: the scatter is measured about a line in λ
        trend = val .+ 3 .* ((λ .- 1e-5) ./ 1e-5)
        s2 = O.BinSpec(10); O._bin_table(s2, λ, Δλ, trend, err, flag; stat = :v2)
        @test 0.9 < O.binning_scatter(s2)[:v2] < 1.1

        # correlated channels (white noise smoothed over 5 channels) agree better than their
        # errors say: χ²/dof well below 1
        w = randn(rng, 124, nrec)
        corr = [sum(w[i:i+4, r]) / sqrt(5) for i in 1:120, r in 1:nrec]
        s3 = O.BinSpec(10); O._bin_table(s3, λ, Δλ, 0.5 .+ σ .* corr, err, flag; stat = :v2)
        @test O.binning_scatter(s3)[:v2] < 0.5
    end

    @testset "block: weights 1/σ², flags, phases" begin
        v = reshape([1.0, 3.0], 2, 1); e = reshape([1.0, 2.0], 2, 1)
        r = O._bin_table(O.BinSpec(2), λ[1:2], Δλ[1:2], v, e, falses(2, 1))
        @test r[3][1] ≈ (1 / 1 + 3 / 4) / (1 + 1 / 4)        # inverse-variance mean
        @test r[4][1] ≈ sqrt(1 / (1 + 1 / 4))
        # one valid channel suffices; none flags the new channel
        fl = reshape([true, false], 2, 1)
        @test O._bin_table(O.BinSpec(2), λ[1:2], Δλ[1:2], v, e, fl)[3][1] == 3.0
        @test O._bin_table(O.BinSpec(2), λ[1:2], Δλ[1:2], v, e, trues(2, 1))[5][1]
        # NaN or non-positive errors do not enter
        @test O._bin_table(O.BinSpec(2), λ[1:2], Δλ[1:2], reshape([NaN, 3.0], 2, 1), e,
                           falses(2, 1))[3][1] == 3.0
        # phases across the ±180° cut average to ±180, not to 0
        p = reshape([179.0, -179.0], 2, 1)
        @test abs(abs(O._bin_table(O.BinSpec(2), λ[1:2], Δλ[1:2], p, ones(2, 1), falses(2, 1);
                                   phase = true)[3][1]) - 180) < 1e-9
        # windows in metres
        s = O.BinSpec([[8e-6, 9e-6], [12e-6, 13e-6]])
        λw, _, _ = O._bin_channels(s, λ, Δλ)
        @test length(λw) == 2 && 8e-6 < λw[1] < 9e-6 && 12e-6 < λw[2] < 13e-6
        @test_throws ArgumentError O.BinSpec(0)
        @test_throws ArgumentError O.BinSpec(10; corr = 0.5)
        @test_throws ArgumentError O.BinSpec(10; mode = :nope)
        @test_throws ArgumentError O.BinSpec([[8e-6, 9e-6]]; mode = :pmoired)
    end

    @testset "pmoired: a regular grid, errors not reduced" begin
        val = fill(0.7, 120, 3); err = fill(0.05, 120, 3)
        λn, Δn, v, e, f = O._bin_table(O.BinSpec(10; mode = :pmoired), λ, Δλ, val, err, falses(120, 3))
        @test length(λn) == 12 && λn[1] ≈ λ[1] && λn[end] ≈ λ[end]
        @test Δn ≈ fill(10Δλ[1], 12)
        @test all(v .≈ 0.7) && all(e .≈ 0.05) && !any(f)
    end

    @testset "readoifits(…; binning)" begin
        file = joinpath(@__DIR__, "..", "demos", "data", "Beauty_Contest_2026", "OBJECT1_N.oifits")
        rd(; kw...) = readoifits(file; T = Float64, warn = false, verbose = false, kw...)
        d0 = rd()[1, 1]
        st = Dict{Symbol,Float64}()
        d1 = rd(binning = 10, binstats = st)[1, 1]
        λ0 = sort(unique(d0.uv_lam)); λ1 = sort(unique(d1.uv_lam))
        @test length(λ0) == 119 && length(λ1) == 12
        @test λ1 ≈ [sum(λ0[i:min(i + 9, 119)]) / length(i:min(i + 9, 119)) for i in 1:10:119]
        @test d1.nv2 < d0.nv2 / 9 && d1.nt3phi < d0.nt3phi / 9
        # uv follows the new wavelength: the baselines in metres are the same set
        m0 = sort(unique(round.(d0.uv[1, :] .* d0.uv_lam; digits = 3)))
        m1 = sort(unique(round.(d1.uv[1, :] .* d1.uv_lam; digits = 3)))
        @test m1 ⊆ m0
        @test haskey(st, :v2) && haskey(st, :t3phi) && all(isfinite, values(st))
        # grouping after averaging: one polychromatic bin per AVERAGED channel
        @test size(rd(binning = 10, polychromatic = true), 1) == 12
        @test length(unique(rd(binning = 10, binmode = :pmoired)[1, 1].uv_lam)) == 11
        # written and read back: the averaged channels survive
        p = tempname() * ".oifits"
        try
            write_oifits(p, rd(binning = 10))
            d2 = readoifits(p; T = Float64, warn = false, verbose = false)[1, 1]
            @test sort(unique(d2.uv_lam)) ≈ λ1 rtol = 1e-6   # EFF_WAVE is single precision in OIFITS
            @test d2.nv2 == d1.nv2
        finally
            rm(p; force = true)
        end
    end
end
