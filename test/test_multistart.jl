using OITOOLS, Test, Random

# `multistart_fit`: PMOIRED's gridFit starts (grid in PMOIRED's count and order, uniform and
# normal draws, `constrain`), the same starts and minima whatever the thread count, the grouping
# into distinct minima, and an orbit's minima merged across whole periods of MJD0.
@testset "multistart_fit" begin
    @testset "starts" begin
        x0 = [0.5, 1.0, 2.0]; free = ["a", "b", "c"]
        X, names = OITOOLS._ms_starts(x0, free, ["a" => (-1, 1, 0.5), "b" => (0, 2, 1)],
                                      Pair{String,Any}[], Pair{String,Any}[], nothing, UInt64(1))
        @test size(X) == (3, 15) && names == ["a", "b"]
        @test X[1, 1:3] == [-1.0, -1.0, -1.0] && X[2, 1:3] == [0.0, 1.0, 2.0]   # the first slowest
        @test all(X[3, :] .== 2.0)                                              # not explored
        @test length(OITOOLS._ms_axis("p", (-12, 12, 1.5))) == 17               # PMOIRED's count
        @test OITOOLS._ms_axis("p", (0, 1, 0.3)) == collect(range(0, 1; length = 4))
        Xr, _ = OITOOLS._ms_starts(x0, free, Pair{String,Any}[], ["a" => (2, 3)], ["b" => (10, 0.1)], 500, UInt64(7))
        @test size(Xr, 2) == 500 && all(2 .<= Xr[1, :] .<= 3)
        @test abs(sum(Xr[2, :]) / 500 - 10) < 0.05 && all(Xr[3, :] .== 2.0)
        Xr2, _ = OITOOLS._ms_starts(x0, free, Pair{String,Any}[], ["a" => (2, 3)], ["b" => (10, 0.1)], 500, UInt64(7))
        @test Xr2 == Xr                                                         # seeded per start
        @test_throws ArgumentError OITOOLS._ms_starts(x0, free, ["d" => (0, 1, 0.5)], Pair{String,Any}[], Pair{String,Any}[], nothing, UInt64(1))
        @test_throws ArgumentError OITOOLS._ms_starts(x0, free, Pair{String,Any}[], ["a" => (0, 1)], ["a" => (0, 1)], 5, UInt64(1))
        @test_throws ArgumentError OITOOLS._ms_starts(x0, free, Pair{String,Any}[], ["a" => (0, 1)], Pair{String,Any}[], nothing, UInt64(1))
    end

    # a binary simulated through V1295 Aql's uv points, noiseless: V² and closure phases
    src = joinpath(@__DIR__, "..", "demos", "data", "V1295_Aql.oifits")
    truth = Dict{String,Any}("A,ud" => 0.4, "B,f" => 0.3, "B,ud" => 0.1, "B,x" => 3.0, "B,y" => -2.0)
    out = joinpath(mktempdir(), "binary.oifits")
    simulate_from_oifits(src, out; mode = "copy_errors", noise = false,
                         flat_model = dict_to_model(truth, String[]), flat_params = Float64[])
    d = readoifits(out; warn = false, verbose = false, T = Float64)[1, 1]
    m = merge(copy(truth), Dict{String,Any}("B,f" => 0.5, "B,x" => 0.0, "B,y" => 0.0))
    free = ["B,f", "B,x", "B,y"]
    kw = (; grid = ["B,x" => (-6, 6, 2), "B,y" => (-6, 6, 2)],
          constrain = [("sqrt(\$B,x^2 + \$B,y^2)", ">=", 1.5)], lb = Dict("B,f" => 0.0), ub = Dict("B,f" => 1.0),
          verbose = false)
    r1 = multistart_fit(m, free, d; kw..., threaded = false)
    @testset "grid, constrain, minima" begin
        @test r1.nrejected == 1                                   # the node at the origin
        @test size(r1.starts, 2) == 48 && length(r1.fits) == 48
        b = r1.minima[1].fit
        @test b.x_opt ≈ [0.3, 3.0, -2.0] rtol = 1e-4              # the truth, from the grid
        @test issorted([mm.fit.chi2r for mm in r1.minima])
        # every start that counts belongs to exactly one minimum
        grouped = reduce(vcat, [mm.starts for mm in r1.minima])
        @test sort(grouped) == findall(==(:minimum), r1.status)
        @test all(s -> s in (:minimum, :failed, :no_uncertainties, :did_not_move, :chi2), r1.status)
        @test r1.minima[1].fit.chi2r <= minimum(f.chi2r for f in r1.fits if f !== nothing) + 1e-12
    end
    @testset "threads change nothing" begin
        rt = multistart_fit(m, free, d; kw..., threaded = true)
        @test rt.starts == r1.starts
        @test [f === nothing ? nothing : f.x_opt for f in rt.fits] == [f === nothing ? nothing : f.x_opt for f in r1.fits]
        @test rt.status == r1.status
        @test [mm.starts for mm in rt.minima] == [mm.starts for mm in r1.minima]
    end
    @testset "deltachi2, bounds" begin
        r = multistart_fit(m, free, d; kw..., deltachi2 = 1.0, threaded = false)
        @test count(==(:chi2), r.status) > 0
        @test all(f -> f.chi2r <= r.minima[1].fit.chi2r + 1.0, (mm.fit for mm in r.minima))
        rc = multistart_fit(m, free, d; grid = ["B,f" => (0.5, 2.0, 0.5)], lb = Dict("B,f" => 0.0),
                            ub = Dict("B,f" => 1.0), verbose = false, threaded = false)
        @test rc.nclamped == 2 && maximum(rc.starts[1, :]) == 1.0
    end

    @testset "an orbit's minima across whole periods" begin
        torb = Dict{String,Any}("1,ud" => 0.3, "2,f" => 0.3, "2,x" => "orbit", "2,y" => "orbit",
                                "2,orb P" => 30.0, "2,orb MJD0" => "\$dMJD0 + 58640", "dMJD0" => 0.0,
                                "2,orb e" => 0.4, "2,orb omega" => 60.0, "2,orb OMEGA" => 120.0,
                                "2,orb incl" => 50.0, "2,orb a" => 4.0)
        o = joinpath(mktempdir(), "orbit.oifits")
        simulate_from_oifits(src, o; mode = "copy_errors", noise = true, seed = 3,
                             flat_model = dict_to_model(torb, String[]), flat_params = Float64[])
        dd = readoifits(o; warn = false, verbose = false, T = Float64)[1, 1]
        # starts near the truth and near the truth one period later: one minimum
        r = multistart_fit(torb, ["dMJD0", "2,orb e", "2,orb a"], dd; grid = ["dMJD0" => (-0.3, 30.3, 30.6 / 3)],
                           lb = Dict("2,orb e" => 0.0, "2,orb a" => 0.0), ub = Dict("2,orb e" => 0.99),
                           verbose = false, threaded = false)
        # with noise (σ as real data give them), the fits from the first and the last start
        # end one period apart, and are one minimum
        f1, f4 = r.fits[1], r.fits[end]
        @test abs(f1.x_opt[1]) < 0.05 && abs(f4.x_opt[1] - 30) < 0.05
        @test f1.x_opt[2:3] ≈ f4.x_opt[2:3] rtol = 1e-4
        best = r.minima[1]
        @test 1 in best.starts && length(r.fits) in best.starts
    end
end
