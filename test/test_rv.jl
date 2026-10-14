# Radial velocities fitted with the interferometry (src/radial_velocity.jl): the CSV reader, the
# residuals against PMOIRED's velocity formulas, and fits with every fitter that takes them —
# Levenberg-Marquardt (FlatModel too, as the bootstrap uses it), NLopt, the bootstrap.
using OITOOLS, Test, Random

@testset "radial velocities" begin
    d = mktempdir()
    @testset "read_rv" begin
        r = read_rv(joinpath(@__DIR__, "..", "demos", "data", "HD210763", "rv_gallenne2023.csv"))
        @test length(r) == 34 && count(==(:a), r.star) == 17 && r.sigma[1] == 0.072
        @test length(read_rv(joinpath(@__DIR__, "..", "demos", "data", "Spica", "rv_spica.csv"))) == 66
        bad = joinpath(d, "bad.csv")
        write(bad, "MJD,RV,sigma\n1,2,3\n")
        @test_throws ArgumentError read_rv(bad)
        write(bad, "MJD,RV,sigma,star\n1,2,3,c\n")
        @test_throws ArgumentError read_rv(bad)
        write(bad, "MJD,RV,sigma,star\n1,2,0,a\n")
        @test_throws ArgumentError read_rv(bad)
    end

    # a binary on an orbit, observed through V1295 Aql's uv points, with velocities of both stars
    src = joinpath(@__DIR__, "..", "demos", "data", "V1295_Aql.oifits")
    truth = Dict{String,Any}("1,ud" => 0.3, "2,f" => 0.3, "2,x" => "orbit", "2,y" => "orbit",
                             "2,orb P" => 30.0, "2,orb MJD0" => 58640.0, "2,orb e" => 0.4,
                             "2,orb omega" => 60.0, "2,orb OMEGA" => 120.0, "2,orb incl" => 50.0,
                             "2,orb a" => 4.0, "2,orb Ka" => 40.0, "2,orb Kb" => 55.0,
                             "2,orb gamma" => -12.0)
    out = joinpath(d, "orbit.oifits")
    simulate_from_oifits(src, out; mode = "copy_errors", noise = false,
                         flat_model = dict_to_model(truth, String[]), flat_params = Float64[])
    data = readoifits(out; warn = false, verbose = false, T = Float64)[1, 1]
    t = collect(range(58600.0, 58700.0; length = 12))
    vel(s) = [orbit_rv(x, s, 30.0, 58640.0, 0.4, 60.0, s === :a ? 40.0 : 55.0, -12.0) for x in t]
    rvf = joinpath(d, "rv.csv")
    open(rvf, "w") do io
        println(io, "# synthetic\nMJD,RV,sigma,star")
        for (x, v) in zip(t, vel(:a)); println(io, "$x,$v,0.1,a"); end
        for (x, v) in zip(t, vel(:b)); println(io, "$x,$v,0.2,b"); end
    end
    rv = read_rv(rvf)
    start = merge(copy(truth), Dict{String,Any}("2,orb Ka" => 35.0, "2,orb Kb" => 60.0,
                                                "2,orb gamma" => -10.0, "2,orb e" => 0.38))
    free = ["2,orb Ka", "2,orb Kb", "2,orb gamma", "2,orb e"]

    @testset "residuals" begin
        fm = dict_to_model(truth, String[])
        @test maximum(abs, OITOOLS.rv_residuals(fm, Float64[], rv, OITOOLS._rv_index(fm, "2"))) < 1e-9
        # M and q in place of Ka and Kb: Kepler III
        Ka, Kb = orbit_semiamplitudes(30.0, 0.4, 50.0, 3.0, 0.8)
        mq = merge(delete!(delete!(copy(truth), "2,orb Ka"), "2,orb Kb"),
                   Dict{String,Any}("2,orb M" => 3.0, "2,orb q" => 0.8))
        fq = dict_to_model(mq, String[])
        rq = OITOOLS.rv_residuals(fq, Float64[], rv, OITOOLS._rv_index(fq, "2"))
        @test rq ≈ [(rv.rv[k] - orbit_rv(rv.mjd[k], rv.star[k], 30.0, 58640.0, 0.4, 60.0,
                                         rv.star[k] === :a ? Ka : Kb, -12.0)) / rv.sigma[k]
                    for k in 1:length(rv)]
        @test Ka / Kb ≈ 0.8                                            # Ka/Kb = Mb/Ma = q
        @test_throws ArgumentError OITOOLS._rv_index(dict_to_model(delete!(copy(truth), "2,orb Kb"), String[]), "2")
        @test_throws ArgumentError OITOOLS._orbit_component(Dict{String,Any}("a,ud" => 1.0))
    end

    @testset "Levenberg-Marquardt and NLopt" begin
        r = fit_model_lsqfit(start, free, data; rv, lb = Dict("2,orb e" => 0.0), ub = Dict("2,orb e" => 0.99))
        p = Dict(zip(r.list_free_params, r.x_opt))
        @test isapprox(p["2,orb Ka"], 40.0; rtol = 1e-3) && isapprox(p["2,orb Kb"], 55.0; rtol = 1e-3)
        @test isapprox(p["2,orb gamma"], -12.0; atol = 0.01) && isapprox(p["2,orb e"], 0.4; rtol = 1e-3)
        @test r.ndof == OITOOLS._ndof(data, [1.0, 1.0, 1.0, 0, 0, 0, 0]) + 24
        @test 0 <= r.chi2_rv <= r.chi2
        # the same fit on the compiled model, as the bootstrap runs it
        fm = dict_to_model(start, free)
        r2 = fit_model_lsqfit(fm, [Float64(start[k]) for k in free], data; rv, rv_component = "2",
                              lb = Dict("2,orb e" => 0.0), ub = Dict("2,orb e" => 0.99))
        @test r2.x_opt ≈ r.x_opt rtol = 1e-6
        n = fit_model(start, free, data; rv, lb = Dict("2,orb e" => 0.0, "2,orb Ka" => 0.0, "2,orb Kb" => 0.0,
                      "2,orb gamma" => -50.0), ub = Dict("2,orb e" => 0.99, "2,orb Ka" => 100.0,
                      "2,orb Kb" => 100.0, "2,orb gamma" => 50.0), method = :LD_LBFGS)
        q = Dict(zip(n.list_free_params, n.x_opt))
        @test isapprox(q["2,orb Ka"], 40.0; rtol = 1e-3) && isapprox(q["2,orb gamma"], -12.0; atol = 0.02)
    end

    @testset "bootstrap resamples the velocities" begin
        rr = OITOOLS.resample_rv(rv, Xoshiro(1))
        @test length(rr) == length(rv) && issubset(rr.mjd, rv.mjd) && rr.mjd != rv.mjd
        b = bootstrap_fit(start, free, data; rv, nboot = 20, seed = 1, verb = false,
                          lb = Dict("2,orb e" => 0.0), ub = Dict("2,orb e" => 0.99))
        @test all(isfinite, b.sigma) && isapprox(b.median[1], 40.0; rtol = 1e-2)
        @test_throws ErrorException bootstrap_fit(start, free, data; rv, fitter = :nlopt, nboot = 2, verb = false)
    end
end
