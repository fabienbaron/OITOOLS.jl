# Binary orbits with PMOIRED's keys and PMOIRED's definition (src/orbit.jl): positions and radial
# velocities against values computed by PMOIRED's own `_orbit` (test/data/pmoired_orbit.txt,
# e < 0.5, where its 20 fixed-point Kepler iterations have converged), Kepler's derivatives, the
# `c,x = c,y = "orbit"` expansion, and an orbit simulated through real uv points and fitted back.
using OITOOLS, Test, DelimitedFiles, ForwardDiff

@testset "orbits (PMOIRED's _orbit)" begin
    @testset "positions and RVs against PMOIRED" begin
        M = readdlm(joinpath(@__DIR__, "data", "pmoired_orbit.txt"))
        # columns: P MJD0 e omega OMEGA incl a Ka Kb gamma t x y VA VB (PMOIRED's output)
        for r in eachrow(M)
            P, MJD0, e, ω, Ω, i, a, Ka, Kb, γ, t, x, y, va, vb = r
            @test OITOOLS.orbit_x(t, P, MJD0, e, ω, Ω, i, a) ≈ x atol = 1e-4
            @test OITOOLS.orbit_y(t, P, MJD0, e, ω, Ω, i, a) ≈ y atol = 1e-4
            @test OITOOLS.orbit_rv(t, :a, P, MJD0, e, ω, Ka, γ) ≈ va atol = 1e-3
            @test OITOOLS.orbit_rv(t, :b, P, MJD0, e, ω, Kb, γ) ≈ vb atol = 1e-3
        end
    end

    @testset "Kepler: solution and derivatives" begin
        for e in (0.0, 0.3, 0.9, 0.999), Mv in (-3.1, -0.01, 0.0, 1e-6, 2.0, 50.0)
            E = OITOOLS.kepler_E(Mv, e)
            @test E - e * sin(E) ≈ Mv atol = 1e-12
        end
        # ROTIR's stall case: |ΔE| held one ulp above the tolerance at e = 0.99
        E = OITOOLS.kepler_E(0.0380132711084365, 0.99)
        @test abs(E - 0.99 * sin(E) - 0.0380132711084365) < 1e-15
        g = ForwardDiff.gradient(p -> OITOOLS.kepler_E(p[1], p[2]), [1.3, 0.7])
        E = OITOOLS.kepler_E(1.3, 0.7)
        @test g ≈ [1, sin(E)] ./ (1 - 0.7 * cos(E)) rtol = 1e-12          # implicit derivative
        el = [42.38, 59785.45, 0.624, 114.08, 77.34, 71.03, 3.743]
        f(p) = OITOOLS.orbit_x(59800.3, p...)
        h = 1e-6 .* max.(abs.(el), 1.0)
        fd = [(f(el .+ h .* (1:7 .== k)) - f(el .- h .* (1:7 .== k))) / 2h[k] for k in 1:7]
        @test ForwardDiff.gradient(f, el) ≈ fd rtol = 1e-5
        @test_throws ArgumentError OITOOLS.kepler_E(1.0, 1.0)
    end

    @testset "x = y = \"orbit\"" begin
        base = Dict{String,Any}("1,ud" => 0.2, "2,f" => 0.4, "2,x" => "orbit", "2,y" => "orbit",
                                "2,orb P" => 42.38, "2,orb MJD0" => 59785.45, "2,orb e" => 0.62,
                                "2,orb omega" => 114.0, "2,orb OMEGA" => 77.3, "2,orb incl" => 71.0,
                                "2,orb a" => 3.74, "2,orb Ka" => 50.3, "2,orb Kb" => 58.8,
                                "2,orb gamma" => 14.9)
        @test isempty(OITOOLS.model_warnings(base, String[]))              # Ka, Kb, γ are not "unused"
        fm = dict_to_model(base, String[])
        uv = [3e7 1e7; -2e7 4e7]; t = [59800.0, 59811.5]
        V = eval_model(fm, Float64[], uv; mjd = t)
        # the same as a component placed at the orbit's position by hand, point by point
        for k in 1:2
            fixed = merge(copy(base), Dict{String,Any}("2,x" => OITOOLS.orbit_x(t[k], 42.38, 59785.45, 0.62, 114.0, 77.3, 71.0, 3.74),
                                                       "2,y" => OITOOLS.orbit_y(t[k], 42.38, 59785.45, 0.62, 114.0, 77.3, 71.0, 3.74)))
            @test eval_model(dict_to_model(fixed, String[]), Float64[], uv[:, k:k])[1] ≈ V[k]
        end
        # a from the total mass and the parallax (Kepler III)
        mp = delete!(merge(copy(base), Dict{String,Any}("2,orb M" => 3.2248, "2,orb plx" => 10.692)), "2,orb a")
        Vm = eval_model(dict_to_model(mp, String[]), Float64[], uv; mjd = t)
        a3 = OITOOLS.orbit_a(42.38, 3.2248, 10.692)
        Va = eval_model(dict_to_model(merge(copy(base), Dict{String,Any}("2,orb a" => a3)), String[]), Float64[], uv; mjd = t)
        @test Vm ≈ Va
        # refusals, each naming what is wrong
        e1 = try dict_to_model(merge(copy(base), Dict{String,Any}("2,y" => 1.0)), String[]); nothing catch e; e end
        @test e1 isa ArgumentError && occursin("both", e1.msg)
        e2 = try dict_to_model(delete!(copy(base), "2,orb P"), String[]); nothing catch e; e end
        @test e2 isa ArgumentError && occursin("2,orb P", e2.msg)
        e3 = try dict_to_model(delete!(copy(base), "2,orb a"), String[]); nothing catch e; e end
        @test e3 isa ArgumentError && occursin("semi-major", e3.msg)
        lb, ub = default_bounds(base, ["2,orb e", "2,orb incl", "2,orb P"])
        @test (lb["2,orb e"], ub["2,orb e"]) == (0.0, 0.99) && ub["2,orb incl"] == 180.0 && lb["2,orb P"] == 0.0
    end

    @testset "expressions of time alone, once per epoch" begin
        m = Dict{String,Any}("1,ud" => "0.2 + 0.01 * (\$MJD - 59800)", "2,f" => "0.4 + 0.1 * \$WL",
                             "2,ud" => "0.1 + 0.001 * (\$MJD - 59800) * \$WL",
                             "2,x" => "orbit", "2,y" => "orbit", "2,orb P" => "\$dP + 42.38", "dP" => 0.0,
                             "2,orb MJD0" => 59785.45, "2,orb e" => 0.62, "2,orb omega" => 114.0,
                             "2,orb OMEGA" => 77.3, "2,orb incl" => 71.0, "2,orb a" => 3.74)
        free = ["dP", "2,orb e", "2,orb a"]
        fm = dict_to_model(m, free); x = [0.01, 0.62, 3.74]
        # runs, a repeated epoch out of order, a run of one
        t  = [59800.0, 59800.0, 59800.0, 59811.5, 59811.5, 59800.0, 59830.25]
        n  = length(t)
        uv = hcat([[3e7 + 1e6k, -2e7 + 3e6k] for k in 1:n]...)
        wl = fill(2.2e-6, n) .+ 1e-8 .* (1:n)
        V, J = eval_model_jac(fm, x, uv; wl, mjd = t)
        @test eval_model(fm, x, uv; wl, mjd = t) ≈ V rtol = 1e-14
        for k in 1:n                                       # each point alone: a run of one
            Vk, Jk = eval_model_jac(fm, x, uv[:, k:k]; wl = wl[k:k], mjd = t[k:k])
            @test Vk[1] == V[k]
            @test Jk[1, :] == J[k, :]
        end
        @test isempty(eval_model(fm, x, zeros(2, 0); wl = Float64[], mjd = Float64[]))
        v = fm.resolver.fn(x, wl, t)                       # the value path, as constraints read it
        i = findfirst(==("1,ud"), fm.resolver.all_names)
        @test v[i] == 0.2 .+ 0.01 .* (t .- 59800)
        @test fm.resolver.fn(x)[i] isa Float64             # MJD = NaN: one value
        tn = [59800.0, NaN, NaN, 59811.5]                  # NaN never equals itself: each its own
        vn = fm.resolver.fn(x, wl[1:4], tn)[i]
        @test vn[1] == 0.2 && isnan(vn[2]) && isnan(vn[3]) && vn[4] ≈ 0.315
        # Float32 parameters and uv with Float64 MJDs (as `readoifits` gives them): the model in Float32
        x32 = Float32.(x); uv32 = Float32.(uv); f32 = (wl = Float32.(wl), mjd = t)
        @test eltype(eval_model(fm, x32, uv32; f32...)) == ComplexF32
        # the same model compiled again generates the same code, so it reuses the compiled methods
        fm2 = dict_to_model(m, free)
        @test typeof(fm2.resolver.fill!) == typeof(fm.resolver.fill!)
        @test typeof(fm2.resolver.fn) == typeof(fm.resolver.fn)
    end

    @testset "simulated through V1295 Aql's uv points, fitted back" begin
        src = joinpath(@__DIR__, "..", "demos", "data", "V1295_Aql.oifits")
        truth = Dict{String,Any}("1,ud" => 0.3, "2,f" => 0.3, "2,x" => "orbit", "2,y" => "orbit",
                                 "2,orb P" => 30.0, "2,orb MJD0" => 58640.0, "2,orb e" => 0.4,
                                 "2,orb omega" => 60.0, "2,orb OMEGA" => 120.0, "2,orb incl" => 50.0,
                                 "2,orb a" => 4.0)
        out = joinpath(mktempdir(), "orbit.oifits")
        simulate_from_oifits(src, out; mode = "copy_errors", noise = false,
                             flat_model = dict_to_model(truth, String[]), flat_params = Float64[])
        d = readoifits(out; warn = false, verbose = false, T = Float64)[1, 1]
        start = merge(copy(truth), Dict{String,Any}("2,orb e" => 0.35, "2,orb omega" => 70.0,
                                                    "2,orb OMEGA" => 115.0, "2,orb incl" => 45.0,
                                                    "2,orb a" => 3.8, "2,orb P" => 30.2))
        free = ["2,orb e", "2,orb omega", "2,orb OMEGA", "2,orb incl", "2,orb a", "2,orb P"]
        r = fit_model_lsqfit(start, free, d; lb = Dict("2,orb e" => 0.0), ub = Dict("2,orb e" => 0.99))
        p = Dict(zip(r.list_free_params, r.x_opt))
        # noiseless, but written in Float32 and fitted to LM's tolerance in units of the noise
        for k in free
            @test p[k] ≈ truth[k] rtol = 1e-3
        end
        @test r.converged
    end
end
