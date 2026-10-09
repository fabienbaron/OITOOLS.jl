# The companion search, on PMOIRED's example 3 (AX Cir, PIONIER): tolerances from PMOIRED's
# printed results (GUI_CompanionSearch.md, phase C), not from this code's own output.
#   UD alone: χ²ᵣ 0.974, ud 0.9309 mas.  UD + companion: ud 0.8235 mas, companion 0.00875 of
#   the total flux at x = +6.228, y = −28.487 mas; ≈ 8σ (CANDID agrees).
@testset "companion search" begin
    ax = readoifits(joinpath(@__DIR__, "..", "demos", "data", "AX_Cir.oifits");
                    T = Float64, warn = false, verbose = false)[1, 1]
    w  = [1.0, 0.0, 1.0]                                  # V² and T3φ, as PMOIRED fitted
    ud = fit_model(Dict{String,Any}("star,ud" => 0.9), ["star,ud"], ax; method = :LD_MMA,
                   weights = w)
    @test ud.x_opt[1] ≈ 0.9309 atol = 0.005
    @test ud.chi2r ≈ 0.974 atol = 0.01
    md = Dict{String,Any}("star,ud" => ud.x_opt[1], "star,f" => "1 - \$comp,f",
                          "comp,f" => 0.01, "comp,x" => 0.0, "comp,y" => 0.0,
                          "smear" => "auto")
    truth = (6.228, -28.487)

    @testset "significance: CANDID's convention" begin
        # With one degree of freedom the χ² tail is a Gaussian's: n_σ = √(χ² ratio).
        @test detection_sigma(9.0, 1.0, 1) ≈ 3.0 rtol = 1e-10
        @test detection_sigma(25.0, 1.0, 1) ≈ 5.0 rtol = 1e-10
        # Past Float64's tail the asymptotic form takes over, continuing upward.
        s = [detection_sigma(r, 1.0, 1500) for r in (1.2, 1.5, 2.0, 3.0)]
        @test issorted(s) && all(isfinite, s)
    end

    # The disc and the companion's flux free at each node: test 1 of the plan.
    g = companion_search(md, ax; weights = w, free = ["comp,f", "star,ud"])
    step = g.x[2] - g.x[1]

    @testset "AX Cir: found where PMOIRED and CANDID find it" begin
        @test g.path === :general
        @test hypot(g.best.x - truth[1], g.best.y - truth[2]) <= step
        @test g.best.f ≈ 0.00875 rtol = 0.15
        @test g.best.params[2] ≈ 0.8235 atol = 0.03         # the disc shrinks with the companion
        @test g.nsigma > 5
        @test g.chi2r_null ≈ 0.974 atol = 0.01
        @test g.fit isa FitResult
        @test g.fit.chi2r ≈ g.best.chi2r rtol = 1e-6
        @test step ≈ 2.02 atol = 0.02                       # λ/(2 B_max)
        @test maximum(abs, g.x) ≤ 36.6 && maximum(abs, g.x) > 34   # R λ/(2 B_max)
        @test occursin("general path", sprint(show, g))
    end

    @testset "the mirror position is decisively worse" begin
        i = findfirst(==(g.best.x), g.x); j = findfirst(==(g.best.y), g.y)
        mirror = g.chi2r[length(g.x) + 1 - i, length(g.y) + 1 - j]
        @test (mirror - g.best.chi2r) * g.ndof > 25          # Δχ² far beyond any 1σ wobble
    end

    @testset "the fast path agrees with the general one where both apply" begin
        md2 = merge(md, Dict{String,Any}("star,ud" => g.best.params[2]))
        fast = companion_search(md2, ax; weights = w, rmin = 26, rmax = 32)
        slow = companion_search(md2, ax; weights = w, rmin = 26, rmax = 32, path = :general)
        @test fast.path === :fast && slow.path === :general
        ok = .!isnan.(fast.chi2r)
        @test ok == .!isnan.(slow.chi2r)
        @test fast.chi2r[ok] ≈ slow.chi2r[ok] rtol = 1e-5
        @test companion_flux_map(fast)[ok] ≈ companion_flux_map(slow)[ok] atol = 1e-5
        # Halley reaches BraDi's answers, falling back where it must.
        hal = companion_search(md2, ax; weights = w, rmin = 26, rmax = 32, solver = :halley)
        @test hal.chi2r[ok] ≈ fast.chi2r[ok] rtol = 1e-6
    end

    @testset "Halley's derivatives are the closed forms' (against finite differences)" begin
        O = OITOOLS
        for mode in ("incoherent", "coherent"),
            ww in ([1.0, 0.0, 1.0], [1.0, 1.0, 1.0, 1.0, 1.0])
            mdd = merge(md, Dict{String,Any}("smear_mode" => mode))
            fp = O._fast_problem(mdd, ax, "comp", O._pad_weights(ww))
            ws = O._fast_coefficients!(O.FastWork(fp), fp, 6.0, -28.5)
            χ(f) = O._fast_chi2(f, ws, fp)
            f, h = 0.012, 1e-4
            c, c1, c2, c3 = O._fast_chi2_d3(f, ws, fp)
            @test c ≈ χ(f)
            @test c1 ≈ (χ(f + h) - χ(f - h)) / 2h rtol = 1e-6
            @test c2 ≈ (χ(f + h) - 2χ(f) + χ(f - h)) / h^2 rtol = 1e-4
            @test c3 ≈ (χ(f + 2h) - 2χ(f + h) + 2χ(f - h) - χ(f - 2h)) / (2h^3) rtol = 1e-2
        end
    end

    @testset "with an image as the primary" begin
        # The disc fitted with the companion, rendered by pixel area into a FITS image in the
        # shared convention and loaded as the primary: the motivating case end to end.
        D, nx, ps = g.best.params[2], 96, 0.02
        c = nx ÷ 2 + 1
        img = zeros(nx, nx)
        for jj in 1:nx, ii in 1:nx, b in 1:6, a in 1:6
            x = (ii - c + (a - 3.5) / 6) * ps; y = (jj - c + (b - 3.5) / 6) * ps
            x^2 + y^2 <= (D / 2)^2 && (img[ii, jj] += 1)
        end
        p = joinpath(mktempdir(), "primary.fits"); writefits(img ./ sum(img), p; pixsize = ps)
        mi = Dict{String,Any}("star,image" => p, "star,f" => "1 - \$comp,f", "comp,f" => 0.01,
                              "comp,x" => 0.0, "comp,y" => 0.0, "smear" => "auto")
        s = companion_search(mi, ax; weights = w)
        @test s.path === :fast
        @test (s.best.x, s.best.y) == (g.best.x, g.best.y)
        # Against the same disc as an analytic component, both fixed: the image path is the
        # analytic one to the rendering's accuracy. (Not against `g`, whose null model REFITS
        # the disc to 0.93 mas; a fixed primary's null cannot be refitted, so its significance
        # answers a different question.)
        a = companion_search(merge(md, Dict{String,Any}("star,ud" => D)), ax; weights = w)
        @test (a.best.x, a.best.y) == (s.best.x, s.best.y)
        @test s.best.f ≈ a.best.f rtol = 0.01
        @test s.nsigma ≈ a.nsigma rtol = 0.02
        @test s.best.f ≈ g.best.f rtol = 0.05
        # The primary's own pixels are never searched.
        s2 = companion_search(mi, ax; weights = w, rmin = 0.0, rmax = 6.0, step = 0.5)
        @test isnan(s2.chi2r[findfirst(==(0.0), s2.x), findfirst(==(0.0), s2.y)])
        @test count(!isnan, s2.chi2r) > 0
    end

    @testset "no companion, no detection" begin
        out = joinpath(mktempdir(), "ud_only.oifits")
        simulate_from_oifits(joinpath(@__DIR__, "..", "demos", "data", "AX_Cir.oifits"), out;
                             flat_model = dict_to_model(Dict{String,Any}("star,ud" => 0.93),
                                                        String[]),
                             flat_params = Float64[], seed = 5)
        nul = readoifits(out; T = Float64, warn = false, verbose = false)[1, 1]
        s = companion_search(merge(md, Dict{String,Any}("star,ud" => 0.93)), nul; weights = w)
        @test s.nsigma < 3
    end

    @testset "Float32 data are searched in Float32" begin
        ax32 = readoifits(joinpath(@__DIR__, "..", "demos", "data", "AX_Cir.oifits");
                          warn = false, verbose = false)[1, 1]
        @test eltype(ax32.uv) === Float32
        fp = OITOOLS._fast_problem(md, ax32, "comp", OITOOLS._pad_weights(w))
        @test fp isa OITOOLS.FastProblem{Float32}
        ws = OITOOLS._fast_coefficients!(OITOOLS.FastWork(fp), fp, 6.0, -28.5)
        @test OITOOLS._fast_chi2(0.008, ws, fp) isa Float32
        @test OITOOLS._fast_chi2_d3(0.008, ws, fp) isa NTuple{4,Float32}
        s32 = companion_search(md, ax32; weights = w, rmin = 26, rmax = 32)
        s64 = companion_search(md, ax; weights = w, rmin = 26, rmax = 32)
        # The same node (its coordinates differ in the 7th digit: the step is computed from
        # Float32 uv). A Float32 χ² summed over ~1500 terms is noisy at ~6e-5 relative, which
        # moves the flux minimum by ~1e-4 here, against a statistical σ of ~4e-4.
        @test hypot(s32.best.x - s64.best.x, s32.best.y - s64.best.y) < 1e-5
        @test s32.best.f ≈ s64.best.f rtol = 0.02
        @test s32.nsigma ≈ s64.nsigma rtol = 0.01
        # The smeared gradient follows the data's precision too.
        m = dict_to_model(merge(md, Dict{String,Any}("smear" => 4)), ["comp,f"])
        c, g = model_to_chi2_fg(m, Float32[0.009], ax32; weights = w)
        @test c isa Float32 && eltype(g) === Float32
    end

    # The hook a window draws the map from while the search runs: called once, before any node
    # is fitted, with the arrays the search then fills and returns. And the grid on its own,
    # which is what a panel shows before running.
    @testset "on_grid hands over the arrays being filled; companion_grid is the same grid" begin
        calls = Any[]
        s = companion_search(md, ax; weights = w, rmin = 26, rmax = 32,
                             on_grid = g -> push!(calls, (g, count(isfinite, g.chi2r))))
        @test length(calls) == 1
        g, finite_then = only(calls)
        @test finite_then == 0                                   # nothing fitted yet
        @test g.chi2r === s.chi2r && g.params === s.params && g.nnodes == s.nnodes
        @test g.x == s.x && g.y == s.y
        # The order nodes are visited in is a permutation, so every node is fitted once.
        @test all(sort(OITOOLS._spread_order(n)) == 1:n for n in (1, 2, 3, 10, 1030, 31112))
        cg = companion_grid(md, ax; rmin = 26, rmax = 32)
        @test cg.x == s.x && length(cg.nodes) == s.nnodes && cg.rmin == 26 && cg.rmax == 32
    end

    @testset "what is refused" begin
        @test_throws ArgumentError companion_search(Dict{String,Any}("star,ud" => 0.9), ax)
        @test_throws ArgumentError companion_search(md, ax; free = ["star,ud"])
        @test_throws ArgumentError companion_search(md, ax; free = ["comp,f", "comp,x"])
        @test_throws ArgumentError companion_search(md, ax; weights = [1, 1, 1, 0, 0, 1])
        nonlin = merge(md, Dict{String,Any}("star,f" => "1 - \$comp,f^2"))
        @test_throws ArgumentError companion_search(nonlin, ax)
    end
end
