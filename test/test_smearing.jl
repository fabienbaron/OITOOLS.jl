# Bandwidth smearing (`"smear" => n` in a model): the observables averaged over each spectral
# channel, judged against averages written out here from the closed forms, not against the
# implementation. AX Cir (PIONIER, R ≈ 18, B_max 51 Mλ) is the data: its smearing radius,
# Rλ/2B_max ≈ 36 mas, is where the companion search of GUI_CompanionSearch.md works.
using Statistics

@testset "bandwidth smearing" begin
    MAS = 2.0626480624709636e8
    ax  = readoifits(joinpath(@__DIR__, "..", "demos", "data", "AX_Cir.oifits");
                     T = Float64, warn = false, verbose = false)[1, 1]
    λ, Δλ = ax.uv_lam, ax.uv_dlam

    # A channel as the model samples it: flat in wavenumber between 1/(λ+Δλ/2) and
    # 1/(λ-Δλ/2). `M` points of a fine midpoint rule, standing in for the integral.
    sub(i, M) = (lo = 1 / (λ[i] + Δλ[i] / 2); hi = 1 / (λ[i] - Δλ[i] / 2);
                 [lo + (k - 0.5) / M * (hi - lo) for k in 1:M])
    # The visibility of a binary at uv point i, sub-wavenumber σ: two points, fluxes 1-f and f.
    binary(i, σ, f, x, y) = (1 - f) + f * cis(-2π * (ax.uv[1, i] * x + ax.uv[2, i] * y) / MAS *
                                              λ[i] * σ)
    model(extra...) = dict_to_model(Dict{String,Any}("p,f" => "1 - \$c,f", "c,f" => 0.3,
                                                     "c,x" => 6.0, "c,y" => -28.0, extra...),
                                    String[])

    @testset "off unless asked for, and then exactly as before" begin
        for d in (Dict{String,Any}(), Dict{String,Any}("smear" => 0), Dict{String,Any}("smear" => 1))
            m = model(d...)
            @test model_to_chi2(m, Float64[], ax) === model_to_chi2(model(), Float64[], ax)
            @test model_to_obs(m, Float64[], ax).v2 == model_to_obs(model(), Float64[], ax).v2
        end
        @test OITOOLS._smeared(model("smear" => 5), ax)
        @test occursin("smear: 5", sprint(show, model("smear" => 5)))
    end

    @testset "a binary's V² against the channel average of |V|²" begin
        for (x, y) in ((6.0, -28.0), (20.0, 30.0))
            obs = model_to_obs(model("smear" => 64, "c,x" => x, "c,y" => y), Float64[], ax).v2
            truth = [mean(abs2(binary(i, σ, 0.3, x, y)) for σ in sub(i, 4000))
                     for i in ax.indx_v2]
            @test maximum(abs.(obs .- truth)) < 2e-4
            plain = model_to_obs(model("c,x" => x, "c,y" => y), Float64[], ax).v2
            @test maximum(abs.(plain .- truth)) > 1e-2       # it matters at this separation
        end
        # The closed form CANDID uses: the cross term carries sinc(B·s·Δσ).
        i = ax.indx_v2[argmax(hypot.(ax.uv[1, ax.indx_v2], ax.uv[2, ax.indx_v2]))]
        lo, hi = 1 / (λ[i] + Δλ[i] / 2), 1 / (λ[i] - Δλ[i] / 2)
        path = (ax.uv[1, i] * 6.0 + ax.uv[2, i] * -28.0) / MAS * λ[i]   # cycles per unit σ
        cf = 0.7^2 + 0.3^2 + 2 * 0.3 * 0.7 * cos(2π * path * (lo + hi) / 2) * sinc(path * (hi - lo))
        j = findfirst(==(i), ax.indx_v2)
        @test model_to_obs(model("smear" => 64), Float64[], ax).v2[j] ≈ cf atol = 2e-4
    end

    @testset "closure phases: the channel average of the triple product" begin
        m = model("smear" => 64, "c,x" => 20.0, "c,y" => 30.0)
        obs = model_to_obs(m, Float64[], ax).t3phi
        t3(i1, i2, i3) = mean(binary(i1, σ, 0.3, 20.0, 30.0) * binary(i2, σ, 0.3, 20.0, 30.0) *
                              binary(i3, σ, 0.3, 20.0, 30.0) for σ in sub(i1, 4000))
        truth = [angle(t3(ax.indx_t3_1[k], ax.indx_t3_2[k], ax.indx_t3_3[k])) * 180 / π
                 for k in eachindex(ax.indx_t3_1)]
        d = mod.(obs .- truth .+ 180, 360) .- 180
        @test maximum(abs.(d)) < 0.05
        # Smearing washes a binary's closure phases out.
        plain = model_to_obs(model("c,x" => 20.0, "c,y" => 30.0), Float64[], ax).t3phi
        @test mean(abs.(obs)) < mean(abs.(plain))
    end

    @testset "a uniform disc's null fills in" begin
        # A disc whose first null falls inside AX Cir's baselines: V² there is ~0 unsmeared
        # and the channel average of the Airy function is not.
        D = 3.8317059702075125 / π / 45e6 * MAS        # first null (J₁ zero) at 45 Mλ
        ud(e...) = dict_to_model(Dict{String,Any}("s,ud" => D, e...), String[])
        Vud(q) = (t = π * q * D / MAS; t == 0 ? 1.0 : 2 * OITOOLS.besselj1(t) / t)
        obs = model_to_obs(ud("smear" => 64), Float64[], ax).v2
        truth = [mean(Vud(hypot(ax.uv[1, i], ax.uv[2, i]) * λ[i] * σ)^2 for σ in sub(i, 4000))
                 for i in ax.indx_v2]
        @test maximum(abs.(obs .- truth)) < 1e-5
        # At the null itself, with AX Cir's channel width: |V|² is 0 at the centre and the
        # channel average is not.
        q = [45e6; 0.0;;]
        Vs = eval_model_sub(ud("smear" => 64), Float64[], q; wl = λ[1], dlam = Δλ[1])
        @test abs2(eval_model(ud(), Float64[], q)[1]) < 1e-12
        lo, hi = 1 / (λ[1] + Δλ[1] / 2), 1 / (λ[1] - Δλ[1] / 2)
        filled = mean(Vud(45e6 * λ[1] * (lo + (k - 0.5) / 4000 * (hi - lo)))^2 for k in 1:4000)
        @test filled > 1e-4
        @test mean(abs2, Vs) ≈ filled rtol = 1e-3
    end

    @testset "χ², its gradient and the residuals carry it through" begin
        m = dict_to_model(Dict{String,Any}("p,f" => "1 - \$c,f", "c,f" => 0.3, "c,x" => 6.0,
                                           "c,y" => -28.0, "smear" => 8),
                          ["c,f", "c,x", "c,y"])
        x = [0.3, 6.0, -28.0]
        c, g = model_to_chi2_fg(m, x, ax)
        @test c ≈ model_to_chi2(m, x, ax)
        h = 1e-6
        fd = [(model_to_chi2(m, x .+ h .* (1:3 .== j), ax) -
               model_to_chi2(m, x .- h .* (1:3 .== j), ax)) / 2h for j in 1:3]
        @test g ≈ fd rtol = 1e-4
        r, J = OITOOLS.residuals_flat_jac(m, x, ax)
        @test sum(abs2, r) ≈ c rtol = 1e-10
        @test size(J) == (length(r), 3)
        @test 2 .* (J' * r) ≈ g rtol = 1e-6
    end

    @testset "a fit of smeared data recovers the companion" begin
        truth = Dict{String,Any}("p,f" => "1 - \$c,f", "c,f" => 0.05, "c,x" => 20.0,
                                 "c,y" => 30.0, "smear" => 8)
        out = joinpath(mktempdir(), "sim.oifits")
        simulate_from_oifits(joinpath(@__DIR__, "..", "demos", "data", "AX_Cir.oifits"), out;
                             flat_model = dict_to_model(truth, String[]),
                             flat_params = Float64[], seed = 3)
        data = readoifits(out; T = Float64, warn = false, verbose = false)[1, 1]
        start = merge(truth, Dict{String,Any}("c,f" => 0.04, "c,x" => 19.8, "c,y" => 30.2))
        res = fit_model(start, ["c,f", "c,x", "c,y"], data; method = :LD_MMA)
        @test 0.5 < res.chi2r < 1.5
        @test abs(res.x_opt[1] - 0.05) < 0.005
        @test hypot(res.x_opt[2] - 20.0, res.x_opt[3] - 30.0) < 0.2
        # Fitted without smearing, the companion comes back fainter: its fringes were washed out.
        nosm = fit_model(delete!(copy(start), "smear"), ["c,f", "c,x", "c,y"], data;
                         method = :LD_MMA)
        @test nosm.x_opt[1] < res.x_opt[1] && nosm.chi2r > res.chi2r
    end

    @testset "an image component keeps its channel's plane" begin
        rot = joinpath(@__DIR__, "data", "rotir")
        f = sort(filter(endswith(".oifits"), readdir(rot; join = true)))[1]
        d = readoifits(f; polychromatic = true, warn = false, verbose = false, T = Float64)[1, 1]
        m = dict_to_model(Dict{String,Any}("sky,image" => joinpath(rot, "zetand_ref_cube.fits"),
                                           "smear" => 4), String[])
        @test isfinite(model_to_chi2(m, Float64[], d))
    end

    @testset "coherent: the visibility averaged, then the observables (CANDID, ABCD)" begin
        obs = model_to_obs(model("smear" => 64, "smear_mode" => "coherent"), Float64[], ax).v2
        truth = [abs2(mean(binary(i, σ, 0.3, 6.0, -28.0) for σ in sub(i, 4000)))
                 for i in ax.indx_v2]
        @test maximum(abs.(obs .- truth)) < 2e-4
        # CANDID's closed form: the sinc on the companion's term, squared with it.
        i = ax.indx_v2[argmax(hypot.(ax.uv[1, ax.indx_v2], ax.uv[2, ax.indx_v2]))]
        lo, hi = 1 / (λ[i] + Δλ[i] / 2), 1 / (λ[i] - Δλ[i] / 2)
        path = (ax.uv[1, i] * 6.0 + ax.uv[2, i] * -28.0) / MAS * λ[i]
        g = sinc(path * (hi - lo))
        cf = 0.7^2 + (0.3 * g)^2 + 2 * 0.3 * 0.7 * cos(2π * path * (lo + hi) / 2) * g
        j = findfirst(==(i), ax.indx_v2)
        @test obs[j] ≈ cf atol = 2e-4
        # The two conventions differ by f²(1 − sinc²) on each point.
        inc = model_to_obs(model("smear" => 64), Float64[], ax).v2
        @test inc[j] - obs[j] ≈ 0.3^2 * (1 - g^2) atol = 2e-4
        cp = model_to_obs(model("smear" => 64, "smear_mode" => "coherent"), Float64[], ax).t3phi
        t3(i1, i2, i3) = (V(k) = mean(binary(k, σ, 0.3, 6.0, -28.0) for σ in sub(i1, 4000));
                          V(i1) * V(i2) * V(i3))
        tr = [angle(t3(ax.indx_t3_1[k], ax.indx_t3_2[k], ax.indx_t3_3[k])) * 180 / π
              for k in eachindex(ax.indx_t3_1)]
        @test maximum(abs.(mod.(cp .- tr .+ 180, 360) .- 180)) < 0.05
        m = dict_to_model(Dict{String,Any}("p,f" => "1 - \$c,f", "c,f" => 0.3, "c,x" => 6.0,
                                           "c,y" => -28.0, "smear" => 8,
                                           "smear_mode" => "coherent"), ["c,x"])
        c, gr = model_to_chi2_fg(m, [6.0], ax)
        h = 1e-6
        @test gr[1] ≈ (model_to_chi2(m, [6.0 + h], ax) - model_to_chi2(m, [6.0 - h], ax)) / 2h rtol = 1e-4
    end

    @testset "auto: CANDID's rule, at the model's extent" begin
        # A companion at 28.6 mas on AX Cir's longest baseline: the phase changes by
        # |uv|·ρ·Δλ/λ cycles across a channel, and four samples per cycle are taken.
        m = model("smear" => "auto")
        ρ = hypot(6.0, -28.0)
        @test m.smear_radius ≈ ρ
        cyc = maximum(hypot(ax.uv[1, i], ax.uv[2, i]) * ρ / MAS * Δλ[i] / λ[i]
                      for i in axes(ax.uv, 2))
        n = smear_samples(m, ax)
        @test n == max(3, ceil(Int, 4cyc))
        @test model_to_chi2(m, Float64[], ax) ≈ model_to_chi2(model("smear" => n), Float64[], ax)
        # A centred disc's extent is its radius.
        @test dict_to_model(Dict{String,Any}("s,ud" => 3.0, "smear" => "auto"),
                            String[]).smear_radius ≈ 1.5
        @test smear_samples(dict_to_model(Dict{String,Any}("s,ud" => 0.5, "smear" => "auto"),
                                          String[]), ax) == 3
        @test occursin("auto", sprint(show, m))
    end

    @testset "what is refused" begin
        @test_throws ErrorException dict_to_model(Dict{String,Any}("s,ud" => 1.0, "smear" => 2.5),
                                                  String[])
        @test_throws ErrorException dict_to_model(Dict{String,Any}("s,ud" => 1.0, "smear" => 3),
                                                  ["smear"])
        @test_logs dict_to_model(Dict{String,Any}("s,ud" => 1.0, "smear" => 3), String[])
        @test_throws ErrorException dict_to_model(Dict{String,Any}("s,ud" => 1.0, "smear" => 3,
                                                                   "smear_mode" => "partial"),
                                                  String[])
        @test_throws ErrorException dict_to_model(Dict{String,Any}("s,ud" => 1.0,
                                                                   "smear" => "sometimes"),
                                                  String[])
        m = dict_to_model(Dict{String,Any}("s,ud" => 1.0, "smear" => 3), String[])
        @test_throws ArgumentError OITOOLS._refuse_smear(m, "SPARCO image + model")
    end
end
