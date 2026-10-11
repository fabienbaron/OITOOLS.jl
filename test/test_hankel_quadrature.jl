# A radial profile's quadrature (PLAN_HankelQuadrature.md): written profiles against closed
# forms and against quadgk given their breakpoints, what the probe finds in a profile, the node
# count, and the derivatives, a step's moving position included, in Float64 and Float32.
using QuadGK, SpecialFunctions, ForwardDiff

@testset "profile quadrature" begin
    O = OITOOLS
    MAS = O.MAS2RAD
    uvof(B) = vcat(permutedims(B .* MAS), zeros(1, length(B)))   # cycles/mas → cycles/rad, on u
                                                                 # (MAS2RAD is mas per radian)
    Bset(L, D) = collect(range(L / 64, L, length = 64)) ./ D      # B_max·D = L lobes
    function quadV(f, pts, B; n = 0)
        N = quadgk(r -> f(r) * r, pts...; rtol = 1e-13)[1]
        [quadgk(r -> f(r) * besselj(n, 2π * b * r) * r, pts...; rtol = 1e-13, maxevals = 10^7)[1] / N
         for b in B]
    end
    vis(md, B; T = Float64) = eval_model(dict_to_model(md, String[]; T), T[], T.(uvof(B)))
    tol(T) = T == Float64 ? 1e-8 : 1e-6

    @testset "the five laws, written as profiles, against their closed forms" begin
        D = 4.0
        for (prof, pars, exact) in (
                ("1.0", Dict{String,Any}(), ζ -> O._ud_v(ζ)),
                ("1 - \$u*(1 - \$MU)", Dict{String,Any}("s,u" => 0.5), ζ -> O._ldlin_v(ζ, 0.5)),
                ("1 - \$qa*(1 - \$MU) - \$qb*(1 - \$MU)^2", Dict{String,Any}("s,qa" => 0.3, "s,qb" => 0.2),
                 ζ -> O._ldquad_v(ζ, 0.3, 0.2)),
                ("\$MU^\$alpha", Dict{String,Any}("s,alpha" => 0.14), ζ -> O._ldpow_v(ζ, 0.14)),
                ("1 - 0.7127*(1 - \$MU^0.5) + 0.0452*(1 - \$MU) + 0.2643*(1 - \$MU^1.5) - 0.1311*(1 - \$MU^2)",
                 Dict{String,Any}(), ζ -> O._ldclaret4_v(ζ, 0.7127, -0.0452, -0.2643, 0.1311))),
            L in (3, 17, 70)
            B = Bset(L, D)
            md = merge(Dict{String,Any}("s,diam" => D, "s,profile" => prof), pars)
            @test maximum(abs.(vis(md, B) .- exact.(π * D .* B))) < 1e-8
        end
        # Float32, against the Float64 closed form, at its own floor
        md = Dict{String,Any}("s,diam" => D, "s,profile" => "\$MU^0.14")
        @test maximum(abs.(vis(md, Bset(3, D); T = Float32) .- O._ldpow_v.(π * D .* Bset(3, D), 0.14))) < 1e-6
    end

    @testset "a ring (r_min > 0) and a Gaussian, against their closed forms" begin
        B = Bset(17, 6.0)
        H(r) = r .* besselj1.(2π .* B .* r) ./ (2π .* B)
        ring = vis(Dict{String,Any}("r,profile" => "1.0", "r,diamin" => 4.0, "r,udout" => 6.0), B)
        @test maximum(abs.(ring .- (H(3.0) .- H(2.0)) ./ ((9 - 4) / 2))) < 1e-8
        # exp(−r²/s²) on r ≤ 6s: the truncation is e^-36
        gauss = vis(Dict{String,Any}("g,profile" => "exp(-(\$R/0.5)^2)", "g,udout" => 6.0), B)
        @test maximum(abs.(gauss .- exp.(-(π * 0.5 .* B) .^ 2))) < 1e-8
    end

    @testset "azimuthal orders 1–3 on a limb-darkened disc, against quadgk" begin
        D = 4.0; B = Bset(17, D)                      # baselines on u: ψ = 0
        f(r) = sqrt(max(0, 1 - (r / 2)^2))^0.14
        H0 = quadV(f, [0.0, 2.0], B)
        for n in 1:3
            md = Dict{String,Any}("s,diam" => D, "s,profile" => "\$MU^0.14",
                                  "s,az amp$n" => 0.3, "s,az projang$n" => 20.0)
            ex = H0 .+ 0.3 .* (-im)^n .* quadV(f, [0.0, 2.0], B; n) .* cos(n * (deg2rad(20) - π / 2))
            @test maximum(abs.(vis(md, B) .- ex)) < 1e-8
        end
    end

    @testset "the node count follows the data and the profile" begin
        md = Dict{String,Any}("s,diam" => 4.0, "s,profile" => "\$MU^0.14")
        for (L, n) in ((3, 32), (70, 104))
            m = dict_to_model(md, String[]); eval_model(m, Float64[], uvof(Bset(L, 4.0)))
            @test collect(keys(only(m.components).kcache.nodes)) == [(n, n)]
        end
        m = dict_to_model(merge(md, Dict{String,Any}("s,nr" => 50)), String[])
        eval_model(m, Float64[], uvof(Bset(3, 4.0)))
        @test first(only(keys(only(m.components).kcache.nodes))) == 50
        @test_throws ErrorException dict_to_model(merge(md, Dict{String,Any}("s,nr" => 1)), String[])
        # a rim sharper than the lobe rule allows for (the α double sigmoid template at its
        # seeds) gets more nodes, and is then right
        md = Dict{String,Any}("d,profile" => "(\$a + (1-\$a)/(1+exp(-(\$R-\$rin)/\$win))) * 1/(1+exp((\$R-\$rout)/\$wout))",
                              "d,a" => 0.17, "d,rin" => 1.25, "d,rout" => 1.24, "d,win" => 0.10,
                              "d,wout" => 0.50, "d,udout" => 8.0)
        sig(x) = 1 / (1 + exp(-x))
        f(r) = (0.17 + 0.83 * sig((r - 1.25) / 0.10)) * sig(-(r - 1.24) / 0.50)
        B = Bset(3, 8.0)
        m = dict_to_model(md, String[]); V = eval_model(m, Float64[], uvof(B))
        @test maximum(abs.(V .- quadV(f, [0.0, 4.0], B))) < 1e-8
        @test only(m.components).kcache.mult > 1
    end

    @testset "steps, kinks and a singular centre: found, and integrated to $(T)'s floor" for T in (Float64, Float32)
        for (name, md, f, pts, found, rmax) in (
                ("step ring", Dict{String,Any}("d,profile" => "(\$R >= 3.0) * (\$R <= 5.0) * 1.0", "d,udout" => 10.0),
                 r -> (3 <= r <= 5) * 1.0, [0.0, 3.0, 5.0], [:jump], 5.0),
                ("docs: (R > 0.25)·R^-0.5", Dict{String,Any}("d,profile" => "(\$R > 0.25) * \$R^(-0.5)", "d,udout" => 4.0),
                 r -> (r > 0.25) * r^-0.5, [0.0, 0.25, 2.0], [:jump], 2.0),
                ("docs: two kinks", Dict{String,Any}("d,profile" => "max(0, 1 - (2*(\$R - 2.0))^2)", "d,udout" => 6.0),
                 r -> max(0, 1 - (2 * (r - 2))^2), [0.0, 1.5, 2.5, 3.0], [:kink, :kink], 3.0),
                ("singular power law", Dict{String,Any}("d,profile" => "\$R^(-1.5)", "d,udout" => 8.0),
                 r -> r^-1.5, [0.0, 4.0], Symbol[], 4.0)),
            L in (3, 17)
            B = Bset(L, 2rmax)
            m = dict_to_model(md, String[]; T)
            V = eval_model(m, T[], T.(uvof(B)))
            @test maximum(abs.(V .- quadV(f, pts, B))) < tol(T)
            lay = only(m.components).kcache.layout
            @test filter(!=(:transition), lay.kinds) == found
            @test maximum(abs.(lay.breaks .- pts[2:end-1] ./ rmax); init = 0.0) < 1e-6
            @test lay.sing == (name == "singular power law")
        end
        # the same step with an azimuthal mode (order 1, odd: the case GL in t alone fails)
        B = Bset(3, 10.0)
        md = Dict{String,Any}("d,profile" => "(\$R >= 3.0) * (\$R <= 5.0) * 1.0", "d,udout" => 10.0,
                              "d,az amp1" => 0.4, "d,az projang1" => 37.0)
        ex = quadV(r -> 1.0, [3.0, 5.0], B) .+ 0.4 .* (-im) .* quadV(r -> 1.0, [3.0, 5.0], B; n = 1) .*
             cos(deg2rad(37) - π / 2)
        @test maximum(abs.(vis(md, B; T) .- ex)) < tol(T)
    end

    @testset "nothing found in smooth profiles" begin
        for md in (Dict{String,Any}("s,diam" => 4.0, "s,profile" => "\$MU^0.14"),
                   Dict{String,Any}("s,diam" => 4.0, "s,profile" => "1 - 0.5*(1 - \$MU)"),
                   Dict{String,Any}("s,udout" => 6.0, "s,profile" => "exp(-(\$R/0.5)^2)"),
                   Dict{String,Any}("s,udout" => 5.0, "s,profile" => "1 - ((\$R - 1.5)/2)^2"),
                   Dict{String,Any}("s,udout" => 8.0, "s,profile" => "1/(1+exp(-(\$R-1.26)/0.26)) * 1/(1+exp((\$R-1.25)/0.43))"),
                   Dict{String,Any}("s,udout" => 8.0, "s,profile" => "exp(-(\$R^2 + 4)/(2*0.09)) * besseli(0, \$R*2/0.09)"))
            m = dict_to_model(md, String[]); eval_model(m, Float64[], uvof(Bset(3, 4.0)))
            lay = only(m.components).kcache.layout
            @test isempty(lay.breaks) && !lay.sing
        end
    end

    @testset "the detector's edge cases" for T in (Float64, Float32)
        h = 0.99 / 4096                                  # the probe's spacing in s
        for (name, md, f, pts) in (
                ("a step on a probe point", Dict{String,Any}("e,profile" => "(\$R >= $(2048h)) * 1.0 + 0.5", "e,udout" => 2.0),
                 r -> (r >= 2048h) * 1.0 + 0.5, [0.0, 2048h, 1.0]),
                ("a step inside the rim", Dict{String,Any}("e,profile" => "(\$R < 0.995) * 1.0 + 0.2", "e,udout" => 2.0),
                 r -> (r < 0.995) * 1.0 + 0.2, [0.0, 0.995, 1.0]),
                ("two steps one probe spacing apart", Dict{String,Any}("e,profile" => "1.0 + (\$R >= 0.5) + (\$R >= $(0.5 + h))", "e,udout" => 2.0),
                 r -> 1.0 + (r >= 0.5) + (r >= 0.5 + h), [0.0, 0.5, 0.5 + h, 1.0]),
                ("a step next to r_min", Dict{String,Any}("e,profile" => "(\$R >= 0.3005) * 1.0 + 0.1", "e,diamin" => 0.6, "e,udout" => 2.0),
                 r -> (r >= 0.3005) * 1.0 + 0.1, [0.3, 0.3005, 1.0]),
                ("a sigmoid narrower than the probe spacing", Dict{String,Any}("e,profile" => "1/(1+exp(-(\$R - 0.4)/1e-5))", "e,udout" => 2.0),
                 r -> 1 / (1 + exp(-(r - 0.4) / 1e-5)), [0.0, 0.4, 1.0]))
            B = Bset(3, 2.0)
            @test maximum(abs.(vis(md, B; T) .- quadV(f, pts, B))) < (T == Float64 ? 1e-8 : 2e-6)
        end
    end

    @testset "a step at a fitted position: ∂V/∂rin through its boundary term" for T in (Float64, Float32)
        B = Bset(3, 10.0); uv = T.(uvof(B))
        m = dict_to_model(Dict{String,Any}("q,profile" => "(\$R >= \$rin) * 1.0", "q,rin" => 3.0,
                                           "q,udout" => 10.0), ["q,rin"]; T)
        for p in (3.0, 3.3)                      # 3.3: the step moved, and is followed
            N = (25 - p^2) / 2
            V = [quadgk(r -> besselj0(2π * b * r) * r, p, 5; rtol = 1e-13)[1] for b in B] ./ N
            dV = p .* (V .- besselj0.(2π .* B .* p)) ./ N
            t = T == Float64 ? 1e-8 : 1e-5
            @test maximum(abs.(eval_model(m, T[p], uv) .- V)) < t
            @test maximum(abs.(eval_model_jac(m, T[p], uv)[2][:, 1] .- dV)) < t
            @test maximum(abs.(eval_model_grad(m, T[p], uv)[2][:, 1] .- dV)) < t
        end
    end

    @testset "the count is judged on 64 baselines and is the count all of them give" for T in (Float64, Float32)
        B = T.(collect(range(0.0, 12 / 8, length = 2000)))                # 0 … 12 lobes
        for (prof, pn, pv) in (("\$MU^\$alpha", ["alpha"], [0.15]),
                               ("(\$R >= \$rin) * 1.0", ["rin"], [2.0]),
                               ("max(0, 1 - (2*(\$R - 2.0))^2)", String[], Float64[]),
                               ("1/(1+exp(-(\$R - 1.5)/1e-5))", String[], Float64[]),
                               ("(\$R/1.0)^(-1.5)", String[], Float64[]))
            pf = O.compile_profile(prof, pn); pp = T.(pv)
            lay = O._hankel_probe(pf, pp, 4.0, 0.0, T)
            n0 = O._hankel_n(8 * Float64(maximum(B)))
            nodes(n) = ((s, μ, g) = O._hankel_nodes(T, 0.0, n; layout = lay, n0); (T(4) .* s, μ, g))
            @test length(O._hankel_subset(B, 64)) == 64 && maximum(O._hankel_subset(B, 64)) == maximum(B)
            @test O._hankel_converge(nodes, pf, pp, B, n0; tol = tol(T))[1] ==
                  O._hankel_converge(nodes, pf, pp, B, n0; tol = tol(T), nsub = length(B))[1]
        end
    end

    @testset "an inclined component is checked once, and is right" begin
        md = Dict{String,Any}("s,diam" => 4.0, "s,profile" => "\$MU^0.14", "s,incl" => 40.0, "s,projang" => 30.0)
        ma = Dict{String,Any}("s,ldpow" => 4.0, "s,alpha" => 0.14, "s,incl" => 40.0, "s,projang" => 30.0)
        uv = uvof(Bset(3, 4.0)); uv[2, :] .= 0.3 .* uv[1, :]
        m = dict_to_model(md, ["s,incl"])
        eval_model(m, [40.0], uv); l1 = only(m.components).kcache.layout
        eval_model(m, [45.0], uv); eval_model_jac(m, [50.0], uv)
        @test only(m.components).kcache.layout === l1
        @test maximum(abs.(eval_model(m, [45.0], uv) .- eval_model(dict_to_model(ma, ["s,incl"]), [45.0], uv))) < 1e-8
    end

    @testset "the direct transform's rrule, with the B, r and g tangents" begin
        s, μ, g = O._hankel_nodes(Float64, 0.0, 32); r = 2 .* s; I = μ .^ 0.14; B = Bset(3, 4.0)
        ȳ = sin.(1:length(B))
        H, pb = O.ChainRulesCore.rrule(O._hankel_wtransform, I, r, g, B; n = 1)
        _, Ī, _, _, B̄ = pb(ȳ)
        @test H ≈ O._hankel_wtransform(I, r, g, B; n = 1)
        @test Ī ≈ ForwardDiff.gradient(x -> sum(ȳ .* O._hankel_wtransform(x, r, g, B; n = 1)), I) rtol = 1e-10
        @test B̄ ≈ ForwardDiff.gradient(x -> sum(ȳ .* O._hankel_wtransform(I, r, g, x; n = 1)), B) rtol = 1e-8
        _, _, r̄, ḡ, _ = pb(ȳ)
        @test r̄ ≈ ForwardDiff.gradient(x -> sum(ȳ .* O._hankel_wtransform(I, x, g, B; n = 1)), r) rtol = 1e-8
        @test ḡ ≈ ForwardDiff.gradient(x -> sum(ȳ .* O._hankel_wtransform(I, r, x, B; n = 1)), g) rtol = 1e-8
    end

    @testset "a free size: diameters, the inner edge and \$RMAX, in V and its derivatives" for T in (Float64, Float32)
        t = tol(T) * 10
        B = Bset(3, 8.0); uv = T.(uvof(B))
        ud(D, b) = (x = π * D * b; x == 0 ? 1.0 : 2besselj1(x) / x)
        # a ring of constant brightness: the difference of two uniform discs
        ring(a, b) = [((b / 2)^2 * ud(b, x) - (a / 2)^2 * ud(a, x)) / ((b / 2)^2 - (a / 2)^2) for x in B]
        m = dict_to_model(Dict{String,Any}("s,diamin" => 3.0, "s,udout" => 8.0, "s,profile" => "1.0"),
                          ["s,diamin", "s,udout"]; T)
        h = 1e-6
        Ja = (ring(3 + h, 8) .- ring(3 - h, 8)) ./ 2h; Jb = (ring(3, 8 + h) .- ring(3, 8 - h)) ./ 2h
        @test maximum(abs.(eval_model(m, T[3, 8], uv) .- ring(3, 8))) < t
        for J in (eval_model_jac(m, T[3, 8], uv)[2], eval_model_grad(m, T[3, 8], uv)[2])
            @test maximum(abs.(J[:, 1] .- Ja)) < t && maximum(abs.(J[:, 2] .- Jb)) < t
        end
        # the size changed: followed, and the matrices are those of the new size
        @test maximum(abs.(eval_model(m, T[2.5, 7], uv) .- ring(2.5, 7))) < t
        # a breakpoint fixed in R moves in s with the size; `$RMAX` is the size
        for (prof, pv, free) in (("(\$R >= 2.0) * exp(-\$R/\$a)", [1.5], ["s,udout", "s,a"]),
                                 ("1 - (\$R/\$RMAX)^2 * \$c", [0.5], ["s,udout", "s,c"]),
                                 ("exp(-(\$R - \$RMIN)/\$c)", [1.5], ["s,diamin", "s,udout", "s,c"]))
            md = Dict{String,Any}("s,diamin" => 2.0, "s,udout" => 8.0, "s,profile" => prof,
                                  "s,$(free[end][3:end])" => pv[1])
            m = dict_to_model(md, free; T); x = T[md[k] for k in free]
            m64 = dict_to_model(md, free)
            Jd = reduce(hcat, [(eval_model(m64, x .+ h .* (1:length(x) .== j), uvof(B)) .-
                                eval_model(m64, x .- h .* (1:length(x) .== j), uvof(B))) ./ 2h
                               for j in eachindex(x)])
            @test maximum(abs.(eval_model_jac(m, x, uv)[2] .- Jd)) < 10t
            @test maximum(abs.(eval_model_grad(m, x, uv)[2] .- Jd)) < 10t
        end
        # a fixed size: the parameters' columns are those of the same model with the size free
        md = Dict{String,Any}("s,udout" => 8.0, "s,profile" => "exp(-(\$R/\$w)^2)", "s,w" => 2.0)
        Jf = eval_model_jac(dict_to_model(md, ["s,w"]; T), T[2], uv)[2]
        Js = eval_model_jac(dict_to_model(md, ["s,udout", "s,w"]; T), T[8, 2], uv)[2]
        @test Jf[:, 1] ≈ Js[:, 2] rtol = 10eps(T)
    end

    @testset "α Cen A, the diameter free: a written power law fits as the closed-form one" begin
        d = readoifits(joinpath(@__DIR__, "..", "demos", "data", "AlphaCenA.oifits");
                       T = Float64, warn = false, verbose = false)[1, 1]
        rc = fit_model_lsqfit(dict_to_model(Dict{String,Any}("s,ldpow" => 8.0, "s,alpha" => 0.3),
                                            ["s,ldpow", "s,alpha"]), [8.0, 0.3], d)
        rp = fit_model_lsqfit(dict_to_model(Dict{String,Any}("s,diam" => 8.0, "s,profile" => "\$MU^\$alpha",
                                                             "s,alpha" => 0.3), ["s,diam", "s,alpha"]), [8.0, 0.3], d)
        @test rp.x_opt ≈ rc.x_opt rtol = 1e-6
        @test rp.chi2 ≈ rc.chi2 rtol = 1e-8
    end

    @testset "the GUI preview's quadrature is the model's" begin
        B = Bset(3, 8.0)
        md = Dict{String,Any}("d,profile" => "(\$R >= 2.0) * exp(-\$R/\$a)", "d,a" => 1.5, "d,udout" => 8.0)
        V = vis(md, B)
        Vp, n, lay, ok, _ = O.hankel_profile_vis(O.compile_profile(md["d,profile"], ["a"]), [1.5], 0.0, 4.0, B)
        @test ok && lay.kinds == [:jump]
        @test maximum(abs.(Vp .- V)) < 1e-12
    end
end
