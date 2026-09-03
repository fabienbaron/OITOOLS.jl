# Absolute-scale check on the analytic components.
#
# Why this exists: nothing else pinned it. The PMOIRED comparison in `demos/` asserts only
# V(0) = 1, |V| <= 1 and monotonicity, all of which a Gaussian of the wrong width satisfies,
# and `test_model_gradients.jl` validates derivatives against the same formula it is testing —
# so a component whose size parameter means the wrong thing passes both.
#
# It did. The Gaussian family fed the `exp(-t^2/ln2)` form a sigma where it takes a HALF-width,
# making every Gaussian narrower than its stated FWHM by sqrt(2 ln 2) = 1.177, and a fitted FWHM
# about 18% high.
#
# The reference here is `hankel_vis` of the brightness distribution the parameter NAMES —
# numerical, independent of the analytic formula, and itself checked below against the uniform
# disc, whose transform is known in closed form.

@testset "component widths mean what they say" begin
    B = collect(range(0.0, 0.6; length = 120))          # cycles/mas
    ρ = B .* 2.0626480624709636e8                        # cycles/rad, as eval_model wants
    uv = permutedims(hcat(ρ, zeros(length(ρ))))

    V(md) = real.(OITOOLS.eval_model(dict_to_model(md, String[]), Float64[], uv))
    "numerical transform of I(r) on a grid fine enough that quadratic convergence has landed"
    function numerical(I_of_r, rmax; n = 40_000)
        r = collect(range(0.0, rmax; length = n))
        OITOOLS.hankel_vis(I_of_r.(r), r, B)
    end

    @testset "the reference itself: uniform disc" begin
        D = 4.0
        num = numerical(r -> r <= D/2 ? 1.0 : 0.0, D/2)
        ana = OITOOLS.analytic_ud_vis(D, B)
        @test maximum(abs.(num .- ana)) < 1e-3          # trapezoid on a hard edge
    end

    @testset "gaussian: fwhm is the full width at half maximum" begin
        F = 3.0
        # I(r) = exp(-4 ln2 r^2 / F^2) is, by definition, a Gaussian of FWHM F.
        num = numerical(r -> exp(-4*log(2) * r^2 / F^2), 8F)
        @test maximum(abs.(V(Dict{String,Any}("g,fwhm" => F, "g,f" => 1.0)) .- num)) < 1e-4
        # and against the closed form mfit uses, exp(-(pi a rho)^2 / (4 ln2)) with a = FWHM
        @test maximum(abs.(V(Dict{String,Any}("g,fwhm" => F, "g,f" => 1.0)) .-
                           (@. exp(-(π*F*B)^2 / (4*log(2)))))) < 1e-9
    end

    @testset "spatial_kernel smooths by its stated FWHM" begin
        F = 2.5
        # A point source smoothed by a kernel of FWHM F IS a Gaussian of FWHM F.
        pt = V(Dict{String,Any}("p,f" => 1.0, "spatial_kernel" => F))
        gs = V(Dict{String,Any}("g,fwhm" => F, "g,f" => 1.0))
        @test maximum(abs.(pt .- gs)) < 1e-12
    end

    @testset "ud: the parameter is the diameter" begin
        D = 5.0
        @test maximum(abs.(V(Dict{String,Any}("s,ud" => D, "s,f" => 1.0)) .-
                           OITOOLS.analytic_ud_vis(D, B))) < 1e-9
    end

    @testset "limb-darkening laws agree where they coincide" begin
        # Exact, and the sharpest check there is on the normalisations: each law reduces to a
        # simpler one when its extra coefficients vanish, and any mistake in the per-term
        # weights of N shows up here as a non-zero difference rather than as a plausible curve.
        θ, u = 7.0, 0.3
        lin = V(Dict{String,Any}("s,ldlin" => θ, "s,u" => u, "s,f" => 1.0))
        ud  = V(Dict{String,Any}("s,ud" => θ, "s,f" => 1.0))

        @test V(Dict{String,Any}("s,ldsqrt" => θ, "s,u" => u, "s,w" => 0.0,
                                 "s,f" => 1.0)) ≈ lin
        @test V(Dict{String,Any}("s,ldsqrt" => θ, "s,u" => 0.0, "s,w" => 0.0,
                                 "s,f" => 1.0)) ≈ ud
        @test V(Dict{String,Any}("s,ldclaret4" => θ, "s,c1" => 0.0, "s,c2" => u,
                                 "s,c3" => 0.0, "s,c4" => 0.0, "s,f" => 1.0)) ≈ lin
        @test V(Dict{String,Any}("s,ldclaret4" => θ, "s,c1" => 0.0, "s,c2" => 0.0,
                                 "s,c3" => 0.0, "s,c4" => 0.0, "s,f" => 1.0)) ≈ ud
    end

    @testset "limb-darkening laws match their own I(mu)" begin
        # Against a direct integration of the profile each law NAMES, so the series expansion
        # and the brightness distribution cannot drift apart.
        function from_profile(Imu, θ; n = 60_001)
            ζ = @. π * θ * ρ / 2.0626480624709636e8
            μ = range(0, 1; length = n); h = step(μ)
            trap(v) = h * (sum(v) - (v[1] + v[end]) / 2)
            den = trap([Imu(m) * m for m in μ])
            [trap([Imu(m) * OITOOLS.besselj0(z * sqrt(max(0, 1 - m^2))) * m for m in μ]) / den
             for z in ζ]
        end
        θ = 6.0
        c, d = 0.35, 0.20
        @test maximum(abs.(V(Dict{String,Any}("s,ldsqrt" => θ, "s,u" => c, "s,w" => d,
                                              "s,f" => 1.0)) .-
                           from_profile(μ -> 1 - c*(1-μ) - d*(1-sqrt(μ)), θ))) < 1e-6
        cc = [0.40, -0.20, 0.50, -0.15]
        @test maximum(abs.(V(Dict{String,Any}("s,ldclaret4" => θ, "s,c1" => cc[1],
                                              "s,c2" => cc[2], "s,c3" => cc[3],
                                              "s,c4" => cc[4], "s,f" => 1.0)) .-
                           from_profile(μ -> 1 - sum(cc[k]*(1-μ^(k/2)) for k in 1:4), θ))) < 1e-6
    end

    @testset "a flat profile is a uniform disc of diameter udout" begin
        D = 6.0
        prof = V(Dict{String,Any}("s,profile" => "1.0", "s,udout" => D,
                                  "s,nr" => 2000, "s,f" => 1.0))
        @test maximum(abs.(prof .- OITOOLS.analytic_ud_vis(D, B))) < 1e-4
    end

    # A visibility written out must agree with the analytic kind it reproduces, in the units
    # the expression is written in. `$B` is Mλ, so the conversion to the cycles/rad the model
    # works in is part of what this pins: an expression right in one unit and evaluated in the
    # other still looks like a plausible curve.
    @testset "a custom visibility function is the model it spells out" begin
        D = 5.0
        vf = Dict{String,Any}(
            "c,visfunc" => "2*besselj1(pi*\$dd*\$B/206.264806) / (pi*\$dd*\$B/206.264806)",
            "c,dd" => D, "c,f" => 1.0)
        @test OITOOLS._identify_kind("c", Dict{String,Any}(vf)) === :visfunc

        # away from B = 0, where 2J₁(x)/x is 0/0; the guarded form is checked below
        Bv = collect(range(0.02, 0.6; length = 80))          # cycles/mas
        uv = vcat((Bv .* OITOOLS._MAS2RAD_PM)', zeros(1, length(Bv)))
        Vv = eval_model(dict_to_model(vf, String[]), Float64[], uv)
        @test maximum(abs.(real.(Vv) .- OITOOLS.analytic_ud_vis(D, Bv))) < 1e-8

        # V IS evaluated at B = 0 -- that point is the total flux, and `chi2_flat`'s flux term
        # and `model_to_sed` both ask for it. An unguarded 0/0 there returns NaN and poisons
        # the whole chi2, so the guard is part of writing one of these.
        z = vcat([0.0], [0.0])
        @test isnan(real(eval_model(dict_to_model(vf, String[]), Float64[],
                                    reshape(z, 2, 1))[1]))
        vg = merge(vf, Dict{String,Any}("c,visfunc" =>
                 "ifelse(\$B < 1e-8, 1.0, 2*besselj1(pi*\$dd*\$B/206.264806) / " *
                 "(pi*\$dd*\$B/206.264806))"))
        @test real(eval_model(dict_to_model(vg, String[]), Float64[],
                              reshape(z, 2, 1))[1]) ≈ 1.0

        # the expression's own parameter is discovered and fittable, and `$B` is not one
        fm = dict_to_model(vf, ["c,dd"])
        @test fm.list_free_params == ["c,dd"]
        @test "c,B" ∉ fm.all_names && "B" ∉ fm.all_names

        # and it is AD-transparent, which is what makes it usable by a gradient fitter
        f(x) = sum(abs2, eval_model(fm, x, uv))   # uv is the B > 0 grid built above
        # `OITOOLS.ForwardDiff`, not a bare `using`: this file is included into Main by
        # runtests.jl, where a `using` at file scope does not reach. OITOOLS depends on
        # ForwardDiff already, so nothing new is needed in the test environment.
        g  = OITOOLS.ForwardDiff.gradient(f, [D])[1]
        h  = 1e-6
        @test g ≈ (f([D+h]) - f([D-h]))/(2h) rtol=1e-6
    end

    # Azimuthal modes carry a sign convention, and a wrong one is invisible in |V|: it rotates
    # the asymmetry rather than changing its magnitude, so only a COMPLEX comparison catches it.
    #
    # The ground truth is an image built with PMOIRED's own image-domain formula --
    # `I(r)·(1 + Σ ampₖ cos(nₖ(PA + φₖ − π/2)))` with `PA = atan2(y, x)`, from `_Vazvar` in its
    # oimodels.py -- numerically Fourier transformed. That pins BOTH that the analytic Hankel
    # sum is right and that `projang` means in OITOOLS what it means in PMOIRED.
    @testset "azimuthal modes follow PMOIRED's projang convention" begin
        rin, rout = 3.0, 5.0
        amp, phi, order = 0.4, 37.0, 1
        mas = π/180/3600/1000

        # numerical FT of the image, on a grid fine enough to resolve the ring edges
        nx, fov = 768, 24.0
        pix = fov/nx
        cen = ((1:nx) .- (nx+1)/2) .* pix
        Bm  = collect(4.0:8.0:52.0)                      # Mλ
        th  = deg2rad(23.0)
        u, v = (Bm .* 1e6) .* cos(th), (Bm .* 1e6) .* sin(th)

        Vimg = zeros(ComplexF64, length(u)); tot = 0.0
        for j in 1:nx, i in 1:nx
            x, y = cen[i], cen[j]
            r = hypot(x, y)
            (rin <= r <= rout) || continue
            I = 1 + amp*cos(order*(atan(y, x) + deg2rad(phi) - π/2))
            tot += I
            @. Vimg += I * cis(-2π*(u*x + v*y)*mas)
        end
        Vimg ./= tot

        md = Dict{String,Any}("d,profile" => "(\$R >= $rin) * (\$R <= $rout) * 1.0",
                              "d,udout" => 2rout, "d,nr" => 400.0, "d,f" => 1.0,
                              "d,az amp1" => amp, "d,az projang1" => phi)
        Vmod = eval_model(dict_to_model(md, String[]), Float64[], vcat(u', v'))

        # 1e-3 is the pixel grid's own accuracy, not the model's; the wrong sign would put an
        # m=1 term of amplitude 0.4 in the wrong place and miss by two orders of magnitude more.
        @test maximum(abs.(Vmod .- Vimg)) < 1e-3
        @test maximum(abs.(angle.(Vmod ./ Vimg))) < deg2rad(0.5)

        # and the discriminating half: the opposite sign does NOT reproduce the same image
        md_flip = merge(md, Dict{String,Any}("d,az projang1" => phi + 180.0))
        Vflip = eval_model(dict_to_model(md_flip, String[]), Float64[], vcat(u', v'))
        @test maximum(abs.(Vflip .- Vimg)) > 0.05
    end
end
