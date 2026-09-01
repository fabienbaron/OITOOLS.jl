# test_model_gradients.jl — finite-difference validation of the analytic gradients
# used by every gradient-based parametric fit.
#
# Why this exists: `fit_model` (NLopt LD_*) and `fit_model_lsqfit` (Levenberg-Marquardt)
# are driven entirely by `chi2_flat_fg` / `residuals_flat_jac`, whose derivatives are
# hand-written ChainRules `rrule`s in model_chainrules.jl plus a hand-written ForwardDiff
# specialisation for `_vis_ldpow` in parse_model.jl. A wrong derivative there does not
# throw and does not produce a NaN — the optimiser simply stops somewhere that is not the
# minimum and reports success, and `stderror` (built from the same Jacobian) comes out
# confidently wrong too.
#
# That is not hypothetical. `∂V/∂α` of the power-law disc carried a `log(ζ)` where the
# derivation needs `log(ζ/2)`, making the derivative too small by exactly `log(2)·V/2`.
# On the Kervella+2017 α Cen A PIONIER data that stopped LD_LBFGS and LsqFit at
# α = 0.121 instead of 0.136 (Δχ² = 34) while Nelder-Mead, COBYLA and UltraNest — none of
# which use the gradient — all found the true minimum. Nothing in the suite noticed:
# test_grad.jl and test_gradient_sign.jl cover the *image* (DFT/NFFT) path only, and
# neither is included in runtests.jl.
#
# Differences come from FiniteDifferences.jl rather than a hand-rolled stencil: it picks
# the step adaptively per point, so these assertions do not silently become a test of a
# badly chosen `h` on some future model whose χ² has a different scale.

using OITOOLS, Test, Printf, FiniteDifferences

# 5th-order central stencil: accurate enough that a tight rtol below is testing the
# analytic derivative rather than the truncation error of the comparison.
const FDM = central_fdm(5, 1)

@testset "parametric model gradients" begin

    data = readoifits(joinpath(@__DIR__, "..", "demos", "data", "AlphaCenA.oifits");
                      warn = false, verbose = false)[1,1]
    W = [1.0, 1.0, 1.0]      # V2 + T3amp + T3phi, so every observable's chain rule is exercised

    # (name, model dict, free parameters, evaluation point)
    #
    # The evaluation point is deliberately NOT the best fit: at a minimum the gradient is
    # zero and any error would be masked. Diameters are ~8 mas against baselines that
    # resolve α Cen A into the third lobe, so ζ = πθρ/C reaches ~11 here — well past the
    # first Bessel zero, which is where a derivative error starts to bias a diameter.
    cases = [
        ("ud",            Dict{String,Any}("s,ud"=>8.0, "s,f"=>1.0),
                          ["s,ud"], [8.2]),
        ("gaussian",      Dict{String,Any}("s,fwhm"=>6.0, "s,f"=>1.0),
                          ["s,fwhm"], [6.3]),
        ("ldlin",         Dict{String,Any}("s,ldlin"=>8.0, "s,u"=>0.2, "s,f"=>1.0),
                          ["s,ldlin","s,u"], [8.2, 0.25]),
        ("ldquad",        Dict{String,Any}("s,ldquad"=>8.0, "s,u"=>0.2, "s,w"=>0.05, "s,f"=>1.0),
                          ["s,ldquad","s,u","s,w"], [8.2, 0.25, -0.04]),
        ("ldpow",         Dict{String,Any}("s,ldpow"=>8.0, "s,alpha"=>0.15, "s,f"=>1.0),
                          ["s,ldpow","s,alpha"], [8.2, 0.13]),
        # Square root and Claret four-parameter: both are sums of half-integer powers of mu,
        # so their rrules are hand-written like ldquad's. Coefficients away from zero and of
        # BOTH signs, because each enters the normalisation N as well as the numerator and a
        # sign error there cancels at c = 0.
        ("ldsqrt",        Dict{String,Any}("s,ldsqrt"=>8.0, "s,u"=>0.30, "s,w"=>0.20,
                                           "s,f"=>1.0),
                          ["s,ldsqrt","s,u","s,w"], [8.2, 0.35, -0.18]),
        ("ldclaret4",     Dict{String,Any}("s,ldclaret4"=>8.0, "s,c1"=>0.40, "s,c2"=>-0.20,
                                           "s,c3"=>0.50, "s,c4"=>-0.15, "s,f"=>1.0),
                          ["s,ldclaret4","s,c1","s,c2","s,c3","s,c4"],
                          [8.2, 0.45, -0.25, 0.44, -0.11]),
        ("ring",          Dict{String,Any}("s,diamin"=>4.0, "s,diamout"=>8.0, "s,f"=>1.0),
                          ["s,diamin","s,diamout"], [4.3, 8.4]),
        ("gaussian_ring", Dict{String,Any}("s,fwhmin"=>3.0, "s,fwhmout"=>7.0, "s,f"=>1.0),
                          ["s,fwhmin","s,fwhmout"], [3.2, 7.3]),
        # two components: exercises flux weighting and the x/y phase factors
        ("ud + point",    Dict{String,Any}("s,ud"=>8.0, "s,f"=>0.9,
                                           "p,x"=>2.0, "p,y"=>-1.5, "p,f"=>0.1),
                          ["s,ud","s,f","p,x","p,y"], [8.2, 0.85, 2.3, -1.7]),
        # radial profile: the Hankel branch rather than the analytic one
        ("hankel profile", Dict{String,Any}("r,profile"=>"exp(-(\$R/\$r,s)^2)",
                                            "r,udout"=>9.0, "r,s"=>2.0, "r,f"=>1.0),
                          ["r,s","r,f"], [2.2, 0.9]),
    ]

    for (name, md, free, x0) in cases
        @testset "$name" begin
            fm = OITOOLS.parse_model(md, free; nB_workspace = size(data.uv, 2))
            f(p) = OITOOLS.chi2_flat(fm, p, data; weights = W)

            _, g   = OITOOLS.chi2_flat_fg(fm, collect(x0), data; weights = W)
            g_fd,  = FiniteDifferences.grad(FDM, f, collect(x0))

            @test all(isfinite, g)
            @test g ≈ g_fd rtol = 1e-5
        end
    end

    # ── the ForwardDiff specialisation of _vis_ldpow ──────────────────────────
    # parse_model.jl carries a hand-written Dual-number method for the power law
    # (besselj(ν, x) returns NaN when ν is a Dual). It duplicates the rrule's algebra,
    # so it can — and did — carry the same mistake independently. Both are compared to
    # finite differences of the primal, never to each other, so a shared error in the two
    # hand derivations cannot cancel.
    @testset "_vis_ldpow ForwardDiff specialisation and rrule" begin
        ρ = [5e6, 1.5e7, 3e7, 5e7, 8e7]
        p0 = [8.45, 0.136]

        J_fd, = FiniteDifferences.jacobian(FDM, p -> OITOOLS.vis_ldpow(p[1], p[2], ρ), p0)
        J_ad  = OITOOLS.ForwardDiff.jacobian(p -> OITOOLS._vis_ldpow(p, ρ), p0)
        @test J_ad ≈ J_fd rtol = 1e-6

        _, pb = OITOOLS.ChainRulesCore.rrule(OITOOLS.vis_ldpow, p0[1], p0[2], ρ)
        J_rr = reduce(vcat, begin
                   seed = zeros(length(ρ)); seed[j] = 1.0
                   _, gθ, gα = pb(seed)
                   [gθ gα]
               end for j in eachindex(ρ))
        @test J_rr ≈ J_fd rtol = 1e-6
    end

    # ── ∂J_ν/∂ν across the full range a fit can reach ────────────────────────
    # The ascending series `_dbesselj_dnu` uses is cancellation-limited and returns
    # garbage above x ≈ 20 (1.6e5 at x = 30, 1e24 at x = 60) — see the note on
    # _DBESSEL_SERIES_XMAX. α Cen A only reaches ζ ≈ 11, but 8.5 mas on a 330 m CHARA
    # baseline reaches ζ ≈ 26 and a 43 mas supergiant on the VLTI reaches ζ ≈ 39, so the
    # range below deliberately runs well past where the series alone stops working. It
    # also straddles the crossover so a botched threshold shows up as a discontinuity.
    @testset "_dbesselj_dnu accuracy" begin
        worst = 0.0
        for ζ in range(0.05, 300.0, length = 400), ν in (1.0, 1.07, 1.15, 1.5, 2.0)
            fd = FDM(n -> OITOOLS.SpecialFunctions.besselj(n, ζ), ν)
            worst = max(worst, abs(OITOOLS._dbesselj_dnu(ν, ζ) - fd))
        end
        @test worst < 1e-8

        # continuity across the series → stencil crossover
        x0 = OITOOLS._DBESSEL_SERIES_XMAX
        for ν in (1.0, 1.07, 1.5)
            @test OITOOLS._dbesselj_dnu(ν, x0 - 1e-9) ≈
                  OITOOLS._dbesselj_dnu(ν, x0 + 1e-9) atol = 1e-9
        end
    end

    # ── end to end: gradient-based and derivative-free fitters must agree ─────
    # The real symptom of a bad derivative is not a NaN, it is two optimisers quietly
    # disagreeing. This asserts they do not.
    @testset "optimiser agreement (power law, α Cen A)" begin
        md   = Dict{String,Any}("s,ldpow"=>8.0, "s,alpha"=>0.15, "s,f"=>1.0)
        free = ["s,ldpow", "s,alpha"]
        lb   = Dict("s,ldpow"=>4.0,  "s,alpha"=>0.0)
        ub   = Dict("s,ldpow"=>12.0, "s,alpha"=>0.6)
        w    = [1.0, 0, 0]

        r_grad = fit_model(md, free, data; lb, ub, weights = w)                       # LD_LBFGS
        r_free = fit_model(md, free, data; lb, ub, weights = w,
                           method = :LN_NELDERMEAD, maxeval = 20000)
        r_lm   = fit_model_lsqfit(md, free, data; lb, ub, weights = w)

        @test r_grad.chi2  ≈ r_free.chi2  rtol = 1e-4
        @test r_lm.chi2    ≈ r_free.chi2  rtol = 1e-4
        @test r_grad.x_opt ≈ r_free.x_opt rtol = 1e-3
        @test r_lm.x_opt   ≈ r_free.x_opt rtol = 1e-3

        # Kervella et al. 2017 (A&A 597, A137) Table 3, power law, α Cen A:
        # θ_LD = 8.502 ± 0.006 mas, α = 0.1404 ± 0.0050, after multiplying the PIONIER
        # wavelength scale by γ = 1.00481 (their Eq. 1). The OIFITS here is uncorrected,
        # so the published diameter corresponds to 8.502/γ = 8.4613 mas.
        @test r_free.x_opt[1] ≈ 8.4613 atol = 0.02
        @test r_free.x_opt[2] ≈ 0.1404 atol = 0.02
    end
end
