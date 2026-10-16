# reconstruct_sparco on V1295 Aql: the star + background + grey image model.
using Test, OITOOLS

@testset "reconstruct_sparco" begin
    file = joinpath(@__DIR__, "..", "demos", "data", "V1295_Aql.oifits")
    data = readoifits(file; filter_bad_data = true, verbose = false)[1, 1]
    nx, pixsize = 64, 0.125
    ft = setup_nfft(data, nx, pixsize)
    x0 = gaussian2d(nx, nx, nx / 6)
    ndata = data.nv2 + data.nt3phi
    for first in (true, false)
        r = reconstruct_sparco(x0, data, ft; lambda_ref = 1.6e-6, star_flux = 0.48, d_env = -0.1,
                               weights = [1.0, 0.0, 1.0], regularizers = [["tv", 1e2]],
                               maxiter = 100, rounds = 2, fit_model_first = first)
        fstar, fbg = r.param_names["star,f0"], r.param_names["bg,f0"]
        # the image's flux at λ0, 1 − f_star − f_bg, stays non-negative
        @test 0 <= fstar <= 1 && 0 <= fbg && fstar + fbg <= 1
        @test fbg ≈ (1 - fstar) * r.param_names["bg,frac"]
        @test r.free_params == ["star,f0", "bg,frac", "d_env"]
        @test r.chi2 / ndata < 1.5
    end
    # the starting fluxes, converted to bg,frac, are the model's
    r = reconstruct_sparco(x0, data, ft; lambda_ref = 1.6e-6, star_flux = 0.4, bg_flux = 0.3,
                           weights = [1.0, 0.0, 1.0], maxiter = 1, rounds = 1)
    @test r.model.list_free_params == ["star,f0", "bg,frac", "d_env"]
    @test_throws ArgumentError reconstruct_sparco(x0, data, ft; lambda_ref = 1.6e-6,
                                                  star_flux = 0.8, bg_flux = 0.3)
    # a background start at the bound: the image weight at λ0 is zero, not negative
    m = r.model
    vals = m.resolver.fn([0.6, 1.0, 0.0], [1.6, 2.2], NaN)
    @test only(unique(vals[findfirst(==("bg,f0"), m.all_names)])) ≈ 0.4
    @test all(abs.(vals[OITOOLS._resolve_w_idx(m, "W")]) .< 1e-12)
end
