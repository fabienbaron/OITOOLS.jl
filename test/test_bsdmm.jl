# test_bsdmm.jl — the ADMM reconstructor, monochromatic and as a wavelength cube.
#
# BSDMM had no tests at all. These pin what the engine promises rather than what it prints: the
# shapes, positivity, that it reduces chi2 from its start, that the history carries the
# residuals, and that the cross-channel group block actually changes a cube. Iteration counts
# are deliberately small — these are checks of the machinery, not of convergence. The cube
# criterion BSDMM minimises is pinned here too, on a file that carries OI_FLUX.

using OITOOLS, Test

@testset "BSDMM" begin
    datadir = joinpath(@__DIR__, "..", "demos", "data")
    chi2(img, ft, data) = Float64(OITOOLS.image_to_chi2(OITOOLS.to_ft_precision(Float64.(img), ft),
                                                        ft, data; weights = [1.0, 1.0, 1.0],
                                                        verb = false))
    gauss(nx, s) = (c = (nx + 1) / 2;
                    g = [exp(-((i - c)^2 + (j - c)^2) / (2s^2)) for i in 1:nx, j in 1:nx];
                    g ./ sum(g))

    @testset "one image" begin
        data = readoifits(joinpath(datadir, "BC2004", "2004-data1.oifits");
                          warn = false, verbose = false)
        nx = 24
        ft = setup_ft(data, nx, 0.4)
        x0 = gauss(nx, 3.0)
        img, hist = reconstruct_bsdmm(x0, data, ft; mu_reg = 1e-3, mu_cen = 1e-3, maxit = 15,
                                      x_maxiter = 5, verb = false, history = true)
        @test size(img) == (nx, nx)
        @test all(isfinite, img) && all(img .>= 0)
        @test chi2(img, ft, data) < chi2(x0, ft, data)
        # The residuals are what say whether ADMM is converging; a history without them would
        # leave the panel nothing to judge a run by.
        @test haskey(hist, :chi2) && haskey(hist, :r_norm) && haskey(hist, :s_norm)
        @test !isempty(hist.chi2)
    end

    @testset "a wavelength cube" begin
        data = readoifits(joinpath(datadir, "2019_v1295Aql.WL_SMOOTH.A.oifits");
                          warn = false, verbose = false, polychromatic = true,
                          merge_oi_wavelength = true, use_vis = false)
        nwav = size(data, 1)
        @test nwav > 1
        nx = 24
        ft = setup_ft(data, nx, 0.4)
        x0 = repeat(gauss(nx, 3.0), 1, 1, nwav, 1)

        plain = reconstruct_bsdmm(x0, data, ft; mu_tv = 1e-3, mu_cen = 1e-3, maxit = 15,
                                  x_maxiter = 5, verb = false)
        @test size(plain) == (nx, nx, nwav, 1)
        @test all(isfinite, plain) && all(plain .>= 0)
        @test chi2(plain, ft, data) < chi2(x0, ft, data)
        # A cube is not one image repeated: the channels see different uv points.
        @test !all(plain[:, :, 1, 1] ./ sum(plain[:, :, 1, 1]) ≈
                   plain[:, :, w, 1] ./ sum(plain[:, :, w, 1]) for w in 2:nwav)

        # The group block is the engine's cross-channel coupling. A weight that changed nothing
        # would be a control on the panel with no effect on the result.
        grouped = reconstruct_bsdmm(x0, data, ft; mu_tv = 1e-3, mu_cen = 1e-3, mu_group = 1e-2,
                                    group_type = :sparsity, maxit = 15, x_maxiter = 5,
                                    verb = false)
        @test size(grouped) == size(plain)
        @test !(grouped ≈ plain)
    end

    @testset "the cube criterion scores OI_FLUX once, under weights[6]" begin
        # Each channel's kernel fits OI_FLUX when weights[6] asks for it. Nothing else may: a
        # second, ungated copy puts chi2 into the cube's total that the ndof does not count and
        # no weight can switch off. On OBJECT2_K such a copy is 1e13 against a V2+T3phi chi2 of
        # 1e6, and the flux normalisation projects its gradient out, so value and gradient
        # disagree and BSDMM's image does not leave its start.
        file = joinpath(datadir, "BC2026", "OBJECT1_N.oifits")
        grey = readoifits(file; warn = false, verbose = false, use_vis = false)[1, 1]
        lo, hi = extrema(Float64.(grey.uv_lam))
        edges = range(lo - 1e-12, hi + 1e-12; length = 4)
        data = readoifits(file; spectralbin = [[edges[i], edges[i+1]] for i in 1:3],
                          use_vis = false, warn = false, verbose = false)
        @test all(d -> d.nflux > 0, data[:, 1])
        nx = 32
        ft = setup_ft(data, nx, 1.0)
        plane = gauss(nx, 5.0)
        cube = repeat(plane, 1, 1, 3, 1)
        g = similar(cube)
        for w in ([1.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0, 0.0, 0.0, 1.0])
            channels = sum(OITOOLS.image_to_chi2(plane, ft[k, 1], data[k, 1]; weights = w)
                           for k in 1:3)
            @test OITOOLS.image_to_chi2(cube, ft, data; weights = w) ≈ channels rtol = 1e-5
            @test OITOOLS.image_to_chi2_fg(cube, g, ft, data; weights = w) ≈ channels rtol = 1e-5
        end
        @test OITOOLS.image_to_chi2(cube, ft, data; weights = [1.0, 1.0, 1.0, 0.0, 0.0, 1.0]) >
              OITOOLS.image_to_chi2(cube, ft, data; weights = [1.0, 1.0, 1.0])
    end
end
