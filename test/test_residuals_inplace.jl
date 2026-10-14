using OITOOLS, Test

# `residuals_flat!` / `residuals_flat_jac!` against `residuals_flat` / `residuals_flat_jac`:
# the same rows, bit for bit, for every observable, both precisions, the per-point paths
# (`$WL`, an orbit on `$MJD`), smearing (the allocating path, copied in) and the flux term.
@testset "residuals written in place" begin
    f = joinpath(@__DIR__, "..", "demos", "data", "V1295_Aql.oifits")
    models = [
        ("binary", Dict{String,Any}("1,ud" => 0.6, "2,f" => 0.3, "2,ud" => 0.2, "2,x" => 2.0, "2,y" => -1.5),
         ["1,ud", "2,f", "2,x", "2,y"]),
        ("chromatic", Dict{String,Any}("1,ud" => "\$d0 * (1 + 0.1 * (\$WL - 1.6))", "d0" => 0.6,
                                      "2,f" => "\$f0 * (\$WL / 1.6)^-2", "f0" => 0.3,
                                      "2,x" => 2.0, "2,y" => -1.5),
         ["d0", "f0", "2,x"]),
        ("orbit", Dict{String,Any}("1,ud" => 0.6, "2,f" => 0.3, "2,x" => "orbit", "2,y" => "orbit",
                                  "2,orb P" => 30.0, "2,orb MJD0" => 58640.0, "2,orb e" => 0.4,
                                  "2,orb omega" => 60.0, "2,orb OMEGA" => 120.0, "2,orb incl" => 50.0,
                                  "2,orb a" => 4.0),
         ["1,ud", "2,f", "2,orb e", "2,orb a", "2,orb MJD0"]),
        ("smeared", Dict{String,Any}("1,ud" => 0.6, "2,f" => 0.3, "2,x" => 8.0, "2,y" => -5.0, "smear" => 3),
         ["1,ud", "2,f", "2,x"]),
    ]
    weightsets = ([1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0],
                  [1.0, 0.5, 2.0, 1.0, 1.0, 0.0, 0.0])
    for T in (Float64, Float32)
        d = readoifits(f; warn = false, verbose = false, T)[1, 1]
        for (name, m, free) in models, w in weightsets
            fm = dict_to_model(m, free)
            x = T[m[k] for k in free]
            r0 = OITOOLS.residuals_flat(fm, x, d; weights = w)
            r = fill(NaN, length(r0))
            @test OITOOLS.residuals_flat!(r, fm, x, d; weights = w) === r
            @test r == r0
            rj0, J0 = OITOOLS.residuals_flat_jac(fm, x, d; weights = w)
            rj, J = fill(NaN, length(rj0)), fill(NaN, size(J0))
            OITOOLS.residuals_flat_jac!(rj, J, fm, x, d; weights = w, jac_cache = OITOOLS.ModelJacCache())
            @test rj == rj0
            @test J == J0
        end
    end

    # the flux term, calibrated (C = 1) and not (C fitted), on a dataset given three fluxes
    d = readoifits(f; warn = false, verbose = false, T = Float64)[1, 1]
    d.flux = [1.0, 0.9, 1.1]; d.flux_err = [0.05, 0.04, 0.06]; d.nflux = 3
    m, free = models[2][2], models[2][3]
    fm = dict_to_model(m, free); x = Float64[m[k] for k in free]
    wf = [1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0]
    for cal in (true, false)
        d.flux_calibrated = cal
        r0 = OITOOLS.residuals_flat(fm, x, d; weights = wf)
        @test OITOOLS.residuals_flat!(similar(r0), fm, x, d; weights = wf) == r0
        rj0, J0 = OITOOLS.residuals_flat_jac(fm, x, d; weights = wf)
        rj, J = OITOOLS.residuals_flat_jac!(similar(rj0), similar(J0), fm, x, d; weights = wf)
        @test rj == rj0
        @test J == J0
    end

    # unfiltered data: an index holds more entries than the valid-point count, and the rows of
    # the invalid ones are NaN in both (`isequal`: NaN equals NaN)
    du = readoifits(f; warn = false, verbose = false, filter_bad_data = false, T = Float64)[1, 1]
    @test length(du.indx_vis) > du.nvisamp
    fm = dict_to_model(models[1][2], models[1][3]); x = Float64[models[1][2][k] for k in models[1][3]]
    w3 = weightsets[3]
    r0 = OITOOLS.residuals_flat(fm, x, du; weights = w3)
    @test isequal(OITOOLS.residuals_flat!(similar(r0), fm, x, du; weights = w3), r0)
    rj0, J0 = OITOOLS.residuals_flat_jac(fm, x, du; weights = w3)
    rj, J = OITOOLS.residuals_flat_jac!(similar(rj0), similar(J0), fm, x, du; weights = w3)
    @test isequal(rj, rj0)
    @test isequal(J, J0)

    # a buffer of the wrong length is an error, not a silent truncation
    fm = dict_to_model(models[1][2], models[1][3]); x = Float64[models[1][2][k] for k in models[1][3]]
    n = length(OITOOLS.residuals_flat(fm, x, d))
    @test_throws DimensionMismatch OITOOLS.residuals_flat!(zeros(n + 1), fm, x, d)
    @test_throws DimensionMismatch OITOOLS.residuals_flat_jac!(zeros(n - 1), zeros(n - 1, 4), fm, x, d)
end
