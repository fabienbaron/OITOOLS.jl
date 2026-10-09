# write_oifits: what is written reads back to the same points -- values, errors, wavelengths and
# uv -- on a v1 file (2004), a file with VIS tables and several instruments (V1295 Aql), one
# with OI_FLUX and differential phases (BC2026 OBJECT1_LM), a spectrally binned read and a
# filtered dataset.
@testset "write_oifits round trip" begin
    dir = mktempdir()
    data = joinpath(@__DIR__, "..", "demos", "data")
    srt(x) = sort(Float64.(x))
    same(a, b, f) = srt(getfield(a, f)) == srt(getfield(b, f))
    rd(f; kw...) = readoifits(f; T = Float64, warn = false, verbose = false, kw...)

    for f in (joinpath(data, "Beauty_Contest_2004", "2004-data1.oifits"),
              joinpath(data, "V1295_Aql.oifits"),
              joinpath(data, "Beauty_Contest_2026", "OBJECT1_LM.oifits"))
        d = rd(f)[1, 1]
        out = write_oifits(joinpath(dir, "out.oifits"), d; overwrite = true)
        e = rd(out)[1, 1]
        for fld in (:v2, :v2_err, :v2_lam, :t3amp, :t3phi, :t3phi_err, :t3_lam,
                    :visamp, :visamp_err, :visphi, :visphi_err, :flux, :flux_err, :flux_lam)
            @test same(d, e, fld)
        end
        @test srt(d.uv[1, d.indx_v2]) ≈ srt(e.uv[1, e.indx_v2]) rtol = 1e-12
        @test srt(d.uv[2, d.indx_t3_1]) ≈ srt(e.uv[2, e.indx_t3_1]) rtol = 1e-12
        @test d.amptyp == e.amptyp && d.phityp == e.phityp
        @test d.flux_calibrated == e.flux_calibrated
        @test sort(d.sta_name) == sort(e.sta_name)
    end

    # Bins: every point of every bin, the channels one OI_WAVELENGTH table.
    b = rd(joinpath(data, "V1295_Aql.oifits"); polychromatic = true)
    c = rd(write_oifits(joinpath(dir, "poly.oifits"), b; overwrite = true); polychromatic = true)
    count(bins) = Dict(l => sum(x.nv2 for x in vec(bins) if x.nv2 > 0 && x.v2_lam[1] == l)
                       for l in unique(x.v2_lam[1] for x in vec(bins) if x.nv2 > 0))
    @test count(b) == count(c)
    @test sum(x.nt3phi for x in b) == sum(x.nt3phi for x in c)

    # A filtered dataset writes what was kept.
    g = rd(joinpath(data, "Beauty_Contest_2004", "2004-data1.oifits"))[1, 1]
    h = filter_data(g, set_data_filter(g; baseline_range = [0.0, 60e6]))
    @test 0 < h.nv2 < g.nv2
    k = rd(write_oifits(joinpath(dir, "filt.oifits"), h; overwrite = true))[1, 1]
    @test same(h, k, :v2) && same(h, k, :t3phi)

    # An existing file is not replaced unless asked.
    @test_throws ArgumentError write_oifits(joinpath(dir, "filt.oifits"), h)
end
