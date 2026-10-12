# write_oifits: what is written reads back to the same points -- values, errors, wavelengths and
# uv -- on a v1 file (2004), a file with VIS tables and several instruments (V1295 Aql), one
# with OI_FLUX and differential phases (BC2026 OBJECT1_LM), a spectrally binned read and a
# filtered dataset. Each point keeps the table, so the instrument, it was read from.
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

    # MIRC-X and PIONIER, merged in one file: written back under their own names.
    d = rd(joinpath(data, "V1295_Aql.oifits"))[1, 1]
    w = read(OITOOLS.OIFITS.OIDataSet, write_oifits(joinpath(dir, "v1295.oifits"), d; overwrite = true))
    names = sort([String(i.insname) for i in w.instr])
    @test "MIRCX" in names && any(startswith("PIONIER"), names)
    @test names == sort(unique([d.tables[t].insname for t in vcat(d.v2_tab, d.t3_tab, d.vis_tab)]))
    for db in vcat(w.vis2, w.t3)        # each table's channels are its instrument's
        @test size(db isa OITOOLS.OIFITS.OI_VIS2 ? db.vis2data : db.t3phi, 1) ==
              length(w.instr[findfirst(i -> i.insname == db.insname, w.instr)].eff_wave)
    end
    e = rd(joinpath(dir, "v1295.oifits"); insname = "MIRCX")[1, 1]
    @test e.nv2 == sum(t -> d.tables[t].insname == "MIRCX", d.v2_tab)

    # Bins: every point of every bin; every bin's points tagged with their tables.
    b = rd(joinpath(data, "V1295_Aql.oifits"); polychromatic = true, merge_oi_wavelength = true)
    c = rd(write_oifits(joinpath(dir, "poly.oifits"), b; overwrite = true); polychromatic = true,
           merge_oi_wavelength = true)
    count(bins) = Dict(l => sum(x.nv2 for x in vec(bins) if x.nv2 > 0 && x.v2_lam[1] == l)
                       for l in unique(x.v2_lam[1] for x in vec(bins) if x.nv2 > 0))
    @test count(b) == count(c)
    @test sum(x.nt3phi for x in b) == sum(x.nt3phi for x in c)
    for x in b
        @test length(x.v2_tab) == x.nv2 && length(x.t3_tab) == length(x.t3phi)
        @test all(t -> x.tables[t].extname == "OI_VIS2", x.v2_tab)
    end

    # Averaged channels keep the tag: one per averaged point, its table's.
    a = rd(joinpath(data, "V1295_Aql.oifits"); binning = 3, redundance_remove = false)[1, 1]
    @test length(a.v2_tab) == length(a.v2) && all(t -> a.tables[t].extname == "OI_VIS2", a.v2_tab)
    @test sort(unique(a.tables[t].insname for t in a.v2_tab)) ==
          sort(unique(d.tables[t].insname for t in d.v2_tab))

    # Merged files: every file's tables, each point still pointing at its own.
    f1, f2 = joinpath(data, "AZ_Cyg_2011.oifits"), joinpath(data, "AZ_Cyg_2014.oifits")
    m = rd([f1, f2])[1, 1]; m1 = rd(f1)[1, 1]; m2 = rd(f2)[1, 1]
    @test m.tables == vcat(m1.tables, m2.tables)
    @test m.v2_tab == vcat(m1.v2_tab, m2.v2_tab .+ length(m1.tables))
    @test m.t3_tab == vcat(m1.t3_tab, m2.t3_tab .+ length(m1.tables))

    # A filtered dataset writes what was kept.
    g = rd(joinpath(data, "Beauty_Contest_2004", "2004-data1.oifits"))[1, 1]
    h = filter_data(g, set_data_filter(g; baseline_range = [0.0, 60e6]))
    @test 0 < h.nv2 < g.nv2
    k = rd(write_oifits(joinpath(dir, "filt.oifits"), h; overwrite = true))[1, 1]
    @test same(h, k, :v2) && same(h, k, :t3phi)
    @test length(h.v2_tab) == h.nv2 && length(h.t3_tab) == length(h.t3phi)
    @test h.v2_tab == g.v2_tab[setdiff(1:g.nv2, set_data_filter(g; baseline_range = [0.0, 60e6])[3])]

    # An existing file is not replaced unless asked.
    @test_throws ArgumentError write_oifits(joinpath(dir, "filt.oifits"), h)
end
