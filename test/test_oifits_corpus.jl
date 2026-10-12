# Every OIFITS file shipped in demos/data: read at the revision it declares, read by readoifits,
# and written back by write_oifits to the same points. Between them they cover PIONIER, MIRC-X
# (OIFITS v1 and v2), the old MIRC (H_PRISM), NPOI (2004), MATISSE and GRAVITY (2026) and the
# 6-way combiner, so a change in the FITS layer that loses a column, a revision or a point on any
# instrument fails here rather than in a user's file.
using OITOOLS, Test, OIFITS
import AstroFITS

@testset "the shipped OIFITS corpus" begin
    root = joinpath(@__DIR__, "..", "demos", "data")
    files = sort([joinpath(r, f) for (r, _, fs) in walkdir(root) for f in fs
                  if occursin(r"\.(oifits|fits)$"i, f) && !occursin("images", r)])
    files = filter(files) do f            # OIFITS files only: an image has no OI_ table
        AstroFITS.FitsFile(ff -> any(h -> startswith(something(h.extname, ""), "OI_"), ff), f)
    end
    @test length(files) >= 19
    rd(f; kw...) = readoifits(f; T = Float64, warn = false, verbose = false, kw...)
    srt(x) = sort(Float64.(x))
    dir = mktempdir()

    @testset "$(relpath(f, root))" for f in files
        # Each table is read at the revision its header declares. A file may break its own
        # revision (polaris' OI_ARRAY claims revision 2 without the FOV column revision 2
        # requires); then, and only then, readoifits falls back to reading it as revision 1.
        direct = try
            read(OIDataSet, f)
        catch e
            @test e isa OIFITS.MissingColumn
            nothing
        end
        if direct === nothing
            @test read(OIDataSet, f; hack_revn = 1) isa OIDataSet
        else
            declared = AstroFITS.FitsFile(f) do ff
                sort([(something(h.extname, ""), get(Int, h, "OI_REVN", 0)) for h in ff
                      if startswith(something(h.extname, ""), "OI_")])
            end
            blocks = vcat([("OI_TARGET", direct.target.revn)],
                          [("OI_ARRAY", b.revn) for b in direct.array],
                          [("OI_WAVELENGTH", b.revn) for b in direct.instr],
                          [("OI_CORR", b.revn) for b in direct.correl],
                          [("OI_VIS", b.revn) for b in direct.vis],
                          [("OI_VIS2", b.revn) for b in direct.vis2],
                          [("OI_T3", b.revn) for b in direct.t3],
                          [("OI_FLUX", b.revn) for b in direct.flux])
            @test sort(blocks) == filter(x -> x[1] != "OI_INSPOL", declared)
        end

        d = rd(f)[1, 1]
        @test d.nv2 + d.nt3phi + d.nvisamp + d.nflux > 0

        # AMPTYP / PHITYP: what the tables say is what the raw header says.
        at, pt = OITOOLS._read_vis_header_keywords(f)
        isempty(at) || @test d.amptyp == OITOOLS._canon_obstype(at)
        isempty(pt) || @test d.phityp == OITOOLS._canon_obstype(pt)

        # Written and read back: the same points. Both reads keep every uv point: merging the
        # ones within `uvtol` of each other picks a representative, and the written file need
        # not merge them the same way, which would move a point by up to `uvtol` and hide a real
        # difference behind a tolerance.
        d = rd(f; redundance_remove = false)[1, 1]
        e = rd(write_oifits(joinpath(dir, "out.oifits"), d; overwrite = true);
               redundance_remove = false)[1, 1]
        @test (e.nv2, e.nt3amp, e.nt3phi, e.nvisamp, e.nvisphi, e.nflux) ==
              (d.nv2, d.nt3amp, d.nt3phi, d.nvisamp, d.nvisphi, d.nflux)
        for fld in (:v2, :v2_err, :v2_lam, :v2_mjd, :t3amp, :t3amp_err, :t3phi, :t3phi_err,
                    :t3_lam, :visamp, :visamp_err, :visphi, :visphi_err, :flux, :flux_err)
            @test isequal(srt(getfield(d, fld)), srt(getfield(e, fld)))
        end
        # uv at the observables' own indices
        for idx in (:indx_v2, :indx_vis, :indx_t3_1, :indx_t3_2, :indx_t3_3), row in 1:2
            @test srt(d.uv[row, getfield(d, idx)]) ≈ srt(e.uv[row, getfield(e, idx)]) rtol = 1e-12
        end
        @test (d.amptyp, d.phityp, d.flux_calibrated) == (e.amptyp, e.phityp, e.flux_calibrated)
        @test sort(d.sta_name) == sort(e.sta_name)

        # Each point names the table it came from, and that table's kind.
        for (tab, n, ext) in ((:v2_tab, :v2, "OI_VIS2"), (:vis_tab, :visamp, "OI_VIS"),
                              (:t3_tab, :t3phi, "OI_T3"), (:flux_tab, :flux, "OI_FLUX")), x in (d, e)
            @test length(getfield(x, tab)) == length(getfield(x, n))
            @test all(t -> x.tables[t].extname == ext, getfield(x, tab))
        end
        # The written file keeps the instruments: the same names, each with its own channels and
        # its own points, which `insname` reads alone; and an OI_VIS table's AMPTYP / PHITYP.
        pts(x) = [(x.tables[t].insname, Float64(l)) for (tab, lam) in
                  ((x.v2_tab, x.v2_lam), (x.vis_tab, x.vis_lam), (x.t3_tab, x.t3_lam),
                   (x.flux_tab, x.flux_lam)) for (t, l) in zip(tab, lam)]
        inst(x) = sort(unique(first.(pts(x))))
        @test inst(d) == inst(e)
        for X in inst(d)
            @test sort(unique(last.(filter(p -> p[1] == X, pts(d))))) ≈
                  sort(unique(last.(filter(p -> p[1] == X, pts(e))))) rtol = 1e-6
            a = rd(f; insname = X, redundance_remove = false)[1, 1]
            b = rd(joinpath(dir, "out.oifits"); insname = X, redundance_remove = false)[1, 1]
            @test (a.nv2, a.nt3amp, a.nt3phi, a.nvisamp, a.nvisphi, a.nflux) ==
                  (b.nv2, b.nt3amp, b.nt3phi, b.nvisamp, b.nvisphi, b.nflux)
        end
        vtyp(x) = sort(unique([(x.tables[t].insname, x.tables[t].amptyp, x.tables[t].phityp)
                               for t in x.vis_tab]))
        @test vtyp(d) == vtyp(e)
        w = read(OIDataSet, joinpath(dir, "out.oifits"))
        @test w.target.revn == 2
        @test all(b.revn == 2 for b in vcat(w.array, w.instr, w.vis, w.vis2, w.t3))
    end
end
