# load_oifits on files that break the standard: each is read table by table, the
# failing table only worked around, and every workaround said. The shipped examples are one
# GRAVITY file of PMOIRED's o Leo (OI_ARRAY without FOV, OI_FLUX with FLUX) and one of HD 210763
# (a ZIP archive; revisions 5 and 3; OI_VIS without VISPHI; TARGET_IDs not in OI_TARGET).
import AstroFITS
const p7zip_jll = OITOOLS.p7zip_jll

@testset "load_oifits: non-standard files" begin
    data = joinpath(@__DIR__, "..", "demos", "data")
    oleo = joinpath(data, "o_Leo", "GRAVI.2022-02-28T042658.272_singlescivis_singlesciviscalibrated.fits")
    hd   = joinpath(data, "HD210763", "2022-08-03-calibrated_cleaned.fits")
    rd(f; kw...) = readoifits(f; T = Float64, warn = false, verbose = false, kw...)[1, 1]
    notes_of(f) = (n = String[]; OITOOLS.load_oifits(f; notes = n); n)

    @testset "GRAVITY (o Leo)" begin
        n = notes_of(oleo)
        @test n == ["OI_ARRAY: read as revision 1 (declared 2: no FOV column)",
                    "OI_FLUX GRAVITY_FT: read as revision 1, FLUX read as FLUXDATA",
                    "OI_FLUX GRAVITY_SC: read as revision 1, FLUX read as FLUXDATA"]
        d = rd(oleo)
        @test d.nv2 > 9000 && d.nt3phi > 6000 && d.nvisphi > 9000 && d.nflux == 6536
        @test Set(t.insname for t in d.tables) == Set(["GRAVITY_SC", "GRAVITY_FT"])
        @test rd(oleo; insname = "GRAVITY_FT").nv2 == 36
        # the flux is the file's FLUX column, every finite value of it
        raw = rd(oleo; filter_bad_data = false, redundance_remove = false)
        col = AstroFITS.FitsFile(f -> vcat([vec(read(h, "FLUX")) for h in f if h.extname == "OI_FLUX"]...), oleo)
        @test sort(filter(isfinite, raw.flux)) == sort(filter(isfinite, col))
        # said once, when warnings are on
        @test_logs (:warn, r"does not follow the OIFITS standard") match_mode = :any readoifits(oleo; verbose = false)
    end

    @testset "OICLEAN, zipped (HD 210763)" begin
        n = notes_of(hd)
        @test n[1] == "a ZIP archive; read its member HD210763/2022-08-03-calibrated_cleaned.fits"
        @test "OI_ARRAY: read as revision 1 (declared 5: a revision OIFITS does not define)" in n
        @test count(contains("no VISPHI/VISPHIERR column (left missing)"), n) == 2
        @test count(contains("rows dropped, their TARGET_ID is not in OI_TARGET"), n) == 6
        @test count(x -> startswith(x, "OI_FLUX") && occursin("skipped", x), n) == 2
        d = rd(hd)
        @test d.nv2 > 3000 && d.nt3phi > 1500 && d.nvisamp > 3000
        @test d.nvisphi == 0 && d.nflux == 0
        @test length(oifits_headers(hd)) == 13
        # a repaired copy is a plain, standard file holding the same points
        out = repair_oifits(hd, joinpath(mktempdir(), "hd.oifits"))
        @test !OITOOLS._is_zip(out) && isempty(notes_of(out))
        @test (rd(out; redundance_remove = false).nv2, rd(out; redundance_remove = false).nvisamp) ==
              (rd(hd; redundance_remove = false).nv2, rd(hd; redundance_remove = false).nvisamp)
    end

    @testset "standard and old files need nothing" begin
        # 0-based TARGET_ID and STA_INDEX (old MIRC, NPOI 2004) are consistent within the file
        for f in (joinpath(data, "Bet_Lyr6T.oifits"), joinpath(data, "Iota_Peg6T.oifits"),
                  joinpath(data, "Beauty_Contest_2004", "2004-data1.oifits"), joinpath(data, "V1295_Aql.oifits"))
            @test isempty(notes_of(f))
        end
        # polaris: its OI_ARRAY tables only, every other table at the revision it declares
        n = notes_of(joinpath(data, "polaris.oifits"))
        @test !isempty(n) && all(startswith("OI_ARRAY: read as revision 1"), n)
    end

    @testset "ZIP archives that are not one FITS file" begin
        dir = mktempdir()
        two = joinpath(dir, "two.fits")
        run(pipeline(`$(p7zip_jll.p7zip()) a -tzip $two $(joinpath(data, "AX_Cir.oifits")) $(joinpath(data, "Pi1_Gru.oifits"))`;
                     stdout = devnull))
        @test OITOOLS._is_zip(two)
        @test_throws ArgumentError OITOOLS.load_oifits(two)
        bad = joinpath(dir, "bad.fits"); write(bad, UInt8['P', 'K', 3, 4, 0, 0, 0, 0])
        @test_throws ArgumentError OITOOLS.load_oifits(bad)
    end
end
