# The Tools tab's file operations: the header inspector (oifits_headers), the validator and its
# repaired copy (validate_oifits, repair_oifits), the splits (split_by_time,
# split_by_instrument, split_by_target) and a FITS image's keywords (fits_set_keywords).
import AstroFITS

@testset "OIFITS tools" begin
    data = joinpath(@__DIR__, "..", "demos", "data")
    rd(f; kw...) = readoifits(f; T = Float64, warn = false, verbose = false, kw...)
    dir = mktempdir()
    sev(r, s) = [x for x in r if x.severity == s]

    @testset "oifits_headers" begin
        h = oifits_headers(joinpath(data, "V1295_Aql.oifits"))
        @test h[1].kind == :primary && all(x -> x.kind == :table, h[2:end])
        wav = [x for x in h if x.extname == "OI_WAVELENGTH"]
        @test "MIRCX" in [x.insname for x in wav] && all(x -> x.revn >= 1, wav)
        v2 = first(x for x in h if x.extname == "OI_VIS2")
        @test v2.nrows > 0 && "VIS2DATA" in [c.name for c in v2.columns]
        @test first(c for c in h[1].cards).name == "SIMPLE"
        # an image: its dimensions
        img = joinpath(dir, "img.fits"); AstroFITS.writefits!(img, ["OBJECT" => "x"], rand(Float32, 8, 6))
        hi = oifits_headers(img)
        @test hi[1].dims == [8, 6] && any(c -> c.name == "OBJECT" && c.value == "x", hi[1].cards)
        # every HDU of a GRAVITY file, the ones OITOOLS does not read included
        oleo = joinpath(homedir(), "SOFTWARE", "PMOIRED_examples", "DATA", "o_Leo")
        if isdir(oleo)
            g = oifits_headers(first(filter(endswith(".fits"), readdir(oleo; join = true))))
            @test length(g) == 14 && g[end].extname == "TELLURICS"
            @test Set(x.insname for x in g if x.extname == "OI_WAVELENGTH") == Set(["GRAVITY_SC", "GRAVITY_FT"])
        end
        zip = joinpath(dir, "zipped.fits"); write(zip, UInt8['P', 'K', 3, 4, 0, 0, 0, 0])
        @test_throws ArgumentError oifits_headers(zip)
    end

    @testset "validate_oifits" begin
        @test isempty(sev(validate_oifits(joinpath(data, "AX_Cir.oifits")), :error))
        r = validate_oifits(joinpath(data, "AlphaCenA.oifits"))
        @test any(x -> occursin("no telescope or station name", x.message), sev(r, :warning))
        # 0-based TARGET_ID and STA_INDEX (old MIRC): OIFITS v2 numbers both from 1
        r = validate_oifits(joinpath(data, "Bet_Lyr6T.oifits"))
        @test any(x -> startswith(x.message, "TARGET_ID 0"), r) && any(x -> startswith(x.message, "STA_INDEX 0"), r)
        @test !any(x -> startswith(x.message, "TARGET_ID"), validate_oifits(joinpath(data, "V1295_Aql.oifits")))
        r = validate_oifits(joinpath(data, "polaris.oifits"))
        @test any(x -> occursin("read as revision 1", x.message), sev(r, :warning))
        # its OI_TARGET is revision 2 and numbers its target 0: the one error
        @test [x.message for x in sev(r, :error)] == ["TARGET_ID 0: OIFITS v2 numbers targets from 1"]
        @test any(x -> occursin("ZIP", x.message), validate_oifits(joinpath(dir, "zipped.fits")))
        # OI_VIS tables that disagree on PHITYP
        p = joinpath(dir, "phityp.oifits")
        write_oifits(p, rd(joinpath(data, "V1295_Aql.oifits"))[1, 1]; overwrite = true)
        AstroFITS.FitsFile(p, "rw") do f
            h = first(h for h in f if h.extname == "OI_VIS"); h["PHITYP"] = "differential"
        end
        @test any(x -> occursin("PHITYP", x.message), sev(validate_oifits(p), :warning))
        oleo = joinpath(homedir(), "SOFTWARE", "PMOIRED_examples", "DATA", "o_Leo")
        if isdir(oleo)
            r = validate_oifits(first(filter(endswith(".fits"), readdir(oleo; join = true))))
            @test count(x -> occursin("FLUX read as FLUXDATA", x.message), sev(r, :warning)) == 2
            @test isempty(sev(r, :error))
            @test any(x -> startswith(x.table, "TELLURICS"), sev(r, :info))
        end
    end

    @testset "repair_oifits" begin
        src = joinpath(data, "polaris.oifits")
        out = repair_oifits(src, joinpath(dir, "polaris.oifits"))
        @test isempty(sev(validate_oifits(out), :error)) && isempty(sev(validate_oifits(out), :warning))
        a = rd(src; redundance_remove = false)[1, 1]; b = rd(out; redundance_remove = false)[1, 1]
        @test (a.nv2, a.nt3amp, a.nt3phi) == (b.nv2, b.nt3amp, b.nt3phi)
        @test sort(a.v2) == sort(b.v2) && sort(a.t3phi) == sort(b.t3phi)
        @test_throws ArgumentError repair_oifits(src, src)
        @test_throws ArgumentError repair_oifits(joinpath(data, "HD140573.oifits"), joinpath(dir, "x.oifits"))
        @test isfile(repair_oifits(joinpath(data, "MWC275_4T.oifits"), joinpath(dir, "mwc.oifits");
                                   targetname = "MWC 275"))
    end

    @testset "splits" begin
        for T in (Float32, Float64)
            d = readoifits(joinpath(data, "V1295_Aql.oifits"); T, warn = false, verbose = false)[1, 1]
            parts, mjd = split_by_time(d)
            @test length(parts) == 6 && issorted(mjd)
            @test sum(p.nv2 for p in parts) == d.nv2 && sum(length(p.t3phi) for p in parts) == length(d.t3phi)
            @test all(p -> eltype(p.v2) == T, parts)
            # fixed windows: the same parts at either precision (points go by their Float64 MJDs)
            @test length(split_by_time(d; binwidth = 1 / 24)[1]) == 15
            w = split_by_time(d; edges = [(58612.0, 58613.0), (58639.0, 58644.0)])[1]
            @test length(w) == 2 && all(p -> all(58612.0 .<= p.v2_mjd .<= 58644.0), w)
        end
        # epochs come from every observable: closure phases with no V² at their epoch are kept
        d = rd(joinpath(data, "V1295_Aql.oifits"))[1, 1]
        late = findall(>(58650.0), d.v2_mjd)
        e = filter_data(d, [Int64[], Int64[], late, Int64[], Int64[]])
        parts, _ = split_by_time(e)
        @test length(parts) == 6 && sum(length(p.t3phi) for p in parts) == length(d.t3phi)
        @test any(p -> p.nv2 == 0 && length(p.t3phi) > 0, parts)
        # by instrument, from the table tags
        s = split_by_instrument(rd(joinpath(data, "MWC275_4T.oifits"))[1, 1])
        @test first.(s)[1] == "MIRCX" && length(s) == 3
        @test all(all(t -> x.tables[t].insname == n, x.v2_tab) for (n, x) in s)
        @test sum(x.nv2 for x in last.(s)) == rd(joinpath(data, "MWC275_4T.oifits"))[1, 1].nv2
        # by target, re-read
        t = split_by_target(joinpath(data, "MWC275_4T.oifits"); T = Float64, warn = false, verbose = false)
        @test first.(t) == ["MWC 275", "HD163296"] && t[1][2][1, 1].nv2 > 0
    end

    @testset "fits_set_keywords" begin
        # a cube with no axis keywords: a disc brightening with wavelength, moving with time
        cube = Float32[hypot(x - 8 - t, y - 8) < 5 ? w + t : 0.01 for x in 1:16, y in 1:16, w in 1:3, t in 1:2]
        bare = joinpath(dir, "bare.fits")
        AstroFITS.writefits!(bare, ["OBJECT" => ("disc", "a test cube")], cube)
        @test_throws ArgumentError read_model_image(bare)                 # no axis keywords
        wl = collect(range(1.5e-6, step = 0.1e-6, length = 3))
        mj = collect(range(58000.5, step = 2.0, length = 2))
        out = fits_set_keywords(bare, joinpath(dir, "keys.fits"); pixsize = 0.25,
                                wavelengths = wl, mjds = mj)
        m = read_model_image(out)
        @test m.wavelengths ≈ wl && m.mjds ≈ mj && m.pixsize ≈ 0.25
        @test m.image ≈ cube ./ sum(cube; dims = (1, 2))                 # planes, normalised
        h = first(OITOOLS._fits_header(out))
        @test h["OBJECT"] == "disc" && h["CTYPE3"] == "WAVE" && h["CTYPE4"] == "TIME"
        @test eltype(OITOOLS.readfits(out)) == Float32                   # the file's own type
        # time on axis 3: the axes are swapped into (wavelength, time)
        swapped = joinpath(dir, "swapped.fits")
        AstroFITS.writefits!(swapped, ["OBJECT" => "x"], permutedims(cube, (1, 2, 4, 3)))
        m2 = read_model_image(fits_set_keywords(swapped, joinpath(dir, "s2.fits"); pixsize = m.pixsize,
                                                axis3 = :time, wavelengths = wl, mjds = mj))
        @test m2.image ≈ m.image && m2.mjds ≈ m.mjds
        # East on the left in the file: flipped East-first
        img = rand(Float64, 8, 6)
        el = joinpath(dir, "eastleft.fits"); AstroFITS.writefits!(el, ["OBJECT" => "x"], reverse(img; dims = 1))
        @test fits_keyword_cube(el; east_left = true)[:, :, 1, 1] == img
        # a 3-D time series; overwriting the file itself
        ts = joinpath(dir, "series.fits"); AstroFITS.writefits!(ts, ["OBJECT" => "x"], rand(8, 8, 3) .+ 0.1)
        fits_set_keywords(ts, ts; pixsize = 0.3, axis3 = :time, mjds = [1.0, 2.0, 3.0])
        @test read_model_image(ts).mjds == [1.0, 2.0, 3.0] && size(read_model_image(ts).image) == (8, 8, 1, 3)
        # the file's reference pixel is kept, and read back as the image's offset
        rp = joinpath(dir, "refpix.fits")
        AstroFITS.writefits!(rp, ["CRPIX1" => 8.0, "CRPIX2" => 9.0, "CDELT1" => 1e-9], Float32.(cube[:, :, 1, 1]))
        @test fits_refpix(rp) == (8.0, 9.0) && fits_refpix(bare) === nothing
        k = joinpath(dir, "rp_keys.fits")
        h(f) = first(OITOOLS._fits_header(f))
        fits_set_keywords(rp, k; pixsize = 0.25)
        @test (h(k)["CRPIX1"], h(k)["CRPIX2"]) == (8.0, 9.0)
        @test all(read_model_image(k).offset .≈ (0.25, 0.0))       # one pixel West of the centre 9
        fits_set_keywords(rp, k; pixsize = 0.25, refpix = nothing)
        @test (h(k)["CRPIX1"], h(k)["CRPIX2"]) == (9.0, 9.0) && read_model_image(k).offset == (0.0, 0.0)
        fits_set_keywords(rp, k; pixsize = 0.25, refpix = (5, 6))
        @test (h(k)["CRPIX1"], h(k)["CRPIX2"]) == (5.0, 6.0)
        # East on the left: the pixel numbered in the file, flipped with the array
        fits_set_keywords(rp, k; pixsize = 0.25, east_left = true)
        @test h(k)["CRPIX1"] == 9.0 && h(k)["CRPIX2"] == 9.0
        # refusals
        @test_throws ArgumentError fits_set_keywords(bare, joinpath(dir, "x.fits"); pixsize = 0.1, mjds = mj)
        @test_throws ArgumentError fits_set_keywords(bare, joinpath(dir, "x.fits"); pixsize = 0.1,
                                                     wavelengths = [1e-6], mjds = mj)
        @test_throws ArgumentError fits_set_keywords(bare, joinpath(dir, "x.fits"); pixsize = 0,
                                                     wavelengths = wl, mjds = mj)
        @test_throws ArgumentError fits_keyword_cube(bare; axis3 = :space)
    end
end
