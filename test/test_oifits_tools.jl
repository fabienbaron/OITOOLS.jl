# The Tools tab's file operations: the header inspector (oifits_headers), the validator and its
# repaired copy (validate_oifits, repair_oifits), and the splits (split_by_time,
# split_by_instrument, split_by_target).
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
        @test any(x -> occursin("not in OI_ARRAY", x.message), sev(r, :warning))
        r = validate_oifits(joinpath(data, "polaris.oifits"))
        @test any(x -> occursin("read as revision 1", x.message), sev(r, :warning))
        @test isempty(sev(r, :error))
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
end
