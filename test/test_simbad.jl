# SIMBAD, offline.
#
# The network is not a test dependency. Everything below runs against a response recorded from
# SIMBAD's TAP service (`test/references/simbad_vega_tap.csv`), which is what makes it possible
# to assert on real column names and real values without a request. The live path is exercised
# only under OITOOLS_SIMBAD_TESTS=1, and then with exactly one query — the whole point of the
# TAP rewrite is that one target costs one request.
#
# Recapture the fixture with:
#
#     julia --project=. -e 'using OITOOLS, Downloads;
#         q = replace(OITOOLS._TARGET_QUERY, "__NAME__" => "Vega");
#         Downloads.download(string(OITOOLS.SIMBAD_TAP_URL,
#             "?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=",
#             OITOOLS.simbad_percent_encode(q)), "test/references/simbad_vega_tap.csv")'

VEGA_TAP = joinpath(@__DIR__, "references", "simbad_vega_tap.csv")

@testset "SIMBAD" begin

@testset "percent encoding" begin
    # The unreserved set passes through untouched; everything else is escaped. A query that
    # escaped too little would be truncated at the first '&' or '#' by the server.
    @test OITOOLS.simbad_percent_encode("abcXYZ019-_.~") == "abcXYZ019-_.~"
    @test OITOOLS.simbad_percent_encode("a b") == "a%20b"
    @test OITOOLS.simbad_percent_encode("SELECT * FROM basic") == "SELECT%20%2A%20FROM%20basic"
    @test OITOOLS.simbad_percent_encode("id='M 31'") == "id%3D%27M%2031%27"
    @test OITOOLS.simbad_percent_encode("a&b#c+d/e") == "a%26b%23c%2Bd%2Fe"
    # Non-ASCII is encoded byte by byte, UTF-8, as RFC 3986 requires
    @test OITOOLS.simbad_percent_encode("α") == "%CE%B1"
end

@testset "CSV parsing" begin
    p = OITOOLS.simbad_parse_csv

    @test p("a,b,c\n1,2,3\n") == [["a","b","c"], ["1","2","3"]]
    @test p("a,b\r\n1,2\r\n")  == [["a","b"], ["1","2"]]      # CRLF
    @test p("a,b\n1,2")        == [["a","b"], ["1","2"]]      # no trailing newline

    # Empty fields survive as empty strings and keep the columns aligned. SIMBAD marks "no
    # value" exactly this way, so collapsing them would shift every column after the gap.
    @test p("a,b,c\n1,,3\n") == [["a","b","c"], ["1","","3"]]
    @test p("a,b,c\n,,\n")   == [["a","b","c"], ["","",""]]

    # Quoted fields: this is why the parser exists. Spectral types and identifiers contain
    # commas ("A0Va_1, B") and SIMBAD quotes them.
    @test p("a,b\n\"x,y\",2\n")        == [["a","b"], ["x,y","2"]]
    @test p("a\n\"he said \"\"hi\"\"\"\n") == [["a"], ["he said \"hi\""]]
    @test p("a,b\n\"line1\nline2\",2\n") == [["a","b"], ["line1\nline2","2"]]
    @test p("a,b\n\"\",2\n")           == [["a","b"], ["","2"]]

    # Blank lines are not rows
    @test p("a,b\n\n1,2\n\n") == [["a","b"], ["1","2"]]
    @test p("") == Vector{String}[]
end

@testset "one request covers astrometry and photometry" begin
    q = OITOOLS._TARGET_QUERY
    # A single statement: no semicolon-separated second query, and both tables in it.
    @test count("SELECT", q) == 1
    @test occursin("basic", q) && occursin("flux", q) && occursin("ident", q)
    # LEFT OUTER, so a target with no published photometry still returns its coordinates
    # instead of reading as "no such target".
    @test occursin(r"LEFT\s+OUTER\s+JOIN\s+flux"i, q)
    @test occursin("__NAME__", q)
end

@testset "record assembly" begin
    cols, rows = OITOOLS.simbad_tap_parse(read(VEGA_TAP, String))
    t = OITOOLS.simbad_record("Vega", cols, rows)

    @test t.name == "Vega"
    @test t.main_id == "* alf Lyr"
    @test t.ra  ≈ 279.234735  atol = 1e-3       # degrees, not hours
    @test t.dec ≈  38.783689  atol = 1e-3
    @test t.plx > 100                            # 130 mas
    @test startswith(t.sptype, "A0")
    @test t.otype != ""
    @test isfinite(t.pmra) && isfinite(t.pmdec)

    # Every band is a key, always. A panel iterating SIMBAD_BANDS must not have to test for
    # presence as well as for NaN.
    @test sort(collect(keys(t.mags))) == sort(SIMBAD_BANDS)
    @test t.mags["V"] ≈ 0.03  atol = 0.01
    @test t.mags["K"] ≈ 0.129 atol = 0.01
    @test t.mags["J"] < 0                        # Vega straddles zero in the near infrared
    # SIMBAD also returns U for this target. Bands outside SIMBAD_BANDS are dropped rather
    # than added as keys, so a panel iterating SIMBAD_BANDS sees exactly what it expects.
    @test !haskey(t.mags, "U")
    @test isnan(t.mags["L"]) && isnan(t.mags["M"]) && isnan(t.mags["N"])
    # Vega is the definition of zero, so an unmeasured band arriving as 0.0 would be
    # indistinguishable from a real measurement. It has to be NaN.
    @test all(v -> isnan(v) || isfinite(v), values(t.mags))
end

@testset "record assembly is order- and case-independent" begin
    cols, rows = OITOOLS.simbad_tap_parse(read(VEGA_TAP, String))
    ref = OITOOLS.simbad_record("Vega", cols, rows)

    # Columns are looked up by name, so permuting them changes nothing. This is what stops a
    # new column in _TARGET_QUERY from silently shifting ra into dec.
    perm  = reverse(1:length(cols))
    cols2 = [uppercase(cols[i]) for i in perm]
    rows2 = [[r[i] for i in perm] for r in rows]
    t = OITOOLS.simbad_record("Vega", cols2, rows2)
    @test t.main_id == ref.main_id
    @test t.ra == ref.ra && t.dec == ref.dec
    @test all(b -> isequal(t.mags[b], ref.mags[b]), SIMBAD_BANDS)
end

@testset "missing values and no match" begin
    cols = ["main_id","ra_deg","dec_deg","pmra","pmdec","plx","rv","sptype","otype","band","mag"]

    # A target with coordinates but no photometry at all: the LEFT OUTER JOIN gives one row
    # with empty flux columns.
    r = [["* tst","10.5","-46.5","","","","","M2III","","",""]]
    t = OITOOLS.simbad_record("Test", cols, r)
    @test t.ra == 10.5 && t.dec == -46.5
    @test isnan(t.pmra) && isnan(t.plx) && isnan(t.rv)
    @test t.sptype == "M2III"
    @test all(isnan, values(t.mags))

    # Southern declinations keep their sign. This is P0-2: TAP returns decimal degrees, so
    # there is no sexagesimal reduction left to get wrong, and the test says so.
    @test OITOOLS.simbad_record("Test", cols,
              [["x","0.0","-0.5","","","","","","","",""]]).dec == -0.5

    # An unknown band in the answer is ignored rather than added as a key
    r2 = [["x","1","2","","","","","","","u_","12.3"],
          ["x","1","2","","","","","","","K","1.5"]]
    t2 = OITOOLS.simbad_record("Test", cols, r2)
    @test !haskey(t2.mags, "u_")
    @test t2.mags["K"] == 1.5

    # No rows means SIMBAD does not know the name — a data answer, not a failure of the query
    @test_throws ErrorException OITOOLS.simbad_record("NoSuchStarXYZ", cols, Vector{String}[])
end

@testset "argument checking happens before the network" begin
    @test_throws ArgumentError simbad_target("")
    @test_throws ArgumentError simbad_target("   ")
end

# The live path. One request, and only when explicitly asked for.
if get(ENV, "OITOOLS_SIMBAD_TESTS", "0") == "1"
    @testset "live query" begin
        t = simbad_target("Vega")
        @test t.main_id == "* alf Lyr"
        @test t.ra  ≈ 279.2347 atol = 1e-3
        @test t.dec ≈  38.7837 atol = 1e-3
        @test isfinite(t.mags["V"])
        @test sort(collect(keys(t.mags))) == sort(SIMBAD_BANDS)
    end
else
    @info "SIMBAD live query skipped; set OITOOLS_SIMBAD_TESTS=1 to run it"
end

end   # SIMBAD
