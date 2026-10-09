# set_data_filter's station, baseline and triangle cuts, on the 2004 beauty-contest file
# (NPOI, six stations): what is removed is exactly what uses the named station, baseline or
# triangle, and everything else is kept.
@testset "set_data_filter: stations, baselines, triangles" begin
    d = readoifits(joinpath(@__DIR__, "..", "demos", "data", "Beauty_Contest_2004",
                            "2004-data1.oifits"); warn = false, verbose = false)[1, 1]
    id = Dict(n => i for (n, i) in zip(d.sta_name, d.sta_index))
    uses(c, s) = id[s] in c
    pair(c, a, b) = sort(collect(c)) == sort([id[a], id[b]])

    # A station: every V² and T3 point using it, and nothing else.
    idx = set_data_filter(d; exclude_stations = ["E02"])
    v2_e02 = findall(c -> uses(c, "E02"), collect(eachcol(d.v2_sta_index)))
    t3_e02 = findall(c -> uses(c, "E02"), collect(eachcol(d.t3_sta_index)))
    @test !isempty(v2_e02) && !isempty(t3_e02)
    @test sort(idx[3]) == v2_e02 && sort(idx[4]) == t3_e02
    f = filter_data(d, idx)
    @test f.nv2 == d.nv2 - length(v2_e02) && f.nt3phi == length(d.t3phi) - length(t3_e02)
    @test !any(c -> uses(c, "E02"), eachcol(f.v2_sta_index))

    # A baseline, in either order: its V² points, and the triangles that contain it.
    for b in ("AC0-AE0", "AE0-AC0")
        idx = set_data_filter(d; exclude_baselines = [b])
        @test sort(idx[3]) == findall(c -> pair(c, "AC0", "AE0"), collect(eachcol(d.v2_sta_index)))
        @test sort(idx[4]) == findall(c -> uses(c, "AC0") && uses(c, "AE0"),
                                      collect(eachcol(d.t3_sta_index)))
    end

    # A triangle, in any order: those T3 points and no V².
    c1 = d.t3_sta_index[:, 1]
    names = [d.sta_name[findfirst(==(k), d.sta_index)] for k in c1]
    idx = set_data_filter(d; exclude_triangles = [join(reverse(names), "-")])
    @test isempty(idx[3])
    @test sort(idx[4]) == findall(c -> sort(collect(c)) == sort(c1), collect(eachcol(d.t3_sta_index)))

    # Combined with the existing cuts, and the refusals.
    idx = set_data_filter(d; exclude_stations = ["E02"], baseline_range = [0.0, 60e6])
    @test issubset(v2_e02, idx[3])
    @test_throws ArgumentError set_data_filter(d; exclude_stations = ["Z9"])
    @test_throws ArgumentError set_data_filter(d; exclude_baselines = ["AC0"])
    @test_throws ArgumentError set_data_filter(d; exclude_triangles = ["AC0-AE0"])
end

# readoifits(...; insname): one instrument's tables only. On a merged MIRC-X + PIONIER file the
# instruments read one at a time account for every V² and T3 point of the full read.
@testset "readoifits: insname" begin
    f = joinpath(@__DIR__, "..", "demos", "data", "V1295_Aql.oifits")
    full = readoifits(f; warn = false, verbose = false)[1, 1]
    names = unique([db.insname for db in OITOOLS.load_oifits(f).instr])
    @test length(names) > 1
    parts = [readoifits(f; insname = n, warn = false, verbose = false)[1, 1] for n in names]
    @test sum(p.nv2 for p in parts) == full.nv2
    @test sum(p.nt3phi for p in parts) == full.nt3phi
    @test all(p -> p.nv2 < full.nv2, parts)
    mircx = readoifits(f; insname = "MIRCX", warn = false, verbose = false)[1, 1]
    @test 0 < mircx.nv2 < full.nv2
    @test_throws ErrorException readoifits(f; insname = "AMBER", warn = false, verbose = false)
end

# Resetting the file's flags. MWC 275 (CHARA 4T) flags 30 triangles whose values are finite, so
# only their flags remove them (most files flag NaN points, which other cuts remove anyway);
# V² flags are set by hand on five valid points.
@testset "ignore_flags" begin
    f = joinpath(@__DIR__, "..", "demos", "data", "MWC275_4T.oifits")
    raw = readoifits(f; filter_bad_data = false, warn = false, verbose = false)[1, 1]
    @test count(raw.t3_flag) == 30
    raw.v2_flag[1:5] .= true
    # set_data_filter: with the flags honoured the cut is a superset of the one without, and
    # what the flags add is exactly the flagged points.
    with    = set_data_filter(raw; filter_bad_data = true)
    without = set_data_filter(raw; filter_bad_data = true, ignore_flags = true)
    for (k, fl, n) in ((3, raw.v2_flag, 5), (4, raw.t3_flag, 30))
        @test issubset(without[k], with[k])
        extra = setdiff(with[k], without[k])
        @test length(extra) == n && all(fl[extra])
    end
    @test count(raw.v2_flag) == 5                  # the data it was given is unchanged
    # readoifits: the flagged triangles come back, and the returned flags are cleared.
    d0 = readoifits(f; warn = false, verbose = false)[1, 1]
    d1 = readoifits(f; ignore_flags = true, warn = false, verbose = false)[1, 1]
    @test length(d1.t3phi) == length(d0.t3phi) + 30
    @test !any(d1.t3_flag) && d1.nv2 == d0.nv2
end

# Error bars: oifits_prep's floors and multipliers, now for VIS and flux as well, and
# set_data_filter's caps.
@testset "error-bar floors and caps" begin
    f = joinpath(@__DIR__, "..", "demos", "data", "V1295_Aql.oifits")
    d = readoifits(f; warn = false, verbose = false)[1, 1]
    @test d.nvisamp > 0

    # V², T3amp and T3φ exactly as before (the formulas written out here), linear and quadrature.
    for quad in (false, true)
        e = deepcopy(d)
        oifits_prep(e; min_v2_err_add = 0.01, min_v2_err_rel = 0.05, v2_err_mult = 1.5,
                    min_t3amp_err_add = 0.02, min_t3amp_err_rel = 0.1, t3amp_err_mult = 2,
                    min_t3phi_err_add = 1.0, t3phi_err_mult = 1.2, quad)
        T = Float32
        old(err, val, add, rel, mult) = quad ?
            sqrt.((val .* rel).^2 .+ (mult .* err).^2 .+ T(add)^2) :
            max.(mult .* err, abs.(val .* rel) .+ add)
        @test e.v2_err ≈ old(d.v2_err, d.v2, 0.01, 0.05, 1.5) nans = true
        @test e.t3amp_err ≈ old(d.t3amp_err, d.t3amp, 0.02, 0.1, 2) nans = true
        @test e.t3phi_err ≈ (quad ? sqrt.((1.2 .* d.t3phi_err).^2 .+ 1) :
                                    max.(1.2 .* d.t3phi_err, 1.0)) nans = true
    end
    # VIS and flux: the new rows, and the defaults change nothing.
    e = deepcopy(d)
    oifits_prep(e)
    @test isequal(e.visamp_err, d.visamp_err) && isequal(e.visphi_err, d.visphi_err)
    oifits_prep(e; min_visamp_err_rel = 0.1, visphi_err_mult = 3, min_visphi_err_add = 5.0)
    ok = isfinite.(d.visamp) .& isfinite.(d.visamp_err)
    @test all(e.visamp_err[ok] .>= 0.1 .* abs.(d.visamp[ok]) .- 1f-6)
    okp = isfinite.(d.visphi_err)
    @test e.visphi_err[okp] ≈ max.(3 .* d.visphi_err[okp], 5.0)

    # Caps: exactly the points above the cap are removed.
    cap = sort(d.v2_err)[end ÷ 2]                   # the median, without Statistics
    idx = set_data_filter(d; max_v2_err = cap)
    @test sort(idx[3]) == findall(>(cap), d.v2_err)
    idx = set_data_filter(d; max_v2_err_rel = 0.2)
    @test sort(idx[3]) == findall(i -> d.v2_err[i] > 0.2 * abs(d.v2[i]), eachindex(d.v2))
    idx = set_data_filter(d; max_t3phi_err = 5.0)
    @test sort(idx[4]) == findall(>(5.0), d.t3phi_err)
    idx = set_data_filter(d; max_visphi_err = 10.0)
    @test sort(idx[2]) == findall(>(10.0), d.visphi_err)
    @test all(isempty, set_data_filter(d))           # no cap, nothing removed
end
