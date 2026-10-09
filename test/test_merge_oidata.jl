# Several files as one dataset: merge_oidata and readoifits(::Vector). The same target at two
# epochs (AZ Cyg 2011 and 2014; ι Peg with 4 and 6 telescopes): every point of every file, each
# with its own stations, and a model's χ² on the merge equal to the sum over the files.
@testset "merge_oidata: several files as one dataset" begin
    data = joinpath(@__DIR__, "..", "demos", "data")
    rd(f; kw...) = readoifits(f; T = Float64, warn = false, verbose = false, kw...)
    stations(d, c) = sort([d.sta_name[findfirst(==(k), d.sta_index)] for k in c])
    model = dict_to_model(Dict{String,Any}("star,ud" => 1.0), String[])

    for pair in (["AZ_Cyg_2011.oifits", "AZ_Cyg_2014.oifits"],
                 ["Iota_Peg4T.oifits", "Iota_Peg6T.oifits"])
        fs = joinpath.(data, pair)
        each = [rd(f)[1, 1] for f in fs]
        m = rd(fs)[1, 1]
        for f in (:nv2, :nt3amp, :nt3phi, :nvisamp, :nvisphi, :nflux, :nuv)
            @test getfield(m, f) == sum(getfield(e, f) for e in each)
        end
        @test m.v2 == reduce(vcat, [e.v2 for e in each])
        @test m.t3phi == reduce(vcat, [e.t3phi for e in each])
        @test allunique(m.sta_name) && m.sta_index == 1:length(m.sta_name)

        # Each point keeps its stations and its uv.
        o2 = 0; o3 = 0
        for e in each
            @test all(stations(e, e.v2_sta_index[:, i]) == stations(m, m.v2_sta_index[:, o2 + i])
                      for i in eachindex(e.v2))
            @test all(stations(e, e.t3_sta_index[:, i]) == stations(m, m.t3_sta_index[:, o3 + i])
                      for i in eachindex(e.t3phi))
            @test m.uv[:, m.indx_v2[o2 .+ eachindex(e.v2)]] == e.uv[:, e.indx_v2]
            @test m.uv[:, m.indx_t3_3[o3 .+ eachindex(e.t3phi)]] == e.uv[:, e.indx_t3_3]
            o2 += length(e.v2); o3 += length(e.t3phi)
        end

        w = [1, 1, 1, 0, 0, 0, 0]
        @test model_to_chi2(model, Float64[], m; weights = w) ≈
              sum(model_to_chi2(model, Float64[], e; weights = w) for e in each) rtol = 1e-12

        # The merge writes and reads back to the same points.
        out = write_oifits(joinpath(mktempdir(), "merged.oifits"), m)
        back = rd(out)[1, 1]
        @test sort(back.v2) == sort(m.v2) && sort(back.t3phi) == sort(m.t3phi)
    end

    # One file is that file; the refusals.
    f1 = joinpath(data, "AZ_Cyg_2011.oifits")
    @test rd([f1])[1, 1].v2 == rd(f1)[1, 1].v2
    a = rd(f1)[1, 1]; b = deepcopy(a)
    @test merge_oidata([a, b]).nv2 == 2a.nv2
    @test_throws ArgumentError merge_oidata(typeof(a)[])
    if !isempty(a.visamp)
        b.phityp = a.phityp == "differential" ? "absolute" : "differential"
        @test_throws ArgumentError merge_oidata([a, b])
    else
        c = rd(joinpath(data, "V1295_Aql.oifits"))[1, 1]; e = deepcopy(c)
        e.phityp = "differential"
        @test_throws ArgumentError merge_oidata([c, e])
    end
    @test_throws ArgumentError rd([joinpath(data, "V1295_Aql.oifits"), f1]; polychromatic = true)
end
