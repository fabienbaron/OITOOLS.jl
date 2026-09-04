# Interoperability: files written through OIFITS.write carry a redundant
# TDIMn card on every scalar column, which makes astropy (and therefore
# PMOIRED) read those columns as (nrow, 1) instead of (nrow,).
# oifits_fix_tdim removes them.

using OITOOLS, Test

@testset "oifits_fix_tdim" begin
    src = joinpath(@__DIR__, "..", "demos", "data", "BC2004", "2004-data1.oifits")
    out = joinpath(mktempdir(), "written.oifits")

    model = dict_to_model(Dict{String,Any}("star,ud" => 3.0, "star,f" => 1.0), String[])
    simulate_from_oifits(src, out; flat_model=model, flat_params=Float64[])
    @test isfile(out)

    n = oifits_fix_tdim(out)
    @test n > 0                      # cfitsio wrote redundant TDIMs
    @test oifits_fix_tdim(out) == 0  # idempotent

    # the file is still a valid OIFITS for us
    d = readoifits(out; T=Float64)[1,1]
    @test d.nv2 > 0

    # and multi-element columns kept their TDIM (i.e. we did not flatten the
    # spectral axis away)
    d2 = readoifits(src; T=Float64, filter_bad_data=false)[1,1]
    @test length(readoifits(out; T=Float64, filter_bad_data=false)[1,1].v2) == length(d2.v2)

    @test_throws ErrorException oifits_fix_tdim(joinpath(@__DIR__, "no_such_file.oifits"))
end
