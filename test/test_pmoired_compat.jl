# The PMOIRED import and export: what they rewrite where the two expression languages differ,
# against PMOIRED's own definitions (its oimodels.py), and a round trip.
using QuadGK, SpecialFunctions

@testset "PMOIRED import and export" begin
    @testset "** is ^, with Python's precedence" begin
        d = pmoired_to_dict("{'a,f': '-2**2', 'b,f': '2**3**2', 'c,f': '4**-0.5'}")
        @test d["a,f"] == "-2^2" && d["b,f"] == "2^3^2" && d["c,f"] == "4^-0.5"
        # Python: -2**2 = -4, 2**3**2 = 512, 4**-0.5 = 0.5
        @test eval(Meta.parse(d["a,f"])) == -4
        @test eval(Meta.parse(d["b,f"])) == 512
        @test eval(Meta.parse(d["c,f"])) == 0.5
        @test occursin("'2**3**2'", dict_to_pmoired(Dict{String,Any}("b,f" => "2^3^2"); check = false))
    end

    @testset "\$WL means the same (microns) in both" begin
        d = pmoired_to_dict("{'s,ud': 1.0, 's,f': '(\$WL/1.67)**-3.8'}")
        @test d["s,f"] == "(\$WL/1.67)^-3.8"
        m = dict_to_model(d, String[])
        @test real(only(eval_model(m, Float64[], zeros(2, 1); wl = [1.67e-6]))) ≈ 1
    end

    @testset "PMOIRED's named ring profiles" begin
        d = pmoired_to_dict("{'r,profile': 'uniform', 's,profile': 'doughnut', 't,profile': 'doughnut3'}")
        x = "(2*(\$R - (\$RMIN + \$RMAX)/2)/(\$RMAX - \$RMIN))"
        @test d["r,profile"] == "1"
        @test d["s,profile"] == "1 - $x^2"
        @test d["t,profile"] == "1 - abs($x)^3"
        # an expression is left as it is
        @test pmoired_to_dict("{'r,profile': '\$R**2'}")["r,profile"] == "\$R^2"

        # the doughnut's V against PMOIRED's formula integrated over the ring
        a, b = 1.0, 3.0                                     # radii, mas
        I(r) = 1 - ((r - (a + b) / 2) / (b - a) * 2)^2
        B = collect(range(0.02, 0.6, length = 12))         # cycles/mas
        V = [quadgk(r -> I(r) * besselj0(2π * β * r) * r, a, b; rtol = 1e-12)[1] for β in B] ./
            quadgk(r -> I(r) * r, a, b; rtol = 1e-12)[1]
        md = merge(pmoired_to_dict("{'d,profile': 'doughnut'}"),
                   Dict{String,Any}("d,diamin" => 2a, "d,diamout" => 2b, "d,f" => 1.0))
        uv = vcat(permutedims(B .* OITOOLS.MAS2RAD), zeros(1, length(B)))
        @test maximum(abs.(eval_model(dict_to_model(md, String[]), Float64[], uv) .- V)) < 1e-9
    end
end
