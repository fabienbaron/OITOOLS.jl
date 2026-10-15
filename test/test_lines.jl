using OITOOLS, Test, ForwardDiff

# Spectral lines on a component's flux (`line_<n>_*`, PMOIRED's): the four profiles against
# PMOIRED's `Ssingle`, the flux the model uses, refusals, and the Jacobian with lines free.
@testset "spectral lines" begin
    ref = Dict{String,Any}()
    for l in eachline(joinpath(@__DIR__, "data", "pmoired_lines.txt"))
        startswith(l, "#") && continue
        f = split(l)
        if f[1] == "wl"
            ref["wl"] = parse.(Float64, f[2:end])
        elseif f[1] == "param"
            get!(ref, f[2], Dict{String,Any}())["param"] =
                Dict(String(split(x, '=')[1]) => parse(Float64, split(x, '=')[2]) for x in f[3:end])
        elseif f[1] == "flux"
            get!(ref, f[2], Dict{String,Any}())["flux"] = parse.(Float64, f[3:end])
        end
    end
    wl = ref["wl"]                                    # µm

    @testset "profiles against PMOIRED" begin
        for case in ("lor", "gau", "gaup", "texp", "tlorp", "tlorm", "two")
            p = ref[case]["param"]
            m = Dict{String,Any}("s,ud" => 0.5)
            for (k, v) in p; m["s," * k] = v; end
            # the flux the model uses (the expansion into `s,ftotal`), at the same wavelengths
            fm = dict_to_model(m, String[])
            i = findfirst(==("s,ftotal"), fm.all_names)
            got = [fm.resolver.fn(Float64[], w * 1e-6, NaN)[i] for w in wl]
            # a one-sided profile is discontinuous at λ₀: the grid's point there is left out,
            # since µm → m → µm can put it on either side
            keep = [all(abs(w - v) > 1e-9 for (k, v) in p if endswith(k, "_wl0")) for w in wl]
            @test got[keep] ≈ ref[case]["flux"][keep] atol = 1e-12
        end
        # the functions themselves, in µm, where nothing is converted
        p = ref["tlorm"]["param"]
        @test [1 + line_trunclor(w, p["line_a_wl0"], p["line_a_f"], p["line_a_trunclor"]) for w in wl] ≈
              ref["tlorm"]["flux"] atol = 1e-12
    end

    @testset "the model's flux, keys, refusals" begin
        m = Dict{String,Any}("A,ud" => 1.2, "B,f" => 0.25, "B,ud" => 0.6, "B,x" => -0.9, "B,y" => 3.5,
                             "A,line_brg_wl0" => 2.166, "A,line_brg_f" => -0.5, "A,line_brg_lorentzian" => 1.0,
                             "B,line_brg_wl0" => "2.16612 * (1 + \$vB / 2.998e5)", "vB" => -20.0,
                             "B,line_brg_f" => -0.3, "B,line_brg_gaussian" => 1.0)
        @test isempty(OITOOLS.model_warnings(m, String[]))
        fm = dict_to_model(m, String[])
        tot, comps = model_to_sed(fm, Float64[], [2.0e-6, 2.166e-6])
        @test comps["A"] ≈ [1 + line_lorentzian(2.0, 2.166, -0.5, 1.0), 0.5] rtol = 1e-12
        wB = 2.16612 * (1 - 20 / 2.998e5)
        @test comps["B"][2] ≈ 0.25 + line_gaussian(2.166, wB, -0.3, 1.0) rtol = 1e-12
        @test tot ≈ comps["A"] .+ comps["B"]
        # the continuum stays a parameter of its own
        fmf = dict_to_model(m, ["B,f"])
        @test fmf.list_free_params == ["B,f"]
        bad(d, msg) = (e = try dict_to_model(d, String[]); nothing catch err; err end;
                       e isa ArgumentError && occursin(msg, e.msg))
        @test bad(delete!(copy(m), "A,line_brg_wl0"), "A,line_brg_wl0")
        @test bad(delete!(copy(m), "A,line_brg_lorentzian"), "one profile")
        @test bad(merge(copy(m), Dict{String,Any}("A,line_brg_gaussian" => 1.0)), "lorentzian and gaussian")
        @test OITOOLS.reserved_suffix("line_brg_wl0") && OITOOLS.reserved_suffix("ftotal")
        # PMOIRED's dict imports as it is
        pm = OITOOLS.pmoired_to_dict("{'A,ud': 1.2, 'A,line_brg_wl0': '2.16612*(1+\$vA/2.998e5)', 'vA': 20, " *
                                     "'A,line_brg_f': -0.5, 'A,line_brg_lorentzian': 1.0}")
        @test isempty(OITOOLS.model_warnings(pm, String[])) && dict_to_model(pm, String[]) isa OITOOLS.FlatModel
    end

    @testset "Jacobian with the lines free" begin
        m = Dict{String,Any}("A,ud" => 1.2, "B,f" => 0.25, "B,ud" => 0.6, "B,x" => -0.9, "B,y" => 3.5,
                             "A,line_brg_wl0" => "2.16612 * (1 + \$vA / 2.998e5)", "vA" => 15.0,
                             "A,line_brg_f" => -0.5, "A,line_brg_lorentzian" => 1.0,
                             "B,line_brg_wl0" => 2.1662, "B,line_brg_f" => -0.3, "B,line_brg_gaussian" => 0.8)
        free = ["vA", "A,line_brg_f", "A,line_brg_lorentzian", "B,f", "B,line_brg_f", "B,line_brg_gaussian"]
        fm = dict_to_model(m, free)
        x = Float64[m[k] for k in free]
        n = 40
        uv = hcat([[3e7 + 1e5k, -2e7 + 2e5k] for k in 1:n]...)
        wl = collect(range(2.160e-6, 2.172e-6; length = n))
        V, J = eval_model_jac(fm, x, uv; wl)
        @test V ≈ eval_model(fm, x, uv; wl) rtol = 1e-12
        Jad = ForwardDiff.jacobian(p -> (v = eval_model(fm, p, uv; wl); [real.(v); imag.(v)]), x)
        @test [real.(J); imag.(J)] ≈ Jad rtol = 1e-8
    end

    @testset "a companion with lines takes the general path" begin
        @test !isempty(OITOOLS._line_names(Dict{String,Any}("B,line_a_f" => 1.0, "B,line_a_wl0" => 2.0,
                                                             "B,line_a_lorentzian" => 1.0), "B"))
    end
end
