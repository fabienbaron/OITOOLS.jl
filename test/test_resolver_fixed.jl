using OITOOLS, Test

# Fixed parameters are data, not code: a model compiled again with other fixed values has the
# same generated functions (nothing new to compile) and the new values.
@testset "fixed values as data" begin
    for m in (Dict{String,Any}("1,ud" => 0.5, "2,f" => 0.3, "2,ud" => 0.2, "2,x" => 1.0, "2,y" => -2.0),
              Dict{String,Any}("1,ud" => "\$d0 * (1 + 0.1 * (\$WL - 1.6))", "d0" => 0.6, "2,f" => 0.3,
                               "2,x" => 1.0, "2,y" => -2.0),
              Dict{String,Any}("1,ud" => 0.3, "2,f" => 0.3, "2,x" => "orbit", "2,y" => "orbit",
                               "2,orb P" => 30.0, "2,orb MJD0" => 58640.0, "2,orb e" => 0.4,
                               "2,orb omega" => 60.0, "2,orb OMEGA" => 120.0, "2,orb incl" => 50.0,
                               "2,orb a" => 4.0))
        free = ["2,f"]
        a = dict_to_model(m, free)
        k = first(sort([k for (k, v) in m if v isa Real && k ∉ free]))
        m2 = merge(copy(m), Dict{String,Any}(k => m[k] + 0.25))
        b = dict_to_model(m2, free)
        @test typeof(a.resolver.fill!) == typeof(b.resolver.fill!)
        @test typeof(a.resolver.fn) == typeof(b.resolver.fn)
        i = findfirst(==(k), a.all_names)
        @test a.resolver.fn([0.3])[i] == m[k] && b.resolver.fn([0.3])[i] == m[k] + 0.25
        # evaluations follow the values: the second model at its own fixed value equals the first
        # model freed and evaluated there
        uv = [3e7 1e7; -2e7 4e7]; t = [58650.0, 58655.0]
        c = dict_to_model(m2, [free; k])
        @test eval_model(b, [0.3], uv; wl = [1.6e-6, 1.7e-6], mjd = t) ==
              eval_model(c, [0.3, m[k] + 0.25], uv; wl = [1.6e-6, 1.7e-6], mjd = t)
    end

    # globals named like the generated functions' arguments do not overwrite them
    for g in ("x", "s", "c", "cols")
        m = Dict{String,Any}("1,ud" => "\$$g * 0.5", g => 1.2, "2,f" => 0.3, "2,x" => 1.0, "2,y" => -2.0)
        uv = [3e7 1e7; -2e7 4e7]
        for free in ([g, "2,f"], ["2,f"])
            fm = dict_to_model(m, free)
            x = Float64[m[k] for k in free]
            ref = eval_model(dict_to_model(merge(copy(m), Dict{String,Any}("1,ud" => 0.6)), free), x, uv)
            @test eval_model(fm, x, uv) ≈ ref
            @test eltype(eval_model(fm, Float32.(x), Float32.(uv))) == ComplexF32
        end
    end

    # the fixed values in the parameters' type: Float32 parameters, Float32 values throughout
    m = Dict{String,Any}("1,ud" => "\$d0 * \$k", "d0" => 0.6, "k" => 1.0, "2,f" => 0.3, "2,x" => 1.0, "2,y" => -2.0)
    fm = dict_to_model(m, ["k"])
    n2i = Dict(n => i for (i, n) in enumerate(fm.all_names))
    v32 = fm.resolver.fn(Float32[1.1]); v64 = fm.resolver.fn([1.1])
    @test v32[n2i["d0"]] isa Float32 && v32[n2i["1,ud"]] isa Float32 && v32[n2i["2,x"]] isa Float32
    @test v64[n2i["d0"]] isa Float64 && v64[n2i["1,ud"]] === 0.6 * 1.1
    uv = Float32[3e7 1e7; -2e7 4e7]
    V, J = eval_model_jac(fm, Float32[1.1], uv)
    @test eltype(V) == ComplexF32 && eltype(J) == ComplexF32
    @test V ≈ eval_model(fm, [1.1], Float64.(uv)) rtol = 1e-5
end
