# Type-stable model evaluation (PLAN_ModelJacobian.md, J1): the resolver is the one dynamic
# call, its output concretely typed, everything after it inferred. Speed and JET counts are
# measured by test/bench/model_jacobian.jl; this pins the properties they depend on.
@testset "typed model evaluation" begin
    d = readoifits(joinpath(@__DIR__, "..", "demos", "data", "AX_Cir.oifits");
                   T = Float64, warn = false, verbose = false)[1, 1]
    uv, wl, mjd = d.uv, d.uv_lam, d.uv_mjd

    @testset "the return type is inferred" begin
        for (md, free) in (
                (Dict{String,Any}("s,ud" => 0.83, "s,f" => "1 - \$c,f", "c,f" => 0.01,
                                  "c,x" => 6.2, "c,y" => -28.5), ["s,ud", "c,f"]),
                (Dict{String,Any}("s,ud" => 0.9, "s,f" => "\$fs * (\$WL/1.6e-6)^(-4)", "fs" => 0.6,
                                  "d,fwhm" => 3.0, "d,f" => "1 - \$fs"), ["fs"]),
                (Dict{String,Any}("r,diamin" => 2.0, "r,diamout" => 3.0, "r,incl" => 40.0,
                                  "r,pa" => 30.0, "spatial_kernel" => 0.3), ["r,incl"]))
            m = dict_to_model(md, free)
            for T in (Float64, Float32)
                x = T.([md[k] for k in free])
                u = T.(uv)
                @test (@inferred eval_model(m, x, u; wl = T.(wl), mjd)) isa Vector{Complex{T}}
            end
        end
    end

    @testset "only what depends on \$WL/\$MJD is per point" begin
        md = Dict{String,Any}("s,ud" => "2 * \$b", "b" => 0.45,
                              "s,f" => "\$fs * (\$WL/1.6e-6)^(-4)", "fs" => 0.6,
                              "g,fwhm" => 3.0, "g,f" => "1 - \$s,f")
        m = dict_to_model(md, ["b", "fs"])
        r = m.resolver
        slot(k) = r.slot[findfirst(==(k), r.all_names)]
        @test slot("s,f") < 0 && slot("g,f") < 0          # chromatic, and what refers to it
        @test slot("s,ud") > 0 && slot("b") > 0 && slot("g,fwhm") > 0
        @test m.chromatic == [true, true]                  # both components read a chromatic flux
        # The typed values are the ones the vector-returning resolver gives.
        x = [0.45, 0.6]
        s = zeros(r.nscalar); cols = [zeros(size(uv, 2)) for _ in 1:r.nchromatic]
        r.fill!(s, cols, x, wl, mjd)
        v = r.fn(x, wl, mjd)
        for (i, k) in enumerate(r.all_names)
            @test (slot(k) > 0 ? s[slot(k)] : cols[-slot(k)]) == v[i]
        end
    end

    @testset "a radial profile's parameters must be scalars" begin
        md = Dict{String,Any}("d,profile" => "exp(-(\$R/\$s)^2)", "d,s" => "1.0 * \$WL / 1.6e-6",
                              "d,udout" => 6.0)
        err = try dict_to_model(md, String[]); nothing catch e; e end
        @test err isa ErrorException && occursin("radial profile", err.msg)
    end
end
