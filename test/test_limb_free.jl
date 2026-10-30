# test_limb_free.jl — the free-form limb-darkened disc (`ldspline`, `ldpiecewise`): its
# visibility against the laws it can represent exactly and against direct quadrature, its
# closed-form derivatives against finite differences, and `ld_free!` / `ld_profile`.

using OITOOLS, Test, FiniteDifferences, QuadGK, SpecialFunctions
const LF = OITOOLS

@testset "free-form limb darkening" begin
    ρ = collect(range(1e5, 2.5e8; length = 300))      # ζ up to ~30 at 8 mas
    linear_nodes(n, u) = [1 - u * (1 - m) for m in LF.ld_nodes(n)][1:end-1]

    @testset "a linear law is represented exactly" begin
        V = LF._vis_ldlin([8.0, 0.4], ρ)
        for K in (:ldspline, :ldpiecewise), n in (2, 3, 5, 9)
            @test maximum(abs.(LF._vis_ldfree(K, vcat(8.0, linear_nodes(n, 0.4)), ρ) .- V)) < 1e-9
        end
    end

    @testset "against direct quadrature" begin
        # V(ζ) = ∫ I J₀(ζ√(1-μ²)) μ dμ / ∫ I μ dμ, with I from `ld_profile`
        for (gk, nodes) in (("ldspline", [0.3, 0.55, 0.7, 0.92]), ("ldpiecewise", [0.2, 0.8, 0.7]))
            md = Dict{String,Any}("s,$gk" => 8.0, "s,f" => 1.0)
            for (k, v) in enumerate(nodes); md["s,ld$k"] = v; end
            I(μ) = LF.ld_profile(md, "s", [μ])[1]
            F = quadgk(μ -> I(μ) * μ, LF.ld_nodes(length(nodes) + 1)...; rtol = 1e-12)[1]
            V = LF._vis_ldfree(Symbol(gk), vcat(8.0, nodes), ρ)
            for i in (1, 60, 150, 300)
                ζ = π * 8.0 * ρ[i] / LF.MAS2RAD
                ref = quadgk(μ -> I(μ) * besselj0(ζ * sqrt(1 - μ^2)) * μ,
                             LF.ld_nodes(length(nodes) + 1)...; rtol = 1e-12)[1] / F
                @test V[i] ≈ ref atol = 1e-10
            end
        end
    end

    @testset "past the table (ζ ≈ 900)" begin
        ρb = collect(range(1.0e8, 2.0e8; length = 5))
        V = LF._vis_ldfree(:ldspline, vcat(300.0, linear_nodes(5, 0.4)), ρb)
        @test maximum(abs.(V .- LF._vis_ldlin([300.0, 0.4], ρb))) < 1e-12
    end

    @testset "derivatives" begin
        for K in (:ldspline, :ldpiecewise)
            p0 = [8.1, 0.45, 0.62, 0.8, 0.93]
            n = length(ρ); W = zeros(n); D = zeros(n, 5); Dρ = zeros(n)
            LF._vis_jac!(W, D, Dρ, Val(K), p0, ρ)
            Jfd, = FiniteDifferences.jacobian(central_fdm(5, 1), p -> LF._vis_ldfree(K, p, ρ), p0)
            @test W ≈ LF._vis_ldfree(K, p0, ρ)
            @test maximum(abs.(D .- Jfd)) < 1e-8
            Jρ = FiniteDifferences.jacobian(central_fdm(5, 1),
                                            r -> LF._vis_ldfree(K, p0, r .* ρ), [1.0])[1]
            @test maximum(abs.(Dρ .* ρ .- vec(Jρ))) < 1e-8
            @test LF.ForwardDiff.jacobian(p -> LF._vis_ldfree(K, p, ρ), p0) ≈ Jfd atol = 1e-8
            # a chromatic parameter (one value per point) reads the same as a scalar
            pc = Any[fill(8.1, n), 0.45, fill(0.62, n), 0.8, 0.93]
            @test LF._vis_ldfree(K, pc, ρ) ≈ LF._vis_ldfree(K, p0, ρ)
        end
        # Float32 stays Float32
        @test eltype(LF._vis_ldfree(:ldspline, Float32[8, 0.4, 0.6, 0.8, 0.9], Float32.(ρ))) == Float32
    end

    @testset "ld_profile" begin
        μ = [0.0, 0.3, 1.0]
        for (md, want) in (
                (Dict{String,Any}("s,ud" => 5.0), [1, 1, 1]),
                (Dict{String,Any}("s,ldlin" => 5.0, "s,u" => 0.4), [0.6, 0.72, 1]),
                (Dict{String,Any}("s,ldquad" => 5.0, "s,u" => 0.4, "s,w" => 0.1), [0.5, 0.671, 1]),
                (Dict{String,Any}("s,ldpow" => 5.0, "s,alpha" => 0.5), [0, sqrt(0.3), 1]),
                (Dict{String,Any}("s,ldpiecewise" => 5.0, "s,ld1" => 0.5, "s,ld2" => 0.9),
                 [0.5, 0.74, 1]))
            @test LF.ld_profile(md, "s", μ) ≈ want atol = 1e-12
        end
        # expressions are evaluated
        md = Dict{String,Any}("s,ldlin" => 5.0, "u0" => 0.2, "s,u" => "\$u0 * 2")
        @test LF.ld_profile(md, "s", [0.0])[1] ≈ 0.6
    end

    @testset "ld_free!" begin
        md = Dict{String,Any}("s,ldlin" => 8.4, "s,u" => 0.3, "s,f" => 1.0)
        free = ["s,ldlin", "s,u"]; lb = Dict("s,ldlin" => 1.0, "s,u" => -1.0)
        ub = Dict("s,ldlin" => 20.0, "s,u" => 1.0)
        ld_free!(md, free, "s"; lb, ub)
        @test md["s,ldspline"] == 8.4 && !haskey(md, "s,ldlin") && !haskey(md, "s,u")
        @test free == ["s,ldspline", "s,ld1", "s,ld2", "s,ld3", "s,ld4"]
        @test lb["s,ldspline"] == 1.0 && ub["s,ldspline"] == 20.0 && !haskey(lb, "s,u")
        @test [md["s,ld$k"] for k in 1:4] ≈ [0.7, 0.775, 0.85, 0.925]
        # the visibility is the linear law's
        fm = LF.parse_model(md, free)
        uv = [ρ'; zeros(1, length(ρ))]
        @test abs.(eval_model(fm, [md[k] for k in free], uv)) ≈ abs.(LF._vis_ldlin([8.4, 0.3], ρ)) atol = 1e-9

        # another node count, the piecewise basis, monotonic: the profile carried over
        ld_free!(md, free, "s"; n = 3, basis = :piecewise, monotonic = true, lb, ub)
        @test haskey(md, "s,ldpiecewise") && !haskey(md, "s,ldspline") && !haskey(md, "s,ld3")
        @test free == ["s,ldpiecewise", "s,ldq1", "s,ldq2"]
        @test (lb["s,ldq1"], ub["s,ldq1"]) == (0.0, 1.0) && !haskey(lb, "s,ld1")
        @test LF.ld_profile(md, "s", [0.0, 0.5, 1.0]) ≈ [0.7, 0.85, 1.0]
        # and back, the nodes numeric again
        ld_free!(md, free, "s"; n = 3, lb, ub)
        @test md["s,ld1"] ≈ 0.7 && !haskey(md, "s,ldq1") && free == ["s,ldspline", "s,ld1", "s,ld2"]
        @test_throws ArgumentError ld_free!(Dict{String,Any}("g,fwhm" => 2.0), String[], "g")
    end

    @testset "keys: no warning, bounds, a model file" begin
        md = Dict{String,Any}("s,ldspline" => 8.0, "s,ldq1" => 0.6, "s,ldq2" => 0.9,
                              "s,ld1" => "\$s,ldq1 * \$s,ld2", "s,ld2" => "\$s,ldq2", "s,f" => 1.0)
        @test isempty(LF._warn_unrecognised_keys(md))
        lb, ub = default_bounds(md, ["s,ldspline", "s,ldq1", "s,ld2"])
        @test (lb["s,ldq1"], ub["s,ldq1"]) == (0.0, 1.0)
        @test (lb["s,ld2"], ub["s,ld2"]) == (0.0, 2.0)
        @test LF.reserved_suffix("ld7") && LF.reserved_suffix("ldq2")
        path = tempname() * ".toml"
        write_model_file(path, md; free = ["s,ldspline", "s,ldq1", "s,ldq2"])
        m = read_model_file(path)
        @test m.model == md
        fm = LF.parse_model(m.model, m.free)
        @test length(fm.components[1].param_idx) == 3
    end
end
