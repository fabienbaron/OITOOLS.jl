# The closed-form model Jacobian (PLAN_ModelJacobian.md, J3) against ForwardDiff through
# `eval_model` (`eval_model_grad_ad`), on the benchmark zoo and on the geometry it does not
# cover, with finite differences on a few columns so the reference is not trusted blindly.
using FiniteDifferences

@testset "model Jacobian" begin
    O = OITOOLS
    include(joinpath(@__DIR__, "bench", "zoo.jl"))

    # Column by column: |J − J_ad| ≤ rtol · max|J_ad[:, k]| (Float32: or 1e-6, if looser).
    function agrees(J, Ja, T)
        rtol = T == Float64 ? 1e-10 : 1e-4
        all(maximum(abs.(J[:, k] .- Ja[:, k])) ≤ max(rtol * maximum(abs, Ja[:, k]),
                                                     T == Float64 ? 0.0 : 1e-6)
            for k in axes(J, 2))
    end

    extra = [
        ("crescent incl/pa", Dict{String,Any}("c,crin" => 2.0, "c,crout" => 3.0, "c,croff" => 0.5,
            "c,crprojang" => 20.0, "c,incl" => 35.0, "c,pa" => 50.0, "c,x" => 0.4, "c,y" => -0.2),
            ["c,crin", "c,croff", "c,incl", "c,pa", "c,x", "c,y"]),
        ("profile + mode, incl/pa", Dict{String,Any}("d,profile" => "exp(-(\$R/\$s)^2)", "d,s" => 1.0,
            "d,udout" => 6.0, "d,az amp1" => 0.3, "d,az projang1" => 40.0, "d,incl" => 30.0,
            "d,pa" => 70.0), ["d,s", "d,az amp1", "d,az projang1", "d,incl", "d,pa"]),
        ("visfunc incl/pa", Dict{String,Any}("v,visfunc" => "exp(-(\$BL*\$a)^2)", "v,a" => 0.01,
            "v,incl" => 50.0, "v,pa" => 10.0), ["v,a", "v,incl", "v,pa"]),
        ("ldpow incl/pa, shifted", Dict{String,Any}("s,ldpow" => 1.5, "s,alpha" => 0.3,
            "s,incl" => 25.0, "s,pa" => -40.0, "s,f" => "1 - \$g,f", "g,fwhm" => 4.0, "g,f" => 0.2,
            "g,x" => 1.0, "g,y" => 0.5), ["s,ldpow", "s,alpha", "s,incl", "s,pa", "g,f", "g,x"]),
        ("chromatic shape and kernel", Dict{String,Any}("s,ud" => "\$d0 * (\$WL/1.6)^\$p", "d0" => 0.9,
            "p" => 0.3, "s,incl" => "20 + 10 * (\$WL - 1.6)", "s,pa" => 30.0,
            "spatial_kernel" => "\$k * \$WL / 1.6", "k" => 0.2), ["d0", "p", "s,pa", "k"]),
        ("image, flux and shift free", Dict{String,Any}("i,image" => IMG, "i,f" => 0.7, "i,x" => 0.3,
            "i,y" => -0.1, "p,f" => "1 - \$i,f"), ["i,f", "i,x", "i,y"]),
    ]

    for T in (Float64, Float32)
        d = readoifits(AXCIR; T, warn = false, verbose = false)[1, 1]
        uv, wl, mjd = d.uv, d.uv_lam, d.uv_mjd
        @testset "$name, $T" for (name, md, free) in [ZOO; extra]
            fm = dict_to_model(md, free); x = T.([O._resolve_numeric(k, md) for k in free])
            V, J = eval_model_jac(fm, x, uv; wl, mjd)
            Va, Ja = eval_model_grad_ad(fm, x, uv; wl, mjd)
            @test eltype(J) == Complex{T} && eltype(V) == Complex{T}
            @test maximum(abs.(V .- eval_model(fm, x, uv; wl, mjd))) ≤ (T == Float64 ? 1e-13 : 1e-6)
            @test agrees(J, Ja, T)
        end
    end

    d = readoifits(AXCIR; T = Float64, warn = false, verbose = false)[1, 1]
    uv, wl, mjd = d.uv, d.uv_lam, d.uv_mjd

    @testset "finite differences, $name" for (name, md, free) in extra
        fm = dict_to_model(md, free); x = [O._resolve_numeric(k, md) for k in free]
        _, J = eval_model_jac(fm, x, uv; wl, mjd)
        fdm = central_fdm(5, 1)
        for k in eachindex(x), i in (7, 400, 1100)
            g = fdm(t -> (y = copy(x); y[k] = t; collect(reim(eval_model(fm, y, uv; wl, mjd)[i]))), x[k])
            @test complex(g...) ≈ J[i, k] rtol = 1e-6 atol = 1e-9
        end
    end

    @testset "zero baseline" begin
        # The flux point: the discs' branch at ζ = 0, the shift's phase 0.
        md = Dict{String,Any}("s,ud" => 0.9, "s,f" => 0.6, "g,fwhm" => 2.0, "g,f" => "1 - \$s,f",
                              "g,x" => 1.0, "g,y" => 2.0)
        fm = dict_to_model(md, ["s,ud", "s,f", "g,fwhm", "g,x"])
        x = [0.9, 0.6, 2.0, 1.0]
        uv0 = [zeros(2) uv[:, 1:20]]
        V, J = eval_model_jac(fm, x, uv0)
        Va, Ja = eval_model_grad_ad(fm, x, uv0)
        @test V[1] ≈ 1 && agrees(J, Ja, Float64)
        @test J[1, 1] == 0 && J[1, 3] == 0 && J[1, 4] == 0 && abs(J[1, 2]) < 1e-15
    end

    # Measured inside a function: a testset's variables are not concretely typed, and the
    # keyword call and returned tuple would be counted against the caller.
    function jac_allocs(V, J, ws, fm, x, uv, wl, mjd)
        eval_model_jac!(V, J, ws, fm, x, uv; wl, mjd)
        return @allocated eval_model_jac!(V, J, ws, fm, x, uv; wl, mjd)
    end

    @testset "workspace: reused, in place, no allocations" begin
        md, free = ZOO[8][2], ZOO[8][3]                     # binary
        fm = dict_to_model(md, free); x = [md[k] for k in free]
        ws = ModelJacWorkspace(fm, x, uv)
        N = size(uv, 2)
        V = zeros(ComplexF64, N); J = zeros(ComplexF64, N, length(x))
        eval_model_jac!(V, J, ws, fm, 1.1 .* x, uv; wl, mjd)
        eval_model_jac!(V, J, ws, fm, x, uv; wl, mjd)
        V2, J2 = eval_model_jac(fm, x, uv; wl, mjd)
        @test V == V2 && J == J2
        for (name, md, free) in ZOO[[1, 3, 7, 8, 9, 11, 21]]
            fm = dict_to_model(md, free); x = [md[k] for k in free]
            ws = ModelJacWorkspace(fm, x, uv)
            J = zeros(ComplexF64, N, length(x))
            @test jac_allocs(V, J, ws, fm, x, uv, wl, mjd) == 0
        end
        @test_throws DimensionMismatch eval_model_jac!(V, zeros(ComplexF64, N, 9), ws, fm, x, uv)
        @test_throws ArgumentError eval_model_jac!(V, J, ws, fm, x, uv; n = 1)
    end
end

# The consumers of the Jacobian (J4): `eval_model_grad` routed to it, the cache the χ² gradient
# and the residual Jacobian keep between calls, and the allocation-free χ² kernels behind them,
# each against the path it replaces.
@testset "model Jacobian consumers" begin
    O = OITOOLS
    include(joinpath(@__DIR__, "bench", "zoo.jl"))
    V1295 = joinpath(@__DIR__, "..", "demos", "data", "V1295_Aql.oifits")
    binary = Dict{String,Any}("s,ud" => 0.83, "s,f" => "1 - \$c,f", "c,f" => 0.05, "c,x" => 6.0,
                              "c,y" => -28.0)
    free = ["s,ud", "c,f", "c,x", "c,y"]
    weightsets = ([1.0, 0, 1], [1.0, 1, 1, 1, 1], [2.0, 0.5, 3, 0.7, 1.3])

    @testset "eval_model_grad is the closed-form Jacobian" begin
        d = readoifits(AXCIR; T = Float64, warn = false, verbose = false)[1, 1]
        fm = dict_to_model(binary, free); x = [0.83, 0.05, 6.0, -28.0]
        @test eval_model_grad(fm, x, d.uv; wl = d.uv_lam) == eval_model_jac(fm, x, d.uv; wl = d.uv_lam)
        # Dual parameters (a Hessian through it) and n ≠ 0 take the AD path.
        H = O.ForwardDiff.jacobian(x_ -> real.(eval_model_grad(fm, x_, d.uv[:, 1:30])[2][:, 1]), x)
        Hfd = FiniteDifferences.jacobian(central_fdm(5, 1),
                  x_ -> real.(eval_model_grad(fm, x_, d.uv[:, 1:30])[2][:, 1]), x)[1]
        @test isapprox(H, Hfd; rtol = 1e-6, atol = 1e-9)
    end

    @testset "ModelJacCache" begin
        d = readoifits(AXCIR; T = Float64, warn = false, verbose = false)[1, 1]
        fm = dict_to_model(binary, free); x = [0.83, 0.05, 6.0, -28.0]
        c = ModelJacCache()
        V, J = eval_model_jac!(c, fm, x, d.uv; wl = d.uv_lam)
        @test (copy(V), copy(J)) == eval_model_jac(fm, x, d.uv; wl = d.uv_lam)
        ws = c.ws
        eval_model_jac!(c, fm, 1.1 .* x, d.uv; wl = d.uv_lam)
        @test c.ws === ws                                       # reused
        V3, J3 = eval_model_jac!(c, fm, x, d.uv[:, 1:40]; wl = d.uv_lam[1:40])
        @test c.ws !== ws && size(J3) == (40, 4)                # rebuilt for other uv
        fm2 = dict_to_model(ZOO[1][2], ZOO[1][3])
        V4, J4 = eval_model_jac!(c, fm2, [0.9], d.uv)
        @test (V4, J4) == eval_model_jac(fm2, [0.9], d.uv)      # rebuilt for another model
        for w in weightsets
            @test model_to_chi2_fg(fm, x, d; weights = w, jac_cache = c) ==
                  model_to_chi2_fg(fm, x, d; weights = w)
            @test O.residuals_flat_jac(fm, x, d; weights = w, jac_cache = c) ==
                  O.residuals_flat_jac(fm, x, d; weights = w)
        end
    end

    @testset "χ² kernel against cvis_to_chi2_fg, $T" for T in (Float64, Float32)
        d = readoifits(V1295; T, warn = false, verbose = false)[1, 1]
        @test d.nvisamp > 0 && d.nt3amp > 0
        fm = dict_to_model(binary, free)
        V = eval_model(fm, T[0.83, 0.05, 6.0, -28.0], d.uv)
        g = similar(V)
        for w in weightsets
            c1, g1 = O.cvis_to_chi2_fg(V, d; weights = O._pad_weights(w))
            c2 = O._cvis_to_chi2_fg!(g, V, d, ntuple(i -> T(O._pad_weights(w)[i]), 7))
            @test c2 isa T
            @test isapprox(c2, c1; rtol = T == Float64 ? 1e-13 : 1e-6)
            @test maximum(abs.(g .- g1)) ≤ (T == Float64 ? 1e-13 : 1e-6) * maximum(abs, g1)
        end
    end

    # Through the module's own name: `O` is a captured local, and a call through it dispatches
    # dynamically, which is the caller's allocation, not the kernel's.
    function kernel_allocs(g, V, d, w)
        OITOOLS._cvis_to_chi2_fg!(g, V, d, w)
        return @allocated OITOOLS._cvis_to_chi2_fg!(g, V, d, w)
    end
    @testset "χ² kernel allocates nothing" begin
        d = readoifits(V1295; T = Float64, warn = false, verbose = false)[1, 1]
        V = eval_model(dict_to_model(binary, free), [0.83, 0.05, 6.0, -28.0], d.uv)
        @test kernel_allocs(similar(V), V, d, (1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0)) == 0
    end

    @testset "smeared χ² gradient, $(basename(f)), $T, $mode" for f in (AXCIR, V1295),
            T in (Float64, Float32), mode in ("incoherent", "coherent")
        d = readoifits(f; T, warn = false, verbose = false)[1, 1]
        fm = dict_to_model(merge(binary, Dict{String,Any}("smear" => 3, "smear_mode" => mode)), free)
        x = T[0.83, 0.05, 6.0, -28.0]
        c = ModelJacCache()
        tol = T == Float64 ? 1e-12 : 1e-5
        for w in weightsets
            # Against χ² = Σr², ∇χ² = 2Jᵀr from the smeared residual Jacobian.
            r, J = O._smeared_residuals_jac(fm, x, d, O._pad_weights(w))
            c2, g2 = model_to_chi2_fg(fm, x, d; weights = w, jac_cache = c)
            @test c2 isa T
            @test isapprox(c2, sum(abs2, r); rtol = tol)
            @test maximum(abs.(g2 .- 2 .* (J' * r))) ≤ tol * maximum(abs, 2 .* (J' * r))
            @test (c2, g2) == model_to_chi2_fg(fm, x, d; weights = w, jac_cache = c) ==
                             model_to_chi2_fg(fm, x, d; weights = w)
        end
    end
end
