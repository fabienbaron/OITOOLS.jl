# The component kernels the model Jacobian is built from (PLAN_ModelJacobian.md, J2), judged
# against references that do not share their code: B_ν against SpecialFunctions' `besselj`,
# ldpow against its series summed in 512-bit arithmetic,
# each law's value against the original `vis_*` / `_vis_*` implementations, every derivative
# against ForwardDiff through those originals and against finite differences.
using ForwardDiff, FiniteDifferences, SpecialFunctions

@testset "model kernels" begin
    O = OITOOLS
    MAS = O.MAS2RAD

    @testset "B_ν against besselj" begin
        xs = [1e-6; 1e-3; 0.05; collect(range(0.1, 60; length = 400))]
        for ν in (1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0)
            got = [O._bnu(Val(ν), x) for x in xs]
            ref = [besselj(ν, x) / x^ν for x in xs]
            scale = maximum(abs, ref)               # 1/(2^ν Γ(ν+1)): B_ν's largest value
            @test maximum(abs.(got .- ref)) / scale < 1e-12
            @test O._bnu(Val(ν), 0.0) ≈ 1 / (2^ν * gamma(ν + 1))
        end
    end

    ρ = collect(range(0, 2.5e8; length = 300))       # cycles/rad, through several lobes

    # The laws as originally written, with `besselj` at every order: an independent reference,
    # since `vis_ld*` themselves now run on the kernels under test.
    B(ν, ζ) = ζ < 1e-8 ? 1 / (2^ν * gamma(ν + 1)) : besselj(ν, ζ) / ζ^ν
    ζs(θ, ρ) = @. π * θ * ρ / MAS
    K32, K54, K74 = sqrt(π / 2), gamma(2.25) * 2^1.25 / 2.5, gamma(2.75) * 2^1.75 / 3.5
    ref_ldlin(p, ρ) = [ζ < 1e-8 ? one(ζ) :
        ((1 - p[2]) * B(1, ζ) + p[2] * K32 * B(1.5, ζ)) / (0.5 - p[2] / 6) for ζ in ζs(p[1], ρ)]
    ref_ldquad(p, ρ) = [ζ < 1e-8 ? one(ζ) :
        ((1 - p[2] - p[3]) * B(1, ζ) + (p[2] + 2p[3]) * K32 * B(1.5, ζ) - 2p[3] * B(2, ζ)) /
        (0.5 - p[2] / 6 - p[3] / 12) for ζ in ζs(p[1], ρ)]
    ref_ldsqrt(p, ρ) = [ζ < 1e-8 ? one(ζ) :
        ((1 - p[2] - p[3]) * B(1, ζ) + p[2] * K32 * B(1.5, ζ) + p[3] * K54 * B(1.25, ζ)) /
        (0.5 - p[2] / 6 - p[3] / 10) for ζ in ζs(p[1], ρ)]
    ref_ldclaret4(p, ρ) = [ζ < 1e-8 ? one(ζ) :
        ((1 - p[2] - p[3] - p[4] - p[5]) * B(1, ζ) + p[2] * K54 * B(1.25, ζ) +
         p[3] * K32 * B(1.5, ζ) + p[4] * K74 * B(1.75, ζ) + p[5] * 2 * B(2, ζ)) /
        (0.5 - p[2] / 10 - p[3] / 6 - 3p[4] / 14 - p[5] / 4) for ζ in ζs(p[1], ρ)]
    # (kind, parameter values, original kernel)
    cases = [
        (:ud, [3.0], O._vis_ud),
        (:gaussian, [2.0], O._vis_gaussian),
        (:ldlin, [3.0, 0.4], ref_ldlin),
        (:ldquad, [3.0, 0.3, 0.2], ref_ldquad),
        (:ldsqrt, [3.0, 0.3, 0.2], ref_ldsqrt),
        (:ldclaret4, [3.0, 0.5, -0.2, 0.3, -0.1], ref_ldclaret4),
        (:ring, [2.0, 3.0], O._vis_ring),
        (:gaussian_ring, [2.0, 3.5], O._vis_gaussian_ring),
        (:point, Float64[], O._vis_point),
        (:resolved, Float64[], O._vis_resolved),
    ]

    @testset "$kind: value and derivatives" for (kind, p, orig) in cases
        n = length(p); N = length(ρ)
        W = zeros(N); D = zeros(N, max(n, 1)); Dρ = zeros(N)
        O._vis_jac!(W, D, Dρ, Val(kind), p, ρ)
        @test maximum(abs.(W .- orig(p, ρ))) < 1e-12
        # ∂/∂params: ForwardDiff through the ORIGINAL kernel, parameter by parameter.
        for j in 1:n
            ref = ForwardDiff.derivative(t -> orig([k == j ? t : p[k] for k in 1:n], ρ), p[j])
            @test maximum(abs.(D[:, j] .- ref)) <= 1e-10 * max(1.0, maximum(abs, ref))
        end
        # ∂/∂ρ, point by point.
        refρ = [ForwardDiff.derivative(r -> orig(p, [r])[1], ρi) for ρi in ρ]
        @test maximum(abs.(Dρ .- refρ)) <= 1e-10 * max(1.0, maximum(abs, refρ))
        # Finite differences at a few points, so the reference is not trusted blindly.
        fdm = central_fdm(5, 1)
        for i in (40, 150, 260), j in 1:n
            fd = fdm(t -> orig([k == j ? t : p[k] for k in 1:n], [ρ[i]])[1], p[j])
            @test D[i, j] ≈ fd rtol = 1e-7 atol = 1e-10
        end
        # No allocation with the buffers given.
        @test (@allocated O._vis_jac!(W, D, Dρ, Val(kind), p, ρ)) == 0
    end

    @testset "ldpow: the order moves with α" begin
        # Reference: V = Σ_k q^k / (k! (ν+1)…(ν+k)), q = −ζ²/4, and its ζ and ν derivatives, summed
        # in 512-bit arithmetic, where the alternating terms cancel harmlessly at every ζ here.
        function ref_pow(ζ, α)
            setprecision(BigFloat, 512) do
                ζb = big(ζ); ν = big(α) / 2 + 1; q = -ζb^2 / 4
                t = big(1); s = big(0); V = big(1); dν = big(0); dq = big(0); qk = big(1)
                for k in 1:400
                    t /= k * (ν + k); s += 1 / (ν + k); dq += k * t * qk
                    qk *= q; V += t * qk; dν -= t * s * qk
                end
                Float64.((V, -dq * ζb / 2, dν / 2))
            end
        end
        # ζ across the series (< 2, < 12), Hankel's expansion and the boundaries between them.
        zs = [0.0; 1e-9; 0.3; 1.999; 2.0; 2.001; collect(2.5:0.75:11.5); 11.999; 12.0; 12.001;
              collect(12.5:1.5:60)]
        # α = 1 ends Hankel's series for J (ν = 3/2) but not for ∂J/∂ν; α = 16 (ν + 1 > 7) takes
        # the library path.
        for α in (0.0, 0.15, 0.6, 1.0, 2.0, 4.0, 10.0, 16.0)
            tol = α > 10 ? 1e-7 : 1e-12
            G, ψ = O._ldpow_consts(α / 2 + 1)
            # below `_ZETA_MIN` the kernel returns V = 1 and zero derivatives, as every law does
            err = maximum(maximum(abs.(collect(O._ldpow_d(ζ, α, G, ψ)) .- collect(ref_pow(ζ, α))))
                          for ζ in zs if ζ ≥ O._ZETA_MIN)
            @test err < tol
            @test maximum(abs(O._ldpow_v(ζ, α) - ref_pow(ζ, α)[1]) for ζ in zs) < tol
            # the value against the library, which shares no code with either branch
            ν = α / 2 + 1
            @test maximum(abs(O._ldpow_v(ζ, α) - gamma(ν + 1) * 2^ν * besselj(ν, ζ) / ζ^ν)
                          for ζ in zs if ζ > 1e-8) < 1e-11
            # ForwardDiff through the value kernel against the hand-written derivatives
            for ζ in (0.7, 5.0, 11.0, 13.0, 40.0)
                d = O._ldpow_d(ζ, α)
                @test ForwardDiff.derivative(z -> O._ldpow_v(z, α), ζ) ≈ d[2] rtol = 1e-8 atol = 1e-12
                α > 10 || @test ForwardDiff.derivative(a -> O._ldpow_v(ζ, a), α) ≈ d[3] rtol = 1e-8 atol = 1e-12
            end
        end

        # The Jacobian kernel and everything routed through it, out to ζ ≈ 25.
        ρw = collect(range(0, 1.7e9; length = 400)); N = length(ρw)
        for α in (0.15, 0.6)
            p = [3.0, α]
            W = zeros(N); D = zeros(N, 2); Dρ = zeros(N)
            O._vis_jac!(W, D, Dρ, Val(:ldpow), p, ρw)
            r = [ref_pow(π * 3.0 * x / MAS, α) for x in ρw]
            @test maximum(abs.(W .- first.(r))) < 2e-10
            @test maximum(abs.(W .- O.vis_ldpow(3.0, α, ρw))) < 1e-12
            @test maximum(abs.(D[:, 1] .- getindex.(r, 2) .* π .* ρw ./ MAS)) < 2e-10 * π * 1.7e9 / MAS
            @test maximum(abs.(D[:, 2] .- last.(r))) < 2e-10
            @test maximum(abs.(Dρ .- getindex.(r, 2) .* π .* 3.0 ./ MAS)) < 2e-10 * π * 3.0 / MAS
            fdm = central_fdm(5, 1)
            for i in (60, 200, 380)
                @test D[i, 2] ≈ fdm(a -> O.vis_ldpow(3.0, a, [ρw[i]])[1], α) rtol = 1e-6 atol = 1e-10
            end
            # the ForwardDiff specialisation and the rrule carry the same numbers
            refα = ForwardDiff.derivative(a -> O._vis_ldpow([3.0, a], ρw), α)
            @test refα ≈ D[:, 2] rtol = 1e-12
            @test (@allocated O._vis_jac!(W, D, Dρ, Val(:ldpow), p, ρw)) == 0
        end
        # a chromatic α, one value per point
        αv = collect(range(0.1, 0.9; length = N)); W = zeros(N); D = zeros(N, 2); Dρ = zeros(N)
        O._vis_jac!(W, D, Dρ, Val(:ldpow), [fill(3.0, N), αv], ρw)
        @test maximum(abs(D[i, 2] - ref_pow(π * 3.0 * ρw[i] / MAS, αv[i])[3]) for i in 1:N) < 2e-10
    end

    @testset "crescent: in (u, v), complex" begin
        t = range(0, 2π; length = 300)
        u = collect(2.0e8 .* sin.(3t) .* cos.(t)); v = collect(2.0e8 .* sin.(3t) .* sin.(t))
        p = [2.0, 3.0, 0.6, 35.0]; N = length(u)
        W = zeros(ComplexF64, N); D = zeros(ComplexF64, N, 4)
        Du = zeros(ComplexF64, N); Dv = zeros(ComplexF64, N)
        O._crescent_jac!(W, D, Du, Dv, p, u, v)
        @test maximum(abs.(W .- O._vis_crescent(p, u, v))) < 1e-12
        cd(f, x) = ForwardDiff.derivative(t -> real(f(t)), x) + im * ForwardDiff.derivative(t -> imag(f(t)), x)
        for j in 1:4
            ref = [cd(s -> O._vis_crescent([k == j ? s : p[k] for k in 1:4], [u[i]], [v[i]])[1], p[j]) for i in 1:N]
            @test maximum(abs.(D[:, j] .- ref)) < 1e-10 * max(1, maximum(abs, ref))
        end
        refu = [cd(s -> O._vis_crescent(p, [s], [v[i]])[1], u[i]) for i in 1:N]
        refv = [cd(s -> O._vis_crescent(p, [u[i]], [s])[1], v[i]) for i in 1:N]
        @test maximum(abs.(Du .- refu)) < 1e-10 * maximum(abs, refu)
        @test maximum(abs.(Dv .- refv)) < 1e-10 * maximum(abs, refv)
        @test (@allocated O._crescent_jac!(W, D, Du, Dv, p, u, v)) == 0
    end

    @testset "visfunc: dual numbers through the user's expression" begin
        fn = O.compile_visfunc("exp(-(\$BL * \$a)^2) * cis(\$c * \$BL) * (\$WL / 1.6)", ["a", "c"])
        Bm = collect(range(0, 120; length = 200)); N = length(Bm)
        wl = collect(range(1.5e-6, 1.8e-6; length = N))
        a = 0.012; c = collect(range(0.01, 0.02; length = N))      # c per point (chromatic)
        W = zeros(ComplexF64, N); D = zeros(ComplexF64, N, 2); DB = zeros(ComplexF64, N)
        O._visfunc_jac!(W, D, DB, fn, (a, c), Bm, wl, NaN)
        @test W ≈ fn(Bm, wl, NaN, a, c)
        h = 1e-7
        @test D[:, 1] ≈ (fn(Bm, wl, NaN, a + h, c) .- fn(Bm, wl, NaN, a - h, c)) ./ 2h rtol = 1e-6
        @test D[:, 2] ≈ (fn(Bm, wl, NaN, a, c .+ h) .- fn(Bm, wl, NaN, a, c .- h)) ./ 2h rtol = 1e-6
        @test DB ≈ (fn(Bm .+ h, wl, NaN, a, c) .- fn(Bm .- h, wl, NaN, a, c)) ./ 2h rtol = 1e-5 atol = 1e-8
    end

    @testset "radial profile with azimuthal modes" begin
        d = readoifits(joinpath(@__DIR__, "..", "demos", "data", "AX_Cir.oifits");
                       T = Float64, warn = false, verbose = false)[1, 1]
        uv = d.uv; u = uv[1, :]; v = uv[2, :]; N = size(uv, 2)
        md = Dict{String,Any}("d,profile" => "exp(-(\$R/\$s)^\$q)", "d,s" => 1.2, "d,q" => 2.0,
                              "d,udout" => 6.0, "d,az amp1" => 0.3, "d,az projang1" => 40.0,
                              "d,az amp2" => 0.1, "d,az projang2" => -20.0)
        free = ["d,s", "d,q", "d,az amp1", "d,az projang1", "d,az amp2", "d,az projang2"]
        m = dict_to_model(md, free); x = [md[k] for k in free]
        comp = only(m.components)
        ρh  = hypot.(u, v); Bh = ρh ./ O._MAS2RAD_PM; ψ = atan.(v, u)
        # the profile's parameters follow its diameters, [d_min, d_max, s, q]
        st = O._hankel_setup(comp, Bh, [0.0, 6.0, 1.2, 2.0]); s, r, μ, g = st.s, st.r, st.μ, st.g
        prof(p) = comp.profile_fn(s, μ, 0.0, 6.0, p...)
        I  = prof([1.2, 2.0]); dI = ForwardDiff.jacobian(prof, [1.2, 2.0])
        W = zeros(ComplexF64, N); D = zeros(ComplexF64, N, 6)
        DB = zeros(ComplexF64, N); Dψ = zeros(ComplexF64, N)
        O._hankel_jac!(W, D, DB, Dψ, r, g, I, dI, Bh, ψ, [1, 2], (0.3, 0.1), (40.0, -20.0), comp.kcache)
        @test W ≈ eval_model(m, x, uv) rtol = 1e-12
        J = eval_model_grad(m, x, uv)[2]                 # ForwardDiff through eval_model
        @test maximum(abs.(D .- J)) < 1e-10 * maximum(abs, J)
        # ∂/∂u, ∂/∂v through Bh and ψ, against finite differences in the uv coordinates.
        dBdu = u ./ (ρh .* O._MAS2RAD_PM); dψdu = -v ./ ρh .^ 2
        h = 1.0                                           # cycles/rad; |uv| ~ 1e7
        fdu = (eval_model(m, x, [u' .+ h; v']) .- eval_model(m, x, [u' .- h; v'])) ./ 2h
        @test DB .* dBdu .+ Dψ .* dψdu ≈ fdu rtol = 1e-5
        # The transform matrices are built once per uv set.
        nb = comp.kcache.nbuilt
        O._hankel_jac!(W, D, DB, Dψ, r, g, I, dI, Bh, ψ, [1, 2], (0.3, 0.1), (40.0, -20.0), comp.kcache)
        eval_model(m, x, uv)
        @test comp.kcache.nbuilt == nb
    end

    @testset "per-point (chromatic) parameters" begin
        N = length(ρ)
        θv = collect(range(2.5, 3.5; length = N)); uv_ = collect(range(0.2, 0.5; length = N))
        W = zeros(N); D = zeros(N, 2); Dρ = zeros(N)
        O._vis_jac!(W, D, Dρ, Val(:ldlin), (θv, uv_), ρ)
        @test W ≈ [ref_ldlin([θv[i], uv_[i]], [ρ[i]])[1] for i in 1:N] rtol = 1e-12
        h = 1e-6
        @test D[:, 1] ≈ (O.vis_ldlin(θv .+ h, uv_, ρ) .- O.vis_ldlin(θv .- h, uv_, ρ)) ./ 2h rtol = 1e-6 atol = 1e-9
    end

    @testset "Float32 in, Float32 out" begin
        ρ32 = Float32.(ρ); N = length(ρ)
        W = zeros(Float32, N); D = zeros(Float32, N, 3); Dρ = zeros(Float32, N)
        O._vis_jac!(W, D, Dρ, Val(:ldquad), Float32[3.0, 0.3, 0.2], ρ32)
        W64 = zeros(N); D64 = zeros(N, 3); Dρ64 = zeros(N)
        O._vis_jac!(W64, D64, Dρ64, Val(:ldquad), [3.0, 0.3, 0.2], ρ)
        @test maximum(abs.(W .- W64)) < 1e-5
        @test maximum(abs.(D .- D64)) < 1e-4 * maximum(abs, D64)
        # ldpow on every branch, the library one (α > 10, `_dbesselj_dnu`'s stencil) included:
        # a Float64 on any branch makes the broadcast box every point
        @test Core.Compiler.return_type(O._ldpow_d, Tuple{Float32,Float32,Float32,Float32}) ==
              NTuple{3,Float32}
        for α in (0.14f0, 14.0f0)
            G, ψ = O._ldpow_consts(α / 2 + 1)
            ζ = collect(range(0.1f0, 110.0f0, length = 200))
            @test eltype(O._ldpow_d.(ζ, α, G, ψ)) == NTuple{3,Float32}
        end
    end
end
