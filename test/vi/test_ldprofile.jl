# The free-form limb-darkening profile by variational inference (`fit_ld_profile_vi`): the
# problem's metric products and gradient against finite differences, for both priors; the
# clamps that keep a wild Newton step cheap; and a simulated star recovered within its band.

using OITOOLS, VarInf, LinearAlgebra, Statistics, Random
const LDVI = Base.get_extension(OITOOLS, :OITOOLSVarInfExt)

const SRC = joinpath(pkgdir(OITOOLS), "demos", "data", "AlphaCenA.oifits")
truth = Dict{String,Any}("star,ldpow" => 8.4, "star,alpha" => 0.2, "star,f" => 1.0)
file = tempname() * ".oifits"
simulate_from_oifits(SRC, file; flat_model = parse_model(truth, ["star,ldpow", "star,alpha"]),
                     flat_params = [8.4, 0.2], seed = 7)
data = readoifits(file; warn = false, verbose = false)[1, 1]
start = Dict{String,Any}("star,ldlin" => 8.4, "star,u" => 0.4, "star,f" => 1.0)
r0 = fit_model_lsqfit(start, ["star,ldlin", "star,u"], data; weights = [1.0, 1, 0, 0, 0, 0, 0])
start["star,ldlin"], start["star,u"] = r0.x_opt

@testset "the problem's derivatives ($prior prior)" for prior in (:fixed, :learned)
    p, md, free = LDVI._ld_vi_problem(start, "star", data; ngrid = 12, prior)
    @test free[1] == "star,ldpiecewise" && length(free) == 12
    nz = VarInf.latent_size(p)
    z = 0.3 .* randn(MersenneTwister(1), nz)
    T(z) = VarInf.transformation(p, z)
    h = 1e-6
    J = reduce(hcat, [(T(z .+ h .* (1:nz .== j)) .- T(z .- h .* (1:nz .== j))) ./ 2h for j in 1:nz])
    v = randn(MersenneTwister(2), nz); w = randn(MersenneTwister(3), length(T(z)))
    @test VarInf.right_sqrt_metric(p, z, v) ≈ J * v rtol = 1e-5
    @test VarInf.left_sqrt_metric(p, z, w) ≈ J' * w rtol = 1e-5
    E, g = VarInf.energy_and_gradient(p, z)
    @test E ≈ 0.5 * (sum(abs2, T(z)) + sum(abs2, z))
    gfd = [(VarInf.energy_and_gradient(p, z .+ h .* (1:nz .== j))[1] -
            VarInf.energy_and_gradient(p, z .- h .* (1:nz .== j))[1]) / 2h for j in 1:nz]
    @test g ≈ gfd rtol = 1e-5
    # the transformation is the whitened residual vector of the model it maps to
    x, _ = LDVI._ld_params(p, z)
    @test T(z) ≈ OITOOLS.residuals_flat(p.fm, x, data; weights = p.weights)
end

@testset "a wild step stays finite and cheap" begin
    p, _, _ = LDVI._ld_vi_problem(start, "star", data; ngrid = 12)
    z = zeros(VarInf.latent_size(p)); z[1] = 60.0; z[2:end] .= 200.0
    t = @elapsed E, g = VarInf.energy_and_gradient(p, z)
    @test isfinite(E) && all(isfinite, g) && t < 1.0
    x, _ = LDVI._ld_params(p, z)
    @test x[1] ≈ 8.4 * exp(0.8) rtol = 0.02
end

@testset "a μ^0.2 star at 8.4 mas" begin
    r = fit_ld_profile_vi(start, "star", data; ngrid = 24, seed = 1)
    @test r isa LDProfileVIResult && r.method === :geovi && r.prior === :fixed
    @test size(r.I_samples) == (24, 100) && r.I_mean[end] == 1 && r.I_std[end] == 0
    @test r.chi2r < 1.3
    @test abs(r.diameter - 8.4) < 2.5 * r.diameter_std
    # the true profile within 2σ at the nodes the data constrain (the limb's μ^0.2 has an
    # infinite slope no grid follows, so μ < 0.1 is left out)
    sel = r.μ .>= 0.1 .&& r.μ .< 1
    z = (r.I_mean[sel] .- r.μ[sel] .^ 0.2) ./ r.I_std[sel]
    @test mean(abs.(z) .< 2) >= 0.9
    @test all(r.I_lo .<= r.I_hi) && all(r.F_lo .<= r.F_hi)
    # I/⟨I⟩ integrates to 1 over the disc
    g0 = OITOOLS._ld_table(:piecewise, 24, 0.0).g0
    @test 2 * dot(g0, r.F_mean) ≈ 1 rtol = 1e-2
    # the model at the variational mean is the grid disc
    @test r.free[1] == "star,ldpiecewise" && haskey(r.model, "star,ld23")
    @test_throws ArgumentError fit_ld_profile_vi(start, "star", data; prior = :other)
    @test_throws ArgumentError fit_ld_profile_vi(start, "star", data; method = :other)
end

@testset "error bars scaled when too small (α Cen A)" begin
    real = readoifits(SRC; warn = false, verbose = false)[1, 1]
    s = Dict{String,Any}("star,ldlin" => 8.4, "star,u" => 0.3, "star,f" => 1.0)
    rl = fit_model_lsqfit(s, ["star,ldlin", "star,u"], real; weights = [1.0, 1, 0, 0, 0, 0, 0])
    s["star,ldlin"], s["star,u"] = rl.x_opt
    r = fit_ld_profile_vi(s, "star", real; scale_errors = true, seed = 1, n_draws = 20)
    @test r.error_scale > 1.5 && 0.8 < r.chi2r < 1.3 && r.chi2r_map < 1.1
    # simulated data, whose error bars are right, are left as they are
    q = fit_ld_profile_vi(start, "star", data; scale_errors = true, seed = 1, n_draws = 20)
    @test q.error_scale < 1.1
end

@testset "bagged VI" begin
    r = fit_ld_profile_vi(start, "star", data; ngrid = 12, seed = 1, n_draws = 10, bootstrap = 4)
    @test r.bootstrap == 4 && r.single isa LDProfileVIResult && r.single.bootstrap == 0
    # four replicates of four antithetic pairs each
    @test size(r.I_samples) == (12, 32) && length(r.diameter_samples) == 32
    @test r.model == r.single.model && r.chi2r == r.single.chi2r
    @test abs(r.diameter - 8.4) < 3 * r.diameter_std
    @test all(r.I_lo .<= r.I_hi) && r.I_mean[end] == 1
    # the same seed, the same pool
    q = fit_ld_profile_vi(start, "star", data; ngrid = 12, seed = 1, n_draws = 10, bootstrap = 4)
    @test q.I_samples == r.I_samples
    @test fit_ld_profile_vi(start, "star", data; ngrid = 12, seed = 1, n_draws = 10).bootstrap == 0
    @test_throws ArgumentError fit_ld_profile_vi(start, "star", data; bootstrap = 1)
end
