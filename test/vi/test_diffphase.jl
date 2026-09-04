using OITOOLS, OIVI, LinearAlgebra, FiniteDifferences, Test
using OITOOLS, VarInf
const OIVI = Base.get_extension(OITOOLS, :OITOOLSVarInfExt)
using .OIVI
using Statistics: mean

oifits = "/home/baron/SOFTWARE/OITOOLS.jl/demos/data/BC2026/OBJECT1_N.oifits"
data = readoifits(oifits; polychromatic=true, filter_bad_data=true, verbose=false, warn=false)

nf_test = 3
npix = 16; pixsize = 0.3
data_sub = data[1:nf_test, :]
ft = setup_ft(data_sub, npix, pixsize)
dp = build_diffphase_config(ft, data_sub, nf_test)

# Build a realistic image (sky model)
freq = [3e8 / mean(data_sub[c,1].uv_lam) for c in 1:nf_test]
p = SkyModelParams(npix, pixsize, freq;
    R_mas=3.0, u=0.2,
    spatial_slope_prior=(-4.0, 1.0),
    spatial_fluct_prior=(1.6487, 2.1612))
xi = 0.1 * randn(_latent_size(p))
image = sky_forward(xi, p)
println("Image size: ", size(image), " range: ", extrema(image))
println("Image flux per channel: ", [sum(image[:,:,c]) for c in 1:nf_test])

# 1. Test diffphi_forward runs
dphi = diffphi_forward(image, dp)
println("\ndiffphi_forward output: size=", size(dphi), " range=", extrema(dphi))

@testset "diffphi_adjoint (finite difference)" begin
    g_dphi = randn(size(dphi))
    g_image = diffphi_adjoint(g_dphi, image, dp)

    v_image = randn(size(image)) .* 1e-3
    fdm = central_fdm(5, 1)
    Jv = jvp(fdm, img -> diffphi_forward(img, dp), (image, v_image))

    lhs = dot(g_image, v_image)
    rhs = dot(g_dphi, Jv)
    rel_err = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-20)
    println("  diffphi_adjoint FD: rel_err = ", rel_err)
    @test rel_err < 1e-4
end

@testset "diffphi_jvp (finite difference)" begin
    d_image = randn(size(image)) .* 1e-3
    fdm = central_fdm(5, 1)
    d_dphi_jvp = diffphi_jvp(image, d_image, dp)
    d_dphi_fd = jvp(fdm, img -> diffphi_forward(img, dp), (image, d_image))

    rel_err = norm(d_dphi_jvp .- d_dphi_fd) / (norm(d_dphi_fd) + 1e-20)
    println("  diffphi_jvp FD: rel_err = ", rel_err)
    @test rel_err < 1e-4
end

@testset "diffphi adjoint-JVP consistency" begin
    v_img = randn(size(image)) .* 1e-3
    v_dphi = randn(size(dphi))
    g_img = diffphi_adjoint(v_dphi, image, dp)
    d_dphi = diffphi_jvp(image, v_img, dp)
    lhs = dot(g_img, v_img)
    rhs = dot(v_dphi, d_dphi)
    rel_err = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-20)
    println("  diffphi adjoint-JVP consistency: rel_err = ", rel_err)
    @test rel_err < 1e-8
end

@testset "energy_fg with diffphi (finite difference)" begin
    using .OIVI: _split_latent

    # Build ObsContext with diffphi
    obs_vec = ObservationConfig[]
    for c in 1:nf_test
        push!(obs_vec, ObservationConfig(ft[c,1], data_sub[c,1]))
    end
    sigma_vec = vcat([vcat(o.v2_err, o.t3amp_err, o.t3phi_err) for o in obs_vec]...)
    sigma_vec = vcat(sigma_vec, vec(dp.diffphi_err))
    ctx = ObsContext(obs_vec, sigma_vec, dp)

    e, g = energy_fg(xi, p, ft, data_sub, ctx)
    println("  energy with diffphi: ", e)

    fdm = central_fdm(5, 1)
    g_fd = grad(fdm, z0 -> energy_fg(z0, p, ft, data_sub, ctx)[1], xi)[1]

    rel_err = norm(g .- g_fd) / (norm(g_fd) + 1e-20)
    println("  energy_fg (with diffphi) FD: rel_err = ", rel_err)
    @test rel_err < 1e-4
end

println("\nAll diffphase FD tests passed!")
