# test_vis_functions.jl — the legacy (param, uv) visibility function API.
#
# This was demos/example_model_fitting_functions.jl, a gallery that computed each visibility
# function and inverted it to an image for visual inspection. It was a test wearing an
# example's clothes: it asserted nothing, so when it broke it broke silently. By the time it
# moved here it had three independent faults — `using PyPlot` (dead since the PythonCall
# migration), an NFFT grid whose nodes fell outside [-1/2, 1/2) because it paired a 3000 m
# baseline with a 3.2 mas field of view, and `color=:red` Symbols that PyCall used to coerce
# for matplotlib and PythonCall does not.
#
# These functions live in vis_functions.jl and are exported, but they are a *parallel* API to
# the flat model dict: only ud, ldlin, ldquad and ldpow are reachable from `dict_to_model`.
# Nothing else covered them.

using OITOOLS, Test, FFTW, NFFT

@testset "legacy visibility functions" begin

    λ    = 1.6e-6
    Bmax = 600.0
    # Odd count so the grid contains exactly (0,0): that is the only sample where V must be
    # 1, and it also exercises each function's own guard against the rho -> 0 singularity.
    N    = 65
    x    = collect(range(-Bmax, Bmax, length = N)) ./ λ
    uv   = Array{Float64}(undef, 2, N * N)
    uv[1, :] = vec(repeat(x, 1, N))
    uv[2, :] = vec(repeat(x, 1, N)')
    ρ    = sqrt.(uv[1, :].^2 .+ uv[2, :].^2)
    zero_baseline = argmin(ρ)
    @test ρ[zero_baseline] == 0.0      # the grid really does contain the origin

    # (name, function, params). Every one of these is exported.
    cases = [
        ("ud",                 visibility_ud,                 [10.0]),
        ("ldlin",              visibility_ldlin,              [10.0, 0.7]),
        ("ldquad",             visibility_ldquad,             [10.0, 0.5, 0.5]),
        ("ldpow",              visibility_ldpow,              [10.0, 0.3]),
        ("ldsquareroot",       visibility_ldsquareroot,       [10.0, 0.5, 0.3]),
        ("ellipse_uniform",    visibility_ellipse_uniform,    [10.0, 2.0, -45.0]),
        ("ellipse_quad",       visibility_ellipse_quad,       [10.0, 0.5, 0.5, 2.0, -45.0]),
        ("thin_ring",          visibility_thin_ring,          [10.0, 0.0, 0.0]),
        ("Gaussian_ring",      visibility_Gaussian_ring,      [10.0, 90.0, 60.0, 0.1]),
        ("Lorentzian_ring",    visibility_Lorentzian_ring,    [10.0, 90.0, 60.0, 0.1]),
        ("Gaussian_ring_az",   visibility_Gaussian_ring_az,   [10.0, -90.0, 60.0, 0.15, 0.5, 0.0, 0.0, 0.0]),
    ]

    for (name, f, p) in cases
        @testset "$name" begin
            V = f(p, uv)
            @test length(V) == size(uv, 2)
            @test all(isfinite, V)
            # A non-negative brightness distribution normalised to unit flux cannot exceed
            # |V| = 1 anywhere, and must reach 1 at zero baseline.
            @test maximum(abs.(V)) <= 1.0 + 1e-6
            @test abs(V[zero_baseline]) ≈ 1.0 atol = 1e-3
        end
    end

    # visibility_annulus carries a `# bugged` annotation at vis_functions.jl:108. Its values
    # are in fact fine; what is wrong is the return type. The degenerate branch
    # (outer - inner < tol) returns a Complex vector while the normal branch returns a real
    # one, so the function is type-unstable and a caller that hits the degenerate case gets a
    # different element type than it was given for every other parameter value. Recorded as
    # broken so that fixing it surfaces here.
    @testset "annulus" begin
        V = visibility_annulus([1.0, 2.0], uv)
        @test length(V) == size(uv, 2)
        @test all(isfinite, V)
        @test maximum(abs.(V)) <= 1.0 + 1e-6
        @test abs(V[zero_baseline]) ≈ 1.0 atol = 1e-3
        @test_broken eltype(visibility_annulus([1.0, 1.0], uv)) == eltype(V)
    end

    # The NFFT geometry the old demo got wrong: the image grid and the uv grid have to be
    # consistent, or plan_nfft rejects the nodes outright.
    @testset "image inversion round trip" begin
        nx      = 64
        pixsize = 0.25                      # mas -> 16 mas field, larger than the 10 mas model
        scale   = pixsize * (pi / 180.0) / 3600000.0 .* [-1; 1] .* uv
        @test maximum(abs.(scale)) < 0.5     # NFFT node range; this is what used to fail

        plan = plan_nfft(scale, (nx, nx), m = 4, σ = 2.0)
        V    = Complex.(visibility_ud([10.0], uv))
        img  = real.(adjoint(plan) * V)
        @test size(img) == (nx, nx)
        @test all(isfinite, img)
        @test maximum(img) > 0
        # A centred uniform disc must invert to something centred and positive-dominated.
        pos = max.(img, 0)
        cy, cx = OITOOLS.cdg(pos)
        @test cx ≈ nx/2 atol = 2.0
        @test cy ≈ nx/2 atol = 2.0
    end

    # The flat-dict equivalent of the same model, which is the recommended API.
    @testset "flat dict model image" begin
        model_dict = Dict{String,Any}(
            "ring,fwhmin"  => 8.5,
            "ring,fwhmout" => 11.5,
            "ring,incl"    => 60.0,
            "ring,pa"      => -90.0,
            "ring,f"       => 0.5,
            "star,ud"      => 0.0,
            "star,f"       => "1.0 - \$ring,f",
        )
        model = dict_to_model(model_dict, String[])
        img   = model_to_image(model, Float64[]; nx = 128, pixsize = 0.25)
        @test size(img) == (128, 128)
        @test all(isfinite, img)
        @test all(img .>= 0)
        @test sum(img) ≈ 1.0 rtol = 1e-6     # model_to_image normalises
    end
end
