# The Image perspective's data layer.
#
# The assertion that matters is §"a reconstruction with nothing but positivity still fits":
# a run that returns an image proves only that it returned. What proves the criterion, the
# Fourier plan and the weights are all connected is χ² going down a long way and the image
# staying non-negative — the second being a property of VMLMB's `lower = 0` rather than of any
# regulariser, which is why it holds with no regularisers at all.

@testset "Imaging data layer" begin

    IMGFILE = joinpath(@__DIR__, "data", "2004-data1.oifits")
    data = readoifits(IMGFILE; warn = false, verbose = false)

    @testset "geometry comes from the data, not from a constant" begin
        s = imaging_defaults(data; nx = 32)
        @test s.nx == 32
        @test s.pixsize ≈ auto_pixsize(data[1, 1])
        @test fov(s) ≈ 32 * s.pixsize
    end

    @testset "the starting image matches the plan's precision" begin
        # image_to_vis dispatches on the element type matching the plan's, so a Float64 start
        # against a Float32 plan is a MethodError from inside the criterion — one that says
        # nothing about precision. The cast belongs here, once.
        s  = ImagingSetup(; nx = 32, pixsize = 0.3)
        ft = setup_ft(data, s.nx, s.pixsize)
        @test ft isa OIft{Float32}
        for kind in (:dirac, :gaussian)
            x = start_image(ImagingSetup(; nx = 32, pixsize = 0.3, startkind = kind), ft)
            @test eltype(x) === Float32
            @test size(x) == (32, 32)
            @test sum(x) ≈ 1.0 rtol = 1e-5      # unit flux, both kinds
            @test all(x .>= 0)
        end
        # a Dirac really is one pixel
        d = start_image(ImagingSetup(; nx = 32, pixsize = 0.3, startkind = :dirac), ft)
        @test count(!iszero, d) == 1
        # and a Gaussian is not
        g = start_image(ImagingSetup(; nx = 32, pixsize = 0.3, startkind = :gaussian), ft)
        @test count(!iszero, g) > 1
        @test_throws ErrorException start_image(
            ImagingSetup(; nx = 32, pixsize = 0.3, startkind = :spiral), ft)
    end

    @testset "weights are three long, and say which three" begin
        # reconstruct takes [V², T3amp, T3φ] — three, not the seven model fitting uses.
        @test imaging_weights() == [1.0, 1.0, 1.0]
        @test length(imaging_weights()) == 3
        @test imaging_weights(t3amp = false) == [1.0, 0.0, 1.0]
        @test imaging_weights(v2 = false, t3amp = false, t3phi = false) == [0.0, 0.0, 0.0]
    end

    @testset "a reconstruction with nothing but positivity still fits" begin
        s = ImagingSetup(; nx = 48, pixsize = 0.3, startkind = :gaussian)
        r = reconstruct_image(data, s; maxiter = 120)

        @test r isa ImagingResult
        @test size(r.image) == (48, 48)
        # positivity is VMLMB's lower bound, in force with no regulariser asking for it
        @test all(r.image .>= 0)
        # and the fit actually descended, by a lot
        @test r.chi2r < r.chi2r_start
        @test r.chi2r < r.chi2r_start / 100
        @test isfinite(r.chi2r)
        @test r.seconds > 0
        @test r.setup === s
        @test r.weights == [1.0, 1.0, 1.0]
        @test occursin("χ²ᵣ", sprint(show, r))
    end

    @testset "more iterations do not make the fit worse" begin
        s = ImagingSetup(; nx = 32, pixsize = 0.35)
        short = reconstruct_image(data, s; maxiter = 30)
        long  = reconstruct_image(data, s; maxiter = 150)
        # VMLMB is a descent method from the same start, so the longer run cannot end higher
        @test long.chi2r <= short.chi2r * (1 + 1e-6)
        @test long.chi2r_start ≈ short.chi2r_start
    end

    @testset "switching an observable off changes the problem" begin
        s = ImagingSetup(; nx = 32, pixsize = 0.35)
        all3 = reconstruct_image(data, s; maxiter = 60)
        v2only = reconstruct_image(data, s; maxiter = 60,
                                   weights = imaging_weights(t3amp = false, t3phi = false))
        # Dropping the closure phases drops their contribution from chi2, so the two numbers
        # are not comparable — what is checked is that the weights reached the criterion at all.
        @test v2only.weights == [1.0, 0.0, 0.0]
        @test v2only.chi2r != all3.chi2r
        @test v2only.chi2r < v2only.chi2r_start
    end

    @testset "no observables is refused rather than fitted" begin
        # An all-zero weight vector gives a criterion with nothing in it, an ndof of zero, and
        # a meaningless chi2. Better to say so than to return one.
        s = ImagingSetup(; nx = 16, pixsize = 0.5)
        @test_throws ErrorException reconstruct_image(data, s;
                                        weights = imaging_weights(v2 = false, t3amp = false,
                                                                  t3phi = false))
    end

    @testset "total flux is reported because nothing constrains it" begin
        # V² and closure phase are both invariant under a global scaling, so an image fitted
        # from them alone is free to drift away from unit flux. Reporting it stops that reading
        # as a bug.
        s = ImagingSetup(; nx = 32, pixsize = 0.35)
        r = reconstruct_image(data, s; maxiter = 60)
        @test r.flux ≈ sum(r.image)
        @test r.flux > 0
    end
end
