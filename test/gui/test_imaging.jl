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
    # Five H-band channels, which is enough to be a cube and small enough to reconstruct in a
    # test. `POLY` in the harness is BC2026's 119-channel N-band file, far too slow for this.
    POLYFILE = joinpath(@__DIR__, "..", "..", "demos", "data",
                        "2019_v1295Aql.WL_SMOOTH.A.oifits")

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
            # `(nx, nx, nwav, nepoch)`, one channel here: the engines see one shape whether the
            # dataset has bins or not, so no reader has to ask which it was handed.
            @test size(x) == (32, 32, 1, 1)
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
        @test size(r.image) == (48, 48, 1, 1)
        @test result_channels(r) == 1
        @test size(result_plane(r)) == (48, 48)
        # positivity is VMLMB's lower bound, in force with no regulariser asking for it
        @test all(r.image .>= 0)
        # and the fit actually descended, by a lot
        @test chi2r(r) < chi2r_start(r)
        @test chi2r(r) < chi2r_start(r) / 100
        @test isfinite(chi2r(r))
        # `image_to_chi2` returns the RAW chi2; the reduced one divides by the points actually
        # fitted. Confusing the two turns a good fit (0.7) into an apparent disaster (315).
        @test r.chi2 ≈ chi2r(r) * r.ndof
        @test r.ndof == sum(b -> b.used ? b.n : 0, r.breakdown)
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
        @test chi2r(long) <= chi2r(short) * (1 + 1e-6)
        @test chi2r_start(long) ≈ chi2r_start(short)
    end

    @testset "switching an observable off changes the problem" begin
        s = ImagingSetup(; nx = 32, pixsize = 0.35)
        all3 = reconstruct_image(data, s; maxiter = 60)
        v2only = reconstruct_image(data, s; maxiter = 60,
                                   weights = imaging_weights(t3amp = false, t3phi = false))
        # Dropping the closure phases drops their contribution from chi2, so the two numbers
        # are not comparable — what is checked is that the weights reached the criterion at all.
        @test v2only.weights == [1.0, 0.0, 0.0]
        @test chi2r(v2only) != chi2r(all3)
        @test chi2r(v2only) < chi2r_start(v2only)
        # only the fitted observable counts toward ndof
        @test v2only.ndof == only(filter(b -> b.name == "V²", v2only.breakdown)).n
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

    @testset "the breakdown says which observable is fitted and which is predicted" begin
        # One aggregate number cannot say that a V²-only fit leaves the closure phases wild,
        # and that is exactly the thing worth seeing.
        s = ImagingSetup(; nx = 32, pixsize = 0.35)
        r = reconstruct_image(data, s; maxiter = 60,
                              weights = imaging_weights(t3amp = false, t3phi = false))
        by = Dict(b.name => b for b in r.breakdown)
        @test by["V²"].used
        @test !by["T3φ"].used
        @test all(isfinite(b.chi2r) for b in r.breakdown)
        # the fitted one is fitted; the unfitted one is merely predicted, and much worse
        @test by["V²"].chi2r < by["T3φ"].chi2r
    end

    # The sampling engines return a distribution, and what they return is worth naming
    # precisely: `reconstruct_squeeze` gives one posterior mean PER CHAIN and returns the best
    # of them, so the ensemble mean and the reconstruction are different images. The panel
    # offers both, and this pins that they really are different.
    @testset "an ensemble is recognised, and only where there is one" begin
        # A point estimate has no ensemble, whatever else the engine attached.
        @test result_ensemble(nothing) === nothing
        @test result_ensemble((; history = NamedTuple[])) === nothing          # BSMEM/BSDMM
        @test result_ensemble((; params = [1.0], param_names = Dict())) === nothing  # SPARCO

        # Per-chain images become mean, spread and members.
        imgs = [fill(Float64(i), 4, 4) for i in 1:4]
        e = result_ensemble((; images = imgs, best_chain = 2))
        @test length(e.samples) == 4
        @test e.source == "4 chains"
        @test e.mean[1, 1] ≈ 2.5
        @test e.sigma[1, 1] ≈ sqrt(sum(abs2, [1.0,2.0,3.0,4.0] .- 2.5) / 3)   # sample std

        # One chain has a mean but no spread; saying "0" would claim a precision that was
        # never measured.
        @test result_ensemble((; images = [imgs[1]])).sigma === nothing

        # Ragged input is refused rather than broadcast into a wrong answer.
        @test result_ensemble((; images = [zeros(4, 4), zeros(5, 5)])) === nothing

        # The engine names its members. SQUEEZE returns one mean per CHAIN and VI returns
        # posterior DRAWS around one centre; the spread means a different thing in each, so
        # the panel must not call both "chains".
        @test result_ensemble((; images = imgs,
                                 ensemble_noun = "posterior sample")).source ==
              "4 posterior samples"
        @test result_ensemble((; images = [imgs[1]],
                                 ensemble_noun = "posterior sample")).source ==
              "1 posterior sample"
    end

    # VI is four algorithms behind one engine entry, chosen by an option. The console echoes
    # the call that ran, so the name has to follow the option rather than the engine.
    @testset "the console names the VI algorithm that ran" begin
        @test engine_call_name(:vmlmb) == "reconstruct"
        @test engine_call_name(:vi) == "reconstruct_hybrid"          # the default
        for v in ("map", "mgvi", "geovi", "hybrid")
            @test engine_call_name(:vi, Dict("vi_engine" => v)) == "reconstruct_" * v
        end
        # A name from nowhere falls back to the entry rather than inventing a function.
        @test engine_call_name(:vi, Dict("vi_engine" => "nonsense")) == "reconstruct_hybrid"
        # Options belonging to another engine change nothing.
        @test engine_call_name(:squeeze, Dict("vi_engine" => "map")) == "reconstruct_squeeze"
    end

    @testset "tempering is offered, and says what it needs" begin
        # The engine is nameable whether or not Pigeons is installed: the panel greys it with a
        # reason, which it can only do if the entry exists.
        @test haskey(IMAGING_ENGINES, :tempering)
        @test engine_call_name(:tempering) == "reconstruct_squeeze_tempered"
        # `shell_optional_engines` is what the panel binds that greying to.
        @test occursin("tempering=", shell_optional_engines())

        # No run, nothing to report — and specifically "" rather than a row of zeros, because
        # a swap acceptance of 0 means a dead rung and printing that for "no data" would
        # invent a fault.
        sh = SHELL[]
        if sh !== nothing
            keep = sh.imaging
            sh.imaging = nothing
            @test shell_tempering_diagnostics() == ""
            sh.imaging = keep
        end
    end

    # ── the wavelength cube ──────────────────────────────────────────────────
    #
    # VMLMB is the one engine that reconstructs one image per spectral bin. These assert the
    # SHAPE and the BOOKKEEPING of that path; whether the pictures are any good is a question
    # for a ground-truth test, and the identity below is the cheap half of one.
    @testset "a cube goes in and a cube comes out" begin
        poly = readoifits(POLYFILE; warn = false, verbose = false, polychromatic = true,
                          merge_oi_wavelength = true, use_vis = false)
        nwav = size(poly, 1)
        @test nwav > 1
        s  = ImagingSetup(; engine = :vmlmb, nx = 24, pixsize = 0.5, startkind = :gaussian)
        ft = setup_ft(poly, s.nx, s.pixsize)

        x0 = start_image(s, ft)
        @test size(x0) == (24, 24, nwav, 1)
        # EVERY channel, not the total: the criterion normalises per cell, so one empty plane
        # divides by zero and the whole run comes back NaN with nothing saying which channel.
        @test all(sum(x0[:, :, w, 1]) ≈ 1 for w in 1:nwav)

        r = reconstruct_image(poly, s; maxiter = 20)
        @test size(r.image) == (24, 24, nwav, 1)
        @test result_channels(r) == nwav
        @test size(result_plane(r, 2)) == (24, 24)
        @test result_plane(r, 99) == result_plane(r, nwav)      # clamped, not thrown
        @test length(r.wavelengths) == nwav
        @test issorted(r.wavelengths)

        # The regression for a χ² that was wrong by roughly the bin count: the reduced χ²
        # must not scale with how the SAME data was split.
        mono = readoifits(POLYFILE; warn = false, verbose = false, use_vis = false)
        rm = reconstruct_image(mono, ImagingSetup(; engine = :vmlmb, nx = 24, pixsize = 0.5,
                                                    startkind = :gaussian); maxiter = 20)
        @test r.ndof == rm.ndof                       # the same points, however they are binned
        @test all(isfinite(b.chi2r) for b in r.breakdown)
    end

    @testset "cross-channel regularisers reach the engine, or are refused" begin
        spatial, trans = split_regularizers(parse_regularizers(
            "l1l2,1e-3,1e-6;transspectral_tv,1e-2;tv,1e-3"))
        @test length(spatial) == 2 && length(trans) == 1
        @test first(trans[1]) == "transspectral_tv"

        poly = readoifits(POLYFILE; warn = false, verbose = false, polychromatic = true,
                          merge_oi_wavelength = true, use_vis = false)
        s = ImagingSetup(; engine = :vmlmb, nx = 24, pixsize = 0.5, startkind = :gaussian)
        a = reconstruct_image(poly, s; maxiter = 30, regularizers = "l1l2,1e-1,1e-6")
        b = reconstruct_image(poly, s; maxiter = 30,
                              regularizers = "l1l2,1e-1,1e-6;transspectral_tv,1e2")
        @test !(a.image ≈ b.image)                    # the cross-channel term did something

        # Refused rather than dropped: an engine that cannot apply one must say so, or the user
        # believes they constrained the spectrum and nothing did.
        @test_throws ErrorException reconstruct_image(
            poly[1:1, :], ImagingSetup(; engine = :bsmem, nx = 24, pixsize = 0.5);
            maxiter = 2, regularizers = "transspectral_tv,1e-2")
    end

    @testset "a saved cube says what its third axis means" begin
        poly = readoifits(POLYFILE; warn = false, verbose = false, polychromatic = true,
                          merge_oi_wavelength = true, use_vis = false)
        lams = bin_wavelengths(poly)
        @test length(lams) == size(poly, 1)
        @test all(1e-7 .< lams .< 1e-4)               # metres, not µm and not m⁻¹

        cube = rand(8, 8, length(lams), 1) .+ 0.1
        path = joinpath(mktempdir(), "cube.fits")
        writefits(cube, path; pixsize = 0.5, wavelengths = lams)
        h = OITOOLS.FITSIO.read_header(OITOOLS.FITSIO.FITS(path)[1])
        @test h["NAXIS3"] == length(lams)
        @test h["CTYPE3"] == "WAVE"
        @test h["CUNIT3"] == "m"
        # the axis reconstructs the centres it was given
        got = h["CRVAL3"] .+ (0:length(lams)-1) .* h["CDELT3"]
        @test all(isapprox.(got, lams; rtol = 1e-6))

        # A grey image gains no spectral axis: a third axis nothing can interpret is worse
        # than none.
        grey = joinpath(mktempdir(), "grey.fits")
        writefits(rand(8, 8), grey; pixsize = 0.5)
        @test !haskey(OITOOLS.FITSIO.read_header(OITOOLS.FITSIO.FITS(grey)[1]), "CTYPE3")
    end
end
