# Images and cubes as model components: the FITS convention shared with ROTIR, its writer and its
# reader. The reader is what a model component will be built from, so what is asserted here is
# the CONVENTION — East-first, phase centre at nx÷2+1, planes matched by wavelength and epoch —
# and not merely that a file comes back as it went in.
using FITSIO

@testset "model images: the FITS convention shared with ROTIR" begin
    d = mktempdir()
    hdr(path) = FITSIO.read_header(FITS(path)[1])
    nx = 16; c = nx ÷ 2 + 1
    unit_planes(a) = (b = copy(a); for k in CartesianIndices(size(b)[3:end])
                          b[:, :, k] ./= sum(b[:, :, k]); end; b)

    @testset "writefits writes the shared header" begin
        img = unit_planes(rand(nx, nx))
        p = joinpath(d, "grey.fits"); writefits(img, p; pixsize = 0.25)
        h = hdr(p)
        @test h["CTYPE1"] == "RA---SIN" && h["CTYPE2"] == "DEC--SIN"
        @test h["CUNIT1"] == "deg" && h["CUNIT2"] == "deg"
        @test h["CDELT1"] > 0 && h["CDELT2"] > 0              # dimension 1 runs East
        @test h["CDELT1"] * 3.6e6 ≈ 0.25
        @test h["CRPIX1"] == c && h["CRPIX2"] == c            # the NFFT's phase centre
        @test h["PIXSIZE"] == 0.25
        @test h["NORMALIZ"] == "PLANE"
        @test fits_pixsize(p) ≈ 0.25 rtol = 1e-10
        @test !haskey(h, "CTYPE3")                            # a grey image has no third axis

        raw = rand(nx, nx) .+ 1
        q = joinpath(d, "unnormalised.fits"); writefits(raw, q; pixsize = 0.25)
        @test !haskey(hdr(q), "NORMALIZ")                     # only said when true
    end

    @testset "round trips: grey, wavelength, epoch, both" begin
        lam_lin = [1.5e-6, 1.6e-6, 1.7e-6]
        lam_irr = [1.5e-6, 1.55e-6, 1.7e-6]
        mjd_irr = [58000.1, 58001.2, 58009.0]

        g = unit_planes(rand(nx, nx))
        p = joinpath(d, "g.fits"); writefits(g, p; pixsize = 0.3)
        mi = read_model_image(p)
        @test size(mi.image) == (nx, nx, 1, 1)
        @test mi.image[:, :, 1, 1] ≈ g
        @test mi.pixsize ≈ 0.3 && mi.normalization === :plane && mi.offset == (0.0, 0.0)
        @test all(isnan, mi.wavelengths) && all(isnan, mi.mjds)

        cube = unit_planes(rand(nx, nx, 3, 1))
        p = joinpath(d, "lin.fits"); writefits(cube, p; pixsize = 0.3, wavelengths = lam_lin)
        h = hdr(p)
        @test h["CTYPE3"] == "WAVE" && h["CUNIT3"] == "m" && haskey(h, "CDELT3")
        @test read_model_image(p).wavelengths == lam_lin

        p = joinpath(d, "irr.fits")
        writefits(cube, p; pixsize = 0.3, wavelengths = lam_irr, bandwidths = fill(2e-8, 3))
        h = hdr(p)
        @test h["WAVEIRR"] == true && !haskey(h, "CDELT3")  # no linear axis that misplaces planes
        mi = read_model_image(p)
        @test mi.wavelengths == lam_irr && mi.bandwidths == fill(2e-8, 3)

        both = unit_planes(rand(nx, nx, 2, 3))
        p = joinpath(d, "both.fits")
        writefits(both, p; pixsize = 0.3, wavelengths = lam_lin[1:2], mjds = mjd_irr)
        h = hdr(p)
        @test h["CTYPE4"] == "TIME" && h["CUNIT4"] == "d" && h["TIMESYS"] == "UTC"
        @test h["TIMEIRR"] == true
        mi = read_model_image(p)
        @test size(mi.image) == (nx, nx, 2, 3)
        @test mi.image ≈ both
        @test mi.wavelengths == lam_lin[1:2] && mi.mjds == mjd_irr

        # Not normalised in the file: each plane by its own SUM, which keeps a slightly
        # negative band-limited plane's shape and flux.
        raw = rand(nx, nx, 2, 1) .+ 0.5; raw[1, 1, :, 1] .= -0.01
        p = joinpath(d, "raw.fits"); writefits(raw, p; pixsize = 0.3, wavelengths = lam_lin[1:2])
        mi = read_model_image(p)
        @test [sum(mi.image[:, :, k, 1]) for k in 1:2] ≈ [1.0, 1.0]
        @test mi.image[1, 1, 1, 1] < 0                      # never clipped
        @test mi.image[:, :, 1, 1] ≈ raw[:, :, 1, 1] ./ sum(raw[:, :, 1, 1])
    end

    @testset "orientation: East stays East" begin
        k = 3
        # Our own file: a point k pixels EAST of the phase centre.
        img = zeros(nx, nx); img[c + k, c] = 1.0
        p = joinpath(d, "east.fits"); writefits(img, p; pixsize = 0.3)
        @test Tuple(argmax(read_model_image(p).image[:, :, 1, 1])) == (c + k, c)

        # A third-party image in the usual astronomical layout: CDELT1 < 0, so the array index
        # increases WEST, and the same EAST point sits at c - k. It must come back at c + k.
        wf(name, keys, vals, a; comments = fill("", length(keys))) =
            (path = joinpath(d, name); f = FITS(path, "w");
             write(f, a; header = FITSHeader(keys, Any[vals...], comments)); close(f); path)
        # Where a pixel IS, in mas East of the origin: its index from the phase centre, plus
        # the offset the reader carries for a reference pixel off that centre.
        east_of(mi) = (Tuple(argmax(mi.image[:, :, 1, 1]))[1] - c) * mi.pixsize + mi.offset[1]
        west_first = zeros(nx, nx); west_first[c - k, c] = 1.0
        q = wf("astro.fits", ["CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2", "CDELT1", "CDELT2",
                              "CRPIX1", "CRPIX2"],
               ["RA---TAN", "DEC--TAN", "deg", "deg", -0.3 / 3.6e6, 0.3 / 3.6e6, c, c],
               west_first)
        # Flipping an even-sized array moves its phase-centre pixel by one, which the reader
        # carries as an offset: the point is k pixels East either way.
        @test east_of(read_model_image(q)) ≈ k * 0.3
    end

    @testset "a reference pixel elsewhere becomes an offset" begin
        img = unit_planes(rand(nx, nx))
        p = joinpath(d, "crpix.fits"); writefits(img, p; pixsize = 0.5)
        f = FITS(p, "r+"); write_key(f[1], "CRPIX1", Float64(c + 2)); close(f)
        mi = read_model_image(p)
        # The origin is the file's reference pixel, two pixels EAST of the phase centre, so the
        # image content sits two pixels WEST of the origin.
        @test mi.offset[1] ≈ -2 * 0.5 && mi.offset[2] == 0.0
    end

    @testset "what is refused, and why" begin
        cube3 = rand(nx, nx, 3) .+ 0.5
        f = FITS(joinpath(d, "nax3.fits"), "w"); write(f, cube3); close(f)
        err = try read_model_image(joinpath(d, "nax3.fits"); pixsize = 0.3); nothing
              catch e; e end
        @test err isa ArgumentError && occursin("CTYPE3", err.msg)

        f = FITS(joinpath(d, "time3.fits"), "w")
        write(f, cube3; header = FITSHeader(["CTYPE3"], Any["TIME"], [""])); close(f)
        @test size(read_model_image(joinpath(d, "time3.fits"); pixsize = 0.3).image) ==
              (nx, nx, 1, 3)

        dark = zeros(nx, nx); dark[1] = -1.0
        p = joinpath(d, "dark.fits"); writefits(dark, p; pixsize = 0.3)
        @test_throws ArgumentError read_model_image(p)

        nopix = joinpath(d, "nopix.fits"); writefits(rand(nx, nx), nopix)
        @test_throws ArgumentError read_model_image(nopix)
        @test read_model_image(nopix; pixsize = 0.2).pixsize == 0.2
    end

    @testset "ROTIR's reference cube, and matching planes to data" begin
        dir  = joinpath(@__DIR__, "data", "rotir")
        mi   = read_model_image(joinpath(dir, "zetand_ref_cube.fits"))
        @test size(mi.image) == (144, 144, 2, 2)
        @test mi.pixsize == 0.10699553372976                # the header's PIXSIZE, exactly
        @test mi.wavelengths == [1.515382109573693e-6, 1.6227540982072242e-6]
        @test mi.mjds == [58696.37925102794, 58701.39097318102]
        @test mi.normalization === :plane && mi.offset == (0.0, 0.0)
        @test all(abs(sum(mi.image[:, :, i, j]) - 1) < 1e-6 for i in 1:2, j in 1:2)
        @test minimum(mi.image) < 0                          # the limb's ringing is kept

        files = sort(filter(endswith(".oifits"), readdir(dir; join = true)))
        for (it, f) in enumerate(files)
            data = readoifits(f; polychromatic = true, warn = false, verbose = false)
            @test size(data, 1) == 6
            for ch in 1:6
                x  = data[ch, 1]
                λ  = Float64(x.uv_lam[1]); t = Float64(x.uv_mjd[1])
                if ch in (1, 4)
                    # The cube's channels 1 and 4 — on the second night too, although that
                    # night's channels sit 1.2e-4 away from the first's.
                    @test model_image_plane(mi, λ, t) == (ch == 1 ? 1 : 2, it)
                else
                    @test_throws ArgumentError model_image_plane(mi, λ, t)
                end
            end
        end
        @test_throws ArgumentError model_image_plane(mi, 1.515382e-6, 58650.0)  # no such night
        # A grey image serves every wavelength and every time.
        grey = joinpath(d, "anyλ.fits"); writefits(unit_planes(rand(nx, nx)), grey; pixsize = 0.3)
        @test model_image_plane(read_model_image(grey), 2.2e-6, 50000.0) == (1, 1)
    end
end
