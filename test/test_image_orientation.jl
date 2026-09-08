# Sky orientation of an image array, and of every display that draws one.
#
# The anchor is not a convention anyone chose in the plotting code: `setup_nfft` hands u to
# the NFFT as node row 1, and NFFT indexes the FIRST array dimension by that row, so a pixel
# at index c+k along dimension 1 IS a source `+k*pixsize` mas EAST. Dimension 2 is North the
# same way. Everything downstream has to agree with that, and for a while nothing checked it:
# both plotters paired a DESCENDING x coordinate vector with reversed axis limits, which drew
# eastward pixels on the western side and mirrored every image the package produced.

# `FITS`/`FITSHeader` for the pixel-size testset at the bottom: this file is included by
# runtests.jl, which does not load FITSIO itself.
using FITSIO

@testset "image orientation" begin
    MAS = 2.0626480624709636e8
    nx, ps = 32, 0.4
    c = nx ÷ 2 + 1                       # the NFFT's zero-frequency pixel
    data = readoifits(joinpath(@__DIR__, "..", "demos", "data", "BC2004", "2004-data1.oifits");
                       filter_bad_data = true, T = Float64)
    ft = setup_ft(data, nx, ps)
    cell = ft[1, 1]
    uv = data[1,1].uv

    @testset "array dimension 1 is East, dimension 2 is North" begin
        for k in (3, -5)
            off = k * ps
            img = zeros(nx, nx); img[c + k, c] = 1.0
            V = OITOOLS.image_to_vis(img, cell.uv)
            east = cis.(-2π / MAS .* (uv[1, :] .* (+off)))    # parse_model's shift, ra = East
            west = cis.(-2π / MAS .* (uv[1, :] .* (-off)))
            @test maximum(abs.(V .- east)) < 1e-6
            @test maximum(abs.(V .- west)) > 1.0

            img2 = zeros(nx, nx); img2[c, c + k] = 1.0
            V2 = OITOOLS.image_to_vis(img2, cell.uv)
            north = cis.(-2π / MAS .* (uv[2, :] .* (+off)))
            @test maximum(abs.(V2 .- north)) < 1e-6
        end
    end

    @testset "imdisp places an eastward pixel at positive x" begin
        # `imdisp` draws `reverse(rotl90(image), dims=2)` under extent [+X,-X,-Y,+Y], with
        # matplotlib's default upper origin. Replicate that placement arithmetic: a cell at
        # column j of an m-column array sits at x = left + (j-0.5)/m * (right-left).
        img = zeros(nx, nx); img[c + 4, c] = 1.0        # 1.6 mas EAST
        M = reverse(rotl90(img), dims = 2)
        i, j = Tuple(findfirst(==(1.0), M))
        X = 0.5 * nx * ps
        x = X + (j - 0.5) / size(M, 2) * (-2X)
        @test x > 0                                      # positive x is East on that axis

        # North is on the other axis and must not have moved with the fix.
        imgn = zeros(nx, nx); imgn[c, c + 4] = 1.0       # 1.6 mas NORTH
        Mn = reverse(rotl90(imgn), dims = 2)
        i2, _ = Tuple(findfirst(==(1.0), Mn))
        y = X + (i2 - 0.5) / size(Mn, 1) * (-2X)
        @test y > 0
        @test Tuple(findfirst(==(1.0), rotl90(imgn)))[1] == i2   # unchanged by the reverse
    end

    @testset "coordinate vectors ascend with the array index" begin
        # What `imdisp_makie`, the GUI's `show_image!` and the snapshot figure all build. The
        # axis limits are reversed separately, which is what puts East on the left; the vector
        # itself must ascend, or the two cancel and the image comes out mirrored.
        half = nx * ps / 2
        xs = range(-half, half; length = nx)
        ys = range(-half, half; length = nx)
        @test xs[c + 4] > xs[c]         # further East  -> larger x
        @test ys[c + 4] > ys[c]         # further North -> larger y
        @test issorted(xs) && issorted(ys)
    end

    # The other half of the same mapping: how large a pixel IS. The Observing panel fills its
    # pixel-size box from the header rather than making the user retype a number the file
    # already carries, so the reader has to agree with the writer — and with everyone else's
    # writer, which uses different units.
    @testset "fits_pixsize reads what writefits wrote, and what others write" begin
        d = mktempdir()
        img = rand(16, 16)
        wf(name, keys, vals, comments) = begin
            path = joinpath(d, name)
            f = FITS(path, "w")
            write(f, img; header = FITSHeader(keys, Any[vals...], comments))
            close(f)
            path
        end

        # This package's own file: radians, and since this change also a CUNIT saying so.
        own = joinpath(d, "own.fits"); writefits(img, own; pixsize = 0.25)
        @test fits_pixsize(own) ≈ 0.25 rtol = 1e-5

        # An OLDER file of ours has the unit only in CDELT's comment. Still readable, or every
        # image written before today would come back 5.7e4 times too small.
        old_ours = wf("old_ours.fits", ["CDELT1", "CDELT2"], [-1.2120e-9, 1.2120e-9],
                      ["Radians per Pixel", "Radians per Pixel"])
        @test fits_pixsize(old_ours) ≈ 0.25 rtol = 1e-3

        # Everyone else: degrees, which is what CUNIT says and what the FITS standard assumes
        # when it is absent.
        @test fits_pixsize(wf("deg.fits", ["CDELT1", "CDELT2", "CUNIT1", "CUNIT2"],
                              [-1.0e-7, 1.0e-7, "deg", "deg"],
                              ["", "", "", ""])) ≈ 0.36
        @test fits_pixsize(wf("nounit.fits", ["CDELT1", "CDELT2"], [-1.0e-7, 1.0e-7],
                              ["", ""])) ≈ 0.36
        @test fits_pixsize(wf("as.fits", ["CDELT2", "CUNIT2"], [5.0e-4, "arcsec"],
                              ["", ""])) ≈ 0.5
        @test fits_pixsize(wf("mas.fits", ["CDELT2", "CUNIT2"], [0.3, "mas"],
                              ["", ""])) ≈ 0.3
        # The CD matrix form, which is what a WCS-aware writer often produces instead.
        @test fits_pixsize(wf("cd.fits", ["CD1_1", "CD2_2", "CUNIT2"],
                              [-2.0e-8, 2.0e-8, "deg"], ["", "", ""])) ≈ 0.072

        # And nothing rather than a guess: no WCS at all, an absurd value, a missing file.
        bare = joinpath(d, "bare.fits"); writefits(img, bare)
        @test fits_pixsize(bare) === nothing
        @test fits_pixsize(wf("absurd.fits", ["CDELT2"], [1.0e30], [""])) === nothing
        @test fits_pixsize(joinpath(d, "not_here.fits")) === nothing
    end
end
