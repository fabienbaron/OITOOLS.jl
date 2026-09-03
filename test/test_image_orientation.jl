# Sky orientation of an image array, and of every display that draws one.
#
# The anchor is not a convention anyone chose in the plotting code: `setup_nfft` hands u to
# the NFFT as node row 1, and NFFT indexes the FIRST array dimension by that row, so a pixel
# at index c+k along dimension 1 IS a source `+k*pixsize` mas EAST. Dimension 2 is North the
# same way. Everything downstream has to agree with that, and for a while nothing checked it:
# both plotters paired a DESCENDING x coordinate vector with reversed axis limits, which drew
# eastward pixels on the western side and mirrored every image the package produced.

@testset "image orientation" begin
    MAS = 2.0626480624709636e8
    nx, ps = 32, 0.4
    c = nx ÷ 2 + 1                       # the NFFT's zero-frequency pixel
    data = readoifits(joinpath(@__DIR__, "..", "demos", "data", "2004-data1.oifits");
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
end
