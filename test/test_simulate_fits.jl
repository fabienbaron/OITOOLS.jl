# Simulating from a FITS file: `image = path` is read by its keywords (`read_model_image`), so a
# 2-D image, a wavelength cube, a time series and a 4-D cube are each observed as what they are.
# Judged against the same file as a model component (`sky,image`) through the same uv points,
# which picks planes by wavelength and date independently of the simulator.
using OITOOLS, Test, Dates, Statistics
import AstroFITS

@testset "simulating from a FITS image or cube" begin
    d = mktempdir()
    src = joinpath(@__DIR__, "..", "demos", "data", "V1295_Aql.oifits")
    data = readoifits(src; warn = false, verbose = false, T = Float64)[1, 1]
    chans = sort(unique(data.uv_lam))                     # 5 channels
    _, nights = split_by_time(data)                       # 6 nights
    nx, ps = 65, 0.05                                     # odd: a flip keeps the centre pixel

    # a disc of diameter D (mas), East-first, centred on nx÷2+1, with a dot East of it so the
    # closure phases see the orientation
    function scene(D; dot = 0.3)
        c = nx ÷ 2 + 1
        img = [hypot(i - c, j - c) * ps < D / 2 ? 1.0 : 0.0 for i in 1:nx, j in 1:nx]
        img[c + 12, c + 5] += dot * sum(img)
        return img
    end
    sim(path; kw...) = (out = joinpath(d, "s$(rand(UInt32)).oifits");
                        simulate_from_oifits(src, out; mode = "copy_errors", noise = false, kw...);
                        readoifits(out; warn = false, verbose = false, T = Float64, filter_bad_data = false)[1, 1])
    asmodel(path) = (; flat_model = dict_to_model(Dict{String,Any}("sky,image" => path), String[]),
                       flat_params = Float64[])
    same(a, b) = isapprox(a.v2, b.v2; rtol = 1e-6, atol = 1e-9) &&
                 maximum(abs.(mod.(a.t3phi .- b.t3phi .+ 180, 360) .- 180)) < 1e-4

    @testset "2-D" begin
        p = joinpath(d, "grey.fits"); writefits(scene(1.2), p; pixsize = ps)
        @test same(sim(p; image = p), sim(p; asmodel(p)...))
    end

    @testset "3-D wavelength cube: planes by wavelength, not by position" begin
        D = [0.8, 1.0, 1.2, 1.4, 1.6]
        up = joinpath(d, "up.fits")
        writefits(cat((scene(x) for x in D)...; dims = 3), up; pixsize = ps, wavelengths = chans)
        down = joinpath(d, "down.fits")
        writefits(cat((scene(x) for x in reverse(D))...; dims = 3), down; pixsize = ps,
                  wavelengths = reverse(chans))
        a = sim(up; image = up)
        @test same(a, sim(up; asmodel(up)...))
        @test same(a, sim(down; image = down))
        # one disc for every channel would not match: the planes really differ
        g = joinpath(d, "flat.fits"); writefits(scene(1.2), g; pixsize = ps)
        @test !same(a, sim(g; image = g))
    end

    @testset "3-D time series and 4-D cube" begin
        D = range(0.8, 1.6; length = length(nights))
        ts = joinpath(d, "series.fits")
        AstroFITS.writefits!(ts, ["OBJECT" => "x"], cat((scene(x) for x in D)...; dims = 3))
        fits_set_keywords(ts, ts; pixsize = ps, axis3 = :time, mjds = nights)
        a = sim(ts; image = ts)
        @test same(a, sim(ts; asmodel(ts)...))
        cube = Array{Float64}(undef, nx, nx, length(chans), length(nights))
        for w in eachindex(chans), t in eachindex(nights)
            cube[:, :, w, t] = scene(0.6 + 0.1w + 0.1t)
        end
        q = joinpath(d, "four.fits"); writefits(cube, q; pixsize = ps, wavelengths = chans, mjds = nights)
        @test same(sim(q; image = q), sim(q; asmodel(q)...))
    end

    @testset "a standard file (East left) is not mirrored" begin
        img = scene(1.2)
        east = joinpath(d, "east.fits"); writefits(img, east; pixsize = ps)
        west = joinpath(d, "westleft.fits")
        AstroFITS.writefits!(west, Pair{String,Any}["CTYPE1" => "RA---SIN", "CDELT1" => -ps / 3.6e6,
                                                    "CUNIT1" => "deg", "CRPIX1" => Float64(nx ÷ 2 + 1),
                                                    "CTYPE2" => "DEC--SIN", "CDELT2" => ps / 3.6e6,
                                                    "CUNIT2" => "deg", "CRPIX2" => Float64(nx ÷ 2 + 1)],
                             reverse(img; dims = 1))
        a, b = sim(east; image = east), sim(west; image = west)
        @test same(a, b)
        @test maximum(abs.(a.t3phi)) > 1                   # the dot makes the phases non-zero
    end

    @testset "pixsize overrides the header" begin
        p = joinpath(d, "grey.fits")
        q = joinpath(d, "grey2.fits"); writefits(scene(1.2), q; pixsize = 2ps)
        @test same(sim(p; image = p, pixsize = 2ps), sim(q; image = q))
    end

    @testset "what is refused" begin
        bare = joinpath(d, "bare.fits"); AstroFITS.writefits!(bare, ["OBJECT" => "x"], rand(nx, nx, 3) .+ 0.1)
        err = try sim(bare; image = bare); nothing catch e; e end
        @test err isa ArgumentError && occursin("CTYPE3", err.msg)
        two = joinpath(d, "two.fits")
        writefits(cat(scene(1.0), scene(1.2); dims = 3), two; pixsize = ps, wavelengths = chans[[1, 4]])
        err = try sim(two; image = two); nothing catch e; e end
        @test err isa ArgumentError && occursin("not one of the image's channels", err.msg)
        @test_throws ArgumentError sim(two; image = rand(nx, nx, 2, 2))
    end

    @testset "model + image: the image as one more component" begin
        p = joinpath(d, "grey.fits")
        both = dict_to_model(Dict{String,Any}("star,ud" => 0.4, "star,f" => 1.0,
                                              "env,image" => p, "env,f" => 0.5), String[])
        a = sim(p; flat_model = both, flat_params = Float64[])
        star = dict_to_model(Dict{String,Any}("star,ud" => 0.4), String[])
        env  = dict_to_model(Dict{String,Any}("env,image" => p), String[])
        uv = a.uv; wl = a.uv_lam; mjd = a.uv_mjd
        V = (eval_model(star, Float64[], uv; wl, mjd) .+ 0.5 .* eval_model(env, Float64[], uv; wl, mjd)) ./ 1.5
        @test isapprox(a.v2, abs2.(V[a.indx_v2]); rtol = 1e-5, atol = 1e-8)
    end

    @testset "the SPARCO star + image template: simulated, then fitted back" begin
        f = read_model_file(joinpath(@__DIR__, "..", "demos", "models", "SPARCO_Star_Image.toml"))
        @test isfile(f.model["env,image"])                  # named relative to the model file
        @test sort(f.free) == ["env,slope", "fs"]
        truth = merge(copy(f.model), Dict{String,Any}("fs" => 0.45, "env,slope" => 1.2))
        out = joinpath(d, "sparco.oifits")
        simulate_from_oifits(src, out; mode = "copy_errors", noise = false,
                             flat_model = dict_to_model(truth, String[]), flat_params = Float64[])
        sdata = readoifits(out; warn = false, verbose = false, T = Float64)[1, 1]
        r = fit_model_lsqfit(f.model, f.free, sdata; lb = f.lb, ub = f.ub)
        i = Dict(zip(r.list_free_params, r.x_opt))
        @test i["fs"] ≈ 0.45 atol = 1e-4
        @test i["env,slope"] ≈ 1.2 atol = 1e-3
    end

    @testset "simulate: a wavelength cube on a night of its own" begin
        dates = collect(DateTime(2018, 8, 13, 3, 0, 0):Minute(30):DateTime(2018, 8, 13, 4, 0, 0))
        t = read_obs_file("default_obs"); t.target = "TESTSTAR"
        t.raep0 = 15 * (20 + 57 / 60); t.decep0 = 46 + 28 / 60
        f, c, w = read_facility_file("CHARA"), read_comb_file("MIRCX"), read_wave_file("MIRCX_LOWH")
        cube = cat((scene(0.6 + 0.1k) for k in eachindex(w.λ))...; dims = 3)
        p = joinpath(d, "mircx.fits"); writefits(cube, p; pixsize = ps, wavelengths = w.λ)
        simrun(; kw...) = (out = joinpath(d, "n$(rand(UInt32)).oifits");
                      simulate(f, t, c, w, dates, out; noise = false, kw...);
                      readoifits(out; warn = false, verbose = false, T = Float64, filter_bad_data = false)[1, 1])
        @test same(simrun(image = p), simrun(; asmodel(p)...))
        @test_throws ArgumentError simrun(image = rand(nx, nx, 2, 2))

        # the image's own channels: a cube on wavelengths no setup file has
        λc = [1.52e-6, 1.58e-6, 1.66e-6, 1.71e-6]
        q = joinpath(d, "own.fits")
        writefits(cat((scene(0.8 + 0.2k) for k in eachindex(λc))...; dims = 3), q; pixsize = ps, wavelengths = λc)
        wq = wave_from_image(q)
        @test wq.λ == λc && wq.δλ ≈ [0.06e-6, 0.07e-6, 0.065e-6, 0.05e-6]   # spacing to neighbours
        wb = joinpath(d, "own_band.fits")
        writefits(cat((scene(1.0) for _ in λc)...; dims = 3), wb; pixsize = ps, wavelengths = λc,
                  bandwidths = fill(0.03e-6, 4))
        @test wave_from_image(wb).δλ == fill(0.03e-6, 4)                   # EFF_BAND when given
        sim_img(; kw...) = (out = joinpath(d, "i$(rand(UInt32)).oifits");
                            simulate(f, t, c, :image, dates, out; noise = false, kw...);
                            readoifits(out; warn = false, verbose = false, T = Float64)[1, 1])
        a = sim_img(image = q)
        @test isapprox(sort(unique(a.uv_lam)), λc; rtol = 1e-6)          # OIFITS keeps Float32
        # a model holding the cube (a uniform-disc star + the cube): the same channels
        star_cube = dict_to_model(Dict{String,Any}("star,ud" => 0.5, "star,f" => 0.6,
                                                   "env,image" => q, "env,f" => 0.4), String[])
        @test wave_from_image(star_cube).λ == λc
        b = sim_img(flat_model = star_cube, flat_params = Float64[])
        @test isapprox(sort(unique(b.uv_lam)), λc; rtol = 1e-6) && b.v2 != a.v2
        # a model with no cube has no channels to give
        @test_throws ArgumentError wave_from_image(dict_to_model(Dict{String,Any}("star,ud" => 0.5), String[]))
        @test_throws ArgumentError sim_img(flat_model = dict_to_model(Dict{String,Any}("star,ud" => 0.5), String[]),
                                           flat_params = Float64[])
    end
end
