# An image or a cube as a model component (`c,image`), judged against answers computed without
# it: the analytic uniform disc (alone, and as the planes of a cube at real data's uv points), an
# analytic point source for the position convention, and a fit recovering injected values.
using LinearAlgebra, Statistics

@testset "the image component" begin
    d   = mktempdir()
    MAS = 2.0626480624709636e8
    v1295 = joinpath(@__DIR__, "..", "demos", "data", "V1295_Aql.oifits")

    # A uniform disc rendered by pixel AREA (ss × ss subsamples per pixel), centred (cx, cy)
    # mas from the phase-centre pixel nx÷2+1; East is dimension 1.
    function ud_image(nx, ps, diam; cx = 0.0, cy = 0.0, ss = 8)
        c = nx ÷ 2 + 1; img = zeros(nx, nx); r2 = (diam / 2)^2
        sub = ((1:ss) .- (ss + 1) / 2) ./ ss
        for j in 1:nx, i in 1:nx, b in sub, a in sub
            x = (i - c + a) * ps - cx; y = (j - c + b) * ps - cy
            x^2 + y^2 <= r2 && (img[i, j] += 1)
        end
        return img ./ sum(img)
    end
    image_model(path, extra::Pair...) =
        dict_to_model(Dict{String,Any}("sky,image" => path, extra...), String[])
    eval0(m, uv; kw...) = eval_model(m, Float64[], uv; kw...)

    # Baselines along a few position angles out to 150 Mλ, where a 1.5 mas disc is well
    # past its first null.
    B  = range(5e6, 1.5e8; length = 40)
    uv = reduce(hcat, [[b * cos(t), b * sin(t)] for b in B, t in (0.0, 0.7, 2.1)])

    @testset "a rendered cube against analytic discs, at V1295 Aql's uv points" begin
        # Two nights of V1295 Aql (epochs 1 and 3) and two of its five channels (1.50, 1.65 µm):
        # a 2 × 2 cube of centred discs, a different diameter in each plane, whose visibilities
        # are known in closed form.
        data = readoifits(v1295; warn = false, verbose = false, T = Float64)[1, 1]
        parts, mjd = split_by_time(data)
        nights = parts[[1, 3]]
        chans = sort(unique(data.uv_lam))
        λ = chans[[1, 4]]
        D = [1.0 1.3; 1.6 1.9]                                # D[wavelength, night], mas
        ps, nx = 0.02, 128
        cube = Array{Float64}(undef, nx, nx, 2, 2)
        for w in 1:2, t in 1:2; cube[:, :, w, t] = ud_image(nx, ps, D[w, t]); end
        p = joinpath(d, "discs.fits")
        writefits(cube, p; pixsize = ps, wavelengths = λ, mjds = mjd[[1, 3]])
        m = image_model(p)
        ud(Dm, uvs) = eval0(dict_to_model(Dict{String,Any}("ud,ud" => Dm), String[]), uvs)
        for (t, x) in enumerate(nights)
            for (w, l) in enumerate(λ)
                i = findall(==(l), x.uv_lam)
                V = eval0(m, x.uv[:, i]; wl = x.uv_lam[i], mjd = x.uv_mjd[i])
                # the area rendering's edge and the pixel's transfer function, ~1e-3
                @test maximum(abs.(V .- ud(D[w, t], x.uv[:, i]))) < 2e-3
            end
            # the three channels the cube does not have: an error naming its channels
            i = findall(==(chans[2]), x.uv_lam)
            err = try eval0(m, x.uv[:, i]; wl = x.uv_lam[i], mjd = x.uv_mjd[i]); nothing
                  catch e; e end
            @test err isa ArgumentError && occursin("1.5", err.msg)
        end

        # Planes are chosen by each point's wavelength and epoch, not by its position in the
        # list: reversed points give reversed visibilities, and the same uv points on the other
        # night are the other night's plane.
        x, y = nights
        ix = findall(==(λ[1]), x.uv_lam); iy = findall(==(λ[2]), y.uv_lam)
        uvs = hcat(x.uv[:, ix], y.uv[:, iy])
        wls = vcat(x.uv_lam[ix], y.uv_lam[iy]); ts = vcat(x.uv_mjd[ix], y.uv_mjd[iy])
        V = eval0(m, uvs; wl = wls, mjd = ts)
        @test eval0(m, reverse(uvs; dims = 2); wl = reverse(wls), mjd = reverse(ts)) ≈ reverse(V)
        other = eval0(m, x.uv[:, ix]; wl = x.uv_lam[ix], mjd = fill(y.uv_mjd[1], length(ix)))
        @test maximum(abs.(other .- ud(D[1, 2], x.uv[:, ix]))) < 2e-3   # night 2's disc
        @test_throws ArgumentError eval0(m, x.uv[:, ix])                  # a cube needs wavelengths
    end

    @testset "a rendered disc is the analytic disc" begin
        p = joinpath(d, "ud.fits"); writefits(ud_image(128, 0.02, 1.5), p; pixsize = 0.02)
        V  = eval0(image_model(p), uv)
        Vt = eval0(dict_to_model(Dict{String,Any}("ud,ud" => 1.5), String[]), uv)
        # The area rendering's edge and the pixel's own transfer function, both ~1e-3 here.
        @test maximum(abs.(V .- Vt)) < 2e-3
        @test maximum(abs.(imag.(V))) < 1e-6                 # centred: no phase
    end

    @testset "x and y shift it, and +x is East" begin
        ps, nx, c = 0.05, 64, 33
        p0 = joinpath(d, "disc0.fits"); writefits(ud_image(nx, ps, 1.0), p0; pixsize = ps)
        p1 = joinpath(d, "disc1.fits")
        writefits(ud_image(nx, ps, 1.0; cx = 4ps, cy = -3ps), p1; pixsize = ps)
        shifted = dict_to_model(Dict{String,Any}("sky,image" => p0, "sky,x" => 4ps,
                                                 "sky,y" => -3ps), String[])
        @test maximum(abs.(eval0(shifted, uv) .- eval0(image_model(p1), uv))) < 1e-6  # NFFT

        dot = zeros(nx, nx); dot[c + 5, c] = 1.0             # 5 pixels EAST of centre
        pd = joinpath(d, "dot.fits"); writefits(dot, pd; pixsize = ps)
        east = dict_to_model(Dict{String,Any}("pt,f" => 1.0, "pt,x" => 5ps, "pt,y" => 0.0),
                             String[])
        @test maximum(abs.(eval0(image_model(pd), uv) .- eval0(east, uv))) < 1e-6
    end

    @testset "a cube gives each channel its own plane" begin
        ps, nx = 0.02, 128
        λ = [1.6e-6, 2.2e-6]
        cube = cat(ud_image(nx, ps, 1.0), ud_image(nx, ps, 1.8); dims = 3)
        p = joinpath(d, "cube.fits"); writefits(cube, p; pixsize = ps, wavelengths = λ)
        n  = size(uv, 2)
        wl = [isodd(i) ? λ[1] : λ[2] for i in 1:n]           # the channels interleaved
        V  = eval0(image_model(p), uv; wl)
        ud(D) = eval0(dict_to_model(Dict{String,Any}("ud,ud" => D), String[]), uv)
        Vt = [isodd(i) ? ud(1.0)[i] : ud(1.8)[i] for i in 1:n]
        @test maximum(abs.(V .- Vt)) < 2e-3

        # `relative` keeps the cube's spectrum: the brighter plane sums to 1, the other to its
        # share, so V(0) is the plane's flux.
        p2 = joinpath(d, "sed.fits")
        writefits(cat(ud_image(nx, ps, 1.0), 0.5 .* ud_image(nx, ps, 1.8); dims = 3), p2;
                  pixsize = ps, wavelengths = λ)
        zero2 = zeros(2, 2)
        @test real.(eval0(image_model(p2, "sky,norm" => "relative"), zero2; wl = λ)) ≈ [1.0, 0.5] atol = 1e-6
        @test real.(eval0(image_model(p2), zero2; wl = λ)) ≈ [1.0, 1.0] atol = 1e-6
        @test_throws ErrorException image_model(p2, "sky,norm" => "peak")
    end

    @testset "the visibilities are computed once per uv set" begin
        p = joinpath(d, "cache.fits"); writefits(ud_image(64, 0.05, 1.0), p; pixsize = 0.05)
        m = dict_to_model(Dict{String,Any}("sky,image" => p, "sky,f" => 0.7, "sky,x" => 0.1,
                                           "pt,f" => 0.3), ["sky,f", "sky,x", "pt,f"])
        cache = only(c for c in m.components if c isa OITOOLS.ImageSpec).cache
        V1 = eval_model(m, [0.7, 0.1, 0.3], uv)
        @test cache.nfft == 1
        V2 = eval_model(m, [0.6, 0.2, 0.4], uv)               # a fit's next step
        @test cache.nfft == 1 && V2 != V1
        _, J = eval_model_grad(m, [0.6, 0.2, 0.4], uv)       # the Jacobian, through the cache
        @test cache.nfft == 1 && all(isfinite, J)
        # V is normalised by the total flux F = 0.6 + 0.4 = 1: ∂V/∂f = (shifted image − V)/F
        @test J[:, 1] ≈ eval_model(m, [1.0, 0.2, 0.0], uv) .- eval_model(m, [0.6, 0.2, 0.4], uv)
        eval_model(m, [0.6, 0.2, 0.4], uv[:, 1:10])          # a new uv set is a new transform
        @test cache.nfft == 2
        eval_model(m, [0.6, 0.2, 0.4], Float32.(uv))         # the same points at Float32
        @test cache.nfft == 3
    end

    @testset "beyond the image's Nyquist frequency" begin
        p = joinpath(d, "coarse.fits"); writefits(ud_image(32, 0.5, 4.0), p; pixsize = 0.5)
        m = image_model(p)
        far = [1e6 3e8; 0.0 0.0]                              # 3e8 > 1/(2 × 0.5 mas)
        V = @test_logs (:warn, r"Nyquist") eval0(m, far)
        @test V[2] == 0 && abs(V[1]) > 0.9
        @test_logs eval0(m, far .* 1.01)                      # one warning per component
    end

    @testset "what is refused, and what the checkers accept" begin
        p = joinpath(d, "chk.fits"); writefits(ud_image(64, 0.05, 1.0), p; pixsize = 0.05)
        md = Dict{String,Any}("sky,image" => p, "sky,pixsize" => 0.05, "sky,norm" => "plane",
                              "sky,f" => "1 - \$pt,f", "pt,f" => 0.2, "pt,x" => 1.0,
                              "pt,y" => 0.0)
        @test_logs dict_to_model(md, ["pt,f", "pt,x"])      # no "keys not used" warning
        @test isempty(OITOOLS.missing_refs(md))
        @test isempty(model_warnings(md, ["pt,f"]))
        @test occursin("[image]", sprint(show, dict_to_model(md, String[])))
        out = mktemp() do path, io
            redirect_stdout(io) do; display_model(md, ["pt,f"]); end
            close(io); read(path, String)
        end
        @test occursin("[image]", out)
        @test all(OITOOLS.reserved_suffix, ("image", "pixsize", "norm"))

        # A pixel scale given in the model wins over the header.
        m = dict_to_model(merge(md, Dict{String,Any}("sky,pixsize" => 0.1)), String[])
        @test only(c for c in m.components if c isa OITOOLS.ImageSpec).image.pixsize == 0.1

        @test_throws ErrorException dict_to_model(merge(md, Dict{String,Any}("sky,pa" => 30.0)),
                                                  String[])
        @test_throws ErrorException dict_to_model(merge(md, Dict{String,Any}("sky,incl" => 30.0)),
                                                  String[])
        @test_throws ArgumentError dict_to_model(
            merge(md, Dict{String,Any}("sky,image" => joinpath(d, "nope.fits"))), String[])

        # A model file names its image relative to itself.
        sub = mkpath(joinpath(d, "models")); cp(p, joinpath(sub, "star.fits"))
        tf = joinpath(sub, "m.toml")
        write(tf, "[model]\n\"sky,image\" = \"star.fits\"\n\"sky,f\" = 1.0\n")
        r = cd(() -> read_model_file(tf), d)
        @test r.model["sky,image"] == joinpath(sub, "star.fits")
        @test dict_to_model(r.model, String[]) isa FlatModel
    end

    @testset "a fit recovers an injected companion" begin
        # Plan Part 12.3's pattern: simulate WITH noise from image + point, fit the point back,
        # and judge the error against σ from the χ² curvature.
        src = joinpath(@__DIR__, "..", "demos", "data", "Beauty_Contest_2004", "2004-data1.oifits")
        p = joinpath(d, "primary.fits"); writefits(ud_image(64, 0.05, 1.5), p; pixsize = 0.05)
        truth = Dict{String,Any}("star,image" => p, "star,f" => "1 - \$pt,f",
                                 "pt,f" => 0.2, "pt,x" => 3.0, "pt,y" => -2.0)
        out = joinpath(d, "sim.oifits")
        simulate_from_oifits(src, out; flat_model = dict_to_model(truth, String[]),
                             flat_params = Float64[], seed = 11)
        data = readoifits(out; T = Float64, filter_bad_data = true)[1, 1]
        free = ["pt,f", "pt,x", "pt,y"]
        start = merge(truth, Dict{String,Any}("pt,f" => 0.15, "pt,x" => 2.8, "pt,y" => -1.8))
        # LD_MMA: the default LD_LBFGS returns FAILURE without a step on this sharp χ²
        # (a companion position free at high SNR), with an analytic primary as well.
        res = fit_model(start, free, data; method = :LD_MMA)
        # The χ² Hessian by central differences (the χ² is not twice differentiable by
        # ForwardDiff), with steps well inside each parameter's σ.
        χ²(x) = model_to_chi2(res.model, x, data)
        h = [1e-4, 1e-3, 1e-3]; e(i) = (v = zeros(3); v[i] = h[i]; v)
        H = [(χ²(res.x_opt .+ e(i) .+ e(j)) - χ²(res.x_opt .+ e(i) .- e(j)) -
              χ²(res.x_opt .- e(i) .+ e(j)) + χ²(res.x_opt .- e(i) .- e(j))) / (4h[i] * h[j])
             for i in 1:3, j in 1:3]
        σ = sqrt.(diag(inv(H ./ 2)))
        @test all(abs.(res.x_opt .- [0.2, 3.0, -2.0]) .< 4 .* σ)
        @test 0.5 < res.chi2r < 1.5
    end
end
