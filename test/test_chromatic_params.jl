# Per-point (chromatic / time-variable) model parameters.
#
# An expression referencing `$WL`, `$MJD` or `$BL` resolves to one value per uv
# point rather than a scalar.  This must work for *geometric* parameters
# (a diameter that varies with wavelength, a limb darkening that varies with
# time, ...), not only for fluxes, and the gradient must follow.

using OITOOLS, Test, Random

@testset "per-point model parameters" begin

# a small, deterministic uv set spanning a few wavelengths
rng  = Xoshiro(42)
npt  = 24
uv   = 3e7 .* (2 .* rand(rng, 2, npt) .- 1)      # cycles/rad
wl   = repeat([1.50e-6, 1.60e-6, 1.70e-6], inner = npt ÷ 3)
mjd  = repeat([55753.0, 55754.0, 55755.0], inner = npt ÷ 3)

WL   = raw"$WL"
MJD  = raw"$MJD"
chrom(v0, slope) = "$v0 * (1 + $slope * ($WL/1.6 - 1))"

@testset "geometric parameters accept expressions" begin
    # every component kind whose geometry used to require a scalar
    cases = Dict(
        "ud"          => Dict{String,Any}("s,ud" => chrom(0.8, 0.05), "s,f" => 1.0),
        "ldlin"       => Dict{String,Any}("s,ldlin" => chrom(0.8, 0.05), "s,u" => 0.3, "s,f" => 1.0),
        "ldlin (u)"   => Dict{String,Any}("s,ldlin" => 0.8, "s,u" => chrom(0.3, 0.2), "s,f" => 1.0),
        "ldquad"      => Dict{String,Any}("s,ldquad" => chrom(0.8, 0.05), "s,u" => 0.2, "s,w" => 0.1, "s,f" => 1.0),
        "ldpow"       => Dict{String,Any}("s,ldpow" => chrom(0.8, 0.05), "s,alpha" => 0.5, "s,f" => 1.0),
        "ldquad (w)"  => Dict{String,Any}("s,ldquad" => 0.8, "s,u" => 0.2, "s,w" => chrom(0.1, 0.3), "s,f" => 1.0),
        "ldpow (α)"   => Dict{String,Any}("s,ldpow" => 0.8, "s,alpha" => chrom(0.5, 0.2), "s,f" => 1.0),
        "ring"        => Dict{String,Any}("r,diamin" => 0.5, "r,diamout" => chrom(1.0, 0.05), "r,f" => 1.0),
        "crescent"    => Dict{String,Any}("c,crin" => 0.6, "c,crout" => chrom(1.0, 0.05),
                                          "c,croff" => 0.5, "c,crprojang" => 30.0, "c,f" => 1.0),
        "gaussian"    => Dict{String,Any}("g,fwhm" => chrom(0.8, 0.05), "g,f" => 1.0),
        "gauss. ring" => Dict{String,Any}("g,fwhmin" => 0.2, "g,fwhmout" => chrom(0.5, 0.05), "g,f" => 1.0),
    )
    for (tag, md) in cases
        m = parse_model(md, String[])
        V = eval_model(m, Float64[], uv; wl=wl, mjd=mjd)
        @test length(V) == npt
        @test all(isfinite, V)
        @test all(abs.(V) .<= 1 + 1e-8)
    end
end

@testset "chromatic == per-wavelength scalar model" begin
    # the reference: evaluate a *scalar* model at each point's own wavelength
    m_c = parse_model(Dict{String,Any}("s,ud" => chrom(0.8, 0.05), "s,f" => 1.0), String[])
    V_c = eval_model(m_c, Float64[], uv; wl=wl, mjd=mjd)
    V_r = similar(V_c)
    for i in 1:npt
        d_i = 0.8 * (1 + 0.05 * (wl[i]/1.6e-6 - 1))
        m_s = parse_model(Dict{String,Any}("s,ud" => d_i, "s,f" => 1.0), String[])
        V_r[i] = eval_model(m_s, Float64[], uv[:, i:i]; wl=[wl[i]], mjd=[mjd[i]])[1]
    end
    @test V_c ≈ V_r
    # and a constant "chromatic" expression reproduces the plain scalar model
    m_0 = parse_model(Dict{String,Any}("s,ud" => 0.8, "s,f" => 1.0), String[])
    m_1 = parse_model(Dict{String,Any}("s,ud" => chrom(0.8, 0.0), "s,f" => 1.0), String[])
    @test eval_model(m_0, Float64[], uv; wl=wl, mjd=mjd) ≈
          eval_model(m_1, Float64[], uv; wl=wl, mjd=mjd)
end

@testset "time-variable parameters" begin
    m = parse_model(Dict{String,Any}(
            "s,ud" => "0.8 + 0.1 * ($MJD - 55753.0)", "s,f" => 1.0), String[])
    V = eval_model(m, Float64[], uv; wl=wl, mjd=mjd)
    @test all(isfinite, V)
    # the three epochs see three different diameters, so their visibilities differ
    @test !(V[1:8] ≈ V[9:16])
end

@testset "gradients follow (ForwardDiff vs finite differences)" begin
    md = Dict{String,Any}(
        "s,d0"    => 0.80,
        "s,slope" => 0.05,
        "s,ud"    => "\$s,d0 * (1 + \$s,slope * ($WL/1.6 - 1))",
        "s,f"     => 1.0,
    )
    free = ["s,d0", "s,slope"]
    m  = parse_model(md, free)
    x0 = Float64[md[p] for p in free]

    V, J = eval_model_grad(m, x0, uv; wl=wl, mjd=mjd)
    @test size(J) == (npt, length(x0))
    for j in eachindex(x0)
        h  = 1e-6
        xp = copy(x0); xp[j] += h
        xm = copy(x0); xm[j] -= h
        Jfd = (eval_model(m, xp, uv; wl=wl, mjd=mjd) .-
               eval_model(m, xm, uv; wl=wl, mjd=mjd)) ./ (2h)
        @test J[:, j] ≈ Jfd rtol=1e-5
    end
end

# The entry points that take an `OIdata` rather than an explicit `wl=`. These are where a
# chromatic model is actually used -- `model_to_residuals` and every `plot_*_residuals` go
# through `model_to_obs`, and `simulate_from_oifits` writes a file from `eval_model` -- and
# each has to take the wavelengths from the data itself. Handing the resolver no wavelength
# returns NaN for every `$WL` expression, which is not an error and shows up only as a plot of
# nothing; `model_to_obs` did exactly that, while `chi2_flat` passed `wl` and was correct.
#
# The assertion is agreement between the two paths, not a value written here: chi2 IS the sum
# of the squared normalised residuals over the weighted observables, so if `model_to_obs` and
# `chi2_flat` see the same model the two numbers are the same number.
@testset "chromatic models resolve against the data's own wavelengths" begin
    datadir = joinpath(@__DIR__, "..", "demos", "data")
    chromatic = Dict{String,Any}(
        "star,ud"   => 2.0, "star,f"   => "0.7 * ($WL/1.6)^-4",
        "disk,fwhm" => 6.0, "disk,f"   => "0.3 * ($WL/1.6)^-1.5")

    for f in (joinpath("Beauty_Contest_2004", "2004-data1.oifits"), "V1295_Aql.oifits")
        d  = readoifits(joinpath(datadir, f); warn = false)[1, 1]
        m  = dict_to_model(chromatic, String[])
        o  = model_to_obs(m, Float64[], d)
        r  = model_to_residuals(m, Float64[], d)

        @test d.nv2 == 0 || all(isfinite, o.v2)
        @test d.nv2 == 0 || all(isfinite, r.v2)
        @test d.nvisamp == 0 || all(isfinite, r.visamp)

        # The default weights are [1,1,1,0,0,0,0], so V2, T3amp and T3phi are the three that
        # enter chi2; visamp and visphi are computed but not weighted.
        res2 = sum(sum(abs2, getproperty(r, k)) for k in (:v2, :t3amp, :t3phi))
        @test res2 ≈ model_to_chi2(m, Float64[], d) rtol = 1e-8
    end

    # A grey model was never affected, and must not have moved.
    grey = Dict{String,Any}("star,ud" => 2.0, "star,f" => 0.7,
                            "disk,fwhm" => 6.0, "disk,f" => 0.3)
    d = readoifits(joinpath(datadir, "Beauty_Contest_2004", "2004-data1.oifits"); warn = false)[1, 1]
    m = dict_to_model(grey, String[])
    r = model_to_residuals(m, Float64[], d)
    @test sum(sum(abs2, getproperty(r, k)) for k in (:v2, :t3amp, :t3phi)) ≈
          model_to_chi2(m, Float64[], d) rtol = 1e-8
end

@testset "\$WL is in microns; every wl argument in metres" begin
    # A flux of (λ/1.6 µm)^-2: 1 at 1.6 µm, 1/4 at 3.2 µm, whichever path evaluates it.
    md = Dict{String,Any}("s,ud" => 1.0, "s,f" => "(\$WL/1.6)^-2")
    m = dict_to_model(md, String[])
    z = zeros(2, 2); w2 = [1.6e-6, 3.2e-6]
    @test model_flux(m, Float64[]; wl = w2) ≈ [1.0, 0.25]
    @test real.(eval_model(m, Float64[], z; wl = w2, normalize = false)) ≈ [1.0, 0.25]
    @test real.(eval_model_jac(m, Float64[], z; wl = w2, normalize = false)[1]) ≈ [1.0, 0.25]
    @test model_to_sed(m, Float64[], w2)[1] ≈ [1.0, 0.25]
    @test eltype(eval_model(dict_to_model(md, String[]; T = Float32), Float32[], Float32.(z);
                            wl = Float32.(w2))) === ComplexF32
    # a visibility function's $WL too
    mv = dict_to_model(Dict{String,Any}("c,visfunc" => "\$WL/1.6 + 0*\$BL", "c,f" => 1.0), String[])
    @test real.(eval_model(mv, Float64[], z; wl = w2)) ≈ [1.0, 2.0]

    # A model file written before $WL was in microns (no wl_unit) means what it meant
    f = tempname() * ".toml"
    write(f, "free = []\n[model]\n\"s,ud\" = 1.0\n\"s,f\" = \"(\$WL/1.6e-6)^-2\"\n")
    r = read_model_file(f)
    @test r.model["s,f"] == "((\$WL*1e-6)/1.6e-6)^-2"
    @test eval_model(dict_to_model(r.model, String[]), Float64[], z; wl = w2) ≈
          eval_model(m, Float64[], z; wl = w2)
    # and a file written now carries the marker, and reads back unchanged
    write_model_file(f, md)
    @test occursin("wl_unit = \"um\"", read(f, String))
    @test read_model_file(f).model["s,f"] == "(\$WL/1.6)^-2"

    # a constant in metres set against $WL is reported
    @test any(occursin("looks like metres", x)
              for x in model_warnings(Dict{String,Any}("s,ud" => 1.0, "s,f" => "(\$WL/1.6e-6)^-2"), String[]))
    @test !any(occursin("looks like metres", x) for x in model_warnings(md, String[]))
end

@testset "visibilities are normalised by the total flux" begin
    uvn = [0.0 1e7 3e7; 0.0 2e7 -1e7]
    wn  = [1.5e-6, 1.6e-6, 1.7e-6]
    fluxes = Dict{String,Any}("s,ud" => 1.0, "s,f" => 1.0, "d,fwhm" => 4.0, "d,f" => 1.5)
    fracs  = Dict{String,Any}("s,ud" => 1.0, "s,f" => 1 / 2.5, "d,fwhm" => 4.0, "d,f" => 1.5 / 2.5)
    V = eval_model(dict_to_model(fluxes, String[]), Float64[], uvn)
    @test V[1] ≈ 1
    @test V ≈ eval_model(dict_to_model(fracs, String[]), Float64[], uvn)
    @test eval_model(dict_to_model(fluxes, String[]), Float64[], uvn; normalize = false) ≈ 2.5 .* V
    @test model_flux(dict_to_model(fluxes, String[]), Float64[]) ≈ [2.5]

    # a chromatic flux: normalised at each point, by that point's total
    chromf = Dict{String,Any}("s,ud" => 1.0, "s,f" => "(\$WL/1.6)^-2", "d,fwhm" => 4.0, "d,f" => 1.0)
    m = dict_to_model(chromf, String[])
    fs = (wn ./ 1.6e-6) .^ -2
    Vs = eval_model(dict_to_model(Dict{String,Any}("s,ud" => 1.0), String[]), Float64[], uvn)
    Vd = eval_model(dict_to_model(Dict{String,Any}("d,fwhm" => 4.0), String[]), Float64[], uvn)
    @test eval_model(m, Float64[], uvn; wl = wn) ≈ (fs .* Vs .+ Vd) ./ (fs .+ 1)
    @test model_flux(m, Float64[]; wl = wn) ≈ fs .+ 1

    # the closed-form Jacobian against ForwardDiff: free fluxes, a derived one, a chromatic one
    for (md, free, x) in ((fluxes, ["s,f", "d,f", "d,fwhm"], [1.0, 1.5, 4.0]),
                          (merge(chromf, Dict{String,Any}("d,f" => "\$k * 2", "k" => 0.4)), ["k", "d,fwhm"], [0.4, 4.0]))
        mm = dict_to_model(md, free)
        V1, J1 = eval_model_jac(mm, x, uvn; wl = wn)
        V2, J2 = eval_model_grad_ad(mm, x, uvn; wl = wn)
        @test V1 ≈ V2 && J1 ≈ J2
    end

    # a fully resolved background: no visibility, its flux in the total, a uniform level in
    # an image (64² pixels; the star is a small Gaussian, without the ringing of a disc's
    # edge, so the corner is the background alone)
    bgm = dict_to_model(Dict{String,Any}("s,fwhm" => 0.5, "s,f" => 1.0, "bg,resolved" => true,
                                         "bg,f" => 0.25), String[])
    Vbg = eval_model(bgm, Float64[], uvn)
    @test Vbg[1] ≈ 1                                    # zero baseline: the background's own flux too
    @test Vbg[2:end] ≈ eval_model(dict_to_model(Dict{String,Any}("s,fwhm" => 0.5), String[]),
                                  Float64[], uvn[:, 2:end]) ./ 1.25
    img = model_to_image(bgm, Float64[]; nx = 64, pixsize = 0.1, normalize = false)
    @test sum(img) ≈ 1.25
    @test img[1, 1] ≈ 0.25 / 64^2 rtol = 1e-6
    @test sum(model_to_image(bgm, Float64[]; nx = 64, pixsize = 0.1)) ≈ 1

    # what the warnings say instead of "the fractions do not sum to one"
    @test any(occursin("every flux is free", x) for x in model_warnings(fluxes, ["s,f", "d,f"]))
    @test !any(occursin("every flux is free", x) for x in model_warnings(fluxes, ["d,f"]))
    @test any(occursin("must be positive", x)
              for x in model_warnings(merge(fluxes, Dict{String,Any}("d,f" => -1.0)), String[]))
end

end # testset
