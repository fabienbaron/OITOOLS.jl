# The two diagnostic panels under the fits table: residuals and the SED.
#
# Both take numbers straight from OITOOLS -- `model_to_residuals` and `model_to_sed` -- so what
# is checked here is the mapping from those into the plots, which is where the panel can lie
# without erroring: the wrong baseline field against an observable, a wavelength axis in the
# wrong unit, or a legend entry left over from the previous model.

@testset "Model diagnostics" begin

    datadir = joinpath(@__DIR__, "..", "..", "demos", "data")
    data = readoifits(joinpath(datadir, "2004-data1.oifits"); warn = false)[1, 1]

    md = Dict{String,Any}("star,ud" => 6.5, "star,f" => 0.7,
                          "disk,fwhm" => 12.0, "disk,f" => 0.3)
    model = dict_to_model(md, String[])
    res   = model_to_residuals(model, Float64[], data)

    @testset "each residual is plotted against its own baseline field" begin
        series = Dict(k => (pts, rms) for (k, pts, rms) in residual_series(data, res))

        # V² against `v2_baseline`, closures against `t3_baseline` -- the geometric mean, not
        # the longest leg. Taking the field from OBS_SPECS is what keeps this agreeing with
        # the observable plots, where the same choice was easy to get wrong.
        v2pts = series[:v2][1]
        @test length(v2pts) == data.nv2
        @test [p[1] for p in v2pts] ≈ data.v2_baseline .* 1e-6
        @test [p[2] for p in v2pts] ≈ res.v2

        t3pts = series[:t3phi][1]
        @test length(t3pts) == data.nt3phi
        @test [p[1] for p in t3pts] ≈ data.t3_baseline .* 1e-6

        # rms is the number the caption reports, and the ±1/±3 lines are drawn in the same
        # units, so the two have to be the same quantity.
        @test series[:v2][2] ≈ sqrt(sum(abs2, res.v2) / length(res.v2))
    end

    @testset "an observable the file does not carry draws nothing" begin
        series = Dict(k => pts for (k, pts, _) in residual_series(data, res))
        @test data.nvisamp == 0
        @test isempty(series[:visamp])
        @test isempty(series[:visphi])
    end

    @testset "mismatched lengths draw nothing rather than mismatched pairs" begin
        # A residual vector and its baselines can only disagree if the data changed under the
        # model. Pairing them anyway would plot points that belong to different observations.
        short = merge(res, (; v2 = res.v2[1:end-1]))
        series = Dict(k => pts for (k, pts, _) in residual_series(data, short))
        @test isempty(series[:v2])
        @test length(series[:t3phi]) == data.nt3phi
    end

    @testset "non-finite residuals are dropped, not drawn at zero" begin
        y = copy(res.v2); y[1] = NaN; y[2] = Inf
        series = Dict(k => (pts, rms) for (k, pts, rms) in residual_series(data, merge(res, (; v2 = y))))
        @test length(series[:v2][1]) == data.nv2 - 2
        @test all(isfinite(p[2]) for p in series[:v2][1])
        @test isfinite(series[:v2][2])
    end

    @testset "the SED axis is in microns while the model is evaluated in metres" begin
        # `model_to_sed` takes metres, as the resolver and `uv_lam` do. A grid in microns
        # evaluates every $WL expression a million times too far out and returns numbers that
        # are wrong without erroring, so the conversion belongs at the drawing.
        chromatic = Dict{String,Any}("star,ud" => 2.0,
                                     "star,f"  => "0.7 * (\$WL/1.6e-6)^-4",
                                     "disk,fwhm" => 6.0,
                                     "disk,f"  => "0.3 * (\$WL/1.6e-6)^-1.5")
        m  = dict_to_model(chromatic, String[])
        wl = [1.5e-6, 1.6e-6, 1.7e-6]
        total, comps = model_to_sed(m, Float64[], wl)
        @test total ≈ [0.7 * (1.5/1.6)^-4 + 0.3 * (1.5/1.6)^-1.5,
                       1.0,
                       0.7 * (1.7/1.6)^-4 + 0.3 * (1.7/1.6)^-1.5]

        fig = Makie.Figure()
        panel = build_sed(fig, Makie.Axis(fig[1, 1]))
        line = update_sed!(panel, wl, total, comps)

        @test [p[1] for p in panel.total[]] ≈ [1.5, 1.6, 1.7]
        @test [p[2] for p in panel.total[]] ≈ total
        @test occursin("2 components", line)
        @test panel.labels[1][] == "disk"       # name order, so colours do not reshuffle
        @test panel.labels[2][] == "star"
    end

    @testset "unused SED slots carry no label and no colour" begin
        # The Legend is built once with every entry -- it cannot be created after the window
        # exists -- so an unused slot is emptied rather than removed. A leftover swatch would
        # read as a component that failed to draw.
        fig = Makie.Figure()
        panel = build_sed(fig, Makie.Axis(fig[1, 1]))
        wl = [1.5e-6, 1.6e-6]
        update_sed!(panel, wl, [1.0, 1.0], Dict("a" => [0.5, 0.5], "b" => [0.5, 0.5]))
        for i in 3:GUI.SED_MAX_COMPONENTS
            @test isempty(panel.comps[i][])
            @test strip(panel.labels[i][]) == ""
            @test panel.colors[i][].alpha == 0
        end
        @test panel.colors[1][].alpha > 0

        # A second model with fewer components must not leave the first one's entries behind.
        update_sed!(panel, wl, [1.0, 1.0], Dict("only" => [1.0, 1.0]))
        @test panel.labels[1][] == "only"
        @test strip(panel.labels[2][]) == ""
        @test panel.colors[2][].alpha == 0
    end
end
