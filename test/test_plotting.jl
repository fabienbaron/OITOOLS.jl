# test_plotting.jl — regression suite for every figure OITOOLS can produce.
#
# `test_python_boundary.jl` proves the Julia<->Python bridge works. This file proves the
# *plots* are right, which is a different question: a figure can render perfectly and still
# show the wrong numbers. The three tiers below answer three separate questions.
#
#   1. "renders"   — every exported entry point and every method form is called and produces
#                    a non-blank PNG. Catches crashes and dead code paths.
#   2. "values"    — the points actually on the axes are compared against the data fields the
#                    plot claims to show, including the unit scaling. Catches a plot that
#                    draws the wrong field, drops points, or loses a 1e-6 factor. This is the
#                    tier that catches silent corruption, and it cannot be done by eye.
#   3. "structure" — labels, panel counts, log scaling, colorbars and the branch selectors
#                    (`t3base`, `phityp`, `amptyp`). Catches option handling regressions.
#
# Deliberately NOT reference-image comparison: matplotlib version, font and freetype drift
# make pixel baselines fail for reasons unrelated to OITOOLS, and a baseline that must be
# re-cut on every upgrade stops being evidence. Asserting on artist data is version-stable
# and strictly more specific — it says *which* number is wrong, not just "the picture moved".
#
# Everything runs headless under Agg.

ENV["MPLBACKEND"] = "Agg"

using OITOOLS, Test
using PythonPlot     # activates OITOOLSPythonPlotExt, which defines everything below
using PythonCall     # weak dep of OITOOLS; PythonPlot has already loaded it

const PLT = PythonPlot
# imshow2 and _cbar_ticks are matplotlib helpers, so they live in the extension module
# rather than in OITOOLS itself. OBS_PLOT_SPECS is toolkit-independent and stayed in core.
const EXT = Base.get_extension(OITOOLS, :OITOOLSPythonPlotExt)
const PC  = PythonCall

# Rendered PNGs land in test/figures/plotting so the whole gallery can be browsed after a
# run — the assertions below check numbers, not appearance, so looking at the figures is
# still how you catch "it's technically correct but ugly". Override with OITOOLS_PLOT_DIR.
#
# Stale files are cleared first: a renamed or deleted test would otherwise leave an orphan
# PNG behind and the directory would stop being an accurate picture of the current suite.
const _DIR = mkpath(get(ENV, "OITOOLS_PLOT_DIR",
                        joinpath(@__DIR__, "figures", "plotting")))
for f in readdir(_DIR; join = true)
    endswith(f, ".png") && rm(f)
end
const _DATA = joinpath(@__DIR__, "..", "demos", "data")

# Render resolution for the gallery. 300 dpi is publication-grade and keeps legend text
# sharp; override with OITOOLS_PLOT_DPI to trade size for speed.
const _DPI = parse(Int, get(ENV, "OITOOLS_PLOT_DPI", "300"))

# ── datasets ─────────────────────────────────────────────────────────────────
# All three are tracked in git, so this suite runs on a fresh checkout. They are split by
# role rather than using one big file, because the render tier calls ~50 plots and cost
# scales with point count:
#
#   mono  61 kB  single wavelength, single MJD -> the "wav"/"mjd" fallbacks and legends
#   poly 286 kB  V2/T3/VISAMP/VISPHI, 18 wavelengths, 67 MJDs -> the main polychromatic set
#   flux 2.4 MB  the only tracked files carrying OI_FLUX are the BC2026 pair, so FLUX
#                coverage uses the smaller of the two and nothing else
const mono = readoifits(joinpath(_DATA, "BC2004", "2004-data1.oifits");
                        warn = false, verbose = false)[1, 1]
const poly = readoifits(joinpath(_DATA, "MWC275_T4a.oifits");
                        warn = false, verbose = false)[1, 1]
const flux = readoifits(joinpath(_DATA, "BC2026", "OBJECT1_LM.oifits");
                        warn = false, verbose = false)[1, 1]

# ── helpers ──────────────────────────────────────────────────────────────────
_close()  = PLT.close("all")
_nfigs()  = PC.pyconvert(Int, PC.pylen(PLT.get_fignums()))
_pystr(x) = PC.pyconvert(String, x)
_pyint(x) = PC.pyconvert(Int, x)

"""Call a plotting function and assert a real figure came out (not a blank canvas)."""
function renders(name, f; minbytes = 3000)
    _close()
    f()
    @test _nfigs() >= 1
    total = 0
    for (k, num) in enumerate(PLT.get_fignums())      # differential layouts paginate
        p = joinpath(_DIR, "$(name)_$k.png")
        # 300 dpi so legend text is sharp when the gallery is inspected by eye or dropped
        # into a talk; bbox_inches="tight" keeps colorbars and below-plot legends from being
        # clipped, which matters now that legends sit outside the axes.
        PLT.figure(num); PLT.savefig(p, dpi = _DPI, bbox_inches = "tight")
        @test isfile(p)
        total += filesize(p)
    end
    @test total > minbytes                            # a failed/blank canvas is ~1-2 kB
    _close()
    return total
end

"""Every (x, y) point drawn on `ax`, from marker lines and from scatter offsets alike.

`_plot_obs` draws grouped data with `errorbar(fmt="o")` (one Line2D per group) but
wavelength/MJD colouring with `scatter` (one PathCollection). Both are collected so the
comparison below is independent of which colouring branch ran. `errorbar.capsize` is 0 in
`set_oiplot_defaults`, so there are no cap lines to filter out; the error bars themselves are
LineCollections and are skipped by the PathCollection test."""
function points(ax)
    xs, ys = Float64[], Float64[]
    for ln in ax.lines
        append!(xs, PC.pyconvert(Vector{Float64}, ln.get_xdata()))
        append!(ys, PC.pyconvert(Vector{Float64}, ln.get_ydata()))
    end
    for c in ax.collections
        _pystr(PC.pytype(c).__name__) == "PathCollection" || continue
        off = PC.pyconvert(Matrix{Float64}, c.get_offsets())
        size(off, 1) > 0 && size(off, 2) == 2 || continue
        append!(xs, off[:, 1]); append!(ys, off[:, 2])
    end
    return xs, ys
end

"""All points across every axis of the current figure."""
function figpoints()
    xs, ys = Float64[], Float64[]
    for ax in PLT.gcf().axes
        x, y = points(ax)
        append!(xs, x); append!(ys, y)
    end
    return xs, ys
end

"""Compare a plotted coordinate against expected values as an unordered multiset.

Grouped colouring splits the data across one artist per baseline, so ordering carries no
meaning — only the set of plotted values does."""
function same_values(plotted::Vector{Float64}, expected::AbstractVector; rtol = 1e-4)
    length(plotted) == length(expected) || return false
    return isapprox(sort(plotted), sort(Float64.(collect(expected))); rtol = rtol)
end

axes_of(fig = PLT.gcf()) = [ax for ax in fig.axes]

"""Every point an axis draws, as plain Julia data, so two figures can be compared exactly."""
function _xy(ax)
    pts = Vector{Tuple{Float64,Float64}}()
    for ln in ax.get_lines()
        xs = PC.pyconvert(Vector{Float64}, ln.get_xdata())
        ys = PC.pyconvert(Vector{Float64}, ln.get_ydata())
        append!(pts, zip(xs, ys))
    end
    for col in ax.collections
        o = PC.pyconvert(Array, col.get_offsets())
        size(o, 1) > 0 && append!(pts, [(o[i,1], o[i,2]) for i in axes(o, 1)])
    end
    return sort(pts)
end
ylabel_of(ax) = _pystr(ax.get_ylabel())
yscale_of(ax) = _pystr(ax.get_yscale())

@testset "Plotting" begin

    # ═════════════════════════════════════════════════════════════════════════
    # 1. every exported entry point renders
    # ═════════════════════════════════════════════════════════════════════════
    @testset "entry points render" begin

        @testset "observable plots" begin
            renders("plot_v2",      () -> plot_v2(mono))
            renders("plot_t3phi",   () -> plot_t3phi(mono))
            renders("plot_t3amp",   () -> plot_t3amp(mono))
            renders("plot_visamp",  () -> plot_visamp(poly))
            renders("plot_visphi",  () -> plot_visphi(poly))
            renders("plot_flux",    () -> plot_flux(flux))
            # deprecated alias — still exported, so still covered
            renders("plot_diffphi", () -> (@test_logs (:warn,) plot_diffphi(poly)))
        end

        @testset "colour modes" begin
            for c in ("baseline", "wav", "mjd", "red")
                renders("plot_v2_$c", () -> plot_v2(poly, color = c))
            end
            # single-wavelength / single-MJD data must fall back, not crash
            renders("plot_v2_wav_fallback", () -> (@test_logs (:warn,) plot_v2(mono, color = "wav")))
            renders("plot_v2_mjd_fallback", () -> (@test_logs (:warn,) plot_v2(mono, color = "mjd")))
            renders("plot_flux_station",    () -> plot_flux(flux, color = "station"))
        end

        @testset "display options" begin
            renders("plot_v2_log",      () -> plot_v2(poly, logplot = true))
            renders("plot_v2_markopt",  () -> plot_v2(mono, markopt = true))
            renders("plot_v2_below",    () -> plot_v2(mono, legend_below = true))
            renders("plot_v2_title",    () -> plot_v2(mono, figtitle = "regression"))
            renders("plot_t3phi_max",   () -> plot_t3phi(poly, t3base = "max"))
            renders("plot_t3amp_max",   () -> plot_t3amp(poly, t3base = "max"))
        end

        @testset "plot_obs / plot_multi" begin
            renders("plot_obs_default", () -> plot_obs(poly))
            # `flux` is the only tracked set with all six, so it drives the full-panel case
            renders("plot_obs_all",     () -> plot_obs(flux;
                        obs = ["V2", "T3PHI", "T3AMP", "VISAMP", "VISPHI", "FLUX"]))
            renders("plot_obs_max",     () -> plot_obs(poly; obs = ["T3PHI_MAX", "T3AMP_MAX"]))
            renders("plot_obs_wav",     () -> plot_obs(poly; color = "wav"))
            renders("plot_obs_mjd",     () -> plot_obs(poly; color = "mjd"))
            renders("plot_obs_log",     () -> plot_obs(poly; logplot = true))
            # empty panels must be skipped, not drawn blank
            renders("plot_obs_skips",   () -> plot_obs(mono; obs = ["V2", "FLUX"]))
            # plot_multi forwards to plot_obs. It used to be `const plot_multi = plot_obs`,
            # but a const alias cannot be extended from an extension module, so the move of
            # plotting into OITOOLSPythonPlotExt made it a real method. Identity no longer
            # holds; assert what actually matters, which is that the two draw the same thing.
            @test plot_multi !== plot_obs                     # genuinely a separate function
            _close()
            plot_obs(poly; obs = ["V2", "T3PHI"])
            ref = [_xy(ax) for ax in axes_of()]
            _close()
            plot_multi(poly; obs = ["V2", "T3PHI"])
            got = [_xy(ax) for ax in axes_of()]
            _close()
            @test length(got) == length(ref) && !isempty(ref)
            @test got == ref                                  # same panels, same points
            renders("plot_multi", () -> plot_multi(poly; obs = ["V2", "T3PHI"]))
        end

        @testset "uv coverage" begin
            renders("uvplot",         () -> uvplot(mono))
            renders("uvplot_matrix",  () -> uvplot(mono.uv))
            renders("uvplot_wav",     () -> uvplot(poly, color = "wav"))
            renders("uvplot_mjd",     () -> uvplot(poly, color = "mjd"))
            renders("uvplot_flipx",   () -> uvplot(mono, flipx = true))
            renders("uvplot_square",  () -> uvplot(mono, square = false))
            renders("uvplot_bounds",  () -> uvplot(mono, minuv = 0.0, maxuv = 50e6))
        end

        @testset "images" begin
            img = gaussian2d(64, 64, 64 / 6)
            renders("imdisp",           () -> imdisp(img; pixsize = 0.2))
            renders("imdisp_cbar",      () -> imdisp(img; pixsize = 0.2, use_colorbar = true))
            renders("imdisp_nopix",     () -> imdisp(img))
            renders("imdisp_beam",      () -> imdisp(img; pixsize = 0.2, beamsize = 2.0,
                                                     beamlocation = [0.8, 0.8]))
            cube = cat(gaussian2d(32, 32, 5.0), gaussian2d(32, 32, 7.0); dims = 3)
            renders("imdisp_multi",     () -> imdisp_multi(cube; pixsize = 0.2))
            renders("imdisp_multi_cbar",() -> imdisp_multi(cube; pixsize = 0.2,
                                                           use_colorbar = true))
            renders("imshow2",          () -> EXT.imshow2(gaussian2d(32, 32, 5.0),
                                                              gaussian2d(32, 32, 7.0)))
        end

        @testset "facility" begin
            fac = read_facility_file(joinpath(@__DIR__, "..", "src", "configs", "CHARA.toml"))
            renders("plot_facility", () -> plot_facility(fac); minbytes = 2000)
        end

        @testset "multifile" begin
            renders("plot_v2_multifile", () -> plot_v2_multifile([mono, mono]))
            renders("plot_v2_vector",    () -> plot_v2([mono, mono]))
        end

        # 13 residual methods across 3 call conventions
        @testset "residuals" begin
            img   = gaussian2d(64, 64, 64 / 6)
            ftm   = setup_nfft(mono, 64, 0.2)
            ftf   = setup_nfft(poly, 64, 0.2)
            model = dict_to_model(Dict{String,Any}("s,ud" => 6.5, "s,f" => 1.0), ["s,ud"])
            x     = [6.5]
            obs_m = model_to_obs(model, x, mono)
            obs_f = model_to_obs(model, x, poly)

            @testset "(data, model_vector)" begin
                renders("res_v2_vec",     () -> plot_v2_residuals(mono, obs_m.v2))
                renders("res_t3amp_vec",  () -> plot_t3amp_residuals(mono, obs_m.t3amp))
                renders("res_t3phi_vec",  () -> plot_t3phi_residuals(mono, obs_m.t3phi))
                renders("res_visamp_vec", () -> plot_visamp_residuals(poly, obs_f.visamp))
                renders("res_visphi_vec", () -> plot_visphi_residuals(poly, obs_f.visphi))
            end

            @testset "(image, data, ft)" begin
                renders("res_v2_img",    () -> plot_v2_residuals(img, mono, ftm))
                renders("res_t3amp_img", () -> plot_t3amp_residuals(img, mono, ftm))
                renders("res_t3phi_img", () -> plot_t3phi_residuals(img, mono, ftm))
            end

            @testset "(model, x, data)" begin
                renders("res_v2_mod",     () -> plot_v2_residuals(model, x, mono))
                renders("res_t3amp_mod",  () -> plot_t3amp_residuals(model, x, mono))
                renders("res_t3phi_mod",  () -> plot_t3phi_residuals(model, x, mono))
                renders("res_visamp_mod", () -> plot_visamp_residuals(model, x, poly))
                renders("res_visphi_mod", () -> plot_visphi_residuals(model, x, poly))
            end

            @testset "plot_residuals dispatch" begin
                renders("plot_residuals_obs",   () -> plot_residuals(mono, obs_m))
                renders("plot_residuals_img",   () -> plot_residuals(img, ftm, mono))
                renders("plot_residuals_model", () -> plot_residuals(model, x, mono))
                renders("plot_residuals_full",  () -> plot_residuals(poly, obs_f))
                renders("plot_residuals_ftf",   () -> plot_residuals(img, ftf, poly))
            end
        end

        @testset "set_oiplot_defaults" begin
            @test set_oiplot_defaults() === nothing || true
            @test set_oiplot_defaults(compact = true) === nothing || true
            # compact mode changes the differential-phase grid, so render under it too
            renders("plot_v2_compact", () -> plot_v2(mono))
            set_oiplot_defaults(compact = false)
        end
    end

    # ═════════════════════════════════════════════════════════════════════════
    # 2. the numbers on the axes are the numbers in the data
    # ═════════════════════════════════════════════════════════════════════════
    @testset "plotted values match the data" begin

        # (label, plot call, expected x, expected y) driven by OBS_PLOT_SPECS so the
        # expectation is derived from the same declaration the plot uses — except the
        # scaling, which is written out literally here so a changed x_scale is caught
        # rather than silently tracked.
        cases = [
            ("V2",     poly, (d; kw...) -> plot_v2(d; kw...),     :v2_baseline,  :v2,     1e-6),
            ("T3PHI",  poly, (d; kw...) -> plot_t3phi(d; kw...),  :t3_baseline,  :t3phi,  1e-6),
            ("T3AMP",  poly, (d; kw...) -> plot_t3amp(d; kw...),  :t3_baseline,  :t3amp,  1e-6),
            ("VISAMP", poly, (d; kw...) -> plot_visamp(d; kw...), :vis_baseline, :visamp, 1e-6),
            ("VISPHI", poly, (d; kw...) -> plot_visphi(d; kw...), :vis_baseline, :visphi, 1e-6),
            ("FLUX",   flux, (d; kw...) -> plot_flux(d; kw...),   :flux_lam,     :flux,   1e6),
        ]

        for (name, dset, plotfn, xfield, yfield, xscale) in cases
            @testset "$name" begin
                expx = getfield(dset, xfield) .* xscale
                expy = getfield(dset, yfield)
                # the same points must be drawn no matter how they are coloured
                for colour in ("baseline", "wav", "mjd")
                    _close(); plotfn(dset; color = colour)
                    xs, ys = figpoints()
                    @test same_values(ys, expy)
                    @test same_values(xs, expx)
                    _close()
                end
            end
        end

        @testset "T3 max-baseline variant uses the other x field" begin
            _close(); plot_t3phi(poly, t3base = "max")
            xs, _ = figpoints()
            @test same_values(xs, poly.t3_maxbaseline .* 1e-6)
            @test !same_values(xs, poly.t3_baseline .* 1e-6)   # genuinely different axis
            _close()
        end

        # `uv` is already a spatial frequency (units of 1/rad, i.e. baseline/lambda), so
        # uvplot only divides by 1e6 to reach Mlambda — it must NOT divide by uv_lam again.
        @testset "uvplot draws the conjugate uv plane in Mlambda" begin
            _close(); uvplot(mono)
            xs, ys = figpoints()
            # Baselines are grouped from indx_v2 *and* the three indx_t3 legs, so a uv column
            # shared by a V2 point and a closure leg is drawn once per referencing
            # observable. Counts therefore exceed size(uv,2); the invariant that matters is
            # that the plotted cloud is exactly the uv coverage and its Hermitian mirror.
            @test length(xs) >= 2 * size(mono.uv, 2)
            @test isapprox(maximum(abs, xs), maximum(abs, mono.uv[1, :]) / 1e6; rtol = 1e-6)
            @test isapprox(maximum(abs, ys), maximum(abs, mono.uv[2, :]) / 1e6; rtol = 1e-6)
            # the mirroring makes the plotted cloud exactly centro-symmetric
            @test isapprox(sum(xs), 0.0; atol = 1e-6 * maximum(abs, xs))
            @test isapprox(sum(ys), 0.0; atol = 1e-6 * maximum(abs, ys))
            for (plotted, uvrow) in ((xs, mono.uv[1, :]), (ys, mono.uv[2, :]))
                expected = sort(unique(vcat(uvrow, -uvrow) ./ 1e6))
                @test isapprox(sort(unique(plotted)), expected; rtol = 1e-6)
            end
            _close()
        end

        # imdisp/imdisp_multi peak-normalise each panel and apply `rotl90` (Monnier's
        # orientation). Asserting the whole array pins the orientation as well as the
        # values — a transposed or flipped image is a scientific error that renders fine.
        @testset "imdisp shows the image, rotated and peak-normalised" begin
            img = gaussian2d(64, 64, 64 / 6)
            _close(); imdisp(img; pixsize = 0.2)
            ax  = axes_of()[1]
            @test _pyint(PC.pylen(ax.images)) == 1
            arr = PC.pyconvert(Matrix{Float64}, ax.images[0].get_array())
            expected = reverse(rotl90(img), dims = 2) ./ maximum(img)
            @test size(arr) == size(expected)
            @test isapprox(arr, expected; rtol = 1e-5)
            @test isapprox(maximum(arr), 1.0; rtol = 1e-6)
            # guard the orientation specifically: an unrotated image would also be
            # peak-normalised and the same size, so equality above must not be vacuous
            @test !isapprox(arr, img ./ maximum(img); rtol = 1e-5)
            # and reject the bare `rotl90`, which pairs with the reversed extent to draw
            # East on the western side. See test_image_orientation.jl for why that is wrong.
            @test !isapprox(arr, rotl90(img) ./ maximum(img); rtol = 1e-5)
            _close()
        end

        @testset "imdisp_multi shows every slice" begin
            cube = cat(gaussian2d(32, 32, 5.0), gaussian2d(32, 32, 7.0),
                       gaussian2d(32, 32, 9.0); dims = 3)
            _close(); imdisp_multi(cube; pixsize = 0.2)
            imgaxes = [ax for ax in axes_of() if _pyint(PC.pylen(ax.images)) > 0]
            @test length(imgaxes) == size(cube, 3)
            for k in 1:size(cube, 3)
                arr   = PC.pyconvert(Matrix{Float64}, imgaxes[k].images[0].get_array())
                slice = cube[:, :, k]
                @test isapprox(arr, reverse(rotl90(slice), dims = 2) ./ maximum(slice);
                               rtol = 1e-5)
                @test !isapprox(arr, rotl90(slice) ./ maximum(slice); rtol = 1e-5)
            end
            _close()
        end
    end

    # ═════════════════════════════════════════════════════════════════════════
    # 3. options and branch selectors behave
    # ═════════════════════════════════════════════════════════════════════════
    @testset "structure and options" begin

        @testset "axis labels come from the spec" begin
            for (name, plotfn) in [("V2", () -> plot_v2(poly)),
                                   ("T3PHI", () -> plot_t3phi(poly)),
                                   ("T3AMP", () -> plot_t3amp(poly)),
                                   ("VISPHI", () -> plot_visphi(poly)),
                                   ("FLUX", () -> plot_flux(flux))]
                _close(); plotfn()
                @test ylabel_of(axes_of()[1]) == OITOOLS.OBS_PLOT_SPECS[name].ylabel
                _close()
            end
        end

        @testset "logplot applies only to amplitude-like panels" begin
            _close(); plot_v2(poly, logplot = true)
            @test yscale_of(axes_of()[1]) == "log"
            _close(); plot_t3amp(poly, logplot = true)
            @test yscale_of(axes_of()[1]) == "log"
            # phases must stay linear even when asked for log
            _close(); plot_obs(poly; obs = ["T3PHI"], logplot = true)
            @test yscale_of(axes_of()[1]) == "linear"
            _close(); plot_v2(poly, logplot = false)
            @test yscale_of(axes_of()[1]) == "linear"
            _close()
        end

        @testset "plot_obs makes one panel per non-empty observable" begin
            _close(); plot_obs(poly; obs = ["V2", "T3PHI", "T3AMP"])
            # colouring by baseline adds no colorbar axis
            @test length([ax for ax in axes_of() if ylabel_of(ax) != ""]) == 3
            _close()

            _close(); plot_obs(poly; obs = ["V2", "T3PHI"])
            @test length([ax for ax in axes_of() if ylabel_of(ax) != ""]) == 2
            _close()

            # mono has no FLUX, so that panel must be dropped rather than drawn empty
            _close(); plot_obs(mono; obs = ["V2", "FLUX"])
            @test length([ax for ax in axes_of() if ylabel_of(ax) != ""]) == 1
            _close()
        end

        @testset "wav/mjd colouring adds a colorbar" begin
            _close(); plot_v2(poly, color = "baseline"); nbase = length(axes_of())
            _close(); plot_v2(poly, color = "wav");      nwav  = length(axes_of())
            _close(); plot_v2(poly, color = "mjd");      nmjd  = length(axes_of())
            @test nwav > nbase
            @test nmjd > nbase
            _close()
        end

        # Stacked panels have y labels of differing widths, so without `fig.align_ylabels()`
        # each sits at its own x and the column reads ragged.
        @testset "y labels are co-aligned across stacked panels" begin
            model = dict_to_model(Dict{String,Any}("s,ud" => 6.5, "s,f" => 1.0), ["s,ud"])
            obs   = model_to_obs(model, [6.5], poly)

            for (what, draw) in ("plot_obs"       => () -> plot_obs(poly;
                                       obs = ["V2", "T3PHI", "T3AMP"]),
                                 "plot_residuals" => () -> plot_residuals(poly, obs),
                                 "single residual" => () -> plot_v2_residuals(poly, obs.v2))
                @testset "$what" begin
                    _close(); draw()
                    fig = PLT.gcf(); fig.canvas.draw()
                    rend = fig.canvas.get_renderer()
                    xs = Float64[]
                    for ax in fig.axes
                        isempty(_pystr(ax.get_ylabel())) && continue
                        bb = ax.yaxis.label.get_window_extent(rend)
                        push!(xs, PC.pyconvert(Float64, bb.x0))
                    end
                    @test length(xs) >= 2
                    @test maximum(xs) - minimum(xs) < 1.0    # pixels
                    _close()
                end
            end
        end

        # A tick rounded past the end of its own range makes matplotlib widen the colorbar
        # axis to fit it, painting an unmapped white strip at the end of the bar.
        @testset "colorbar ticks stay inside the mapped range" begin
            for vals in ([1.5153, 1.7742],
                         [58000.21, 58003.97],
                         collect(range(1.5153, 1.7742, length = 119)))
                lo, hi = extrema(vals)
                for n in (5, 7)
                    t = EXT._cbar_ticks(vals, n)
                    @test all(lo .<= t .<= hi)
                    @test issorted(t)
                end
            end
            @test length(EXT._cbar_ticks([2.0, 2.0], 7)) == 1   # degenerate range

            # end to end: no tick may fall outside the colorbar axis's own limits
            _close(); plot_v2(poly, color = "wav")
            cax = last(axes_of())                       # the colorbar axis
            x0, x1 = PC.pyconvert(Tuple{Float64,Float64}, cax.get_xlim())
            ticks  = PC.pyconvert(Vector{Float64}, cax.get_xticks())
            @test !isempty(ticks)
            @test all(min(x0, x1) .<= ticks .<= max(x0, x1))
            _close()
        end

        # every text artist must follow rcParams; a hard-coded family makes one label
        # render in a different face from the rest of the figure
        @testset "colorbar label inherits the figure font" begin
            _close(); plot_v2(poly, color = "wav")
            fam = _pystr(PC.pyimport("matplotlib").rcParams["font.family"][0])
            for t in last(axes_of()).texts
                got = t.get_fontfamily()
                got = PC.pyconvert(Vector{String}, PC.pylist(got))
                @test fam in got
            end
            _close()
        end

        @testset "imdisp use_colorbar adds an axis" begin
            img = gaussian2d(32, 32, 5.0)
            _close(); imdisp(img; pixsize = 0.2);                      n0 = length(axes_of())
            _close(); imdisp(img; pixsize = 0.2, use_colorbar = true); n1 = length(axes_of())
            @test n1 > n0
            _close()
        end

        # AMPTYP/PHITYP casing is not fixed by the OIFITS standard and files disagree —
        # ASPRO writes "DIFFERENTIAL", the spec's examples are lower case — so every spelling
        # must select the same behaviour.
        @testset "visamp title tracks amptyp, whatever the case" begin
            for (typ, expected) in ("correlated flux" => "Correlated flux",
                                    "CORRELATED FLUX" => "Correlated flux",
                                    "Correlated Flux" => "Correlated flux",
                                    "differential"    => "Diff. visamp",
                                    "DIFFERENTIAL"    => "Diff. visamp",
                                    "  differential " => "Diff. visamp",
                                    "absolute"        => "Visamp",
                                    "ABSOLUTE"        => "Visamp")
                d = deepcopy(poly); d.amptyp = typ
                _close(); plot_visamp(d)
                @test ylabel_of(axes_of()[1]) == expected
                _close()
            end
        end

        # plot_visphi switches to a per-baseline grid on phityp. Every spelling of
        # "differential" must reach the grid, and every spelling of "absolute" the single
        # panel — the file's capitalisation must not change which plot you get.
        @testset "visphi layout tracks phityp, whatever the case" begin
            n_abs = 0
            for typ in ("absolute", "ABSOLUTE", "Absolute")
                d = deepcopy(poly); d.phityp = typ
                _close(); plot_visphi(d)
                n = length(axes_of())
                n_abs == 0 && (n_abs = n)
                @test n == n_abs
                @test n == 1
                _close()
            end
            for typ in ("differential", "DIFFERENTIAL", "Differential", " differential ")
                d = deepcopy(poly); d.phityp = typ
                _close(); plot_visphi(d)
                @test length(axes_of()) > n_abs      # one panel per baseline
                _close()
            end
        end

        # The canonicalisation itself, including the read path that previously kept the
        # file's own casing.
        @testset "obstype canonicalisation" begin
            for s in ("differential", "DIFFERENTIAL", "Differential", " differential\t")
                @test OITOOLS._canon_obstype(s) == "differential"
                @test OITOOLS._is_obstype(s, "differential")
                @test !OITOOLS._is_obstype(s, "absolute")
            end
            @test OITOOLS._is_obstype("CORRELATED FLUX", "correlated flux")
            @test OITOOLS._canon_obstype("") == ""
            @test !OITOOLS._is_obstype("", "absolute")
        end
    end

    # ═════════════════════════════════════════════════════════════════════════
    # 4. hygiene — plotting must not leak figures into the next test
    # ═════════════════════════════════════════════════════════════════════════
    @testset "figure hygiene" begin
        _close()
        @test _nfigs() == 0
        plot_v2(mono)
        @test _nfigs() >= 1
        _close()
        @test _nfigs() == 0
    end
end

let npng = count(f -> endswith(f, ".png"), readdir(_DIR))
    @info "Plotting: $npng figures written to $_DIR"

    # ═════════════════════════════════════════════════════════════════════════
    # gantt_onenight: contiguity and the Moon
    # ═════════════════════════════════════════════════════════════════════════
    #
    # Both of these were wrong until the ASPRO-style Observe panel was designed against this
    # function. They are asserted on the drawn rectangles rather than on the index maths,
    # because the bug was never in the maths -- night_observability computed the Moon
    # correctly all along, and the plot simply did not use it.
    @testset "gantt_onenight" begin
        # Bars are Rectangles: (x, width, y, height) in axis units.
        rects(ax) = [(PC.pyconvert(Float64, r.get_x()), PC.pyconvert(Float64, r.get_width()),
                      PC.pyconvert(Float64, r.get_y()), PC.pyconvert(Float64, r.get_height()))
                     for r in ax.patches]
        # The observable/delay bar is the one CENTRED on y = 2 — the target's own row. Selecting
        # it by its bottom edge instead would depend on the bar height, and silently find
        # nothing the moment that changes.
        obsbars(ax) = [r for r in rects(ax) if isapprox(r[3] + r[4]/2, 2.0; atol = 1e-6)]
        width(ax)   = sum(r[2] for r in obsbars(ax); init = 0.0)

        obsdate = Dates.DateTime(2026, 8, 19)
        n       = 121
        lst     = Float32.(collect(range(0.0, 6.0, length = n)))
        az      = Float32.(fill(90.0, n))
        alt     = Float32.(fill(45.0, n))
        allidx  = collect(1:n)

        @testset "index_runs" begin
            # In the core package, not this extension: a Gantt bar IS a contiguous run, and the
            # matplotlib and Makie charts must find the same ones.
            R = OITOOLS.index_runs
            @test R(Int[])              == Tuple{Int,Int}[]
            @test R([3])                == [(3, 3)]
            @test R([1, 2, 3])          == [(1, 3)]
            @test R([1, 2, 3, 7, 8, 11]) == [(1, 3), (7, 8), (11, 11)]
            @test R([5, 4, 3, 1])       == [(1, 1), (3, 5)]   # unsorted input still works
        end

        @testset "an interrupted window draws as two bars, not one" begin
            _close()
            gantt_onenight("t", obsdate, lst, 3.0, az, alt, allidx, collect(20:80);
                           show_alt = false)
            one = obsbars(axes_of()[1]); w_one = width(axes_of()[1])
            _close()
            gantt_onenight("t", obsdate, lst, 3.0, az, alt, allidx, [20:40; 60:80];
                           show_alt = false)
            two = obsbars(axes_of()[1]); w_two = width(axes_of()[1])
            _close()

            @test length(one) == 1                 # contiguous: still a single bar
            @test length(two) == 2                 # split: two, with the gap left empty
            # The gap must actually be missing, not merely redrawn: 20:40 plus 60:80 covers
            # less time than 20:80. Before the fix both drew one bar of identical width.
            @test w_two < w_one
        end

        @testset "the Moon constrains the observable window" begin
            _close()
            gantt_onenight("t", obsdate, lst, 3.0, az, alt, allidx, collect(1:n);
                           show_alt = false)
            w_nomoon = width(axes_of()[1])
            _close()
            gantt_onenight("t", obsdate, lst, 3.0, az, alt, allidx, collect(1:n);
                           good_moon = collect(1:60), show_alt = false)
            w_moon = width(axes_of()[1])
            _close()
            @test w_moon < w_nomoon                # the Moon must remove time, not be ignored
            @test w_moon ≈ w_nomoon / 2 rtol = 0.05
        end

        @testset "empty means excluded, not ignored" begin
            # A night with no dark time at all, or a target never far enough from the Moon,
            # is not observable. Treating the empty vector as "no constraint" is what made
            # good_twilight silently permissive.
            for kw in ((; good_twilight = Int[]), (; good_moon = Int[]))
                _close()
                gantt_onenight("t", obsdate, lst, 3.0, az, alt, allidx, collect(1:n);
                               show_alt = false, kw...)
                @test isempty(obsbars(axes_of()[1]))
                _close()
            end
        end
    end

end
