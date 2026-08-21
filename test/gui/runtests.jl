# LOAD ORDER IS LOAD-BEARING: matplotlib before the Makie stack, always.
#
# Both stacks ship their own native harfbuzz. Makie pulls Julia's HarfBuzz_jll, and
# conda's `libraqm.so.0` -- which matplotlib imports for text shaping -- is linked against
# conda's. Whichever `libharfbuzz.so.0` is already mapped wins for both, and the JLL build
# does not export `hb_ft_font_get_ft_face`, so importing matplotlib second dies with
#
#     libraqm.so.0: undefined symbol: hb_ft_font_get_ft_face
#
# Measured both ways in one process: Makie first fails, matplotlib first works. Loading
# matplotlib here, at the top, is what makes the rest of this file order-independent -- the
# Makie loads further down are then no-ops as far as the linker is concerned.
#
# This is the same shape of bug as the Qt clash `check_qt_conflict` reports: two builds of
# one native library in a single process. It is a test-only concern, because the application
# never puts matplotlib and Makie together -- OITOOLS keeps matplotlib in an extension, and
# the GUI draws with Makie alone. Only this harness deliberately loads both, to compare them.
ENV["MPLBACKEND"] = "Agg"
using PythonPlot
PythonPlot.matplotlib          # force the import NOW, while conda's harfbuzz can still win

# GLMakie, QMLMakie and QML are what activate OITOOLSGUIExt, where the whole GUI lives.
using OITOOLS, GLMakie, QMLMakie, QML, Test, Dates
const GUI = Base.get_extension(OITOOLS, :OITOOLSGUIExt)
GUI === nothing && error("OITOOLSGUIExt did not load; the GUI tests cannot run")
using .GUI

const OIF = joinpath(@__DIR__, "..", "..", "demos", "data")
const MONO = joinpath(OIF, "AlphaCenA.oifits")
const POLY = joinpath(OIF, "BC2026", "OBJECT1_N.oifits")

@testset "OITOOLSGUI" begin

    @testset "session basics" begin
        s = Session()
        @test isempty(s.datasets) && isempty(s.models) && isempty(s.images)
        @test isempty(s.log)
        # the config catalogue comes from OITOOLS and must be populated
        @test !isempty(s.configs.facilities)
        @test "CHARA" in s.configs.facilities
        @test haskey(s.configs.wave_combiner, "MIRCX_LOWH")
        @test occursin("0 dataset", sprint(show, s))
    end

    @testset "literal rendering round-trips" begin
        L = GUI._literal
        for v in (1, -3, 2.5, 1e-17, 0.1 + 0.2, true, false, nothing, :abc,
                  "plain", "with \"quotes\"", "with\\backslash",
                  [1, 2, 3], [1.5, 2.5], ("a", 2))
            @test eval(Meta.parse(L(v))) == v || (v === nothing && eval(Meta.parse(L(v))) === nothing)
        end
        # Floats must survive exactly — a rounded literal silently changes a fit.
        x = 0.1 + 0.2
        @test eval(Meta.parse(L(x))) === x
        # Dicts render in sorted key order so exported scripts are diffable
        d = Dict("b" => 2.0, "a" => 1.0)
        @test L(d) == "Dict(\"a\" => 1.0, \"b\" => 2.0)"
        @test eval(Meta.parse(L(d))) == d
    end

    @testset "load_dataset! logs a replayable call" begin
        s = Session()
        e = load_dataset!(s, MONO; warn = false, verbose = false)
        @test e.name == "data1"
        @test e.data isa AbstractArray
        @test length(s.log) == 1
        @test s.log[1].binding == "data1"
        @test occursin("readoifits", s.log[1].code)
        # defaults must not be echoed, chosen values must be
        @test occursin("warn = false", s.log[1].code)
        @test !occursin("use_v2", s.log[1].code)     # left at its default

        # a second load gets a distinct binding
        load_dataset!(s, MONO; warn = false, verbose = false)
        @test s.datasets[2].name == "data2"

        @test_throws ArgumentError load_dataset!(s, joinpath(OIF, "does_not_exist.oifits"))
    end

    @testset "exported script is runnable and reproduces the session" begin
        s = Session()
        load_dataset!(s, MONO; warn = false, verbose = false)
        script = export_script(s)
        @test occursin("using OITOOLS", script)

        # Run it in a fresh module — the real test of the reproducibility claim.
        m = Module(:ReplayMono)
        Core.eval(m, :(using OITOOLS))
        Core.eval(m, Meta.parseall(script))
        replayed = Core.eval(m, :data1)

        orig = s.datasets[1].data[1, 1]
        rep  = replayed[1, 1]
        @test rep.nv2 == orig.nv2
        @test rep.nuv == orig.nuv
        @test rep.v2 == orig.v2
        @test rep.uv == orig.uv
    end

    @testset "polychromatic read options survive the round trip" begin
        if isfile(POLY)
            s = Session()
            load_dataset!(s, POLY; polychromatic = true, merge_oi_wavelength = true,
                          warn = false, verbose = false)
            code = s.log[1].code
            @test occursin("polychromatic = true", code)
            @test occursin("merge_oi_wavelength = true", code)

            m = Module(:ReplayPoly)
            Core.eval(m, :(using OITOOLS))
            Core.eval(m, Meta.parseall(export_script(s)))
            rep = Core.eval(m, :data1)
            # the binning is what makes this worth checking: shape must match exactly
            @test size(rep) == size(s.datasets[1].data)
            @test size(rep, 1) > 1                     # really did split by wavelength
        else
            @test_skip "polychromatic test data not present"
        end
    end

    @testset "filter_dataset! derives without mutating the source" begin
        s = Session()
        load_dataset!(s, MONO; warn = false, verbose = false)
        n0 = s.datasets[1].data[1, 1].nv2
        f = filter_dataset!(s, "data1"; filter_bad_data = true, baseline_range = [5e6, 100e6])
        @test f.name == "data2"
        @test f.derived_from == "data1"
        @test s.datasets[1].data[1, 1].nv2 == n0        # source untouched — this is the undo story
        @test f.data[1, 1].nv2 <= n0                    # filtering can only remove

        m = Module(:ReplayFilter)
        Core.eval(m, :(using OITOOLS))
        Core.eval(m, Meta.parseall(export_script(s)))
        @test Core.eval(m, :(data2[1,1].nv2)) == f.data[1, 1].nv2
    end

    @testset "export_script writes a file" begin
        s = Session()
        load_dataset!(s, MONO; warn = false, verbose = false)
        p = tempname() * ".jl"
        @test export_script(s, p) == p
        @test isfile(p) && filesize(p) > 0
        @test occursin("readoifits", read(p, String))
    end
end


@testset "Makie plots" begin
    # CairoMakie so the figures render with no GPU and no display: the GUI uses GLMakie
    # through QMLMakie, but the *content* of a figure is identical either way.
    using .GUI: baseline_names, PlotData

    s = Session()
    load_dataset!(s, MONO; warn = false, verbose = false)
    d = s.datasets[1].data[1, 1]

    @testset "zoom is bounded at both ends" begin
        # The wheel arithmetic, without a mouse. QMLMakie hands Makie the raw Qt delta when a
        # device reports pixelDelta -- ±120 for a discrete wheel under libinput -- and Makie
        # raises (1-speed) to it, so an unbounded axis reached 1e49 in a handful of clicks and
        # then threw `AssertionError: vmin != vmax` out of the tick locator. These bounds are
        # what stops that, so they are worth asserting rather than trusting.
        fig = Makie.Figure(); ax = Makie.Axis(fig[1, 1])
        c = build_canvas(fig, ax)
        update_canvas!(c, d, :uv; color = :baseline)
        hx, hy = c.homespan[]
        @test hx > 0 && hy > 0

        span() = (Float64(ax.finallimits[].widths[1]), Float64(ax.finallimits[].widths[2]))

        # Far past the stop in each direction, then assert where it came to rest. Bounds are
        # asserted to within 10%: the stop cannot be landed on exactly through a Float32 rect
        # under DataAspect, and a tenth of a stop is far below anything a reader would see.
        for _ in 1:200; zoom_step!(c, 1); end
        sx, sy = span()
        # Never crosses the bound, and comes to rest within one detent of it.
        @test ZOOM_MIN_SPAN * hx ≤ sx ≤ ZOOM_MIN_SPAN * hx * ZOOM_PER_DETENT
        @test ZOOM_MIN_SPAN * hy ≤ sy ≤ ZOOM_MIN_SPAN * hy * ZOOM_PER_DETENT
        @test zoom_step!(c, 1) == false        # overzooming does nothing at all
        @test span() == (sx, sy)               # and leaves the view untouched

        for _ in 1:400; zoom_step!(c, -1); end
        sx, sy = span()
        @test ZOOM_MAX_SPAN * hx / ZOOM_PER_DETENT ≤ sx ≤ ZOOM_MAX_SPAN * hx
        @test ZOOM_MAX_SPAN * hy / ZOOM_PER_DETENT ≤ sy ≤ ZOOM_MAX_SPAN * hy
        @test zoom_step!(c, -1) == false
        @test span() == (sx, sy)

        # Every limit stays finite: the assertion in the tick locator fires on the way to Inf.
        @test all(isfinite, (ax.finallimits[].origin..., ax.finallimits[].widths...))

        # Aspect is preserved, which DataAspect makes a correctness claim about the sky.
        zoom_step!(c, 0)                          # no-op
        update_canvas!(c, d, :uv; color = :baseline)
        r0 = (/)(span()...)
        for _ in 1:5; zoom_step!(c, 1); end
        @test (/)(span()...) ≈ r0 rtol = 1e-6
    end

    @testset "baseline names and point info" begin
        n = baseline_names(d)
        @test length(n) == d.nv2
        @test all(occursin("-", x) for x in n)
        i = point_info(d, 1)
        # U+03BC GREEK SMALL LETTER MU, not U+00B5 MICRO SIGN. They look alike in an editor
        # and are not interchangeable here: Makie's default font has a glyph for the Greek
        # letter and none for the micro sign, so a label carrying the wrong one renders as a
        # blank box next to a λ that comes out fine. Comparing by codepoint is the only way
        # this is visible from a test.
        @test occursin("Mλ", i) && occursin("μm", i) && occursin("V²", i)
        @test !occursin("\u00b5", i)
        # the units must be the corrected ones: B/λ, not metres
        b = Float64(d.v2_baseline[1]) / 1e6
        @test occursin(string(round(b, digits = 2)), i)
    end

    @testset "uv_figure" begin
        pd = uv_figure(d)
        @test pd isa PlotData
        # every uv point is drawn (V² baselines AND closure legs), not just the V² subset
        @test length(pd.info) == 2 * d.nuv
        @test d.nuv > d.nv2                      # closure legs really do add coverage
        # conjugate half describes the same measurements
        @test pd.info[1] == pd.info[d.nuv + 1]
        # equal aspect is mandatory or the coverage is misread
        @test pd.axis.aspect[] isa Makie.DataAspect

        p = tempname() * ".png"
        Makie.save(p, pd.figure)
        @test isfile(p) && filesize(p) > 5000

        # colour modes
        for m in (:baseline, :wav, :mjd, :none)
            @test uv_figure(d; color = m) isa PlotData
        end
        @test length(uv_figure(d; conjugate = false).info) == d.nuv
    end

    @testset "observable_figure" begin
        for w in (:v2, :t3phi, :t3amp)
            pd = observable_figure(d, w)
            @test pd isa PlotData
            p = tempname() * ".png"
            Makie.save(p, pd.figure)
            @test filesize(p) > 5000
        end
        @test_throws ArgumentError observable_figure(d, :nonsense)
        # log scale is an axis option, not a reimplementation
        pd = observable_figure(d, :v2; logscale = true)
        @test pd.axis.yscale[] === log10
    end

    @testset "figures are built from the data, not hardcoded" begin
        pd = observable_figure(d, :v2)
        @test length(pd.info) == d.nv2
        @test occursin(basename(d.filename), pd.axis.title[])
    end
end

@testset "ports match their oiplot originals" begin
    # Every assertion is against what oiplot.jl actually drew (test/plotport.jl). No
    # assertion encodes an expectation about the data — that is how the first uv port shipped
    # drawing only V² baselines while uvplot draws the closure legs too.
    ENV["MPLBACKEND"] = "Agg"
    include(joinpath(@__DIR__, "plotport.jl"))

    sess = Session()
    load_dataset!(sess, MONO; warn = false, verbose = false)
    d = sess.datasets[1].data[1, 1]

    # Tick COUNT may differ by a few: matplotlib and Makie use different nice-number
    # algorithms, so identical limits can give 8 vs 10. A large difference still fails.
    TICKSLACK = 3

    @testset "uvplot" begin
        names = GUI.uv_point_labels(d)[1]
        ref = oiplot_ref(() -> uvplot(d))
        got = makie_ref(uv_figure(d); colors = makie_color_map(names),
                        legend = (haslegend = true, nlegend = length(unique(names)),
                                  legendsize = 10.0, legendncol = 5))
        diffs = compare_plots(ref, got; tickslack = TICKSLACK)
        isempty(diffs) || @info "uvplot port differences" diffs
        @test isempty(diffs)
        @test !isempty(got.points) && !isempty(got.colors)
    end

    # Every observable oiplot draws from this dataset, each against its own original.
    # t3phi/t3amp are grouped by TRIPLET and default to the geometric-mean baseline; both
    # were wrong in the first port and are asserted here rather than declared away.
    for (kind, mplfn) in ((:v2,        () -> plot_v2(d)),
                          (:t3phi,     () -> plot_t3phi(d)),
                          (:t3amp,     () -> plot_t3amp(d)),
                          (:t3phi_max, () -> plot_t3phi(d; t3base = "max")),
                          (:t3amp_max, () -> plot_t3amp(d; t3base = "max")))
        @testset "$kind" begin
            spec  = GUI.OBS_SPECS[kind]
            names = GUI.group_names(d, spec)
            ref = oiplot_ref(mplfn)
            got = makie_ref(observable_figure(d, kind); colors = makie_color_map(names),
                            legend = (haslegend = true, nlegend = length(unique(names)),
                                      legendsize = 8.0, legendncol = 4))
            diffs = compare_plots(ref, got; tickslack = TICKSLACK)
            isempty(diffs) || @info "$kind port differences" diffs
            @test isempty(diffs)
            @test !isempty(got.colors)
        end
    end

    # VIS and FLUX need a file that has those tables; AlphaCenA has neither.
    # VIS and FLUX need a file that has those tables; AlphaCenA has neither.
    @testset "visamp and flux" begin
        if isfile(POLY)
            sp = Session()
            load_dataset!(sp, POLY; warn = false, verbose = false)
            dp = sp.datasets[1].data[1, 1]
            for (kind, mplfn, n) in ((:visamp, () -> plot_visamp(dp), dp.nvisamp),
                                     (:flux,   () -> plot_flux(dp),   dp.nflux))
                n > 0 || continue
                spec  = GUI.OBS_SPECS[kind]
                names = GUI.group_names(dp, spec)
                # Each oiplot function has its own default colouring: plot_flux uses "wav",
                # so it draws a colorbar and NO legend. Mirrored by ObsSpec.default_color.
                leg = spec.default_color in (:wav, :mjd) ?
                      (haslegend = false, nlegend = 0, legendsize = NaN, legendncol = 0) :
                      (haslegend = true, nlegend = length(unique(names)),
                       legendsize = 8.0, legendncol = 4)
                ref = oiplot_ref(mplfn)
                got = makie_ref(observable_figure(dp, kind);
                                colors = makie_color_map(names), legend = leg)
                diffs = compare_plots(ref, got; tickslack = TICKSLACK)
                isempty(diffs) || @info "$kind port differences" diffs
                @test isempty(diffs)
            end

            # KNOWN GAP, asserted so it cannot be forgotten: with phityp="differential"
            # oiplot draws visphi as a paginated per-baseline grid against wavelength. That
            # layout is not ported, and drawing the absolute one instead would look right and
            # be a different figure — so the port refuses.
            if occursin("differential", lowercase(String(dp.phityp)))
                @test_throws ArgumentError observable_figure(dp, :visphi)
            end
        else
            @test_skip "polychromatic test data not present"
        end
    end
end


@testset "the shell's drawing path" begin
    # These three bugs all shipped because the shell called plot_into! while the port harness
    # tested uv_figure/observable_figure — two implementations, only one of them checked:
    #   * the uv plot was not isotropic (no DataAspect)
    #   * it had no legend
    #   * choosing a colour mode threw a MethodError and killed the application
    # They are now one implementation; these tests assert that, and exercise every view and
    # colour mode the dropdowns offer.
    using .GUI: draw!, plot_into!

    sess = Session()
    load_dataset!(sess, MONO; warn = false, verbose = false)
    d = sess.datasets[1].data[1, 1]

    @testset "uv keeps equal aspect and a legend" begin
        fig = Makie.Figure(); ax = Makie.Axis(fig[1, 1]); extras = Any[]
        p, info = plot_into!(fig, ax, d, :uv; extras = extras)
        @test ax.aspect[] isa Makie.DataAspect      # isotropic, or the coverage is misread
        @test !isempty(extras)                      # legend present
        @test length(info) == 2 * d.nuv
    end

    @testset "every dropdown combination draws without throwing" begin
        fig = Makie.Figure(); ax = Makie.Axis(fig[1, 1]); extras = Any[]
        for kind in (:uv, :v2, :t3phi, :t3amp), color in (:baseline, :wav, :mjd, :none)
            p, info = plot_into!(fig, ax, d, kind; color = color, extras = extras)
            @test !isempty(info)
            @test length(extras) == 1               # exactly one legend/colorbar, not stacked
        end
    end

    @testset "redrawing does not stack legends" begin
        fig = Makie.Figure(); ax = Makie.Axis(fig[1, 1]); extras = Any[]
        for _ in 1:5
            plot_into!(fig, ax, d, :uv; extras = extras)
        end
        @test length(extras) == 1
    end

    # ═════════════════════════════════════════════════════════════════════════
    # graphics environment (src/graphics.jl)
    # ═════════════════════════════════════════════════════════════════════════
    #
    # configure_graphics! writes to ENV, so every test here saves and restores the variables
    # it touches. The WSL-only branch cannot be exercised on a normal Linux box, so what is
    # tested is the guard logic — each condition that switches the workaround OFF — plus the
    # driver search, which is made reachable by pointing LIBGL_DRIVERS_PATH at a temp dir.
    @testset "graphics environment" begin
        # Core, not the extension: `configure_graphics!` has to be callable before GLMakie is
        # loaded, so it cannot live somewhere GLMakie's arrival creates.
        G = OITOOLS
        touched = ["OITOOLSGUI_NO_GPU_SETUP", "GALLIUM_DRIVER", "LIBGL_ALWAYS_SOFTWARE",
                   "LIBGL_DRIVERS_PATH", "MESA_D3D12_DEFAULT_ADAPTER_NAME"]
        saved = Dict(k => get(ENV, k, nothing) for k in touched)
        restore() = for (k, v) in saved
            v === nothing ? delete!(ENV, k) : (ENV[k] = v)
        end
        for k in touched; delete!(ENV, k); end

        try
            @testset "_isset treats off-values as unset" begin
                for off in ("", "0", "false", "no", "NO", " False ")
                    ENV["OITOOLSGUI_NO_GPU_SETUP"] = off
                    @test !G._isset("OITOOLSGUI_NO_GPU_SETUP")
                end
                for on in ("1", "true", "yes", "d3d12")
                    ENV["OITOOLSGUI_NO_GPU_SETUP"] = on
                    @test G._isset("OITOOLSGUI_NO_GPU_SETUP")
                end
                delete!(ENV, "OITOOLSGUI_NO_GPU_SETUP")
                @test !G._isset("OITOOLSGUI_NO_GPU_SETUP")   # absent is also unset
            end

            @testset "the kill switch wins over everything" begin
                ENV[G.NO_GPU_SETUP] = "1"
                r = configure_graphics!(verbose = false)
                @test r.applied == false
                @test occursin(G.NO_GPU_SETUP, r.reason)
                @test isempty(r.vars)
                @test !haskey(ENV, "GALLIUM_DRIVER")         # and it set nothing
                delete!(ENV, G.NO_GPU_SETUP)
            end

            @testset "an explicit driver choice is left alone" begin
                # Both of these mean "I am choosing the renderer myself", which is exactly
                # what you do when checking whether a bug belongs to the driver.
                ENV["LIBGL_ALWAYS_SOFTWARE"] = "1"
                r = configure_graphics!(verbose = false)
                @test r.applied == false
                @test !haskey(ENV, "GALLIUM_DRIVER")
                delete!(ENV, "LIBGL_ALWAYS_SOFTWARE")

                ENV["GALLIUM_DRIVER"] = "llvmpipe"
                r = configure_graphics!(verbose = false)
                @test r.applied == false
                @test ENV["GALLIUM_DRIVER"] == "llvmpipe"    # unchanged
                delete!(ENV, "GALLIUM_DRIVER")
            end

            @testset "driver search honours LIBGL_DRIVERS_PATH first" begin
                mktempdir() do dir
                    ENV["LIBGL_DRIVERS_PATH"] = dir
                    @test first(G.dri_driver_dirs()) == dir  # searched before the defaults
                    @test "/usr/lib/dri" in G.dri_driver_dirs()

                    # NOT `=== nothing`: a real Mesa install ships /usr/lib/dri/d3d12_dri.so
                    # (this machine has one), and asserting its absence would encode a fact
                    # about the test machine rather than about the search. What must hold is
                    # the ORDER -- nothing is picked up from the empty dir, and once the dir
                    # has a driver it wins over whatever the system provides.
                    @test G.d3d12_driver() != joinpath(dir, "d3d12_dri.so")
                    fake = joinpath(dir, "d3d12_dri.so")
                    write(fake, "")
                    @test G.d3d12_driver() == fake           # searched before the defaults

                    # multiple entries, colon-separated, empties skipped
                    ENV["LIBGL_DRIVERS_PATH"] = "::" * dir
                    @test G.d3d12_driver() == fake
                end
                delete!(ENV, "LIBGL_DRIVERS_PATH")
            end

            @testset "predicates agree with the platform" begin
                @test G.is_wsl() isa Bool
                @test G.has_wsl_nvidia() isa Bool
                Sys.islinux() || @test G.is_wsl() == false   # WSL is a Linux
                # is_wsl is the gate for everything else, so off WSL the answer is fixed
                if !G.is_wsl()
                    r = configure_graphics!(verbose = false)
                    @test r.applied == false
                    @test occursin("not WSL", r.reason)
                end
            end

            @testset "idempotent: calling twice changes nothing" begin
                r1 = configure_graphics!(verbose = false)
                snapshot = Dict(k => get(ENV, k, nothing) for k in touched)
                r2 = configure_graphics!(verbose = false)
                @test r2.applied == false                    # second call never re-applies
                @test all(get(ENV, k, nothing) == snapshot[k] for k in touched)
                # a NamedTuple with the documented shape, whichever branch was taken
                for r in (r1, r2)
                    @test r.reason isa AbstractString && !isempty(r.reason)
                    @test r.vars isa AbstractDict
                end
            end
        finally
            restore()
        end
    end

    # ═════════════════════════════════════════════════════════════════════════
    # UI scale policy (src/scaling.jl)
    # ═════════════════════════════════════════════════════════════════════════
    #
    # Only the policy is testable here, and that is the point of the split: detection happens
    # in QML, from a live Screen property, so this file can cover every branch of "what does
    # the environment variable mean" without a display.
    @testset "UI scale policy" begin
        G = GUI
        saved = get(ENV, G.UI_SCALE_VAR, nothing)
        restore() = saved === nothing ? delete!(ENV, G.UI_SCALE_VAR) : (ENV[G.UI_SCALE_VAR] = saved)
        delete!(ENV, G.UI_SCALE_VAR)

        try
            @testset "unset and the spellings of automatic" begin
                for v in (nothing, "", "  ", "auto", "AUTO", "0")
                    v === nothing ? delete!(ENV, G.UI_SCALE_VAR) : (ENV[G.UI_SCALE_VAR] = v)
                    r = ui_scale_override(verbose = false)
                    @test r.auto
                    @test r.scale == G.UI_SCALE_AUTO == 0.0   # the sentinel QML tests with
                    @test !isempty(r.reason)
                end
                delete!(ENV, G.UI_SCALE_VAR)
            end

            @testset "an explicit scale is taken as given" begin
                for v in ("1", "1.25", "1.5", "2", "0.5", "4")
                    ENV[G.UI_SCALE_VAR] = v
                    r = ui_scale_override(verbose = false)
                    @test !r.auto
                    @test r.scale == parse(Float64, v)
                    @test occursin(v, r.reason) || occursin(string(parse(Float64, v)), r.reason)
                end
            end

            @testset "out of range is clamped, not obeyed" begin
                # 150 is the classic typo: a percentage where a factor was wanted. Obeying it
                # would open a window whose title bar is off-screen, i.e. unrecoverable
                # without editing the environment -- so it is clamped and said out loud.
                for (req, want) in ("150" => G.UI_SCALE_MAX, "8" => G.UI_SCALE_MAX,
                                    "0.1" => G.UI_SCALE_MIN, "-2" => G.UI_SCALE_MIN)
                    ENV[G.UI_SCALE_VAR] = req
                    r = ui_scale_override(verbose = false)
                    @test !r.auto
                    @test r.scale == want
                    @test G.UI_SCALE_MIN <= r.scale <= G.UI_SCALE_MAX
                    @test occursin("clamped", r.reason)
                end
            end

            @testset "nonsense falls back to automatic, never throws" begin
                for v in ("big", "1.2.3", "NaN", "Inf", "1,5", "%150")
                    ENV[G.UI_SCALE_VAR] = v
                    r = ui_scale_override(verbose = false)
                    @test r.auto                      # never a bogus number
                    @test r.scale == G.UI_SCALE_AUTO
                end
            end

            @testset "the result always has the shape QML is handed" begin
                for v in ("", "1.5", "999", "rubbish")
                    ENV[G.UI_SCALE_VAR] = v
                    r = ui_scale_override(verbose = false)
                    @test r.scale isa Float64 && isfinite(r.scale)
                    @test r.auto isa Bool
                    @test r.reason isa AbstractString && !isempty(r.reason)
                    # QML branches on `> 0`, so auto and explicit must be distinguishable
                    @test (r.scale > 0) == !r.auto
                end
            end
        finally
            restore()
        end
    end

    @testset "console pane" begin
        G = GUI
        sh = G.ShellState(Session(), nothing, nothing, nothing, String[], Any[], 0,
                          :uv, :baseline, false, false, "", String[], nothing, "", nothing, nothing, nothing,
                          nothing, nothing, nothing, nothing, Any[], nothing)

        G.console!(sh, "load_dataset!(session, \"a.oifits\")"; kind = :cmd)
        G.console!(sh, "1520 points")
        G.console!(sh, "boom"; kind = :err)
        @test length(sh.console) == 3
        # The markers are what make the pane readable, and what tells a replayable command
        # apart from commentary.
        @test occursin("> load_dataset!", sh.console[1])
        @test occursin("  1520 points",   sh.console[2])
        @test occursin("! boom",          sh.console[3])
        @test all(l -> occursin(r"^\d\d:\d\d:\d\d ", l), sh.console)   # timestamped

        @testset "bounded so a long session cannot grow without bound" begin
            for i in 1:(G.CONSOLE_MAX + 50)
                G.console!(sh, "line $i")
            end
            @test length(sh.console) == G.CONSOLE_MAX
            @test occursin("line $(G.CONSOLE_MAX + 50)", sh.console[end])  # newest kept
            @test !any(l -> occursin("> load_dataset!", l), sh.console)    # oldest dropped
        end
    end

    # ═════════════════════════════════════════════════════════════════════════
    # the live canvas agrees with draw!
    # ═════════════════════════════════════════════════════════════════════════
    #
    # These are two implementations of the same picture: `draw!` builds a fresh figure (and is
    # what the port harness checks against oiplot), while `update_canvas!` pushes into
    # pre-built plots because the live GL context forbids creating any. Two implementations
    # that drift is a mistake this project has already made once -- the shell's first copy of
    # the uv plot lost the equal aspect and the legend, and only the other one was tested. So
    # the agreement is asserted rather than assumed.
    @testset "live canvas matches draw!" begin
        for (kind, color) in ((:uv, :baseline), (:uv, :wav), (:uv, :mjd),
                              (:v2, :baseline), (:v2, :wav), (:t3phi, :baseline))
            @testset "$kind / $color" begin
                pd = GUI.canvas_data(d, kind; color = color)

                f1 = Makie.Figure(); a1 = Makie.Axis(f1[1, 1]); extras = Any[]
                p1, info1 = plot_into!(f1, a1, d, kind; color = color, extras = extras)
                ref = sort([(round(Float64(q[1]), digits = 6), round(Float64(q[2]), digits = 6))
                            for q in p1[1][]])

                got = sort([(round(Float64(x), digits = 6), round(Float64(y), digits = 6))
                            for (x, y) in zip(pd.x, pd.y)])

                @test got == ref                       # same points, same order-insensitive set
                @test length(pd.info) == length(info1)  # same per-point identification
                @test pd.info == info1
                @test length(pd.colors) == length(pd.x) # one colour per point, always resolved
                @test eltype(pd.colors) == Makie.RGBAf
            end
        end

        @testset "axis properties match" begin
            for kind in (:uv, :v2)
                pd = GUI.canvas_data(d, kind; color = :baseline)
                f = Makie.Figure(); a = Makie.Axis(f[1, 1])
                plot_into!(f, a, d, kind; color = :baseline, extras = Any[])
                @test pd.xlabel == a.xlabel[]
                @test pd.ylabel == a.ylabel[]
                # uv coverage must stay isotropic; that is the property whose loss made the
                # first port draw sheared baselines.
                @test pd.isotropic == (a.aspect[] isa Makie.DataAspect)
            end
        end

        @testset "the legend has one entry per group, and none when colour is continuous" begin
            pd_base = GUI.canvas_data(d, :uv; color = :baseline)
            names, _ = uv_point_labels(d)
            @test length(pd_base.legend) == length(unique(names))
            @test first.(pd_base.legend) == sort(unique(names))   # oiplot's ordering
            # Every uv point must be named. uv_point_labels has to read indx_vis as well as
            # indx_v2 and the three indx_t3_* maps: on a dataset with OI_VIS, VIS-only points
            # are reachable from none of the others (2190 of 4932 on v1295Aql) and would keep
            # their initial "", drawn in one colour under a blank legend row. A dataset
            # without OI_VIS cannot catch this, hence both files below.
            @test !any(isempty, first.(pd_base.legend))
            for c in (:wav, :mjd)
                pd = GUI.canvas_data(d, :uv; color = c)
                @test isempty(pd.legend)          # a colorbar replaces it
                @test !isempty(pd.cvals)
                @test !isempty(pd.clabel)
            end
        end

        @testset "every uv point gets a baseline name" begin
            for f in ("2004-data1.oifits", "v1295Aql.oifits")
                dd = readoifits(joinpath(@__DIR__, "data", f);
                                warn = false, verbose = false)[1, 1]
                nm, inf = uv_point_labels(dd)
                @testset "$f" begin
                    @test length(nm) == dd.nuv
                    @test !any(isempty, nm)          # the OI_VIS gap would show up here
                    @test !any(isempty, inf)
                end
            end
        end

        @testset "ramp_colors" begin
            R = GUI.ramp_colors
            @test isempty(R(Float64[]))
            @test length(R([1.0, 2.0, 3.0])) == 3
            @test eltype(R([1.0, 2.0, 3.0])) == Makie.RGBAf
            @test all(c -> c == R([5.0, 5.0])[1], R([5.0, 5.0]))   # degenerate range is flat
            lo, hi = R([0.0, 1.0])[1], R([0.0, 1.0])[2]
            @test lo != hi                                          # and a real range is not
        end
    end

    include(joinpath(@__DIR__, "test_model.jl"))   # Model perspective data layer
    include(joinpath(@__DIR__, "test_imaging.jl")) # Image perspective data layer
    include(joinpath(@__DIR__, "test_observing.jl")) # Observe perspective data layer
    include(joinpath(@__DIR__, "test_gantt.jl"))     # the Gantt port vs its oiplot original

    @testset "the shell path and the figure builders agree" begin
        # Same implementation, so the harness result for uv_figure applies to the shell too.
        fig = Makie.Figure(); ax = Makie.Axis(fig[1, 1]); extras = Any[]
        p_shell, info_shell = plot_into!(fig, ax, d, :v2; extras = extras)
        pd = observable_figure(d, :v2)
        @test info_shell == pd.info
        @test ax.xlabel[] == pd.axis.xlabel[]
        @test ax.ylabel[] == pd.axis.ylabel[]
        shell_pts = sort([(Float64(q[1]), Float64(q[2])) for q in p_shell[1][]])
        fig_pts   = sort([(Float64(q[1]), Float64(q[2])) for q in pd.points[1][]])
        @test shell_pts == fig_pts
    end
end
