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
using OITOOLS, GLMakie, QMLMakie, QML, Test, Dates, TOML
using PairPlots                # the Modeling tab's Corner view (OITOOLSPairPlotsExt)
using Statistics
# Pigeons, for parallel tempering in the Modeling tab (OITOOLSPigeonsExt), when installed
try
    @eval using Pigeons
catch
end
const GUI = Base.get_extension(OITOOLS, :OITOOLSGUIExt)
GUI === nothing && error("OITOOLSGUIExt did not load; the GUI tests cannot run")
using .GUI

const OIF = joinpath(@__DIR__, "..", "..", "demos", "data")
const MONO = joinpath(OIF, "AlphaCenA.oifits")
const POLY = joinpath(OIF, "Beauty_Contest_2026", "OBJECT1_N.oifits")

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

    # Explore's Dataset section (GUI_Explore.md, Phase 2): every row equals the number read
    # directly from the OIdata or the file's tables, for a single-channel file and a
    # polychromatic merged one -- the panel reports the library's numbers, not its own.
    @testset "the dataset summary is the data's own numbers" begin
        G = GUI
        gd = joinpath(@__DIR__, "data")
        for (f, kw) in (("2004-data1.oifits", ()), ("v1295Aql.oifits", ()),
                        ("v1295Aql.oifits", (; polychromatic = true)))
            s = Session()
            e = load_dataset!(s, joinpath(gd, f); warn = false, verbose = false, kw...)
            rows = Dict(G.dataset_summary(e))
            data = vec(e.data)
            hdr = OITOOLS.load_oifits(e.path)
            @test rows["target"] == G._first_few(unique(String[t.target for t in hdr.target]))
            @test rows["instrument"] == G._first_few(unique(String[i.insname for i in hdr.instr]))
            used = unique(vcat((vcat(vec(d.v2_sta_index), vec(d.t3_sta_index), vec(d.vis_sta_index))
                                for d in data)...))
            @test startswith(rows["stations"], string(length(used), ": "))
            mjd = vcat((Float64.(d.uv_mjd) for d in data)...)
            f4(x) = GUI.Printf.@sprintf("%.4f", x)
            @test occursin(f4(minimum(mjd)), rows["MJD"]) && occursin(f4(maximum(mjd)), rows["MJD"])
            lam = vcat((Float64.(d.uv_lam) for d in data)...)
            @test startswith(rows["channels"], string(length(unique(round.(lam; sigdigits = 10)))))
            nv2 = sum(d.nv2 for d in data); nt3 = sum(d.nt3phi for d in data)
            @test occursin("V² $nv2", rows["points"]) && occursin("T3φ $nt3", rows["points"])
            @test haskey(rows, "bins") == (size(e.data) != (1, 1))
        end
        # The 2004 file by name: one target, one channel, six stations, no OI_VIS.
        s = Session()
        e = load_dataset!(s, joinpath(gd, "2004-data1.oifits"); warn = false, verbose = false)
        r = Dict(G.dataset_summary(e))
        @test r["target"] == "FKV0497" && r["instrument"] == "NPOI_2004-04-12"
        @test r["stations"] == "6: E02 AC0 AE0 AW0 W07 AN0"
        @test startswith(r["channels"], "1 at ") && !occursin("VIS", r["points"])
        f = filter_dataset!(s, "data1"; baseline_range = [0.0, 60e6])
        @test Dict(G.dataset_summary(f))["filtered from"] == "data1"
        # a file of two instruments read with one: the summary names that one
        grav = joinpath(pkgdir(OITOOLS), "demos", "data", "o_Leo",
                        "GRAVI.2022-02-28T042658.272_singlescivis_singlesciviscalibrated.fits")
        g = load_dataset!(s, grav; insname = "GRAVITY_SC", warn = false, verbose = false)
        @test Dict(G.dataset_summary(g))["instrument"] == "GRAVITY_SC"
        g2 = load_dataset!(s, grav; warn = false, verbose = false)
        @test Set(split(Dict(G.dataset_summary(g2))["instrument"], ", ")) == Set(["GRAVITY_SC", "GRAVITY_FT"])
    end

    # Explore's Filter section (GUI_Explore.md, Phase 4), through the shell calls the panel
    # makes: each kind of setting applied, the result compared with the library called
    # directly on the parent, the settings kept with the result, the exported script
    # replaying every cut, and the Model panel then evaluating the filtered data.
    @testset "the Filter section: each kind, against the library, and replayed" begin
        G = GUI
        prev = G.SHELL[]
        gd = joinpath(@__DIR__, "data")
        sess = Session()
        load_dataset!(sess, joinpath(gd, "2004-data1.oifits"); warn = false, verbose = false)
        load_dataset!(sess, joinpath(gd, "v1295Aql.oifits"); warn = false, verbose = false)
        load_dataset!(sess, joinpath(OIF, "MWC275_4T.oifits"); warn = false, verbose = false)
        sfig = Makie.Figure(); sax = Makie.Axis(sfig[1, 1])
        sh = G.ShellState(sess, sfig, sax, nothing, String[], Any[], 1,
                          :v2, :baseline, false, :none, :none, false, "", String[],
                          G.build_canvas(sfig, sax), "", nothing, nothing, nothing, nothing,
                          nothing, nothing, nothing, Any[], nothing, nothing, nothing,
                          nothing, "", "", nothing)
        G.SHELL[] = sh
        try
            sel(name) = (sh.current = findfirst(e -> e.name == name, sess.datasets))
            direct(name; kw...) = (d = G.dataset(sess, name).data[1, 1];
                                   OITOOLS.filter_data(d, OITOOLS.set_data_filter(d; kw...)))
            function apply(name, spec)
                sel(name)
                @test !startswith(G.shell_apply_filter(spec), "!")
                e = G.current_dataset(sh)
                @test e.derived_from == name && e.filterspec == spec
                @test G.shell_filter_spec() == spec
                return e.data[1, 1]
            end

            # The choices come from the file: 2004-data1's six stations and fifteen baselines.
            sel("data1")
            ch = split(G.shell_filter_choices(), '\n')
            @test count(startswith("station\t"), ch) == 6
            @test count(startswith("baseline\t"), ch) == 15
            @test any(startswith("night\t"), ch) && any(==("target\tFKV0497"), ch)

            # The preview counts what Apply removes, and draws it on the V² plot.
            sel("data1"); sh.kind = :v2
            pv = G.shell_filter_preview("bl_max\t60")
            want = direct("data1"; baseline_range = [0.0, 60e6])
            d0 = G.dataset(sess, "data1").data[1, 1]
            @test occursin("Filtered $(d0.nv2 - want.nv2)/$(d0.nv2) V²", pv)
            @test length(sh.canvas.prevpoints[]) == d0.nv2 - want.nv2
            @test G.shell_filter_preview_clear() == "" && isempty(sh.canvas.prevpoints[])
            @test startswith(G.shell_filter_preview("bl_min\t9\nbl_max\t1"), "! the baseline range is empty")
            @test startswith(G.shell_apply_filter(""), "! nothing to filter")

            # Ranges: spatial frequency, baseline length, wavelength, nights.
            @test apply("data1", "bl_max\t60").nv2 == want.nv2
            @test apply("data1", "bl_unit\tm\nbl_min\t15").nv2 ==
                  direct("data1"; baseline_length_range = [15.0, 1e99]).nv2
            @test apply("data2", "wav_max\t1.6").nv2 ==
                  direct("data2"; wav_range = [-1.0, 1.6e-6]).nv2
            n1 = minimum(floor.(Int, G.dataset(sess, "data2").data[1, 1].uv_mjd))
            @test apply("data2", "nights\t$n1").nv2 ==
                  direct("data2"; mjd_range = [[Float64(n1), n1 + 1.0]]).nv2
            # The choices offer a night as its MJD window, and a night is a run of epochs, not a
            # UT day: one crossing 0h UT stays one.
            sel("data2")
            w1 = split(first(filter(startswith("night\t"), split(G.shell_filter_choices(), '\n'))), '\t')[2]
            a1, b1 = parse.(Float64, split(w1, ':'))
            @test apply("data2", "nights\t$w1").nv2 == direct("data2"; mjd_range = [[a1, b1]]).nv2
            @test [n.n for n in G._nights([100.9, 101.1, 101.2, 102.95])] == [3, 1]
            @test G._night_window(G._nights([100.9, 101.1, 102.95]), 1) == (100.4, (101.1 + 102.95) / 2)
            @test startswith(G.shell_apply_filter("nights\tyesterday"), "! nights")
            # Telescopes, baselines, triangles.
            @test apply("data1", "exclude_stations\tE02").nv2 ==
                  direct("data1"; exclude_stations = ["E02"]).nv2
            @test apply("data1", "exclude_baselines\tAC0-AE0").nv2 ==
                  direct("data1"; exclude_baselines = ["AC0-AE0"]).nv2
            # Error bars: a cap removes, a floor inflates.
            @test apply("data1", "v2_cap\t0.005").nv2 == direct("data1"; max_v2_err = 0.005).nv2
            f = apply("data1", "v2_add\t0.05")
            @test f.nv2 == d0.nv2 && all(f.v2_err .>= 0.05f0)
            # Quality.
            @test apply("data1", "bad\t1\nsnr\t5").nv2 ==
                  direct("data1"; filter_bad_data = true, filter_v2_snr_threshold = 5.0).nv2
            # Read options: an instrument, and the flags reset (MWC 275 flags 30 finite triangles).
            @test apply("data2", "insname\tMIRCX").nv2 ==
                  readoifits(joinpath(gd, "v1295Aql.oifits"); insname = "MIRCX",
                             warn = false, verbose = false)[1, 1].nv2
            m0 = G.dataset(sess, "data3").data[1, 1]
            @test length(apply("data3", "ignore_flags\t1").t3phi) == length(m0.t3phi) + 30

            # Applying from a filtered dataset filters its unfiltered source afresh.
            g = apply("data1", "bl_max\t60")
            nm = G.current_dataset(sh).name
            @test !startswith(G.shell_apply_filter("bl_min\t60"), "!")
            @test G.current_dataset(sh).derived_from == "data1"
            @test G.current_dataset(sh).data[1, 1].nv2 == direct("data1"; baseline_range = [60e6, 1e99]).nv2
            @test occursin("data1 › filtered", G.shell_dataset_names())

            # The exported script makes every one of these datasets.
            m = Module(:ReplayFilterSection)
            Core.eval(m, :(using OITOOLS))
            Core.eval(m, Meta.parseall(export_script(sess)))
            for e in sess.datasets
                e.derived_from === nothing && continue
                @test Core.eval(m, :($(Symbol(e.name))[1, 1].nv2)) == e.data[1, 1].nv2
                @test Core.eval(m, :(length($(Symbol(e.name))[1, 1].t3phi))) == length(e.data[1, 1].t3phi)
            end

            # The Model panel evaluates the dataset on screen, the filtered one included.
            G.shell_open_model(joinpath(pkgdir(OITOOLS), "demos", "models", "Uniform_Disk.toml"))
            sel("data1"); c_all = G.shell_model_chi2()
            sel(nm);      c_cut = G.shell_model_chi2()
            @test c_all != c_cut && !isempty(c_cut)
        finally
            G.SHELL[] = prev
        end
    end

    # Explore's selection (GUI_Explore.md, Phase 5), without a mouse: the polygon the mouse
    # handler hands over selects known points, "Add to filter" turns them into indices of the
    # unfiltered dataset, and Apply removes exactly those -- also when the selection was made
    # on a filtered copy, and for uv points. The exported script replays it.
    @testset "selection: a known set flagged, on the source and through a filtered copy" begin
        G = GUI
        prev = G.SHELL[]
        sess = Session()
        load_dataset!(sess, joinpath(@__DIR__, "data", "2004-data1.oifits"); warn = false, verbose = false)
        sfig = Makie.Figure(); sax = Makie.Axis(sfig[1, 1])
        sh = G.ShellState(sess, sfig, sax, nothing, String[], Any[], 1,
                          :v2, :baseline, false, :none, :none, false, "", String[],
                          G.build_canvas(sfig, sax), "", nothing, nothing, nothing, nothing,
                          nothing, nothing, nothing, Any[], nothing, nothing, nothing,
                          nothing, "", "", nothing)
        G.SHELL[] = sh
        try
            d = sess.datasets[1].data[1, 1]
            x = Float64.(d.v2_baseline) ./ 1e6; y = Float64.(d.v2)
            box(x0, x1, y0, y1) = [Makie.Point2f(x0, y0), Makie.Point2f(x1, y0),
                                   Makie.Point2f(x1, y1), Makie.Point2f(x0, y1)]
            want = findall(i -> 20 < x[i] < 30 && 0.08 < y[i] < 0.35, eachindex(x))
            @test !isempty(want)

            # The tools switch Makie's rectangle zoom off and back on.
            @test G.shell_set_select_tool("box") == ""
            @test startswith(G.shell_set_select_tool("magic"), "!")
            @test G.select_polygon!(sh, box(20, 30, 0.08, 0.35)) == length(want)
            @test G.SELECTION[].idx == want && length(sh.canvas.selpoints[]) == length(want)
            @test occursin("$(length(want)) V² points selected on data1", G.shell_selection_text())
            spec = G.shell_selection_spec()
            @test spec == "remove_v2\t" * join(want, ",")

            # Through the Filter section: exactly those points go.
            @test !startswith(G.shell_apply_filter(spec), "!")
            f = G.current_dataset(sh)
            @test f.derived_from == "data1" && f.data[1, 1].nv2 == d.nv2 - length(want)
            @test f.data[1, 1].v2 == d.v2[setdiff(1:d.nv2, want)]

            # On the filtered copy, the same box names the same SOURCE points (none left there,
            # so pick a different region and check it maps back by identity).
            fd = f.data[1, 1]
            fx = Float64.(fd.v2_baseline) ./ 1e6; fy = Float64.(fd.v2)
            onf = findall(i -> 45 < fx[i] < 70, eachindex(fx))
            @test G.select_polygon!(sh, box(45, 70, -1, 2)) == length(onf)
            spec2 = G.shell_selection_spec()
            src = parse.(Int, split(split(spec2, '\t')[2], ','))
            @test d.v2[src] == fd.v2[onf] && d.v2_mjd[src] == fd.v2_mjd[onf]

            # uv points: a point or its conjugate inside; removing them removes what uses them.
            sh.current = 1; sh.kind = :uv
            u = Float64.(d.uv[1, :]) ./ 1e6; v = Float64.(d.uv[2, :]) ./ 1e6
            inuv = findall(i -> (40 < u[i] < 80 && 0 < v[i] < 60) || (40 < -u[i] < 80 && 0 < -v[i] < 60), eachindex(u))
            @test G.select_polygon!(sh, box(40, 80, 0, 60)) == length(inuv)
            uspec = G.shell_selection_spec()
            @test startswith(uspec, "remove_uv\t")
            @test !startswith(G.shell_apply_filter(uspec), "!")
            g = G.current_dataset(sh).data[1, 1]
            @test g.nuv == d.nuv - length(inuv) || g.nuv < d.nuv

            @test G.shell_clear_selection() == "" && G.SELECTION[] === nothing
            @test G.shell_set_select_tool("none") == ""

            m = Module(:ReplaySelection)
            Core.eval(m, :(using OITOOLS))
            Core.eval(m, Meta.parseall(export_script(sess)))
            for e in sess.datasets
                e.derived_from === nothing && continue
                @test Core.eval(m, :($(Symbol(e.name))[1, 1].nv2)) == e.data[1, 1].nv2
            end
        finally
            G.SEL_TOOL[] = :none; G.SELECTION[] = nothing
            G.SHELL[] = prev
        end
    end

    # Explore's Export (GUI_Explore.md, Phase 6): a filtered dataset to OIFITS reads back to the
    # points kept, and the source file is never the target.
    @testset "export: filtered dataset to OIFITS" begin
        G = GUI
        prev = G.SHELL[]
        sess = Session()
        load_dataset!(sess, joinpath(@__DIR__, "data", "2004-data1.oifits"); warn = false, verbose = false)
        sfig = Makie.Figure(); sax = Makie.Axis(sfig[1, 1])
        sh = G.ShellState(sess, sfig, sax, nothing, String[], Any[], 1,
                          :v2, :baseline, false, :none, :none, false, "", String[],
                          G.build_canvas(sfig, sax), "", nothing, nothing, nothing, nothing,
                          nothing, nothing, nothing, Any[], nothing, nothing, nothing,
                          nothing, "", "", nothing)
        G.SHELL[] = sh
        dir = mktempdir()
        try
            d = sess.datasets[1].data[1, 1]
            @test startswith(G.shell_export_oifits(sess.datasets[1].path), "!")
            @test !startswith(G.shell_apply_filter("bl_max\t60"), "!")
            f = G.current_dataset(sh).data[1, 1]
            @test 0 < f.nv2 < d.nv2

            out = joinpath(dir, "filtered.oifits")
            @test G.shell_export_oifits("file://" * out) == "saved filtered.oifits"
            back = readoifits(out; warn = false, verbose = false)[1, 1]
            @test sort(back.v2) == sort(f.v2) && sort(back.t3phi) == sort(f.t3phi)
            @test occursin("write_oifits(\"$out\"", G.shell_console())
        finally
            G.SHELL[] = prev
        end
    end

    # Several files as one dataset (GUI_TODO.md, row 4): the open dialog hands over one path per
    # line; the dataset is every point of both files, names them in its summary, re-reads both
    # when a filter needs a read option, refuses both as export targets, and the exported
    # script rebuilds it.
    @testset "several files opened as one dataset" begin
        G = GUI
        prev = G.SHELL[]
        sess = Session()
        sfig = Makie.Figure(); sax = Makie.Axis(sfig[1, 1])
        sh = G.ShellState(sess, sfig, sax, nothing, String[], Any[], 1,
                          :v2, :baseline, false, :none, :none, false, "", String[],
                          G.build_canvas(sfig, sax), "", nothing, nothing, nothing, nothing,
                          nothing, nothing, nothing, Any[], nothing, nothing, nothing,
                          nothing, "", "", nothing)
        G.SHELL[] = sh
        demo = joinpath(pkgdir(OITOOLS), "demos", "data")
        fs = joinpath.(demo, ["AZ_Cyg_2011.oifits", "AZ_Cyg_2014.oifits"])
        try
            each = [readoifits(f; warn = false, verbose = false)[1, 1] for f in fs]
            st = G.shell_open("file://" * fs[1] * "\n" * fs[2] * "\n")
            @test !startswith(st, "could not")
            e = G.current_dataset(sh)
            @test e.files == fs && e.path == fs[1]
            @test e.data[1, 1].nv2 == sum(x.nv2 for x in each)
            @test e.data[1, 1].nt3phi == sum(x.nt3phi for x in each)
            @test G.source_label(e) == "AZ_Cyg_2011.oifits + 1 more"
            @test occursin("AZ_Cyg_2011.oifits + 1 more", G.shell_dataset_names())
            @test occursin("files\t2: ", G.shell_dataset_summary())

            # One file through the same call is the plain single-file dataset.
            G.shell_open(fs[1])
            @test isempty(G.current_dataset(sh).files)
            sh.current = 1

            # A filter with a read option re-reads both files.
            @test !startswith(G.shell_apply_filter("ignore_flags\t1"), "!")
            f = G.current_dataset(sh)
            @test f.files == fs && f.derived_from == e.name
            @test f.data[1, 1].nv2 >= e.data[1, 1].nv2

            # Neither source file is a target, for export or for a simulation.
            @test startswith(G.shell_export_oifits(fs[2]), "!")
            @test occursin("readoifits([", G.shell_script())

            m = Module(:ReplaySeveral)
            Core.eval(m, :(using OITOOLS))
            Core.eval(m, Meta.parseall(export_script(sess)))
            @test Core.eval(m, :($(Symbol(e.name))[1, 1].nv2)) == e.data[1, 1].nv2
            @test Core.eval(m, :($(Symbol(f.name))[1, 1].nv2)) == f.data[1, 1].nv2
        finally
            G.SHELL[] = prev
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

        # A binned dataset is filtered bin by bin and stays binned (it used to come back as its
        # first bin alone, every other channel silently dropped).
        sp = Session()
        load_dataset!(sp, joinpath(@__DIR__, "data", "v1295Aql.oifits"); polychromatic = true,
                      warn = false, verbose = false)
        src = sp.datasets[1].data
        @test length(src) > 1
        fp = filter_dataset!(sp, "data1"; baseline_range = [0.0, 100e6])
        @test size(fp.data) == size(src)
        for k in eachindex(src)
            want = OITOOLS.filter_data(src[k], OITOOLS.set_data_filter(src[k]; baseline_range = [0.0, 100e6]))
            @test fp.data[k].nv2 == want.nv2 && fp.data[k].nt3phi == want.nt3phi
            @test fp.data[k].nv2 < src[k].nv2                     # really cut
        end
        @test fp.readopts == sp.datasets[1].readopts              # read options, not the filter
        mp = Module(:ReplayFilterBinned)
        Core.eval(mp, :(using OITOOLS))
        Core.eval(mp, Meta.parseall(export_script(sp)))
        @test Core.eval(mp, :(size(data2))) == size(src)
        @test Core.eval(mp, :([d.nv2 for d in data2])) == [d.nv2 for d in fp.data]
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
    using .GUI: baseline_names
    # The Makie API is the core package's now: OITOOLSMakieExt fills stubs declared there.
    using OITOOLS: PlotData, imdisp_makie, imdisp_multi_makie, plot_residuals_makie,
                   plot_obs_makie, plot_v2_multifile_makie, plot_facility_makie

    s = Session()
    load_dataset!(s, MONO; warn = false, verbose = false)
    d = s.datasets[1].data[1, 1]

    # The circles behind a uv plot, and a colorbar that stays beside the plot it describes.
    @testset "uv circles, and a colorbar column the plot's own width" begin
        @test GUI.uv_ring_radii(125.0) == [25.0, 50.0, 75.0, 100.0, 125.0]
        @test GUI.uv_ring_radii(60.0)  == [20.0, 40.0, 60.0]
        @test isempty(GUI.uv_ring_radii(0.0)) && isempty(GUI.uv_ring_radii(NaN))
        for rmax in (0.7, 9.0, 31.0, 125.0, 310.0, 4000.0)
            r = GUI.uv_ring_radii(rmax)
            @test 2 <= length(r) <= 5 && r[end] <= rmax
        end

        fig = Makie.Figure(); ax = Makie.Axis(fig[1, 1])
        c = build_canvas(fig, ax)
        update_canvas!(c, d, :uv; color = :baseline)
        rmax = maximum(hypot.(Float64.(d.uv[1, :]), Float64.(d.uv[2, :]))) / 1e6
        @test count(p -> isnan(p[1]), c.uvrings[]) == length(GUI.uv_ring_radii(rmax)) >= 2
        update_canvas!(c, d, :v2)
        @test isempty(c.uvrings[])

        # Baseline colours: no colorbar, so the column is free for the legend below it.
        update_canvas!(c, d, :uv; color = :baseline)
        @test fig.layout.colsizes[1] isa Makie.Auto
        # A colorbar beside a fixed-shape plot: the column becomes the plot's own width.
        update_canvas!(c, d, :uv; color = :mjd)
        @test fig.layout.colsizes[1] isa Makie.Aspect
        GUI.show_image!(c, rand(16, 16), 0.5)
        @test fig.layout.colsizes[1] isa Makie.Aspect && fig.layout.colsizes[1].ratio == 1
        @test isempty(c.uvrings[])
        GUI.hide_image!(c)
        update_canvas!(c, d, :v2)
        @test fig.layout.colsizes[1] isa Makie.Auto
    end

    # Explore §2.6 D3: per baseline on, then off, left the panels' scatter and error bars drawing
    # in the collapsed grid row, as coloured dashes around the legend.
    @testset "turning the panel grid off leaves none of its plots drawing" begin
        G = GUI
        fig = Makie.Figure(); ax = Makie.Axis(fig[1, 1])
        c = build_canvas(fig, ax)
        G.show_panels!(c, true)
        n = G.update_panels!(c, d, :v2)
        @test n > 0 && all(!isempty(c.panels.points[k][]) for k in 1:n)
        G.show_panels!(c, false)
        update_canvas!(c, d, :v2; color = :baseline)
        @test all(isempty(p[]) for p in c.panels.points)
        @test all(isempty(p[]) for p in c.panels.errpoints)
        G.show_panels!(c, true)
        @test G.update_panels!(c, d, :v2) == n && !isempty(c.panels.points[1][])
    end

    # Log axes are labelled at whole decades and simple multiples, never 10^-0.5.
    @testset "log ticks are never fractional powers of ten" begin
        G = GUI
        # At most two significant figures: 0.02, 0.5, 0.31 -- never 10^-0.5 = 0.316…
        isnice(v) = (m = 10 * v / 10.0^floor(log10(v)); isapprox(m, round(m); atol = 1e-6))
        for (lo, hi) in ((0.01, 0.9), (1e-6, 1.0), (0.2, 0.8), (3e-4, 2e-3), (1e-12, 10.0), (0.31, 0.36))
            v = G.log_tick_values(lo, hi)
            @test length(v) >= 2 && all(lo .<= v .<= hi)
            @test all(isnice, v)
            @test all(l -> !occursin('^', l), G.log_tick_label.(v))
        end
        @test G.log_tick_values(0.01, 0.9) == [0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
        @test G.log_tick_values(1e-6, 1.0) == [1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.1, 1.0]
        @test length(G.log_tick_values(1e-12, 10.0)) <= 9
        fig = Makie.Figure(); ax = Makie.Axis(fig[1, 1])
        c = build_canvas(fig, ax)
        update_canvas!(c, d, :v2; logscale = true)
        @test ax.yticks[] === G.log_tick_values
        update_canvas!(c, d, :v2)
        @test ax.ytickformat[] === Makie.automatic
    end

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

        # Right-click reset, without a mouse. Makie's own `:limitreset` returns to
        # `ax.limits[]`, which `zoom_step!` overwrites on every step, so the stock reset came
        # back to the last zoom rather than to the whole picture.
        reset_view!(c)
        @test all(isapprox.(span(), (hx, hy); rtol = 1e-6))
        # and it is idempotent, rather than creeping outwards on each press
        reset_view!(c); @test all(isapprox.(span(), (hx, hy); rtol = 1e-6))

        # An image reverses x to put East on the left. A reset must not quietly clear that:
        # `limits!` takes the reversal as an ORDERING, so handing it an ascending pair mirrors
        # the reconstruction east-for-west.
        show_image!(c, [Float64(i) for i in 1:16, _ in 1:16], 0.25)
        @test ax.xreversed[]
        for _ in 1:5; zoom_step!(c, 1); end
        reset_view!(c)
        @test ax.xreversed[]
        @test isapprox(Float64(ax.finallimits[].widths[1]), 16 * 0.25; rtol = 1e-6)

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

    @testset "uvplot_makie" begin
        pd = uvplot_makie(d)
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
            @test uvplot_makie(d; color = m) isa PlotData
        end
        @test length(uvplot_makie(d; conjugate = false).info) == d.nuv
    end

    @testset "plot_observable_makie" begin
        for w in (:v2, :t3phi, :t3amp)
            pd = plot_observable_makie(d, w)
            @test pd isa PlotData
            p = tempname() * ".png"
            Makie.save(p, pd.figure)
            @test filesize(p) > 5000
        end
        @test_throws ArgumentError plot_observable_makie(d, :nonsense)
        # log scale is an axis option, not a reimplementation
        pd = plot_observable_makie(d, :v2; logscale = true)
        @test pd.axis.yscale[] === log10
    end

    @testset "figures are built from the data, not hardcoded" begin
        pd = plot_observable_makie(d, :v2)
        @test length(pd.info) == d.nv2
        @test occursin(basename(d.filename), pd.axis.title[])
    end

    # ── the ports of oiplot's remaining figures ──────────────────────────────
    #
    # These are the Makie counterparts of `imdisp`, `imdisp_multi` and `plot_residuals`. They
    # exist so that a build without PythonPlot -- a PackageCompiler sysimage, or anyone who
    # simply has not installed matplotlib -- can still draw them.

    @testset "imdisp_makie" begin
        img = [exp(-((i - 32)^2 + (j - 32)^2) / 50) for i in 1:64, j in 1:64]
        pd = imdisp_makie(img; pixsize = 0.25)
        @test pd isa PlotData
        @test pd.axis.aspect[] == 1

        # East LEFT: the x limits run from positive to negative, which is Monnier's
        # orientation and what `imdisp` draws. Getting this backwards mirrors every
        # reconstruction, and nothing else in the figure would look wrong.
        lo, hi = pd.axis.finallimits[].origin[1], pd.axis.finallimits[].widths[1]
        @test occursin("← E", pd.axis.xlabel[])
        @test occursin("→ N", pd.axis.ylabel[])

        # A flat vector of nx² pixels is the other shape `imdisp` takes
        @test imdisp_makie(vec(img); pixsize = 0.25) isa PlotData
        @test_throws ErrorException imdisp_makie(rand(63))     # not square

        # pixsize <= 0 labels in pixels rather than inventing a scale
        @test occursin("pixels", imdisp_makie(img).axis.xlabel[])

        # An all-zero image is drawn rather than divided by roughly zero
        @test imdisp_makie(zeros(8, 8); pixsize = 0.1) isa PlotData
    end

    @testset "imdisp_multi_makie" begin
        img  = [exp(-((i - 16)^2 + (j - 16)^2) / 20) for i in 1:32, j in 1:32]
        cube = cat(img, img .* 0.4, img .* 0.1; dims = 3)
        pd = imdisp_multi_makie(cube; pixsize = 0.25)
        @test length(pd.axis) == 3

        # Per-panel normalisation is `imdisp_multi`'s behaviour: every panel peaks at 1, which
        # is why the panels are NOT comparable in brightness.
        @test imdisp_multi_makie(cube; labels = ["a", "b", "c"]) isa PlotData
        @test_throws ArgumentError imdisp_multi_makie(cube; labels = ["a"])
        # shared_scale is the opt-out, for when the comparison is the point
        @test imdisp_multi_makie(cube; shared_scale = true) isa PlotData
    end

    @testset "plot_residuals_makie" begin
        md = Dict{String,Any}("s,ud" => 6.8, "s,f" => 1.0)
        fm = dict_to_model(md, ["s,ud"])
        pd = plot_residuals_makie(fm, [6.8], d)
        @test pd isa PlotData
        # Two panels per observable present, data over residual
        nobs = (d.nv2 > 0) + (d.nt3amp > 0) + (d.nt3phi > 0) +
               (d.nvisamp > 0) + (d.nvisphi > 0)
        @test length(pd.axis) == 2 * nobs
        @test all(c -> plot_residuals_makie(fm, [6.8], d; color = c) isa PlotData,
                  (:baseline, :wav, :mjd))

        # A phase residual is wrapped before it is divided by σ. Without that a model 359°
        # from the data reads as a huge outlier instead of the 1° agreement it is.
        # `RESIDUAL_KINDS`/`PHASE_KINDS` belong to OITOOLSMakieExt, so they are reached the
        # way this file reaches any extension internal.
        MK = Base.get_extension(OITOOLS, :OITOOLSMakieExt)
        @test OITOOLS.mod360([359.0])[1] ≈ -1.0 atol = 1e-9
        @test :t3phi_max in MK.PHASE_KINDS
        @test :visphi in MK.PHASE_KINDS
        @test !(:v2 in MK.PHASE_KINDS)

        # Closure residuals go against the LONGEST leg, not the geometric mean — the same
        # choice `plot_residuals` makes, and a visibly different x axis if it is missed.
        @test OITOOLS.OBS_SPECS[:t3phi_max].x === :t3_maxbaseline
        @test any(kf -> kf[1] === :t3phi_max, MK.RESIDUAL_KINDS)
    end

    @testset "plot_obs_makie, plot_v2_multifile_makie, plot_facility_makie" begin
        pd = plot_obs_makie(d)
        nobs = (d.nv2 > 0) + (d.nt3amp > 0) + (d.nt3phi > 0) +
               (d.nvisamp > 0) + (d.nvisphi > 0) + (d.nflux > 0)
        @test length(pd.axis) == nobs
        # Only the bottom panel keeps its x label; repeating it per panel is what makes a
        # stacked figure unreadable.
        @test isempty(pd.axis[1].xlabel[])
        @test !isempty(pd.axis[end].xlabel[])
        @test plot_obs_makie(d; kinds = [:v2]) isa PlotData
        @test_throws ArgumentError plot_obs_makie(d; kinds = Symbol[])

        # One legend entry per baseline across every file, not one per file
        m = plot_v2_multifile_makie([d, d])
        @test m isa PlotData
        @test occursin("2 files", m.axis.title[])

        f = plot_facility_makie(read_facility_file("CHARA"))
        @test f isa PlotData
        # Equal aspect: an array on unequal axes misreports which baselines are long.
        @test f.axis.aspect[] isa Makie.DataAspect
    end

    @testset "a log y axis survives non-positive data" begin
        # Real V² contains negative values from noise. matplotlib omits them silently; Makie
        # throws a DomainError out of the RENDER, which surfaces at `save` rather than at the
        # call. Every figure with a log option therefore has to drop those points AND clamp the
        # lower error whisker, which can go negative even for a positive point. Both were
        # missed on the first pass here, and neither shows up until the figure is saved.
        # A fixture that actually contains non-positive V². The one this file loads by default
        # does not, and a log-axis test on all-positive data proves nothing.
        dn = readoifits(joinpath(@__DIR__, "..", "..", "demos", "data", "polaris.oifits");
                        warn = false, verbose = false)[1, 1]
        @test any(v -> v <= 0, dn.v2)
        out = joinpath(mktempdir(), "log.png")
        for p in (plot_observable_makie(dn, :v2; logscale = true),
                  plot_obs_makie(dn; logscale = true),
                  plot_v2_multifile_makie([dn, dn]; logscale = true))
            @test (Makie.save(out, p.figure); filesize(out) > 1000)
        end
        md = Dict{String,Any}("s,ud" => 1.5, "s,f" => 1.0)
        pr = plot_residuals_makie(dict_to_model(md, ["s,ud"]), [1.5], dn; logscale = true)
        @test (Makie.save(out, pr.figure); filesize(out) > 1000)
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
        got = makie_ref(uvplot_makie(d); colors = makie_color_map(names),
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
            got = makie_ref(plot_observable_makie(d, kind); colors = makie_color_map(names),
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
                got = makie_ref(plot_observable_makie(dp, kind);
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
                @test_throws ArgumentError plot_observable_makie(dp, :visphi)
            end
        else
            @test_skip "polychromatic test data not present"
        end
    end
end


@testset "the shell's drawing path" begin
    # These three bugs all shipped because the shell called plot_into! while the port harness
    # tested uvplot_makie/plot_observable_makie — two implementations, only one of them checked:
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

    @testset "panel density from the EDID" begin
        G = GUI
        # a 128-byte block whose first detailed timing is w x h on a wmm x hmm image
        function edid_block(w, h, wmm, hmm)
            e = zeros(UInt8, 128)
            e[55:56] .= (0x02, 0x3a)
            e[57] = w & 0xff; e[59] = (w >> 8) << 4
            e[60] = h & 0xff; e[62] = (h >> 8) << 4
            e[67] = wmm & 0xff; e[68] = hmm & 0xff
            e[69] = ((wmm >> 8) << 4) | (hmm >> 8)
            return e
        end
        @test G.edid_native_mode(edid_block(1920, 1080, 527, 296)) ==
              (; width = 1920, height = 1080, width_mm = 527, height_mm = 296)
        @test G.edid_native_mode(zeros(UInt8, 128)) === nothing      # no timing descriptor
        @test G.edid_native_mode(edid_block(1920, 1080, 0, 0)) === nothing  # no size
        @test G.edid_native_mode(UInt8[]) === nothing

        drm_root = mktempdir()
        function drm_connector(name, status, e)
            dir = mkpath(joinpath(drm_root, name))
            write(joinpath(dir, "status"), status * "\n"); write(joinpath(dir, "edid"), e)
        end
        drm_connector("card1-eDP-1", "connected", edid_block(2560, 1600, 344, 215))
        drm_connector("card1-HDMI-A-1", "connected", edid_block(1920, 1080, 527, 296))
        drm_connector("card1-DP-1", "disconnected", edid_block(2560, 1600, 600, 375))
        drm_connector("card1-DP-2", "connected", UInt8[])
        if Sys.islinux()
            @test G.panel_dpi(2560, 1600; root = drm_root) ≈ 2560 / (344 / 25.4)    # 189
            @test G.panel_dpi(1920, 1080; root = drm_root) ≈ 1920 / (527 / 25.4)    # 92.5
            @test G.panel_dpi(3840, 2160; root = drm_root) == 0.0                   # no such display
            drm_connector("card2-DP-3", "connected", edid_block(1920, 1080, 344, 194))
            @test G.panel_dpi(1920, 1080; root = drm_root) == 0.0                   # two that disagree
        end
        @test G.panel_dpi(1920, 1080; root = joinpath(drm_root, "absent")) == 0.0
        if Sys.islinux()
            @test G.panel_dpi(2561, 1600; root = drm_root) ≈ 2560 / (344 / 25.4)    # a size rounded through 1.5
        end

        # the whole-window scale for a density
        @test OITOOLS.window_scale_for(189.0) == 1.25
        @test OITOOLS.window_scale_for(141.6) == 1.0
        @test OITOOLS.window_scale_for(92.5) == 1.0
        @test OITOOLS.window_scale_for(282.0) == 2.0
    end

    @testset "Qt's scale before it starts" begin
        qt_env = (OITOOLS.QT_SCALE_VARS..., "XDG_CONFIG_HOME", "WSL_DISTRO_NAME", "PATH")
        qt_keep = Dict(k => ENV[k] for k in qt_env if haskey(ENV, k))
        try
            foreach(k -> delete!(ENV, k), OITOOLS.QT_SCALE_VARS)
            qt_cfg = mktempdir(); ENV["XDG_CONFIG_HOME"] = qt_cfg
            # one set by hand is left alone
            ENV["QT_SCALE_FACTOR"] = "2"
            rq = configure_qt_scale!(verbose = false)
            @test !rq.applied && ENV["QT_SCALE_FACTOR"] == "2"
            delete!(ENV, "QT_SCALE_FACTOR")
            # a saved window_scale wins, on any session
            mkpath(joinpath(qt_cfg, "oitools"))
            write(joinpath(qt_cfg, "oitools", "gui.toml"), "ui_scale = 0.7\nwindow_scale = 1.25\n")
            rq = configure_qt_scale!(verbose = false)
            if Sys.islinux()
                @test rq.applied && ENV["QT_SCALE_FACTOR"] == "1.25"
            else
                @test !rq.applied
            end
            # under WSL, from Windows' display scale in the registry: 150% stands for 180 dpi
            if Sys.islinux()
                delete!(ENV, "QT_SCALE_FACTOR")
                rm(joinpath(qt_cfg, "oitools", "gui.toml"))
                ENV["WSL_DISTRO_NAME"] = "test"
                bin = mktempdir()
                write(joinpath(bin, "reg.exe"), "#!/bin/sh\nprintf '\\r\\n    AppliedDPI    REG_DWORD    0x90\\r\\n'\n")
                chmod(joinpath(bin, "reg.exe"), 0o755)
                ENV["PATH"] = bin * ":" * get(qt_keep, "PATH", "/usr/bin")
                rq = configure_qt_scale!(verbose = false)
                @test rq.applied && ENV["QT_SCALE_FACTOR"] == "1.25"
            end
        finally
            foreach(k -> delete!(ENV, k), qt_env)
            merge!(ENV, qt_keep)
        end
    end

    @testset "console pane" begin
        G = GUI
        sh = G.ShellState(Session(), nothing, nothing, nothing, String[], Any[], 0,
                          :uv, :baseline, false, :none, :none, false, "", String[], nothing, "", nothing, nothing, nothing,
                          nothing, nothing, nothing, nothing, Any[], nothing, nothing, nothing,
                          nothing, "", "", nothing)

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
    include(joinpath(@__DIR__, "test_diagnostics.jl")) # residuals and SED panels
    include(joinpath(@__DIR__, "test_imaging.jl")) # Image perspective data layer
    include(joinpath(@__DIR__, "test_observing.jl")) # Observe perspective data layer
    include(joinpath(@__DIR__, "test_tools.jl"))     # Tools perspective data layer
    include(joinpath(@__DIR__, "test_observables.jl")) # observables kept in Explore, fitted in Modeling
    include(joinpath(@__DIR__, "test_batch_a.jl"))   # Filter averaging, merge, cube viewer, fit each
    include(joinpath(@__DIR__, "test_gantt.jl"))     # the Gantt port vs its oiplot original

    # ═════════════════════════════════════════════════════════════════════
    # Save PNG must rebuild the panel that is on SCREEN
    # ═════════════════════════════════════════════════════════════════════
    #
    # The snapshot cannot read the framebuffer (src/gui/snapshot.jl says why), so it re-runs
    # the builders — and every piece of view state the canvas carries has to be handed to them
    # again. Each one that is not is a file that quietly disagrees with the panel it was taken
    # from: residuals and the model/image overplot were both absent, the figure was rebuilt at
    # twice the scene's size (which re-lays out ticks, legend and labels rather than enlarging
    # them), and a zoomed view saved the whole dataset.
    @testset "a saved figure is the picture on the canvas" begin
        G = GUI
        sess = Session()
        load_dataset!(sess, MONO; warn = false, verbose = false)
        sfig = Makie.Figure(); sax = Makie.Axis(sfig[1, 1])
        sh = G.ShellState(sess, sfig, sax, nothing, String[], Any[], 1,
                          :v2, :baseline, false, :none, :none, false, "", String[],
                          G.build_canvas(sfig, sax), "", nothing, nothing, nothing, nothing,
                          nothing, nothing, nothing, Any[], nothing, nothing, nothing,
                          nothing, "", "", nothing)

        # A stand-in reconstruction, so the image residuals and the image overplot both have a
        # producer. What it looks like does not matter; that both paths reach the file does.
        nx, px = 32, 0.5
        img = [exp(-((i - nx / 2)^2 + (j - nx / 2)^2) / (2 * 4.0^2)) for i in 1:nx, j in 1:nx]
        img ./= sum(img)
        # Reshaped: a result holds `(nx, nx, nwav, nepoch)`, one channel for a grey run.
        sh.imaging = G.ImagingResult(reshape(img, nx, nx, 1, 1), 0.0, 0.0, 1, 1.0, 0, 0.0,
                                     G.ImagingSetup(; nx = nx, pixsize = px),
                                     [1.0, 1.0, 1.0], NamedTuple[])

        function snap(sh, which = "explore")
            f = G._snapshot_figure(sh, which, (900, 650), false)
            @test f isa Makie.Figure
            a = f.content[findfirst(b -> b isa Makie.Axis, f.content)]
            return (fig = f, axis = a,
                    scatters = [p for p in a.scene.plots if p isa Makie.Scatter])
        end

        for (kind, res, ov) in ((:v2, :none, :none), (:v2, :image, :none),
                                (:v2, :none, :image), (:t3phi, :image, :none),
                                (:uv, :none, :image))
            @testset "$kind residual=$res overlay=$ov" begin
                sh.kind, sh.residual, sh.overlay = kind, res, ov
                G.refresh_plot!(sh)
                s = snap(sh)
                # The data scatter first, the overplot second: `build_canvas` creates them in
                # that order and the snapshot goes through the same builder.
                @test s.scatters[1][1][] == sh.canvas.points[]
                @test s.scatters[2][1][] == sh.canvas.overpoints[]
                @test s.axis.title[]  == sh.canvas.axis.title[]
                @test s.axis.ylabel[] == sh.canvas.axis.ylabel[]
                @test s.axis.yscale[] === sh.canvas.axis.yscale[]
                # Framed like the scene, not like the size QML asked for.
                @test Makie.size(s.fig.scene) == Makie.size(sh.canvas.figure.scene)
            end
        end

        @testset "a log scale is saved as a log scale" begin
            sh.kind, sh.residual, sh.overlay, sh.logy = :v2, :none, :none, true
            G.refresh_plot!(sh)
            s = snap(sh)
            @test s.axis.yscale[] === log10
            @test s.scatters[1][1][] == sh.canvas.points[]
            sh.logy = false
        end

        @testset "the file follows the zoom" begin
            sh.kind, sh.residual, sh.overlay = :v2, :none, :none
            G.refresh_plot!(sh)
            @test G.zoom_step!(sh.canvas, 3.0)
            fl = sh.canvas.axis.finallimits[]
            lim = snap(sh).axis.limits[]
            @test lim[1] !== nothing && lim[2] !== nothing
            @test minimum(lim[1]) ≈ Float64(fl.origin[1])
            @test maximum(lim[1]) ≈ Float64(fl.origin[1] + fl.widths[1])
            @test minimum(lim[2]) ≈ Float64(fl.origin[2])
            @test maximum(lim[2]) ≈ Float64(fl.origin[2] + fl.widths[2])
        end

        @testset "an image panel keeps its colorbar label and its field of view" begin
            ifig = Makie.Figure(); iax = Makie.Axis(ifig[1, 1]; title = "reconstruction")
            sh.imcanvas = G.build_canvas(ifig, iax)
            G.show_image!(sh.imcanvas, img, px; label = "mean flux / pixel")
            s = snap(sh, "image")
            @test s.axis.title[] == "reconstruction"
            cb = s.fig.content[findfirst(b -> b isa Makie.Colorbar, s.fig.content)]
            @test cb.label[] == "mean flux / pixel"
            # East still on the left, and the same field of view: `show_image!` reverses x.
            @test s.axis.xreversed[] == sh.imcanvas.axis.xreversed[]
            lim = s.axis.limits[]
            @test maximum(lim[1]) ≈ nx * px / 2
            @test minimum(lim[1]) ≈ -nx * px / 2
        end

        @testset "the starting image is drawn as one plane" begin
            # `start_image` hands out a cube and the canvas draws a plane. Without a live
            # canvas `show_image!` is never reached, which is how a cube got to it unnoticed.
            ifig = Makie.Figure(); iax = Makie.Axis(ifig[1, 1])
            sh.imcanvas = G.build_canvas(ifig, iax)
            mode = String(G.ImagingSetup(; nx = nx, pixsize = px).mode)
            prev = G.SHELL[]
            G.SHELL[] = sh
            try
                @test G.shell_show_start_image(nx, px, mode, "gaussian") ==
                      "showing the starting image"
                im = G._canvas_image(sh.imcanvas)
                @test im !== nothing && size(im.image) == (nx, nx)
                @test sh.imcanvas.cbarlabel[] == "starting image"
                # A slider left on a channel a grey start does not have.
                sh.imchannel = 7
                @test G.shell_show_start_image(nx, px, mode, "dirac") ==
                      "showing the starting image"
                @test size(G._canvas_image(sh.imcanvas).image) == (nx, nx)
            finally
                G.SHELL[] = prev
                sh.imchannel = 1
            end
        end

        @testset "right click returns the image to exactly its field of view" begin
            ifig = Makie.Figure(); iax = Makie.Axis(ifig[1, 1])
            sh.imcanvas = G.build_canvas(ifig, iax)
            G.show_image!(sh.imcanvas, img, px)
            view() = (fl = sh.imcanvas.axis.finallimits[];
                      (Float64(fl.widths[1]), Float64(fl.widths[2])))
            field = (nx * px, nx * px)
            rev = sh.imcanvas.axis.xreversed[]
            @test all(isapprox.(view(), field; rtol = 1e-5))
            @test G.zoom_step!(sh.imcanvas, -2.0)
            @test view()[1] > 1.2 * field[1]
            prev = G.SHELL[]
            G.SHELL[] = sh
            try
                G.shell_reset_image_zoom()
                # The field itself, not the field padded by the axis's 5% autolimit margin.
                @test all(isapprox.(view(), field; rtol = 1e-5))
                @test sh.imcanvas.axis.xreversed[] == rev      # East still on the left
            finally
                G.SHELL[] = prev
            end
        end
    end

    @testset "the shell path and the figure builders agree" begin
        # Same implementation, so the harness result for uvplot_makie applies to the shell too.
        fig = Makie.Figure(); ax = Makie.Axis(fig[1, 1]); extras = Any[]
        p_shell, info_shell = plot_into!(fig, ax, d, :v2; extras = extras)
        pd = plot_observable_makie(d, :v2)
        @test info_shell == pd.info
        @test ax.xlabel[] == pd.axis.xlabel[]
        @test ax.ylabel[] == pd.axis.ylabel[]
        shell_pts = sort([(Float64(q[1]), Float64(q[2])) for q in p_shell[1][]])
        fig_pts   = sort([(Float64(q[1]), Float64(q[2])) for q in pd.points[1][]])
        @test shell_pts == fig_pts
    end
    @testset "file:// URLs survive the round trip" begin
        # QML speaks URLs and Julia speaks paths, and every crossing used to hand-roll the
        # conversion — `"file://" * path` out, seven characters off on the way back. That pair
        # is wrong on Windows, where the URL needs three slashes and the drive inside the path,
        # so a round trip is the assertion that matters rather than either half alone.
        for p in ("/home/user/data", "/tmp/a b/c#d.oifits", "/tmp/100%/x?y", "/plain/path")
            @test GUI.strip_file_url(GUI.file_url(p)) == p
        end
        @test GUI.strip_file_url("/not/a/url")   == "/not/a/url"    # a typed path is left alone
        @test GUI.strip_file_url("file:///a/b")  == "/a/b"
        @test startswith(GUI.file_url("/a b"), "file://")
        @test !occursin(" ", GUI.file_url("/a b"))                  # escaped, not raw
    end

    @testset "the appearance config, saved and reset" begin
        # The settings file is per-user and read at startup, so the reset has to remove the
        # FILE: a window restored to its built-in look while the file still says otherwise
        # comes back tweaked at the next launch, which is the confusion the button ends.
        # Point the config at a scratch directory rather than writing the developer's own.
        mktempdir() do dir
            withenv("XDG_CONFIG_HOME" => dir, "APPDATA" => dir) do
                path = GUI.gui_settings_file()
                @test startswith(path, dir)
                @test !isfile(path)
                @test GUI.shell_reset_settings() == ""          # nothing saved yet

                saved = GUI.shell_save_settings(join(["ui_scale\t1.5",
                                                      "controls_style\tFusion",
                                                      "ui_font_pt\t13"], "\n"))
                @test saved == path
                @test isfile(path)
                d = TOML.parsefile(path)
                @test d["ui_scale"] == 1.5
                @test d["controls_style"] == "Fusion"

                @test GUI.shell_reset_settings() == path
                @test !isfile(path)
                @test GUI.shell_reset_settings() == ""          # and again is harmless
            end
        end
    end

    @testset "the controls styles offered" begin
        # Line 1 is the style in force, the rest are the choices. QML slices it exactly there,
        # so an extra header line would silently offer the current style as a choice.
        lines = split(GUI.shell_controls_styles(), '\n')
        @test length(lines) == length(GUI.CONTROLS_STYLES) + 1
        @test lines[2:end] == collect(GUI.CONTROLS_STYLES)
        @test "Basic" in lines[2:end]
        @test !isempty(lines[1])
    end

    @testset "the file picker's places and volumes" begin
        # Every other place is somewhere INSIDE one tree. On Windows that is not enough: `D:`
        # cannot be reached by walking up from `C:`, so without volumes a file on another drive
        # is only openable by typing its path — confirmed missing on Windows before this.
        vols = GUI.picker_volumes()
        @test !isempty(vols)
        @test all(isdir(last(v)) for v in vols)
        @test all(!isempty(first(v)) for v in vols)

        if Sys.iswindows()
            @test any(v -> occursin(r"^[A-Z]:$", first(v)), vols)      # at least one drive
            @test any(v -> uppercase(last(v)) == "C:\\", vols)
        else
            @test any(v -> last(v) == "/", vols)                        # the root is reachable
        end

        rows = [split(r, '\t') for r in split(GUI.picker_places(), '\n') if !isempty(r)]
        @test all(length(r) == 2 for r in rows)
        @test all(isdir(r[2]) for r in rows)
        @test any(r -> r[1] == "Home", rows)
        # the volumes are carried through into the places listing
        for (label, path) in vols
            @test any(r -> r[1] == label && r[2] == path, rows)
        end

        # A path typed into the File field replaces the folder when it is absolute or starts
        # with `~`, and is a name in it otherwise.
        d = mktempdir()
        @test GUI.picker_join(d, "a.py") == joinpath(d, "a.py")
        @test GUI.picker_join(d, "~") == homedir()
        @test GUI.picker_join(d, joinpath("~", "x.py")) == joinpath(homedir(), "x.py")
        @test GUI.picker_join(d, d) == d
    end

    # The companion search, through the shell calls the Model panel makes: AX Cir and its shipped
    # model, the grid preview, the search on a worker, the map in the χ² map chart, Adopt.
    # PMOIRED finds the companion at (+6.23, −28.49) mas; the best node must be within a step.
    @testset "the companion search runs from the Model panel and finds AX Cir's companion" begin
        G = GUI
        prev = G.SHELL[]
        sess = Session()
        load_dataset!(sess, joinpath(pkgdir(OITOOLS), "demos", "data", "AX_Cir.oifits");
                      warn = false, verbose = false)
        sh = G.ShellState(sess, nothing, nothing, nothing, String[], Any[], 1,
                          :uv, :baseline, false, :none, :none, false, "", String[], nothing,
                          "", nothing, nothing, nothing, nothing, nothing, nothing, nothing,
                          Any[], nothing, nothing, nothing, nothing, "", "", nothing)
        cfig = Makie.Figure(); cax = Makie.Axis(cfig[1, 1])
        sh.chi2map = G.build_chi2_map(cfig, cax)
        G.SHELL[] = sh
        try
            G.shell_open_model(joinpath(pkgdir(OITOOLS), "demos", "models", "AX_Cir_Binary.toml"))
            @test G.shell_companion_candidates() == "companion"

            # The header χ²r is the library's, over the ticked observables: V² and T3φ give the
            # fit's 0.7305, and ticking T3amp too changes the number rather than leaving it.
            d1 = sess.datasets[1].data[1, 1]
            fm = parse_model(G._model(sh).dict, G._model(sh).free)
            x1 = Float64[G._model(sh).dict[k] for k in G._model(sh).free]
            w13 = [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
            @test G.shell_model_chi2(true, false, true, false, false, false, false) ==
                  GUI.Printf.@sprintf("%.4f", OITOOLS.chi2_flat(fm, x1, d1; weights = w13) / OITOOLS._ndof(d1, w13))
            @test parse(Float64, G.shell_model_chi2(true, false, true, false, false, false, false)) ≈ 0.7305 atol = 1e-3
            @test G.shell_model_chi2(true, true, true, false, false, false, false) !=
                  G.shell_model_chi2(true, false, true, false, false, false, false)
            @test G.shell_model_chi2(false, false, false, false, false, false, false) == ""

            # Smearing is a menu, not a row: the model file's "auto" shows there, and the menu
            # writes the keys back (off removes them).
            @test G.shell_smearing() == "auto\tincoherent"
            @test !occursin("smear", G.shell_model_rows())
            m = G._model(sh)
            @test G.shell_set_smearing("3", "coherent") == "" && G.shell_smearing() == "3\tcoherent"
            @test m.dict["smear"] == 3 && m.dict["smear_mode"] == "coherent"
            @test G.shell_set_smearing("off", "coherent") == "" && !haskey(m.dict, "smear") &&
                  !haskey(m.dict, "smear_mode") && G.shell_smearing() == "off\tincoherent"
            @test startswith(G.shell_set_smearing("1", "incoherent"), "!")
            @test startswith(G.shell_set_smearing("auto", "sideways"), "!")
            @test G.shell_set_smearing("auto", "incoherent") == "" && m.dict["smear"] == "auto" &&
                  !haskey(m.dict, "smear_mode")
            # The dialog's prefill, to two significant figures: the step rounded down (never
            # coarser than the default), the annulus widened.
            dg = companion_grid(G._model(sh).dict, sess.datasets[1].data[1, 1]; companion = "companion")
            df = parse.(Float64, split(G.shell_companion_defaults("companion"), '\t'))
            @test df[1] <= dg.step
            @test df[2] <= dg.rmin && df[3] >= dg.rmax
            @test all(isapprox.(df, [dg.step, dg.rmin, dg.rmax]; rtol = 0.1))
            @test G.shell_companion_defaults("companion") == "2.0\t2.0\t37"     # AX Cir
            pv = G.shell_companion_grid("companion", NaN, NaN, NaN, NaN)
            @test occursin(r"^\d+ nodes · step ", pv)
            @test startswith(G.shell_companion_grid("companion", NaN, 50.0, 10.0, NaN), "!")

            # The panel's free rows: the companion's position is dropped (the grid sets it), the
            # flux has to be free, and OI_FLUX is refused.
            free = "star,ud\t0\t3\ncompanion,f\t0\t0.5\ncompanion,x\t-60\t60\ncompanion,y\t-60\t60"
            @test startswith(G.shell_companion_search("star,ud\t0\t3", true, false, true,
                                                      false, false, false, false, "companion", NaN, NaN,
                                                      NaN, NaN, NaN, "bradi"), "! free companion,f")
            @test startswith(G.shell_companion_search(free, true, false, true, false, false, true,
                                                      false, "companion", NaN, NaN, NaN, NaN, NaN,
                                                      "vmlmb"), "! the search fits")
            @test startswith(G.shell_companion_search(free, true, false, true, false, false, false, false,
                                                      "companion", NaN, NaN, NaN, NaN, NaN, "lbfgs"),
                             "! unknown per-node fitter")
            # The primary held: only the companion's flux is fitted, by the fast path.
            @test G.shell_companion_search(free, true, false, true, false, false, false, false, "companion",
                                           NaN, NaN, NaN, NaN, NaN, "bradi", false) == "searching…"
            held = ""
            for _ in 1:1200
                r = G.shell_job_poll()
                startswith(r, "running") || (held = r; break)
                sleep(0.25)
            end
            @test startswith(held, "done\tbest node")
            @test last(sh.fits).map.free == ["companion,f"] && last(sh.fits).map.path === :fast

            @test G.shell_companion_search(free, true, false, true, false, false, false, false, "companion",
                                           NaN, NaN, NaN, NaN, NaN, "vmlmb", true) == "searching…"
            line = ""
            for _ in 1:1200
                r = G.shell_job_poll()
                startswith(r, "running") || (line = r; break)
                sleep(0.25)
            end
            @test startswith(line, "done\tbest node")
            f = last(sh.fits)
            s = f.map
            @test s isa CompanionSearch && f.optimiser == "companion search"
            @test s.free == ["star,ud", "companion,f"]
            @test hypot(s.best.x - 6.23, s.best.y + 28.49) <= s.x[2] - s.x[1]   # one grid step
            @test s.nsigma > 5
            # Nine fields per row: `refreshFits()` in ModelTab.qml appends only rows of that length.
            @test all(length(split(l, '\t')) == 9 for l in split(G.shell_fit_rows(), '\n'))
            info = split(G.shell_chi2_map_info(), '\t')
            @test length(info) == 6 && info[6] == "companion"
            @test sh.chi2map.axis.xreversed[] && sh.chi2map.best[] == [Makie.Point2f(s.best.x, s.best.y)]
            @test G.shell_companion_progress() == ""             # cleared once it finished

            # Drawn while it runs: half the nodes fitted, the rest blank, the axes already the
            # grid's, East to the left, the best node so far ringed.
            part = copy(s.chi2r); part[:, 1:end÷2] .= NaN
            g = (; x = s.x, y = s.y, chi2r = part, params = s.params, nnodes = s.nnodes)
            done, total = G.update_companion_progress!(sh.chi2map, g, "companion")
            @test done == count(isfinite, part) && total == s.nnodes && 0 < done < total
            @test sh.chi2map.axis.xreversed[] && count(isnan, sh.chi2map.zs[]) == count(isnan, part)
            @test occursin("of $(s.nnodes) nodes", sh.chi2map.axis.title[])
            empty = (; g..., chi2r = fill(NaN, size(part)))
            @test G.update_companion_progress!(sh.chi2map, empty, "companion") == (0, s.nnodes)
            @test isempty(sh.chi2map.best[])
            G.update_chi2_map!(sh.chi2map, s)

            G.shell_companion_map_view("flux")
            @test sh.chi2map.zs[] ≈ Float32.(companion_flux_map(s)) nans = true

            # Adopt: the best node's fitted values and the position the grid set.
            vals = Dict(split(l, '=')[1] => parse(Float64, split(l, '=')[2])
                        for l in split(G.shell_fit_values(), '\n'))
            @test vals["companion,x"] == s.best.x && vals["companion,y"] == s.best.y
            @test vals["companion,f"] == s.best.f && haskey(vals, "star,ud")

            # Refined: the best minima fitted again with the position free. Adopt takes the
            # refined position, once, and the row carries its σ.
            @test G.shell_companion_search(free, true, false, true, false, false, false, false, "companion",
                                           NaN, NaN, NaN, NaN, NaN, "vmlmb", true, 3) == "searching…"
            line = ""
            for _ in 1:1200
                r = G.shell_job_poll()
                startswith(r, "running") || (line = r; break)
                sleep(0.25)
            end
            @test startswith(line, "done\tbest refined minimum")
            sr = last(sh.fits).map
            @test !isempty(sr.refined) && sr.best.x == sr.refined[1].x
            @test hypot(sr.best.x - 6.228, sr.best.y + 28.487) < 0.2
            fv = split(G.shell_fit_values(), '\n')
            @test count(startswith("companion,x="), fv) == 1 && count(startswith("companion,y="), fv) == 1
            vals = Dict(split(l, '=')[1] => parse(Float64, split(l, '=')[2]) for l in fv)
            @test vals["companion,x"] == sr.best.x && vals["companion,y"] == sr.best.y
            @test occursin("companion,x = ", split(G.shell_fit_rows(), '\n')[1]) &&
                  occursin("±", split(G.shell_fit_rows(), '\n')[1])

            # Multi-start (`multistart_fit`), as the dialog calls it: a grid on the companion's
            # position around PMOIRED's answer, a condition dropping one node, the minima in the
            # map panel, a minimum selected for Adopt.
            mfree = "star,ud\t0\t3\ncompanion,f\t0\t0.5\ncompanion,x\t-60\t60\ncompanion,y\t-60\t60"
            explore = "companion,x\tgrid\t0\t12\t4\ncompanion,y\tgrid\t-32\t-24\t4\nstar,ud\tstart\t\t\t"
            margs(ex, conds) = (mfree, "", ex, conds, true, false, true, false, false, false, false,
                                20, 1, NaN, 100, false)
            @test startswith(G.shell_multistart(margs("star,ud\tstart\t\t\t", "")...), "! explore")
            @test startswith(G.shell_multistart(margs("companion,x\tgrid\t0\t12\t", "")...), "! ")
            @test startswith(G.shell_multistart(margs(explore, "\$companion,x >")...), "! ")
            @test G.shell_multistart(margs(explore, "\$companion,x + \$companion,y < -13")...) == "searching…"
            ml = ""
            for _ in 1:1200
                r = G.shell_job_poll()
                startswith(r, "running") || (ml = r; break)
                sleep(0.25)
            end
            @test startswith(ml, "done\t")
            fm = last(sh.fits)
            ms = fm.map
            @test ms isa MultiStartResult && fm.optimiser == "multi-start"
            @test size(ms.starts, 2) == 11 && ms.nrejected == 1          # (12, −24) breaks the condition
            @test hypot(fm.result.x_opt[3] - 6.228, fm.result.x_opt[4] + 28.487) < 0.2
            @test any(l -> occursin("minima; best χ²r", l), sh.console)
            info = split(G.shell_chi2_map_info(), '\t')
            @test length(info) == 6 && info[6] == "minima" && info[1] == "companion,x"
            mins = split(G.shell_multistart_minima(), '\n')
            @test length(mins) == length(ms.minima) && startswith(mins[1], GUI.Printf.@sprintf("%.4f", fm.result.chi2r))
            @test split(G.shell_multistart_names(), '\n')[1:2] == ["companion,x", "companion,y"]
            @test length(sh.chi2map.mmin[1][]) == length(ms.minima) && !sh.chi2map.axis.xreversed[]
            @test isempty(sh.chi2map.mseg[]) == false && length(sh.chi2map.mseg[]) == 2 * count(==(:minimum), ms.status)
            @test G.shell_multistart_view("star,ud", "companion,f") == "" && sh.chi2map.axis.xlabel[] == "star,ud"
            if length(ms.minima) > 1
                @test G.shell_multistart_select(2) == ""
                @test last(sh.fits).result === ms.minima[2].fit
                vals = Dict(split(l, '=')[1] => parse(Float64, split(l, '=')[2])
                            for l in split(G.shell_fit_values(), '\n'))
                @test vals["companion,x"] == ms.minima[2].fit.x_opt[3]
            end
            @test startswith(G.shell_multistart_select(999), "!")
            # PMOIRED's expl dict as the dialog's rows; a global named as the model import names it
            @test G.shell_multistart_pmoired("{'grid':{'c,x':(-40, 40, 4)}, 'rand':{'MJD0-59800':(-20, 20)}, 'randn':{'c,f':(0.1, 0.01)}}") ==
                  "c,x\tgrid\t-40\t40\t4\nMJD0_59800\tuniform\t-20\t20\nc,f\tnormal\t0.1\t0.01"
            @test startswith(G.shell_multistart_pmoired("not a dict"), "!")
            # the companion's position fixed and gridded: held at each node, the rest fitted
            hfree = "star,ud\t0\t3\ncompanion,f\t0\t0.5"
            hex = "companion,x\tgrid\t4\t8\t2\ncompanion,y\tgrid\t-30\t-26\t2"
            G.shell_set_param("companion,x", "mode", "fixed"); G.shell_set_param("companion,y", "mode", "fixed")
            @test G.shell_multistart(hfree, "", hex, "", true, false, true, false, false, false, false,
                                     20, 1, NaN, 100, false) == "searching…"
            hl = ""
            for _ in 1:1200
                r = G.shell_job_poll()
                startswith(r, "running") || (hl = r; break)
                sleep(0.25)
            end
            @test startswith(hl, "done\t")
            mh = last(sh.fits).map
            @test mh.held == ["companion,x", "companion,y"] && size(mh.starts, 2) == 9
            hv = Dict(split(l, '=')[1] => parse(Float64, split(l, '=')[2])
                      for l in split(G.shell_fit_values(), '\n'))
            @test (hv["companion,x"], hv["companion,y"]) == (6.0, -28.0)      # the node nearest PMOIRED's
            @test occursin("companion,x=", split(G.shell_multistart_minima(), '\n')[1])
            @test "companion,x" in split(G.shell_multistart_names(), '\n')

            # a companion search after it empties the minima
            G.update_chi2_map!(sh.chi2map, sr)
            @test isempty(sh.chi2map.mmin[1][]) && isempty(sh.chi2map.mseg[])
        finally
            G.SHELL[] = prev
        end
    end

    # The Uncertainty selector decides the ± in the fits table: LM's own σ, none, or a
    # bootstrap run after the fit -- whose σ must be the library's for the same call.
    @testset "the Uncertainty selector reaches the fit" begin
        G = GUI
        prev = G.SHELL[]
        sess = Session()
        load_dataset!(sess, joinpath(pkgdir(OITOOLS), "demos", "data", "AX_Cir.oifits");
                      warn = false, verbose = false)
        sh = G.ShellState(sess, nothing, nothing, nothing, String[], Any[], 1,
                          :uv, :baseline, false, :none, :none, false, "", String[], nothing,
                          "", nothing, nothing, nothing, nothing, nothing, nothing, nothing,
                          Any[], nothing, nothing, nothing, nothing, "", "", nothing)
        G.SHELL[] = sh
        function wait_fit()
            for _ in 1:2400
                r = G.shell_job_poll()
                startswith(r, "running") || return r
                sleep(0.25)
            end
            return ""
        end
        free = "star,ud\t0\t3\ncompanion,f\t0\t0.5\ncompanion,x\t-60\t60\ncompanion,y\t-60\t60"
        fit(unc; cons = "", opt = "lsqfit", nboot = 12) =
            G.shell_fit_model(free, cons, "", true, false, true, false, false, false, false, opt, 200,
                              "", "", 60, 0, unc, nboot, "replacement", "config", 7, false)
        try
            G.shell_open_model(joinpath(pkgdir(OITOOLS), "demos", "models", "AX_Cir_Binary.toml"))
            @test fit("jacobian") == "fitting…" && startswith(wait_fit(), "done\tchi2r")
            f = last(sh.fits)
            @test f.sigma == f.result.stderror && f.boot === nothing && f.optimiser == "lsqfit"
            @test occursin("±", split(G.shell_fit_rows(), '\n')[1])
            # The correlations follow the same source: LM's covariance here. PMOIRED: ud–f −65%.
            corr(t) = (l = split(t, '\n'); (split(l[1], '\t'),
                       [parse.(Float64, split(x, '\t')) for x in l[2:end]]))
            nm, R = corr(G.shell_fit_correlation())
            @test nm == f.result.list_free_params
            s = f.result.stderror
            @test all(isapprox.(reduce(vcat, permutedims.(R)), f.result.covar ./ (s * s'); atol = 0.005))
            iud, ifl = findfirst(==("star,ud"), nm), findfirst(==("companion,f"), nm)
            @test R[iud][ifl] ≈ -0.65 atol = 0.03

            @test fit("none") == "fitting…" && startswith(wait_fit(), "done\tchi2r")
            @test last(sh.fits).sigma === nothing
            @test !occursin("±", split(G.shell_fit_rows(), '\n')[1])
            @test G.shell_fit_correlation() == ""
            @test startswith(G.shell_corner_plot(600, 600), "! the latest fit has no samples")

            @test fit("bootstrap") == "fitting…" && startswith(wait_fit(), "done\tchi2r")
            f = last(sh.fits)
            @test f.optimiser == "lsqfit + bootstrap 12" && f.boot isa OITOOLS.BootstrapResult
            mb = copy(G._model(sh).dict)
            for (k, v) in zip(f.result.list_free_params, f.result.x_opt); mb[k] = v; end
            _, lb, ub = G.parse_free_lines(free)
            ref = bootstrap_fit(mb, f.free, sess.datasets[1].data[1, 1]; nboot = 12,
                                mode = :replacement, granularity = :config, fitter = :lsqfit,
                                lb, ub, weights = G.fitting_weights(; t3amp = false),
                                maxeval = 200, seed = 7, threaded = false, verb = false)
            @test f.sigma ≈ ref.sigma rtol = 1e-8
            nb, Rb = corr(G.shell_fit_correlation())
            @test all(isapprox.(reduce(vcat, permutedims.(Rb)), f.boot.correlation; atol = 0.005))

            # The corner plot of the bootstrap's kept replicates, drawn offscreen to a PNG.
            @test occursin("corner=1", G.shell_optional_engines())
            pd = plot_corner_makie(f.boot; truth = f.boot.x_opt)
            @test count(c -> c isa Makie.Axis, pd.figure.content) == 4 * 5 ÷ 2   # lower triangle
            p1 = G.shell_corner_plot(600, 600)
            @test isfile(p1) && filesize(p1) > 10_000
            p2 = G.shell_corner_plot(600, 600)
            @test p2 != p1 && isfile(p2) && !isfile(p1)    # a new name each time, the old removed
            out = joinpath(mktempdir(), "corner")
            @test G.shell_save_corner(out) == "" && filesize(out * ".png") == filesize(p2)
            @test occursin("±", split(G.shell_fit_rows(), '\n')[1])

            # The two samplers beside nested sampling: σ from the posterior, and its corner plot
            sfit(opt, sopts) = G.shell_fit_model(free, "", "", true, false, true, false, false, false,
                                                 false, opt, 200, "", "", 60, 0, "posterior", 12,
                                                 "replacement", "config", 7, false, false, sopts)
            @test sfit("ensemble", "sweeps\t300") == "fitting…" && startswith(wait_fit(), "done\tchi2r")
            f = last(sh.fits)
            @test f.result isa SamplerResult && f.result.sampler === :ensemble && f.optimiser == "ensemble"
            @test f.sigma ≈ vec(std(f.result.posterior; dims = 1))
            @test occursin("fit_model_ensemble(", join(sh.console, "\n"))
            pe = G.shell_corner_plot(600, 600)
            @test isfile(pe) && filesize(pe) > 10_000
            if Base.get_extension(OITOOLS, :OITOOLSPigeonsExt) !== nothing
                @test sfit("pigeons", "rounds\t6\nchains\t6\nexplorer\tslice") == "fitting…" &&
                      startswith(wait_fit(), "done\tchi2r")
                f = last(sh.fits)
                @test f.result isa SamplerResult && f.result.sampler === :pigeons
                @test f.sigma ≈ vec(std(f.result.posterior; dims = 1))
                @test occursin("fit_model_pigeons(", join(sh.console, "\n"))
                pp = G.shell_corner_plot(600, 600)
                @test isfile(pp) && filesize(pp) > 10_000
            end
            @test startswith(fit("bootstrap"; opt = "ensemble"), "! bootstrap refits")

            @test startswith(fit("bootstrap"; opt = "grid"), "! bootstrap refits")
            @test startswith(fit("bootstrap"; cons = "star,ud\t>\t0.1"), "! bootstrap replicates")
            @test startswith(fit("sideways"), "! unknown uncertainty")
        finally
            G.SHELL[] = prev
        end
    end

    # A global belongs to no component, so no component card can remove one. Open the SPARCO
    # model, delete the star and the disc, and `fs` was left behind -- in the fit vector, with
    # its bounds, and unreachable: `shell_remove_component` only deletes keys prefixed
    # `name,`, and a global has no comma. The model could not be cleared to start another.
    @testset "a global survives its components, and can still be removed" begin
        G = GUI
        prev = G.SHELL[]
        sh = G.ShellState(Session(), nothing, nothing, nothing, String[], Any[], 0,
                          :uv, :baseline, false, :none, :none, false, "", String[], nothing,
                          "", nothing, nothing, nothing, nothing, nothing, nothing, nothing,
                          Any[], nothing, nothing, nothing, nothing, "", "", nothing)
        G.SHELL[] = sh
        try
            G.shell_open_model(joinpath(pkgdir(OITOOLS), "demos", "models",
                                        "SPARCO_Star_Disk.toml"))
            m = G._model(sh)
            @test haskey(m.dict, "fs") && "fs" in m.free

            @test G.shell_remove_component("star") == ""
            @test G.shell_remove_component("disk") == ""
            @test collect(keys(m.dict)) == ["fs"]        # the bug: the global outlives both

            # Globals only. A component's keys are what the parser identifies it BY, so a
            # single-key deletion there leaves a component it can no longer place.
            @test startswith(G.shell_remove_param("disk,fwhm"), "!")
            @test startswith(G.shell_remove_param("nope"), "!")

            @test G.shell_remove_param("fs") == ""
            @test isempty(m.dict) && isempty(m.free) && isempty(m.lb) && isempty(m.ub)
            @test G.shell_model_rows() == ""

            # Removed while something still refers to it, the result is reported rather than
            # rewritten -- `model_warnings` names every expression left pointing at nothing.
            m.dict = Dict{String,Any}("star,ud" => 3.0, "star,f" => "\$fs * 1.0",
                                      "disk,fwhm" => 5.0, "disk,f" => "1 - \$fs",
                                      "fs" => 0.4)
            m.free = ["fs"]
            @test G.shell_remove_param("fs") == ""
            w = split(G.shell_model_warnings(), "\n")
            @test length(w) == 2
            @test all(occursin("fs", x) for x in w)
        finally
            G.SHELL[] = prev
        end
    end

    # The expression dialog's checker and keyword list answer the question the parser answers
    # later, so every verdict is held against whether the model actually evaluates.
    @testset "the expression checker agrees with the parser" begin
        G = GUI
        prev = G.SHELL[]
        sh = G.ShellState(Session(), nothing, nothing, nothing, String[], Any[], 0,
                          :uv, :baseline, false, :none, :none, false, "", String[], nothing,
                          "", nothing, nothing, nothing, nothing, nothing, nothing, nothing,
                          Any[], nothing, nothing, nothing, nothing, "", "", nothing)
        G.SHELL[] = sh
        uv = [1e7 2e7; 0.0 0.0]
        runs(d) = try
            eval_model(dict_to_model(Dict{String,Any}(d), String[]), Float64[], uv;
                       wl = fill(1.6e-6, 2), mjd = fill(57000.0, 2))
            true
        catch
            false
        end
        verdicts(key, expr) = Dict(f[1] => f[2] for f in
                                   split.(split(G.shell_check_expression(key, expr), "\n"), '\t'))
        enabled(key) = [f[1] for f in split.(split(G.shell_expression_keywords(key), "\n"), '\t')
                        if f[3] == "1"]
        try
            m = G._model(sh)

            # An ordinary derived value looks a reference up literally: `$ud` in `star,f` is
            # not `star,ud`, and the checker says what to write instead.
            @test G.shell_add_component("star", "ud") == ""
            @test G.shell_set_param("star,f", "expr", "\$ud * 0.1") == ""
            f = split(G.shell_check_expression("star,f", "\$ud * 0.1"), '\t')
            @test f[2] == "unknown" && occursin("\$star,ud", f[3])
            @test !runs(m.dict)
            @test verdicts("star,f", "\$star,ud * 0.1")["star,ud"] == "ok"
            @test G.shell_set_param("star,f", "expr", "\$star,ud * 0.1") == ""
            @test runs(m.dict)
            @test G.shell_remove_component("star") == ""

            # A uniform disc changes law like any other disc: the diameter, its free state and
            # its bounds carry over; the coefficients are seeded, and dropped going back.
            @test G.shell_add_component("star", "ud") == ""
            @test G.shell_ld_law("star") == "ud\tuniform"
            m.dict["star,ud"] = 8.3
            @test G.shell_set_param("star,ud", "mode", "free") == ""
            G.shell_set_param("star,ud", "lb", 0.0); G.shell_set_param("star,ud", "ub", 20.0)
            @test G.shell_set_ld_law("star", "ldclaret4") == ""
            @test m.dict["star,ldclaret4"] == 8.3 && "star,ldclaret4" in m.free
            @test (m.lb["star,ldclaret4"], m.ub["star,ldclaret4"]) == (0.0, 20.0)
            @test all(haskey(m.dict, "star,c$i") for i in 1:4) && !haskey(m.dict, "star,ud")
            @test G.shell_set_ld_law("star", "ud") == ""
            @test m.dict["star,ud"] == 8.3 && "star,ud" in m.free
            @test !any(haskey(m.dict, "star,c$i") for i in 1:4)
            # the add-component dialog offers the same laws, in the same order
            kinds = split.(split(G.shell_component_kinds(), '\n'), '\t')
            @test [k[3] for k in kinds if k[1] == "disc"] == [l.key for l in G.DISC_LAWS]
            @test G.shell_remove_component("star") == ""

            # The visibility function's own seed checks clean, `$BL` included, and its palette
            # offers the baseline beside the wavelength and the epoch.
            @test G.shell_add_component("comp", "vis_func") == ""
            @test all(==("ok"), values(verdicts("comp,visfunc", m.dict["comp,visfunc"])))
            @test runs(m.dict)
            @test sort(enabled("comp,visfunc")) == ["\$BL", "\$MJD", "\$WL"]
            @test verdicts("comp,visfunc", "\$BL * \$WL")["WL"] == "ok"
            @test verdicts("comp,visfunc", "\$BL * \$R")["R"] == "scope"
            @test G.shell_remove_component("comp") == ""

            @test all(length(split(l, '\t')) == 5
                      for l in split(G.shell_expression_keywords("star,f"), "\n"))
            @test sort(enabled("ring,profile")) ==
                  sort(["\$R", "\$MU", "\$D", "\$RMIN", "\$RMAX", "\$DMIN", "\$DMAX"])
            @test sort(enabled("star,f")) == ["\$MJD", "\$WL"]
        finally
            G.SHELL[] = prev
        end
    end

    # The two models the click-through walkthrough could not build -- SPARCO, and anything
    # sharing a parameter -- built here from an empty panel with the panel's own calls.
    @testset "+ parameter and Link build what clicks could not" begin
        G = GUI
        prev = G.SHELL[]
        sh = G.ShellState(Session(), nothing, nothing, nothing, String[], Any[], 0,
                          :uv, :baseline, false, :none, :none, false, "", String[], nothing,
                          "", nothing, nothing, nothing, nothing, nothing, nothing, nothing,
                          Any[], nothing, nothing, nothing, nothing, "", "", nothing)
        G.SHELL[] = sh
        uv = [1e7 2e7; 0.0 0.0]
        runs(d) = try
            eval_model(dict_to_model(Dict{String,Any}(d), String[]), Float64[], uv;
                       wl = fill(1.6e-6, 2), mjd = fill(57000.0, 2))
            true
        catch
            false
        end
        try
            m = G._model(sh)

            # SPARCO: a global flux fraction and a spectral slope on the disc.
            @test G.shell_add_component("star", "point") == ""
            @test G.shell_add_component("disk", "gaussian") == ""
            @test G.shell_add_param("__global__", "fs", 0.4) == ""
            @test G.shell_add_param("disk", "slope", 0.5) == ""
            @test m.dict["fs"] == 0.4 && m.dict["disk,slope"] == 0.5
            @test G.shell_set_param("fs", "mode", "free") == ""
            @test G.shell_set_param("disk,slope", "mode", "free") == ""
            den = "(\$fs * (\$WL/1.6)^(-4) + (1 - \$fs) * (\$WL/1.6)^\$disk,slope)"
            @test G.shell_set_param("star,f", "expr", "\$fs * (\$WL/1.6)^(-4) / " * den) == ""
            @test G.shell_set_param("disk,f", "expr",
                                    "(1 - \$fs) * (\$WL/1.6)^\$disk,slope / " * den) == ""
            demo = read_model_file(joinpath(pkgdir(OITOOLS), "demos", "models",
                                            "SPARCO_Star_Disk.toml"))
            @test sort(collect(keys(m.dict))) == sort(collect(keys(demo.model)))
            @test sort(m.free) == sort(demo.free)
            @test isempty(G.shell_model_warnings())
            @test runs(m.dict)

            # A name the parser already reads, a built-in, a bad identifier, a duplicate, a
            # component that is not there: each refused, and the dialog's check says so first.
            for (c, n) in (("disk", "ud"), ("disk", "u"), ("disk", "az amp1"), ("__global__", "WL"),
                           ("__global__", "BL"), ("__global__", "2x"), ("__global__", "a,b"),
                           ("__global__", "fs"), ("disk", "slope"), ("nope", "a"),
                           ("__global__", "spatial_kernel"), ("__global__", ""))
                @test !isempty(G.shell_param_name_problem(c, n))
                @test startswith(G.shell_add_param(c, n, 1.0), "!")
            end
            @test m.dict["disk,fwhm"] isa Real                 # still a Gaussian

            # The baseline is `$BL`, so `B` is an ordinary name: PMOIRED's own for a spectral
            # slope (`$disk,B`), and a global too.
            @test G.shell_param_name_problem("disk", "B") == ""
            @test G.shell_add_param("disk", "B", 2.0) == "" && m.dict["disk,B"] == 2.0
            @test G.shell_remove_param("disk,B") == ""
            @test G.shell_param_name_problem("__global__", "B") == ""

            # What "+ parameter" adds, Remove takes away; what the parser reads, it does not.
            @test G.shell_param_removable("fs") == "1"
            @test G.shell_param_removable("disk,slope") == "1"
            @test G.shell_param_removable("disk,fwhm") == "0"
            @test G.shell_param_removable("star,f") == "0"
            @test startswith(G.shell_remove_param("disk,fwhm"), "!")
            @test G.shell_add_param("disk", "spare", 2.0) == ""
            @test G.shell_remove_param("disk,spare") == ""
            @test !haskey(m.dict, "disk,spare")

            # What nothing uses, Validation reports once per key, and its "Remove unused"
            # clears in one go.
            @test G.shell_add_param("disk", "spare", 2.0) == ""
            @test G.shell_add_param("disk", "other", 3.0) == ""
            @test occursin("disk,spare", G.shell_model_inspection())
            wl = filter(l -> occursin("disk,spare", l), split(G.shell_model_warnings(), '\n'))
            @test wl == ["disk,spare is not used: no expression refers to \$disk,spare"]
            @test G.shell_remove_unused() == ""
            @test !haskey(m.dict, "disk,spare") && !haskey(m.dict, "disk,other")
            @test haskey(m.dict, "disk,slope")                 # used, so kept
            @test startswith(G.shell_remove_unused(), "!")

            # Link: two orientations through one global, which takes the first's value and is
            # free because one of them was.
            empty!(m.dict); empty!(m.free); empty!(m.lb); empty!(m.ub)
            @test G.shell_add_component("ring", "ring") == ""
            @test G.shell_add_component("disk", "gaussian") == ""
            @test G.shell_set_component_geometry("ring", "orientation", true) == ""
            @test G.shell_set_component_geometry("disk", "orientation", true) == ""
            G.shell_set_param("ring,incl", "value", 30.0); G.shell_set_param("ring,incl", "mode", "free")
            G.shell_set_param("disk,incl", "value", 45.0)
            @test startswith(G.shell_link_params("incl", "ring,incl"), "!")         # one is not a link
            @test G.shell_link_params("incl", "ring,incl\tdisk,incl") == ""
            @test m.dict["incl"] == 30.0 && "incl" in m.free
            @test m.dict["ring,incl"] == "\$incl" && m.dict["disk,incl"] == "\$incl"
            @test !("ring,incl" in m.free) && !haskey(m.lb, "ring,incl")
            @test haskey(m.lb, "incl") && haskey(m.ub, "incl")
            @test runs(m.dict)
            @test startswith(G.shell_link_params("incl", "ring,pa\tdisk,pa"), "!")  # name taken
            @test startswith(G.shell_link_params("pa", "ring,incl\tdisk,pa"), "!")  # an expression

            # Through a global, so removing either side leaves the other intact.
            @test G.shell_remove_component("ring") == ""
            @test !occursin("does not define", G.shell_model_warnings())
            @test runs(m.dict)
        finally
            G.SHELL[] = prev
        end
    end

    @testset "an error in a Qt callback comes back as a string, not through Qt" begin
        # Every callback is registered through `_guarded`; an exception reaching Qt freezes the
        # window, so it has to come back as the panels' "! ..." string and a console line.
        G = GUI
        prev = G.SHELL[]
        sh = G.ShellState(Session(), nothing, nothing, nothing, String[], Any[], 0,
                          :uv, :baseline, false, :none, :none, false, "", String[], nothing,
                          "", nothing, nothing, nothing, nothing, nothing, nothing, nothing,
                          Any[], nothing, nothing, nothing, nothing, "", "", nothing)
        G.SHELL[] = sh
        try
            n0 = length(sh.console)
            msg = redirect_stderr(devnull) do
                G._guarded(x -> error("broken on purpose: $x"), 3)
            end
            @test startswith(msg, "! ") && occursin("broken on purpose: 3", msg)
            @test length(sh.console) == n0 + 1 && occursin("broken on purpose: 3", sh.console[end])
            @test G._guarded(+, 1, 2) == 3
            @test G._guarded(G.shell_console) == G.shell_console()
            # No shell yet, as at start-up: still a string, and nothing thrown.
            G.SHELL[] = nothing
            @test startswith(redirect_stderr(() -> G._guarded(G.shell_console), devnull), "! ")
        finally
            G.SHELL[] = prev
        end
    end

    @testset "every kind the parser reads has a name in the table's header" begin
        # The header reads `disk · Gaussian` from `KIND_LABELS` (sent to the QML at startup), or
        # from the law selector's names for a disc; the model's description reads the same. A
        # kind with neither shows its raw symbol, so a kind added to the parser must be named.
        laws = Set(l.key for l in GUI.DISC_LAWS)
        for k in keys(OITOOLS._KIND_GEOMETRY)
            @test haskey(GUI.KIND_LABELS, string(k)) || string(k) in laws
        end
        @test length(split(GUI.shell_kind_labels(), "\n")) == length(GUI.KIND_LABELS)
        @test GUI._kind_label("ldpow") == "power law disc"
        @test GUI._kind_label("hankel") == "written profile"
    end

    @testset "the model is described by its components, and Save writes that" begin
        # The title is the file; the description is `name: kind` per component, so a law switch
        # or a rename shows at once and a saved file says what it holds.
        G = GUI
        prev = G.SHELL[]
        sh = G.ShellState(Session(), nothing, nothing, nothing, String[], Any[], 0,
                          :uv, :baseline, false, :none, :none, false, "", String[], nothing,
                          "", nothing, nothing, nothing, nothing, nothing, nothing, nothing,
                          Any[], nothing, nothing, nothing, nothing, "", "", nothing)
        G.SHELL[] = sh
        out = tempname() * ".toml"
        try
            f = joinpath(pkgdir(OITOOLS), "demos", "models", "Limb_Darkened_Disk_Claret4.toml")
            @test G.shell_open_model(f) == "Limb_Darkened_Disk_Claret4"
            @test occursin("named in the file \"Claret four-parameter limb-darkened disk\"",
                           sh.console[end])
            @test G.shell_model_description() == "star: four-parameter disc"
            @test G.shell_set_ld_law("star", "ldpow") == ""
            @test G.shell_model_description() == "star: power law disc"
            @test G.shell_rename_component("star", "sun") == ""
            @test G.shell_model_description() == "sun: power law disc"
            @test G.shell_save_model(out) == ""
            @test read_model_file(out).name == "sun: power law disc"
            @test G._model(sh).name == "Limb_Darkened_Disk_Claret4"   # a name of its own is kept
            @test G._model(sh).path == out && !G._model(sh).dirty
        finally
            rm(out; force = true)
            G.SHELL[] = prev
        end
    end

    @testset "several models: open, select, rename, duplicate, remove, fits" begin
        G = GUI
        prev = G.SHELL[]
        sess = Session()
        load_dataset!(sess, joinpath(pkgdir(OITOOLS), "demos", "data", "AlphaCenA.oifits");
                      warn = false, verbose = false)
        sh = G.ShellState(sess, nothing, nothing, nothing, String[], Any[], 1,
                          :uv, :baseline, false, :none, :none, false, "", String[], nothing,
                          "", nothing, nothing, nothing, nothing, nothing, nothing, nothing,
                          Any[], nothing, nothing, nothing, nothing, "", "", nothing)
        G.SHELL[] = sh
        names() = [split(l, '\t') for l in split(G.shell_model_names(), '\n')]
        selected() = only(findall(f -> f[4] == "1", names()))
        dir = joinpath(pkgdir(OITOOLS), "demos", "models")
        out = tempname() * ".toml"
        try
            @test names() == [["model 1", "0", "", "1"]]
            # The empty, unedited first model is filled; an edited one is kept beside the next.
            @test G.shell_open_model(joinpath(dir, "Uniform_Disk.toml")) == "Uniform_Disk"
            @test length(sess.models) == 1
            G.shell_set_model_dirty(true)
            @test G.shell_open_model(joinpath(dir, "Limb_Darkened_Disk_Claret4.toml")) ==
                  "Limb_Darkened_Disk_Claret4"
            @test [f[1] for f in names()] == ["Uniform_Disk", "Limb_Darkened_Disk_Claret4"]
            @test selected() == 2 && names()[1][2] == "1" && names()[2][2] == "0"
            @test names()[2][3] == joinpath(dir, "Limb_Darkened_Disk_Claret4.toml")
            @test G.shell_open_model(joinpath(dir, "Uniform_Disk.toml")) == "Uniform_Disk 2"
            @test G.shell_remove_model(3) == ""

            # The fit history is the window's, each fit naming its model; Adopt takes a fit of
            # the selected model only. The fit panel's rows belong to the model.
            ld = copy(G._model(sh).dict)
            d = sess.datasets[1].data[1, 1]
            r = fit_model_lsqfit(Dict{String,Any}("star,ud" => 6.0, "star,f" => 1.0), ["star,ud"], d)
            push!(sh.fits, (; result = r, optimiser = "lsqfit", seconds = 0.0, free = ["star,ud"],
                              map = nothing, chi2_start = NaN, sigma = nothing, boot = nothing,
                              model = G._model(sh)))
            G.shell_set_model_fit_rows("star,ud\t>\t1\t0.001", "star,ud\t8\t0.5")
            @test startswith(G.shell_fit_values(), "star,ud=")
            @test G.shell_select_model(1) == ""
            @test split(G.shell_fit_rows(), '\t')[9] == "Limb_Darkened_Disk_Claret4"
            @test startswith(G.shell_fit_values(), "! the latest fit is of Limb_Darkened")
            @test G.shell_model_constraints() == "" && G.shell_model_priors() == ""
            @test G._model(sh).dirty && G._model(sh).dict != ld
            @test G.shell_select_model(2) == ""
            @test startswith(G.shell_fit_values(), "star,ud=")
            @test G.shell_model_constraints() == "star,ud\t>\t1\t0.001"
            @test G.shell_model_priors() == "star,ud\t8\t0.5"
            @test G._model(sh).dict == ld

            # Rename: not empty, not another model's name.
            @test startswith(G.shell_rename_model("  "), "!")
            @test startswith(G.shell_rename_model("Uniform_Disk"), "!")
            @test G.shell_rename_model("LD") == "" && names()[2][1] == "LD"
            @test split(G.shell_fit_rows(), '\t')[9] == "LD"                # the fit follows

            # Duplicate: the model and its rows; New: an empty `model k`.
            @test G.shell_duplicate_model() == "LD copy"
            @test selected() == 3 && G._model(sh).dict == ld && G._model(sh).dict !== ld
            @test G.shell_model_constraints() == "star,ud\t>\t1\t0.001"
            @test G.shell_new_model() == "model 4"
            @test isempty(G._model(sh).dict) && selected() == 4

            # A `model k` saved takes the file's name.
            G._model(sh).dict = Dict{String,Any}("s,ud" => 1.0, "s,f" => 1.0)
            @test G.shell_save_model(out) == ""
            @test names()[4][1] == first(splitext(basename(out))) && names()[4][2] == "0"

            # Nothing switches while a fit runs.
            sh.job = :running
            @test startswith(G.shell_select_model(1), "!")
            @test startswith(G.shell_new_model(), "!")
            @test startswith(G.shell_remove_model(1), "!")
            @test startswith(G.shell_open_model(joinpath(dir, "Uniform_Disk.toml")), "!")
            sh.job = nothing
            @test selected() == 4 && length(sess.models) == 4

            # Remove: the selected one goes to its neighbour; an earlier one keeps the selection.
            @test G.shell_remove_model(4) == "" && selected() == 3
            @test G.shell_remove_model(1) == "" && selected() == 2 && names()[2][1] == "LD copy"
            @test G.shell_remove_model(1) == "" && G.shell_remove_model(1) |> startswith("!")
            @test names() == [["LD copy", "1", "", "1"]]
            @test split(G.shell_fit_rows(), '\t')[9] == ""                  # one model: no name
        finally
            rm(out; force = true)
            G.SHELL[] = prev
        end
    end

    @testset "every Error bars column has a tooltip" begin
        # One tip per column of `errCols`, each opening with the column's name, in the same order.
        qml = read(joinpath(pkgdir(OITOOLS), "src", "gui", "qml", "FilterPanel.qml"), String)
        cols = match(r"readonly property var errCols: \[([^\]]*)\]", qml)[1]
        ncols = count(==('"'), cols) ÷ 2
        i = findfirst("readonly property var errTips", qml)
        block = qml[last(i):findnext("\n    ]", qml, last(i))[1]]
        heads = [m[1] for m in eachmatch(r"^\s*\"([a-z ]+):"m, block)]
        @test heads == ["floor", "relative floor", "multiplier", "cap", "relative cap"]
        @test length(heads) == ncols == 5
    end

    @testset "the palette offers the model's parameters; a fit cut short says so" begin
        G = GUI
        prev = G.SHELL[]
        sess = Session()
        load_dataset!(sess, joinpath(pkgdir(OITOOLS), "demos", "data", "AlphaCenA.oifits");
                      warn = false, verbose = false)
        sh = G.ShellState(sess, nothing, nothing, nothing, String[], Any[], 1,
                          :uv, :baseline, false, :none, :none, false, "", String[], nothing,
                          "", nothing, nothing, nothing, nothing, nothing, nothing, nothing,
                          Any[], nothing, nothing, nothing, nothing, "", "", nothing)
        G.SHELL[] = sh
        try
            m = G._model(sh)
            m.dict = Dict{String,Any}("star,ud" => 8.0, "star,f" => 1.0, "disk,fwhm" => 20.0,
                                      "disk,spectrum" => "\$disk,A*(\$WL/1.67)^\$disk,B",
                                      "disk,A" => 1.0, "disk,B" => 2.0, "fs" => 0.5)
            # Written in full, own component first, then globals; never the key itself, and a
            # spectrum is the `f` it sets.
            ps = split(G.shell_expression_params("disk,spectrum"), '\n')
            @test ps[1:3] == ["\$disk,A", "\$disk,B", "\$disk,fwhm"]
            @test "\$fs" in ps && "\$star,ud" in ps
            @test "\$disk,f" ∉ ps && "\$disk,spectrum" ∉ ps
            @test findfirst(==("\$fs"), ps) < findfirst(==("\$star,ud"), ps)
            @test !("\$star,ud" in split(G.shell_expression_params("star,ud"), '\n'))
            m.dict["star,profile"] = "\$nope"; m.dict["star,diamout"] = 8.0
            @test G.shell_expression_params("disk,A") == ""   # does not compile: the checker says why
            empty!(m.dict)

            # One iteration of Levenberg-Marquardt does not converge, and the row says so.
            d = sess.datasets[1].data[1, 1]
            ud = Dict{String,Any}("star,ud" => 6.0, "star,f" => 1.0)
            for it in (1, 200)
                r = fit_model_lsqfit(ud, ["star,ud"], d; maxIter = it)
                push!(sh.fits, (; result = r, optimiser = "lsqfit", seconds = 0.0,
                                  free = ["star,ud"], map = nothing, chi2_start = NaN,
                                  sigma = nothing, boot = nothing))
            end
            # The data want 8.3 mas; capped at 6 the minimum holds the diameter on its bound.
            ud5 = merge(ud, Dict{String,Any}("star,ud" => 5.0))
            r = fit_model_lsqfit(ud5, ["star,ud"], d; ub = Dict("star,ud" => 6.0))
            push!(sh.fits, (; result = r, optimiser = "lsqfit", seconds = 0.0,
                              free = ["star,ud"], map = nothing, chi2_start = NaN,
                              sigma = nothing, boot = nothing))
            rows = split(G.shell_fit_rows(), '\n')
            @test all(length(split(l, '\t')) == 9 for l in rows)
            @test endswith(split(rows[1], '\t')[5], " · on a bound: star,ud")
            @test !occursin("not converged", rows[2])              # newest first: 200 iterations
            @test endswith(split(rows[3], '\t')[5], " · not converged")

            # The uv title names the dataset as its row does, not only its first file.
            pd = G.canvas_data(d, :uv; source = "first.fits + 25 more")
            @test pd.title == "uv coverage — first.fits + 25 more"
            @test G.canvas_data(d, :uv).title == "uv coverage — " * basename(d.filename)
        finally
            G.SHELL[] = prev
        end
    end

    @testset "a profile written with \$RMIN and \$RMAX previews" begin
        # PMOIRED's doughnut is written on the grid's own radii; the preview has to compile it
        # on that grid, as the model does.
        G = GUI
        prev = G.SHELL[]
        prfig = Makie.Figure()
        pplot = G.build_profile_plot(prfig, Makie.Axis(prfig[1, 1]), Makie.Axis(prfig[1, 2]))
        sh = G.ShellState(Session(), nothing, nothing, nothing, String[], Any[], 0,
                          :uv, :baseline, false, :none, :none, false, "", String[], nothing,
                          "", nothing, nothing, nothing, nothing, nothing, nothing, nothing,
                          Any[], nothing, pplot, nothing, nothing, "", "", nothing)
        G.SHELL[] = sh
        py = tempname() * ".py"
        try
            write(py, """{'rim,diamin': 10., 'rim,diamout': 20., 'rim,profile': 'doughnut',
                          'rim,f': 1.0}""")
            @test G.shell_import_pmoired(py) == first(splitext(basename(py)))
            @test occursin("\$RMIN", G._model(sh).dict["rim,profile"])
            g = split(G.shell_profile_grid("rim"), '\t')
            @test parse(Float64, g[2]) == 5.0 && parse(Float64, g[3]) == 10.0
            @test parse(Int, g[6]) > 0                       # the quadrature's node count
            @test G.shell_profile_curves("rim") == ""

            # Pasted, the same dict gives the same model; the console line keeps `$WL` as written.
            filed = copy(G._model(sh).dict)
            @test G.shell_import_pmoired_text(read(py, String)) == "pasted PMOIRED model"
            @test G._model(sh).dict == filed
            @test occursin("pmoired_to_dict(raw\"\"\"", sh.console[end - 1] * sh.console[end])
            @test startswith(G.shell_import_pmoired_text("  \n"), "!")
            @test startswith(G.shell_import_pmoired_text("{'a,ud': }"), "!")
            @test G._model(sh).dict == filed                 # a failed paste changes nothing

            # The template menu's Doughnut is the expression the import writes: 1 mid-grid, 0 at
            # both rims, whatever the grid.
            t = split(G.shell_profile_template("doughnut"), '\t')
            @test t[1] == G._model(sh).dict["rim,profile"]
            @test isempty(t[3])                             # no parameters of its own
            for (lo, hi) in ((5.0, 10.0), (0.0, 2.5))
                fn = OITOOLS.compile_profile(String(t[1]), String[]; rmin = lo, rmax = hi)
                r = [lo, (lo + hi) / 2, hi]
                @test fn(r, zero(r)) ≈ [0.0, 1.0, 0.0] atol = 1e-12
            end
        finally
            rm(py; force = true)
            G.SHELL[] = prev
        end
    end

    @testset "an image component: added from its file, read back for its card" begin
        G = GUI
        prev = G.SHELL[]
        sh = G.ShellState(Session(), nothing, nothing, nothing, String[], Any[], 0,
                          :uv, :baseline, false, :none, :none, false, "", String[], nothing,
                          "", nothing, nothing, nothing, nothing, nothing, nothing, nothing,
                          Any[], nothing, nothing, nothing, nothing, "", "", nothing)
        G.SHELL[] = sh
        d = mktempdir()
        disc(nx, ps, D) = (c = nx ÷ 2 + 1;
                           [hypot((i - c) * ps, (j - c) * ps) <= D / 2 ? 1.0 : 0.0
                            for i in 1:nx, j in 1:nx])
        p1 = joinpath(d, "star.fits"); writefits(disc(64, 0.05, 1.0), p1; pixsize = 0.05)
        p2 = joinpath(d, "wider.fits"); writefits(disc(64, 0.05, 2.0), p2; pixsize = 0.05)
        p3 = joinpath(d, "noscale.fits"); writefits(disc(64, 0.05, 1.0), p3)
        uv = [1e7 2e7; 0.0 0.0]
        try
            m = G._model(sh)
            @test occursin("\tfits_image\t", G.shell_component_kinds())
            @test startswith(G.shell_add_component("img", "fits_image"), "!")  # needs its file

            @test G.shell_add_image_component("img", p1) == ""
            @test m.dict["img,image"] == abspath(p1) && m.dict["img,f"] == 1.0
            @test isempty(m.free)
            @test startswith(G.shell_add_image_component("img", p1), "!")       # name taken
            # The file, scale and normalisation are the card's, not rows of the table.
            @test !occursin("img,image", G.shell_model_rows())
            @test occursin("img\timage\timage", G.shell_model_components())

            f = split(G.shell_image_info("img"), '\t')
            @test f[1] == abspath(p1) && occursin("64×64", f[2]) && occursin("grey", f[2])
            @test parse(Float64, f[3]) ≈ 0.05 && f[4] == "header"
            @test f[5] == "" && f[6] == "plane"

            # Scale: an override, then back to the header's.
            @test G.shell_set_image("img", "pixsize", "0.1") == ""
            @test m.dict["img,pixsize"] == 0.1
            @test split(G.shell_image_info("img"), '\t')[4] == "model"
            @test startswith(G.shell_set_image("img", "pixsize", "-1"), "!")
            @test G.shell_set_image("img", "pixsize", "") == ""
            @test !haskey(m.dict, "img,pixsize")

            # Normalisation.
            @test G.shell_set_image("img", "norm", "relative") == ""
            @test m.dict["img,norm"] == "relative"
            @test startswith(G.shell_set_image("img", "norm", "peak"), "!")
            @test G.shell_set_image("img", "norm", "") == "" && !haskey(m.dict, "img,norm")

            # Another file; one that is not there is refused and changes nothing.
            @test startswith(G.shell_set_image("img", "file", joinpath(d, "nope.fits")), "!")
            @test m.dict["img,image"] == abspath(p1)
            @test G.shell_set_image("img", "file", p2) == ""
            @test m.dict["img,image"] == abspath(p2)

            # Position yes, orientation no: the parser refuses incl/pa on an image.
            @test startswith(G.shell_set_component_geometry("img", "orientation", true), "!")
            @test !haskey(m.dict, "img,incl")
            @test G.shell_set_component_geometry("img", "position", true) == ""

            # A companion beside it shares the flux, and the model evaluates.
            @test G.shell_add_component("pt", "point") == ""
            @test m.dict["img,f"] ≈ 0.5 && m.dict["pt,f"] ≈ 0.5
            @test isempty(G.shell_model_warnings())
            V = eval_model(dict_to_model(Dict{String,Any}(m.dict), String[]), Float64[], uv)
            @test abs(V[1]) < 1.0 && all(isfinite, V)

            # A file without a pixel scale is added with a placeholder the card flags.
            @test G.shell_add_image_component("bare", p3) == ""
            @test m.dict["bare,pixsize"] == 1.0
            @test split(G.shell_image_info("bare"), '\t')[4] == "none"
            @test startswith(G.shell_set_image("bare", "pixsize", ""), "!")

            # A file the reader refuses never enters the model.
            cube = joinpath(d, "cube3.fits")
            OITOOLS.AstroFITS.writefits!(cube, Pair{String,Any}[], rand(16, 16, 3))
            n = length(m.dict)
            @test startswith(G.shell_add_image_component("bad", cube), "!")
            @test length(m.dict) == n

            # Rendered from the model Julia holds: the image's keys are not table rows, and a
            # model rebuilt from the rows drew this disc as a point.
            only_img = Dict{String,Any}(k => v for (k, v) in m.dict if startswith(k, "img,"))
            r = G.render_model_image(only_img, String[], Float64[]; nx = 64, pixsize = 0.05)
            @test count(>(0.01 * maximum(r.image)), r.image) > 100
            @test !startswith(G.shell_model_image(64, 0.05, 0, 0), "!")

            # A cube depends on wavelength through its planes, with no $WL anywhere.
            pc = joinpath(d, "cube.fits")
            writefits(cat(disc(64, 0.05, 1.0), disc(64, 0.05, 2.0); dims = 3), pc;
                      pixsize = 0.05, wavelengths = [1.6e-6, 2.2e-6])
            @test G.model_depends_on(Dict{String,Any}("c,image" => pc, "c,f" => 1.0)).wl
            @test !G.model_depends_on(Dict{String,Any}("c,image" => p1, "c,f" => 1.0)).wl
            @test occursin("2 channels 1.6–2.2 µm", G.image_facts(pc).summary)

            # Image › Save FITS… writes what this reads: a reconstruction saved by the Image
            # tab's own call and loaded as a component has the reconstruction's visibilities
            # at the data's uv points — same centre, same orientation, same scale.
            data = readoifits(joinpath(@__DIR__, "data", "2004-data1.oifits");
                              warn = false, verbose = false, T = Float64)
            rs = G.ImagingSetup(; nx = 32, pixsize = auto_pixsize(data[1, 1]),
                                startkind = :gaussian)
            sh.imaging = G.reconstruct_image(data, rs; maxiter = 30)
            saved = joinpath(d, "reconstruction.fits")
            @test startswith(G.shell_save_image(saved), "saved")
            @test G.shell_add_image_component("rec", saved) == ""
            plane = G.result_plane(sh.imaging)
            ft = setup_ft(data, rs.nx, rs.pixsize)
            Vimg = image_to_vis(plane, ft[1, 1].uv)
            Vcmp = eval_model(dict_to_model(Dict{String,Any}("rec,image" => m.dict["rec,image"]),
                                            String[]), Float64[], data[1, 1].uv;
                              wl = data[1, 1].uv_lam)
            @test maximum(abs.(Vcmp .- Vimg)) < 1e-5
            G.shell_remove_component("rec")

            # Saved and read back, the model still names the file.
            toml = joinpath(d, "m.toml")
            @test G.shell_save_model(toml) == ""
            r = read_model_file(toml)
            @test r.model["img,image"] == abspath(p2)
            @test dict_to_model(r.model, String[]) isa FlatModel
        finally
            G.SHELL[] = prev
        end
    end
end

# The QML itself, last: it opens a real window in a subprocess, and nothing may create a
# Makie plot afterwards (see test_qml.jl for why it is a subprocess at all).
include(joinpath(@__DIR__, "test_qml.jl"))
