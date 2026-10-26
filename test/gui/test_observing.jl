# The Observe perspective's data layer.
#
# These assert ASTRONOMY, not self-consistency. A planning tool that agrees with itself and
# disagrees with the sky is worse than none, and the three bugs this perspective was blocked on
# (twilight in the wrong units, southern declinations with the wrong sign, RA in hours where
# degrees were documented) were all of that kind: plausible numbers, wrong answers. So the
# checks below are things independently known — Vega is a northern summer target, Canopus never
# rises from California, a June night is short — rather than comparisons against whatever the
# code happens to return.

@testset "Observe data layer" begin

    # Vega: high overhead at CHARA (dec +38.8 against latitude +34.2), a summer target.
    VEGA_RA, VEGA_DEC = 279.2347, 38.7837
    # Canopus: far south, and so never usefully above the horizon from California.
    CANOPUS_RA, CANOPUS_DEC = 95.988, -52.696
    JUNE, DECEMBER = DateTime(2026, 6, 21), DateTime(2026, 12, 21)

    @testset "the catalogue links a spectral setup to its combiner" begin
        c = config_catalog()
        @test "CHARA" in c.facilities
        @test !isempty(c.combiners)
        @test isempty(c.unknown)          # every shipped config classifies
        # wave_combiner is the ONLY machine-readable link between the two, and without it a
        # panel offers spectral setups the chosen instrument cannot produce.
        for (wave, comb) in c.wave_combiner
            @test wave in c.wavelengths
            @test comb in c.combiners
        end
        @test !isempty(c.wave_combiner)
    end

    @testset "the noise switches reach simulate, and the command log says so" begin
        # `shell_simulate` had no coverage at all, so nothing checked that a control on the
        # panel reaches the core call. These two are worth pinning in particular: both default
        # ON and both change the data, so a silent failure to forward them looks like a working
        # simulation with the wrong uncertainties.
        model = joinpath(@__DIR__, "..", "..", "demos", "models", "Uniform_Disk.toml")
        isfile(model) || (model = joinpath(@__DIR__, "..", "..", "demos", "models",
                                           "Limb_Darkened_Disk_Linear.toml"))
        @test isfile(model)
        tmp = mktempdir()
        # `shell_simulate` reaches the session through the global SHELL, so the test has to
        # install one. ShellState takes every field positionally, exactly as the console-pane
        # testset builds it; if a field is added there, this call has to follow.
        prev = GUI.SHELL[]
        GUI.SHELL[] = GUI.ShellState(Session(), nothing, nothing, nothing, String[], Any[], 0,
                          :uv, :baseline, false, :none, :none, false, "", String[], nothing, "",
                          nothing, nothing, nothing, nothing, nothing, nothing, nothing, Any[],
                          nothing, nothing, nothing, nothing, "", "", nothing)
        opts(extra) = join(vcat(["source\tmodel", "path\t" * model, "mag\t8", "mag_ao\t8",
                                 "seed\t1", "n_samples\t20"], extra), "\n")
        function run_sim(tag, extra, comb, wave)
            out = joinpath(tmp, tag * ".oifits")
            r = GUI.shell_simulate(out, "CHARA", "", comb, wave, "T", VEGA_RA, VEGA_DEC,
                                   "2026-06-21", -1.0, 1.0, 30.0, opts(extra))
            (r, out)
        end

        # SPICA is the only shipped combiner with a fringe tracker, so it is the only one where
        # the switch can be observed at all.
        r1, o1 = run_sim("ft_on",  ["fringe_tracker\t1"], "SPICA", "SPICA_LR")
        r2, o2 = run_sim("ft_off", ["fringe_tracker\t0"], "SPICA", "SPICA_LR")
        @test startswith(r1, "ok\t") && startswith(r2, "ok\t")
        e1 = readoifits(o1; T = Float64, filter_bad_data = false)[1, 1]
        e2 = readoifits(o2; T = Float64, filter_bad_data = false)[1, 1]
        # Tracking integrates longer, so the error bars must come out SMALLER at this magnitude.
        mean_err(d) = sum(Float64.(d.v2_err)) / length(d.v2_err)
        @test mean_err(e1) < mean_err(e2)

        # The systematic is drawn per baseline per night, so switching it off leaves the data
        # with less scatter about the same truth while the error bars are unchanged.
        r3, o3 = run_sim("sys_on",  ["systematics\t1"], "MIRCX", "MIRCX_LOWH")
        r4, o4 = run_sim("sys_off", ["systematics\t0"], "MIRCX", "MIRCX_LOWH")
        @test startswith(r3, "ok\t") && startswith(r4, "ok\t")
        e3 = readoifits(o3; T = Float64, filter_bad_data = false)[1, 1]
        e4 = readoifits(o4; T = Float64, filter_bad_data = false)[1, 1]
        @test e3.v2 != e4.v2                       # the draw reaches the data
        @test e3.v2_err ≈ e4.v2_err                # and not the error bars

        # The command log has to be runnable: a non-default switch that does not appear in it
        # produces a script that does not reproduce the run it claims to.
        sh  = GUI._shell()
        log = join([string(l) for l in sh.console], "\n")
        @test occursin("systematics = false", log)
        @test occursin("fringe_tracker = false", log)
        GUI.SHELL[] = prev
        rm(tmp; recursive = true, force = true)
    end

    @testset "an image source of any dimension, read by its keywords" begin
        tmp = mktempdir()
        img(nx) = [hypot(i - nx ÷ 2 - 1, j - nx ÷ 2 - 1) < 8 ? 1.0 : 0.0 for i in 1:nx, j in 1:nx]
        g = joinpath(tmp, "grey.fits"); writefits(img(33), g; pixsize = 0.1)
        i = GUI.simulate_source_info("image", g)
        @test i.ok && i.ndims == 2 && occursin("33×33 · 1 λ · 1 epoch · 0.1 mas", i.summary)
        w = read_wave_file("MIRCX_LOWH")
        c = joinpath(tmp, "cube.fits")
        writefits(cat((img(33) for _ in w.λ)...; dims = 3), c; pixsize = 0.1, wavelengths = w.λ)
        i = GUI.simulate_source_info("image", c)
        @test i.ok && i.ndims == 3 && i.nwav == length(w.λ) && occursin("$(length(w.λ)) λ (", i.summary)
        q = joinpath(tmp, "four.fits")
        writefits(cat((img(33) for _ in 1:6)...; dims = 3) |> x -> reshape(x, 33, 33, 2, 3), q;
                  pixsize = 0.1, wavelengths = [1.6e-6, 2.2e-6], mjds = [60000.0, 60001.0, 60002.0])
        i = GUI.simulate_source_info("image", q)
        @test i.ok && i.ndims == 4 && occursin("2 λ (1.6–2.2 µm) · 3 epochs (MJD 60000.0–60002.0)", i.summary)
        @test GUI.simulate_source_info("cube", q).ok                 # an older session's kind
        bare = joinpath(tmp, "bare.fits")
        OITOOLS.AstroFITS.writefits!(bare, ["OBJECT" => "x"], rand(33, 33, 3) .+ 0.1)
        i = GUI.simulate_source_info("image", bare)
        @test !i.ok && occursin("CTYPE3", i.summary) && occursin("FITS Image Handling", i.summary)

        # through the panel: the wavelength cube, simulated by path, logged with it
        prev = GUI.SHELL[]
        GUI.SHELL[] = GUI.ShellState(Session(), nothing, nothing, nothing, String[], Any[], 0,
                          :uv, :baseline, false, :none, :none, false, "", String[], nothing, "",
                          nothing, nothing, nothing, nothing, nothing, nothing, nothing, Any[],
                          nothing, nothing, nothing, nothing, "", "", nothing)
        out = joinpath(tmp, "sim.oifits")
        r = GUI.shell_simulate(out, "CHARA", "", "MIRCX", "MIRCX_LOWH", "T", VEGA_RA, VEGA_DEC,
                               "2026-06-21", -1.0, 1.0, 30.0,
                               join(["source\timage", "path\t" * c, "pixsize\t0.1", "mag\t5",
                                     "mag_ao\t5", "seed\t1", "n_samples\t20"], "\n"))
        @test startswith(r, "ok\t")
        log = join(GUI._shell().console, "\n")
        @test occursin("image = \"$(c)\", pixsize = 0.1", log)
        r = GUI.shell_simulate(joinpath(tmp, "no.oifits"), "CHARA", "", "MIRCX", "MIRCX_LOWH", "T",
                               VEGA_RA, VEGA_DEC, "2026-06-21", -1.0, 1.0, 30.0,
                               join(["source\timage", "path\t" * bare, "pixsize\t0.1"], "\n"))
        @test startswith(r, "!") && occursin("CTYPE3", r)

        # the image's own channels, for an image and for a model holding the cube
        own = joinpath(tmp, "own.fits")
        λc = [1.52e-6, 1.6e-6, 1.68e-6]
        writefits(cat((img(33) for _ in λc)...; dims = 3), own; pixsize = 0.1, wavelengths = λc)
        mfile = joinpath(tmp, "star_cube.toml")
        write(mfile, """
            [model]
            "star,ud" = 0.5
            "star,f" = 0.6
            "env,image" = "own.fits"
            "env,f" = 0.4
            """)
        for (src, path) in (("image", own), ("model", mfile))
            o = joinpath(tmp, "own_$(src).oifits")
            r = GUI.shell_simulate(o, "CHARA", "", "MIRCX", "MIRCX_LOWH", "T", VEGA_RA, VEGA_DEC,
                                   "2026-06-21", -1.0, 1.0, 30.0,
                                   join(["source\t" * src, "path\t" * path, "pixsize\t0.1",
                                         "image_channels\t1", "mag\t5", "mag_ao\t5", "seed\t1",
                                         "n_samples\t20"], "\n"))
            @test startswith(r, "ok\t")
            x = readoifits(o; T = Float64, warn = false, verbose = false)[1, 1]
            @test isapprox(sort(unique(x.uv_lam)), λc; rtol = 1e-6)
        end
        @test occursin("simulate(facility, target, combiner, :image, dates", join(GUI._shell().console, "\n"))
        # fractional noise: the error bars are the fractions asked, and the call is logged
        o = joinpath(tmp, "frac.oifits")
        r = GUI.shell_simulate(o, "CHARA", "", "MIRCX", "MIRCX_LOWH", "T", VEGA_RA, VEGA_DEC,
                               "2026-06-21", -1.0, 1.0, 30.0,
                               join(["source\timage", "path\t" * own, "pixsize\t0.1", "image_channels\t1",
                                     "noise\tfractional", "frac_v2\t2", "frac_visamp\t1", "frac_t3amp\t1",
                                     "frac_visphi\t0.5", "frac_t3phi\t0.3", "frac_flux\t1", "seed\t1",
                                     "observability\t1", "alt_limit\t25", "alt_max\t85"], "\n"))
        @test startswith(r, "ok\t")
        x = readoifits(o; T = Float64, warn = false, verbose = false, filter_bad_data = false)[1, 1]
        @test all(x.t3phi_err .≈ 0.3) && all(x.v2_err .≥ 0.02 * 0.01 - 1e-12)
        @test occursin("noise = FractionalNoise(v2 = 0.02, t3phi = 0.3)", join(GUI._shell().console, "\n"))
        @test occursin("observability = (min_elevation = 25.0, max_elevation = 85.0)", join(GUI._shell().console, "\n"))
        GUI.SHELL[] = prev
        rm(tmp; recursive = true, force = true)
    end

    @testset "telescope order is the facility's own" begin
        tels = facility_telescopes("CHARA")
        @test tels == ["S1", "S2", "E1", "E2", "W1", "W2"]
        # compute_delays indexes into exactly this order, so a config vector built from names
        # has to line up with it positionally.
        @test telescope_config(["S1", "E1"], tels) == [1, 0, 1, 0, 0, 0]
        @test telescope_config(String[], tels)     == zeros(Int, 6)
        @test telescope_config(tels, tels)         == ones(Int, 6)
        # the delay reference, marked as the panel sends it
        @test telescope_config(["S1", "E1", "W2*"], tels) == [1, 0, 1, 0, 0, 2]
    end

    @testset "the delay reference: the most central telescope, on the equator's plane" begin
        f = read_facility_file("CHARA")
        tels = facility_telescopes("CHARA")
        @test tels[central_telescope(f, ones(Int, 6))] == "W2"
        @test tels[central_telescope(f, [1, 1, 1, 1, 1, 0])] == "S2"
        @test tels[central_telescope(f, telescope_config(["E1", "E2"], tels))] in ("E1", "E2")
        @test central_telescope(f, zeros(Int, 6)) == 0
        @test GUI.shell_central_telescope("CHARA", "S1 S2 E1 E2 W1") == "S2"
        @test GUI.shell_central_telescope("CHARA", "") == ""
        # the shell's single-setting search takes the reference from the marked list: Capella's
        # best setting with W2 as reference (test_pop_plan.jl)
        rows = GUI.shell_best_pops("CHARA", 45.998, "2026-09-15T00:00:00", 79.172,
                                   "S1 S2 E1 E2 W1 W2*", 1)
        @test startswith(rows, "1 1 4 4 1 5\t")
    end

    @testset "Vega is a summer target from CHARA, and not a winter one" begin
        june = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE)
        dec  = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, DECEMBER)
        # dec +38.8 against latitude +34.2 transits within 5 degrees of the zenith
        @test maximum(june.alt) > 80
        @test observable_hours(june) > 3
        # In December it barely scrapes the 25-degree floor -- 28.5 at transit -- so it is
        # nominally "up" for a few minutes and useless. The check is the ratio, not zero:
        # asserting zero would have been an accident of the old 30-degree default.
        @test maximum(dec.alt) < 30
        @test length(dec.good_alt) < length(june.good_alt) / 10
        @test observable_hours(dec) < 0.5
    end

    @testset "a far-southern target never rises from a northern site" begin
        p = night_plan("CHARA", "Canopus", CANOPUS_RA, CANOPUS_DEC, DateTime(2026, 2, 1))
        # dec -52.7 from latitude +34.2 cannot exceed 90 - 34.2 - 52.7 ≈ 3 degrees
        @test maximum(p.alt) < 5
        @test observable_hours(p) == 0
    end

    @testset "the June night is short, and it is what limits the target" begin
        # alt_max = 90 on purpose: with the default 85 the zenith cap ALSO bites (Vega transits
        # at 85.4), and this test is about twilight being computed rather than assumed. The cap
        # gets its own test below.
        p = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE; alt_max = 90.0)
        @test length(p.good_alt) > length(p.good_twilight)
        @test observable_hours(p) ≈ length(p.good_twilight) / 60 rtol = 1e-6
        # and a June night at mid-latitude is short
        @test 3 < length(p.good_twilight) / 60 < 8
    end

    @testset "the defaults are ASPRO's, and the zenith cap is not decorative" begin
        @test DEFAULT_ALT_LIMIT == 25.0
        @test DEFAULT_ALT_MAX   == 85.0
        d = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE)
        @test d.alt_limit == 25.0 && d.alt_max == 85.0
        # Vega transits at 85.4 from CHARA, so the cap removes real time rather than nothing
        open_ = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE; alt_max = 90.0)
        @test observable_hours(d) < observable_hours(open_)
    end

    @testset "the delay check is opt-in, and says when it did not run" begin
        # Applying it with an unsearched POP reports far less time than is available, so it
        # must not be the default — and "not checked" must not look like "passed".
        p = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE)
        @test p.delay_applied == false
        @test p.good_delay == collect(1:length(p.ha))   # covers the night: not a constraint

        cfg = telescope_config(["S1", "S2", "E1", "E2"], facility_telescopes("CHARA"))
        q = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE; config = cfg, use_delay = true)
        @test q.delay_applied == true
        @test observable_hours(q) <= observable_hours(p)   # a constraint can only remove time
    end

    @testset "searching the POPs buys real time" begin
        cfg = telescope_config(["S1", "S2", "E1", "E2"], facility_telescopes("CHARA"))
        base = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE;
                          config = cfg, use_delay = true)               # POP 1 throughout
        pops = best_pops("CHARA", VEGA_DEC, base.ha, cfg; n = 3)
        @test length(pops) == 3
        @test issorted([p.score for p in pops]; rev = true)             # best first
        best = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE;
                          config = cfg, pop = pops[1].pop, use_delay = true)
        # This is the whole point of the POP search: an arbitrary configuration wastes the night
        @test observable_hours(best) > observable_hours(base)
        @test occursin("\t", pop_rows(pops))
    end

    @testset "observable hours do not depend on how finely the night was sampled" begin
        # Every duration on a plan is a COUNT of grid samples, so it means nothing without the
        # step. When the step had to be supplied by the caller, three of the four callers did
        # not: the panel samples every ten minutes and the summary divided by one, so Vega came
        # back at 0.43 h instead of 4.33 and read as a target that cannot be observed.
        fine   = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE; step_minutes = 1)
        coarse = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE; step_minutes = 10)
        @test observable_hours(coarse) > 3
        @test observable_hours(fine) ≈ observable_hours(coarse) rtol = 0.05
        # and the table agrees with the plan it was built from
        row = split(split(plan_rows([coarse]), "\n")[1], "\t")
        @test parse(Float64, row[4]) ≈ observable_hours(coarse) atol = 0.01
    end

    @testset "the now marker belongs to the night in progress" begin
        # LST repeats every day, so the current LST lands inside almost any night's span. A red
        # line saying "you are here" on a chart of a night months away is simply false.
        # A year out, rather than a fixed past date: run the suite on the day that date names
        # and a fixed one would be the night in progress, once a year.
        next_june = DateTime(Dates.year(Dates.now()) + 1, 6, 21)
        far = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, next_june)
        @test !far.is_tonight
        @test !isfinite(gantt_geometry(far).now)
    end

    @testset "the summary table reports a window, not just a total" begin
        rows = split(plan_rows([night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE),
                                night_plan("CHARA", "Canopus", CANOPUS_RA, CANOPUS_DEC, JUNE)]), "\n")
        @test length(rows) == 2
        vega = split(rows[1], "\t")
        @test vega[1] == "Vega"
        @test parse(Float64, vega[4]) > 3          # hours
        @test parse(Float64, vega[2]) != parse(Float64, vega[3])   # a real window
        # a target that never rises still gets a row, with zero hours rather than no row
        @test split(rows[2], "\t")[1] == "Canopus"
        @test parse(Float64, split(rows[2], "\t")[4]) == 0
    end

    @testset "SIMBAD resolution answers in the panel's own format" begin
        # The reply is parsed by ObserveTab.qml's applySimbad, which splits on tabs and
        # branches on field 1. Both branches are asserted here because a malformed reply
        # would be read by QML as a successful resolve and write NaNs into the row.
        @test startswith(GUI.shell_simbad(""), "bad\t")
        @test startswith(GUI.shell_simbad("   "), "bad\t")

        # The network is not a test dependency. This runs only when asked for, because a
        # suite that fails on a SIMBAD outage is a suite people learn to ignore.
        if get(ENV, "OITOOLS_SIMBAD_TESTS", "0") == "1"
            f = split(GUI.shell_simbad("Vega"), "\t")
            @test f[1] == "ok"
            @test length(f) == 9
            @test parse(Float64, f[2]) ≈ VEGA_RA  atol = 1e-3   # degrees, not hours
            @test parse(Float64, f[3]) ≈ VEGA_DEC atol = 1e-3
            @test f[8] == "* alf Lyr"
            @test startswith(f[9], "A0")

            # Southern declinations keep their sign: this is P0-2, and a sexagesimal parser
            # that applies the minus to the degrees field alone gets it wrong by 2×(′+″).
            g = split(GUI.shell_simbad("Alpha Cen A"), "\t")
            @test g[1] == "ok"
            @test parse(Float64, g[3]) < -60
            # α Cen A has no published R or I magnitude. That must arrive as an empty field,
            # not as 0.0 — magnitude zero is Vega-bright, and writing it for "unknown" would
            # silently make an unmeasured band the brightest in the row.
            t = simbad_target("Alpha Cen A")
            @test isnan(t.mags["R"])
            @test t.mags["K"] < 0
            @test t.plx > 700                       # the nearest star system

            @test_throws Exception simbad_target("NoSuchStarXYZ")
        end

        # Rejected before the network, so this costs nothing and runs everywhere
        @test_throws ArgumentError simbad_target("  ")
        @test issorted(indexin(["V", "J", "H", "K"], SIMBAD_BANDS))   # panel order
    end

    @testset "the observing band comes from the setup's wavelengths, not its name" begin
        band(x) = split(GUI.shell_sim_band(x), "\t")

        h = band("MIRCX_LOWH")
        @test h[1] == "H"
        @test parse(Float64, h[2]) ≈ 1.58 atol = 0.05      # µm

        k = band("MYSTIC_LOWK")
        @test k[1] == "K"
        @test parse(Float64, k[2]) ≈ 2.19 atol = 0.05

        # The name is not the band. MIRCX_LOWJ's channels sit at 1.44 µm, which is nearer H
        # than J -- so a panel that read the letter out of the setup name would feed the noise
        # model a J magnitude for an H-band observation.
        j = band("MIRCX_LOWJ")
        @test parse(Float64, j[2]) < 1.5
        @test j[1] != "J"

        # An unknown setup is an ordinary answer, not an error: the panel shows no band and
        # leaves whatever magnitude is already there alone.
        @test GUI.shell_sim_band("no such setup") == "\t"
        @test GUI.shell_sim_band("") == "\t"
        @test GUI.shell_sim_band("   ") == "\t"
    end

    @testset "the moon is reported even when it does not exclude anything" begin
        p = night_plan("CHARA", "Vega", VEGA_RA, VEGA_DEC, JUNE)
        @test 0 <= p.moon_fli <= 1
        @test length(p.moon_sep) == length(p.alt)
        @test all(0 .<= p.moon_sep .<= 180)
        @test occursin("moon", sprint(show, p))
    end

    @testset "a multi-POP plan: segments, changes, and the charts drawn from them" begin
        ra, dec, d = 279.2347, 38.7837, DateTime(2026, 6, 21)
        pp = night_pop_plan("CHARA", ra, dec, d; max_changes = 2, max_drop = 1)
        @test pp.changes == 1 && pp.minutes[1] > 0 && pp.minutes[2] > 0
        p = night_plan("CHARA", "Vega", ra, dec, d; use_delay = true, popplan = pp)
        @test is_multi_pop(p)
        @test length(p.changes) == pp.changes
        # the minutes observed are the plan's, on the night's grid
        @test abs(length(observable_indices(p)) - sum(pp.minutes)) <= 2 * length(pp.segments)
        # nothing is observed while the POPs change ...
        for (c0, c1) in p.changes
            @test c1 - c0 == 10
            @test isempty(intersect(observable_indices(p), c0:c1-1))
        end
        # E1 dropped in the last segment, after the last change: no second change for it
        @test p.segments[end].tier == 2 && p.segments[end].config[3] == 0
        # E1 W2 W1 S2 S1 E2, the dropped E1 a cross; "−E1" where it is lost
        @test segment_label(p, p.segments[end]) == "× " *
              join(string.(p.segments[end].pop[[6, 5, 2, 1, 4]]), " ")
        @test transition_label(p, length(p.segments)) == "−E1"
        @test transition_label(p, 1) == ""

        geo = gantt_geometry(p)
        @test geo.multi
        # the first setting where it starts, the change as "→" where it starts
        @test length(geo.poplabels) == 2 && !geo.poplabels[1][5] && geo.poplabels[2][5]
        @test geo.poplabels[2][1] ≈ unwrap_lst(p.lst)[first(only(p.changes))]
        # "POP" where it changes and the time it starts just before, "DROP" before E1 is lost
        # and the time it is lost
        @test count(m -> m[2] == "POP", geo.marks) == 1
        @test count(m -> m[2] == "DROP", geo.marks) == 1
        x = unwrap_lst(p.lst)
        c0 = first(only(p.changes))
        @test any(m -> m[3] === :end && m[1] ≈ x[c0] && occursin(r"^\d\d:\d\d$", m[2]), geo.marks)
        @test any(l -> l.kind === :time && l.x ≈ x[p.segments[end].i0], geo.labels)
        @test any(geo.poplabels[2][4]) && !any(geo.poplabels[1][4])
        # Every character these labels draw is pre-warmed in the glyph atlas: one inserted
        # once Qt owns the GL context breaks the plot that needed it.
        drawn = vcat([l[3] for l in geo.poplabels], [l[3] for l in geo.seglabels],
                     [m[2] for m in geo.marks], [geo.subtitle, "→ ", "EL", "×", "+E1"])
        @test all(c -> c in GUI.PLOT_GLYPHS, join(drawn))
        # Vega passes above the soft limit either side of transit
        @test length(geo.ellabels) == count(b -> b.x1 - b.x0 >= 0.3, geo.softbars) >= 1
        @test startswith(geo.subtitle, string(length(p.segments), " segments · 1 change"))
        @test any(l -> l[3] == "−E1", geo.seglabels)
        det = gantt_geometry(p; detailed = true)
        @test any(r -> r[2] == "POPs", det.rows) && length(det.switches) == 1
        # ... but the Gantt draws no gap for it: the bar before runs on to the one after
        c0, c1 = only(p.changes)
        x = unwrap_lst(p.lst)
        @test any(b -> b.x1 ≈ x[c1], geo.bars) && any(b -> b.x0 ≈ x[c1], geo.bars)
        # the delay curves break while the carts move
        dg = delay_plot_geometry(p)
        @test all(c -> all(isnan, c.y[c0:c1-1]) && all(isfinite, c.y[1:c0-1]), dg.curves)
        @test length(dg.switches) == 1 && all(isfinite, dg.ylim)

        # one setting: the plan and its chart are what they were
        p1 = night_plan("CHARA", "Vega", ra, dec, d; use_delay = true, pop = [1, 1, 2, 3, 1, 5])
        @test !is_multi_pop(p1) && isempty(p1.changes)
        @test pop_label(p1) == "E1:P2 W2:P5 W1:P1 S2:P1 S1:P1 E2:P3"     # the panel's order
        g1geo = gantt_geometry(p1)
        @test !g1geo.multi && isempty(g1geo.poplabels) && isempty(g1geo.seglabels)
    end

    @testset "POP plans as text, as a summary, and ranked for the POP strategies tab" begin
        ra, dec, d = 279.2347, 38.7837, DateTime(2026, 6, 21)
        f = read_facility_file("CHARA")
        pp = night_pop_plan(f, ra, dec, d; max_changes = 2, max_drop = 1)
        q = parse_pop_plan_text(pop_plan_text(pp))
        @test q.changes == pp.changes && q.minutes == pp.minutes
        @test [sg.pop for sg in q.segments] == [sg.pop for sg in pp.segments]
        @test [sg.config for sg in q.segments] == [sg.config for sg in pp.segments]
        sm = pop_plan_summary(pp, f, ra, ones(Int, 6))
        @test occursin(" → ", sm) && occursin(r"−E1 \d+:\d\d", sm)
        # a telescope missing from the start is a cross, not a drop
        hd = (345.067176, -22.524347, DateTime(2026, 10, 26))    # SIMBAD
        php = night_pop_plan(f, hd...; max_changes = 2, max_drop = 1)
        shd = pop_plan_summary(php, f, hd[1], ones(Int, 6))
        @test startswith(shd, "× ") && !occursin("−", shd)
        p_hd = night_plan("CHARA", "HD 217357", hd...; use_delay = true, popplan = php)
        g_hd = gantt_geometry(p_hd)
        @test startswith(only(g_hd.poplabels)[3], "× ") && isempty(g_hd.seglabels)
        @test !any(m -> m[2] == "DROP", g_hd.marks)
        rows = split(GUI.shell_pop_plans("CHARA", dec, "2026-06-21T00:00:00", ra,
                                     "S1 S2 E1 E2 W1 W2", "W2", 2, 1, 10), "\n")
        @test length(rows) >= 2
        f6 = [parse(Float64, split(r, "\t")[4]) for r in rows]
        tot = [parse(Float64, split(r, "\t")[5]) for r in rows]
        # best first: full-array hours never increase down the table
        @test issorted(f6; rev = true)
        @test all(r -> length(split(r, "\t")) == 6, rows)
        @test first(tot) == maximum(tot)
        # the plan column is what shell_gantt takes
        @test parse_pop_plan_text(split(first(rows), "\t")[6]).minutes == pp.minutes
    end

    @testset "a night's POP timeline" begin
        ra, dec, d = 279.2347, 38.7837, DateTime(2026, 6, 21)
        f = read_facility_file("CHARA")
        pp = night_pop_plan(f, ra, dec, d; max_changes = 2, max_drop = 1)
        p = night_plan("CHARA", "Vega", ra, dec, d; use_delay = true, popplan = pp)
        p1 = night_plan("CHARA", "Vega", ra, dec, d; use_delay = true, pop = [1, 1, 2, 3, 1, 5])
        obs(q) = Set(observable_indices(q))
        spans(rows) = [r for r in rows if !(r.kind in (:drop, :readd))]
        # The rows cover the dark hours once, in order; the settings are the observed samples.
        function tiles(q, rows)
            sp = spans(rows)
            @test first(sp).i0 == minimum(q.good_twilight)
            @test last(sp).i1 == maximum(q.good_twilight)
            @test all(sp[k].i0 == sp[k-1].i1 + 1 for k in 2:length(sp))
            @test sum(r.hours for r in rows if r.kind === :setting) ≈ observable_hours(q)
            @test all(r -> r.kind !== :setting || all(in(obs(q)), r.i0:r.i1), rows)
            @test all(r -> r.kind !== :gap || !any(in(obs(q)), r.i0:r.i1), rows)
            # each gap says what excludes it: elevation first, then the Moon, then the delay
            for r in filter(r -> r.kind === :gap, rows), i in r.i0:r.i1
                want = !(i in q.good_alt) ? :elevation : !(i in q.good_moon) ? :moon : :delay
                @test r.why === want
            end
            # each span ends where the next begins, on the clock
            @test all(sp[k].t0 ≈ sp[k-1].t1 for k in 2:length(sp))
        end

        rows = pop_plan_timeline(p)
        tiles(p, rows)
        kinds = [r.kind for r in rows]
        # one POP change, 10 minutes, listing only the carts that move
        ch = only(filter(r -> r.kind === :change, rows))
        @test ch.hours ≈ 10 / 60
        k = findfirst(==(:change), kinds)
        before, after = rows[k-1], rows[findnext(==(:setting), kinds, k)]
        @test !isempty(ch.moved)
        @test all(((t, a, b),) -> before.pop[t] == a && after.pop[t] == b && a != b, ch.moved)
        @test length(ch.moved) == count(t -> before.pop[t] != after.pop[t], findall(>(0), p.config))
        # E1 dropped with the POPs unchanged: a drop row, no change row next to it
        kd = only(findall(==(:drop), kinds))
        @test rows[kd].tels == [3] && rows[kd+1].kind === :setting && rows[kd+1].config[3] == 0
        @test rows[kd-1].kind === :setting && rows[kd-1].pop == rows[kd+1].pop
        @test isnan(rows[kd].t1) && rows[kd].t0 ≈ rows[kd+1].t0

        # one setting: settings and gaps only
        r1 = pop_plan_timeline(p1)
        tiles(p1, r1)
        @test all(r -> r.kind in (:setting, :gap), r1)
        @test all(r -> r.kind !== :setting || r.pop == p1.pop, r1)

        # the clocks move the times and nothing else
        x = unwrap_lst(p.lst)
        for sys in (:utc, :local)
            rs = pop_plan_timeline(p; system = sys)
            @test [r.kind for r in rs] == kinds
            @test [r.hours for r in rs] == [r.hours for r in rows]
            @test isequal([r.ha0 for r in rs], [r.ha0 for r in rows])
            @test all(r -> r.t0 ≈ GUI._gantt_clock(p, x[r.i0], sys), rs)
        end
        u = pop_plan_timeline(p; system = :utc); l = pop_plan_timeline(p; system = :local)
        off = l[1].t0 - u[1].t0
        @test off ≈ GUI.local_utc_offset(f, d) && all(k -> l[k].t0 - u[k].t0 ≈ off, eachindex(u))
        # local midnight falls inside this night: the times increase through it
        t0s = [r.t0 for r in spans(l)]
        @test issorted(t0s) && first(t0s) < 24 < last(t0s)
        txt = timeline_text(l, p; system = :local)
        @test occursin(r"\n\s+23:\d\d", txt) || occursin(r"\n\s+0:\d\d", txt)
        @test startswith(txt, "LOCAL") && occursin("POP change", txt) && occursin("drop E1", txt)

        # synthetic plans: the same night with other segments
        remake(q; kw...) = NightPlan((get(kw, fn, getfield(q, fn)) for fn in fieldnames(NightPlan))...)
        s1 = p.segments[1]
        a = s1.i0; b = s1.i1; m = (a + b) ÷ 2; q3 = div(b - m, 2) + m
        drop = copy(s1.config); drop[3] = 0
        seg(i0, i1, cfg, pop) = GUI.PlanSegment((i0, i1, cfg, pop, 1))
        # a telescope dropped, then added back, POPs unchanged
        pr = remake(p; segments = [seg(a, m, s1.config, s1.pop), seg(m + 1, q3, drop, s1.pop),
                                   seg(q3 + 1, b, s1.config, s1.pop), p.segments[2:end]...])
        rr = pop_plan_timeline(pr)
        kr = [r.kind for r in rr]
        i = findfirst(==(:drop), kr)
        @test rr[i].tels == [3] && kr[i+1] === :setting && kr[i+2] === :readd && rr[i+2].tels == [3]
        tiles(pr, rr)
        # a POP change and a drop at the same boundary: the change, then the drop
        s2 = p.segments[2]
        d2 = copy(s2.config); d2[3] = 0
        pc = remake(p; segments = [s1, seg(s2.i0, s2.i1, d2, s2.pop), p.segments[3:end]...])
        kc = [r.kind for r in pop_plan_timeline(pc)]
        j = findfirst(==(:change), kc)
        @test kc[j+1:j+3] == [:drop, :setting, :readd] || kc[j+1:j+2] == [:drop, :setting]
        # Capella through LST 0 h: settings in time order, every POP change a change row
        cra, cdec, cd = 79.172, 45.998, DateTime(2026, 9, 15)
        cp = night_pop_plan(f, cra, cdec, cd; max_changes = 2, max_drop = 1)
        cq = night_plan("CHARA", "Capella", cra, cdec, cd; use_delay = true, popplan = cp)
        crows = pop_plan_timeline(cq)
        tiles(cq, crows)
        cs = filter(r -> r.kind in (:setting, :change), crows)
        @test !any(k -> cs[k-1].kind === :setting && cs[k].kind === :setting &&
                        cs[k-1].pop[cq.config .> 0] != cs[k].pop[cq.config .> 0], 2:length(cs))
        @test count(r -> r.kind === :change, crows) == cp.changes
        @test issorted([r.ha0 for r in crows if r.kind === :setting])
        @test first(crows).kind === :gap && first(crows).why === :elevation
        @test occursin("elevation-restricted", timeline_text(crows, cq))
        # the candidates chart: single settings and multi-POP plans as stacked rows, a plan's
        # row with its own POPs, changes and drops
        g2 = gantt_stack_geometry([p1, p], ["#1", "#2"], "Vega · 2 candidates")
        @test [r[2] for r in g2.rows] == ["#1", "#2"] && g2.rows[1][1] > g2.rows[2][1]
        @test count(m -> m[2] == "POP", g2.marks) == length(p.changes)
        @test count(m -> m[2] == "DROP", g2.marks) == 1
        @test length(g2.poplabels) == 1 + count(l -> true, gantt_geometry(p).poplabels)
        @test all(m -> m[4] ≈ g2.rows[2][1], g2.marks)          # on the plan's row
        @test GUI.shell_plans_gantt("CHARA", "Vega", ra, dec, "2026-06-21T00:00:00",
                                    "S1 S2 E1 E2 W1 W2", "#1\t1 1 2 3 1 5\t\n#2\t\t" * pop_plan_text(pp),
                                    true) |> x -> parse(Int, x) > 0
        @test startswith(GUI.shell_plans_gantt("CHARA", "Vega", ra, dec, "2026-06-21T00:00:00",
                                               "", "", true), "!")
        # the delay reference: `*` after its POP; "Ref → E2" with its time where it changes and
        # the POPs do not; with a POP change, only the new POPs, `*` on the new reference
        w2 = copy(s1.config); w2[6] = 2
        e2 = copy(s1.config); e2[4] = 2                     # E2 the reference, W2 still in
        nw = copy(e2); nw[6] = 0                            # W2 dropped, E2 the reference
        pf = remake(p; config = w2, segments = [seg(a, m, w2, s1.pop), seg(m + 1, q3, e2, s1.pop),
                                                 seg(q3 + 1, b, nw, s1.pop)], changes = Tuple{Int,Int}[])
        @test segment_label(pf, pf.segments[1]) == join([string(s1.pop[i], i == 6 ? "*" : "") for i in [3, 6, 5, 2, 1, 4]], " ")
        @test occursin("*", segment_label(pf, pf.segments[2])) && endswith(segment_label(pf, pf.segments[2]), "*")
        @test transition_label(pf, 2) == "Ref → E2"
        @test transition_label(pf, 3) == "−W2"                 # E2 was the reference already
        @test GUI.reference_change(pf, 2) == 4 && GUI.reference_change(pf, 3) == 0
        gf = gantt_geometry(pf)
        xf = unwrap_lst(pf.lst)
        @test any(l -> l[3] == "Ref → E2" && l[1] ≈ xf[m + 1], gf.seglabels)
        @test any(l -> l.kind === :time && l.x ≈ xf[m + 1], gf.labels)
        @test all(c -> c in GUI.PLOT_GLYPHS, "Ref → E2*")
        rf = pop_plan_timeline(pf)
        kf = [r.kind for r in rf]
        i = findfirst(==(:ref), kf)
        @test rf[i].tels == [4] && kf[i+1] === :setting && isnan(rf[i].t1)
        @test occursin("Ref → E2", timeline_text(rf, pf))
        @test GUI.timeline_fields(rf[i+1], pf)[4] == segment_label(pf, pf.segments[2])
        # reference and POPs both change: no "Ref", the new POPs carry the `*`
        p2 = copy(s1.pop); p2[3] = s1.pop[3] == 1 ? 2 : 1
        pc2 = remake(pf; segments = [seg(a, m, w2, s1.pop), seg(m + 1, b, e2, p2)])
        @test transition_label(pc2, 2) == "" && GUI.reference_change(pc2, 2) == 0
        @test !any(==(:ref), [r.kind for r in pop_plan_timeline(pc2)])
        # the single-setting title and the stacked rows mark it too
        pw = night_plan("CHARA", "Vega", ra, dec, d; use_delay = true, pop = [1, 1, 2, 3, 1, 5],
                        config = [1, 1, 1, 1, 1, 2])
        @test occursin("W2:P5*", pop_label(pw))
        @test any(l -> occursin("5*", l[3]), gantt_stack_geometry([pw], ["#1"]).poplabels)
        # no darkness: no rows
        @test isempty(pop_plan_timeline(remake(p; good_twilight = Int[])))

        # the shell builds the same night as the Gantt and gives the same rows
        txtplan = pop_plan_text(pp)
        sh = split(GUI.shell_pop_timeline("CHARA", "Vega", ra, dec, "2026-06-21T00:00:00",
                                          "S1 S2 E1 E2 W1 W2", "", true, DEFAULT_ALT_LIMIT,
                                          DEFAULT_ALT_MAX, "local", txtplan, 10), "\n")
        q = night_plan("CHARA", "Vega", ra, dec, d; use_delay = true,
                       popplan = parse_pop_plan_text(txtplan))
        @test sh == [join(GUI.timeline_fields(r, q), "\t") for r in pop_plan_timeline(q; system = :local)]
        @test startswith(GUI.shell_pop_timeline("CHARA", "Vega", ra, dec, "2026-06-21T00:00:00",
                                                "S1 S2 E1 E2 W1 W2", "", true, DEFAULT_ALT_LIMIT,
                                                DEFAULT_ALT_MAX, "utc", txtplan, 10, true), "UTC")
        @test startswith(GUI.shell_pop_timeline("CHARA", "Vega", ra, dec, "not a date", "", "", true), "!")

        # hour angles on a clock; summaries on a clock
        h = -2.25
        @test GUI._hhmm(hour_angle_clock(p, h)) == GUI._hhmm(ra / 15 + h)
        @test hour_angle_clock(p, h; system = :utc) ≈ GUI.lst_to_utc(p, hour_angle_clock(p, h))
        cfg = ones(Int, 6)
        @test pop_plan_summary(pp, f, ra, cfg; clock = identity) == pop_plan_summary(pp, f, ra, cfg)
        rows4 = split(GUI.shell_best_pops("CHARA", dec, "2026-06-21T00:00:00", ra,
                                          "S1 S2 E1 E2 W1 W2", 2, DEFAULT_ALT_LIMIT,
                                          DEFAULT_ALT_MAX, "utc"), "\n")
        @test all(r -> length(split(r, "\t")) == 4 && occursin(r"^\d+:\d\d – \d+:\d\d$", split(r, "\t")[4]), rows4)
    end

    @testset "the top single POP settings, one row each" begin
        ra, dec, d = 279.2347, 38.7837, DateTime(2026, 6, 21)
        f = read_facility_file("CHARA")
        ha, ok = GUI._night_sky(f, ra, dec, d, ones(Int, 6))
        res = best_pops(f, dec, ha[ok], ones(Int, 6); n = 3)
        plans = [night_plan("CHARA", "Vega", ra, dec, d; pop = r.pop, use_delay = true) for r in res]
        geo = gantt_stack_geometry(plans, ["#1", "#2", "#3"])
        @test geo.stack && length(geo.rows) == 3
        # the first at the top, rows STACK_ROW_STEP apart
        ys = [r[1] for r in geo.rows]
        @test ys == sort(ys; rev = true) && all(≈(GUI.STACK_ROW_STEP), -diff(ys))
        # each row's POPs above its bar, in the panel's order
        @test length(geo.poplabels) == 3
        @test geo.poplabels[1][3] == join(string.(plans[1].pop[[3, 6, 5, 2, 1, 4]]), " ")
        # every bar sits on one of the rows
        @test all(b -> any(r -> abs(b.y - r[1]) < 1e-9, geo.rows), geo.bars)
        @test length(geo.transits) == 3
    end

    @testset "the settings that cover the night's thirds best" begin
        ra, dec, d = 279.2347, 38.7837, DateTime(2026, 6, 21)
        f = read_facility_file("CHARA")
        ha, ok = GUI._night_sky(f, ra, dec, d, ones(Int, 6))
        cov = coverage_pops(f, dec, ha[ok], ones(Int, 6))
        @test [c.part for c in cov] == ["early", "middle", "late"]
        # each is in delay inside its own third
        h = sort(Float64.(ha[ok])); lo, hi = extrema(h); w = (hi - lo) / 3
        for (k, c) in enumerate(cov)
            a, b = lo + (k - 1) * w, k == 3 ? hi : lo + k * w
            sub = Float32[x for x in h if a <= x <= b]
            @test !isempty(in_delay(f, dec, sub, ones(Int, 6), c.pop).good_delay)
        end
        plans = [night_plan("CHARA", "Vega", ra, dec, d; pop = c.pop, use_delay = true) for c in cov]
        geo = gantt_stack_geometry(plans, ["early", "middle", "late"], "best time coverage POPs")
        @test geo.subtitle == "best time coverage POPs" && [r[2] for r in geo.rows] == ["early", "middle", "late"]
        @test isempty(coverage_pops(f, dec, Float32[], ones(Int, 6)))
        # each pick is a best setting for its third, and among the ties the picks move the
        # fewest POPs from third to third: W2 stays on POP 5 throughout
        for (k, c) in enumerate(cov)
            a, b = lo + (k - 1) * w, k == 3 ? hi : lo + k * w
            sub = Float32[x for x in h if a <= x <= b]
            @test c.score == best_pop(f, dec, sub, ones(Int, 6); n_best = 1)[1].score
        end
        @test sum(count(i -> cov[k-1].pop[i] != cov[k].pop[i], 1:6) for k in 2:3) == 3
        @test all(c -> c.pop[6] == 5, cov)
    end

    @testset "a target list to and from CSV" begin
        tmp = mktempdir()
        rows = join(["Vega\t279.234735\t38.783689\t0.03\t-0.18\t-0.03\t0.13\t1",
                     "HD 217357, \"GJ 884\"\t345.067176\t-22.524347\t7.869\t\t\t4.478\t0"], "\n")
        p = joinpath(tmp, "targets")
        @test GUI.shell_save_targets(p, rows) == ""
        @test isfile(p * ".csv")
        back = split(GUI.shell_load_targets(p * ".csv"), "\n")
        @test length(back) == 2
        f1, f2 = split(back[1], "\t"), split(back[2], "\t")
        @test f1[1] == "Vega" && f1[8] == "1" && parse(Float64, f1[2]) == 279.234735
        # a name with a comma and quotes survives; unknown magnitudes stay unknown
        @test f2[1] == "HD 217357, \"GJ 884\"" && f2[5] == "" && f2[8] == "0"
        # names alone: unresolved targets
        write(joinpath(tmp, "names.txt"), "Altair\n# a comment\nDeneb\n")
        nm = split(GUI.shell_load_targets(joinpath(tmp, "names.txt")), "\n")
        @test [split(r, "\t")[1] for r in nm] == ["Altair", "Deneb"]
        @test all(r -> split(r, "\t")[8] == "0", nm)
        @test startswith(GUI.shell_load_targets(joinpath(tmp, "none.csv")), "!")
    end

    @testset "the facility's display order" begin
        @test telescope_display_order("CHARA") == [3, 6, 5, 2, 1, 4]       # E1 W2 W1 S2 S1 E2
        f = read_facility_file("VLTI_AT_small")
        @test telescope_display_order(f) == collect(1:f.ntel)              # none declared
    end
end
