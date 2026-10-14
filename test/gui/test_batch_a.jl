# What the window does with features the core already has: averaging channels from Explore's
# Filter, merging datasets, the cube viewer and a FITS file's keywords (Tools), and fitting each
# of several datasets (Modeling).

@testset "Filter averaging, merge, cube viewer, fit each" begin
    G = GUI
    prev = G.SHELL[]
    sess = Session()
    sfig = Makie.Figure(); sax = Makie.Axis(sfig[1, 1])
    sh = G.ShellState(sess, sfig, sax, nothing, String[], Any[], 0,
                      :v2, :baseline, false, :none, :none, false, "", String[],
                      G.build_canvas(sfig, sax), "", nothing, nothing, nothing, nothing,
                      nothing, nothing, nothing, Any[], nothing, nothing, nothing,
                      nothing, "", "", nothing)
    G.SHELL[] = sh
    prevcube = G.CUBE_CANVAS[]
    cfig = Makie.Figure(); cax = Makie.Axis(cfig[1, 1])
    G.CUBE_CANVAS[] = G.build_canvas(cfig, cax)
    dir = mktempdir()
    replay(m) = (mod = Module(m); Core.eval(mod, :(using OITOOLS));
                 Core.eval(mod, Meta.parseall(export_script(sess))); mod)
    try
        @testset "Explore › Filter › Spectral: average every N channels" begin
            k = G.filter_spec_kwargs("binning\t10\nbincorr\t2")
            @test k.read == Dict(:binning => 10, :bincorr => 2.0)
            @test G.filter_spec_kwargs("binning\t10\nbinmode\tpmoired").read ==
                  Dict(:binning => 10, :binmode => :pmoired)
            @test isempty(G.filter_spec_kwargs("binning\t1").read)
            @test startswith(G.filter_spec_kwargs("binning\tmany"), "!")
            G.shell_open(POLY)
            src = sess.datasets[end]
            @test occursin("averaging every 10", G.shell_filter_preview("binning\t10"))
            @test !startswith(G.shell_apply_filter("binning\t10"), "!")
            e = G.current_dataset(sh)
            nch(x) = length(unique(x.data[1, 1].v2_lam))
            @test nch(e) == cld(nch(src), 10) && e.readopts[:binning] == 10
            @test G.shell_filter_spec() == "binning\t10"
            m = replay(:ReplayBin)
            @test Core.eval(m, :($(Symbol(e.name))[1, 1].nv2)) == e.data[1, 1].nv2
        end

        @testset "Tools › Spectral binning: group" begin
            @test G.group_windows([1.0, 2.0, 3.0, 4.0, 5.0], 2) == [[0.5, 2.5], [2.5, 4.5], [4.5, 5.5]]
            @test length(G.group_windows([1.0, 1.0, 2.0], 1)) == 2            # repeated channels once
            src = readoifits(POLY; warn = false, verbose = false)[1, 1]
            nch = length(unique(src.v2_lam))
            @test startswith(G.shell_tools_group(POLY, 10), "grouped every 10")
            g = G.current_dataset(sh)
            @test size(g.data, 1) == cld(nch, 10) && g.label == "grouped ×10"
            @test sum(d.nv2 for d in g.data) == src.nv2                          # every channel kept
            @test all(length(unique(d.v2_lam)) <= 10 for d in g.data)
            r = replay(:ReplayGroup)
            @test Core.eval(r, :(size($(Symbol(g.name))))) == size(g.data)
        end

        @testset "Tools › Merge files" begin
            f1, f2 = joinpath(OIF, "AZ_Cyg_2011.oifits"), joinpath(OIF, "AZ_Cyg_2014.oifits")
            a = readoifits(f1; warn = false, verbose = false)[1, 1]
            b = readoifits(f2; warn = false, verbose = false)[1, 1]
            msg = G.shell_tools_merge(f1 * "\n" * f2)
            @test startswith(msg, "merged 2 files")
            m = G.current_dataset(sh)
            @test m.files == [f1, f2] && m.data[1, 1].nv2 == a.nv2 + b.nv2 && m.label == "merged"
            r = replay(:ReplayMerge)
            @test Core.eval(r, :($(Symbol(m.name))[1, 1].nv2)) == m.data[1, 1].nv2
            @test startswith(G.shell_tools_merge(f1), "!")
            @test startswith(G.shell_tools_merge(f1 * "\n" * f1), "!")
        end

        @testset "Tools › Cube viewer" begin
            cube = Float64[x + 10y + 100w for x in 1:8, y in 1:8, w in 1:3]
            p = joinpath(dir, "cube.fits")
            writefits(cube, p; pixsize = 0.5, wavelengths = [1.6e-6, 1.7e-6, 1.8e-6])
            labels = split(G.shell_cube_open(p), "\n")
            @test length(labels) == 3 && occursin("1.7 µm", labels[2])
            @test G.shell_cube_show(2) == labels[2]
            m = G.CUBE[]
            @test Float32.(m.image[:, :, 2, 1]) == G.CUBE_CANVAS[].imagedata[]
            @test G.shell_cube_show(99) == labels[3]            # clamped to the last plane
            img = joinpath(dir, "img.fits"); writefits(rand(6, 6), img; pixsize = 0.3)
            @test G.shell_cube_open(img) == "plane 1"
            @test startswith(G.shell_cube_open(joinpath(OIF, "AX_Cir.oifits")), "!")
            # a third axis without CTYPE3: wavelength from CUNIT3 / CRVAL3 / CDELT3, else unlabelled
            q = joinpath(dir, "nowcs.fits")
            OITOOLS.AstroFITS.writefits!(q, ["CUNIT3" => "m", "CRVAL3" => 2.0e-6, "CDELT3" => 1.0e-9,
                                             "CRPIX3" => 1.0], rand(8, 8, 4))
            lq = split(G.shell_cube_open(q), "\n")
            @test length(lq) == 4 && occursin("2.003 µm", lq[4])
            r = joinpath(dir, "bare.fits"); OITOOLS.AstroFITS.writefits!(r, ["OBJECT" => "x"], rand(8, 8, 3))
            @test split(G.shell_cube_open(r), "\n") == ["plane 1", "plane 2", "plane 3"]
        end

        @testset "Tools › Set keywords, Save FITS" begin
            # a cube with no axis keywords: 3 wavelength channels × 2 epochs
            cube = Float32[hypot(x - 8 - t, y - 8) < 5 ? w + t : 0.01 for x in 1:16, y in 1:16, w in 1:3, t in 1:2]
            bare = joinpath(dir, "bare.fits")
            OITOOLS.AstroFITS.writefits!(bare, ["OBJECT" => "disc"], cube)
            kv = Dict(split(l, "\t"; limit = 2) for l in split(G.shell_fits_keywords(bare), "\n"))
            @test kv["dims"] == "16 16 3 2" && kv["pixsize"] == "" && kv["eastleft"] == "0"
            @test kv["axis3"] == "wave" && kv["wstart"] == ""
            cards = [split(l, "\t") for l in split(G.shell_fits_keywords(bare), "\n") if startswith(l, "card")]
            @test length(cards) == 20 && all(c -> c[3] == "0", cards)     # 4 axes × 5 cards, none there
            # no CTYPE3, but CUNIT3 / CRVAL3 / CDELT3 in microns: the wavelengths are prefilled
            amhra = joinpath(dir, "amhra.fits")
            OITOOLS.AstroFITS.writefits!(amhra, ["CTYPE1" => "RA---TAN", "CDELT1" => -2.7e-7, "CDELT2" => 2.7e-7,
                                                 "CRVAL3" => 1.0, "CDELT3" => 1.5, "CUNIT3" => "MICRON"],
                                         rand(8, 8, 5))
            ka = Dict(split(l, "\t"; limit = 2) for l in split(G.shell_fits_keywords(amhra), "\n"))
            @test ka["wstart"] == "1" && ka["wstep"] == "1.5" && ka["eastleft"] == "1"
            ca = Dict(c[2] => c[3] for c in (split(l, "\t") for l in split(G.shell_fits_keywords(amhra), "\n"))
                      if c[1] == "card")
            @test ca["CTYPE1"] == "1" && ca["CDELT3"] == "1" && ca["CTYPE3"] == "0" && ca["CRPIX3"] == "0"
            G.FITSKW[] = nothing                                        # nothing set yet
            @test startswith(G.shell_fits_save(joinpath(dir, "z.fits")), "!")
            # Apply: values in the dialog's units; problems named
            @test startswith(G.shell_fits_apply(bare, "", "mas", false, "wave", "1.5", "0.1", "µm",
                                                "58000.5", "2", "d"), "! the pixel size")
            @test startswith(G.shell_fits_apply(bare, "0.25", "mas", false, "wave", "1.5", "", "µm",
                                                "58000.5", "2", "d"), "! the wavelength step")
            labs = split(G.shell_fits_apply(bare, "250", "µas", false, "wave", "1.5", "0.1", "µm",
                                            "58000.5", "48", "h"), "\n")
            @test length(labs) == 6 && occursin("1.5 µm", labs[1]) && G.shell_fits_keyword_source() == bare
            # the viewer's two axes, and a channel of an epoch
            ax = split(G.shell_cube_axes(), "\n")
            @test ax[1] == "3\t2" && split(ax[2], "\t") == ["λ=1.5 µm", "λ=1.6 µm", "λ=1.7 µm"]
            @test split(ax[3], "\t") == ["MJD 58000.5000", "MJD 58002.5000"]
            @test G.shell_cube_show_at(2, 2) == "λ 1.6 µm · MJD 58002.5000"
            @test G.CUBE_CANVAS[].imagedata[] == Float32.(cube[:, :, 2, 2])
            @test G.shell_cube_show_at(9, 9) == "λ 1.7 µm · MJD 58002.5000"   # clamped
            out = joinpath(dir, "keys.fits")
            @test G.shell_fits_save(out) == "wrote keys.fits"
            m = read_model_image(out)
            @test m.wavelengths ≈ [1.5e-6, 1.6e-6, 1.7e-6] && m.mjds ≈ [58000.5, 58002.5]
            @test m.pixsize ≈ 0.25 && m.image ≈ cube ./ sum(cube; dims = (1, 2))
            # the logged line writes the same file
            line = last(filter(l -> occursin("fits_set_keywords(", l), sh.console))
            line = line[findfirst("fits_set_keywords(", line)[1]:end]
            again = joinpath(dir, "again.fits")
            mod = Module(:ReplayKeywords); Core.eval(mod, :(using OITOOLS))
            Core.eval(mod, Meta.parse(replace(line, out => again)))
            @test read_model_image(again).wavelengths ≈ m.wavelengths && read_model_image(again).mjds ≈ m.mjds
            # the saved file's own keywords prefill the dialog
            kv = Dict(split(l, "\t"; limit = 2) for l in split(G.shell_fits_keywords(out), "\n"))
            @test parse(Float64, kv["pixsize"]) ≈ 0.25 && parse(Float64, kv["wstart"]) ≈ 1.5
            @test parse(Float64, kv["tstep"]) ≈ 2.0
        end

        @testset "Modeling › SPARCO fluxes" begin
            p = joinpath(dir, "two.toml")
            write(p, """
                free = ["disk,fwhm"]
                [model]
                "star,ud" = 0.3
                "star,f" = 1.0
                "disk,fwhm" = 3.0
                "disk,f" = 1.0
                """)
            G.shell_open_model(p)
            @test startswith(G.shell_sparco_fluxes("star", "star", -4, 1.6), "!")
            @test startswith(G.shell_sparco_fluxes("star", "halo", -4, 1.6), "!")
            @test startswith(G.shell_sparco_fluxes("star", "disk", -4, 0), "!")
            @test G.shell_sparco_fluxes("star", "disk", -4, 1.6) == ""
            m = G._model(sh)
            @test m.dict["star,f"] == "\$fs * (\$WL/1.6)^(-4) / (\$fs * (\$WL/1.6)^(-4) + (1 - \$fs) * (\$WL/1.6)^\$disk,slope)"
            @test occursin("(1 - \$fs) * (\$WL/1.6)^\$disk,slope /", m.dict["disk,f"])
            @test m.dict["fs"] == 0.5 && m.dict["disk,slope"] == 0.0
            @test issubset(["fs", "disk,slope", "disk,fwhm"], m.free) && !("star,f" in m.free)
            @test (m.lb["fs"], m.ub["fs"], m.lb["disk,slope"], m.ub["disk,slope"]) == (0.0, 1.0, -4.0, 4.0)
            # the fluxes sum to 1 at any wavelength, whatever fs and the slope
            fm = dict_to_model(merge(copy(m.dict), Dict{String,Any}("fs" => 0.3, "disk,slope" => 1.7)), String[])
            for λ in (1.5e-6, 1.6e-6, 2.2e-6)
                @test real(eval_model(fm, Float64[], zeros(2, 1); wl = [λ])[1]) ≈ 1.0
            end
            # applied again: fs and the slope keep their values
            m.dict["fs"] = 0.42
            @test G.shell_sparco_fluxes("star", "disk", -4, 1.6) == "" && m.dict["fs"] == 0.42
        end

        @testset "Modeling › orbit and radial velocities" begin
            prevorb = G.ORBIT_CANVAS[]
            ofig = Makie.Figure()
            G.ORBIT_CANVAS[] = G.build_orbit_view(ofig, Makie.Axis(ofig[1, 1]), Makie.Axis(ofig[1, 2]))
            # a binary on an orbit, observed through V1295 Aql's uv points, with velocities
            truth = Dict{String,Any}("1,ud" => 0.3, "2,f" => 0.3, "2,x" => "orbit", "2,y" => "orbit",
                                     "2,orb P" => 30.0, "2,orb MJD0" => 58640.0, "2,orb e" => 0.4,
                                     "2,orb omega" => 60.0, "2,orb OMEGA" => 120.0, "2,orb incl" => 50.0,
                                     "2,orb a" => 4.0, "2,orb Ka" => 40.0, "2,orb Kb" => 55.0,
                                     "2,orb gamma" => -12.0)
            oif = joinpath(dir, "binary.oifits")
            simulate_from_oifits(joinpath(OIF, "V1295_Aql.oifits"), oif; mode = "copy_errors",
                                 noise = false, flat_model = dict_to_model(truth, String[]), flat_params = Float64[])
            rvf = joinpath(dir, "rv.csv")
            open(rvf, "w") do io
                println(io, "MJD,RV,sigma,star")
                for t in 58600.0:9.0:58700.0, s in (:a, :b)
                    println(io, t, ",", orbit_rv(t, s, 30.0, 58640.0, 0.4, 60.0, s === :a ? 40.0 : 55.0, -12.0),
                            ",0.1,", s)
                end
            end
            G.shell_open(oif)
            mf = joinpath(dir, "binary.toml")
            write(mf, """
                [model]
                "1,ud" = 0.3
                "2,f" = 0.3
                """)
            G.shell_open_model(mf)
            # Geometry › orbit: x = y = "orbit", the elements as rows, MJD0 at the data's mean date
            @test G.shell_component_geometry("2") == "0\t0\t0"
            @test G.shell_set_component_geometry("2", "orbit", true) == ""
            md = G._model(sh).dict
            @test md["2,x"] == "orbit" && haskey(md, "2,orb P") && 58600 < md["2,orb MJD0"] < 58700
            @test G.shell_component_geometry("2") == "1\t0\t1" && G.shell_orbit_component() == "2"
            for (k, v) in truth; startswith(k, "2,orb ") && k ∉ ("2,orb Ka", "2,orb Kb", "2,orb gamma") && (md[k] = v); end
            # the RV table adds the velocity elements to the orbit's component
            @test startswith(G.shell_rv_open(rvf), "rv.csv: 24 velocities (12 a, 12 b)")
            @test md["2,orb Ka"] == 10.0 && md["2,orb gamma"] == 0.0
            @test startswith(G.shell_rv_open(joinpath(dir, "nope.csv")), "!")
            G.shell_rv_open(rvf)
            # the fit, with the velocities
            free = "2,orb Ka\t0\t100\n2,orb Kb\t0\t100\n2,orb gamma\t-50\t50"
            flags = (true, true, true, false, false, false, false)
            fitargs(opt) = (free, "", "", flags..., opt, 100, "", "", 60, 0, "jacobian", 200,
                            "replacement", "config", 1, true, true)          # the last: fit RVs
            @test startswith(G.shell_fit_model(fitargs("grid")...), "!")
            @test G.shell_fit_model(fitargs("lsqfit")...) == "fitting…"
            line = ""
            for _ in 1:600
                f = split(G.shell_job_poll(), "\t"; limit = 2)
                f[1] == "done" && (line = f[2]; break)
                sleep(0.1)
            end
            r = sh.fits[end].result
            p = Dict(zip(r.list_free_params, r.x_opt))
            @test isapprox(p["2,orb Ka"], 40.0; rtol = 1e-3) && isapprox(p["2,orb gamma"], -12.0; atol = 0.01)
            @test any(l -> occursin("radial velocities", l), sh.console)
            # the Orbit view: the track, the nights, both stars' velocities
            @test startswith(G.shell_orbit_view(), "2: P 30 d, e 0.400")
            ov = G.ORBIT_CANVAS[]
            @test length(ov.track[]) == 400 && !isempty(ov.epochs[])
            @test length(ov.dots[1][]) == 12 && length(ov.dots[2][]) == 12 && length(ov.curve[1][]) == 300
            # off its orbit
            @test G.shell_set_component_geometry("2", "orbit", false) == ""
            @test !haskey(md, "2,orb P") && md["2,x"] == 0.0 && startswith(G.shell_orbit_view(), "!")
            G.shell_rv_close(); G.ORBIT_CANVAS[] = prevorb
        end

        @testset "Modeling › Fit each" begin
            G.shell_tools_split(joinpath(OIF, "V1295_Aql.oifits"), "gap", 0.5, "")
            parts = sess.datasets[end-1:end]
            G.shell_open_model(joinpath(pkgdir(OITOOLS), "demos", "models", "Uniform_Disk.toml"))
            G.shell_set_param("star,ud", "value", 0.5)
            free = "star,ud\t0.05\t3"
            flags = (true, false, true, false, false, false, false)
            names = join([e.name for e in parts], "\n")
            @test startswith(G.shell_fit_each("", free, "", "", flags..., "lsqfit", 100, "jacobian"), "!")
            @test startswith(G.shell_fit_each(names, free, "", "", flags..., "grid", 100, "jacobian"), "!")
            @test startswith(G.shell_fit_each(names, free, "", "", flags..., "lsqfit", 100, "jacobian"), "fitting 2")
            line = ""
            for _ in 1:600
                f = split(G.shell_job_poll(), "\t"; limit = 2)
                f[1] == "done" && (line = f[2]; break)
                sleep(0.1)
            end
            @test line == "fitted 2 of 2 datasets"
            t = split(G.shell_epochs_table(), "\n")
            @test t[1] == "dataset\tMJD\tchi2r\tndof\tstar,ud\tσ star,ud" && length(t) == 3
            md = G._model(sh).dict
            for (row, e) in zip(t[2:end], parts)
                c = split(row, "\t")
                @test c[1] == e.name
                r = fit_model_lsqfit(md, ["star,ud"], e.data[1, 1];
                                     lb = Dict("star,ud" => 0.05), ub = Dict("star,ud" => 3.0),
                                     weights = G.fitting_weights(; t3amp = false), maxIter = 100)
                @test parse(Float64, c[5]) ≈ r.x_opt[1] rtol = 1e-4
                @test parse(Float64, c[2]) ≈ e.data[1, 1].mean_mjd atol = 1e-4
            end
            out = joinpath(dir, "epochs.tsv")
            @test startswith(G.shell_save_epochs(out), "saved") && readlines(out)[1] == t[1]
        end
    finally
        G.SHELL[] = prev
        G.CUBE_CANVAS[] = prevcube
    end
end
