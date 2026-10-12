# What the window does with features the core already has: averaging channels from Explore's
# Filter, merging datasets and the cube viewer (Tools), and fitting each of several datasets
# (Modeling).

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
