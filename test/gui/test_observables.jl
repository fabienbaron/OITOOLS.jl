# Choosing observables: Explore's Filter keeps or leaves out each
# observable of the DATASET, one half of an amplitude/phase pair through `filter_data`'s `drop`;
# Modeling chooses again, per fit, among what is kept, with VIS amplitude and phase as two boxes.

@testset "Observables: Explore's drop, Modeling's visamp / visφ" begin
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
    has(t) = Set(l[5:end] for l in split(t, "\n") if startswith(l, "has\t"))
    flags(s) = Dict(split(kv, "=")[1] => split(kv, "=")[2] == "1" for kv in split(s, ","))
    try
        @testset "the spec" begin
            k = G.filter_spec_kwargs("drop\tt3amp,visphi")
            @test k.drop == [:t3amp, :visphi] && isempty(k.read) && isempty(k.filter)
            @test startswith(G.filter_spec_kwargs("drop\tclosure"), "!")
            @test G.filter_spec_kwargs("use_t3\t0").read == Dict(:use_t3 => false)
        end

        @testset "what the data hold" begin
            G.shell_open(joinpath(OIF, "Beauty_Contest_2026", "OBJECT2_K.oifits"))
            @test has(G.shell_filter_choices()) == Set(["v2", "t3phi", "visamp", "flux"])
            G.shell_open(joinpath(OIF, "MWC275_4T.oifits"))
            @test has(G.shell_filter_choices()) == Set(["v2", "t3amp", "t3phi", "visamp", "visphi"])
        end

        @testset "dropping one half" begin
            src = sess.datasets[end]; d = src.data[1, 1]
            @test !startswith(G.shell_apply_filter("drop\tt3amp"), "!")
            e = G.current_dataset(sh); f = e.data[1, 1]
            @test f.nt3amp == 0 && f.nt3phi == d.nt3phi && f.nv2 == d.nv2 && f.nvisphi == d.nvisphi
            @test G.shell_filter_spec() == "drop\tt3amp"
            o = flags(G.shell_observables())
            @test !o["t3amp"] && o["t3phi"] && o["visamp"] && o["visphi"]
            # the script rebuilds it, drop included
            @test occursin("drop = [:t3amp]", G.shell_script())
            m = Module(:ReplayDrop)
            Core.eval(m, :(using OITOOLS))
            Core.eval(m, Meta.parseall(export_script(sess)))
            @test Core.eval(m, :($(Symbol(e.name))[1, 1].nt3amp)) == 0
            @test Core.eval(m, :($(Symbol(e.name))[1, 1].nt3phi)) == f.nt3phi
            # both halves: the table is left out when read
            G.shell_select_dataset(findfirst(==(src), sess.datasets))
            @test !startswith(G.shell_apply_filter("use_t3\t0"), "!")
            @test isempty(G.current_dataset(sh).data[1, 1].t3phi)
        end

        @testset "Modeling: visamp and visφ apart" begin
            G.shell_select_dataset(findfirst(e -> isempty(e.filterspec) && e.path == joinpath(OIF, "MWC275_4T.oifits"),
                                             sess.datasets))
            @test G.fitting_weights(; visamp = true) == [1, 1, 1, 1, 0, 0, 0]
            @test G.fitting_weights(; v2 = false, t3amp = false, t3phi = false, visphi = true) == [0, 0, 0, 0, 1, 0, 0]
            G.shell_open_model(joinpath(pkgdir(OITOOLS), "demos", "models", "Uniform_Disk.toml"))
            fm = parse_model(G._model(sh).dict, G._model(sh).free)
            d = G.current_dataset(sh).data[1, 1]
            x = Float64[G._model(sh).dict[k] for k in G._model(sh).free]
            for w in ([0.0, 0, 0, 1, 0, 0, 0], [0.0, 0, 0, 0, 1, 0, 0])
                c = G.shell_model_chi2(false, false, false, w[4] > 0, w[5] > 0, false, false)
                @test c == GUI.Printf.@sprintf("%.4f", OITOOLS.chi2_flat(fm, x, d; weights = w) / OITOOLS._ndof(d, w))
            end
            @test G.shell_model_chi2(false, false, false, true, false, false, false) !=
                  G.shell_model_chi2(false, false, false, false, true, false, false)
        end
    finally
        G.SHELL[] = prev
    end
end
