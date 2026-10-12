# The Tools perspective's data layer: headers, validation, the repaired copy, the splits and
# channel averaging, through the shell calls the tab makes. A split or an average is a logged
# re-read, so the exported script must rebuild every part.

@testset "Tools data layer" begin
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
    dir = mktempdir()
    v1295 = joinpath(OIF, "V1295_Aql.oifits")
    try
        @testset "headers" begin
            lines = split(G.shell_tools_headers(v1295), "\n")
            @test length(lines) == length(oifits_headers(v1295))
            @test startswith(lines[1], "1\tprimary")
            @test any(l -> occursin("\tOI_WAVELENGTH\tMIRCX\t", l), lines)
            txt = G.shell_tools_hdu(v1295, 2)
            @test occursin("Columns (", txt) && occursin("Header (", txt) && occursin("XTENSION", txt)
            @test startswith(G.shell_tools_hdu(v1295, 999), "!")
            @test occursin("oifits_headers(", G.shell_console())
        end

        @testset "validate and repair" begin
            r = G.shell_tools_validate(MONO)
            @test any(l -> startswith(l, "warning\t") && occursin("not in OI_ARRAY", l), split(r, "\n"))
            pol = joinpath(OIF, "polaris.oifits")
            @test any(l -> startswith(l, "warning\t") && occursin("revision 1", l),
                      split(G.shell_tools_validate(pol), "\n"))
            out = joinpath(dir, "polaris_repaired.oifits")
            @test startswith(G.shell_tools_repair(pol, out, ""), "wrote")
            @test !any(l -> startswith(l, "error\t") || startswith(l, "warning\t"),
                       split(G.shell_tools_validate(out), "\n"))
            @test startswith(G.shell_tools_repair(pol, pol, ""), "!")
            mwc = joinpath(OIF, "MWC275_4T.oifits")
            @test split(G.shell_tools_targets(mwc), "\n") == ["MWC 275", "HD163296"]
            @test startswith(G.shell_tools_repair(mwc, joinpath(dir, "m.oifits"), ""), "!")
            @test startswith(G.shell_tools_repair(mwc, joinpath(dir, "m.oifits"), "MWC 275"), "wrote")
        end

        @testset "split" begin
            @test startswith(G.shell_tools_split("gap", 0.5, ""), "!")       # nothing loaded
            @test !startswith(G.shell_open(v1295), "!")
            src = sess.datasets[end]; d = src.data[1, 1]
            @test length(split(G.shell_tools_split_preview("instrument", 0), "\n")) == 7
            pv = split(G.shell_tools_split_preview("gap", 0.5), "\n")
            @test length(pv) == 6
            @test sum(parse(Int, match(r"\t(\d+) V²", l)[1]) for l in pv) == d.nv2
            n0 = length(sess.datasets)
            msg = G.shell_tools_split("gap", 0.5, dir)
            @test occursin("into 6 datasets, 6 written", msg)
            parts = sess.datasets[n0+1:end]
            @test length(parts) == 6 && all(e -> e.derived_from === nothing && !isempty(e.label), parts)
            @test sum(e.data[1, 1].nv2 for e in parts) == d.nv2
            @test sum(length(e.data[1, 1].t3phi) for e in parts) == length(d.t3phi)
            @test length(filter(endswith(".oifits"), readdir(dir))) >= 6
            @test G.current_dataset(sh) === parts[end]
            @test occursin(parts[1].label, G.shell_dataset_names())
            # by instrument: each part holds its instrument's points only
            G.shell_select_dataset(n0)
            G.shell_tools_split("instrument", 0, "")
            ins = sess.datasets[end-6:end]
            @test [e.label for e in ins][1] == "MIRCX"
            @test sum(e.data[1, 1].nv2 for e in ins) == d.nv2
            @test all(all(t -> x.data[1, 1].tables[t].insname == x.label, x.data[1, 1].v2_tab) for x in ins)
            # a filtered dataset is not re-read: its cut would be lost
            G.shell_select_dataset(n0)
            @test !startswith(G.shell_apply_filter("bl_max\t200"), "!")
            @test startswith(G.shell_tools_split("gap", 0.5, ""), "!")
        end

        @testset "average channels" begin
            G.shell_open(POLY)
            src = sess.datasets[end]
            msg = G.shell_tools_bin(10, "block", 1)
            @test startswith(msg, "averaged ×10") && occursin("scatter", msg)
            b = sess.datasets[end]
            @test b.label == "binned ×10" && b.readopts[:binning] == 10
            nch(e) = length(unique(e.data[1, 1].v2_lam))
            @test nch(b) == cld(nch(src), 10)
            @test startswith(G.shell_tools_bin(4, "pmoired", 1), "averaged ×4")
            @test sess.datasets[end].readopts[:binmode] == :pmoired
            @test startswith(G.shell_tools_bin(1, "block", 1), "!")
        end

        @testset "the script rebuilds every part" begin
            m = Module(:ReplayTools)
            Core.eval(m, :(using OITOOLS))
            Core.eval(m, Meta.parseall(export_script(sess)))
            for e in sess.datasets
                isempty(e.label) && continue
                @test Core.eval(m, :($(Symbol(e.name))[1, 1].nv2)) == e.data[1, 1].nv2
            end
        end
    finally
        G.SHELL[] = prev
    end
end
