# Observing › a list simulated, each target from its own sky: two models opened in Modeling, a
# list of four targets loaded, three models set through the Model popup (a Modeling model, a FITS
# image, a uniform disc for the calibrator), the fourth left without one (○, not simulated),
# Simulate into one file per target, the files and the calibrator's CATEGORY checked, a model
# renamed and the reference following, the list saved and read back, Simulate into one file,
# Open in Explore: one dataset per target.
#
#     julia --project=. test/gui/walks/simtargets.jl   (the checks read the files with OITOOLS)

include(joinpath(@__DIR__, "walk.jl"))
using .Walk
using OITOOLS

const MODELS = joinpath(Walk.REPO, "demos", "models")
const TMP = mktempdir()
const DISC = joinpath(TMP, "disc.fits")
writefits([hypot(i - 17, j - 17) < 4 ? 1.0 : 0.0 for i in 1:33, j in 1:33], DISC; pixsize = 0.1)
list_file = joinpath(TMP, "list.csv")
write(list_file, """
name,ra,dec,V,J,H,K,simbad,cal,for
A,279.2347,38.7837,,,3.0,,0,0,
B,280.0,39.0,,,4.0,,0,0,
C,281.0,37.0,,,4.0,,0,1,A
D,282.0,38.0,,,4.0,,0,0,
""")

shown(w, q) = try locate(w, q; timeout = 2).w > 0 catch; false end
status(w) = value(w, "obsStatus")

w = start("simtargets"; data_dir = TMP)
const OUT = joinpath(w.out, "night.oifits")
const BOTH = joinpath(w.out, "night_all.oifits")
file(n) = joinpath(w.out, "night_$n.oifits")
try
    stage(w, "two models in Modeling") do
        click(w, "text:Modeling"; pause = 1.5)
        click(w, "text:Open…"; pause = 1)
        pickfile(w, joinpath(MODELS, "Uniform_Disk.toml"))
        click(w, "text:Open…"; pause = 1)
        pickfile(w, joinpath(MODELS, "AX_Cir_Binary.toml"))
        check(w, "the binary selected (1=yes)", locate(w, "modelSelect").current == "AX_Cir_Binary", 1; tol = 0)
    end
    stage(w, "a list, three sources set") do
        click(w, "text:Observing"; pause = 2)
        click(w, "targetsLoad"; pause = 1)
        pickfile(w, list_file)
        check(w, "four targets (1=yes)", value(w, "obsTargetName:3") == "D", 1; tol = 0)
        check(w, "no model line yet (1=yes)", !shown(w, "obsTargetSimText:0"), 1; tol = 0)
        check(w, "A: ○ (1=yes)", value(w, "obsTargetModelMark:0") == "○", 1; tol = 0)
        check(w, "Simulate waits for a model (1=yes)", value(w, "simBlockReason") == "no target has a model", 1; tol = 0)
        typein(w, "obsDate", "2026-06-21")
        typein(w, "obsHaMin", "-1")
        typein(w, "obsHaMax", "1")

        click(w, "obsTargetSim:0"; pause = 1)
        click(w, "simDlgKind"; pause = 0.8)
        shot(w, "01a-popup")
        click(w, "popup:Uniform_Disk"; pause = 0.8)
        check(w, "the dialog finds it (1=yes)", value(w, "simDlgStatus") == "Uniform_Disk", 1; tol = 0)
        click(w, "simDlgOk"; pause = 1)
        check(w, "A: Uniform_Disk (1=yes)", value(w, "obsTargetSimText:0") == "Uniform_Disk", 1; tol = 0)
        check(w, "A: ✓ (1=yes)", value(w, "obsTargetModelMark:0") == "✓", 1; tol = 0)

        click(w, "obsTargetSim:1"; pause = 1)
        choose(w, "simDlgKind", "FITS image…")
        click(w, "simDlgFile"; pause = 1)
        pickfile(w, DISC)
        check(w, "the header's pixel size (1=yes)", value(w, "simDlgPix") == "0.1000", 1; tol = 0)
        click(w, "simDlgOk"; pause = 1)
        check(w, "B: disc.fits (1=yes)", value(w, "obsTargetSimText:1") == "disc.fits", 1; tol = 0)

        click(w, "obsTargetSim:2"; pause = 1)
        choose(w, "simDlgKind", "Uniform disc…")
        typein(w, "simDlgDiam", "0.6"; enter = "Tab")
        click(w, "simDlgOk"; pause = 1)
        check(w, "C: UD 0.6 mas (1=yes)", value(w, "obsTargetSimText:2") == "UD 0.6 mas", 1; tol = 0)
        check(w, "D: ○ (1=yes)", value(w, "obsTargetModelMark:3") == "○", 1; tol = 0)
        check(w, "3 of 4 targets (1=yes)", value(w, "simModelCount") == "3 of 4 targets", 1; tol = 0)
        shot(w, "01-sources")
    end
    stage(w, "Simulate: one file per target, D not simulated") do
        click(w, "text:Output file…"; pause = 1)
        typein(w, "pickerName", OUT; pause = 2)
        check(w, "the files named per target (1=yes)",
              value(w, "simOutputText") == joinpath(w.out, "night_<target>.oifits"), 1; tol = 0)
        click(w, "text:Simulate"; pause = 2)
        waitfor(() -> startswith(status(w), "3 simulated"), "the simulations"; timeout = 600)
        check(w, "three, no failure (1=yes)", status(w) == "3 simulated", 1; tol = 0)
        check(w, "three files (1=yes)", all(isfile, file.(["A", "B", "C"])) && !isfile(file("D")), 1; tol = 0)
        check(w, "C written as a calibrator (1=yes)",
              OITOOLS.load_oifits(file("C")).target.list[1].category == "CAL", 1; tol = 0)
        check(w, "A as a science target (1=yes)",
              OITOOLS.load_oifits(file("A")).target.list[1].category == "SCI", 1; tol = 0)
        shot(w, "02-simulated")
    end
    stage(w, "a model renamed: A follows it") do
        click(w, "text:Modeling"; pause = 1.5)
        l = locate(w, "modelSelect")
        Walk._x(w, "mousemove", "--window", w.wid, string(l.x + l.w - 8), string(l.y + l.h ÷ 2), "click", "1")
        sleep(0.5)
        click(w, "popup:Uniform_Disk"; pause = 1)
        typein(w, "modelSelect", "UD star"; pause = 1)
        click(w, "text:Observing"; pause = 2)
        check(w, "A: UD star (1=yes)", value(w, "obsTargetSimText:0") == "UD star", 1; tol = 0)
    end
    stage(w, "the list saved and read back with its sources") do
        click(w, "targetsSave"; pause = 1)
        typein(w, "pickerName", joinpath(TMP, "saved.csv"); pause = 2)
        click(w, "targetsLoad"; pause = 1)
        pickfile(w, joinpath(TMP, "saved.csv"))
        check(w, "A (1=yes)", value(w, "obsTargetSimText:0") == "UD star", 1; tol = 0)
        check(w, "B (1=yes)", value(w, "obsTargetSimText:1") == "disc.fits", 1; tol = 0)
        check(w, "C (1=yes)", value(w, "obsTargetSimText:2") == "UD 0.6 mas", 1; tol = 0)
        check(w, "D has none (1=yes)", !shown(w, "obsTargetSimText:3"), 1; tol = 0)
    end
    stage(w, "Simulate: every target in one file") do
        tick(w, "simFilePerTarget", false)
        click(w, "text:Output file…"; pause = 1)
        typein(w, "pickerName", BOTH; pause = 2)
        check(w, "one file (1=yes)", value(w, "simOutputText") == BOTH, 1; tol = 0)
        click(w, "text:Simulate"; pause = 2)
        waitfor(() -> startswith(status(w), "3 simulated"), "the simulations"; timeout = 600)
        check(w, "into night_all (1=yes)", status(w) == "3 simulated into " * BOTH, 1; tol = 0)
        check(w, "three targets in it (1=yes)", list_oifits_targets(BOTH) == ["A", "B", "C"], 1; tol = 0)
        check(w, "C a calibrator (1=yes)",
              [t.category for t in OITOOLS.load_oifits(BOTH).target.list] == ["SCI", "SCI", "CAL"], 1; tol = 0)
        shot(w, "03-one-file")
    end
    stage(w, "Open in Explore: one dataset per target") do
        click(w, "text:Open in Explore"; pause = 3)
        waitfor(() -> startswith(locate(w, "datasetBox").current, "data"), "the dataset"; timeout = 120)
        check(w, "A selected (1=yes)", endswith(locate(w, "datasetBox").current, "night_all.oifits · A"), 1; tol = 0)
        click(w, "datasetBox"; pause = 0.5)
        check(w, "C listed (1=yes)", shown(w, "popup:data3 — night_all.oifits · C"), 1; tol = 0)
        press(w, "Escape")
        shot(w, "04-explore")
    end
catch e
    println("walk stopped: ", sprint(showerror, e))
finally
    finish(w)
end
