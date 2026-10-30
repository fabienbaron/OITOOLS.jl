# Modeling › free-form limb darkening: α Cen A with the linear-law model, the law switched to
# free-form (5 nodes from the linear profile), the node count, the basis and the monotonic
# option changed from the Limb darkening box, a fit, the profile by variational inference and
# adopted, and back to a law.
#
#     julia --project=. test/gui/walks/ldfree.jl

include(joinpath(@__DIR__, "walk.jl"))
using .Walk

const DATA = joinpath(Walk.REPO, "demos", "data")
const MODELS = joinpath(Walk.REPO, "demos", "models")

shown(w, q) = try locate(w, q; timeout = 2).w > 0 catch; false end
num(w, q) = parse(Float64, value(w, q))

w = start("ldfree"; data_dir = DATA)
try
    stage(w, "α Cen A and the linear law") do
        click(w, "text:Open OIFITS…")
        click(w, "pickerFiles/text:AlphaCenA.oifits")
        click(w, "pickerOk"; pause = 1)
        waitfor(() -> startswith(locate(w, "datasetBox").current, "data1"), "the dataset"; timeout = 120)
        click(w, "text:Modeling"; pause = 1.5)
        click(w, "text:Open…"; pause = 1)
        pickfile(w, joinpath(MODELS, "Limb_Darkened_Disk_Linear.toml"))
        typein(w, "param:star,u:value", "0.4")
        check(w, "the law menu (1=yes)", locate(w, "ldLaw").current == "linear", 1; tol = 0)
        check(w, "no node settings for a law (1=yes)", !shown(w, "ldNodes"), 1; tol = 0)
    end
    stage(w, "free-form: 5 nodes from the linear profile") do
        choose(w, "ldLaw", "free-form")
        check(w, "5 nodes (1=yes)", value(w, "ldNodes") == "5", 1; tol = 0)
        check(w, "spline (1=yes)", locate(w, "ldBasis").current == "spline", 1; tol = 0)
        check(w, "not monotonic (1=yes)", !locate(w, "ldMonotonic").checked, 1; tol = 0)
        check(w, "the diameter kept", num(w, "param:star,ldspline:value"), 8.0; tol = 1e-6)
        check(w, "I at μ = 0 is 1 − u", num(w, "param:star,ld1:value"), 0.6; tol = 1e-6)
        check(w, "I at μ = 0.75", num(w, "param:star,ld4:value"), 0.9; tol = 1e-6)
        check(w, "no fifth free node (1=yes)", !shown(w, "param:star,ld5:value"), 1; tol = 0)
        shot(w, "01-free-form")
    end
    stage(w, "nodes, basis, monotonic") do
        spin(w, "ldNodes", 6); sleep(1)
        check(w, "6 nodes: ld5 (1=yes)", shown(w, "param:star,ld5:value"), 1; tol = 0)
        check(w, "I at μ = 0.8", num(w, "param:star,ld5:value"), 0.92; tol = 1e-6)
        choose(w, "ldBasis", "straight lines")
        check(w, "the piecewise disc (1=yes)", shown(w, "param:star,ldpiecewise:value"), 1; tol = 0)
        tick(w, "ldMonotonic", true); sleep(1)
        check(w, "ratios (1=yes)", shown(w, "param:star,ldq1:value"), 1; tol = 0)
        check(w, "the first ratio", num(w, "param:star,ldq1:value"), 0.6 / 0.68; tol = 1e-4)
        check(w, "still 6 nodes (1=yes)", value(w, "ldNodes") == "6", 1; tol = 0)
        shot(w, "02-monotonic")
    end
    stage(w, "a fit") do
        click(w, "text:Model Optimization"; pause = 1)
        fitnow(w)
        chi2r, p = fitresult(w, 1)
        check(w, "fitted (1=yes)", isfinite(chi2r), 1; tol = 0)
        check(w, "the ratios in [0, 1] (1=yes)",
              all(0 <= p["star,ldq$k"][1] <= 1 for k in 1:5), 1; tol = 0)
        shot(w, "03-fit")
    end
    stage(w, "the free-form profile by variational inference") do
        check(w, "the disc offered (1=yes)", locate(w, "ldviDisc").current == "star", 1; tol = 0)
        spin(w, "ldviNodes", 12)
        spin(w, "ldviDraws", 10)
        spin(w, "ldviBoot", 4)
        click(w, "ldviRun"; pause = 2)
        waitfor(() -> value(w, "ldviRun") == "Run", "the VI run"; timeout = 600)
        check(w, "the view shown (1=yes)", startswith(value(w, "ldviText"), "star: θ = "), 1; tol = 0)
        check(w, "bagged (1=yes)", endswith(value(w, "ldviText"), "bagged × 4"), 1; tol = 0)
        shot(w, "04-vi")
        choose(w, "ldviNorm", "⟨I⟩ = 1"); sleep(1)
        shot(w, "05-vi-mean")
        click(w, "ldviAdopt"; pause = 1)
        click(w, "text:Model Creation"; pause = 1)
        check(w, "12 nodes adopted (1=yes)", value(w, "ldNodes") == "12", 1; tol = 0)
        check(w, "straight lines (1=yes)", locate(w, "ldBasis").current == "straight lines", 1; tol = 0)
        check(w, "nodes fixed (1=yes)", locate(w, "param:star,ld1:mode").current == "fixed", 1; tol = 0)
    end
    stage(w, "back to a law") do
        choose(w, "ldLaw", "quadratic")
        check(w, "no node settings (1=yes)", !shown(w, "ldNodes"), 1; tol = 0)
        check(w, "no nodes left (1=yes)", !shown(w, "param:star,ld1:value") &&
                                           !shown(w, "param:star,ldq1:value"), 1; tol = 0)
        check(w, "the quadratic law (1=yes)", shown(w, "param:star,w:value"), 1; tol = 0)
        shot(w, "06-law")
    end
catch e
    println("walk stopped: ", sprint(showerror, e))
finally
    finish(w)
end
