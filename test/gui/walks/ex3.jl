# EX3 in the window — AX Cir's faint companion, PMOIRED's PIONIER file: a uniform disc, then the
# binary, Companion search… with the dialog's defaults (refinement of the best 5), Adopt, Fit,
# a bootstrap. V² and T3φ, smearing auto. Numbers against the script's (GUI_TODO, "EX3
# reproduced") within 0.05σ — the fit stops when a step would gain less than 0.001·χ²r, about
# 0.03σ per parameter, so the window's start (the search's minimum) sets the last digits — and
# PMOIRED's.
#
#     julia test/gui/walks/ex3.jl

include(joinpath(@__DIR__, "walk.jl"))
using .Walk

const DATA = joinpath(homedir(), "SOFTWARE", "PMOIRED_examples", "DATA", "AXCir")
const MODELS = joinpath(Walk.REPO, "demos", "models")
isfile(joinpath(DATA, "AXCir.oifits")) || error("needs the PMOIRED example data in $DATA")

# within `nσ` of the script's value, σ its printed uncertainty
near(w, label, got, want, σ; nσ = 0.05) = check(w, label, got, want; tol = nσ * σ)

w = start("ex3n"; data_dir = DATA)
try
    stage(w, "open AXCir.oifits") do
        click(w, "text:Open OIFITS…")
        click(w, "pickerFiles/text:AXCir.oifits")
        click(w, "pickerOk"; pause = 1)
        waitfor(() -> startswith(locate(w, "datasetBox").current, "data1"), "the dataset"; timeout = 120)
        shot(w, "01-open")
    end
    stage(w, "uniform disc, V² + T3φ, smearing auto") do
        click(w, "text:Modeling"; pause = 1.5)
        click(w, "text:Open…"; pause = 1)
        pickfile(w, joinpath(MODELS, "Uniform_Disk.toml"))
        typein(w, "param:star,ud:value", "0.8")
        click(w, "text:Model Optimization"; pause = 1)
        tick(w, "t3amp", false)
        choose(w, "smearing", "auto")
        choose(w, "param:star,ud:mode", "free")
        fitnow(w)
        chi2r, p = fitresult(w, 1)
        near(w, "UD θ", p["star,ud"][1], 0.93074, 0.0061)
        check(w, "UD σθ", p["star,ud"][2], 0.0061; tol = 1e-4)
        check(w, "UD χ²r", chi2r, 0.9753; tol = 1e-4)
        shot(w, "02-ud")
    end
    stage(w, "binary, companion search with the defaults") do
        click(w, "text:Open…"; pause = 1)
        pickfile(w, joinpath(MODELS, "AX_Cir_Binary.toml"))
        click(w, "companionButton"; pause = 1.5)
        click(w, "text:Search"; pause = 2)
        waitfor(() -> value(w, "fitButton") == "Fit", "the search"; timeout = 600)
        sleep(1)
        _, p = fitresult(w, 2)
        # refined best against PMOIRED's 6.228 ± 0.058, −28.487 ± 0.079
        check(w, "search x", p["companion,x"][1], 6.228; tol = 0.058)
        check(w, "search y", p["companion,y"][1], -28.487; tol = 0.079)
        shot(w, "03-search")
    end
    stage(w, "adopt, fit") do
        click(w, "text:Adopt"; pause = 1.5)
        fitnow(w)
        chi2r, p = fitresult(w, 3)
        check(w, "binary χ²r", chi2r, 0.7305; tol = 1e-4)
        near(w, "star,ud", p["star,ud"][1], 0.82478, 0.0078)
        check(w, "σ star,ud", p["star,ud"][2], 0.00782; tol = 1e-5)
        near(w, "companion,f", p["companion,f"][1], 0.0086221, 0.0004)
        near(w, "companion,x", p["companion,x"][1], 6.2250, 0.059)
        near(w, "companion,y", p["companion,y"][1], -28.4939, 0.079)
        shot(w, "04-fit")
    end
    stage(w, "bootstrap, PMOIRED's resampling") do
        click(w, "text:Adopt"; pause = 1.5)
        choose(w, "uncertainty", "bootstrap")
        choose(w, "bootMode", "pmoired")
        spin(w, "nboot", 200)
        fitnow(w; timeout = 900)
        _, p = fitresult(w, 4)
        # PMOIRED: ud 0.0080, f 0.00036, x 0.071, y 0.067 — σ compared within a factor 2
        for (k, s) in (("star,ud", 0.0080), ("companion,f", 0.00036), ("companion,x", 0.071),
                       ("companion,y", 0.067))
            check(w, "bootstrap σ $k (×½–×2 of PMOIRED's)", p[k][2] / s, 1.0; tol = 1.0)
        end
        shot(w, "05-bootstrap")
    end
catch e
    println("walk stopped: ", sprint(showerror, e))
finally
    exit(finish(w) ? 0 : 1)
end
