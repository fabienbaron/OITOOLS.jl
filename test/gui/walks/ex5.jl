# EX5 in the window, companion half — o Leo, one GRAVITY file: instrument GRAVITY_SC averaged
# ×20, |V| and T3φ, Companion search… over ±12 mas in 1.5 mas steps, the sizes freed and
# refitted, a bootstrap. Numbers against PMOIRED's printed results (the notebook).
#
#     julia test/gui/walks/ex5.jl

include(joinpath(@__DIR__, "walk.jl"))
using .Walk

const DATA = joinpath(homedir(), "SOFTWARE", "PMOIRED_examples", "DATA", "o_Leo")
const FILE = "GRAVI.2022-02-28T042658.272_singlescivis_singlesciviscalibrated.fits"
isfile(joinpath(DATA, FILE)) || error("needs the PMOIRED example data in $DATA")

# within `nσ` of PMOIRED's σ
near(w, label, got, want, σ; nσ = 1.0) = check(w, label, got, want; tol = nσ * σ)

w = start("ex5n"; data_dir = DATA)
try
    stage(w, "open, GRAVITY_SC, ×20 averaging") do
        click(w, "text:Open OIFITS…")
        click(w, "pickerFiles/text:" * FILE)
        click(w, "pickerOk"; pause = 1)
        waitfor(() -> startswith(locate(w, "datasetBox").current, "data1"), "the dataset"; timeout = 180)
        click(w, "text:Tables"; pause = 1)
        choose(w, "filterInstrument", "GRAVITY_SC")
        click(w, "text:Spectral"; pause = 1)
        spin(w, "filterAvgN", 20)
        sleep(1)
        click(w, "text:Apply filter")
        waitfor(() -> startswith(locate(w, "datasetBox").current, "data2"), "the averaged dataset"; timeout = 180)
        shot(w, "01-data")
    end
    stage(w, "binary model, |V| + T3φ") do
        click(w, "text:Modeling"; pause = 1.5)
        pastepmoired(w, "{'A,ud':0.1, 'B,ud':0.1, 'B,x':3, 'B,y':3, 'B,f':0.6}")
        click(w, "text:Model Optimization"; pause = 1)
        tick(w, "v2", false); tick(w, "t3amp", false); tick(w, "visφ", false)
        tick(w, "visamp", true); tick(w, "t3phi", true)
        choose(w, "param:B,f:mode", "free")
    end
    stage(w, "companion search ±12 mas, 1.5 mas steps") do
        click(w, "companionButton"; pause = 1.5)
        check(w, "default step 1.5 mas", parse(Float64, value(w, "compStep")), 1.5; tol = 0)
        typein(w, "compRmax", "12"; enter = "Tab")
        typein(w, "compFmax", "1"; enter = "Tab")
        click(w, "text:Search"; pause = 2)
        waitfor(() -> value(w, "fitButton") == "Fit", "the search"; timeout = 600)
        sleep(1)
        _, p = fitresult(w, 1)
        near(w, "grid B,f", p["B,f"][1], 0.2879, 0.0047)
        near(w, "grid B,x", p["B,x"][1], -1.035, 0.013)
        near(w, "grid B,y", p["B,y"][1], 3.663, 0.037)
        shot(w, "02-search")
    end
    stage(w, "refit, sizes free") do
        click(w, "text:Adopt"; pause = 1.5)
        for k in ("A,ud", "B,ud", "B,x", "B,y")
            choose(w, "param:$k:mode", "free")
        end
        fitnow(w)
        _, p = fitresult(w, 2)
        near(w, "refit A,ud", p["A,ud"][1], 1.2564, 0.0022)
        near(w, "refit B,f", p["B,f"][1], 0.25364, 0.00043)
        near(w, "refit B,x", p["B,x"][1], -0.9187, 0.0011)
        near(w, "refit B,y", p["B,y"][1], 3.5519, 0.0031)
        shot(w, "03-refit")
    end
    stage(w, "bootstrap, PMOIRED's resampling") do
        click(w, "text:Adopt"; pause = 1.5)
        choose(w, "uncertainty", "bootstrap")
        choose(w, "bootMode", "pmoired")
        spin(w, "nboot", 200)
        fitnow(w; timeout = 900)
        _, p = fitresult(w, 3)
        # PMOIRED, 500 replicates: A,ud 1.2570 +0.0044 −0.0106; B,f 0.25374 +0.00129 −0.00058;
        # B,x −0.9185 ± 0.0024; B,y 3.5505 +0.0045 −0.0095 — σ compared within a factor 2
        for (k, s) in (("A,ud", 0.0075), ("B,f", 0.0009), ("B,x", 0.0024), ("B,y", 0.007))
            check(w, "bootstrap σ $k (×½–×2 of PMOIRED's)", p[k][2] / s, 1.0; tol = 1.0)
        end
        shot(w, "04-bootstrap")
    end
catch e
    println("walk stopped: ", sprint(showerror, e))
finally
    exit(finish(w) ? 0 : 1)
end
