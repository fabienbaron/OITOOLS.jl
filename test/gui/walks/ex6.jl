# EX6 in the window — HD 210763, one orbit through all of PMOIRED's 11 GRAVITY files: instrument
# GRAVITY_SC averaged ×100 (PMOIRED's averaging), floors 1° on T3φ and 1 % on |V|, the notebook's
# dict pasted into Import PMOIRED, |V| + T3φ; the orbit fitted, then with Gallenne et al.'s radial
# velocities; the Orbit view. Numbers against PMOIRED's printed results.
#
#     julia test/gui/walks/ex6.jl

include(joinpath(@__DIR__, "walk.jl"))
using .Walk

const DATA = joinpath(homedir(), "SOFTWARE", "PMOIRED_examples", "DATA", "HD210763")
const RV = joinpath(Walk.REPO, "demos", "data", "HD210763", "rv_gallenne2023.csv")
isdir(DATA) || error("needs the PMOIRED example data in $DATA")
files = sort(filter(endswith(".fits"), readdir(DATA)))

near(w, label, got, want, σ; nσ = 1.0) = check(w, label, got, want; tol = nσ * σ)

# the notebook's dict (Gallenne et al. 2023 table 6, ω and Ω − 180°), with the semi-major axis
# and the velocity elements of its grid-search cell
const DICT = """{'2,ud': 0.1, '1,ud': 0.2, '2,f': 0.371, '2,x': 'orbit', '2,y': 'orbit',
 '2,orb P': '\$P-42.38 + 42.38', 'P-42.38': 0.0011,
 '2,orb MJD0': '\$MJD0-59800 + 59800', 'MJD0-59800': 54275.93+42.3811*130-59800,
 '2,orb e': 0.6228, '2,orb omega': 293.965-180, '2,orb OMEGA': 257.6-180, '2,orb incl': 71.0,
 '2,orb gamma': 14.94, '2,orb Ka': 50.28, '2,orb Kb': 58.76, '2,orb a': 3.76}"""

# the notebook's exploration (EX6, "grid search orbit"), its globals as PMOIRED names them
const EXPL = "{'rand':{'2,orb incl':(30, 150), '2,orb OMEGA':(0, 180), '2,orb omega':(-180, 180), " *
             "'MJD0-59800': (-20, 20), 'P-42.38': (-1, 1), '2,orb Ka':(45, 60), '2,orb Kb':(45, 60)}}"

w = start("ex6n"; data_dir = DATA)
try
    stage(w, "open 11 files, GRAVITY_SC, ×100 (PMOIRED's averaging), floors") do
        click(w, "text:Open OIFITS…")
        click(w, "pickerFiles/text:" * files[1])
        shiftclick(w, "pickerFiles/text:" * files[end])
        click(w, "pickerOk"; pause = 1)
        waitfor(() -> startswith(locate(w, "datasetBox").current, "data1"), "the dataset"; timeout = 300)
        click(w, "text:Tables"; pause = 1)
        choose(w, "filterInstrument", "GRAVITY_SC")
        click(w, "text:Spectral"; pause = 1)
        spin(w, "filterAvgN", 100)
        sleep(1)
        choose(w, "filterAvgMode", "pmoired")
        click(w, "text:Error bars"; pause = 1)
        typein(w, "err:t3phi_add", "1")
        typein(w, "err:visamp_rel", "0.01")
        click(w, "text:Apply filter")
        waitfor(() -> startswith(locate(w, "datasetBox").current, "data2"), "the filtered dataset"; timeout = 300)
        shot(w, "01-data")
    end
    stage(w, "the notebook's dict, |V| + T3φ") do
        click(w, "text:Modeling"; pause = 1.5)
        pastepmoired(w, DICT)
        click(w, "text:Model Optimization"; pause = 1)
        tick(w, "v2", false); tick(w, "t3amp", false); tick(w, "visφ", false)
        tick(w, "visamp", true); tick(w, "t3phi", true)
        choose(w, "smearing", "auto")
        for k in ("1,ud", "2,f", "2,orb OMEGA", "2,orb a", "2,orb e", "2,orb incl", "2,orb omega",
                  "MJD0_59800", "P_42_38")
            choose(w, "param:$(occursin(',', k) ? k : "__global__," * k):mode", "free")
        end
    end
    stage(w, "the orbit, interferometry only") do
        fitnow(w)
        c2, p = fitresult(w, 1)
        # PMOIRED: χ²r 2.167 (dof 1660)
        check(w, "χ²r (PMOIRED 2.167)", c2, 2.167; tol = 0.05)
        for (k, v, s) in (("2,f", 0.37299, 0.0007), ("2,orb OMEGA", 77.344, 0.066), ("2,orb a", 3.7432, 0.0035),
                          ("2,orb e", 0.62428, 0.00091), ("2,orb incl", 71.027, 0.033),
                          ("2,orb omega", 114.083, 0.060), ("MJD0_59800", -14.549, 0.012))
            near(w, "orbit $k", p[k][1], v, s)
        end
        shot(w, "02-orbit")
    end
    stage(w, "with the radial velocities") do
        click(w, "text:Adopt"; pause = 1.5)
        click(w, "rvOpen"; pause = 1)
        pickfile(w, RV)
        check(w, "RVs ticked once the table is open", locate(w, "fitRV").checked, true; tol = 0)
        for k in ("2,orb Ka", "2,orb Kb", "2,orb gamma")
            choose(w, "param:$k:mode", "free")
        end
        fitnow(w)
        _, p = fitresult(w, 2)
        # ours (script, the same data): Ka 50.325, Kb 58.782, γ 14.931, e 0.62320; PMOIRED (RVs
        # weighted per data block): 50.309, 58.76, 14.929, 0.62295
        for (k, v, s) in (("2,orb Ka", 50.309, 0.054), ("2,orb Kb", 58.76, 0.12), ("2,orb gamma", 14.929, 0.028),
                          ("2,orb e", 0.62295, 0.00052), ("2,orb omega", 114.147, 0.042))
            near(w, "with RVs $k (within 2σ of PMOIRED)", p[k][1], v, s; nσ = 2.0)
        end
        shot(w, "03-rv")
    end
    stage(w, "multi-start: the notebook's 100 random starts") do
        _, p2 = fitresult(w, 2)
        c2 = parse(Float64, value(w, "fit:fit 2:chi2r"))
        click(w, "multiButton"; pause = 1.5)
        click(w, "multiPaste"; pause = 0.2)
        Walk._x(w, "type", "--window", w.wid, "--delay", "2", "--", EXPL)
        click(w, "multiPasteApply"; pause = 1)
        check(w, "expl applied: incl uniform", locate(w, "multi:2,orb incl:kind").current == "uniform", true; tol = 0)
        check(w, "expl applied: MJD0_59800 renamed", locate(w, "multi:MJD0_59800:kind").current == "uniform", true; tol = 0)
        click(w, "text:Search"; pause = 1)
        waitfor(() -> value(w, "fitButton") == "Fit", "the multi-start search"; timeout = 600)
        sleep(1)
        c3, p3 = fitresult(w, 3)
        # the global minimum of 100 random starts is the fit from the published solution
        check(w, "global minimum χ²r = the fit from the published orbit", c3, c2; tol = 1e-3)
        for k in ("2,orb e", "2,orb Ka", "2,orb OMEGA")
            near(w, "global minimum $k", p3[k][1], p2[k][1], p2[k][2])
        end
        shot(w, "04-minima")
        click(w, "multiMinimum"; pause = 0.5)
        press(w, "Home"); press(w, "Down"); press(w, "Return"; pause = 1)
        println("minimum combo shows: ", repr(locate(w, "multiMinimum").current))
        check(w, "second minimum chosen", startswith(locate(w, "multiMinimum").current, "2:"), true; tol = 0)
        c3b = parse(Float64, value(w, "fit:fit 3:chi2r"))
        check(w, "second minimum selected: higher χ²r", c3b > c3, true; tol = 0)
    end
    stage(w, "Orbit view") do
        click(w, "orbitButton"; pause = 2)
        cap = value(w, "diagCaption")
        println("caption: ", cap)
        check(w, "the view's caption names 34 RVs", occursin("34 RVs", cap), true; tol = 0)
        shot(w, "05-orbit-view")
    end
catch e
    println("walk stopped: ", sprint(showerror, e))
finally
    exit(finish(w) ? 0 : 1)
end
