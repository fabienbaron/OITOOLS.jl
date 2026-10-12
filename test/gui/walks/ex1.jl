# EX1 in the window — α Cen A's diameter from PMOIRED's 8 PIONIER files (V² only): uniform disc,
# Claret-4 with fixed coefficients, power law, bootstrap. Each number is read from the window and
# compared with the script's (GUI_TODO, "EX1 — script"), at the printed digits.
#
#     julia test/gui/walks/ex1.jl

include(joinpath(@__DIR__, "walk.jl"))
using .Walk

const DATA = joinpath(homedir(), "SOFTWARE", "PMOIRED_examples", "DATA", "alphaCenA")
const MODELS = joinpath(Walk.REPO, "demos", "models")
isdir(DATA) || error("needs the PMOIRED example data in $DATA")

w = start("ex1n"; data_dir = DATA)
try
    stage(w, "open the 8 files as one dataset") do
        click(w, "text:Open OIFITS…")
        click(w, "pickerFiles/text:20160523_AlphaCenA_1.fits")
        shiftclick(w, "pickerFiles/text:PIONI.2016-05-31T04_39_53.560_oidataCalibrated.fits")
        click(w, "pickerOk"; pause = 1)
        waitfor(() -> startswith(locate(w, "datasetBox").current, "data1"), "the dataset"; timeout = 120)
        shot(w, "01-open")
    end
    stage(w, "V² relative error floor 0.01") do
        click(w, "text:Error bars")
        typein(w, "err:v2_rel", "0.01")
        click(w, "text:Apply filter")
        waitfor(() -> startswith(locate(w, "datasetBox").current, "data2"), "the filtered dataset")
        shot(w, "02-filter")
    end
    stage(w, "uniform disc") do
        click(w, "text:Modeling"; pause = 1.5)
        click(w, "text:Open…"; pause = 1)
        pickfile(w, joinpath(MODELS, "Uniform_Disk.toml"))
        click(w, "text:Model Optimization"; pause = 1)
        tick(w, "t3amp", false); tick(w, "t3phi", false)
        choose(w, "smearing", "auto")
        choose(w, "param:star,ud:mode", "free")
        fitnow(w)
        chi2r, p = fitresult(w, 1)
        check(w, "UD θ", p["star,ud"][1], 8.2989; tol = 1e-4)
        check(w, "UD σθ", p["star,ud"][2], 0.00383; tol = 1e-5)
        check(w, "UD χ²r", chi2r, 18.9908; tol = 1e-4)
        shot(w, "03-ud")
    end
    stage(w, "Claret-4, coefficients fixed") do
        click(w, "text:Open…"; pause = 1)
        pickfile(w, joinpath(MODELS, "Limb_Darkened_Disk_Claret4.toml"))
        check(w, "smearing kept on auto (1 = yes)", locate(w, "smearing").current == "auto", 1; tol = 0)
        for (k, v) in zip(("c1", "c2", "c3", "c4"), ("0.7127", "-0.0452", "-0.2643", "0.1311"))
            typein(w, "param:star,$k:value", v)
            choose(w, "param:star,$k:mode", "fixed")
        end
        fitnow(w)
        chi2r, p = fitresult(w, 2)
        check(w, "Claret-4 θ", p["star,ldclaret4"][1], 8.4945; tol = 1e-4)
        check(w, "Claret-4 σθ", p["star,ldclaret4"][2], 0.00222; tol = 1e-5)
        check(w, "Claret-4 χ²r", chi2r, 4.4487; tol = 1e-4)
        shot(w, "04-claret")
    end
    stage(w, "power law") do
        click(w, "text:Adopt"; pause = 1)
        click(w, "text:Model Creation"; pause = 1)
        choose(w, "ldLaw", "power law")
        choose(w, "param:star,alpha:mode", "free")
        check(w, "α upper bound", parse(Float64, value(w, "param:star,alpha:ub")), 3.0; tol = 0)
        click(w, "text:Model Optimization"; pause = 1)
        fitnow(w)
        chi2r, p = fitresult(w, 3)
        check(w, "power θ", p["star,ldpow"][1], 8.4644; tol = 1e-4)
        check(w, "power σθ", p["star,ldpow"][2], 0.00522; tol = 1e-5)
        check(w, "power α", p["star,alpha"][1], 0.14665; tol = 2e-5)
        check(w, "power σα", p["star,alpha"][2], 0.00414; tol = 1e-5)
        check(w, "power χ²r", chi2r, 4.0533; tol = 1e-4)
        shot(w, "05-power")
    end
    stage(w, "bootstrap 100, seed 1") do
        click(w, "text:Adopt"; pause = 1)
        choose(w, "uncertainty", "bootstrap")
        spin(w, "nboot", 100)
        fitnow(w)
        _, p = fitresult(w, 4)
        check(w, "bootstrap σθ", p["star,ldpow"][2], 0.0105; tol = 1e-4)
        check(w, "bootstrap σα", p["star,alpha"][2], 0.00924; tol = 1e-5)
        click(w, "text:Correlations"; pause = 2)
        check(w, "bootstrap corr(θ, α)", parse(Float64, value(w, "corr:1:0")), 0.84; tol = 0.005)
        shot(w, "06-bootstrap")
    end
catch e
    println("walk stopped: ", sprint(showerror, e))
finally
    exit(finish(w) ? 0 : 1)
end
