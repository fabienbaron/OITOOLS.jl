# EX2 in the window — FS CMa from PMOIRED's 26 PIONIER files: the notebook's error-bar caps and
# floors, a star + Gaussian built by hand (step 1), then the notebook's param dicts pasted into
# Import PMOIRED (`ex2/ex2_step*.py`, steps 2–6), each with the notebook's free set. Numbers read
# from the window against the script's (GUI_TODO, "EX2 — script" and "EX2 in the window").
#
#     julia test/gui/walks/ex2.jl

include(joinpath(@__DIR__, "walk.jl"))
using .Walk

const DATA = joinpath(homedir(), "SOFTWARE", "PMOIRED_examples", "DATA", "FS_CMa")
isdir(DATA) || error("needs the PMOIRED example data in $DATA")
stepdict(k) = read(joinpath(@__DIR__, "ex2", "ex2_step$k.py"), String)

function freeall(w, comp, keys)
    for k in keys
        choose(w, "param:$comp,$k:mode", "free")
    end
end

# value and σ of `name` in fit `n`, checked against the script within `nσ` of the window's σ
function checkfit(w, n, chi2r_want, want; nσ = 0.05)
    chi2r, p = fitresult(w, n)
    check(w, "fit $n χ²r", chi2r, chi2r_want; tol = 1e-4)
    for (name, v) in want
        check(w, "fit $n $name", p[name][1], v; tol = nσ * p[name][2] + 1e-4)
    end
end

w = start("ex2n"; data_dir = DATA)
try
    stage(w, "open the 26 files as one dataset") do
        click(w, "text:Open OIFITS…")
        click(w, "pickerFiles/text:PIONI.2017-10-21T08_37_13.075_oidataCalibrated.fits")
        shiftclick(w, "pickerFiles/text:PIONIER.2013-02-22T04_41_51_oidataCalibrated.fits")
        click(w, "pickerOk"; pause = 1)
        waitfor(() -> startswith(locate(w, "datasetBox").current, "data1"), "the dataset"; timeout = 180)
        shot(w, "01-open")
    end
    stage(w, "caps and floors") do
        click(w, "text:Error bars")
        typein(w, "err:v2_rel", "0.05"); typein(w, "err:v2_caprel", "0.17")
        typein(w, "err:t3phi_add", "1"); typein(w, "err:t3phi_cap", "5")
        click(w, "text:Apply filter")
        waitfor(() -> startswith(locate(w, "datasetBox").current, "data2"), "the filtered dataset")
        shot(w, "02-filter")
    end
    stage(w, "1. star + Gaussian, by hand") do
        click(w, "text:Modeling"; pause = 1.5)
        for (name, kind) in (("star", "disc"), ("disk", "Gaussian"))
            click(w, "text:+ component"; pause = 1)
            typein(w, "newCompName", name; pause = 0.2, enter = "Tab")
            choose(w, "newCompKind", kind)
            click(w, "text:OK"; pause = 1.5)
        end
        typein(w, "param:star,ud:value", "0.1"); choose(w, "param:star,ud:mode", "fixed")
        typein(w, "param:star,f:value", "1")
        typein(w, "param:disk,fwhm:value", "18"); typein(w, "param:disk,f:value", "1")
        choose(w, "param:disk,f:mode", "free")
        # a flux is a relative weight: free, it has no upper bound (blank)
        check(w, "a free flux unbounded above (1 = yes)", isempty(value(w, "param:disk,f:ub")), 1; tol = 0)
        click(w, "text:Model Optimization"; pause = 1)
        tick(w, "t3amp", false); tick(w, "t3phi", false)
        choose(w, "smearing", "auto")
        fitnow(w)
        check(w, "V² kept", parse(Int, split(value(w, "fit:fit 1:points"))[1]), 612; tol = 0)
        checkfit(w, 1, 35.8255, Dict("disk,fwhm" => 17.841, "disk,f" => 1.5512))
        shot(w, "03-step1")
    end
    stage(w, "2. spectral index") do
        pastepmoired(w, stepdict(2))
        freeall(w, "disk", ("fwhm", "A", "B"))
        fitnow(w)
        checkfit(w, 2, 7.8616, Dict("disk,fwhm" => 17.978, "disk,A" => 1.3981, "disk,B" => 2.3468))
        shot(w, "04-step2")
    end
    stage(w, "3. inclination and PA") do
        pastepmoired(w, stepdict(3))
        freeall(w, "disk", ("fwhm", "A", "B", "incl", "projang"))
        fitnow(w)
        checkfit(w, 3, 6.9306, Dict("disk,fwhm" => 20.785, "disk,incl" => 56.50, "disk,projang" => 66.26))
        shot(w, "05-step3")
    end
    stage(w, "4. doughnut rim + az mode") do
        pastepmoired(w, stepdict(4))
        freeall(w, "rim", ("diamin", "diamout", "az amp1", "az projang1", "incl", "projang", "A", "B"))
        tick(w, "t3phi", true)
        fitnow(w)
        check(w, "V² + T3φ kept", parse(Int, split(value(w, "fit:fit 4:points"))[1]), 975; tol = 0)
        checkfit(w, 4, 8.7949, Dict("rim,diamout" => 32.324, "rim,incl" => 55.449, "rim,projang" => 64.981,
                                    "rim,az amp1" => 0.68668))
        shot(w, "06-step4")
    end
    stage(w, "5. + 1/R² outer disc") do
        pastepmoired(w, stepdict(5))
        freeall(w, "rim", ("A", "B", "az amp1", "az projang1", "diamout", "incl", "projang"))
        freeall(w, "disk", ("A", "B", "az amp1", "az projang1"))
        fitnow(w)
        checkfit(w, 5, 5.6826, Dict("rim,diamout" => 28.877, "rim,incl" => 52.957, "rim,projang" => 68.125))
        shot(w, "07-step5")
    end
    stage(w, "6. off-centre rim and disc") do
        pastepmoired(w, stepdict(6))
        freeall(w, "rim", ("A", "B", "az projang1", "diamout", "incl", "projang", "x", "y"))
        freeall(w, "disk", ("A", "B", "az projang1"))
        fitnow(w)
        checkfit(w, 6, 5.2745, Dict("rim,x" => 1.54, "rim,y" => 1.15, "rim,diamout" => 28.62,
                                    "rim,incl" => 58.31, "rim,projang" => 78.11); nσ = 0.1)
        shot(w, "08-step6")
    end
catch e
    println("walk stopped: ", sprint(showerror, e))
finally
    exit(finish(w) ? 0 : 1)
end
