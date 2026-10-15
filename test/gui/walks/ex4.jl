# EX4 part 1 in the window — AMHRA's Be star + Keplerian disc cube (PMOIRED's example 4) observed
# with CHARA: Tools (Set keywords, its reference pixel kept; Save FITS; four channels side by
# side, stretched as PMOIRED shows them), Observing (MYSTIC on the cube's own channels, (0h, −24°)
# at transit, PMOIRED's fractional noise), Exploring. Checks: the saved cube's keywords, the
# simulated file's error bars, and its χ² against the noiseless truth.
#
#     julia --project=. test/gui/walks/ex4.jl      (the checks read the files with OITOOLS)

include(joinpath(@__DIR__, "walk.jl"))
using .Walk
using OITOOLS

mean(x) = sum(x) / length(x)

const DATA = joinpath(homedir(), "SOFTWARE", "PMOIRED_examples", "DATA", "BeModelsAMHRA")
const CUBE = joinpath(DATA, "output_Be_1675260915739.fits")
isfile(CUBE) || error("needs the PMOIRED example data in $DATA")

w = start("ex4n"; data_dir = DATA)
const KEYED = joinpath(w.out, "amhra_keys.fits")
const OUT = joinpath(w.out, "amhra_chara.oifits")
try
    stage(w, "Tools: Set keywords, Save FITS") do
        click(w, "text:Tools"; pause = 1.5)
        typein(w, "toolsFitsPath", CUBE; pause = 3)
        click(w, "toolsSetKeywords"; pause = 1)
        # prefilled from the header: rad pixels, CUNIT3 m, CRPIX 128
        check(w, "pixel prefilled (mas)", parse(Float64, value(w, "kwPix")), 0.0548345; tol = 1e-6)
        check(w, "λ start prefilled (µm)", parse(Float64, value(w, "kwWStart")), 2.15962; tol = 1e-6)
        check(w, "reference pixel x prefilled", parse(Float64, value(w, "kwRefX")), 128; tol = 0)
        shot(w, "01-keywords")
        click(w, "text:Apply"; pause = 2)
        click(w, "toolsSaveFits"; pause = 1)
        typein(w, "pickerName", KEYED; pause = 3)
        m = read_model_image(KEYED)
        check(w, "saved: 51 channels", length(m.wavelengths), 51; tol = 0)
        check(w, "saved: λ 25 (µm)", m.wavelengths[25] * 1e6, 2.1661; tol = 1e-9)
        check(w, "saved: offset (mas), CRPIX 128 kept", m.offset[1], 0.0548345; tol = 1e-6)
    end
    stage(w, "Tools: four channels side by side, power 0.2") do
        typein(w, "toolsFitsPath", KEYED; pause = 3)
        tick(w, "toolsTiles", true)
        typein(w, "toolsTileList", "2.1650 2.1658 2.1664 2.1701"; pause = 2)
        typein(w, "toolsPower", "0.2"; pause = 2)
        check(w, "side by side shown (1 = yes)", locate(w, "toolsTilesArea").w > 0, 1; tol = 0)
        shot(w, "02-channels")
        tick(w, "toolsTiles", false)
    end
    stage(w, "Observing: CHARA + MYSTIC, the cube's channels, fractional noise") do
        click(w, "text:Observing"; pause = 2)
        typein(w, "obsTargetName:0", "AMHRA_Be")
        typein(w, "obsTargetRA:0", "0")
        typein(w, "obsTargetDec:0", "-24")
        choose(w, "obsCombiner", "MYSTIC")
        typein(w, "obsDate", "2024-09-22")
        typein(w, "obsHaMin", "0")
        typein(w, "obsHaMax", "0")
        choose(w, "simSourceBox", "Image")
        click(w, "text:FITS…"; pause = 1)
        typein(w, "pickerName", KEYED; pause = 3)
        tick(w, "simImageChannels", true)
        choose(w, "simNoiseMode", "fractional")
        for (k, v) in (("v2", "0.01"), ("visamp", "0.01"), ("t3amp", "0.01"), ("visphi", "0.01"),
                       ("t3phi", "0.01"), ("flux", "0.01"))
            typein(w, "simFrac:$k", v; enter = "Tab")
        end
        click(w, "text:Output file…"; pause = 1)
        typein(w, "pickerName", OUT; pause = 2)
        shot(w, "03-observe")
        click(w, "text:Simulate"; pause = 2)
        waitfor(() -> isfile(OUT), "the simulated file"; timeout = 300)
        sleep(3)
        d = readoifits(OUT; T = Float64, filter_bad_data = false, warn = false, verbose = false)[1, 1]
        check(w, "V² points (15 baselines × 51)", d.nv2, 765; tol = 0)
        check(w, "T3φ error (°)", maximum(d.t3phi_err), 0.01; tol = 1e-9)
        check(w, "V² error / V² (floored)", maximum(d.v2_err ./ max.(abs.(d.v2), 0.01)), 1e-4; tol = 2e-5)
        # the truth on the same uv, no noise: χ²r of the data against it is 1 within its scatter
        truth = joinpath(w.out, "truth.oifits")
        simulate_from_oifits(OUT, truth; mode = "copy_errors", image = KEYED, noise = false)
        q = readoifits(truth; T = Float64, filter_bad_data = false, warn = false, verbose = false)[1, 1]
        check(w, "χ²r V² against the truth", mean(((d.v2 .- q.v2) ./ d.v2_err) .^ 2), 1.0; tol = 0.15)
        check(w, "χ²r T3φ against the truth", mean(((d.t3phi .- q.t3phi) ./ d.t3phi_err) .^ 2), 1.0; tol = 0.15)
    end
    stage(w, "Exploring: the simulated data") do
        click(w, "text:Open in Explore"; pause = 3)
        waitfor(() -> startswith(locate(w, "datasetBox").current, "data"), "the dataset"; timeout = 120)
        shot(w, "04-explore")
    end
catch e
    println("walk stopped: ", sprint(showerror, e))
finally
    exit(finish(w) ? 0 : 1)
end
