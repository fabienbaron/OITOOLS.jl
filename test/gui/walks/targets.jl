# Observing › Targets in the window: a list of names loaded and resolved, a list loaded over
# another, a list read by column name, a calibrator ticked, saved and read back, both file
# dialogs cancelled, and a save over an existing file declined and accepted.
# Needs the network for SIMBAD.
#
#     julia test/gui/walks/targets.jl
#
# The last stages save the shipped list's contents (twelve targets) over a scratch file.

include(joinpath(@__DIR__, "walk.jl"))
using .Walk

const TMP = mktempdir()
names_file = joinpath(TMP, "names.txt")
write(names_file, "Altair\n# a comment\nDeneb\n")
one_file = joinpath(TMP, "polaris.csv")
write(one_file, "name,ra,dec,V,J,H,K,simbad\nPolaris,37.95456067,89.26410897,1.98,,,1.39,1\n")
# columns by name, in another order; K = 0 exactly, H unknown, a column of the user's own
byname_file = joinpath(TMP, "byname.csv")
write(byname_file, "Dec,Name,RA,K,H,notes,simbad\n89.26410897,Polaris,37.95456067,0.0,,\"pole star, F7\",1\n")
existing = joinpath(TMP, "existing.csv")
write(existing, "old\n")

cached(w) = value(w, "obsCacheCount")
name(w, i) = try value(w, "obsTargetName:$i") catch; nothing end

w = start("targets"; data_dir = TMP)
try
    stage(w, "a file of names loads as unresolved targets") do
        click(w, "text:Observing"; pause = 2)
        click(w, "targetsLoad"; pause = 1)
        pickfile(w, names_file)
        check(w, "two targets, none resolved (1=yes)", cached(w) == "0/2 cached", 1; tol = 0)
        check(w, "in the file's order (1=yes)", name(w, 0) == "Altair" && name(w, 1) == "Deneb", 1; tol = 0)
        click(w, "text:Resolve all"; pause = 1)
        waitfor(() -> cached(w) == "2/2 cached", "both resolved"; timeout = 60)
        check(w, "Altair's RA from SIMBAD", parse(Float64, value(w, "obsTargetRA:0")), 297.6958; tol = 0.01)
        shot(w, "01-names")
    end
    stage(w, "loading replaces the list") do
        click(w, "targetsLoad"; pause = 1)
        pickfile(w, one_file)
        check(w, "one target, resolved (1=yes)", cached(w) == "1/1 cached", 1; tol = 0)
        check(w, "it is Polaris, and the second row is gone (1=yes)",
              name(w, 0) == "Polaris" && name(w, 1) === nothing, 1; tol = 0)
    end
    stage(w, "a file read by column name: unknown magnitudes empty, other columns kept") do
        click(w, "targetsLoad"; pause = 1)
        pickfile(w, byname_file)
        check(w, "Polaris, resolved (1=yes)", name(w, 0) == "Polaris" && cached(w) == "1/1 cached", 1; tol = 0)
        check(w, "RA from its column", parse(Float64, value(w, "obsTargetRA:0")), 37.9546; tol = 1e-3)
        check(w, "K = 0 shown (1=yes)", value(w, "obsTargetMagK:0") == "0.00", 1; tol = 0)
        check(w, "H and V empty (1=yes)",
              value(w, "obsTargetMagH:0") == "" && value(w, "obsTargetMagV:0") == "", 1; tol = 0)
        choose(w, "obsCombiner", "MYSTIC"; pause = 2)
        check(w, "a K setup takes K = 0 (1=yes)", value(w, "simMag") == "0.00", 1; tol = 0)
        click(w, "targetsSave"; pause = 1)
        pickfile(w, joinpath(TMP, "back.csv"))
        txt = read(joinpath(TMP, "back.csv"), String)
        check(w, "saved with its own column, H empty (1=yes)",
              startswith(txt, "name,ra,dec,V,J,H,K,simbad,cal,for,notes\nPolaris,37.95456067,89.26410897,,,,0,1,0,,\"pole star, F7\""), 1; tol = 0)
        choose(w, "obsCombiner", "MIRCX"; pause = 2)
    end
    stage(w, "Cal: the bar's colour, saved and read back") do
        check(w, "not a calibrator on load (1=yes)", !locate(w, "obsTargetCal:0").checked, 1; tol = 0)
        tick(w, "obsTargetCal:0", true)
        click(w, "obsCompute"; pause = 4)
        shot(w, "02-cal")
        click(w, "targetsSave"; pause = 1)
        pickfile(w, joinpath(TMP, "cal.csv"))
        check(w, "saved with cal = 1 (1=yes)",
              occursin("\nPolaris,37.95456067,89.26410897,,,,0,1,1,", read(joinpath(TMP, "cal.csv"), String)), 1; tol = 0)
        tick(w, "obsTargetCal:0", false)
        click(w, "targetsLoad"; pause = 1)
        pickfile(w, joinpath(TMP, "cal.csv"))
        check(w, "read back as a calibrator (1=yes)", locate(w, "obsTargetCal:0").checked, 1; tol = 0)
    end
    stage(w, "cancelling either dialog changes nothing") do
        click(w, "targetsLoad"; pause = 1)
        click(w, "text:Cancel"; pause = 1)
        check(w, "the list is as it was (1=yes)", cached(w) == "1/1 cached" && name(w, 0) == "Polaris", 1; tol = 0)
        nfiles = length(readdir(TMP))
        click(w, "targetsSave"; pause = 1)
        click(w, "text:Cancel"; pause = 1)
        check(w, "no file written (1=yes)", length(readdir(TMP)) == nfiles, 1; tol = 0)
    end
    stage(w, "Load… opens on the shipped lists") do
        click(w, "targetsLoad"; pause = 1)
        click(w, "pickerFiles/text:bright_stars.csv")
        click(w, "pickerOk"; pause = 1)
        check(w, "twelve targets, all resolved (1=yes)", cached(w) == "12/12 cached", 1; tol = 0)
        check(w, "Polaris first (1=yes)", name(w, 0) == "Polaris", 1; tol = 0)
    end
    stage(w, "saving over a file asks first") do
        click(w, "targetsSave"; pause = 1)
        pickfile(w, existing)
        click(w, "text:No"; pause = 1)
        check(w, "declined: the file is untouched (1=yes)", read(existing, String) == "old\n", 1; tol = 0)
        click(w, "text:Cancel"; pause = 1)
        click(w, "targetsSave"; pause = 1)
        pickfile(w, existing)
        click(w, "text:Yes"; pause = 1)
        txt = read(existing, String)
        check(w, "accepted: the list replaces it (1=yes)",
              startswith(txt, "name,ra,dec,V,J,H,K,simbad,cal,for\nPolaris,") && count("\n", txt) == 13, 1; tol = 0)
        shot(w, "03-saved")
    end
catch e
    println(sprint(showerror, e))
finally
    finish(w)
end
