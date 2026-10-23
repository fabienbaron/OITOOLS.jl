# Observing › Targets in the window: a list of names loaded and resolved, a list loaded over
# another, both file dialogs cancelled, and a save over an existing file declined and accepted.
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
              startswith(txt, "name,ra,dec,V,J,H,K,simbad\nPolaris,") && count("\n", txt) == 13, 1; tol = 0)
        shot(w, "02-saved")
    end
catch e
    println(sprint(showerror, e))
finally
    finish(w)
end
