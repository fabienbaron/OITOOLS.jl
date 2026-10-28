# Observing › a list of targets under one shared POP plan in the window: a list with two
# calibrators loaded, the list's strategy computed, a row picked, the shared plans searched and
# one used, the settings of the POP strategies tab, the POPs set by hand, and back to one target.
#
#     julia test/gui/walks/multitarget.jl

include(joinpath(@__DIR__, "walk.jl"))
using .Walk

const TMP = mktempdir()
list_file = joinpath(TMP, "list.csv")
write(list_file, """
name,ra,dec,V,J,H,K,simbad,cal,for
Alpheratz,2.0965,29.0904,,,,,0,0,
Mirach,17.4330,35.6206,,,,,0,0,
Hamal,31.7934,23.4624,,,,,0,0,
cal1,4.0,28.0,,,,,0,1,Alpheratz
cal2,33.5,22.0,,,,,0,1,Hamal;Mirach
""")

status(w) = value(w, "obsStatus")
# shown: found and of some width (the locator finds no hidden item)
shown(w, q) = try locate(w, q; timeout = 2).w > 0 catch; false end

w = start("multitarget"; data_dir = TMP)
try
    stage(w, "a list with two calibrators") do
        click(w, "text:Observing"; pause = 2)
        click(w, "targetsLoad"; pause = 1)
        pickfile(w, list_file)
        check(w, "five targets (1=yes)", value(w, "obsTargetName:4") == "cal2", 1; tol = 0)
        check(w, "cal2 is for Hamal and Mirach (1=yes)", value(w, "obsTargetFor:4") == "Hamal;Mirach", 1; tol = 0)
        check(w, "a science target has no for field (1=yes)", !shown(w, "obsTargetFor:0"), 1; tol = 0)
    end
    stage(w, "the list under one shared plan") do
        typein(w, "obsDate", "2026-11-09")
        choose(w, "popStrategy", "POP changes on a single line, all targets")
        click(w, "obsCompute"; pause = 5)
        check(w, "five targets charted (1=yes)", startswith(status(w), "5 targets · "), 1; tol = 0)
        check(w, "a score (1=yes)", occursin("score", status(w)), 1; tol = 0)
        shot(w, "01-list")
    end
    stage(w, "a row picked: the chart's timeline is that target's") do
        click(w, "obsTargetDot:2"; pause = 4)
        click(w, "text:POP strategies"; pause = 2)
        check(w, "timeline names Hamal (1=yes)", value(w, "popTimelineTitle") == "Hamal · chart", 1; tol = 0)
        check(w, "no target menu for a list (1=yes)", !shown(w, "popTarget"), 1; tol = 0)
        check(w, "the calibrator settings shown (1=yes)",
              shown(w, "popBracket") && shown(w, "popNeedsCal"), 1; tol = 0)
    end
    stage(w, "Search POPs: the shared plans, then Use") do
        click(w, "popSearch"; pause = 6)
        check(w, "a first plan listed (1=yes)", locate(w, "planRank:all:1").w > 0, 1; tol = 0)
        shot(w, "02-plans")
        click(w, "planUse:all:1"; pause = 5)
        check(w, "in use (1=yes)", value(w, "planUse:all:1") == "in use", 1; tol = 0)
        click(w, "planUse:all:1"; pause = 5)
        check(w, "back to AutoPOPs (1=yes)", value(w, "planUse:all:1") == "Use", 1; tol = 0)
    end
    stage(w, "the tab's settings: no calibrators needed, half a night, a longer stretch") do
        before = status(w)
        tick(w, "needs cal", false)
        sleep(4)
        check(w, "still five targets (1=yes)", startswith(status(w), "5 targets · "), 1; tol = 0)
        typein(w, "popNightFrom", "7"; pause = 4)
        check(w, "half a night computed (1=yes)", startswith(status(w), "5 targets · "), 1; tol = 0)
        shot(w, "03-half-night")
        typein(w, "popNightFrom", ""; enter = "BackSpace")       # emptied: the whole night
        press(w, "Return"; pause = 4)
        spin(w, "popMinStretch", 60)
        sleep(4)
        check(w, "with a 60-minute stretch (1=yes)", startswith(status(w), "5 targets · "), 1; tol = 0)
        tick(w, "needs cal", true)
        sleep(4)
        spin(w, "popMinStretch", 15)
        sleep(4)
        check(w, "back to the first chart (1=yes)", status(w) == before, 1; tol = 0)
    end
    stage(w, "the POPs set by hand: one setting for the list") do
        tick(w, "AutoPOPs", false)
        click(w, "obsCompute"; pause = 5)
        check(w, "no change (1=yes)", startswith(status(w), "5 targets · 0 changes"), 1; tol = 0)
        shot(w, "04-by-hand")
        tick(w, "AutoPOPs", true)
    end
    stage(w, "back to one target") do
        choose(w, "popStrategy", "POP changes on a single line")
        click(w, "obsCompute"; pause = 4)
        check(w, "one target's chart (1=yes)", !startswith(status(w), "5 targets"), 1; tol = 0)
        click(w, "text:POP strategies"; pause = 2)
        check(w, "the target menu is back (1=yes)", shown(w, "popTarget"), 1; tol = 0)
        shot(w, "05-single")
    end
catch e
    println(sprint(showerror, e))
finally
    finish(w)
end
