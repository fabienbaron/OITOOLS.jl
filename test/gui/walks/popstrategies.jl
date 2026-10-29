# Observing › POP strategies in the window: the reference telescope's default, the timeline of the
# chart, of a multi-POP plan and of a single setting, on two clocks, a plan adopted with Use; the
# strategy and target menus, and the candidates chart for Capella on a night through LST 0 h.
#
#     julia test/gui/walks/popstrategies.jl

include(joinpath(@__DIR__, "walk.jl"))
using .Walk

# The timeline's "what" column, row by row, until the first row that is not there.
function timeline(w)
    out = String[]
    for k in 0:40
        t = try value(w, "tlWhat:$k") catch; break end
        push!(out, t)
    end
    return out
end
clockrow(w, k) = value(w, "tlFrom:$k")
# the rows of the plans table
count_plans(w) = count(k -> (try locate(w, "planRank:Vega:$k"; timeout = 1).w > 0 catch; false end), 1:20)

w = start("popstrategies")
try
    stage(w, "Observing opens with W2 as the reference") do
        click(w, "text:Observing"; pause = 2)
        check(w, "reference is W2 (1=yes)", locate(w, "obsReference").current == "W2", 1; tol = 0)
    end
    stage(w, "the chart's timeline: one setting") do
        typein(w, "obsDate", "2026-06-21")
        click(w, "obsCompute"; pause = 3)
        click(w, "text:POP strategies"; pause = 2)
        check(w, "title names the chart (1=yes)", value(w, "popTimelineTitle") == "Vega · chart", 1; tol = 0)
        rows = timeline(w)
        check(w, "rows shown (≥ 2)", length(rows) >= 2, 1; tol = 0)
        check(w, "no POP change with one setting (1=yes)", !any(startswith("POP change"), rows), 1; tol = 0)
        check(w, "times as H:MM (1=yes)", occursin(r"^\d+:\d\d$", clockrow(w, 0)), 1; tol = 0)
        shot(w, "01-single")
    end
    stage(w, "POP changes on a single line: a change and a drop in the timeline") do
        choose(w, "popStrategy", "POP changes on a single line")
        click(w, "obsCompute"; pause = 4)
        rows = timeline(w)
        check(w, "a POP change row (1=yes)", any(startswith("POP change"), rows), 1; tol = 0)
        shot(w, "02-multi-lst")
    end
    stage(w, "the clock: Local moves the times") do
        lst0 = clockrow(w, 0)
        choose(w, "obsTimeSystem", "Local"; pause = 3)
        loc0 = clockrow(w, 0)
        check(w, "first time changed with the clock (1=yes)", loc0 != lst0, 1; tol = 0)
        check(w, "same rows (1=yes)", length(timeline(w)) >= 3, 1; tol = 0)
        shot(w, "03-multi-local")
    end
    stage(w, "Search POPs: a plan's timeline, then Use") do
        click(w, "text:Search POPs"; pause = 6)
        click(w, "planRank:Vega:2"; pause = 2)
        check(w, "title names plan 2 (1=yes)", value(w, "popTimelineTitle") == "Vega · plan 2", 1; tol = 0)
        picked = timeline(w)
        shot(w, "04-plan2")
        click(w, "planUse:Vega:2"; pause = 4)
        check(w, "after Use, the chart's (1=yes)", value(w, "popTimelineTitle") == "Vega · chart", 1; tol = 0)
        check(w, "the chart's timeline is plan 2's (1=yes)", timeline(w) == picked, 1; tol = 0)
        click(w, "popTimelineCopy"; pause = 1)
        check(w, "copied (1=yes)", value(w, "obsStatus") == "timeline copied", 1; tol = 0)
        shot(w, "05-used")
    end
    stage(w, "Single POPs: a setting's timeline, the table on the same clock") do
        click(w, "planUse:Vega:2"; pause = 3)           # back to AutoPOPs
        choose(w, "popStrategy", "Single POPs")
        click(w, "text:Search POPs"; pause = 4)
        click(w, "popRank:1"; pause = 2)
        check(w, "title names the setting (1=yes)",
              startswith(value(w, "popTimelineTitle"), "Vega · ") &&
              value(w, "popTimelineTitle") != "Vega · chart", 1; tol = 0)
        check(w, "no POP change (1=yes)", !any(startswith("POP change"), timeline(w)), 1; tol = 0)
        shot(w, "06-setting")
    end
    stage(w, "the strategy menu: section titles cannot be chosen") do
        before = locate(w, "popStrategyTab").current
        click(w, "popStrategyTab"; pause = 0.5)
        for title in ("— single target", "— multiple targets")
            l = locate(w, "popup:" * title)
            check(w, "\"$title\" listed, not choosable (1=yes)", l.w > 0 && !l.enabled, 1; tol = 0)
        end
        l = locate(w, "popup:POP changes on a single line, all targets")
        check(w, "the list's strategy choosable (1=yes)", l.w > 0 && l.enabled, 1; tol = 0)
        press(w, "Escape")
        check(w, "unchanged (1=yes)", locate(w, "popStrategyTab").current == before, 1; tol = 0)
        choose(w, "popStrategyTab", "Three-row POP (time coverage)")
        check(w, "the left menu follows (1=yes)",
              locate(w, "popStrategy").current == "Three-row POP (time coverage)", 1; tol = 0)
    end
    stage(w, "Capella through LST 0 h: target menu, candidates, Use") do
        click(w, "text:Add target"; pause = 1)
        typein(w, "obsTargetName:1", "Capella")
        typein(w, "obsTargetRA:1", "79.172")
        typein(w, "obsTargetDec:1", "45.998")
        typein(w, "obsDate", "2026-09-15")
        choose(w, "popStrategyTab", "POP changes && Telescope Drops on a single line")
        choose(w, "popTarget", "Capella"; pause = 3)
        check(w, "the list follows the menu (1=yes)",
              startswith(value(w, "popTimelineTitle"), "Capella"), 1; tol = 0)
        click(w, "popSearch"; pause = 8)
        check(w, "candidates drawn (1=yes)", locate(w, "plansGanttArea").w > 0, 1; tol = 0)
        click(w, "planRank:Capella:1"; pause = 2)
        rows = timeline(w)
        check(w, "first gap is elevation-restricted (1=yes)", first(rows) == "elevation-restricted", 1; tol = 0)
        shot(w, "07-capella-candidates")
        click(w, "popTimelineUse"; pause = 4)
        check(w, "after Use, the chart's (1=yes)", value(w, "popTimelineTitle") == "Capella · chart", 1; tol = 0)
        check(w, "the chart's timeline is the candidate's (1=yes)", timeline(w) == rows, 1; tol = 0)
        shot(w, "08-capella-used")
    end
    stage(w, "the POP change time, in the tab and in the POPs panel") do
        choose(w, "popStrategyTab", "POP changes on a single line"; pause = 3)
        spin(w, "popChangeMinutesTab", 20)
        check(w, "the panel's box follows (1=yes)", value(w, "popChangeMinutes") == "20", 1; tol = 0)
        spin(w, "popChangeMinutes", 10)
        check(w, "the tab's box follows (1=yes)", value(w, "popChangeMinutesTab") == "10", 1; tol = 0)
        choose(w, "popStrategyTab", "POP changes && Telescope Drops on a single line"; pause = 3)
    end
    stage(w, "Vega on 2026-09-15: W2 dropped, the reference moves to S2 then E2") do
        choose(w, "popTarget", "Vega"; pause = 3)
        click(w, "obsCompute"; pause = 4)
        rows = timeline(w)
        check(w, "Ref → S2 then Ref → E2 in the timeline (1=yes)",
              (a = findfirst(==("Ref → S2"), rows); b = findfirst(==("Ref → E2"), rows);
               a !== nothing && b !== nothing && a < b), 1; tol = 0)
        shot(w, "09-vega-reference-timeline")
        click(w, "text:Gantt"; pause = 3)
        shot(w, "10-vega-reference-gantt")
    end
    stage(w, "Precise: change times per POP and of the carts") do
        click(w, "text:POP strategies"; pause = 2)
        choose(w, "popStrategyTab", "POP changes on a single line"; pause = 3)
        click(w, "popSearch"; pause = 6)
        spin(w, "popPerPop", 5)
        typein(w, "popCartSpeed", "0.05"; pause = 1)
        click(w, "popPrecise"; pause = 8)
        check(w, "the tables say precise (1=yes)", endswith(value(w, "obsStatus"), "plans, precise"), 1; tol = 0)
        check(w, "a first plan (1=yes)", locate(w, "planRank:Vega:1").w > 0, 1; tol = 0)
        shot(w, "11-precise")
        click(w, "planUse:Vega:1"; pause = 4)
        check(w, "used (1=yes)", value(w, "planUse:Vega:1") == "in use", 1; tol = 0)
        click(w, "planUse:Vega:1"; pause = 3)
        click(w, "popSearch"; pause = 6)
        check(w, "Search POPs is the fast search again (1=yes)",
              value(w, "obsStatus") == "$(count_plans(w)) plans", 1; tol = 0)
    end
catch e
    println("walk stopped: ", sprint(showerror, e))
finally
    exit(finish(w) ? 0 : 1)
end
