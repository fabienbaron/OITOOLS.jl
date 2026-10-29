# Modeling › several models in the window: a model opened, edited, a second opened beside it,
# a fit on the first, back and forth with the edit and the fit kept, Adopt refused for a fit of
# the other model, New, Duplicate, rename, Remove.
#
#     julia test/gui/walks/models.jl

include(joinpath(@__DIR__, "walk.jl"))
using .Walk

const DATA = joinpath(Walk.REPO, "demos", "data")
const MODELS = joinpath(Walk.REPO, "demos", "models")

# shown: found and of some width (the locator finds no hidden item)
shown(w, q) = try locate(w, q; timeout = 2).w > 0 catch; false end
current(w) = locate(w, "modelSelect").current
# The box is editable: a click on its text edits the name, so open its list from the arrow.
function pick(w, name)
    l = locate(w, "modelSelect")
    Walk._x(w, "mousemove", "--window", w.wid, string(l.x + l.w - 8), string(l.y + l.h ÷ 2), "click", "1")
    sleep(0.5)
    click(w, "popup:" * name; pause = 1)
    current(w) == name || error("modelSelect shows \"$(current(w))\" after picking \"$name\"")
end

w = start("models"; data_dir = DATA)
try
    stage(w, "a dataset and a first model") do
        click(w, "text:Open OIFITS…")
        click(w, "pickerFiles/text:AlphaCenA.oifits")
        click(w, "pickerOk"; pause = 1)
        waitfor(() -> startswith(locate(w, "datasetBox").current, "data1"), "the dataset"; timeout = 120)
        click(w, "text:Modeling"; pause = 1.5)
        check(w, "an empty model 1 (1=yes)", current(w) == "model 1", 1; tol = 0)
        click(w, "text:Open…"; pause = 1)
        pickfile(w, joinpath(MODELS, "Uniform_Disk.toml"))
        check(w, "named after the file, in place of model 1 (1=yes)", current(w) == "Uniform_Disk", 1; tol = 0)
        check(w, "not marked unsaved (1=yes)", !shown(w, "modelDirtyMark"), 1; tol = 0)
        typein(w, "param:star,ud:value", "8")
        check(w, "an edit marks it (1=yes)", shown(w, "modelDirtyMark"), 1; tol = 0)
    end
    stage(w, "a fit, then a second model beside the first") do
        click(w, "text:Model Optimization"; pause = 1)
        choose(w, "param:star,ud:mode", "free")
        fitnow(w)
        check(w, "fit 1 (1=yes)", value(w, "fit:fit 1:label") == "fit 1", 1; tol = 0)
        click(w, "text:Open…"; pause = 1)
        pickfile(w, joinpath(MODELS, "AX_Cir_Binary.toml"))
        check(w, "the binary selected (1=yes)", current(w) == "AX_Cir_Binary", 1; tol = 0)
        check(w, "the fit names its model (1=yes)",
              value(w, "fit:fit 1:label") == "fit 1 · Uniform_Disk", 1; tol = 0)
        tick(w, "Console", true)
        click(w, "text:Adopt"; pause = 1)
        check(w, "Adopt refused for the other model's fit (1=yes)",
              occursin("the latest fit is of Uniform_Disk", value(w, "consoleText")), 1; tol = 0)
        tick(w, "Console", false)
        shot(w, "01-two")
    end
    stage(w, "back to the first: its edit and its mark kept") do
        pick(w, "Uniform_Disk")
        check(w, "the disc's parameter (1=yes)", shown(w, "param:star,ud:value"), 1; tol = 0)
        check(w, "still unsaved (1=yes)", shown(w, "modelDirtyMark"), 1; tol = 0)
        check(w, "the free parameter kept (1=yes)",
              locate(w, "param:star,ud:mode").current == "free", 1; tol = 0)
        pick(w, "AX_Cir_Binary")
        check(w, "the binary's companion (1=yes)", shown(w, "param:companion,x:value"), 1; tol = 0)
        check(w, "the binary not marked (1=yes)", !shown(w, "modelDirtyMark"), 1; tol = 0)
        shot(w, "02-back")
    end
    stage(w, "New, Duplicate, rename, Remove") do
        click(w, "modelNew"; pause = 1)
        check(w, "model 3 (1=yes)", current(w) == "model 3", 1; tol = 0)
        check(w, "empty (1=yes)", !shown(w, "param:star,ud:value"), 1; tol = 0)
        shot(w, "03-new")
        pick(w, "AX_Cir_Binary")
        click(w, "modelDuplicate"; pause = 1)
        check(w, "a copy (1=yes)", current(w) == "AX_Cir_Binary copy", 1; tol = 0)
        check(w, "with the companion (1=yes)", shown(w, "param:companion,x:value"), 1; tol = 0)
        typein(w, "modelSelect", "binary B"; pause = 1)
        check(w, "renamed (1=yes)", current(w) == "binary B", 1; tol = 0)
        typein(w, "modelSelect", "Uniform_Disk"; pause = 1)
        check(w, "a taken name refused, and nothing else selected (1=yes)", current(w) == "binary B", 1; tol = 0)
        sleep(1)
        check(w, "still binary B (1=yes)", current(w) == "binary B", 1; tol = 0)
        click(w, "modelRemove"; pause = 1)
        try
            locate(w, "text:Discard"; timeout = 2)
            click(w, "text:Discard"; pause = 1)
        catch
        end
        check(w, "removed, the one before selected (1=yes)", current(w) == "model 3", 1; tol = 0)
        shot(w, "04-list")
    end
catch e
    println(sprint(showerror, e))
finally
    finish(w)
end
