// OITOOLS — the Model perspective: one model dict, edited three ways and explained once.
//
// The central design problem is §3.2 of the plan, and it is not "how do I show a checkbox".
// In OITOOLS a parameter is in exactly one of three states, and they are NOT independent: a
// name in `list_free_params` must be numeric, because the resolver raises on a string there.
// Free and derived are therefore mutually exclusive, and the table has to make that
// unrepresentable rather than report it after a failed fit. Hence a three-way mode selector
// per row, not a "free" tick — a tick has no way to say "this is an expression".
//
// The second thing this perspective exists for is the inspector (§3.3). An unrecognised
// geometry key does not error: `"ring,gaussian_ring"` alongside a `fwhm` silently yields a
// plain Gaussian and one inert parameter. Two shipped demos did exactly that. The parser can
// already name the keys it ignored, so the inspector shows them instead of leaving the user to
// wonder why a ring fits like a blob.
//
// QML holds no model of its own. Every list below is filled by `refreshModel` from the shell
// and every edit goes back through it, so the table cannot drift from what a fit or an
// optimiser will see. The data layer is `src/gui/model.jl` — `model_rows`, `model_inspection`
// and `free_parameter_vector`, whose fields are exactly the roles of `paramModel` below.

import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml            // supplies `Julia`, through which the model is rendered
import Makie            // MakieArea, for this perspective's own canvas

Item {
    id: root

    // ── scaling, supplied by the window ───────────────────────────────────────
    //
    // Same contract as Main.qml: the window has already worked out the screen's scale, and a
    // tab that re-derived it would drift away from the chrome around it.
    property real uiScale:   1.0
    property real fontScale: 1.0
    property real baseFontPt: 11

    function dp(px)     { return Math.round(px * uiScale) }
    function pt(points) { return points * fontScale }

    // ── palette ──────────────────────────────────────────────────────────────
    // A paragraph of explanation costs vertical space on every visit and is read once. The
    // glyph keeps it one hover away without spending the room, which is what lets the panel
    // show the model rather than a manual about the model.
    // Drawn rather than typed: U+24D8 (CIRCLED LATIN SMALL LETTER I) is missing from the UI
    // font here and renders as a tofu box. A circle with an `i` in it always renders.
    component InfoTip: Rectangle {
        property string tip: ""
        implicitWidth: root.dp(14)
        implicitHeight: root.dp(14)
        radius: width / 2
        color: infoHover.hovered ? "#5b7fb0" : "#a8a8a8"
        Label {
            anchors.centerIn: parent
            text: "i"
            color: "white"
            font.bold: true
            font.pointSize: root.pt(root.baseFontPt - 3)
        }
        HoverHandler { id: infoHover; cursorShape: Qt.WhatsThisCursor }
        ToolTip {
            visible: infoHover.hovered && tip.length > 0
            delay: 150
            // A fixed width, because a ToolTip sized by its text runs the full width of the
            // screen for a paragraph and becomes unreadable.
            contentItem: Text {
                text: tip
                wrapMode: Text.WordWrap
                width: root.dp(360)
                font.pointSize: root.pt(root.baseFontPt - 1)
            }
        }
    }

    readonly property color cFixed:  "#8a8a8a"
    readonly property color cFree:   "#1a7f37"
    readonly property color cExpr:   "#7048b6"
    readonly property color cBad:    "#c02020"
    readonly property color cWarn:   "#b06000"
    readonly property color cPanel:  "#f6f6f6"
    readonly property color cLine:   "#dcdcdc"

    function modeColor(m) {
        return m === "free" ? cFree : m === "expr" ? cExpr : cFixed
    }

    // ── the model, as the parser sees it ─────────────────────────────────────
    //
    // `paramModel` mirrors `ParamRow` one-for-one, so Julia fills it straight from
    // `model_rows(model_dict, list_free_params; lb, ub)` with no reshaping. `fitindex` is the
    // position in `list_free_params`, and therefore the index into every `x` an optimiser
    // touches and every `x_opt` a result reports — showing it is what makes an exported
    // script legible.
    property string modelName: "untitled"
    property bool   modelDirty: false
    property string modelPath: ""

    // Starts empty. A model invented by the interface is a model the user did not write, and
    // there is no way to tell it apart from one they did once it is on screen -- the fitted
    // numbers, the χ² and the exported script would all describe a star nobody asked about.
    // Open a model file (demos/models has a uniform disc and a limb-darkened one) or add
    // components.
    ListModel {
        id: paramModel
        // roles: comp, param, key, mode, value, expr, lb, ub, fitindex, atbound, kind
        // (`comp` is ParamRow.component renamed: `component` is a QML keyword)
        // Roles are defined by the first append, so every append must carry all eleven.
        // Filled from `model_rows(model_dict, list_free_params; lb, ub)`, in its order.
    }

    // Components, as the parser identified them (§3.3 (i)). `geometryKey` is the key that
    // DECIDED `kind`: `_GEOMETRY_KEYS` is searched in order and the first match wins, so
    // naming the deciding key is worth more than naming the kind alone.
    ListModel {
        id: componentModel
        // roles: name, kind, geometryKey, nparams, nunrecognised
        // Filled from `model_inspection(model_dict).components`.
    }

    // Keys the parser ignored. Not an error, and that is the whole problem: they change what
    // the component IS, silently.
    property var unrecognisedKeys: []          // model_inspection(...).unrecognised
    property bool broadcasting: false          // model_inspection(...).broadcasting
    property var globalKeys: []                // model_inspection(...).globals

    // Kinds the "+ component" dialog offers, as {key, label}. Read from Julia so the list is
    // the parser's own rather than a copy that can fall behind it.
    // [{ key, label, subs: [{ key, label, base }] }], in menu order.
    property var componentKinds: []
    property var ldLaws: []
    property var profileTemplates: []

    // Carried from a chosen template until its parameters are created: the seeds that describe
    // the shape, and the grid radius it needs. `_component_seed`'s generic 2.0 describes
    // nothing, and a grid too small truncates the profile without saying so.
    property string profileSeeds: ""
    property real   profileSuggestedRmax: 0

    // I(r_max)/max(I), and only when the profile is still falling there. A shape cut off
    // mid-decline is being modelled as something narrower than what was written.
    property real   profileEdgeFrac: 0

    // Which optional geometry the selected component carries.
    property bool hasPosition: false
    property bool hasOrientation: false

    // Row index of a parameter key, or -1. `commitParam` addresses rows; the profile
    // panel knows keys.
    function rowOfKey(key) {
        for (var i = 0; i < paramModel.count; ++i)
            if (paramModel.get(i).key === key) return i
        return -1
    }

    function refreshGeometry() {
        root.hasPosition = false; root.hasOrientation = false
        if (root.selectedComponent.length === 0) return
        var f = Julia.shell_component_geometry(root.selectedComponent).split("\t")
        if (f.length === 2) {
            root.hasPosition    = f[0] === "1"
            root.hasOrientation = f[1] === "1"
        }
    }

    function setGeometry(which, on) {
        var e = Julia.shell_set_component_geometry(root.selectedComponent, which, on)
        root.consoleChanged()
        if (e.length > 0 && e.charAt(0) === "!") { root.fitText = e; return }
        root.modelDirty = true
        refreshModel()
    }

    // The selected component's azimuthal modes, read back from the model. The panel keeps no
    // list of its own: a mode is two keys in the dict, and anything the panel remembered
    // separately could disagree with what the parser will build.
    function refreshAzModes() {
        azModel.clear()
        if (root.selectedComponent.length === 0) return
        var r = Julia.shell_az_modes(root.selectedComponent)
        if (r.length === 0) return
        var lines = r.split("\n")
        for (var i = 0; i < lines.length; ++i) {
            var f = lines[i].split("\t")
            if (f.length !== 3) continue
            azModel.append({ n: parseInt(f[0]),
                             amp: parseFloat(f[1]),
                             projang: parseFloat(f[2]) })
        }
    }

    // Add and remove go through Julia so the amp/projang pair stays a pair.
    function addAzMode() {
        var e = Julia.shell_add_az_mode(root.selectedComponent)
        root.consoleChanged()
        if (e.length > 0 && e.charAt(0) === "!") { root.fitText = e; return }
        root.modelDirty = true
        refreshModel()
    }

    function removeAzMode(order) {
        var e = Julia.shell_remove_az_mode(root.selectedComponent, order)
        root.consoleChanged()
        if (e.length > 0 && e.charAt(0) === "!") { root.fitText = e; return }
        root.modelDirty = true
        refreshModel()
    }

    // Edited here or edited in the parameter table is the same edit: both commit the same key
    // through commitParam, so a mode freed for fitting keeps its bounds and its fit index.
    function commitAzField(order, suffix, text) {
        var row = root.rowOfKey(root.selectedComponent + ",az " + suffix + order)
        if (row < 0) { refreshModel(); return }
        commitParam(row, "value", text)
    }

    // The selected component's limb-darkening law, or "" when it is not a limb-darkened disc.
    property string ldLawKey: ""
    property string ldLawLabel: ""

    function refreshLdLaw() {
        root.ldLawKey = ""; root.ldLawLabel = ""
        if (root.selectedComponent.length === 0) return
        var r = Julia.shell_ld_law(root.selectedComponent)
        if (r.length === 0) return
        var f = r.split("\t")
        root.ldLawKey = f[0]; root.ldLawLabel = f.length > 1 ? f[1] : f[0]
    }

    // Rows in `expr` mode. The Resolution view draws a DAG of $-references, and with no derived
    // parameter there are no references and no DAG -- so the count is what decides whether that
    // panel is worth the height it takes.
    property int nDerived: 0

    // The selected component's profile, when it has one. Empty for every other kind, which is
    // what the preview and the grid readout key off.
    property string profileKey: ""
    property string profileText: ""
    property string profileError: ""

    // One route for the profile string, whether it arrives from the editor's Apply, from focus
    // leaving it, or from the expression dialog on the `profile` row: `shell_set_param`'s "expr"
    // field writes the string and drops the row from the fit vector in one call.
    function commitProfile(text) {
        if (root.profileKey.length === 0 || text === root.profileText) return
        var err = Julia.shell_set_param(root.profileKey, "expr", text)
        if (err.length > 0 && err.charAt(0) === "!") { root.fitText = err; return }
        root.modelDirty = true
        refreshModel()
        root.consoleChanged()
    }

    function refreshProfile() {
        root.profileKey = ""
        root.profileError = ""
        root.profileText = ""
        root.profileParams = []
        if (root.selectedComponent.length === 0) return
        for (var i = 0; i < paramModel.count; ++i) {
            if (paramModel.get(i).key === root.selectedComponent + ",profile") {
                root.profileKey = paramModel.get(i).key
                break
            }
        }
        if (root.profileKey.length === 0) return

        for (var j = 0; j < paramModel.count; ++j)
            if (paramModel.get(j).key === root.profileKey)
                root.profileText = paramModel.get(j).expr

        var g = Julia.shell_profile_grid(root.selectedComponent)
        if (g.charAt(0) === "!") { root.profileError = g; return }
        var f = g.split("\t")
        if (f.length >= 4) {
            // The KEY, and whether it was halved: three of the four are diameters, so a grid
            // twice the expected size is otherwise a puzzle with nothing on screen to solve it.
            root.profileRMaxKey = f[0] + (f[0] === "r_max" ? "" : "/2")
            root.profileRMin = parseFloat(f[1])
            root.profileRMax = parseFloat(f[2])
            root.profileNr   = parseInt(f[3])
            root.profileEdgeFrac = f.length > 4 ? parseFloat(f[4]) : 0
        }

        // What `compile_profile` will treat as a parameter, split into what exists and what
        // writing the expression has asked for but not yet created.
        var found = [], want = [], miss = []
        for (var r = 0; r < paramModel.count; ++r) {
            var k = paramModel.get(r).key
            if (k.indexOf(root.selectedComponent + ",") === 0) {
                var suf = k.substring(root.selectedComponent.length + 1)
                if (root.profileText.indexOf("$" + suf) >= 0) found.push(suf)
            }
        }
        var mt = Julia.shell_profile_params(root.profileKey, root.profileText)
        if (mt.length > 0) {
            var ml = mt.split("\n")
            for (var q = 0; q < ml.length; ++q) {
                var mf = ml[q].split("\t")
                if (mf.length === 2) { want.push(mf[0] + " (missing)"); miss.push("$" + mf[0]) }
            }
        }
        root.profileParams = found.concat(want)
        root.profileMissing = miss

        // Drawing is Julia's: the curves go straight into the preview's Observables.
        var err = Julia.shell_profile_curves(root.selectedComponent)
        root.profileError = (err.length > 0 && err.charAt(0) === "!") ? err : ""
        root.profileDirty()
    }

    signal profileDirty()

    // Free-parameter names, for the grid's two axis pickers. Refreshed with the model, because
    // freeing or fixing a parameter changes what a grid can be run over.
    property var freeNames: []

    // The grid's step size on each axis, from the bounds and the step count. A user thinks in
    // "how fine is this", and the answer is the bound span divided by the steps — which the
    // panel can work out and they should not have to.
    function _spanOf(key) {
        for (var i = 0; i < paramModel.count; ++i) {
            var r = paramModel.get(i)
            if (r.key === key) return (isFinite(r.ub) && isFinite(r.lb)) ? (r.ub - r.lb) : NaN
        }
        return NaN
    }
    readonly property real gridStep1: _spanOf(gridP1.currentText) / Math.max(1, gridN.value - 1)
    readonly property real gridStep2: _spanOf(gridP2.currentText) / Math.max(1, gridN2.value - 1)

    // What the last grid fit mapped, or empty when the last fit was not one. The map panel
    // shows the chart only when there is a chart to show.
    // Everything the last fit printed, shown verbatim in the console under the results table.
    property string fitOutput: ""

    property string chi2MapP1: ""
    property string chi2MapP2: ""

    // Which nested sampler is available, and what to call it. Empty when neither extension is
    // loaded, which is what the entry is gated on.
    property string nestedBackend: ""
    property string nestedLabel: ""

    Component.onCompleted: {
        var nb = Julia.shell_nested_backend().split("\t")
        root.nestedBackend = nb[0]
        root.nestedLabel   = nb.length > 1 ? nb[1] : ""

        // Consecutive lines sharing a category become one menu entry with a sub-choice;
        // a category of one is offered as itself, since a menu of one is not a choice.
        var cats = []
        var txt = Julia.shell_component_kinds()
        if (txt.length > 0) {
            var lines = txt.split("\n")
            for (var i = 0; i < lines.length; ++i) {
                var f = lines[i].split("\t")
                if (f.length !== 5) continue
                var last = cats.length > 0 ? cats[cats.length - 1] : null
                if (last === null || last.key !== f[0])
                    cats.push({ key: f[0], label: f[1], subs: [] })
                cats[cats.length - 1].subs.push({ key: f[2], label: f[3], base: f[4] })
            }
        }
        root.componentKinds = cats

        var lw = []
        var lt = Julia.shell_ld_laws()
        if (lt.length > 0) {
            var ll = lt.split("\n")
            for (var m = 0; m < ll.length; ++m) {
                var lf = ll[m].split("\t")
                if (lf.length === 3) lw.push({ key: lf[0], label: lf[1], note: lf[2] })
            }
        }
        root.ldLaws = lw

        var tp = []
        var tt = Julia.shell_profile_templates()
        if (tt.length > 0) {
            var tl = tt.split("\n")
            for (var q = 0; q < tl.length; ++q) {
                var tf = tl[q].split("\t")
                if (tf.length === 3) tp.push({ key: tf[0], label: tf[1], note: tf[2] })
            }
        }
        root.profileTemplates = tp
    }

    readonly property int nFree: {
        var n = 0
        for (var i = 0; i < paramModel.count; ++i)
            if (paramModel.get(i).mode === "free") n++
        return n
    }

    // ── constraints (relations bounds cannot express) ────────────────────────
    ListModel {
        id: constraintModel
        // roles: lhs, op, rhs, tol, satisfied
        // `satisfied` is set by `shell_check_constraints`; `constraintsAsLines` sends them
        // to the fitter.
    }

    // ── Gaussian priors ──────────────────────────────────────────────────────
    ListModel {
        id: priorModel
        // roles: expr, target, sigma
        // Sent to the fitter by `priorsAsLines` as the `priors` vector of tuples.
    }

    // ── model → image ────────────────────────────────────────────────────────
    property int    imgNx: 128
    property real   imgPixsize: 0.1
    property real   imgWl: 1.6e-6            // metres, matching $WL
    property real   imgMjd: 0.0
    property bool   hasRender: false
    property string renderText: ""

    // Whether λ and MJD mean anything for THIS model, so the controls can be dead when they
    // would change nothing. Julia decides, by looking for $WL and $MJD in the expressions.
    property bool dependsWl: false
    property bool dependsMjd: false
    // The wavelength a chromatic value in the table is shown at. Named in the banner: a
    // number with no wavelength beside it is one arbitrary entry of a vector presented as
    // though it were the parameter.
    property real displayWl: NaN

    function modelAsLines() {
        var out = []
        for (var i = 0; i < paramModel.count; ++i) {
            var r = paramModel.get(i)
            out.push(r.key + "=" + (r.mode === "expr" ? r.expr : r.value))
        }
        return out.join("\n")
    }

    function renderModel() {
        renderText = Julia.shell_model_image(modelAsLines(), imgNx, imgPixsize,
                                             dependsWl ? imgWl : 0, dependsMjd ? imgMjd : 0)
        hasRender = renderText.indexOf("!") !== 0
        modelArea.update()
        consoleChanged()
    }


    // One picker for every plot area on this tab; `which` says which one asked for it.
    FilePicker {
        id: savePngDialog
        property string which: ""
        property real pxw: 1200
        property real pxh: 900
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Save plot as PNG"
        saveMode: true
        defaultSuffix: "png"
        filters: [{ label: "PNG images (*.png)", patterns: "*.png" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) {
            Julia.shell_save_figure(savePngDialog.which, path,
                                    savePngDialog.pxw, savePngDialog.pxh, false)
            root.consoleChanged()
        }
    }

    // Twice the size it has on screen, so the file is usable in a talk rather than being a
    // screenshot of a panel. Capped because the figure is rendered into a real framebuffer and
    // a software GL stack refuses the very large ones.
    function savePng(which, area) {
        savePngDialog.which = which
        savePngDialog.pxw = Math.min(2400, Math.max(640, area.width * 2))
        savePngDialog.pxh = Math.min(1800, Math.max(480, area.height * 2))
        savePngDialog.openAt("")
    }

    signal consoleChanged()

    // ── running a fit ────────────────────────────────────────────────────────
    //
    // The panel holds strings and the fitters take dicts, bound dicts, constraint objects and
    // a seven-element weight vector; the conversion lives in Julia so this file never has to
    // know which fitter wants which shape.
    //
    // Deferred by a tick, as the reconstruction is: the call holds Qt's thread, so the status
    // has to be on screen before it starts.
    property string fitText: ""

    function freeAsLines() {
        var rows = []
        for (var i = 0; i < paramModel.count; ++i) {
            var r = paramModel.get(i)
            if (r.mode !== "free") continue
            rows.push({ key: r.key, lb: r.lb, ub: r.ub, idx: r.fitindex })
        }
        // x[i] ↔ free[i], so the rows go in fit-vector order and not in table order.
        rows.sort(function (a, b) { return a.idx - b.idx })
        var out = []
        for (var k = 0; k < rows.length; ++k)
            out.push(rows[k].key + "\t" + rows[k].lb + "\t" + rows[k].ub)
        return out.join("\n")
    }

    function constraintsAsLines() {
        var out = []
        for (var i = 0; i < constraintModel.count; ++i) {
            var c = constraintModel.get(i)
            if (!c.lhs || c.lhs.length === 0) continue
            out.push(c.lhs + "\t" + c.op + "\t" + c.rhs + "\t" + c.tol)
        }
        return out.join("\n")
    }

    function priorsAsLines() {
        var out = []
        for (var i = 0; i < priorModel.count; ++i) {
            var p = priorModel.get(i)
            if (!p.expr || p.expr.length === 0) continue
            out.push(p.expr + "\t" + p.target + "\t" + p.sigma)
        }
        return out.join("\n")
    }

    function startFit() {
        if (!canRun) return
        running = true
        fitText = "fitting…"
        fitTimer.restart()
    }

    Timer {
        id: fitTimer
        interval: 120; repeat: false
        onTriggered: {
            root.fitText = Julia.shell_fit_model(
                root.modelAsLines(), root.freeAsLines(),
                root.constraintsAsLines(), root.priorsAsLines(),
                root.useV2, root.useT3amp, root.useT3phi,
                root.useCvis, root.useFlux, root.useDiffvis,
                root.optimiser, root.maxeval,
                gridP1.currentText, gridP2.currentText, gridN.value, gridN2.value)
            root.fitOutput = ""
            if (root.fitText.length > 0 && root.fitText.charAt(0) === "!") {
                root.running = false
                root.consoleChanged()
                return
            }
            // The call starts the fit on a worker and returns; this follows it, so the
            // optimiser's trace appears while it runs rather than after.
            fitJobTimer.start()
        }
    }

    Timer {
        id: fitJobTimer
        interval: 200
        repeat: true
        running: false
        onTriggered: {
            var f = Julia.shell_job_poll().split("\t")
            if (f[0] === "running") { root.fitOutput = f[1]; return }
            stop()
            root.running = false
            root.fitOutput = Julia.shell_fit_output()
            if (f[0] === "done") root.fitText = f[1]
            root.refreshFits()
            root.refreshChi2Map()
            root.consoleChanged()
        }
    }

    function refreshFits() {
        fitModel.clear()
        var rows = Julia.shell_fit_rows()
        if (rows.length === 0) return
        var lines = rows.split("\n")
        for (var i = 0; i < lines.length; ++i) {
            var f = lines[i].split("\t")
            if (f.length === 10)
                fitModel.append({ label: f[0], optimiser: f[1],
                                  chi2rStart: parseFloat(f[2]), chi2r: parseFloat(f[3]),
                                  work: f[4],
                                  ndof: parseInt(f[5]), nfree: parseInt(f[6]),
                                  ret: "", chi2: 0.0,
                                  aic: parseFloat(f[7]), bic: parseFloat(f[8]),
                                  params: f[9] })
        }
    }

    // Write the fitted values back into the table, so what is shown is what was fitted and the
    // model can then be saved or rendered as one.
    // Each value goes to Julia, not into the table. Writing the row alone would leave the
    // model still holding the starting guess while the table showed the fitted number -- and
    // the chi2 readout, the render and the next fit would all be of the value on screen.
    function adoptFit() {
        var vals = Julia.shell_fit_values()
        if (vals.length === 0) return
        var lines = vals.split("\n")
        for (var i = 0; i < lines.length; ++i) {
            var j = lines[i].lastIndexOf("=")
            if (j < 0) continue
            var key = lines[i].substring(0, j)
            var v   = parseFloat(lines[i].substring(j + 1))
            if (isNaN(v)) continue
            var err = Julia.shell_set_param(key, "value", v)
            if (err.length > 0 && err.charAt(0) === "!") root.fitText = err
        }
        modelDirty = true
        refreshModel()
        consoleChanged()
    }

    // ── live feedback (§3.6) ─────────────────────────────────────────────────
    property real  currentChi2r: NaN            // chi2_flat, re-run after every edit
    property int   currentNdof: 0
    property var   validationWarnings: []       // display_model's own checks, rendered inline

    // ── observables (§1.8) ───────────────────────────────────────────────────
    //
    // Model fitting takes a SEVEN-element weight vector,
    // [V2, T3amp, T3phi, visamp, visphi, flux, diffphase] — unlike imaging's three. `cvis`
    // drives visamp AND visphi together, which is why it is one box and not two.
    property bool haveCvis:    false
    property bool haveV2:      false
    property bool haveT3amp:   false
    property bool haveT3phi:   false
    property bool haveDiffvis: false
    property bool haveFlux:    false

    readonly property bool useCvis:    cvisBox.checked
    readonly property bool useV2:      v2Box.checked
    readonly property bool useT3amp:   t3ampBox.checked
    readonly property bool useT3phi:   t3phiBox.checked
    readonly property bool useDiffvis: diffvisBox.checked
    readonly property bool useFlux:    fluxBox.checked
    readonly property bool anyObservable:
        useCvis || useV2 || useT3amp || useT3phi || useDiffvis || useFlux

    // ── optimiser and uncertainty: two separate, composable selectors (§3.7, §3.8) ──
    //
    // Bootstrap is not a fifth optimiser. It is a wrapper that runs an inner fitter many
    // times, so presenting it beside NLopt and LM would mis-model what it is.
    // lsqfit by default: it is the one reached for most often, and it is the only fitter that
    // returns a covariance, so the default path is also the one that can report an uncertainty.
    property string optimiser: "lsqfit"
    property string uncertainty: "none"         // none | jacobian | bootstrap | posterior

    property int  maxeval: 2000
    property real ftolRel: 1e-8
    property real xtolRel: 1e-6
    property int  maxIter: 200

    property int    nboot: 200
    property string bootMode: "replacement"     // replacement | halfsample | weights | pmoired
    property string bootGranularity: "config"
    property int    bootSeed: 1                 // set by default: an unseeded run is not reproducible
    property bool   bootThreaded: true

    property int  nlive: 400
    property bool useStepsampler: false
    property int  nsteps: 400
    property real fracRemain: 0.001

    property bool running: false
    property string statusText: "idle"

    // Nested sampling makes the bounds the prior, so every free parameter needs finite ones
    // and both samplers error without them. The entry is therefore disabled with
    // the reason shown rather than after a long run. `atbound` is a different thing and does
    // not gate anything.
    readonly property bool allBoundsFinite: {
        for (var i = 0; i < paramModel.count; ++i) {
            var r = paramModel.get(i)
            if (r.mode !== "free") continue
            if (!isFinite(r.lb) || !isFinite(r.ub)) return false
        }
        return true
    }

    // Julia pushes this after a load. Without it "no observables selected" is what an empty
    // session reports, which sends the user looking at tick boxes instead of at the File menu.
    property bool hasDataset: false

    readonly property string blockedReason:
          !hasDataset        ? "no dataset — open an OIFITS first"
        : nFree === 0        ? "no free parameters"
        : !anyObservable     ? "no observables selected"
        : running            ? "a fit is running"
        : ""
    readonly property bool canRun: blockedReason === ""

    // ── fit history (§3.9) ───────────────────────────────────────────────────
    ListModel {
        id: fitModel
        // roles: label, optimiser, chi2, chi2rStart, chi2r, work, ndof, nfree, ret,
        //        aic, bic, params. Roles are defined by the first append, so every
        //        append must carry all of them.
        // One row per completed fit, so competing models compare side by side.
    }

    function fmt(v, d) { return (v === undefined || v === null || isNaN(v)) ? "—" : v.toFixed(d) }

    // chi2 spans orders of magnitude: a starting model can score millions where the fitted one
    // scores single figures, and four decimals of 2283398.5050 is eleven characters of column
    // for two that matter. Fixed where it is readable, exponential where it is not.
    function fmtChi2(v) {
        if (v === undefined || v === null || isNaN(v)) return "—"
        return Math.abs(v) >= 1e5 ? v.toExponential(3) : v.toFixed(4)
    }

    // The `$`-references an expression makes, in the order they appear and without repeats.
    // A reference is `$name,suffix` for a component parameter or `$name` for a global, and a
    // component name cannot contain a comma, so nothing more is needed to know where one ends.
    function refsOf(expr) {
        if (!expr || expr.length === 0) return []
        var m = expr.match(/\$[A-Za-z_][A-Za-z0-9_]*(?:,[A-Za-z_][A-Za-z0-9_ ]*)?/g)
        if (m === null) return []
        var out = []
        for (var i = 0; i < m.length; ++i) {
            var r = m[i].replace(/\s+$/, "")
            if (out.indexOf(r) < 0) out.push(r)
        }
        return out
    }

    // The component a reference or a key belongs to, or "" for a global, for $WL/$MJD, and for
    // a name no component has. The caller checks before selecting: an expression is free to
    // name something that is not there, and the unrecognised list is what reports that.
    function componentOf(ref) {
        var t = ref.charAt(0) === "$" ? ref.substring(1) : ref
        var i = t.indexOf(",")
        if (i < 0) return ""
        var name = t.substring(0, i)
        for (var k = 0; k < componentModel.count; ++k)
            if (componentModel.get(k).name === name) return name
        return ""
    }

    // ═════════════════════════════════════════════════════════════════════════
    ColumnLayout {
        anchors.fill: parent
        anchors.margins: dp(8)
        spacing: dp(8)

        // ── toolbar: the model as a document ─────────────────────────────────
        RowLayout {
            Layout.fillWidth: true
            spacing: dp(6)

            Label { text: "Model"; font.bold: true }
            Label {
                text: root.modelName + (root.modelDirty ? " •" : "")
                color: root.modelDirty ? root.cWarn : "#444"
                elide: Text.ElideMiddle
                Layout.maximumWidth: dp(220)
            }

            ToolSeparator {}

            // A model dict does not describe a fit on its own: the free list, the bounds, the
            // constraints and the priors all change the answer. The TOML file carries all five.
            Button {
                text: "Open…"
                // Open where the models are, not where the last dialog happened to be: the
                // shipped starting points live in demos/models, and a model file is what this
                // button is for. A path already in use wins, so re-opening returns to it.
                onClicked: openModelDialog.openAt(root.modelPath.length > 0 ? root.modelPath
                                                  : Julia.picker_examples("model"))
            }
            Button {
                text: "Save…"
                enabled: paramModel.count > 0
                onClicked: saveModelDialog.openAt(root.modelPath.length > 0 ? root.modelPath
                                                  : Julia.picker_examples("model"))
            }

            ToolSeparator {}

            Button {
                text: "Import PMOIRED…"
                onClicked: importPmoiredDialog.openAt(Julia.picker_examples("pmoired"))
            }
            Button {
                text: "Export PMOIRED…"
                enabled: paramModel.count > 0
                // Export warns rather than failing when the model uses something PMOIRED
                // cannot represent — the OITOOLS-only geometries. Writing a model that means
                // something else on the other side is the failure worth preventing.
                onClicked: exportPmoiredDialog.openAt("")
            }

            Item { Layout.fillWidth: true }

            // χ² updated on every edit (§3.6) — cheap for analytic components, still
            // interactive for Hankel ones at nr = 100.
            Label {
                text: "χ²ᵣ " + root.fmt(root.currentChi2r, 3)
                font.bold: true
                color: isNaN(root.currentChi2r) ? "#999" : "#222"
            }
            Label {
                text: root.nFree + " free"
                color: "#666"
                font.pointSize: root.pt(root.baseFontPt - 1)
            }
        }

        // ── broadcasting and unrecognised-key banners ────────────────────────
        //
        // Both are silent in the library and both change what is being fitted, so they are
        // banners rather than something to find in a panel.
        Rectangle {
            Layout.fillWidth: true
            visible: root.unrecognisedKeys.length > 0
            color: "#fdeaea"
            border.color: root.cBad
            radius: dp(3)
            implicitHeight: unrecLabel.implicitHeight + dp(10)
            Label {
                id: unrecLabel
                anchors.fill: parent
                anchors.margins: dp(5)
                wrapMode: Text.WordWrap
                color: root.cBad
                font.pointSize: root.pt(root.baseFontPt - 1)
                text: "Keys the parser ignored: " + root.unrecognisedKeys.join(", ") +
                      " — these do not error, they change what the component is."
            }
        }

        Rectangle {
            Layout.fillWidth: true
            visible: root.broadcasting
            color: "#fff6e0"
            border.color: root.cWarn
            radius: dp(3)
            implicitHeight: bcLabel.implicitHeight + dp(10)
            Label {
                id: bcLabel
                anchors.fill: parent
                anchors.margins: dp(5)
                wrapMode: Text.WordWrap
                color: root.cWarn
                font.pointSize: root.pt(root.baseFontPt - 1)
                // If ANY derived expression references $WL or $MJD the resolver broadcasts
                // ALL of them, so every parameter becomes a per-uv-point vector. That is a
                // behavioural switch and the user should see it flip.
                // The table has to show a single number per row, so it shows the value at one
                // wavelength; naming it here is what stops that number being read as the
                // parameter itself.
                text: "Broadcasting is on: one chromatic expression makes every parameter a " +
                      "per-uv-point vector, not just the chromatic one." +
                      (isNaN(root.displayWl)
                        ? "  Load a dataset to see them at a wavelength."
                        : "  Values in the table are shown at λ = " +
                          (root.displayWl * 1e6).toFixed(3) + " µm.")
            }
        }

        // ── the three editors, side by side ──────────────────────────────────
        // ── the three things this perspective does ───────────────────────────
        //
        // Creating a model, seeing it, and fitting it are separate activities that share
        // one model, not one activity with three panels. Splitting them keeps each page
        // legible; the toolbar and the banners above stay common, because they describe
        // the model itself rather than what is being done to it.
        SplitView {
            Layout.fillWidth: true
            Layout.fillHeight: true

            ColumnLayout {
                SplitView.preferredWidth: root.dp(620)
                SplitView.minimumWidth: root.dp(420)
                spacing: root.dp(6)

                // Component cards. Kind is inferred from the geometry key exactly as
                // `_identify_kind` does; the card names the deciding key so a surprising
                // kind is traceable to the key that caused it.
                //
                // "+ component" sits outside the Flow, and first. Inside it, the button
                // trailed the cards and moved every time one was added or removed, so the
                // place you click to add a second component is not where you clicked to add
                // the first.
                RowLayout {
                    Layout.fillWidth: true
                    spacing: root.dp(6)

                    Button {
                        text: "+ component"
                        Layout.alignment: Qt.AlignTop
                        onClicked: {
                            addComponentDialog.suggestName()
                            addComponentDialog.open()
                        }
                    }

                    Flow {
                        Layout.fillWidth: true
                        spacing: root.dp(6)

                        Repeater {
                            model: componentModel
                            delegate: Rectangle {
                                width: cardCol.implicitWidth + root.dp(16)
                                height: cardCol.implicitHeight + root.dp(10)
                                radius: root.dp(4)
                                color: model.name === root.selectedComponent ? "#e8f0fe" : root.cPanel
                                border.color: model.nunrecognised > 0 ? root.cBad : root.cLine

                                // Declared BEFORE the content: later siblings sit on top and are
                                // offered the click first, so the − button is reachable instead of
                                // being swallowed by the card's own select-on-click.
                                MouseArea {
                                    anchors.fill: parent
                                    onClicked: root.selectedComponent = model.name
                                }
                                RowLayout {
                                    id: cardCol
                                    anchors.centerIn: parent
                                    spacing: root.dp(6)

                                    ColumnLayout {
                                        spacing: 0
                                        Label {
                                            text: model.name
                                            font.bold: true
                                            visible: !nameEdit.visible
                                        }
                                        // The name is the prefix of every key the component owns
                                        // and of every `$name,suffix` an expression refers to, so
                                        // it is worth being able to change without rebuilding the
                                        // component around it.
                                        TextField {
                                            id: nameEdit
                                            visible: false
                                            selectByMouse: true
                                            font.bold: true
                                            implicitWidth: root.dp(110)
                                            // visible is cleared BEFORE the rename: the rename
                                            // refreshes the model, which destroys this delegate,
                                            // and assigning to it afterwards would be a write to
                                            // an object that no longer exists.
                                            onAccepted: {
                                                var was = model.name, want = text
                                                visible = false
                                                root.renameComponent(was, want)
                                            }
                                            onActiveFocusChanged: if (!activeFocus) visible = false
                                            Keys.onEscapePressed: visible = false
                                        }
                                        Label {
                                            text: model.kind + " ← " + model.geometryKey
                                            color: "#666"
                                            font.pointSize: root.pt(root.baseFontPt - 2)
                                        }
                                    }
                                    ToolButton {
                                        text: "✎"
                                        implicitWidth: root.dp(22)
                                        implicitHeight: root.dp(22)
                                        ToolTip.visible: hovered
                                        ToolTip.text: "rename " + model.name
                                        onClicked: {
                                            nameEdit.text = model.name
                                            nameEdit.visible = true
                                            nameEdit.forceActiveFocus()
                                            nameEdit.selectAll()
                                        }
                                    }
                                    ToolButton {
                                        text: "−"
                                        implicitWidth: root.dp(22)
                                        implicitHeight: root.dp(22)
                                        ToolTip.visible: hovered
                                        ToolTip.text: "remove " + model.name + " and its parameters"
                                        onClicked: root.removeComponent(model.name)
                                    }
                                }
                            }
                        }
                        }
                }

                // ── the parameter table ──────────────────────────────────────
                Rectangle {
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    color: "white"
                    border.color: root.cLine
                    radius: root.dp(3)

                    ColumnLayout {
                        anchors.fill: parent
                        anchors.margins: root.dp(1)
                        spacing: 0

                        // header
                        Rectangle {
                            Layout.fillWidth: true
                            implicitHeight: root.dp(24)
                            color: root.cPanel
                            RowLayout {
                                anchors.fill: parent
                                anchors.leftMargin: root.dp(6)
                                anchors.rightMargin: root.dp(6)
                                spacing: root.dp(4)
                                Label { text: "parameter"; Layout.preferredWidth: root.dp(150)
                                        font.pointSize: root.pt(root.baseFontPt - 2); color: "#555" }
                                Label { text: "mode"; Layout.preferredWidth: root.dp(78)
                                        font.pointSize: root.pt(root.baseFontPt - 2); color: "#555" }
                                Label { text: "value"; Layout.preferredWidth: root.dp(96)
                                        font.pointSize: root.pt(root.baseFontPt - 2); color: "#555" }
                                Label { text: "lower"; Layout.preferredWidth: root.dp(76)
                                        font.pointSize: root.pt(root.baseFontPt - 2); color: "#555" }
                                Label { text: "upper"; Layout.preferredWidth: root.dp(76)
                                        font.pointSize: root.pt(root.baseFontPt - 2); color: "#555" }
                                Label { text: "range"; Layout.fillWidth: true
                                        font.pointSize: root.pt(root.baseFontPt - 2); color: "#555" }
                                Label { text: "x[i]"; Layout.preferredWidth: root.dp(34)
                                        horizontalAlignment: Text.AlignRight
                                        font.pointSize: root.pt(root.baseFontPt - 2); color: "#555" }
                            }
                        }

                        ListView {
                            id: paramList
                            Layout.fillWidth: true
                            Layout.fillHeight: true
                            clip: true
                            model: paramModel
                            ScrollBar.vertical: ScrollBar {}

                            // Section headers group by component, with globals first — the
                            // parser already groups bare keys under `__global__`, and globals
                            // are where linked parameters live.
                            section.property: "comp"
                            section.criteria: ViewSection.FullString
                            section.delegate: Rectangle {
                                width: paramList.width
                                implicitHeight: root.dp(22)
                                color: "#eef1f4"
                                Label {
                                    anchors.verticalCenter: parent.verticalCenter
                                    x: root.dp(6)
                                    text: section === "__global__" ? "globals" : section
                                    font.bold: true
                                    font.pointSize: root.pt(root.baseFontPt - 1)
                                    color: "#333"
                                }
                            }

                            delegate: Rectangle {
                                id: rowItem
                                width: paramList.width
                                implicitHeight: root.dp(26)
                                color: index === paramList.currentIndex ? "#eaf2ff"
                                     : (index % 2 ? "#fbfbfb" : "white")

                                // Captured here, not read as `model.x` further down: inside a
                                // ComboBox `model` means the ComboBox's own list, which
                                // silently shadows the delegate's row.
                                readonly property int    rIndex:   index
                                readonly property string rParam:   model.param
                                readonly property string rMode:    model.mode
                                readonly property real   rValue:   model.value
                                readonly property real   rLb:      model.lb
                                readonly property real   rUb:      model.ub
                                readonly property int    rFitIndex: model.fitindex
                                readonly property bool   rAtBound: model.atbound

                                MouseArea {
                                    anchors.fill: parent
                                    onClicked: paramList.currentIndex = rowItem.rIndex
                                }

                                RowLayout {
                                    anchors.fill: parent
                                    anchors.leftMargin: root.dp(6)
                                    anchors.rightMargin: root.dp(6)
                                    spacing: root.dp(4)

                                    Label {
                                        Layout.preferredWidth: root.dp(150)
                                        Layout.maximumWidth: root.dp(150)
                                        text: rowItem.rParam
                                        elide: Text.ElideRight
                                        // Free in green, derived in violet, fixed in plain --
                                        // the same three colours the mode selector uses, so the
                                        // state of a row is readable without reading its mode.
                                        // At-bound still wins: it is a warning about a free
                                        // parameter and would be invisible in the free colour.
                                        color: rowItem.rAtBound ? root.cWarn
                                             : rowItem.rMode === "free" ? root.cFree
                                             : rowItem.rMode === "expr" ? root.cExpr : "#222"
                                        font.bold: rowItem.rMode === "free"
                                    }

                                    // The three-way selector. Not a checkbox: a tick cannot
                                    // express "this is an expression", and free/derived are
                                    // mutually exclusive rather than independent.
                                    ComboBox {
                                        id: modeBox
                                        Layout.preferredWidth: root.dp(78)
                                        Layout.maximumWidth: root.dp(78)
                                        implicitHeight: root.dp(22)
                                        font.pointSize: root.pt(root.baseFontPt - 2)
                                        // "derived" is what this panel, the inspector and the
                                        // documentation call the state; "expr" is what the shell
                                        // and `model_rows` call it. The label follows the prose
                                        // and the wire value follows the parser, so the one word
                                        // the user has to find is the one they have been reading.
                                        model: ["fixed", "free", "derived"]
                                        currentIndex: rowItem.rMode === "free" ? 1
                                                    : rowItem.rMode === "expr" ? 2 : 0
                                        onActivated: root.requestModeChange(rowItem.rIndex,
                                                         currentIndex === 2 ? "expr" : currentText)
                                    }

                                    // Expr rows show the expression; its value is a computed
                                    // read-only display, since editing it would mean nothing.
                                    TextField {
                                        Layout.preferredWidth: root.dp(96)
                                        Layout.maximumWidth: root.dp(96)
                                        implicitHeight: root.dp(22)
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                        readOnly: rowItem.rMode === "expr"
                                        color: rowItem.rMode === "expr" ? root.cExpr
                                             : rowItem.rMode === "free" ? root.cFree : "#222"
                                        text: isNaN(rowItem.rValue) ? "—" : rowItem.rValue.toFixed(4)
                                        horizontalAlignment: Text.AlignRight
                                        onEditingFinished: root.commitParam(rowItem.rIndex, "value", text)
                                    }

                                    TextField {
                                        Layout.preferredWidth: root.dp(76)
                                        Layout.maximumWidth: root.dp(76)
                                        implicitHeight: root.dp(22)
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                        enabled: rowItem.rMode === "free"
                                        opacity: enabled ? 1 : 0.35
                                        text: rowItem.rMode === "free" ? root.fmt(rowItem.rLb, 3) : ""
                                        horizontalAlignment: Text.AlignRight
                                        onEditingFinished: root.commitParam(rowItem.rIndex, "lb", text)
                                    }
                                    TextField {
                                        Layout.preferredWidth: root.dp(76)
                                        Layout.maximumWidth: root.dp(76)
                                        implicitHeight: root.dp(22)
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                        enabled: rowItem.rMode === "free"
                                        opacity: enabled ? 1 : 0.35
                                        text: rowItem.rMode === "free" ? root.fmt(rowItem.rUb, 3) : ""
                                        horizontalAlignment: Text.AlignRight
                                        onEditingFinished: root.commitParam(rowItem.rIndex, "ub", text)
                                    }

                                    // Where the value sits inside its box. `default_bounds` is
                                    // deliberately generous, so a span-relative test would call
                                    // every small diameter pinned; what is worth showing is the
                                    // value landing ON a bound, which is what `atbound` means.
                                    // Always present, even on a fixed row that has no range to
                                    // show. A RowLayout ignores INVISIBLE items, so hiding this
                                    // one removed the only `fillWidth` item from the row and the
                                    // leftover width was redistributed across the fields --
                                    // which is why fixed rows sat wider and further right than
                                    // free ones. It stays in the layout and draws nothing.
                                    Rectangle {
                                        readonly property bool showRange: rowItem.rMode === "free"
                                        Layout.fillWidth: true
                                        implicitHeight: root.dp(6)
                                        color: showRange ? "#e9e9e9" : "transparent"
                                        radius: root.dp(3)

                                        Rectangle {
                                            visible: parent.showRange
                                                     && isFinite(rowItem.rLb) && isFinite(rowItem.rUb)
                                                     && rowItem.rUb > rowItem.rLb
                                            width: root.dp(6)
                                            height: parent.height
                                            radius: width / 2
                                            color: rowItem.rAtBound ? root.cWarn : root.cFree
                                            x: {
                                                if (!(rowItem.rUb > rowItem.rLb)) return 0
                                                var t = (rowItem.rValue - rowItem.rLb)
                                                        / (rowItem.rUb - rowItem.rLb)
                                                t = Math.max(0, Math.min(1, t))
                                                return t * (parent.width - width)
                                            }
                                        }
                                        // An unbounded free parameter has no position to show,
                                        // and saying so beats drawing a marker that means nothing.
                                        Label {
                                            anchors.centerIn: parent
                                            visible: parent.showRange
                                                     && !(isFinite(rowItem.rLb) && isFinite(rowItem.rUb))
                                            text: "unbounded"
                                            color: "#999"
                                            font.pointSize: root.pt(root.baseFontPt - 3)
                                        }
                                    }

                                    Label {
                                        Layout.preferredWidth: root.dp(34)
                                        horizontalAlignment: Text.AlignRight
                                        // x[i] ↔ list_free_params[i] holds throughout OITOOLS,
                                        // so the free rows in order ARE the parameter vector
                                        // every optimiser sees. Showing i makes an exported
                                        // script legible.
                                        text: rowItem.rFitIndex > 0 ? String(rowItem.rFitIndex) : ""
                                        color: root.cFree
                                        font.pointSize: root.pt(root.baseFontPt - 2)
                                    }
                                }
                            }
                        }

                        // ── row actions ──────────────────────────────────────
                        Rectangle {
                            Layout.fillWidth: true
                            implicitHeight: root.dp(32)
                            color: root.cPanel
                            RowLayout {
                                anchors.fill: parent
                                anchors.leftMargin: root.dp(6)
                                anchors.rightMargin: root.dp(6)
                                spacing: root.dp(6)

                                Button {
                                    text: "Default bounds"
                                    implicitHeight: root.dp(24)
                                    enabled: root.nFree > 0
                                    // With a dataset loaded the angular-size ceiling comes from
                                    // the coverage itself — 2 λ/B_min, the largest scale the
                                    // shortest baseline senses — rather than from a constant.
                                    ToolTip.visible: hovered
                                    ToolTip.text: root.hasDataset
                                                  ? "bounds for every free parameter, sizes scaled from the uv coverage"
                                                  : "bounds for every free parameter (load data to scale sizes from the coverage)"
                                    onClicked: {
                                        var err = Julia.shell_default_bounds()
                                        root.consoleChanged()
                                        if (err.length > 0 && err.charAt(0) === "!") {
                                            root.fitText = err; return
                                        }
                                        root.modelDirty = true
                                        root.refreshModel()
                                    }
                                }
                                Button {
                                    text: "Link"
                                    implicitHeight: root.dp(24)
                                    // Selecting two parameters and linking creates a global and
                                    // rewrites both rows to reference it ("$PA"). Tying two
                                    // position angles together becomes one gesture instead of
                                    // an expression-editing exercise. Unlinking inlines the
                                    // current value back into each.
                                    // NOT WIRED: needs multi-select in the parameter table, which
                                    // the table does not have, plus a shell call to create the
                                    // global and rewrite both rows to reference it.
                                }
                                Item { Layout.fillWidth: true }
                                Label {
                                    text: root.nFree + " in the fit vector"
                                    color: "#666"
                                    font.pointSize: root.pt(root.baseFontPt - 2)
                                }
                            }
                        }
                    }
                }
            }

            ColumnLayout {
                SplitView.fillWidth: true
                SplitView.minimumWidth: root.dp(340)
                spacing: root.dp(6)

                // The three views of a model, at the level the inspector tabs used to sit.
                // The parameter table beside them is NOT part of the stack: it is the model,
                // and it stays put whichever view is showing.
                TabBar {
                    id: modelTabs
                    Layout.fillWidth: true
                    TabButton { text: "Model Creation" }
                    TabButton { text: "Model visualization" }
                    TabButton { text: "Model Optimization" }
                }

                StackLayout {
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    currentIndex: modelTabs.currentIndex

                        // ── Model Creation ───────────────────────────────────────
                        //
                        // Inspector, Components and Constraints on ONE surface rather than
                        // behind three tabs. They describe the same model from different
                        // angles, and the misparse the inspector exists to catch is only
                        // visible when the component that caused it is in the same view.
                        //
                        // Each section keeps its own ColumnLayout: flattening them into the
                        // shared one made every `spacing:` a second assignment to the same
                        // property, which QML rejects outright.
                        ScrollView {
                            clip: true
                            ColumnLayout {
                                width: root.width * 0.4
                                spacing: root.dp(8)

                                ColumnLayout {
                                    spacing: root.dp(8)

                                    // Warnings first. They are the reason to open this panel at
                                    // all, and under three views of a model that parsed cleanly
                                    // they are read last or not at all.
                                    //
                                    // `display_model` performs exactly these checks -- value<lb,
                                    // value>ub, lb>=ub, flux fractions not summing to 1 -- so what
                                    // is rendered here are its results, not a second implementation
                                    // of them that could disagree with the fit.
                                    Label {
                                        visible: root.validationWarnings.length > 0
                                        text: "Validation"
                                        font.bold: true
                                    }
                                    Repeater {
                                        model: root.validationWarnings
                                        delegate: Label {
                                            Layout.fillWidth: true
                                            wrapMode: Text.WordWrap
                                            color: root.cWarn
                                            font.pointSize: root.pt(root.baseFontPt - 1)
                                            text: "• " + modelData
                                        }
                                    }
                                    // The checks above are `display_model`'s own — value<lb,
                                    // value>ub, lb>=ub, flux fractions not summing to 1 — rendered
                                    // through `shell_model_warnings` rather than reimplemented.

                                    // ── Dependencies: what is computed from what ──────
                                    //
                                    // Hidden until something derives from something else: with
                                    // no `$`-reference there is nothing to resolve, and an empty
                                    // box is height taken from the panels that do have something
                                    // to say.
                                    //
                                    // One row per derived parameter, in the resolver's own order
                                    // -- globals, then fixed, then derived, topologically sorted
                                    // -- with each `$`-reference as a chip. The chips are the
                                    // point: a reference is a parameter somewhere else in this
                                    // model, and clicking one selects the component it belongs to
                                    // instead of leaving you to find it in the table.
                                    // "Dependencies", not "Resolution". This is the resolver's
                                    // evaluation order -- globals, then fixed, then derived,
                                    // topologically sorted -- but on a panel read by people who
                                    // measure angular resolution for a living, that word names
                                    // the wrong thing entirely.
                                    Label {
                                        text: "Dependencies"; font.bold: true
                                        visible: root.nDerived > 0
                                    }
                                    InfoTip {
                                        visible: root.nDerived > 0 && root.broadcasting
                                        tip: "One expression uses $WL or $MJD, so the resolver broadcasts every parameter: each is a vector over the uv points, not a number."
                                    }
                                    Repeater {
                                        model: paramModel
                                        delegate: Rectangle {
                                            id: resRow
                                            property string keyText: model.key
                                            property string exprText: model.expr
                                            visible: model.mode === "expr"
                                            Layout.fillWidth: true
                                            implicitHeight: visible ? rCol.implicitHeight + root.dp(10) : 0
                                            color: "white"
                                            border.color: root.cLine
                                            radius: root.dp(3)
                                            ColumnLayout {
                                                id: rCol
                                                x: root.dp(6); y: root.dp(5)
                                                width: parent.width - root.dp(12)
                                                spacing: root.dp(2)
                                                Label {
                                                    text: resRow.keyText + " = " + resRow.exprText
                                                    color: root.cExpr
                                                    font.family: "monospace"
                                                    font.pointSize: root.pt(root.baseFontPt - 1)
                                                    elide: Text.ElideRight
                                                    Layout.fillWidth: true
                                                }
                                                Flow {
                                                    Layout.fillWidth: true
                                                    spacing: root.dp(4)
                                                    Repeater {
                                                        model: root.refsOf(resRow.exprText)
                                                        delegate: Rectangle {
                                                            width: refLabel.implicitWidth + root.dp(10)
                                                            height: refLabel.implicitHeight + root.dp(4)
                                                            radius: root.dp(3)
                                                            color: root.cPanel
                                                            border.color: root.cLine
                                                            Label {
                                                                id: refLabel
                                                                anchors.centerIn: parent
                                                                text: modelData
                                                                font.family: "monospace"
                                                                font.pointSize: root.pt(root.baseFontPt - 2)
                                                                color: root.componentOf(modelData).length > 0
                                                                       ? root.cExpr : "#666"
                                                            }
                                                            MouseArea {
                                                                anchors.fill: parent
                                                                cursorShape: Qt.PointingHandCursor
                                                                onClicked: {
                                                                    var c = root.componentOf(modelData)
                                                                    if (c.length > 0) root.selectedComponent = c
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }

                                }

                                ColumnLayout {
                                    spacing: root.dp(8)

                                    RowLayout {
                                        Layout.fillWidth: true
                                        spacing: root.dp(8)
                                        // Which component is selected is shown by the highlighted
                                        // card in the left panel; repeating its name here said
                                        // nothing the panel had not already said.
                                        Label {
                                            visible: root.selectedComponent === ""
                                            text: "Select a component"
                                            font.bold: true
                                        }
                                        Label {
                                            visible: root.selectedKind.length > 0
                                            text: root.selectedKind + " ← " + root.selectedGeometryKey
                                            color: "#666"
                                            font.pointSize: root.pt(root.baseFontPt - 2)
                                        }
                                        Item { Layout.fillWidth: true }
                                    }

                                    // ── optional geometry, for any kind ──────────────
                                    //
                                    // Presence toggles rather than a name field: the optional
                                    // keys are a closed set of four, so what a component has
                                    // should be visible rather than remembered. Added in pairs
                                    // for the reason the azimuthal modes are -- an inclination
                                    // with no position angle inclines about whichever axis the
                                    // sky happens to give, and a position angle with no
                                    // inclination rotates a circle.
                                    GroupBox {
                                        Layout.fillWidth: true
                                        visible: root.selectedComponent !== ""
                                        title: "Geometry"

                                        RowLayout {
                                            spacing: root.dp(16)
                                            CheckBox {
                                                text: "position (x, y)"
                                                checked: root.hasPosition
                                                ToolTip.visible: hovered
                                                ToolTip.text: "offset from the phase centre, in mas"
                                                onToggled: root.setGeometry("position", checked)
                                            }
                                            CheckBox {
                                                text: "inclination & PA"
                                                checked: root.hasOrientation
                                                ToolTip.visible: hovered
                                                ToolTip.text: "project the component on sky: " +
                                                              "incl in degrees (0 = face-on), " +
                                                              "PA north through east"
                                                onToggled: root.setGeometry("orientation", checked)
                                            }
                                            Item { Layout.fillWidth: true }
                                            Label {
                                                text: "seeded at 0, which leaves the component " +
                                                      "exactly as it was"
                                                color: "#666"
                                                font.pointSize: root.pt(root.baseFontPt - 2)
                                            }
                                        }
                                    }

                                    // Only a Hankel component HAS a radial profile. Offering the
                                    // I(r) box for a Gaussian invites writing one, which would add
                                    // a `profile` key and silently re-identify the component as
                                    // :hankel -- the same class of quiet misparse the Structure
                                    // panel exists to catch.
                                    Label {
                                        Layout.fillWidth: true
                                        visible: root.selectedComponent !== "" &&
                                                 root.selectedGeometryKey !== "profile"
                                        wrapMode: Text.WordWrap
                                        color: "#666"
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                        text: root.selectedComponent + " is a " + root.selectedKind +
                                              ", which has no radial profile and no azimuthal " +
                                              "modes. Both belong to a component built from a " +
                                              "`profile` key."
                                    }

                                    // Each picture beside the editor it illustrates, rather than both in a column
                                    // of their own. The two columns were independent, so a preview lined up with
                                    // its section only by accident and drifted the moment the editor grew -- which
                                    // it does as soon as a profile needs a parameter created.
                                    ColumnLayout {
                                        Layout.fillWidth: true
                                        spacing: root.dp(8)
                                        visible: root.selectedGeometryKey === "profile"

                                        RowLayout {
                                            Layout.fillWidth: true
                                            Layout.alignment: Qt.AlignTop
                                            spacing: root.dp(8)
                                        // ── freeform radial profile (§3.4) ───────────────
                                        GroupBox {
                                            Layout.fillWidth: true
                                            // AlignTop, like the preview beside it: without it the
                                            // RowLayout stretches this box to the row height the
                                            // taller preview sets, and the frames stop lining up.
                                            Layout.alignment: Qt.AlignTop
                                            visible: root.selectedGeometryKey === "profile"
                                            title: "Radial profile"

                                            ColumnLayout {
                                                spacing: root.dp(6)

                                                // Templates write into the EDITOR, not the model: choosing a shape is a
                                                // draft, and Apply/Revert is already here to arbitrate it, so nothing is
                                                // replaced invisibly.
                                                RowLayout {
                                                    Layout.fillWidth: true
                                                    spacing: root.dp(6)
                                                    Label { text: "start from"; color: "#666" }
                                                    ComboBox {
                                                        id: templateBox
                                                        Layout.fillWidth: true
                                                        textRole: "label"
                                                        valueRole: "key"
                                                        model: root.profileTemplates
                                                        currentIndex: -1
                                                        displayText: currentIndex < 0 ? "choose a profile…" : currentText
                                                        onActivated: {
                                                            var t = Julia.shell_profile_template(
                                                                        root.profileTemplates[currentIndex].key)
                                                            if (t.charAt(0) === "!") { root.fitText = t; return }
                                                            var f = t.split("\t")
                                                            profileEdit.text = f[0]
                                                            root.profileSuggestedRmax = parseFloat(f[1])
                                                            root.profileSeeds = f.length > 2 ? f[2] : ""
                                                        }
                                                    }
                                                }
                                                Label {
                                                    Layout.fillWidth: true
                                                    visible: templateBox.currentIndex >= 0 && root.profileTemplates.length > 0
                                                    wrapMode: Text.WordWrap
                                                    color: "#666"
                                                    font.pointSize: root.pt(root.baseFontPt - 2)
                                                    text: (templateBox.currentIndex >= 0 && root.profileTemplates.length > 0)
                                                          ? root.profileTemplates[templateBox.currentIndex].note : ""
                                                }
                                                // The profile IS a parameter of the model -- it shows
                                                // in the table as the component's `profile` row -- but
                                                // it is the one worth editing beside its own picture,
                                                // so it is here too. Committed on focus loss rather
                                                // than per keystroke: every edit re-evaluates the
                                                // profile and its transform.
                                                TextArea {
                                                    id: profileEdit
                                                    Layout.fillWidth: true
                                                    Layout.preferredHeight: root.dp(56)
                                                    placeholderText: "I(r) in $R and $MU, e.g. (1-$s)*(1-$MU^$a)"
                                                    font.family: "monospace"
                                                    font.pointSize: root.pt(root.baseFontPt - 1)
                                                    text: root.profileText
                                                    onActiveFocusChanged: if (!activeFocus) root.commitProfile(text)
                                                }

                                                // An explicit commit as well as focus-out. Focus-out
                                                // alone is not enough: clicking dead space does not
                                                // move focus, so an edit could sit on screen looking
                                                // applied while the model still held the old profile.
                                                RowLayout {
                                                    Layout.fillWidth: true
                                                    spacing: root.dp(6)
                                                    Button {
                                                        text: "Apply"
                                                        enabled: profileEdit.text !== root.profileText
                                                        onClicked: root.commitProfile(profileEdit.text)
                                                    }
                                                    Button {
                                                        text: "Revert"
                                                        enabled: profileEdit.text !== root.profileText
                                                        onClicked: profileEdit.text = root.profileText
                                                    }
                                                    Label {
                                                        Layout.fillWidth: true
                                                        visible: profileEdit.text !== root.profileText
                                                        text: "edited, not applied"
                                                        color: root.cWarn
                                                        font.pointSize: root.pt(root.baseFontPt - 2)
                                                    }
                                                }

                                                // Writing a name in a profile CREATES the need for a
                                                // parameter rather than referring to one, so the offer
                                                // to make them belongs beside the expression.
                                                Button {
                                                    Layout.fillWidth: true
                                                    visible: root.profileMissing.length > 0
                                                    text: "Create " + root.profileMissing.join(", ")
                                                    onClicked: {
                                                        var e = Julia.shell_add_profile_params(
                                                                    root.profileKey, root.profileText,
                                                                    root.profileSeeds)
                                                        root.consoleChanged()
                                                        if (e.length > 0 && e.charAt(0) === "!") {
                                                            root.fitText = e; return
                                                        }
                                                        root.modelDirty = true
                                                        root.refreshModel()
                                                    }
                                                }

                                                // Which key set r_max is a genuine source of confusion:
                                                // udout/2 → diamout/2 → diam/2 → r_max, first one present
                                                // wins. So the grid is stated rather than assumed.
                                                RowLayout {
                                                    Layout.fillWidth: true
                                                    spacing: root.dp(6)
                                                    Label {
                                                        text: "grid: r ∈ [" + root.fmt(root.profileRMin, 3) +
                                                              ", " + root.fmt(root.profileRMax, 3) + "] mas"
                                                        font.pointSize: root.pt(root.baseFontPt - 2)
                                                        color: "#666"
                                                    }
                                                    Label {
                                                        text: "from " + root.profileRMaxKey
                                                        font.pointSize: root.pt(root.baseFontPt - 2)
                                                        color: root.cExpr
                                                    }
                                                    Item { Layout.fillWidth: true }
                                                    Label {
                                                        text: root.profileNr + " points"
                                                        font.pointSize: root.pt(root.baseFontPt - 2)
                                                        color: "#666"
                                                    }
                                                }
                                                // The grid the parser built: r_min from diamin/2, r_max
                                                // and the key that set it, and nr. Read back from
                                                // `shell_profile_grid` rather than recomputed here, so it
                                                // is the grid the model will actually be evaluated on.

                                                // The grid is cutting the profile off, rather than the profile having
                                                // died away by the edge. Actionable rather than informational: the
                                                // number on its own invites "and?", so the fix sits beside it.
                                                RowLayout {
                                                    Layout.fillWidth: true
                                                    visible: root.profileEdgeFrac > 0.01
                                                    spacing: root.dp(6)
                                                    Label {
                                                        Layout.fillWidth: true
                                                        wrapMode: Text.WordWrap
                                                        color: root.cBad
                                                        font.pointSize: root.pt(root.baseFontPt - 2)
                                                        text: "truncated: I is still " +
                                                              Math.round(root.profileEdgeFrac * 100) +
                                                              "% of its peak at r_max, so the grid cuts the profile short"
                                                    }
                                                    Button {
                                                        text: "Widen"
                                                        implicitHeight: root.dp(22)
                                                        ToolTip.visible: hovered
                                                        ToolTip.text: "grow the grid until the profile has room to fall away"
                                                        onClicked: {
                                                            // The template's own radius when there is one, else double.
                                                            var want = root.profileSuggestedRmax > root.profileRMax
                                                                       ? root.profileSuggestedRmax : root.profileRMax * 2
                                                            var r = root.rowOfKey(root.selectedComponent + ",udout")
                                                            if (r >= 0) root.commitParam(r, "value", String(2 * want))
                                                        }
                                                    }
                                                }
                                                // compile_profile treats every name that is not $R/$MU as
                                                // a parameter. Listing what it found makes a typo surface
                                                // as a spurious new parameter instead of a confusing error.
                                                Label {
                                                    Layout.fillWidth: true
                                                    wrapMode: Text.WordWrap
                                                    font.pointSize: root.pt(root.baseFontPt - 2)
                                                    color: root.profileParams.length > 0 ? "#333" : "#999"
                                                    text: root.profileParams.length > 0
                                                          ? "parameters discovered: " + root.profileParams.join(", ")
                                                          : "no parameters discovered yet"
                                                }
                                                // `compile_profile`'s discovered names, resolved
                                                // component-qualified first and then global — so a typo
                                                // shows up as a spurious new parameter rather than as a
                                                // confusing error at fit time.
                                            }
                                        }

                                            // A GroupBox, not a bare Rectangle, so its frame starts where
                                            // the editor's frame starts. A Rectangle beside a GroupBox
                                            // begins at the row's top and therefore sits a title-height
                                            // above the thing it illustrates.
                                            GroupBox {
                                                // A fixed width beside a fillWidth editor. Sharing the row by
                                                // fillWidth collapses the picture to nothing, which is what
                                                // happened when these moved out of their own column.
                                                Layout.preferredWidth: root.dp(300)
                                                Layout.alignment: Qt.AlignTop
                                                title: "I(r) and V(B)"

                                            // I(r) and V(B) share one figure and one mount: they are
                                            // read together, and two MakieAreas would mean two GL
                                            // surfaces for what is a single pair of axes.
                                            Rectangle {
                                                anchors.fill: parent
                                                implicitHeight: root.dp(240)
                                                color: "white"; border.color: root.cLine
                                                MakieArea {
                                                    id: profileArea
                                                    anchors.fill: parent
                                                    anchors.margins: root.dp(1)
                                                    visible: root.profileError.length === 0
                                                    scene: profilePlot
                                                }
                                                Label {
                                                    anchors.fill: parent
                                                    anchors.margins: root.dp(6)
                                                    visible: root.profileError.length > 0
                                                    wrapMode: Text.WordWrap
                                                    horizontalAlignment: Text.AlignHCenter
                                                    verticalAlignment: Text.AlignVCenter
                                                    color: root.cBad
                                                    font.pointSize: root.pt(root.baseFontPt - 2)
                                                    text: root.profileError
                                                }
                                                // An Observable assignment does not by itself ask Qt
                                                // for a frame; Julia has already written the curves.
                                                Connections {
                                                    target: root
                                                    function onProfileDirty() { profileArea.update() }
                                                }
                                            }
                                            }
                                        }

                                        RowLayout {
                                            Layout.fillWidth: true
                                            Layout.alignment: Qt.AlignTop
                                            spacing: root.dp(8)
                                        // ── the limb-darkening law ───────────────────────
                                    //
                                    // Beside the component rather than in the kind list: every
                                    // one of these is a disc, and which law it wears is a
                                    // question about the atmosphere, not the geometry.
                                    GroupBox {
                                        Layout.fillWidth: true
                                        visible: root.ldLawKey.length > 0
                                        title: "Limb darkening"

                                        ColumnLayout {
                                            spacing: root.dp(4)
                                            RowLayout {
                                                Layout.fillWidth: true
                                                spacing: root.dp(6)
                                                Label { text: "law"; color: "#666" }
                                                ComboBox {
                                                    Layout.fillWidth: true
                                                    textRole: "label"
                                                    valueRole: "key"
                                                    model: root.ldLaws
                                                    currentIndex: {
                                                        for (var i = 0; i < root.ldLaws.length; ++i)
                                                            if (root.ldLaws[i].key === root.ldLawKey)
                                                                return i
                                                        return -1
                                                    }
                                                    onActivated: {
                                                        var e = Julia.shell_set_ld_law(
                                                                    root.selectedComponent,
                                                                    root.ldLaws[currentIndex].key)
                                                        root.consoleChanged()
                                                        if (e.length > 0 && e.charAt(0) === "!") {
                                                            root.fitText = e; return
                                                        }
                                                        root.modelDirty = true
                                                        root.refreshModel()
                                                    }
                                                }
                                            }
                                            Label {
                                                Layout.fillWidth: true
                                                wrapMode: Text.WordWrap
                                                color: "#666"
                                                font.pointSize: root.pt(root.baseFontPt - 2)
                                                text: {
                                                    for (var i = 0; i < root.ldLaws.length; ++i)
                                                        if (root.ldLaws[i].key === root.ldLawKey)
                                                            return root.ldLaws[i].note
                                                    return ""
                                                }
                                            }
                                            // The diameter is the same quantity in every law and
                                            // is kept; the coefficients are not and are reseeded.
                                            InfoTip { tip: "Changing the law keeps the diameter and reseeds the coefficients: a value fitted under one law is not a value under another." }
                                        }
                                    }

                                    // ── azimuthal modes (§3.5) ───────────────────────
                                        GroupBox {
                                            id: azBox
                                            Layout.fillWidth: true
                                            // AlignTop, like the preview beside it: without it the
                                            // RowLayout stretches this box to the row height the
                                            // taller preview sets, and the frames stop lining up.
                                            Layout.alignment: Qt.AlignTop
                                            visible: root.selectedGeometryKey === "profile"
                                            title: "Azimuthal modes"

                                            ColumnLayout {
                                                spacing: root.dp(4)

                                                InfoTip { tip: "amp and projang are added and removed together — the parser errors if one is missing." }

                                                Repeater {
                                                    model: azModel
                                                    delegate: RowLayout {
                                                        Layout.fillWidth: true
                                                        spacing: root.dp(4)
                                                        Label {
                                                            text: "mode " + model.n
                                                            Layout.preferredWidth: root.dp(56)
                                                            font.pointSize: root.pt(root.baseFontPt - 1)
                                                        }
                                                        Label { text: "amp"; color: "#666"
                                                                font.pointSize: root.pt(root.baseFontPt - 2) }
                                                        TextField {
                                                            Layout.preferredWidth: root.dp(70)
                                                            implicitHeight: root.dp(22)
                                                            text: model.amp.toFixed(3)
                                                            font.pointSize: root.pt(root.baseFontPt - 1)
                                                            onEditingFinished:
                                                                root.commitAzField(model.n, "amp", text)
                                                        }
                                                        Label { text: "projang"; color: "#666"
                                                                font.pointSize: root.pt(root.baseFontPt - 2) }
                                                        TextField {
                                                            Layout.preferredWidth: root.dp(70)
                                                            implicitHeight: root.dp(22)
                                                            text: model.projang.toFixed(2)
                                                            font.pointSize: root.pt(root.baseFontPt - 1)
                                                            onEditingFinished:
                                                                root.commitAzField(model.n, "projang", text)
                                                        }
                                                        Item { Layout.fillWidth: true }
                                                        Button {
                                                            text: "−"
                                                            implicitWidth: root.dp(26)
                                                            implicitHeight: root.dp(22)
                                                            onClicked: root.removeAzMode(model.n)
                                                        }
                                                    }
                                                }

                                                RowLayout {
                                                    Layout.fillWidth: true
                                                    Button {
                                                        text: "+ mode"
                                                        implicitHeight: root.dp(24)
                                                        onClicked: root.addAzMode()
                                                    }
                                                    Item { Layout.fillWidth: true }
                                                }

                                                InfoTip { tip: "projang follows PMOIRED's convention, so a model moves between the two unchanged." }

                                            }
                                        }

                                            GroupBox {
                                                Layout.preferredWidth: root.dp(300)
                                                Layout.alignment: Qt.AlignTop
                                                title: "asymmetry"
                                                Rectangle {
                                                    anchors.fill: parent
                                                    implicitHeight: root.dp(110)
                                                    color: "white"; border.color: root.cLine
                                                    Label {
                                                        anchors.centerIn: parent; color: "#999"
                                                        text: "brightness preview of the asymmetry"
                                                        font.pointSize: root.pt(root.baseFontPt - 2)
                                                    }
                                                }
                                                // NOT WIRED: needs a MakieArea and a figure of its own,
                                                // built before the window like every other panel, showing
                                                // `model_to_image` of this component alone. The Model
                                                // visualization tab already renders the whole model, which
                                                // is where the asymmetry can be seen today.
                                            }
                                        }
                                    }
                                }

                                RowLayout {
                                    Layout.topMargin: root.dp(6)
                                    spacing: root.dp(5)
                                    Label {
                                        text: "Constraints"
                                        font.bold: true
                                        font.pointSize: root.pt(root.baseFontPt + 1)
                                    }
                                    InfoTip { tip: "A bound is a box on one parameter. A relation between two is not, and needs a constraint. Either side may be an expression; the right may be a number." }
                                    Item { Layout.fillWidth: true }
                                }
                                ColumnLayout {
                                    spacing: root.dp(8)

                                    Repeater {
                                        model: constraintModel
                                        delegate: RowLayout {
                                            Layout.fillWidth: true
                                            spacing: root.dp(4)
                                            TextField {
                                                Layout.fillWidth: true
                                                implicitHeight: root.dp(22)
                                                text: model.lhs
                                                font.pointSize: root.pt(root.baseFontPt - 1)
                                            }
                                            ComboBox {
                                                Layout.preferredWidth: root.dp(62)
                                                implicitHeight: root.dp(22)
                                                model: ["<", "<=", ">", ">=", "="]
                                                font.pointSize: root.pt(root.baseFontPt - 2)
                                            }
                                            TextField {
                                                Layout.fillWidth: true
                                                implicitHeight: root.dp(22)
                                                text: model.rhs
                                                font.pointSize: root.pt(root.baseFontPt - 1)
                                            }
                                            Label {
                                                text: "tol"; color: "#666"
                                                font.pointSize: root.pt(root.baseFontPt - 2)
                                            }
                                            TextField {
                                                Layout.preferredWidth: root.dp(64)
                                                implicitHeight: root.dp(22)
                                                text: String(model.tol)
                                                font.pointSize: root.pt(root.baseFontPt - 1)
                                            }
                                            // A constraint the starting model already breaks is not an
                                            // error — the fit will push it back — but it just as often
                                            // means the relation was written backwards.
                                            Label {
                                                text: model.satisfied ? "✓" : "✗"
                                                color: model.satisfied ? root.cFree : root.cBad
                                            }
                                            Button {
                                                text: "−"
                                                implicitWidth: root.dp(26)
                                                implicitHeight: root.dp(22)
                                                onClicked: constraintModel.remove(index)
                                            }
                                        }
                                    }

                                    RowLayout {
                                        Layout.fillWidth: true
                                        Button {
                                            text: "+ constraint"
                                            implicitHeight: root.dp(24)
                                            onClicked: constraintModel.append({
                                                lhs: "", op: ">", rhs: "", tol: 0.001, satisfied: true })
                                        }
                                        // The one constraint almost every multi-component model
                                        // wants, built from the components that are actually
                                        // there. The expression is one term per component, so by
                                        // three components it is long enough to mistype -- and a
                                        // normalisation that quietly leaves one out is worse than
                                        // not having one.
                                        Button {
                                            text: "Fluxes sum to 1"
                                            implicitHeight: root.dp(24)
                                            ToolTip.visible: hovered
                                            ToolTip.text: "Deriving one flux from the others is " +
                                                          "exact and needs no constraint; this is " +
                                                          "for when you want every flux free."
                                            onClicked: {
                                                var c = Julia.shell_flux_constraint()
                                                if (c.charAt(0) === "!") { root.fitText = c; return }
                                                var f = c.split("\t")
                                                for (var i = 0; i < constraintModel.count; ++i)
                                                    if (constraintModel.get(i).lhs === f[0]) return
                                                constraintModel.append({ lhs: f[0], op: f[1],
                                                                         rhs: f[2],
                                                                         tol: parseFloat(f[3]),
                                                                         satisfied: true })
                                            }
                                        }
                                        // Checks what is ON SCREEN, including a constraint
                                        // being edited and not yet used in a fit -- the panel's
                                        // own lines go to Julia rather than anything stored.
                                        Button {
                                            text: "Check"
                                            implicitHeight: root.dp(24)
                                            enabled: constraintModel.count > 0
                                            onClicked: {
                                                var r = Julia.shell_check_constraints(
                                                            root.constraintsAsLines())
                                                if (r.length > 0 && r.charAt(0) === "!") {
                                                    root.fitText = r; return
                                                }
                                                var f = r.split("\n")
                                                for (var i = 0; i < f.length &&
                                                                i < constraintModel.count; ++i)
                                                    constraintModel.setProperty(i, "satisfied",
                                                                                f[i] === "1")
                                            }
                                        }
                                        Item { Layout.fillWidth: true }
                                    }

                                    Label {
                                        Layout.fillWidth: true
                                        wrapMode: Text.WordWrap
                                        color: root.optimiser === "lsqfit" || root.optimiser === "nested"
                                               ? root.cWarn : "#666"
                                        font.pointSize: root.pt(root.baseFontPt - 2)
                                        // The difference is not a detail: a soft penalty can be
                                        // overruled by a steep enough χ², and then the fit reports a
                                        // number that is neither the constrained nor the free answer.
                                        text: (root.optimiser === "lsqfit" || root.optimiser === "nested")
                                              ? "This fitter enforces constraints with a SOFT penalty — a " +
                                                "steep χ² can overrule them. Use an NLopt optimiser for " +
                                                "constraints that must hold."
                                              : "NLopt takes these as real nonlinear constraints, so they " +
                                                "hold at the optimum. An algorithm that cannot is wrapped " +
                                                "in AUGLAG rather than replaced."
                                    }

                                    RowLayout {
                                        spacing: root.dp(5)
                                        Label { text: "Priors"; font.bold: true }
                                        InfoTip { tip: "Gaussian pulls: (expression, target, σ). The expression may be a derived quantity, not just a parameter name.\n\nPriors and constraints are supported by the Dict form of the fitters only; the FlatModel form rejects them rather than dropping them silently." }
                                        Item { Layout.fillWidth: true }
                                    }
                                    Repeater {
                                        model: priorModel
                                        delegate: RowLayout {
                                            Layout.fillWidth: true
                                            spacing: root.dp(4)
                                            TextField {
                                                Layout.fillWidth: true
                                                implicitHeight: root.dp(22)
                                                text: model.expr
                                                font.pointSize: root.pt(root.baseFontPt - 1)
                                            }
                                            TextField {
                                                Layout.preferredWidth: root.dp(70)
                                                implicitHeight: root.dp(22)
                                                text: model.target.toFixed(4)
                                                font.pointSize: root.pt(root.baseFontPt - 1)
                                            }
                                            Label { text: "±"; color: "#666" }
                                            TextField {
                                                Layout.preferredWidth: root.dp(64)
                                                implicitHeight: root.dp(22)
                                                text: model.sigma.toFixed(4)
                                                font.pointSize: root.pt(root.baseFontPt - 1)
                                            }
                                            Button {
                                                text: "−"
                                                implicitWidth: root.dp(26)
                                                implicitHeight: root.dp(22)
                                                onClicked: priorModel.remove(index)
                                            }
                                        }
                                    }
                                    Button {
                                        text: "+ prior"
                                        implicitHeight: root.dp(24)
                                        onClicked: priorModel.append({ expr: "", target: 0.0, sigma: 1.0 })
                                    }
                                }

                            }
                        }

                    ColumnLayout {
                        spacing: dp(8)

                        RowLayout {
                            Layout.fillWidth: true
                            spacing: dp(6)

                            Label { text: "nx" }
                            SpinBox {
                                from: 16; to: 2048; stepSize: 16; editable: true
                                value: root.imgNx
                                onValueModified: root.imgNx = value
                            }
                            Label { text: "pixel size" }
                            TextField {
                                Layout.preferredWidth: dp(80)
                                text: root.imgPixsize.toFixed(4)
                                horizontalAlignment: Text.AlignRight
                                onEditingFinished: {
                                    var v = parseFloat(text)
                                    if (!isNaN(v) && v > 0) root.imgPixsize = v; else text = root.imgPixsize.toFixed(4)
                                }
                            }
                            Label { text: "mas   FOV " + (root.imgNx * root.imgPixsize).toFixed(2) + " mas"
                                    color: "#666" }

                            ToolSeparator {}

                            // Live only when the model actually references them. A control that
                            // changes nothing invites the conclusion that the render is broken.
                            Label { text: "λ"; enabled: root.dependsWl }
                            TextField {
                                Layout.preferredWidth: dp(90)
                                enabled: root.dependsWl
                                text: root.imgWl.toExponential(3)
                                horizontalAlignment: Text.AlignRight
                                onEditingFinished: { var v = parseFloat(text); if (!isNaN(v)) root.imgWl = v }
                            }
                            Label { text: "m"; color: "#666"; enabled: root.dependsWl }

                            Label { text: "MJD"; enabled: root.dependsMjd }
                            TextField {
                                Layout.preferredWidth: dp(90)
                                enabled: root.dependsMjd
                                text: root.imgMjd.toFixed(3)
                                horizontalAlignment: Text.AlignRight
                                onEditingFinished: { var v = parseFloat(text); if (!isNaN(v)) root.imgMjd = v }
                            }

                            Item { Layout.fillWidth: true }
                            Button { text: "Render"; onClicked: root.renderModel() }
                        }

                        Label {
                            Layout.fillWidth: true
                            visible: !root.dependsWl && !root.dependsMjd
                            text: "this model references neither $WL nor $MJD, so it is the same at every " +
                                  "wavelength and epoch"
                            color: "#888"
                            font.pointSize: pt(baseFontPt - 1)
                            wrapMode: Text.WordWrap
                        }

                        Rectangle {
                            Layout.fillWidth: true
                            Layout.fillHeight: true
                            Layout.minimumHeight: dp(240)
                            color: "#f4f4f4"
                            border.color: "#ddd"

                            MakieArea {
                                id: modelArea
                                anchors.fill: parent
                                scene: modelPlot
                            }

                                // ── Save PNG ──────────────────────────────────
                                //
                                // Overlaid on the plot rather than placed in a toolbar, so every
                                // plot area in the application carries the same control in the
                                // same corner. Dimmed until hovered: it sits over the figure, and
                                // a button competing with the data for attention is worse than one
                                // that has to be looked for.
                                Button {
                                    anchors.right: parent.right
                                    anchors.top: parent.top
                                    anchors.margins: root.dp(6)
                                    text: "Save PNG"
                                    enabled: root.hasRender
                                    opacity: hovered ? 1.0 : 0.5
                                    ToolTip.visible: hovered
                                    ToolTip.text: "write this plot to a PNG file"
                                    onClicked: root.savePng("model", modelArea)
                                }
                            Rectangle {
                                anchors.fill: parent
                                visible: !root.hasRender
                                color: "#f4f4f4"
                                Label {
                                    anchors.centerIn: parent
                                    color: "#888"
                                    text: "press Render to see the model"
                                }
                            }
                        }

                        Label {
                            Layout.fillWidth: true
                            text: root.renderText
                            color: root.renderText.indexOf("!") === 0 ? "#c62828" : "#444"
                            elide: Text.ElideRight
                        }
                    }

                    ColumnLayout {
                        spacing: dp(8)
                        Rectangle {
                            Layout.fillWidth: true
                            implicitHeight: fitCol.implicitHeight + root.dp(12)
                            color: root.cPanel
                            border.color: root.cLine
                            radius: root.dp(3)

                            ColumnLayout {
                                id: fitCol
                                x: root.dp(8); y: root.dp(6)
                                width: parent.width - root.dp(16)
                                spacing: root.dp(6)

                                RowLayout {
                                    Layout.fillWidth: true
                                    spacing: root.dp(10)

                                    Label { text: "Observables"; font.bold: true }

                                    // Ticked-but-absent is the state that must be unrepresentable. Unticked
                                    // but present is legitimate and common — fitting V² only, or dropping a
                                    // suspect T3amp — so every enabled box stays the user's to clear.
                                    CheckBox {
                                        id: cvisBox; text: "cvis"
                                        enabled: root.haveCvis; checked: root.haveCvis
                                        ToolTip.visible: hovered && !enabled
                                        ToolTip.text: "the file has no OI_VIS"
                                    }
                                    CheckBox {
                                        id: v2Box; text: "v2"
                                        enabled: root.haveV2; checked: root.haveV2
                                        ToolTip.visible: hovered && !enabled
                                        ToolTip.text: "the file has no OI_VIS2"
                                    }
                                    CheckBox {
                                        id: t3ampBox; text: "t3amp"
                                        enabled: root.haveT3amp; checked: root.haveT3amp
                                        ToolTip.visible: hovered && !enabled
                                        ToolTip.text: "the file has no T3AMP"
                                    }
                                    CheckBox {
                                        id: t3phiBox; text: "t3phi"
                                        enabled: root.haveT3phi; checked: root.haveT3phi
                                        ToolTip.visible: hovered && !enabled
                                        ToolTip.text: "the file has no T3PHI"
                                    }
                                    CheckBox {
                                        id: diffvisBox; text: "diffvis"
                                        enabled: root.haveDiffvis; checked: false
                                        ToolTip.visible: hovered && !enabled
                                        ToolTip.text: "needs a polychromatic file with differential phase"
                                    }
                                    CheckBox {
                                        id: fluxBox; text: "flux"
                                        enabled: root.haveFlux; checked: false
                                        ToolTip.visible: hovered && !enabled
                                        ToolTip.text: "the file has no OI_FLUX"
                                    }

                                    Item { Layout.fillWidth: true }

                                }

                                // The grid's own settings, on their own row: five more controls
                                // beside the optimiser pushed the Fit button off the edge of the
                                // window, and a row that only exists for one optimiser should not
                                // cost the others their width.
                                RowLayout {
                                    Layout.fillWidth: true
                                    spacing: root.dp(8)
                                    visible: root.optimiser === "grid"

                                    // Two parameters only: the cost is exponential, and three axes
                                    // is a picture nobody can read. Every other free parameter is
                                    // held at its current value, which makes this a SLICE through
                                    // the χ² surface rather than a profile.
                                    Label { text: "Grid over"; font.bold: true }
                                    ComboBox {
                                        id: gridP1
                                        Layout.preferredWidth: root.dp(170)
                                        model: root.freeNames
                                    }
                                    Label { text: "×"; color: "#666" }
                                    ComboBox {
                                        id: gridP2
                                        Layout.preferredWidth: root.dp(170)
                                        model: root.freeNames
                                        currentIndex: Math.min(1, Math.max(0, root.freeNames.length - 1))
                                    }
                                    // One count per axis. The two parameters rarely deserve the
                                    // same resolution — a diameter may need 200 steps where its
                                    // darkening coefficient needs 40 — and tying them together
                                    // spends the budget on whichever axis did not need it.
                                    Label { text: "steps"; color: "#666" }
                                    SpinBox {
                                        id: gridN
                                        from: 4; to: 500; stepSize: 4
                                        value: 60
                                        editable: true
                                        Layout.preferredWidth: root.dp(110)
                                        ToolTip.visible: hovered
                                        ToolTip.text: "steps along " + gridP1.currentText
                                    }
                                    Label { text: "×"; color: "#666" }
                                    SpinBox {
                                        id: gridN2
                                        from: 4; to: 500; stepSize: 4
                                        value: 60
                                        editable: true
                                        Layout.preferredWidth: root.dp(110)
                                        ToolTip.visible: hovered
                                        ToolTip.text: "steps along " + gridP2.currentText
                                    }
                                    // fillWidth + elide, NOT a natural width. At its full length
                                    // this label made the row wider than the panel, which widened
                                    // the whole column and pushed the Fit button off the right
                                    // edge of the window. The step sizes move to the tooltip.
                                    Label {
                                        Layout.fillWidth: true
                                        elide: Text.ElideRight
                                        text: "= " + (gridN.value * gridN2.value) + " evaluations"
                                        color: "#888"
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                        ToolTip.visible: gridStepHover.hovered
                                        ToolTip.text: "step size: " + root.fmt(root.gridStep1, 4) +
                                                      " in " + gridP1.currentText + ", " +
                                                      root.fmt(root.gridStep2, 4) +
                                                      " in " + gridP2.currentText
                                        HoverHandler { id: gridStepHover }
                                    }
                                    Label {
                                        visible: root.nFree < 2
                                        text: "needs two free parameters"
                                        color: root.cBad
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                    }
                                    Label {
                                        visible: root.nFree >= 2 && gridP1.currentText === gridP2.currentText
                                        text: "pick two different parameters"
                                        color: root.cBad
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                    }
                                }

                                RowLayout {
                                    Layout.fillWidth: true
                                    spacing: root.dp(8)

                                    // How the point estimate is found.
                                    Label { text: "Optimiser"; font.bold: true }
                                    ComboBox {
                                        id: optBox
                                        // Wide enough for the longest entry: a dropdown that elides its
                                        // own descriptions hides the part distinguishing the choices.
                                        Layout.preferredWidth: root.dp(330)
                                        // `root.optimisers` carries the names, the keys and the
                                        // reason each unavailable entry cannot run. The three
                                        // that get used are first; the rest keep their place
                                        // below a separator rather than competing for the top
                                        // of the list.
                                        model: root.optimisers
                                        textRole: "name"
                                        valueRole: "key"
                                        onActivated: root.optimiser = root.optimiserKey(index)

                                        // The Image tab's treatment of Pigeons and OIVI, applied
                                        // to a sampler that is not loaded: dimmed, still readable,
                                        // and named alongside what would make it run. An entry
                                        // that looks like every other one and then refuses to be
                                        // chosen reads as a broken list.
                                        delegate: ItemDelegate {
                                            id: optDelegate
                                            required property int index
                                            required property var model
                                            width: optBox.width
                                            enabled: root.optimiserUsable(optDelegate.model)
                                            highlighted: optBox.highlightedIndex === optDelegate.index
                                            contentItem: RowLayout {
                                                spacing: root.dp(8)
                                                Label {
                                                    Layout.fillWidth: true
                                                    text: optDelegate.model.name
                                                    elide: Text.ElideRight
                                                    opacity: optDelegate.enabled ? 1.0 : 0.5
                                                }
                                                Label {
                                                    text: optDelegate.model.reason
                                                    visible: text.length > 0
                                                    color: "#c62828"
                                                    font.pointSize: root.pt(root.baseFontPt - 2)
                                                }
                                            }
                                        }
                                    }

                                    Label {
                                        visible: root.optimiser === "nested" && !root.allBoundsFinite
                                        text: "needs finite bounds on every free parameter"
                                        color: root.cBad
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                    }
                                    // Same treatment the Image tab gives Pigeons and OIVI: an
                                    // engine that cannot run is named along with what would
                                    // make it run, rather than silently missing.
                                    InfoTip {
                                        visible: root.nestedBackend.length === 0
                                        tip: "Nested sampling needs a sampler: start with Julia by adding NestedSamplers, or use UltraNest through PythonCall."
                                    }


                                    ToolSeparator {}

                                    // How the uncertainties are obtained. Separate, because they compose:
                                    // bootstrap is a wrapper that runs an inner fitter many times, not a
                                    // fifth optimiser.
                                    Label { text: "Uncertainty"; font.bold: true }
                                    ComboBox {
                                        id: uncBox
                                        Layout.preferredWidth: root.dp(160)
                                        // "skip", not "none": this picks whether uncertainties are
                                        // estimated at all, and "none" reads like a result the fit
                                        // came back with rather than a step being left out.
                                        model: ["skip", "Jacobian", "bootstrap", "posterior"]
                                        onActivated: root.uncertainty =
                                            ["none", "jacobian", "bootstrap", "posterior"][index]
                                    }
                                    Label {
                                        visible: root.uncertainty === "jacobian" && root.optimiser !== "lsqfit"
                                        text: "Jacobian covariance comes from Levenberg-Marquardt only"
                                        color: root.cWarn
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                    }
                                    Label {
                                        visible: root.uncertainty === "posterior" && root.optimiser !== "nested"
                                        text: "a posterior needs nested sampling"
                                        color: root.cWarn
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                    }

                                    Item { Layout.fillWidth: true }

                                    Label {
                                        visible: !root.canRun
                                        text: root.blockedReason
                                        color: root.cBad
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                    }
                                    Button {
                                        text: root.running ? "…" : "Fit"
                                        enabled: !root.running && root.canRun
                                        onClicked: root.startFit()
                                    }
                                }

                                // Bootstrap settings, shown only when bootstrap is the choice.
                                RowLayout {
                                    Layout.fillWidth: true
                                    visible: root.uncertainty === "bootstrap"
                                    spacing: root.dp(6)

                                    Label { text: "nboot"; color: "#666" }
                                    SpinBox {
                                        from: 10; to: 10000; stepSize: 10; value: root.nboot
                                        implicitHeight: root.dp(26)
                                        onValueModified: root.nboot = value
                                    }
                                    Label { text: "mode"; color: "#666" }
                                    ComboBox {
                                        Layout.preferredWidth: root.dp(140)
                                        model: ["replacement", "halfsample", "weights", "pmoired"]
                                        onActivated: root.bootMode =
                                            ["replacement", "halfsample", "weights", "pmoired"][index]
                                    }
                                    Label {
                                        visible: root.bootMode === "pmoired"
                                        text: "biased low by ≈√2"
                                        color: root.cWarn
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                    }
                                    Label { text: "granularity"; color: "#666" }
                                    ComboBox {
                                        Layout.preferredWidth: root.dp(110)
                                        // :config is the default and the right one for calibrated
                                        // interferometric data — a block is one (MJD, configuration).
                                        model: ["config", "epoch", "point"]
                                        onActivated: root.bootGranularity = ["config", "epoch", "point"][index]
                                    }
                                    Label { text: "seed"; color: "#666" }
                                    SpinBox {
                                        from: 0; to: 999999; value: root.bootSeed
                                        implicitHeight: root.dp(26)
                                        onValueModified: root.bootSeed = value
                                        // Set by default. Replicate i uses Xoshiro(seed+i), so results are
                                        // scheduling-independent — and without it the exported script does not
                                        // reproduce the numbers, which breaks the command log's guarantee.
                                    }
                                    CheckBox {
                                        text: "threaded"; checked: root.bootThreaded
                                        onToggled: root.bootThreaded = checked
                                    }
                                    Item { Layout.fillWidth: true }
                                    Label {
                                        text: "inner fitter: " + (root.optimiser === "lsqfit" ? "lsqfit" : "nlopt")
                                        color: "#666"
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                        // Defaulted to whatever is selected above, and said out loud rather
                                        // than silently using :lsqfit.
                                    }
                                }
                            }
                        }
                        Rectangle {
                            Layout.fillWidth: true
                            Layout.preferredHeight: root.dp(120)
                            color: "white"
                            border.color: root.cLine
                            radius: root.dp(3)

                            ColumnLayout {
                                anchors.fill: parent
                                anchors.margins: root.dp(1)
                                spacing: 0

                                Rectangle {
                                    Layout.fillWidth: true
                                    implicitHeight: root.dp(24)
                                    color: root.cPanel
                                    RowLayout {
                                        anchors.fill: parent
                                        anchors.leftMargin: root.dp(6)
                                        anchors.rightMargin: root.dp(6)
                                        spacing: root.dp(6)
                                        Label { text: "Fits"; font.bold: true }
                                        Item { Layout.fillWidth: true }
                                        // The model as it stands, so this works before a fit as
                                        // well as after one: whether a starting point is anywhere
                                        // near the data is worth knowing before spending an
                                        // optimiser on it.
                                        Button {
                                            text: "Residuals"; implicitHeight: root.dp(20)
                                            enabled: paramModel.count > 0
                                            ToolTip.visible: hovered
                                            ToolTip.text: "(model − data)/σ against baseline, " +
                                                          "for the values in the table"
                                            onClicked: root.showResiduals()
                                        }
                                        // Shows the map the grid search made; it does not run a
                                        // second one. Running a grid search IS the grid-search
                                        // optimiser, and a button that quietly ran another would
                                        // be a different fit under the same name.
                                        Button {
                                            text: "χ² map"; implicitHeight: root.dp(20)
                                            enabled: root.chi2MapP1.length > 0
                                            ToolTip.visible: hovered
                                            ToolTip.text: root.chi2MapP1.length > 0
                                                ? "show the grid search over " + root.chi2MapP1 +
                                                  " and " + root.chi2MapP2
                                                : "run the Grid search optimiser to produce one"
                                            onClicked: root.diagView = "chi2map"
                                        }
                                        // Not a second image surface. The Model visualization tab
                                        // already renders `model_to_image` of the current values,
                                        // and "Adopt" writes a fit's values into the table — so
                                        // this shows the fitted model by going there, rather than
                                        // by drawing the same picture somewhere else.
                                        Button {
                                            text: "Model image"; implicitHeight: root.dp(20)
                                            enabled: paramModel.count > 0
                                            ToolTip.visible: hovered
                                            ToolTip.text: "render the model as it stands, in Model visualization"
                                            onClicked: { modelTabs.currentIndex = 1; root.renderModel() }
                                        }
                                        // Only for a model that varies with wavelength. One that
                                        // does not has an SED of flat lines, which is a control
                                        // that looks broken rather than a spectrum.
                                        Button {
                                            text: "SED"; implicitHeight: root.dp(20)
                                            enabled: paramModel.count > 0 && root.dependsWl
                                            ToolTip.visible: hovered
                                            ToolTip.text: root.dependsWl
                                                ? "flux against wavelength, total and per component"
                                                : "this model references no $WL, so its SED is flat"
                                            onClicked: root.showSed()
                                        }
                                        Button {
                                            text: "Adopt"; implicitHeight: root.dp(20)
                                            enabled: fitModel.count > 0
                                            ToolTip.visible: hovered
                                            ToolTip.text: "write the fitted values into the table"
                                            onClicked: root.adoptFit()
                                        }
                                    }
                                }

                                // Column headers. Six numbers in a row are unreadable without
                                // them, and `ndof` in particular is a DATA-POINT count rather
                                // than a degrees-of-freedom, which the label has to say.
                                Rectangle {
                                    Layout.fillWidth: true
                                    implicitHeight: root.dp(18)
                                    color: "#f7f7f7"
                                    RowLayout {
                                        anchors.fill: parent
                                        spacing: root.dp(6)
                                        Label { Layout.preferredWidth: root.dp(140); text: "fit"
                                                leftPadding: root.dp(6); color: "#777"
                                                font.pointSize: root.pt(root.baseFontPt - 2) }
                                        Label { Layout.preferredWidth: root.dp(130); text: "optimiser"
                                                color: "#777"; font.pointSize: root.pt(root.baseFontPt - 2) }
                                        Label { Layout.preferredWidth: root.dp(92); text: "χ²r start"
                                                color: "#777"; font.pointSize: root.pt(root.baseFontPt - 2) }
                                        Label { Layout.preferredWidth: root.dp(92); text: "χ²r final"
                                                color: "#777"; font.pointSize: root.pt(root.baseFontPt - 2) }
                                        Label { Layout.preferredWidth: root.dp(80); text: "work"
                                                color: "#777"; font.pointSize: root.pt(root.baseFontPt - 2) }
                                        Label { Layout.preferredWidth: root.dp(90); text: "points"
                                                color: "#777"; font.pointSize: root.pt(root.baseFontPt - 2) }
                                        Label { Layout.preferredWidth: root.dp(80); text: "AIC"
                                                color: "#777"; font.pointSize: root.pt(root.baseFontPt - 2) }
                                        Label { Layout.preferredWidth: root.dp(80); text: "BIC"
                                                color: "#777"; font.pointSize: root.pt(root.baseFontPt - 2) }
                                        Label { Layout.fillWidth: true; text: "optimised parameters"
                                                color: "#777"; font.pointSize: root.pt(root.baseFontPt - 2) }
                                    }
                                }

                                ListView {
                                    id: fitList
                                    Layout.fillWidth: true
                                    Layout.fillHeight: true
                                    clip: true
                                    model: fitModel
                                    ScrollBar.vertical: ScrollBar {}

                                    delegate: RowLayout {
                                        width: fitList.width
                                        height: root.dp(22)
                                        spacing: root.dp(6)
                                        Label { Layout.preferredWidth: root.dp(140); text: model.label
                                                leftPadding: root.dp(6); elide: Text.ElideRight }
                                        Label { Layout.preferredWidth: root.dp(130); text: model.optimiser
                                                color: "#666" }
                                        // Where the fit began and where it ended, in the same
                                        // units and against the same ndof. A fit whose two
                                        // numbers are equal did not move.
                                        Label { Layout.preferredWidth: root.dp(92)
                                                text: root.fmtChi2(model.chi2rStart); color: "#666" }
                                        Label { Layout.preferredWidth: root.dp(92)
                                                text: root.fmtChi2(model.chi2r) }
                                        // Iterations or chi2 evaluations -- whichever the
                                        // optimiser reports; they are not the same quantity, so
                                        // the unit is in the cell rather than in the header.
                                        Label {
                                            Layout.preferredWidth: root.dp(80)
                                            text: model.work; color: "#666"
                                            ToolTip.visible: fitWorkHover.hovered && model.work !== "—"
                                            ToolTip.text: model.work.indexOf("it") >= 0
                                                ? "Levenberg–Marquardt iterations, each costing a residual and a Jacobian"
                                                : "χ² evaluations; NLopt reports no iteration count"
                                            HoverHandler { id: fitWorkHover }
                                        }
                                        // ndof counts data points and is NOT reduced by the number of free
                                        // parameters, so it is labelled as such and AIC/BIC are computed
                                        // properly rather than read off chi2r.
                                        Label { Layout.preferredWidth: root.dp(90);  text: model.ndof + " pts"
                                                color: "#666" }
                                        Label { Layout.preferredWidth: root.dp(80);  text: root.fmt(model.aic, 1)
                                                color: "#666" }
                                        Label { Layout.preferredWidth: root.dp(80);  text: root.fmt(model.bic, 1)
                                                color: "#666" }
                                        // Where the parameters landed. Two fits that differ only
                                        // in χ² are two numbers; this is what makes them two
                                        // models. `± σ` appears only where the fitter actually
                                        // produced one — LM's Jacobian, or a posterior.
                                        Label {
                                            Layout.fillWidth: true
                                            text: model.params
                                            font.family: "monospace"
                                            font.pointSize: root.pt(root.baseFontPt - 2)
                                            elide: Text.ElideRight
                                            ToolTip.visible: fitParamHover.hovered && text.length > 0
                                            ToolTip.text: model.params
                                            HoverHandler { id: fitParamHover }
                                        }
                                    }
                                }
                            }
                        }
                        // What the optimiser printed, verbatim. The table says where the fit
                        // landed; this says how it got there — NLopt's evaluation trace, the
                        // sampler's progress, whatever the fitter chose to report.
                        GroupBox {
                            title: "Optimiser output"
                            Layout.fillWidth: true
                            Layout.preferredHeight: root.dp(150)
                            visible: root.fitOutput.length > 0

                            OutputConsole {
                                anchors.fill: parent
                                text: root.fitOutput
                                fontPointSize: root.pt(root.baseFontPt - 2)
                            }
                        }

                        // What the optimiser had to say, beyond one row in the table. Grid
                        // search is the one that produces a picture today; the others each
                        // have one worth drawing here (a corner plot for nested sampling,
                        // residuals for the rest) and this is the space for them.
                        Rectangle {
                            Layout.fillWidth: true
                            Layout.fillHeight: true
                            Layout.minimumHeight: root.dp(220)
                            color: "#f4f4f4"
                            border.color: root.cLine
                            radius: root.dp(3)

                            MakieArea {
                                id: chi2MapArea
                                anchors.fill: parent
                                anchors.margins: root.dp(1)
                                visible: root.diagView === "chi2map"
                                scene: chi2MapPlot
                            }

                            MakieArea {
                                id: residArea
                                anchors.fill: parent
                                anchors.margins: root.dp(1)
                                visible: root.diagView === "residuals"
                                scene: residPlot
                            }

                            MakieArea {
                                id: sedArea
                                anchors.fill: parent
                                anchors.margins: root.dp(1)
                                visible: root.diagView === "sed"
                                scene: sedPlot
                            }

                                // ── Save PNG ──────────────────────────────────
                                //
                                // Overlaid on the plot rather than placed in a toolbar, so every
                                // plot area in the application carries the same control in the
                                // same corner. Dimmed until hovered: it sits over the figure, and
                                // a button competing with the data for attention is worse than one
                                // that has to be looked for.
                                Button {
                                    anchors.right: parent.right
                                    anchors.top: parent.top
                                    anchors.margins: root.dp(6)
                                    text: "Save PNG"
                                    // Saves whichever of the two is on screen. One button in
                                    // the corner every other plot area puts one in, rather than
                                    // a second that appears and disappears with the view.
                                    enabled: root.diagView === "residuals" ? root.residText.length > 0
                                           : root.diagView === "sed"       ? root.sedText.length > 0
                                           : root.chi2MapP1.length > 0
                                    opacity: hovered ? 1.0 : 0.5
                                    ToolTip.visible: hovered
                                    ToolTip.text: "write this plot to a PNG file"
                                    onClicked: root.diagView === "residuals" ? root.savePng("residuals", residArea)
                                             : root.diagView === "sed"       ? root.savePng("sed", sedArea)
                                             : root.savePng("chi2map", chi2MapArea)
                                }

                            // Covers the chart until there is one. The Makie figure exists from
                            // startup -- it has to, the GL context belongs to Qt afterwards --
                            // so without this an empty axis reads as a finished, featureless map.
                            Rectangle {
                                anchors.fill: parent
                                anchors.margins: root.dp(1)
                                visible: root.diagView === "chi2map" && root.chi2MapP1.length === 0
                                color: "#f4f4f4"
                                ColumnLayout {
                                    anchors.centerIn: parent
                                    spacing: root.dp(4)
                                    Label {
                                        Layout.alignment: Qt.AlignHCenter
                                        text: "χ² map"
                                        color: "#888"
                                    }
                                    Label {
                                        Layout.alignment: Qt.AlignHCenter
                                        text: "run a grid search to see the surface"
                                        color: "#aaa"
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                    }
                                }
                            }

                            // The contours are the quantitative part, so they are named.
                            Label {
                                visible: root.diagView === "chi2map" && root.chi2MapP1.length > 0
                                anchors.left: parent.left
                                anchors.bottom: parent.bottom
                                anchors.margins: root.dp(6)
                                text: "contours: Δχ² = 2.30, 6.17, 11.8  (68.3, 95.4, 99.7% for two parameters)"
                                color: "#666"
                                font.pointSize: root.pt(root.baseFontPt - 2)
                            }

                            Label {
                                visible: root.diagView === "sed" && root.sedText.length > 0
                                anchors.left: parent.left
                                anchors.bottom: parent.bottom
                                anchors.margins: root.dp(6)
                                text: root.sedText
                                color: "#666"
                                font.pointSize: root.pt(root.baseFontPt - 2)
                            }

                            // rms per observable, which is what the ±1 and ±3 lines are drawn
                            // at -- so the caption and the picture say the same thing.
                            Label {
                                visible: root.diagView === "residuals" && root.residText.length > 0
                                anchors.left: parent.left
                                anchors.bottom: parent.bottom
                                anchors.margins: root.dp(6)
                                text: root.residText
                                color: "#666"
                                font.pointSize: root.pt(root.baseFontPt - 2)
                            }
                        }
                    }
                }
            }
        }
    }

    // ═════════════════════════════════════════════════════════════════════════
    // state the panels above read
    // ═════════════════════════════════════════════════════════════════════════

    property string selectedComponent: ""
    onSelectedComponentChanged: { refreshProfile(); refreshLdLaw(); refreshGeometry() }

    // The geometry key that DECIDED the selected component's kind, or "" when nothing is
    // selected. `_identity_kind` searches the keys in order and the first match wins, so this
    // is the key that made the component what it is.
    readonly property string selectedGeometryKey: {
        for (var i = 0; i < componentModel.count; ++i)
            if (componentModel.get(i).name === selectedComponent)
                return componentModel.get(i).geometryKey
        return ""
    }
    readonly property string selectedKind: {
        for (var i = 0; i < componentModel.count; ++i)
            if (componentModel.get(i).name === selectedComponent)
                return componentModel.get(i).kind
        return ""
    }

    property real   profileRMin: 0.0
    property real   profileRMax: 1.0
    property string profileRMaxKey: "r_max"
    property int    profileNr: 100
    property var    profileParams: []
    property var    profileMissing: []      // names the expression needs and the model lacks

    ListModel {
        id: azModel
        // roles: n, amp, projang — always both, never one
    }

    // Send one edited field to Julia and re-read what came back. Not a local assignment: the
    // resolver may change other rows -- a derived parameter referencing this one moves with it --
    // so the table is refreshed from the model rather than patched in place.
    //
    // A field that will not parse is refused and the row redrawn with the old value, which is
    // how the box shows that nothing was accepted.
    function commitParam(row, field, text) {
        if (row < 0 || row >= paramModel.count) return
        var v = parseFloat(text)
        if (isNaN(v)) { refreshModel(); return }
        var err = Julia.shell_set_param(paramModel.get(row).key, field, v)
        if (err.length > 0 && err.charAt(0) === "!") root.fitText = err
        root.modelDirty = true
        refreshModel()
        root.consoleChanged()
    }

    // Re-read the whole model from Julia. QML keeps no model of its own -- every edit goes to
    // the shell and comes back through here, so the table cannot drift from what a fit or an
    // exported script would use.
    function refreshModel() {
        paramModel.clear()
        var rows = Julia.shell_model_rows()
        if (rows.length > 0) {
            var lines = rows.split("\n")
            for (var i = 0; i < lines.length; ++i) {
                var f = lines[i].split("\t")
                if (f.length !== 10) continue
                paramModel.append({ comp: f[0], param: f[1], key: f[2], mode: f[3],
                                    value: parseFloat(f[4]), expr: f[5],
                                    lb: parseFloat(f[6]), ub: parseFloat(f[7]),
                                    fitindex: parseInt(f[8]), atbound: f[9] === "true",
                                    kind: "" })
            }
        }

        componentModel.clear()
        var comps = Julia.shell_model_components()
        if (comps.length > 0) {
            var clines = comps.split("\n")
            for (var c = 0; c < clines.length; ++c) {
                var g = clines[c].split("\t")
                if (g.length !== 5) continue
                componentModel.append({ name: g[0], kind: g[1], geometryKey: g[2],
                                        nparams: parseInt(g[3]),
                                        nunrecognised: parseInt(g[4]) })
            }
        }

        // Three facts the model does not otherwise show: keys the parser ignored, whether the
        // resolver broadcasts every parameter over wavelength, and the globals.
        var insp = Julia.shell_model_inspection().split("\n")
        var nd = 0
        for (var d = 0; d < paramModel.count; ++d)
            if (paramModel.get(d).mode === "expr") nd++
        root.nDerived = nd

        root.unrecognisedKeys = (insp.length > 0 && insp[0].length > 0) ? insp[0].split(" ") : []
        root.broadcasting     = (insp.length > 1 && insp[1] === "true")
        root.globalKeys       = (insp.length > 2 && insp[2].length > 0) ? insp[2].split(" ") : []

        if (root.selectedComponent.length > 0) {
            var still = false
            for (var k = 0; k < componentModel.count; ++k)
                if (componentModel.get(k).name === root.selectedComponent) still = true
            if (!still) root.selectedComponent = ""
        }
        if (root.selectedComponent.length === 0 && componentModel.count > 0)
            root.selectedComponent = componentModel.get(0).name
        refreshProfile()
        refreshLdLaw()
        refreshGeometry()
        refreshAzModes()

        var chi2 = Julia.shell_model_chi2()
        root.currentChi2r = chi2.length > 0 ? parseFloat(chi2) : NaN

        var fn = Julia.shell_free_names()
        root.freeNames = fn.length > 0 ? fn.split("\n") : []

        // What is wrong with the model before anyone fits it. The flux-sum check is the one
        // that matters most: a second component added with f = 1 doubles the model's flux, and
        // the only other symptom is a χ² in the hundreds of millions.
        var w = Julia.shell_model_warnings()
        root.validationWarnings = w.length > 0 ? w.split("\n") : []

        // Residuals of the model as it was two keystrokes ago would be the wrong picture, so
        // they follow the edits -- quietly, since one transcript line per keystroke is not a
        // transcript. Only while the panel is actually showing them.
        if (root.diagView === "residuals") {
            root.residText = Julia.shell_model_residuals(true)
            residArea.update()
        }

        // Which of the render panel's two controls are live, and whether the SED is a curve
        // or a set of flat lines.
        var dep = Julia.shell_model_depends().split("\t")
        root.dependsWl  = dep.length > 0 && dep[0] === "1"
        root.dependsMjd = dep.length > 1 && dep[1] === "1"
        root.displayWl  = dep.length > 2 && dep[2].length > 0 ? parseFloat(dep[2]) : NaN
    }

    // Which diagnostic the surface under the fits table is showing. One rectangle, two
    // figures: only one is ever on screen, and a plot that appeared somewhere else would be
    // missed. "chi2map" is the default because the grid search is the optimiser that fills it
    // on its own; residuals are drawn when asked for.
    property string diagView: "chi2map"
    property string residText: ""
    property string sedText: ""

    function showSed() {
        root.sedText = Julia.shell_model_sed()
        root.consoleChanged()
        if (root.sedText.length > 0 && root.sedText.charAt(0) === "!") {
            root.fitText = root.sedText; return
        }
        root.diagView = "sed"
        sedArea.update()
    }

    function showResiduals() {
        root.residText = Julia.shell_model_residuals()
        root.consoleChanged()
        if (root.residText.length > 0 && root.residText.charAt(0) === "!") {
            root.fitText = root.residText; return
        }
        root.diagView = "residuals"
        residArea.update()
    }

    // Whether the most recent fit left a χ² map behind, and over what.
    function refreshChi2Map() {
        var info = Julia.shell_chi2_map_info()
        if (info.length === 0) { root.chi2MapP1 = ""; root.chi2MapP2 = ""; return }
        var f = info.split("\t")
        if (f.length !== 5) { root.chi2MapP1 = ""; root.chi2MapP2 = ""; return }
        root.chi2MapP1 = f[0]
        root.chi2MapP2 = f[1]
        chi2MapArea.update()
    }

    // `x[i] ↔ list_free_params[i]` holds throughout OITOOLS, so the free rows IN ORDER are the
    // parameter vector every optimiser sees and every `x_opt` reports against. Anything that
    // adds, frees or removes a row has to renumber, or the index shown beside a row stops
    // being the index the fitter uses -- which is the one thing the column is there to say.
    function renumberFree() {
        var n = 0
        for (var i = 0; i < paramModel.count; ++i) {
            var free = paramModel.get(i).mode === "free"
            paramModel.setProperty(i, "fitindex", free ? ++n : 0)
        }
    }

    // Removing a component takes its parameters with it. A row keyed "disk,fwhm" with no disk
    // component is not a model, and leaving it behind would keep fitting something the cards
    // no longer show.
    function removeComponent(name) {
        var err = Julia.shell_remove_component(name)
        if (err.length > 0) { root.fitText = err; return }
        if (root.selectedComponent === name) root.selectedComponent = ""
        root.modelDirty = true
        refreshModel()
        root.consoleChanged()
    }

    // Adding goes through Julia too, which builds the key set the parser identifies the kind
    // by. A component assembled in QML could name keys OITOOLS does not recognise, and the
    // failure would be a silently different component rather than an error.
    function addComponent(name, kind) {
        var err = Julia.shell_add_component(name, kind)
        if (err.length > 0) { root.fitText = err; return false }
        root.selectedComponent = name
        root.modelDirty = true
        refreshModel()
        root.consoleChanged()
        return true
    }

    // Renaming is Julia's job for the same reason adding is: the name is a prefix shared by
    // the component's keys, its bounds, its place in the free list and any expression that
    // refers to it, and moving some of those without the others is not a rename.
    function renameComponent(oldName, newName) {
        var err = Julia.shell_rename_component(oldName, newName)
        if (err.length > 0) { root.fitText = err; return }
        if (root.selectedComponent === oldName) root.selectedComponent = newName
        root.modelDirty = true
        refreshModel()
        root.consoleChanged()
    }

    // ── mode changes, where the cleverness lives (§3.2) ──────────────────────
    //
    // Each conversion has a consequence the user should see rather than discover:
    // expr → free needs a numeric seed, free → expr silently shortens the fit vector, and
    // fixed → free needs bounds from somewhere.
    function requestModeChange(row, newMode) {
        var r = paramModel.get(row)
        if (r.mode === newMode) return

        if (r.mode === "expr" && newMode === "free") {
            // The resolver has already evaluated the expression, so the current value is the
            // obvious seed. Offering it beats asking the user to retype a number they can see.
            seedDialog.row = row
            seedDialog.seedValue = r.value
            seedDialog.open()
            return
        }
        if (newMode === "expr") {
            // From EITHER of the other two. Going expr means supplying an expression, and
            // `shell_set_param(key, "mode", "expr")` cannot: all it does is drop the row from
            // the free list, so a fixed row asked to become derived kept its number, came back
            // from `model_rows` as fixed again, and the click did nothing at all.
            exprDialog.row = row
            exprDialog.open()
            return
        }
        root.applyMode(row, newMode)
    }

    // The expression ITSELF, which is the part `applyMode` cannot carry: `shell_set_param`'s
    // "mode" field only decides membership of the fit vector, so sending `mode = expr` and
    // nothing else left the typed string on the floor and the value still a number -- and
    // `model_rows`, reading a number, reported the row as fixed again.
    //
    // The "expr" field writes the string into the dict AND drops the row from the free list,
    // which is the whole conversion in one call.
    function commitExpr(row, text) {
        if (row < 0 || row >= paramModel.count) return
        if (text.length === 0) { refreshModel(); return }
        var err = Julia.shell_set_param(paramModel.get(row).key, "expr", text)
        if (err.length > 0 && err.charAt(0) === "!") { root.fitText = err; refreshModel(); return }
        root.modelDirty = true
        refreshModel()
        refreshProfile()
        root.consoleChanged()
    }

    // Through the shell, not a local assignment: `list_free_params` lives in Julia, and a row
    // that reads "free" while the fit vector says otherwise is the one thing the mode column
    // exists to prevent. `shell_set_param` also holds the exclusion -- free drops an expression,
    // derived drops the row from the fit vector -- so the rule stays in one place.
    function applyMode(row, newMode) {
        if (row < 0 || row >= paramModel.count) return
        var err = Julia.shell_set_param(paramModel.get(row).key, "mode", newMode)
        if (err.length > 0 && err.charAt(0) === "!") { root.fitText = err; refreshModel(); return }
        root.modelDirty = true
        refreshModel()
        root.consoleChanged()
    }

    // The optimiser list, with the reason each unavailable entry cannot run. `reason` drives
    // both the greying and the red note beside the row, exactly as the Image tab's engine list
    // does -- the two panels ask the same question and should not answer it differently.
    //
    // NLopt's names carry an LD_/LN_ prefix encoding exactly one thing --
    // `use_grad = startswith(method, "LD_")`. That is said in words in `name`; `key` is still
    // the symbol, because that is what reaches `Opt(method, n)`.
    readonly property var optimisers: [
        { key: "lsqfit",        name: "Levenberg-Marquardt  ·  lsqfit",                reason: "" },
        { key: "LN_NELDERMEAD", name: "Nelder-Mead  ·  derivative-free",               reason: "" },
        { key: "nested",
          name:   "Nested sampling  ·  " + (root.nestedBackend.length > 0
                                            ? root.nestedLabel : "no sampler loaded"),
          reason: root.nestedBackend.length > 0 ? "" : "needs Nautilus or UltraNest" },
        { key: "grid",          name: "Grid search  ·  two parameters, χ² map",        reason: "" },
        { key: "",              name: "──────────",                                    reason: "" },
        { key: "LD_LBFGS",      name: "L-BFGS  ·  gradient",                           reason: "" },
        { key: "LD_MMA",        name: "MMA  ·  gradient, takes constraints",           reason: "" },
        { key: "LD_SLSQP",      name: "SLSQP  ·  gradient, takes constraints",         reason: "" },
        { key: "LN_COBYLA",     name: "COBYLA  ·  derivative-free, takes constraints", reason: "" }
    ]

    // The separator has no key and is not an optimiser; an entry with a reason cannot run.
    function optimiserUsable(m) { return m.key.length > 0 && m.reason.length === 0 }

    // Selecting something unusable keeps whatever was chosen before rather than silently
    // switching to it -- the red note beside the row says why it was refused.
    function optimiserKey(i) {
        if (i < 0 || i >= root.optimisers.length) return optimiser
        var m = root.optimisers[i]
        return optimiserUsable(m) ? m.key : optimiser
    }

    // ═════════════════════════════════════════════════════════════════════════
    // dialogs
    // ═════════════════════════════════════════════════════════════════════════

    Dialog {
        id: seedDialog
        property int row: -1
        property real seedValue: 0
        title: "Free this parameter"
        modal: true
        anchors.centerIn: Overlay.overlay
        standardButtons: Dialog.Ok | Dialog.Cancel

        ColumnLayout {
            spacing: root.dp(6)
            Label {
                text: "A free parameter must be numeric — the resolver raises on an expression."
                wrapMode: Text.WordWrap
                Layout.maximumWidth: root.dp(360)
            }
            RowLayout {
                Label { text: "seed at" }
                TextField {
                    id: seedField
                    text: seedDialog.seedValue.toFixed(6)
                    Layout.preferredWidth: root.dp(120)
                }
                Label { text: "(current resolved value)"; color: "#666" }
            }
        }
        onAccepted: root.applyMode(row, "free")
    }

    Dialog {
        id: exprDialog
        property int row: -1
        title: "Make this an expression"
        modal: true
        anchors.centerIn: Overlay.overlay
        standardButtons: Dialog.Ok | Dialog.Cancel

        // The key being edited decides which implicit variables are in scope, so both the
        // palette and the checker are re-read whenever the dialog opens on a different row.
        property string editKey: row >= 0 && row < paramModel.count ? paramModel.get(row).key : ""
        property var keywords: []
        property var refs: []

        function refresh() {
            var kw = []
            var txt = Julia.shell_expression_keywords(editKey)
            if (txt.length > 0) {
                var lines = txt.split("\n")
                for (var i = 0; i < lines.length; ++i) {
                    var f = lines[i].split("\t")
                    if (f.length === 4)
                        kw.push({ name: f[0], meaning: f[1], enabled: f[2] === "1", reason: f[3] })
                }
            }
            keywords = kw
            check()
        }

        // Names a profile expression needs that the model does not have yet. Empty for any
        // other key: only a profile discovers its parameters from what is written.
        property var missing: []

        function check() {
            var miss = []
            var mt = Julia.shell_profile_params(editKey, exprField.text)
            if (mt.length > 0) {
                var ml = mt.split("\n")
                for (var j = 0; j < ml.length; ++j) {
                    var mf = ml[j].split("\t")
                    if (mf.length === 2) miss.push(mf[0])
                }
            }
            missing = miss

            var out = []
            var txt = Julia.shell_check_expression(editKey, exprField.text)
            if (txt.length > 0) {
                var lines = txt.split("\n")
                for (var i = 0; i < lines.length; ++i) {
                    var f = lines[i].split("\t")
                    if (f.length === 3) out.push({ ref: f[0], cls: f[1], msg: f[2] })
                }
            }
            refs = out
        }

        function insert(kw) {
            exprField.insert(exprField.cursorPosition, kw)
            exprField.forceActiveFocus()
            check()
        }

        onOpened: refresh()

        ColumnLayout {
            spacing: root.dp(6)
            Label {
                // The count changing behind the user's back is the thing to prevent.
                text: "This removes the parameter from the fit vector: " +
                      Math.max(0, root.nFree - 1) + " free parameters will remain."
                wrapMode: Text.WordWrap
                Layout.maximumWidth: root.dp(360)
                visible: !exprDialog.editKey.endsWith(",profile")
            }
            InfoTip { tip: "The radial profile I(r). Every name here that is not $R or $MU is a parameter, discovered from the expression and resolved component-qualified first, then global." }
            TextField {
                id: exprField
                Layout.fillWidth: true
                Layout.preferredWidth: root.dp(360)
                placeholderText: exprDialog.editKey.endsWith(",profile")
                                 ? "I(r), e.g. exp(-($R / $scale)^2 / 2)"
                                 : "expression, e.g. $star,ud * 3"
                font.family: "monospace"
                onTextEdited: exprDialog.check()
            }

            // ── the implicit variables, and which of them apply here ──────────
            //
            // All four always listed. The two that are out of scope are disabled with the
            // reason rather than hidden: a keyword that is simply absent reads as one that
            // does not exist, and the whole difficulty with these is that nothing anywhere
            // names them.
            Label {
                text: "Keywords"
                font.bold: true
                font.pointSize: root.pt(root.baseFontPt - 1)
            }
            Repeater {
                model: exprDialog.keywords
                delegate: RowLayout {
                    Layout.maximumWidth: root.dp(360)
                    spacing: root.dp(6)
                    Button {
                        text: modelData.name
                        enabled: modelData.enabled
                        opacity: enabled ? 1.0 : 0.5
                        font.family: "monospace"
                        implicitWidth: root.dp(52)
                        onClicked: exprDialog.insert(modelData.name)
                    }
                    Label {
                        Layout.fillWidth: true
                        wrapMode: Text.WordWrap
                        text: modelData.enabled ? modelData.meaning : modelData.reason
                        color: modelData.enabled ? "#666" : "#c62828"
                        font.pointSize: root.pt(root.baseFontPt - 2)
                    }
                }
            }

            // A profile is the one place where writing a name CREATES the need for a
            // parameter rather than referring to one, so the offer to make them belongs here
            // and not in the component card.
            Button {
                Layout.fillWidth: true
                visible: exprDialog.missing.length > 0
                text: exprDialog.missing.length === 1
                      ? "Create $" + exprDialog.missing[0]
                      : "Create " + exprDialog.missing.length + " missing parameters"
                onClicked: {
                    var err = Julia.shell_add_profile_params(exprDialog.editKey, exprField.text)
                    root.consoleChanged()
                    if (err.length > 0 && err.charAt(0) === "!") { root.fitText = err; return }
                    root.modelDirty = true
                    root.refreshModel()
                    exprDialog.check()
                }
            }

            // ── what the references resolve to, as they are typed ─────────────
            //
            // `dict_to_model` accepts an out-of-scope keyword and a misspelt parameter alike,
            // and both surface much later as the same `UndefVarError` naming an internal
            // module. Telling them apart here is the point.
            Repeater {
                model: exprDialog.refs
                delegate: Label {
                    Layout.maximumWidth: root.dp(360)
                    wrapMode: Text.WordWrap
                    visible: modelData.cls !== "ok" || modelData.msg.length > 0
                    font.pointSize: root.pt(root.baseFontPt - 2)
                    color: modelData.cls === "unknown" ? root.cBad
                         : modelData.cls === "scope"   ? root.cWarn : "#666"
                    text: "$" + modelData.ref + " — " + modelData.msg
                }
            }
        }
        onAccepted: root.commitExpr(row, exprField.text)
    }

    // Name and kind together: the kind decides which keys are written, and the parser reads the
    // kind back OFF those keys, so choosing it afterwards is not a thing that can be done.
    Dialog {
        id: addComponentDialog
        title: "Add a component"
        modal: true
        anchors.centerIn: Overlay.overlay
        standardButtons: Dialog.Ok | Dialog.Cancel
        width: root.dp(320)

        // Named after the kind -- `star`, `ring`, `gauss` -- and numbered only when that name
        // is already taken, so the first uniform disk is `star` rather than `star1`. Unique
        // against what is already there, so OK is never rejected for a reason the dialog could
        // have avoided.
        //
        // Set once the user has typed over the suggestion, so changing the kind afterwards
        // does not discard a name they chose.
        property bool nameEdited: false

        function taken(n) {
            for (var k = 0; k < componentModel.count; ++k)
                if (componentModel.get(k).name === n) return true
            return false
        }

        // The subtypes of the selected category, or [] before one is selected.
        readonly property var subs:
            (kindField.currentIndex >= 0 && kindField.currentIndex < root.componentKinds.length)
                ? root.componentKinds[kindField.currentIndex].subs : []

        // What OK sends: the subtype's id, which is what `shell_add_component` takes.
        readonly property string chosenId:
            (subs.length === 0) ? ""
          : (subs.length === 1) ? subs[0].key
          : (subField.currentIndex >= 0 && subField.currentIndex < subs.length)
                ? subs[subField.currentIndex].key : subs[0].key

        function suggestName() {
            var base = "c"
            for (var j = 0; j < subs.length; ++j)
                if (subs[j].key === chosenId) base = subs[j].base
            var name = base, i = 1
            while (taken(name)) { i++; name = base + i }
            nameField.text = name
            nameEdited = false
        }

        ColumnLayout {
            anchors.fill: parent
            spacing: root.dp(6)
            RowLayout {
                Layout.fillWidth: true
                Label { text: "name"; color: "#666"; Layout.preferredWidth: root.dp(50) }
                TextField {
                    id: nameField
                    Layout.fillWidth: true
                    selectByMouse: true
                    onTextEdited: addComponentDialog.nameEdited = true
                    onAccepted: addComponentDialog.accept()
                }
            }
            RowLayout {
                Layout.fillWidth: true
                Label { text: "kind"; color: "#666"; Layout.preferredWidth: root.dp(50) }
                ComboBox {
                    id: kindField
                    Layout.fillWidth: true
                    textRole: "label"
                    valueRole: "key"
                    model: root.componentKinds
                    // `onActivated`, not `onCurrentIndexChanged`: it fires only when the user
                    // picks, so populating the box does not overwrite a typed name.
                    onActivated: {
                        subField.currentIndex = 0
                        if (!addComponentDialog.nameEdited) addComponentDialog.suggestName()
                    }
                }
            }

            // Shown only where the shape HAS a sub-choice: a disc's law across the face, a
            // ring's cross-section. A category of one is already fully chosen.
            RowLayout {
                Layout.fillWidth: true
                visible: addComponentDialog.subs.length > 1
                Label {
                    text: kindField.currentValue === "disc" ? "law" : "shape"
                    color: "#666"; Layout.preferredWidth: root.dp(50)
                }
                ComboBox {
                    id: subField
                    Layout.fillWidth: true
                    textRole: "label"
                    valueRole: "key"
                    model: addComponentDialog.subs
                    onActivated: if (!addComponentDialog.nameEdited) addComponentDialog.suggestName()
                }
            }
            // The law's own note, for the discs, which is where the choice is least obvious.
            Label {
                Layout.fillWidth: true
                visible: text.length > 0
                wrapMode: Text.WordWrap
                color: "#666"
                font.pointSize: root.pt(root.baseFontPt - 2)
                text: {
                    for (var i = 0; i < root.ldLaws.length; ++i)
                        if (root.ldLaws[i].key === addComponentDialog.chosenId)
                            return root.ldLaws[i].note
                    if (addComponentDialog.chosenId === "ring_profile")
                        return "the ring's radial cross-section, written in $R and evaluated " +
                               "between diamin/2 and diamout/2"
                    if (addComponentDialog.chosenId === "image_func")
                        return "I(r) written in $R and $MU, Hankel-transformed to a visibility"
                    if (addComponentDialog.chosenId === "vis_func")
                        return "V written directly in $B, the baseline in Mλ"
                    return ""
                }
            }
            InfoTip { tip: "Values start somewhere sane rather than at anything measured, the geometry parameter starts free, and the flux is shared with whatever is already there." }
        }

        onAccepted: root.addComponent(nameField.text, chosenId)
    }

    // FilePicker, not QtQuick.Dialogs.FileDialog. That dialog sets `visible = false` and
    // emits `accepted` while leaving its window on screen on some systems -- measured with
    // `test/gui/filepicker_min.jl`, and see FilePicker.qml for what else was ruled out.
    // A Popup has no window of its own to leave behind.
    FilePicker {
        id: openModelDialog
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Open model file"
        filters: [{ label: "Model files (*.toml)", patterns: "*.toml" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) {
            var name = Julia.shell_open_model(path)
            root.consoleChanged()
            if (name.length === 0 || name.charAt(0) === "!") { root.fitText = name; return }
            root.modelPath = path
            root.modelName = name
            root.modelDirty = false
            root.selectedComponent = ""
            root.refreshModel()
        }
    }

    FilePicker {
        id: saveModelDialog
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Save model file"
        saveMode: true
        defaultSuffix: "toml"
        filters: [{ label: "Model files (*.toml)", patterns: "*.toml" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) {
            var err = Julia.shell_save_model(path)
            root.consoleChanged()
            if (err.length > 0) { root.fitText = err; return }
            root.modelPath = path
            root.modelDirty = false
        }
    }

    FilePicker {
        id: importPmoiredDialog
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Import a PMOIRED model"
        filters: [{ label: "Python files (*.py)", patterns: "*.py" },
                  { label: "All files", patterns: "*" }]
        // Import is a transpiler, not a validator: `shell_import_pmoired` runs the inspector
        // afterwards and puts anything the parser did not recognise in the console, which is
        // where a model that transpiled cleanly but means something else shows up.
        //
        // Replaces the model, so the free list, the bounds and the selection go with it --
        // keeping a selection that named a component of the previous model would leave the
        // panel describing something no longer there.
        onAccepted: function (path) {
            var name = Julia.shell_import_pmoired(path)
            root.consoleChanged()
            if (name.length === 0 || name.charAt(0) === "!") { root.fitText = name; return }
            root.modelPath = ""          // a .py is not a model file we can save back to
            root.modelName = name
            root.modelDirty = true       // nothing on disk holds this as a TOML yet
            root.selectedComponent = ""
            root.refreshModel()
        }
    }

    FilePicker {
        id: exportPmoiredDialog
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Export as a PMOIRED model"
        saveMode: true
        defaultSuffix: "py"
        filters: [{ label: "Python files (*.py)", patterns: "*.py" },
                  { label: "All files", patterns: "*" }]
        // `dict_to_pmoired_file` warns rather than failing on what PMOIRED cannot represent --
        // the OITOOLS-only geometries (ldlin, ldquad, ldpow, resolved). Those warnings go to
        // the console: silently writing a model that means something else on the other side is
        // the failure worth preventing.
        onAccepted: function (path) {
            var err = Julia.shell_export_pmoired(path)
            root.consoleChanged()
            if (err.length > 0 && err.charAt(0) === "!") root.fitText = err
        }
    }
}
