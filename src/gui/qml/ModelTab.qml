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
// This file is the view only. Every point where Julia will be called is marked `// wire:`.
// The data layer behind it already exists and is tested — `src/gui/model.jl` supplies
// `model_rows`, `model_inspection` and `free_parameter_vector`, whose fields are exactly the
// roles of `paramModel` below.

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
        // wire: model_rows(model_dict, list_free_params; lb, ub) → these rows, in order
    }

    // Components, as the parser identified them (§3.3 (i)). `geometryKey` is the key that
    // DECIDED `kind`: `_GEOMETRY_KEYS` is searched in order and the first match wins, so
    // naming the deciding key is worth more than naming the kind alone.
    ListModel {
        id: componentModel
        // roles: name, kind, geometryKey, nparams, nunrecognised
        // wire: model_inspection(model_dict).components
    }

    // Keys the parser ignored. Not an error, and that is the whole problem: they change what
    // the component IS, silently.
    property var unrecognisedKeys: []          // wire: model_inspection(...).unrecognised
    property bool broadcasting: false          // wire: model_inspection(...).broadcasting
    property var globalKeys: []                // wire: model_inspection(...).globals

    // Kinds the "+ component" dialog offers, as {key, label}. Read from Julia so the list is
    // the parser's own rather than a copy that can fall behind it.
    property var componentKinds: []

    // Rows in `expr` mode. The Resolution view draws a DAG of $-references, and with no derived
    // parameter there are no references and no DAG -- so the count is what decides whether that
    // panel is worth the height it takes.
    property int nDerived: 0

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

        var out = []
        var txt = Julia.shell_component_kinds()
        if (txt.length > 0) {
            var lines = txt.split("\n")
            for (var i = 0; i < lines.length; ++i) {
                var f = lines[i].split("\t")
                if (f.length === 3) out.push({ key: f[0], label: f[1], base: f[2] })
            }
        }
        root.componentKinds = out
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
        // wire: read_model_file(...).constraints; check_constraints(...) sets `satisfied`
    }

    // ── Gaussian priors ──────────────────────────────────────────────────────
    ListModel {
        id: priorModel
        // roles: expr, target, sigma
        // wire: the `priors` vector of (expr, target, sigma) tuples
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
            if (f.length === 8)
                fitModel.append({ label: f[0], optimiser: f[1], chi2r: parseFloat(f[2]),
                                  ndof: parseInt(f[3]), nevals: parseInt(f[4]),
                                  ret: "", chi2: 0.0,
                                  aic: parseFloat(f[5]), bic: parseFloat(f[6]),
                                  params: f[7] })
        }
    }

    // Write the fitted values back into the table, so what is shown is what was fitted and the
    // model can then be saved or rendered as one.
    function adoptFit() {
        var vals = Julia.shell_fit_values()
        if (vals.length === 0) return
        var lines = vals.split("\n")
        for (var i = 0; i < lines.length; ++i) {
            var j = lines[i].lastIndexOf("=")
            if (j < 0) continue
            var key = lines[i].substring(0, j)
            var v   = parseFloat(lines[i].substring(j + 1))
            for (var r = 0; r < paramModel.count; ++r)
                if (paramModel.get(r).key === key) paramModel.setProperty(r, "value", v)
        }
        modelDirty = true
    }

    // ── live feedback (§3.6) ─────────────────────────────────────────────────
    property real  currentChi2r: NaN            // wire: chi2_flat after every edit
    property int   currentNdof: 0
    property var   validationWarnings: []       // wire: display_model's own checks, rendered inline

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
        // roles: label, optimiser, chi2, chi2r, ndof, nevals, ret, aic, bic, params
        // wire: one row per completed fit, so competing models compare side by side
    }

    function fmt(v, d) { return (v === undefined || v === null || isNaN(v)) ? "—" : v.toFixed(d) }

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
                // cannot represent — the OITOOLS-only geometries, and azimuthal modes, whose
                // `az projang` differs by π/2 between the two packages. Writing a model that
                // means something else on the other side is the failure worth preventing.
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
                text: "Broadcasting is on: one chromatic expression makes every parameter a " +
                      "per-uv-point vector, not just the chromatic one."
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
                                        color: rowItem.rAtBound ? root.cWarn : "#222"
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
                                        model: ["fixed", "free", "expr"]
                                        currentIndex: rowItem.rMode === "free" ? 1
                                                    : rowItem.rMode === "expr" ? 2 : 0
                                        onActivated: root.requestModeChange(rowItem.rIndex, currentText)
                                    }

                                    // Expr rows show the expression; its value is a computed
                                    // read-only display, since editing it would mean nothing.
                                    TextField {
                                        Layout.preferredWidth: root.dp(96)
                                        Layout.maximumWidth: root.dp(96)
                                        implicitHeight: root.dp(22)
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                        readOnly: rowItem.rMode === "expr"
                                        color: rowItem.rMode === "expr" ? root.cExpr : "#222"
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
                                    // Pass the dataset and the angular-size ceiling comes from
                                    // the coverage itself — 2 λ/B_min, the largest scale the
                                    // shortest baseline senses — rather than from a constant.
                                    // wire: default_bounds(model_dict, free; data)
                                }
                                Button {
                                    text: "Link"
                                    implicitHeight: root.dp(24)
                                    // Selecting two parameters and linking creates a global and
                                    // rewrites both rows to reference it ("$PA"). Tying two
                                    // position angles together becomes one gesture instead of
                                    // an expression-editing exercise. Unlinking inlines the
                                    // current value back into each.
                                    // wire: rewrite model_dict, add the global, re-run model_rows
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

                                Label {
                                    text: "What the parser understood"
                                    font.bold: true
                                    font.pointSize: root.pt(root.baseFontPt + 1)
                                    Layout.topMargin: root.dp(6)
                                }
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
                                    // wire: display_model already performs exactly these checks —
                                    // value<lb, value>ub, lb>=ub, flux fractions not summing to 1.
                                    // Render its results rather than reimplementing them.

                                    Label { text: "Structure"; font.bold: true }
                                    Repeater {
                                        model: componentModel
                                        delegate: Rectangle {
                                            Layout.fillWidth: true
                                            implicitHeight: sCol.implicitHeight + root.dp(10)
                                            color: root.cPanel
                                            border.color: root.cLine
                                            radius: root.dp(3)
                                            ColumnLayout {
                                                id: sCol
                                                x: root.dp(6); y: root.dp(5)
                                                width: parent.width - root.dp(12)
                                                spacing: root.dp(1)
                                                Label {
                                                    text: model.name + " — " + model.kind
                                                    font.bold: true
                                                }
                                                Label {
                                                    text: "decided by: " + model.geometryKey
                                                    color: "#666"
                                                    font.pointSize: root.pt(root.baseFontPt - 2)
                                                }
                                                Label {
                                                    visible: model.nunrecognised > 0
                                                    text: model.nunrecognised + " unrecognised key(s), ignored"
                                                    color: root.cBad
                                                    font.pointSize: root.pt(root.baseFontPt - 2)
                                                }
                                            }
                                        }
                                    }

                                    // ── Resolution: what depends on what ──────────────
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
                                    Label {
                                        text: "Resolution"; font.bold: true
                                        visible: root.nDerived > 0
                                    }
                                    Label {
                                        visible: root.nDerived > 0 && root.broadcasting
                                        Layout.fillWidth: true
                                        wrapMode: Text.WordWrap
                                        color: root.cWarn
                                        font.pointSize: root.pt(root.baseFontPt - 2)
                                        text: "One expression uses $WL or $MJD, so the resolver " +
                                              "broadcasts every parameter: each is a vector over " +
                                              "the uv points, not a number."
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

                                    // ── Numeric: every parameter at its current value ─────
                                    //
                                    // `all_names` as the resolver leaves it, which is the vector
                                    // the vis functions are actually handed -- free entries in
                                    // green, derived in violet, so the three modes are the same
                                    // colours here as in the table. Clicking one selects its
                                    // component. Empty for an empty model, and hidden then.
                                    Label {
                                        text: "Numeric"; font.bold: true
                                        visible: componentModel.count > 0
                                    }
                                    Label {
                                        visible: componentModel.count > 0 && root.broadcasting
                                        Layout.fillWidth: true
                                        wrapMode: Text.WordWrap
                                        color: "#666"
                                        font.pointSize: root.pt(root.baseFontPt - 2)
                                        text: "Broadcasting: each value below is the first entry " +
                                              "of a per-λ vector."
                                    }
                                    Flow {
                                        visible: componentModel.count > 0
                                        Layout.fillWidth: true
                                        spacing: root.dp(4)
                                        Repeater {
                                            model: paramModel
                                            delegate: Rectangle {
                                                id: numChip
                                                property string keyText: model.key
                                                width: numLabel.implicitWidth + root.dp(10)
                                                height: numLabel.implicitHeight + root.dp(4)
                                                radius: root.dp(3)
                                                color: "white"
                                                border.color: root.cLine
                                                Label {
                                                    id: numLabel
                                                    anchors.centerIn: parent
                                                    text: numChip.keyText + " = " + root.fmt(model.value, 4)
                                                    font.family: "monospace"
                                                    font.pointSize: root.pt(root.baseFontPt - 2)
                                                    color: model.mode === "free" ? root.cFree
                                                         : model.mode === "expr" ? root.cExpr
                                                         : root.cFixed
                                                }
                                                MouseArea {
                                                    anchors.fill: parent
                                                    cursorShape: Qt.PointingHandCursor
                                                    onClicked: {
                                                        var c = root.componentOf(numChip.keyText)
                                                        if (c.length > 0) root.selectedComponent = c
                                                    }
                                                }
                                            }
                                        }
                                    }

                                }

                                Label {
                                    text: "Components"
                                    font.bold: true
                                    font.pointSize: root.pt(root.baseFontPt + 1)
                                    Layout.topMargin: root.dp(6)
                                }
                                ColumnLayout {
                                    spacing: root.dp(8)

                                    RowLayout {
                                        Layout.fillWidth: true
                                        spacing: root.dp(8)
                                        Label {
                                            text: root.selectedComponent === ""
                                                  ? "Select a component"
                                                  : "Component: " + root.selectedComponent
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

                                    // Editors left, pictures right. The three previews are the same
                                    // kind of thing -- what the numbers on the left look like -- so
                                    // they belong in one column where they can be compared, rather
                                    // than buried one inside each editor.
                                    RowLayout {
                                        Layout.fillWidth: true
                                        spacing: root.dp(8)
                                        visible: root.selectedGeometryKey === "profile"

                                    ColumnLayout {
                                        Layout.fillWidth: true
                                        Layout.alignment: Qt.AlignTop
                                        spacing: root.dp(8)

                                    // ── freeform radial profile (§3.4) ───────────────
                                    GroupBox {
                                        Layout.fillWidth: true
                                        visible: root.selectedGeometryKey === "profile"
                                        title: "Radial profile"

                                        ColumnLayout {
                                            spacing: root.dp(6)

                                            TextArea {
                                                Layout.fillWidth: true
                                                Layout.preferredHeight: root.dp(56)
                                                placeholderText: "I(r) in $R and $MU, e.g. (1-$s)*(1-$MU^$a)"
                                                font.family: "monospace"
                                                font.pointSize: root.pt(root.baseFontPt - 1)
                                                // wire: "cn,profile"; compile_profile on every edit
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
                                            // wire: r_min from diamin/2; r_max and the key that set it;
                                            // nr (default 100)

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
                                            // wire: compile_profile's discovered names, resolved
                                            // component-qualified first, then global
                                        }
                                    }

                                    // ── azimuthal modes (§3.5) ───────────────────────
                                    GroupBox {
                                        Layout.fillWidth: true
                                        visible: root.selectedGeometryKey === "profile"
                                        title: "Azimuthal modes"

                                        ColumnLayout {
                                            spacing: root.dp(4)

                                            Label {
                                                Layout.fillWidth: true
                                                wrapMode: Text.WordWrap
                                                color: "#666"
                                                font.pointSize: root.pt(root.baseFontPt - 2)
                                                // The parser errors unless both keys exist, so the UI
                                                // adds and removes them strictly as a pair and never
                                                // exposes one without the other.
                                                text: "amp and projang are added and removed together — " +
                                                      "the parser errors if one is missing."
                                            }

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
                                                    }
                                                    Label { text: "projang"; color: "#666"
                                                            font.pointSize: root.pt(root.baseFontPt - 2) }
                                                    TextField {
                                                        Layout.preferredWidth: root.dp(70)
                                                        implicitHeight: root.dp(22)
                                                        text: model.projang.toFixed(2)
                                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                                    }
                                                    Item { Layout.fillWidth: true }
                                                    Button {
                                                        text: "−"
                                                        implicitWidth: root.dp(26)
                                                        implicitHeight: root.dp(22)
                                                        onClicked: azModel.remove(index)   // wire: remove BOTH keys
                                                    }
                                                }
                                            }

                                            RowLayout {
                                                Layout.fillWidth: true
                                                Button {
                                                    text: "+ mode"
                                                    implicitHeight: root.dp(24)
                                                    onClicked: azModel.append({
                                                        n: azModel.count + 1, amp: 0.0, projang: 0.0 })
                                                    // wire: add "cn,az amp<N>" AND "cn,az projang<N>"
                                                }
                                                Item { Layout.fillWidth: true }
                                            }

                                            Label {
                                                Layout.fillWidth: true
                                                wrapMode: Text.WordWrap
                                                color: root.cWarn
                                                font.pointSize: root.pt(root.baseFontPt - 2)
                                                // Matters the moment someone imports a PMOIRED model.
                                                text: "Convention: OITOOLS uses +π/2 where PMOIRED uses −π/2."
                                            }

                                        }
                                    }

                                    }   // end of the left column

                                    // The pictures. Editing a profile and watching the brightness
                                    // distribution, its visibility signature and the asymmetry move
                                    // together is the whole point of doing this in a GUI.
                                    ColumnLayout {
                                        Layout.preferredWidth: root.dp(260)
                                        Layout.maximumWidth: root.dp(320)
                                        Layout.alignment: Qt.AlignTop
                                        spacing: root.dp(6)

                                        Rectangle {
                                            Layout.fillWidth: true
                                            Layout.preferredHeight: root.dp(120)
                                            color: "white"; border.color: root.cLine
                                            Label {
                                                anchors.centerIn: parent; color: "#999"
                                                text: "I(r), with μ axis"
                                                font.pointSize: root.pt(root.baseFontPt - 2)
                                            }
                                            // wire: MAKIE MOUNT. $MU = sqrt(1-(r/r_max)²), so both
                                            // axes are shown — limb-darkening profiles are written
                                            // in μ and users need both readings.
                                        }
                                        Rectangle {
                                            Layout.fillWidth: true
                                            Layout.preferredHeight: root.dp(120)
                                            color: "white"; border.color: root.cLine
                                            Label {
                                                anchors.centerIn: parent; color: "#999"
                                                text: "V(B)"
                                                font.pointSize: root.pt(root.baseFontPt - 2)
                                            }
                                            // wire: MAKIE MOUNT — the Hankel transform of what was typed
                                        }
                                        Rectangle {
                                            Layout.fillWidth: true
                                            Layout.preferredHeight: root.dp(110)
                                            color: "white"; border.color: root.cLine
                                            Label {
                                                anchors.centerIn: parent; color: "#999"
                                                text: "brightness preview of the asymmetry"
                                                font.pointSize: root.pt(root.baseFontPt - 2)
                                            }
                                            // wire: MAKIE MOUNT — model_to_image of this component
                                        }
                                    }
                                    }   // end of the two-column row
                                }

                                Label {
                                    text: "Constraints"
                                    font.bold: true
                                    font.pointSize: root.pt(root.baseFontPt + 1)
                                    Layout.topMargin: root.dp(6)
                                }
                                ColumnLayout {
                                    spacing: root.dp(8)

                                    Label {
                                        Layout.fillWidth: true
                                        wrapMode: Text.WordWrap
                                        color: "#666"
                                        font.pointSize: root.pt(root.baseFontPt - 2)
                                        // Bounds are a box, one parameter at a time. A relation between
                                        // parameters is not a box, which is the whole reason this list
                                        // exists next to the lower/upper columns.
                                        text: "A bound is a box on one parameter. A relation between two " +
                                              "is not, and needs a constraint. Either side may be an " +
                                              "expression; the right may be a number."
                                    }

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
                                        Button {
                                            text: "Check"
                                            implicitHeight: root.dp(24)
                                            // wire: check_constraints(constraints, model_dict) → satisfied
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

                                    Label { text: "Priors"; font.bold: true }
                                    Label {
                                        Layout.fillWidth: true
                                        wrapMode: Text.WordWrap
                                        color: "#666"
                                        font.pointSize: root.pt(root.baseFontPt - 2)
                                        text: "Gaussian pulls: (expression, target, σ). The expression may " +
                                              "be a derived quantity, not just a parameter name."
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
                                    Label {
                                        Layout.fillWidth: true
                                        wrapMode: Text.WordWrap
                                        color: "#666"
                                        font.pointSize: root.pt(root.baseFontPt - 3)
                                        // Both are compiled into the model, so both need the dict.
                                        text: "Priors and constraints are supported by the Dict form of the " +
                                              "fitters only; the FlatModel form rejects them rather than " +
                                              "dropping them silently."
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
                                    Label {
                                        visible: root.nestedBackend.length === 0
                                        Layout.fillWidth: true
                                        wrapMode: Text.WordWrap
                                        text: "Nested sampling needs a sampler: start with " +
                                              "`using Nautilus` (pure Julia) or " +
                                              "`using PythonCall` (UltraNest)"
                                        color: root.cWarn
                                        font.pointSize: root.pt(root.baseFontPt - 1)
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
                                        Button { text: "Residuals";  implicitHeight: root.dp(20)
                                                 enabled: fitModel.count > 0 }  // wire: plot_residuals
                                        Button { text: "χ² map";     implicitHeight: root.dp(20)
                                                 enabled: root.nFree >= 2 }     // wire: grid search over two free params
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
                                        Button { text: "SED";        implicitHeight: root.dp(20)
                                                 enabled: fitModel.count > 0 }  // wire: model_to_sed
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
                                        Label { Layout.preferredWidth: root.dp(80); text: "χ²r"
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
                                        Label { Layout.preferredWidth: root.dp(80);  text: root.fmt(model.chi2r, 4) }
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
                                scene: chi2MapPlot
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
                                    enabled: root.chi2MapP1.length > 0
                                    opacity: hovered ? 1.0 : 0.5
                                    ToolTip.visible: hovered
                                    ToolTip.text: "write this plot to a PNG file"
                                    onClicked: root.savePng("chi2map", chi2MapArea)
                                }

                            // Covers the chart until there is one. The Makie figure exists from
                            // startup -- it has to, the GL context belongs to Qt afterwards --
                            // so without this an empty axis reads as a finished, featureless map.
                            Rectangle {
                                anchors.fill: parent
                                anchors.margins: root.dp(1)
                                visible: root.chi2MapP1.length === 0
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
                                visible: root.chi2MapP1.length > 0
                                anchors.left: parent.left
                                anchors.bottom: parent.bottom
                                anchors.margins: root.dp(6)
                                text: "contours: Δχ² = 2.30, 6.17, 11.8  (68.3, 95.4, 99.7% for two parameters)"
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

        var chi2 = Julia.shell_model_chi2()
        root.currentChi2r = chi2.length > 0 ? parseFloat(chi2) : NaN

        var fn = Julia.shell_free_names()
        root.freeNames = fn.length > 0 ? fn.split("\n") : []

        // What is wrong with the model before anyone fits it. The flux-sum check is the one
        // that matters most: a second component added with f = 1 doubles the model's flux, and
        // the only other symptom is a χ² in the hundreds of millions.
        var w = Julia.shell_model_warnings()
        root.validationWarnings = w.length > 0 ? w.split("\n") : []
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
        if (r.mode === "free" && newMode === "expr") {
            exprDialog.row = row
            exprDialog.open()
            return
        }
        root.applyMode(row, newMode)
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

        ColumnLayout {
            spacing: root.dp(6)
            Label {
                // The count changing behind the user's back is the thing to prevent.
                text: "This removes the parameter from the fit vector: " +
                      Math.max(0, root.nFree - 1) + " free parameters will remain."
                wrapMode: Text.WordWrap
                Layout.maximumWidth: root.dp(360)
            }
            TextField {
                id: exprField
                Layout.fillWidth: true
                Layout.preferredWidth: root.dp(360)
                placeholderText: "expression, e.g. $star,ud * 3"
                font.family: "monospace"
            }
        }
        onAccepted: root.applyMode(row, "expr")
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

        function suggestName() {
            var base = "c"
            for (var j = 0; j < root.componentKinds.length; ++j)
                if (root.componentKinds[j].key === kindField.currentValue)
                    base = root.componentKinds[j].base
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
                    onActivated: if (!addComponentDialog.nameEdited) addComponentDialog.suggestName()
                }
            }
            Label {
                Layout.fillWidth: true
                wrapMode: Text.WordWrap
                color: "#666"
                font.pointSize: root.pt(root.baseFontPt - 2)
                text: "Values start somewhere sane rather than at anything measured, the " +
                      "geometry parameter starts free, and the flux is shared with whatever " +
                      "is already there."
            }
        }

        onAccepted: root.addComponent(nameField.text, kindField.currentValue)
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
        // wire: pmoired_to_julia_file. Import is a transpiler, not a validator, so run the
        // inspector's unrecognised-key check afterwards.
    }

    FilePicker {
        id: exportPmoiredDialog
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Export as a PMOIRED model"
        saveMode: true
        defaultSuffix: "py"
        filters: [{ label: "Python files (*.py)", patterns: "*.py" },
                  { label: "All files", patterns: "*" }]
        // wire: dict_to_pmoired_file, surfacing its warnings -- the OITOOLS-only geometries
        // (ldlin, ldquad, ldpow, resolved) and azimuthal modes cannot be represented, and
        // writing them anyway would mean something different on the other side.
    }
}
