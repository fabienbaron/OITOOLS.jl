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

    // Free-parameter names, for the grid's two axis pickers. Refreshed with the model, because
    // freeing or fixing a parameter changes what a grid can be run over.
    property var freeNames: []

    // What the last grid fit mapped, or empty when the last fit was not one. The map panel
    // shows the chart only when there is a chart to show.
    property string chi2MapP1: ""
    property string chi2MapP2: ""

    Component.onCompleted: {
        var out = []
        var txt = Julia.shell_component_kinds()
        if (txt.length > 0) {
            var lines = txt.split("\n")
            for (var i = 0; i < lines.length; ++i) {
                var f = lines[i].split("\t")
                if (f.length === 2) out.push({ key: f[0], label: f[1] })
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
                gridP1.currentText, gridP2.currentText, gridN.value)
            root.running = false
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
            if (f.length === 7)
                fitModel.append({ label: f[0], optimiser: f[1], chi2r: parseFloat(f[2]),
                                  ndof: parseInt(f[3]), nevals: parseInt(f[4]),
                                  ret: "", chi2: 0.0,
                                  aic: parseFloat(f[5]), bic: parseFloat(f[6]) })
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
    property bool vonmises: false

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

    // UltraNest errors unless EVERY free parameter has finite bounds, so it is disabled with
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
        // roles: label, optimiser, chi2, chi2r, ndof, nevals, ret, aic, bic
        // wire: one row per completed fit, so competing models compare side by side
    }

    function fmt(v, d) { return (v === undefined || v === null || isNaN(v)) ? "—" : v.toFixed(d) }

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
                onClicked: openModelDialog.openAt("")
            }
            Button {
                text: "Save…"
                enabled: paramModel.count > 0
                onClicked: saveModelDialog.openAt("")
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
                                    Label { text: model.name; font.bold: true }
                                    Label {
                                        text: model.kind + " ← " + model.geometryKey
                                        color: "#666"
                                        font.pointSize: root.pt(root.baseFontPt - 2)
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

                    Button {
                        text: "+ component"
                        onClicked: {
                            addComponentDialog.suggestName()
                            addComponentDialog.open()
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
                                        implicitHeight: root.dp(22)
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                        readOnly: rowItem.rMode === "expr"
                                        color: rowItem.rMode === "expr" ? root.cExpr : "#222"
                                        text: isNaN(rowItem.rValue) ? "—" : rowItem.rValue.toFixed(4)
                                        horizontalAlignment: Text.AlignRight
                                        // wire: commit the edited value into model_dict, then
                                        // re-resolve and refresh chi2
                                        onEditingFinished: root.modelDirty = true
                                    }

                                    TextField {
                                        Layout.preferredWidth: root.dp(76)
                                        implicitHeight: root.dp(22)
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                        enabled: rowItem.rMode === "free"
                                        opacity: enabled ? 1 : 0.35
                                        text: rowItem.rMode === "free" ? root.fmt(rowItem.rLb, 3) : ""
                                        horizontalAlignment: Text.AlignRight
                                        onEditingFinished: root.modelDirty = true   // wire: lb
                                    }
                                    TextField {
                                        Layout.preferredWidth: root.dp(76)
                                        implicitHeight: root.dp(22)
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                        enabled: rowItem.rMode === "free"
                                        opacity: enabled ? 1 : 0.35
                                        text: rowItem.rMode === "free" ? root.fmt(rowItem.rUb, 3) : ""
                                        horizontalAlignment: Text.AlignRight
                                        onEditingFinished: root.modelDirty = true   // wire: ub
                                    }

                                    // Where the value sits inside its box. `default_bounds` is
                                    // deliberately generous, so a span-relative test would call
                                    // every small diameter pinned; what is worth showing is the
                                    // value landing ON a bound, which is what `atbound` means.
                                    Rectangle {
                                        Layout.fillWidth: true
                                        implicitHeight: root.dp(6)
                                        visible: rowItem.rMode === "free"
                                        color: "#e9e9e9"
                                        radius: root.dp(3)

                                        Rectangle {
                                            visible: isFinite(rowItem.rLb) && isFinite(rowItem.rUb)
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
                                            visible: !(isFinite(rowItem.rLb) && isFinite(rowItem.rUb))
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

                                    Label { text: "Resolution"; font.bold: true }
                                    Rectangle {
                                        Layout.fillWidth: true
                                        implicitHeight: root.dp(140)
                                        color: "white"
                                        border.color: root.cLine
                                        radius: root.dp(3)
                                        Label {
                                            anchors.centerIn: parent
                                            color: "#999"
                                            horizontalAlignment: Text.AlignHCenter
                                            // The resolver's own order: globals → fixed → derived,
                                            // topologically sorted. Drawn as a DAG of $-references so
                                            // dependencies are visible rather than inferred.
                                            text: "dependency DAG of $-references\n(globals → fixed → derived)"
                                        }
                                        // wire: the resolver's topological order → a small DAG
                                    }

                                    Label { text: "Numeric"; font.bold: true }
                                    Rectangle {
                                        Layout.fillWidth: true
                                        implicitHeight: root.dp(120)
                                        color: "white"
                                        border.color: root.cLine
                                        radius: root.dp(3)
                                        Label {
                                            anchors.centerIn: parent
                                            color: "#999"
                                            horizontalAlignment: Text.AlignHCenter
                                            // Chromatic parameters are vectors, so they render as a
                                            // curve versus λ rather than as a number.
                                            text: root.broadcasting
                                                  ? "all_names at current values — per-λ curves"
                                                  : "all_names at current values"
                                        }
                                        // wire: the evaluated all_names vector, free entries highlighted
                                    }

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
                                }

                                Label {
                                    text: "Components"
                                    font.bold: true
                                    font.pointSize: root.pt(root.baseFontPt + 1)
                                    Layout.topMargin: root.dp(6)
                                }
                                ColumnLayout {
                                    spacing: root.dp(8)

                                    Label {
                                        text: root.selectedComponent === ""
                                              ? "Select a component"
                                              : "Component: " + root.selectedComponent
                                        font.bold: true
                                    }

                                    // ── freeform radial profile (§3.4) ───────────────
                                    GroupBox {
                                        Layout.fillWidth: true
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

                                            RowLayout {
                                                Layout.fillWidth: true
                                                spacing: root.dp(6)
                                                // I(r) and V(B) side by side: editing the profile and
                                                // seeing both the brightness distribution and its
                                                // visibility signature update together is the whole
                                                // point of a GUI here.
                                                Rectangle {
                                                    Layout.fillWidth: true
                                                    Layout.preferredHeight: root.dp(120)
                                                    color: "white"; border.color: root.cLine
                                                    Label {
                                                        anchors.centerIn: parent; color: "#999"
                                                        text: "I(r), with μ axis"
                                                        font.pointSize: root.pt(root.baseFontPt - 2)
                                                    }
                                                    // wire: MAKIE MOUNT. $MU = sqrt(1-(r/r_max)²), so
                                                    // both axes are shown — limb-darkening profiles are
                                                    // written in μ and users need both readings.
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
                                                    // wire: MAKIE MOUNT — the Hankel transform of what
                                                    // was typed
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
                                            // wire: compile_profile's discovered names, resolved
                                            // component-qualified first, then global
                                        }
                                    }

                                    // ── azimuthal modes (§3.5) ───────────────────────
                                    GroupBox {
                                        Layout.fillWidth: true
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
                                    }
                                }

                                Label {
                                    text: "Constraints"
                                    font.bold: true
                                    font.pointSize: root.pt(root.baseFontPt + 1)
                                    Layout.topMargin: root.dp(6)
                                }
                                ColumnLayout {
                                    spacing: root.dp(8)

                                    Label { text: "Constraints"; font.bold: true }
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
                                        color: root.optimiser === "lsqfit" || root.optimiser === "ultranest"
                                               ? root.cWarn : "#666"
                                        font.pointSize: root.pt(root.baseFontPt - 2)
                                        // The difference is not a detail: a soft penalty can be
                                        // overruled by a steep enough χ², and then the fit reports a
                                        // number that is neither the constrained nor the free answer.
                                        text: (root.optimiser === "lsqfit" || root.optimiser === "ultranest")
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

                                    CheckBox {
                                        text: "von Mises"
                                        checked: root.vonmises
                                        onToggled: root.vonmises = checked
                                        ToolTip.visible: hovered
                                        ToolTip.text: "circular statistic for T3phi instead of a wrapped difference"
                                    }
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
                                    Label { text: "points/axis"; color: "#666" }
                                    SpinBox {
                                        id: gridN
                                        from: 8; to: 400; stepSize: 4
                                        value: 60
                                        editable: true
                                        Layout.preferredWidth: root.dp(120)
                                        ToolTip.visible: hovered
                                        ToolTip.text: "the grid costs this squared: 60 is 3600 χ² evaluations"
                                    }
                                    Label {
                                        text: gridN.value + "² = " + (gridN.value * gridN.value) + " evaluations"
                                        color: "#888"
                                        font.pointSize: root.pt(root.baseFontPt - 1)
                                    }
                                    Item { Layout.fillWidth: true }
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
                                        // NLopt's names carry an LD_/LN_ prefix encoding exactly one
                                        // thing — `use_grad = startswith(method, "LD_")`. That is said
                                        // in words here instead, so the list reads as algorithms rather
                                        // than as symbols. The SYMBOL is unchanged: `optimiserKey` still
                                        // returns :LD_LBFGS, because that is what reaches Opt(method, n).
                                        // The three that get used are first; the rest keep
                                        // their place below a separator rather than competing
                                        // for the top of the list.
                                        model: ["Levenberg-Marquardt  ·  lsqfit",
                                                "Nelder-Mead  ·  derivative-free",
                                                "Nested sampling  ·  UltraNest",
                                                "Grid search  ·  two parameters, χ² map",
                                                "──────────",
                                                "L-BFGS  ·  gradient",
                                                "MMA  ·  gradient, takes constraints",
                                                "SLSQP  ·  gradient, takes constraints",
                                                "COBYLA  ·  derivative-free, takes constraints"]
                                        onActivated: root.optimiser = root.optimiserKey(index)
                                    }

                                    Label {
                                        visible: root.optimiser === "ultranest" && !root.allBoundsFinite
                                        text: "needs finite bounds on every free parameter"
                                        color: root.cBad
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
                                        model: ["none", "Jacobian", "bootstrap", "posterior"]
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
                                        visible: root.uncertainty === "posterior" && root.optimiser !== "ultranest"
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
                                        Button { text: "Model image"; implicitHeight: root.dp(20)
                                                 enabled: fitModel.count > 0 }  // wire: model_to_image
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
                                        Label { Layout.fillWidth: true; text: model.ret; color: "#666" }
                                    }
                                }
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

    property real   profileRMin: 0.0
    property real   profileRMax: 1.0
    property string profileRMaxKey: "r_max"
    property int    profileNr: 100
    property var    profileParams: []

    ListModel {
        id: azModel
        // roles: n, amp, projang — always both, never one
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

    function applyMode(row, newMode) {
        paramModel.setProperty(row, "mode", newMode)
        root.modelDirty = true
        // wire: rebuild list_free_params from the rows in order, refresh fitindex on every
        // row, pull default bounds for a newly free one, and re-run model_rows
    }

    // Index 3 is the separator, which is not an optimiser: selecting it keeps whatever was
    // chosen before rather than silently switching to something else.
    function optimiserKey(i) {
        var k = ["lsqfit", "LN_NELDERMEAD", "ultranest", "grid", "",
                 "LD_LBFGS", "LD_MMA", "LD_SLSQP", "LN_COBYLA"][i]
        return (k === undefined || k.length === 0) ? optimiser : k
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

        // "c1", "c2", ... — unique against what is already there, so OK is never rejected for
        // a reason the dialog could have avoided.
        function suggestName() {
            var i = 1, taken = true
            while (taken) {
                taken = false
                for (var k = 0; k < componentModel.count; ++k)
                    if (componentModel.get(k).name === "c" + i) { taken = true; i++; break }
            }
            nameField.text = "c" + i
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
                }
            }
            Label {
                Layout.fillWidth: true
                wrapMode: Text.WordWrap
                color: "#666"
                font.pointSize: root.pt(root.baseFontPt - 2)
                text: "Values start somewhere sane rather than at anything measured; set them " +
                      "in the table, then free what should be fitted."
            }
        }

        onAccepted: root.addComponent(nameField.text, kindField.currentValue)
    }

    // FilePicker, not QtQuick.Dialogs.FileDialog. That dialog sets `visible = false` and
    // emits `accepted` while leaving its window on screen on some systems -- measured with
    // `test/gui/filepicker_min.jl`. A Popup has no window of its own to leave behind.
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
