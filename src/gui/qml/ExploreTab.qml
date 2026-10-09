// OITOOLS — the Explore perspective: the loaded data, plotted.
//
// The plot kind, its colouring, the compare views (residuals, model and image observables)
// and the per-group panels. The window owns the dataset, the status line and the PNG picker,
// and is told what happens here through the signals below; everything this tab shows is
// drawn by Julia into the `plot` scene (`src/gui/livecanvas.jl`).

import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml            // supplies `Julia`
import Makie            // MakieArea

Item {
    id: root

    // ── scaling, supplied by the window ───────────────────────────────────────
    property real uiScale:   1.0
    property real fontScale: 1.0
    property real baseFontPt: 11

    function dp(px)     { return Math.round(px * uiScale) }
    function pt(points) { return points * fontScale }

    // ── what the window tells this tab, and what it is told back ─────────────
    property bool hasDataset: false
    // The dataset's name as the window's Dataset box shows it, for the frame's Dataset section,
    // and what is in it, from `shell_dataset_summary`.
    property string datasetName: ""
    property var summaryRows: []
    function refreshSummary() {
        var t = Julia.shell_dataset_summary()
        var out = []
        if (t.length > 0) {
            var lines = t.split("\n")
            for (var i = 0; i < lines.length; ++i) {
                var f = lines[i].split("\t")
                if (f.length === 2) out.push({ label: f[0], value: f[1] })
            }
        }
        root.summaryRows = out
    }
    // Whether the frame of controls is open; the plot keeps its width when it is not.
    property bool frameOpen: true
    // The plot kind as the menu shows it, for the window's PNG file name.
    readonly property string currentKind: kindBox.currentText
    property string pickText: ""

    signal statusReported(string text)      // a status line for the window's header
    signal actionDone()                     // the window refreshes everything that follows the data
    signal savePngRequested(var area)       // the window owns the picker

    // Repaint the plot. This is NOT optional and it is easy to forget: the shell updates
    // Makie Observables rather than inserting plots, and an Observable assignment does not by
    // itself ask Qt for a new frame.
    function redraw() { makieArea.update() }

    function refreshGroupingNoun() { root.groupingNoun = Julia.shell_grouping_noun(effectiveKind()) }

    // What one panel of the per-group view would hold, for whichever observable is selected.
    property string groupingNoun: ""

    // The kind actually plotted. `t3phi`/`t3amp` each have a `_max` twin in `OBS_SPECS` that
    // differs only in which baseline it plots against, so the tick selects the twin rather
    // than adding a second axis convention anywhere.
    function effectiveKind() {
        return (maxBlBox.enabled && maxBlBox.checked) ? kindBox.currentText + "_max"
                                                      : kindBox.currentText
    }

    // Both are "", "model" or "image", and at most one of the four ticks is ever on: a
    // residual has the model subtracted already, so overplotting it there would draw the
    // prediction against its own difference.
    property string residualView: ""
    property string overlayView: ""
    // The Compare section's choice: "data", "model", "image", "residmodel" or "residimage".
    // One radio group, so "at most one" is the widget's own rule rather than code.
    property string compareMode: "data"

    // uv coverage is geometry, not an observable — there is nothing for a model to predict and
    // nothing to subtract. The ticks are disabled there rather than hidden, like every other
    // control in this bar.
    // Whether there is anything to compare against. Refreshed with the observables, so it
    // follows a fit, a reconstruction and a dataset change without a signal of its own.
    property bool haveModelObs: false
    property bool haveImageObs: false

    function refreshCompareAvailability() {
        var parts = Julia.shell_compare_available().split(",")
        for (var i = 0; i < parts.length; ++i) {
            var kv = parts[i].split("=")
            if (kv[0] === "model") root.haveModelObs = (kv[1] === "1")
            if (kv[0] === "image") root.haveImageObs = (kv[1] === "1")
        }
        // A choice that has just lost its source stops claiming to show it.
        var needModel = root.compareMode === "model" || root.compareMode === "residmodel"
        var needImage = root.compareMode === "image" || root.compareMode === "residimage"
        if ((needModel && !root.haveModelObs) || (needImage && !root.haveImageObs)) {
            root.compareMode = "data"
            root.setCompareView()
        }
    }

    function comparableKind() {
        return ["v2", "t3amp", "t3phi", "visamp", "visphi"].indexOf(kindBox.currentText) >= 0
    }

    // Why a compare box is disabled. BOTH reasons can hold at once -- a uv plot with no model
    // built -- and a tooltip naming only the first sends the user to fix the wrong one, so this
    // reports every reason it finds. `needModel` picks which of the two artefacts is wanted.
    function compareWhy(needModel) {
        var why = []
        if (!comparableKind())
            why.push("'" + kindBox.currentText + "' is not a comparable observable; " +
                     "plot v2, t3amp, t3phi, visamp or visphi")
        if (needModel  && !haveModelObs) why.push("no model defined")
        if (!needModel && !haveImageObs) why.push("no image reconstructed")
        return why.join(" — ")
    }

    function setCompareView() {
        root.residualView = root.compareMode === "residmodel" ? "model"
                          : root.compareMode === "residimage" ? "image" : ""
        root.overlayView  = root.compareMode === "model" ? "model"
                          : root.compareMode === "image" ? "image" : ""
        // Both modes are set every time: setting only the overlay when no residual is chosen
        // left the residual mode on, and the residuals stayed on screen.
        if (root.residualView.length > 0) {
            root.statusReported(Julia.shell_set_residual_mode(root.residualView))
        } else {
            Julia.shell_set_residual_mode("")
            root.statusReported(Julia.shell_set_overlay_mode(root.overlayView))
        }
        root.actionDone()
    }

    function setView() {
        root.groupingNoun = Julia.shell_grouping_noun(effectiveKind())
        root.statusReported(Julia.shell_set_view(effectiveKind(), colorBox.currentText,
                                                 logyBox.enabled && logyBox.checked,
                                                 panelsBox.enabled && panelsBox.checked))
        root.actionDone()
    }

    // Plot kinds, with whether the loaded file can supply each. Julia decides availability so
    // the menu and `canvas_data` cannot disagree; see `shell_plot_kinds`.
    ListModel {
        id: kindModel
        ListElement { name: "uv";         available: true }
        ListElement { name: "v2";         available: true }
        ListElement { name: "t3phi";      available: true }
        ListElement { name: "t3amp";      available: true }
        ListElement { name: "visamp";     available: true }
        ListElement { name: "visphi";     available: true }
        ListElement { name: "diffphi";    available: true }
        ListElement { name: "diffvisamp"; available: true }
        ListElement { name: "flux";       available: true }
    }

    function refreshPlotKinds() {
        var s = Julia.shell_plot_kinds()
        if (s.length === 0) return                    // no dataset: leave everything selectable
        var have = {}
        var parts = s.split(",")
        for (var i = 0; i < parts.length; ++i) {
            var kv = parts[i].split("=")
            if (kv.length === 2) have[kv[0]] = (kv[1] === "1")
        }
        for (var r = 0; r < kindModel.count; ++r) {
            var n = kindModel.get(r).name
            kindModel.setProperty(r, "available", have[n] === true)
        }
        // If the current view just became unavailable — a new file without that table — fall
        // back to uv rather than leaving a selection that cannot be drawn.
        if (!kindModel.get(kindBox.currentIndex).available) {
            kindBox.currentIndex = 0
            root.setView()
        }
    }

    // Poll the last identified point. Julia produces it from Makie's event stream, so there
    // is no QML signal to connect to; 150 ms is imperceptible for a label and costs nothing.
    Timer {
        interval: 150; running: true; repeat: true
        onTriggered: {
            var t = Julia.shell_pick_text()
            if (t !== root.pickText) root.pickText = t
        }
    }

    // A frame of controls on the left and the plot on the right, as in Modeling and Imaging.
    SplitView {
        anchors.fill: parent
        orientation: Qt.Horizontal

        // ── the frame ────────────────────────────────────────────────────────────
        ScrollView {
            id: frame
            visible: root.frameOpen
            SplitView.preferredWidth: root.dp(290)
            SplitView.minimumWidth: root.dp(220)
            contentWidth: availableWidth
            clip: true

            ColumnLayout {
                width: frame.availableWidth
                spacing: dp(8)

                GroupBox {
                    title: "Plot"
                    Layout.fillWidth: true
                    Layout.leftMargin: dp(6); Layout.rightMargin: dp(6); Layout.topMargin: dp(6)

                    GridLayout {
                        anchors.fill: parent
                        columns: 2
                        columnSpacing: dp(8)
                        rowSpacing: dp(4)

                        Label { text: "Plot" }
                        ComboBox {
                            id: kindBox
                            Layout.fillWidth: true
                            textRole: "name"
                            // `diffphi` is `visphi` against WAVELENGTH — OIFITS carries
                            // differential phase in OI_VIS with PHITYP=differential, so there is
                            // no separate field, only a more useful axis.
                            //
                            // Entries the loaded file cannot supply are shown but disabled,
                            // rather than hidden: which observables a file carries is worth
                            // seeing, and a list that changes length between files is hard to
                            // aim at.
                            model: kindModel
                            delegate: ItemDelegate {
                                id: kindItem
                                required property var model
                                required property int index
                                width: kindBox.width
                                enabled: model.available
                                highlighted: kindBox.highlightedIndex === index
                                // Dimmed by hand: this style draws a disabled ItemDelegate
                                // exactly like an enabled one (GUI_Explore.md D2).
                                contentItem: Label {
                                    text: kindItem.model.name
                                    opacity: kindItem.enabled ? 1.0 : 0.4
                                }
                                // The greyed rows still need a reason, or an unavailable kind
                                // looks like a broken menu rather than a file without that table.
                                ToolTip.visible: hovered && !model.available
                                ToolTip.text: "this file has no " + model.name + " data"
                            }
                            onActivated: root.setView()
                        }

                        Label { text: "Colour by" }
                        ComboBox {
                            id: colorBox
                            Layout.fillWidth: true
                            model: ["baseline", "wav", "mjd", "none"]
                            onActivated: root.setView()
                        }

                        CheckBox {
                            id: logyBox
                            Layout.columnSpan: 2
                            text: "log y"
                            // Only where it means something. A phase takes both signs, and a
                            // log axis on one drops half the points rather than rescaling them.
                            enabled: ["v2", "visamp", "t3amp", "flux", "diffvisamp"]
                                     .indexOf(kindBox.currentText) >= 0
                            onToggled: root.setView()
                        }
                        CheckBox {
                            id: panelsBox
                            Layout.columnSpan: 2
                            // Labelled from Julia, because the grouping is not always a
                            // baseline: closure quantities group by triplet and flux by station,
                            // and "per baseline" over a grid of triangles is simply wrong.
                            text: "per " + (root.groupingNoun.length > 0 ? root.groupingNoun : "group")
                            // uv coverage has no groups to split by — it is geometry, not an
                            // observable — and a single-wavelength file gives one point per panel.
                            enabled: root.groupingNoun.length > 0
                            onToggled: root.setView()
                        }
                        CheckBox {
                            id: maxBlBox
                            Layout.columnSpan: 2
                            text: "max baseline"
                            // Closure quantities have three baselines and one x position. The
                            // default is the geometric mean; this plots them against the LONGEST
                            // leg instead, which is what sets the resolution the triangle
                            // actually reaches. On the 2004 contest data the two differ by a
                            // factor 1.6, so it is not a cosmetic choice.
                            enabled: ["t3phi", "t3amp"].indexOf(kindBox.currentText) >= 0
                            ToolTip.visible: hovered
                            ToolTip.text: "closure quantities against their longest baseline, " +
                                          "not the geometric mean of the three"
                            onToggled: root.setView()
                        }
                    }
                }

                // Residuals and the model/image observables live here rather than inside
                // Modeling and Imaging: they are views OF the data, so they belong beside the
                // other views of it. Each follows the plot kind — `v2` becomes v2 residuals, on
                // the same axes and in the same groups. One choice at a time: a residual has the
                // model subtracted already, so overplotting it there would draw the prediction
                // against its own difference.
                GroupBox {
                    title: "Compare"
                    Layout.fillWidth: true
                    Layout.leftMargin: dp(6); Layout.rightMargin: dp(6)

                    ButtonGroup { id: compareGroup }
                    ColumnLayout {
                        anchors.fill: parent
                        spacing: 0
                        Repeater {
                            model: [
                                { key: "data",       label: "data only",          needModel: false },
                                { key: "model",      label: "model observables",  needModel: true  },
                                { key: "image",      label: "image observables",  needModel: false },
                                { key: "residmodel", label: "residuals (model)",  needModel: true  },
                                { key: "residimage", label: "residuals (image)",  needModel: false }
                            ]
                            RadioButton {
                                required property var modelData
                                text: modelData.label
                                ButtonGroup.group: compareGroup
                                checked: root.compareMode === modelData.key
                                enabled: modelData.key === "data" ||
                                         (root.comparableKind() &&
                                          (modelData.needModel ? root.haveModelObs : root.haveImageObs))
                                ToolTip.visible: hovered && !enabled
                                ToolTip.text: root.compareWhy(modelData.needModel)
                                onClicked: { root.compareMode = modelData.key; root.setCompareView() }
                            }
                        }
                    }
                }

                // What is plotted: read-only, from the file and the loaded data. Below the
                // controls, so its length -- which changes with the dataset -- never moves them.
                GroupBox {
                    title: "Dataset"
                    Layout.fillWidth: true
                    Layout.leftMargin: dp(6); Layout.rightMargin: dp(6)

                    GridLayout {
                        anchors.fill: parent
                        columns: 2
                        columnSpacing: dp(8)
                        rowSpacing: dp(2)

                        Label {
                            Layout.columnSpan: 2
                            Layout.fillWidth: true
                            text: root.datasetName.length > 0 ? root.datasetName : "no dataset loaded"
                            color: root.datasetName.length > 0 ? "#000" : "#888"
                            font.bold: root.datasetName.length > 0
                            elide: Text.ElideMiddle
                        }
                        Repeater {
                            model: root.summaryRows
                            delegate: Label {
                                required property var modelData
                                required property int index
                                Layout.row: 1 + index
                                Layout.column: 0
                                Layout.alignment: Qt.AlignTop
                                text: modelData.label
                                color: "#666"
                                font.pointSize: pt(baseFontPt - 1)
                            }
                        }
                        Repeater {
                            model: root.summaryRows
                            delegate: Label {
                                required property var modelData
                                required property int index
                                Layout.row: 1 + index
                                Layout.column: 1
                                Layout.fillWidth: true
                                Layout.preferredWidth: 1
                                text: modelData.value
                                wrapMode: Text.Wrap
                                font.pointSize: pt(baseFontPt - 1)
                            }
                        }
                    }
                }

                Item { Layout.fillHeight: true }
            }
        }

        // ── the plot ─────────────────────────────────────────────────────────────
        ColumnLayout {
            SplitView.fillWidth: true
            spacing: 0

            Item {
                Layout.fillWidth: true
                Layout.fillHeight: true

                MakieArea {
                    id: makieArea
                    anchors.fill: parent
                    scene: plot

                    // Picking uses a TapHandler, NOT a MouseArea laid over the scene.
                    //
                    // MakieArea's own input handling lives in a MouseArea whose `visible` is
                    // bound to `parent.activeFocus`, and focus is granted by a TapHandler in
                    // MakieArea that fires on the first press. A MouseArea of ours stacked on
                    // top consumed that press, so MakieArea never took focus, its MouseArea
                    // never became visible, and Makie received nothing at all: no zoom, no
                    // pan, no picking. `propagateComposedEvents` does not help -- it forwards
                    // only composed events (clicked/doubleClicked/pressAndHold), never the
                    // press, the drag, or the wheel that zoom and pan are built from.
                    //
                    // Pointer handlers cooperate instead of stealing: this one coexists with
                    // MakieArea's TapHandler, and its default DragThreshold gesture policy
                    // bows out as soon as the pointer moves, so a drag still pans.
                    // No pointer handlers here on purpose. MakieArea fills itself with a
                    // MouseArea accepting all buttons, and a MouseArea takes an exclusive
                    // grab on press, which cancels any TapHandler of ours the moment the
                    // plot has focus. Clicks are read from Makie's event stream instead --
                    // see install_interactions! in src/shell.jl -- and the result is polled
                    // below.
                }

                // Over the plot, as the other perspectives place theirs: the frame's toggle in
                // one corner, Save PNG in the other. Dimmed until hovered, so neither competes
                // with the data for attention.
                Button {
                    anchors.left: parent.left
                    anchors.top: parent.top
                    anchors.margins: dp(6)
                    text: root.frameOpen ? "◀" : "▶"
                    opacity: hovered ? 1.0 : 0.5
                    ToolTip.visible: hovered
                    ToolTip.text: root.frameOpen ? "hide the controls" : "show the controls"
                    onClicked: root.frameOpen = !root.frameOpen
                }
                Button {
                    anchors.right: parent.right
                    anchors.top: parent.top
                    anchors.margins: dp(6)
                    text: "Save PNG"
                    enabled: root.hasDataset
                    opacity: hovered ? 1.0 : 0.5
                    ToolTip.visible: hovered
                    ToolTip.text: "write this plot to a PNG file"
                    onClicked: root.savePngRequested(makieArea)
                }
            }

            Label {
                Layout.fillWidth: true
                Layout.margins: dp(6)
                text: root.pickText.length > 0 ? root.pickText : "click a point to identify it"
                color: root.pickText.length > 0 ? "#000" : "#888"
                font.family: "monospace"
                elide: Text.ElideRight
            }
        }
    }
}
