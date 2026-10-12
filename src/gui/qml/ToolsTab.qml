// OITOOLS — the Tools perspective: operations on files and datasets that are none of the other
// four.
//
// Files: the header inspector and the validator work on any FITS file, chosen here or taken
// from the current dataset, and the validator's "repaired copy" writes a new file. Current
// dataset: channel averaging and the splits, each a re-read of the dataset's file that becomes
// a new dataset (so Explore, Modeling and Imaging use it at once), and Export. Everything is
// done by `src/gui/tools.jl`; this file only lays it out.

import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml            // supplies `Julia`
import Makie            // MakieArea, for the cube viewer's canvas (`cubePlot`, built in gui())

Item {
    id: root

    // ── scaling, supplied by the window ───────────────────────────────────────
    property real uiScale:   1.0
    property real fontScale: 1.0
    property real baseFontPt: 11
    function dp(px)     { return Math.round(px * uiScale) }
    function pt(points) { return points * fontScale }

    signal consoleChanged()                 // Julia wrote to the console
    signal datasetsChanged()                // new datasets were made; the window lists them
    signal statusReported(string text)

    // What the result pane shows: "", "headers", "validate", "split" or "cube".
    property string mode: ""
    property string headersPath: ""

    function failed(t) { return t.length > 0 && t[0] === "!" }
    function previewSplit(by, param) {
        var t = Julia.shell_tools_split_preview(filePath.text, by, param)
        root.mode = "split"
        resultArea.text = failed(t) ? t : t.split("\n").length + " parts:\n" + t.split("\t").join("    ")
    }
    function runSplit(by, param, folder) {
        report(Julia.shell_tools_split(filePath.text, by, param, folder))
        root.datasetsChanged()
    }
    function report(t) { root.statusReported(t); root.consoleChanged() }

    function refreshTargets() {
        var t = filePath.text.length > 0 ? Julia.shell_tools_targets(filePath.text) : ""
        var names = (t.length > 0 && !failed(t)) ? t.split("\n") : []
        targetBox.model = names
        targetBox.visible = names.length > 1
    }

    function useCurrentDataset(field) {
        var p = Julia.shell_tools_source()
        if (p.length === 0) { root.statusReported("! no single-file dataset is selected"); return }
        field.text = p
    }

    function showHeaders() {
        if (fitsPath.text.length === 0) return
        var t = Julia.shell_tools_headers(fitsPath.text)
        root.consoleChanged()
        if (failed(t)) { root.statusReported(t); return }
        var rows = t.split("\n")
        hduModel.clear()
        for (var i = 0; i < rows.length; ++i) {
            var f = rows[i].split("\t")
            if (f.length < 7) continue
            hduModel.append({ idx: parseInt(f[0]), kind: f[1], extname: f[2] || "(primary)",
                              insname: f[3], revn: f[4], nrows: f[5], extra: f[6] })
        }
        root.headersPath = fitsPath.text
        root.mode = "headers"
        hduList.currentIndex = 0
        showHdu(1)
    }

    function showHdu(i) { resultArea.text = Julia.shell_tools_hdu(root.headersPath, i) }

    function validate() {
        if (filePath.text.length === 0) return
        var t = Julia.shell_tools_validate(filePath.text)
        root.consoleChanged()
        root.mode = "validate"
        if (failed(t)) { resultArea.text = t; return }
        if (t.length === 0) { resultArea.text = "Nothing found."; return }
        var rows = t.split("\n"), out = []
        var tag = { error: "ERROR  ", warning: "warning", info: "note   " }
        for (var i = 0; i < rows.length; ++i) {
            var f = rows[i].split("\t")
            out.push((tag[f[0]] || f[0]) + "  " + (f[1].length > 0 ? f[1] + ": " : "") + f[2])
        }
        resultArea.text = out.join("\n")
    }

    // ── pickers ───────────────────────────────────────────────────────────────
    // One picker per file field: OIFITS handling takes OIFITS files, FITS handling FITS files.
    FilePicker {
        id: openPicker
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Choose an OIFITS file"
        filters: [{ label: "OIFITS files (*.oifits *.fits *.fit)", patterns: "*.oifits,*.fits,*.fit" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) { filePath.text = path }
    }
    FilePicker {
        id: fitsPicker
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Choose a FITS file"
        filters: [{ label: "FITS files (*.fits)", patterns: "*.fits" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) { fitsPath.text = path }
    }
    FilePicker {
        id: mergePicker
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Merge OIFITS files  (Ctrl- or Shift-click: several)"
        multiSelect: true
        filters: [{ label: "OIFITS files (*.oifits *.fits)", patterns: "*.oifits,*.fits" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (paths) {
            report(Julia.shell_tools_merge(paths))
            root.datasetsChanged()
        }
    }
    FilePicker {
        id: exportPicker
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Export the file as OIFITS v2"
        saveMode: true
        defaultSuffix: "oifits"
        filters: [{ label: "OIFITS files (*.oifits *.fits)", patterns: "*.oifits *.fits" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) {
            report(Julia.shell_tools_export(filePath.text, path,
                                            targetBox.visible ? targetBox.currentText : ""))
        }
    }
    FilePicker {
        id: cubeSavePicker
        property real pxw: 1000
        property real pxh: 900
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Save the plane as PNG"
        saveMode: true
        defaultSuffix: "png"
        filters: [{ label: "PNG images (*.png)", patterns: "*.png" }, { label: "All files", patterns: "*" }]
        onAccepted: function (path) {
            Julia.shell_save_figure("cube", path, cubeSavePicker.pxw, cubeSavePicker.pxh, false)
            root.consoleChanged()
        }
    }

    // The cube on show: one label per plane, from `shell_cube_open`.
    property var cubePlanes: []
    function openCube(path) {
        var t = Julia.shell_cube_open(path)
        root.consoleChanged()
        if (failed(t)) { root.statusReported(t); return }
        root.cubePlanes = t.split("\n")
        cubePlane.value = 1
        cubeLabel.text = root.cubePlanes[0]
        root.mode = "cube"
        cubeArea.update()
    }
    function showPlane(k) {
        cubeLabel.text = Julia.shell_cube_show(k)
        cubeArea.update()
    }

    FilePicker {
        id: folderPicker
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Folder for the split files"
        chooseFolder: true
        onAccepted: function (path) { folderField.text = path }
    }

    // Controls on the left, results on the right; the split is the user's to move.
    SplitView {
        anchors.fill: parent
        anchors.margins: dp(6)
        orientation: Qt.Horizontal

        // ── controls ──────────────────────────────────────────────────────────
        ScrollView {
            // Wide enough for the widest row on opening; narrower is the user's choice.
            SplitView.preferredWidth: Math.max(dp(420), controls.implicitWidth + dp(24))
            SplitView.minimumWidth: dp(260)
            contentWidth: availableWidth
            clip: true

            ColumnLayout {
                id: controls
                width: parent.width
                spacing: dp(8)

                // ── OIFITS handling ──────────────────────────────────────────────
                GroupBox {
                    Layout.fillWidth: true
                    label: Label { text: "OIFITS handling"; font.bold: true; font.pointSize: pt(baseFontPt + 2) }
                    ColumnLayout {
                        anchors.fill: parent
                        spacing: dp(6)

                        RowLayout {
                            Layout.fillWidth: true
                            TextField {
                                id: filePath
                                objectName: "toolsFilePath"
                                Layout.fillWidth: true
                                Layout.minimumWidth: dp(120)
                                placeholderText: "an OIFITS file"
                                selectByMouse: true
                                onTextChanged: root.refreshTargets()
                            }
                            Button { text: "…"; objectName: "toolsFilePick"; onClicked: openPicker.openAt("")
                                     ToolTip.visible: hovered; ToolTip.text: "Choose a file" }
                            Button { text: "Current"; onClicked: root.useCurrentDataset(filePath)
                                     ToolTip.visible: hovered; ToolTip.text: "The current dataset's file" }
                        }
                        RowLayout {
                            visible: targetBox.visible
                            Label { text: "target"; color: "#666" }
                            ComboBox { id: targetBox; visible: false; Layout.fillWidth: true }
                        }
                        Button {
                            objectName: "toolsValidate"
                            text: "Validate"; enabled: filePath.text.length > 0
                            onClicked: validate()
                            ToolTip.visible: hovered
                            ToolTip.text: "Check against the OIFITS standard"
                        }

                        Label { text: "Spectral binning"; font.bold: true; Layout.topMargin: dp(4) }
                        RowLayout {
                            Label { text: "every"; color: "#666" }
                            SpinBox { id: binN; from: 1; to: 1000; value: 10; editable: true }
                            Label { text: "channels"; color: "#666" }
                            ComboBox {
                                id: binKind
                                model: ["average", "group"]
                                ToolTip.visible: hovered
                                ToolTip.text: "average: fewer channels · group: bins, every channel kept"
                            }
                        }
                        RowLayout {
                            visible: binKind.currentIndex === 0
                            ComboBox {
                                id: binMode
                                model: ["block", "pmoired"]
                                ToolTip.visible: hovered
                                ToolTip.text: "block: weighted means · pmoired: PMOIRED's binning"
                            }
                            Label { text: "c"; color: "#666"; visible: binMode.currentIndex === 0 }
                            TextField {
                                id: binC; text: "1"; visible: binMode.currentIndex === 0
                                Layout.preferredWidth: dp(50)
                                validator: DoubleValidator { bottom: 1 }
                                ToolTip.visible: hovered
                                ToolTip.text: "Adjacent-channel correlation (errors × √c)"
                            }
                        }
                        Button {
                            objectName: "toolsAverage"
                            text: "Bin into a new dataset"
                            enabled: filePath.text.length > 0
                            onClicked: {
                                report(binKind.currentIndex === 0
                                       ? Julia.shell_tools_bin(filePath.text, binN.value, binMode.currentText,
                                                               parseFloat(binC.text) || 1)
                                       : Julia.shell_tools_group(filePath.text, binN.value))
                                root.datasetsChanged()
                            }
                        }

                        Label { text: "Temporal binning"; font.bold: true; Layout.topMargin: dp(4) }
                        RowLayout {
                            Layout.fillWidth: true
                            ComboBox {
                                id: timeBy
                                Layout.fillWidth: true
                                textRole: "text"; valueRole: "value"
                                model: [{ text: "gaps longer than (days)", value: "gap" },
                                        { text: "windows of (days)", value: "width" }]
                            }
                            TextField {
                                id: timeParam; text: "0.5"
                                Layout.preferredWidth: dp(60)
                                validator: DoubleValidator { bottom: 0 }
                            }
                        }
                        RowLayout {
                            Layout.fillWidth: true
                            CheckBox { id: writeFiles; text: "also write each epoch to" }
                            TextField {
                                id: folderField; Layout.fillWidth: true; enabled: writeFiles.checked
                                Layout.minimumWidth: dp(60)
                                placeholderText: "folder"
                            }
                            Button { text: "…"; enabled: writeFiles.checked; onClicked: folderPicker.openAt("") }
                        }
                        RowLayout {
                            Button {
                                text: "Preview"; enabled: filePath.text.length > 0
                                onClicked: root.previewSplit(timeBy.currentValue, parseFloat(timeParam.text) || 0)
                            }
                            Button {
                                objectName: "toolsSplit"
                                text: "One dataset per epoch"; enabled: filePath.text.length > 0
                                onClicked: root.runSplit(timeBy.currentValue, parseFloat(timeParam.text) || 0,
                                                         writeFiles.checked ? folderField.text : "")
                            }
                        }

                        Label { text: "Split"; font.bold: true; Layout.topMargin: dp(4) }
                        RowLayout {
                            ComboBox {
                                id: splitBy
                                textRole: "text"; valueRole: "value"
                                model: [{ text: "by instrument", value: "instrument" },
                                        { text: "by target", value: "target" }]
                            }
                            Button { text: "Preview"; enabled: filePath.text.length > 0
                                     onClicked: root.previewSplit(splitBy.currentValue, 0) }
                            Button {
                                objectName: "toolsSplitBy"
                                text: "Split"; enabled: filePath.text.length > 0
                                onClicked: root.runSplit(splitBy.currentValue, 0, "")
                            }
                        }

                        Label { text: "Merge"; font.bold: true; Layout.topMargin: dp(4) }
                        Button {
                            objectName: "toolsMerge"
                            text: "Merge files…"
                            onClicked: mergePicker.openAt("")
                            ToolTip.visible: hovered
                            ToolTip.text: "Several OIFITS files as one dataset"
                        }

                        Button {
                            objectName: "toolsExport"
                            Layout.topMargin: dp(6)
                            text: "Export OIFITS…"
                            enabled: filePath.text.length > 0
                            onClicked: {
                                exportPicker.suggestedName = Julia.picker_suggest(
                                    filePath.text.replace(/^.*\//, "").replace(/\.(oi)?fits$/i, "") + "_v2", "oifits")
                                exportPicker.openAt("")
                            }
                            ToolTip.visible: hovered
                            ToolTip.text: "Write the file as clean OIFITS v2"
                        }
                    }
                }

                // ── FITS handling ────────────────────────────────────────────────
                GroupBox {
                    Layout.fillWidth: true
                    label: Label { text: "FITS handling"; font.bold: true; font.pointSize: pt(baseFontPt + 2) }
                    ColumnLayout {
                        anchors.fill: parent
                        spacing: dp(6)

                        RowLayout {
                            Layout.fillWidth: true
                            TextField {
                                id: fitsPath
                                objectName: "toolsFitsPath"
                                Layout.fillWidth: true
                                Layout.minimumWidth: dp(120)
                                placeholderText: "a FITS file"
                                selectByMouse: true
                            }
                            Button { text: "…"; objectName: "toolsFitsPick"; onClicked: fitsPicker.openAt("")
                                     ToolTip.visible: hovered; ToolTip.text: "Choose a file" }
                            Button { text: "Current"; onClicked: root.useCurrentDataset(fitsPath)
                                     ToolTip.visible: hovered; ToolTip.text: "The current dataset's file" }
                        }
                        RowLayout {
                            Button {
                                objectName: "toolsHeaders"
                                text: "Headers"; enabled: fitsPath.text.length > 0
                                onClicked: showHeaders()
                                ToolTip.visible: hovered
                                ToolTip.text: "Every HDU, its columns and header cards"
                            }
                            Button {
                                objectName: "toolsCubeOpen"
                                text: "Cube"; enabled: fitsPath.text.length > 0
                                onClicked: root.openCube(fitsPath.text)
                                ToolTip.visible: hovered
                                ToolTip.text: "The image, one plane at a time"
                            }
                        }
                        RowLayout {
                            visible: root.cubePlanes.length > 1
                            Label { text: "plane"; color: "#666" }
                            SpinBox {
                                id: cubePlane
                                from: 1; to: Math.max(1, root.cubePlanes.length); value: 1; editable: true
                                onValueModified: root.showPlane(value)
                            }
                            Label { text: "of " + root.cubePlanes.length; color: "#666" }
                        }
                        Slider {
                            Layout.fillWidth: true
                            visible: root.cubePlanes.length > 1
                            from: 1; to: Math.max(1, root.cubePlanes.length); stepSize: 1
                            value: cubePlane.value
                            onMoved: { cubePlane.value = value; root.showPlane(value) }
                        }
                    }
                }
                Item { Layout.fillHeight: true }
            }
        }

        // ── results ───────────────────────────────────────────────────────────
        ColumnLayout {
            SplitView.fillWidth: true
            SplitView.minimumWidth: dp(300)
            spacing: dp(4)

            Label {
                font.bold: true
                text: root.mode === "headers" ? "Headers — " + root.headersPath.replace(/^.*\//, "")
                    : root.mode === "validate" ? "Validation — " + filePath.text.replace(/^.*\//, "")
                    : root.mode === "split" ? "Split preview"
                    : root.mode === "cube" ? "Cube — " + cubeLabel.text : "Results"
            }
            Label { id: cubeLabel; visible: false }
            // The HDU list over the selected HDU's cards (Headers), the split the user's to move;
            // otherwise the report or the cube alone.
            SplitView {
                Layout.fillWidth: true
                Layout.fillHeight: true
                orientation: Qt.Vertical
            ListView {
                id: hduList
                objectName: "toolsHduList"
                visible: root.mode === "headers"
                SplitView.preferredHeight: dp(220)
                SplitView.minimumHeight: dp(60)
                clip: true
                model: ListModel { id: hduModel }
                delegate: ItemDelegate {
                    width: ListView.view.width
                    highlighted: ListView.isCurrentItem
                    text: model.idx + "   " + model.extname +
                          (model.insname.length > 0 ? "   " + model.insname : "") +
                          (model.revn !== "0" ? "   rev " + model.revn : "") +
                          (model.kind === "table" ? "   " + model.nrows + " rows, " : "   ") + model.extra
                    onClicked: { hduList.currentIndex = index; root.showHdu(model.idx) }
                }
                ScrollBar.vertical: ScrollBar {}
            }
            Rectangle {
                visible: root.mode === "cube"
                SplitView.fillHeight: root.mode === "cube"
                color: "#f4f4f4"
                border.color: "#ddd"
                MakieArea {
                    id: cubeArea
                    objectName: "toolsCubeArea"
                    anchors.fill: parent
                    scene: cubePlot
                }
                Button {
                    anchors.right: parent.right
                    anchors.top: parent.top
                    anchors.margins: dp(6)
                    text: "Save PNG"
                    opacity: hovered ? 1.0 : 0.6
                    onClicked: {
                        cubeSavePicker.pxw = Math.min(2400, Math.max(640, cubeArea.width))
                        cubeSavePicker.pxh = Math.min(1800, Math.max(480, cubeArea.height))
                        cubeSavePicker.suggestedName = Julia.picker_suggest("cube\t" + cubeLabel.text, "png")
                        cubeSavePicker.openAt("")
                    }
                }
            }
            ScrollView {
                visible: root.mode !== "cube"
                SplitView.fillHeight: root.mode !== "cube"
                SplitView.minimumHeight: dp(60)
                TextArea {
                    id: resultArea
                    objectName: "toolsResult"
                    readOnly: true
                    selectByMouse: true
                    wrapMode: TextEdit.NoWrap
                    font.family: "monospace"
                    font.pointSize: pt(baseFontPt - 1)
                }
            }
            }
        }
    }
}
