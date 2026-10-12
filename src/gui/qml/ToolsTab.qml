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
    signal exportRequested()                // the window owns the Export picker
    signal statusReported(string text)

    // What the result pane shows: "", "headers" or "validate".
    property string mode: ""
    property string headersPath: ""

    function failed(t) { return t.length > 0 && t[0] === "!" }
    function report(t) { root.statusReported(t); root.consoleChanged() }

    function useCurrentDataset() {
        var p = Julia.shell_tools_source()
        if (p.length === 0) { root.statusReported("! no single-file dataset is selected"); return }
        filePath.text = p
    }

    function refreshTargets() {
        var t = filePath.text.length > 0 ? Julia.shell_tools_targets(filePath.text) : ""
        var names = (t.length > 0 && !failed(t)) ? t.split("\n") : []
        targetBox.model = names
        targetBox.visible = names.length > 1
    }

    function showHeaders() {
        if (filePath.text.length === 0) return
        var t = Julia.shell_tools_headers(filePath.text)
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
        root.headersPath = filePath.text
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
    FilePicker {
        id: openPicker
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Choose a FITS / OIFITS file"
        filters: [{ label: "FITS files (*.oifits *.fits *.fit)", patterns: "*.oifits,*.fits,*.fit" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) { filePath.text = path }
    }
    FilePicker {
        id: repairPicker
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Write the repaired copy as"
        saveMode: true
        defaultSuffix: "oifits"
        filters: [{ label: "OIFITS files (*.oifits *.fits)", patterns: "*.oifits *.fits" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) {
            report(Julia.shell_tools_repair(filePath.text, path,
                                            targetBox.visible ? targetBox.currentText : ""))
        }
    }
    FilePicker {
        id: folderPicker
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Folder for the split files"
        chooseFolder: true
        onAccepted: function (path) { folderField.text = path }
    }

    RowLayout {
        anchors.fill: parent
        anchors.margins: dp(6)
        spacing: dp(8)

        // ── controls ──────────────────────────────────────────────────────────
        ScrollView {
            Layout.preferredWidth: dp(400)
            Layout.fillHeight: true
            contentWidth: availableWidth
            clip: true

            ColumnLayout {
                width: parent.width
                spacing: dp(6)

                Label { text: "File"; font.bold: true; font.pointSize: pt(baseFontPt + 1) }
                Label {
                    Layout.fillWidth: true; wrapMode: Text.WordWrap
                    color: "#666"; font.pointSize: pt(baseFontPt - 1)
                    text: "Any FITS file: its headers, or a check of its OIFITS tables."
                }
                RowLayout {
                    Layout.fillWidth: true
                    TextField {
                        id: filePath
                        objectName: "toolsFilePath"
                        Layout.fillWidth: true
                        placeholderText: "a FITS / OIFITS file"
                        selectByMouse: true
                        onTextChanged: refreshTargets()
                    }
                    Button { text: "…"; onClicked: openPicker.openAt("")
                             ToolTip.visible: hovered; ToolTip.text: "Choose a file" }
                }
                Button {
                    text: "Use the current dataset's file"
                    onClicked: useCurrentDataset()
                }
                RowLayout {
                    Button {
                        objectName: "toolsHeaders"
                        text: "Headers"; enabled: filePath.text.length > 0
                        onClicked: showHeaders()
                        ToolTip.visible: hovered
                        ToolTip.text: "Every HDU, the ones OITOOLS does not read included; select one for its columns and cards"
                    }
                    Button {
                        objectName: "toolsValidate"
                        text: "Validate"; enabled: filePath.text.length > 0
                        onClicked: validate()
                        ToolTip.visible: hovered
                        ToolTip.text: "Check the file against the OIFITS standard and its data for suspect values"
                    }
                    Button {
                        text: "Repaired copy…"; enabled: filePath.text.length > 0
                        onClicked: {
                            repairPicker.suggestedName = Julia.picker_suggest(
                                filePath.text.replace(/^.*\//, "").replace(/\.(oi)?fits$/i, "") + "_repaired", "oifits")
                            repairPicker.openAt("")
                        }
                        ToolTip.visible: hovered
                        ToolTip.text: "Read the file and write it again as clean OIFITS v2 (instruments kept). " +
                                      "Not kept: non-OIFITS HDUs, columns OITOOLS does not read, flagged points, OI_CORR."
                    }
                }
                RowLayout {
                    visible: targetBox.visible
                    Label { text: "Target"; color: "#666" }
                    ComboBox { id: targetBox; visible: false; Layout.fillWidth: true }
                }

                Rectangle { Layout.fillWidth: true; height: 1; color: "#ddd"; Layout.topMargin: dp(6) }

                Label { text: "Current dataset"; font.bold: true; font.pointSize: pt(baseFontPt + 1) }
                Label {
                    Layout.fillWidth: true; wrapMode: Text.WordWrap
                    color: "#666"; font.pointSize: pt(baseFontPt - 1)
                    text: "Averaging and splitting read the selected dataset's file again; each " +
                          "result is a new dataset, listed in the Dataset box."
                }

                Label { text: "Average channels"; font.bold: true }
                RowLayout {
                    Label { text: "every"; color: "#666" }
                    SpinBox { id: binN; from: 2; to: 1000; value: 10; editable: true }
                    Label { text: "channels"; color: "#666" }
                }
                RowLayout {
                    ComboBox {
                        id: binMode
                        model: ["block", "pmoired"]
                        ToolTip.visible: hovered
                        ToolTip.text: "block: 1/σ²-weighted means, circular-mean phases, σ·√c errors and a scatter check. " +
                                      "pmoired: PMOIRED's resampling, to reproduce its results."
                    }
                    Label { text: "c"; color: "#666"; visible: binMode.currentIndex === 0 }
                    TextField {
                        id: binC; text: "1"; visible: binMode.currentIndex === 0
                        Layout.preferredWidth: dp(50)
                        validator: DoubleValidator { bottom: 1 }
                        ToolTip.visible: hovered
                        ToolTip.text: "Correlation between adjacent channels: errors are multiplied by √c (1 = independent)"
                    }
                }
                Button {
                    objectName: "toolsAverage"
                    text: "Average into a new dataset"
                    onClicked: {
                        report(Julia.shell_tools_bin(binN.value, binMode.currentText, parseFloat(binC.text) || 1))
                        root.datasetsChanged()
                    }
                    ToolTip.visible: hovered
                    ToolTip.text: "Averaging makes fewer channels, each the mean of the ones it replaces. " +
                                  "Not the same as grouping (Imaging › Setup), which keeps every channel in fewer bins."
                }

                Label { text: "Split"; font.bold: true; Layout.topMargin: dp(4) }
                RowLayout {
                    ComboBox {
                        id: splitBy
                        Layout.fillWidth: true
                        textRole: "text"; valueRole: "value"
                        model: [{ text: "by time: gaps longer than (days)", value: "gap" },
                                { text: "by time: windows of (days)", value: "width" },
                                { text: "by instrument", value: "instrument" },
                                { text: "by target", value: "target" }]
                    }
                    TextField {
                        id: splitParam; text: "0.5"
                        Layout.preferredWidth: dp(60)
                        enabled: splitBy.currentValue === "gap" || splitBy.currentValue === "width"
                        validator: DoubleValidator { bottom: 0 }
                    }
                }
                RowLayout {
                    CheckBox { id: writeFiles; text: "also write each part to" }
                    TextField {
                        id: folderField; Layout.fillWidth: true; enabled: writeFiles.checked
                        placeholderText: "folder"
                    }
                    Button { text: "…"; enabled: writeFiles.checked; onClicked: folderPicker.openAt("") }
                }
                RowLayout {
                    Button {
                        text: "Preview"
                        onClicked: {
                            var t = Julia.shell_tools_split_preview(splitBy.currentValue, parseFloat(splitParam.text) || 0)
                            root.mode = "split"
                            resultArea.text = failed(t) ? t : t.split("\n").length + " parts:\n" + t.split("\t").join("    ")
                        }
                    }
                    Button {
                        objectName: "toolsSplit"
                        text: "Split into datasets"
                        onClicked: {
                            report(Julia.shell_tools_split(splitBy.currentValue, parseFloat(splitParam.text) || 0,
                                                           writeFiles.checked ? folderField.text : ""))
                            root.datasetsChanged()
                        }
                    }
                }

                Rectangle { Layout.fillWidth: true; height: 1; color: "#ddd"; Layout.topMargin: dp(6) }
                Button {
                    text: "Export the current dataset as OIFITS…"
                    onClicked: root.exportRequested()
                }
                Item { Layout.fillHeight: true }
            }
        }

        // ── results ───────────────────────────────────────────────────────────
        ColumnLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: dp(4)

            Label {
                font.bold: true
                text: root.mode === "headers" ? "Headers — " + root.headersPath.replace(/^.*\//, "")
                    : root.mode === "validate" ? "Validation — " + filePath.text.replace(/^.*\//, "")
                    : root.mode === "split" ? "Split preview" : "Results"
            }
            ListView {
                id: hduList
                objectName: "toolsHduList"
                visible: root.mode === "headers"
                Layout.fillWidth: true
                Layout.preferredHeight: Math.min(contentHeight, dp(220))
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
            ScrollView {
                Layout.fillWidth: true
                Layout.fillHeight: true
                TextArea {
                    id: resultArea
                    objectName: "toolsResult"
                    readOnly: true
                    selectByMouse: true
                    wrapMode: TextEdit.NoWrap
                    font.family: "monospace"
                    font.pointSize: pt(baseFontPt - 1)
                    placeholderText: "Headers, validation reports and split previews appear here."
                }
            }
        }
    }
}
