// OITOOLS — the Tools perspective: operations on files that are none of the other four.
//
// OIFITS Handling works on the file in its field: Headers, Validate, spectral and temporal
// binning, the splits and Merge (each a re-read that becomes a new dataset, so Explore,
// Modeling and Imaging use it at once), and Export OIFITS. FITS Image Handling works on the file
// in its own field: Headers, the cube viewer, Set keywords and Save FITS. Everything is done by
// `src/gui/tools.jl`; this file only lays it out.

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
        showHeaders(p)
    }

    function showHeaders(path) {
        if (path.length === 0) return
        var t = Julia.shell_tools_headers(path)
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
        root.headersPath = path
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
    // One picker per file field: OIFITS Handling takes OIFITS files and shows their headers, FITS
    // Image Handling takes FITS files and views them.
    FilePicker {
        id: openPicker
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Choose an OIFITS file"
        filters: [{ label: "OIFITS files (*.oifits *.fits *.fit)", patterns: "*.oifits,*.fits,*.fit" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) { filePath.text = path; root.showHeaders(path) }
    }
    FilePicker {
        id: fitsPicker
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Choose a FITS file"
        filters: [{ label: "FITS files (*.fits)", patterns: "*.fits" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) { fitsPath.text = path; root.openCube(path) }
    }
    FilePicker {
        id: fitsSavePicker
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Save FITS"
        saveMode: true
        defaultSuffix: "fits"
        filters: [{ label: "FITS files (*.fits)", patterns: "*.fits" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) {
            report(Julia.shell_fits_save(path))
            root.keywordSource = Julia.shell_fits_keyword_source()
            if (path === fitsPath.text) root.openCube(path)
        }
    }

    // The file Set keywords was applied to; Save FITS writes it.
    property string keywordSource: ""

    // ── Set keywords ──────────────────────────────────────────────────────────
    Dialog {
        id: keywordsDialog
        objectName: "keywordsDialog"
        title: "Set keywords"
        modal: true
        anchors.centerIn: Overlay.overlay
        width: root.dp(460)
        property string path: ""
        property var dims: []
        property string problem: ""
        // the file's WCS cards: { name, ok, value }
        property var cards: []
        readonly property int n3: dims.length > 2 ? dims[2] : 1
        readonly property int n4: dims.length > 3 ? dims[3] : 1
        // planes along wavelength and time, after the axis-3 choice
        readonly property int nWave: kwAxis3.currentIndex === 0 ? n3 : n4
        readonly property int nTime: kwAxis3.currentIndex === 0 ? n4 : n3
        function openFor(p) {
            var t = Julia.shell_fits_keywords(p)
            if (root.failed(t)) { root.statusReported(t); return }
            var kv = {}, cs = []
            var lines = t.split("\n")
            for (var i = 0; i < lines.length; ++i) {
                var f = lines[i].split("\t")
                if (f[0] === "card") cs.push({ name: f[1], ok: f[2] === "1", value: f.length > 3 ? f[3] : "" })
                else kv[f[0]] = f.length > 1 ? f[1] : ""
            }
            cards = cs
            path = p
            dims = kv["dims"].split(" ").map(Number)
            kwPix.text = kv["pixsize"]; kwPixUnit.currentIndex = 0
            kwEast.checked = kv["eastleft"] === "1"
            kwAxis3.currentIndex = kv["axis3"] === "time" ? 1 : 0
            kwWStart.text = kv["wstart"]; kwWStep.text = kv["wstep"]; kwWUnit.currentIndex = 1
            kwTStart.text = kv["tstart"]; kwTStep.text = kv["tstep"]; kwTUnit.currentIndex = 0
            kwRefX.text = kv["refx"]; kwRefY.text = kv["refy"]
            problem = ""
            open()
        }
        GridLayout {
            width: root.dp(420)
            columns: 5
            columnSpacing: root.dp(6)
            rowSpacing: root.dp(6)

            Label { text: "size"; color: "#666" }
            Label { Layout.columnSpan: 4; text: keywordsDialog.dims.join(" × ") }

            Label { text: "pixel"; color: "#666" }
            TextField { id: kwPix; objectName: "kwPix"; Layout.preferredWidth: root.dp(110) }
            ComboBox { id: kwPixUnit; objectName: "kwPixUnit"; Layout.columnSpan: 2
                       Layout.preferredWidth: root.dp(90); model: ["mas", "µas", "arcsec", "deg", "rad"] }
            CheckBox {
                id: kwEast; objectName: "kwEast"; text: "East left"
                ToolTip.visible: hovered; ToolTip.text: "The file stores East on the left"
            }

            Label { text: "ref. pixel"; color: "#666" }
            TextField {
                id: kwRefX; objectName: "kwRefX"; Layout.preferredWidth: root.dp(110)
                placeholderText: "centre"
                ToolTip.visible: hovered; ToolTip.text: "Origin, x (CRPIX1)"
            }
            TextField {
                id: kwRefY; objectName: "kwRefY"; Layout.columnSpan: 3; Layout.preferredWidth: root.dp(90)
                placeholderText: "centre"
                ToolTip.visible: hovered; ToolTip.text: "Origin, y (CRPIX2)"
            }

            Label { text: "axis 3"; color: "#666"; visible: keywordsDialog.dims.length > 2 }
            ComboBox {
                id: kwAxis3; objectName: "kwAxis3"
                visible: keywordsDialog.dims.length > 2
                Layout.columnSpan: 4
                Layout.fillWidth: true
                model: keywordsDialog.dims.length > 3 ? ["wavelength (axis 4: time)", "time (axis 4: wavelength)"]
                                                      : ["wavelength", "time"]
            }

            Label { text: "λ start"; color: "#666"; visible: keywordsDialog.nWave > 1 }
            TextField { id: kwWStart; objectName: "kwWStart"; Layout.preferredWidth: root.dp(110)
                        visible: keywordsDialog.nWave > 1 }
            Label { text: "step"; color: "#666"; visible: keywordsDialog.nWave > 1 }
            TextField { id: kwWStep; objectName: "kwWStep"; Layout.preferredWidth: root.dp(90)
                        visible: keywordsDialog.nWave > 1 }
            ComboBox { id: kwWUnit; objectName: "kwWUnit"; Layout.preferredWidth: root.dp(80)
                       visible: keywordsDialog.nWave > 1; model: ["m", "µm", "nm", "Å"] }

            Label { text: "MJD start"; color: "#666"; visible: keywordsDialog.nTime > 1 }
            TextField { id: kwTStart; objectName: "kwTStart"; Layout.preferredWidth: root.dp(110)
                        visible: keywordsDialog.nTime > 1 }
            Label { text: "step"; color: "#666"; visible: keywordsDialog.nTime > 1 }
            TextField { id: kwTStep; objectName: "kwTStep"; Layout.preferredWidth: root.dp(90)
                        visible: keywordsDialog.nTime > 1 }
            ComboBox { id: kwTUnit; objectName: "kwTUnit"; Layout.preferredWidth: root.dp(80)
                       visible: keywordsDialog.nTime > 1; model: ["d", "h", "min", "s"] }

            // What the file's header has now: ✓ present (its value on hover), ✗ missing.
            Label { text: "WCS in the file"; color: "#666"; Layout.columnSpan: 5; Layout.topMargin: root.dp(4) }
            // one row per axis: CTYPE CUNIT CRVAL CDELT CRPIX
            Grid {
                objectName: "kwCards"
                Layout.columnSpan: 5
                columns: 5
                columnSpacing: root.dp(12)
                rowSpacing: root.dp(2)
                Repeater {
                    model: keywordsDialog.cards
                    Label {
                        required property var modelData
                        text: (modelData.ok ? "\u2713 " : "\u2717 ") + modelData.name
                        color: modelData.ok ? "#1a7f37" : "#b00020"
                        font.pointSize: root.pt(root.baseFontPt - 1)
                        HoverHandler { id: cardHover }
                        ToolTip.visible: cardHover.hovered && modelData.ok
                        ToolTip.text: modelData.value
                    }
                }
            }

            Label {
                Layout.columnSpan: 5
                Layout.preferredWidth: root.dp(420)
                visible: keywordsDialog.problem.length > 0
                text: keywordsDialog.problem
                color: "#b00020"
                wrapMode: Text.WordWrap
            }
        }
        // Apply is not an accept role: the dialog closes only when the values are taken, and
        // otherwise stays open saying what is wrong.
        footer: DialogButtonBox {
            Button { text: "Apply"; onClicked: keywordsDialog.apply() }
            Button { text: "Cancel"; DialogButtonBox.buttonRole: DialogButtonBox.RejectRole }
        }
        function apply() {
            var t = Julia.shell_fits_apply(path, kwPix.text, kwPixUnit.currentText, kwEast.checked,
                                           kwAxis3.currentIndex === 1 ? "time" : "wave",
                                           kwWStart.text, kwWStep.text, kwWUnit.currentText,
                                           kwTStart.text, kwTStep.text, kwTUnit.currentText,
                                           kwRefX.text, kwRefY.text)
            root.consoleChanged()
            if (root.failed(t)) { problem = t.substring(2); return }
            problem = ""
            close()
            root.keywordSource = Julia.shell_fits_keyword_source()
            root.refreshCubeAxes()
        }
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
            Julia.shell_save_figure(root.cubeTiles ? "cubetiles" : "cube", path, cubeSavePicker.pxw,
                                    cubeSavePicker.pxh, false)
            root.consoleChanged()
        }
    }

    // The cube on show: its wavelength channels and epochs (`shell_cube_axes`), each with its
    // value where the file gives one (µm, MJD), and the channel and epoch on the canvas.
    property int cubeNW: 0
    property int cubeNT: 0
    property var cubeWL: []
    property var cubeMJ: []
    property int cubeW: 1
    property int cubeT: 1
    // channels side by side instead of one; the list typed (channels or µm), the labels shown
    property bool cubeTiles: false
    property string cubeTileText: ""
    function openCube(path) {
        var t = Julia.shell_cube_open(path)
        root.consoleChanged()
        if (failed(t)) { root.statusReported(t); return }
        refreshCubeAxes()
    }
    // After Julia has opened a cube (or set its keywords): read its axes, show the first plane.
    function refreshCubeAxes() {
        var t = Julia.shell_cube_axes()
        if (failed(t)) return
        var l = t.split("\n"), n = l[0].split("\t")
        cubeNW = parseInt(n[0]); cubeNT = parseInt(n[1])
        cubeWL = l[1].split("\t"); cubeMJ = l[2].split("\t")
        root.mode = "cube"
        showAt(1, 1)
    }
    function showAt(w, t) {
        cubeW = Math.max(1, Math.min(cubeNW, w)); cubeT = Math.max(1, Math.min(cubeNT, t))
        cubeLabel.text = Julia.shell_cube_show_at(cubeW, cubeT)
        cubeArea.update()
        if (cubeTiles) showTiles()
    }
    function showTiles() {
        var t = Julia.shell_cube_tiles(tileList.text, cubeT)
        if (failed(t)) { root.statusReported(t); return }
        cubeTileText = t.split("\t").join("  ·  ")
        tilesArea.update()
    }

    // One channel axis: ◀ slider ▶, the index and the channel's value (µm or MJD) when known.
    component AxisSlider: RowLayout {
        id: ax
        property string name: "Channel"
        property int n: 1
        property int value: 1
        property string info: ""
        signal stepped(int k)
        spacing: root.dp(2)
        // as wide as the longer name, so the two sliders line up
        Label { text: ax.name; color: "#666"; Layout.preferredWidth: nameWidth.advanceWidth("Channel") + root.dp(4) }
        FontMetrics { id: nameWidth }
        ToolButton { text: "◀"; enabled: ax.value > 1; onClicked: ax.stepped(ax.value - 1) }
        Slider {
            Layout.fillWidth: true
            Layout.minimumWidth: root.dp(60)
            from: 1; to: Math.max(1, ax.n); stepSize: 1; snapMode: Slider.SnapAlways
            value: ax.value
            onMoved: ax.stepped(Math.round(value))
        }
        ToolButton { text: "▶"; enabled: ax.value < ax.n; onClicked: ax.stepped(ax.value + 1) }
        Label { text: ax.value + "/" + ax.n; color: "#666" }
        Label { text: ax.info; Layout.preferredWidth: root.dp(116); Layout.leftMargin: root.dp(6) }
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

                // ── OIFITS Handling ──────────────────────────────────────────────
                GroupBox {
                    Layout.fillWidth: true
                    label: Label { text: "OIFITS Handling"; font.bold: true; font.pointSize: pt(baseFontPt + 2) }
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
                                onAccepted: root.showHeaders(text)
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
                        RowLayout {
                            Button {
                                objectName: "toolsOifitsHeaders"
                                text: "Headers"; enabled: filePath.text.length > 0
                                onClicked: showHeaders(filePath.text)
                                ToolTip.visible: hovered
                                ToolTip.text: "Every HDU, its columns and header cards"
                            }
                            Button {
                                objectName: "toolsValidate"
                                text: "Validate"; enabled: filePath.text.length > 0
                                onClicked: validate()
                                ToolTip.visible: hovered
                                ToolTip.text: "Check against the OIFITS standard"
                            }
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

                // ── FITS Image Handling ────────────────────────────────────────────────
                GroupBox {
                    Layout.fillWidth: true
                    label: Label { text: "FITS Image Handling"; font.bold: true; font.pointSize: pt(baseFontPt + 2) }
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
                                onAccepted: root.openCube(text)
                            }
                            Button { text: "…"; objectName: "toolsFitsPick"; onClicked: fitsPicker.openAt("")
                                     ToolTip.visible: hovered; ToolTip.text: "Choose a file" }
                            Button {
                                objectName: "toolsHeaders"
                                text: "Headers"; enabled: fitsPath.text.length > 0
                                onClicked: showHeaders(fitsPath.text)
                                ToolTip.visible: hovered
                                ToolTip.text: "Every HDU, its columns and header cards"
                            }
                        }
                        // View, beside the channel slider: wavelength channels when the cube
                        // has several, else its epochs; a cube with both gets an Epoch slider too.
                        RowLayout {
                            Button {
                                id: viewButton
                                objectName: "toolsView"
                                text: "View"; enabled: fitsPath.text.length > 0
                                onClicked: root.openCube(fitsPath.text)
                                ToolTip.visible: hovered
                                ToolTip.text: "The image, one channel at a time"
                            }
                            AxisSlider {
                                objectName: "toolsChannel"
                                Layout.fillWidth: true
                                visible: root.cubeNW > 1 || root.cubeNT > 1
                                readonly property bool wave: root.cubeNW > 1
                                n: wave ? root.cubeNW : root.cubeNT
                                value: wave ? root.cubeW : root.cubeT
                                info: (wave ? root.cubeWL : root.cubeMJ)[value - 1] || ""
                                onStepped: function (k) { wave ? root.showAt(k, root.cubeT) : root.showAt(1, k) }
                            }
                        }
                        // channels side by side, and the stretch of every plane shown
                        RowLayout {
                            visible: root.mode === "cube"
                            CheckBox {
                                id: tilesBox
                                objectName: "toolsTiles"
                                text: "Side by side"
                                checked: root.cubeTiles
                                onToggled: { root.cubeTiles = checked; if (checked) root.showTiles() }
                                ToolTip.visible: hovered; ToolTip.text: "Several channels at once"
                            }
                            TextField {
                                id: tileList
                                objectName: "toolsTileList"
                                visible: root.cubeTiles
                                Layout.fillWidth: true
                                placeholderText: "4 channels"
                                onEditingFinished: root.showTiles()
                                ToolTip.visible: hovered; ToolTip.text: "Channels, or λ in µm (up to 6)"
                            }
                            Label { text: "power"; color: "#666" }
                            TextField {
                                objectName: "toolsPower"
                                Layout.preferredWidth: root.dp(60)
                                text: "1"
                                onEditingFinished: {
                                    var t = Julia.shell_cube_power(text)
                                    if (root.failed(t)) { root.statusReported(t); return }
                                    root.showAt(root.cubeW, root.cubeT)
                                }
                                ToolTip.visible: hovered; ToolTip.text: "Brightness stretch, I^p"
                            }
                        }
                        AxisSlider {
                            objectName: "toolsEpoch"
                            Layout.fillWidth: true
                            // under the Channel slider, past the View button
                            Layout.leftMargin: viewButton.width + root.dp(5)
                            visible: root.cubeNW > 1 && root.cubeNT > 1
                            name: "Epoch"
                            n: root.cubeNT
                            value: root.cubeT
                            info: root.cubeMJ[value - 1] || ""
                            onStepped: function (k) { root.showAt(root.cubeW, k) }
                        }
                        RowLayout {
                            Button {
                                objectName: "toolsSetKeywords"
                                text: "Set keywords…"; enabled: fitsPath.text.length > 0
                                onClicked: keywordsDialog.openFor(fitsPath.text)
                                ToolTip.visible: hovered
                                ToolTip.text: "Pixel size, orientation, wavelength and time axes"
                            }
                            Button {
                                objectName: "toolsSaveFits"
                                text: "Save FITS…"
                                enabled: root.keywordSource.length > 0 && root.keywordSource === fitsPath.text
                                onClicked: {
                                    fitsSavePicker.suggestedName = fitsPath.text.replace(/^.*\//, "")
                                    fitsSavePicker.openAt(fitsPath.text)
                                }
                                ToolTip.visible: hovered
                                ToolTip.text: "Write the file with the keywords set"
                            }
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
                    : root.mode === "cube" ? "Cube — " + (root.cubeTiles ? root.cubeTileText : cubeLabel.text)
                    : "Results"
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
                    visible: !root.cubeTiles
                    scene: cubePlot
                }
                MakieArea {
                    id: tilesArea
                    objectName: "toolsTilesArea"
                    anchors.fill: parent
                    visible: root.cubeTiles
                    scene: cubeTilesPlot
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
