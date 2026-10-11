// OITOOLS — Explore's Filter section: what to cut from the dataset on screen.
//
// Nothing here is applied until Apply. While the settings are edited the plot shows the points
// that would go, hollow and grey, with a running count; Apply makes a NEW dataset ("data1 ›
// filtered") from the unfiltered one, logged as the OITOOLS calls that make it, and selects it.
// The source is never changed, so choosing it again in the Dataset box is the undo. Reset
// clears the settings.
//
// The section only gathers settings. They go to Julia as `key\tvalue` lines
// (`filter_spec_kwargs` in src/gui/actions_data.jl says what each key becomes), and every
// count and every cut is computed there.

import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml

ColumnLayout {
    id: root
    spacing: dp(4)

    property real uiScale:   1.0
    property real fontScale: 1.0
    property real baseFontPt: 11
    function dp(px)     { return Math.round(px * uiScale) }
    function pt(points) { return points * fontScale }

    signal statusReported(string text)
    signal applied()                 // a new dataset was made and selected
    signal previewChanged()          // the plot needs a repaint

    // ── what the dataset offers, from `shell_filter_choices` ──────────────────
    property var extentWav: ["", ""]
    property var extentBl: ["", ""]
    property var extentBlm: ["", ""]
    property var nights: []          // [{ mjd, label }]
    property var instruments: []
    property var targets: []
    property var stations: []
    property var baselines: []
    property var triangles: []
    // Ticked = kept. Keyed by name, so a list rebuilt for another dataset starts all kept.
    property var keptNights: ({})
    property var keptStations: ({})
    property var keptBaselines: ({})
    property var keptTriangles: ({})

    property string previewText: ""

    // Points picked by hand on the plot (Explore's Select section), by key: "remove_v2" ->
    // [indices into the unfiltered dataset].
    property var handPicked: ({})
    readonly property var handLabels: ({ remove_uv: "uv", remove_vis: "VIS", remove_v2: "V²",
                                         remove_t3: "T3", remove_flux: "flux" })
    function addHandPicked(line) {
        if (line.length === 0) return
        if (line.charAt(0) === "!") { root.previewText = line; return }
        var f = line.split("\t")
        if (f.length !== 2) return
        var m = Object.assign({}, root.handPicked)
        var have = m[f[0]] ? m[f[0]].slice() : []
        var add = f[1].split(",")
        for (var i = 0; i < add.length; ++i)
            if (add[i].length > 0 && have.indexOf(add[i]) < 0) have.push(add[i])
        m[f[0]] = have
        root.handPicked = m
        root.edited()
    }
    readonly property string handText: {
        var parts = []
        for (var k in root.handPicked)
            if (root.handPicked[k].length > 0) parts.push(root.handPicked[k].length + " " + root.handLabels[k])
        return parts.join(" · ")
    }
    // Set while the fields are filled programmatically, so that does not count as an edit.
    property bool loading: false

    function refreshChoices() {
        var t = Julia.shell_filter_choices()
        var ni = [], ins = [], tg = [], st = [], bl = [], tr = []
        var lines = t.length > 0 ? t.split("\n") : []
        for (var i = 0; i < lines.length; ++i) {
            var f = lines[i].split("\t")
            switch (f[0]) {
            case "extent_wav": root.extentWav = [f[1], f[2]]; break
            case "extent_bl":  root.extentBl  = [f[1], f[2]]; break
            case "extent_blm": root.extentBlm = [f[1], f[2]]; break
            case "night":      ni.push({ mjd: f[1], label: f[2] }); break
            case "insname":    ins.push(f[1]); break
            case "target":     tg.push(f[1]); break
            case "station":    st.push(f[1]); break
            case "baseline":   bl.push(f[1]); break
            case "triangle":   tr.push(f[1]); break
            }
        }
        root.nights = ni; root.instruments = ins; root.targets = tg
        root.stations = st; root.baselines = bl; root.triangles = tr
        root.keptNights = ({}); root.keptStations = ({})
        root.keptBaselines = ({}); root.keptTriangles = ({})
    }

    function isKept(map, key) { return map[key] !== false }
    function setKept(which, key, on) {
        var m = Object.assign({}, root[which]); m[key] = on; root[which] = m
        root.edited()
    }

    // ── the settings as `key\tvalue` lines ───────────────────────────────────
    readonly property var errRows: ["v2", "t3amp", "t3phi", "visamp", "visphi", "flux"]
    readonly property var errLabels: ({ v2: "V²", t3amp: "T3amp", t3phi: "T3φ (°)",
                                        visamp: "VISamp", visphi: "VISφ (°)", flux: "flux" })
    readonly property var errCols: ["add", "rel", "mult", "cap", "caprel"]
    // What each column does to an error bar, by column (errCols' order). x is the value, σ the
    // file's error bar, σ′ the one used: σ′ = max(×·σ, |x|·rel. + floor), as oifits_prep does.
    readonly property var errTips: [
        "floor: the smallest error bar allowed, in the observable's own units (degrees for a phase).\n" +
        "σ′ = max(×·σ, |x|·rel. + floor), where x is the value and σ the file's error bar",
        "relative floor: a fraction of the value, added to the floor (0.05 = 5 % of |x|).\n" +
        "σ′ = max(×·σ, |x|·rel. + floor). A phase has none",
        "multiplier: every error bar is scaled by this before the floors are compared (1 = unchanged).\n" +
        "σ′ = max(×·σ, |x|·rel. + floor)",
        "cap: points whose error bar is above this are REMOVED, in the observable's own units.\n" +
        "Compared with the file's σ, before any floor or multiplier",
        "relative cap: points whose error bar is above this fraction of |x| are REMOVED " +
        "(0.2 = σ above 20 % of the value).\nCompared with the file's σ, before any floor or " +
        "multiplier. A phase has none"
    ]
    // Field text by key ("v2_add", …), so the grid below is generated rather than written out.
    property var errText: ({})
    // Typical MIRC-X error-bar floors (V² max(5 %, 2e-4), T3amp max(10 %, 1e-5), T3φ 1°):
    // shown greyed in the empty fields, and filled in by the "MIRC-X typical" button.
    readonly property var errTypical: ({ v2_add: "2e-4", v2_rel: "0.05",
                                         t3amp_add: "1e-5", t3amp_rel: "0.1",
                                         t3phi_add: "1" })
    function fillTypical() {
        var m = Object.assign({}, root.errText)
        for (var k in root.errTypical) m[k] = root.errTypical[k]
        root.errText = m
        root.edited()
    }
    function setErr(key, text) {
        var m = Object.assign({}, root.errText); m[key] = text; root.errText = m
        root.edited()
    }

    function spec() {
        var out = []
        function put(k, v) { if (v !== undefined && v !== null && String(v).length > 0) out.push(k + "\t" + v) }
        put("wav_min", wavMin.text); put("wav_max", wavMax.text)
        put("bl_unit", blUnit.currentIndex === 1 ? "m" : "Mλ")
        put("bl_min", blMin.text); put("bl_max", blMax.text)
        var allNights = true, kept = []
        for (var i = 0; i < root.nights.length; ++i) {
            if (root.isKept(root.keptNights, root.nights[i].mjd)) kept.push(root.nights[i].mjd)
            else allNights = false
        }
        if (!allNights) out.push("nights\t" + kept.join(","))
        if (insBox.currentIndex > 0) put("insname", insBox.currentText)
        if (targetBox.currentIndex > 0) put("target", targetBox.currentText)
        if (!useV2.checked)   out.push("use_v2\t0")
        if (!useT3.checked)   out.push("use_t3\t0")
        if (!useVis.checked)  out.push("use_vis\t0")
        if (!useFlux.checked) out.push("use_flux\t0")
        function dropped(list, map) {
            var d = []
            for (var j = 0; j < list.length; ++j) if (!root.isKept(map, list[j])) d.push(list[j])
            return d.join(",")
        }
        put("exclude_stations",  dropped(root.stations,  root.keptStations))
        put("exclude_baselines", dropped(root.baselines, root.keptBaselines))
        put("exclude_triangles", dropped(root.triangles, root.keptTriangles))
        for (var r = 0; r < root.errRows.length; ++r)
            for (var c = 0; c < root.errCols.length; ++c) {
                var key = root.errRows[r] + "_" + root.errCols[c]
                put(key, root.errText[key])
            }
        if (badBox.checked) {
            out.push("bad\t1")
            put("snr", snrField.text)
            put("v2_min", v2Min.text); put("v2_max", v2Max.text)
            put("t3amp_min", t3Min.text); put("t3amp_max", t3Max.text)
        }
        if (!flagsBox.checked) out.push("ignore_flags\t1")
        for (var hk in root.handPicked)
            if (root.handPicked[hk].length > 0) out.push(hk + "\t" + root.handPicked[hk].join(","))
        return out.join("\n")
    }

    // Whether anything differs from "keep everything".
    readonly property bool hasEdits: spec().split("\n").filter(function (l) {
        return l.length > 0 && l.indexOf("bl_unit\t") !== 0 }).length > 0

    // The settings of a filtered dataset, shown again when it is selected.
    function load(text) {
        root.loading = true
        clearFields()
        var lines = text.length > 0 ? text.split("\n") : []
        var e = {}, hp = {}
        for (var i = 0; i < lines.length; ++i) {
            var f = lines[i].split("\t")
            if (f.length !== 2) continue
            var k = f[0], v = f[1]
            switch (k) {
            case "wav_min": wavMin.text = v; break
            case "wav_max": wavMax.text = v; break
            case "bl_unit": blUnit.currentIndex = v === "m" ? 1 : 0; break
            case "bl_min": blMin.text = v; break
            case "bl_max": blMax.text = v; break
            case "nights":
                var keep = v.split(","), m = {}
                for (var j = 0; j < root.nights.length; ++j)
                    m[root.nights[j].mjd] = keep.indexOf(root.nights[j].mjd) >= 0
                root.keptNights = m; break
            case "insname": insBox.currentIndex = Math.max(0, insBox.model.indexOf(v)); break
            case "target": targetBox.currentIndex = Math.max(0, targetBox.model.indexOf(v)); break
            case "use_v2": useV2.checked = v !== "0"; break
            case "use_t3": useT3.checked = v !== "0"; break
            case "use_vis": useVis.checked = v !== "0"; break
            case "use_flux": useFlux.checked = v !== "0"; break
            case "exclude_stations":  root.keptStations  = dropMap(v); break
            case "exclude_baselines": root.keptBaselines = dropMap(v); break
            case "exclude_triangles": root.keptTriangles = dropMap(v); break
            case "bad": badBox.checked = v === "1"; break
            case "snr": snrField.text = v; break
            case "v2_min": v2Min.text = v; break
            case "v2_max": v2Max.text = v; break
            case "t3amp_min": t3Min.text = v; break
            case "t3amp_max": t3Max.text = v; break
            case "ignore_flags": flagsBox.checked = v !== "1"; break
            default:
                if (k.indexOf("remove_") === 0) hp[k] = v.split(",")
                else e[k] = v
            }
        }
        root.errText = e
        root.handPicked = hp
        root.loading = false
    }
    function dropMap(v) { var m = {}; var l = v.split(","); for (var i = 0; i < l.length; ++i) m[l[i]] = false; return m }

    function clearFields() {
        wavMin.text = ""; wavMax.text = ""; blMin.text = ""; blMax.text = ""; blUnit.currentIndex = 0
        root.keptNights = ({}); root.keptStations = ({}); root.keptBaselines = ({}); root.keptTriangles = ({})
        insBox.currentIndex = 0; targetBox.currentIndex = 0
        useV2.checked = true; useT3.checked = true; useVis.checked = true; useFlux.checked = true
        root.errText = ({})
        badBox.checked = false; snrField.text = ""; v2Min.text = ""; v2Max.text = ""
        t3Min.text = ""; t3Max.text = ""; flagsBox.checked = true
        root.handPicked = ({})
    }

    // Called by the window when the dataset changes: new choices, and a filtered dataset's
    // own settings shown again.
    function refresh() {
        refreshChoices()
        load(Julia.shell_filter_spec())
        preview()
    }

    // Every edit restarts a short timer, so typing a number previews once, not per keystroke.
    function edited() { if (!root.loading) previewTimer.restart() }
    Timer { id: previewTimer; interval: 250; repeat: false; onTriggered: root.preview() }
    function preview() {
        if (!root.hasEdits) {
            root.previewText = ""
            Julia.shell_filter_preview_clear()
        } else {
            root.previewText = Julia.shell_filter_preview(root.spec())
        }
        root.previewChanged()
    }

    // ── layout ──────────────────────────────────────────────────────────────────
    component Section: Button {
        property bool open: false
        flat: true
        Layout.fillWidth: true
        contentItem: Label { text: (parent.open ? "▾ " : "▸ ") + parent.text; font.bold: true }
        onClicked: open = !open
    }

    // Ranges
    Section { id: secRanges; text: "Ranges" }
    GridLayout {
        visible: secRanges.open
        columns: 4; columnSpacing: dp(4); rowSpacing: dp(2)
        Layout.leftMargin: dp(8)
        Label { text: "λ (µm)" }
        TextField { Layout.preferredWidth: dp(64); onTextEdited: root.edited(); id: wavMin; placeholderText: root.extentWav[0] }
        Label { text: "–" }
        TextField { Layout.preferredWidth: dp(64); onTextEdited: root.edited(); id: wavMax; placeholderText: root.extentWav[1] }
        ComboBox {
            id: blUnit
            Layout.preferredWidth: dp(92)
            model: ["B (Mλ)", "B (m)"]
            onActivated: root.edited()
        }
        TextField { Layout.preferredWidth: dp(64); onTextEdited: root.edited(); id: blMin; placeholderText: blUnit.currentIndex === 1 ? root.extentBlm[0] : root.extentBl[0] }
        Label { text: "–" }
        TextField { Layout.preferredWidth: dp(64); onTextEdited: root.edited(); id: blMax; placeholderText: blUnit.currentIndex === 1 ? root.extentBlm[1] : root.extentBl[1] }
    }
    Label {
        visible: secRanges.open && root.nights.length > 1
        Layout.leftMargin: dp(8)
        text: "nights"
    }
    Repeater {
        model: secRanges.open && root.nights.length > 1 ? root.nights : []
        CheckBox {
            required property var modelData
            Layout.leftMargin: dp(16)
            text: modelData.label
            checked: root.isKept(root.keptNights, modelData.mjd)
            onToggled: root.setKept("keptNights", modelData.mjd, checked)
        }
    }

    // Tables: instrument, target and observable tables, decided when the file is read.
    Section { id: secWhat; text: "Tables" }
    GridLayout {
        visible: secWhat.open
        columns: 2; columnSpacing: dp(6); rowSpacing: dp(2)
        Layout.leftMargin: dp(8)
        Label { text: "instrument" }
        ComboBox {
            id: insBox
            Layout.fillWidth: true
            model: ["all"].concat(root.instruments)
            onActivated: root.edited()
        }
        Label { text: "target" }
        ComboBox {
            id: targetBox
            Layout.fillWidth: true
            model: ["all"].concat(root.targets)
            onActivated: root.edited()
        }
        Label { text: "tables" }
        Flow {
            Layout.fillWidth: true
            spacing: dp(2)
            CheckBox { id: useV2;   text: "V²";   checked: true; onToggled: root.edited() }
            CheckBox { id: useT3;   text: "T3";   checked: true; onToggled: root.edited() }
            CheckBox { id: useVis;  text: "VIS";  checked: true; onToggled: root.edited() }
            CheckBox { id: useFlux; text: "flux"; checked: true; onToggled: root.edited() }
        }
    }

    // Telescopes and baselines, from the file: unticking a station removes everything using it.
    Section { id: secTel; text: "Telescopes & baselines" }
    ColumnLayout {
        visible: secTel.open
        Layout.leftMargin: dp(8)
        spacing: 0
        Label { text: "stations"; color: "#666" }
        Flow {
            Layout.fillWidth: true
            Repeater {
                model: root.stations
                CheckBox {
                    required property string modelData
                    text: modelData
                    checked: root.isKept(root.keptStations, modelData)
                    onToggled: root.setKept("keptStations", modelData, checked)
                }
            }
        }
        Label { text: "baselines"; color: "#666" }
        Flow {
            Layout.fillWidth: true
            Repeater {
                model: root.baselines
                CheckBox {
                    required property string modelData
                    text: modelData
                    checked: root.isKept(root.keptBaselines, modelData)
                    onToggled: root.setKept("keptBaselines", modelData, checked)
                }
            }
        }
        Label { text: "triangles"; color: "#666"; visible: root.triangles.length > 0 }
        Flow {
            Layout.fillWidth: true
            Repeater {
                model: root.triangles
                CheckBox {
                    required property string modelData
                    text: modelData
                    checked: root.isKept(root.keptTriangles, modelData)
                    onToggled: root.setKept("keptTriangles", modelData, checked)
                }
            }
        }
    }

    // Error bars: these CHANGE the data rather than remove it (floors, multipliers), apart
    // from the caps, which remove points whose error is above them.
    Section { id: secErr; text: "Error bars" }
    GridLayout {
        visible: secErr.open
        columns: 6; columnSpacing: dp(3); rowSpacing: dp(1)
        Layout.leftMargin: dp(8)
        Label { text: "" }
        Repeater {
            model: ["floor", "rel.", "×", "cap", "cap rel."]
            Label {
                required property string modelData
                required property int index
                text: modelData; color: "#666"; font.pointSize: pt(baseFontPt - 2)
                HoverHandler { id: headHover }
                ToolTip.visible: headHover.hovered
                ToolTip.text: root.errTips[index]
            }
        }
        Repeater {
            model: root.errRows.length * 6
            Item {
                id: cell
                required property int index
                readonly property int r: Math.floor(index / 6)
                readonly property int c: index % 6
                readonly property string obs: root.errRows[r]
                readonly property bool phase: obs === "t3phi" || obs === "visphi"
                implicitWidth: c === 0 ? lab.implicitWidth : fld.implicitWidth
                implicitHeight: fld.implicitHeight
                Label { id: lab; visible: cell.c === 0; text: root.errLabels[cell.obs]; anchors.verticalCenter: parent.verticalCenter }
                TextField {
                    id: fld
                    visible: cell.c > 0
                    width: root.dp(44)
                    implicitWidth: root.dp(44)
                    readonly property string key: cell.obs + "_" + root.errCols[Math.max(0, cell.c - 1)]
                    // A phase has no relative floor or relative cap.
                    enabled: !(cell.phase && (cell.c === 2 || cell.c === 5))
                    text: root.errText[key] !== undefined ? root.errText[key] : ""
                    placeholderText: root.errTypical[key] !== undefined ? root.errTypical[key] : ""
                    onTextEdited: root.setErr(key, text)
                    ToolTip.visible: hovered && enabled
                    ToolTip.text: root.errTips[Math.max(0, cell.c - 1)]
                }
            }
        }
    }

    Button {
        visible: secErr.open
        Layout.leftMargin: dp(8)
        text: "MIRC-X typical"
        ToolTip.visible: hovered
        ToolTip.text: "fill in the floors shown greyed: V² max(5 %, 0.0002), T3amp max(10 %, 1e-5), T3φ 1°"
        onClicked: root.fillTypical()
    }

    // Quality: the file's flags and the usual cuts.
    Section { id: secQual; text: "Quality" }
    GridLayout {
        visible: secQual.open
        columns: 4; columnSpacing: dp(4); rowSpacing: dp(2)
        Layout.leftMargin: dp(8)
        CheckBox {
            id: flagsBox
            Layout.columnSpan: 4
            text: "honour the file's flags"
            checked: true
            onToggled: root.edited()
            ToolTip.visible: hovered
            ToolTip.text: "unticked: flagged points are kept unless another cut removes them (re-reads the file)"
        }
        CheckBox {
            id: badBox
            Layout.columnSpan: 4
            text: "drop bad points"
            onToggled: root.edited()
            ToolTip.visible: hovered
            ToolTip.text: "NaN or non-positive errors, and the cuts below"
        }
        Label { text: "V² SNR ≥"; enabled: badBox.checked }
        TextField { Layout.preferredWidth: dp(64); onTextEdited: root.edited(); id: snrField; enabled: badBox.checked; placeholderText: "0.01" }
        Item { Layout.columnSpan: 2 }
        Label { text: "V² range"; enabled: badBox.checked }
        TextField { Layout.preferredWidth: dp(64); onTextEdited: root.edited(); id: v2Min; enabled: badBox.checked; placeholderText: "-1" }
        Label { text: "–"; enabled: badBox.checked }
        TextField { Layout.preferredWidth: dp(64); onTextEdited: root.edited(); id: v2Max; enabled: badBox.checked; placeholderText: "2" }
        Label { text: "T3amp range"; enabled: badBox.checked }
        TextField { Layout.preferredWidth: dp(64); onTextEdited: root.edited(); id: t3Min; enabled: badBox.checked; placeholderText: "-1" }
        Label { text: "–"; enabled: badBox.checked }
        TextField { Layout.preferredWidth: dp(64); onTextEdited: root.edited(); id: t3Max; enabled: badBox.checked; placeholderText: "1.5" }
    }

    // Spectral binning stays where the polychromatic reconstruction is set up.
    Section { id: secSpec; text: "Spectral" }
    Label {
        visible: secSpec.open
        Layout.leftMargin: dp(8)
        Layout.fillWidth: true
        wrapMode: Text.WordWrap
        text: "Binning by wavelength is in Imaging › Setup."
        color: "#666"
    }

    // Points picked by hand on the plot (Select, below), shown here because they are part of
    // what Apply removes.
    RowLayout {
        visible: root.handText.length > 0
        Layout.fillWidth: true
        Label {
            Layout.fillWidth: true
            wrapMode: Text.WordWrap
            text: "hand-picked: " + root.handText
            font.pointSize: pt(baseFontPt - 1)
        }
        Button {
            text: "clear"
            implicitHeight: dp(22)
            onClicked: { root.handPicked = ({}); root.edited() }
        }
    }

    // ── the count, Apply and Reset ─────────────────────────────────────────────
    Label {
        Layout.fillWidth: true
        Layout.topMargin: dp(4)
        wrapMode: Text.WordWrap
        text: root.previewText.length > 0 ? root.previewText : "nothing filtered"
        color: root.previewText.indexOf("!") === 0 ? "#b00020" : "#444"
        font.pointSize: pt(baseFontPt - 1)
    }
    RowLayout {
        Button {
            text: "Apply filter"
            enabled: root.hasEdits && root.previewText.indexOf("!") !== 0
            onClicked: {
                root.statusReported(Julia.shell_apply_filter(root.spec()))
                root.applied()
            }
        }
        Button {
            text: "Reset"
            enabled: root.hasEdits
            ToolTip.visible: hovered
            ToolTip.text: "clear these settings; choosing the unfiltered dataset in the Dataset box is the undo"
            onClicked: { root.loading = true; root.clearFields(); root.loading = false; root.preview() }
        }
    }
}
