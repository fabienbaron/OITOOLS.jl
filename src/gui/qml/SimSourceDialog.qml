// What one target is simulated from, in one menu: none (not simulated), a model held in the
// Modeling tab, a uniform disc, a FITS image, or a model file.

import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml

Dialog {
    id: dlg

    property real uiScale: 1.0
    property real fontScale: 1.0
    property real baseFontPt: 10
    function dp(px) { return Math.round(px * uiScale) }

    property int    row: -1
    property string targetName: ""
    property var    modelNames: []
    property string kind: ""
    property string ref: ""
    property real   pixsize: 0
    property bool   channels: false
    property real   diam: 0.5
    property bool   usable: true
    property string statusText: ""

    signal chosen(int row, string kind, string ref, real pixsize, bool channels, real diam)

    // One menu: none, the Modeling tab's models by name, then the other kinds.
    readonly property var choices: {
        var c = [{ label: "None", kind: "", ref: "" }]
        for (var i = 0; i < modelNames.length; i++)
            c.push({ label: modelNames[i], kind: "session", ref: modelNames[i] })
        c.push({ label: "Uniform disc…", kind: "ud", ref: "" })
        c.push({ label: "FITS image…", kind: "image", ref: "" })
        c.push({ label: "Model file…", kind: "modelfile", ref: "" })
        return c
    }
    function choiceIndex() {
        for (var i = 0; i < choices.length; i++)
            if (choices[i].kind === kind && (kind !== "session" || choices[i].ref === ref)) return i
        return 0
    }

    title: "Model for " + targetName
    modal: true
    anchors.centerIn: Overlay.overlay
    width: dp(480)

    function openFor(i, name, k, r, px, ch, d, names) {
        row = i; targetName = name; modelNames = names
        kind = k; ref = r; pixsize = px; channels = ch; diam = d > 0 ? d : 0.5
        check()
        open()
    }
    // The source as Simulate would read it: green when usable, red with the reason.
    function check() {
        if (kind.length === 0) { usable = true; statusText = "not simulated"; return }
        var f = Julia.shell_sim_source(kind, kind === "ud" ? String(diam) : ref).split("\t")
        usable = f[0] === "ok"
        statusText = f.length > 1 ? f[1] : f[0]
        var px = f.length > 2 ? parseFloat(f[2]) : 0
        if (kind === "image" && usable && px > 0 && pixsize <= 0) pixsize = px
    }
    function choose(i) {
        var c = choices[i]
        if (c.kind === kind && c.ref === ref) return
        if (c.kind !== kind || c.kind === "session") ref = c.ref
        kind = c.kind
        pixsize = 0
        check()
    }

    ColumnLayout {
        width: parent.width
        spacing: dlg.dp(8)

        RowLayout {
            spacing: dlg.dp(8)
            Label { text: "model" }
            ComboBox {
                objectName: "simDlgKind"
                Layout.preferredWidth: dlg.dp(260)
                model: dlg.choices.map(function (c) { return c.label })
                currentIndex: dlg.choiceIndex()
                onActivated: function (i) { dlg.choose(i) }
            }
        }
        RowLayout {
            visible: dlg.kind === "modelfile" || dlg.kind === "image"
            spacing: dlg.dp(8)
            Button {
                objectName: "simDlgFile"
                text: "File…"
                onClicked: dlg.kind === "image" ? imagePicker.openAt(dlg.ref) : modelPicker.openAt(dlg.ref)
            }
            Label {
                Layout.fillWidth: true
                text: dlg.ref.length > 0 ? dlg.ref.replace(/^.*[\\/]/, "") : "no file chosen"
                color: dlg.ref.length > 0 ? "#000" : "#888"
                elide: Text.ElideMiddle
            }
        }
        RowLayout {
            visible: dlg.kind === "image"
            spacing: dlg.dp(8)
            Label { text: "pixel size" }
            TextField {
                objectName: "simDlgPix"
                Layout.preferredWidth: dlg.dp(80)
                text: dlg.pixsize > 0 ? dlg.pixsize.toFixed(4) : ""
                placeholderText: "header"
                onEditingFinished: { var v = parseFloat(text); dlg.pixsize = isNaN(v) || v <= 0 ? 0 : v }
            }
            Label { text: "mas"; color: "#666" }
            CheckBox {
                objectName: "simDlgChannels"
                text: "image's channels"
                checked: dlg.channels
                onToggled: dlg.channels = checked
            }
        }
        RowLayout {
            visible: dlg.kind === "ud"
            spacing: dlg.dp(8)
            Label { text: "diameter" }
            TextField {
                objectName: "simDlgDiam"
                Layout.preferredWidth: dlg.dp(80)
                text: String(dlg.diam)
                validator: DoubleValidator { bottom: 0.0 }
                onEditingFinished: { var v = parseFloat(text); if (!isNaN(v)) { dlg.diam = v; dlg.check() } }
            }
            Label { text: "mas"; color: "#666" }
        }
        Label {
            objectName: "simDlgStatus"
            Layout.fillWidth: true
            text: dlg.statusText
            color: dlg.kind.length === 0 ? "#888" : dlg.usable ? "#2e7d32" : "#b00000"
            elide: Text.ElideMiddle
        }
    }

    footer: DialogButtonBox {
        Button {
            objectName: "simDlgOk"
            text: "OK"
            DialogButtonBox.buttonRole: DialogButtonBox.AcceptRole
        }
        Button { text: "Cancel"; DialogButtonBox.buttonRole: DialogButtonBox.RejectRole }
    }
    onAccepted: chosen(row, kind, kind === "ud" ? "" : ref, kind === "image" ? pixsize : 0,
                       kind === "image" && channels, kind === "ud" ? diam : 0)

    FilePicker {
        id: imagePicker
        uiScale: dlg.uiScale; fontScale: dlg.fontScale; baseFontPt: dlg.baseFontPt
        title: "Image or cube for " + dlg.targetName
        filters: [{ label: "FITS files (*.fits *.fit)", patterns: "*.fits,*.fit" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) { dlg.ref = path; dlg.pixsize = 0; dlg.check() }
    }
    FilePicker {
        id: modelPicker
        uiScale: dlg.uiScale; fontScale: dlg.fontScale; baseFontPt: dlg.baseFontPt
        title: "Model file for " + dlg.targetName
        filters: [{ label: "Model files (*.toml)", patterns: "*.toml" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) { dlg.ref = path; dlg.check() }
    }
}
