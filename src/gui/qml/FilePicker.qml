// The file picker, drawn inside this window.
//
// `QtQuick.Dialogs.FileDialog` is not used, and the reason is measured rather than assumed:
// `test/gui/filepicker_min.jl` is that dialog and nothing else, and on some systems it reports
//
//     signal: visibleChanged -> false
//     signal: accepted; visible is now false
//
// while its window is still on screen. The dialog is a separate top-level window; where Qt
// cannot get a usable GL surface for it — X11 forwarding, a failing EGL, a software fallback —
// hiding it does not unmap it, and destroying the QML object does not either. Nothing in QML
// can close a window Qt will not close.
//
// A Popup is drawn on the parent window's own scene graph. No second window is created, so
// none can be left behind, and the failure mode is gone by construction rather than by
// workaround. The cost is that the listing, the places and the filtering come from Julia
// instead of from the toolkit — `src/gui/filepicker.jl` is that side.
//
// Usage:
//
//     FilePicker {
//         id: picker
//         title: "Open OIFITS"
//         filters: [{ label: "OIFITS files", patterns: "*.oifits,*.fits" },
//                   { label: "All files",    patterns: "*" }]
//         onAccepted: (path) => doSomethingWith(path)
//     }
//     ...
//     picker.openAt("")        // "" means "wherever Julia thinks is sensible"

import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml            // supplies `Julia`, which every listing call below goes through

Popup {
    id: root

    // ── contract ─────────────────────────────────────────────────────────────
    property real uiScale:    1.0
    property real fontScale:  1.0
    property real baseFontPt: 11
    function dp(px)     { return Math.round(px * uiScale) }
    function pt(points) { return points * fontScale }

    property string title: "Open a file"
    property bool   saveMode: false            // show a filename field and confirm overwrites
    property string defaultSuffix: ""
    // Each entry is { label, patterns } with patterns a comma-separated glob list.
    property var    filters: [{ label: "All files", patterns: "*" }]
    property string folder: ""
    property string selectedName: ""
    property bool   showHidden: false

    signal accepted(string path)
    signal rejected()

    modal: true
    focus: true
    closePolicy: Popup.CloseOnEscape
    anchors.centerIn: Overlay.overlay
    width:  Math.min(parent ? parent.width  - dp(80) : dp(900), dp(900))
    height: Math.min(parent ? parent.height - dp(80) : dp(560), dp(560))
    padding: dp(10)

    readonly property string currentPatterns:
        (filters && filters.length > 0 && filterBox.currentIndex >= 0
         && filterBox.currentIndex < filters.length)
            ? filters[filterBox.currentIndex].patterns : "*"

    // ── opening ──────────────────────────────────────────────────────────────
    function openAt(hint) {
        folder = Julia.picker_start(hint === undefined ? "" : hint)
        selectedName = ""
        nameField.text = ""
        refresh()
        open()
    }

    function refresh() {
        entryModel.clear()
        var rows = Julia.picker_list(folder, currentPatterns, showHidden ? true : false)
        if (rows.length > 0) {
            var lines = rows.split("\n")
            for (var i = 0; i < lines.length; ++i) {
                var f = lines[i].split("\t")
                if (f.length < 5) continue
                entryModel.append({ kind: f[0], name: f[1], size: f[2],
                                    modified: f[3], sortkey: f[4] })
            }
        }
        fileList.currentIndex = entryModel.count > 0 ? 0 : -1
        statusLabel.text = entryModel.count + (entryModel.count === 1 ? " item" : " items")
    }

    onFolderChanged: if (visible) refresh()
    onOpened: fileList.forceActiveFocus()
    onClosed: root.rejected()

    ListModel { id: entryModel }

    // The path as clickable segments. A long path in one label is unreadable and, worse,
    // gives no way back up except a button pressed repeatedly.
    function segments() {
        var p = folder.replace(/\/+$/, "")
        if (p.length === 0) return [{ label: "/", path: "/" }]
        var parts = p.split("/")
        var out = [{ label: "/", path: "/" }]
        var acc = ""
        for (var i = 1; i < parts.length; ++i) {
            if (parts[i].length === 0) continue
            acc = acc + "/" + parts[i]
            out.push({ label: parts[i], path: acc })
        }
        return out
    }

    function activate(i) {
        if (i < 0 || i >= entryModel.count) return
        var e = entryModel.get(i)
        if (e.kind === "dir") {
            folder = Julia.picker_join(folder, e.name)
        } else if (e.kind === "file") {
            selectedName = e.name
            nameField.text = e.name
            root.confirm()
        }
    }

    function confirm() {
        var name = saveMode ? nameField.text.trim() : selectedName
        if (name.length === 0) { statusLabel.text = "nothing selected"; return }
        if (saveMode && defaultSuffix.length > 0 && name.indexOf(".") < 0)
            name = name + "." + defaultSuffix
        var full = Julia.picker_join(folder, name)

        // A typed name that turns out to be a directory should walk into it, not be "opened".
        if (Julia.picker_kind(full) === "dir") { folder = full; nameField.text = ""; return }

        if (saveMode && Julia.picker_would_overwrite(full) === "1") {
            overwriteDialog.target = full
            overwriteDialog.open()
            return
        }
        root.finish(full)
    }

    // `close()` fires onClosed, which reports a rejection; suppress that for a real accept.
    property bool _accepting: false
    function finish(path) {
        _accepting = true
        close()
        _accepting = false
        root.accepted(path)
    }

    // ── layout ───────────────────────────────────────────────────────────────
    contentItem: ColumnLayout {
        spacing: root.dp(8)

        RowLayout {
            Layout.fillWidth: true
            Label { text: root.title; font.bold: true
                    font.pointSize: root.pt(root.baseFontPt + 1) }
            Item { Layout.fillWidth: true }
            ToolButton {
                text: "↑"
                ToolTip.visible: hovered; ToolTip.text: "parent folder"
                implicitHeight: root.dp(26); implicitWidth: root.dp(30)
                onClicked: root.folder = Julia.picker_parent(root.folder)
            }
        }

        // breadcrumb
        Flickable {
            Layout.fillWidth: true
            implicitHeight: root.dp(26)
            contentWidth: crumbRow.width
            flickableDirection: Flickable.HorizontalFlick
            clip: true
            Row {
                id: crumbRow
                spacing: 0
                Repeater {
                    model: root.segments()
                    delegate: Row {
                        spacing: 0
                        Label {
                            text: "›"; color: "#999"
                            visible: index > 0
                            padding: root.dp(3)
                            anchors.verticalCenter: parent.verticalCenter
                        }
                        Button {
                            flat: true
                            text: modelData.label
                            implicitHeight: root.dp(24)
                            font.pointSize: root.pt(root.baseFontPt - 1)
                            onClicked: root.folder = modelData.path
                        }
                    }
                }
            }
        }

        RowLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: root.dp(8)

            // places
            Rectangle {
                Layout.preferredWidth: root.dp(170)
                Layout.fillHeight: true
                color: "#f6f6f6"
                border.color: "#dcdcdc"
                radius: root.dp(3)
                ListView {
                    id: placeList
                    anchors.fill: parent
                    anchors.margins: root.dp(2)
                    clip: true
                    model: placeModel
                    delegate: ItemDelegate {
                        width: placeList.width
                        implicitHeight: root.dp(26)
                        text: model.label
                        font.pointSize: root.pt(root.baseFontPt - 1)
                        onClicked: root.folder = model.path
                    }
                }
            }

            // listing
            Rectangle {
                Layout.fillWidth: true
                Layout.fillHeight: true
                color: "white"
                border.color: "#dcdcdc"
                radius: root.dp(3)

                ColumnLayout {
                    anchors.fill: parent
                    anchors.margins: root.dp(1)
                    spacing: 0

                    Rectangle {
                        Layout.fillWidth: true
                        implicitHeight: root.dp(22)
                        color: "#f0f0f0"
                        RowLayout {
                            anchors.fill: parent
                            anchors.leftMargin: root.dp(8)
                            anchors.rightMargin: root.dp(8)
                            Label { text: "Name"; Layout.fillWidth: true; color: "#555"
                                    font.pointSize: root.pt(root.baseFontPt - 2) }
                            Label { text: "Size"; Layout.preferredWidth: root.dp(80)
                                    horizontalAlignment: Text.AlignRight; color: "#555"
                                    font.pointSize: root.pt(root.baseFontPt - 2) }
                            Label { text: "Modified"; Layout.preferredWidth: root.dp(120)
                                    horizontalAlignment: Text.AlignRight; color: "#555"
                                    font.pointSize: root.pt(root.baseFontPt - 2) }
                        }
                    }

                    ListView {
                        id: fileList
                        Layout.fillWidth: true
                        Layout.fillHeight: true
                        clip: true
                        model: entryModel
                        focus: true
                        keyNavigationEnabled: true
                        ScrollBar.vertical: ScrollBar {}
                        highlightMoveDuration: 0

                        // Backspace goes up, Return opens, Escape cancels — the three gestures
                        // every file list has had for thirty years.
                        Keys.onReturnPressed: root.activate(currentIndex)
                        Keys.onEnterPressed:  root.activate(currentIndex)
                        Keys.onPressed: function (event) {
                            if (event.key === Qt.Key_Backspace) {
                                root.folder = Julia.picker_parent(root.folder)
                                event.accepted = true
                            }
                        }
                        onCurrentIndexChanged: {
                            if (currentIndex < 0 || currentIndex >= entryModel.count) return
                            var e = entryModel.get(currentIndex)
                            if (e.kind === "file") {
                                root.selectedName = e.name
                                if (!root.saveMode || nameField.text.length === 0)
                                    nameField.text = e.name
                            }
                        }

                        delegate: ItemDelegate {
                            id: rowItem
                            width: fileList.width
                            implicitHeight: root.dp(24)
                            highlighted: ListView.isCurrentItem
                            // Captured, because `model` inside a nested Control means that
                            // control's own model rather than this row.
                            readonly property string rKind: model.kind
                            readonly property string rName: model.name
                            readonly property string rSize: model.size
                            readonly property string rMod:  model.modified
                            readonly property int    rIndex: index

                            onClicked: fileList.currentIndex = rIndex
                            onDoubleClicked: root.activate(rIndex)

                            contentItem: RowLayout {
                                spacing: root.dp(6)
                                // Geometric glyphs, not emoji: the emoji faces are not in
                                // the fonts Qt falls back to here and render as nothing at
                                // all, which is worse than no icon.
                                Label {
                                    text: rowItem.rKind === "dir" ? "▸" : "·"
                                    color: rowItem.rKind === "dir" ? "#1a4d8f" : "#aaa"
                                    Layout.preferredWidth: root.dp(10)
                                    horizontalAlignment: Text.AlignHCenter
                                    font.pointSize: root.pt(root.baseFontPt - 1)
                                }
                                Label {
                                    text: rowItem.rName
                                    Layout.fillWidth: true
                                    elide: Text.ElideMiddle
                                    color: rowItem.rKind === "dir" ? "#1a4d8f" : "#222"
                                    font.bold: rowItem.rKind === "dir"
                                }
                                Label {
                                    text: rowItem.rSize
                                    Layout.preferredWidth: root.dp(80)
                                    horizontalAlignment: Text.AlignRight
                                    color: "#666"
                                    font.pointSize: root.pt(root.baseFontPt - 2)
                                }
                                Label {
                                    text: rowItem.rMod
                                    Layout.preferredWidth: root.dp(120)
                                    horizontalAlignment: Text.AlignRight
                                    color: "#666"
                                    font.pointSize: root.pt(root.baseFontPt - 2)
                                }
                            }
                        }
                    }
                }
            }
        }

        RowLayout {
            Layout.fillWidth: true
            spacing: root.dp(6)
            Label { text: root.saveMode ? "Save as" : "File"; color: "#555" }
            TextField {
                id: nameField
                Layout.fillWidth: true
                implicitHeight: root.dp(28)
                readOnly: !root.saveMode
                placeholderText: root.saveMode ? "filename" : ""
                onAccepted: root.confirm()
            }
            ComboBox {
                id: filterBox
                Layout.preferredWidth: root.dp(200)
                implicitHeight: root.dp(28)
                textRole: "label"
                model: root.filters
                onActivated: root.refresh()
            }
        }

        RowLayout {
            Layout.fillWidth: true
            spacing: root.dp(6)
            CheckBox {
                text: "hidden files"
                checked: root.showHidden
                font.pointSize: root.pt(root.baseFontPt - 1)
                onToggled: { root.showHidden = checked; root.refresh() }
            }
            Label { id: statusLabel; color: "#666"
                    font.pointSize: root.pt(root.baseFontPt - 1) }
            Item { Layout.fillWidth: true }
            Button { text: "Cancel"; implicitHeight: root.dp(28); onClicked: root.close() }
            Button {
                text: root.saveMode ? "Save" : "Open"
                implicitHeight: root.dp(28)
                highlighted: true
                enabled: root.saveMode ? nameField.text.trim().length > 0
                                       : fileList.currentIndex >= 0
                onClicked: root.activate(fileList.currentIndex)
            }
        }
    }

    ListModel { id: placeModel }

    Component.onCompleted: {
        var rows = Julia.picker_places()
        if (rows.length === 0) return
        var lines = rows.split("\n")
        for (var i = 0; i < lines.length; ++i) {
            var f = lines[i].split("\t")
            if (f.length === 2) placeModel.append({ label: f[0], path: f[1] })
        }
    }

    Dialog {
        id: overwriteDialog
        property string target: ""
        title: "Replace the file?"
        modal: true
        anchors.centerIn: Overlay.overlay
        standardButtons: Dialog.Yes | Dialog.No
        Label {
            text: overwriteDialog.target + "\nalready exists."
            wrapMode: Text.WordWrap
        }
        onAccepted: root.finish(overwriteDialog.target)
    }
}
