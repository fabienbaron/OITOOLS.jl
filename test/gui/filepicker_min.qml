// The file picker, and nothing else.
//
// No Makie, no GL, no session, no OITOOLS: just a window, a button and a FileDialog. If the
// picker misbehaves here, the cause is QtQuick.Dialogs or QML.jl and nothing this package
// does; if it behaves here but not in the GUI, the cause is something the GUI adds. That
// split is the whole point, and it is not answerable from inside the full application.
//
// Every lifecycle signal is logged through Julia so the transcript survives the window.

import QtQuick
import QtQuick.Window
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs
import jlqml

ApplicationWindow {
    id: win
    visible: true
    width: 700; height: 420
    title: "file picker only"

    // `mode` selects which of the candidate handlers runs, so the three can be compared
    // without editing anything:
    //   none    -- let the dialog manage its own lifecycle, touch nothing
    //   close   -- call close() inside onAccepted, which is what the GUI does today
    //   timer   -- close(), then defer the "work" by 250 ms, also what the GUI does today
    //   destroy -- build the dialog at click time and destroy() the OBJECT afterwards. If the
    //              window outlives `visible = false` because the QML object still owns it,
    //              destroying the object is what takes the window with it.
    //   inline  -- no FileDialog at all: an in-window Popup listing the directory. There is no
    //              second window, so there is nothing that can be left behind. This is the
    //              candidate that cannot fail the observed way, at the cost of being our own
    //              picker rather than the platform's.
    property string mode: pickerMode

    function log(s) { Julia.picker_log(s) }

    Component.onCompleted: {
        log("window up, mode=" + mode)
        log("dialog visible at start: " + openDialog.visible)
    }

    ColumnLayout {
        anchors.fill: parent
        anchors.margins: 12
        spacing: 8

        Label { text: "mode: " + win.mode; font.bold: true }
        Label { id: stateLabel; text: "dialog.visible = " + openDialog.visible }
        Label { id: pickedLabel; text: "picked: (nothing yet)"; elide: Text.ElideMiddle
                Layout.fillWidth: true }

        Button {
            id: openButton
            text: "Open a file…"
            onClicked: {
                win.log("button clicked; mode=" + win.mode)
                if (win.mode === "inline") {
                    inlinePicker.folder = Julia.picker_home()
                    inlinePicker.open()
                } else if (win.mode === "destroy") {
                    win.spawnDialog()
                } else {
                    openDialog.open()
                }
            }
        }

        Item { Layout.fillHeight: true }
    }

    FileDialog {
        id: openDialog
        title: "Pick a file"
        currentFolder: initialFolder
        nameFilters: ["All files (*)"]

        // The claim in Main.qml is that this never fires. Test it rather than repeat it.
        onVisibleChanged: win.log("  signal: visibleChanged -> " + visible)

        onAccepted: {
            win.log("  signal: accepted; visible is now " + visible)
            var picked = selectedFile.toString()
            pickedLabel.text = "picked: " + picked

            if (win.mode === "close") {
                win.log("  calling close()")
                openDialog.close()
                win.log("  after close(), visible = " + visible)
            } else if (win.mode === "timer") {
                win.log("  calling close()")
                openDialog.close()
                win.log("  after close(), visible = " + visible)
                workTimer.restart()
            } else {
                win.log("  mode=none: not touching the dialog")
            }
            checkTimer.restart()
        }

        onRejected: win.log("  signal: rejected; visible is now " + visible)
    }

    // ── candidate: build it, use it, destroy it ──────────────────────────────
    //
    // A QML object that still exists still owns whatever window it created. If hiding the
    // dialog leaves its window mapped, deleting the object is the next thing to try, because
    // it takes the window down through a different path than `visible = false`.
    property var spawned: null

    function spawnDialog() {
        if (spawned !== null) { win.log("  a spawned dialog already exists"); return }
        spawned = dialogFactory.createObject(win)
        if (spawned === null) { win.log("  ! could not create the dialog"); return }
        win.log("  created dialog object; calling open()")
        spawned.open()
    }

    Component {
        id: dialogFactory
        FileDialog {
            title: "Pick a file"
            parentWindow: win               // imperative creation gives it no parent otherwise
            currentFolder: initialFolder
            nameFilters: ["All files (*)"]
            onVisibleChanged: win.log("  signal: visibleChanged -> " + visible)
            onAccepted: {
                win.log("  signal: accepted; visible is now " + visible)
                pickedLabel.text = "picked: " + selectedFile.toString()
                win.log("  destroying the dialog object")
                win.spawned.destroy()
                win.spawned = null
                checkTimer.restart()
            }
            onRejected: {
                win.log("  signal: rejected")
                win.spawned.destroy(); win.spawned = null
            }
        }
    }

    // ── candidate: our own picker, inside this window ────────────────────────
    //
    // A Popup is drawn on the parent window's own scene graph, so no second window is created
    // and none can survive. Julia lists the directory, which is the one thing QML cannot do.
    Popup {
        id: inlinePicker
        anchors.centerIn: Overlay.overlay
        width: Math.min(win.width - 60, 560)
        height: Math.min(win.height - 60, 340)
        modal: true
        property string folder: ""

        onFolderChanged: {
            entryModel.clear()
            var rows = Julia.picker_listdir(folder)
            if (rows.length === 0) return
            var lines = rows.split("\n")
            for (var i = 0; i < lines.length; ++i) {
                var kv = lines[i].split("\t")
                if (kv.length === 2) entryModel.append({ kind: kv[0], name: kv[1] })
            }
        }
        onOpened: { entryList.forceActiveFocus(); win.log("  inline picker opened at " + folder) }
        onClosed:  win.log("  inline picker closed; visible = " + visible)

        ListModel { id: entryModel }

        ColumnLayout {
            anchors.fill: parent
            spacing: 6
            Label { text: inlinePicker.folder; elide: Text.ElideMiddle; Layout.fillWidth: true }
            ListView {
                id: entryList
                Layout.fillWidth: true
                Layout.fillHeight: true
                clip: true
                model: entryModel
                focus: true
                keyNavigationEnabled: true
                currentIndex: 0
                highlight: Rectangle { color: "#cfe3ff" }
                ScrollBar.vertical: ScrollBar {}

                function activate(i) {
                    if (i < 0 || i >= entryModel.count) return
                    var e = entryModel.get(i)
                    if (e.kind === "file") {
                        pickedLabel.text = "picked: " +
                            Julia.picker_join(inlinePicker.folder, e.name)
                        win.log("  inline pick: " + e.name)
                        inlinePicker.close()
                        checkTimer.restart()
                    } else {
                        inlinePicker.folder = Julia.picker_join(inlinePicker.folder, e.name)
                        entryList.currentIndex = 0
                    }
                }
                Keys.onReturnPressed: activate(currentIndex)
                Keys.onEnterPressed:  activate(currentIndex)

                delegate: ItemDelegate {
                    width: entryList.width
                    text: (model.kind === "file" ? "  " : "▸ ") + model.name
                    highlighted: ListView.isCurrentItem
                    onClicked: { entryList.currentIndex = index; entryList.activate(index) }
                }
            }
            Button { text: "Cancel"; onClicked: inlinePicker.close() }
        }
    }

    // Work standing in for the seconds-long load the real GUI does here, so the question
    // "does blocking the GUI thread keep the dialog on screen" gets an answer.
    Timer {
        id: workTimer
        interval: 250; repeat: false
        onTriggered: {
            win.log("  deferred work starting (blocking 2 s)")
            Julia.picker_block(2)
            win.log("  deferred work done; dialog visible = " + openDialog.visible)
        }
    }

    // The state some time after everything settled: this is what the user actually sees.
    Timer {
        id: checkTimer
        interval: 4000; repeat: false
        onTriggered: {
            var v = win.mode === "inline"  ? inlinePicker.visible
                  : win.mode === "destroy" ? (win.spawned === null ? false : win.spawned.visible)
                                           : openDialog.visible
            win.log("VERDICT after 4 s: picker visible (per QML) = " + v)
            win.log("        >>> now LOOK: is the picker still on screen? <<<")
            stateLabel.text = "picker visible (per QML) = " + v
        }
    }
}
