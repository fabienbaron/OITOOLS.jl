// OITOOLS — one window, four perspectives over one session.
//
// The context bar, task tray and command log sit OUTSIDE the tab stack on purpose: the
// dataset is what the perspectives share, background work must survive a tab switch, and the
// log records across modes. Putting any of them inside a tab would make this four programs
// sharing a title bar.
//
// QML holds no state. Everything is read from, or pushed through, Julia.

import QtQuick
import QtQuick.Window
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs
import jlqml
import Makie

ApplicationWindow {
    id: win

    // ── one scale factor for all the chrome ───────────────────────────────────
    //
    // Qt already scales `font.pointSize` by the font DPI, so text grows on a HiDPI screen by
    // itself. Pixel quantities do not, so without `dp()` the text outgrows the containers
    // holding it.
    //
    // autoScale is exactly the ratio Qt is applying to text, so at the default the fonts are
    // left completely alone (fontScale == 1) and only the pixels move. An explicit
    // OITOOLSGUI_SCALE is a departure from what the screen asked for, so it has to move both
    // -- hence fontScale, which is that departure and nothing else. Scaling pointSize by the
    // full uiScale would apply the screen's DPI twice.
    //
    // Detection is here rather than in Julia because Screen.logicalPixelDensity is live: it
    // is right before any window exists, and it follows the window to a monitor with a
    // different scale. Julia supplies only uiScaleOverride (0 = decide for yourself).
    readonly property real autoScale:
        Math.max(0.5, Math.min(4.0, Screen.logicalPixelDensity * 25.4 / 96.0))
    // 1.25 by default, not 1.0. Measured on the development display: at autoScale == 1.0 the
    // window is legibly too small, and every screen we have tried reports 96 dpi and a device
    // pixel ratio of 1, so autoScale alone never grows anything. OITOOLSGUI_SCALE overrides it
    // outright -- set it to 1.0 to get the old size back.
    readonly property real defaultDensity: 1.25
    readonly property real uiScale:
        uiScaleOverride > 0 ? uiScaleOverride : autoScale * defaultDensity
    readonly property real fontScale: uiScale / autoScale

    function dp(px)     { return Math.round(px * uiScale) }   // pixel lengths
    function pt(points) { return points * fontScale }         // explicit point sizes

    // Base font for the whole window. Qt's default is whatever the platform theme says, which
    // under WSL is 9 pt at 96 dpi -- small on a large panel. Everything inherits this unless
    // it sets its own size, so raising it here raises the whole UI, and OITOOLSGUI_SCALE
    // multiplies it along with the spacing.
    property real baseFontPt: 11
    font.pointSize: pt(baseFontPt)

    // Never open larger than the screen: dp(1280) on a 2x panel is 2560 logical pixels, which
    // would put the window's own edges off the desktop before the user has done anything.
    width:  Math.min(Screen.desktopAvailableWidth  * 0.95, dp(1280))
    // The plot is one of five stacked regions -- context bar, tabs, plot, pick line, console --
    // so the window needs the height to leave it more than a band.
    height: Math.min(Screen.desktopAvailableHeight * 0.95, dp(1040))
    visible: true
    title: "OITOOLS"
    // Test hook: fullscreen makes screen coordinates map 1:1 onto client coordinates inside a
    // nested compositor, so automated clicking needs no per-window offset. Off unless asked.
    visibility: fullscreenOnStart ? Window.FullScreen : Window.Windowed

    property string status: initialStatus
    property string pickText: ""

    // Every action funnels through here so the console cannot go stale: the pane is the
    // record of what happened, and a record that updates only sometimes is worse than none.
    function refreshConsole() { if (logToggle.checked) logArea.text = Julia.shell_console() }

    // Repaint the plot. This is NOT optional and it is easy to forget: the shell now updates
    // Makie Observables rather than inserting plots, and an Observable assignment does not by
    // itself ask Qt for a new frame. The old code repainted only as a side effect of creating
    // plots -- the same side effect that was corrupting the GL context. MakieArea's own mouse
    // handlers call update() for exactly this reason.
    function redrawPlot() { if (typeof makieArea !== "undefined") makieArea.update() }

    function afterAction() { refreshConsole(); redrawPlot() }

    // Printed once at startup so the scale is diagnosable on a machine nobody can inspect
    // remotely. If the window looks wrong, this line says whether the screen asked for it
    // (autoScale) or an override did, and OITOOLSGUI_SCALE is how you argue with it.
    Component.onCompleted: console.log(
        "OITOOLSGUI ui scale: screen=" + (Screen.logicalPixelDensity * 25.4).toFixed(1) + " dpi"
        + "  devicePixelRatio=" + Screen.devicePixelRatio.toFixed(2)
        + "  autoScale=" + autoScale.toFixed(3)
        + (uiScaleOverride > 0 ? "  override=" + uiScaleOverride.toFixed(3) : "  (no override)")
        + "  -> uiScale=" + uiScale.toFixed(3) + " fontScale=" + fontScale.toFixed(3))

    // ── context bar: what every perspective is looking at ─────────────────────
    header: ToolBar {
        RowLayout {
            anchors.fill: parent
            anchors.leftMargin: dp(8); anchors.rightMargin: dp(8)
            spacing: dp(10)

            Label { text: "Dataset:"; font.bold: true }
            ComboBox {
                id: datasetBox
                Layout.preferredWidth: dp(320)
                model: []
                enabled: model.length > 0
                onActivated: { win.status = Julia.shell_select_dataset(currentIndex + 1); win.afterAction() }
            }
            Button {
                text: "Open OIFITS…"
                onClicked: openDialog.open()
            }
            Label { text: win.status; elide: Text.ElideMiddle; Layout.fillWidth: true }

            // ── task tray ────────────────────────────────────────────────────
            BusyIndicator { id: busy; running: false; implicitWidth: dp(20); implicitHeight: dp(20) }
            Label { id: taskLabel; text: "idle"; color: "#666" }
        }
    }

    FileDialog {
        id: openDialog
        title: "Open OIFITS"
        // Start somewhere with data in it rather than wherever Qt defaults to. Julia works out
        // where: beside the last dataset opened, else the bundled data/ directory.
        currentFolder: initialFolder
        nameFilters: ["OIFITS files (*.oifits *.fits)", "All files (*)"]

        // Close, then work — with enough delay for the close to actually happen.
        //
        // The load takes seconds and runs synchronously on the GUI thread, and Qt cannot
        // finish tearing the dialog down while one of our handlers is running, so the dialog
        // stays on screen for the whole load unless the work is deferred.
        //
        // Two tempting alternatives do not work here: `options: DontUseNativeDialog` stops the
        // dialog accepting a selection, and `onVisibleChanged` never fires. Raise the interval
        // if a platform needs longer.
        onAccepted: {
            var picked = selectedFile.toString()
            openDialog.close()
            busy.running = true; taskLabel.text = "loading…"
            loadTimer.path = picked
            loadTimer.restart()
        }
    }

    // Poll the last identified point. Julia produces it from Makie's event stream, so there
    // is no QML signal to connect to; 150 ms is imperceptible for a label and costs nothing.
    Timer {
        interval: 150; running: true; repeat: true
        onTriggered: {
            var t = Julia.shell_pick_text()
            if (t !== win.pickText) win.pickText = t
        }
    }

    Timer {
        id: loadTimer
        property string path: ""
        interval: 250; repeat: false
        onTriggered: {
            win.status = Julia.shell_open(path)
            datasetBox.model = Julia.shell_dataset_names().split("\n").filter(function (s) { return s.length > 0 })
            datasetBox.currentIndex = datasetBox.model.length - 1
            win.afterAction()
            busy.running = false; taskLabel.text = "idle"
        }
    }

    FileDialog {
        id: saveDialog
        title: "Export script"
        fileMode: FileDialog.SaveFile
        nameFilters: ["Julia files (*.jl)"]
        // Same treatment as openDialog: writing the script is fast, but the dialog still
        // cannot finish closing while this handler runs, and the symptom is identical.
        onAccepted: {
            var target = selectedFile.toString()
            saveDialog.close()
            exportTimer.path = target
            exportTimer.restart()
        }
    }

    Timer {
        id: exportTimer
        property string path: ""
        interval: 250; repeat: false
        onTriggered: { win.status = Julia.shell_export(path); win.refreshConsole() }
    }

    ColumnLayout {
        anchors.fill: parent
        spacing: 0

        TabBar {
            id: tabs
            Layout.fillWidth: true
            // Explore first, and so the default. It is the only perspective that does
            // anything yet, and opening on a placeholder makes the application look broken
            // before it has done a thing.
            currentIndex: 0
            TabButton { text: "Explore" }
            TabButton { text: "Observe" }
            TabButton { text: "Model" }
            TabButton { text: "Image" }
        }

        StackLayout {
            id: stack
            Layout.fillWidth: true
            Layout.fillHeight: true
            currentIndex: tabs.currentIndex

            // ── Explore ──────────────────────────────────────────────────────
            ColumnLayout {
                spacing: dp(6)
                RowLayout {
                    Layout.fillWidth: true
                    Layout.margins: dp(6)
                    spacing: dp(8)
                    Label { text: "Plot:" }
                    ComboBox {
                        id: kindBox
                        model: ["uv", "v2", "t3phi", "t3amp"]
                        onActivated: { win.status = Julia.shell_set_view(currentText, colorBox.currentText); win.afterAction() }
                    }
                    Label { text: "Colour by:" }
                    ComboBox {
                        id: colorBox
                        model: ["baseline", "wav", "mjd", "none"]
                        onActivated: { win.status = Julia.shell_set_view(kindBox.currentText, currentText); win.afterAction() }
                    }
                    Item { Layout.fillWidth: true }
                    Label {
                        text: "scroll = zoom · drag = pan · click = identify · right click = reset view"
                        color: "#888"; font.pointSize: pt(9)
                    }
                }

                MakieArea {
                    id: makieArea
                    Layout.fillWidth: true
                    Layout.fillHeight: true
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

                Label {
                    Layout.fillWidth: true
                    Layout.margins: dp(6)
                    text: win.pickText.length > 0 ? win.pickText : "click a point to identify it"
                    color: win.pickText.length > 0 ? "#000" : "#888"
                    font.family: "monospace"
                    elide: Text.ElideRight
                }
            }

            // ── Observe ──────────────────────────────────────────────────────
            Item {
                Label {
                    anchors.centerIn: parent
                    horizontalAlignment: Text.AlignHCenter
                    color: "#666"
                    text: "Observe — target resolution, array and instrument, observability.\n" +
                          "Blocked on the twilight and RA-convention fixes; see the plan."
                }
            }

            // ── Model ────────────────────────────────────────────────────────
            Item {
                Label {
                    anchors.centerIn: parent
                    horizontalAlignment: Text.AlignHCenter
                    color: "#666"
                    text: "Model — component tree, parameter table (fixed / free / expression),\n" +
                          "sky canvas, optimiser and uncertainty selectors."
                }
            }

            // ── Image ────────────────────────────────────────────────────────
            Item {
                Label {
                    anchors.centerIn: parent
                    horizontalAlignment: Text.AlignHCenter
                    color: "#666"
                    text: "Image — engine (VMLMB / MaxEnt / ADMM / SPARCO / SQUEEZE / VI),\n" +
                          "engine-dependent regularisation, live reconstruction."
                }
            }
        }

        // ── command log: spans every perspective ──────────────────────────────
        Rectangle {
            Layout.fillWidth: true
            Layout.preferredHeight: logToggle.checked ? dp(200) : dp(32)
            color: "#f4f4f4"
            border.color: "#ddd"

            ColumnLayout {
                anchors.fill: parent
                spacing: 0
                RowLayout {
                    Layout.fillWidth: true
                    Layout.margins: dp(4)
                    CheckBox {
                        id: logToggle
                        text: "Console"
                        checked: true
                        onToggled: if (checked) win.afterAction()
                    }
                    Label {
                        text: "·  > lines are the script; ! lines are errors"
                        color: "#888"; font.pointSize: pt(baseFontPt - 2)
                        visible: logToggle.checked
                    }
                    Item { Layout.fillWidth: true }
                    Button {
                        text: "Copy script"
                        enabled: logToggle.checked
                        onClicked: { logArea.text = Julia.shell_script(); logArea.selectAll(); logArea.copy() }
                    }
                    Button {
                        text: "Export script…"
                        enabled: logToggle.checked
                        onClicked: saveDialog.open()
                    }
                }
                ScrollView {
                    id: logScroll
                    visible: logToggle.checked
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    TextArea {
                        id: logArea
                        readOnly: true
                        font.family: "monospace"
                        font.pointSize: pt(baseFontPt - 2)
                        text: ""
                        // Follow the tail. Without this the newest line -- which is the error
                        // you just caused -- is the one line you cannot see.
                        onTextChanged: cursorPosition = length
                    }
                }
            }
        }
    }

    // Automation hooks, both inert unless the caller asks for them.
    //
    // The ready timer fires once the window is up, so anything it does runs on the GUI thread
    // against a live GL screen -- the state a headless unit test cannot reach. gui() installs
    // the Julia side; with no hook installed this is a no-op call.
    Timer {
        interval: 700
        running: true
        repeat: false
        onTriggered: { win.status = Julia.shell_ready(); win.afterAction() }
    }

    // Closes the window after a delay so the shell can be exercised under a virtual display
    // with nobody to close it.
    Timer {
        interval: autoQuitMs > 0 ? autoQuitMs : 1000
        running: autoQuitMs > 0
        repeat: false
        onTriggered: Qt.quit()
    }
}
