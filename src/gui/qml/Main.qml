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

    // Never open larger than the screen: dp(1600) on a 2x panel is 3200 logical pixels, which
    // would put the window's own edges off the desktop before the user has done anything.
    //
    // 1600 rather than something squarer because every perspective is a settings column beside
    // a plot: Observing puts the Gantt next to the target list, Imaging the reconstruction next
    // to the regularisers. Width is what those layouts spend, and too little of it pushes the
    // plot into a strip.
    width:  Math.min(Screen.desktopAvailableWidth  * 0.95, dp(1600))
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

    function afterAction() { refreshConsole(); redrawPlot(); refreshObservables() }

    function setView() {
        win.status = Julia.shell_set_view(kindBox.currentText, colorBox.currentText,
                                          logyBox.enabled && logyBox.checked)
        win.afterAction()
    }

    // Which observables the loaded file actually holds. Julia counts VALID points, not tables,
    // and pushes the answer to every panel that gates a tick box on it. Without this the boxes
    // stay off after a load and the run buttons never enable — they have no other source of
    // truth, and defaulting them to on would offer observables the file does not contain.
    function refreshObservables() {
        var have = { v2: false, t3amp: false, t3phi: false,
                     cvis: false, flux: false, diffvis: false }
        var s = Julia.shell_observables()
        if (s.length > 0) {
            var parts = s.split(",")
            for (var i = 0; i < parts.length; ++i) {
                var kv = parts[i].split("=")
                if (kv.length === 2) have[kv[0]] = (kv[1] === "1")
            }
        }
        imageTab.haveV2      = have.v2
        imageTab.haveT3amp   = have.t3amp
        imageTab.haveT3phi   = have.t3phi
        imageTab.haveCvis    = have.cvis
        imageTab.haveFlux    = have.flux
        imageTab.haveDiffvis = have.diffvis
        imageTab.resetObservables()

        // Geometry the data suggests. A constant pixel size is not a neutral default: too
        // coarse and the image cannot represent what the data resolves, too fine and it is
        // unconstrained at the pixel scale.
        var g = Julia.shell_image_defaults()
        if (g.length > 0) {
            var gp = g.split(",")
            for (var k = 0; k < gp.length; ++k) {
                var kv2 = gp[k].split("=")
                if (kv2.length !== 2) continue
                if (kv2[0] === "pixsize") imageTab.pixsize = parseFloat(kv2[1])
                if (kv2[0] === "nx")      imageTab.nx      = parseInt(kv2[1])
            }
        }
        // The plan is keyed on the dataset as well as the geometry, so a new file needs one
        // even when nx and the pixel size happen to be unchanged.
        imageTab.refreshPlan()

        modelTab.haveV2      = have.v2
        modelTab.haveT3amp   = have.t3amp
        modelTab.haveT3phi   = have.t3phi
        modelTab.haveCvis    = have.cvis
        modelTab.haveFlux    = have.flux
        modelTab.haveDiffvis = have.diffvis
        modelTab.hasDataset  = s.length > 0
    }

    // Printed once at startup so the scale is diagnosable on a machine nobody can inspect
    // remotely. If the window looks wrong, this line says whether the screen asked for it
    // (autoScale) or an override did, and OITOOLSGUI_SCALE is how you argue with it.
    Component.onCompleted: console.log(
        "OITOOLSGUI ui scale: screen=" + (Screen.logicalPixelDensity * 25.4).toFixed(1) + " dpi"
        + "  devicePixelRatio=" + Screen.devicePixelRatio.toFixed(2)
        + "  autoScale=" + autoScale.toFixed(3)
        + (uiScaleOverride > 0 ? "  override=" + uiScaleOverride.toFixed(3) : "  (no override)")
        + "  -> uiScale=" + uiScale.toFixed(3) + " fontScale=" + fontScale.toFixed(3))

    // Reads the dataset list back from the session and selects the newest entry. Called both
    // after a load and at startup, because a Session may already hold datasets before the
    // window exists -- `oitoolsgui.jl file.oifits` loads them and then calls gui().
    function refreshDatasets() {
        var names = Julia.shell_dataset_names().split("\n").filter(function (s) { return s.length > 0 })
        datasetBox.model = names
        datasetBox.currentIndex = names.length - 1
    }

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
                // A dialog already on screen must not be opened a second time. Under the XDG
                // portal the dialog arrives over D-Bus and takes a visible moment, which is
                // precisely when a second click happens; the extra dialog then outlives the
                // one that was answered and reads as a picker that would not close.
                enabled: !openDialog.visible
                onClicked: openDialog.openAt(initialFolder)
            }
            Label { text: win.status; elide: Text.ElideMiddle; Layout.fillWidth: true }

            // ── task tray ────────────────────────────────────────────────────
            BusyIndicator { id: busy; running: false; implicitWidth: dp(20); implicitHeight: dp(20) }
            Label { id: taskLabel; text: "idle"; color: "#666" }
        }
    }

    // The pickers are drawn INSIDE this window. QtQuick.Dialogs.FileDialog leaves its own
    // window mapped on some systems — measured, see FilePicker.qml — and a window Qt will not
    // close cannot be closed from QML. A Popup has no window of its own to leave behind.
    FilePicker {
        id: openDialog
        uiScale: win.uiScale; fontScale: win.fontScale; baseFontPt: win.baseFontPt
        title: "Open OIFITS"
        filters: [{ label: "OIFITS files (*.oifits *.fits)", patterns: "*.oifits,*.fits" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) {
            busy.running = true; taskLabel.text = "loading…"
            // Still deferred: the load holds the GUI thread for seconds, and nothing repaints
            // while one of our handlers runs.
            loadTimer.path = path
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
            win.refreshDatasets()
            win.afterAction()
            busy.running = false; taskLabel.text = "idle"
            // The load holds the GUI thread for seconds, so by the time it returns the window
            // manager has had every chance to put focus somewhere else. Claim it back rather
            // than leaving the user to hunt for the window.
            win.raise(); win.requestActivate()
        }
    }

    FilePicker {
        id: saveDialog
        uiScale: win.uiScale; fontScale: win.fontScale; baseFontPt: win.baseFontPt
        title: "Export script"
        saveMode: true
        defaultSuffix: "jl"
        filters: [{ label: "Julia files (*.jl)", patterns: "*.jl" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) { exportTimer.path = path; exportTimer.restart() }
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
            //
            // `initialTab` overrides it, which is how an automated run reaches a perspective
            // it would otherwise never render: a tab that is never current is constructed but
            // never laid out, so its layout warnings never appear.
            currentIndex: initialTab
            TabButton { text: "Exploring" }
            TabButton { text: "Observing" }
            TabButton { text: "Modeling" }
            TabButton { text: "Imaging" }
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
                        Layout.preferredWidth: dp(130)
                        // `diffphi` is `visphi` against WAVELENGTH — OIFITS carries differential
                        // phase in OI_VIS with PHITYP=differential, so there is no separate
                        // field, only a more useful axis.
                        model: ["uv", "v2", "t3phi", "t3amp", "visamp", "visphi",
                                "diffphi", "diffvisamp", "flux"]
                        onActivated: win.setView()
                    }
                    Label { text: "Colour by:" }
                    ComboBox {
                        id: colorBox
                        model: ["baseline", "wav", "mjd", "none"]
                        onActivated: win.setView()
                    }
                    CheckBox {
                        id: logyBox
                        text: "log y"
                        // Only where it means something. A phase takes both signs, and a log
                        // axis on one drops half the points rather than rescaling them.
                        enabled: ["v2", "visamp", "t3amp", "flux", "diffvisamp"]
                                 .indexOf(kindBox.currentText) >= 0
                        onToggled: win.setView()
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
            // Sibling .qml files in this directory resolve as types without an import.
            // The three scale values are passed rather than reached for: a tab is a plain
            // component, so it cannot see `win`.
            ObserveTab {
                id: observeTab
                uiScale:    win.uiScale
                fontScale:  win.fontScale
                baseFontPt: win.baseFontPt
                onConsoleChanged: win.refreshConsole()
            }

            // ── Model ────────────────────────────────────────────────────────
            ModelTab {
                id: modelTab
                uiScale:    win.uiScale
                fontScale:  win.fontScale
                baseFontPt: win.baseFontPt
                onConsoleChanged: win.refreshConsole()
            }

            // ── Image ────────────────────────────────────────────────────────
            ImageTab {
                id: imageTab
                uiScale:    win.uiScale
                fontScale:  win.fontScale
                baseFontPt: win.baseFontPt
                onConsoleChanged: win.refreshConsole()
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
                        enabled: logToggle.checked && !saveDialog.visible
                        onClicked: saveDialog.openAt(initialFolder)
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
        onTriggered: {
            win.refreshDatasets()
            win.status = Julia.shell_ready()
            win.afterAction()
        }
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
