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
    // qtTextScale is exactly the ratio Qt is applying to text by itself, so fontScale carries
    // only the DEPARTURE from it and the pointSize never gets the screen's DPI twice.
    //
    // It reads logicalPixelDensity, which is 96 dpi on essentially every Linux screen -- that
    // is the point: it says what Qt did, not what the screen is.
    readonly property real qtTextScale:
        Math.max(0.5, Math.min(4.0, Screen.logicalPixelDensity * 25.4 / 96.0))

    // What the screen actually is, from PHYSICAL pixel density.
    //
    // `Screen.pixelDensity` comes from the EDID and genuinely differs between machines, where
    // `logicalPixelDensity` does not. Measured: a 24" 1920x1080 desktop reports 92.6 dpi and a
    // 16" laptop panel around 189, while BOTH report logicalPixelDensity 96 and
    // devicePixelRatio 1. Scaling off the logical figure is why one hardcoded factor could fit
    // one machine and be far too large on the other.
    readonly property real physicalDpi: Screen.pixelDensity * 25.4

    // Anchored on two points judged by eye, not derived: 0.875 on this 92.6 dpi desktop and
    // 1.25 on a laptop panel of roughly 189 dpi. `refDpi`/`refScale`/`dpiExponent` ARE that
    // fit -- move them if a screen reads wrong, rather than trying to re-derive a rule.
    //
    // The exponent falls out as 0.5. Matching physical SIZE would be 1.0 and give 0.61 here,
    // too small: a desktop monitor sits further away than a laptop screen and wants larger
    // elements to subtend the same angle. Viewing distance is not knowable, so the root splits
    // the difference, and it happens to pass through both judgements.
    //
    // This governs the CHROME only. The window is sized from the screen (below) and plot text
    // from the density (`live_plot_scale` in Julia) -- tying either to this factor was tried
    // and both were wrong: the window shrank when the widgets did, leaving less room for the
    // data rather than more.
    //
    // 189 is an estimate for the laptop, not a reading. The startup line below prints the
    // physical dpi, so it can be replaced with the real one and the exponent refitted.
    readonly property real refDpi:      92.6
    readonly property real refScale:    0.875
    readonly property real dpiExponent: 0.5
    readonly property real autoScale:
        Math.max(0.5, Math.min(4.0,
            physicalDpi > 0 ? refScale * Math.pow(physicalDpi / refDpi, dpiExponent)
                            : refScale))

    // Turned by hand in the settings panel; 0 means "no override". It wins over both the
    // startup variable and the screen, because it is the most recent thing the user said.
    property real uiScaleUser: 0
    property string uiFontFamily: ""

    // What the plot layer was last TOLD, as opposed to what it is drawing: zero means the
    // value is computed from the screen. "Save defaults" stores these rather than the numbers
    // in force, so a scale worked out from this monitor's DPI is not pinned onto the next one.
    property real plotScaleUser: 0
    property real markerSizeUser: 0

    readonly property real uiScale: uiScaleUser > 0     ? uiScaleUser
                                  : uiScaleOverride > 0 ? uiScaleOverride
                                                        : autoScale
    readonly property real fontScale: uiScale / qtTextScale

    function dp(px)     { return Math.round(px * uiScale) }   // pixel lengths
    function pt(points) { return points * fontScale }         // explicit point sizes

    // Base font for the whole window. Qt's default is whatever the platform theme says, which
    // under WSL is 9 pt at 96 dpi -- small on a large panel. Everything inherits this unless
    // it sets its own size, so raising it here raises the whole UI, and OITOOLSGUI_SCALE
    // multiplies it along with the spacing.
    property real baseFontPt: 11
    font.pointSize: pt(baseFontPt)
    // Empty means whatever the platform theme chose, which is the sane default.
    font.family: uiFontFamily.length > 0 ? uiFontFamily : Qt.application.font.family

    // Sized from the SCREEN, not through `dp()`.
    //
    // Passing it through the UI scale ties the window to the widgets, and the two want opposite
    // things: asking for smaller chrome shrank the window to 1050x682 on a 1920x1080 desktop,
    // so there was LESS room for the plot and the tick labels collided. Smaller widgets should
    // buy more space for the data, not less.
    //
    // A wide fraction rather than a square one because every perspective is a settings column
    // beside a plot: Observing puts the Gantt next to the target list, Imaging the
    // reconstruction next to the regularisers. Width is what those layouts spend.
    //
    // The minimum keeps a small screen usable; the fraction keeps a large one from opening a
    // window with its edges off the desktop.
    width:  Math.max(900, Math.round(Screen.desktopAvailableWidth  * 0.83))
    // The plot is one of five stacked regions -- context bar, tabs, plot, pick line, console --
    // so the window needs the height to leave it more than a band.
    height: Math.max(600, Math.round(Screen.desktopAvailableHeight * 0.95))
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

    // What one panel of the per-group view would hold, for whichever observable is selected.
    property string groupingNoun: ""

    function setView() {
        win.groupingNoun = Julia.shell_grouping_noun(kindBox.currentText)
        win.status = Julia.shell_set_view(kindBox.currentText, colorBox.currentText,
                                          logyBox.enabled && logyBox.checked,
                                          panelsBox.enabled && panelsBox.checked)
        win.afterAction()
    }

    // Which observables the loaded file actually holds. Julia counts VALID points, not tables,
    // and pushes the answer to every panel that gates a tick box on it. Without this the boxes
    // stay off after a load and the run buttons never enable — they have no other source of
    // truth, and defaulting them to on would offer observables the file does not contain.
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
            win.setView()
        }
    }

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
        win.refreshPlotKinds()
        win.groupingNoun = Julia.shell_grouping_noun(kindBox.currentText)

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
    Component.onCompleted: {
        // Saved defaults first: they can override the scale, and `shell_ui_scale` below has to
        // be told the value that actually wins.
        applySavedSettings()
        // Hand the scale to Julia before anything is drawn: Makie font and marker sizes are
        // computed there, and they have to follow the same factor as the chrome around them.
        Julia.shell_ui_scale(uiScale, physicalDpi)
        console.log(
        "OITOOLSGUI ui scale: screen=" + physicalDpi.toFixed(1) + " dpi physical, "
        + (Screen.logicalPixelDensity * 25.4).toFixed(1) + " logical"
        + "  devicePixelRatio=" + Screen.devicePixelRatio.toFixed(2)
        + "  autoScale=" + autoScale.toFixed(3)
        + (uiScaleOverride > 0 ? "  override=" + uiScaleOverride.toFixed(3) : "  (no override)")
        + "  -> uiScale=" + uiScale.toFixed(3) + " fontScale=" + fontScale.toFixed(3))
    }

    // Applies whatever "Save defaults" last wrote. Silent when nothing has been saved, which is
    // the normal case -- the built-in defaults are a perfectly good answer, and a settings file
    // that cannot be read must not stop the window from opening.
    function applySavedSettings() {
        var txt = Julia.shell_load_settings()
        if (txt.length === 0) return
        var lines = txt.split("\n")
        for (var i = 0; i < lines.length; ++i) {
            var f = lines[i].split("\t")
            if (f.length !== 2) continue
            var v = parseFloat(f[1])
            switch (f[0]) {
            case "ui_scale":    if (v > 0) win.uiScaleUser = v; break
            case "ui_font":     win.uiFontFamily = f[1]; break
            case "ui_font_pt":  if (v > 0) win.baseFontPt = v; break
            case "plot_scale":  win.plotScaleUser = v > 0 ? v : 0
                                Julia.shell_set_plot_scale(win.plotScaleUser); break
            case "marker_size": win.markerSizeUser = v > 0 ? v : 0
                                Julia.shell_set_marker_size(win.markerSizeUser); break
            }
        }
        console.log("OITOOLSGUI applied saved appearance defaults")
    }

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

            // Settings. Left of the dataset because it governs the whole window rather than
            // any one perspective, and because the dataset is what the eye should land on.
            ToolButton {
                id: settingsButton
                text: "\u2699"                    // GEAR; see settingsPanel for the fallback
                font.pointSize: pt(baseFontPt + 3)
                ToolTip.visible: hovered
                ToolTip.text: "Appearance settings"
                onClicked: {
            if (settingsPanel.opened) { settingsPanel.close(); return }
            // Read the live values back, so the panel opens showing what is in force rather
            // than what it last displayed.
            uiScaleSpin.value = Math.round(win.uiScale * 100)
            uiFontSizeSpin.value = Math.round(win.baseFontPt)
            // The combo only reports what was last picked FROM it, so a family that arrived
            // from the saved defaults has to be found in the list by hand. -1 for a name that
            // is not on the list, and 0 for none at all, both mean "(system default)".
            uiFontBox.currentIndex = Math.max(0, uiFontBox.model.indexOf(win.uiFontFamily))
            // Plot scale and marker size live in Julia, so ask it rather than trusting the
            // number the spin box was built with.
            var f = Julia.shell_plot_scale().split("\t")
            if (f.length === 3) {
                plotScaleSpin.value  = Math.round(parseFloat(f[0]) * 100)
                win.plotScaleUser    = parseFloat(f[1])
                win.markerSizeUser   = parseFloat(f[2])
                markerSpin.value     = Math.round(win.markerSizeUser)
            }
            settingsPanel.open()
        }
            }

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

    // ── appearance settings ───────────────────────────────────────────────────
    //
    // The controls show what is CURRENTLY in force, and the scale, font and plot controls
    // change it live. Plot font and theme are shown but disabled, each for a stated reason.
    // "Save defaults" writes the lot to a per-user file that the window reads at startup.
    Popup {
        id: settingsPanel
        // Centred. Anchored to the top-left it landed on the Observing settings column, which
        // is also a stack of labelled rows in a light panel -- close enough in appearance that
        // it was not obvious which one had focus. The middle of the window belongs to nothing
        // else, so a panel there reads as a dialog.
        x: Math.round((win.width  - width)  / 2)
        y: Math.round((win.height - height) / 2)
        width: dp(430)
        // Dimmed behind, for the same reason: it separates the panel from whatever it covers.
        modal: true
        dim: true
        focus: true
        closePolicy: Popup.CloseOnEscape | Popup.CloseOnPressOutsideParent
        padding: dp(12)

        ColumnLayout {
            width: parent.width
            spacing: dp(10)

            RowLayout {
                Layout.fillWidth: true
                Label { text: "Appearance"; font.bold: true }
                Item { Layout.fillWidth: true }
                Label {
                    text: "plot font and theme need a restart"
                    color: "#888"
                    font.pointSize: pt(baseFontPt - 2)
                }
            }

            GridLayout {
                columns: 3
                columnSpacing: dp(8)
                rowSpacing: dp(6)
                Layout.fillWidth: true

                // ── UI ────────────────────────────────────────────────────────
                Label { text: "UI scale"; color: "#666" }
                SpinBox {
                    id: uiScaleSpin
                    from: 50; to: 400; stepSize: 5
                    value: Math.round(win.uiScale * 100)
                    editable: true
                    textFromValue: function (v) { return (v / 100).toFixed(2) }
                    valueFromText: function (t) { return Math.round(parseFloat(t) * 100) }
                    Layout.fillWidth: true
                    // Everything sized through dp()/pt() rebinds, so the window resizes as the
                    // number changes rather than on close.
                    onValueModified: win.uiScaleUser = value / 100
                }
                Button {
                    text: "auto"
                    enabled: win.uiScaleUser > 0
                    ToolTip.visible: hovered
                    ToolTip.text: uiScaleOverride > 0
                        ? "back to OITOOLSGUI_SCALE=" + uiScaleOverride.toFixed(2)
                        : "back to the value computed from " + physicalDpi.toFixed(0) + " dpi"
                    onClicked: win.uiScaleUser = 0
                }

                Label { text: "UI font"; color: "#666" }
                ComboBox {
                    id: uiFontBox
                    model: ["(system default)", "DejaVu Sans", "Noto Sans", "Liberation Sans",
                            "JuliaMono"]
                    Layout.fillWidth: true
                    onActivated: win.uiFontFamily = currentIndex === 0 ? "" : currentText
                }
                SpinBox {
                    id: uiFontSizeSpin
                    from: 6; to: 24; stepSize: 1
                    value: Math.round(baseFontPt)
                    editable: true
                    onValueModified: win.baseFontPt = value
                }

                // ── plots ─────────────────────────────────────────────────────
                Label { text: "Plot font"; color: "#666" }
                ComboBox {
                    id: plotFontBox
                    // Shown, not settable. Makie takes the font at Figure CONSTRUCTION; every
                    // route to change it afterwards throws `Failed to resolve
                    // data_boundingbox` out of its compute graph -- tried on the axis label,
                    // tick label and colorbar attributes alike. Changing it means restarting.
                    model: ["DejaVu Sans", "JuliaMono", "Noto Sans", "Liberation Sans",
                            "(Makie default)"]
                    enabled: false
                    Layout.fillWidth: true
                    ToolTip.visible: hovered
                    ToolTip.text: "Set at startup; Makie cannot change a figure's font once it " +
                                  "exists. Edit PLOT_FONT in src/gui/plots.jl."
                }
                SpinBox {
                    id: plotScaleSpin
                    from: 50; to: 500; stepSize: 5
                    value: 119                     // replaced on open by the live value
                    editable: true
                    textFromValue: function (v) { return (v / 100).toFixed(2) }
                    valueFromText: function (t) { return Math.round(parseFloat(t) * 100) }
                    onValueModified: {
                        win.plotScaleUser = value / 100
                        Julia.shell_set_plot_scale(win.plotScaleUser)
                    }
                }

                Label { text: "Plot symbols"; color: "#666" }
                SpinBox {
                    id: markerSpin
                    from: 0; to: 30; stepSize: 1
                    value: 0
                    editable: true
                    Layout.fillWidth: true
                    // 0 means "whatever the plot chooses", which differs per view -- uv
                    // coverage draws smaller points than an observable plot.
                    textFromValue: function (v) { return v === 0 ? "auto" : String(v) }
                    valueFromText: function (t) { return t === "auto" ? 0 : parseInt(t) }
                    onValueModified: {
                        win.markerSizeUser = value
                        Julia.shell_set_marker_size(value)
                    }
                }
                Label { text: "px"; color: "#888"; font.pointSize: pt(baseFontPt - 2) }

                // ── theme ─────────────────────────────────────────────────────
                Label { text: "Theme"; color: "#666" }
                ComboBox {
                    id: themeBox
                    model: ["Follow the desktop", "Light", "Dark"]
                    enabled: false
                    Layout.fillWidth: true
                    ToolTip.visible: hovered
                    ToolTip.text: "The panels hardcode light colours, so choosing Dark would " +
                                  "leave parts of the window unreadable. Needs a theming pass."
                }
                Label {
                    // Honest about the state of it: the panels hardcode light colours, so a
                    // dark desktop session leaves parts of the window unreadable today.
                    text: "dark is unstyled"
                    color: "#888"
                    font.pointSize: pt(baseFontPt - 2)
                }
            }

            RowLayout {
                Layout.fillWidth: true
                Label {
                    id: savedLabel
                    Layout.fillWidth: true
                    elide: Text.ElideMiddle
                    color: "#888"
                    font.pointSize: pt(baseFontPt - 2)
                }
                Button {
                    text: "Save defaults"
                    ToolTip.visible: hovered
                    ToolTip.text: "write these as the startup defaults for every project"
                    onClicked: {
                        // The OVERRIDES, not the values in force: a scale of 0 means "work it
                        // out from the screen", which is the right thing to carry to a machine
                        // with a different one.
                        var path = Julia.shell_save_settings(
                            [ "ui_scale\t"    + win.uiScaleUser,
                              "ui_font\t"     + win.uiFontFamily,
                              "ui_font_pt\t"  + win.baseFontPt,
                              "plot_scale\t"  + win.plotScaleUser,
                              "marker_size\t" + win.markerSizeUser ].join("\n"))
                        savedLabel.text = path.length > 0 ? "saved to " + path
                                                          : "could not save — see the console"
                    }
                }
                Button { text: "Close"; onClicked: settingsPanel.close() }
            }
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
                        textRole: "name"
                        // `diffphi` is `visphi` against WAVELENGTH — OIFITS carries differential
                        // phase in OI_VIS with PHITYP=differential, so there is no separate
                        // field, only a more useful axis.
                        //
                        // Entries the loaded file cannot supply are shown but disabled, rather
                        // than hidden: which observables a file carries is worth seeing, and a
                        // list that changes length between files is hard to aim at.
                        model: kindModel
                        delegate: ItemDelegate {
                            required property var model
                            required property int index
                            width: kindBox.width
                            text: model.name
                            enabled: model.available
                            highlighted: kindBox.highlightedIndex === index
                            // The greyed rows still need a reason, or an unavailable kind looks
                            // like a broken menu rather than a file without that table.
                            ToolTip.visible: hovered && !model.available
                            ToolTip.text: "this file has no " + model.name + " data"
                        }
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
                    CheckBox {
                        id: panelsBox
                        // Labelled from Julia, because the grouping is not always a baseline:
                        // closure quantities group by triplet and flux by station, and a tick
                        // reading "per baseline" over a grid of triangles is simply wrong.
                        text: "per " + (win.groupingNoun.length > 0 ? win.groupingNoun : "group")
                        // uv coverage has no groups to split by — it is geometry, not an
                        // observable — and a single-wavelength file gives one point per panel.
                        enabled: win.groupingNoun.length > 0
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
                        // Off by default: the transcript is still recorded either way -- this
                        // only decides whether the pane is filled -- and the plot is what the
                        // window should open on.
                        checked: false
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
