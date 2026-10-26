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

    // ── a fixed light palette ─────────────────────────────────────────────────
    //
    // Every panel, table and label below paints itself with a fixed light colour, and the
    // Makie canvases are light too. Qt's controls, left alone, follow the SYSTEM theme -- so on
    // a dark desktop the window came out mixed: dark chrome around white tables and white
    // plots, with grey-on-grey text that was hard to read.
    //
    // `prefer_light_controls!` pins the style to Basic for the same reason; this pins the
    // colours the style derives from, which is the half a style choice does not cover. The
    // values are the ones already used throughout this file, so nothing shifts on a light
    // desktop -- they only stop a dark one leaking in.
    //
    // If a dark theme is ever wanted, this block plus `style_axis!` for the figures are the two
    // places that would have to change, and the 71 literals below would have to become palette
    // roles first.
    palette {
        window:          "#f4f4f4"
        windowText:      "#222222"
        base:            "#ffffff"
        alternateBase:   "#fbfbfb"
        text:            "#222222"
        button:          "#efefef"
        buttonText:      "#222222"
        placeholderText: "#888888"
        light:           "#ffffff"
        midlight:        "#eeeeee"
        mid:             "#dddddd"
        dark:            "#999999"
        shadow:          "#666666"
        highlight:       "#3874d8"
        highlightedText: "#ffffff"
        toolTipBase:     "#ffffe1"
        toolTipText:     "#222222"
        brightText:      "#ffffff"
        link:            "#0645ad"
        // The roles above apply to every colour group, so a disabled control would keep
        // the same dark text as an enabled one.
        disabled {
            windowText:  "#a0a0a0"
            text:        "#a0a0a0"
            buttonText:  "#a0a0a0"
        }
    }

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

    // What the screen actually is: its density in DEVICE pixels per inch.
    //
    // `Screen.pixelDensity` is per logical pixel, so it is multiplied by the device pixel
    // ratio: a 16" 2560x1600 panel at 150% reports 126 and is 189. On Linux the panel's own
    // EDID wins when one matches the screen's mode, since under X11 Qt takes the size in mm
    // from the X server, which can report a nominal 96-100 dpi for any panel. Under WSL, which
    // has no EDID and the same nominal density, Windows' display scale stands for it.
    readonly property real screenDpr: Screen.devicePixelRatio > 0 ? Screen.devicePixelRatio : 1.0
    readonly property real panelDpi: Julia.shell_panel_dpi(Math.round(Screen.width * screenDpr),
                                                           Math.round(Screen.height * screenDpr))
    readonly property real deviceDpi: panelDpi > 0 ? panelDpi
                                                   : Screen.pixelDensity * 25.4 * screenDpr
    // Per logical pixel, the unit `dp()` and the plot scale work in.
    readonly property real physicalDpi: deviceDpi / screenDpr

    // The size on screen is uiScale x screenDpr, set by the device density alone:
    // `refScale` at `refDpi`, growing as density^`dpiExponent`. Judged by eye: 0.73 on a
    // 92.6 dpi 1920x1080 desktop at 100%, and 0.62 on a 189 dpi 16" 2560x1600 laptop at 150%
    // (0.93 on screen), which the exponent 0.34 fits. The OS scaling is then divided out, so
    // one panel gets the same size whatever scaling or platform it runs under; on Linux the
    // launcher sets that scaling from the panel when the desktop does not (`configure_qt_scale!`).
    //
    // This governs the CHROME only. The window is sized from the screen (below) and plot text
    // from the density (`live_plot_scale` in Julia).
    readonly property real refDpi:      92.6
    readonly property real refScale:    0.73
    readonly property real dpiExponent: 0.34
    readonly property real autoScale:
        Math.max(0.5, Math.min(4.0,
            (deviceDpi > 0 ? refScale * Math.pow(deviceDpi / refDpi, dpiExponent) : refScale)
            / screenDpr))

    // Turned by hand in the settings panel; 0 means "no override". It wins over both the
    // startup variable and the screen, because it is the most recent thing the user said.
    // The out-of-the-box appearance, in one place. Both the initial values below and the
    // settings panel's "Reset to defaults" read it, so the button cannot drift from the
    // declarations — a reset that restored the wrong number would be worse than none.
    //
    // Zero is not a size anywhere here: it is "unset", and the shipped value behind it lives in
    // Julia (`PLOT_SCALE_AT_REF_DPI`, `ZOOM_PER_DETENT`, each plot's own marker size) or is
    // computed from the screen. So the reset sends zero rather than a number, and the boxes are
    // filled from what Julia then reports — repeating 1.25 here would pin today's constant as a
    // user override and survive a change to it.
    // ── the window's own font ─────────────────────────────────────────────────
    //
    // Qt resolves a family name through the SYSTEM font database, so naming one only works
    // where the machine has it -- stock Windows has no Noto. These two faces come out of the
    // MakieAssets artifact, which the bundle already carries for the plots, and a FontLoader
    // puts them in Qt's database at startup: the same UI font on every platform, nothing extra
    // shipped. Empty when the assets cannot be found, and then the platform's font is used.
    readonly property var shippedFontFiles: {
        var t = Julia.shell_ui_font_files()
        return t.length > 0 ? t.split("\n") : []
    }
    FontLoader { id: shippedFontRegular
                 source: win.shippedFontFiles.length > 0 ? win.shippedFontFiles[0] : "" }
    FontLoader { id: shippedFontBold
                 source: win.shippedFontFiles.length > 1 ? win.shippedFontFiles[1] : "" }
    // Bold is loaded into the same family, so `font.bold` picks the real face rather than
    // having Qt smear the regular one.
    readonly property string shippedFontFamily:
        shippedFontRegular.status === FontLoader.Ready ? shippedFontRegular.name : ""

    readonly property var appearanceDefaults: ({
        uiScaleUser: 0, windowScaleUser: 0, uiFontFamily: "", baseFontPt: 11, plotScaleUser: 0, markerSizeUser: 0,
        // The shipped style is Julia's to name (`DEFAULT_CONTROLS_STYLE`), so it is asked for
        // rather than repeated here: a reset has to restore what the window would ship with,
        // not a literal that was true when this line was written.
        zoomStepUser: 0, controlsStyle: Julia.shell_default_controls_style(),
        startTab: "exploring", startFullscreen: false, rememberGeometry: false,
        obsTimeSystem: "LST", obsElevMin: 25.0, obsElevMax: 85.0
    })

    // A multiple of the size computed for the screen (`autoScale`, or OITOOLSGUI_SCALE when
    // set): 0 means none, so `uiSize` is 1 by default on every screen.
    property real uiScaleUser: appearanceDefaults.uiScaleUser
    readonly property real uiSize: uiScaleUser > 0 ? uiScaleUser : 1.0
    // Read by Julia from the saved settings before the window exists (`_initial_tab`,
    // `fullscreenOnStart`), so these take effect at the next launch.
    property string startTab: appearanceDefaults.startTab
    property bool startFullscreen: appearanceDefaults.startFullscreen
    // Written at close, with the window's geometry, when set.
    property bool rememberGeometry: appearanceDefaults.rememberGeometry
    // Qt's own scale (QT_SCALE_FACTOR), read by the launcher before Qt starts, so a change
    // takes effect at the next launch. 0 is automatic; `windowScaleRunning` is the one in force.
    property real windowScaleUser: appearanceDefaults.windowScaleUser
    readonly property real windowScaleRunning: Julia.shell_window_scale()
    // What "auto" gives at the next launch, and what that launch will use.
    readonly property real windowScaleAuto: Julia.shell_auto_window_scale()
    readonly property real windowScaleTarget: windowScaleUser > 0 ? windowScaleUser : windowScaleAuto
    property string uiFontFamily: appearanceDefaults.uiFontFamily
    // Qt Quick Controls style. Read by `apply_controls_style!` from the settings file BEFORE
    // the QML is loaded, so a change here takes effect at the next launch, not this one.
    property string controlsStyle: appearanceDefaults.controlsStyle

    // What the plot layer was last TOLD, as opposed to what it is drawing: zero means the
    // value is computed from the screen. "Save config" stores these rather than the numbers
    // in force, so a scale worked out from this monitor's DPI is not pinned onto the next one.
    property real plotScaleUser: appearanceDefaults.plotScaleUser
    property real markerSizeUser: appearanceDefaults.markerSizeUser
    // How far one wheel detent zooms a plot. A setting because it depends on the pointing
    // device as much as on taste: a detented wheel, a free-spinning one and a touchpad all
    // deliver different amounts of scroll for the same gesture. Filled from Julia when the
    // panel opens, since it is Julia that zooms.
    property real zoomStepUser: appearanceDefaults.zoomStepUser
    // Where this user keeps their own data. Saved in the config, and consulted by `picker_start`
    // and the window's opening folder ahead of the shipped demo data -- so after setting it once,
    // every Open and Save starts where the work is rather than inside the package.
    property string workDirectory: ""
    // Applied at once, not at "Save config": that button is for the NEXT launch, and a folder
    // that appears to do nothing until you press an unrelated button reads as broken.
    onWorkDirectoryChanged: Julia.shell_set_work_directory(workDirectory)

    readonly property real uiScale: (uiScaleOverride > 0 ? uiScaleOverride : autoScale) * uiSize
    readonly property real fontScale: uiScale / qtTextScale

    function dp(px)     { return Math.round(px * uiScale) }   // pixel lengths
    function pt(points) { return points * fontScale }         // explicit point sizes

    // Base font for the whole window. Qt's default is whatever the platform theme says, which
    // under WSL is 9 pt at 96 dpi -- small on a large panel. Everything inherits this unless
    // it sets its own size, so raising it here raises the whole UI, and OITOOLSGUI_SCALE
    // multiplies it along with the spacing.
    property real baseFontPt: appearanceDefaults.baseFontPt
    font.pointSize: pt(baseFontPt)
    // Empty means the window's own default: the shipped Noto face, or the platform's font if
    // the assets could not be found.
    font.family: uiFontFamily.length > 0     ? uiFontFamily
               : shippedFontFamily.length > 0 ? shippedFontFamily
                                              : Qt.application.font.family

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


    // Saves the Exploring plot. Each perspective carries its own, because the picker's title,
    // suffix and target all belong to the plot that asked for it.
    FilePicker {
        id: savePngDialog
        property real pxw: 1200
        property real pxh: 900
        uiScale: win.uiScale; fontScale: win.fontScale; baseFontPt: win.baseFontPt
        title: "Save plot as PNG"
        saveMode: true
        defaultSuffix: "png"
        filters: [{ label: "PNG images (*.png)", patterns: "*.png" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) {
            Julia.shell_save_figure("explore", path, savePngDialog.pxw, savePngDialog.pxh, false)
            win.refreshConsole()
        }
    }

    // Explore's Export: the dataset on screen, as filtered, to a new OIFITS file.
    FilePicker {
        id: exportDataDialog
        uiScale: win.uiScale; fontScale: win.fontScale; baseFontPt: win.baseFontPt
        title: "Save the dataset as OIFITS"
        saveMode: true
        defaultSuffix: "oifits"
        filters: [{ label: "OIFITS files (*.oifits *.fits)", patterns: "*.oifits *.fits" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) {
            win.status = Julia.shell_export_oifits(path)
            win.refreshConsole()
        }
    }
    function exportData() {
        exportDataDialog.suggestedName = Julia.picker_suggest(
            (datasetBox.currentText || "data").replace(/\.(oi)?fits$/i, ""), "oifits")
        exportDataDialog.openAt("")
    }

    // The panel's own size, which Julia uses only when it cannot measure the live scene: the
    // file has to be laid out like the picture on screen, and Makie sizes text in points, so
    // asking for a bigger figure would re-lay it out rather than enlarge it. Resolution comes
    // from `px_per_unit` at save time instead. Capped because the figure is rendered into a
    // real framebuffer and a software GL stack refuses the very large ones.
    function savePng(which, area) {
        savePngDialog.pxw = Math.min(2400, Math.max(640, area.width))
        savePngDialog.pxh = Math.min(1800, Math.max(480, area.height))
        // The dataset and the view are what distinguish one of these files from the next.
        savePngDialog.suggestedName = Julia.picker_suggest(
            (datasetBox.currentText || "plot") + "\t" + (exploreTab.currentKind || ""), "png")
        savePngDialog.openAt("")
    }

    // Every action funnels through here so the console cannot go stale: the pane is the
    // record of what happened, and a record that updates only sometimes is worse than none.
    function refreshConsole() { if (logToggle.checked) logArea.text = Julia.shell_console() }

    // Repaint the plot. This is NOT optional and it is easy to forget: the shell now updates
    // Makie Observables rather than inserting plots, and an Observable assignment does not by
    // itself ask Qt for a new frame. The old code repainted only as a side effect of creating
    // plots -- the same side effect that was corrupting the GL context. MakieArea's own mouse
    // handlers call update() for exactly this reason.
    function redrawPlot() { exploreTab.redraw() }

    function afterAction() { refreshConsole(); redrawPlot(); refreshObservables()
                             exploreTab.refreshCompareAvailability()
                             exploreTab.datasetName = datasetBox.currentText
                             exploreTab.refreshSummary()
                             exploreTab.refreshFilter() }

    // Plot scale, marker size and zoom step live in Julia, so ask it rather than trusting the
    // numbers the spin boxes were built with. The first field is the scale in FORCE, which is
    // what the box shows: `plotScaleUser` may be 0 for "computed from the screen", and a box
    // reading 0 would say nothing about what the plots are doing.
    function readPlotSettings() {
        var f = Julia.shell_plot_scale().split("\t")
        if (f.length !== 5) return
        plotScaleSpin.value  = Math.round(parseFloat(f[0]) * 100)
        win.plotScaleUser    = parseFloat(f[1])
        win.markerSizeUser   = parseFloat(f[2])
        markerSpin.value     = Math.round(win.markerSizeUser)
        zoomSpin.value       = Math.round(parseFloat(f[3]) * 100)
        win.zoomStepUser     = parseFloat(f[4])
    }

    // The settings boxes whose values live elsewhere -- the start-up choices here, the
    // Observing ones in its tab -- set from those, since a box the user changed no longer
    // follows its binding.
    function refreshSettingsBoxes() {
        startTabBox.currentIndex = Math.max(0, startTabBox.model.map(function (t) { return t.toLowerCase() })
                                                              .indexOf(win.startTab))
        startFullscreenBox.checked  = win.startFullscreen
        rememberGeometryBox.checked = win.rememberGeometry
        controlsStyleBox.currentIndex = Math.max(0, controlsStyleBox.model.indexOf(win.controlsStyle))
        obsTimeBox.currentIndex  = Math.max(0, observeTab.timeSystems.indexOf(observeTab.timeSystem))
        obsElevMinField.text     = observeTab.altLimit.toFixed(1)
        obsElevMaxField.text     = observeTab.altMax.toFixed(1)
    }

    // The Gantt's clock, from the settings panel: as the Observing tab's own box does, the
    // chart on screen is redrawn in it.
    function setObserveTimeSystem(name) {
        if (observeTab.timeSystem === name) return
        observeTab.timeSystem = name
        if (observeTab.hasPlan) observeTab.computePlan()
    }

    // The geometry goes to the settings file at close when the user asked for it to be kept;
    // a full-screen window has none worth keeping.
    onClosing: function (close) {
        if (!rememberGeometry || visibility !== Window.Windowed) return
        Julia.shell_save_settings(["remember_geometry\t1", "window_x\t" + x, "window_y\t" + y,
                                   "window_width\t" + width, "window_height\t" + height].join("\n"))
    }

    // Which observables the loaded file actually holds. Julia counts VALID points, not tables,
    // and pushes the answer to every panel that gates a tick box on it. Without this the boxes
    // stay off after a load and the run buttons never enable — they have no other source of
    // truth, and defaulting them to on would offer observables the file does not contain.
    function refreshObservables() {
        var have = { v2: false, t3amp: false, t3phi: false, cvis: false,
                     visamp: false, visphi: false, flux: false, diffvis: false }
        var s = Julia.shell_observables()
        if (s.length > 0) {
            let parts = s.split(",")
            for (let i = 0; i < parts.length; ++i) {
                let kv = parts[i].split("=")
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
        exploreTab.refreshPlotKinds()
        exploreTab.refreshGroupingNoun()

        // Geometry the data suggests. A constant pixel size is not a neutral default: too
        // coarse and the image cannot represent what the data resolves, too fine and it is
        // unconstrained at the pixel scale.
        var g = Julia.shell_image_defaults()
        if (g.length > 0) {
            let gp = g.split(",")
            for (let k = 0; k < gp.length; ++k) {
                let kv2 = gp[k].split("=")
                if (kv2.length !== 2) continue
                if (kv2[0] === "pixsize") imageTab.pixsize = parseFloat(kv2[1])
                if (kv2[0] === "nx")      imageTab.nx      = parseInt(kv2[1])
            }
        }
        // A new file arrives unbinned, and its channel count is what the binning control has
        // to offer -- so this comes before the plan, which is built per bin.
        imageTab.binMode = "single"
        imageTab.refreshBins()
        // The plan is keyed on the dataset as well as the geometry, so a new file needs one
        // even when nx and the pixel size happen to be unchanged.
        imageTab.refreshPlan()

        modelTab.haveV2      = have.v2
        modelTab.haveT3amp   = have.t3amp
        modelTab.haveT3phi   = have.t3phi
        modelTab.haveVisamp  = have.visamp
        modelTab.haveVisphi  = have.visphi
        modelTab.haveFlux    = have.flux
        modelTab.haveDiffvis = have.diffvis
        modelTab.hasDataset  = s.length > 0
        exploreTab.hasDataset = s.length > 0
    }

    // Printed once at startup so the scale is diagnosable on a machine nobody can inspect
    // remotely. If the window looks wrong, this line says whether the screen asked for it
    // (autoScale) or an override did, and OITOOLSGUI_SCALE is how you argue with it.
    Component.onCompleted: {
        // Saved defaults first: they can override the scale, and `shell_ui_scale` below has to
        // be told the value that actually wins.
        applySavedSettings()
        // Which optional engines this build can run. A session property, not a dataset one:
        // it depends on whether the package extension loaded, so it is read once here rather
        // than left hardcoded false in ImageTab as both of these were.
        var eng = Julia.shell_optional_engines().split(",")
        for (let i = 0; i < eng.length; ++i) {
            let kv = eng[i].split("=")
            if (kv[0] === "vi")        imageTab.oiviAvailable    = (kv[1] === "1")
            if (kv[0] === "tempering") { imageTab.pigeonsAvailable = (kv[1] === "1")
                                         modelTab.pigeonsAvailable = (kv[1] === "1") }
            if (kv[0] === "corner")    modelTab.cornerAvailable  = (kv[1] === "1")
        }
        // Hand the scale to Julia before anything is drawn: Makie font and marker sizes are
        // computed there, and they have to follow the same factor as the chrome around them.
        Julia.shell_ui_scale(uiScale, physicalDpi)
        if (verboseStartup) console.log(
        "OITOOLSGUI ui scale: screen=" + (Screen.pixelDensity * 25.4).toFixed(1) + " dpi physical, "
        + (Screen.logicalPixelDensity * 25.4).toFixed(1) + " logical"
        + "  devicePixelRatio=" + Screen.devicePixelRatio.toFixed(2)
        + "  panel=" + (panelDpi > 0 ? panelDpi.toFixed(1) + " dpi" : "unknown")
        + "  device=" + deviceDpi.toFixed(1) + " dpi"
        + "  autoScale=" + autoScale.toFixed(3)
        + (uiScaleOverride > 0 ? "  override=" + uiScaleOverride.toFixed(3) : "  (no override)")
        + "  -> uiScale=" + uiScale.toFixed(3) + " fontScale=" + fontScale.toFixed(3))
    }

    // Applies whatever "Save config" last wrote. Silent when nothing has been saved, which is
    // the normal case -- the built-in defaults are a perfectly good answer, and a settings file
    // that cannot be read must not stop the window from opening.
    function applySavedSettings() {
        var txt = Julia.shell_load_settings()
        if (txt.length === 0) return
        var lines = txt.split("\n")
        var geo = {}
        for (let i = 0; i < lines.length; ++i) {
            let f = lines[i].split("\t")
            if (f.length !== 2) continue
            let v = parseFloat(f[1])
            switch (f[0]) {
            // "ui_size", not the older "ui_scale", which held the factor itself
            case "ui_size":     if (v > 0) win.uiScaleUser = v; break
            case "start_tab":   win.startTab = f[1]; break
            case "start_fullscreen":  win.startFullscreen  = (v === 1); break
            case "remember_geometry": win.rememberGeometry = (v === 1); break
            case "window_x": case "window_y": case "window_width": case "window_height":
                                geo[f[0]] = v; break
            case "obs_time_system": observeTab.timeSystem = f[1]; break
            case "obs_elev_min": if (v >= 0 && v <= 90) observeTab.altLimit = v; break
            case "obs_elev_max": if (v >= 0 && v <= 90) observeTab.altMax = v; break
            case "window_scale": win.windowScaleUser = v > 0 ? v : 0; break
            case "ui_font":     win.uiFontFamily = f[1]; break
            case "controls_style": win.controlsStyle = f[1]; break
            case "ui_font_pt":  if (v > 0) win.baseFontPt = v; break
            case "plot_scale":  win.plotScaleUser = v > 0 ? v : 0
                                Julia.shell_set_plot_scale(win.plotScaleUser); break
            case "marker_size": win.markerSizeUser = v > 0 ? v : 0
                                Julia.shell_set_marker_size(win.markerSizeUser); break
            case "zoom_step":   if (v > 1) { win.zoomStepUser = v
                                             Julia.shell_set_zoom_step(v) } break
            case "work_directory": win.workDirectory = f[1]; break
            }
        }
        // The saved geometry, kept on the desktop: a window last closed on a monitor that is
        // no longer attached must not open off screen.
        if (rememberGeometry && !fullscreenOnStart && geo.window_width > 0 && geo.window_height > 0) {
            win.width  = Math.min(Math.max(900, geo.window_width),  Screen.desktopAvailableWidth)
            win.height = Math.min(Math.max(600, geo.window_height), Screen.desktopAvailableHeight)
            if (isFinite(geo.window_x) && isFinite(geo.window_y)) {
                win.x = Math.max(0, Math.min(geo.window_x, Screen.desktopAvailableWidth  - win.width))
                win.y = Math.max(0, Math.min(geo.window_y, Screen.desktopAvailableHeight - win.height))
            }
        }
        if (verboseStartup) console.log("OITOOLSGUI applied saved appearance defaults")
    }

    // Reads the dataset list back from the session and selects the newest entry. Called both
    // after a load and at startup, because a Session may already hold datasets before the
    // window exists -- `oitoolsgui.jl file.oifits` loads them and then calls gui().
    // After Julia has moved to another dataset (a filter was applied), the box follows it.
    function selectCurrentDataset() {
        refreshDatasets()
        afterAction()
    }

    function refreshDatasets() {
        var names = Julia.shell_dataset_names().split("\n").filter(function (s) { return s.length > 0 })
        datasetBox.model = names
        // The dataset Julia is showing, not the newest: with several loaded before the window
        // opened, the newest is not the one every perspective is looking at.
        datasetBox.currentIndex = Math.max(0, Julia.shell_current_dataset() - 1)
        // Observing's "Use existing OIFITS as template" simulates through the loaded file's uv
        // points, so it needs to know there IS one — and to say so beside the button rather
        // than after Simulate is pressed.
        observeTab.datasetCount = names.length
    }

    // ── context bar: what every perspective is looking at ─────────────────────
    header: ToolBar {
        RowLayout {
            anchors.fill: parent
            anchors.leftMargin: dp(8); anchors.rightMargin: dp(8)
            spacing: dp(10)

            // About: the version, licence and credits, styled as the settings button beside it.
            ToolButton {
                id: aboutButton
                objectName: "aboutButton"
                text: "?"
                font.pointSize: pt(baseFontPt + 8)
                font.bold: true
                implicitWidth: dp(38)
                implicitHeight: dp(38)
                background: Rectangle {
                    radius: dp(3)
                    border.color: "#000000"
                    border.width: 1
                    color: aboutButton.down ? "#d4d4d4"
                           : aboutButton.hovered ? "#ececec" : "transparent"
                }
                ToolTip.visible: hovered
                ToolTip.text: "About OITOOLS"
                onClicked: aboutPanel.opened ? aboutPanel.close() : aboutPanel.openWith(Julia.shell_about())
            }

            // Settings. Left of the dataset because it governs the whole window rather than
            // any one perspective, and because the dataset is what the eye should land on.
            ToolButton {
                id: settingsButton
                text: "\u2699"                    // GEAR; see settingsPanel for the fallback
                // Sized well above the toolbar's text. This is the only control in the bar that
                // is a symbol rather than a word, so at text size it read as punctuation next
                // to "Dataset:" instead of as the button it is.
                font.pointSize: pt(baseFontPt + 11)
                // The glyph is taller than a ToolButton's default content box, so give it room
                // rather than letting the bar clip its top and bottom.
                implicitWidth: dp(38)
                implicitHeight: dp(38)
                // A ToolButton draws no frame until it is hovered, so the one symbol in the
                // bar had nothing marking it as pressable. The outline says so standing
                // still; the fill is what still answers the pointer.
                background: Rectangle {
                    radius: dp(3)
                    border.color: "#000000"
                    border.width: 1
                    color: settingsButton.down ? "#d4d4d4"
                           : settingsButton.hovered ? "#ececec" : "transparent"
                }
                ToolTip.visible: hovered
                ToolTip.text: "Appearance settings"
                onClicked: {
            if (settingsPanel.opened) { settingsPanel.close(); return }
            // Read the live values back, so the panel opens showing what is in force rather
            // than what it last displayed.
            uiScaleSpin.value = Math.round(win.uiSize * 100)
            windowScaleSpin.value = Math.round(win.windowScaleTarget * 100)
            win.refreshSettingsBoxes()
            uiFontSizeSpin.value = Math.round(win.baseFontPt)
            // The combo only reports what was last picked FROM it, so a family that arrived
            // from the saved defaults has to be found in the list by hand. -1 for a name that
            // is not on the list, and 0 for none at all, both mean "(system default)".
            uiFontBox.currentIndex = Math.max(0, uiFontBox.model.indexOf(win.uiFontFamily))
            readPlotSettings()
            settingsPanel.open()
        }
            }

            Label { text: "Dataset:"; font.bold: true }
            ComboBox {
                id: datasetBox
                objectName: "datasetBox"
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
                onClicked: openDialog.openAt(Julia.shell_initial_folder())
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
        title: "Open OIFITS  (Ctrl- or Shift-click: several files as one dataset)"
        multiSelect: true
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

    // A settings group's title, with what it applies to or when on the right.
    component SettingsHeader: RowLayout {
        id: header
        property string title
        property string note
        Layout.fillWidth: true
        Layout.topMargin: dp(4)
        Label { text: header.title; font.bold: true }
        Item { Layout.fillWidth: true }
        Label { text: header.note; color: "#888"; font.pointSize: pt(baseFontPt - 2) }
    }

    // ── appearance settings ───────────────────────────────────────────────────
    //
    // The controls show what is CURRENTLY in force, and the scale, font and plot controls
    // change it live. The plot font is shown but disabled, for a stated reason, and the
    // controls style takes effect at the next launch because Qt reads it once, at load.
    // "Save config" writes the lot to a per-user file that the window reads at startup.
    // About: version, licence, author and acknowledgements. Fixed widths on the panel and its
    // wrapped labels, neither read from the other (see the settings panel on binding loops).
    Popup {
        id: aboutPanel
        objectName: "aboutPanel"
        x: Math.max(dp(12), Math.round((parent.width  - width)  / 2))
        y: Math.max(dp(12), Math.round((parent.height - height) / 2))
        width: dp(460)
        modal: true
        dim: true
        closePolicy: Popup.CloseOnEscape | Popup.CloseOnPressOutsideParent
        padding: dp(16)
        property var info: ({})
        function openWith(t) {
            var kv = {}
            var lines = t.split("\n")
            for (let i = 0; i < lines.length; ++i) {
                let f = lines[i].split("\t"); kv[f[0]] = f.length > 1 ? f[1] : ""
            }
            info = kv
            open()
        }
        ColumnLayout {
            width: dp(428)
            spacing: dp(8)
            Label { text: "OITOOLS"; font.bold: true; font.pointSize: pt(baseFontPt + 6) }
            Label {
                objectName: "aboutVersion"
                text: "Version " + (aboutPanel.info.version || "?") + "  ·  Julia " + (aboutPanel.info.julia || "?")
            }
            Label {
                objectName: "aboutLicense"
                Layout.preferredWidth: dp(428)
                wrapMode: Text.WordWrap
                visible: (aboutPanel.info.license || "").length > 0
                text: "License: " + (aboutPanel.info.license || "")
            }
            Label {
                Layout.preferredWidth: dp(428)
                wrapMode: Text.WordWrap
                text: "Developed by Prof. Fabien Baron, Georgia State University"
            }
            Label { text: "Acknowledgements"; font.bold: true; Layout.topMargin: dp(4) }
            Label {
                Layout.preferredWidth: dp(428)
                wrapMode: Text.WordWrap
                // each credit kept on one line; the lines break between them
                text: ["Jean-Marie Mariotti Center", "Prof. Ryan Norris (NMT)", "Dr. Arturo Martinez",
                       "PMOIRED by Dr. Antoine Mérand", "VMLMB by Dr. Éric Thiébaut"].map(function (c) {
                          return c.replace(/ /g, "\u00a0") }).join(" · ")
            }
            Button {
                Layout.alignment: Qt.AlignRight
                text: "Close"
                onClicked: aboutPanel.close()
            }
        }
    }

    Popup {
        id: settingsPanel
        // Centred. Anchored to the top-left it landed on the Observing settings column, which
        // is also a stack of labelled rows in a light panel -- close enough in appearance that
        // it was not obvious which one had focus. The middle of the window belongs to nothing
        // else, so a panel there reads as a dialog.
        //
        // Capped, and clamped, because it can outgrow the window. With no height set a Popup
        // takes its content's, and the content is font-sized: at a large enough UI font the
        // panel ran off the bottom and took the buttons with it -- the settings panel being
        // exactly where one goes to undo a setting like that. The content scrolls instead.
        //
        // Against `parent`, NOT `win`: a Popup is positioned inside the window's CONTENT item,
        // which is the window less its header. Measured against `win.height` the panel sat a
        // toolbar's height too low and ran off the bottom even when it would otherwise fit --
        // photographed under Xvfb at 26 pt before this line was written.
        //
        // The width grows with the FONT, not with the content: 860 at 11 pt, enough for a work
        // directory path, and the rows are text, so a larger font needs proportionally more room or the spin boxes get
        // squeezed. Deriving it from the content's implicit width instead is what Qt reports as
        // a binding loop -- the content's width comes back from the panel's, so asking the
        // content how wide the panel should be closes the circle.
        x: Math.max(dp(12), Math.round((parent.width  - width)  / 2))
        y: Math.max(dp(12), Math.round((parent.height - height) / 2))
        width:  Math.min(dp(860) * Math.max(1, baseFontPt / 11), parent.width - dp(24))
        height: Math.min(implicitHeight, parent.height - dp(24))
        // Dimmed behind, for the same reason: it separates the panel from whatever it covers.
        modal: true
        dim: true
        focus: true
        closePolicy: Popup.CloseOnEscape | Popup.CloseOnPressOutsideParent
        padding: dp(12)

        // The panel is capped at the window's height, so its content has to be able to scroll --
        // and to clip, since a control drawn outside the panel is worse than one below the fold.
        ScrollView {
            id: settingsScroll
            anchors.fill: parent
            clip: true
            rightPadding: dp(12)                   // room for the vertical scrollbar
            // The content is exactly as wide as the view: the rows are a label/control grid that
            // stretches, so there is nothing to scroll to sideways, and any binding that lets
            // the content decide the width closes a loop back through the panel's own width.
            contentWidth: availableWidth
            ScrollBar.horizontal.policy: ScrollBar.AlwaysOff

            ColumnLayout {
                id: settingsColumn
                width: settingsScroll.availableWidth
                spacing: dp(10)

                // Grouped by what a setting changes and when it takes effect: the window now,
                // the next launch, the plots, the Observing tab, and where files open.

                // ── Window ────────────────────────────────────────────────────
                SettingsHeader { title: "Window"; note: "applies now" }
                GridLayout {
                    columns: 2
                    columnSpacing: dp(8)
                    rowSpacing: dp(6)
                    Layout.fillWidth: true

                    Label { text: "UI scale"; color: "#666"; Layout.preferredWidth: dp(120) }
                    RowLayout {
                        spacing: dp(8)
                        // A multiple of the size computed for this screen, so 1.00 is the default
                        // on every screen. Everything sized through dp()/pt() rebinds, so the
                        // window resizes as the number changes.
                        SpinBox {
                            id: uiScaleSpin
                            objectName: "uiScaleSpin"
                            Layout.preferredWidth: dp(120)
                            from: 50; to: 300; stepSize: 5
                            value: Math.round(win.uiSize * 100)
                            editable: true
                            textFromValue: function (v) { return (v / 100).toFixed(2) }
                            valueFromText: function (t) { return Math.round(parseFloat(t.replace(",", ".")) * 100) }
                            validator: RegularExpressionValidator { regularExpression: /[0-9]*[.,]?[0-9]*×?/ }
                            onValueModified: win.uiScaleUser = value / 100
                        }
                        Button {
                            text: "auto"
                            enabled: win.uiScaleUser > 0
                            ToolTip.visible: hovered
                            ToolTip.text: uiScaleOverride > 0
                                ? "1.00 = OITOOLSGUI_SCALE=" + uiScaleOverride.toFixed(2)
                                : "1.00 = sized for " + deviceDpi.toFixed(0) + " dpi"
                            onClicked: { win.uiScaleUser = 0; uiScaleSpin.value = 100 }
                        }
                    }

                    Label { text: "Font"; color: "#666"; Layout.preferredWidth: dp(120) }
                    RowLayout {
                        spacing: dp(8)
                        ComboBox {
                            id: uiFontBox
                            Layout.preferredWidth: dp(240)
                            // Index 0 is "whatever the window would use if you had never set
                            // this", and it NAMES that font: which one it is depends on whether
                            // the shipped face loaded. When it did, the platform's own font goes
                            // at the end instead.
                            model: win.shippedFontFamily.length > 0
                                 ? ["Noto Sans", "DejaVu Sans", "Liberation Sans", "JuliaMono",
                                    "(system font)"]
                                 : ["(system default)", "DejaVu Sans", "Noto Sans",
                                    "Liberation Sans", "JuliaMono"]
                            onActivated: win.uiFontFamily =
                                currentIndex === 0             ? ""
                              : currentText === "(system font)" ? Qt.application.font.family
                                                                : currentText
                        }
                        Label { text: "size"; color: "#888"; font.pointSize: pt(baseFontPt - 2) }
                        SpinBox {
                            id: uiFontSizeSpin
                            Layout.preferredWidth: dp(120)
                            from: 6; to: 24; stepSize: 1
                            value: Math.round(baseFontPt)
                            editable: true
                            onValueModified: win.baseFontPt = value
                        }
                    }
                }

                // ── Start-up ──────────────────────────────────────────────────
                SettingsHeader { title: "Start-up"; note: "applies at next launch" }
                GridLayout {
                    columns: 2
                    columnSpacing: dp(8)
                    rowSpacing: dp(6)
                    Layout.fillWidth: true

                    Label { text: "Controls style"; color: "#666"; Layout.preferredWidth: dp(120) }
                    ComboBox {
                        id: controlsStyleBox
                        Layout.preferredWidth: dp(240)
                        // The list comes from Julia (CONTROLS_STYLES), so offering a style the
                        // bundled Qt does not carry cannot happen. Line 1 is the style the window
                        // is RUNNING; the rest are the choices.
                        property var styleInfo: Julia.shell_controls_styles().split("\n")
                        model: styleInfo.slice(1)
                        currentIndex: Math.max(0, model.indexOf(win.controlsStyle))
                        onActivated: win.controlsStyle = currentText
                        ToolTip.visible: hovered
                        ToolTip.text: "running: " + styleInfo[0]
                    }

                    // Linux only: elsewhere nothing reads it before Qt starts.
                    Label { text: "Linux UI scaler"; color: "#666"; Layout.preferredWidth: dp(120); visible: Qt.platform.os === "linux" }
                    RowLayout {
                        spacing: dp(8)
                        visible: Qt.platform.os === "linux"
                        SpinBox {
                            id: windowScaleSpin
                            objectName: "windowScaleSpin"
                            Layout.preferredWidth: dp(120)
                            from: 50; to: 400; stepSize: 25
                            value: Math.round(win.windowScaleTarget * 100)
                            editable: true
                            textFromValue: function (v) { return (v / 100).toFixed(2) }
                            valueFromText: function (t) { return Math.round(parseFloat(t.replace(",", ".")) * 100) }
                            validator: RegularExpressionValidator { regularExpression: /[0-9]*[.,]?[0-9]*×?/ }
                            onValueModified: win.windowScaleUser = value / 100
                            ToolTip.visible: hovered
                            ToolTip.text: "running: " + win.windowScaleRunning.toFixed(2)
                        }
                        Button {
                            text: "auto"
                            enabled: win.windowScaleUser > 0
                            ToolTip.visible: hovered
                            ToolTip.text: "auto: " + win.windowScaleAuto.toFixed(2)
                            onClicked: {
                                win.windowScaleUser = 0
                                windowScaleSpin.value = Math.round(win.windowScaleAuto * 100)
                            }
                        }
                    }

                    Label { text: "Open on"; color: "#666"; Layout.preferredWidth: dp(120) }
                    ComboBox {
                        id: startTabBox
                        objectName: "startTabBox"
                        Layout.preferredWidth: dp(240)
                        model: ["Exploring", "Observing", "Modeling", "Imaging", "Tools"]
                        currentIndex: Math.max(0, model.map(function (t) { return t.toLowerCase() })
                                                       .indexOf(win.startTab))
                        onActivated: win.startTab = currentText.toLowerCase()
                    }

                    Label { text: "Window"; color: "#666"; Layout.preferredWidth: dp(120) }
                    RowLayout {
                        spacing: dp(12)
                        CheckBox {
                            id: startFullscreenBox
                            objectName: "startFullscreenBox"
                            text: "full screen"
                            checked: win.startFullscreen
                            onToggled: win.startFullscreen = checked
                        }
                        CheckBox {
                            id: rememberGeometryBox
                            objectName: "rememberGeometryBox"
                            text: "remember size and position"
                            checked: win.rememberGeometry
                            onToggled: win.rememberGeometry = checked
                        }
                    }
                }

                // ── Plots ─────────────────────────────────────────────────────
                SettingsHeader { title: "Plots"; note: "applies now" }
                GridLayout {
                    columns: 2
                    columnSpacing: dp(8)
                    rowSpacing: dp(6)
                    Layout.fillWidth: true

                    Label { text: "Text scale"; color: "#666"; Layout.preferredWidth: dp(120) }
                    RowLayout {
                        spacing: dp(8)
                        SpinBox {
                            id: plotScaleSpin
                            objectName: "plotScaleSpin"
                            Layout.preferredWidth: dp(120)
                            from: 50; to: 500; stepSize: 5
                            value: 119                     // replaced on open by the live value
                            editable: true
                            textFromValue: function (v) { return (v / 100).toFixed(2) }
                            valueFromText: function (t) { return Math.round(parseFloat(t.replace(",", ".")) * 100) }
                            validator: RegularExpressionValidator { regularExpression: /[0-9]*[.,]?[0-9]*×?/ }
                            onValueModified: {
                                win.plotScaleUser = value / 100
                                Julia.shell_set_plot_scale(win.plotScaleUser)
                            }
                        }
                        Button {
                            text: "auto"
                            enabled: win.plotScaleUser > 0
                            ToolTip.visible: hovered
                            ToolTip.text: "from the screen's density"
                            onClicked: { Julia.shell_set_plot_scale(0); readPlotSettings() }
                        }
                    }

                    Label { text: "Symbols"; color: "#666"; Layout.preferredWidth: dp(120) }
                    RowLayout {
                        spacing: dp(8)
                        SpinBox {
                            id: markerSpin
                            Layout.preferredWidth: dp(120)
                            from: 0; to: 30; stepSize: 1
                            value: 0
                            editable: true
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
                    }

                    Label { text: "Wheel zoom"; color: "#666"; Layout.preferredWidth: dp(120) }
                    RowLayout {
                        spacing: dp(8)
                        SpinBox {
                            id: zoomSpin
                            Layout.preferredWidth: dp(120)
                            from: 102; to: 300; stepSize: 5
                            value: 125
                            editable: true
                            textFromValue: function (v) { return (v / 100).toFixed(2) + "×" }
                            valueFromText: function (t) { return Math.round(parseFloat(t.replace(",", ".")) * 100) }
                            validator: RegularExpressionValidator { regularExpression: /[0-9]*[.,]?[0-9]*×?/ }
                            onValueModified: {
                                win.zoomStepUser = value / 100
                                Julia.shell_set_zoom_step(win.zoomStepUser)
                            }
                        }
                        Button {
                            text: "auto"
                            enabled: win.zoomStepUser > 0
                            onClicked: { Julia.shell_set_zoom_step(0); readPlotSettings() }
                        }
                        Label { text: "per detent"; color: "#888"; font.pointSize: pt(baseFontPt - 2) }
                    }

                    // Shown, not settable: Makie takes the font when a Figure is built, and every
                    // route to change it afterwards throws out of its compute graph.
                    Label { text: "Font"; color: "#666"; Layout.preferredWidth: dp(120) }
                    Label { text: Julia.shell_plot_font() + "  (fixed)" }
                }

                // ── Observing ─────────────────────────────────────────────────
                SettingsHeader { title: "Observing"; note: "applies now" }
                GridLayout {
                    columns: 2
                    columnSpacing: dp(8)
                    rowSpacing: dp(6)
                    Layout.fillWidth: true

                    Label { text: "Time axis"; color: "#666"; Layout.preferredWidth: dp(120) }
                    ComboBox {
                        id: obsTimeBox
                        objectName: "obsTimeBox"
                        Layout.preferredWidth: dp(120)
                        model: observeTab.timeSystems
                        currentIndex: Math.max(0, observeTab.timeSystems.indexOf(observeTab.timeSystem))
                        onActivated: win.setObserveTimeSystem(currentText)
                    }

                    Label { text: "Elevation"; color: "#666"; Layout.preferredWidth: dp(120) }
                    RowLayout {
                        spacing: dp(8)
                        TextField {
                            id: obsElevMinField
                            objectName: "obsElevMin"
                            Layout.preferredWidth: dp(70)
                            text: observeTab.altLimit.toFixed(1)
                            validator: DoubleValidator { bottom: 0.0; top: 90.0; decimals: 1
                                                         notation: DoubleValidator.StandardNotation }
                            onEditingFinished: observeTab.altLimit = parseFloat(text)
                        }
                        Label { text: "to"; color: "#888"; font.pointSize: pt(baseFontPt - 2) }
                        TextField {
                            id: obsElevMaxField
                            objectName: "obsElevMax"
                            Layout.preferredWidth: dp(70)
                            text: observeTab.altMax.toFixed(1)
                            validator: DoubleValidator { bottom: 0.0; top: 90.0; decimals: 1
                                                         notation: DoubleValidator.StandardNotation }
                            onEditingFinished: observeTab.altMax = parseFloat(text)
                        }
                        Label { text: "°"; color: "#888" }
                    }

                    // No light/dark selector: the window pins a light palette (see the top of
                    // this file) because every panel, table and Makie canvas in it is drawn light.
                    // A dark mode is a real piece of work — the colour literals would have to
                    // become palette roles first — not a dropdown, so there is no dropdown.
                }

                // ── Files ─────────────────────────────────────────────────────
                SettingsHeader { title: "Files"; note: "where Open and Save start" }

                RowLayout {
                    Layout.fillWidth: true
                    spacing: dp(8)
                    Label { text: "Work directory"; color: "#666"; Layout.preferredWidth: dp(120) }
                    TextField {
                        id: workDirField
                        Layout.fillWidth: true
                        implicitHeight: dp(28)
                        text: win.workDirectory
                        placeholderText: "unset — the shipped demo data is used"
                        font.pointSize: pt(baseFontPt - 1)
                        ToolTip.visible: hovered && text.length > 0
                        ToolTip.text: text
                        // Typed as well as browsed: a path pasted from a terminal is often the
                        // fastest way in, and on a machine where the folder is not mounted yet
                        // it is the only way.
                        onEditingFinished: win.workDirectory = text.trim()
                    }
                    Button {
                        text: "Browse…"
                        implicitHeight: dp(28)
                        onClicked: workDirDialog.openAt(win.workDirectory)
                    }
                    Button {
                        text: "Clear"
                        implicitHeight: dp(28)
                        enabled: win.workDirectory.length > 0
                        ToolTip.visible: hovered
                        ToolTip.text: "forget it, and go back to opening in the demo data"
                        onClicked: win.workDirectory = ""
                    }
                }

                RowLayout {
                    Layout.fillWidth: true
                    // The buttons get the row to themselves. What happened is reported on the line
                    // BELOW: sharing the row, the message was elided to a middle-truncated path
                    // sitting beside the buttons, which reads as a stray fragment rather than as
                    // the answer to what the button just did — and it squeezed the buttons besides.
                    Item { Layout.fillWidth: true }
                    Button {
                        text: "Save config"
                        ToolTip.visible: hovered
                        ToolTip.text: "write these settings to the per-user config file, which " +
                                      "the window reads at every launch"
                        onClicked: {
                            // The OVERRIDES, not the values in force: a scale of 0 means "work it
                            // out from the screen", which is the right thing to carry to a machine
                            // with a different one.
                            var path = Julia.shell_save_settings(
                                [ "ui_size\t"     + win.uiScaleUser,
                                  "start_tab\t"   + win.startTab,
                                  "start_fullscreen\t"  + (win.startFullscreen ? 1 : 0),
                                  "remember_geometry\t" + (win.rememberGeometry ? 1 : 0),
                                  "obs_time_system\t"   + observeTab.timeSystem,
                                  "obs_elev_min\t"      + observeTab.altLimit,
                                  "obs_elev_max\t"      + observeTab.altMax,
                                  "window_scale\t" + win.windowScaleUser,
                                  "ui_font\t"     + win.uiFontFamily,
                                  "controls_style\t" + win.controlsStyle,
                                  "ui_font_pt\t"  + win.baseFontPt,
                                  "plot_scale\t"  + win.plotScaleUser,
                                  "marker_size\t" + win.markerSizeUser,
                                  "zoom_step\t"   + win.zoomStepUser,
                                  "work_directory\t" + win.workDirectory ].join("\n"))
                            savedLabel.text = path.length > 0 ? "saved · " + path
                                                              : "could not save — see the console"
                        }
                    }
                    Button {
                        text: "Reset to defaults"
                        ToolTip.visible: hovered
                        ToolTip.text: "back to the out-of-the-box appearance, and delete the saved config"
                        onClicked: {
                            // Both halves, and the file matters more than the window: settings are
                            // applied at startup, so restoring the look while leaving the file in
                            // place would come back tweaked at the next launch.
                            var d = win.appearanceDefaults
                            // The work directory goes with them: "reset to defaults" deletes the
                            // config file, so leaving the path on screen would show a setting that
                            // no longer exists anywhere.
                            win.workDirectory  = ""
                            win.uiScaleUser    = d.uiScaleUser
                            win.windowScaleUser = d.windowScaleUser
                            win.uiFontFamily   = d.uiFontFamily
                            win.baseFontPt     = d.baseFontPt
                            win.controlsStyle  = d.controlsStyle
                            uiFontBox.currentIndex        = 0
                            uiFontSizeSpin.value          = d.baseFontPt
                            win.startTab         = d.startTab
                            win.startFullscreen  = d.startFullscreen
                            win.rememberGeometry = d.rememberGeometry
                            win.setObserveTimeSystem(d.obsTimeSystem)
                            observeTab.altLimit  = d.obsElevMin
                            observeTab.altMax    = d.obsElevMax
                            uiScaleSpin.value             = 100
                            windowScaleSpin.value         = Math.round(win.windowScaleTarget * 100)
                            win.refreshSettingsBoxes()
                            controlsStyleBox.currentIndex = Math.max(0, controlsStyleBox.model.indexOf(d.controlsStyle))
                            // The plot side is Julia's, and zero means "work it out from the
                            // screen" rather than "zero" -- so the setters go first and the boxes
                            // are read back from what Julia then computed, not from d.
                            Julia.shell_set_plot_scale(d.plotScaleUser)
                            Julia.shell_set_marker_size(d.markerSizeUser)
                            Julia.shell_set_zoom_step(d.zoomStepUser)
                            readPlotSettings()
                            win.redrawPlot()
                            var removed = Julia.shell_reset_settings()
                            savedLabel.text = removed.length > 0 ? "reset · removed " + removed
                                                                 : "reset · nothing had been saved"
                        }
                    }
                    Button { text: "Close"; onClicked: settingsPanel.close() }
                }

                // Always exactly one line. Empty until a button was pressed, it made the panel
                // grow the first time anything was saved; wrapped, it grew again on a long
                // path. And it says what it IS, because a truncated path on its own reads as a
                // stray fragment rather than as the answer to what the button just did.
                Label {
                    id: savedLabel
                    Layout.fillWidth: true
                    elide: Text.ElideMiddle
                    maximumLineCount: 1
                    text: "config: " + Julia.shell_settings_path()
                    color: "#888"
                    font.pointSize: pt(baseFontPt - 2)
                    HoverHandler { id: savedHover }
                    ToolTip.visible: savedHover.hovered && savedLabel.truncated
                    ToolTip.text: savedLabel.text
                }

                // Which code is actually running. The first thing to establish about any bug
                // report, and the settings panel is where a user already comes to look at what
                // the window is doing rather than at their data.
                Label {
                    Layout.fillWidth: true
                    horizontalAlignment: Text.AlignRight
                    // The graphics API Qt Quick runs on, as Qt reports it.
                    readonly property string gfx: GraphicsInfo.api === GraphicsInfo.OpenGL
                        ? "OpenGL " + GraphicsInfo.majorVersion + "." + GraphicsInfo.minorVersion
                        : ""
                    text: Julia.shell_version() + "  ·  " + Julia.shell_diagnostics()
                          + (gfx.length > 0 ? "  ·  " + gfx : "")
                    color: "#888"
                    font.pointSize: pt(baseFontPt - 2)
                }
            }
        }
    }

    FilePicker {
        id: workDirDialog
        uiScale: win.uiScale; fontScale: win.fontScale; baseFontPt: win.baseFontPt
        title: "Choose work directory"
        chooseFolder: true
        onAccepted: function (path) { win.workDirectory = path }
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
            // A model fitted or an image reconstructed in another perspective is something
            // Explore can compare against, and a model removed is something it no longer can:
            // its compare boxes are re-read whenever it is shown.
            onCurrentIndexChanged: if (currentIndex === 0) exploreTab.refreshCompareAvailability()
            TabButton { text: "Exploring" }
            TabButton { text: "Observing" }
            TabButton { text: "Modeling" }
            TabButton { text: "Imaging" }
            TabButton { text: "Tools" }
        }

        StackLayout {
            id: stack
            Layout.fillWidth: true
            Layout.fillHeight: true
            currentIndex: tabs.currentIndex

            // ── Explore ──────────────────────────────────────────────────────
            ExploreTab {
                id: exploreTab
                uiScale:    win.uiScale
                fontScale:  win.fontScale
                baseFontPt: win.baseFontPt
                onStatusReported: function (text) { win.status = text }
                onActionDone: win.afterAction()
                onSavePngRequested: function (area) { win.savePng("explore", area) }
                onExportRequested: win.exportData()
                // A new dataset: the box lists and selects it, and everything follows it.
                onFilterApplied: win.selectCurrentDataset()
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
                // A simulated file becomes an ordinary dataset: the same load path a file
                // chosen by hand takes, then the Exploring tab, because the point of the
                // button is to LOOK at what was written.
                onOpenRequested: function (path) {
                    win.status = Julia.shell_open(path)
                    win.refreshDatasets()
                    win.afterAction()
                    tabs.currentIndex = 0
                }
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

            // ── Tools ────────────────────────────────────────────────────────
            // Splits and averages make new datasets, which the window lists and selects.
            ToolsTab {
                id: toolsTab
                uiScale:    win.uiScale
                fontScale:  win.fontScale
                baseFontPt: win.baseFontPt
                onConsoleChanged: win.refreshConsole()
                onStatusReported: function (text) { win.status = text }
                onDatasetsChanged: win.selectCurrentDataset()
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
                        onClicked: saveDialog.openAt(Julia.shell_initial_folder())
                    }
                }
                ScrollView {
                    id: logScroll
                    visible: logToggle.checked
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    TextArea {
                        id: logArea
                        objectName: "consoleText"
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

    // Window walks (`OITOOLSGUI_WALK`, see `shell_walk_request`): a request names a control,
    // the reply says where it is in window coordinates. Inert when `walkDir` is empty.
    //
    //   locate q   q = objectName | text:T (a visible item whose text is T) | popup:T (the same,
    //              inside an open popup) | a/b (b searched inside a); one match, or an error
    //   reveal q   the same, after scrolling it into view
    //   dump       every visible named item
    function walkKids(o) {
        const k = []
        if (o.modal !== undefined && o.contentItem !== undefined && o.opened !== undefined) {
            if (!o.opened) return k
            for (const x of [o.contentItem, o.header, o.footer]) if (x) k.push(x)
            return k
        }
        if (o.data && typeof o.data !== "function" && o.data.length !== undefined)
            for (let i = 0; i < o.data.length; i++) k.push(o.data[i])
        if (o.popup && typeof o.popup === "object") k.push(o.popup)
        return k
    }
    function walkSearch(root, test, all) {
        const hits = [], seen = new Set()
        const go = (o, inPopup) => {
            if (!o || seen.has(o)) return
            seen.add(o)
            if (o.visible === false) return
            const pop = inPopup || (o.modal !== undefined && o.opened === true)
            if (o.mapToItem && test(o, pop)) hits.push(o)
            for (const c of walkKids(o)) go(c, pop)
        }
        go(root, false)
        // a Label and a Button may share a text: the one that can be clicked is meant
        const live = hits.filter(o => typeof o.clicked === "function" || typeof o.toggled === "function")
        return all || live.length === 0 ? hits : live
    }
    function walkFind(q) {
        let root = win
        const parts = q.split("/")
        for (let i = 0; i < parts.length; i++) {
            const p = parts[i]
            const test = p.startsWith("text:") ? ((o, pop) => o.text === p.slice(5))
                       : p.startsWith("popup:") ? ((o, pop) => pop && o.text === p.slice(6))
                       : ((o, pop) => o.objectName === p)
            const hits = walkSearch(root, test)
            if (hits.length !== 1)
                throw new Error(hits.length + " matches for " + p +
                                hits.map(o => " [" + String(o) + "]").join(""))
            root = hits[0]
        }
        return root
    }
    // One tab-separated line: name x y w h enabled checked current text (newlines as spaces).
    // A control with no text (a SpinBox) reports its value there.
    function walkDescribe(o) {
        const pt = o.mapToItem(null, 0, 0)
        const flat = v => (v === undefined || v === null) ? "" : String(v).replace(/[\t\n]/g, " ")
        return [o.objectName, Math.round(pt.x), Math.round(pt.y), Math.round(o.width),
                Math.round(o.height), o.enabled, flat(o.checked), flat(o.currentText),
                flat(o.text !== undefined ? o.text : o.value)].join("\t")
    }
    // Scroll every Flickable around `o` (a ScrollView's, a ListView) until `o` is in view.
    function walkReveal(o) {
        for (let p = o.parent; p; p = p.parent) {
            if (p.contentY === undefined || p.contentItem === undefined || !(p.height > 0)) continue
            const r = o.mapToItem(p.contentItem, 0, 0)
            const most = Math.max(0, p.contentHeight - p.height)
            if (r.y < p.contentY)
                p.contentY = Math.max(0, r.y - 10)
            else if (r.y + o.height > p.contentY + p.height)
                p.contentY = Math.min(most, r.y + o.height - p.height + 10)
        }
        return o
    }
    function walkAnswer(req) {
        const sp = req.indexOf(" ")
        const cmd = sp < 0 ? req : req.slice(0, sp), arg = sp < 0 ? "" : req.slice(sp + 1)
        if (cmd === "locate") return walkDescribe(walkFind(arg))
        if (cmd === "reveal") return walkDescribe(walkReveal(walkFind(arg)))
        if (cmd === "dump")
            return walkSearch(win, (o, pop) => o.objectName !== "" && o.objectName !== "label", true).map(walkDescribe).join("\n")
        throw new Error("unknown walk request: " + cmd)
    }
    Timer {
        interval: 150
        running: walkDir !== ""
        repeat: true
        onTriggered: {
            const req = Julia.shell_walk_request()
            if (req === "") return
            let out
            try { out = win.walkAnswer(req) } catch (e) { out = "!\t" + String(e) }
            Julia.shell_walk_reply(out)
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
