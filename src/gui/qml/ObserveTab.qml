// Observe — planning a night before it happens.
//
// This is the only perspective that PRODUCES a dataset rather than consuming one, so it ends in
// a "Simulate" action and an explicit hand-off to Explore rather than in the shared dataset
// selector. The context bar, task tray and command log belong to the window, not to this tab,
// and are deliberately not repeated here.
//
// CHARA only, by decision rather than by omission. The delay-line and POP machinery in plan.jl
// is CHARA's unless the caller supplies a pop_array and an airpath, and the VLTI facility TOMLs
// carry no delay_lengths at all, so a VLTI feasibility answer here would be confidently wrong.
// One array whose answers are true beats two where half of them lie.
//
// QML holds no state that Julia owns. Every list below is a placeholder standing in for a real
// call, and each seam carries a `// wire:` comment naming the function that will fill it.

import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs
import jlqml
import Makie            // MakieArea, for this perspective's own Gantt canvas

Item {
    id: root

    // ── scale contract, identical to Main.qml's ───────────────────────────────
    //
    // The window owns the detection and passes the result down, so a tab can be laid out at a
    // different density for testing without touching the shell.
    property real uiScale:   1.0
    property real fontScale: 1.0
    property real baseFontPt: 11
    function dp(px)     { return Math.round(px * uiScale) }   // pixel lengths
    function pt(points) { return points * fontScale }         // explicit point sizes

    // ── instrument state ──────────────────────────────────────────────────────
    property var    facilities: ["CHARA"]
    property string facility: "CHARA"
    property var    telescopeNames: []
    // 0 = unused, 1 = used, 2 = reference — the encoding get_baselines() reads, so this vector
    // goes to best_pop / in_delay / simulate untranslated.
    property var    telescopeConfig: []
    property string referenceTelescope: ""
    property var    combiners: []
    property string combiner: ""
    property var    spectralSetups: []
    property string spectralSetup: ""

    // ── target state ──────────────────────────────────────────────────────────
    //
    // targetModel is what the rows edit; `targets` is the flat snapshot Julia reads, so the list
    // stays editable in QML without Julia having to walk a ListModel.
    property var targets: []
    property int currentTargetIndex: 0
    property int cacheCount: 0

    // ── conditions ────────────────────────────────────────────────────────────
    //
    // seeing / r0 / t0 are mutable fields on FacilityConfig; the AO guide magnitude is not —
    // AOConfig is immutable, so changing it means constructing a new one and reassigning
    // facility.ao.
    property real seeing: 1.0
    property real r0:     0.078
    property real t0:     0.0032
    property real magAO:  2.0

    // ── time ──────────────────────────────────────────────────────────────────
    property string dateISO: "2026-08-19"
    property real   haMin: -4.0
    property real   haMax:  4.0
    property int    stepMinutes: 10
    // ASPRO's defaults. 25 rather than 30 because the last five degrees are a real part of a
    // short night; 85 rather than 90 because near the zenith the delay lines and the azimuth
    // axis both have to track fast, and the cap genuinely bites -- Vega transits at 85.4 from
    // CHARA.
    property real   altLimit:   25.0
    property real   altMax:     85.0
    property real   moonMinSep: 30.0
    property real   darkOffset:  0.0
    // Dark-window readout. The two hours are decimal local time and are properties rather than
    // text because the Gantt axis does arithmetic on them.
    property real   darkStartHour: 20.0
    property real   darkEndHour:    5.0
    property string darkWindowText: "—"

    // ── views ─────────────────────────────────────────────────────────────────
    // Two charts and the POP table, and nothing else. The Moon is reported as a number on the
    // plan line rather than given a view of its own, and uv coverage and SNR belong to a
    // dataset — so the way to see them is to Simulate and open the result in Exploring, which
    // already draws both, rather than to build a second uv plot here.
    property var viewNames: ["Gantt", "Delay carts", "POPs"]
    // Readonly and driven by the TabBar: a ComboBox or TabBar drops a two-way currentIndex
    // binding the instant the user clicks, so the control is the source of truth and this is a
    // read-only mirror of it.
    readonly property int currentView: viewTabs.currentIndex
    // ── simulation output ─────────────────────────────────────────────────────
    property string outputFile: ""
    property string lastSimulatedFile: ""
    property real   simMag:      2.0
    property real   simMagAO:    2.0
    property bool   simNoise:    true
    property bool   simDebias:   true
    property int    simNSamples: 100
    // Seeded by default. An unseeded run cannot be reproduced by the exported script, which is
    // the one guarantee the command log makes.
    property int    simSeed: 1
    property bool   simObservability: true

    // ── what the simulation is OF ────────────────────────────────────────────
    //
    // `simulate` takes a sky, and the two kinds are genuinely different arguments rather than
    // one argument with a flag: an image (2-D grey, or 3-D one plane per channel) goes in as
    // `image` with a `pixsize`, a model goes in as `flat_model` plus its current `flat_params`.
    // So the panel has to know which it holds before it can say Simulate is possible.
    property string simSource: "image"          // image | cube | model
    property string simImagePath: ""
    property string simModelPath: ""            // empty ⇒ the Modeling perspective's model
    property real   simPixsize: 0.1             // mas, images only; a model carries its own scale
    property string simSourceText: "no file chosen"
    property bool   simSourceOk: false

    function refreshSimSource() {
        var path = simSource === "model" ? simModelPath : simImagePath
        var r = Julia.shell_sim_source(simSource, path)
        var f = r.split("\t")
        simSourceOk   = f[0] === "ok"
        simSourceText = f.length > 1 ? f[1] : r
    }
    onSimSourceChanged:    refreshSimSource()
    onSimImagePathChanged: refreshSimSource()
    onSimModelPathChanged: refreshSimSource()

    property string statusText: "ready"

    // ── derived state the controls enable on ──────────────────────────────────
    readonly property int nSelectedTelescopes: {
        var n = 0
        for (var i = 0; i < telescopeConfig.length; ++i) if (telescopeConfig[i] > 0) n += 1
        return n
    }
    readonly property var selectedTelescopes: {
        var out = []
        for (var i = 0; i < telescopeConfig.length; ++i)
            if (telescopeConfig[i] > 0) out.push(telescopeNames[i])
        return out
    }
    readonly property int nNamedTargets: {
        var n = 0
        for (var i = 0; i < targets.length; ++i)
            if (targets[i].name && targets[i].name.length > 0) n += 1
        return n
    }
    readonly property int nBaselines:
        referenceTelescope.length > 0 ? Math.max(0, nSelectedTelescopes - 1)
                                      : nSelectedTelescopes * (nSelectedTelescopes - 1) / 2
    // One string rather than a set of booleans: the reason a run is impossible has to be
    // readable next to the button, not deduced from which control happens to be grey.
    readonly property string blockReason:
        nNamedTargets === 0        ? "no target" :
        nSelectedTelescopes < 2    ? "select at least 2 telescopes" :
        combiner.length === 0      ? "no combiner" :
        spectralSetup.length === 0 ? "no spectral setup" :
        // A simulation is OF something. Without a usable sky there is nothing to observe, and
        // that is as much a blocker as having no combiner.
        !simSourceOk               ? "no source: " + simSourceText :
        outputFile.length === 0    ? "no output file" : ""
    readonly property bool canSimulate: blockReason.length === 0

    // ── placeholder catalogue ─────────────────────────────────────────────────
    //
    // wire: list_configs() -> (facilities, combiners, wavelengths, targets, wave_combiner).
    // wave_combiner is the only machine-readable link from a spectral setup back to its combiner,
    // so it is what setupsFor() will consult instead of matching on filename.
    readonly property var telescopeCatalog: ({
        "CHARA": ["S1", "S2", "E1", "E2", "W1", "W2"]
    })
    readonly property var combinerCatalog: ({
        "CHARA": ["MIRCX", "MYSTIC", "SPICA"]
    })
    readonly property var setupCatalog: ({
        "MIRCX":  ["MIRCX_LOWH", "MIRCX_LOWJ"],
        "MYSTIC": ["MYSTIC_LOWK"],
        "SPICA":  ["SPICA_LR"]
    })

    function telescopesFor(f) { return telescopeCatalog[f] !== undefined ? telescopeCatalog[f] : [] }
    function combinersFor(f)  { return combinerCatalog[f]  !== undefined ? combinerCatalog[f]  : [] }
    function setupsFor(c)     { return setupCatalog[c]     !== undefined ? setupCatalog[c]     : [] }

    // The three instrument controls cascade because their contents are not independent: a
    // combiner belongs to one facility and a spectral setup belongs to one combiner. Letting a
    // stale pair survive — MYSTIC_LOWK under SPICA — would produce a simulation nobody could
    // have observed, so every level rebuilds the ones below it.
    function selectFacility(f) {
        root.facility = f
        root.telescopeNames = telescopesFor(f)
        var cfg = []
        for (var i = 0; i < root.telescopeNames.length; ++i) cfg.push(1)
        root.telescopeConfig = cfg
        root.referenceTelescope = ""
        root.combiners = combinersFor(f)
        selectCombiner(root.combiners.length > 0 ? root.combiners[0] : "")
    }

    function selectCombiner(c) {
        root.combiner = c
        root.spectralSetups = setupsFor(c)
        root.spectralSetup = root.spectralSetups.length > 0 ? root.spectralSetups[0] : ""
    }

    function setTelescopeUsed(i, used) {
        var cfg = root.telescopeConfig.slice()
        cfg[i] = used ? 1 : 0
        root.telescopeConfig = cfg
        // A reference that is no longer in the array is not a reference.
        if (!used && root.telescopeNames[i] === root.referenceTelescope) setReference("")
    }

    // A reference telescope reduces the baselines to reference-to-each instead of every pair
    // (get_baselines, astrometry.jl). That is a real loss of uv coverage, so it is an explicit
    // choice and defaults to none.
    function setReference(name) {
        var cfg = root.telescopeConfig.slice()
        var ok = false
        for (var i = 0; i < cfg.length; ++i) {
            if (cfg[i] === 2) cfg[i] = 1
            if (root.telescopeNames[i] === name && cfg[i] > 0) { cfg[i] = 2; ok = true }
        }
        root.telescopeConfig = cfg
        root.referenceTelescope = ok ? name : ""
    }

    function configString() {
        var parts = []
        for (var i = 0; i < root.telescopeConfig.length; ++i)
            if (root.telescopeConfig[i] > 0)
                parts.push(root.telescopeNames[i] + (root.telescopeConfig[i] === 2 ? "*" : ""))
        return parts.length > 0 ? parts.join(" ") : "none"
    }

    // ── target list plumbing ──────────────────────────────────────────────────
    ListModel {
        id: targetModel
        // Two northern seeds so the panel is never empty on a first run.
        ListElement { name: "Vega";   ra: 279.234735; dec: 38.783689; magV: 0.03; magJ: -0.18; magH: -0.03; magK: 0.13; cached: false }
        ListElement { name: "Altair"; ra: 297.695827; dec:  8.868322; magV: 0.76; magJ: 0.35; magH: 0.24; magK: 0.24; cached: false }
    }

    function refreshTargets() {
        var out = []
        var nc = 0
        for (var i = 0; i < targetModel.count; ++i) {
            var t = targetModel.get(i)
            out.push({ name: t.name, ra: t.ra, dec: t.dec,
                       magV: t.magV, magJ: t.magJ, magH: t.magH, magK: t.magK,
                       cached: t.cached })
            if (t.cached) nc += 1
        }
        root.targets = out
        root.cacheCount = nc
    }

    function addTarget() {
        targetModel.append({ name: "", ra: 0.0, dec: 0.0,
                             magV: 0.0, magJ: 0.0, magH: 0.0, magK: 0.0, cached: false })
        root.currentTargetIndex = targetModel.count - 1
        refreshTargets()
    }

    function removeTarget(i) {
        if (targetModel.count <= 1) return          // the Gantt needs at least one row to draw
        targetModel.remove(i)
        if (root.currentTargetIndex >= targetModel.count) root.currentTargetIndex = targetModel.count - 1
        refreshTargets()
    }

    function setTargetField(i, key, value) {
        var patch = {}
        patch[key] = value
        // Editing name, RA or Dec by hand makes the row no longer what SIMBAD returned, so the
        // cache flag drops — otherwise the indicator would vouch for hand-typed coordinates.
        if (key === "name" || key === "ra" || key === "dec") patch["cached"] = false
        targetModel.set(i, patch)
        refreshTargets()
    }

    // Applies one shell_simbad reply to one row. A magnitude SIMBAD does not publish comes
    // back as an empty field, and the row keeps whatever was there rather than being zeroed:
    // magnitude 0 is a real and very bright value, so writing it for "unknown" would put a
    // first-magnitude star in every empty band.
    function applySimbad(i, reply) {
        var f = String(reply).split("\t")
        if (f[0] !== "ok") {
            root.statusText = "SIMBAD: " + (f.length > 1 ? f[1] : "failed")
            return false
        }
        var patch = { ra: parseFloat(f[1]), dec: parseFloat(f[2]), cached: true }
        var bands = ["magV", "magJ", "magH", "magK"]
        for (var b = 0; b < bands.length; ++b)
            if (f[3 + b] !== "" && f[3 + b] !== undefined) patch[bands[b]] = parseFloat(f[3 + b])
        targetModel.set(i, patch)
        refreshTargets()
        return true
    }

    // The query is a network round trip on Qt's thread, so the status line is set and the call
    // deferred by a tick — otherwise the window repaints only after the answer arrives and the
    // user sees nothing at all until it does.
    function resolveTarget(i) {
        var nm = targetModel.get(i).name
        if (!nm || nm.length === 0) return
        root.statusText = "resolving " + nm + "…"
        resolveTimer.queue = [i]
        resolveTimer.restart()
    }

    function resolveAll() {
        var q = []
        for (var i = 0; i < targetModel.count; ++i) {
            var t = targetModel.get(i)
            // Rows already cached are left alone: they were resolved this session and asking
            // again would spend a round trip to be told the same thing.
            if (t.name && t.name.length > 0 && !t.cached) q.push(i)
        }
        if (q.length === 0) { root.statusText = "nothing to resolve"; return }
        root.statusText = "resolving " + q.length + " target" + (q.length === 1 ? "" : "s") + "…"
        resolveTimer.queue = q
        resolveTimer.restart()
    }

    // One row per tick, so a list of ten targets repaints between queries instead of freezing
    // for the whole batch.
    Timer {
        id: resolveTimer
        property var queue: []
        property int failures: 0
        interval: 20
        repeat: false
        onTriggered: {
            if (queue.length === 0) return
            var i = queue.shift()
            if (!root.applySimbad(i, Julia.shell_simbad(targetModel.get(i).name))) failures += 1
            root.consoleChanged()       // shell_simbad logs what it resolved; show it
            if (queue.length > 0) {
                restart()
            } else if (failures > 0) {
                failures = 0
            } else {
                root.statusText = "resolved"
            }
        }
    }

    function clearCache() {
        for (var i = 0; i < targetModel.count; ++i) targetModel.set(i, { cached: false })
        refreshTargets()
        // Clears only the per-row flags. There is no Julia-side cache to drop: shell_simbad
        // queries SIMBAD every time, so a cleared row genuinely re-resolves.
        root.statusText = "SIMBAD flags cleared"
    }

    Component.onCompleted: { selectFacility(root.facility); refreshTargets() }

    // ── Gantt time axis arithmetic ────────────────────────────────────────────
    //
    // The night wraps past midnight, so the span is a modular difference and not a subtraction.
    readonly property real ganttSpanHours: {
        var s = darkEndHour - darkStartHour
        while (s <= 0) s += 24
        return s
    }
    readonly property int ganttTickCount:
        Math.max(2, Math.floor(darkStartHour + ganttSpanHours) - Math.ceil(darkStartHour) + 1)
    function ganttTickHour(k) { return (Math.ceil(darkStartHour) + k) % 24 }
    function ganttTickFrac(k) { return (Math.ceil(darkStartHour) + k - darkStartHour) / ganttSpanHours }

    // ── layout ────────────────────────────────────────────────────────────────
    ColumnLayout {
        anchors.fill: parent
        anchors.margins: dp(6)
        spacing: dp(6)

        RowLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            spacing: dp(8)

            // ── settings column ───────────────────────────────────────────────
            // ── settings, and the action they feed ───────────────────────────
            //
            // The panels scroll; Compute and the POPs it uses do not. An action button that
            // scrolls out of view is one the user has to hunt for, and this one sits where
            // Run sits on the Imaging tab.
            ColumnLayout {
                // One width for the column and everything in it. Two different numbers left a
                // blank strip: the panels sized to the ScrollView while the row below sized to
                // its own content, so the column ended up wider than the panels it held and the
                // difference showed as dead space against the night panel.
                Layout.preferredWidth: dp(440)
                Layout.maximumWidth: dp(440)
                Layout.fillHeight: true
                spacing: dp(6)

                ScrollView {
                    id: settingsScroll
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    contentWidth: availableWidth

                    ColumnLayout {
                        width: settingsScroll.availableWidth
                        spacing: dp(8)

                        // ── targets ───────────────────────────────────────────────
                        Rectangle {
                            Layout.fillWidth: true
                            implicitHeight: targetPanel.implicitHeight + dp(12)
                            color: "#fbfbfb"
                            border.color: "#ddd"

                            ColumnLayout {
                                id: targetPanel
                                anchors.fill: parent
                                anchors.margins: dp(6)
                                spacing: dp(6)

                                RowLayout {
                                    Layout.fillWidth: true
                                    Label { text: "Targets"; font.bold: true }
                                    Item { Layout.fillWidth: true }
                                    Label {
                                        // The cache is per session and worth showing: a resolve
                                        // answered from cache looks exactly like one that reached the
                                        // network, right up until the network is down.
                                        text: root.cacheCount + "/" + targetModel.count + " cached"
                                        color: "#888"; font.pointSize: pt(baseFontPt - 2)
                                    }
                                    Button {
                                        text: "Clear cache"
                                        enabled: root.cacheCount > 0
                                        onClicked: root.clearCache()
                                    }
                                }

                                Repeater {
                                    id: targetRepeater
                                    model: targetModel

                                    Rectangle {
                                        id: targetRow
                                        required property int index
                                        required property string name
                                        required property real ra
                                        required property real dec
                                        required property real magV
                                        required property real magJ
                                        required property real magH
                                        required property real magK
                                        required property bool cached

                                        Layout.fillWidth: true
                                        implicitHeight: targetRowColumn.implicitHeight + dp(10)
                                        color: index === root.currentTargetIndex ? "#eef4fb" : "#ffffff"
                                        border.color: index === root.currentTargetIndex ? "#9bbfe0" : "#e4e4e4"

                                        MouseArea {
                                            // Selecting a row must not eat the click that lands in one
                                            // of its editors, so the press is passed on rather than
                                            // consumed.
                                            anchors.fill: parent
                                            onPressed: function (mouse) {
                                                root.currentTargetIndex = targetRow.index
                                                mouse.accepted = false
                                            }
                                        }

                                        ColumnLayout {
                                            id: targetRowColumn
                                            anchors.fill: parent
                                            anchors.margins: dp(5)
                                            spacing: dp(4)

                                            RowLayout {
                                                Layout.fillWidth: true
                                                spacing: dp(4)
                                                Rectangle {
                                                    // Per-row cache indicator: filled means these
                                                    // numbers came from a resolve this session, hollow
                                                    // means they were typed or have been edited since.
                                                    Layout.alignment: Qt.AlignVCenter
                                                    implicitWidth: dp(9); implicitHeight: dp(9)
                                                    radius: implicitWidth / 2
                                                    color: targetRow.cached ? "#59a14f" : "transparent"
                                                    border.color: targetRow.cached ? "#59a14f" : "#bbb"
                                                }
                                                TextField {
                                                    Layout.fillWidth: true
                                                    text: targetRow.name
                                                    placeholderText: "target name"
                                                    onEditingFinished: root.setTargetField(targetRow.index, "name", text)
                                                }
                                                Button {
                                                    text: "SIMBAD"
                                                    enabled: targetRow.name.length > 0
                                                    onClicked: root.resolveTarget(targetRow.index)
                                                }
                                                Button {
                                                    text: "−"
                                                    enabled: targetModel.count > 1
                                                    onClicked: root.removeTarget(targetRow.index)
                                                }
                                            }

                                            RowLayout {
                                                Layout.fillWidth: true
                                                spacing: dp(4)
                                                // Decimal degrees throughout, matching
                                                // TargetConfig.raep0 and night_observability.
                                                // Sexagesimal is where the sign of a "-00 30 00"
                                                // declination goes missing.
                                                Label { text: "RA°"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                                TextField {
                                                    Layout.fillWidth: true
                                                    text: targetRow.ra.toFixed(6)
                                                    validator: DoubleValidator {
                                                        bottom: 0.0; top: 360.0; decimals: 6
                                                        notation: DoubleValidator.StandardNotation
                                                    }
                                                    onEditingFinished: root.setTargetField(targetRow.index, "ra", parseFloat(text))
                                                }
                                                Label { text: "Dec°"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                                TextField {
                                                    Layout.fillWidth: true
                                                    text: targetRow.dec.toFixed(6)
                                                    validator: DoubleValidator {
                                                        bottom: -90.0; top: 90.0; decimals: 6
                                                        notation: DoubleValidator.StandardNotation
                                                    }
                                                    onEditingFinished: root.setTargetField(targetRow.index, "dec", parseFloat(text))
                                                }
                                            }

                                            RowLayout {
                                                Layout.fillWidth: true
                                                spacing: dp(4)
                                                Label { text: "V"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                                TextField {
                                                    Layout.fillWidth: true
                                                    text: targetRow.magV.toFixed(2)
                                                    validator: DoubleValidator { decimals: 3; notation: DoubleValidator.StandardNotation }
                                                    onEditingFinished: root.setTargetField(targetRow.index, "magV", parseFloat(text))
                                                }
                                                Label { text: "J"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                                TextField {
                                                    Layout.fillWidth: true
                                                    text: targetRow.magJ.toFixed(2)
                                                    validator: DoubleValidator { decimals: 3; notation: DoubleValidator.StandardNotation }
                                                    onEditingFinished: root.setTargetField(targetRow.index, "magJ", parseFloat(text))
                                                }
                                                Label { text: "H"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                                TextField {
                                                    Layout.fillWidth: true
                                                    text: targetRow.magH.toFixed(2)
                                                    validator: DoubleValidator { decimals: 3; notation: DoubleValidator.StandardNotation }
                                                    onEditingFinished: root.setTargetField(targetRow.index, "magH", parseFloat(text))
                                                }
                                                Label { text: "K"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                                TextField {
                                                    Layout.fillWidth: true
                                                    text: targetRow.magK.toFixed(2)
                                                    validator: DoubleValidator { decimals: 3; notation: DoubleValidator.StandardNotation }
                                                    onEditingFinished: root.setTargetField(targetRow.index, "magK", parseFloat(text))
                                                }
                                            }
                                        }
                                    }
                                }

                                RowLayout {
                                    Layout.fillWidth: true
                                    Button { text: "Add target"; onClicked: root.addTarget() }
                                    Button {
                                        text: "Resolve all"
                                        enabled: root.nNamedTargets > 0
                                        onClicked: root.resolveAll()
                                    }
                                    Item { Layout.fillWidth: true }
                                }
                            }
                        }

                        // ── instrument ────────────────────────────────────────────
                        Rectangle {
                            Layout.fillWidth: true
                            implicitHeight: instrumentPanel.implicitHeight + dp(12)
                            color: "#fbfbfb"
                            border.color: "#ddd"

                            ColumnLayout {
                                id: instrumentPanel
                                anchors.fill: parent
                                anchors.margins: dp(6)
                                spacing: dp(6)

                                RowLayout {
                                    Layout.fillWidth: true
                                    Label { text: "Instrument"; font.bold: true }
                                    Item { Layout.fillWidth: true }
                                    Label { text: "CHARA only"; color: "#888"; font.pointSize: pt(baseFontPt - 2) }
                                }

                                GridLayout {
                                    Layout.fillWidth: true
                                    columns: 2
                                    columnSpacing: dp(6)
                                    rowSpacing: dp(4)

                                    Label { text: "Facility"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                    ComboBox {
                                        id: facilityBox
                                        Layout.fillWidth: true
                                        model: root.facilities
                                        // One entry today, but still a ComboBox: the cascade below is
                                        // written against a facility choice, and hiding the control
                                        // would hide what everything else is derived from.
                                        enabled: root.facilities.length > 1
                                        onActivated: root.selectFacility(currentText)
                                    }

                                    Label { text: "Combiner"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                    ComboBox {
                                        id: combinerBox
                                        Layout.fillWidth: true
                                        model: root.combiners
                                        enabled: root.combiners.length > 0
                                        currentIndex: Math.max(0, root.combiners.indexOf(root.combiner))
                                        // A ComboBox drops its currentIndex binding the moment the user
                                        // picks something, so after the cascade replaces the model the
                                        // index is re-asserted explicitly rather than left to a binding
                                        // that may no longer exist.
                                        onModelChanged: currentIndex = Math.max(0, root.combiners.indexOf(root.combiner))
                                        onActivated: root.selectCombiner(currentText)
                                    }

                                    Label { text: "Spectral setup"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                    ComboBox {
                                        id: setupBox
                                        Layout.fillWidth: true
                                        model: root.spectralSetups
                                        // Filtered by combiner, not merely sorted: the `combiner` key
                                        // inside each wavelength TOML is the only machine-readable link
                                        // between the two, and a setup belonging to another combiner
                                        // would be accepted by simulate() and quietly produce nonsense.
                                        enabled: root.spectralSetups.length > 0
                                        currentIndex: Math.max(0, root.spectralSetups.indexOf(root.spectralSetup))
                                        onModelChanged: currentIndex = Math.max(0, root.spectralSetups.indexOf(root.spectralSetup))
                                        onActivated: root.spectralSetup = currentText
                                    }

                                    Label { text: "Reference"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                    ComboBox {
                                        id: referenceBox
                                        Layout.fillWidth: true
                                        // Choosing one reduces the baselines to reference-to-each
                                        // instead of every pair, so "(none)" is the default and the
                                        // reduction is always deliberate.
                                        model: ["(none)"].concat(root.selectedTelescopes)
                                        enabled: root.nSelectedTelescopes > 1
                                        currentIndex: Math.max(0, model.indexOf(root.referenceTelescope))
                                        onModelChanged: currentIndex = Math.max(0, model.indexOf(root.referenceTelescope))
                                        onActivated: root.setReference(currentIndex === 0 ? "" : currentText)
                                    }
                                }

                                Label { text: "Telescopes"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                Flow {
                                    Layout.fillWidth: true
                                    spacing: dp(8)
                                    Repeater {
                                        id: telescopeRepeater
                                        model: root.telescopeNames
                                        CheckBox {
                                            required property int index
                                            required property string modelData
                                            text: modelData + (root.telescopeConfig[index] === 2 ? " ★" : "")
                                            checked: root.telescopeConfig[index] > 0
                                            onToggled: root.setTelescopeUsed(index, checked)
                                        }
                                    }
                                }
                                Label {
                                    Layout.fillWidth: true
                                    text: root.nSelectedTelescopes < 2
                                          ? "select at least 2 telescopes"
                                          : root.configString() + "  ·  " + root.nBaselines + " baselines"
                                    color: root.nSelectedTelescopes < 2 ? "#b00000" : "#666"
                                    font.pointSize: pt(baseFontPt - 2)
                                    elide: Text.ElideRight
                                }
                            }
                        }

                        // ── POPs, by hand ────────────────────────────────────────────────
                        //
                        // Beside the Instrument panel that determines them: a POP is a beam
                        // path through the array, so it belongs with the telescopes it routes
                        // rather than with the button that consumes it.
                        //
                        // Search POPs writes here too, so the searched configuration and a
                        // hand-dialled one are the same field rather than two competing
                        // sources of truth.
                        Rectangle {
                            Layout.fillWidth: true
                            implicitHeight: popCol.implicitHeight + dp(12)
                            color: "#f6f6f6"
                            border.color: "#dcdcdc"
                            radius: dp(3)

                            ColumnLayout {
                                id: popCol
                                x: dp(6); y: dp(6)
                                width: parent.width - dp(12)
                                spacing: dp(4)

                                RowLayout {
                                    Layout.fillWidth: true
                                    Label { text: "POPs"; font.bold: true }
                                    Item { Layout.fillWidth: true }
                                    CheckBox {
                                        text: "AutoPOPs"
                                        checked: root.autoPops
                                        onToggled: root.autoPops = checked
                                        ToolTip.visible: hovered
                                        ToolTip.text: "Search for the best POPs each time the " +
                                                      "plan is computed, and adopt them below. " +
                                                      "Untick to use the numbers as set here."
                                    }
                                }

                                // One dropdown per telescope, labelled with the telescope. CHARA has five
                                // POPs; the number IS the beam path, so a free text box invites 0 and 9.
                                Flow {
                                    Layout.fillWidth: true
                                    spacing: dp(8)
                                    Repeater {
                                        model: root.telescopeNames
                                        delegate: RowLayout {
                                            required property int index
                                            required property string modelData
                                            spacing: dp(3)
                                            Label {
                                                text: modelData
                                                font.bold: true
                                                color: root.telescopeConfig[index] > 0 ? "#222" : "#aaa"
                                            }
                                            ComboBox {
                                                implicitWidth: dp(58)
                                                implicitHeight: dp(26)
                                                model: ["1", "2", "3", "4", "5"]
                                                // Output, not input, while AutoPOPs is on: the
                                                // next Compute would overwrite anything set here.
                                                enabled: root.telescopeConfig[index] > 0 && !root.autoPops
                                                currentIndex: Math.max(0, root.popAt(index) - 1)
                                                onActivated: root.setPop(index, currentIndex + 1)
                                            }
                                        }
                                    }
                                }

                                // The same numbers as text, so a row can be pasted in from elsewhere or
                                // copied out — including from the POPs table, whose Use button writes here.
                                RowLayout {
                                    Layout.fillWidth: true
                                    spacing: dp(6)
                                    TextField {
                                        id: popField
                                        Layout.fillWidth: true
                                        enabled: !root.autoPops
                                        placeholderText: root.autoPops
                                            ? "set by AutoPOPs"
                                            : "paste six POP numbers, e.g. 1 3 3 4 1 1"
                                        text: root.popString
                                        font.family: "monospace"
                                        onEditingFinished: root.setPopString(text)
                                    }
                                    Button {
                                        text: "Copy"
                                        enabled: root.popString.length > 0
                                        onClicked: { popField.selectAll(); popField.copy(); popField.deselect() }
                                    }
                                    Button {
                                        text: "Clear"
                                        enabled: root.popString.length > 0
                                        onClicked: root.popString = ""
                                    }
                                }
                            }
                        }

                        // ── conditions ────────────────────────────────────────────
                        Rectangle {
                            Layout.fillWidth: true
                            implicitHeight: conditionsPanel.implicitHeight + dp(12)
                            color: "#fbfbfb"
                            border.color: "#ddd"

                            ColumnLayout {
                                id: conditionsPanel
                                anchors.fill: parent
                                anchors.margins: dp(6)
                                spacing: dp(6)

                                Label { text: "Conditions"; font.bold: true }

                                GridLayout {
                                    Layout.fillWidth: true
                                    // One pair per row. At four columns the two labels and their two fields
                                    // set a minimum width the panel could not meet, so it pushed
                                    // into the Gantt beside it. Taller is free here; wider is not.
                                    columns: 2
                                    columnSpacing: dp(6)
                                    rowSpacing: dp(4)

                                    Label { text: "Seeing (\")"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                    TextField {
                                        id: seeingField
                                        Layout.fillWidth: true
                                        text: root.seeing.toFixed(2)
                                        validator: DoubleValidator {
                                            bottom: 0.1; top: 5.0; decimals: 2
                                            notation: DoubleValidator.StandardNotation
                                        }
                                        // wire: facility.seeing = value (mutable field on FacilityConfig)
                                        onEditingFinished: root.seeing = parseFloat(text)
                                    }
                                    Label { text: "r0 (m)"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                    TextField {
                                        id: r0Field
                                        Layout.fillWidth: true
                                        text: root.r0.toFixed(4)
                                        validator: DoubleValidator {
                                            bottom: 0.001; top: 2.0; decimals: 4
                                            notation: DoubleValidator.StandardNotation
                                        }
                                        // wire: facility.r0 = value. seeing and r0 describe the same
                                        // atmosphere, so whichever of the two Julia derives from the
                                        // other has to be written back here or the readouts disagree.
                                        onEditingFinished: root.r0 = parseFloat(text)
                                    }

                                    Label { text: "t0 (s)"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                    TextField {
                                        id: t0Field
                                        Layout.fillWidth: true
                                        text: root.t0.toFixed(5)
                                        validator: DoubleValidator {
                                            bottom: 0.0001; top: 1.0; decimals: 6
                                            notation: DoubleValidator.StandardNotation
                                        }
                                        // wire: facility.t0 = value
                                        onEditingFinished: root.t0 = parseFloat(text)
                                    }
                                    Label { text: "AO guide mag"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                    TextField {
                                        id: magAOField
                                        Layout.fillWidth: true
                                        text: root.magAO.toFixed(2)
                                        validator: DoubleValidator {
                                            bottom: -5.0; top: 20.0; decimals: 2
                                            notation: DoubleValidator.StandardNotation
                                        }
                                        // wire: AOConfig is immutable, so this constructs a new one and
                                        // reassigns facility.ao rather than assigning a field.
                                        onEditingFinished: root.magAO = parseFloat(text)
                                    }
                                }
                            }
                        }

                        // ── time ──────────────────────────────────────────────────
                        Rectangle {
                            Layout.fillWidth: true
                            implicitHeight: timePanel.implicitHeight + dp(12)
                            color: "#fbfbfb"
                            border.color: "#ddd"

                            ColumnLayout {
                                id: timePanel
                                anchors.fill: parent
                                anchors.margins: dp(6)
                                spacing: dp(6)

                                Label { text: "Time"; font.bold: true }

                                RowLayout {
                                    Layout.fillWidth: true
                                    spacing: dp(4)
                                    Label { text: "Date"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                    TextField {
                                        id: dateField
                                        Layout.preferredWidth: dp(120)
                                        text: root.dateISO
                                        inputMask: "9999-99-99"
                                        // wire: recompute the dark window and every view for this date
                                        onEditingFinished: root.dateISO = text
                                    }
                                    // Stepping is done in Julia: adding a month to 31 January has to
                                    // land on 28 February, which is calendar arithmetic and not
                                    // something to attempt on a string here.
                                    //
                                    // Double arrows move a whole month, which is the step that
                                    // matters for observability — a target's window shifts by about
                                    // two hours of LST per month, so a month is the unit in which
                                    // "is it up yet this season" gets answered.
                                    Button {
                                        text: "◀◀"
                                        // Sized to the glyph. A default-width button here is
                                        // mostly padding, and five of them pushed the panel
                                        // wider than the column that holds it.
                                        implicitWidth: dp(30)
                                        ToolTip.visible: hovered; ToolTip.text: "one month earlier"
                                        onClicked: root.shiftDate(0, -1)
                                    }
                                    Button {
                                        text: "◀"
                                        // Sized to the glyph. A default-width button here is
                                        // mostly padding, and five of them pushed the panel
                                        // wider than the column that holds it.
                                        implicitWidth: dp(30)
                                        ToolTip.visible: hovered; ToolTip.text: "previous night"
                                        onClicked: root.shiftDate(-1, 0)
                                    }
                                    // Sized to its own label. Only the arrows needed shrinking;
                                    // this one carries a word.
                                    Button { text: "Today"; onClicked: root.shiftDate(0, 0) }
                                    Button {
                                        text: "▶"
                                        // Sized to the glyph. A default-width button here is
                                        // mostly padding, and five of them pushed the panel
                                        // wider than the column that holds it.
                                        implicitWidth: dp(30)
                                        ToolTip.visible: hovered; ToolTip.text: "next night"
                                        onClicked: root.shiftDate(1, 0)
                                    }
                                    Button {
                                        text: "▶▶"
                                        // Sized to the glyph. A default-width button here is
                                        // mostly padding, and five of them pushed the panel
                                        // wider than the column that holds it.
                                        implicitWidth: dp(30)
                                        ToolTip.visible: hovered; ToolTip.text: "one month later"
                                        onClicked: root.shiftDate(0, 1)
                                    }
                                    Item { Layout.fillWidth: true }
                                }

                                RowLayout {
                                    Layout.fillWidth: true
                                    spacing: dp(4)
                                    Label { text: "Dark window"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                    Label {
                                        Layout.fillWidth: true
                                        // wire: sunrise_sunset(date, lat, lon) plus dark_offset — the
                                        // same window night_observability uses. Reads "—" until it has
                                        // been computed rather than showing a plausible default nobody
                                        // asked for.
                                        text: root.darkWindowText
                                        font.family: "monospace"
                                        elide: Text.ElideRight
                                        ToolTip.visible: hovered
                                        ToolTip.text: "Astronomical night: the hours between " +
                                            "evening and morning twilight, which is the window " +
                                            "observability is computed in and the grey band " +
                                            "behind the Gantt."
                                        HoverHandler { id: dwHover }
                                        property bool hovered: dwHover.hovered
                                    }
                                }

                                GridLayout {
                                    Layout.fillWidth: true
                                    columns: 4
                                    columnSpacing: dp(6)
                                    rowSpacing: dp(4)

                                    Label { text: "HA min (h)"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                    TextField {
                                        id: haMinField
                                        Layout.fillWidth: true
                                        text: root.haMin.toFixed(2)
                                        validator: DoubleValidator {
                                            bottom: -12.0; top: 12.0; decimals: 2
                                            notation: DoubleValidator.StandardNotation
                                        }
                                        onEditingFinished: root.haMin = parseFloat(text)
                                    }
                                    Label { text: "HA max (h)"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                    TextField {
                                        id: haMaxField
                                        Layout.fillWidth: true
                                        text: root.haMax.toFixed(2)
                                        validator: DoubleValidator {
                                            bottom: -12.0; top: 12.0; decimals: 2
                                            notation: DoubleValidator.StandardNotation
                                        }
                                        onEditingFinished: root.haMax = parseFloat(text)
                                    }

                                    Label { text: "Step (min)"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                    SpinBox {
                                        id: stepBox
                                        Layout.fillWidth: true
                                        from: 1; to: 120; value: root.stepMinutes
                                        onValueModified: root.stepMinutes = value
                                    }
                                    Label {
                                        Layout.columnSpan: 2
                                        Layout.fillWidth: true
                                        // best_pop scores a count of steps, so it only reads as minutes
                                        // at a 1-minute step. Saying so beside the control is cheaper
                                        // than explaining a POP table whose units quietly changed.
                                        text: root.stepMinutes === 1
                                              ? "POP score is in minutes"
                                              : "POP score counts steps of " + root.stepMinutes + " min"
                                        color: "#888"; font.pointSize: pt(baseFontPt - 2)
                                        elide: Text.ElideRight
                                    }

                                    Label { text: "Elev min (°)"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                    TextField {
                                        id: altLimitField
                                        Layout.fillWidth: true
                                        text: root.altLimit.toFixed(1)
                                        validator: DoubleValidator {
                                            bottom: 0.0; top: 90.0; decimals: 1
                                            notation: DoubleValidator.StandardNotation
                                        }
                                        onEditingFinished: root.altLimit = parseFloat(text)
                                    }
                                    Label { text: "Elev max (°)"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                    TextField {
                                        id: altMaxField
                                        Layout.fillWidth: true
                                        text: root.altMax.toFixed(1)
                                        validator: DoubleValidator {
                                            bottom: 0.0; top: 90.0; decimals: 1
                                            notation: DoubleValidator.StandardNotation
                                        }
                                        onEditingFinished: root.altMax = parseFloat(text)
                                    }

                                    Label { text: "Moon sep (°)"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                    TextField {
                                        id: moonSepField
                                        Layout.fillWidth: true
                                        text: root.moonMinSep.toFixed(1)
                                        validator: DoubleValidator {
                                            bottom: 0.0; top: 180.0; decimals: 1
                                            notation: DoubleValidator.StandardNotation
                                        }
                                        onEditingFinished: root.moonMinSep = parseFloat(text)
                                    }
                                    Label { text: "Dark offset (h)"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                                    TextField {
                                        id: darkOffsetField
                                        Layout.fillWidth: true
                                        text: root.darkOffset.toFixed(2)
                                        validator: DoubleValidator {
                                            bottom: 0.0; top: 6.0; decimals: 2
                                            notation: DoubleValidator.StandardNotation
                                        }
                                        onEditingFinished: root.darkOffset = parseFloat(text)
                                    }
                                }
                            }
                        }
                    }
                }

                        // ── compute ──────────────────────────────────────────────────
                        //
                        // Explicit, not automatic on every edit. A night is cheap to compute once the
                        // code is warm (0.1 ms), but recomputing on each keystroke while a coordinate
                        // is being typed would draw a chart for every half-entered number.
                        RowLayout {
                            Layout.fillWidth: true
                            spacing: dp(8)

                            Button {
                                text: "Compute"
                                enabled: root.nNamedTargets > 0 && root.dateISO.length > 0
                                ToolTip.visible: hovered && !enabled
                                ToolTip.text: "needs a named target and a date"
                                onClicked: root.computePlan()
                            }
                            CheckBox {
                                text: "Show details"
                                enabled: root.useDelay
                                checked: root.detailed
                                onToggled: { root.detailed = checked; if (root.hasPlan) root.computePlan() }
                                ToolTip.visible: hovered
                                ToolTip.text: "ASPRO's detailed view: one row per BASELINE, so the " +
                                              "baseline that closes the night is the short bar. Needs " +
                                              "the delay check on. Off shows only the answer."
                            }
                            // Errors only. The summary that used to sit here was longer than the
                            // strip and elided to nothing useful; it is in the console, and the
                            // chart answers the question better than a sentence about it does.
                            Label {
                                Layout.fillWidth: true
                                visible: root.planText.indexOf("!") === 0
                                text: root.planText
                                elide: Text.ElideRight
                                color: "#c62828"
                            }
                        }

            }


            // ── views ─────────────────────────────────────────────────────────
            ColumnLayout {
                Layout.fillWidth: true
                Layout.fillHeight: true
                spacing: dp(4)


                TabBar {
                    id: viewTabs
                    Layout.fillWidth: true
                    // Gantt first, because it is the view this perspective exists for.
                    currentIndex: 0
                    Repeater {
                        model: root.viewNames
                        TabButton {
                            required property string modelData
                            text: modelData
                        }
                    }
                }

                StackLayout {
                    id: viewStack
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    currentIndex: viewTabs.currentIndex

                    // ── elevation / observability ─────────────────────────────
                    ColumnLayout {
                        spacing: dp(4)

                        RowLayout {
                            Layout.fillWidth: true
                            spacing: dp(10)
                            Label { text: "Night of " + root.dateISO; font.bold: true }
                            Label {
                                text: root.darkWindowText
                                color: "#888"; font.pointSize: pt(baseFontPt - 2)
                            }
                            Item { Layout.fillWidth: true }

                            // Legend. The colours are the contract between this component and
                            // whatever Makie draws into the mount below.
                            RowLayout {
                                spacing: dp(4)
                                Rectangle { implicitWidth: dp(10); implicitHeight: dp(10); radius: dp(2); color: "#4c78a8"; border.color: "#999" }
                                Label { text: "observable"; color: "#666"; font.pointSize: pt(baseFontPt - 3) }
                            }
                            RowLayout {
                                spacing: dp(4)
                                Rectangle { implicitWidth: dp(10); implicitHeight: dp(10); radius: dp(2); color: "#59a14f"; border.color: "#999" }
                                Label { text: "in delay"; color: "#666"; font.pointSize: pt(baseFontPt - 3) }
                            }
                            RowLayout {
                                spacing: dp(4)
                                Rectangle { implicitWidth: dp(10); implicitHeight: dp(10); radius: dp(2); color: "#f2c14e"; border.color: "#999" }
                                Label { text: "twilight"; color: "#666"; font.pointSize: pt(baseFontPt - 3) }
                            }
                            RowLayout {
                                spacing: dp(4)
                                Rectangle { implicitWidth: dp(10); implicitHeight: dp(10); radius: dp(2); color: "#e15759"; border.color: "#999" }
                                Label {
                                    text: "moon < " + root.moonMinSep.toFixed(0) + "°"
                                    color: "#666"; font.pointSize: pt(baseFontPt - 3)
                                }
                            }
                            RowLayout {
                                spacing: dp(4)
                                Rectangle { implicitWidth: dp(10); implicitHeight: dp(10); radius: dp(2); color: "#cccccc"; border.color: "#999" }
                                Label { text: "below limit"; color: "#666"; font.pointSize: pt(baseFontPt - 3) }
                            }
                        }

                        // Time axis, running the full width of the mount below it.
                        RowLayout {
                            Layout.fillWidth: true
                            spacing: 0
                            Item {
                                id: ganttAxis
                                Layout.fillWidth: true
                                Layout.preferredHeight: dp(20)
                                Repeater {
                                    model: root.ganttTickCount
                                    Item {
                                        required property int index
                                        x: ganttAxis.width * root.ganttTickFrac(index)
                                        y: 0
                                        width: 0
                                        height: ganttAxis.height
                                        Label {
                                            anchors.horizontalCenter: parent.left
                                            anchors.top: parent.top
                                            text: root.ganttTickHour(index) + "h"
                                            color: "#666"; font.pointSize: pt(baseFontPt - 3)
                                        }
                                        Rectangle {
                                            anchors.bottom: parent.bottom
                                            anchors.horizontalCenter: parent.left
                                            width: dp(1); height: dp(4)
                                            color: "#bbb"
                                        }
                                    }
                                }
                            }
                        }

                        RowLayout {
                            Layout.fillWidth: true
                            Layout.fillHeight: true
                            spacing: 0

                            // No row-label column: the chart names its own rows on the y axis,
                            // so a second list of the same names beside it said nothing and took
                            // width the night could have used. Targets are selected in the
                            // Targets panel, which is where the rest of their settings are.
                            Rectangle {
                                id: ganttMount
                                // This perspective's own figure, `ganttPlot`, built beside the
                                // other two before the window exists. The chart is static:
                                // a Gantt is read, not manipulated, so there is no picking and
                                // no second input path through MakieArea.
                                Layout.fillWidth: true
                                Layout.fillHeight: true
                                color: "#ffffff"
                                border.color: "#ddd"

                                MakieArea {
                                    id: ganttArea
                                    anchors.fill: parent
                                    scene: ganttPlot

                                    // What the chart says under the pointer. Julia reads the
                                    // position from Makie's own event stream -- the route
                                    // picking already uses -- so there is no second coordinate
                                    // system to keep in step with this one.
                                    HoverHandler { id: ganttHover }
                                    ToolTip.visible: ganttHover.hovered &&
                                                     root.hoverText.length > 0
                                    ToolTip.text: root.hoverText
                                    ToolTip.delay: 250

                                    // Polled rather than pushed: Makie has no QML-visible signal
                                    // for a mouse move, and 120 ms is imperceptible for a
                                    // readout while costing nothing when the pointer is still.
                                    Timer {
                                        interval: 120
                                        repeat: true
                                        running: ganttHover.hovered && root.hasPlan
                                        onTriggered: {
                                            var t = Julia.shell_gantt_hover()
                                            if (t.length > 0) t += root.magnitudeLine(t)
                                            if (t !== root.hoverText) root.hoverText = t
                                        }
                                    }
                                }

                                // Over the canvas until a night has been computed: an empty
                                // axis reads as a failed plan rather than as one not yet run.
                                Rectangle {
                                    anchors.fill: parent
                                    visible: !root.hasPlan
                                    color: "#ffffff"
                                    Label {
                                        anchors.centerIn: parent
                                        horizontalAlignment: Text.AlignHCenter
                                        color: "#888"
                                        text: "choose a target and a date, then Compute"
                                    }
                                }
                            }
                        }
                    }

                    // ── delay carts (the `chara_plan` view) ──────────────────
                    //
                    // Where the Gantt reduces each baseline to in-or-out, this keeps the
                    // number: a baseline running at 44 m of a 45.7 m limit is feasible and
                    // about to stop being so, and only the cart position says that.
                    Rectangle {
                        color: "#ffffff"
                        border.color: "#ddd"

                        MakieArea {
                            id: delayArea
                            anchors.fill: parent
                            scene: delayPlot
                        }

                        Rectangle {
                            anchors.fill: parent
                            visible: !root.hasPlan
                            color: "#ffffff"
                            Label {
                                anchors.centerIn: parent
                                horizontalAlignment: Text.AlignHCenter
                                color: "#888"
                                text: root.hasPlan
                                      ? "choose a target and a date, then Compute"
                                      : "choose a target and a date, then Compute"
                            }
                        }
                    }


                    ColumnLayout {
                        spacing: dp(4)

                        RowLayout {
                            Layout.fillWidth: true
                            Layout.margins: dp(4)
                            spacing: dp(8)
                            Button {
                                text: "Search POPs"
                                enabled: root.nSelectedTelescopes >= 2 && root.nNamedTargets > 0
                                // CHARA_POP_ARRAY and CHARA_AIRPATH are the defaults inside
                                // best_pop, which is the other half of why this tab is CHARA-only.
                                // The best result is adopted rather than merely listed: leaving
                                // the user to copy it is how the searched POPs get forgotten and
                                // the delay check reports a third of the real window.
                                onClicked: { root.findPops(); root.statusText =
                                             popModel.count > 0 ? "POPs " + root.popString
                                                                : "no POP configuration found" }
                            }
                            Label { text: "n_best"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                            SpinBox { id: nBestBox; from: 1; to: 20; value: 5 }
                            Label { text: "min minutes"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                            SpinBox { id: minMinutesBox; from: 0; to: 600; value: 10 }
                            Item { Layout.fillWidth: true }
                            Label {
                                text: "config: " + root.configString()
                                color: "#888"; font.pointSize: pt(baseFontPt - 2)
                                elide: Text.ElideRight
                            }
                        }

                        Rectangle {
                            Layout.fillWidth: true
                            Layout.preferredHeight: dp(150)
                            color: "#fafafa"
                            border.color: "#ddd"
                            // wire: compute_delays(facility, dec, ha, config) and in_delay(...) —
                            // OPD per baseline against the per-telescope delay_length, drawn by
                            // Makie.
                            Label {
                                anchors.centerIn: parent
                                horizontalAlignment: Text.AlignHCenter
                                color: "#666"
                                text: "delay vs hour angle\ncompute_delays() / in_delay()"
                                font.family: "monospace"
                            }
                        }

                        // POP results are a table and not a plot because the answer is a
                        // configuration to dial in, and a configuration has to be read off exactly.
                        Rectangle {
                            Layout.fillWidth: true
                            Layout.fillHeight: true
                            color: "#ffffff"
                            border.color: "#ddd"

                            ColumnLayout {
                                anchors.fill: parent
                                anchors.margins: dp(4)
                                spacing: dp(2)

                                RowLayout {
                                    Layout.fillWidth: true
                                    spacing: dp(6)
                                    Label { Layout.preferredWidth: dp(40);  text: "#";        font.bold: true; font.pointSize: pt(baseFontPt - 2) }
                                    // The telescopes, in the order the POP numbers below are
                                    // written. "POPs" alone leaves the reader counting columns
                                    // to work out which number belongs to which telescope, and
                                    // the whole point of the table is dialling those in exactly.
                                    Label { Layout.fillWidth: true
                                            text: root.telescopeNames.length > 0
                                                  ? "POPs   (" + root.telescopeNames.join("  ") + ")"
                                                  : "POPs"
                                            font.bold: true; font.pointSize: pt(baseFontPt - 2) }
                                    Label { Layout.preferredWidth: dp(90);  text: "score";    font.bold: true; font.pointSize: pt(baseFontPt - 2) }
                                    Label { Layout.preferredWidth: dp(130); text: "HA range"; font.bold: true; font.pointSize: pt(baseFontPt - 2) }
                                }
                                Rectangle { Layout.fillWidth: true; Layout.preferredHeight: dp(1); color: "#ddd" }

                                ListView {
                                    id: popList
                                    Layout.fillWidth: true
                                    Layout.fillHeight: true
                                    clip: true
                                    model: popModel
                                    delegate: RowLayout {
                                        id: popRow
                                        required property int index
                                        required property string pops
                                        required property string score
                                        required property string haRange
                                        width: popList.width
                                        // A ListView sizes delegates horizontally only, and a
                                        // RowLayout outside a layout parent keeps height 0, so
                                        // the rows exist but draw nothing without this.
                                        height: implicitHeight
                                        spacing: dp(6)

                                        Label { Layout.preferredWidth: dp(40);  text: popRow.index + 1; font.family: "monospace" }
                                        Label { Layout.fillWidth: true;         text: popRow.pops;      font.family: "monospace"; elide: Text.ElideRight
                                                color: popRow.pops === root.popString ? "#1a7f37" : "#222"
                                                font.bold: popRow.pops === root.popString }
                                        Label { Layout.preferredWidth: dp(90);  text: popRow.score;     font.family: "monospace" }
                                        Label { Layout.preferredWidth: dp(130); text: popRow.haRange;   font.family: "monospace" }
                                        Button {
                                            text: popRow.pops === root.popString ? "in use" : "Use"
                                            enabled: popRow.pops !== root.popString
                                            implicitHeight: dp(20)
                                            ToolTip.visible: hovered
                                            ToolTip.text: "use these POPs, and stop AutoPOPs from replacing them"
                                            onClicked: root.adoptPops(popRow.pops)
                                        }
                                    }
                                }

                                Label {
                                    Layout.fillWidth: true
                                    visible: popModel.count === 0
                                    text: "no POP search run yet"
                                    color: "#aaa"
                                    horizontalAlignment: Text.AlignHCenter
                                }
                            }
                        }
                    }
                }
            }
        }

        // ── output ────────────────────────────────────────────────────────────
        Rectangle {
            Layout.fillWidth: true
            implicitHeight: outputColumn.implicitHeight + dp(12)
            color: "#f4f4f4"
            border.color: "#ddd"

            ColumnLayout {
                id: outputColumn
                anchors.fill: parent
                anchors.margins: dp(6)
                spacing: dp(6)

                RowLayout {
                    Layout.fillWidth: true
                    spacing: dp(8)

                    Label { text: "mag"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                    TextField {
                        id: simMagField
                        Layout.preferredWidth: dp(70)
                        text: root.simMag.toFixed(2)
                        validator: DoubleValidator {
                            bottom: -5.0; top: 20.0; decimals: 2
                            notation: DoubleValidator.StandardNotation
                        }
                        onEditingFinished: root.simMag = parseFloat(text)
                    }
                    Label { text: "mag_ao"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                    TextField {
                        id: simMagAoField
                        Layout.preferredWidth: dp(70)
                        text: root.simMagAO.toFixed(2)
                        validator: DoubleValidator {
                            bottom: -5.0; top: 20.0; decimals: 2
                            notation: DoubleValidator.StandardNotation
                        }
                        onEditingFinished: root.simMagAO = parseFloat(text)
                    }
                    CheckBox {
                        id: noiseBox
                        text: "noise"
                        checked: root.simNoise
                        onToggled: root.simNoise = checked
                    }
                    CheckBox {
                        id: debiasBox
                        text: "debias"
                        // Debiasing only means something once there is noise to debias, and
                        // n_samples only sizes the noise realisation, so both follow the box above.
                        enabled: root.simNoise
                        checked: root.simDebias && root.simNoise
                        onToggled: root.simDebias = checked
                    }
                    Label { text: "n_samples"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                    SpinBox {
                        id: nSamplesBox
                        Layout.preferredWidth: dp(120)
                        from: 1; to: 100000; value: root.simNSamples
                        enabled: root.simNoise
                        onValueModified: root.simNSamples = value
                    }
                    Label { text: "seed"; color: "#666"; font.pointSize: pt(baseFontPt - 2) }
                    SpinBox {
                        id: seedBox
                        Layout.preferredWidth: dp(120)
                        from: 0; to: 999999; value: root.simSeed
                        onValueModified: root.simSeed = value
                    }
                    CheckBox {
                        id: observabilityBox
                        text: "observability filter"
                        checked: root.simObservability
                        // Opt-in inside simulate() as well: with it off every epoch is used whether
                        // or not the target was up, which is the right default for a pure uv
                        // simulator and the wrong one for a night plan.
                        onToggled: root.simObservability = checked
                    }
                    Item { Layout.fillWidth: true }
                }

                // ── source ───────────────────────────────────────────────────
                RowLayout {
                    Layout.fillWidth: true
                    spacing: dp(8)

                    Label { text: "Source"; font.bold: true }
                    ComboBox {
                        Layout.preferredWidth: dp(140)
                        model: ["Image", "Image cube", "Model"]
                        currentIndex: root.simSource === "cube" ? 1
                                    : root.simSource === "model" ? 2 : 0
                        onActivated: root.simSource = ["image", "cube", "model"][currentIndex]
                    }

                    Button {
                        text: root.simSource === "model" ? "Model file…" : "FITS…"
                        onClicked: root.simSource === "model" ? simModelDialog.openAt("")
                                                             : simImageDialog.openAt("")
                    }
                    Button {
                        visible: root.simSource === "model"
                        text: "From Modeling tab"
                        // The model being edited next door, at its current values — no file, no
                        // round trip. `simulate` wants `flat_model` and `flat_params`, which is
                        // exactly what that perspective already holds.
                        onClicked: { root.simModelPath = ""; root.refreshSimSource() }
                    }

                    Label { text: "pixel size"; visible: root.simSource !== "model"; color: "#666" }
                    // A plain TextField, as everywhere else in this file: `NumField` is an
                    // inline component of ImageTab and inline components are private to the
                    // file that declares them.
                    TextField {
                        id: simPixField
                        visible: root.simSource !== "model"
                        Layout.preferredWidth: dp(80)
                        text: root.simPixsize.toFixed(4)
                        horizontalAlignment: Text.AlignRight
                        onEditingFinished: {
                            var v = parseFloat(text)
                            if (!isNaN(v) && v > 0) root.simPixsize = v
                            else text = root.simPixsize.toFixed(4)
                        }
                    }
                    Label { text: "mas"; visible: root.simSource !== "model"; color: "#666"
                            font.pointSize: pt(baseFontPt - 2) }

                    Label {
                        Layout.fillWidth: true
                        text: root.simSourceText
                        color: root.simSourceOk ? "#2e7d32" : "#b00000"
                        elide: Text.ElideMiddle
                        font.pointSize: pt(baseFontPt - 1)
                    }
                }

                RowLayout {
                    Layout.fillWidth: true
                    spacing: dp(8)

                    Button {
                        text: "Output file…"
                        onClicked: outputDialog.openAt("")
                    }
                    Label {
                        Layout.fillWidth: true
                        text: root.outputFile.length > 0 ? root.outputFile : "no output file chosen"
                        color: root.outputFile.length > 0 ? "#000" : "#888"
                        font.family: "monospace"
                        elide: Text.ElideMiddle
                    }
                    Label {
                        visible: !root.canSimulate
                        text: root.blockReason
                        color: "#b00000"
                        font.pointSize: pt(baseFontPt - 2)
                    }
                    Button {
                        id: simulateButton
                        text: "Simulate"
                        enabled: root.canSimulate
                        onClicked: simulateTimer.restart()
                    }
                    Button {
                        id: openInExploreButton
                        text: "Open in Explore"
                        // Cross-tab moves are explicit verbs, never implicit state changes, so this
                        // stays a button and stays disabled until there is a file to carry across.
                        enabled: root.lastSimulatedFile.length > 0
                        // wire: hand lastSimulatedFile to the shell's open path and switch the
                        // window to the Explore tab.
                        onClicked: root.statusText = "open " + root.lastSimulatedFile
                    }
                }

                Label {
                    Layout.fillWidth: true
                    text: root.statusText
                    color: "#666"
                    font.pointSize: pt(baseFontPt - 2)
                    elide: Text.ElideRight
                }
            }
        }
    }

    // POP search results. Empty until a search runs; the first append fixes the roles, which are
    // pops, score and haRange.
    FilePicker {
        id: simImageDialog
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Image or cube for the simulation"
        filters: [{ label: "FITS files (*.fits *.fit)", patterns: "*.fits,*.fit" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) { root.simImagePath = path }
    }

    FilePicker {
        id: simModelDialog
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Model file for the simulation"
        filters: [{ label: "Model files (*.toml)", patterns: "*.toml" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) { root.simModelPath = path }
    }

    ListModel { id: popModel }

    FilePicker {
        id: outputDialog
        uiScale: root.uiScale; fontScale: root.fontScale; baseFontPt: root.baseFontPt
        title: "Simulated OIFITS"
        saveMode: true
        defaultSuffix: "oifits"
        filters: [{ label: "OIFITS files (*.oifits *.fits)", patterns: "*.oifits,*.fits" },
                  { label: "All files", patterns: "*" }]
        onAccepted: function (path) { root.outputFile = path }
    }

    // The same deferred idiom Main.qml uses for its file dialogs: simulate() runs synchronously on
    // the GUI thread and Qt cannot repaint while a handler is running, so the work is pushed one
    // tick out and the interface gets to show that it started.
    Timer {
        id: simulateTimer
        interval: 250; repeat: false
        onTriggered: {
            // wire: simulate(facility, target, combiner, wavelength, dates, out_file; mag, mag_ao,
            // noise, debias, n_samples, seed, observability) — dates come from the HA range and
            // step, and observability is the named-tuple form when the filter box is ticked. Then
            // set root.lastSimulatedFile and register the file in the session as a new dataset.
            root.statusText = "simulate -> " + root.outputFile
        }
    }

    // ── computing a night ────────────────────────────────────────────────────
    //
    // The delay-line check is OPT-IN. Applying it with an unsearched POP configuration reports
    // far less time than is really available -- measured on Vega from CHARA, 0 minutes with all
    // six telescopes at POP 1, against 380 with the POPs `best_pop` finds -- so a panel that
    // applied it by default would call every target unobservable.
    property bool hasPlan: false
    // Always on. The delay lines are what makes a Gantt an OBSERVABILITY chart rather than a
    // horizon chart: a target above the horizon that no cart can reach is not observable, and
    // showing it as though it were is the one answer this panel must not give. Measured at
    // 0.1 ms on top of a 0.1 ms plan, so there is nothing to opt out of.
    property bool useDelay: true

    // AutoPOPs: search for the POPs whenever the plan is REGENERATED, not when it is redrawn.
    // Measured at 0.5 ms for four telescopes and 11.6 ms for six, which is inside a frame --
    // so the beam paths can simply always be the best ones, the way ASPRO treats them, rather
    // than something the user has to remember to search for and then apply.
    property bool autoPops: true
    // Summary by default, as in ASPRO: "when can I observe this" is the question asked first.
    // Detailed is the delay view -- one row per baseline -- and answers "why not".
    property bool detailed: false
    property string planText: ""

    // The Gantt hover readout: what was computed at the row and instant under the pointer.
    // Magnitudes are added here because they belong to the target list, which QML owns.
    property string hoverText: ""

    // Stepping the date recomputes if a plan is already showing: the chart is what the date
    // control is FOR, and leaving it stale after a step is how a user reads last week's night.
    function shiftDate(days, months) {
        root.dateISO = Julia.shell_shift_date(root.dateISO, days, months)
        if (root.hasPlan) root.computePlan()
    }

    function computePlan() {
        var t = targetModel.count > 0 ? targetModel.get(Math.max(0, currentTargetIndex)) : null
        if (t === null) { planText = "no target"; return }
        var tels = root.selectedTelescopes.join(" ")

        // Search BEFORE planning, so the night is computed with the POPs it reports rather
        // than with the previous ones. Assigning popString is enough for the dropdowns and the
        // text field to follow: they read `popAt`, which reads popString.
        if (root.autoPops) root.searchAndAdoptPops(t, tels)

        var reply = Julia.shell_gantt(root.facility, t.name, t.ra, t.dec, root.dateISO,
                                      tels, root.popString, root.useDelay, root.detailed,
                                      root.altLimit, root.altMax)
        // `summary \t dark window`. The summary is not shown here: it is long, it overflowed
        // the strip it sat in, and the console already keeps it. The chart itself is the
        // answer, and hovering a bar gives the numbers.
        var f = reply.split("\t")
        planText = f[0]
        if (f.length > 1 && f[1].length > 0) root.darkWindowText = f[1]
        hasPlan = planText.indexOf("!") !== 0
        ganttArea.update()
        delayArea.update()
        root.consoleChanged()
    }

    // ── POP helpers ──────────────────────────────────────────────────────────
    //
    // `popString` is the single source of truth — the dropdowns, the text field and the search
    // results all read and write it. Two representations of one setting that can disagree is
    // how a user ends up computing a night with POPs they can see are not the ones shown.
    // Runs the POP search and takes the best row. Returns quietly on failure -- an empty
    // result means no configuration reaches the target, which the plan itself then reports as
    // zero observable hours rather than as an error here.
    function searchAndAdoptPops(t, tels) {
        var rows = Julia.shell_best_pops(root.facility, t.dec, root.dateISO, t.ra, tels, 3)
        if (!rows || rows.length === 0) return false
        var best = rows.split("\n")[0].split("\t")[0]
        if (!best || best.length === 0) return false
        root.popString = best
        return true
    }

    // The magnitudes for whichever target the hover names. They live in `targetModel`, so
    // Julia has never seen them and cannot put them in the readout itself.
    function magnitudeLine(hover) {
        var name = hover.split("\n")[0]
        for (var i = 0; i < targetModel.count; ++i) {
            var t = targetModel.get(i)
            if (t.name !== name) continue
            var parts = []
            var bands = [["V", t.magV], ["J", t.magJ], ["H", t.magH], ["K", t.magK]]
            for (var b = 0; b < bands.length; ++b)
                if (!isNaN(bands[b][1])) parts.push(bands[b][0] + " " + bands[b][1].toFixed(2))
            return parts.length > 0 ? "\n" + parts.join("   ") : ""
        }
        return ""
    }

    function popAt(i) {
        var f = popString.trim().split(/\s+/)
        if (i < 0 || i >= f.length) return 1
        var v = parseInt(f[i])
        return (isNaN(v) || v < 1 || v > 5) ? 1 : v
    }

    function setPop(i, v) {
        var n = telescopeNames.length
        var f = popString.trim().length > 0 ? popString.trim().split(/\s+/) : []
        while (f.length < n) f.push("1")
        f[i] = String(v)
        popString = f.slice(0, n).join(" ")
        if (hasPlan && useDelay) computePlan()
    }

    // Accepts whatever a paste looks like — spaces, commas, tabs — because the numbers get
    // copied out of tables and messages, not typed in a fixed format.
    function setPopString(text) {
        var f = text.trim().split(/[\s,]+/).filter(function (t) { return t.length > 0 })
        var out = []
        for (var i = 0; i < f.length && i < telescopeNames.length; ++i) {
            var v = parseInt(f[i])
            out.push(String(isNaN(v) || v < 1 || v > 5 ? 1 : v))
        }
        popString = out.join(" ")
        if (hasPlan && useDelay) computePlan()
    }

    // Taking a row from the search table by hand. AutoPOPs has to stop, or the next plan runs
    // the search again and replaces what was just chosen -- and the left-panel dropdowns, which
    // are disabled while AutoPOPs is on, would show a setting the user cannot touch.
    function adoptPops(text) {
        root.autoPops = false
        setPopString(text)
    }

    function findPops() {
        var t = targetModel.count > 0 ? targetModel.get(Math.max(0, currentTargetIndex)) : null
        if (t === null) return
        popModel.clear()
        var rows = Julia.shell_best_pops(root.facility, t.dec, root.dateISO, t.ra,
                                         root.selectedTelescopes.join(" "), nBestBox.value)
        root.consoleChanged()
        if (rows.length === 0) return
        var lines = rows.split("\n")
        for (var i = 0; i < lines.length; ++i) {
            var f = lines[i].split("\t")
            // The delegate declares these three as `required`, and a required role the model
            // does not supply means the delegate is never created — no rows, no error.
            if (f.length === 3)
                popModel.append({ pops: f[0], score: f[1], haRange: f[2] })
        }
        // The best one is the one worth using, so take it rather than making the user copy it.
        if (popModel.count > 0) root.popString = popModel.get(0).pops
    }

    property string popString: ""

    // Raised when Julia has written to the shared console and the window should re-read it.
    signal consoleChanged()

}
