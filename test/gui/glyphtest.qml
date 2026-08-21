// Window for test/gui/glyphtest.jl. Deliberately almost empty: the point is to put a Makie
// figure inside a Qt GL context with as little else as possible, so anything seen here cannot
// be blamed on the rest of the GUI.
//
// The decoy figures matter. The application builds six Makie figures into one context, and
// this instantiates the first `nfigs` of them as MakieAreas — only the first is visible, but
// the others exist and have screens, exactly as the perspectives do.

import QtQuick
import QtQuick.Window
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import Makie

ApplicationWindow {
    id: win
    visible: true
    width: 1150
    height: 900
    title: "glyph test"

    // Sets the labels once the window is up, which is when the GUI sets its own — and when
    // the GL context already belongs to Qt's render thread.
    Timer {
        interval: 2000
        running: applyAfter
        repeat: false
        onTriggered: { Julia.apply_labels(); area1.update() }
    }

    // Sets the labels once the window is up, which is when the GUI sets its own — and when the
    // GL context already belongs to Qt's render thread.
    Timer {
        interval: 2000
        running: applyAfter
        repeat: false
        onTriggered: { Julia.apply_labels(); area1.update() }
    }

    // Closes itself so a harness can screenshot and move on. 0 means stay open.
    Timer {
        interval: autoQuitMs > 0 ? autoQuitMs : 1
        running: autoQuitMs > 0
        repeat: false
        onTriggered: Qt.quit()
    }

    ColumnLayout {
        anchors.fill: parent
        anchors.margins: 6
        spacing: 6

        // Qt's own rendering of the same characters, as a side-by-side control: if these are
        // clean while the Makie ones are not, Qt's text stack is fine and Makie's atlas is not.
        Label {
            Layout.fillWidth: true
            text: caption
            font.family: "monospace"
            color: "#444"
        }
        Label {
            Layout.fillWidth: true
            text: "Qt draws: v2 — AZCYG2011_FINAL2018.oifits   Baseline (Mλ)   V²   E1-E2"
            font.pointSize: 11
        }
        Rectangle { Layout.fillWidth: true; Layout.preferredHeight: 1; color: "#ccc" }

        Rectangle {
            Layout.fillWidth: true
            Layout.fillHeight: true
            color: "white"
            border.color: "#ddd"

            MakieArea {
                id: area1
                anchors.fill: parent
                anchors.margins: 1
                scene: plot1
            }
        }

        // The decoys. Zero height: they are here to be CREATED, not to be looked at.
        Item {
            Layout.fillWidth: true
            Layout.preferredHeight: 0
            clip: true
            MakieArea { width: 2; height: 2; scene: plot2; visible: nfigs >= 2 }
            MakieArea { width: 2; height: 2; scene: plot3; visible: nfigs >= 3 }
            MakieArea { width: 2; height: 2; scene: plot4; visible: nfigs >= 4 }
            MakieArea { width: 2; height: 2; scene: plot5; visible: nfigs >= 5 }
            MakieArea { width: 2; height: 2; scene: plot6; visible: nfigs >= 6 }
        }
    }
}
