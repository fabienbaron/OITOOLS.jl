// A scrolling console for whatever an engine or optimiser printed.
//
// One component rather than one per tab: Imaging and Modeling both show the raw output of a
// long-running Julia call, and the two had begun to diverge in scroll behaviour and font.
//
// The text arrives as plain text and is coloured HERE, not in Julia. What comes out of an
// optimiser is its own business; which parts of it a reader needs to find quickly is the
// view's. Keeping the split that way also means the shell can go on returning something that
// can be diffed and tested as text.

import QtQuick
import QtQuick.Controls

Item {
    id: root

    // What the engine printed, verbatim.
    property string text: ""
    property real   fontPointSize: 9
    // Follow the tail as new output arrives. Off while the user is reading back through it.
    property bool   follow: true

    // A long run can print megabytes, and colouring all of it on every update is wasted work
    // nobody sees. Only the tail is rendered; the rest is still in the model.
    property int    maxLines: 2000

    implicitHeight: 160

    // Qt renders HTML in a RichText TextArea, so the payload has to be escaped before any
    // markup is added -- an optimiser that prints `<` would otherwise eat the rest of the line.
    function _escape(s) {
        return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;")
    }

    // OITOOLS' own colours, not ours. `chi2_flat` prints V2 red, T3amp blue and T3phi green;
    // `display_model` warns in yellow and passes in green. Those choices arrive as ANSI escape
    // sequences — the shell captures with `IOContext(:color => true)` so `printstyled` emits
    // them — and this turns them into the spans Qt understands. Guessing at meaning with
    // regexes was the alternative and it would have disagreed with the library.
    //
    // Tuned for a light background: ANSI "yellow" and "white" are illegible on white, so those
    // two are darkened rather than reproduced literally.
    readonly property var _ansi: ({
        "30": "#000000", "31": "#c62828", "32": "#1a7f37", "33": "#a06000",
        "34": "#1f4fd8", "35": "#a4288a", "36": "#0b7a8a", "37": "#666666",
        "90": "#888888", "91": "#e04b4b", "92": "#3aa35a", "93": "#b8860b",
        "94": "#4a72e8", "95": "#c05ac0", "96": "#2f9bab", "97": "#444444"
    })

    function _ansiToHtml(line) {
        var out = "", open = 0, i = 0
        var re = /\x1b\[([0-9;]*)m/g, m, last = 0
        while ((m = re.exec(line)) !== null) {
            out += _escape(line.substring(last, m.index))
            last = re.lastIndex
            var codes = m[1].length === 0 ? ["0"] : m[1].split(";")
            for (var c = 0; c < codes.length; ++c) {
                var code = codes[c]
                if (code === "0" || code === "39" || code === "22") {
                    while (open > 0) { out += "</span>"; open-- }
                } else if (code === "1") {
                    out += '<span style="font-weight:bold">'; open++
                } else if (_ansi[code] !== undefined) {
                    out += '<span style="color:' + _ansi[code] + '">'; open++
                }
            }
        }
        out += _escape(line.substring(last))
        while (open > 0) { out += "</span>"; open-- }
        return out
    }

    readonly property string _html: {
        if (text.length === 0) return ""
        var lines = text.split("\n")
        var from = Math.max(0, lines.length - maxLines)
        var out = []
        if (from > 0)
            out.push('<span style="color:#999">… ' + from + ' earlier lines not shown …</span>')
        for (var i = from; i < lines.length; ++i) out.push(_ansiToHtml(lines[i]))
        // <pre> keeps the column alignment the optimiser laid out. Without it Qt collapses the
        // runs of spaces that make an iteration table a table.
        return "<pre style=\"margin:0\">" + out.join("\n") + "</pre>"
    }

    ScrollView {
        id: scroller
        anchors.fill: parent
        clip: true
        ScrollBar.vertical.policy: ScrollBar.AsNeeded

        TextArea {
            id: area
            readOnly: true
            // Selectable: the useful thing to do with a trace that did not converge is paste it
            // somewhere else.
            selectByMouse: true
            textFormat: TextEdit.RichText
            wrapMode: TextArea.NoWrap
            font.family: "monospace"
            font.pointSize: root.fontPointSize
            text: root._html
            background: null
            // Jump to the newest output, which is the part worth seeing while a run is going.
            onTextChanged: if (root.follow) cursorPosition = length
        }
    }
}
