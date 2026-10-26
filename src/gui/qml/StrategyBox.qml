// The POP search strategy, in two sections: single target and multiple targets. A section's
// title is a row that cannot be chosen; an entry with `off` is listed but cannot be chosen.

import QtQuick
import QtQuick.Controls

ComboBox {
    id: box

    // [{ text, strategy }]: `strategy` < 0 for a section title.
    property var rows: []
    property int strategy: 0
    property real uiScale: 1.0
    signal chosen(int strategy)

    function rowOf(st) {
        for (var i = 0; i < rows.length; ++i) if (rows[i].strategy === st) return i
        return -1
    }
    function sync() { currentIndex = rowOf(strategy) }

    model: rows
    textRole: "text"
    implicitHeight: Math.max(Math.round(26 * uiScale), implicitContentHeight + topPadding + bottomPadding)
    // "&&" is a literal "&" in the items, which read "&" as a mnemonic.
    displayText: currentText.replace("&&", "&")
    Component.onCompleted: sync()
    onStrategyChanged: sync()
    onRowsChanged: sync()

    delegate: ItemDelegate {
        required property var modelData
        required property int index
        width: box.width
        text: modelData.text
        enabled: modelData.strategy >= 0 && !modelData.off
        font.bold: modelData.strategy < 0
        leftPadding: Math.round((modelData.strategy >= 0 ? 22 : 6) * box.uiScale)
        highlighted: box.highlightedIndex === index
    }

    onActivated: function (index) {
        var r = rows[index]
        if (r.strategy < 0 || r.off) { sync(); return }
        chosen(r.strategy)
    }
}
