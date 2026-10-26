#!/bin/sh
#
# Start the application, having settled the windowing system first.
#
# THIS IS NOT OPTIONAL ON LINUX. A bundle is a sysimage, and every package `__init__` in a
# sysimage runs before any user code: GLFW has called `glfwInit` and QML has constructed its
# `QGuiApplication` before `julia_main` gets control. Neither choice can be revised from
# inside the process, so the only moment left is here, before it starts.
#
# Left alone on a Wayland session, GLFW takes its hardcoded X11 and Qt follows the session to
# Wayland. That is two windowing systems and two EGL display connections in one process, and
# it is the configuration OITOOLS' `prefer_native_wayland!` exists to avoid.
#
# Anything set by hand is honoured -- this only fills in a blank, so
#
#     QT_QPA_PLATFORM=xcb ./OITOOLS          # both halves on X11 instead
#
# still does what it says.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)

if [ -n "${WAYLAND_DISPLAY:-}" ] && [ -z "${JULIA_GLFW_PLATFORM:-}" ] && [ -z "${QT_QPA_PLATFORM:-}" ]; then
    JULIA_GLFW_PLATFORM=wayland
    export JULIA_GLFW_PLATFORM
fi

# The bundled libxkbcommon has its BUILD sandbox path compiled in as the X11 locale directory
# (/workspace/destdir/share/X11/locale), so on a Wayland session it cannot find the compose
# tables and every launch reports "No Compose file for locale" while dead keys stop working.
# It does read XLOCALEDIR. Bundle copy first, host second, and never override a value the user
# set -- someone who has pointed this at their own tables meant it.
if [ -z "${XLOCALEDIR:-}" ]; then
    for d in "$here/share/X11/locale" /usr/share/X11/locale /usr/local/share/X11/locale; do
        if [ -f "$d/compose.dir" ]; then XLOCALEDIR=$d; export XLOCALEDIR; break; fi
    done
fi

# Qt's scale, as `configure_qt_scale!` (src/graphics.jl) sets it in a checkout: the saved
# `window_scale` when positive, else the panel's density / 150 to the nearest 0.25 (at least 1):
# on X11 the EDID density, divided by what Xft.dpi already applies; on WSL, whose WSLg gives Qt
# a scale of 1, 120 dpi per 100% of Windows' display scale (AppliedDPI in the registry). Left alone when any Qt scale
# variable is set, on a Wayland session other than WSLg (its compositor scales Qt), and when the
# connected displays differ in density.
qt_scale() {
    for v in QT_SCALE_FACTOR QT_SCREEN_SCALE_FACTORS QT_ENABLE_HIGHDPI_SCALING QT_AUTO_SCREEN_SCALE_FACTOR; do
        eval "[ -n \"\${$v+x}\" ]" && return 0
    done
    cfg=${XDG_CONFIG_HOME:-$HOME/.config}/oitools/gui.toml
    saved=$(sed -n 's/^window_scale *= *\([0-9.]*\).*/\1/p' "$cfg" 2>/dev/null | head -n 1)
    if [ -n "$saved" ] && awk "BEGIN { exit !($saved > 0) }"; then
        f=$saved
    elif [ -n "${WSL_DISTRO_NAME:-}" ]; then
        reg=$(command -v reg.exe || echo /mnt/c/Windows/System32/reg.exe)
        [ -x "$reg" ] || return 0
        hex=$("$reg" query 'HKCU\Control Panel\Desktop\WindowMetrics' /v AppliedDPI 2>/dev/null |
              sed -n 's/.*REG_DWORD[[:space:]]*0x\([0-9a-fA-F]*\).*/\1/p' | head -n 1)
        [ -n "$hex" ] || return 0
        f=$(awk -v d="$((0x$hex))" 'BEGIN { if (d > 0) { s = int(4 * (120 * d / 96) / 150 + 0.5) / 4
                                                           if (s < 1) s = 1; printf "%.3f", s } }')
    else
        if [ -n "${WAYLAND_DISPLAY:-}" ] && [ "${QT_QPA_PLATFORM:-}" != xcb ]; then return 0; fi
        dpis=""
        for d in /sys/class/drm/card*-*; do
            [ "$(cat "$d/status" 2>/dev/null)" = connected ] || continue
            dpis="$dpis $(od -An -tu1 -j54 -N18 "$d/edid" 2>/dev/null | tr '\n' ' ' | awk '
                NF >= 15 && !($1 == 0 && $2 == 0) {
                    w = $3 + int($5 / 16) * 256; wm = $13 + int($15 / 16) * 256
                    if (w > 0 && wm > 0) printf "%.2f", w / (wm / 25.4) }')"
        done
        xft=$(xrdb -query 2>/dev/null | sed -n 's/^Xft\.dpi:[[:space:]]*\([0-9.]*\).*/\1/p' | head -n 1)
        f=$(echo $dpis | awk -v xft="${xft:-96}" '
            { for (i = 1; i <= NF; i++) { if (i == 1 || $i < lo) lo = $i; if (i == 1 || $i > hi) hi = $i } n = NF }
            END { if (!n || hi - lo >= 1) exit
                  s = int(4 * lo / 150 + 0.5) / 4; if (s < 1) s = 1; if (xft <= 0) xft = 96
                  printf "%.3f", s / (xft / 96) }')
    fi
    [ -n "$f" ] || return 0
    awk "BEGIN { exit !($f > 1.01 || $f < 0.99) }" || return 0
    QT_SCALE_FACTOR=$f
    export QT_SCALE_FACTOR
}
qt_scale || true

exec "$here/bin/OITOOLSApp" "$@"
