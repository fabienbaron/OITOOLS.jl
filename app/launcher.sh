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

exec "$here/bin/OITOOLSApp" "$@"
