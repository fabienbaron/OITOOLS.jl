#!/usr/bin/env bash
# Which file-dialog backend will Qt use, and does its window outlive the selection?
#
#     test/gui/dialog_backend.sh
#
# The file picker is not one thing. Qt6 picks one of three implementations depending on the
# session, they look different, they are different processes, and they close by different
# routes:
#
#   Qt Quick     drawn by Qt itself, Qt-styled, its own top-level window. Used when nothing
#                better is available -- no D-Bus session, no portal.
#   XDG portal   drawn by xdg-desktop-portal-gtk/-kde in ANOTHER PROCESS, GTK- or KDE-styled.
#                Used on a normal desktop session. Qt asks for it over D-Bus and hands over a
#                parent-window handle; if that handle does not survive, the window manager has
#                no reason to keep the dialog above the main window.
#   Wayland      either of the above, but as an xdg_toplevel whose create/destroy sequence is
#                not the X11 one.
#
# A bug report of "the picker will not close" means something different in each case, so this
# prints the facts that decide it rather than guessing.
set -uo pipefail

say() { printf '  %-26s %s\n' "$1" "$2"; }

echo "── session ─────────────────────────────────────────────────────────────"
say "DISPLAY"             "${DISPLAY:-<unset>}"
say "WAYLAND_DISPLAY"     "${WAYLAND_DISPLAY:-<unset>}"
say "XDG_SESSION_TYPE"    "${XDG_SESSION_TYPE:-<unset>}"
say "XDG_CURRENT_DESKTOP" "${XDG_CURRENT_DESKTOP:-<unset>}"
say "QT_QPA_PLATFORM"     "${QT_QPA_PLATFORM:-<unset, Qt decides>}"
say "QT_QPA_PLATFORMTHEME" "${QT_QPA_PLATFORMTHEME:-<unset, Qt decides>}"

echo
echo "── portal ──────────────────────────────────────────────────────────────"
if [ -n "${DBUS_SESSION_BUS_ADDRESS:-}" ]; then
    say "session bus" "present"
    if command -v busctl >/dev/null 2>&1 &&
       busctl --user status org.freedesktop.portal.Desktop >/dev/null 2>&1; then
        say "portal service" "RUNNING — Qt will very likely use the portal dialog"
    else
        say "portal service" "not running (it can still start on demand)"
    fi
else
    say "session bus" "ABSENT — Qt falls back to its own Qt Quick dialog"
fi
for impl in gtk kde gnome; do
    [ -e "/usr/share/xdg-desktop-portal/portals/$impl.portal" ] && say "portal backend" "$impl installed"
done

echo
echo "── what to look at ─────────────────────────────────────────────────────"
cat <<'TXT'
  Open the GUI, press "Open OIFITS…", and answer two questions.

  1. Does the picker look like Qt (a plain list, an "Open" button) or like your desktop's
     own file chooser (GTK/KDE styling, a "Select" or "Open" button, a Recent/Home sidebar)?
     That names the backend.

  2. With the picker on screen, run this in another terminal and keep the output:

         xdotool search --onlyvisible --name . | while read w; do
             printf '%s  %s\n' "$w" "$(xdotool getwindowname $w)"
         done

     Then pick a file, press Open, wait for the plot, and run it AGAIN.
     If a window named after the dialog is still listed the second time, the dialog outlived
     the selection; if it is gone, what remains on screen is a SECOND dialog opened earlier.

  On Wayland xdotool sees nothing, which is itself the answer: the dialog is a Wayland
  surface, and this is the case the X11 reasoning above does not cover.
TXT
