#!/usr/bin/env bash
# The click test, on real GPU hardware and native Wayland.
#
#     test/gui_click_wayland.sh
#
# gui_click.sh runs under Xvfb, which means llvmpipe and X11. This runs the same interactions
# under a nested GPU Weston (x11 backend, GL renderer -> EGL/GBM -> renderD128 -> radeonsi),
# with GLMakie forced onto native Wayland rather than XWayland. Two things differ from the
# Xvfb path, and both are exactly where the plot surface has already had one bug:
#
#   * Qt Quick popups are xdg_popup SURFACES here, not X11 child windows, so the
#     create/destroy sequence around a ComboBox or FileDialog is genuinely different.
#   * GLMakie is on Wayland too. Without OITOOLSGUI's prefer_native_wayland!, GLFW.jl's
#     hard-coded X11 would put Qt on Wayland and GLMakie on XWayland -- two windowing
#     systems, two EGL connections, in one process.
#
# Input goes through XTEST against the compositor's own X11 window, which feeds its wl_seat.
# Two requirements, both verified before starting:
#   * the physical X session must be UNLOCKED (kscreenlocker grabs input and swallows XTEST)
#   * plain `xdotool click`, never `xdotool key --window` (XSendEvent, which Weston ignores)
set -uo pipefail
# Repo root: this script lives in test/gui/, and --project=bin needs the root.
cd "$(dirname "$0")/../.."

WL_SOCK="wl-guitest"
RUNTIME_DIR="${XDG_RUNTIME_DIR:-/run/user/$(id -u)}"
LOGDIR="/tmp/wayland-gui-$(id -u)"; mkdir -p "$LOGDIR"
WLOG="$LOGDIR/weston.log"
DUMP="$(mktemp)"
READY="$(mktemp)"; : > "$READY"
JLOG="$LOGDIR/julia.log"
SHOTDIR="${OITOOLSGUI_SHOTS:-$LOGDIR/shots}"; mkdir -p "$SHOTDIR"
FAIL=0
WESTON_PID=""; GUI_PID=""

export DISPLAY="${HOST_DISPLAY:-:0}"
XA=$(tr '\0' '\n' < /proc/"$(pgrep -u "$(id -u)" -x plasmashell | head -1)"/environ 2>/dev/null \
     | sed -n 's/^XAUTHORITY=//p' | head -1)
export XAUTHORITY="${XA:-$HOME/.Xauthority}"

cleanup() {
    [ -n "$GUI_PID" ] && kill "$GUI_PID" 2>/dev/null
    [ -n "$WESTON_PID" ] && kill "$WESTON_PID" 2>/dev/null
    rm -f "$RUNTIME_DIR/$WL_SOCK" "$RUNTIME_DIR/$WL_SOCK.lock" 2>/dev/null
    return 0
}
trap cleanup EXIT

# -- preflight ---------------------------------------------------------------
if ! xdpyinfo >/dev/null 2>&1; then echo "  FAIL  cannot reach X display $DISPLAY"; exit 1; fi
LOCKED=$(loginctl show-session 2 -p LockedHint 2>/dev/null | cut -d= -f2)
if [ "$LOCKED" = "yes" ]; then
    echo "  FAIL  the X session is LOCKED; kscreenlocker swallows XTEST input."
    echo "        Unlock the physical session and re-run."
    exit 1
fi

# -- precompile FIRST, while an X display is still reachable -------------------
# Julia precompiles in a subprocess, and that subprocess runs neither prefer_native_wayland!
# nor with DISPLAY set (the client is launched with DISPLAY unset, deliberately). GLFW.jl's
# hard-coded X11 init then fails there and the extension never builds. So do it out here,
# where $DISPLAY still points at the real session.
echo "  precompiling with DISPLAY=$DISPLAY ..."
if ! julia --project=bin -e 'using Pkg; Pkg.precompile()' > "$JLOG" 2>&1; then
    echo "  FAIL  precompilation failed"; tail -20 "$JLOG"; exit 1
fi

# -- nested GPU compositor ---------------------------------------------------
rm -f "$RUNTIME_DIR/$WL_SOCK" "$RUNTIME_DIR/$WL_SOCK.lock" 2>/dev/null
XDG_RUNTIME_DIR="$RUNTIME_DIR" weston --backend=x11 --renderer=gl \
    --width=1400 --height=1000 --socket="$WL_SOCK" --no-config > "$WLOG" 2>&1 &
WESTON_PID=$!
for _ in $(seq 1 100); do [ -S "$RUNTIME_DIR/$WL_SOCK" ] && break; sleep 0.2; done
if [ ! -S "$RUNTIME_DIR/$WL_SOCK" ]; then echo "  FAIL  compositor did not start"; tail -20 "$WLOG"; exit 1; fi

if grep -qiE 'llvmpipe|softpipe|swrast' "$WLOG"; then
    echo "  FAIL  compositor fell back to software rendering"; grep -i "GL renderer" "$WLOG"; exit 1
fi
echo "  $(grep -m1 'GL renderer' "$WLOG" | sed 's/.*GL renderer/GL renderer/')"
WESTON_WID=$(grep -o 'window id [0-9]*' "$WLOG" | tail -1 | awk '{print $3}')
echo "  compositor window $WESTON_WID"

# -- the GUI, inside the compositor, on native Wayland ------------------------
XDG_RUNTIME_DIR="$RUNTIME_DIR" WAYLAND_DISPLAY="$WL_SOCK" XDG_SESSION_TYPE=wayland \
QT_QPA_PLATFORM=wayland QT_ASSUME_STDERR_HAS_CONSOLE=1 \
MESA_LOADER_DRIVER_OVERRIDE=radeonsi \
XCOMPOSEFILE=/usr/share/X11/locale/en_US.UTF-8/Compose \
OITOOLSGUI_SCALE=1.0 OITOOLSGUI_CONSOLE_DUMP="$DUMP" OITOOLSGUI_READY_FILE="$READY" \
env -u DISPLAY -u XAUTHORITY \
  julia --project=bin -e '
      # prefer_native_wayland! must run BEFORE `using GLMakie`: GLFW.jl hard-codes X11 at
      # module init, so with DISPLAY unset it dies with "X11: The DISPLAY environment variable
      # is missing" before the hint can be applied. GLFW_jll alone activates the extension
      # that defines it, which is the reason that extension is separate from the GUI one.
      using OITOOLS
      configure_graphics!()
      using GLFW_jll
      prefer_native_wayland!()
      using GLMakie, QMLMakie, QML
      GUI = Base.get_extension(OITOOLS, :OITOOLSGUIExt)
      # Fullscreen so screen coordinates map 1:1 onto client coordinates, and the click
      # arithmetic is just the compositor window origin.
      ENV["OITOOLSGUI_FULLSCREEN"] = "1"
      OITOOLS.gui(; autoquit_ms = 120000,
          on_ready = () -> write(ENV["OITOOLSGUI_READY_FILE"], "ready"))' > "$JLOG" 2>&1 &
GUI_PID=$!

# Readiness is the ready hook firing, NOT the console dump: the dump is written when the
# window CLOSES, so waiting on it would wait for the whole run to finish and then click into
# an empty compositor.
for _ in $(seq 1 180); do [ -s "$READY" ] && break; sleep 1; done
if [ ! -s "$READY" ]; then
    echo "  FAIL  the GUI never became ready inside the compositor"; tail -25 "$JLOG"; exit 1
fi
sleep 3

eval "$(xdotool getwindowgeometry --shell "$WESTON_WID")"     # X, Y, WIDTH, HEIGHT
echo "  compositor at ${X},${Y} ${WIDTH}x${HEIGHT}"
xdotool windowactivate --sync "$WESTON_WID" 2>/dev/null
sleep 1

# XTEST, in ROOT coordinates: compositor origin + position inside the compositor.
CL()   { xdotool mousemove $((X + $1)) $((Y + $2)) click 1; sleep "${3:-1.5}"; }
SHOT() { import -window "$WESTON_WID" "$SHOTDIR/$1.png" 2>/dev/null; }

SHOT 00_start
CL 480 42 2           # Explore tab
CL 470 14 3           # "Open OIFITS..."
SHOT 01_dialog
CL 700 342 1          # v1295Aql.oifits
CL 801 603 8          # the dialog's Open button
SHOT 02_after_open
CL 105 73 2 ; xdotool key Down; sleep 1; xdotool key Return; sleep 3
SHOT 03_kind
CL 313 73 2 ; xdotool key Down; sleep 1; xdotool key Return; sleep 3
SHOT 04_colour

wait "$GUI_PID" 2>/dev/null

echo "-- console transcript --"; sed 's/^/  /' "$DUMP"
echo "-- checks --"
grep -qF "OITOOLSGUI ready" "$DUMP" && echo "  PASS  GUI started on Wayland" || { echo "  FAIL  no start"; FAIL=1; }
grep -qF "load_dataset!"    "$DUMP" && echo "  PASS  file dialog picked a file" || { echo "  FAIL  dialog picked nothing"; FAIL=1; }
if grep -qE "glGenBuffers|gl_renderobject" "$DUMP"; then
    echo "  FAIL  GL context lost after the dialog"; FAIL=1
else
    echo "  PASS  no GL context loss after the dialog"
fi
grep -qE "^[0-9:]+ +!" "$DUMP" && { echo "  FAIL  error lines in the transcript"; FAIL=1; } || echo "  PASS  no error lines"
grep -qiE 'llvmpipe|softpipe|swrast' "$JLOG" && { echo "  FAIL  client fell back to software"; FAIL=1; } \
                                             || echo "  PASS  no software fallback"
echo "  screenshots: $SHOTDIR"
exit $FAIL
