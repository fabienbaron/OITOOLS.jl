#!/usr/bin/env bash
# Click-driven GUI test: drives the real widgets with xdotool on a virtual screen.
#
#     test/gui_click.sh              # headless, via Xvfb
#     test/gui_click.sh --display    # against $DISPLAY, e.g. your own machine
#
# WHY THIS EXISTS SEPARATELY FROM gui_xvfb.jl.
#
# gui_xvfb.jl calls the shell callbacks directly. That covers everything reachable from Julia,
# and it PASSES even on a machine where the GUI is visibly broken -- because it never opens a
# POPUP or a DIALOG. This script does, and that turns out to be the whole bug:
#
#     ! could not load v1295Aql.oifits:
#           gl_renderobject: glGenBuffers returned invalid id. OpenGL Context active?
#
# Opening the file dialog leaves a different GL context current. `draw!` then calls
# `empty!(ax)` and builds fresh plots, and building a plot allocates GL buffers -- which fails
# with no context bound. It reproduces here on llvmpipe, so the fix is verifiable without a
# WSL machine.
#
# Coordinates below were read off screenshots (import -window), not guessed. If the layout
# changes they must be re-read; a missed click shows up as a missing console line, not as a
# silent pass.
set -uo pipefail
# Repo root: this script lives in test/gui/, and --project=bin needs the root.
cd "$(dirname "$0")/../.."

DUMP="$(mktemp)"
JLOG="$(mktemp)"
SHOTDIR="${OITOOLSGUI_SHOTS:-$(mktemp -d)}"
FAIL=0
GUI_PID=""
XVFB_PID=""

if [ "${1:-}" = "--display" ]; then
    DISP="${DISPLAY:?--display needs DISPLAY set}"
else
    DISP=":99"
fi

cleanup() {
    [ -n "$GUI_PID" ] && kill "$GUI_PID" 2>/dev/null
    [ -n "$XVFB_PID" ] && kill "$XVFB_PID" 2>/dev/null
    return 0
}
trap cleanup EXIT

say() { printf '  [%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

say "display      $DISP"
say "julia log    $JLOG"
say "console dump $DUMP"
say "screenshots  $SHOTDIR"
if [ "$DISP" != ":99" ]; then
    say "NOTE: on a real display this moves your mouse pointer and clicks. Keep hands off"
    say "      for ~90 s, and do not let another window take focus."
fi

if [ "$DISP" = ":99" ]; then
    Xvfb "$DISP" -screen 0 1400x1000x24 >/dev/null 2>&1 &
    XVFB_PID=$!
    sleep 2
fi

# Pin the scale: the coordinates below were read off screenshots at 1.0, and the default is
# deliberately 1.25. A cosmetic default must not decide whether the test can find a button.
say "precompiling (silent first run otherwise: this is the slow part) ..."
if ! DISPLAY="$DISP" julia --project=bin -e 'using Pkg; Pkg.precompile()' > "$JLOG" 2>&1; then
    say "FAIL  precompilation failed:"; tail -25 "$JLOG"; exit 1
fi
say "precompiled; starting the GUI"

# Pin the picker to this directory: the GUI now opens in demos/data, and this script
# selects its file by row.
OITOOLSGUI_DATA_DIR="$PWD/test/gui/data" OITOOLSGUI_DEBUG_PICK="${OITOOLSGUI_DEBUG_PICK:-}" OITOOLSGUI_SCALE=1.0 OITOOLSGUI_CONSOLE_DUMP="$DUMP" DISPLAY="$DISP" \
${QT_QPA_PLATFORM:+QT_QPA_PLATFORM="$QT_QPA_PLATFORM"} \
  julia --project=bin -e '
      # The order is forced: Mesa reads its driver variables and GLFW its platform hint when
      # the first OpenGL context is created, so both calls must happen ABOVE `using GLMakie`.
      # Loading GLMakie first is why GLFW then dies with "X11: The DISPLAY environment
      # variable is missing" in a Wayland session -- GLFW.jl hard-codes X11.
      using OITOOLS
      configure_graphics!()
      using GLFW_jll
      prefer_native_wayland!()
      using GLMakie, QMLMakie, QML
      GUI = Base.get_extension(OITOOLS, :OITOOLSGUIExt)
      @info "packages loaded; opening the window"
      OITOOLS.gui(; autoquit_ms = 70000)' > "$JLOG" 2>&1 &
GUI_PID=$!

WID=""
for i in $(seq 1 240); do
    WID=$(DISPLAY="$DISP" xdotool search --name "^OITOOLS$" 2>/dev/null | head -1)
    [ -n "$WID" ] && break
    if ! kill -0 "$GUI_PID" 2>/dev/null; then
        say "FAIL  the GUI process exited before a window appeared:"; tail -25 "$JLOG"; exit 1
    fi
    if [ $((i % 15)) -eq 0 ]; then
        say "still waiting for the window (${i}s) ...  last log line: $(tail -1 "$JLOG" 2>/dev/null | cut -c1-90)"
    fi
    sleep 1
done
if [ -z "$WID" ]; then
    say "FAIL  no window titled OITOOLS after 240s."
    # A window that is on screen but invisible to xdotool means the client is not an X11
    # client at all. Under WSLg (and any Wayland session) Qt prefers the wayland plugin
    # whenever WAYLAND_DISPLAY is set, and a Wayland surface has no X window for XTEST to
    # find. That is not a hang: the GUI is running, just not reachable this way.
    if [ -n "${WAYLAND_DISPLAY:-}" ]; then
        say ""
        say "      WAYLAND_DISPLAY=$WAYLAND_DISPLAY is set, so Qt very likely opened a"
        say "      WAYLAND window. xdotool drives X11 only and cannot see or click it."
        say "      Re-run forcing Qt onto X11:"
        say ""
        say "        QT_QPA_PLATFORM=xcb ./test/gui_click.sh --display"
        say ""
        say "      If that works, the GUI itself is fine and only the automation was blind."
    fi
    say "      xdotool sees $(DISPLAY="$DISP" xdotool search --name . 2>/dev/null | wc -l) X windows in total"
    say "--- julia log ---"
    if [ -s "$JLOG" ]; then tail -25 "$JLOG"; else say "(the julia log is empty)"; fi
    exit 1
fi
say "window $WID found; driving it now"
sleep 3

CL()   { DISPLAY="$DISP" xdotool mousemove --window "$WID" "$1" "$2" click 1; sleep "${3:-1.5}"; }
SHOT() { DISPLAY="$DISP" import -window "$WID" "$SHOTDIR/$1.png" 2>/dev/null; }

# Explore is the first tab and selected on start, so there is nothing to click. Click it
# anyway: if the tab order is ever changed again, this is what notices.
CL 160 42 2          # Explore tab (first of four)
SHOT 01_explore

# -- the file picker: opening it is what breaks the GL context ----------------
#
# Driven by keyboard, not by coordinates. The picker gives its list focus on open, the entries
# are sorted with directories first, and Return opens the current row -- so row 0 is
# 2004-data1.oifits and the plot-pick coordinates further down stay valid. Picking by pixel
# silently changed WHICH file was loaded when the picker changed, and every later coordinate
# went stale with it.
CL 470 14 3          # "Open OIFITS..."
SHOT 02_dialog
DISPLAY="$DISP" xdotool key Down; sleep 1     # row 1: v1295Aql.oifits (polychromatic)
DISPLAY="$DISP" xdotool key Return; sleep 8
# The polychromatic file on purpose: it carries 9864 uv points against the mono file's 494, and
# the point-picking check further down clicks a fixed spot. On a sparse plot that spot lands on
# empty space, and a miss is indistinguishable from a broken pick.
SHOT 03_after_open

# -- redraws after the dialog has been open ----------------------------------
CL 105 73 2 ; DISPLAY="$DISP" xdotool key Down; sleep 1; DISPLAY="$DISP" xdotool key Return; sleep 3
SHOT 04_plot_kind
CL 313 73 2 ; DISPLAY="$DISP" xdotool key Down; sleep 1; DISPLAY="$DISP" xdotool key Return; sleep 3
SHOT 05_colour

# -- click a data point, and reset the view ----------------------------------
# Picking is read from Makie's event stream, not from a QML handler, so the only way to know
# it works is to click the plot and look at the label underneath. Two clicks: the first gives
# the plot keyboard focus (MakieArea gates its input on that), the second actually picks.
# Click where the points actually are, not an arbitrary spot: the search returns nothing over
# empty space, so a miss is indistinguishable from a broken pick. By this point the view is V²
# against baseline, and the dense band sits a little below the middle of the axis.
#
# As a FRACTION of the window, not a fixed pixel. A fixed pixel silently stops landing on data
# the moment the window's default size changes -- which is what happened when the default width
# went from 1280 to 1600: 700,440 had been on the band and became empty sky above it, and the
# run then reported a broken pick rather than a stale coordinate. Fractions ride the resize,
# because the axis fills its share of the window either way.
eval "$(DISPLAY="$DISP" xdotool getwindowgeometry --shell "$WID")"
PICK_X=$((WIDTH * 40 / 100)); PICK_Y=$((HEIGHT * 49 / 100))
say "picking at ${PICK_X},${PICK_Y} of ${WIDTH}x${HEIGHT}"
CL "$PICK_X" "$PICK_Y" 2
CL "$PICK_X" "$PICK_Y" 3
SHOT 06_picked
# right click = reset view; anywhere inside the axis will do
DISPLAY="$DISP" xdotool mousemove --window "$WID" $((WIDTH * 43 / 100)) $((HEIGHT * 28 / 100)) click 3
sleep 3
SHOT 07_reset

wait "$GUI_PID" 2>/dev/null

echo "-- console transcript --"
sed 's/^/  /' "$DUMP"

echo "-- checks --"
if grep -qF "OITOOLSGUI ready" "$DUMP"; then echo "  PASS  GUI started"
else echo "  FAIL  GUI never started"; FAIL=1; fi

if grep -qF "load_dataset!" "$DUMP"; then echo "  PASS  the file dialog picked a file"
else echo "  FAIL  the file dialog picked nothing (clicks missed?)"; FAIL=1; fi

# THE regression check. Anything drawn after a dialog or popup must still draw.
if grep -qF "glGenBuffers" "$DUMP" || grep -qF "gl_renderobject" "$DUMP"; then
    echo "  FAIL  GL context lost after the dialog (this is the known bug)"; FAIL=1
else
    echo "  PASS  no GL context loss after the dialog"
fi

if grep -qE "^[0-9:]+ +!" "$DUMP"; then
    echo "  FAIL  error lines in the transcript"; FAIL=1
else
    echo "  PASS  no error lines"
fi

# Picking reads Makie's event stream and searches the drawn points, because Makie.pick needs
# a buffer QMLMakie never renders. Clicking the plot is the only way to check it.
PICK=$(grep "^PICKTEXT: " "$DUMP" | sed 's/^PICKTEXT: //')
if [ -n "$PICK" ]; then
    echo "  PASS  clicking identified a point: $PICK"
else
    echo "  FAIL  clicking a point identified nothing"
    echo "        Look at 06_picked.png first: if the click landed on empty axis rather than"
    echo "        on the data, the fractions above need re-reading, not the picking code."
    FAIL=1
fi

echo "  screenshots: $SHOTDIR"
exit $FAIL
