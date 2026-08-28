#!/usr/bin/env bash
# Capture the four perspectives for README.md and the docs.
#
#     test/gui/screenshots.sh                 # headless, via Xvfb
#     DATA=demos/data/v1295Aql.oifits test/gui/screenshots.sh
#
# Writes docs/src/assets/gui/{exploring,observing,modeling,imaging}.png.
#
# Xvfb rather than a real screen, so the images are the same size and free of whatever
# window decorations, themes and notification bubbles the developer's desktop happens to have.
# The window is 83% x 95% of the screen (Main.qml:112), so the screen size below sets the
# shot size.
#
# xdotool drives X11 only. Under Xvfb there is no WAYLAND_DISPLAY, so Qt takes its X11 path and
# the window is a real X window -- which is exactly why this cannot be run against a Wayland
# session directly.
set -uo pipefail
cd "$(dirname "$0")/../.."

DISP=${SHOT_DISPLAY:-:97}
GEOM=${GEOM:-1920x1200x24}
DATA=${DATA:-demos/data/2004-data1.oifits}
OUT=${OUT:-docs/src/assets/gui}
LOG=$(mktemp /tmp/shots-XXXX.log)
say() { printf '  %s\n' "$*"; }

mkdir -p "$OUT"
[ -f "$DATA" ] || { say "no such data file: $DATA"; exit 1; }

Xvfb "$DISP" -screen 0 "$GEOM" >/dev/null 2>&1 &
XPID=$!
trap 'kill $XPID 2>/dev/null' EXIT
sleep 2

say "launching the GUI on $DISP ($GEOM) with $DATA"
DISPLAY="$DISP" julia --project=bin bin/oitoolsgui.jl "$DATA" > "$LOG" 2>&1 &
GPID=$!

WID=""
for i in $(seq 1 300); do
    kill -0 "$GPID" 2>/dev/null || { say "the GUI exited before showing a window"; tail -15 "$LOG"; exit 1; }
    WID=$(DISPLAY="$DISP" xdotool search --name "^OITOOLS$" 2>/dev/null | head -1)
    [ -n "$WID" ] && break
    sleep 1
done
[ -n "$WID" ] || { say "no window after 300 s"; tail -15 "$LOG"; exit 1; }
say "window $WID up; letting it settle"
sleep 6

CL()   { DISPLAY="$DISP" xdotool mousemove --window "$WID" "$1" "$2" click 1; sleep "${3:-3}"; }
# Escape before every tab switch. A modal dialog swallows the click otherwise, and the capture
# then shows the previous tab under a new name -- which is exactly what happened on the run
# that produced an "imaging.png" of the Modeling tab.
ESC()  { DISPLAY="$DISP" xdotool key --window "$WID" Escape; sleep 1; }
SHOT() { DISPLAY="$DISP" import -window "$WID" "$OUT/$1.png" 2>/dev/null && \
         say "$(printf '%-10s %s' "$1" "$(identify -format '%wx%h  %b' "$OUT/$1.png")")"; }

# One click per tab, left to right. The TabBar fills the window width and the four buttons
# share it equally, so the centres are at 1/8, 3/8, 5/8 and 7/8 of the width -- computed from
# the window rather than hardcoded, so changing GEOM does not silently capture the same tab
# four times. (It did exactly that on the first run: the pixel guesses all landed on Exploring
# and three of the four files came out byte-identical.)
WW=$(DISPLAY="$DISP" xdotool getwindowgeometry --shell "$WID" | sed -n 's/^WIDTH=//p')
say "window width $WW; tab centres at 1/8, 3/8, 5/8, 7/8"
TAB() { ESC; CL $(( WW * (2 * $1 + 1) / 8 )) "${TABY:-46}"; }

# Each tab is driven to something worth looking at first. An idle tab photographs as a page of
# placeholders -- "no reconstruction yet", "choose a target and a date" -- which documents the
# layout and nothing else.
#
# The coordinates below are window-relative and were read off a capture at GEOM=1920x1200,
# which puts the window at 1594x1140. Change GEOM and they must be re-read; the tab centres
# above are computed, these are not.

TAB 0                                  # Exploring: the uv plot is drawn on load
SHOT exploring

TAB 1                                  # Observing: pick the first target, then Compute
CL 18 116 1
CL 45 963 6
SHOT observing

TAB 2                                  # Modeling: a model with something in it
# "+ component" opens a modal dialog (name / kind / OK); OK is at 807,643. Twice, so the
# parameter table has more than one component in it.
# The button MOVES: each component adds a chip to its left, so the second add is not at the
# same x as the first. Clicking 51 twice adds one component and then selects its chip.
for x in 51 142; do
    CL "$x" 110 2                      # + component, at its shifted position
    CL 807 643 2                       # OK
done
SHOT modeling

TAB 3                                  # Imaging: an actual reconstruction
CL 45 1068 2                           # Run
say "waiting for the reconstruction"
sleep "${RUNWAIT:-75}"
SHOT imaging

kill "$GPID" 2>/dev/null; wait "$GPID" 2>/dev/null
say "log kept at $LOG"
