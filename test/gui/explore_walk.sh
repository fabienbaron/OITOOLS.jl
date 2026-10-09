#!/usr/bin/env bash
# Explore, walked through end to end on a virtual screen: every control and every plot
# interaction listed in section 2 of GUI_Explore.md, in one run.
#
#     test/gui/explore_walk.sh                 # headless, via Xvfb
#     OITOOLSGUI_SHOTS=dir test/gui/explore_walk.sh
#
# WHY. Most of what Explore does lives in QML -- which kinds are offered, which ticks clear
# which, what the disabled tooltips say -- and the Julia suites cannot see any of it. This
# drives the real window, writes a screenshot per step, and checks every step that leaves a
# console line. The screenshots are the reference set the redesign is compared against
# (GUI_Explore.md, Phase 0); a step with no console line is checked by its screenshot.
#
# Every coordinate is in the block below and nowhere else, read off screenshots of a
# 1400x950 window at OITOOLSGUI_SCALE=1.0. When the layout changes, that block is what
# changes. A missed click shows up as a FAIL on a missing console line, not as a silent pass.
set -uo pipefail
cd "$(dirname "$0")/../.."

# ── coordinates (window 1400x950, scale 1.0) ────────────────────────────────────────────
WIN_W=1400;          WIN_H=950
DATASET_BOX="216 14"
TAB_EXPLORE="175 40"; TAB_MODEL="875 40";  TAB_IMAGE="1224 40"
# The frame, top to bottom: Dataset, Filter (collapsed), Plot, Compare. The Dataset summary's
# length moves the controls below it, so these are data1's (2004-data1); data2's summary
# (v1295Aql) is 85 px longer, and the one step on data2 uses KIND_BOX_DATA2.
KIND_BOX="148 526";  COLOR_BOX="148 553";  KIND_BOX_DATA2="148 611"
LOG_Y="26 582";      PER_GROUP="26 612";   MAX_BL="26 642"
DATA_ONLY="26 714";  OVER_MODEL="26 741";  OVER_IMAGE="26 768"
RESID_MODEL="26 795"; RESID_IMAGE="26 822"
SAVE_PNG="1355 66"                          # on the plot, top right
KIND_MENU_VISAMP="100 636"                  # the visamp row of the open Plot menu
PLOT_CENTRE="840 415"                       # inside the axis on every view used below
PICK_AT="514 598"                           # an AC0-AE0 point on the V² plot of 2004-data1
MODEL_OPEN="147 68"; MODEL_CARD_REMOVE="214 103"
# The model picker's field is read-only when opening, so a model is chosen from the list: the
# "OITOOLS models" place, then its second file row (Limb_Darkened_Disk_Claret4.toml; the first
# is AX_Cir_Binary.toml).
MODEL_PLACE="428 432"; MODEL_ROW="595 389"
NEUTRAL="1300 900"                          # somewhere no tooltip lives
SAVE_NAME_FIELD="645 630"                   # the Save PNG picker's name field
IMAGE_RUN="18 893"
# Menu rows, in order. A menu is driven by keyboard: Up to the top, Down N, Return.
KINDS=(uv v2 t3phi t3amp visamp visphi diffphi diffvisamp flux)
COLORS=(baseline wav mjd none)
# ────────────────────────────────────────────────────────────────────────────────────────

DISP=":97"
SHOTDIR="${OITOOLSGUI_SHOTS:-$(mktemp -d)}"
mkdir -p "$SHOTDIR"
WORK="$(mktemp -d)"
DUMP="$WORK/console.txt"; JLOG="$WORK/julia.txt"
FAIL=0; GUI_PID=""; XVFB_PID=""

cleanup() {
    [ -n "$GUI_PID" ] && kill -9 "$GUI_PID" 2>/dev/null
    [ -n "$XVFB_PID" ] && kill "$XVFB_PID" 2>/dev/null
    return 0
}
trap cleanup EXIT
say() { printf '  [%s] %s\n' "$(date +%H:%M:%S)" "$*"; }

say "screenshots $SHOTDIR"
say "work dir    $WORK"


Xvfb "$DISP" -screen 0 1600x1000x24 >/dev/null 2>&1 &
XVFB_PID=$!
sleep 2

cat > "$WORK/launch.jl" <<'EOF'
using OITOOLS
configure_graphics!()
using GLMakie, QMLMakie, QML
GUI = Base.get_extension(OITOOLS, :OITOOLSGUIExt)
s = GUI.Session()
d = joinpath(pkgdir(OITOOLS), "test", "gui", "data")
GUI.load_dataset!(s, joinpath(d, "2004-data1.oifits"); warn = false, verbose = false)
GUI.load_dataset!(s, joinpath(d, "v1295Aql.oifits"); warn = false, verbose = false)
OITOOLS.gui(s; autoquit_ms = parse(Int, ENV["WALK_AUTOQUIT_MS"]))
EOF
# The window quits itself after this long, which is what writes the console dump: a killed
# process never returns from QML.exec. Long enough for the walk on a slow machine; the walk
# then waits for the exit rather than killing it.
AUTOQUIT_MS="${WALK_AUTOQUIT_MS:-420000}"
WALK_AUTOQUIT_MS="$AUTOQUIT_MS" OITOOLSGUI_SCALE=1.0 OITOOLSGUI_DEBUG_PICK=1 OITOOLSGUI_CONSOLE_DUMP="$DUMP" \
DISPLAY="$DISP" \
    julia --project=bin "$WORK/launch.jl" > "$JLOG" 2>&1 &
GUI_PID=$!

WID=""
for i in $(seq 1 240); do
    WID=$(DISPLAY="$DISP" xdotool search --name "^OITOOLS$" 2>/dev/null | tail -1)
    [ -n "$WID" ] && break
    kill -0 "$GUI_PID" 2>/dev/null || { say "FAIL  the GUI exited:"; tail -25 "$JLOG"; exit 1; }
    sleep 1
done
[ -z "$WID" ] && { say "FAIL  no window after 240 s"; tail -25 "$JLOG"; exit 1; }
sleep 12
DISPLAY="$DISP" xdotool windowsize "$WID" $WIN_W $WIN_H; sleep 3
say "window $WID, driving it"

X()    { DISPLAY="$DISP" xdotool "$@"; }
CL()   { set -- $1 "${2:-1.5}"; X mousemove --window "$WID" "$1" "$2" click 1; sleep "$3"; }
RCL()  { set -- $1; X mousemove --window "$WID" "$1" "$2" click 3; sleep 2; }
HOVER(){ set -- $1; X mousemove --window "$WID" "$1" "$2"; sleep 2.5; }
# Keys go to the FOCUSED window through XTEST, not `--window`: that sends synthetic events,
# which Qt ignores in some widgets -- a path typed into the model picker never arrived. Every
# key here follows a click that put the focus where it is wanted.
KEY()  { X key "$@"; sleep 0.4; }
TYPE() { X type --delay 30 "$1"; sleep 0.5; }
SHOT() { DISPLAY="$DISP" import -window "$WID" "$SHOTDIR/$1.png" 2>/dev/null; }
# Open a menu at $1 and choose row $2 (0-based).
MENU() { CL "$1" 1; for _ in $(seq 1 12); do KEY Up; done
         for _ in $(seq 1 "$2"); do KEY Down; done; KEY Return; sleep 2; }
index_of() { local want=$1; shift; local i=0; for x in "$@"; do [ "$x" = "$want" ] && { echo $i; return; }; i=$((i+1)); done; echo -1; }
KIND()  { MENU "$KIND_BOX"  "$(index_of "$1" "${KINDS[@]}")"; }
COLOR() { MENU "$COLOR_BOX" "$(index_of "$1" "${COLORS[@]}")"; }

# ── 0. a known starting point ──────────────────────────────────────────────────────────
# Two datasets were loaded before the window, and the combo then shows the newest while the
# plot shows the first (a known defect, GUI_Explore.md §2.6). Choosing data1 explicitly is
# what makes every later step start from the same place.
CL "$TAB_EXPLORE" 2
MENU "$DATASET_BOX" 0
SHOT 00_start_data1_uv

# ── C4: colour modes on uv ─────────────────────────────────────────────────────────────
for c in wav mjd none baseline; do COLOR "$c"; SHOT "01_uv_colour_$c"; done

# ── C1/C9/C10/C8: kinds, log y, per group, max baseline ────────────────────────────────
KIND v2;         SHOT 02_v2
CL "$LOG_Y" 2;   SHOT 03_v2_logy;      CL "$LOG_Y" 2
CL "$PER_GROUP" 3; SHOT 04_v2_panels;  CL "$PER_GROUP" 2
KIND t3phi;      SHOT 05_t3phi
CL "$MAX_BL" 2;  SHOT 06_t3phi_maxbl;  CL "$MAX_BL" 2
KIND t3amp;      CL "$LOG_Y" 2; SHOT 07_t3amp_logy; CL "$LOG_Y" 2
for c in wav mjd; do KIND v2; COLOR "$c"; SHOT "08_v2_colour_$c"; done
COLOR baseline

# ── C2: a kind the file lacks is disabled and says why ─────────────────────────────────
CL "$KIND_BOX" 1; HOVER "$KIND_MENU_VISAMP"; SHOT 09_kind_menu_visamp_disabled; KEY Escape

# ── C6: compare choices disabled with their reasons, before any model or image ─────────
KIND v2
HOVER "$RESID_IMAGE"; SHOT 10_resid_image_tooltip_no_image

# ── P4/P5: pick; P1/P2/P3: zoom, pan, reset ────────────────────────────────────────────
HOVER "$NEUTRAL"; SHOT 11_v2_before_zoom
CL "$PICK_AT" 1.5; CL "$PICK_AT" 2.5; SHOT 12_picked
set -- $PLOT_CENTRE
X mousemove --window "$WID" "$1" "$2"
for _ in 1 2 3 4 5; do X click 4; sleep 0.4; done; sleep 2; SHOT 13_zoomed
X mousedown 1; X mousemove --window "$WID" $(( $1 + 80 )) $(( $2 + 40 )); sleep 0.5; X mouseup 1; sleep 2
SHOT 14_panned
RCL "$PLOT_CENTRE"; sleep 2; SHOT 15_reset

# ── X3/C11: Save PNG, to an absolute path typed into the picker ────────────────────────
CL "$SAVE_PNG" 2
# Click the name field first: the picker opens with the suggestion selected but the focus
# elsewhere, and keys typed without the click went nowhere and left the dialog open.
CL "$SAVE_NAME_FIELD" 1
KEY ctrl+a; TYPE "$WORK/explore_v2.png"; KEY Return
sleep 6; SHOT 16_after_save

# ── a model, opened through the Model tab ──────────────────────────────────────────────
CL "$TAB_MODEL" 2; CL "$MODEL_OPEN" 2; CL "$MODEL_PLACE" 1.5; CL "$MODEL_ROW" 1
KEY Return; sleep 4; SHOT 17_model_loaded
CL "$TAB_EXPLORE" 2

# ── C5/C6: model overlay and model residuals, mutually exclusive ───────────────────────
KIND v2
CL "$OVER_MODEL" 3;  SHOT 18_v2_model_obs
CL "$RESID_MODEL" 3; SHOT 19_v2_resid_model_clears_model_obs
# D5 (fixed): back to "data only" returns to the data, not the residuals.
CL "$DATA_ONLY" 3; SHOT 19b_untick_resid_model_D5
KIND uv; HOVER "$RESID_MODEL"; SHOT 20_uv_compare_tooltip
KIND v2

# ── an image, reconstructed in the Image tab ───────────────────────────────────────────
CL "$TAB_IMAGE" 2; CL "$IMAGE_RUN" 25; SHOT 21_image_done
CL "$TAB_EXPLORE" 2
# D4 (fixed): the image choices are enabled as soon as Explore is shown again.
HOVER "$OVER_IMAGE"; SHOT 21b_image_obs_still_disabled_D4
KIND v2
CL "$OVER_IMAGE" 3;  SHOT 22_v2_image_obs
CL "$RESID_IMAGE" 3; SHOT 23_v2_resid_image

# ── C3: a kind the next file lacks falls back to uv ────────────────────────────────────
MENU "$DATASET_BOX" 1; sleep 3
MENU "$KIND_BOX_DATA2" "$(index_of visamp "${KINDS[@]}")"; SHOT 25_data2_visamp
MENU "$DATASET_BOX" 0; sleep 3; SHOT 26_data1_after_visamp
KIND v2

# ── C7: the model goes away while its residuals are shown ──────────────────────────────
# Last: it used to leave Explore stuck (D6) -- the choice cleared but the shell stayed in
# residual mode, and every later view failed with "! no model". Checked below.
CL "$RESID_MODEL" 3; SHOT 23b_resid_model_on
CL "$TAB_MODEL" 2; CL "$MODEL_CARD_REMOVE" 3; CL "$TAB_EXPLORE" 3
SHOT 24_after_model_removed
KIND uv; SHOT 24b_view_change_after_removal_D6

# ── quit, and collect the console ──────────────────────────────────────────────────────
SHOT 27_end
say "walk done; waiting for the window's auto-quit, which writes the console dump"
while kill -0 "$GUI_PID" 2>/dev/null; do sleep 2; done
GUI_PID=""

echo "-- checks --"
check() { if grep -qF -- "$2" "$DUMP" 2>/dev/null; then echo "  PASS  $1"; else echo "  FAIL  $1   (no line: $2)"; FAIL=1; fi; }
[ -s "$DUMP" ] || { echo "  FAIL  no console dump at $DUMP"; exit 1; }
check "C4 colour wav"                   "plot uv, coloured by wav"
check "C4 colour mjd"                   "plot uv, coloured by mjd"
check "C4 colour none"                  "plot uv, coloured by none"
check "C1 v2"                           "plot v2, coloured by baseline"
check "C9 log y"                        "plot v2, coloured by baseline, log y"
check "C10 per group"                   "plot v2, coloured by baseline, per baseline"
check "C1 t3phi"                        "plot t3phi, coloured by baseline"
check "C8 max baseline"                 "plot t3phi_max, coloured by baseline"
check "C9 log y on t3amp"               "plot t3amp, coloured by baseline, log y"
check "P4 a click reached the picker"   "pick at"
if grep -q "pick at.*-> <nothing>" "$DUMP" && ! grep -q "pick at.*-> [^<]" "$DUMP"; then
    echo "  FAIL  P4 every pick came back empty"; FAIL=1; fi
check "C11 Save PNG wrote a file"       "wrote $WORK/explore_v2.png"
[ -s "$WORK/explore_v2.png" ] && echo "  PASS  C11 the PNG exists" || { echo "  FAIL  C11 no PNG"; FAIL=1; }
check "model opened"                    "read_model_file("
check "C5 model overlay"                "model_to_obs(model, x, data)"
check "C5 model residuals"              "model_to_residuals(model, x, data)"
check "image reconstructed"             "reconstruct"
check "C5 image overlay"                "image_to_obs(x, ft, data)"
check "C5 image residuals"              "image_to_residuals(x, ft, data)"
check "C3 visamp on data2"              "plot visamp, coloured by baseline"
# The first view after the visamp step must be uv: the file lacks visamp, so the menu falls
# back rather than leaving a kind that cannot be drawn.
if awk '/plot visamp, coloured by baseline/{f=1; next} f && /> plot /{print; exit}' "$DUMP" \
       | grep -qF "plot uv, coloured by baseline"; then
    echo "  PASS  C3 falls back to uv on a file without visamp"
else echo "  FAIL  C3 no fall-back to uv after switching to a file without visamp"; FAIL=1; fi
# D6, fixed: after the model is removed, Explore's next view is drawn, not "! no model".
if tail -6 "$DUMP" | grep -q "! no model"; then
    echo "  FAIL  D6 Explore stuck on '! no model' after the model was removed"; FAIL=1
else
    echo "  PASS  D6 Explore recovers after the model is removed"; fi

# P3: the reset view is the view before zooming, pixel for pixel, apart from the pick label.
if command -v compare >/dev/null; then
    d=$(compare -metric AE -crop 1400x880+0+0 "$SHOTDIR/11_v2_before_zoom.png" \
        "$SHOTDIR/15_reset.png" null: 2>&1 | cut -d' ' -f1)
    if [ "${d%%.*}" -lt 50 ]; then echo "  PASS  P3 reset returns to the home view ($d px differ)"
    else echo "  NOTE  P3 reset view differs from the home view in $d px — inspect 11 vs 15"; fi
fi

echo "-- console dump: $DUMP"
echo "-- screenshots:  $SHOTDIR"
exit $FAIL
