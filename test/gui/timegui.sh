#!/usr/bin/env bash
# Time launch -> window for the GUI, under whatever environment is given.
#
#     test/gui/timegui.sh QT_QPA_PLATFORM=xcb
#     test/gui/timegui.sh OITOOLSGUI_GLFW_X11=1
#
# Set SYSIMAGE to time a PackageCompiler build instead of a cold session, which is the
# comparison bin/compile.jl exists to be judged on:
#
#     SYSIMAGE=bin/oitoolsgui.so test/gui/timegui.sh
#
# Reports the time, or says why there is none: a crash, or a window that never appeared.
# Run each twice and take the second — the first pays the glyph-atlas cache and any
# extension precompilation.
set -uo pipefail
cd "$(dirname "$0")/../.."
LOG=$(mktemp /tmp/timegui-XXXX.log)
DATA=${DATA:-demos/data/BC2004/2004-data1.oifits}

t0=$(date +%s.%N)
SYSARG=()
[ -n "${SYSIMAGE:-}" ] && SYSARG=(--sysimage "$SYSIMAGE")
env "$@" julia --project=bin "${SYSARG[@]}" bin/oitoolsgui.jl "$DATA" > "$LOG" 2>&1 &
pid=$!

win=""
for i in $(seq 1 600); do                      # 60 s ceiling, not forever
    if ! kill -0 "$pid" 2>/dev/null; then
        wait "$pid"; rc=$?
        echo "[$*] the GUI exited before showing a window (status $rc)"
        echo "  last lines of $LOG:"; tail -6 "$LOG" | sed 's/^/    /'
        exit 1
    fi
    win=$(xdotool search --name '^OITOOLS$' 2>/dev/null | tail -1)
    [ -n "$win" ] && break
    sleep 0.1
done
t1=$(date +%s.%N)

if [ -z "$win" ]; then
    echo "[$*] no window after 60 s; the process is still running"
    echo "  last lines of $LOG:"; tail -6 "$LOG" | sed 's/^/    /'
else
    python3 -c "print('[%s] launch -> window: %.1f s' % ('$*', $t1 - $t0))"
fi
kill "$pid" 2>/dev/null
wait "$pid" 2>/dev/null
echo "  log kept at $LOG"
