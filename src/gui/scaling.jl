# UI scaling policy.
#
# The window mixes two sizing systems, and before this file they scaled by different rules:
#
#   font.pointSize   points      Qt converts using the font DPI, so text already grows on a
#                                HiDPI screen without anyone asking it to
#   width/spacing    logical px  scaled by the device pixel ratio only, so containers did NOT
#                                grow with the text they hold
#
# Left unscaled, text outgrows the things that contain it: a clipped dataset combo, a cramped
# log panel, buttons whose labels do not fit. So the pixel quantities are scaled by the same
# factor Qt applies to the text.
#
# DETECTION IS NOT HERE. `Screen.logicalPixelDensity` is a live QML property -- it is correct
# before any window exists and it updates when the window is dragged to a monitor with a
# different scale, neither of which a value computed once in Julia could manage. So QML works
# the factor out and this file supplies only the override, which keeps the policy (parsing,
# range, what counts as "auto") in plain Julia where it can be tested without a display.
#
# Plot text is not scaled on the STATIC path. `draw!` leaves Makie font sizes at the 12 pt /
# 8 pt / 10 pt measured from oiplot, so a figure saved from a script is the figure in the
# paper, and the port harness keeps asserting against those numbers unchanged.
#
# The live canvas is a different question, because it is read on screen rather than on paper,
# and it follows the UI scale -- see `live_plot_scale`. It has to: a window whose chrome
# tracks the screen while its plot labels do not ends up with 8 pt tick labels beside 14 pt
# buttons, which is how the same code came to look right on a 189 dpi laptop and too large on
# a 92 dpi desktop.

"Environment variable that overrides the detected UI scale."
const UI_SCALE_VAR = "OITOOLSGUI_SCALE"

"Sentinel handed to QML meaning \"work it out from the screen\"."
const UI_SCALE_AUTO = 0.0

# Below 0.5 the window is unreadable; above 4.0 it cannot fit a 4K panel. A value outside the
# range is far more likely to be a typo (a stray percentage, say) than an intention, so it is
# clamped and reported rather than obeyed.
const UI_SCALE_MIN = 0.5
const UI_SCALE_MAX = 4.0

"""
    ui_scale_override(; verbose = true) -> (; scale, auto, reason)

Read `OITOOLSGUI_SCALE` and turn it into the value QML expects.

`scale` is `UI_SCALE_AUTO` (0.0) when QML should detect the scale itself, which is the
default and the case for every input that does not name a usable number. `auto` says which
of the two happened, and `reason` is the one-line explanation, shown at startup.

    OITOOLSGUI_SCALE=1.5     # everything 50% larger than the screen asks for
    OITOOLSGUI_SCALE=auto    # or unset, or empty: let QML decide

Out-of-range values are clamped to [\$UI_SCALE_MIN, \$UI_SCALE_MAX] with a warning; anything
unparseable warns and falls back to automatic. Neither is an error: a bad scale should not
stop the GUI from opening, because the GUI is how you would fix it.
"""
function ui_scale_override(; verbose::Bool = true)
    raw = strip(get(ENV, UI_SCALE_VAR, ""))
    auto(reason) = (scale = UI_SCALE_AUTO, auto = true, reason = reason)

    isempty(raw)                    && return auto("unset; QML scales from the screen")
    lowercase(raw) in ("auto", "0") && return auto("$UI_SCALE_VAR=$raw; QML scales from the screen")

    parsed = tryparse(Float64, raw)
    if parsed === nothing || !isfinite(parsed)
        verbose && @warn "$UI_SCALE_VAR is not a number, ignoring it" value = raw
        return auto("$UI_SCALE_VAR=$raw is not a number; falling back to automatic")
    end

    if parsed < UI_SCALE_MIN || parsed > UI_SCALE_MAX
        clamped = clamp(parsed, UI_SCALE_MIN, UI_SCALE_MAX)
        verbose && @warn "$UI_SCALE_VAR is outside the usable range, clamping" requested =
            parsed used = clamped range = (UI_SCALE_MIN, UI_SCALE_MAX)
        return (scale = clamped, auto = false,
                reason = "$UI_SCALE_VAR=$parsed clamped to $clamped")
    end

    return (scale = parsed, auto = false, reason = "$UI_SCALE_VAR=$parsed")
end

"""
The UI scale QML settled on, pushed back from the window.

QML owns the detection (`Screen` is live and Julia's copy could not follow a move to another
monitor), so the value has to travel back for anything drawn in Julia to match it. Until the
window reports in, this holds `UI_SCALE_REFERENCE` -- `build_canvas` runs before any window
exists, and every later `update_canvas!` restyles from the live value anyway.
"""
const UI_SCALE_LIVE = Ref(1.25)

"""
Physical density of the screen, in dpi, as QML read it. Zero until the window reports in.

Plot typography is sized from THIS rather than from the UI scale, so the two can be judged
and retuned separately -- they were tied together at first, and every adjustment to one then
moved the other. Both happen to follow the same square root of density; they simply carry
their own anchors.
"""
const SCREEN_DPI_LIVE = Ref(0.0)

"""
Plot scale, anchored on the same desktop the UI scale is: 1.19 at `PLOT_REF_DPI`.

The exponent is 0.5, which carries it to 1.7 at 189 dpi. Chrome uses ~0.9 -- see
`src/gui/qml/Main.qml`, which owns that half.
"""
const PLOT_REF_DPI          = 92.6
const PLOT_SCALE_AT_REF_DPI = 1.19

"Fallback ratio for a screen whose density is unknown, so plots still track the chrome."
const UI_SCALE_REFERENCE       = 0.656
const PLOT_SCALE_AT_REFERENCE  = 1.19

"""
Record what QML measured. Called once, from the window's `Component.onCompleted`.

`dpi` is the physical density and may be 0 where the screen does not report one -- a virtual
display, say -- in which case plot sizing falls back to the UI scale.
"""
function set_ui_scale!(x::Real, dpi::Real = 0.0)
    isfinite(dpi) && dpi > 0 && (SCREEN_DPI_LIVE[] = Float64(dpi))
    (isfinite(x) && x > 0) || return UI_SCALE_LIVE[]
    return UI_SCALE_LIVE[] = clamp(Float64(x), UI_SCALE_MIN, UI_SCALE_MAX)
end

"""
Live overrides from the settings panel. Zero means "no override, compute it".

Separate from the environment variables, which are read once at startup: these are turned by
hand while the window is open, and have to win over both the detected value and the variable.
"""
const PLOT_SCALE_USER  = Ref(0.0)
const MARKER_SIZE_USER = Ref(0.0)

"Set the plot scale by hand. Zero restores the value computed from the screen."
function set_plot_scale!(x::Real)
    PLOT_SCALE_USER[] = (isfinite(x) && x > 0) ? clamp(Float64(x), 0.2, 5.0) : 0.0
    return PLOT_SCALE_USER[]
end

"Set the data-point size by hand, in points. Zero restores each plot's own default."
function set_marker_size!(x::Real)
    MARKER_SIZE_USER[] = (isfinite(x) && x > 0) ? clamp(Float64(x), 1.0, 40.0) : 0.0
    return MARKER_SIZE_USER[]
end
