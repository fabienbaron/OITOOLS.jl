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
# Plot text is deliberately NOT scaled. Makie font sizes stay at the 12 pt / 8 pt / 10 pt
# measured from oiplot, so a figure on screen is the figure in the paper, and the port harness
# keeps asserting against those numbers unchanged. The cost is that on a HiDPI panel the plot
# labels read smaller than the chrome around them.

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
