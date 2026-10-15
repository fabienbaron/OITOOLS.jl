# Spectral lines on a component's flux, as PMOIRED defines them (`pmoired/oimodels.py`,
# `_Fspectrum`), so PMOIRED's keys mean here what they mean there:
#
#   c,line_<n>_wl0     the line's centre (µm)
#   c,line_<n>_f       its peak, added to the continuum (< 0 for an absorption line)
#   and one profile, its width in nm:
#   c,line_<n>_lorentzian   FWHM w: f (w/2)² / ((λ − λ₀)² + (w/2)²)
#   c,line_<n>_gaussian     FWHM w: f exp(−|λ − λ₀|^p / (2 (w/2.35482)^p)), p = c,line_<n>_power (2)
#   c,line_<n>_truncexp     w: f exp(−|λ − λ₀|^q / w^q) redward, exp(−|λ − λ₀|^2q / (w/10)^2q)
#                           blueward, q = c,line_<n>_pow (1)
#   c,line_<n>_trunclor     w: the Lorentzian on one side of λ₀ (red for w > 0, blue for w < 0)
#
# The component's flux is its continuum (`c,f`, or 1 without one) plus its lines, at each
# point's wavelength: `_expand_lines!` writes that sum as the expression `c,ftotal`, which the
# component's flux then is. Each line key is a parameter like any other — fixed, free, or an
# expression (a centre written as a velocity, `2.16612*(1 + $v/2.998e5)`).

"""
    line_lorentzian(wl, wl0, f, fwhm)

A Lorentzian line of peak `f` at `wl0` (µm), FWHM `fwhm` (nm), at wavelength `wl` (µm):
PMOIRED's `line_<n>_lorentzian`.
"""
line_lorentzian(wl, wl0, f, fwhm) = (h = fwhm / 2000; f * h^2 / ((wl - wl0)^2 + h^2))

"""
    line_gaussian(wl, wl0, f, fwhm, power = 2)

A Gaussian line (generalised by `power`) of peak `f` at `wl0` (µm), FWHM `fwhm` (nm):
PMOIRED's `line_<n>_gaussian` with `line_<n>_power`.
"""
line_gaussian(wl, wl0, f, fwhm, power = 2) =
    f * exp(-abs(wl - wl0)^power / (2 * (fwhm / 1000 / 2.35482)^power))

"""
    line_truncexp(wl, wl0, f, width, pow = 1)

PMOIRED's one-sided exponential line (`line_<n>_truncexp`, `line_<n>_pow`): exp(−|Δλ|^pow /
w^pow) redward of `wl0`, and a ten times narrower exp(−|Δλ|^2pow / (w/10)^2pow) blueward; `width`
in nm.
"""
function line_truncexp(wl, wl0, f, width, pow = 1)
    w = width / 1000
    d = abs(wl - wl0)
    return f * (wl >= wl0 ? exp(-d^pow / w^pow) : exp(-d^(2pow) / (w / 10)^(2pow)))
end

"""
    line_trunclor(wl, wl0, f, width)

PMOIRED's one-sided Lorentzian (`line_<n>_trunclor`): the Lorentzian of FWHM |`width`| (nm)
redward of `wl0` for `width > 0`, blueward for `width < 0`.
"""
function line_trunclor(wl, wl0, f, width)
    side = width > 0 ? wl >= wl0 : wl <= wl0
    return side ? line_lorentzian(wl, wl0, f, width) : zero(line_lorentzian(wl, wl0, f, width))
end

# The keys of the lines: `line_<name>_<part>`.
const _LINE_PARTS = ("wl0", "f", "lorentzian", "gaussian", "power", "truncexp", "pow", "trunclor")
const _LINE_KEY = r"^line_(.+)_(wl0|f|lorentzian|gaussian|power|truncexp|pow|trunclor)$"

# The line names of component `cn` in `d`, sorted.
_line_names(d, cn) = sort(unique([m[1] for k in keys(d) if startswith(k, cn * ",")
                                  for m in (match(_LINE_KEY, k[(length(cn) + 2):end]),)
                                  if m !== nothing]))

# Write a component's continuum plus its lines as `cn,ftotal`. Errors name what is missing.
function _expand_lines!(d::Dict, cn::AbstractString)
    names = _line_names(d, cn)
    isempty(names) && return d
    k(n, p) = "$cn,line_$(n)_$p"
    ref(n, p) = "\$" * k(n, p)
    terms = String[]
    for n in names
        for p in ("wl0", "f")
            haskey(d, k(n, p)) || throw(ArgumentError("$cn's line $n has no $(k(n, p))"))
        end
        profiles = [p for p in ("lorentzian", "gaussian", "truncexp", "trunclor") if haskey(d, k(n, p))]
        length(profiles) == 1 || throw(ArgumentError(
            "$cn's line $n needs one profile — $(k(n, "lorentzian")), _gaussian, _truncexp or " *
            "_trunclor" * (isempty(profiles) ? "" : ", not $(join(profiles, " and "))")))
        p = profiles[1]
        args = "\$WL, $(ref(n, "wl0")), $(ref(n, "f")), $(ref(n, p))"
        push!(terms, p == "gaussian" ?
                     "line_gaussian($args, $(haskey(d, k(n, "power")) ? ref(n, "power") : "2"))" :
                     p == "truncexp" ?
                     "line_truncexp($args, $(haskey(d, k(n, "pow")) ? ref(n, "pow") : "1"))" :
                     "line_$p($args)")
    end
    cont = haskey(d, "$cn,f") ? "\$$cn,f" : "1"
    d["$cn,ftotal"] = cont * " + " * join(terms, " + ")
    return d
end
