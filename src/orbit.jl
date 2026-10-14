# Binary orbits, as PMOIRED defines them (`pmoired/oimodels.py`, `_orbit`), so PMOIRED's keys
# (`c,x = c,y = "orbit"`, `c,orb P`, `MJD0`, `e`, `omega`, `OMEGA`, `incl`, `a`, …) mean here
# what they mean there:
#
#   M = 2π(t − MJD0)/P, Kepler E − e sin E = M, ν the true anomaly, r = a(1 − e²)/(1 + e cos ν);
#   the position (x East, y North, mas) is (r cos ν, r sin ν) rotated by −ω, by i + π about the
#   x axis, and by Ω − π/2 — Ω from North; ω is the secondary's, relative to the primary
#   (orbitize!'s convention, 180° from Gallenne et al. 2023's);
#   VA = γ − Ka[cos(ω + ν) + e cos ω], VB = γ + Kb[cos(ω + ν) + e cos ω] (km/s, positive receding).
#
# `orbit_x` / `orbit_y` are called from the model's expressions (`_preprocess_sugar` writes
# them for a component on an orbit), point by point at each point's MJD, and differentiated
# there on dual numbers: Kepler's equation is solved on the values, and one Newton step taken on
# the duals gives its exact derivatives (the implicit dE = (dM + sin E de)/(1 − e cos E)).

# The plain value of a number, through any nesting of dual numbers.
_plain(x::ForwardDiff.Dual) = _plain(ForwardDiff.value(x))
_plain(x::Real) = x

"""
    kepler_E(M, e) -> E

The eccentric anomaly solving E − e sin E = M (radians, 0 ≤ e < 1), by Newton's method on the
values (in Float64), then one more Newton step in the arguments' own type: on dual numbers it
carries the derivatives exactly, with the same value as on plain numbers.
"""
function kepler_E(M::Real, e::Real)
    m, ε = Float64(_plain(M)), Float64(_plain(e))
    0 <= ε < 1 || throw(ArgumentError("an orbit's eccentricity is in [0, 1), not $ε"))
    # M folded into [0, π] (E is odd in M, periodic in 2π) and E₀ = min(π, M/(1 − e)): f is
    # convex there and f(E₀) ≥ 0, so Newton descends onto the root without overshoot for any
    # e < 1 (ROTIR's kepler_E). Stopped at 4 eps, or when the step stops shrinking — near e = 1
    # rounding in f, amplified by 1/(1 − e cos E), can hold |ΔE| one ulp above any tolerance.
    mr = rem2pi(m, RoundNearest)
    mf = abs(mr)
    E = ε == 0 ? mf : min(Float64(π), mf / (1 - ε))
    prev = Inf
    for _ in 1:60
        ε == 0 && break
        dE = (E - ε * sin(E) - mf) / (1 - ε * cos(E))
        E -= dE
        (abs(dE) < 4eps() || abs(dE) >= prev) && break
        prev = abs(dE)
    end
    # One more Newton step in the arguments' own type: on dual numbers it carries the derivatives,
    # with the same value as on plain numbers.
    E = oftype(float(_plain(M) + _plain(e)), copysign(E, mr) + (m - mr))
    return E - (E - e * sin(E) - M) / (1 - e * cos(E))
end

# The true anomaly and the separation (in the units of `a`) at time t (MJD).
function _orbit_nu_r(t, P, MJD0, e, a)
    E = kepler_E(π * (2 * (t - MJD0) / P), e)
    ν = 2 * atan(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2))
    return ν, a * (1 - e * cos(E))
end

# PMOIRED's rotations: −ω, then i + π about x, then Ω − π/2 ("ref is North").
function _orbit_xy(t, P, MJD0, e, omega, OMEGA, incl, a)
    ν, r = _orbit_nu_r(t, P, MJD0, e, a)
    ω, Ω, i = deg2rad(omega), deg2rad(OMEGA), deg2rad(incl)
    x0, y0 = r * cos(ν), r * sin(ν)
    x1 = x0 * cos(ω) - y0 * sin(ω)
    y1 = x0 * sin(ω) + y0 * cos(ω)
    y2 = y1 * cos(i + π)
    θ = (2Ω - π) / 2
    return x1 * cos(θ) + y2 * sin(θ), -x1 * sin(θ) + y2 * cos(θ)
end

"""
    orbit_x(t, P, MJD0, e, omega, OMEGA, incl, a)
    orbit_y(t, P, MJD0, e, omega, OMEGA, incl, a)

The secondary's offset from the primary (mas, x East, y North) at time `t` (MJD), for a period
`P` (days), periastron `MJD0`, eccentricity `e`, `omega`, `OMEGA`, `incl` (degrees) and
semi-major axis `a` (mas) — PMOIRED's `_orbit`. What a component on an orbit (`c,x = c,y =
"orbit"`) has as its position.
"""
orbit_x(t, P, MJD0, e, omega, OMEGA, incl, a) = _orbit_xy(t, P, MJD0, e, omega, OMEGA, incl, a)[1]
"The northward offset of [`orbit_x`](@ref)'s orbit (mas)."
orbit_y(t, P, MJD0, e, omega, OMEGA, incl, a) = _orbit_xy(t, P, MJD0, e, omega, OMEGA, incl, a)[2]

"""
    orbit_rv(t, star, P, MJD0, e, omega, K, gamma; gamma_d = 0) -> km/s

The radial velocity of the primary (`star = :a`, semi-amplitude `K` = Ka) or the secondary
(`:b`, Kb) at time `t` (MJD), PMOIRED's: VA = γ − Ka[cos(ω + ν) + e cos ω],
VB = γ + Kb[cos(ω + ν) + e cos ω], plus `gamma_d`·(t − MJD0).
"""
function orbit_rv(t, star::Symbol, P, MJD0, e, omega, K, gamma; gamma_d = 0)
    ν, _ = _orbit_nu_r(t, P, MJD0, e, one(e))
    ω = deg2rad(omega)
    s = star === :a ? -1 : star === :b ? 1 :
        throw(ArgumentError("star is :a (primary) or :b (secondary), not $(repr(star))"))
    return gamma + gamma_d * (t - MJD0) + s * K * (cos(ω + ν) + e * cos(ω))
end

"""
    orbit_semiamplitudes(P, e, incl, M, q) -> (Ka, Kb)

The velocity semi-amplitudes (km/s) from the total mass `M` (M☉) and the mass ratio
`q = Mb/Ma`, by Kepler's third law (PMOIRED's): K = 2π a sin i / (P √(1 − e²)) with a in km,
Ka = K/(1 + 1/q), Kb = K/(1 + q).
"""
function orbit_semiamplitudes(P, e, incl, M, q)
    c(x) = oftype(float(_plain(P)), x)
    a_au = cbrt(M * (P / c(365.25))^2)
    K = π * (2 * a_au * c(1.495978707e8) * sin(deg2rad(incl)) / (P * 86400 * sqrt(1 - e^2)))
    return K / (1 + 1 / q), K / (1 + q)
end

"""
    orbit_a(P, M, plx) -> mas

The apparent semi-major axis from the total mass `M` (M☉), the parallax `plx` (mas) and the
period `P` (days), by Kepler's third law: a = plx · (M (P/365.25)²)^(1/3).
"""
orbit_a(P, M, plx) = plx * cbrt(M * (P / oftype(float(_plain(P)), 365.25))^2)

# The element names a component on an orbit may carry (`c,orb <name>`).
const _ORBIT_POSITION_ELEMENTS = ("P", "MJD0", "e", "omega", "OMEGA", "incl")
const _ORBIT_ELEMENTS = ("P", "MJD0", "e", "omega", "OMEGA", "incl", "a", "M", "plx", "q",
                         "Ka", "Kb", "gamma")

# Default bounds of the elements a fit frees (`default_bounds`).
const _ORBIT_BOUNDS = Dict("P" => (0.0, Inf), "e" => (0.0, 0.99), "omega" => (-360.0, 360.0),
                           "OMEGA" => (-360.0, 360.0), "incl" => (0.0, 180.0), "a" => (0.0, Inf),
                           "M" => (0.0, Inf), "plx" => (0.0, Inf), "q" => (0.0, Inf),
                           "Ka" => (0.0, Inf), "Kb" => (0.0, Inf))

# Write `c,x` / `c,y` = "orbit" as the expressions that compute them. Errors name what is missing.
function _expand_orbit!(d::Dict, cn::AbstractString)
    kx, ky = "$cn,x", "$cn,y"
    ox, oy = get(d, kx, nothing) == "orbit", get(d, ky, nothing) == "orbit"
    (ox || oy) || return d
    ox && oy || throw(ArgumentError("$cn: x and y are both \"orbit\" or neither"))
    el(n) = "\$$cn,orb $n"
    missing_ = [n for n in _ORBIT_POSITION_ELEMENTS if !haskey(d, "$cn,orb $n")]
    isempty(missing_) || throw(ArgumentError(
        "$cn is on an orbit but has no " * join(("$cn,orb $n" for n in missing_), ", ")))
    a = haskey(d, "$cn,orb a") ? el("a") :
        haskey(d, "$cn,orb M") && haskey(d, "$cn,orb plx") ?
            "orbit_a($(el("P")), $(el("M")), $(el("plx")))" :
        throw(ArgumentError("$cn is on an orbit but has no semi-major axis: give $cn,orb a, or " *
                            "$cn,orb M and $cn,orb plx"))
    args = join(["\$MJD", el("P"), el("MJD0"), el("e"), el("omega"), el("OMEGA"), el("incl"), a], ", ")
    d[kx] = "orbit_x($args)"
    d[ky] = "orbit_y($args)"
    return d
end
