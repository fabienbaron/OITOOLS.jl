# The Orbit view (Modeling, under the fits table): the relative orbit on the sky, East left, with
# the model's positions at the dataset's epochs, and the radial velocities against phase with the
# model's curves when a table is open. Drawn from the model's own elements (`c,orb …`), resolved
# as the fit resolves them, so an expression such as `$dP + 42.38` shows its value.

# Its canvas, built by `gui()` before the window exists (a plot created afterwards allocates GPU
# buffers with no GL context bound).
const ORBIT_CANVAS = Ref{Any}(nothing)

"""
    build_orbit_view(fig, axsky, axrv) -> NamedTuple

Every plot the Orbit view draws into, empty. **Call before the window is shown.**
"""
function build_orbit_view(fig, axsky, axrv)
    for ax in (axsky, axrv), it in (:scrollzoom, :dragpan, :rectanglezoom, :limitreset)
        Makie.deregister_interaction!(ax, it)
    end
    track  = Makie.Observable(Makie.Point2f[])
    epochs = Makie.Observable(Makie.Point2f[])
    Makie.lines!(axsky, track; color = :gray40, linewidth = 1.4)
    Makie.scatter!(axsky, [Makie.Point2f(0, 0)]; marker = :star5, markersize = 16, color = :black)
    Makie.scatter!(axsky, epochs; color = :firebrick, markersize = 9)
    axsky.xreversed = true
    axsky.aspect = Makie.DataAspect()
    axsky.xlabel = "Δα (mas, East left)"; axsky.ylabel = "Δδ (mas)"; axsky.title = "orbit"

    curve = [Makie.Observable(Makie.Point2f[]) for _ in 1:2]
    dots  = [Makie.Observable(Makie.Point2f[]) for _ in 1:2]
    for (k, c) in enumerate((:royalblue, :darkorange))
        Makie.lines!(axrv, curve[k]; color = c, linewidth = 1.6)
        Makie.scatter!(axrv, dots[k]; color = c, markersize = 7)
    end
    axrv.xlabel = "phase"; axrv.ylabel = "RV (km/s)"; axrv.title = "radial velocities (a blue, b orange)"
    return (; figure = fig, sky = axsky, rv = axrv, track, epochs, curve, dots)
end

"""
    update_orbit_view!(panel, el, mjds, rv, K) -> nothing

Draw the orbit of elements `el` (P, MJD0, e, omega, OMEGA, incl, a), its positions at `mjds`,
and — with a table `rv` and semi-amplitudes `K = (Ka, Kb, gamma)` — the velocities against phase.
"""
function update_orbit_view!(panel, el, mjds, rv, K)
    panel === nothing && return nothing
    P, MJD0 = el.P, el.MJD0
    pos(t) = Makie.Point2f(orbit_x(t, P, MJD0, el.e, el.omega, el.OMEGA, el.incl, el.a),
                           orbit_y(t, P, MJD0, el.e, el.omega, el.OMEGA, el.incl, el.a))
    ts = range(MJD0, MJD0 + P; length = 400)
    panel.track[]  = [pos(t) for t in ts]
    panel.epochs[] = [pos(t) for t in mjds]
    Makie.autolimits!(panel.sky)
    phase(t) = mod((t - MJD0) / P, 1.0)
    if rv === nothing || K === nothing
        for k in 1:2; panel.curve[k][] = Makie.Point2f[]; panel.dots[k][] = Makie.Point2f[]; end
    else
        ph = range(0, 1; length = 300)
        for (k, s) in enumerate((:a, :b))
            Kk = s === :a ? K.Ka : K.Kb
            panel.curve[k][] = [Makie.Point2f(p, orbit_rv(MJD0 + p * P, s, P, MJD0, el.e, el.omega,
                                                          Kk, K.gamma)) for p in ph]
            sel = findall(==(s), rv.star)
            panel.dots[k][] = [Makie.Point2f(phase(rv.mjd[i]), rv.rv[i]) for i in sel]
        end
    end
    Makie.autolimits!(panel.rv)
    return nothing
end
