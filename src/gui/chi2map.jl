# The χ² map, in Makie.
#
# Same arrangement as the Gantt: OITOOLS' `chi2_map` computes the surface and this only draws
# it, so what is on screen is what an exported script reproduces.
#
# What is drawn, and why it is not the raw χ²
# -------------------------------------------
# The heatmap shows log10(Δχ² + 1), not χ². A χ² surface spans orders of magnitude between the
# minimum and the edge of a generous grid -- measured on α Cen A, 1444 at the best point
# against 3.4e5 at the corner -- so a linear map of the raw values is one dark region and one
# bright corner, with every feature near the minimum crushed into a single colour.
# `demos/example_grid_search_limb_darkening.jl` reaches the same conclusion and the same fix.
#
# The contours are drawn on the RAW Δχ², at the two-parameter confidence levels. Those are the
# quantitative part of the figure; the colours are only there to show the shape around them.

"Marker for the grid minimum. Drawn over the contours, so it stays visible inside the 1σ region."
const CHI2_BEST_COLOR = Makie.RGBAf(1.0, 1.0, 1.0, 1.0)

"""
    build_chi2_map(figure, axis) -> NamedTuple

Create every plot the χ² map will ever need. **Call before the window is shown**, for the
reason `build_canvas` documents: after QML is up the GL context belongs to Qt's render thread
and inserting a plot allocates buffers with none bound.

`contour!` is created here with placeholder data. Its levels are an attribute rather than
geometry, so they can be re-set later; what cannot happen later is the plot's creation.
"""
function build_chi2_map(fig, ax)
    # Not zoomable, for the Gantt's reason: the axes are the grid, and panning off it shows
    # nothing, because nothing was computed there.
    for it in (:scrollzoom, :dragpan, :rectanglezoom, :limitreset)
        Makie.deregister_interaction!(ax, it)
    end

    xs   = Makie.Observable(Float32[0, 1])
    ys   = Makie.Observable(Float32[0, 1])
    zs   = Makie.Observable(zeros(Float32, 2, 2))     # log10(Δχ² + 1), what the colours show
    raw  = Makie.Observable(zeros(Float32, 2, 2))     # Δχ², what the contours are drawn on
    best = Makie.Observable(Makie.Point2f[])

    hm = Makie.heatmap!(ax, xs, ys, zs; colormap = Makie.Reverse(:gist_heat))
    ct = Makie.contour!(ax, xs, ys, raw;
                        levels = collect(delta_chi2_levels), color = :black, linewidth = 1.2)
    # The minimum, as a ring rather than a dot: a filled marker hides the very cell it reports.
    Makie.scatter!(ax, best; color = (:white, 0.0), strokecolor = CHI2_BEST_COLOR,
                   strokewidth = 2, markersize = 14, marker = :circle)
    Makie.scatter!(ax, best; color = CHI2_BEST_COLOR, markersize = 3)

    cblim = Makie.Observable((0.0f0, 1.0f0))
    cbar  = Makie.Colorbar(fig[1, 2]; colormap = Makie.Reverse(:gist_heat), limits = cblim,
                           label = "log₁₀(Δχ² + 1)")
    Makie.colgap!(fig.layout, 1, 8)

    ax.title  = ""
    ax.xlabel = ""
    ax.ylabel = ""

    return (; figure = fig, axis = ax, xs, ys, zs, raw, best, heatmap = hm, contour = ct,
              cbarlimits = cblim, colorbar = cbar)
end

"""
    update_chi2_map!(c, m::Chi2Map)

Draw one map. Allocates only the arrays handed to the Observables.

The axes carry the parameter names, and the title the best point, because a χ² map with
unlabelled axes is unreadable and a map without its minimum stated is unquotable.
"""
function update_chi2_map!(c, m::Chi2Map)
    lo = minimum(m.chi2)
    dchi2 = Float32.(m.chi2 .- lo)

    c.xs[]  = Float32.(m.v1)
    c.ys[]  = Float32.(m.v2)
    c.raw[] = dchi2
    z = log10.(dchi2 .+ 1.0f0)
    c.zs[]  = z
    c.cbarlimits[] = (0.0f0, Float32(maximum(z) > 0 ? maximum(z) : 1))

    b = argmin(m.chi2)
    c.best[] = [Makie.Point2f(m.v1[b[1]], m.v2[b[2]])]

    c.axis.xlabel = m.p1
    c.axis.ylabel = m.p2
    c.axis.title  = @sprintf("χ²r = %.4g at %s = %.4g, %s = %.4g",
                             lo / m.ndof, m.p1, m.v1[b[1]], m.p2, m.v2[b[2]])
    # The grid IS the plot: there is nothing outside it to pan to, and Makie's default margin
    # would leave a border of blank axis that looks like unexplored parameter space.
    Makie.xlims!(c.axis, first(m.v1), last(m.v1))
    Makie.ylims!(c.axis, first(m.v2), last(m.v2))
    return c
end
