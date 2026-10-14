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
    ring = Makie.scatter!(ax, best; color = (:white, 0.0), strokecolor = CHI2_BEST_COLOR,
                          strokewidth = 2, markersize = 14, marker = :circle)
    dot  = Makie.scatter!(ax, best; color = CHI2_BEST_COLOR, markersize = 3)

    # A multi-start search: each start joined to the minimum it reached, the starts as `+`
    # (set aside ones as red `×`), the minima coloured by χ²ᵣ. Empty for the other maps.
    mseg   = Makie.Observable(Makie.Point2f[])
    mstart = Makie.Observable(Makie.Point2f[])
    mbad   = Makie.Observable(Makie.Point2f[])
    Makie.linesegments!(ax, mseg; color = (:gray40, 0.35), linewidth = 0.8)
    Makie.scatter!(ax, mstart; marker = '+', markersize = 9, color = (:gray30, 0.6))
    Makie.scatter!(ax, mbad; marker = :xcross, markersize = 7, color = (:firebrick, 0.6))
    # positions and colours set together (`Makie.update!`), never one length ahead of the other
    mmin = Makie.scatter!(ax, Makie.Point2f[]; color = Float32[], colormap = :viridis,
                          colorrange = (0.0f0, 1.0f0), markersize = 9, strokecolor = :black,
                          strokewidth = 0.5)

    cblim = Makie.Observable((0.0f0, 1.0f0))
    cbar  = Makie.Colorbar(fig[1, 2]; colormap = Makie.Reverse(:gist_heat), limits = cblim,
                           label = "log₁₀(Δχ² + 1)")
    Makie.colgap!(fig.layout, 1, 8)

    ax.title  = ""
    ax.xlabel = ""
    ax.ylabel = ""

    return (; figure = fig, axis = ax, xs, ys, zs, raw, best, heatmap = hm, contour = ct,
              cbarlimits = cblim, colorbar = cbar, mseg, mstart, mbad, mmin, ring, dot)
end

"""
    update_chi2_map!(c, m::Chi2Map)

Draw one map. Allocates only the arrays handed to the Observables.

The axes carry the parameter names, and the title the best point, because a χ² map with
unlabelled axes is unreadable and a map without its minimum stated is unquotable.
"""
function update_chi2_map!(c, m::Chi2Map)
    _clear_minima!(c)
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

    # A companion map before this one leaves the sky's aspect behind; a parameter grid has no
    # equal units on its two axes. (`xlims!` below, low to high, unreverses the x axis.)
    c.axis.aspect = nothing
    c.heatmap.colormap = Makie.Reverse(:gist_heat)
    c.colorbar.colormap = Makie.Reverse(:gist_heat)
    c.colorbar.label = "log₁₀(Δχ² + 1)"
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

"""
    update_chi2_map!(c, s::CompanionSearch; show = :chi2r)

Draw a companion search into the same chart: χ²ᵣ at every node (`show = :chi2r`), or the
companion flux fitted there (`show = :flux`), with the nodes that were not searched (inside
`rmin`, outside `rmax`, masked) left blank and the best node ringed.

On the sky, so East is to the LEFT (`x` grows leftwards) and the axes share a scale, the
convention `imdisp` and `test_image_orientation.jl` pin. The map shows χ²ᵣ itself rather than
log Δχ²: CANDID and PMOIRED show it that way, and the colours are compared against the χ²ᵣ of
no companion, which the title gives. No contours: their levels are for a parameter grid's Δχ².
"""
function update_chi2_map!(c, s::CompanionSearch; show::Symbol = :chi2r)
    z = show === :flux ? companion_flux_map(s) : s.chi2r
    b = s.best
    title = @sprintf("best %.2f mas at PA %.1f°, %s,f = %.4g, χ²ᵣ %.4g (%.4g without): %.1fσ",
                     b.sep, b.pa, s.companion, b.f, b.chi2r, s.chi2r_null, s.nsigma)
    return _draw_companion_map!(c, s.x, s.y, z, [Makie.Point2f(b.x, b.y)], title,
                                show === :flux ? s.companion * ",f" : "χ²ᵣ")
end

"""
    update_companion_progress!(c, g, companion) -> (done, total)

Draw a companion search while it runs: `g` is what [`companion_search`](@ref) hands its
`on_grid` hook, whose χ²ᵣ array its threads are filling. Nodes not yet fitted are blank, so the
axes and the grid's extent appear as soon as the search starts and the map fills in. The best
node so far is ringed.
"""
function update_companion_progress!(c, g, companion::AbstractString)
    z = copy(g.chi2r)                     # a snapshot: the search's threads keep writing
    done = count(isfinite, z)
    best = Makie.Point2f[]
    if done > 0
        ib = argmin(v -> isfinite(z[v]) ? z[v] : Inf, CartesianIndices(z))
        best = [Makie.Point2f(g.x[ib[1]], g.y[ib[2]])]
    end
    _draw_companion_map!(c, g.x, g.y, z, best,
                         @sprintf("searching for %s: %d of %d nodes", companion, done, g.nnodes),
                         "χ²ᵣ")
    return done, g.nnodes
end

# The sky's conventions for a companion map: East to the LEFT, equal scales on the two axes,
# blank (NaN) where no node was fitted.
function _draw_companion_map!(c, x, y, z, best, title, label)
    _clear_minima!(c)
    zz = Float32.(z)
    c.xs[]  = Float32.(x)
    c.ys[]  = Float32.(y)
    c.raw[] = zeros(Float32, size(zz))
    c.zs[]  = zz
    ok = filter(isfinite, zz)
    lo, hi = isempty(ok) ? (0.0f0, 1.0f0) : extrema(ok)
    c.cbarlimits[] = (lo, hi > lo ? hi : lo + 1)
    c.colorbar.label = label
    # Neither end white: the nodes not searched are blank (NaN), and with the grid map's
    # colours the best node was white too, indistinguishable from the blank centre.
    c.heatmap.colormap = :viridis
    c.colorbar.colormap = :viridis
    c.best[] = best
    c.axis.xlabel = "x (mas, East)"
    c.axis.ylabel = "y (mas, North)"
    c.axis.aspect = Makie.DataAspect()
    c.axis.title = title
    # High to low: `xlims!` sets `xreversed` from the order it is given, so this is what puts
    # East on the left.
    Makie.xlims!(c.axis, last(x), first(x))
    Makie.ylims!(c.axis, first(y), last(y))
    return c
end


# The multi-start plots emptied, for the maps that do not use them.
function _clear_minima!(c)
    for o in (c.mseg, c.mstart, c.mbad)
        o[] = Makie.Point2f[]
    end
    Makie.update!(c.mmin, Makie.Point2f[]; color = Float32[])
    c.ring.strokecolor = CHI2_BEST_COLOR; c.dot.color = CHI2_BEST_COLOR
    return c
end

"""
    update_chi2_map!(c, r::MultiStartResult; p1, p2, selected = 1)

Draw a multi-start search over the parameters `p1` and `p2` (default the first two explored,
or free): every start joined to the minimum its fit reached, the starts of fits set aside as
red crosses, the minima coloured by log₁₀ χ²ᵣ, and minimum `selected` ringed. No surface: the
search samples the parameter space at its starts, it does not map it.
"""
function update_chi2_map!(c, r::MultiStartResult; p1::AbstractString = "",
                          p2::AbstractString = "", selected::Integer = 1)
    # the free parameters, then the held ones (a fit's held values are its start's)
    names = [r.list_free_params; r.held]
    pick = [r.explored; filter(n -> n ∉ r.explored, names)]
    a = isempty(p1) || p1 ∉ names ? pick[1] : String(p1)
    b = isempty(p2) || p2 ∉ names ? pick[min(2, length(pick))] : String(p2)
    ia, ib = findfirst(==(a), names), findfirst(==(b), names)
    nf = length(r.list_free_params)
    pt(v) = Makie.Point2f(v[ia], v[ib])
    ends(s) = [r.fits[s].x_opt; r.starts[(nf + 1):end, s]]
    owner = Dict{Int,Int}()
    for (k, m) in enumerate(r.minima), s in m.starts
        owner[s] = k
    end
    seg = Makie.Point2f[]; st = Makie.Point2f[]; bad = Makie.Point2f[]
    for s in axes(r.starts, 2)
        p0 = pt(view(r.starts, :, s))
        if haskey(owner, s)
            push!(st, p0)
            push!(seg, p0, pt(ends(s)))
        else
            push!(bad, p0)
        end
    end
    mins = [pt(ends(m.best)) for m in r.minima]
    col = Float32[log10(max(m.fit.chi2r, 1e-12)) for m in r.minima]
    lo, hi = isempty(col) ? (0.0f0, 1.0f0) : extrema(col)
    # nothing behind the points: the surface's plots get one blank cell
    c.xs[] = Float32[0, 1]; c.ys[] = Float32[0, 1]
    c.raw[] = zeros(Float32, 2, 2); c.zs[] = fill(NaN32, 2, 2)
    crange = (lo, hi > lo ? hi : lo + 1)
    c.mseg[] = seg; c.mstart[] = st; c.mbad[] = bad
    Makie.update!(c.mmin, mins; color = col, colorrange = crange)
    c.best[] = isempty(mins) ? Makie.Point2f[] : [mins[clamp(selected, 1, length(mins))]]
    c.ring.strokecolor = :black; c.dot.color = :black    # white is for the maps' dark colours
    c.cbarlimits[] = crange
    c.colorbar.colormap = :viridis
    c.colorbar.label = "log₁₀ χ²ᵣ of the minimum"
    c.axis.aspect = nothing
    c.axis.xreversed = false              # a companion map before this one leaves East left
    c.axis.xlabel = a; c.axis.ylabel = b
    k = clamp(selected, 1, max(length(r.minima), 1))
    c.axis.title = isempty(r.minima) ? "no minimum" :
        @sprintf("%d starts, %d minima; minimum %d: χ²ᵣ %.4g from %d start(s)",
                 size(r.starts, 2), length(r.minima), k, r.minima[k].fit.chi2r,
                 length(r.minima[k].starts))
    Makie.autolimits!(c.axis)
    return (a, b)
end
