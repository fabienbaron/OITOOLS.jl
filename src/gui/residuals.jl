# Normalised model residuals, in Makie.
#
# What is drawn, and why it is one axis rather than ten
# ----------------------------------------------------
# `plot_residuals` in oiplot.jl stacks a (data + model, residual) row pair per observable --
# up to ten axes for a dataset carrying all five. That is the right figure to publish and the
# wrong one for a 220 dp panel under the fits table, where ten rows leaves each about 20 px.
#
# So this draws the residuals ONLY, all observables in one axis against baseline, one colour
# each. That is the question the panel is here to answer: are the residuals centred on zero,
# is one observable systematically worse than the others, and is there structure against
# baseline that a χ² of 3.2 does not show. Reading a model curve against its data is what the
# Model visualization tab and `plot_residuals` are for.
#
# The y axis is (model - data)/σ, so the ±1 and ±3 lines mean the same thing for every
# observable and the five series are directly comparable -- which is the whole reason they can
# share an axis at all.

"Colour per observable. Fixed rather than taken from the group palette: these five are always
the same five, and a legend whose colours move between datasets cannot be learned."
const RESIDUAL_COLORS = (v2     = Makie.RGBAf(0.12, 0.35, 0.72, 0.75),
                         t3amp  = Makie.RGBAf(0.85, 0.42, 0.06, 0.75),
                         t3phi  = Makie.RGBAf(0.15, 0.60, 0.25, 0.75),
                         visamp = Makie.RGBAf(0.65, 0.20, 0.65, 0.75),
                         visphi = Makie.RGBAf(0.75, 0.15, 0.20, 0.75))

"Order the series are built and reported in: V² first, closure quantities next, complex last."
const RESIDUAL_KINDS = (:v2, :t3amp, :t3phi, :visamp, :visphi)

"The `OBS_SPECS` entry each residual series takes its x axis from."
const RESIDUAL_SPECS = (v2 = :v2, t3amp = :t3amp, t3phi = :t3phi,
                        visamp = :visamp, visphi = :visphi)

"""
    build_residuals(figure, axis) -> NamedTuple

Create every plot the residual panel will ever need. **Call before the window is shown**, for
the reason `build_canvas` documents: after QML is up the GL context belongs to Qt's render
thread and inserting a plot allocates buffers with none bound.

One scatter per observable, all empty. A dataset without T3amp leaves that series empty rather
than absent, so no plot is created or destroyed when the data changes.
"""
function build_residuals(fig, ax)
    # Not zoomable, for the χ² map's reason: the axis is the whole dataset, and there is
    # nothing outside it to pan to.
    for it in (:scrollzoom, :dragpan, :rectanglezoom, :limitreset)
        Makie.deregister_interaction!(ax, it)
    end

    pts = NamedTuple(k => Makie.Observable(Makie.Point2f[]) for k in RESIDUAL_KINDS)

    # The bands first, so points draw over them.
    Makie.hlines!(ax, [-3.0, 3.0]; color = (:black, 0.18), linewidth = 1, linestyle = :dot)
    Makie.hlines!(ax, [-1.0, 1.0]; color = (:black, 0.28), linewidth = 1, linestyle = :dash)
    Makie.hlines!(ax, [0.0];       color = (:black, 0.55), linewidth = 1.2)

    scatters = NamedTuple(k => Makie.scatter!(ax, pts[k];
                                              color = RESIDUAL_COLORS[k],
                                              markersize = 5, marker = OIPLOT_MARKER)
                          for k in RESIDUAL_KINDS)

    ax.xlabel = "Baseline (Mλ)"
    ax.ylabel = "(model − data) / σ"
    ax.title  = ""

    # The legend is built here with every entry, for the same reason the scatters are: a
    # `Legend` created after the window is up allocates. Entries for absent observables are
    # hidden by emptying their series, which leaves the label in place -- so the panel says
    # "t3amp 0" rather than silently dropping it, and a missing observable is visible.
    leg = Makie.Legend(fig[1, 2],
                       [scatters[k] for k in RESIDUAL_KINDS],
                       [string(k) for k in RESIDUAL_KINDS];
                       framevisible = false, labelsize = OIPLOT_LEGEND_SIZE + 1,
                       patchsize = (10.0f0, 10.0f0), rowgap = 1, padding = (2, 2, 2, 2))
    Makie.colgap!(fig.layout, 1, 4)

    return (; figure = fig, axis = ax, points = pts, scatters, legend = leg)
end

"""
    residual_series(data, res) -> Vector{Tuple{Symbol,Vector{Point2f},Float64}}

Pair each residual vector with the baselines it belongs to, and its rms.

`model_to_residuals` returns one vector per observable in the data's own point order, so the
x values are that observable's own baseline field -- `t3_baseline` for the closure quantities,
`vis_baseline` for the complex ones. Taking them from `OBS_SPECS` rather than naming the
fields here is what keeps this agreeing with the observable plots, where the same choice of
geometric-mean baseline was easy to get wrong and was.
"""
function residual_series(data, res)
    out = Tuple{Symbol,Vector{Makie.Point2f},Float64}[]
    for k in RESIDUAL_KINDS
        y = getproperty(res, k)
        spec = OBS_SPECS[RESIDUAL_SPECS[k]]
        x = getfield(data, spec.x)
        # A residual vector and its baselines can disagree in length only if the data changed
        # under the model; draw nothing rather than a scatter of mismatched pairs.
        if isempty(y) || length(x) != length(y)
            push!(out, (k, Makie.Point2f[], NaN))
            continue
        end
        good = findall(i -> isfinite(y[i]) && isfinite(x[i]), eachindex(y))
        pts  = [Makie.Point2f(x[i] * spec.xscale, y[i]) for i in good]
        rms  = isempty(good) ? NaN : sqrt(sum(abs2, @view y[good]) / length(good))
        push!(out, (k, pts, rms))
    end
    return out
end

"""
    update_residuals!(panel, data, res) -> String

Draw one model's residuals and return the one-line summary the console shows.

The summary reports rms per observable, not χ². rms is the same number the ±1 and ±3 lines
are drawn at, so the text and the picture say the same thing; a reduced χ² is reported by the
fits table and would mean a second convention to reconcile.
"""
function update_residuals!(panel, data, res)
    panel === nothing && return ""
    series = residual_series(data, res)

    lo, hi = Inf, -Inf
    for (k, pts, _) in series
        panel.points[k][] = pts
        for p in pts
            lo = min(lo, p[2]); hi = max(hi, p[2])
        end
    end

    # x from the data, y symmetric about zero -- in that order, because `autolimits!` resets
    # both. The axis is a signed residual, and a y range that is not symmetric makes a centred
    # scatter look biased.
    Makie.autolimits!(panel.axis)
    m = (isfinite(lo) && isfinite(hi)) ? max(abs(lo), abs(hi), 3.5) : 3.5   # never inside ±3
    Makie.ylims!(panel.axis, -1.05m, 1.05m)

    parts = String[]
    for (k, pts, rms) in series
        isempty(pts) && continue
        push!(parts, string(k, " ", length(pts), " pts rms ", round(rms; digits = 2)))
    end
    return isempty(parts) ? "no residuals to draw" : join(parts, ", ")
end
