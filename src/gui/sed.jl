# The model's spectral energy distribution, in Makie.
#
# `model_to_sed` evaluates V(0,0) -- the zero-baseline visibility, which IS the total flux --
# at each wavelength, plus each component's own `f` as the resolver evaluates it there. So the
# curves answer one question: how the model divides its light between components across the
# band, and whether the total is the constant the flux convention expects.
#
# It is drawn only for a model that references `$WL`. A model that does not is the same model
# at every wavelength, and a panel of flat lines is not a spectrum -- it is a control that
# looks broken. `shell_model_sed` says so rather than drawing them.

"How many component curves the panel can draw. Beyond this the total is still right and the
per-component curves are dropped, which the summary says."
const SED_MAX_COMPONENTS = 8

"""
    build_sed(figure, axis) -> NamedTuple

Create every plot the SED panel will ever need. **Call before the window is shown**, for the
reason `build_canvas` documents.

The component count is not known until a model exists, so a fixed set of lines is built and
the unused ones left empty. A line created when a component is added would allocate GL buffers
after the window is up, which is the one thing that cannot happen.
"""
function build_sed(fig, ax)
    for it in (:scrollzoom, :dragpan, :rectanglezoom, :limitreset)
        Makie.deregister_interaction!(ax, it)
    end

    comps = [Makie.Observable(Makie.Point2f[]) for _ in 1:SED_MAX_COMPONENTS]
    total = Makie.Observable(Makie.Point2f[])

    # Colour is an Observable so an unused slot can be made transparent: the Legend draws a
    # swatch from the plot whether or not the plot holds any points, and a column of coloured
    # stubs against blank labels reads as components that failed to draw.
    colors = [Makie.Observable(Makie.RGBAf(0, 0, 0, 0)) for _ in 1:SED_MAX_COMPONENTS]
    lines  = [Makie.lines!(ax, comps[i]; color = colors[i], linewidth = 1.6)
              for i in 1:SED_MAX_COMPONENTS]
    # The total last, so it draws over the components it is the sum of.
    tl = Makie.lines!(ax, total; color = :black, linewidth = 2.4)

    ax.xlabel = "Wavelength (μm)"
    ax.ylabel = "Flux"
    ax.title  = ""

    # Built here with every entry for the reason the lines are, and relabelled per model:
    # `Legend` entries take their text from an Observable, but the Legend itself cannot be
    # created once the window is up.
    labels = [Makie.Observable(" ") for _ in 1:SED_MAX_COMPONENTS]
    leg = Makie.Legend(fig[1, 2], [tl; lines], ["total"; labels];
                       framevisible = false, labelsize = OIPLOT_LEGEND_SIZE + 1,
                       patchsize = (14.0f0, 10.0f0), rowgap = 1, padding = (2, 2, 2, 2))
    Makie.colgap!(fig.layout, 1, 4)

    return (; figure = fig, axis = ax, comps, total, lines, totalline = tl, labels, colors,
              legend = leg)
end

"""
    update_sed!(panel, wl, total, comp_fluxes) -> String

Draw one model's SED and return the one-line summary.

`wl` is in metres, as `model_to_sed` takes it and as the data carries it; the axis is in
microns, so the conversion is here rather than in the caller. Components are drawn in name
order so the colours do not reshuffle when the model is edited.
"""
function update_sed!(panel, wl, total, comp_fluxes)
    panel === nothing && return ""
    names = sort(collect(keys(comp_fluxes)))
    ndrawn = min(length(names), SED_MAX_COMPONENTS)

    um = 1e6 .* wl
    panel.total[] = [Makie.Point2f(um[i], total[i]) for i in eachindex(um)]
    for i in 1:SED_MAX_COMPONENTS
        if i <= ndrawn
            f = comp_fluxes[names[i]]
            panel.comps[i][]  = [Makie.Point2f(um[j], f[j]) for j in eachindex(um)]
            panel.labels[i][] = names[i]
            panel.colors[i][] = Makie.to_color(OIPLOT_COLORS[mod1(i, length(OIPLOT_COLORS))])
        else
            panel.comps[i][]  = Makie.Point2f[]
            panel.labels[i][] = " "
            panel.colors[i][] = Makie.RGBAf(0, 0, 0, 0)
        end
    end
    Makie.autolimits!(panel.axis)

    dropped = length(names) - ndrawn
    parts = ["$(length(wl)) wavelength" * (length(wl) == 1 ? "" : "s"),
             "$(ndrawn) component" * (ndrawn == 1 ? "" : "s")]
    dropped > 0 && push!(parts, "$(dropped) not drawn")
    push!(parts, "total " * string(round(minimum(total); digits = 4)) * "–" *
                 string(round(maximum(total); digits = 4)))
    return join(parts, ", ")
end
