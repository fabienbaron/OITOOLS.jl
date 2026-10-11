# Corner plots, via PairPlots.jl.
#
# `plot_corner_makie` is declared in the core package and left unimplemented there because drawing
# a good corner plot is a package's worth of work — kernel density estimates, credible-interval
# contours, marginal histograms and a triangular layout that stays readable at ten parameters.
# PairPlots.jl does all of it and is Makie-native, so it composes with OITOOLSMakieExt rather
# than duplicating it.
#
# This is the Python-free counterpart of `plot_ultranest_corner`, which wraps
# `ultranest.plot.cornerplot` and can therefore only ever serve that one sampler. This takes a
# plain sample matrix, so it serves both nested samplers, the bootstrap and anything else that
# produces a posterior.
module OITOOLSPairPlotsExt

using OITOOLS
using PairPlots
using Makie

import OITOOLS: plot_corner_makie

using OITOOLS: PlotData, PLOT_FONTS, BootstrapResult

"""
    plot_corner_makie(posterior, names; kwargs...) -> PlotData

Corner plot of an `(n_samples, n_params)` sample matrix, one label per column.

```julia
using OITOOLS, Nautilus, PairPlots, CairoMakie
r = fit_model_nested(model_dict, free, data; lb, ub)
Makie.save("corner.png", plot_corner_makie(r.posterior, r.list_free_params).figure)
```

`truth`, one value per column, marks those values (the best fit, say) on every panel.
`kwargs` are passed to `PairPlots.pairplot`. A parameter whose samples are all identical is
dropped with a warning rather than passed on: a zero-width column has no density to estimate,
and PairPlots fails inside the KDE rather than at the call.
"""
function plot_corner_makie(posterior::AbstractMatrix{<:Real},
                       names::AbstractVector{<:AbstractString}; truth = nothing, kwargs...)
    size(posterior, 2) == length(names) || throw(ArgumentError(
        "got $(length(names)) names for $(size(posterior, 2)) columns"))
    size(posterior, 1) > 1 || throw(ArgumentError(
        "a corner plot needs more than one sample; got $(size(posterior, 1))"))

    keep = [j for j in axes(posterior, 2)
            if !isapprox(minimum(@view posterior[:, j]), maximum(@view posterior[:, j]))]
    dropped = setdiff(axes(posterior, 2), keep)
    isempty(dropped) ||
        @warn "plot_corner_makie: dropping parameters with no spread" params = names[dropped]
    isempty(keep) && throw(ArgumentError("every parameter has zero spread; nothing to draw"))

    # PairPlots takes a Tables.jl table; a NamedTuple of columns is the least ceremony, and
    # the Symbol keys are what label the axes.
    tbl = NamedTuple{Tuple(Symbol.(names[keep]))}(
              Tuple(Vector{Float64}(@view posterior[:, j]) for j in keep))

    fig = Makie.Figure(fonts = PLOT_FONTS)
    if truth === nothing
        PairPlots.pairplot(fig[1, 1], tbl; kwargs...)
    else
        length(truth) == length(names) || throw(ArgumentError(
            "got $(length(truth)) truth values for $(length(names)) columns"))
        t = NamedTuple{keys(tbl)}(Tuple(Float64(truth[j]) for j in keep))
        PairPlots.pairplot(fig[1, 1], tbl, PairPlots.Truth(t; label = "best fit"); kwargs...)
    end
    return PlotData(fig, nothing, nothing, String[])
end

# A bootstrap: the replicates its clipping kept, which are what its σ and correlations are of.
plot_corner_makie(r::BootstrapResult; kwargs...) =
    plot_corner_makie(r.samples[r.mask, :], r.list_free_params; kwargs...)

# A `NestedResult` — anything carrying `posterior` and the names beside it.
plot_corner_makie(r; kwargs...) =
    plot_corner_makie(r.posterior, r.list_free_params; kwargs...)

# Building the figure is most of a first corner plot: measured, 10.9 s of a 17.6 s first draw
# (the other 6.3 s is GLMakie rendering it, which needs a GL context and so belongs to the
# sysimage trace, not here). Three columns with a marked best fit, as the Modeling tab draws
# them, so the 2-D contour panels and the truth lines are compiled too.
using PrecompileTools

@compile_workload begin
    X = [sin(0.7i + j) + 0.1j * cos(1.3i) for i in 1:200, j in 1:3]
    plot_corner_makie(X, ["a", "b", "c"]; truth = [0.0, 0.0, 0.0])
end

end # module
