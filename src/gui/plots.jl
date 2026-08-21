# Interactive plots, as Makie figures.
#
# These return a `Figure` and the per-point metadata a click should report. Makie supplies
# everything the design would otherwise hand-roll — axes, ticks, legends, colorbars, log scales,
# scroll-zoom and drag-pan — so this file is only about turning an `OIdata` into the right
# marks. Hit-testing is not among them: `Makie.pick` needs a pick buffer that QMLMakie never
# renders, so `shell_pick` searches the points directly.
#
# Backend is chosen by the caller, and that is what keeps them testable: the GUI activates
# GLMakie (through QMLMakie, whose `MakieArea` forwards mouse events into the scene), while
# the tests activate CairoMakie and render to PNG with no GPU or display at all.
#
# Static publication figures still go through oiplot.jl's matplotlib functions — an exported
# script must reproduce the exact figure a paper used, and that means matplotlib, not Makie.

"""
    PlotData(figure, axis, points, info)

A built figure plus what a click needs. `points` is the Makie plot object `Makie.pick` will
return; `info[i]` describes point `i` of it.
"""
struct PlotData
    figure :: Any
    axis   :: Any
    points :: Any
    info   :: Vector{String}
end

# oiplot.jl's palette, copied exactly (src/oiplot.jl `oiplot_colors`). NOT tab20 — an
# earlier version of this file assumed tab20 and every baseline came out the wrong colour,
# so the same data read differently on screen and in a paper figure.
const OIPLOT_COLORS = [
    "black", "gold", "chartreuse", "blue", "red", "pink", "lightgray", "darkorange",
    "darkgreen", "aqua", "fuchsia", "saddlebrown", "dimgray", "darkslateblue", "violet",
    "indigo", "navy", "dodgerblue", "sienna", "olive", "purple", "darkorchid", "tomato",
    "darkturquoise", "steelblue", "seagreen", "darkgoldenrod", "darkseagreen",
]

# oiplot's error-bar colour and default marker, so the encoding matches end to end.
const OIPLOT_ECOLOR = "gainsboro"
const OIPLOT_MARKER = :circle

# Style taken from oiplot.jl's globals and the rcParams `set_oiplot_defaults` installs, so a
# figure looks the same on screen as in the paper. Measured from a live matplotlib axes
# rather than assumed: label and tick labels 12 pt, legend 8 pt in 4 columns.
const OIPLOT_LABELSIZE     = 12.0
const OIPLOT_TICKLABELSIZE = 12.0
const OIPLOT_LEGEND_SIZE   = 8.0    # oiplot_legend_fontsize, used by _plot_obs
const OIPLOT_LEGEND_NCOL   = 4      # oiplot_legend_ncol
# uvplot does NOT use those globals: it hardcodes fontsize=10 with ncol=5 below the axes
# (ncol=3 upper-right when legend_below=false). Matching _plot_obs there would be wrong.
const UVPLOT_LEGEND_SIZE   = 10.0
const UVPLOT_LEGEND_NCOL   = 5
const OIPLOT_XTICKS        = 10     # matplotlib's default density; Makie's is far sparser

"""
Font for plot text.

Chosen for coverage, not for looks: it has to carry `μ` and `°` as well as `λ`. The default
face does not, so Makie fetches them from a fallback -- and under QMLMakie those come out as
overlapping strokes, while `λ`, which the default face does carry, is clean in the same
string. The identical strings render correctly through GLMakie's offscreen `save()`, so it is
the fallback path and not the glyphs.

Applied through `Figure(fonts = ...)` at construction, NOT by assigning to `*font` attributes.
Every per-attribute route was tried and all of them throw `Failed to resolve data_boundingbox`
out of Makie's compute graph -- `Axis.xlabelfont`, the ticklabel fonts, and
`Colorbar.ticklabelfont` alike -- whether set before or after the first plot.
"""
const PLOT_FONT = "DejaVu Sans"

"""
Every character the GUI's Makie labels can contain.

Makie's cached texture atlas pre-renders **only** `a-z`, `A-Z`, `0-9`, `.`, `-` and the Unicode
minus (`render_default_glyphs!`), and only for the fonts in its own default theme. Every other
character -- the space included -- is inserted into the atlas the first time it is drawn.

Under QMLMakie that runtime insertion is what corrupts: the glyph arrives with strokes of a
neighbouring atlas cell drawn through it. The same string rendered by plain GLMakie to a PNG is
clean, so it is the shared-context upload rather than the font or the text itself. And because
this GUI pins `PLOT_FONT` rather than Makie's default, NONE of the pre-rendered glyphs apply --
every character is a runtime insertion, which is why the corruption is widespread rather than
occasional, and why it varies by driver.

`prewarm_glyphs!` renders all of these into the atlas before the window exists, so the initial
texture upload carries the complete set and nothing has to be added afterwards.
"""
const PLOT_GLYPHS = Char[
    # printable ASCII: letters, digits, and every punctuation mark a label or tick can carry
    (' ':'~')...,
    # the Greek and typographic characters the axis labels, titles and readouts use
    'λ', 'μ', 'α', 'δ', 'χ', 'Δ', 'θ', 'σ', 'π',
    '°', '²', '³', '±', '×', '·', '—', '–', '…', '≤', '≥', '≈', '∈', '→', '←',
    'ᵣ', '₀', '₁', '₂', '\u2212',                    # subscripts and the Unicode MINUS SIGN
]

"""
    prewarm_glyphs!(fonts = PLOT_FONTS)

Insert every glyph in [`PLOT_GLYPHS`](@ref) into Makie's texture atlas. **Call before the
window exists.**

Once QML is up the GL context belongs to Qt's render thread, and a glyph met after that point
is uploaded into a texture that is not ours to touch. Filling the atlas first is what makes the
labels render as themselves — see `PLOT_GLYPHS` for the measurements behind this.

Missing glyphs are skipped rather than thrown: a font without `≥` should cost that one
character, not the window.
"""
function prewarm_glyphs!(fonts = PLOT_FONTS)
    atlas = Makie.get_texture_atlas()
    n = 0
    for name in unique(values(fonts))
        font = try
            Makie.to_font(name)
        catch
            continue
        end
        for c in PLOT_GLYPHS
            try
                Makie.insert_glyph!(atlas, c, font)
                n += 1
            catch
                # A glyph this face does not have. Makie will fall back when it is drawn.
            end
        end
    end
    return n
end

"""
    style_axis!(ax; scale = 1.0)

Apply oiplot's typography to a Makie axis.

`scale` multiplies every size. It defaults to 1.0, which is oiplot's measured 12 pt and what
the port harness compares against; the live GUI passes something larger, because a figure
sized for a paper column is not sized for a window across the room.
"""
function style_axis!(ax; scale::Real = 1.0)
    ax.xlabelsize[] = OIPLOT_LABELSIZE * scale
    ax.ylabelsize[] = OIPLOT_LABELSIZE * scale
    ax.titlesize[]  = OIPLOT_LABELSIZE * scale
    ax.xticklabelsize[] = OIPLOT_TICKLABELSIZE * scale
    ax.yticklabelsize[] = OIPLOT_TICKLABELSIZE * scale
    # Makie defaults to ~3 ticks where matplotlib gives ~10; without this the axes read very
    # differently even though the data is identical.
    ax.xticks[] = Makie.LinearTicks(OIPLOT_XTICKS)
    ax.yticks[] = Makie.LinearTicks(OIPLOT_XTICKS)
    # Restore Makie's own 5% headroom, which something else on this axis has taken away.
    #
    # `build_canvas` pre-creates a hidden heatmap for the Image perspective, and image-like
    # plots declare a zero autolimit margin -- an image wants no padding around it. The axis
    # adopts that for EVERY plot on it, so the scatter views ended up with limits sitting
    # exactly on the data and the outermost points bisected by the frame. Worse, the limit is
    # stored as Float32 and can round just inside a Float64 datum, which put two uv points
    # outside the axis entirely. Set on every redraw, because that is when it matters.
    ax.xautolimitmargin[] = (0.05f0, 0.05f0)
    ax.yautolimitmargin[] = (0.05f0, 0.05f0)
    return ax
end

"""
    add_baseline_legend!(figure, axis, names; col = 2)

Legend with one entry per baseline group, in oiplot's order, colours, font size and column
count. `uvplot` and the observable plots both draw one; a port without it looks materially
different from the published figure.
"""
function add_baseline_legend!(fig, ax, names;
                             size::Real = OIPLOT_LEGEND_SIZE,
                             ncol::Integer = OIPLOT_LEGEND_NCOL,
                             below::Bool = false)
    cmap   = baseline_color_map(names)
    groups = sort(unique(names))
    elems  = [Makie.MarkerElement(color = cmap[g], marker = OIPLOT_MARKER, markersize = 8)
              for g in groups]
    slot = below ? fig[2, 1] : fig[1, 2]
    block = Makie.Legend(slot, elems, groups;
                         labelsize = size, nbanks = ncol, framevisible = true,
                         orientation = below ? :horizontal : :vertical)
    # `block` is returned so a redraw can remove the previous legend; without it every
    # replot stacks another one into the layout.
    return (haslegend = true, nlegend = length(groups),
            legendsize = Float64(size), legendncol = Int(ncol), block = block)
end

"""
    baseline_color_map(names) -> Dict{String,String}

Group-to-colour mapping, following `_plot_obs` exactly: groups are
`sort(unique(names))` — **sorted**, not first-appearance — and the colour is
`oiplot_colors[mod1(i, 28)]` for the group's position `i`.

Both halves matter. Getting the palette right but the ordering wrong still recolours every
baseline.
"""
function baseline_color_map(names)
    groups = sort(unique(names))
    Dict(g => OIPLOT_COLORS[mod1(i, length(OIPLOT_COLORS))] for (i, g) in enumerate(groups))
end

_station_lut(d) = Dict(si => (k <= length(d.sta_name) ? d.sta_name[k] : string(si))
                       for (k, si) in enumerate(d.sta_index))

"""
    ObsSpec

Mirror of oiplot.jl's `ObsPlotSpec`: which arrays feed a plot, how points are grouped for
colouring, and the axis labels. Driving from a table rather than an `if` chain is deliberate —
it is how oiplot does it, so the two stay aligned when either gains an observable.

`xscale` converts to display units (1e-6 for Mλ, 1e6 for μm), exactly as oiplot's does.
"""
struct ObsSpec
    y :: Symbol; yerr :: Symbol; x :: Symbol; lam :: Symbol; mjd :: Symbol; sta :: Symbol
    grouping :: Symbol      # :baseline, :triplet or :station
    ylabel :: String; xlabel :: String; xscale :: Float64
    default_color :: Symbol # each oiplot function has its own default; plot_flux uses "wav"
end

# Values copied from OBS_PLOT_SPECS. Two things here are easy to get wrong and both were:
# closure quantities default to the GEOMETRIC-MEAN baseline (`t3_baseline`), not the longest
# leg, and flux is plotted against WAVELENGTH, not a baseline.
const OBS_SPECS = Dict{Symbol,ObsSpec}(
    :v2        => ObsSpec(:v2, :v2_err, :v2_baseline, :v2_lam, :v2_mjd,
                          :v2_sta_index, :baseline, "V²", "Baseline (Mλ)", 1e-6, :baseline),
    :t3phi     => ObsSpec(:t3phi, :t3phi_err, :t3_baseline, :t3_lam, :t3_mjd,
                          :t3_sta_index, :triplet, "T3φ (°)", "Baseline (Mλ)", 1e-6, :baseline),
    :t3phi_max => ObsSpec(:t3phi, :t3phi_err, :t3_maxbaseline, :t3_lam, :t3_mjd,
                          :t3_sta_index, :triplet, "T3φ (°)", "Max. baseline (Mλ)", 1e-6, :baseline),
    :t3amp     => ObsSpec(:t3amp, :t3amp_err, :t3_baseline, :t3_lam, :t3_mjd,
                          :t3_sta_index, :triplet, "T3amp", "Baseline (Mλ)", 1e-6, :baseline),
    :t3amp_max => ObsSpec(:t3amp, :t3amp_err, :t3_maxbaseline, :t3_lam, :t3_mjd,
                          :t3_sta_index, :triplet, "T3amp", "Max. baseline (Mλ)", 1e-6, :baseline),
    :visamp    => ObsSpec(:visamp, :visamp_err, :vis_baseline, :vis_lam, :vis_mjd,
                          :vis_sta_index, :baseline, "Visamp", "Baseline (Mλ)", 1e-6, :baseline),
    :visphi    => ObsSpec(:visphi, :visphi_err, :vis_baseline, :vis_lam, :vis_mjd,
                          :vis_sta_index, :baseline, "Visφ (°)", "Baseline (Mλ)", 1e-6, :baseline),
    :flux      => ObsSpec(:flux, :flux_err, :flux_lam, :flux_lam, :flux_mjd,
                          :flux_sta_index, :station, "Flux", "Wavelength (μm)", 1e6, :wav),
    # Against WAVELENGTH, not baseline. There is no separate differential-phase field: OIFITS
    # carries it in OI_VIS with `PHITYP=differential`, so this is `visphi` plotted the way a
    # differential quantity is read — one curve per baseline across the band. Against baseline
    # it is a scatter of every channel at once, which is the same numbers and no signal.
    :diffphi   => ObsSpec(:visphi, :visphi_err, :vis_lam, :vis_lam, :vis_mjd,
                          :vis_sta_index, :baseline, "Diff. phase (°)",
                          "Wavelength (μm)", 1e6, :baseline),
    :diffvisamp => ObsSpec(:visamp, :visamp_err, :vis_lam, :vis_lam, :vis_mjd,
                          :vis_sta_index, :baseline, "Diff. visamp",
                          "Wavelength (μm)", 1e6, :baseline),
)

"""
Observables a log y axis is meaningful for.

Positive quantities only: a phase crosses zero and takes both signs, so a log axis on one is
not a display choice but a way to silently drop half the data.
"""
const LOG_Y_KINDS = Set([:v2, :visamp, :t3amp, :t3amp_max, :flux, :diffvisamp])

"""Baseline label per point, e.g. `"S1-E1"` (oiplot's `get_baseline_names`)."""
function baseline_names(d, sta::Symbol = :v2_sta_index)
    lut = _station_lut(d)
    idx = getfield(d, sta)
    [string(get(lut, idx[1, i], "?"), "-", get(lut, idx[2, i], "?")) for i in 1:size(idx, 2)]
end

"""Triplet label per closure point, e.g. `"S1-E1-W2"` (oiplot's `get_triplet_names`)."""
function triplet_names(d, sta::Symbol = :t3_sta_index)
    lut = _station_lut(d)
    idx = getfield(d, sta)
    [string(get(lut, idx[1, i], "?"), "-", get(lut, idx[2, i], "?"), "-",
            get(lut, idx[3, i], "?")) for i in 1:size(idx, 2)]
end

"""
Station label per flux point (oiplot's `get_station_names`).

Index 0 means a calibrated source spectrum with no telescope attached, which oiplot labels
"Calibrated"; it must not be looked up in `sta_name`.
"""
station_names(d, sta::Symbol = :flux_sta_index) =
    [si == 0 ? "Calibrated" : get(_station_lut(d), si, string(si)) for si in getfield(d, sta)]

"""Group labels for a spec, following its `grouping` rule."""
group_names(d, spec::ObsSpec) =
    spec.grouping === :triplet ? triplet_names(d, spec.sta) :
    spec.grouping === :station ? station_names(d, spec.sta) :
                                 baseline_names(d, spec.sta)

"""
    grouping_noun(kind) -> String

What one panel of the per-group view contains, for the control that switches it on.

Not always a baseline: closure quantities group by TRIPLET and flux by STATION, and a tick box
reading "per baseline" over a grid of closure triangles is simply wrong. `ObsSpec.grouping`
already knows, so the label follows it.
"""
function grouping_noun(kind::Symbol)
    spec = get(OBS_SPECS, kind, nothing)
    spec === nothing && return "group"
    return spec.grouping === :triplet ? "triplet" :
           spec.grouping === :station ? "station" : "baseline"
end

"""
    panel_data(d, kind; color = nothing) -> (; panels, xlabel, ylabel, grouping)

One entry of `panels` per group, each `(; name, x, y, err, color)` with `x` in μm.

Against WAVELENGTH always, whatever the observable's usual x axis is. That is the point of the
view: a baseline's spectrum is what shows a line, or a slope that says the calibration drifted,
and neither is visible when every baseline is piled onto one axis against its length.

Groups are ordered by name so the panels keep their positions between redraws — a panel that
moves when the colour mode changes is one the eye has to find again.
"""
function panel_data(d, kind::Symbol; color = nothing)
    spec = get(OBS_SPECS, kind, nothing)
    spec === nothing && throw(ArgumentError(
        "no per-group view for $(repr(kind)); it is not an observable"))

    names = group_names(d, spec)
    lam   = Float64.(getfield(d, spec.lam)) .* 1e6      # μm
    y     = Float64.(getfield(d, spec.y))
    e     = Float64.(getfield(d, spec.yerr))
    isempty(y) && throw(ArgumentError("no $(kind) data in $(basename(d.filename))"))

    cmap = baseline_color_map(names)
    panels = NamedTuple[]
    for g in sort(unique(names))
        idx = findall(==(g), names)
        # Sorted by wavelength, so a line through the points is the spectrum rather than the
        # order the rows happen to sit in the table.
        perm = sortperm(lam[idx])
        k = idx[perm]
        push!(panels, (; name = g, x = lam[k], y = y[k], err = e[k],
                         color = Makie.RGBAf(Makie.to_color(cmap[g]))))
    end
    return (; panels, xlabel = "λ (μm)", ylabel = spec.ylabel,
              grouping = grouping_noun(kind))
end

"""Generic click description for any observable."""
obs_info(d, spec::ObsSpec, i::Integer, names) =
    string(names[i],
           "   x = ", round(Float64(getfield(d, spec.x)[i]) * spec.xscale, digits = 3),
           "   λ = ", round(Float64(getfield(d, spec.lam)[i]) * 1e6, digits = 3), " μm",
           "   ", strip(spec.ylabel), " = ",
           round(Float64(getfield(d, spec.y)[i]), digits = 4),
           " ± ", round(Float64(getfield(d, spec.yerr)[i]), digits = 4))

"""
    point_info(data, i) -> String

What a click reports for V² point `i`: station pair, spatial frequency in Mλ, wavelength,
MJD, value and error, and whether the point is flagged.
"""
function point_info(d, i::Integer, names = baseline_names(d))
    string(names[i],
           "   B/λ = ", round(Float64(d.v2_baseline[i]) / 1e6, digits = 2), " Mλ",
           "   λ = ", round(Float64(d.v2_lam[i]) * 1e6, digits = 3), " μm",
           "   MJD = ", round(d.v2_mjd[i], digits = 4),
           "   V² = ", round(Float64(d.v2[i]), digits = 4),
           " ± ", round(Float64(d.v2_err[i]), digits = 4),
           d.v2_flag[i] ? "   [FLAGGED]" : "")
end

_colors_for(d, mode::Symbol, names, spec::ObsSpec) =
    if mode === :baseline || mode === :station || mode === :triplet
        cmap = baseline_color_map(names)
        [cmap[n] for n in names]
    elseif mode === :wav
        Float64.(getfield(d, spec.lam)) .* 1e6
    elseif mode === :mjd
        Float64.(getfield(d, spec.mjd))
    else
        "#1f77b4"
    end

"""
    uv_point_labels(data) -> (names, info)

Per **uv point** (not per V² point) station-pair labels and click descriptions.

`OIdata.uv` holds every sampled spatial frequency: the V² baselines and the three legs of
each closure triangle. `indx_v2` and `indx_t3_1/2/3` map observables onto those columns, and
after filtering or uv de-duplication a column can be shared, so the mapping has to be read
rather than assumed.

Points reached only through a closure triangle have no V² value; they are labelled as legs.
"""
function uv_point_labels(d)
    lut   = _station_lut(d)
    names = fill("", d.nuv)
    info  = fill("", d.nuv)
    bn    = baseline_names(d)

    for i in 1:d.nv2
        k = d.indx_v2[i]
        (1 <= k <= d.nuv) || continue
        names[k] = bn[i]
        info[k]  = point_info(d, i, bn)
    end

    # OI_VIS next. Missing this branch left every VIS-only uv point unnamed: on v1295Aql that
    # was 2190 of 4932 points, all of them in `indx_vis` and none reachable from V2 or T3, so
    # they were drawn in one colour under a blank legend entry. A dataset without OI_VIS
    # (2004-data1) showed nothing wrong, which is why it survived.
    nvis = size(d.vis_sta_index, 2)
    for i in 1:min(nvis, length(d.indx_vis))
        k = d.indx_vis[i]
        (1 <= k <= d.nuv) || continue
        isempty(names[k]) || continue           # a V² label is at least as informative
        nm = string(get(lut, d.vis_sta_index[1, i], "?"), "-",
                    get(lut, d.vis_sta_index[2, i], "?"))
        names[k] = nm
        info[k]  = string(nm, "   OI_VIS",
                          "   B/λ = ", round(hypot(Float64(d.uv[1, k]),
                                                   Float64(d.uv[2, k])) / 1e6, digits = 2),
                          " Mλ   λ = ", round(Float64(d.uv_lam[k]) * 1e6, digits = 3), " μm")
    end

    nt3 = size(d.t3_sta_index, 2)
    legs = ((1, 2), (2, 3), (1, 3))
    for (li, idxv) in enumerate((d.indx_t3_1, d.indx_t3_2, d.indx_t3_3))
        a, b = legs[li]
        for i in 1:min(nt3, length(idxv))
            k = idxv[i]
            (1 <= k <= d.nuv) || continue
            isempty(names[k]) || continue          # a V² label is more informative
            nm = string(get(lut, d.t3_sta_index[a, i], "?"), "-",
                        get(lut, d.t3_sta_index[b, i], "?"))
            names[k] = nm
            info[k]  = string(nm, "   closure leg ", li,
                              "   B/λ = ", round(hypot(Float64(d.uv[1, k]),
                                                       Float64(d.uv[2, k])) / 1e6, digits = 2),
                              " Mλ   λ = ", round(Float64(d.uv_lam[k]) * 1e6, digits = 3), " μm")
        end
    end
    return names, info
end

"""
    draw!(figure, axis, data, kind; color, conjugate, logscale, markersize, extras) -> (plot, info)

**The single drawing implementation.** Everything else here is a wrapper.

`kind` is `:uv` or any key of `OBS_SPECS`. The axis is cleared and redrawn in place; `extras`
is a vector of legend/colorbar blocks belonging to previous draws, which are removed first.

One implementation, so what the harness tests is what ships. A second copy for the shell
would be free to lose the equal-aspect ratio or the legend without any test noticing.
"""
function draw!(fig, ax, d, kind::Symbol;
               color::Union{Nothing,Symbol} = nothing, conjugate::Bool = true,
               logscale::Bool = false, markersize::Union{Nothing,Real} = nothing,
               extras::Vector{Any} = Any[])
    for b in extras
        try; Makie.delete!(b); catch; end
    end
    empty!(extras)
    Makie.empty!(ax)

    if kind === :uv
        color = color === nothing ? :baseline : color
        names, info = uv_point_labels(d)
        u = Float64.(d.uv[1, :]) ./ 1e6
        v = Float64.(d.uv[2, :]) ./ 1e6
        x, y, inf = conjugate ? (vcat(u, -u), vcat(v, -v), vcat(info, info)) : (u, v, info)

        c = if color === :baseline
                cmap = baseline_color_map(names)
                cols = [cmap[n] for n in names]
                conjugate ? vcat(cols, cols) : cols
            elseif color === :wav
                w = Float64.(d.uv_lam) .* 1e6
                conjugate ? vcat(w, w) : w
            elseif color === :mjd
                m = Float64.(d.uv_mjd)
                conjugate ? vcat(m, m) : m
            else
                "#1f77b4"
            end

        ax.xlabel[] = "U (Mλ)"
        ax.ylabel[] = "V (Mλ)"
        ax.title[]  = "uv coverage — " * basename(d.filename)
        # Equal aspect is not cosmetic: without it the baselines are sheared and the coverage
        # is misread. Losing this in the shell's copy is what made the plot non-isotropic.
        ax.aspect[] = Makie.DataAspect()
        ax.yscale[] = identity
        style_axis!(ax)

        p = Makie.scatter!(ax, x, y; color = c, marker = OIPLOT_MARKER,
                           markersize = something(markersize, 6))
        if color === :wav || color === :mjd
            push!(extras, Makie.Colorbar(fig[1, 2], p;
                                         label = color === :wav ? "λ (μm)" : "MJD"))
        else
            push!(extras, add_baseline_legend!(fig, ax, names;
                                               size = UVPLOT_LEGEND_SIZE,
                                               ncol = UVPLOT_LEGEND_NCOL, below = true).block)
        end
        return p, inf
    end

    spec = get(OBS_SPECS, kind, nothing)
    spec === nothing && throw(ArgumentError(
        "unknown plot $(repr(kind)); use :uv or one of $(sort(collect(keys(OBS_SPECS))))"))
    color = color === nothing ? spec.default_color : color

    # plot_visphi has TWO layouts. With phityp="differential" oiplot draws a paginated grid,
    # one panel per baseline, against wavelength — not this plot. Refusing is honest:
    # silently drawing the absolute layout would look right and be a different figure.
    if kind === :visphi && occursin("differential", lowercase(String(d.phityp)))
        throw(ArgumentError(
            "this dataset has phityp=\"differential\": its phase is measured against the " *
            "other channels, so it belongs against WAVELENGTH — use the `diffphi` view."))
    end

    names = group_names(d, spec)
    x = Float64.(getfield(d, spec.x)) .* spec.xscale
    y = Float64.(getfield(d, spec.y))
    e = Float64.(getfield(d, spec.yerr))
    isempty(x) && throw(ArgumentError("no $(kind) data in $(basename(d.filename))"))
    info = [obs_info(d, spec, i, names) for i in eachindex(x)]
    c    = _colors_for(d, color, names, spec)

    ax.xlabel[] = spec.xlabel
    ax.ylabel[] = spec.ylabel
    ax.title[]  = string(kind) * " — " * basename(d.filename)
    ax.aspect[] = nothing                     # only uv coverage is isotropic
    style_axis!(ax)

    if logscale
        # Non-positive points are dropped, matching matplotlib under `set_yscale("log")` —
        # which is what oiplot's `logplot=true` produces, and what the port tests compare
        # against. The lower whisker of a surviving point can still reach below zero, and
        # Makie transforms the bar's bounding box, so that one is truncated at the axis floor.
        keep = [isfinite(v) && v > 0 for v in y]
        any(keep) || throw(ArgumentError(
            "cannot use a log scale: no positive $(kind) values in $(basename(d.filename))"))
        x, y, e, c, info = x[keep], y[keep], e[keep], c[keep], info[keep]
        floorv = minimum(y) / 10
        Makie.errorbars!(ax, x, y, y .- max.(y .- e, floorv), e; color = OIPLOT_ECOLOR)
    else
        Makie.errorbars!(ax, x, y, e; color = OIPLOT_ECOLOR)
    end
    p = Makie.scatter!(ax, x, y; color = c, marker = OIPLOT_MARKER,
                       markersize = something(markersize, 7))
    # Log scale is applied AFTER plotting, and only with limits that are valid for it.
    # V² contains zeros, and Makie validates limits against the scale the moment it is set:
    # `Invalid y-limits (0.0, 10.0) for scale log10`. matplotlib simply omits non-positive
    # points, so do the same rather than refuse a plot the user can legitimately ask for.
    if logscale
        lo, hi = extrema(y)
        Makie.ylims!(ax, lo / 2, hi * 2)
        ax.yscale[] = log10
        # LinearTicks on a log axis bunches every label into the top decade; see livecanvas.jl.
        ax.yticks[] = Makie.LogTicks(Makie.LinearTicks(OIPLOT_XTICKS))
    else
        ax.yscale[] = identity
    end

    if color === :wav || color === :mjd
        push!(extras, Makie.Colorbar(fig[1, 2], p;
                                     label = color === :wav ? "λ (μm)" : "MJD"))
    else
        push!(extras, add_baseline_legend!(fig, ax, names).block)
    end
    return p, info
end

"""
    uv_figure(data; color = :baseline, conjugate = true) -> PlotData

uv coverage in a fresh figure. Draws **every** uv point — V² baselines and closure-triangle
legs alike — which is what `uvplot` shows and what the coverage actually is.
"""
function uv_figure(d; color::Symbol = :baseline, conjugate::Bool = true, markersize::Real = 6)
    fig = Makie.Figure()
    ax  = Makie.Axis(fig[1, 1])
    p, info = draw!(fig, ax, d, :uv; color = color, conjugate = conjugate,
                    markersize = markersize)
    return PlotData(fig, ax, p, info)
end

"""
    observable_figure(data, which; color = nothing, logscale = false) -> PlotData

An observable against its natural x axis, in a fresh figure. `which` is any key of
`OBS_SPECS`. `color = nothing` uses that observable's own oiplot default.
"""
function observable_figure(d, which::Symbol = :v2; color::Union{Nothing,Symbol} = nothing,
                           logscale::Bool = false, markersize::Real = 7)
    fig = Makie.Figure()
    ax  = Makie.Axis(fig[1, 1])
    p, info = draw!(fig, ax, d, which; color = color, logscale = logscale,
                    markersize = markersize)
    return PlotData(fig, ax, p, info)
end

"""
    plot_into!(figure, axis, data, kind; kwargs...) -> (plot, info)

Redraw an existing axis — what the shell uses. A thin alias for [`draw!`](@ref) so the shell
and the figure builders cannot drift apart again.
"""
plot_into!(fig, ax, d, kind::Symbol; kwargs...) = draw!(fig, ax, d, kind; kwargs...)
