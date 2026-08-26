# The Makie drawing layer.
#
# Turns an `OIdata` into marks, using the descriptions in src/plot_data.jl. Makie supplies
# everything the design would otherwise hand-roll -- axes, ticks, legends, colorbars, log
# scales, scroll-zoom and drag-pan -- so this file is only about which marks to place.
#
# Hit-testing is not among them: `Makie.pick` needs a pick buffer that QMLMakie never renders,
# so `shell_pick` searches the points directly.
#
# The Makie BACKEND is the caller's choice, and that is what keeps this testable: the GUI
# activates GLMakie (through QMLMakie, whose `MakieArea` forwards mouse events into the
# scene), a script can activate CairoMakie and get a vector figure with no GPU or display at
# all, and the tests do exactly that.

"""
    prewarm_glyphs!(fonts = PLOT_FONTS) -> Int

Insert every glyph in [`PLOT_GLYPHS`](@ref) into Makie's texture atlas, and return how many
insertions were made. **Call before the first `Figure` exists**, and call it as a statement:

    prewarm_glyphs!()

NOT inside a `@debug`/`@info` string. Julia's logging macros do not evaluate their
interpolations when the level is disabled, so `@debug "\$(prewarm_glyphs!())"` silently never
runs -- which is exactly how the first version of this shipped doing nothing.

Missing glyphs are skipped rather than thrown: a font without `≥` should cost that one
character, not the window.
"""
function prewarm_glyphs!(fonts = PLOT_FONTS)
    # A warmed atlas is worth keeping. Rasterising is a signed-distance field per glyph --
    # measured at ~23 ms each, so 121 of them cost 2.4 s of every launch, before the first
    # plot. Storing the result and loading it back costs 0.02 s.
    #
    # The cache is keyed on the glyph set AND the Makie version, because it is Makie's own
    # binary format: a stale file after an upgrade would be silently wrong, and silently wrong
    # here means the corrupted glyphs this function exists to prevent. Anything unexpected
    # falls through to rasterising, which is always correct if slower.
    cache = _glyph_cache_path()
    if cache !== nothing && isfile(cache)
        try
            Makie.TEXTURE_ATLASES[(ATLAS_RESOLUTION, ATLAS_PIX_PER_GLYPH)] =
                Makie.load_texture_atlas(cache)
            return 0
        catch
            rm(cache; force = true)
        end
    end

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
                # A glyph this face does not have. Makie falls back when it is drawn.
            end
        end
    end
    if cache !== nothing
        try
            mkpath(dirname(cache))
            Makie.store_texture_atlas(cache, atlas)
        catch
            # A cache that cannot be written costs 2.4 s next launch and nothing else.
        end
    end
    return n
end

"The atlas `get_texture_atlas()` returns by default, and the key it is stored under."
const ATLAS_RESOLUTION    = 2048
const ATLAS_PIX_PER_GLYPH = 64

"""
Where the warmed glyph atlas is cached, or `nothing` if there is nowhere to put it.

Keyed on the glyph set and the Makie version: the file is Makie's own binary format, so one
written by a different version — or covering a different set of characters — must not be
loaded. A wrong atlas does not fail loudly; it draws wrong glyphs.
"""
function _glyph_cache_path()
    base = Sys.iswindows() ? get(ENV, "LOCALAPPDATA", get(ENV, "APPDATA", "")) :
           get(ENV, "XDG_CACHE_HOME", joinpath(homedir(), ".cache"))
    isempty(base) && return nothing
    key = string(hash((PLOT_GLYPHS, PLOT_FONT, string(pkgversion(Makie)),
                       ATLAS_RESOLUTION, ATLAS_PIX_PER_GLYPH)); base = 16)
    return joinpath(base, "oitools", "glyph_atlas_$(key).bin")
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
    uvplot_makie(data; color = :baseline, conjugate = true) -> PlotData

uv coverage in a fresh figure. Draws **every** uv point — V² baselines and closure-triangle
legs alike — which is what `uvplot` shows and what the coverage actually is.
"""
function uvplot_makie(d; color::Symbol = :baseline, conjugate::Bool = true, markersize::Real = 6)
    fig = Makie.Figure()
    ax  = Makie.Axis(fig[1, 1])
    p, info = draw!(fig, ax, d, :uv; color = color, conjugate = conjugate,
                    markersize = markersize)
    return PlotData(fig, ax, p, info)
end

"""
    plot_observable_makie(data, which; color = nothing, logscale = false) -> PlotData

An observable against its natural x axis, in a fresh figure. `which` is any key of
`OBS_SPECS`. `color = nothing` uses that observable's own oiplot default.
"""
function plot_observable_makie(d, which::Symbol = :v2; color::Union{Nothing,Symbol} = nothing,
                           logscale::Bool = false, markersize::Real = 7)
    fig = Makie.Figure()
    ax  = Makie.Axis(fig[1, 1])
    p, info = draw!(fig, ax, d, which; color = color, logscale = logscale,
                    markersize = markersize)
    return PlotData(fig, ax, p, info)
end

# ── One name per observable ──────────────────────────────────────────────────
#
# `plot_observable_makie` is the implementation and takes the observable as an argument,
# because the GUI dispatches on a runtime `kind` read from `OBS_SPECS` — a function per
# observable could not serve it. These wrappers exist so that a script does not have to know
# that: `plot_v2` and `plot_v2_makie` are the same call in the two backends, which is the whole
# point of the suffix.
for (fname, kind) in ((:plot_v2_makie,      :v2),
                      (:plot_t3phi_makie,   :t3phi),
                      (:plot_t3amp_makie,   :t3amp),
                      (:plot_visamp_makie,  :visamp),
                      (:plot_visphi_makie,  :visphi),
                      (:plot_flux_makie,    :flux),
                      (:plot_diffphi_makie, :diffphi))
    @eval begin
        """
            $($fname)(data; kwargs...) -> PlotData

        `$($(QuoteNode(kind)))` against its natural x axis, in Makie — the counterpart of
        `$(replace(string($(QuoteNode(fname))), "_makie" => ""))`.

        A thin wrapper over [`plot_observable_makie`](@ref), which takes any key of `OBS_SPECS`
        and is what to reach for when the observable is chosen at run time.
        """
        $fname(d; kwargs...) = plot_observable_makie(d, $(QuoteNode(kind)); kwargs...)
    end
end

"""
    plot_into!(figure, axis, data, kind; kwargs...) -> (plot, info)

Redraw an existing axis — what the shell uses. A thin alias for `draw!` so the shell
and the figure builders cannot drift apart again.
"""
plot_into!(fig, ax, d, kind::Symbol; kwargs...) = draw!(fig, ax, d, kind; kwargs...)

# ─────────────────────────────────────────────────────────────────────────────
# Images
# ─────────────────────────────────────────────────────────────────────────────
#
# The Makie counterparts of oiplot's `imdisp`/`imdisp_multi`. They are separate names rather
# than extra methods because the two backends hand back different things — a matplotlib figure
# and a Makie one — and a single name over both would have to lie about which.
#
# ORIENTATION IS EAST-LEFT, NORTH-UP, which is Monnier's convention and what `imdisp` draws.
# `imdisp` gets there with `rotl90` and a reversed extent; here the x coordinate vector simply
# descends, which is the same picture and one fewer transform to get wrong. The GUI's
# `show_image!` does it this way too, so a reconstruction looks identical in a window and in a
# saved figure.

"""
    imdisp_makie(image; pixsize = -1, kwargs...) -> PlotData

One image, East left and North up, normalised to its own maximum.

`image` may be a matrix or a flat vector of `nx²` pixels. `pixsize` is in mas; `-1` labels the
axes in pixels instead.

| keyword | default | meaning |
|---|---|---|
| `colormap` | `:gist_heat` | any Makie colormap; `imdisp`'s default is the same ramp |
| `title` | `""` | axis title |
| `colorbar` | `false` | draw a colour bar beside the image |
| `beamsize` | `-1` | if positive, draw a beam circle of this diameter (mas) |
| `beamlocation` | `(0.8, 0.8)` | fractional position of that circle |

Normalising to the maximum is `imdisp`'s behaviour and is kept deliberately: it makes two
reconstructions of the same target comparable by eye. An image whose maximum is below `1e-20`
is drawn unnormalised, with a warning, rather than divided by roughly zero.
"""
function imdisp_makie(image::AbstractArray; pixsize::Real = -1.0,
                      colormap = :gist_heat, title::AbstractString = "",
                      colorbar::Bool = false, beamsize::Real = -1.0,
                      beamlocation = (0.8, 0.8), scale::Real = 1.0)
    fig = Makie.Figure(fonts = PLOT_FONTS)
    ax  = Makie.Axis(fig[1, 1]; aspect = 1)
    hm  = image_into!(fig, ax, image; pixsize, colormap, title, colorbar,
                      beamsize, beamlocation, scale)
    return PlotData(fig, ax, hm, String[])
end

"""
    image_into!(figure, axis, image; kwargs...) -> heatmap

Draw one image into an existing axis. [`imdisp_makie`](@ref) is this plus a figure, and
[`imdisp_multi_makie`](@ref) calls it once per panel — one implementation, so a cube panel and a
single image cannot drift apart.

`normalise` is the divisor: `nothing` uses the image's own maximum.
"""
function image_into!(fig, ax, image::AbstractArray; pixsize::Real = -1.0,
                     colormap = :gist_heat, title::AbstractString = "",
                     colorbar::Bool = false, beamsize::Real = -1.0,
                     beamlocation = (0.8, 0.8), scale::Real = 1.0,
                     normalise::Union{Nothing,Real} = nothing)
    img = if ndims(image) == 1
        n = isqrt(length(image))
        n * n == length(image) ||
            error("a flat image must be square; got $(length(image)) pixels")
        reshape(image, n, n)
    else
        image
    end
    nx, ny = size(img)

    pixmode = pixsize <= 0
    ps = pixmode ? 1.0 : Float64(pixsize)

    # `nothing` means this panel's own maximum, which is `imdisp`'s behaviour. A number is
    # what `imdisp_multi_makie` passes to put every panel on one scale.
    norm = normalise === nothing ? maximum(img) : normalise
    if abs(norm) < 1e-20
        @warn "Maximum of image < tol; drawing it unnormalised" maxlog = 1
        norm = one(eltype(img))
    end

    halfx = nx * ps / 2
    halfy = ny * ps / 2
    xs = range(halfx, -halfx; length = nx)      # descending: East is to the LEFT
    ys = range(-halfy, halfy; length = ny)

    hm = Makie.heatmap!(ax, xs, ys, Float32.(img ./ norm);
                        colormap = colormap, interpolate = false)
    style_axis!(ax; scale)
    ax.title  = title
    ax.xlabel = pixmode ? "x ← E (pixels)" : "x ← E (mas)"
    ax.ylabel = pixmode ? "y → N (pixels)" : "y → N (mas)"
    Makie.limits!(ax, halfx, -halfx, -halfy, halfy)

    if beamsize > 0
        # Drawn in data coordinates, at a fraction of the half-width, exactly as `imdisp`
        # places it. White and unstroked, so it reads as a beam and not as a source.
        cx = halfx * beamlocation[1]
        cy = -halfy * beamlocation[2]
        Makie.poly!(ax, Makie.Circle(Makie.Point2f(cx, cy), Float32(beamsize / 2));
                    color = :white, strokewidth = 0)
    end
    colorbar && Makie.Colorbar(fig[1, 2], hm; label = "flux / pixel")
    return hm
end

"""
    imdisp_multi_makie(cube; labels = nothing, kwargs...) -> PlotData

A stack of images as a grid of panels, one per slice of an `(nx, ny, nslices)` array.

Each panel is normalised to **its own** maximum, which is what `imdisp_multi` does and is
worth knowing: it makes the faint channels readable and makes the panels non-comparable in
brightness at the same time. Pass `shared_scale = true` for one normalisation across the cube
when the comparison is the point.

Keywords otherwise as [`imdisp_makie`](@ref); `labels` titles the panels and defaults to the
slice number.
"""
function imdisp_multi_makie(cube::AbstractArray{<:Real,3};
                           labels = nothing, pixsize::Real = -1.0,
                           colormap = :gist_heat, colorbar::Bool = false,
                           beamsize::Real = -1.0, beamlocation = (0.8, 0.8),
                           shared_scale::Bool = false, scale::Real = 1.0)
    n = size(cube, 3)
    n > 0 || throw(ArgumentError("empty cube"))
    ttl = labels === nothing ? string.(1:n) : String.(labels)
    length(ttl) == n ||
        throw(ArgumentError("got $(length(ttl)) labels for $n slices"))

    ncol = ceil(Int, sqrt(n))
    nrow = ceil(Int, n / ncol)
    fig  = Makie.Figure(fonts = PLOT_FONTS)

    # One divisor for the whole cube when the panels are meant to be compared. Applied by
    # pre-scaling each slice, so `image_into!` stays the single drawing implementation.
    gmax = shared_scale ? maximum(cube) : nothing
    axes = Any[]
    for k in 1:n
        r, c = fld(k - 1, ncol) + 1, mod(k - 1, ncol) + 1
        ax = Makie.Axis(fig[r, c]; aspect = 1)
        image_into!(fig, ax, @view(cube[:, :, k]); pixsize, colormap, title = ttl[k],
                    colorbar = false, beamsize, beamlocation, scale, normalise = gmax)
        push!(axes, ax)
    end
    colorbar && n > 0 && Makie.Colorbar(fig[1:nrow, ncol + 1];
                                        colormap = colormap, label = "flux / pixel")
    return PlotData(fig, axes, nothing, String[])
end

# ─────────────────────────────────────────────────────────────────────────────
# Residuals
# ─────────────────────────────────────────────────────────────────────────────

"""
Observables a residual figure draws, and the `model_to_obs` field each one reads.

Closure quantities are keyed to the `_max` specs — residuals are plotted against the LONGEST
leg of the triangle, not its geometric mean. `plot_residuals` does the same, and the two
choices put visibly different points on the x axis, so this is not a detail.
"""
const RESIDUAL_KINDS = ((:v2, :v2), (:t3amp_max, :t3amp), (:t3phi_max, :t3phi),
                        (:visamp, :visamp), (:visphi, :visphi))

"Observables whose residual is an ANGLE, and so has to be wrapped before it is divided by σ."
const PHASE_KINDS = Set([:t3phi_max, :visphi])

"""
    plot_residuals_makie(data, obs; color = :baseline, logscale = false) -> PlotData
    plot_residuals_makie(model::FlatModel, x, data; kwargs...)
    plot_residuals_makie(x, ft, data; kwargs...)

Data against model, with normalised residuals beneath, for every observable present.

`obs` is the NamedTuple `model_to_obs` or `image_to_obs` returns; the two convenience methods
build it. Each observable gets a tall panel of data and model together and a short one of
`(model − data) / σ` under it, sharing an x axis — the 3:1 split `plot_residuals` uses.

A phase residual is wrapped through `mod360` BEFORE dividing by σ. Without that a model 359°
from the data reads as a 359σ outlier rather than the 1° agreement it is.

`color` is `:baseline`, `:wav` or `:mjd`, as everywhere else.
"""
function plot_residuals_makie(data, obs::NamedTuple; color::Symbol = :baseline,
                         logscale::Bool = false, scale::Real = 1.0)
    counts = (; v2 = data.nv2, t3amp_max = data.nt3amp, t3phi_max = data.nt3phi,
                visamp = data.nvisamp, visphi = data.nvisphi)
    present = [(k, f) for (k, f) in RESIDUAL_KINDS
               if getfield(counts, k) > 0 && hasproperty(obs, f) &&
                  length(getfield(obs, f)) > 0]
    isempty(present) &&
        throw(ArgumentError("no observable in $(basename(data.filename)) has both data and a model"))

    fig = Makie.Figure(fonts = PLOT_FONTS)
    axes = Any[]
    for (row, (kind, field)) in enumerate(present)
        spec  = OBS_SPECS[kind]
        ydata = Float64.(getfield(data, spec.y))
        yerr  = Float64.(getfield(data, spec.yerr))
        ymod  = Float64.(getfield(obs, field))
        xs    = Float64.(getfield(data, spec.x)) .* spec.xscale
        names = group_names(data, spec)
        cols  = _colors_for(data, color, names, spec)

        # 3:1, and both panels on one x axis: a residual is read against the point above it.
        ax_d = Makie.Axis(fig[2row - 1, 1])
        ax_r = Makie.Axis(fig[2row, 1])
        Makie.rowsize!(fig.layout, 2row - 1, Makie.Relative(3 / (4 * length(present))))
        Makie.rowsize!(fig.layout, 2row,     Makie.Relative(1 / (4 * length(present))))

        # On a log axis the non-positive points are dropped from the DATA panel only. The
        # residual panel underneath is in units of σ and stays linear, so it keeps every point
        # — including the ones the panel above cannot show.
        uselog = logscale && kind in LOG_Y_KINDS
        dx, dy, de, dm, dc = if uselog
            keep, ky, ke, km = _poslog(ydata, yerr, ymod)
            xs[keep], ky, ke, km, cols isa AbstractVector ? cols[keep] : cols
        else
            xs, ydata, yerr, ymod, cols
        end
        _errorbars!(ax_d, dx, dy, de; logscale = uselog)
        Makie.scatter!(ax_d, dx, dy; color = dc, marker = OIPLOT_MARKER, markersize = 7)
        Makie.scatter!(ax_d, dx, dm; color = :black, marker = :xcross, markersize = 7)
        style_axis!(ax_d; scale)
        ax_d.ylabel = spec.ylabel
        ax_d.xticklabelsvisible = false
        uselog && _log_yaxis!(ax_d, dy)

        resid = kind in PHASE_KINDS ? mod360(ymod .- ydata) ./ yerr : (ymod .- ydata) ./ yerr
        Makie.hlines!(ax_r, [0.0]; color = (:black, 0.4), linewidth = 1)
        Makie.scatter!(ax_r, xs, resid; color = cols, marker = OIPLOT_MARKER, markersize = 6)
        style_axis!(ax_r; scale)
        ax_r.ylabel = "σ"
        # Only the bottom panel carries the x label; the rest would repeat it five times.
        ax_r.xlabel = row == length(present) ? spec.xlabel : ""
        row < length(present) && (ax_r.xticklabelsvisible = false)

        Makie.linkxaxes!(ax_d, ax_r)
        push!(axes, ax_d, ax_r)
    end
    return PlotData(fig, axes, nothing, String[])
end

plot_residuals_makie(model::FlatModel, x::AbstractVector, data; kwargs...) =
    plot_residuals_makie(data, OITOOLS.model_to_obs(model, x, data); kwargs...)

plot_residuals_makie(x, ft, data; kwargs...) =
    plot_residuals_makie(data, OITOOLS.image_to_obs(x, ft, data); kwargs...)

# ─────────────────────────────────────────────────────────────────────────────
# Multi-panel, multi-file and array layout
# ─────────────────────────────────────────────────────────────────────────────

"""
    _log_yaxis!(ax, y) -> Bool

Put `ax` on a log y scale with limits and ticks valid for one, or leave it linear and return
`false` when nothing is positive.

Extracted because `draw!`, [`plot_residuals_makie`](@ref) and [`plot_v2_multifile_makie`](@ref) all need
it and getting it wrong is not a cosmetic error: Makie validates limits against the scale the
moment it is set, so a V² array containing zeros throws `Invalid y-limits (0.0, 10.0) for scale
log10` rather than drawing anything.

Setting the scale is only half of it. Makie also applies `log10` to every datum when the figure
is rendered, so the non-positive POINTS have to be left out of the plot as well — see
[`_poslog`](@ref), and `draw!`, which does the same thing inline.
"""
function _log_yaxis!(ax, y)
    pos = filter(v -> isfinite(v) && v > 0, y)
    isempty(pos) && return false
    lo, hi = extrema(pos)
    Makie.ylims!(ax, lo / 2, hi * 2)
    ax.yscale[] = log10
    # LinearTicks on a log axis bunches every label into the top decade; see livecanvas.jl.
    ax.yticks[] = Makie.LogTicks(Makie.LinearTicks(OIPLOT_XTICKS))
    return true
end

"""
    _poslog(y, others...) -> (mask, y[mask], others[mask]...)

The mask of points a log y axis can draw, applied to `y` and to everything plotted beside it.

matplotlib silently omits non-positive points under `set_yscale("log")`; Makie throws a
`DomainError` out of the render instead, which surfaces at `save` time rather than at the call
that caused it. Dropping them here matches oiplot and keeps the error at the right place.
"""
function _poslog(y, others...)
    keep = [isfinite(v) && v > 0 for v in y]
    return (keep, y[keep], (o[keep] for o in others)...)
end

"""
    _errorbars!(ax, x, y, e; logscale)

Error bars that survive a log axis.

A point can be positive while `y - e` is not, and Makie transforms the whole bar, so an
unclamped lower whisker throws a `DomainError` out of the render even after [`_poslog`](@ref)
has dropped the non-positive points. The whisker is truncated at a decade below the smallest
plotted value, which is what `draw!` does and what a log axis shows anyway.
"""
function _errorbars!(ax, x, y, e; logscale::Bool)
    if logscale && !isempty(y)
        floorv = minimum(y) / 10
        Makie.errorbars!(ax, x, y, y .- max.(y .- e, floorv), e; color = OIPLOT_ECOLOR)
    else
        Makie.errorbars!(ax, x, y, e; color = OIPLOT_ECOLOR)
    end
end

"""
    plot_obs_makie(data; kinds = nothing, color = :baseline, logscale = false) -> PlotData

Several observables stacked in one figure on a shared x axis — the counterpart of `plot_obs`
(and of `plot_multi`, which is an alias for it).

`kinds` defaults to every observable the data actually contains, in the order `OBS_SPECS`
lists them. Only the bottom panel carries the x label and its tick labels; repeating them
five times is what makes a stacked figure unreadable.

Each panel is drawn by `draw!` — the same implementation `plot_observable_makie` uses — so
a panel here and a single-observable figure cannot disagree.
"""
function plot_obs_makie(d; kinds = nothing, color::Union{Nothing,Symbol} = :baseline,
                            logscale::Bool = false, scale::Real = 1.0)
    counts = (; v2 = d.nv2, t3amp = d.nt3amp, t3phi = d.nt3phi,
                visamp = d.nvisamp, visphi = d.nvisphi, flux = d.nflux)
    want = kinds === nothing ?
           [k for k in (:v2, :t3amp, :t3phi, :visamp, :visphi, :flux)
            if hasproperty(counts, k) && getfield(counts, k) > 0] :
           collect(Symbol.(kinds))
    isempty(want) &&
        throw(ArgumentError("no observable to plot in $(basename(d.filename))"))

    fig  = Makie.Figure(fonts = PLOT_FONTS)
    axes = Any[]
    extras = Any[]
    for (row, kind) in enumerate(want)
        ax = Makie.Axis(fig[row, 1])
        draw!(fig, ax, d, kind; color, logscale, extras)
        ax.title = ""                       # one title for the figure, not one per panel
        if row < length(want)
            ax.xlabel = ""
            ax.xticklabelsvisible = false
        end
        push!(axes, ax)
    end
    length(axes) > 1 && Makie.linkxaxes!(axes...)
    Makie.Label(fig[0, 1], basename(d.filename); fontsize = OIPLOT_LABELSIZE * scale)
    return PlotData(fig, axes, nothing, String[])
end

"""
    plot_v2_multifile_makie(datasets; logscale = false) -> PlotData

V² from several datasets on one axis — the counterpart of `plot_v2_multifile`.

**One legend entry per baseline, not per file.** A baseline that appears in five nights is one
colour and one entry; the alternative is a legend longer than the plot, which is what makes
the matplotlib version hard to read with more than a couple of files.

That also means the figure answers "does this baseline agree between nights", which is the
question a multi-file V² plot is for. Colouring by file instead would answer a different one,
and neither `plot_v2_multifile` nor this offers it.
"""
function plot_v2_multifile_makie(datasets::AbstractVector; logscale::Bool = false,
                             scale::Real = 1.0, markersize::Real = 6)
    isempty(datasets) && throw(ArgumentError("no datasets"))
    allnames = reduce(vcat, [baseline_names(d) for d in datasets])
    cmap = baseline_color_map(allnames)      # one palette across every file

    fig = Makie.Figure(fonts = PLOT_FONTS)
    ax  = Makie.Axis(fig[1, 1])
    for d in datasets
        names = baseline_names(d)
        x = Float64.(d.v2_baseline) .* 1e-6
        y = Float64.(d.v2)
        e = Float64.(d.v2_err)
        c = [cmap[n] for n in names]
        if logscale
            keep, y, e = _poslog(y, e)
            x, c = x[keep], c[keep]
        end
        _errorbars!(ax, x, y, e; logscale)
        Makie.scatter!(ax, x, y; color = c, marker = OIPLOT_MARKER, markersize = markersize)
    end
    style_axis!(ax; scale)
    ax.xlabel = "Baseline (Mλ)"
    ax.ylabel = "V²"
    ax.title  = string(length(datasets), " files")
    # Across every file: one axis, so the limits have to cover all of them.
    logscale && _log_yaxis!(ax, reduce(vcat, [Float64.(d.v2) for d in datasets]))
    add_baseline_legend!(fig, ax, allnames)
    return PlotData(fig, ax, nothing, String[])
end

"""
    plot_facility_makie(facility) -> PlotData

Telescope positions of a `FacilityConfig`, labelled and to scale — the counterpart of
`plot_facility`.

Equal aspect, because an interferometric array plotted on unequal axes misreports the very
thing the figure is for: which baselines are long and in which direction.
"""
function plot_facility_makie(facility; scale::Real = 1.0)
    xyz = facility.sta_xyz'                 # (3, ntel)
    x = Float64.(xyz[1, :])
    y = Float64.(xyz[2, :])

    fig = Makie.Figure(fonts = PLOT_FONTS)
    ax  = Makie.Axis(fig[1, 1]; aspect = Makie.DataAspect())
    Makie.scatter!(ax, x, y; color = :black, marker = OIPLOT_MARKER, markersize = 10)
    Makie.text!(ax, x, y; text = String.(facility.sta_names),
                offset = (6, 6), fontsize = OIPLOT_LABELSIZE * scale)
    style_axis!(ax; scale)
    ax.xlabel = "x (m)"
    ax.ylabel = "y (m)"
    ax.title  = string(facility.name, "  —  ", facility.ntel, " telescopes")
    return PlotData(fig, ax, nothing, String[])
end
