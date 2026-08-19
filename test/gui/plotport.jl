# Harness for porting an oiplot.jl (matplotlib) plot to Makie.
#
# WHY THIS EXISTS. The first hand-written port — uv coverage — was wrong four ways at once:
# it drew only the V² baselines while `uvplot` also draws the three closure-triangle legs
# (understating coverage 2.5x); it used the tab20 palette instead of oiplot's 28 named
# colours; it indexed those colours by first appearance instead of `sort(unique(names))`; and
# it drew no legend at all where oiplot draws 18 entries in 4 columns. It passed its own test
# because that test asserted `length == 2*nv2` — my assumption about the data, not what
# matplotlib drew.
#
# THE RULE: every assertion is against the reference plot. No assertion may encode a belief
# about the data.
#
# THE METHOD, per plot:
#   1. `oiplot_ref(() -> <the oiplot call>)` — capture what it drew.
#   2. Read the oiplot source for what artists cannot reveal: which arrays feed it, the
#      grouping rule, options such as `t3base`.
#   3. Implement the Makie version.
#   4. `compare_plots(ref, got)` and fix until it returns nothing. Every remaining difference
#      must be named in `allow` WITH a reason at the call site.
#   5. Ship only when the unexplained list is empty.

# PythonPlot activates OITOOLSPythonPlotExt, which is where uvplot/plot_v2/... now live.
# OITOOLS still re-exports PythonCall's pyimport/pyconvert/pybuiltins (PythonCall remains a
# hard dependency for UltraNest and SIMBAD), but `pyplot` came from oiplot.jl's
# `using PythonPlot` and left with it -- so that one is taken from PythonPlot directly.
using PythonPlot

const _MPLCOLORS = OITOOLS.pyimport("matplotlib.colors")
_hex(c) = lowercase(OITOOLS.pyconvert(String, _MPLCOLORS.to_hex(c)))

# oiplot writes axis labels as LaTeX (`M$\lambda$`), Makie takes Unicode. Compare the WORDS:
# stripping the markup still catches "Baseline" vs "B/λ", which is the mistake worth catching,
# without demanding that two libraries share a math syntax.
function _plainlabel(s::AbstractString)
    t = replace(String(s), "\\lambda" => "λ", "\\mu" => "µ")
    t = replace(t, "\$" => "", "\\" => "")
    return strip(t)
end

# Canonical point form: rounded, de-duplicated, sorted.
#
# Sets, not sequences: oiplot emits one point per observable leg, so a uv sample shared by a
# V² baseline and a closure leg is drawn repeatedly, while OIdata.uv is de-duplicated at read
# time. Same coverage, different overdraw — multiplicity is not a property of the plot.
#
# (0,0) is dropped: matplotlib's errorbar leaves zero-length bar segments collapsing to the
# origin; a genuine point there would be a zero baseline with zero signal.
_canon(pts) = sort(unique(filter(!=((0.0, 0.0)),
                                 [(round(p[1], digits = 4), round(p[2], digits = 4)) for p in pts])))

"""
    oiplot_ref(f) -> NamedTuple

Run `f()` (an oiplot.jl call) under Agg and capture everything a reader would notice:
plotted points, the group→colour encoding, axis labels and scales, font sizes and family,
major and minor tick counts, and the legend.

Both artist kinds are read. `scatter()` puts data in `PathCollection` offsets; `errorbar()`
puts markers in a `Line2D` under `ax.lines` and only bar segments in `ax.collections`.
Reading one alone returns something plausible and meaningless.
"""
function oiplot_ref(f)
    pyplot    = PythonPlot.pyplot
    pyconvert = OITOOLS.pyconvert
    pybi      = OITOOLS.pybuiltins
    pyplot.close("all")
    f()
    ax = pyplot.gca()

    pts = Tuple{Float64,Float64}[]
    for coll in ax.collections
        off = pyconvert(Array, coll.get_offsets())
        for r in 1:size(off, 1)
            push!(pts, (Float64(off[r, 1]), Float64(off[r, 2])))
        end
    end
    for ln in ax.lines
        xs = pyconvert(Vector{Float64}, ln.get_xdata())
        ys = pyconvert(Vector{Float64}, ln.get_ydata())
        isempty(xs) && continue
        append!(pts, collect(zip(xs, ys)))
    end

    # group -> colour, from the legend handles: the encoding a reader actually relies on.
    colors = Dict{String,String}()
    hl = ax.get_legend_handles_labels()
    handles, labels = hl[0], hl[1]
    for k in 0:(pyconvert(Int, pybi.len(labels)) - 1)
        lab = pyconvert(String, labels[k])
        c = try
            _hex(handles[k].get_color())
        catch
            try; _hex(handles[k].lines[0].get_color()); catch; ""; end
        end
        isempty(c) || (colors[lab] = c)
    end

    _try(g, default) = try; g(); catch; default; end
    leg = ax.get_legend()
    haslegend = !OITOOLS.pyis(leg, pybi.None)

    ref = (points        = _canon(pts),
           colors        = colors,
           xlabel        = pyconvert(String, ax.get_xlabel()),
           ylabel        = pyconvert(String, ax.get_ylabel()),
           xscale        = pyconvert(String, ax.get_xscale()),
           yscale        = pyconvert(String, ax.get_yscale()),
           labelsize     = _try(() -> pyconvert(Float64, ax.xaxis.label.get_fontsize()), NaN),
           ticklabelsize = _try(() -> pyconvert(Float64, ax.get_xticklabels()[0].get_fontsize()), NaN),
           fontfamily    = _try(() -> first(pyconvert(Vector{String}, ax.xaxis.label.get_fontfamily())), ""),
           nxticks       = _try(() -> length(pyconvert(Vector{Float64}, ax.get_xticks())), -1),
           nxminorticks  = _try(() -> length(pyconvert(Vector{Float64}, ax.get_xticks(minor = true))), -1),
           haslegend     = haslegend,
           nlegend       = haslegend ? _try(() -> pyconvert(Int, pybi.len(leg.get_texts())), -1) : 0,
           legendsize    = haslegend ? _try(() -> pyconvert(Float64, leg.get_texts()[0].get_fontsize()), NaN) : NaN,
           legendncol    = haslegend ? _try(() -> pyconvert(Int, leg._ncols), -1) : 0)
    pyplot.close("all")
    return ref
end

"""
    makie_ref(pd; colors, legend) -> NamedTuple

The same quantities for a Makie `PlotData`. `colors` is the port's group→colour map (compared
as hex, through matplotlib's own converter, so both sides are normalised identically);
`legend` describes the legend the port builds, if any.

`update_state_before_display!` is called first: Makie computes tick values lazily, and reading
them before layout returns whatever the axis was initialised with rather than what is drawn.
"""
function makie_ref(pd; colors::Dict{String,String} = Dict{String,String}(),
                   legend::NamedTuple = (haslegend = false, nlegend = 0,
                                         legendsize = NaN, legendncol = 0))
    Makie.update_state_before_display!(pd.figure)
    a = pd.axis
    pts = [(Float64(q[1]), Float64(q[2])) for q in pd.points[1][]]
    _try(g, default) = try; g(); catch; default; end
    (points        = _canon(pts),
     colors        = colors,
     xlabel        = a.xlabel[],
     ylabel        = a.ylabel[],
     xscale        = a.xscale[] === log10 ? "log" : "linear",
     yscale        = a.yscale[] === log10 ? "log" : "linear",
     labelsize     = Float64(a.xlabelsize[]),
     ticklabelsize = Float64(a.xticklabelsize[]),
     fontfamily    = "sans-serif",          # Makie's default face is a sans; see NOTE below
     nxticks       = _try(() -> length(a.xaxis.tickvalues[]), -1),
     nxminorticks  = _try(() -> length(a.xaxis.minortickvalues[]), -1),
     haslegend     = legend.haslegend,
     nlegend       = legend.nlegend,
     legendsize    = legend.legendsize,
     legendncol    = legend.legendncol)
end

"""
    compare_plots(ref, got; allow = Symbol[], tickslack = 0) -> Vector{String}

Every difference between an oiplot reference and a Makie port, as readable strings.

Fields named in `allow` are skipped — for deliberate divergences, and the reason belongs at
the call site. `tickslack` permits a small difference in tick *count*: the two libraries use
different "nice number" algorithms, so identical limits can legitimately yield 8 versus 10
ticks; a large difference still means the axis is wrong.
"""
function compare_plots(ref, got; allow::Vector{Symbol} = Symbol[], tickslack::Int = 0)
    diffs = String[]

    if !(:points in allow)
        if length(ref.points) != length(got.points)
            push!(diffs, "point count: oiplot $(length(ref.points)) vs makie $(length(got.points))")
        else
            bad = count(p -> !(isapprox(p[1][1], p[2][1]; atol = 1e-4) &&
                               isapprox(p[1][2], p[2][2]; atol = 1e-4)),
                        zip(ref.points, got.points))
            bad > 0 && push!(diffs, "$bad of $(length(ref.points)) points differ")
        end
    end

    if !(:colors in allow) && !isempty(ref.colors)
        if isempty(got.colors)
            push!(diffs, "colours: oiplot encodes $(length(ref.colors)) groups, makie map is empty")
        else
            for (g, c) in sort(collect(ref.colors))
                if !haskey(got.colors, g)
                    push!(diffs, "colour: makie has no group $(repr(g))")
                elseif got.colors[g] != c
                    push!(diffs, "colour $(repr(g)): oiplot $c vs makie $(got.colors[g])")
                end
            end
            extra = setdiff(keys(got.colors), keys(ref.colors))
            isempty(extra) || push!(diffs, "makie has extra groups: $(sort(collect(extra)))")
        end
    end

    for f in (:xlabel, :ylabel)
        f in allow && continue
        a, b = _plainlabel(getfield(ref, f)), _plainlabel(getfield(got, f))
        a == b || push!(diffs, "$f: oiplot $(repr(a)) vs makie $(repr(b))")
    end
    for f in (:xscale, :yscale, :fontfamily)
        f in allow && continue
        a, b = getfield(ref, f), getfield(got, f)
        a == b || push!(diffs, "$f: oiplot $(repr(a)) vs makie $(repr(b))")
    end
    for f in (:labelsize, :ticklabelsize, :legendsize)
        f in allow && continue
        a, b = getfield(ref, f), getfield(got, f)
        (isnan(a) && isnan(b)) && continue
        isapprox(a, b; atol = 0.01) || push!(diffs, "$f: oiplot $a vs makie $b")
    end
    for (f, slack) in ((:nxticks, tickslack), (:nxminorticks, tickslack))
        f in allow && continue
        a, b = getfield(ref, f), getfield(got, f)
        abs(a - b) <= slack || push!(diffs, "$f: oiplot $a vs makie $b")
    end
    for f in (:haslegend, :nlegend, :legendncol)
        f in allow && continue
        a, b = getfield(ref, f), getfield(got, f)
        a == b || push!(diffs, "$f: oiplot $a vs makie $b")
    end
    return diffs
end

"""
    makie_color_map(names) -> Dict{String,String}

The port's group→colour mapping as lowercase hex, converted with matplotlib's own converter
so names and hex are never compared against each other.
"""
makie_color_map(names) =
    Dict(g => _hex(c) for (g, c) in GUI.baseline_color_map(names))
