#
# TO DO: plot_t3phivsdata, plot_t3phivsmodel, everything should be able to take as input a 1D array of OIdata
#

# gather common display tasks
using PyPlot,PyCall, LaTeXStrings, Statistics

# ── Styling globals (set by set_oiplot_defaults) ──────────────────────────────
global oiplot_markersize    = 3.0
global oiplot_scatter_size  = 6.0
global oiplot_elinewidth    = 1.0
global oiplot_cbar_pad      = 0.18
global oiplot_cbar_fraction = 0.05
global oiplot_legend_fontsize = 8
global oiplot_legend_ncol     = 4
global oiplot_legend_ncol_below = 8
global oiplot_show_title    = true
global oiplot_compact       = false
global oiplot_figsize       = (12, 6)

"""
    set_oiplot_defaults(; compact=oiplot_compact)

Apply consistent matplotlib style settings. When `compact=true`, all sizes are
reduced for stacking multiple plots in a single figure. Calling without arguments
preserves the current compact state.
"""
function set_oiplot_defaults(; compact::Bool=oiplot_compact)
    global oiplot_compact = compact
    rc = PyDict(pyimport("matplotlib")."rcParams")
    rc["font.family"] = ["serif"]
    rc["legend.numpoints"] = [1]
    rc["legend.handletextpad"] = [0.3]
    if compact
        rc["font.size"]              = [9]
        rc["xtick.major.size"]       = [3]
        rc["ytick.major.size"]       = [3]
        rc["xtick.minor.size"]       = [3]
        rc["ytick.minor.size"]       = [3]
        rc["xtick.major.width"]      = [0.5]
        rc["ytick.major.width"]      = [0.5]
        rc["xtick.minor.width"]      = [0.5]
        rc["ytick.minor.width"]      = [0.5]
        rc["lines.markeredgewidth"]  = [0.5]
        global oiplot_markersize       = 1.5
        global oiplot_scatter_size     = 3.0
        global oiplot_elinewidth       = 0.5
        global oiplot_cbar_pad         = 0.10
        global oiplot_cbar_fraction    = 0.03
        global oiplot_legend_fontsize  = 6
        global oiplot_legend_ncol      = 6
        global oiplot_legend_ncol_below = 10
        global oiplot_show_title       = false
        global oiplot_figsize          = (8, 3)
    else
        rc["font.size"]              = [14]
        rc["xtick.major.size"]       = [6]
        rc["ytick.major.size"]       = [6]
        rc["xtick.minor.size"]       = [6]
        rc["ytick.minor.size"]       = [6]
        rc["xtick.major.width"]      = [1]
        rc["ytick.major.width"]      = [1]
        rc["xtick.minor.width"]      = [1]
        rc["ytick.minor.width"]      = [1]
        rc["lines.markeredgewidth"]  = [1]
        global oiplot_markersize       = 3.0
        global oiplot_scatter_size     = 6.0
        global oiplot_elinewidth       = 1.0
        global oiplot_cbar_pad         = 0.18
        global oiplot_cbar_fraction    = 0.05
        global oiplot_legend_fontsize  = 8
        global oiplot_legend_ncol      = 4
        global oiplot_legend_ncol_below = 8
        global oiplot_show_title       = true
        global oiplot_figsize          = (12, 6)
    end
end

const global oiplot_colors=["black", "gold","chartreuse","blue","red", "pink","lightgray","darkorange","darkgreen","aqua",
"fuchsia","saddlebrown","dimgray","darkslateblue","violet","indigo","blue","dodgerblue",
"sienna","olive","purple","darkorchid","tomato","darkturquoise","steelblue","seagreen","darkgoldenrod","darkseagreen","salmon","slategray","lime","coral","maroon","mistyrose","sandybrown","tan","olivedrab"]

const global oiplot_markers=["o","s","v","P","*","x","^","D","p",1,"<","H","X","4",4,"_","1",6,"8","d",9]

# ── Data normalization ────────────────────────────────────────────────────────
"""
    as_datavec(data) -> Vector{OIdata}

Ensure `data` is a flat `Vector{OIdata}`, wrapping a single `OIdata` or
flattening a 2-D array as needed.
"""
function as_datavec(data::Union{OIdata, AbstractArray{<:OIdata}})
    data isa OIdata && return [data]
    ndims(data) == 2 && return vec(data)
    return data
end

# ── Generic observable plotting infrastructure ────────────────────────────────

struct ObsPlotSpec
    plot_title::String
    ylabel::String
    xlabel::String
    y_field::Symbol          # :v2, :visamp, :t3phi, ...
    yerr_field::Symbol       # :v2_err, :visamp_err, ...
    x_field::Symbol          # :v2_baseline, :vis_baseline, :flux_lam, ...
    lam_field::Symbol        # :v2_lam, :vis_lam, :t3_lam, :flux_lam
    mjd_field::Symbol        # :v2_mjd, :vis_mjd, :t3_mjd, :flux_mjd
    sta_index_field::Symbol  # :v2_sta_index, :vis_sta_index, :flux_sta_index
    grouping::Symbol         # :baseline, :triplet, or :station
    logplot_ok::Bool
    x_scale::Float64         # multiplicative factor for x values (1e-6 for Mλ, 1e6 for μm)
end

const OBS_PLOT_SPECS = Dict{String, ObsPlotSpec}(
    "V2" => ObsPlotSpec(
        "V²", "V²",
        L"Baseline (M$\lambda$)",
        :v2, :v2_err, :v2_baseline, :v2_lam, :v2_mjd,
        :v2_sta_index, :baseline, true, 1e-6),
    "T3PHI" => ObsPlotSpec(
        "T3φ", "T3φ (°)",
        L"Baseline (M$\lambda$)",
        :t3phi, :t3phi_err, :t3_baseline, :t3_lam, :t3_mjd,
        :t3_sta_index, :triplet, false, 1e-6),
    "T3PHI_MAX" => ObsPlotSpec(
        "T3φ", "T3φ (°)",
        L"Max. baseline (M$\lambda$)",
        :t3phi, :t3phi_err, :t3_maxbaseline, :t3_lam, :t3_mjd,
        :t3_sta_index, :triplet, false, 1e-6),
    "T3AMP" => ObsPlotSpec(
        "T3amp", "T3amp",
        L"Baseline (M$\lambda$)",
        :t3amp, :t3amp_err, :t3_baseline, :t3_lam, :t3_mjd,
        :t3_sta_index, :triplet, true, 1e-6),
    "T3AMP_MAX" => ObsPlotSpec(
        "T3amp", "T3amp",
        L"Max. baseline (M$\lambda$)",
        :t3amp, :t3amp_err, :t3_maxbaseline, :t3_lam, :t3_mjd,
        :t3_sta_index, :triplet, true, 1e-6),
    "VISAMP" => ObsPlotSpec(
        "Visamp", "Visamp",
        L"Baseline (M$\lambda$)",
        :visamp, :visamp_err, :vis_baseline, :vis_lam, :vis_mjd,
        :vis_sta_index, :baseline, true, 1e-6),
    "VISPHI" => ObsPlotSpec(
        "Visφ", "Visφ (°)",
        L"Baseline (M$\lambda$)",
        :visphi, :visphi_err, :vis_baseline, :vis_lam, :vis_mjd,
        :vis_sta_index, :baseline, false, 1e-6),
    "FLUX" => ObsPlotSpec(
        "Flux", "Flux",
        "Wavelength (μm)",
        :flux, :flux_err, :flux_lam, :flux_lam, :flux_mjd,
        :flux_sta_index, :station, false, 1e6),
)

# Normalize color aliases to a canonical form
function canonical_color(color::String)
    color in ("baseline", "base", "bases", "baselines") && return "baseline"
    color in ("station", "stations")                     && return "station"
    color in ("wavelength", "wav", "wavs", "wavelengths") && return "wav"
    color in ("mjd", "time", "timestamp", "timestamps")  && return "mjd"
    return color
end

# Place colorbar label to the right instead of below
function cbar_label_right!(cbar, label::String)
    fs = oiplot_compact ? 9 : 14
    cbar.ax.text(1.02, 0.5, label, transform=cbar.ax.transAxes,
                 va="center", ha="left", fontsize=fs, family="serif")
end

# Set up MJD colorbar ticks
function setup_mjd_colorbar(sc)
    cbar = colorbar(sc, aspect=50, orientation="horizontal",
                    pad=oiplot_cbar_pad, fraction=oiplot_cbar_fraction)
    cbar_label_right!(cbar, "MJD")
    mjdvals = sc.get_array()
    mjds = sort(unique(mjdvals))
    if length(mjds) < 5
        cbar.set_ticks(round.(mjds*100)/100)
    else
        cbar_range = round.(collect(range(minimum(mjdvals), maximum(mjdvals), length=5))*100)/100
        cbar.set_ticks(cbar_range)
    end
    return cbar
end

# Set up wavelength colorbar ticks
function setup_wav_colorbar(sc)
    cbar = colorbar(sc, aspect=50, orientation="horizontal",
                    pad=oiplot_cbar_pad, fraction=oiplot_cbar_fraction)
    cbar_label_right!(cbar, "λ (μm)")
    wavvals = sc.get_array()
    cbar_range = floor.(collect(range(minimum(wavvals), maximum(wavvals), length=11))*100)/100
    cbar.set_ticks(cbar_range)
    return cbar
end

"""
    plot_obs(data, spec; color="baseline", logplot=false, figsize=oiplot_figsize,
             figtitle="", markopt=false, legend_below=false)

Generic plotting engine for any interferometric observable described by an
`ObsPlotSpec`. All public `plot_*` functions delegate here.
"""
function plot_obs(data::Union{OIdata, AbstractArray{<:OIdata}}, spec::ObsPlotSpec;
                  color::String="baseline", logplot::Bool=false, figsize=oiplot_figsize,
                  figtitle::String="", markopt::Bool=false, legend_below::Bool=false,
                  plot_title::String=spec.plot_title, ylabel_str::String=spec.ylabel,
                  xlabel_str::String=spec.xlabel)
    set_oiplot_defaults()
    data = as_datavec(data)

    fig = figure(string(figtitle, plot_title, " data"), figsize=figsize, facecolor="White")
    clf()
    ax = gca()
    if logplot && spec.logplot_ok
        ax.set_yscale("log")
    end

    cc = canonical_color(color)
    xs = spec.x_scale

    if cc == "baseline" || cc == "station"
        # Group by baseline/triplet/station name
        if spec.grouping == :triplet
            label_lists = [get_triplet_names(data[n].sta_name, getfield(data[n], spec.sta_index_field)) for n in eachindex(data)]
        elseif spec.grouping == :station
            label_lists = [get_station_names(data[n].sta_name, getfield(data[n], spec.sta_index_field)) for n in eachindex(data)]
        else
            label_lists = [get_baseline_names(data[n].sta_name, getfield(data[n], spec.sta_index_field)) for n in eachindex(data)]
        end
        groups = sort(unique(vcat(label_lists...)))
        for i in eachindex(groups)
            loc = [findall(label_lists[n] .== groups[i]) for n in eachindex(data)]
            xvals = vcat([getfield(data[n], spec.x_field)[loc[n]] for n in eachindex(data)]...) * xs
            yvals = vcat([getfield(data[n], spec.y_field)[loc[n]] for n in eachindex(data)]...)
            yerr  = vcat([getfield(data[n], spec.yerr_field)[loc[n]] for n in eachindex(data)]...)
            if markopt
                errorbar(xvals, yvals, yerr=yerr, fmt="o", marker=oiplot_markers[i],
                         markeredgecolor="Black", markersize=oiplot_markersize, color="Black",
                         ecolor="Gainsboro", elinewidth=oiplot_elinewidth, label=groups[i])
            else
                errorbar(xvals, yvals, yerr=yerr, fmt="o",
                         markeredgecolor=oiplot_colors[i], color=oiplot_colors[i],
                         markersize=oiplot_markersize, ecolor="Gainsboro",
                         elinewidth=oiplot_elinewidth, label=groups[i])
            end
        end
        if legend_below
            ax.legend(fontsize=oiplot_legend_fontsize, fancybox=true, shadow=true,
                      ncol=oiplot_legend_ncol_below,
                      loc="upper center", bbox_to_anchor=(0.5, -0.15))
        else
            ax.legend(fontsize=oiplot_legend_fontsize, fancybox=true, shadow=true,
                      ncol=oiplot_legend_ncol, loc="best")
        end
    elseif cc == "wav"
        wavcol = vcat([getfield(data[n], spec.lam_field)*1e6 for n in eachindex(data)]...)
        xvals  = vcat([getfield(data[n], spec.x_field) for n in eachindex(data)]...) * xs
        yvals  = vcat([getfield(data[n], spec.y_field) for n in eachindex(data)]...)
        yerr   = vcat([getfield(data[n], spec.yerr_field) for n in eachindex(data)]...)
        sc = scatter(xvals, yvals, c=wavcol, cmap="Spectral_r", alpha=1.0,
                     s=oiplot_scatter_size, zorder=100)
        errorbar(xvals, yvals, yerr=yerr, fmt="none", marker="none",
                 ecolor="Gainsboro", elinewidth=oiplot_elinewidth, zorder=0)
        setup_wav_colorbar(sc)
    elseif cc == "mjd"
        mjdcol = vcat([getfield(data[n], spec.mjd_field) for n in eachindex(data)]...)
        xvals  = vcat([getfield(data[n], spec.x_field) for n in eachindex(data)]...) * xs
        yvals  = vcat([getfield(data[n], spec.y_field) for n in eachindex(data)]...)
        yerr   = vcat([getfield(data[n], spec.yerr_field) for n in eachindex(data)]...)
        sc = scatter(xvals, yvals, c=mjdcol, cmap="plasma", alpha=1.0,
                     s=oiplot_scatter_size, zorder=100)
        errorbar(xvals, yvals, yerr=yerr, fmt="none", marker="none",
                 ecolor="Gainsboro", elinewidth=oiplot_elinewidth, zorder=0)
        setup_mjd_colorbar(sc)
    else
        xvals = vcat([getfield(data[n], spec.x_field) for n in eachindex(data)]...) * xs
        yvals = vcat([getfield(data[n], spec.y_field) for n in eachindex(data)]...)
        yerr  = vcat([getfield(data[n], spec.yerr_field) for n in eachindex(data)]...)
        errorbar(xvals, yvals, yerr=yerr, fmt="o", markersize=oiplot_markersize,
                 color=color, ecolor="Gainsboro", elinewidth=oiplot_elinewidth)
    end

    if oiplot_show_title
        title(plot_title * " data")
    end
    xlabel(xlabel_str)
    PyPlot.ylabel(ylabel_str)
    ax.grid(true, which="both", color="Grey", linestyle=":")
    tight_layout()
    show(block=false)
end




# Overloaded uvplot functions
function uvplot(uv::Array{Float64,2})
    u = uv[1,:]/1e6
    v = uv[2,:]/1e6
    fig = figure("UV plot",figsize=(8,8),facecolor="White")
    set_oiplot_defaults()
    clf();
    ax = gca()
    markeredgewidth=0.1
    ax.locator_params(axis ="y", nbins=20)
    ax.locator_params(axis ="x", nbins=20)
    ax.set_aspect("equal")
    scatter(u, v,alpha=1.0, s = 12.0,color="Black")
    scatter(-u, -v,alpha=1.0, s = 12.0, color="Black")
    title("UV coverage")
    xlabel(L"U (M$\lambda$)")
    ylabel(L"V (M$\lambda$)")
    ax.grid(true,which="both",color="Grey", linestyle=":");
    tight_layout();
end

"""
    uvplot(data; kwargs...)

Plot the UV coverage for one or more `OIdata` bins. Conjugate baselines (-u, -v) are
always shown. Accepts a single `OIdata`, a vector, or the 2-D array returned by
`readoifits` (flattened automatically).

# Keyword arguments
- `color` — colouring scheme:
  - `"baseline"` (default) — one colour per baseline pair, labelled legend
  - `"wav"` / `"wavelength"` — coloured by λ (μm) with a horizontal colorbar
  - `"mjd"` / `"time"` — coloured by MJD with a horizontal colorbar
  - any other string — all points black
- `figsize` — figure size tuple. Default: `(10, 10)`.
- `minuv`, `maxuv` — axis limits in Mλ when `square=true`.
- `square` — force equal axis limits. Default: `true`.
- `legend_below` — place the baseline legend below the plot. Default: `true`.
- `figtitle` — plot title string.
- `cmap` — matplotlib colormap for `wav` / `mjd` modes. Default: `"Spectral_r"`.
- `flipx` — invert the x-axis (East to the right). Default: `false`.
- `filename` — save to file if non-empty.
"""
function uvplot(data::Union{OIdata, AbstractArray{<:OIdata}};color::String="baseline",filename="", figsize=(10,10), minuv= -1e99, maxuv= 1e99, square = true, legend_below = true, figtitle = "", windowtitle="", cmap="Spectral_r", flipx = false)
    set_oiplot_defaults()
    if data isa OIdata
        data = [data]
    end
    if ndims(data) == 2
        data = vec(data)
    end
    nuv = sum(data[i].nuv for i=1:length(data))
    mean_mjd = mean(data[i].mean_mjd for i=1:length(data))
    if windowtitle==""
        windowtitle = string("Mean MJD: $(round(mean_mjd*100)/100), nuv: $(nuv)")
    end
    fig = figure(windowtitle,figsize=figsize,facecolor="White")
    set_oiplot_defaults()
    clf();
    ax = gca()
    ax.locator_params(axis ="y", nbins=10)
    ax.locator_params(axis ="x", nbins=10)
    ax.set_aspect(1.0)
    if (color == "baseline" || color =="base" || color =="bases" || color =="baselines") # we need to identify corresponding baselines #TBD --> could be offloaded to readoifits
        baseline_list_v2 = [get_baseline_names(data[n].sta_name,data[n].v2_sta_index) for n=1:length(data)];
        baseline_list_t3 = [get_triple_baselines_names(data[n].sta_name,data[n].t3_sta_index) for n=1:length(data)];
        baseline=sort(unique(vcat(vcat(baseline_list_v2...),vec(hcat(baseline_list_t3...)))))
        if length(baseline)>length(oiplot_colors)
            @warn("I ran out of colors!")
        end
        indx_v2 = [data[n].indx_v2 for n=1:length(data)]
        indx_t3 = [hcat(data[n].indx_t3_1,data[n].indx_t3_2, data[n].indx_t3_3)' for n=1:length(data)]
        for i=1:length(baseline)
            loc =  [vcat(indx_v2[n][findall(baseline_list_v2[n] .== baseline[i])], indx_t3[n][findall(baseline_list_t3[n] .== baseline[i])]) for n=1:length(data)]
            u = vcat([data[n].uv[1,loc[n]] for n=1:length(data)]...)/1e6;
            v = vcat([data[n].uv[2,loc[n]] for n=1:length(data)]...)/1e6;
            scatter( u,  v, alpha=1.0, s=12.0, color=oiplot_colors[i],label=baseline[i]) #TBD: handle case where length(baseline)>length(oiplot_colors)
            scatter(-u, -v, alpha=1.0, s=12.0, color=oiplot_colors[i])
        end
               
        if legend_below == false
            ax.legend(fontsize=10, fancybox=true, shadow=true, ncol=3, loc="upper right")
        else
            ax.legend(fontsize=10, fancybox=true, shadow=true, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.10));
            tight_layout();
        end
    elseif (color == "wavelength" || color == "wav" || color =="wavs" || color =="wavelengths")
        u = vcat([data[n].uv[1,:]/1e6 for n=1:length(data)]...)
        v = vcat([data[n].uv[2,:]/1e6 for n=1:length(data)]...)
        wavcol = vcat([data[n].uv_lam*1e6 for n=1:length(data)]...)
        wavs=sort(unique(wavcol))
        scatter([u;-u], [v;-v], alpha=1.0, s = 12.0, c=[wavcol;wavcol], cmap=cmap)
        cbar = colorbar(ax=ax, aspect=50, orientation="horizontal", label="Wavelength (μm)", pad=0.1, fraction=0.02)
        if length(wavs)==2
            cbar_range = round.(collect(range(ceil(minimum(wavcol)*100)/100, floor(maximum(wavcol)*100)/100,length=2))*100)/100
            cbar.set_ticks(cbar_range)
            cbar.set_ticklabels(cbar_range)
        elseif length(wavs)==1
            cbar_range = wavs
            cbar.set_ticks(cbar_range)
            cbar.set_ticklabels(cbar_range)
        else # >= 9
            nwavs = min(length(wavs), 10) # max 10 divisions
            minrange = ceil(minimum(wavcol)*100)/100
            maxrange = floor(maximum(wavcol)*100)/100
            step = round((maxrange-minrange)/nwavs*100)/100
            cbar_range = collect(range(minrange, maxrange, step=step))
            cbar.set_ticks(cbar_range)
        end
    elseif (color == "mjd" || color == "time")
        u = vcat([data[n].uv[1,:]/1e6 for n=1:length(data)]...)
        v = vcat([data[n].uv[2,:]/1e6 for n=1:length(data)]...)
        mjdcol = vcat([data[n].uv_mjd for n=1:length(data)]...)
        scatter([u;-u], [v;-v], alpha=1.0, s = 12.0, c=[mjdcol;mjdcol], cmap=cmap)
        cbar = colorbar(ax=ax, aspect=50, orientation="horizontal", label="MJD", pad=0.1, fraction=0.02)
        mjds=unique(mjdcol)
        if length(mjds)<5
            cbar.set_ticks(round.(sort(mjds)*100)/100)
            cbar.set_ticklabels(round.(sort(mjds)*100)/100)
        else
            cbar_range = round.(collect(range(ceil(minimum(mjdcol)*100)/100, floor(maximum(mjdcol)*100)/100,length=5))*100)/100
            cbar.set_ticks(cbar_range)
            cbar.set_ticklabels(cbar_range)
        end
    else
        u = vcat([data[n].uv[1,:]/1e6 for n=1:length(data)]...)
        v = vcat([data[n].uv[2,:]/1e6 for n=1:length(data)]...)
        scatter(u, v,alpha=1.0, s = 12.0,color="Black")
        scatter(-u, -v,alpha=1.0, s = 12.0, color="Black")
    end
    if square==true
        if minuv<-1e98
            minuv = minimum([ax.get_xlim()[1],ax.get_ylim()[1]] )
        end
        if maxuv>1e98
            maxuv = maximum([ax.get_xlim()[2],ax.get_ylim()[2]])
        end     
        ax.set_xlim((minuv,maxuv))
        ax.set_ylim((minuv,maxuv))
    end
    title(figtitle)
    xlabel(L"U (M$\lambda$)")
    ylabel(L"V (M$\lambda$)")
    ax.grid(true,which="both",color="Grey",linestyle=":");
    ax.set_aspect("equal")
    if flipx == true
        ax.invert_xaxis()
    end
    if filename !=""
        savefig(filename)
    end
    tight_layout();
    show(block=false)
end


"""
    plot_v2(data; kwargs...)

Plot squared visibilities V² vs baseline length (Mλ). Accepts a single `OIdata`,
vector, or 2-D array.

# Keyword arguments
- `color` — `"baseline"` (default), `"wav"`, `"mjd"`, or an explicit colour string.
- `logplot` — logarithmic y-axis. Default: `false`.
- `figsize` — figure size. Default: `(12, 6)`.
- `figtitle` — prefix for the figure window title.
- `markopt` — use distinct marker shapes per baseline. Default: `false`.
- `legend_below` — place legend below the plot. Default: `false`.
"""
function plot_v2(data::Union{OIdata, AbstractArray{<:OIdata}}; figsize=oiplot_figsize,
                 logplot=false, color::String="baseline", markopt=false,
                 legend_below=false, figtitle="")
    plot_obs(data, OBS_PLOT_SPECS["V2"]; color=color, logplot=logplot,
             figsize=figsize, figtitle=figtitle, markopt=markopt,
             legend_below=legend_below)
end


"""
    plot_t3phi(data; kwargs...)

Plot closure phases T3φ (degrees) vs a representative baseline length (Mλ).
Accepts a single `OIdata`, vector, or 2-D array.

# Keyword arguments
- `t3base` — x-axis baseline: `"geom"` (default, geometric mean) or `"max"` (longest side).
- `color` — `"baseline"` (default), `"wav"`, `"mjd"`, or explicit colour string.
- `figsize` — figure size. Default: `(12, 6)`.
- `figtitle` — prefix for the figure window title.
- `markopt` — use distinct marker shapes. Default: `false`.
- `legend_below` — place legend below the plot. Default: `false`.
"""
function plot_t3phi(data::Union{OIdata, AbstractArray{<:OIdata}}; figsize=oiplot_figsize,
                    color::String="baseline", markopt=false, legend_below=false,
                    t3base="geom", figtitle="")
    spec = t3base == "max" ? OBS_PLOT_SPECS["T3PHI_MAX"] : OBS_PLOT_SPECS["T3PHI"]
    plot_obs(data, spec; color=color, figsize=figsize, figtitle=figtitle,
             markopt=markopt, legend_below=legend_below)
end

function plot_v2_residuals(x, data::OIdata, ft::Array{NFFT.NFFTPlan{Float64, 2, 1}, 1}; logplot = false, y_range=[], res_range=[])
    plot_v2_residuals(data, image_to_v2(x, data, ft), logplot = logplot, y_range=y_range, res_range=res_range);
end

function plot_t3phi_residuals(x, data::OIdata, ft::Array{NFFT.NFFTPlan{Float64, 2, 1}, 1}; logplot = false, y_range=[], res_range=[])
    plot_t3phi_residuals(data, image_to_t3phi(x, data, ft), logplot = logplot, y_range=y_range, res_range=res_range);
end

function plot_t3amp_residuals(x, data::OIdata, ft::Array{NFFT.NFFTPlan{Float64, 2, 1}, 1}; logplot = false, y_range=[], res_range=[])
    plot_t3amp_residuals(data, image_to_t3amp(x, data, ft), logplot = logplot, y_range=y_range, res_range=res_range);
end


function plot_v2_residuals(data::OIdata, v2_model::Array{Float64,1}; figsize=(12,12), logplot = false, y_range=[], res_range=[]) #plots V2 data vs v2 model
    set_oiplot_defaults()

    fig = figure("V2 plot - Model vs Data",figsize=figsize,facecolor="White");
    fig.subplots(nrows=2, ncols=1, sharex=true)
    ax = subplot(2, 1, 1)
    if logplot==true
        ax.set_yscale("log")
    end
    errorbar(data.v2_baseline/1e6,data.v2,yerr=data.v2_err,fmt="o", markersize=1, color="Grey", ecolor="Grey")
    plot(data.v2_baseline/1e6, v2_model, color="Red", linestyle="none", marker="o", markersize=2)
    title("Squared Visibility Amplitudes - Model vs data plot")
    if y_range != []
        ylim(y_range)
    end
    ylabel("Squared Visibility Amplitudes")
    ax.grid(true,which="both",color="Grey",linestyle=":");
    ax = subplot(2, 1, 2)
    ax.axhline(color="Grey")
    plot(data.v2_baseline/1e6, (v2_model - data.v2)./data.v2_err, color="Grey", linestyle="none", marker="o", markersize=2)
    xlabel(L"Baseline (M$\lambda$)")
    ylabel("Residuals (number of sigma)")
    if res_range !=[]
        ylim(res_range)
    end
    ax.grid(true,which="both",color="Grey",linestyle=":")
    tight_layout()
    show(block=false)
end

function plot_t3phi_residuals(data::OIdata, t3phi_model::Array{Float64,1}; figsize=(12,12), logplot = false, y_range=[],  res_range=[]) #plots V2 data vs v2 model
    set_oiplot_defaults()

    fig = figure("Closure phase plot - Model vs Data",figsize=figsize,facecolor="White");
    fig.subplots(nrows=2, ncols=1, sharex=true)
    ax = subplot(2, 1, 1)
    if logplot==true
        ax.set_yscale("log")
    end
    errorbar(data.t3_maxbaseline/1e6,data.t3phi,yerr=data.t3phi_err,fmt="o", markersize=1, color="Grey", ecolor="Grey")
    plot(data.t3_maxbaseline/1e6, t3phi_model, color="Red", linestyle="none", marker="o", markersize=2)
    title("Closure phases - Model vs data plot")
    if y_range != []
        ylim(y_range)
    end
    ylabel("Closure phases (degrees)")
    ax.grid(true,which="both",color="Grey",linestyle=":");
    ax = subplot(2, 1, 2)
    ax.axhline(color="Grey")
    plot(data.t3_maxbaseline/1e6, mod360(t3phi_model - data.t3phi)./data.t3phi_err,color="Grey", linestyle="none", marker="o", markersize=1)
    if res_range !=[]
        ylim(res_range)
    end
    xlabel(L"Max Baseline (M$\lambda$)")
    ylabel("Residuals (number of sigma)")
    ax.grid(true,which="both",color="Grey",linestyle=":")
    tight_layout()
    show(block=false)
end

function plot_t3amp_residuals(data::OIdata, t3amp_model::Array{Float64,1}; figsize=(12,12), logplot = false, y_range=[],  res_range=[]) #plots V2 data vs v2 model
    set_oiplot_defaults()

    fig = figure("Triple amplitudes plot - Model vs Data",figsize=figsize,facecolor="White");
    fig.subplots(nrows=2, ncols=1, sharex=true)
    ax = subplot(2, 1, 1)
    if logplot==true
        ax.set_yscale("log")
    end
    errorbar(data.t3_maxbaseline/1e6,data.t3amp,yerr=data.t3amp_err,fmt="o", markersize=1, color="Grey", ecolor="Grey")
    plot(data.t3_maxbaseline/1e6, t3amp_model, color="Red", linestyle="none", marker="o", markersize=2)
    title("Triple amplitudes - Model vs data plot")
    if y_range != []
        ylim(y_range)
    end
    ylabel("Triple amplitudes")
    ax.grid(true,which="both",color="Grey",linestyle=":");
    ax = subplot(2, 1, 2)
    ax.axhline(color="Grey")
    plot(data.t3_maxbaseline/1e6, (t3amp_model - data.t3amp)./data.t3amp_err,color="Grey", linestyle="none", marker="o", markersize=2)
    if res_range !=[]
        ylim(res_range)
    end
    xlabel(L"Max Baseline (M$\lambda$)")
    ylabel("Residuals (number of sigma)")
    ax = gca();
    ax.grid(true,which="both",color="Grey",linestyle=":")
    tight_layout()
    show(block=false)
end


"""
    plot_t3amp(data; kwargs...)

Plot triple amplitudes T3amp vs a representative baseline length (Mλ).
Accepts a single `OIdata`, vector, or 2-D array.

# Keyword arguments
- `t3base` — x-axis baseline: `"geom"` (default, geometric mean) or `"max"` (longest side).
- `color` — `"baseline"` (default), `"wav"`, `"mjd"`, or explicit colour string.
- `logplot` — use logarithmic y-axis. Default: `false`.
- `figsize` — figure size. Default: `(12, 6)`.
- `figtitle` — prefix for the figure window title.
- `markopt` — use distinct marker shapes. Default: `false`.
- `legend_below` — place legend below the plot. Default: `false`.
"""
function plot_t3amp(data::Union{OIdata, AbstractArray{<:OIdata}}; figsize=oiplot_figsize,
                    color::String="baseline", markopt=false, legend_below=false,
                    t3base="geom", logplot=false, figtitle="")
    spec = t3base == "max" ? OBS_PLOT_SPECS["T3AMP_MAX"] : OBS_PLOT_SPECS["T3AMP"]
    plot_obs(data, spec; color=color, logplot=logplot, figsize=figsize,
             figtitle=figtitle, markopt=markopt, legend_below=legend_below)
end

"""
    plot_flux(data; kwargs...)

Plot flux spectra vs wavelength (μm). Accepts a single `OIdata`, vector, or 2-D array.

# Keyword arguments
- `color` — colouring scheme:
  - `"wav"` (default) — coloured by wavelength with a horizontal colorbar.
  - `"station"` — one colour per station; OI_FLUX entries with `CALSTAT=C`
    (i.e. `flux_sta_index == 0`) are labelled `"Calibrated"`.
  - `"mjd"` / `"time"` — scatter coloured by MJD with a horizontal colorbar.
  - any other string — treated as a matplotlib colour applied uniformly.
- `figsize` — figure size. Default: `(12, 6)`.
- `figtitle` — prefix for the figure window title.
- `markopt` — use black markers with distinct shapes. Default: `false`.
- `legend_below` — place legend below the plot. Default: `false`.
"""
function plot_flux(data::Union{OIdata, AbstractArray{<:OIdata}}; figsize=oiplot_figsize,
                   color::String="wav", markopt=false, legend_below=false, figtitle="")
    plot_obs(data, OBS_PLOT_SPECS["FLUX"]; color=color, figsize=figsize,
             figtitle=figtitle, markopt=markopt, legend_below=legend_below)
end

"""
    plot_visphi(data; kwargs...)

Plot visibility phases (degrees). The layout is chosen automatically from the
`PHITYP` header stored in `data.phityp`:

- **differential** (`phityp == "differential"`) — one subplot per baseline with
  wavelength (nm) on the x-axis.
- **absolute** (default) — single panel with baseline (Mλ) on the x-axis, coloured
  by baseline, wavelength, or MJD.

Accepts a single `OIdata`, vector, or 2-D array.

# Keyword arguments
- `color` — `"baseline"` (default), `"wav"`, `"mjd"`, or explicit colour string
  (ignored for differential layout).
- `figsize` — figure size. Default: `(12, 6)`.
- `figtitle` — prefix for the figure window title.
- `markopt` — use distinct marker shapes. Default: `false`.
- `legend_below` — place legend below the plot. Default: `false`.
"""
function plot_visphi(data::Union{OIdata, AbstractArray{<:OIdata}}; figsize=oiplot_figsize,
                     color::String="baseline", markopt=false, legend_below=false, figtitle="")
    dvec = as_datavec(data)
    is_diff = any(d.phityp == "differential" for d in dvec)

    if is_diff
        # ── Differential phase layout: one subplot per baseline, wavelength on x-axis ──
        set_oiplot_defaults()
        baseline_list_vis = [get_baseline_names(dvec[n].sta_name, dvec[n].vis_sta_index) for n in eachindex(dvec)]
        baselines = sort(unique(vcat(baseline_list_vis...)))
        nbase = length(baselines)
        fig, ax = plt.subplots(num=string(figtitle, "Differential phase data"),
                               nrows=nbase, sharex=true, figsize=figsize, facecolor="White")
        suptitle("Differential phase data")
        subplots_adjust(hspace=0.0)
        mx = matplotlib.ticker.MultipleLocator(20)
        if nbase == 1
            ax = [ax]
        end
        for i in 1:nbase
            plt.axes(ax[i])
            ax[i].set_title(baselines[i], x=0.9, y=0.75)
            loc = [findall(baseline_list_vis[n] .== baselines[i]) for n in eachindex(dvec)]
            wavcol  = vcat([dvec[n].vis_lam[loc[n]]*1e9 for n in eachindex(dvec)]...)
            visphi  = vcat([dvec[n].visphi[loc[n]]       for n in eachindex(dvec)]...)
            vp_err  = vcat([dvec[n].visphi_err[loc[n]]   for n in eachindex(dvec)]...)
            errorbar(wavcol, visphi, yerr=vp_err, fmt="o", markersize=0.5,
                     ecolor="Gainsboro", elinewidth=0.5)
            ax[i].xaxis.set_minor_locator(mx)
            if i == nbase
                xlabel("Wavelength (nm)")
            end
            ylabel("Δφ (°)")
            ax[i].grid(true, which="both", color="Grey", linestyle=":")
        end
        ax[nbase].tick_params(axis="x", which="major", length=10.0)
        ax[nbase].tick_params(axis="x", which="minor", length=5.0)
        tight_layout()
        show(block=false)
    else
        # ── Absolute phase layout: delegate to generic engine ──
        plot_obs(data, OBS_PLOT_SPECS["VISPHI"]; color=color, figsize=figsize,
                 figtitle=figtitle, markopt=markopt, legend_below=legend_below)
    end
end

"""
    plot_visamp(data; kwargs...)

Plot visibility amplitudes vs baseline length (Mλ). Accepts a single `OIdata`,
vector, or 2-D array. The title and y-label adapt to the `AMPTYP` header stored
in `data.amptyp` (`"absolute"`, `"differential"`, `"correlated flux"`, or generic).

# Keyword arguments
- `color` — `"baseline"` (default), `"wav"`, `"mjd"`, or explicit colour string.
- `logplot` — use logarithmic y-axis. Default: `false`.
- `figsize` — figure size. Default: `(12, 6)`.
- `figtitle` — prefix for the figure window title.
- `markopt` — use distinct marker shapes. Default: `false`.
- `legend_below` — place legend below the plot. Default: `false`.
"""
function plot_visamp(data::Union{OIdata, AbstractArray{<:OIdata}}; figsize=oiplot_figsize,
                     logplot=false, color::String="baseline", markopt=false,
                     legend_below=false, figtitle="")
    dvec = as_datavec(data)
    amptyp = dvec[1].amptyp
    if amptyp == "correlated flux"
        amp_title  = "Correlated flux"
        amp_ylabel = "Correlated flux"
    elseif amptyp == "differential"
        amp_title  = "Diff. visamp"
        amp_ylabel = "Diff. visamp"
    else
        amp_title  = "Visamp"
        amp_ylabel = "Visamp"
    end
    plot_obs(data, OBS_PLOT_SPECS["VISAMP"]; color=color, logplot=logplot,
             figsize=figsize, figtitle=figtitle, markopt=markopt,
             legend_below=legend_below, plot_title=amp_title, ylabel_str=amp_ylabel)
end

"""
    plot_diffphi(data; kwargs...)

Deprecated — use `plot_visphi` instead, which auto-detects differential phases
from `data.phityp`.
"""
function plot_diffphi(data::Union{OIdata, AbstractArray{<:OIdata}}; kwargs...)
    @warn("plot_diffphi is deprecated, use plot_visphi instead (auto-detects differential phases)")
    plot_visphi(data; kwargs...)
end


function plot_v2_multifile(data::AbstractVector{<:OIdata}; logplot = false, remove = false,idpoint=false,clean=false,filename="",legend_below=true)
    global v2base=[]
    global v2value=[]
    global v2err=[]
    global clickbase=Array{Int64,2}
    global clickname=[]
    global clickmjd=[]
    global clicklam=[]
    global clickdlam=[]
    global clickfile=[]
    axiscount=0
    testaxis=0
    fig = figure("V2 data",figsize=(10,5),facecolor="White");
    if clean == true
        clf();
    end
    ax = gca();
    if logplot==true
        ax.set_yscale("log")
    end
    #if remove == true # No support for removing points across multimple nights just yet, although this could be VERY useful....
    #    fig.canvas.mpl_connect("button_press_event",onclickv2)
    #end
    for i=1:length(data) #plot each night
        nv2=data[i].nv2
        sta_name=data[i].sta_name
        v2_sta_index=data[i].v2_sta_index
        baseline_v2=data[i].v2_baseline
        v2_data=data[i].v2
        v2_data_err=data[i].v2_err
        if idpoint == true
            v2base=vcat(v2base,baseline_v2)
            v2value=vcat(v2value,v2_data)
            v2err=vcat(v2err,v2_data_err)
            if i==1
                clickbase=cat(data[i].v2_sta_index,dims=2)
            else
                clickbase=cat(clickbase,data[i].v2_sta_index,dims=2)
            end
            clickname=vcat(clickname,data[i].sta_name)
            clickmjd=vcat(clickmjd,data[i].v2_mjd)
            clicklam=vcat(clicklam,data[i].v2_lam)
            clickdlam=vcat(clickdlam,data[i].v2_dlam)
            basearray=Array{String,1}(undef,data[i].nv2)
            basearray[1:data[i].nv2].=data[i].filename
            if i==1
                clickfile=cat(basearray,dims=1)
            else
                clickfile=cat(clickfile,basearray,dims=1)
            end
        end
        baseline_list=get_baseline_names(sta_name,v2_sta_index)
        for j=1:length(unique(baseline_list))
            baseline=unique(baseline_list)[j]
            loc=findall(baseline_list->baseline_list==baseline,baseline_list)
            errorbar(baseline_v2[loc]/1e6,v2_data[loc],yerr=v2_data_err[loc],fmt="o",marker=oiplot_markers[j], markeredgecolor=oiplot_colors[i],color=oiplot_colors[i],markersize=3,ecolor="Gainsboro",elinewidth=1.0,label=baseline)
            if axiscount==0
                if (length(unique(baseline_list)))==15
                    ax.legend(fontsize=8, fancybox=true, shadow=true, ncol=8,loc="upper center", bbox_to_anchor=(0.5, -0.10));
                    testaxis=1
                end
            end
        end
        if testaxis ==1
            axiscount=1
        end
    end
    title("Squared Visibility Amplitude Data")
    xlabel(L"Baseline (M$\lambda$)")
    ylabel("Squared Visibility Amplitudes")
    ax.grid(true,which="both",color="Grey",linestyle=":")
    tight_layout()
    
    #if idpoint==true
    #    cid=fig.canvas.mpl_connect("button_press_event",onclickidentify)
    #end
    if filename !=""
        savefig(filename)
    end
    show(block=false)
end



    



"""
    imdisp(image; kwargs...)

Display a 2-D reconstructed image. The image is normalised to its maximum and
oriented with East left / North up (Monnier convention). Accepts a flat vector
(square image assumed) or a 2-D matrix.

# Keyword arguments
- `pixsize` — pixel scale in mas. When `-1` (default) the axes show pixel indices.
- `colormap` — matplotlib colormap. Default: `"gist_heat"`.
- `figtitle` — figure window title. Default: `"OITOOLS image"`.
- `tickinterval` — minor-tick spacing in mas. Default: `0.5`; auto-scaled for large images.
- `use_colorbar` — show a right-side colorbar. Default: `false`.
- `beamsize` — if `> 0`, draw a filled white circle of this diameter (mas) to indicate the PSF.
- `beamlocation` — `[fx, fy]` fractional position of the beam circle. Default: `[0.8, 0.8]`.
"""
function imdisp(image; figtitle="OITOOLS image", colormap = "gist_heat", pixsize = -1.0, tickinterval = 0.5, use_colorbar = false, beamsize = -1, beamlocation = [])
    fig = figure(figtitle,figsize=(6,6),facecolor="White")
    clf();
    nx=ny=-1;
    pixmode = false;
    if pixsize == -1
        pixmode = true;
        pixsize = 1
    end
    scaling_factor = maximum(image);
    if abs.(scaling_factor) <  1e-20
        scaling_factor = 1.0;
        @warn("Maximum of image < tol");
    end

    img = []
    if ndims(image) ==1
        ny=nx=Int64(sqrt(length(image)))
        img = imshow(rotl90(reshape(image,nx,nx))/scaling_factor, ColorMap(colormap), interpolation="none", extent=[0.5*nx*pixsize,-0.5*nx*pixsize, -0.5*ny*pixsize,0.5*ny*pixsize]); # uses Monnier's orientation
    else
        nx,ny = size(image);
        img = imshow(rotl90(image)/scaling_factor, ColorMap(colormap), interpolation="none", extent=[0.5*nx*pixsize,-0.5*nx*pixsize,-0.5*ny*pixsize,0.5*ny*pixsize]); # uses Monnier's orientation
    end
    if pixmode == false
        xlabel("x ← E (mas)")
        ylabel("y → N (mas)")
    end

    ax = gca()
    ax.set_aspect("equal")
    if (size(image,1)*pixsize)>1000
        tickinterval = 50.0
    elseif (size(image,1)*pixsize)>100
        tickinterval = 5.0
        # ax.invert_xaxis();
    end
    mx = matplotlib.ticker.MultipleLocator(tickinterval) # Define interval of minor ticks
    ax.xaxis.set_minor_locator(mx) # Set interval of minor ticks
    ax.yaxis.set_minor_locator(mx) # Set interval of minor ticks
    ax.xaxis.set_tick_params(which="major",length=10,width=2)
    ax.xaxis.set_tick_params(which="minor",length=5,width=1)
    ax.yaxis.set_tick_params(which="major",length=10,width=2)
    ax.yaxis.set_tick_params(which="minor",length=5,width=1)
    if use_colorbar == true
        divider = pyimport("mpl_toolkits.axes_grid1").make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.2)
        colorbar(img, cax=cax, ax=ax)
    end

    if beamsize > 0
        if beamlocation == []
            beamlocation = [.8, .8]
        end
        c = matplotlib.patches.Circle((0.5*nx*pixsize*beamlocation[1],-0.5*ny*pixsize*beamlocation[2]),0.5*beamsize,fc="white",ec="white",linewidth=0)
        ax.add_artist(c)
    end
    tight_layout()
    show(block=false)
end

"""
    imdisp_multi(cube; kwargs...)

Display a 3-D stack of images (e.g. polychromatic channels or multi-temporal epochs)
in a grid of subplots, one panel per slice.  Each slice is normalised independently
and rendered with the same formatting as [`imdisp`](@ref).

`cube` is an `(nx, nx, nslices)` array.

# Keyword arguments
- `labels`        — vector of strings for subplot titles. Default: `"1"`, `"2"`, ….
- `pixsize`       — pixel scale in mas (`-1` = pixel indices).
- `colormap`      — matplotlib colormap name. Default: `"gist_heat"`.
- `figtitle`      — figure window title. Default: `"OITOOLS images"`.
- `tickinterval`  — minor-tick spacing (auto-adjusted for large fields of view).
- `use_colorbar`  — show a colour bar beside each panel.
- `beamsize`      — if `> 0`, draw a beam circle on each panel.
- `beamlocation`  — `[x, y]` fractional position for the beam circle.
"""
function imdisp_multi(cube::Array{Float64,3};
        labels::Vector{String} = String[],
        figtitle::String       = "OITOOLS images",
        colormap::String       = "gist_heat",
        pixsize::Float64       = -1.0,
        tickinterval::Float64  = 0.5,
        use_colorbar::Bool     = false,
        beamsize::Float64      = -1.0,
        beamlocation::Vector{Float64} = Float64[])

    nx, ny, nslices = size(cube)
    nside = ceil(Int64, sqrt(nslices))

    pixmode = pixsize == -1
    pix = pixmode ? 1.0 : pixsize

    # Auto-adjust tick interval to field of view
    fov = nx * pix
    if fov > 1000
        tickinterval = 50.0
    elseif fov > 100
        tickinterval = 5.0
    end

    fig = figure(figtitle, figsize=(4*nside, 4*nside), facecolor="White")
    clf()

    for i in 1:nslices
        fig.add_subplot(nside, nside, i)

        lab = isempty(labels) ? "$i" : labels[i]
        title(lab)

        panel = cube[:, :, i]
        scaling = maximum(panel)
        if abs(scaling) < 1e-20
            scaling = 1.0
        end

        img = imshow(rotl90(panel) / scaling, ColorMap(colormap),
                      interpolation="none",
                      extent=[0.5*nx*pix, -0.5*nx*pix, -0.5*ny*pix, 0.5*ny*pix])

        if !pixmode
            xlabel("x ← E (mas)")
            ylabel("y → N (mas)")
        end

        ax = gca()
        ax.set_aspect("equal")
        mx = matplotlib.ticker.MultipleLocator(tickinterval)
        ax.xaxis.set_minor_locator(mx)
        ax.yaxis.set_minor_locator(mx)
        ax.xaxis.set_tick_params(which="major", length=10, width=2)
        ax.xaxis.set_tick_params(which="minor", length=5, width=1)
        ax.yaxis.set_tick_params(which="major", length=10, width=2)
        ax.yaxis.set_tick_params(which="minor", length=5, width=1)

        if use_colorbar
            divider = pyimport("mpl_toolkits.axes_grid1").make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.2)
            colorbar(img, cax=cax, ax=ax)
        end

        if beamsize > 0
            bloc = isempty(beamlocation) ? [0.8, 0.8] : beamlocation
            c = matplotlib.patches.Circle(
                (0.5*nx*pix*bloc[1], -0.5*ny*pix*bloc[2]),
                0.5*beamsize, fc="white", ec="white", linewidth=0)
            ax.add_artist(c)
        end
    end
    tight_layout()
    show(block=false)
end


function get_station_names(sta_names, sta_indx::AbstractVector)
    return [si == 0 ? "Calibrated" : sta_names[si] for si in sta_indx]
end

function get_baseline_names(sta_names,sta_indx)
    nbaselines = size(sta_indx,2)
    baseline_names=Array{String}(undef,nbaselines)
    for i=1:nbaselines
        baseline_names[i]=string(sta_names[sta_indx[1,i]],"-",sta_names[sta_indx[2,i]])
    end
    return baseline_names
end

function get_triplet_names(sta_names, sta_indx)
    nt3=size(sta_indx,2)
    triplet_names=Array{String}(undef,nt3)
    for i=1:nt3
        triplet_names[i]=string(sta_names[sta_indx[1,i]],"-",sta_names[sta_indx[2,i]],"-",sta_names[sta_indx[3,i]])
    end
    return triplet_names
end

function get_triple_baselines_names(sta_names, sta_indx)
    nt3=size(sta_indx,2)
    baseline_names=Array{String}(undef, 3, nt3)
    for i=1:nt3
        baseline_names[1, i]=string(sta_names[sta_indx[1,i]],"-",sta_names[sta_indx[2,i]])
        baseline_names[2, i]=string(sta_names[sta_indx[2,i]],"-",sta_names[sta_indx[3,i]])
        baseline_names[3, i]=string(sta_names[sta_indx[1,i]],"-",sta_names[sta_indx[3,i]])
    end
    return baseline_names
end


function plot_facility(facility)
    coords = facility.sta_xyz'  # (3, ntel)
    scatter(coords[1,:], coords[2,:])
    axis("equal")
    for i=1:facility.ntel
        text(coords[1,i]+1, coords[2,i]+1, facility.sta_names[i])
    end
end
    
function imshow2(img1, img2)
    fig,ax = subplots(nrows=1, ncols=2)
    ax[1].imshow(img1)
    ax[2].imshow(img2)
end