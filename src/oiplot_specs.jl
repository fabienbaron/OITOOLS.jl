# oiplot_specs.jl — what to plot, described without saying how to draw it.
#
# These definitions were part of oiplot.jl until plotting moved into
# ext/OITOOLSPythonPlotExt.jl. They stayed behind in the core package on purpose: none of
# them touches matplotlib. They are the enumeration of OITOOLS' observables (which field,
# which errors, which x axis, how points are grouped and scaled) plus the categorical
# palette, and every plotting front-end needs them -- the matplotlib extension, the Makie
# layer in OITOOLSGUI, and any future one.
#
# Keeping them here is what lets a GUI drive its plot menu from OBS_PLOT_SPECS and colour
# its baselines identically to oiplot without loading Python. Duplicating the table in a
# front-end instead is how the two drift apart.

using LaTeXStrings   # the axis labels carry their own TeX

# Saturated categorical palette. tab20 was tried instead, but half of it is pastel "light"
# variants that wash out at the small marker sizes used here — with many baselines the plot
# stopped being readable. These are the strongly-saturated named colours used historically.
#
# 28 entries, indexed with mod1, so a dataset with more baselines than colours repeats rather
# than erroring. Note the historical list named "blue" twice, which silently gave two
# baselines the same colour; the second is "navy" here.
const global oiplot_colors = [
    "black", "gold", "chartreuse", "blue", "red", "pink", "lightgray", "darkorange",
    "darkgreen", "aqua", "fuchsia", "saddlebrown", "dimgray", "darkslateblue", "violet",
    "indigo", "navy", "dodgerblue", "sienna", "olive", "purple", "darkorchid", "tomato",
    "darkturquoise", "steelblue", "seagreen", "darkgoldenrod", "darkseagreen",
]

const global oiplot_markers=["o","s","v","P","*","x","^","D","p",1,"<","H","X","4",4,"_","1",6,"8","d",9]

# ── Data normalization ────────────────────────────────────────────────────────
"""
    as_datavec(data) -> Vector{<:OIdata}

Ensure `data` is a flat `Vector{<:OIdata}`, wrapping a single `OIdata` or
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
