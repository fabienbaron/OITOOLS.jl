##############################################################################
#  show_pmoired_models.jl
#
#  Replicates the visual plots from PMOIRED's
#  "Model definitions and examples.ipynb" notebook.
#
#  Covers:
#    1. Basic models (ud, profile disk, ring, gaussian, crescent, gaussian ring)
#    2. Geometric transformations (x/y offset, incl/projang)
#    3. Spatial kernel (crescent with and without kernel)
#    4. Combining components (star + disk)
#    5. Parameter expressions (shared PA/INC)
#
##############################################################################

using OITOOLS
using PythonPlot, Printf

# ─────────────────────────────────────────────────────────────────────────────
# Helper: show model image with parameter overlay
# ─────────────────────────────────────────────────────────────────────────────

"""
    show_model(model; fov=1.5, nx=128, title="", ax=nothing, kwargs...)

Synthesise and display a model image, similar to PMOIRED's `showModel`.
This is our main debugging tool to check if we're compatible!
`fov` is the full field-of-view in mas (axes go from -fov/2 to fov/2). `nx` is the grid size (default 128, matching PMOIRED).
"""
function show_model(model_dict::Dict{String,Any};
                    fov::Real     = 3.0,
                    nx::Int       = 128,
                    title_str::String = "",
                    ax = nothing,
                    colormap::String = "inferno",
                    power::Real = 1.0,
                    wl = nothing,
                    show_legend::Bool = true)

    half_fov = fov / 2        # fov is full width (PMOIRED convention)
    pixsize = fov / nx        # derive pixel size from FOV and grid size

    # Separate list_free_params = [] since we just want to display
    model = dict_to_model(model_dict, String[])
    img   = model_to_image(model, Float64[]; nx, pixsize, wl)

    # Apply power scaling for display (like PMOIRED's imPow)
    img_disp = img .^ power

    extent = [half_fov, -half_fov, -half_fov, half_fov]   # E left, W right (RA convention)

    if isnothing(ax)
        ax = gca()
    end

    ax.imshow(rotl90(img_disp), extent=extent, cmap=colormap,
              origin="lower", interpolation="bilinear")
    ax.set_xlabel("E ← ΔRA (mas) → W")
    ax.set_ylabel("S ← ΔDec (mas) → N")
    if !isempty(title_str)
        ax.set_title(title_str)
    end

    # Overlay parameter text
    if show_legend
        lines = String[]
        for k in sort(collect(keys(model_dict)))
            v = model_dict[k]
            if v isa Number
                push!(lines, @sprintf("%-20s = %.3g", k, v))
            elseif v isa String
                v_safe = replace(v, raw"$" => raw"\$")  # escape $ for matplotlib
                push!(lines, @sprintf("%-20s = %s", k, v_safe))
            end
        end
        txt = join(lines, "\n")
        ax.text(0.98, 0.02, txt, transform=ax.transAxes,
                fontsize=6, color="white", family="monospace",
                verticalalignment="bottom", horizontalalignment="right",
                bbox=Dict("facecolor"=>"black", "alpha"=>0.5, "pad"=>3))
    end

    return img
end


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Basic models
# ═══════════════════════════════════════════════════════════════════════════════

models_basic = [
    (Dict{String,Any}("star,ud" => 1.0, "star,f" => 1.0),
     "Uniform disk"),
    (Dict{String,Any}("star,diam" => 1.0, "star,profile" => "\$MU^0.5", "star,f" => 1.0),
     "Disk with profile (\$MU^0.5)"),
    (Dict{String,Any}("ring,udout" => 1.0,
                       "ring,profile" => "(\$R > 0.25) * \$R^(-0.5)",
                       "ring,f" => 1.0),
     "Ring with profile (\$R^-0.5)"),
    (Dict{String,Any}("g,fwhm" => 0.5, "g,f" => 1.0),
     "Gaussian"),
    (Dict{String,Any}("cr,crin" => 0.8, "cr,crout" => 1.0,
                       "cr,croff" => 0.8, "cr,crprojang" => 120.0,
                       "cr,f" => 1.0),
     "Crescent"),
    (Dict{String,Any}("gr,fwhmin" => 0.2, "gr,fwhmout" => 0.5, "gr,f" => 1.0),
     "Gaussian ring"),
]

fig, axes = subplots(2, 3; figsize=(14, 9))
fig.suptitle("1. Basic models", fontsize=14, fontweight="bold")
for (i, (model_dict, name)) in enumerate(models_basic)
    ax = axes[i]
    show_model(model_dict; fov=1.5, ax, title_str=name)
end
tight_layout()
savefig("pmoired_1_basic_models.png", dpi=150)
println("Saved pmoired_1_basic_models.png")


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Simple geometric transformations
# ═══════════════════════════════════════════════════════════════════════════════

models_transform = [
    (Dict{String,Any}("star,ud" => 1.0, "star,f" => 1.0),
     "UD (reference)"),
    (Dict{String,Any}("star,ud" => 0.5, "star,x" => 0.5, "star,y" => -0.2,
                       "star,f" => 1.0),
     "UD offset (x=0.5, y=-0.2)"),
    (Dict{String,Any}("star,ud" => 1.0, "star,incl" => 60.0, "star,projang" => 30.0,
                       "star,f" => 1.0),
     "UD incl=60, PA=30"),
    (Dict{String,Any}("ring,diamin" => 0.5, "ring,diamout" => 1.0,
                       "ring,incl" => 45.0, "ring,projang" => 60.0,
                       "ring,f" => 1.0),
     "Ring incl=45, PA=60"),
    (Dict{String,Any}("cr,crin" => 0.5, "cr,crout" => 1.0,
                       "cr,croff" => 0.8, "cr,crprojang" => 120.0,
                       "cr,incl" => 45.0, "cr,projang" => 30.0,
                       "cr,f" => 1.0),
     "Crescent incl=45, PA=30"),
    (Dict{String,Any}("g,fwhm" => 0.5, "g,incl" => 60.0, "g,projang" => 0.0,
                       "g,f" => 1.0),
     "Gaussian incl=60, PA=0"),
]

fig, axes = subplots(2, 3; figsize=(14, 9))
fig.suptitle("2. Geometric transformations", fontsize=14, fontweight="bold")
for (i, (model_dict, name)) in enumerate(models_transform)
    ax = axes[i]
    show_model(model_dict; fov=1.5, ax, title_str=name)
end
tight_layout()
savefig("pmoired_2_transforms.png", dpi=150)
println("Saved pmoired_2_transforms.png")


# ═══════════════════════════════════════════════════════════════════════════════
# 2b. Azimuthal variations (PMOIRED notebook Cell 7 & 8)
# ═══════════════════════════════════════════════════════════════════════════════

Is = [1, 2, 3]   # harmonic orders

# Row 1: no az, then individual harmonics 1, 2, 3
# Row 2: all harmonics combined, plus ring+profile with az amp1
fig, axes = subplots(2, 3; figsize=(14, 9))
fig.suptitle("2b. Azimuthal variations", fontsize=14, fontweight="bold")

# Row 1: no az + individual harmonics 1, 2, 3
for (col, i) in enumerate([0; Is])
    model_dict = Dict{String,Any}(
        "disk,diamin" => 1.0, "disk,diamout" => 3.0,
        "disk,profile" => "1",
        "disk,projang" => -20.0, "disk,incl" => 60.0, "disk,f" => 1.0,
    )
    if i > 0
        model_dict["disk,az amp$i"]     = 1.0 / length(Is)
        model_dict["disk,az projang$i"] = Float64(i * 30)
    end
    row = col <= 3 ? 1 : 2
    c   = col <= 3 ? col : col - 3
    show_model(model_dict; fov=3.5, ax=axes[row, c],
               title_str= i == 0 ? "No variation" : "az order $i")
end

# Row 2, col 2: all harmonics combined
model_dict_all = Dict{String,Any}(
    "disk,diamin" => 1.0, "disk,diamout" => 3.0,
    "disk,profile" => "1",
    "disk,projang" => -20.0, "disk,incl" => 60.0, "disk,f" => 1.0,
)
for i in Is
    model_dict_all["disk,az amp$i"]     = 1.0 / length(Is)
    model_dict_all["disk,az projang$i"] = Float64(i * 30)
end
show_model(model_dict_all; fov=3.5, ax=axes[2, 2], title_str="All harmonics combined")

# Row 2, col 3: ring + profile + az (Cell 8)
model_dict_az_prof = Dict{String,Any}(
    "disk,diamin" => 0.5, "disk,diamout" => 1.0,
    "disk,profile" => "\$R^(-2)",
    "disk,az amp1" => 0.8, "disk,az projang1" => 60.0,
    "disk,projang" => 45.0, "disk,incl" => -30.0,
    "disk,x" => -0.2, "disk,y" => 0.1, "disk,f" => 1.0,
)
show_model(model_dict_az_prof; fov=1.5, ax=axes[2, 3], title_str="Ring + profile + az amp1")

tight_layout()
savefig("pmoired_2b_azimuthal.png", dpi=150)
println("Saved pmoired_2b_azimuthal.png")


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Spatial kernel
# ═══════════════════════════════════════════════════════════════════════════════

model_dict_nokernel = Dict{String,Any}(
    "cr,crin" => 0.5, "cr,crout" => 1.0,
    "cr,croff" => 0.8, "cr,crprojang" => 120.0,
    "cr,incl" => 45.0, "cr,projang" => 30.0,
    "cr,f" => 1.0)

model_dict_kernel = copy(model_dict_nokernel)
model_dict_kernel["spatial_kernel"] = 0.2

fig, axes = subplots(1, 2; figsize=(10, 5))
fig.suptitle("3. Spatial kernel (Gaussian smoothing)", fontsize=14, fontweight="bold")
show_model(model_dict_nokernel; fov=1.5, ax=axes[1], title_str="Without kernel")
show_model(model_dict_kernel;   fov=1.5, ax=axes[2], title_str="With spatial_kernel=0.2")
tight_layout()
savefig("pmoired_3_spatial_kernel.png", dpi=150)
println("Saved pmoired_3_spatial_kernel.png")


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Combining components
# ═══════════════════════════════════════════════════════════════════════════════

# Star + inclined disk with profile and az variation (PMOIRED Cell 12)
model_dict_multi = Dict{String,Any}(
    "star,fwhm"          => 0.1,
    "star,spectrum"       => "\$WL^(-3)",
    "disk,diamin"        => 0.5,
    "disk,diamout"       => 1.0,
    "disk,profile"       => "\$R^(-2)",
    "disk,az amp1"       => 1.0,
    "disk,az projang1"   => 60.0,
    "disk,projang"       => 45.0,
    "disk,incl"          => -30.0,
    "disk,x"             => -0.05,
    "disk,y"             => 0.05,
    "disk,spectrum"       => "5 * \$WL^(-2)",
)

# Images at min/max WL + SED (like PMOIRED with WL = linspace(1,2,10))
WL_range = range(1.0, 2.0; length=10)

fig = figure(figsize=(14, 5))
fig.suptitle("4. Combining components", fontsize=14, fontweight="bold")

ax1 = fig.add_subplot(1, 3, 1)
show_model(model_dict_multi; fov=1.5, ax=ax1,
           title_str=@sprintf("λ = %.1f µm", first(WL_range)),
           power=0.5, wl=first(WL_range))

ax2 = fig.add_subplot(1, 3, 2)
show_model(model_dict_multi; fov=1.5, ax=ax2,
           title_str=@sprintf("λ = %.1f µm", last(WL_range)),
           power=0.5, wl=last(WL_range))

# SED panel: computed from model via model_to_sed
ax3 = fig.add_subplot(1, 3, 3)
model_multi = dict_to_model(model_dict_multi, String[])
wl_sed = collect(range(first(WL_range), last(WL_range); length=200))
f_total, f_comps = model_to_sed(model_multi, Float64[], wl_sed)
ax3.plot(wl_sed, f_total, label="total", linewidth=2, color="black")
for (name, f) in sort(collect(f_comps))
    ax3.plot(wl_sed, f, label=name, linestyle="--")
end
ax3.set_xlabel("λ (µm)")
ax3.set_ylabel("Flux (a.u.)")
ax3.set_title("SED")
ax3.legend(fontsize=8)
ax3.set_xlim(first(WL_range), last(WL_range))

tight_layout()
savefig("pmoired_4_combining.png", dpi=150)
println("Saved pmoired_4_combining.png")


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Parameter expressions (shared PA, INC) — PMOIRED Cell 14
# ═══════════════════════════════════════════════════════════════════════════════

WL_expr = range(1.0, 2.0; length=10)

model_dict_expr = Dict{String,Any}(
    "PA"                  => 60.0,
    "INC"                 => 45.0,
    "inner,fwhm"          => 1.0,
    "inner,projang"       => "\$PA",
    "inner,incl"          => "\$INC",
    "inner,spectrum"      => "\$WL^(-3)",
    "outer,diamin"        => "2 * \$inner,fwhm",
    "outer,diamout"       => "3 * \$outer,diamin",
    "outer,projang"       => "\$PA",
    "outer,incl"          => "\$INC",
    "outer,spectrum"      => "\$F0 * \$WL^1",
    "F0"                  => 0.5,
    "outer,az amp1"       => 0.9,
    "outer,az projang1"   => 90.0,
    "outer,profile"       => "\$R^\$PL",
    "PL"                  => -2.0,
)

fig = figure(figsize=(14, 5))
fig.suptitle("5. Parameter expressions", fontsize=14, fontweight="bold")

ax1 = fig.add_subplot(1, 3, 1)
show_model(model_dict_expr; fov=9.0, ax=ax1, power=0.5,
           title_str="λ = $(first(WL_expr)) µm", wl=first(WL_expr))

ax2 = fig.add_subplot(1, 3, 2)
show_model(model_dict_expr; fov=9.0, ax=ax2, power=0.5,
           title_str="λ = $(last(WL_expr)) µm", wl=last(WL_expr))

# SED
ax3 = fig.add_subplot(1, 3, 3)
model_expr = dict_to_model(model_dict_expr, String[])
wl_sed_e = collect(range(first(WL_expr), last(WL_expr); length=200))
f_total_e, f_comps_e = model_to_sed(model_expr, Float64[], wl_sed_e)
ax3.plot(wl_sed_e, f_total_e, label="total", linewidth=2, color="black")
for (name, f) in sort(collect(f_comps_e))
    ax3.plot(wl_sed_e, f, label=name, linestyle="--")
end
ax3.set_xlabel("λ (µm)"); ax3.set_ylabel("Flux (a.u.)")
ax3.set_title("SED"); ax3.legend(fontsize=8)
ax3.set_xlim(first(WL_expr), last(WL_expr))

tight_layout()
savefig("pmoired_5_expressions.png", dpi=150)
println("Saved pmoired_5_expressions.png")


# ═══════════════════════════════════════════════════════════════════════════════
# 5b. Spiral pattern (PMOIRED Cell 16)
# ═══════════════════════════════════════════════════════════════════════════════

# Build N concentric rings whose az projang1 rotates with radius → spiral
N_rings = 5
model_dict_spiral = Dict{String,Any}(
    "Din"     => 1.0,
    "Dout"    => 7.0,
    "INCL"    => 45.0,
    "PROJANG" => 60.0,
    "pitch"   => 60.0,
    "PAin"    => -45.0,
)
for i in 0:(N_rings-1)
    k = "ring$i"
    model_dict_spiral["$k,diamin"]       = "\$Din + $i/$N_rings * (\$Dout - \$Din)"
    model_dict_spiral["$k,diamout"]      = "\$Din + $(i+1)/$N_rings * (\$Dout - \$Din)"
    model_dict_spiral["$k,profile"]      = "1"
    model_dict_spiral["$k,az amp1"]      = 1.0
    model_dict_spiral["$k,az projang1"]  = "\$PAin + \$pitch * \$$k,diamin"
    model_dict_spiral["$k,incl"]         = "\$INCL"
    model_dict_spiral["$k,projang"]      = "\$PROJANG"
end

fig, axes = subplots(1, 2; figsize=(10, 5))
fig.suptitle("5b. Spiral from azimuthal variations", fontsize=14, fontweight="bold")

show_model(model_dict_spiral; fov=8.0, ax=axes[1], title_str="Without kernel", show_legend=false)

# Add spatial kernel for smoothing
model_dict_spiral_k = copy(model_dict_spiral)
model_dict_spiral_k["spatial_kernel"] = "0.5 * (\$Dout - \$Din) / $N_rings"
show_model(model_dict_spiral_k; fov=8.0, ax=axes[2], title_str="With spatial kernel", show_legend=false)

tight_layout()
savefig("pmoired_5b_spiral.png", dpi=150)
println("Saved pmoired_5b_spiral.png")


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Gallery: all shapes side by side at different inclinations
# ═══════════════════════════════════════════════════════════════════════════════

shapes = [
    ("UD",             Dict{String,Any}("s,ud"=>1.0, "s,f"=>1.0)),
    ("Gaussian",       Dict{String,Any}("s,fwhm"=>0.6, "s,f"=>1.0)),
    ("Ring",           Dict{String,Any}("s,diamin"=>0.6, "s,diamout"=>1.0, "s,f"=>1.0)),
    ("Gaussian ring",  Dict{String,Any}("s,fwhmin"=>0.3, "s,fwhmout"=>0.6, "s,f"=>1.0)),
    ("Crescent",       Dict{String,Any}("s,crin"=>0.6, "s,crout"=>1.0,
                                         "s,croff"=>0.7, "s,crprojang"=>90.0, "s,f"=>1.0)),
    ("LD pow (α=0.5)", Dict{String,Any}("s,ldpow"=>1.0, "s,alpha"=>0.5, "s,f"=>1.0)),
]
incls = [0.0, 30.0, 60.0]

fig, axes = subplots(length(shapes), length(incls); figsize=(4*length(incls), 4*length(shapes)))
fig.suptitle("6. All shapes at various inclinations (PA=45)", fontsize=14, fontweight="bold")

for (j, incl) in enumerate(incls)
    for (i, (name, base_model_dict)) in enumerate(shapes)
        model_dict = copy(base_model_dict)
        if incl > 0
            model_dict["s,incl"]    = incl
            model_dict["s,projang"] = 45.0
        end
        ax = axes[i, j]
        label = j == 1 ? "$name, i=$(Int(incl))" : "i=$(Int(incl))"
        show_model(model_dict; fov=1.5, ax, title_str=label)
        if j > 1
            ax.set_ylabel("")
        end
    end
end
tight_layout()
savefig("pmoired_6_gallery.png", dpi=150)
println("Saved pmoired_6_gallery.png")


# ═══════════════════════════════════════════════════════════════════════════════
# 7. diam+thick sugar demonstration
# ═══════════════════════════════════════════════════════════════════════════════

fig, axes = subplots(1, 4; figsize=(16, 4))
fig.suptitle("7. diam + thick ring sugar", fontsize=14, fontweight="bold")
for (i, thick) in enumerate([0.05, 0.3, 0.6, 1.0])
    model_dict = Dict{String,Any}("r,diam" => 2.0, "r,thick" => thick, "r,f" => 1.0)
    ax = axes[i]
    show_model(model_dict; fov=3.0, ax, title_str=@sprintf("thick=%.2f", thick))
end
tight_layout()
savefig("pmoired_7_thick.png", dpi=150)
println("Saved pmoired_7_thick.png")


println("\nAll plots saved.")
