# example_simulate_polychromatic_disk.jl
#
# Simulate CHARA/MIRC-X observations of a chromatic, time-variable disk
# with an off-axis companion that produces non-zero differential phases.
#
# The model has three components:
#   - star:      unresolved central source (V = 1)
#   - disk:      uniform disk whose diameter grows with wavelength
#                (e.g. optically thin dust; diameter ∝ λ^α)
#   - companion: faint off-axis point source at (dx, dy) mas
#                Its phase signature is chromatic because the fringe
#                spacing changes with λ, producing a differential phase
#                signal across the H-band.
#
# The simulation writes OI_VIS2, OI_VIS (with differential phases),
# OI_T3, and OI_FLUX tables via the OIFITS.jl library.

using OITOOLS
using Dates
using OIFITS   # for the read-back verification at the end

# ═══════════════════════════════════════════════════════════════════════════════
# 1. Observation setup
# ═══════════════════════════════════════════════════════════════════════════════

# Six hours of CHARA observation, one snapshot every 20 minutes
dates = collect(DateTime(2025, 6, 15, 3, 0, 0):Minute(20):DateTime(2025, 6, 15, 9, 0, 0))

# Read instrument configuration files
facility   = read_facility_file("CHARA")
combiner   = read_comb_file("MIRCX")
wavelength = read_wave_file("MIRCX_LOWH")

# Target: a fictional bright giant in Cygnus (H ≈ 2.0)
target = read_obs_file("default_obs")
target.target  = "HD 999999"
target.raep0   = 15*[20, 30, 0.0]' * [1.0, 1/60.0, 1/3600.0]   # RA  hms -> degrees
target.decep0  =    [40, 0, 0.0]'  * [1.0, 1/60.0, 1/3600.0]   # Dec dms -> degrees

# ═══════════════════════════════════════════════════════════════════════════════
# 2. Model definition — chromatic disk + off-axis companion
# ═══════════════════════════════════════════════════════════════════════════════

# Reference wavelength (H-band centre)
λ0 = 1.65e-6   # metres

# Chromatic flux fractions:
#   star  = 0.50 at λ0
#   disk  = 0.47 at λ0
#   comp  = 0.03 at λ0 (faint companion)

model = Dict{String,Any}(
    # ── Star (unresolved point source) ────────────────────────────────
    "star,f"     => 0.50,

    # ── Disk (uniform disk, chromatic diameter) ───────────────────────
    # diameter grows as λ^1.2 — dust emission extends at longer wavelengths
    "disk,ud"    => "3.0 * (\$WL / $(λ0))^1.2",
    "disk,f"     => 0.47,

    # ── Companion (off-axis point source) ─────────────────────────────
    # 3 mas East, 1 mas North — produces λ-dependent phase signature
    "comp,f"     => "1.0 - \$star,f - \$disk,f",
    "comp,x"     => 3.0,    # ΔRA  offset [mas]
    "comp,y"     => 1.0,    # ΔDEC offset [mas]
)

# No free parameters for simulation
free_params = String[]

# Compile the model
params = dict_to_model(model, free_params)
x = Float64[]

# ═══════════════════════════════════════════════════════════════════════════════
# 3. Simulate and save
# ═══════════════════════════════════════════════════════════════════════════════

out_file = "./data/simulation_chromatic_disk.oifits"

simulate(facility, target, combiner, wavelength, dates, out_file;
         flat_model  = params,
         flat_params = x,
         mag = 2.0)

println("Wrote: $out_file")

# ═══════════════════════════════════════════════════════════════════════════════
# 4. Quick verification: read back via OIFITS.jl and inspect all tables
# ═══════════════════════════════════════════════════════════════════════════════

ds = OIFITS.OIDataSet(out_file)

println("\nOIFITS tables written:")
println("  OI_VIS2  blocks: ", length(ds.vis2))
println("  OI_VIS   blocks: ", length(ds.vis),  "  (PHITYP = ", ds.vis[1].phityp, ")")
println("  OI_T3    blocks: ", length(ds.t3))
println("  OI_FLUX  blocks: ", length(ds.flux), "  (CALSTAT = ", ds.flux[1].calstat, ")")

nw = length(ds.instr[1].eff_wave)
println("  Wavelengths    : ", nw, " channels  (",
        round(ds.instr[1].eff_wave[end]*1e6, digits=3), " – ",
        round(ds.instr[1].eff_wave[1]*1e6, digits=3), " μm)")

# Show differential phase range — should be non-zero thanks to the companion
dphi = ds.vis[1].visphi
println("\n  Differential phase range: ",
        round(minimum(dphi), digits=2), " to ",
        round(maximum(dphi), digits=2), " deg")
println("  (companion at 3 mas E, 1 mas N produces a λ-dependent phase signal)")

# Read back via OITOOLS for plotting
data = readoifits(out_file)[1,1]
println("\n  V² points : ", length(data.v2))
println("  T3φ points: ", length(data.t3phi))
println("  MJD range : ", extrema(data.v2_mjd))

# Plot UV coverage coloured by wavelength
uvplot(data, color="wav")
