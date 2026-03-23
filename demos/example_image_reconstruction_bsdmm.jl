# BSDMM image reconstruction example
#
# Compares BSDMM (proximal TV + centering) with the standard VMLMB
# reconstruction on simulated interferometric data.

using OITOOLS
using OptimPackNextGen
using LinearAlgebra, Printf

oifitsfile = joinpath(@__DIR__, "data", "2004-data1.oifits")
data = readoifits(oifitsfile)

nx = 64
pixsize = 0.2
ft = setup_ft(data, nx, pixsize)

x_start = gaussian2d(nx, nx, nx / 6)

println("=== Starting chi2 ===")
image_to_chi2(x_start, ft, data; verb=true)

# ── BSDMM reconstruction ──────────────────────────────────────────────────
println("\n=== BSDMM: TV + centering ===")
x_bsdmm = reconstruct_bsdmm(x_start, data, ft;
                              mu_reg=1e5, mu_cen=1e2,
                              reg_type=:tv, tv_niter=50,
                              rho_init=1e4, adaptive=false,
                              maxit=500, x_maxiter=5)
println("\nBSDMM final:")
image_to_chi2(x_bsdmm, ft, data; verb=true)

# ── Standard reconstruct (baseline) ───────────────────────────────────────
println("\n=== Standard reconstruct (baseline) ===")
x_ref = reconstruct(copy(x_start), data, ft;
                     regularizers=[["centering", 1e2], ["tv", 1e5]],
                     verb=true, maxiter=500)
println("\nStandard final:")
image_to_chi2(x_ref, ft, data; verb=true)

# ── Side-by-side comparison ───────────────────────────────────────────────
using PyPlot
cube = cat(x_bsdmm, x_ref; dims=3)
imdisp_multi(cube; labels=["BSDMM", "Standard"], pixsize=pixsize,
             figtitle="BSDMM vs Standard reconstruction")
savefig("bsdmm_comparison.png", dpi=150, bbox_inches="tight")
println("Saved bsdmm_comparison.png")
