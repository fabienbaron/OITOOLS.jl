# test_gradient_sign.jl
#
# Finite-difference gradient validation for the DFT and NFFT chi2 paths.
# Tests V2, T3amp, and T3phi gradient components independently.
#
# Usage:  julia --project=. test/test_gradient_sign.jl
#
# This script is sign-convention-independent: it validates that the analytic
# gradient matches the finite-difference gradient regardless of the DFT kernel
# sign.  Run it before AND after sign convention changes to verify correctness.
#
# Note: chi2_fg computes the gradient WITHOUT the flux-normalization correction
# (that correction lives in crit_fg).  We therefore test crit_fg with empty
# regularizers, which applies the full chain rule through x/sum(x).

using OITOOLS, FITSIO, LinearAlgebra, Printf

# ── Load data ────────────────────────────────────────────────────────────────
datafile = joinpath(@__DIR__, "..", "demos", "data", "2004-data1.oifits")
if !isfile(datafile)
    error("Test data not found: $datafile")
end
data = readoifits(datafile; warn=false)[1,1]

# ── Load ground truth image and resize to small grid for speed ───────────────
truefile = joinpath(@__DIR__, "..", "demos", "data", "2004true.fits")
if isfile(truefile)
    img_full = Float64.(read(FITS(truefile)[1]))
    # Downsample to 32×32 by block averaging
    nx_full = size(img_full, 1)
    nx = 32
    block = nx_full ÷ nx
    x = zeros(nx, nx)
    for j in 1:nx, i in 1:nx
        x[i,j] = sum(img_full[(i-1)*block+1:i*block, (j-1)*block+1:j*block])
    end
    pixsize = 0.101 * (nx_full / nx)  # scale pixsize to match
else
    # Fallback: random positive image
    nx = 32
    x = rand(nx, nx) .+ 0.01
    pixsize = 0.2
end
x .= max.(x, 1e-10)  # ensure strictly positive (needed for crit_fg normalization)
x ./= sum(x)          # normalize

# ── Setup DFT and NFFT ──────────────────────────────────────────────────────
dft = setup_dft(data, nx, pixsize)
ft  = setup_nfft(data, nx, pixsize)

# ── Finite-difference gradient check ────────────────────────────────────────
function fd_gradient(chi2_func, x0; eps=1e-7)
    g_fd = similar(x0)
    f0 = chi2_func(x0)
    for k in eachindex(x0)
        xp = copy(x0); xp[k] += eps
        xm = copy(x0); xm[k] -= eps
        g_fd[k] = (chi2_func(xp) - chi2_func(xm)) / (2eps)
    end
    return g_fd
end

function check_gradient(label, fg_func, f_func, x0; eps=1e-7, tol=5e-4)
    g_analytic = similar(x0)
    f = fg_func(x0, g_analytic)
    g_fd = fd_gradient(f_func, x0; eps)

    # Relative error (avoid division by zero)
    max_g = max(maximum(abs.(g_analytic)), maximum(abs.(g_fd)), 1e-30)
    abs_err = maximum(abs.(g_analytic .- g_fd))
    rel_err = abs_err / max_g

    pass = rel_err < tol
    status = pass ? "PASS" : "FAIL"
    @printf("  %-25s  rel_err=%.2e  abs_err=%.2e  [%s]\n", label, rel_err, abs_err, status)
    if !pass
        idx = argmax(abs.(g_analytic .- g_fd))
        @printf("    worst pixel %d: analytic=%.8e  fd=%.8e\n",
                idx, g_analytic[idx], g_fd[idx])
    end
    return pass
end

# ── Wrapper: crit_fg without regularization (includes normalization fix) ─────
# crit_fg signature: crit_fg(x_2d, g_2d, ft, data; weights, regularizers, verb)
# It calls chi2_fg then applies the normalization correction:
#   g[:] = (g .- sum(x.*g)/flux) / flux

function make_fg(transform, data, nx, weights)
    function fg(x0_vec, g_vec)
        x2 = reshape(x0_vec, nx, nx)
        g2 = zeros(nx, nx)
        f = crit_fg(x2, g2, transform, data; weights=weights, regularizers=[], warn=false)
        g_vec .= vec(g2)
        return f
    end
    function f_only(x0_vec)
        x2 = reshape(x0_vec, nx, nx)
        g2 = zeros(nx, nx)
        return crit_fg(x2, g2, transform, data; weights=weights, regularizers=[], warn=false)
    end
    return fg, f_only
end

# ── Test each observable independently ───────────────────────────────────────
println("\n=== Gradient validation (DFT path) ===")
global all_pass = true

for (label, w) in [("DFT V2",     [1.0, 0.0, 0.0]),
                    ("DFT T3amp",  [0.0, 1.0, 0.0]),
                    ("DFT T3phi",  [0.0, 0.0, 1.0]),
                    ("DFT all",    [1.0, 1.0, 1.0])]
    fg, fo = make_fg(dft, data, nx, w)
    global all_pass &= check_gradient(label, fg, fo, vec(x))
end

println("\n=== Gradient validation (NFFT path) ===")

for (label, w) in [("NFFT V2",     [1.0, 0.0, 0.0]),
                    ("NFFT T3amp",  [0.0, 1.0, 0.0]),
                    ("NFFT T3phi",  [0.0, 0.0, 1.0]),
                    ("NFFT all",    [1.0, 1.0, 1.0])]
    fg, fo = make_fg(ft, data, nx, w)
    global all_pass &= check_gradient(label, fg, fo, vec(x))
end

# ── Cross-check: DFT vs NFFT forward observables ────────────────────────────
println("\n=== Forward observable consistency (DFT vs NFFT) ===")

cvis_dft  = image_to_vis(x, dft)
cvis_nfft = image_to_vis(x, ft)

v2_dft  = abs2.(cvis_dft[data.indx_v2])
v2_nfft = abs2.(cvis_nfft[data.indx_v2])
v2_err  = maximum(abs.(v2_dft .- v2_nfft)) / max(maximum(abs.(v2_dft)), 1e-30)
@printf("  V2 DFT vs NFFT max rel diff:    %.2e\n", v2_err)

_, _, t3phi_dft  = vis_to_t3(cvis_dft,  data.indx_t3_1, data.indx_t3_2, data.indx_t3_3)
_, _, t3phi_nfft = vis_to_t3(cvis_nfft, data.indx_t3_1, data.indx_t3_2, data.indx_t3_3)
t3phi_err = maximum(abs.(mod360(t3phi_dft .- t3phi_nfft)))
@printf("  T3phi DFT vs NFFT max abs diff: %.2e deg\n", t3phi_err)

_, t3amp_dft, _  = vis_to_t3(cvis_dft,  data.indx_t3_1, data.indx_t3_2, data.indx_t3_3)
_, t3amp_nfft, _ = vis_to_t3(cvis_nfft, data.indx_t3_1, data.indx_t3_2, data.indx_t3_3)
t3amp_err = maximum(abs.(t3amp_dft .- t3amp_nfft)) / max(maximum(abs.(t3amp_dft)), 1e-30)
@printf("  T3amp DFT vs NFFT max rel diff: %.2e\n", t3amp_err)

# ── Summary ──────────────────────────────────────────────────────────────────
println()
if all_pass
    println("ALL GRADIENT TESTS PASSED")
else
    println("SOME GRADIENT TESTS FAILED")
    exit(1)
end
