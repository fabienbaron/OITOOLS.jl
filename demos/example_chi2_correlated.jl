# Example: Correlated V² chi2 using OI_CORR from GRAVITY data
#
# Reads demos/data/BC2026/OBJECT2_K.oifits, extracts V² data and the
# associated OI_CORR correlation matrix, builds the full covariance,
# and computes chi2 = r' Sigma^{-1} r  via LDLt factorization.
#
# All code is self-contained — nothing is added to src/.

using OITOOLS
using SparseArrays, LinearAlgebra, SpecialFunctions

# ── 1. Read the data ─────────────────────────────────────────────────────────
oifitsfile = joinpath(@__DIR__, "data", "BC2026", "OBJECT2_K.oifits")
data = readoifits(oifitsfile)[1,1]

println("V² points:  $(data.nv2)")
println("Corr matrix size: $(size(data.v2_corr))")
println("Corr idx range:   $(extrema(data.v2_corr_idx))")

# ── 2. Group V² data by MJD (one timestamp = one correlated block) ───────────
# Each timestamp has nbase × nwave = 6 × 233 = 1398 V² points,
# matching the NDATA=1398 of the OI_CORR block.

unique_mjds = sort(unique(data.v2_mjd))
println("Number of timestamps: $(length(unique_mjds))")

nwave   = 233   # spectral channels per baseline (GRAVITY K-band)
nbase   = 6     # baselines for 4 telescopes
nblock  = nwave * nbase  # 1398

# ── 3. Simple model: point source (V² = 1 everywhere) ────────────────────────
v2_model = ones(data.nv2)

# ── 4. Build covariance and compute chi2 per timestamp ────────────────────────
# OI_CORR stores the correlation matrix C (diagonal = 1, off-diagonal = rho).
# Covariance:  Sigma = diag(sigma) * C * diag(sigma)
# chi2 = r' * Sigma^{-1} * r   via  LDLt(Symmetric(Sigma))

C = Matrix(data.v2_corr)  # dense 1398×1398 — small enough for GRAVITY

chi2_correlated_total = 0.0
chi2_diagonal_total   = 0.0

for mjd in unique_mjds
    # indices of V² points for this timestamp
    mask = findall(data.v2_mjd .== mjd)
    n = length(mask)

    if n != nblock
        @warn "Timestamp MJD=$mjd has $n points, expected $nblock — skipping"
        continue
    end

    # Correlation indices for this block
    cidx = data.v2_corr_idx[mask]

    # Check that indices form a contiguous 1:nblock block
    if minimum(cidx) < 1
        @warn "Timestamp MJD=$mjd has zero correlation indices — skipping"
        continue
    end

    # Extract the sub-block of the correlation matrix
    C_block = C[cidx, cidx]

    # Data and errors for this block
    v2_obs = data.v2[mask]
    v2_err = data.v2_err[mask]
    v2_mod = v2_model[mask]

    # Residuals
    r = v2_obs .- v2_mod

    # Build covariance:  Sigma = diag(sigma) * C * diag(sigma)
    D = Diagonal(v2_err)
    Sigma = D * C_block * D

    # ── Correlated chi2 via LDLt ──
    F = ldlt(Symmetric(Sigma))
    chi2_corr = dot(r, F \ r)
    chi2_correlated_total += chi2_corr

    # ── Diagonal chi2 for comparison ──
    chi2_diag = sum((r ./ v2_err).^2)
    chi2_diagonal_total += chi2_diag
end

ndof = data.nv2  # total number of V² points (no free params for point source)

println("\n── Results (point-source model, V² = 1) ──")
println("Total correlated chi2:  $(round(chi2_correlated_total, digits=2))")
println("Total diagonal chi2:    $(round(chi2_diagonal_total, digits=2))")
println("N data points:          $ndof")
println("Correlated chi2/N:      $(round(chi2_correlated_total/ndof, digits=4))")
println("Diagonal chi2/N:        $(round(chi2_diagonal_total/ndof, digits=4))")
println("Ratio corr/diag:        $(round(chi2_correlated_total/chi2_diagonal_total, digits=4))")

# ── 5. Now try with a uniform-disk model ─────────────────────────────────────
println("\n── Trying UD model fit ──")

# Quick grid search over UD diameter
best_diam = 0.0
best_chi2 = Inf

for diam in 0.1:0.05:5.0
    v2_mod = Float64[]
    for i in 1:data.nv2
        B = data.v2_baseline[i]
        lam = data.v2_lam[i]
        # spatial frequency in rad^{-1}
        x = pi * diam * (pi/180/3600/1000) * B / lam
        if abs(x) < 1e-12
            push!(v2_mod, 1.0)
        else
            V = 2 * besselj1(x) / x
            push!(v2_mod, V^2)
        end
    end

    # Compute correlated chi2
    chi2_total = 0.0
    for mjd in unique_mjds
        mask = findall(data.v2_mjd .== mjd)
        length(mask) != nblock && continue
        cidx = data.v2_corr_idx[mask]
        minimum(cidx) < 1 && continue

        r = data.v2[mask] .- v2_mod[mask]
        D = Diagonal(data.v2_err[mask])
        Sigma = D * C[cidx, cidx] * D
        F = ldlt(Symmetric(Sigma))
        chi2_total += dot(r, F \ r)
    end

    if chi2_total < best_chi2
        best_chi2 = chi2_total
        best_diam = diam
    end
end

println("Best UD diameter:       $(best_diam) mas")
println("Best correlated chi2:   $(round(best_chi2, digits=2))")
println("Best chi2/N:            $(round(best_chi2/ndof, digits=4))")

