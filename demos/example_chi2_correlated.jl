# Example: Correlated V² chi2 using OI_CORR from GRAVITY data
#
# Reads demos/data/BC2026/OBJECT2_K.oifits, extracts V² data and the
# associated OI_CORR correlation matrix, builds the full covariance,
# and computes chi2 = r' Σ⁻¹ r  via Bunch-Kaufman (LDLᵀ with pivoting).
#
# All code is self-contained — nothing is added to src/.
#
# GRAVITY quirk: all rows in OI_VIS2 have CORRINDX_VIS2DATA = 1 even though
# the correlation block is 1398×1398 (6 baselines × 233 channels).
# The per-baseline offset is implicit from row order within each timestamp.
# This example reconstructs the correct global correlation index per point.

using OITOOLS
using SparseArrays, LinearAlgebra, SpecialFunctions

function main()
    # ── 1. Read the data ─────────────────────────────────────────────────────
    oifitsfile = joinpath(@__DIR__, "data", "BC2026", "OBJECT2_K.oifits")
    data = readoifits(oifitsfile)[1,1]

    println("V² points:  $(data.nv2)")
    println("Corr matrix size: $(size(data.v2_corr))")

    ndata_corr = size(data.v2_corr, 1)  # 1398

    # ── 2. Build correct global correlation indices ──────────────────────────
    # Group by MJD. Within each timestamp, assign baselines their proper
    # offset into the correlation block: baseline k → channels (k-1)*nwave+1 : k*nwave
    unique_mjds = unique(data.v2_mjd)
    println("Number of timestamps: $(length(unique_mjds))")

    # Detect nwave from data: within one timestamp, count points on first baseline
    mjd1_mask = findall(data.v2_mjd .== data.v2_mjd[1])
    sta1 = data.v2_sta_index[:, mjd1_mask[1]]  # first baseline's station pair
    nwave = count(i -> data.v2_sta_index[:, mjd1_mask[i]] == sta1, 1:length(mjd1_mask))
    nbase = length(mjd1_mask) ÷ nwave
    println("Detected: nwave=$nwave  nbase=$nbase  nblock=$(nwave*nbase)")

    # Build corrected cidx: within each timestamp, consecutive groups of nwave
    # points get offsets 0, nwave, 2*nwave, ...
    cidx_fixed = zeros(Int, data.nv2)
    for mjd in unique_mjds
        mask = findall(data.v2_mjd .== mjd)
        n = length(mask)
        if n != nwave * nbase
            @warn "MJD=$mjd has $n points, expected $(nwave*nbase) — skipping"
            continue
        end
        for b in 0:nbase-1
            offset = b * nwave
            for w in 1:nwave
                idx = mask[b * nwave + w]
                cidx_fixed[idx] = offset + w   # 1-based global index
            end
        end
    end

    println("Fixed cidx range: $(extrema(cidx_fixed[cidx_fixed .> 0]))")

    # Dense correlation matrix — small enough for GRAVITY
    C = Matrix(data.v2_corr)

    # ── 3. Point-source model (V² = 1) ───────────────────────────────────────
    chi2_corr_total, chi2_diag_total, nskipped =
        _correlated_chi2(ones(data.nv2), data, C, cidx_fixed, unique_mjds,
                         nwave * nbase, ndata_corr)

    ndof = data.nv2
    println("\n── Results (point-source model, V² = 1) ──")
    println("Timestamps processed:   $(length(unique_mjds) - nskipped) / $(length(unique_mjds))")
    println("Total correlated chi2:  $(round(chi2_corr_total, digits=2))")
    println("Total diagonal chi2:    $(round(chi2_diag_total, digits=2))")
    println("N data points:          $ndof")
    println("Correlated chi2/N:      $(round(chi2_corr_total/ndof, digits=4))")
    println("Diagonal chi2/N:        $(round(chi2_diag_total/ndof, digits=4))")
    println("Ratio corr/diag:        $(round(chi2_corr_total/chi2_diag_total, digits=4))")

    # ── 4. UD grid search ────────────────────────────────────────────────────
    println("\n── UD grid search ──")

    best_diam = 0.0
    best_chi2 = Inf
    nblock = nwave * nbase

    B   = data.v2_baseline
    lam = data.v2_lam

    for diam in 0.1:0.05:5.0
        x      = @. pi * diam * (pi/180/3600/1000) * B / lam
        V      = @. ifelse(abs(x) < 1e-12, 1.0, 2 * besselj1(x) / x)
        v2_mod = @. V^2

        chi2_total, _, _ = _correlated_chi2(v2_mod, data, C, cidx_fixed,
                                            unique_mjds, nblock, ndata_corr)

        if chi2_total < best_chi2
            best_chi2 = chi2_total
            best_diam = diam
        end
    end

    println("Best UD diameter:       $(best_diam) mas")
    println("Best correlated chi2:   $(round(best_chi2, digits=2))")
    println("Best chi2/N:            $(round(best_chi2/ndof, digits=4))")
end

# ── Helper: correlated chi2 summed over all timestamps ────────────────────────
function _correlated_chi2(v2_mod, data, C, cidx_fixed, unique_mjds,
                          nblock, ndata_corr)
    chi2_corr_total = 0.0
    chi2_diag_total = 0.0
    nskipped = 0

    for mjd in unique_mjds
        mask = findall(data.v2_mjd .== mjd)
        cidx = cidx_fixed[mask]

        r      = data.v2[mask] .- v2_mod[mask]
        v2_err = data.v2_err[mask]

        # Fall back to diagonal if invalid indices
        if any(cidx .== 0) || maximum(cidx) > ndata_corr || length(mask) != nblock
            chi2_d = sum((r ./ v2_err).^2)
            chi2_corr_total += chi2_d
            chi2_diag_total += chi2_d
            nskipped += 1
            continue
        end

        # Covariance: Σ = diag(σ) · C_block · diag(σ)
        D     = Diagonal(v2_err)
        Sigma = Symmetric(D * C[cidx, cidx] * D)

        chi2_corr_total += dot(r, Sigma \ r)
        chi2_diag_total += sum((r ./ v2_err).^2)
    end

    return chi2_corr_total, chi2_diag_total, nskipped
end

main()
