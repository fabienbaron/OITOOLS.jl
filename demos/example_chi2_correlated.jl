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
    data = readoifits(oifitsfile, filter_bad_data=false)[1,1]

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

    # Build a mapping from (baseline_pair, wavelength_index) → global corr index.
    # We assign baselines a canonical order from the first complete timestamp.
    baseline_order = Dict{Tuple{Int,Int}, Int}()  # (sta1,sta2) → 0-based baseline index
    for b in 0:nbase-1
        idx = mjd1_mask[b * nwave + 1]
        pair = (data.v2_sta_index[1, idx], data.v2_sta_index[2, idx])
        baseline_order[pair] = b
    end
    # Build per-point wavelength index within each baseline
    # Within a timestamp, points for one baseline appear as consecutive nwave entries
    # with wavelengths in order.
    lam_sorted = data.v2_lam[mjd1_mask[1:nwave]]  # reference wavelength grid

    # Build corrected cidx for every point
    cidx_fixed = zeros(Int, data.nv2)
    for mjd in unique_mjds
        mask = findall(data.v2_mjd .== mjd)
        for idx in mask
            pair = (data.v2_sta_index[1, idx], data.v2_sta_index[2, idx])
            b = get(baseline_order, pair, -1)
            b == -1 && continue
            # Find wavelength channel (exact match in the reference grid)
            lam_val = data.v2_lam[idx]
            w = findfirst(==(lam_val), lam_sorted)
            w === nothing && continue
            cidx_fixed[idx] = b * nwave + w  # 1-based global index
        end
    end

    ngood = count(>(0), cidx_fixed)
    println("Points with valid corr index: $ngood / $(data.nv2)")

    # Dense correlation matrix — small enough for GRAVITY
    C = Matrix(data.v2_corr)

    # ── 3. Point-source model (V² = 1) ───────────────────────────────────────
    chi2_corr_total, chi2_diag_total, npts_corr, npts_diag =
        _correlated_chi2(ones(data.nv2), data, C, cidx_fixed, unique_mjds, ndata_corr)

    println("\n── Results (point-source model, V² = 1) ──")
    println("Points with correlation: $npts_corr")
    println("Points diagonal only:    $npts_diag")
    println("Total correlated chi2:  $(round(chi2_corr_total, digits=2))")
    println("Total diagonal chi2:    $(round(chi2_diag_total, digits=2))")
    println("Correlated chi2/N:      $(round(chi2_corr_total/data.nv2, digits=4))")
    println("Diagonal chi2/N:        $(round(chi2_diag_total/data.nv2, digits=4))")
    println("Ratio corr/diag:        $(round(chi2_corr_total/chi2_diag_total, digits=4))")

    # ── 4. UD grid search ────────────────────────────────────────────────────
    println("\n── UD grid search ──")

    best_diam = 0.0
    best_chi2 = Inf

    B   = data.v2_baseline
    lam = data.v2_lam

    for diam in 0.1:0.05:5.0
        x      = @. pi * diam * (pi/180/3600/1000) * B / lam
        V      = @. ifelse(abs(x) < 1e-12, 1.0, 2 * besselj1(x) / x)
        v2_mod = @. V^2

        chi2_total, _, _, _ = _correlated_chi2(v2_mod, data, C, cidx_fixed,
                                               unique_mjds, ndata_corr)

        if chi2_total < best_chi2
            best_chi2 = chi2_total
            best_diam = diam
        end
    end

    println("Best UD diameter:       $(best_diam) mas")
    println("Best correlated chi2:   $(round(best_chi2, digits=2))")
    println("Best chi2/N:            $(round(best_chi2/data.nv2, digits=4))")
end

# ── Helper: correlated chi2 summed over all timestamps ────────────────────────
# For each timestamp, extract points with valid correlation indices,
# build the sub-block of Σ, and compute rᵀ Σ⁻¹ r.
# Points without valid correlation indices use diagonal chi2.
function _correlated_chi2(v2_mod, data, C, cidx_fixed, unique_mjds, ndata_corr)
    chi2_corr_total = 0.0
    chi2_diag_total = 0.0
    npts_corr = 0
    npts_diag = 0

    for mjd in unique_mjds
        mask = findall(data.v2_mjd .== mjd)
        cidx = cidx_fixed[mask]

        # Split into points with and without valid correlation
        good = findall(c -> c > 0 && c <= ndata_corr, cidx)
        bad  = findall(c -> c <= 0 || c > ndata_corr, cidx)

        # Diagonal fallback for points without correlation
        if !isempty(bad)
            r_bad = data.v2[mask[bad]] .- v2_mod[mask[bad]]
            chi2_d = sum((r_bad ./ data.v2_err[mask[bad]]).^2)
            chi2_corr_total += chi2_d
            chi2_diag_total += chi2_d
            npts_diag += length(bad)
        end

        # Correlated chi2 for points with valid indices
        if !isempty(good)
            idx_g  = mask[good]
            cidx_g = cidx[good]
            r      = data.v2[idx_g] .- v2_mod[idx_g]
            v2_err = data.v2_err[idx_g]

            D     = Diagonal(v2_err)
            Sigma = Symmetric(D * C[cidx_g, cidx_g] * D)

            chi2_corr_total += dot(r, Sigma \ r)
            chi2_diag_total += sum((r ./ v2_err).^2)
            npts_corr += length(good)
        end
    end

    return chi2_corr_total, chi2_diag_total, npts_corr, npts_diag
end

main()
