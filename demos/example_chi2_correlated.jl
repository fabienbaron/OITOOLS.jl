# Example: Correlated V² chi2 using OI_CORR from GRAVITY data
#
# Reads demos/data/BC2026/OBJECT2_K.oifits, extracts V² data and the
# associated OI_CORR correlation matrix, whitens the data, and computes
# chi2 = |w|² where w = L⁻¹ r  (L from Cholesky of covariance).
#
# All code is self-contained — nothing is added to src/.

using OITOOLS
using SparseArrays, LinearAlgebra, SpecialFunctions

"""
    CorrBlock

Precomputed factorization for one correlated timestamp block.
- `mask`:  indices into the flat data arrays
- `L_inv`: inverse of the lower-Cholesky factor of Σ = diag(σ)·C·diag(σ),
           so that  w = L_inv * r  gives whitened residuals with chi2 = |w|².
"""
struct CorrBlock
    mask::Vector{Int}
    L_inv::Matrix{Float64}
end

"""
    precompute_blocks(data, C, ndata_corr) -> Vector{CorrBlock}

Factorize the covariance for each timestamp once.  Returns a vector of
`CorrBlock`s; each holds the data-point mask and the precomputed L⁻¹.
Points without valid correlation indices get a diagonal-only block
(L_inv = diag(1/σ)).
"""
function precompute_blocks(data, C, ndata_corr)
    unique_mjds = unique(data.v2_mjd)
    blocks = CorrBlock[]

    for mjd in unique_mjds
        mask = findall(data.v2_mjd .== mjd)
        cidx = data.v2_corr_idx[mask]
        σ    = data.v2_err[mask]

        good = findall(c -> c > 0 && c <= ndata_corr, cidx)
        bad  = findall(c -> c <= 0 || c > ndata_corr, cidx)

        # Diagonal-only block for points without correlation
        if !isempty(bad)
            L_inv_diag = Diagonal(1.0 ./ σ[bad])
            push!(blocks, CorrBlock(mask[bad], Matrix(L_inv_diag)))
        end

        # Correlated block
        if !isempty(good)
            D = Diagonal(σ[good])
            C_block = C[cidx[good], cidx[good]]
            Σ = Symmetric(D * C_block * D)
            L_inv = try
                L = cholesky(Σ).L
                Matrix(inv(L))
            catch
                # Fallback: eigendecomposition with clamped eigenvalues
                F = eigen(Σ)
                λ = max.(F.values, 1e-12)
                Matrix((Diagonal(1.0 ./ sqrt.(λ)) * F.vectors')')
            end
            push!(blocks, CorrBlock(mask[good], L_inv))
        end
    end
    return blocks
end

"""
    chi2_whitened(v2_mod, data, blocks) -> (chi2_corr, chi2_diag)

Compute chi2 using precomputed whitening blocks.  No factorization happens here.
"""
function chi2_whitened(v2_mod, data, blocks)
    chi2_corr = 0.0
    chi2_diag = 0.0
    for blk in blocks
        r = data.v2[blk.mask] .- v2_mod[blk.mask]
        w = blk.L_inv * r
        chi2_corr += dot(w, w)
        chi2_diag += sum((r ./ data.v2_err[blk.mask]).^2)
    end
    return chi2_corr, chi2_diag
end

function main()
    # ── 1. Read the data ─────────────────────────────────────────────────────
    oifitsfile = joinpath(@__DIR__, "data", "BC2026", "OBJECT2_K.oifits")
    data = readoifits(oifitsfile, filter_bad_data=false)[1,1]

    println("V² points:  $(data.nv2)")
    println("Corr matrix size: $(size(data.v2_corr))")
    println("Corr idx range:   $(extrema(data.v2_corr_idx))")

    ndata_corr = size(data.v2_corr, 1)
    C = Matrix(data.v2_corr)

    # ── 2. Precompute whitening blocks (once) ────────────────────────────────
    print("Precomputing whitening blocks... ")
    blocks = precompute_blocks(data, C, ndata_corr)
    println("done ($(length(blocks)) blocks)")

    # ── 3. Point-source model (V² = 1) ───────────────────────────────────────
    chi2_corr, chi2_diag = chi2_whitened(ones(data.nv2), data, blocks)

    println("\n── Results (point-source model, V² = 1) ──")
    println("Total correlated chi2:  $(round(chi2_corr, digits=2))")
    println("Total diagonal chi2:    $(round(chi2_diag, digits=2))")
    println("Correlated chi2/N:      $(round(chi2_corr/data.nv2, digits=4))")
    println("Diagonal chi2/N:        $(round(chi2_diag/data.nv2, digits=4))")
    println("Ratio corr/diag:        $(round(chi2_corr/chi2_diag, digits=4))")

    # ── 4. UD grid search (fast — no factorization in the loop) ──────────────
    println("\n── UD grid search ──")
    best_diam = 0.0
    best_chi2 = Inf

    B   = data.v2_baseline
    lam = data.v2_lam

    for diam in 0.1:0.05:5.0
        x      = @. pi * diam * (pi/180/3600/1000) * B / lam
        V      = @. ifelse(abs(x) < 1e-12, 1.0, 2 * besselj1(x) / x)
        v2_mod = @. V^2

        chi2, _ = chi2_whitened(v2_mod, data, blocks)
        if chi2 < best_chi2
            best_chi2 = chi2
            best_diam = diam
        end
    end

    println("Best UD diameter:       $(best_diam) mas")
    println("Best correlated chi2:   $(round(best_chi2, digits=2))")
    println("Best chi2/N:            $(round(best_chi2/data.nv2, digits=4))")
end

main()
