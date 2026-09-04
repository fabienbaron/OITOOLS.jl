# ============================================================================
# Differential phase: forward, adjoint and JVP.
# ============================================================================

using NFFT

"""
    DiffPhaseConfig

Precomputed cross-channel differential phase configuration.
Holds per-channel VIS NFFT plans, common-baseline indices, and observed data.
"""
struct DiffPhaseConfig
    ft_vis::Vector{NFFTPlan}       # per-channel VIS NFFT plan (ft[2])
    cidx::Vector{Vector{Int}}      # per-channel indices into common baselines
    ncommon::Int                   # number of common baselines
    diffphi_data::Matrix{Float64}  # ncommon × nwavs — observed (degrees)
    diffphi_err::Matrix{Float64}   # ncommon × nwavs — errors (degrees)
end

# `mod360` is OITOOLS', imported by the extension module. This file used to carry its own
# copy — the same wrap, hard-coded to Float64 where the core one follows the array's element
# type. Two definitions of one convention is precisely the drift this merge removes.

"""
    build_diffphase_config(ft, data, nf) -> Union{Nothing, DiffPhaseConfig}

Build a `DiffPhaseConfig` from OITOOLS ft plans and data matrices.
Returns `nothing` if no differential phases are present or `nf == 1`.

Baseline matching uses `(sta1, sta2, round(mjd, 6))` keys, same as OITOOLS.
"""
function build_diffphase_config(ft::AbstractMatrix, data::AbstractMatrix, nf::Int)
    nf <= 1 && return nothing
    data[1,1].phityp == "differential" || return nothing

    # Per-channel VIS plan, by name: `NFFTCell` puts it at index 2, but naming it says which
    # plan this is rather than where it happens to sit.
    ft_vis = NFFTPlan[]
    for c in 1:nf
        ft_ch = ft[c, 1]
        push!(ft_vis, ft_ch isa NFFTCell ? ft_ch.vis : ft_ch[2])
    end

    # Build per-channel baseline keys: (sta1, sta2, round(mjd, 6))
    vis_keys = Vector{Vector{Tuple{Int,Int,Float64}}}(undef, nf)
    for c in 1:nf
        d = data[c, 1]
        nv = length(d.indx_vis)
        vis_keys[c] = [(d.vis_sta_index[1,j], d.vis_sta_index[2,j],
                        round(d.vis_mjd[j], digits=6)) for j in 1:nv]
    end

    # Common baselines present in every channel
    common_set = reduce(intersect, Set.(vis_keys))
    if isempty(common_set)
        return nothing
    end

    # Per-channel index arrays for common baselines
    cidx = Vector{Vector{Int}}(undef, nf)
    for c in 1:nf
        cidx[c] = [j for (j, k) in enumerate(vis_keys[c]) if k in common_set]
    end
    ncommon = length(cidx[1])

    # Extract observed diffphi data and errors
    diffphi_data = hcat([data[c,1].visphi[cidx[c]] for c in 1:nf]...)   # ncommon × nf
    diffphi_err  = hcat([data[c,1].visphi_err[cidx[c]] for c in 1:nf]...) # ncommon × nf

    return DiffPhaseConfig(ft_vis, cidx, ncommon, diffphi_data, diffphi_err)
end


# ============================================================================
# Forward: image cube → differential phases
# ============================================================================

"""
    diffphi_forward(image, dp) -> Matrix{Float64} (ncommon × nwavs)

Compute model differential phases from a polychromatic image cube.

`image` is `npix × npix × nf`. For each channel, complex visibilities are
computed via NFFT, then differential phase is:
    φ_diff_i = arg(V_i / V_ref_i) · 180/π
where V_ref_i = mean(V_j for j ≠ i).
"""
function diffphi_forward(image::AbstractArray{<:Real, 3}, dp::DiffPhaseConfig)
    nf = length(dp.ft_vis)
    ncommon = dp.ncommon

    # Compute per-channel cvis at common baselines
    cvis_all = Matrix{ComplexF64}(undef, ncommon, nf)
    for c in 1:nf
        img_c = image[:, :, c]
        flux = sum(img_c)
        img_norm = Complex{Float64}.(img_c ./ flux)
        cvis_full = dp.ft_vis[c] * img_norm
        cvis_all[:, c] = cvis_full[dp.cidx[c]]
    end

    # Reference visibility and differential phase
    cvis_sum = vec(sum(cvis_all, dims=2))
    diffphi_model = Matrix{Float64}(undef, ncommon, nf)
    for c in 1:nf
        cvis_ref = (cvis_sum .- cvis_all[:, c]) ./ (nf - 1)
        diffphi_model[:, c] = angle.(cvis_all[:, c] ./ cvis_ref) .* (180.0 / π)
    end

    return diffphi_model
end


# ============================================================================
# Adjoint (VJP): diffphi cotangent → image gradient
# ============================================================================

"""
    diffphi_adjoint(g_diffphi, image, dp) -> g_image (npix × npix × nf)

Adjoint (VJP) of differential phase operator.

Maps a cotangent in diffphi-space back to image-space gradients.
Includes both direct and cross-channel terms from the chain rule through
V_ref_i = mean(V_j for j ≠ i).

Uses our observe.jl convention: `g_img_norm = Re(adjoint(ft) * g_cvis)`.
The Wirtinger conjugate gradient of `φ = arg(z)·180/π` is
`g_cvis = i·(180/π)·g_φ·z/|z|²` (where z is the complex visibility ratio).
"""
function diffphi_adjoint(g_diffphi::AbstractMatrix{Float64},
                         image::AbstractArray{<:Real, 3},
                         dp::DiffPhaseConfig)
    nf = length(dp.ft_vis)
    npix = size(image, 1)
    ncommon = dp.ncommon

    # Forward pass to get cvis at common baselines
    cvis_all = Matrix{ComplexF64}(undef, ncommon, nf)
    fluxes = Vector{Float64}(undef, nf)
    for c in 1:nf
        img_c = image[:, :, c]
        fluxes[c] = sum(img_c)
        img_norm = Complex{Float64}.(img_c ./ fluxes[c])
        cvis_full = dp.ft_vis[c] * img_norm
        cvis_all[:, c] = cvis_full[dp.cidx[c]]
    end
    cvis_sum = vec(sum(cvis_all, dims=2))

    # Reference visibilities
    cvis_ref = Matrix{ComplexF64}(undef, ncommon, nf)
    for c in 1:nf
        cvis_ref[:, c] = (cvis_sum .- cvis_all[:, c]) ./ (nf - 1)
    end

    # Wirtinger gradient of φ = arg(z)·180/π w.r.t. complex z:
    #   g_z = i·(180/π)·g_φ / conj(z) = i·(180/π)·g_φ·z/|z|²
    # For φ_c = arg(V_c / V_ref_c), with z_c = V_c/V_ref_c:
    #   d arg(z) = Im(dz/z) = Im(dV/V - dVref/Vref)
    # Direct contribution to g_cvis_k from ∂φ_k/∂V_k:
    #   g_cvis_k += i·(180/π)·g_φ_k·V_k/|V_k|²
    # Cross-channel: V_ref_i depends on V_k for i≠k via mean,
    # contributing -1/(nf-1) of the Vref term.
    scale = im * (180.0 / π)

    g_image = zeros(npix, npix, nf)

    for k in 1:nf
        img_k = image[:, :, k]
        cvis_full_k = dp.ft_vis[k] * Complex{Float64}.(img_k ./ fluxes[k])

        g_cvis = zeros(ComplexF64, length(cvis_full_k))

        # Direct term: ∂φ_k/∂V_k → g_cvis_k = i·(180/π)·g_φ_k·V_k/|V_k|²
        g_cvis[dp.cidx[k]] .= scale .* g_diffphi[:, k] .* cvis_all[:, k] ./ abs2.(cvis_all[:, k])

        # Cross-channel terms: channel k contributes to V_ref_i for all i≠k
        # V_ref_i = (Σ V_j - V_i)/(nf-1), so ∂V_ref_i/∂V_k = 1/(nf-1) for i≠k
        # g_cvis_ref_i = -i·(180/π)·g_φ_i·V_ref_i/|V_ref_i|²  (from the -dVref/Vref term)
        # g_cvis_k += g_cvis_ref_i · 1/(nf-1)
        for i in 1:nf
            i == k && continue
            g_cvis[dp.cidx[k]] .+= (-scale / (nf - 1)) .* g_diffphi[:, i] .* cvis_ref[:, i] ./ abs2.(cvis_ref[:, i])
        end

        # NFFT adjoint → image gradient (Re convention, matching observe.jl)
        g_img_norm = real.(vec(adjoint(dp.ft_vis[k]) * g_cvis))
        g_img_k = reshape(g_img_norm, npix, npix)

        # Flux normalization correction: image was normalized by flux before NFFT
        inv_flux = 1.0 / fluxes[k]
        img_norm_k = img_k .* inv_flux
        correction = sum(g_img_k .* img_norm_k)
        g_image[:, :, k] .= (g_img_k .- correction) .* inv_flux
    end

    return g_image
end


# ============================================================================
# JVP: image tangent → diffphi tangent
# ============================================================================

"""
    diffphi_jvp(image, d_image, dp) -> d_diffphi (ncommon × nwavs)

JVP (Jacobian-vector product) of the differential phase operator.

Maps an image-space tangent to a diffphi-space tangent, accounting for
cross-channel coupling through the reference visibility.
"""
function diffphi_jvp(image::AbstractArray{<:Real, 3},
                     d_image::AbstractArray{<:Real, 3},
                     dp::DiffPhaseConfig)
    nf = length(dp.ft_vis)
    ncommon = dp.ncommon

    # Forward: compute cvis and d_cvis at common baselines
    cvis_all = Matrix{ComplexF64}(undef, ncommon, nf)
    d_cvis_all = Matrix{ComplexF64}(undef, ncommon, nf)
    for c in 1:nf
        img_c = image[:, :, c]
        d_img_c = d_image[:, :, c]
        flux = sum(img_c)
        d_flux = sum(d_img_c)

        # Tangent of img/flux: (d_img·flux - img·d_flux) / flux²
        img_norm = Complex{Float64}.(img_c ./ flux)
        d_img_norm = Complex{Float64}.((d_img_c .* flux .- img_c .* d_flux) ./ flux^2)

        cvis_full = dp.ft_vis[c] * img_norm
        d_cvis_full = dp.ft_vis[c] * d_img_norm

        cvis_all[:, c] = cvis_full[dp.cidx[c]]
        d_cvis_all[:, c] = d_cvis_full[dp.cidx[c]]
    end

    # Tangent of reference visibility
    cvis_sum = vec(sum(cvis_all, dims=2))
    d_cvis_sum = vec(sum(d_cvis_all, dims=2))

    # Tangent of diffphi: d arg(V_i / V_ref_i) · 180/π
    # arg(z) derivative: d arg(z) = Im(dz / z)
    # z_i = V_i / V_ref_i, dz_i = (dV_i · V_ref_i - V_i · dV_ref_i) / V_ref_i²
    # So d arg(z_i) = Im(dz_i / z_i) = Im((dV_i/V_i) - (dV_ref_i/V_ref_i))
    d_diffphi = Matrix{Float64}(undef, ncommon, nf)
    for c in 1:nf
        cvis_ref = (cvis_sum .- cvis_all[:, c]) ./ (nf - 1)
        d_cvis_ref = (d_cvis_sum .- d_cvis_all[:, c]) ./ (nf - 1)

        # d arg(V/Vref) = Im(dV/V - dVref/Vref) · 180/π
        d_diffphi[:, c] = imag.(d_cvis_all[:, c] ./ cvis_all[:, c] .-
                                d_cvis_ref ./ cvis_ref) .* (180.0 / π)
    end

    return d_diffphi
end
