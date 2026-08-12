# chi2_flat.jl
#
# Chi-squared evaluation and gradient for flat-dict parametric model fitting.
# Operates on FlatModel (from parse_model.jl) and OIdata (from readoifits.jl).
#
# Load order:
#   include("test_resolvers.jl")     # SharedUtils, RGF
#   include("test_hankel.jl")        # Hankel + chainrules
#   include("model_chainrules.jl")   # vis_ud, vis_ldlin, ...
#   include("parse_model.jl")        # FlatModel, eval_model, eval_model_grad
#   include("readoifits.jl")         # OIdata struct
#   include("chi2_flat.jl")          # this file
#
# Interface
# ---------
#   chi2_flat(model, x, data; weights, verb, vonmises)
#       → T           (forward only — cheap, for initial checks or gradient-free optimizers)
#
#   chi2_flat_fg(model, x, data; weights, verb, vonmises)
#       → (T, Vector{T})               (chi2 + gradient — for gradient-based optimizers)
#
#   chi2_flat_reduced(model, x, data; weights, vonmises)
#       → T           (chi2 / ndof — for human-readable goodness-of-fit)
#
# Polychromatic interface (data::Vector{<:OIdata}):
#   chi2_flat(model, x, data_array; weights, verb, vonmises)
#   chi2_flat_fg(model, x, data_array; weights, verb, vonmises)
#       Evaluates the model at each channel's UV points and optionally computes
#       differential phase across wavelength channels.
#
# Weights convention (7-element):
#   weights = [w_v2, w_t3amp, w_t3phi, w_visamp, w_visphi, w_flux, w_diffphase]
#   Default: [1, 1, 1, 0, 0, 0, 0]  (V2 + T3 only, matching the imaging code default)
#   Set w_visamp/w_visphi > 0 to include absolute visamp/visphi.
#   Set w_flux > 0 to include OI_FLUX data.
#   Set w_diffphase > 0 to include differential phase (polychromatic only).
#   For backward compatibility, 5-element weights are zero-padded to 7.
#
# Gradient derivation
# -------------------
# Let V ∈ ℂ^{nuv} = eval_model(model, x, data.uv)
# Let J ∈ ℂ^{nuv × nparams} = ∂V/∂x  (from eval_model_grad)
#
# For a scalar observable O(V) and residual r = ∂chi2/∂O, the parameter
# gradient is computed via the Wirtinger chain rule through the real-valued
# loss.  The resulting formulas are identical in structure to the DFT matrix
# formulas in oichi2.jl, with J replacing dft:
#
#   V2:      g += 4·Re( Jᵀ[indx_v2,:]  · (r_v2  .* conj(V[indx_v2]))          )
#   T3amp:   g += Re(   Jᵀ[indx_t3_k,:]· (r_t3a .* conj(V[k]) ./ |V[k]| ·
#                                          |V[other legs]|)  ) summed over k=1,2,3
#   T3phi:   g += Im(   Jᵀ[indx_t3_k,:]· (r_t3p .* conj(V[k]) ./ |V[k]|²)    ) summed k
#   visamp:  g += Re(   Jᵀ[indx_vis,:] · (r_va  .* conj(V[indx_vis]) ./ |V[indx_vis]|) )
#   visphi:  g += Im(   Jᵀ[indx_vis,:] · (r_vp  .* conj(V[indx_vis]) ./ |V[indx_vis]|²))
#   flux:    g += 2·C·Σ(r_flux)·Re(J_zero)  (C=1 if calibrated; C_flux constant if uncalibrated)
#   diffphi: g += Im(   Jᵀ[indx_vis,:] · (r_dp  .* conj(V[indx_vis]) ./ |V[indx_vis]|²))
#                 (per-channel, following oichi2.jl leading-order approximation)
#
# where Jᵀ is the plain (non-conjugate) transpose, matching the DFT convention.

using LinearAlgebra, Printf


# ─────────────────────────────────────────────────────────────────────────────
# _pad_weights: ensure weights is 7-element (backward compat with 5-element)
# ─────────────────────────────────────────────────────────────────────────────

_pad_weights(w::AbstractVector{<:Real}) = length(w) >= 7 ? w[1:7] : vcat(w, zeros(7 - length(w)))

# Working precision for a chi2 evaluation: whatever the visibilities and the data promote to.
# The default weights are a Float64 literal vector, so they are converted rather than allowed
# to drag a Float32 pipeline up to Float64.
_chi2_eltype(V::AbstractVector{<:Complex}, data::OIdata{T}) where {T} =
    promote_type(real(eltype(V)), T)

_pad_weights(w::AbstractVector{<:Real}, ::Type{T}) where {T<:AbstractFloat} =
    convert(Vector{T}, _pad_weights(w))


# ─────────────────────────────────────────────────────────────────────────────
# mod360: wrap phase residuals to [-180, 180]  (same as oichi2.jl)
# ─────────────────────────────────────────────────────────────────────────────

function _mod360(x)
    T = float(eltype(x))
    return mod.(mod.(x .+ T(180), T(360)) .+ T(360), T(360)) .- T(180)
end


# ─────────────────────────────────────────────────────────────────────────────
# _chi2_terms: compute all chi2 components from model cvis
# Returns a NamedTuple with the values needed by both forward and gradient paths.
# ─────────────────────────────────────────────────────────────────────────────

function _chi2_terms(V::AbstractVector{<:Complex},
                     data::OIdata,
                     weights::AbstractVector{<:Real},
                     vonmises::Bool)

    T = _chi2_eltype(V, data)
    w = _pad_weights(weights, T)
    w_v2, w_t3amp, w_t3phi, w_visamp, w_visphi = w[1], w[2], w[3], w[4], w[5]

    # ── V2 ────────────────────────────────────────────────────────────────────
    chi2_v2 = zero(T)
    r_v2    = T[]                   # residual scale: (v2_model - v2) / err²
    v2_model = T[]

    if w_v2 > 0 && data.nv2 > 0
        v2_model = abs2.(V[data.indx_v2])
        r_v2     = (v2_model .- data.v2) ./ data.v2_err.^2
        chi2_v2  = sum(r_v2 .* (v2_model .- data.v2))   # = norm(...)^2
    end

    # ── T3 products ───────────────────────────────────────────────────────────
    chi2_t3amp = zero(T)
    chi2_t3phi = zero(T)
    r_t3amp    = T[]
    r_t3phi    = T[]
    t3amp_model = T[]
    t3phi_model = T[]
    V1 = Complex{T}[]; V2c = Complex{T}[]; V3 = Complex{T}[]  # per-leg cvis

    need_t3 = (w_t3amp > 0 && data.nt3amp > 0) || (w_t3phi > 0 && data.nt3phi > 0)
    if need_t3
        V1  = V[data.indx_t3_1]
        V2c = V[data.indx_t3_2]
        V3  = V[data.indx_t3_3]
        t3  = V1 .* V2c .* V3
        t3amp_model = abs.(t3)
        t3phi_model = angle.(t3) .* T(180 / π)
    end

    if w_t3amp > 0 && data.nt3amp > 0
        r_t3amp   = T(2) .* (t3amp_model .- data.t3amp) ./ data.t3amp_err.^2
        chi2_t3amp = sum(r_t3amp .* (t3amp_model .- data.t3amp)) / 2   # norm(...)^2
    end

    if w_t3phi > 0 && data.nt3phi > 0
        if !vonmises
            dphi       = _mod360(t3phi_model .- data.t3phi)
            r_t3phi    = T(360 / π) .* dphi ./ data.t3phi_err.^2
            chi2_t3phi = sum(dphi .* dphi ./ data.t3phi_err.^2)
        else
            dphi_rad   = (t3phi_model .- data.t3phi) .* T(π / 180)
            r_t3phi    = T(2) .* data.t3phi_vonmises_err .* sin.(dphi_rad)
            chi2_t3phi = sum(T(-2) .* data.t3phi_vonmises_err .* cos.(dphi_rad)
                             .+ data.t3phi_vonmises_chi2_offset)
        end
    end

    # ── Absolute VISAMP / VISPHI ───────────────────────────────────────────────
    chi2_visamp  = zero(T)
    chi2_visphi  = zero(T)
    r_visamp     = T[]
    r_visphi     = T[]
    Vvis         = Complex{T}[]
    visamp_model = T[]
    visphi_model = T[]

    need_vis = (w_visamp > 0 && data.nvisamp > 0) || (w_visphi > 0 && data.nvisphi > 0)
    if need_vis
        Vvis         = V[data.indx_vis]
        visamp_model = abs.(Vvis)
        visphi_model = angle.(Vvis) .* T(180 / π)
    end

    if w_visamp > 0 && data.nvisamp > 0
        r_visamp    = T(2) .* (visamp_model .- data.visamp) ./ data.visamp_err.^2
        chi2_visamp = sum(r_visamp .* (visamp_model .- data.visamp)) / 2
    end

    if w_visphi > 0 && data.nvisphi > 0
        dphi        = _mod360(visphi_model .- data.visphi)
        r_visphi    = T(360 / π) .* dphi ./ data.visphi_err.^2
        chi2_visphi = sum(dphi .* dphi ./ data.visphi_err.^2)
    end

    return (;
        # scalar chi2 components
        chi2_v2, chi2_t3amp, chi2_t3phi, chi2_visamp, chi2_visphi,
        # model observables (needed for verb output)
        v2_model, t3amp_model, t3phi_model, visamp_model, visphi_model,
        # weighted residuals (needed for gradient)
        r_v2, r_t3amp, r_t3phi, r_visamp, r_visphi,
        # per-leg cvis slices (needed for gradient)
        V1, V2=V2c, V3, Vvis,
    )
end


# ─────────────────────────────────────────────────────────────────────────────
# _ndof: number of data points (weighted count)
# ─────────────────────────────────────────────────────────────────────────────

function _ndof(data::OIdata, weights::AbstractVector{<:Real})
    w = _pad_weights(weights)
    n  = (w[1] > 0) * data.nv2
    n += (w[2] > 0) * data.nt3amp
    n += (w[3] > 0) * data.nt3phi
    n += (w[4] > 0) * data.nvisamp
    n += (w[5] > 0) * data.nvisphi
    n += (w[6] > 0) * data.nflux
    return max(n, 1)
end

function _ndof(data::Vector{<:OIdata}, weights::AbstractVector{<:Real})
    n = sum(_ndof(d, weights) for d in data)
    w = _pad_weights(weights)
    if w[7] > 0 && length(data) > 1
        n += sum(d.nvisphi for d in data)
    end
    return max(n, 1)
end


# ─────────────────────────────────────────────────────────────────────────────
# _total_chi2: weighted sum of components
# ─────────────────────────────────────────────────────────────────────────────

function _total_chi2(t, weights; chi2_flux=zero(eltype(weights)), chi2_diffphase=zero(eltype(weights)))
    w = _pad_weights(weights)
    return (w[1] * t.chi2_v2
          + w[2] * t.chi2_t3amp
          + w[3] * t.chi2_t3phi
          + w[4] * t.chi2_visamp
          + w[5] * t.chi2_visphi
          + w[6] * chi2_flux
          + w[7] * chi2_diffphase)
end


# ─────────────────────────────────────────────────────────────────────────────
# _accumulate_g_cvis!: compute the complex adjoint source vector from chi2 terms.
#
# Convention: g_cvis is defined so that for any linear map V = A * z (z ∈ ℝⁿ):
#   ∂chi²/∂z = real(Aᵀ * g_cvis)
#
# This means phase-type observables (T3phi, visphi) contribute via -im factors,
# so that imag(Jᵀ * rhs) becomes real(Jᵀ * (-im * rhs)) inside the combined g_cvis.
#
# Uses scatter-add loops to correctly handle duplicate T3 indices.
# ─────────────────────────────────────────────────────────────────────────────

function _accumulate_g_cvis!(g_cvis::AbstractVector{<:Complex},
                             V::AbstractVector{<:Complex},
                             t::NamedTuple, data::OIdata,
                             weights::AbstractVector{<:Real})
    w = _pad_weights(weights)
    w_v2, w_t3amp, w_t3phi, w_visamp, w_visphi = w[1], w[2], w[3], w[4], w[5]
    # Divide-safety floor at the working precision: eps(Float64) would be nine orders of
    # magnitude below eps(Float32) and so no guard at all for a Float32 pipeline.
    ε = eps(real(eltype(V)))

    # V2: g_cvis[k] += w_v2 * 4 * r_v2 * conj(V[k])
    if w_v2 > 0 && data.nv2 > 0
        for j in eachindex(t.r_v2)
            k = data.indx_v2[j]
            g_cvis[k] += w_v2 * 4 * t.r_v2[j] * conj(V[k])
        end
    end

    # T3amp: scatter-add over three legs
    if w_t3amp > 0 && data.nt3amp > 0
        a1 = abs.(t.V1); a2 = abs.(t.V2); a3 = abs.(t.V3)
        safe1 = max.(a1, ε); safe2 = max.(a2, ε); safe3 = max.(a3, ε)
        for j in eachindex(t.r_t3amp)
            c1 = t.r_t3amp[j] * conj(t.V1[j]) / safe1[j] * a2[j] * a3[j]
            c2 = t.r_t3amp[j] * conj(t.V2[j]) / safe2[j] * a1[j] * a3[j]
            c3 = t.r_t3amp[j] * conj(t.V3[j]) / safe3[j] * a1[j] * a2[j]
            g_cvis[data.indx_t3_1[j]] += w_t3amp * c1
            g_cvis[data.indx_t3_2[j]] += w_t3amp * c2
            g_cvis[data.indx_t3_3[j]] += w_t3amp * c3
        end
    end

    # T3phi: scatter-add with -im factor
    if w_t3phi > 0 && data.nt3phi > 0
        safe1 = max.(abs2.(t.V1), ε)
        safe2 = max.(abs2.(t.V2), ε)
        safe3 = max.(abs2.(t.V3), ε)
        for j in eachindex(t.r_t3phi)
            c1 = t.r_t3phi[j] * conj(t.V1[j]) / safe1[j]
            c2 = t.r_t3phi[j] * conj(t.V2[j]) / safe2[j]
            c3 = t.r_t3phi[j] * conj(t.V3[j]) / safe3[j]
            g_cvis[data.indx_t3_1[j]] += -im * w_t3phi * c1
            g_cvis[data.indx_t3_2[j]] += -im * w_t3phi * c2
            g_cvis[data.indx_t3_3[j]] += -im * w_t3phi * c3
        end
    end

    # Visamp: g_cvis[k] += w_visamp * r_visamp * conj(Vvis) / |Vvis|
    if w_visamp > 0 && data.nvisamp > 0
        safe = max.(abs.(t.Vvis), ε)
        for j in eachindex(t.r_visamp)
            k = data.indx_vis[j]
            g_cvis[k] += w_visamp * t.r_visamp[j] * conj(t.Vvis[j]) / safe[j]
        end
    end

    # Visphi: g_cvis[k] += -im * w_visphi * r_visphi * conj(Vvis) / |Vvis|²
    if w_visphi > 0 && data.nvisphi > 0
        safe = max.(abs2.(t.Vvis), ε)
        for j in eachindex(t.r_visphi)
            k = data.indx_vis[j]
            g_cvis[k] += -im * w_visphi * t.r_visphi[j] * conj(t.Vvis[j]) / safe[j]
        end
    end

    return g_cvis
end


# ─────────────────────────────────────────────────────────────────────────────
# cvis_to_chi2_fg: chi2 + complex adjoint source from complex visibilities.
#
# Building block for all chi2+gradient code paths (model, image, SPARCO, hybrid).
# Handles: V², T3amp, T3phi, absolute visamp, absolute visphi.
# Does NOT handle: flux (separate observable), differential observables (cross-channel).
#
# Returns (chi2, g_cvis) where g_cvis satisfies:
#   g_params = real(Jᵀ · g_cvis)              for any model Jacobian J
#   g_image  = real(DFTᵀ · g_cvis)            for DFT image reconstruction
#   g_image  = real(NFFT^H · conj(g_cvis))    for NFFT image reconstruction
# ─────────────────────────────────────────────────────────────────────────────

"""
    cvis_to_chi2_f(V, data; weights=[1,1,1,0,0,0,0], vonmises=false, model_flux=NaN) -> chi2

Chi-squared from complex visibilities `V` and interferometric data `data`, without the
gradient.

Handles the same per-channel observables as [`cvis_to_chi2_fg`](@ref) — V², T3amp, T3phi,
visamp, visphi — and additionally the **flux** term: pass the model's zero-baseline flux as
`model_flux` and give `weights[6] > 0`. (`cvis_to_chi2_fg` cannot do this, because the flux
*gradient* needs the model Jacobian at zero baseline, which is not derivable from `V` alone.)

Precision follows the inputs: `ComplexF32` visibilities against an `OIdata{Float32}` return a
`Float32` chi². Accumulation is via `sum`, which is pairwise, so Float32 holds ~1e-8 relative
accuracy even for millions of residuals.
"""
function cvis_to_chi2_f(V::AbstractVector{<:Complex}, data::OIdata;
                        weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                        vonmises::Bool = false,
                        model_flux::Real = NaN)
    T = _chi2_eltype(V, data)
    w = _pad_weights(weights, T)
    t = _chi2_terms(V, data, w, vonmises)
    chi2_fl = zero(T)
    if w[6] > 0 && data.nflux > 0 && isfinite(model_flux)
        chi2_fl = T(_chi2_flux(model_flux, data).chi2)
    end
    return _total_chi2(t, w; chi2_flux=chi2_fl)
end

"""
    cvis_to_chi2_fg(V, data; weights=[1,1,1,0,0,0,0], vonmises=false) -> (chi2, g_cvis)

Compute chi-squared and the complex adjoint source `g_cvis` from complex
visibilities `V` and interferometric data `data`.

This is the core building block for all gradient-based chi² code paths
(model fitting, image reconstruction, SPARCO, hybrid model+image).
It handles per-channel observables: V², T3amp, T3phi, visamp, visphi.
It does **not** handle flux or differential phase (cross-channel observables).

Precision follows the inputs rather than being forced: `ComplexF32` visibilities against an
`OIdata{Float32}` return a `Float32` chi² and a `Vector{ComplexF32}` adjoint source, so a
Float32 pipeline stays in Float32 end to end. See [`cvis_to_chi2_f`](@ref) if you do not
need the gradient.

The adjoint source `g_cvis` satisfies:
- `g_params = real(Jᵀ · g_cvis)` for any model Jacobian J
- `g_image  = real(DFTᵀ · g_cvis)` for DFT image reconstruction
- `g_image  = real(NFFT^H · conj(g_cvis))` for NFFT image reconstruction
"""
function cvis_to_chi2_fg(V::AbstractVector{<:Complex}, data::OIdata;
                         weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                         vonmises::Bool = false)
    T = _chi2_eltype(V, data)
    w = _pad_weights(weights, T)
    t = _chi2_terms(V, data, w, vonmises)
    chi2 = _total_chi2(t, w)
    g_cvis = zeros(Complex{T}, length(V))
    _accumulate_g_cvis!(g_cvis, V, t, data, w)
    return chi2, g_cvis
end


# ─────────────────────────────────────────────────────────────────────────────
# _chi2_flux: compute flux chi2.
# If flux is calibrated (CALSTAT=C), compare model directly to data.
# If uncalibrated (CALSTAT=U), fit optimal scaling C = Σ(fm·fd·w) / Σ(fm²·w).
# ─────────────────────────────────────────────────────────────────────────────

function _chi2_flux(model_flux::Real, data::OIdata)
    fm = fill(model_flux, data.nflux)
    if data.flux_calibrated
        residual = (fm .- data.flux) ./ data.flux_err.^2
        chi2 = sum(((fm .- data.flux) ./ data.flux_err).^2)
        C = 1.0
    else
        w  = 1.0 ./ data.flux_err.^2
        C  = sum(fm .* data.flux .* w) / sum(fm.^2 .* w)
        residual = (C .* fm .- data.flux) ./ data.flux_err.^2
        chi2 = sum(((C .* fm .- data.flux) ./ data.flux_err).^2)
    end
    return (; chi2, C, residual)
end


# ─────────────────────────────────────────────────────────────────────────────
# _chi2_flux_polychromatic: flux chi2 across multiple wavelength channels.
# Calibrated: direct comparison (no scaling).
# Uncalibrated: single C_flux across all channels.
# ─────────────────────────────────────────────────────────────────────────────

function _chi2_flux_polychromatic(model_fluxes::AbstractVector{<:Real}, data::Vector{<:OIdata})
    Tf = promote_type(float(eltype(model_fluxes)), eltype(first(data).flux))
    fm = Tf[]; fd = Tf[]; fe = Tf[]; chan_idx = Int[]
    nwavs = length(data)
    for i in 1:nwavs
        d = data[i]; d.nflux > 0 || continue
        append!(fm, fill(model_fluxes[i], d.nflux))
        append!(fd, d.flux)
        append!(fe, d.flux_err)
        append!(chan_idx, fill(i, d.nflux))
    end
    isempty(fm) && return (; chi2=zero(Tf), C=one(Tf), residual=Tf[], chan_idx=Int[])
    calibrated = !isempty(data) && data[1].flux_calibrated
    if calibrated
        residual = (fm .- fd) ./ fe.^2
        chi2 = sum(((fm .- fd) ./ fe).^2)
        C = one(Tf)
    else
        w = one(Tf) ./ fe.^2
        C = sum(fm .* fd .* w) / sum(fm.^2 .* w)
        residual = (C .* fm .- fd) ./ fe.^2
        chi2 = sum(((C .* fm .- fd) ./ fe).^2)
    end
    return (; chi2, C, residual, chan_idx)
end


# ─────────────────────────────────────────────────────────────────────────────
# _verb_print: match the color-coded output format of oichi2.jl
# ─────────────────────────────────────────────────────────────────────────────

function _verb_print(t, data::OIdata, weights::AbstractVector{<:Real}, chi2_total::Real;
                     chi2_flux::Real=0.0, C_flux::Real=NaN,
                     chi2_diffphase::Real=0.0)
    w = _pad_weights(weights)
    ndof = _ndof(data, weights)
    w[1] > 0 && data.nv2    > 0 && printstyled(@sprintf("V2: %.3f ",    t.chi2_v2    / data.nv2),    color=:red)
    w[2] > 0 && data.nt3amp > 0 && printstyled(@sprintf("T3A: %.3f ",   t.chi2_t3amp / data.nt3amp), color=:blue)
    w[3] > 0 && data.nt3phi > 0 && printstyled(@sprintf("T3P: %.3f ",   t.chi2_t3phi / data.nt3phi), color=:green)
    w[4] > 0 && data.nvisamp> 0 && printstyled(@sprintf("VA: %.3f ",    t.chi2_visamp/ data.nvisamp),color=:magenta)
    w[5] > 0 && data.nvisphi> 0 && printstyled(@sprintf("VP: %.3f ",    t.chi2_visphi/ data.nvisphi),color=:cyan)
    if w[6] > 0 && data.nflux > 0
        if C_flux == 1.0
            printstyled(@sprintf("FLUX: %.3f ", chi2_flux / data.nflux), color=:yellow)
        else
            printstyled(@sprintf("FLUX: %.3f (C=%.4f) ", chi2_flux / data.nflux, C_flux), color=:yellow)
        end
    end
    printstyled(@sprintf("χ²r: %.3f\n", chi2_total / ndof), color=:white)
end


# ─────────────────────────────────────────────────────────────────────────────
# chi2_flat  (forward only — monochromatic)
# ─────────────────────────────────────────────────────────────────────────────

"""
    chi2_flat(model, x, data::OIdata; weights=[1,1,1,0,0,0,0], verb=false, vonmises=false)
        -> T

Evaluate the weighted chi-squared for `model` at parameter vector `x` against
interferometric data `data`.

No gradient is computed.  Use `chi2_flat_fg` for gradient-based optimizers.

`weights` selects which observables contribute:
[V2, T3amp, T3phi, visamp, visphi, flux, diffphase].
5-element weights are zero-padded to 7 for backward compatibility.
"""
function chi2_flat(model::FlatModel,
                   x::AbstractVector,
                   data::OIdata;
                   weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                   verb::Bool = false,
                   vonmises::Bool = false)

    w = _pad_weights(weights)
    need_flux = w[6] > 0 && data.nflux > 0

    # Evaluate model (prepend zero baseline if flux needed)
    # Pass per-UV wavelength and MJD so chromatic expressions ($WL, $MJD) work.
    _wl  = hasproperty(data, :uv_lam) ? data.uv_lam : nothing
    _mjd = hasproperty(data, :uv_mjd) ? data.uv_mjd : nothing
    if need_flux
        uv_aug  = hcat(zeros(2, 1), data.uv)
        wl_aug  = isnothing(_wl)  ? nothing : vcat([_wl[1]], _wl)
        mjd_aug = isnothing(_mjd) ? nothing : vcat([_mjd[1]], _mjd)
        V_aug   = eval_model(model, x, uv_aug; wl=wl_aug, mjd=mjd_aug)
        model_flux = real(V_aug[1])
        V = V_aug[2:end]
    else
        V = eval_model(model, x, data.uv; wl=_wl, mjd=_mjd)
        model_flux = NaN
    end

    chi2 = cvis_to_chi2_f(V, data; weights=w, vonmises=vonmises, model_flux=model_flux)

    if verb
        t = _chi2_terms(V, data, w, vonmises)
        chi2_fl = 0.0; C_fl = NaN
        if need_flux
            fl = _chi2_flux(model_flux, data)
            chi2_fl = fl.chi2; C_fl = fl.C
        end
        _verb_print(t, data, w, chi2; chi2_flux=chi2_fl, C_flux=C_fl)
    end
    return chi2
end


# ─────────────────────────────────────────────────────────────────────────────
# chi2_flat_fg  (chi2 + gradient — monochromatic)
# ─────────────────────────────────────────────────────────────────────────────

"""
    chi2_flat_fg(model, x, data::OIdata; weights=[1,1,1,0,0,0,0], verb=false, vonmises=false)
        -> (T, Vector{T})

Evaluate chi-squared and its gradient w.r.t. `x`.

The gradient is computed analytically via the Wirtinger chain rule applied to
the complex visibility Jacobian J = ∂V/∂x from `eval_model_grad`.

Returns `(chi2, grad)`; the element type follows `x` and `data.uv`, so Float32 inputs give a Float32 gradient.
"""
function chi2_flat_fg(model::FlatModel,
                      x::AbstractVector,
                      data::OIdata;
                      weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                      verb::Bool = false,
                      vonmises::Bool = false)

    w = _pad_weights(weights)
    w_flux = w[6]
    need_flux = w_flux > 0 && data.nflux > 0

    # ── Forward pass + Jacobian ─────────────────────────────────────────────
    _wl  = hasproperty(data, :uv_lam) ? data.uv_lam : nothing
    _mjd = hasproperty(data, :uv_mjd) ? data.uv_mjd : nothing
    if need_flux
        uv_aug  = hcat(zeros(2, 1), data.uv)
        wl_aug  = isnothing(_wl)  ? nothing : vcat([_wl[1]], _wl)
        mjd_aug = isnothing(_mjd) ? nothing : vcat([_mjd[1]], _mjd)
        V_aug, J_aug = eval_model_grad(model, x, uv_aug; wl=wl_aug, mjd=mjd_aug)
        model_flux = real(V_aug[1])
        J_zero = J_aug[1:1, :]
        V = V_aug[2:end]
        J = J_aug[2:end, :]
    else
        V, J = eval_model_grad(model, x, data.uv; wl=_wl, mjd=_mjd)
        model_flux = NaN
    end

    # ── Per-channel chi2 + adjoint source ──────────────────────────────────
    chi2_obs, g_cvis = cvis_to_chi2_fg(V, data; weights=w, vonmises=vonmises)

    # ── Flux chi2 ───────────────────────────────────────────────────────────
    Tc = typeof(chi2_obs)
    chi2_fl = zero(Tc); C_fl = Tc(NaN); r_flux = Tc[]
    if need_flux
        fl = _chi2_flux(model_flux, data)
        chi2_fl = fl.chi2; C_fl = fl.C; r_flux = fl.residual
    end

    chi2 = chi2_obs + Tc(w_flux) * chi2_fl

    # ── Gradient via single matrix-vector product ──────────────────────────
    g = real.(transpose(J) * g_cvis)

    # Flux gradient
    if need_flux
        g .+= w_flux .* 2 .* C_fl .* sum(r_flux) .* real.(vec(J_zero))
    end

    if verb
        t = _chi2_terms(V, data, w, vonmises)
        _verb_print(t, data, w, chi2; chi2_flux=chi2_fl, C_flux=C_fl)
    end
    return chi2, g
end


# ─────────────────────────────────────────────────────────────────────────────
# chi2_flat  (forward only — polychromatic)
# ─────────────────────────────────────────────────────────────────────────────

"""
    chi2_flat(model, x, data::Vector{<:OIdata}; weights=[1,1,1,0,0,0,0], verb=false, vonmises=false)
        -> T

Polychromatic chi-squared.  Evaluates the model at each channel's UV points and
computes per-channel V2/T3/visamp/visphi/flux plus cross-channel differential phase.

Differential phase follows the ASPRO2 convention used in oichi2.jl:
  V_ref_i = mean of all other channels' cvis at the same baselines
  diffphi_i = angle(V_i / V_ref_i)
"""
function chi2_flat(model::FlatModel,
                   x::AbstractVector,
                   data::Vector{<:OIdata};
                   weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                   verb::Bool = false,
                   vonmises::Bool = false)

    nwavs = length(data)
    # Precision follows the caller, as in the monochromatic path.
    Tp = promote_type(float(eltype(x)), eltype(first(data).uv))
    w = _pad_weights(weights, Tp)

    # Evaluate model at each channel's UV points
    need_flux = w[6] > 0 && any(d.nflux > 0 for d in data)
    if need_flux
        cvis = Vector{Vector{Complex{Tp}}}(undef, nwavs)
        model_fluxes = Vector{Tp}(undef, nwavs)
        for i in 1:nwavs
            _wl_i  = hasproperty(data[i], :uv_lam) ? data[i].uv_lam : nothing
            _mjd_i = hasproperty(data[i], :uv_mjd) ? data[i].uv_mjd : nothing
            uv_aug  = hcat(zeros(2, 1), data[i].uv)
            wl_aug  = isnothing(_wl_i)  ? nothing : vcat([_wl_i[1]], _wl_i)
            mjd_aug = isnothing(_mjd_i) ? nothing : vcat([_mjd_i[1]], _mjd_i)
            V_aug   = eval_model(model, x, uv_aug; wl=wl_aug, mjd=mjd_aug)
            model_fluxes[i] = real(V_aug[1])
            cvis[i] = V_aug[2:end]
        end
    else
        cvis = [eval_model(model, x, data[i].uv;
                           wl  = hasproperty(data[i], :uv_lam) ? data[i].uv_lam : nothing,
                           mjd = hasproperty(data[i], :uv_mjd) ? data[i].uv_mjd : nothing)
                for i in 1:nwavs]
        model_fluxes = Tp[]
    end

    # Per-channel V2, T3, visamp, visphi
    chi2 = zero(Tp)
    for i in 1:nwavs
        t = _chi2_terms(cvis[i], data[i], w, vonmises)
        chi2 += _total_chi2(t, w)  # flux/diffphase weights handled separately below
    end

    # Flux across all channels (single C_flux)
    chi2_fl = zero(Tp)
    if need_flux
        fl = _chi2_flux_polychromatic(model_fluxes, data)
        chi2_fl = Tp(fl.chi2)
        chi2 += w[6] * chi2_fl
    end

    # Differential phase across channels
    chi2_dp = zero(Tp)
    if w[7] > 0 && nwavs > 1
        chi2_dp = Tp(_diffphase_chi2(cvis, data, nwavs))
        chi2 += w[7] * chi2_dp
    end

    if verb
        ndof = _ndof(data, w)
        printstyled(@sprintf("χ²r: %.3f\n", chi2 / ndof), color=:white)
    end
    return chi2
end


# ─────────────────────────────────────────────────────────────────────────────
# _diffphase_chi2: differential phase chi2 across wavelength channels
# ─────────────────────────────────────────────────────────────────────────────

function _diffphase_chi2(cvis, data, nwavs)
    nvis = data[1].nvisphi
    nvis > 0 || return 0.0
    cvis_vis = hcat([cvis[i][data[i].indx_vis] for i in 1:nwavs]...)
    cvis_sum = vec(sum(cvis_vis, dims=2))
    chi2 = 0.0
    for i in 1:nwavs
        cvis_ref = (cvis_sum .- cvis_vis[:, i]) ./ (nwavs - 1)
        diffphi_model = angle.(cvis_vis[:, i] ./ cvis_ref) .* (180.0 / π)
        dphi = _mod360(diffphi_model .- data[i].visphi)
        chi2 += sum((dphi ./ data[i].visphi_err).^2)
    end
    return chi2
end


# ─────────────────────────────────────────────────────────────────────────────
# chi2_flat_fg  (chi2 + gradient — polychromatic)
# ─────────────────────────────────────────────────────────────────────────────

"""
    chi2_flat_fg(model, x, data::Vector{<:OIdata}; weights=[1,1,1,0,0,0,0], verb=false, vonmises=false)
        -> (T, Vector{T})

Polychromatic chi-squared and gradient.  Computes per-channel V2/T3/visamp/visphi/flux
gradients plus cross-channel differential phase gradient.
"""
function chi2_flat_fg(model::FlatModel,
                      x::AbstractVector,
                      data::Vector{<:OIdata};
                      weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                      verb::Bool = false,
                      vonmises::Bool = false)

    nwavs = length(data)
    nparams = length(x)
    Tp = promote_type(float(eltype(x)), eltype(first(data).uv))
    w = _pad_weights(weights, Tp)
    w_flux, w_dp = w[6], w[7]

    need_flux = w_flux > 0 && any(d.nflux > 0 for d in data)
    need_dp   = w_dp > 0 && nwavs > 1

    # ── Evaluate model + Jacobian at all channels ──────────────────────────
    cvis = Vector{Vector{Complex{Tp}}}(undef, nwavs)
    Jac  = Vector{Matrix{Complex{Tp}}}(undef, nwavs)
    model_fluxes = Vector{Tp}(undef, nwavs)
    J_zeros = Vector{Vector{Complex{Tp}}}(undef, nwavs)

    for i in 1:nwavs
        _wl_i  = hasproperty(data[i], :uv_lam) ? data[i].uv_lam : nothing
        _mjd_i = hasproperty(data[i], :uv_mjd) ? data[i].uv_mjd : nothing
        if need_flux
            uv_aug  = hcat(zeros(2, 1), data[i].uv)
            wl_aug  = isnothing(_wl_i)  ? nothing : vcat([_wl_i[1]], _wl_i)
            mjd_aug = isnothing(_mjd_i) ? nothing : vcat([_mjd_i[1]], _mjd_i)
            V_aug, J_aug = eval_model_grad(model, x, uv_aug; wl=wl_aug, mjd=mjd_aug)
            model_fluxes[i] = real(V_aug[1])
            J_zeros[i] = vec(J_aug[1, :])
            cvis[i] = V_aug[2:end]
            Jac[i]  = J_aug[2:end, :]
        else
            V_i, J_i = eval_model_grad(model, x, data[i].uv; wl=_wl_i, mjd=_mjd_i)
            cvis[i] = V_i
            Jac[i]  = J_i
        end
    end

    # ── Per-channel chi2 + gradient ────────────────────────────────────────
    chi2 = zero(Tp)
    g = zeros(Tp, nparams)

    for i in 1:nwavs
        chi2_i, g_cvis_i = cvis_to_chi2_fg(cvis[i], data[i]; weights=w, vonmises=vonmises)
        chi2 += chi2_i
        g .+= real.(transpose(Jac[i]) * g_cvis_i)
    end

    # ── Flux chi2 + gradient (across all channels) ─────────────────────────
    chi2_fl = zero(Tp)
    if need_flux
        fl = _chi2_flux_polychromatic(model_fluxes, data)
        chi2_fl = Tp(fl.chi2)
        chi2 += w_flux * chi2_fl
        # Gradient: for each channel, g += 2·C·Σ(r_flux for that channel)·Re(J_zero_i)
        for i in 1:nwavs
            idx = findall(fl.chan_idx .== i)
            isempty(idx) && continue
            g .+= w_flux .* 2.0 .* fl.C .* sum(fl.residual[idx]) .* real.(J_zeros[i])
        end
    end

    # ── Differential phase chi2 + gradient ─────────────────────────────────
    if need_dp
        nvis = data[1].nvisphi
        if nvis > 0
            cvis_vis = hcat([cvis[i][data[i].indx_vis] for i in 1:nwavs]...)
            cvis_sum = vec(sum(cvis_vis, dims=2))
            for i in 1:nwavs
                cvis_ref = (cvis_sum .- cvis_vis[:, i]) ./ (nwavs - 1)
                diffphi_model = angle.(cvis_vis[:, i] ./ cvis_ref) .* (180.0 / π)
                dphi = _mod360(diffphi_model .- data[i].visphi)
                chi2 += w_dp * sum((dphi ./ data[i].visphi_err).^2)
                # Gradient: leading-order through V_i only (following oichi2.jl)
                ε = eps(real(eltype(cvis_vis)))
                dT3 = -360 / π .* dphi ./ data[i].visphi_err.^2
                safe = max.(abs2.(cvis_vis[:, i]), ε)
                rhs = dT3 .* conj.(cvis_vis[:, i]) ./ safe
                g .+= w_dp .* imag.(transpose(Jac[i][data[i].indx_vis, :]) * rhs)
            end
        end
    end

    if verb
        ndof = _ndof(data, w)
        printstyled(@sprintf("χ²r: %.3f\n", chi2 / ndof), color=:white)
    end
    return chi2, g
end


# ─────────────────────────────────────────────────────────────────────────────
# residuals_flat / residuals_flat_jac — residual vector and Jacobian for LsqFit
#
# LsqFit minimises Σ r_i².  We define r_i = (model_i - data_i) / σ_i
# for each enabled observable, concatenated into a single flat vector.
# The Jacobian J[i,j] = ∂r_i/∂x_j is computed analytically from the
# complex visibility Jacobian ∂V/∂x via the Wirtinger chain rule.
# ─────────────────────────────────────────────────────────────────────────────

"""
    residuals_flat(model, x, data; weights, vonmises) -> Vector{Float64}

Compute the weighted residual vector `r` such that `chi2 = dot(r, r)`.

Each element is `(model_observable - data_observable) / σ`, concatenated
in order: [V2; T3amp; T3phi; visamp; visphi; flux].
"""
function residuals_flat(model::FlatModel,
                        x::AbstractVector,
                        data::OIdata;
                        weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                        vonmises::Bool = false)

    w = _pad_weights(weights)
    need_flux = w[6] > 0 && data.nflux > 0

    _wl  = hasproperty(data, :uv_lam) ? data.uv_lam : nothing
    _mjd = hasproperty(data, :uv_mjd) ? data.uv_mjd : nothing
    if need_flux
        uv_aug  = hcat(zeros(2, 1), data.uv)
        wl_aug  = isnothing(_wl)  ? nothing : vcat([_wl[1]], _wl)
        mjd_aug = isnothing(_mjd) ? nothing : vcat([_mjd[1]], _mjd)
        V_aug   = eval_model(model, x, uv_aug; wl=wl_aug, mjd=mjd_aug)
        model_flux = real(V_aug[1])
        V = V_aug[2:end]
    else
        V = eval_model(model, x, data.uv; wl=_wl, mjd=_mjd)
        model_flux = NaN
    end

    r = Float64[]

    # V2
    if w[1] > 0 && data.nv2 > 0
        v2_model = abs2.(V[data.indx_v2])
        append!(r, sqrt(w[1]) .* (v2_model .- data.v2) ./ data.v2_err)
    end

    # T3amp
    if w[2] > 0 && data.nt3amp > 0
        t3 = V[data.indx_t3_1] .* V[data.indx_t3_2] .* V[data.indx_t3_3]
        append!(r, sqrt(w[2]) .* (abs.(t3) .- data.t3amp) ./ data.t3amp_err)
    end

    # T3phi
    if w[3] > 0 && data.nt3phi > 0
        t3 = V[data.indx_t3_1] .* V[data.indx_t3_2] .* V[data.indx_t3_3]
        t3phi_model = angle.(t3) .* T(180 / π)
        append!(r, sqrt(w[3]) .* _mod360(t3phi_model .- data.t3phi) ./ data.t3phi_err)
    end

    # Visamp
    if w[4] > 0 && data.nvisamp > 0
        visamp_model = abs.(V[data.indx_vis])
        append!(r, sqrt(w[4]) .* (visamp_model .- data.visamp) ./ data.visamp_err)
    end

    # Visphi
    if w[5] > 0 && data.nvisphi > 0
        visphi_model = angle.(V[data.indx_vis]) .* (180.0 / π)
        append!(r, sqrt(w[5]) .* _mod360(visphi_model .- data.visphi) ./ data.visphi_err)
    end

    # Flux
    if need_flux
        fl = _chi2_flux(model_flux, data)
        fm = fill(fl.C * model_flux, data.nflux)
        append!(r, sqrt(w[6]) .* (fm .- data.flux) ./ data.flux_err)
    end

    return r
end


"""
    residuals_flat_jac(model, x, data; weights, vonmises) -> (Vector{Float64}, Matrix{Float64})

Compute the weighted residual vector and its Jacobian `J[i,j] = ∂r_i/∂x_j`.

Returns `(r, J)` where `chi2 = dot(r, r)` and `∂chi2/∂x = 2 Jᵀ r`.
"""
function residuals_flat_jac(model::FlatModel,
                            x::AbstractVector,
                            data::OIdata;
                            weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                            vonmises::Bool = false)

    w = _pad_weights(weights)
    w_v2, w_t3amp, w_t3phi, w_visamp, w_visphi, w_flux = w[1], w[2], w[3], w[4], w[5], w[6]
    need_flux = w_flux > 0 && data.nflux > 0
    nparams = length(x)

    _wl  = hasproperty(data, :uv_lam) ? data.uv_lam : nothing
    _mjd = hasproperty(data, :uv_mjd) ? data.uv_mjd : nothing
    if need_flux
        uv_aug  = hcat(zeros(2, 1), data.uv)
        wl_aug  = isnothing(_wl)  ? nothing : vcat([_wl[1]], _wl)
        mjd_aug = isnothing(_mjd) ? nothing : vcat([_mjd[1]], _mjd)
        V_aug, dV_aug = eval_model_grad(model, x, uv_aug; wl=wl_aug, mjd=mjd_aug)
        model_flux = real(V_aug[1])
        V  = V_aug[2:end]
        dV = dV_aug[2:end, :]
    else
        V, dV = eval_model_grad(model, x, data.uv; wl=_wl, mjd=_mjd)
        model_flux = NaN
    end

    # LsqFit consumes these, so they stay Float64 (see fit_model.jl).
    r_parts = Vector{Float64}[]
    J_parts = Matrix{Float64}[]
    ε = eps(Float64)

    # ── V2: r = (|V|² - d) / σ,  ∂r/∂x = 2 Re(V* · ∂V/∂x) / σ ──────────
    if w_v2 > 0 && data.nv2 > 0
        idx = data.indx_v2
        Vv2 = V[idx]
        v2_model = abs2.(Vv2)
        r_v2 = sqrt(w_v2) .* (v2_model .- data.v2) ./ data.v2_err
        # J_v2[i,j] = sqrt(w) * 2 Re(conj(V_i) * dV[i,j]) / σ_i
        J_v2 = sqrt(w_v2) .* 2.0 .* real.(conj.(Vv2) .* dV[idx, :]) ./ data.v2_err
        push!(r_parts, r_v2); push!(J_parts, J_v2)
    end

    # ── T3amp: r = (|T3| - d) / σ ────────────────────────────────────────────
    if w_t3amp > 0 && data.nt3amp > 0
        V1 = V[data.indx_t3_1]; V2c = V[data.indx_t3_2]; V3 = V[data.indx_t3_3]
        t3 = V1 .* V2c .* V3
        t3amp_model = abs.(t3)
        r_t3a = sqrt(w_t3amp) .* (t3amp_model .- data.t3amp) ./ data.t3amp_err
        # ∂|T3|/∂x = Re(conj(T3)/|T3| · ∂T3/∂x)
        # ∂T3/∂x = V2·V3·∂V1/∂x + V1·V3·∂V2/∂x + V1·V2·∂V3/∂x
        safe_amp = max.(t3amp_model, ε)
        phase_factor = conj.(t3) ./ safe_amp
        dT3 = (V2c .* V3) .* dV[data.indx_t3_1, :] .+
               (V1 .* V3) .* dV[data.indx_t3_2, :] .+
               (V1 .* V2c) .* dV[data.indx_t3_3, :]
        J_t3a = sqrt(w_t3amp) .* real.(phase_factor .* dT3) ./ data.t3amp_err
        push!(r_parts, r_t3a); push!(J_parts, J_t3a)
    end

    # ── T3phi: r = mod360(φ_model - φ_data) / σ ──────────────────────────────
    if w_t3phi > 0 && data.nt3phi > 0
        V1 = V[data.indx_t3_1]; V2c = V[data.indx_t3_2]; V3 = V[data.indx_t3_3]
        t3 = V1 .* V2c .* V3
        t3phi_model = angle.(t3) .* T(180 / π)
        dphi = _mod360(t3phi_model .- data.t3phi)
        r_t3p = sqrt(w_t3phi) .* dphi ./ data.t3phi_err
        # ∂angle(T3)/∂x = Im(conj(T3)/|T3|² · ∂T3/∂x) · 180/π
        safe_amp2 = max.(abs2.(t3), ε)
        phase_factor = conj.(t3) ./ safe_amp2
        dT3 = (V2c .* V3) .* dV[data.indx_t3_1, :] .+
               (V1 .* V3) .* dV[data.indx_t3_2, :] .+
               (V1 .* V2c) .* dV[data.indx_t3_3, :]
        J_t3p = sqrt(w_t3phi) .* (180.0 / π) .* imag.(phase_factor .* dT3) ./ data.t3phi_err
        push!(r_parts, r_t3p); push!(J_parts, J_t3p)
    end

    # ── Visamp: r = (|V| - d) / σ ────────────────────────────────────────────
    if w_visamp > 0 && data.nvisamp > 0
        idx = data.indx_vis
        Vvis = V[idx]
        visamp_model = abs.(Vvis)
        r_va = sqrt(w_visamp) .* (visamp_model .- data.visamp) ./ data.visamp_err
        safe_amp = max.(visamp_model, ε)
        J_va = sqrt(w_visamp) .* real.(conj.(Vvis) ./ safe_amp .* dV[idx, :]) ./ data.visamp_err
        push!(r_parts, r_va); push!(J_parts, J_va)
    end

    # ── Visphi: r = mod360(φ_model - φ_data) / σ ─────────────────────────────
    if w_visphi > 0 && data.nvisphi > 0
        idx = data.indx_vis
        Vvis = V[idx]
        visphi_model = angle.(Vvis) .* T(180 / π)
        dphi = _mod360(visphi_model .- data.visphi)
        r_vp = sqrt(w_visphi) .* dphi ./ data.visphi_err
        safe_amp2 = max.(abs2.(Vvis), ε)
        J_vp = sqrt(w_visphi) .* (180.0 / π) .* imag.(conj.(Vvis) ./ safe_amp2 .* dV[idx, :]) ./ data.visphi_err
        push!(r_parts, r_vp); push!(J_parts, J_vp)
    end

    # ── Flux (C=1 for calibrated; for uncalibrated, C treated as constant) ────
    if need_flux
        fl = _chi2_flux(model_flux, data)
        fm = fill(fl.C * model_flux, data.nflux)
        r_fl = sqrt(w_flux) .* (fm .- data.flux) ./ data.flux_err
        # J_flux[i,j] = sqrt(w) * C * Re(∂V_zero/∂x_j) / σ_i
        J_zero = real.(dV_aug[1:1, :])  # 1 × nparams
        J_fl = sqrt(w_flux) .* fl.C .* repeat(J_zero, data.nflux, 1) ./ data.flux_err
        push!(r_parts, r_fl); push!(J_parts, J_fl)
    end

    r = vcat(r_parts...)
    J = vcat(J_parts...)
    return r, J
end


# ─────────────────────────────────────────────────────────────────────────────
# chi2_flat_reduced  (chi2 / ndof — for display)
# ─────────────────────────────────────────────────────────────────────────────

"""
    chi2_flat_reduced(model, x, data; weights=[1,1,1,0,0,0,0], vonmises=false)
        -> T

Return chi2 / ndof (reduced chi-squared).  Useful for assessing goodness of fit.
Works with both monochromatic (OIdata) and polychromatic (Vector{<:OIdata}) data.
"""
function chi2_flat_reduced(model::FlatModel,
                           x::AbstractVector,
                           data::Union{OIdata, Vector{<:OIdata}};
                           weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                           vonmises::Bool = false)
    chi2 = chi2_flat(model, x, data; weights, vonmises)
    return chi2 / _ndof(data, weights)
end


# ─────────────────────────────────────────────────────────────────────────────
# chi2_flat_print  (verbose summary without running a separate forward pass)
# ─────────────────────────────────────────────────────────────────────────────

"""
    chi2_flat_print(model, x, data; weights=[1,1,1,0,0,0,0], vonmises=false)

Print a per-observable chi2 breakdown.  Does a single forward pass.
"""
function chi2_flat_print(model::FlatModel,
                         x::AbstractVector,
                         data::Union{OIdata, Vector{<:OIdata}};
                         weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                         vonmises::Bool = false)
    chi2_flat(model, x, data; weights, verb=true, vonmises)
    nothing
end

"""
    model_to_chi2(model, x, data; weights=[1,1,1,0,0,0,0], verb=false, vonmises=false)

Alias for `chi2_flat`.  Compute the weighted chi-squared of a parametric
model against data, for consistency with the `model_to_obs` / `model_to_vis` family.
"""
const model_to_chi2 = chi2_flat

"""
    model_to_chi2_fg(model, x, data; weights=[1,1,1,0,0,0,0], verb=false, vonmises=false)

Alias for `chi2_flat_fg`.  Compute chi-squared and its gradient w.r.t. `x`,
for consistency with the `model_to_obs` / `model_to_vis` family.
"""
const model_to_chi2_fg = chi2_flat_fg
