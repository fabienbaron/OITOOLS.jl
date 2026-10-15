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
#   chi2_flat(model, x, data; weights, verb)
#       → T           (forward only — cheap, for initial checks or gradient-free optimizers)
#
#   chi2_flat_fg(model, x, data; weights, verb)
#       → (T, Vector{T})               (chi2 + gradient — for gradient-based optimizers)
#
#   chi2_flat_reduced(model, x, data; weights)
#       → T           (chi2 / ndof — for human-readable goodness-of-fit)
#
# Polychromatic interface (data::Vector{<:OIdata}):
#   chi2_flat(model, x, data_array; weights, verb)
#   chi2_flat_fg(model, x, data_array; weights, verb)
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
_chi2_eltype(V::AbstractVecOrMat{<:Complex}, data::OIdata{T}) where {T} =
    promote_type(real(eltype(V)), T)

# `Real`, not `AbstractFloat`: a ForwardDiff Dual is the working type when a gradient is taken
# through the χ² itself, as it is for a smeared model.
_pad_weights(w::AbstractVector{<:Real}, ::Type{T}) where {T<:Real} =
    convert(Vector{T}, _pad_weights(w))


# ─────────────────────────────────────────────────────────────────────────────
# Bandwidth smearing: observables from visibilities across each channel
# ─────────────────────────────────────────────────────────────────────────────
#
# A smeared model gives an N × nsub matrix (`eval_model_sub`), one column per sub-wavelength.
# Every observable is formed per column and averaged over the columns: V² as the mean |V|², a
# triple product as the mean product, a visibility as the mean visibility. For an N-vector the
# same expressions reduce to the plain ones exactly, so one code path serves both.

_rows(V::AbstractVector, i) = V[i]
_rows(V::AbstractMatrix, i) = V[i, :]
_rowmean(A::AbstractVector) = A
_rowmean(A::AbstractMatrix) = vec(sum(A; dims = 2)) ./ size(A, 2)

"""Whether `model` smears over `data`'s channels: `smear` ≥ 2 or auto, and channel widths."""
_smeared(model::FlatModel, data::OIdata) =
    (model.smear >= 2 || model.smear == -1) && hasproperty(data, :uv_dlam) &&
    any(d -> isfinite(d) && d > 0, data.uv_dlam)

"""
    smear_samples(model, data) -> Int

How many sub-wavelengths per channel `model` uses on `data`: `smear` as given, or for
`"smear" => "auto"` CANDID's rule (`_estimateNsmear`), four samples per cycle of the largest
phase change across a channel, at least 3. The phase change is that of a point at the model's
extent (`FlatModel.smear_radius`, fixed when the model is compiled, so a fit never sees the
sampling change) on the data's longest projected baseline.
"""
function smear_samples(model::FlatModel, data::OIdata)
    model.smear == -1 || return model.smear
    cycles = maximum(hypot(data.uv[1, i], data.uv[2, i]) * model.smear_radius / _MAS2RAD_PM *
                     data.uv_dlam[i] / data.uv_lam[i] for i in axes(data.uv, 2); init = 0.0)
    return max(3, ceil(Int, 4 * cycles))
end

# For a criterion that forms its observables from one visibility per point: a smeared model is
# refused rather than evaluated at the channels' centres without saying so.
_refuse_smear(model::FlatModel, what) = (model.smear >= 2 || model.smear == -1) && throw(ArgumentError(
    "$what does not apply bandwidth smearing; remove \"smear\" from the model or set it to 0"))

# The model's total flux at the first point's wavelength and epoch, what OI_FLUX is compared
# with, and its gradient. (`eval_model` divides by it, so V at zero baseline is 1.)
function _flux_at(model::FlatModel, x::AbstractVector, wl, mjd)
    return only(model_flux(model, x; wl = isnothing(wl) ? nothing : wl[1:1],
                           mjd = isnothing(mjd) ? nothing : mjd[1:1]))
end
function _flux_grad_at(model::FlatModel, x::AbstractVector, wl, mjd)
    f(y) = _flux_at(model, y, wl, mjd)
    return f(x), ForwardDiff.gradient(f, x)
end

# The smeared visibilities of `data`'s uv points, and the model's flux when asked for.
function _smeared_cvis(model::FlatModel, x::AbstractVector, data::OIdata, need_flux::Bool)
    wl  = data.uv_lam
    mjd = hasproperty(data, :uv_mjd) ? data.uv_mjd : nothing
    V = eval_model_sub(model, x, data.uv; wl, mjd, dlam = data.uv_dlam,
                       nsub = smear_samples(model, data))
    # Coherent: the visibility averaged over the channel and the observables formed from that,
    # as an ABCD / coherent-flux measurement gives them (and as CANDID models it).
    model.smear_coherent && (V = _rowmean(V))
    flux = need_flux ? _flux_at(model, x, wl, mjd) : NaN
    return V, flux
end


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
                     weights::AbstractVector{<:Real})

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
        dphi       = _mod360(t3phi_model .- data.t3phi)
        r_t3phi    = T(360 / π) .* dphi ./ data.t3phi_err.^2
        chi2_t3phi = sum(dphi .* dphi ./ data.t3phi_err.^2)
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


# ── the smeared χ², kept apart from the unsmeared path ──────────────────────────────────
#
# Separate functions rather than branches inside the unsmeared ones: a branch that assigns
# the visibilities a matrix in one arm widens their type in the other, and the unsmeared path
# then pays for dispatch it never needed.

# The averaged observables of an N × nsub matrix: V² = mean |V|², T3 = mean of the products,
# VIS = mean V.
function _smeared_observables(V::AbstractMatrix, data::OIdata)
    v2  = _rowmean(abs2.(_rows(V, data.indx_v2)))
    t3  = _rowmean(_rows(V, data.indx_t3_1) .* _rows(V, data.indx_t3_2) .* _rows(V, data.indx_t3_3))
    vis = _rowmean(_rows(V, data.indx_vis))
    return v2, t3, vis
end

"""`_chi2_terms` for smeared visibilities: the χ² components of the averaged observables."""
function _chi2_terms(V::AbstractMatrix{<:Complex}, data::OIdata, weights::AbstractVector{<:Real})
    T = real(eltype(V))
    w = _pad_weights(weights, T)
    v2, t3, vis = _smeared_observables(V, data)
    sq(x) = sum(abs2, x)
    t3amp, t3phi = abs.(t3), angle.(t3) .* T(180 / π)
    visamp, visphi = abs.(vis), angle.(vis) .* T(180 / π)
    return (;
        chi2_v2     = w[1] > 0 && data.nv2 > 0    ? sq((v2 .- data.v2) ./ data.v2_err) : zero(T),
        chi2_t3amp  = w[2] > 0 && data.nt3amp > 0 ? sq((t3amp .- data.t3amp) ./ data.t3amp_err) : zero(T),
        chi2_t3phi  = w[3] > 0 && data.nt3phi > 0 ? sq(_mod360(t3phi .- data.t3phi) ./ data.t3phi_err) : zero(T),
        chi2_visamp = w[4] > 0 && data.nvisamp > 0 ? sq((visamp .- data.visamp) ./ data.visamp_err) : zero(T),
        chi2_visphi = w[5] > 0 && data.nvisphi > 0 ? sq(_mod360(visphi .- data.visphi) ./ data.visphi_err) : zero(T),
        v2_model = v2, t3amp_model = t3amp, t3phi_model = t3phi,
        visamp_model = visamp, visphi_model = visphi)
end

function _chi2_flat_smeared(model, x, data::OIdata, weights, verb)
    w = _pad_weights(weights)
    need_flux = w[6] > 0 && data.nflux > 0
    V, model_flux = _smeared_cvis(model, x, data, need_flux)
    t = _chi2_terms(V, data, w)
    chi2_fl = need_flux ? _chi2_flux(model_flux, data).chi2 : zero(real(eltype(V)))
    chi2 = _total_chi2(t, w; chi2_flux = chi2_fl)
    verb && _verb_print(t, data, w, chi2; chi2_flux = chi2_fl,
                        C_flux = need_flux ? _chi2_flux(model_flux, data).C : NaN)
    return chi2
end

function _chi2_flat_smeared(model, x, data::Vector{<:OIdata}, w, need_flux, Tp, verb)
    nwavs = length(data)
    chi2 = zero(Tp)
    model_fluxes = need_flux ? Vector{Tp}(undef, nwavs) : Tp[]
    vmean = Vector{Any}(undef, nwavs)                  # each channel's mean V, for diffphase
    for (i, d) in enumerate(data)
        if _smeared(model, d)
            V, fl = _smeared_cvis(model, x, d, need_flux)
        else
            V = eval_model(model, x, d.uv; wl = d.uv_lam,
                           mjd = hasproperty(d, :uv_mjd) ? d.uv_mjd : nothing)
            fl = need_flux ? _flux_at(model, x, d.uv_lam,
                                      hasproperty(d, :uv_mjd) ? d.uv_mjd : nothing) : NaN
        end
        chi2 += _total_chi2(_chi2_terms(V, d, w), w)
        need_flux && (model_fluxes[i] = fl)
        vmean[i] = _rowmean(V)
    end
    need_flux && (chi2 += w[6] * Tp(_chi2_flux_polychromatic(model_fluxes, data).chi2))
    # A differential phase is a phase of a visibility, so a smeared channel gives its mean.
    w[7] > 0 && nwavs > 1 && (chi2 += w[7] * Tp(_diffphase_chi2(vmean, data, nwavs)))
    verb && printstyled(@sprintf("χ²r: %.3f\n", chi2 / _ndof(data, w)), color = :white)
    return chi2
end

function _residuals_smeared(model, x, data::OIdata, w, need_flux)
    V, model_flux = _smeared_cvis(model, x, data, need_flux)
    Vm = V isa AbstractMatrix ? V : reshape(V, :, 1)
    v2, t3, vis = _smeared_observables(Vm, data)
    r = promote_type(Float64, real(eltype(V)))[]
    w[1] > 0 && data.nv2 > 0 && append!(r, sqrt(w[1]) .* (v2 .- data.v2) ./ data.v2_err)
    w[2] > 0 && data.nt3amp > 0 && append!(r, sqrt(w[2]) .* (abs.(t3) .- data.t3amp) ./ data.t3amp_err)
    w[3] > 0 && data.nt3phi > 0 &&
        append!(r, sqrt(w[3]) .* _mod360(angle.(t3) .* (180 / π) .- data.t3phi) ./ data.t3phi_err)
    w[4] > 0 && data.nvisamp > 0 && append!(r, sqrt(w[4]) .* (abs.(vis) .- data.visamp) ./ data.visamp_err)
    w[5] > 0 && data.nvisphi > 0 &&
        append!(r, sqrt(w[5]) .* _mod360(angle.(vis) .* (180 / π) .- data.visphi) ./ data.visphi_err)
    if need_flux
        fl = _chi2_flux(model_flux, data)
        append!(r, sqrt(w[6]) .* (fill(fl.C * model_flux, data.nflux) .- data.flux) ./ data.flux_err)
    end
    return r
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


"""
    scatter_obs_cotangent!(g_cvis, V, data; g_v2, g_t3amp, g_t3phi, g_visamp, g_visphi,
                           scale_v2 = 1, …, V1 = nothing, V2 = nothing, V3 = nothing,
                           Vvis = nothing) -> g_cvis

Adjoint of the map from complex visibilities to observables: given `dL/d(observable)` for **any**
scalar `L`, accumulate `dL/dV` into `g_cvis`.

Each observable contributes `(dL/dobs) * (dobs/dV)`, scatter-added at its own `indx_*` so that
a uv point shared between observables — and the repeated legs of a closure triangle — sum
rather than overwrite.

chi² is the special case `L = chi²`, which is what [`_accumulate_g_cvis!`](@ref) supplies. It is
not the only one: variational inference needs the same operator applied to an arbitrary
direction, because the Fisher metric is the forward map linearised about a point rather than
the gradient at the data. Having one definition is the point — two would be free to disagree
about a sign or a factor, and only one of them would be under test.

`g_t3phi` and `g_visphi` are cotangents **per radian**, not per degree, even though the stored
observables are in degrees.

`data` is duck-typed rather than an `OIdata`: it needs `indx_v2` and `indx_t3_1/2/3`, and
`indx_vis` only when a vis cotangent is given. The VI extension's `ObservationConfig` carries
the same index vectors without being an `OIdata`, and requiring one there would mean either a
second copy of this scatter or a conversion that exists only to satisfy a signature.

Pass `V1`/`V2`/`V3`/`Vvis` when the caller already holds the per-leg slices; they are taken
from `V` otherwise. `scale_*` folds a per-observable weight in without materialising a scaled
copy of the cotangent, which matters because this runs once per criterion evaluation.
"""
function scatter_obs_cotangent!(g_cvis::AbstractVector{<:Complex},
                                V::AbstractVector{<:Complex},
                                data;
                                g_v2 = nothing, g_t3amp = nothing, g_t3phi = nothing,
                                g_visamp = nothing, g_visphi = nothing,
                                scale_v2 = 1, scale_t3amp = 1, scale_t3phi = 1,
                                scale_visamp = 1, scale_visphi = 1,
                                V1 = nothing, V2 = nothing, V3 = nothing, Vvis = nothing)
    # Divide-safety floor at the working precision: eps(Float64) would be nine orders of
    # magnitude below eps(Float32) and so no guard at all for a Float32 pipeline.
    ε = eps(real(eltype(V)))

    # V²  =  |V|²,      d(V²)/dV = 2 conj(V)
    if g_v2 !== nothing
        @inbounds for j in eachindex(g_v2)
            k = data.indx_v2[j]
            g_cvis[k] += scale_v2 * 2 * g_v2[j] * conj(V[k])
        end
    end

    need_t3 = g_t3amp !== nothing || g_t3phi !== nothing
    if need_t3
        z1 = V1 === nothing ? view(V, data.indx_t3_1) : V1
        z2 = V2 === nothing ? view(V, data.indx_t3_2) : V2
        z3 = V3 === nothing ? view(V, data.indx_t3_3) : V3

        # T3amp = |V1||V2||V3|,   d/dV1 = conj(V1)/|V1| * |V2| * |V3|
        if g_t3amp !== nothing
            @inbounds for j in eachindex(g_t3amp)
                a1 = abs(z1[j]); a2 = abs(z2[j]); a3 = abs(z3[j])
                c = scale_t3amp * g_t3amp[j]
                g_cvis[data.indx_t3_1[j]] += c * conj(z1[j]) / max(a1, ε) * a2 * a3
                g_cvis[data.indx_t3_2[j]] += c * conj(z2[j]) / max(a2, ε) * a1 * a3
                g_cvis[data.indx_t3_3[j]] += c * conj(z3[j]) / max(a3, ε) * a1 * a2
            end
        end

        # T3phi = arg(V1 V2 V3),  d/dV1 = -im conj(V1)/|V1|²  in the g_cvis convention above
        if g_t3phi !== nothing
            @inbounds for j in eachindex(g_t3phi)
                c = -im * scale_t3phi * g_t3phi[j]
                g_cvis[data.indx_t3_1[j]] += c * conj(z1[j]) / max(abs2(z1[j]), ε)
                g_cvis[data.indx_t3_2[j]] += c * conj(z2[j]) / max(abs2(z2[j]), ε)
                g_cvis[data.indx_t3_3[j]] += c * conj(z3[j]) / max(abs2(z3[j]), ε)
            end
        end
    end

    if g_visamp !== nothing || g_visphi !== nothing
        zv = Vvis === nothing ? view(V, data.indx_vis) : Vvis

        # visamp = |V|,   d/dV = conj(V)/|V|
        if g_visamp !== nothing
            @inbounds for j in eachindex(g_visamp)
                k = data.indx_vis[j]
                g_cvis[k] += scale_visamp * g_visamp[j] * conj(zv[j]) / max(abs(zv[j]), ε)
            end
        end

        # visphi = arg(V), d/dV = -im conj(V)/|V|²
        if g_visphi !== nothing
            @inbounds for j in eachindex(g_visphi)
                k = data.indx_vis[j]
                g_cvis[k] += -im * scale_visphi * g_visphi[j] * conj(zv[j]) / max(abs2(zv[j]), ε)
            end
        end
    end

    return g_cvis
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
    # The chi2 cotangents. `_chi2_terms` returns each `r_*` already equal to d(chi2)/d(obs),
    # except V², whose residual is defined without the factor 2 that differentiating a square
    # introduces — hence `2 * w[1]` there and a bare weight everywhere else. Phase residuals
    # carry the degree-to-radian conversion, which is why the scatter takes radians.
    scatter_obs_cotangent!(g_cvis, V, data;
        g_v2      = (w[1] > 0 && data.nv2 > 0)     ? t.r_v2     : nothing, scale_v2     = 2 * w[1],
        g_t3amp   = (w[2] > 0 && data.nt3amp > 0)  ? t.r_t3amp  : nothing, scale_t3amp  = w[2],
        g_t3phi   = (w[3] > 0 && data.nt3phi > 0)  ? t.r_t3phi  : nothing, scale_t3phi  = w[3],
        g_visamp  = (w[4] > 0 && data.nvisamp > 0) ? t.r_visamp : nothing, scale_visamp = w[4],
        g_visphi  = (w[5] > 0 && data.nvisphi > 0) ? t.r_visphi : nothing, scale_visphi = w[5],
        V1 = t.V1, V2 = t.V2, V3 = t.V3, Vvis = t.Vvis)
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
    cvis_to_chi2_f(V, data; weights=[1,1,1,0,0,0,0], model_flux=NaN) -> chi2

Chi-squared from complex visibilities `V` and interferometric data `data`, without the
gradient.

Handles the same per-channel observables as [`cvis_to_chi2_fg`](@ref) — V², T3amp, T3phi,
visamp, visphi — and additionally the **flux** term: pass the model's total flux
([`model_flux`](@ref)) as `model_flux` and give `weights[6] > 0`. (`cvis_to_chi2_fg` cannot do
this, because the flux *gradient* needs ∂F/∂x, which is not derivable from `V` alone.)

Precision follows the inputs: `ComplexF32` visibilities against an `OIdata{Float32}` return a
`Float32` chi². Accumulation is via `sum`, which is pairwise, so Float32 holds ~1e-8 relative
accuracy even for millions of residuals.
"""
function cvis_to_chi2_f(V::AbstractVector{<:Complex}, data::OIdata;
                        weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                        model_flux::Real = NaN)
    T = _chi2_eltype(V, data)
    w = _pad_weights(weights, T)
    t = _chi2_terms(V, data, w)
    chi2_fl = zero(T)
    if w[6] > 0 && data.nflux > 0 && isfinite(model_flux)
        chi2_fl = T(_chi2_flux(model_flux, data).chi2)
    end
    return _total_chi2(t, w; chi2_flux=chi2_fl)
end

"""
    cvis_to_chi2_fg(V, data; weights=[1,1,1,0,0,0,0]) -> (chi2, g_cvis)

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
                         weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
    T = _chi2_eltype(V, data)
    w = _pad_weights(weights, T)
    t = _chi2_terms(V, data, w)
    chi2 = _total_chi2(t, w)
    g_cvis = zeros(Complex{T}, length(V))
    _accumulate_g_cvis!(g_cvis, V, t, data, w)
    return chi2, g_cvis
end


# `cvis_to_chi2_fg` as loops into a caller's `g_cvis`, allocating nothing: the same terms as
# `_chi2_terms` and the same cotangents as `scatter_obs_cotangent!`, accumulated point by point
# rather than summed pairwise, so the two agree to rounding. `w` holds the seven weights in
# the working precision `T`; the χ² sums run in at least Float64, since a sequential Float32
# sum over thousands of terms would lose the accuracy `sum`'s pairwise order keeps. Model fitting calls this once per gradient, from many threads at
# once in the companion search, where the garbage collector is what stops them scaling.
function _cvis_to_chi2_fg!(g_cvis::AbstractVector{Complex{T}}, V::AbstractVector{<:Complex},
                           data::OIdata, w::NTuple{7,T}) where {T}
    fill!(g_cvis, zero(Complex{T}))
    ε = eps(real(eltype(V)))
    d2r = T(180 / π); r2 = T(360 / π)
    wrap(x) = mod(mod(x + T(180), T(360)) + T(360), T(360)) - T(180)
    A = promote_type(T, Float64)
    c2v2 = zero(A); c2ta = zero(A); c2tp = zero(A); c2va = zero(A); c2vp = zero(A)
    @inbounds begin
        if w[1] > 0 && data.nv2 > 0
            for j in eachindex(data.indx_v2)
                k = data.indx_v2[j]
                m = T(abs2(V[k])) - T(data.v2[j]); e2 = T(data.v2_err[j])^2
                c2v2 += m * m / e2
                g_cvis[k] += 2 * w[1] * 2 * (m / e2) * conj(V[k])
            end
        end
        doa = w[2] > 0 && data.nt3amp > 0
        dop = w[3] > 0 && data.nt3phi > 0
        if doa || dop
            for j in eachindex(data.indx_t3_1)
                k1 = data.indx_t3_1[j]; k2 = data.indx_t3_2[j]; k3 = data.indx_t3_3[j]
                z1 = V[k1]; z2 = V[k2]; z3 = V[k3]
                t3 = z1 * z2 * z3
                if doa
                    m = T(abs(t3)) - T(data.t3amp[j]); e2 = T(data.t3amp_err[j])^2
                    c2ta += m * m / e2
                    a1 = abs(z1); a2 = abs(z2); a3 = abs(z3)
                    c = w[2] * (2 * m / e2)
                    g_cvis[k1] += c * conj(z1) / max(a1, ε) * a2 * a3
                    g_cvis[k2] += c * conj(z2) / max(a2, ε) * a1 * a3
                    g_cvis[k3] += c * conj(z3) / max(a3, ε) * a1 * a2
                end
                if dop
                    dφ = wrap(T(angle(t3)) * d2r - T(data.t3phi[j])); e2 = T(data.t3phi_err[j])^2
                    c2tp += dφ * dφ / e2
                    c = -im * (w[3] * (r2 * dφ / e2))
                    g_cvis[k1] += c * conj(z1) / max(abs2(z1), ε)
                    g_cvis[k2] += c * conj(z2) / max(abs2(z2), ε)
                    g_cvis[k3] += c * conj(z3) / max(abs2(z3), ε)
                end
            end
        end
        doa = w[4] > 0 && data.nvisamp > 0
        dop = w[5] > 0 && data.nvisphi > 0
        if doa || dop
            for j in eachindex(data.indx_vis)
                k = data.indx_vis[j]; z = V[k]
                if doa
                    m = T(abs(z)) - T(data.visamp[j]); e2 = T(data.visamp_err[j])^2
                    c2va += m * m / e2
                    g_cvis[k] += w[4] * (2 * m / e2) * conj(z) / max(abs(z), ε)
                end
                if dop
                    dφ = wrap(T(angle(z)) * d2r - T(data.visphi[j])); e2 = T(data.visphi_err[j])^2
                    c2vp += dφ * dφ / e2
                    g_cvis[k] += -im * (w[5] * (r2 * dφ / e2)) * conj(z) / max(abs2(z), ε)
                end
            end
        end
    end
    return T(w[1] * c2v2 + w[2] * c2ta + w[3] * c2tp + w[4] * c2va + w[5] * c2vp)
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
    chi2_flat(model, x, data::OIdata; weights=[1,1,1,0,0,0,0], verb=false)
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
                   verb::Bool = false)

    _smeared(model, data) && return _chi2_flat_smeared(model, x, data, weights, verb)
    w = _pad_weights(weights)
    need_flux = w[6] > 0 && data.nflux > 0

    # Evaluate the model, and its total flux if the flux is fitted.
    # Pass per-UV wavelength and MJD so chromatic expressions ($WL, $MJD) work.
    _wl  = hasproperty(data, :uv_lam) ? data.uv_lam : nothing
    _mjd = hasproperty(data, :uv_mjd) ? data.uv_mjd : nothing
    V = eval_model(model, x, data.uv; wl=_wl, mjd=_mjd)
    model_flux = need_flux ? _flux_at(model, x, _wl, _mjd) : NaN

    chi2 = cvis_to_chi2_f(V, data; weights=w, model_flux=model_flux)

    if verb
        t = _chi2_terms(V, data, w)
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
    chi2_flat_fg(model, x, data::OIdata; weights=[1,1,1,0,0,0,0], verb=false,
                 jac_cache=nothing) -> (T, Vector{T})

Evaluate chi-squared and its gradient w.r.t. `x`.

The gradient is computed analytically via the Wirtinger chain rule applied to
the complex visibility Jacobian J = ∂V/∂x from `eval_model_grad`.

Returns `(chi2, grad)`; the element type follows `x` and `data.uv`, so Float32 inputs give a Float32 gradient.
A [`ModelJacCache`](@ref) as `jac_cache` keeps the Jacobian's buffers between calls.
"""
function chi2_flat_fg(model::FlatModel,
                      x::AbstractVector,
                      data::OIdata;
                      weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                      verb::Bool = false,
                      jac_cache::Union{Nothing,ModelJacCache} = nothing)

    _smeared(model, data) && return _smeared_chi2_fg(model, x, data, _pad_weights(weights), jac_cache)
    w = _pad_weights(weights)
    w_flux = w[6]
    need_flux = w_flux > 0 && data.nflux > 0

    # ── Forward pass + Jacobian ─────────────────────────────────────────────
    _wl  = hasproperty(data, :uv_lam) ? data.uv_lam : nothing
    _mjd = hasproperty(data, :uv_mjd) ? data.uv_mjd : nothing
    V, J = _model_jac(jac_cache, model, x, data.uv; wl=_wl, mjd=_mjd)
    model_flux, J_zero = need_flux ? _flux_grad_at(model, x, _wl, _mjd) : (NaN, nothing)

    # ── Per-channel chi2 + adjoint source ──────────────────────────────────
    T = _chi2_eltype(V, data)
    g_cvis = Vector{Complex{T}}(undef, length(V))
    chi2_obs = _cvis_to_chi2_fg!(g_cvis, V, data, ntuple(i -> T(w[i]), Val(7)))

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
        t = _chi2_terms(V, data, w)
        _verb_print(t, data, w, chi2; chi2_flux=chi2_fl, C_flux=C_fl)
    end
    return chi2, g
end


# ─────────────────────────────────────────────────────────────────────────────
# chi2_flat  (forward only — polychromatic)
# ─────────────────────────────────────────────────────────────────────────────

"""
    chi2_flat(model, x, data::Vector{<:OIdata}; weights=[1,1,1,0,0,0,0], verb=false)
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
                   verb::Bool = false)

    nwavs = length(data)
    # Precision follows the caller, as in the monochromatic path.
    Tp = promote_type(float(eltype(x)), eltype(first(data).uv))
    w = _pad_weights(weights, Tp)

    # Evaluate model at each channel's UV points
    need_flux = w[6] > 0 && any(d.nflux > 0 for d in data)
    any(d -> _smeared(model, d), data) &&
        return _chi2_flat_smeared(model, x, data, w, need_flux, Tp, verb)
    if need_flux
        cvis = Vector{Vector{Complex{Tp}}}(undef, nwavs)
        model_fluxes = Vector{Tp}(undef, nwavs)
        for i in 1:nwavs
            _wl_i  = hasproperty(data[i], :uv_lam) ? data[i].uv_lam : nothing
            _mjd_i = hasproperty(data[i], :uv_mjd) ? data[i].uv_mjd : nothing
            model_fluxes[i] = _flux_at(model, x, _wl_i, _mjd_i)
            cvis[i] = eval_model(model, x, data[i].uv; wl=_wl_i, mjd=_mjd_i)
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
        t = _chi2_terms(cvis[i], data[i], w)
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
        # A differential phase is a phase of a visibility, so a smeared channel contributes its
        # mean visibility.
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
    chi2_flat_fg(model, x, data::Vector{<:OIdata}; weights=[1,1,1,0,0,0,0], verb=false)
        -> (T, Vector{T})

Polychromatic chi-squared and gradient.  Computes per-channel V2/T3/visamp/visphi/flux
gradients plus cross-channel differential phase gradient.
"""
function chi2_flat_fg(model::FlatModel,
                      x::AbstractVector,
                      data::Vector{<:OIdata};
                      weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                      verb::Bool = false)

    if any(d -> _smeared(model, d), data)
        wp = _pad_weights(weights)
        # Channel by channel through the smeared residuals, unless a term couples the channels
        # (one flux scale for all of them, differential phase).
        if wp[6] == 0 && wp[7] == 0
            c = zero(Float64); g = zeros(length(x))
            for d in data
                cd, gd = chi2_flat_fg(model, x, d; weights = wp)
                c += cd; g .+= gd
            end
            return c, g
        end
        f = p -> chi2_flat(model, p, data; weights)
        return f(x), ForwardDiff.gradient(f, x)
    end
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
            fl_i, gf_i = _flux_grad_at(model, x, _wl_i, _mjd_i)
            model_fluxes[i] = fl_i
            J_zeros[i] = gf_i
            cvis[i], Jac[i] = eval_model_grad(model, x, data[i].uv; wl=_wl_i, mjd=_mjd_i)
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
        chi2_i, g_cvis_i = cvis_to_chi2_fg(cvis[i], data[i]; weights=w)
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
    residuals_flat(model, x, data; weights) -> Vector{Float64}

Compute the weighted residual vector `r` such that `chi2 = dot(r, r)`.

Each element is `(model_observable - data_observable) / σ`, concatenated
in order: [V2; T3amp; T3phi; visamp; visphi; flux].
"""
function residuals_flat(model::FlatModel,
                        x::AbstractVector,
                        data::OIdata;
                        weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0])

    w = _pad_weights(weights)
    need_flux = w[6] > 0 && data.nflux > 0

    _wl  = hasproperty(data, :uv_lam) ? data.uv_lam : nothing
    _mjd = hasproperty(data, :uv_mjd) ? data.uv_mjd : nothing
    _smeared(model, data) && return _residuals_smeared(model, x, data, w, need_flux)
    V = eval_model(model, x, data.uv; wl=_wl, mjd=_mjd)
    model_flux = need_flux ? _flux_at(model, x, _wl, _mjd) : NaN

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
        t3phi_model = angle.(t3) .* (180.0 / π)
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
    residuals_flat_jac(model, x, data; weights, jac_cache = nothing)
        -> (Vector{Float64}, Matrix{Float64})

Compute the weighted residual vector and its Jacobian `J[i,j] = ∂r_i/∂x_j`.

Returns `(r, J)` where `chi2 = dot(r, r)` and `∂chi2/∂x = 2 Jᵀ r`. A
[`ModelJacCache`](@ref) as `jac_cache` keeps the model Jacobian's buffers between calls.
"""
function residuals_flat_jac(model::FlatModel,
                            x::AbstractVector,
                            data::OIdata;
                            weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                            jac_cache::Union{Nothing,ModelJacCache} = nothing)

    if _smeared(model, data)
        r, J = _smeared_residuals_jac(model, x, data, _pad_weights(weights), jac_cache)
        return Float64.(r), Float64.(J)                # LsqFit consumes these
    end
    w = _pad_weights(weights)
    w_v2, w_t3amp, w_t3phi, w_visamp, w_visphi, w_flux = w[1], w[2], w[3], w[4], w[5], w[6]
    need_flux = w_flux > 0 && data.nflux > 0
    nparams = length(x)

    _wl  = hasproperty(data, :uv_lam) ? data.uv_lam : nothing
    _mjd = hasproperty(data, :uv_mjd) ? data.uv_mjd : nothing
    V, dV = _model_jac(jac_cache, model, x, data.uv; wl=_wl, mjd=_mjd)
    model_flux, dflux = need_flux ? _flux_grad_at(model, x, _wl, _mjd) : (NaN, nothing)

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
        t3phi_model = angle.(t3) .* (180.0 / π)
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
        visphi_model = angle.(Vvis) .* (180.0 / π)
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
        # J_flux[i,j] = sqrt(w) * C * ∂F/∂x_j / σ_i
        J_zero = permutedims(dflux)  # 1 × nparams
        J_fl = sqrt(w_flux) .* fl.C .* repeat(J_zero, data.nflux, 1) ./ data.flux_err
        push!(r_parts, r_fl); push!(J_parts, J_fl)
    end

    r = vcat(r_parts...)
    J = vcat(J_parts...)
    return r, J
end


# ─────────────────────────────────────────────────────────────────────────────
# residuals_flat! / residuals_flat_jac! — the same rows, written into the caller's buffers
#
# Element by element the arithmetic of `residuals_flat` / `residuals_flat_jac`, in the same
# order, so the values are identical; without the copies of `V[idx]` and `dV[idx, :]`, the
# temporaries and the final `vcat`. A smeared model takes the allocating path and is copied in.
# ─────────────────────────────────────────────────────────────────────────────

# The rows `residuals_flat` returns for `data` under the padded weights `w`: an observable
# counts when weighted and measured at all, with one row per entry of its index, as
# `residuals_flat` loops over them.
_residual_nrows(data, w, need_flux) =
    (w[1] > 0 && data.nv2 > 0 ? length(data.indx_v2) : 0) +
    (w[2] > 0 && data.nt3amp > 0 ? length(data.indx_t3_1) : 0) +
    (w[3] > 0 && data.nt3phi > 0 ? length(data.indx_t3_1) : 0) +
    (w[4] > 0 && data.nvisamp > 0 ? length(data.indx_vis) : 0) +
    (w[5] > 0 && data.nvisphi > 0 ? length(data.indx_vis) : 0) +
    (need_flux ? data.nflux : 0)

# The rows of the VISPHI block among a dataset's residuals (empty when not weighted).
function _visphi_rows(data, w)
    k = (w[1] > 0 && data.nv2 > 0 ? length(data.indx_v2) : 0) +
        (w[2] > 0 && data.nt3amp > 0 ? length(data.indx_t3_1) : 0) +
        (w[3] > 0 && data.nt3phi > 0 ? length(data.indx_t3_1) : 0) +
        (w[4] > 0 && data.nvisamp > 0 ? length(data.indx_vis) : 0)
    return w[5] > 0 && data.nvisphi > 0 ? ((k + 1):(k + length(data.indx_vis))) : (1:0)
end

# One wrapped phase residual, as `_mod360` computes it for a vector.
@inline function _mod360_1(x)
    T = float(typeof(x))
    return mod(mod(x + T(180), T(360)) + T(360), T(360)) - T(180)
end

"""
    residuals_flat!(r, model, x, data; weights) -> r

[`residuals_flat`](@ref) written into `r`, whose length is the number of rows it returns.
"""
function residuals_flat!(r::AbstractVector, model::FlatModel, x::AbstractVector, data::OIdata;
                         weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                         ref = nothing)
    if _smooths(ref)
        _smoothed_rows!(r, nothing, model, x, data, _pad_weights(weights), ref, nothing)
    else
        _residuals_flat!(r, model, x, data, weights)
    end
    ref === nothing || _apply_reference!(r, nothing, ref, model, x, data, _pad_weights(weights))
    return r
end

# The VISPHI rows of `r` (and `J`) referenced by the measurement model's context `ref`.
function _apply_reference!(r, J, ref, model, x, data, w)
    ref.mm.visphi === :continuum || return 0
    rows = _visphi_rows(data, w)
    isempty(rows) && return 0
    vals = model.resolver.fn(_plain.(x))
    return _reference_phases!(view(r, rows), J === nothing ? nothing : view(J, rows, :), ref, vals,
                              data.visphi_err, sqrt(w[5]), data.visphi)
end

function _residuals_flat!(r::AbstractVector, model::FlatModel, x::AbstractVector, data::OIdata,
                          weights::AbstractVector{<:Real})
    w = _pad_weights(weights)
    need_flux = w[6] > 0 && data.nflux > 0
    _wl  = hasproperty(data, :uv_lam) ? data.uv_lam : nothing
    _mjd = hasproperty(data, :uv_mjd) ? data.uv_mjd : nothing
    n = _residual_nrows(data, w, need_flux)
    length(r) == n || throw(DimensionMismatch("residuals_flat!: $n rows for a buffer of $(length(r))"))
    if _smeared(model, data)
        copyto!(r, _residuals_smeared(model, x, data, w, need_flux))
        return r
    end
    V = eval_model(model, x, data.uv; wl=_wl, mjd=_mjd)
    model_flux = need_flux ? _flux_at(model, x, _wl, _mjd) : NaN
    return _residual_rows!(r, V, data, w, need_flux, model_flux)
end

# The rows of `residuals_flat!` from the model's visibilities: a function of its own, so the
# loops are compiled for `V`'s concrete type.
function _residual_rows!(r, V, data, w, need_flux, model_flux)
    k = 0
    @inbounds begin
        if w[1] > 0 && data.nv2 > 0
            s = sqrt(w[1]); idx = data.indx_v2
            for i in eachindex(idx)
                r[k += 1] = s * (abs2(V[idx[i]]) - data.v2[i]) / data.v2_err[i]
            end
        end
        if w[2] > 0 && data.nt3amp > 0
            s = sqrt(w[2]); i1, i2, i3 = data.indx_t3_1, data.indx_t3_2, data.indx_t3_3
            for i in eachindex(i1)
                t3 = V[i1[i]] * V[i2[i]] * V[i3[i]]
                r[k += 1] = s * (abs(t3) - data.t3amp[i]) / data.t3amp_err[i]
            end
        end
        if w[3] > 0 && data.nt3phi > 0
            s = sqrt(w[3]); i1, i2, i3 = data.indx_t3_1, data.indx_t3_2, data.indx_t3_3
            for i in eachindex(i1)
                t3 = V[i1[i]] * V[i2[i]] * V[i3[i]]
                r[k += 1] = s * _mod360_1(angle(t3) * (180.0 / π) - data.t3phi[i]) / data.t3phi_err[i]
            end
        end
        if w[4] > 0 && data.nvisamp > 0
            s = sqrt(w[4]); idx = data.indx_vis
            for i in eachindex(idx)
                r[k += 1] = s * (abs(V[idx[i]]) - data.visamp[i]) / data.visamp_err[i]
            end
        end
        if w[5] > 0 && data.nvisphi > 0
            s = sqrt(w[5]); idx = data.indx_vis
            for i in eachindex(idx)
                r[k += 1] = s * _mod360_1(angle(V[idx[i]]) * (180.0 / π) - data.visphi[i]) / data.visphi_err[i]
            end
        end
    end
    if need_flux
        fl = _chi2_flux(model_flux, data)
        s = sqrt(w[6]); fm = fl.C * model_flux
        @inbounds for i in 1:data.nflux
            r[k += 1] = s * (fm - data.flux[i]) / data.flux_err[i]
        end
    end
    return r
end

"""
    residuals_flat_jac!(r, J, model, x, data; weights, jac_cache = nothing) -> (r, J)

[`residuals_flat_jac`](@ref) written into `r` and `J` (rows × parameters).
"""
function residuals_flat_jac!(r::AbstractVector, J::AbstractMatrix, model::FlatModel,
                             x::AbstractVector, data::OIdata;
                             weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                             jac_cache::Union{Nothing,ModelJacCache} = nothing, ref = nothing)
    if _smooths(ref)
        _smoothed_rows!(r, J, model, x, data, _pad_weights(weights), ref, jac_cache)
    else
        _residuals_flat_jac!(r, J, model, x, data, weights, jac_cache)
    end
    ref === nothing || _apply_reference!(r, J, ref, model, x, data, _pad_weights(weights))
    return r, J
end

function _residuals_flat_jac!(r::AbstractVector, J::AbstractMatrix, model::FlatModel,
                              x::AbstractVector, data::OIdata, weights, jac_cache)
    w = _pad_weights(weights)
    need_flux = w[6] > 0 && data.nflux > 0
    np = length(x)
    n = _residual_nrows(data, w, need_flux)
    (length(r) == n && size(J) == (n, np)) ||
        throw(DimensionMismatch("residuals_flat_jac!: $n rows × $np parameters for buffers " *
                                "$(length(r)), $(size(J))"))
    if _smeared(model, data)
        r0, J0 = residuals_flat_jac(model, x, data; weights, jac_cache)
        copyto!(r, r0); copyto!(J, J0)
        return r, J
    end
    _wl  = hasproperty(data, :uv_lam) ? data.uv_lam : nothing
    _mjd = hasproperty(data, :uv_mjd) ? data.uv_mjd : nothing
    V, dV = _model_jac(jac_cache, model, x, data.uv; wl=_wl, mjd=_mjd)
    model_flux, dflux = need_flux ? _flux_grad_at(model, x, _wl, _mjd) : (NaN, nothing)
    return _residual_jac_rows!(r, J, V, dV, data, w, need_flux, model_flux, dflux, np)
end

# The rows of `residuals_flat_jac!` from the model's visibilities and their Jacobian (`V`, `dV`
# come out of an untyped cache): a function of its own, compiled for their concrete types.
function _residual_jac_rows!(r, J, V, dV, data, w, need_flux, model_flux, dflux, np)
    ε = eps(Float64)
    k = 0
    @inbounds begin
        if w[1] > 0 && data.nv2 > 0
            s = sqrt(w[1]); idx = data.indx_v2
            for i in eachindex(idx)
                k += 1; Vi = V[idx[i]]; e = data.v2_err[i]
                r[k] = s * (abs2(Vi) - data.v2[i]) / e
                for j in 1:np
                    J[k, j] = s * 2.0 * real(conj(Vi) * dV[idx[i], j]) / e
                end
            end
        end
        if w[2] > 0 && data.nt3amp > 0
            s = sqrt(w[2]); i1, i2, i3 = data.indx_t3_1, data.indx_t3_2, data.indx_t3_3
            for i in eachindex(i1)
                k += 1; V1 = V[i1[i]]; V2c = V[i2[i]]; V3 = V[i3[i]]; e = data.t3amp_err[i]
                t3 = V1 * V2c * V3
                a = abs(t3)
                r[k] = s * (a - data.t3amp[i]) / e
                pf = conj(t3) / max(a, ε)
                c1, c2, c3 = V2c * V3, V1 * V3, V1 * V2c
                for j in 1:np
                    dT3 = c1 * dV[i1[i], j] + c2 * dV[i2[i], j] + c3 * dV[i3[i], j]
                    J[k, j] = s * real(pf * dT3) / e
                end
            end
        end
        if w[3] > 0 && data.nt3phi > 0
            s = sqrt(w[3]); i1, i2, i3 = data.indx_t3_1, data.indx_t3_2, data.indx_t3_3
            for i in eachindex(i1)
                k += 1; V1 = V[i1[i]]; V2c = V[i2[i]]; V3 = V[i3[i]]; e = data.t3phi_err[i]
                t3 = V1 * V2c * V3
                r[k] = s * _mod360_1(angle(t3) * (180.0 / π) - data.t3phi[i]) / e
                pf = conj(t3) / max(abs2(t3), ε)
                c1, c2, c3 = V2c * V3, V1 * V3, V1 * V2c
                for j in 1:np
                    dT3 = c1 * dV[i1[i], j] + c2 * dV[i2[i], j] + c3 * dV[i3[i], j]
                    J[k, j] = s * (180.0 / π) * imag(pf * dT3) / e
                end
            end
        end
        if w[4] > 0 && data.nvisamp > 0
            s = sqrt(w[4]); idx = data.indx_vis
            for i in eachindex(idx)
                k += 1; Vi = V[idx[i]]; e = data.visamp_err[i]
                a = abs(Vi)
                r[k] = s * (a - data.visamp[i]) / e
                pf = conj(Vi) / max(a, ε)
                for j in 1:np
                    J[k, j] = s * real(pf * dV[idx[i], j]) / e
                end
            end
        end
        if w[5] > 0 && data.nvisphi > 0
            s = sqrt(w[5]); idx = data.indx_vis
            for i in eachindex(idx)
                k += 1; Vi = V[idx[i]]; e = data.visphi_err[i]
                r[k] = s * _mod360_1(angle(Vi) * (180.0 / π) - data.visphi[i]) / e
                pf = conj(Vi) / max(abs2(Vi), ε)
                for j in 1:np
                    J[k, j] = s * (180.0 / π) * imag(pf * dV[idx[i], j]) / e
                end
            end
        end
    end
    if need_flux
        fl = _chi2_flux(model_flux, data)
        s = sqrt(w[6]); fm = fl.C * model_flux
        @inbounds for i in 1:data.nflux
            k += 1; e = data.flux_err[i]
            r[k] = s * (fm - data.flux[i]) / e
            for j in 1:np
                J[k, j] = s * fl.C * dflux[j] / e
            end
        end
    end
    return r, J
end

# The sub-wavelength points of a smeared evaluation, and the adjoint buffer that goes with them:
# kept in the cache, since they depend only on the data and the number of samples.
function _smeared_points(::Nothing, data, ns, mjd, ::Type{C}) where {C}
    return (_sub_points(data.uv, data.uv_lam, data.uv_dlam, mjd, ns)...,
            Vector{C}(undef, size(data.uv, 2) * ns))
end
function _smeared_points(c::ModelJacCache, data, ns, mjd, ::Type{C}) where {C}
    sp = c.sub
    if !(sp isa Tuple && sp[1] === data && sp[2] == ns && sp[7] isa Vector{C})
        sp = c.sub = (data, ns, _sub_points(data.uv, data.uv_lam, data.uv_dlam, mjd, ns)...,
                      Vector{C}(undef, size(data.uv, 2) * ns))
    end
    return sp[3], sp[4], sp[5], sp[6], sp[7]
end

"""
    _smeared_chi2_fg(model, x, data, w, jac_cache = nothing) -> (chi2, grad)

`chi2_flat_fg` for a smeared model. The chain rule of `_smeared_residuals_jac`,
contracted before it is materialised: each averaged observable's cotangent goes back to every
sub-wavelength visibility it averages (`G`, one entry per sub-point), and the gradient is
`real(Jᵀ G)` over the sub-points' model Jacobian. χ² is `Σ r²` over the same residuals.
"""
function _smeared_chi2_fg(model::FlatModel, x::AbstractVector, data::OIdata, w,
                          jac_cache = nothing)
    N  = size(data.uv, 2)
    ns = smear_samples(model, data)
    mjd = hasproperty(data, :uv_mjd) ? data.uv_mjd : nothing
    # The precision `eval_model_grad` works in: the parameters' and the uv coordinates'.
    C  = Complex{promote_type(float(eltype(x)), eltype(data.uv))}
    uvs, wls, mjds, plane, G = _smeared_points(jac_cache, data, ns, mjd, C)
    Vf, Jf = _model_jac(jac_cache, model, x, uvs; wl = wls, mjd = mjds, plane_wl = plane)
    T  = real(C)
    chi2, g = _smeared_chi2_g(G, reshape(convert(Vector{C}, Vf), N, ns),
                              convert(Matrix{C}, Jf), data, ntuple(i -> T(w[i]), Val(7)),
                              model.smear_coherent)
    if w[6] > 0 && data.nflux > 0
        # As in `_smeared_residuals_jac`: the model's total flux, the scale held fixed.
        Fz, gz = _flux_grad_at(model, x, data.uv_lam, mjd)
        fl = _chi2_flux(Fz, data)
        sw = sqrt(w[6])
        r  = sw .* (fl.C * Fz .- data.flux) ./ data.flux_err
        chi2 += sum(abs2, r)
        g .+= 2 * sum(r .* sw .* fl.C ./ data.flux_err) .* gz
    end
    return chi2, g
end

# The averaged observables' χ² and the cotangent of every sub-point visibility, into `G`,
# then the gradient. Typed by `T` behind this barrier, since the cache's buffers are not.
function _smeared_chi2_g(G::Vector{Complex{T}}, V::Matrix{Complex{T}}, Jf::Matrix{Complex{T}},
                         data::OIdata, w::NTuple{7,T}, coherent::Bool) where {T}
    N, ns = size(V)
    fill!(G, zero(Complex{T}))
    Gm = reshape(G, N, ns)
    ε = eps(T); d2r = T(180 / π); r2 = T(360 / π); inv_ns = inv(T(ns))
    wrap(x) = mod(mod(x + T(180), T(360)) + T(360), T(360)) - T(180)
    vc(q) = (s = zero(Complex{T}); @inbounds(for k in 1:ns; s += V[q, k]; end); s * inv_ns)
    A = promote_type(T, Float64)
    chi2 = zero(A)
    @inbounds begin
        if w[1] > 0 && data.nv2 > 0
            for (j, q) in enumerate(data.indx_v2)
                v2 = zero(T); z = zero(Complex{T})
                if coherent
                    z = vc(q); v2 = abs2(z)
                else
                    for k in 1:ns; v2 += abs2(V[q, k]); end
                    v2 *= inv_ns
                end
                m = v2 - T(data.v2[j]); e2 = T(data.v2_err[j])^2
                chi2 += w[1] * m * m / e2
                c = w[1] * 2 * m / e2 * 2 * inv_ns
                for k in 1:ns
                    Gm[q, k] += c * conj(coherent ? z : V[q, k])
                end
            end
        end
        doa = w[2] > 0 && data.nt3amp > 0
        dop = w[3] > 0 && data.nt3phi > 0
        if doa || dop
            for j in eachindex(data.indx_t3_1)
                q1 = data.indx_t3_1[j]; q2 = data.indx_t3_2[j]; q3 = data.indx_t3_3[j]
                a = b = c3 = zero(Complex{T})
                if coherent
                    a = vc(q1); b = vc(q2); c3 = vc(q3); t3 = a * b * c3
                else
                    t3 = zero(Complex{T})
                    for k in 1:ns; t3 += V[q1, k] * V[q2, k] * V[q3, k]; end
                    t3 *= inv_ns
                end
                γ = zero(Complex{T})
                if doa
                    amp = abs(t3)
                    m = amp - T(data.t3amp[j]); e2 = T(data.t3amp_err[j])^2
                    chi2 += w[2] * m * m / e2
                    γ += w[2] * 2 * m / e2 * conj(t3) / max(amp, ε)
                end
                if dop
                    dφ = wrap(angle(t3) * d2r - T(data.t3phi[j])); e2 = T(data.t3phi_err[j])^2
                    chi2 += w[3] * dφ * dφ / e2
                    γ += -im * (w[3] * r2 * dφ / e2) * conj(t3) / max(abs2(t3), ε)
                end
                γ *= inv_ns
                for k in 1:ns
                    v1, v2, v3 = coherent ? (a, b, c3) : (V[q1, k], V[q2, k], V[q3, k])
                    Gm[q1, k] += γ * v2 * v3
                    Gm[q2, k] += γ * v1 * v3
                    Gm[q3, k] += γ * v1 * v2
                end
            end
        end
        doa = w[4] > 0 && data.nvisamp > 0
        dop = w[5] > 0 && data.nvisphi > 0
        if doa || dop
            for (j, q) in enumerate(data.indx_vis)
                z = vc(q)
                γ = zero(Complex{T})
                if doa
                    amp = abs(z)
                    m = amp - T(data.visamp[j]); e2 = T(data.visamp_err[j])^2
                    chi2 += w[4] * m * m / e2
                    γ += w[4] * 2 * m / e2 * conj(z) / max(amp, ε)
                end
                if dop
                    dφ = wrap(angle(z) * d2r - T(data.visphi[j])); e2 = T(data.visphi_err[j])^2
                    chi2 += w[5] * dφ * dφ / e2
                    γ += -im * (w[5] * r2 * dφ / e2) * conj(z) / max(abs2(z), ε)
                end
                γ *= inv_ns
                for k in 1:ns
                    Gm[q, k] += γ
                end
            end
        end
    end
    return T(chi2), real.(transpose(Jf) * G)
end

"""
    _smeared_residuals_jac(model, x, data, w, jac_cache = nothing) -> (r, J)

`residuals_flat_jac` for a smeared model, by the chain rule rather than by differentiating the
χ²: `eval_model_grad` at every sub-wavelength point gives each `V_k` and `∂V_k/∂x`, and the
averaged observables differentiate term by term,

    ∂⟨|V_k|²⟩ = ⟨2 Re(V̄_k ∂V_k)⟩,   ∂⟨V₁V₂V₃⟩ = ⟨V₂V₃ ∂V₁ + V₁V₃ ∂V₂ + V₁V₂ ∂V₃⟩,   ∂⟨V⟩ = ⟨∂V⟩;

in the coherent mode `V_k` and `∂V_k` are averaged first. Amplitudes and phases then follow
as in the unsmeared Jacobian.
"""
function _smeared_residuals_jac(model::FlatModel, x::AbstractVector, data::OIdata, w,
                                jac_cache = nothing)
    N  = size(data.uv, 2)
    ns = smear_samples(model, data)
    mjd = hasproperty(data, :uv_mjd) ? data.uv_mjd : nothing
    uvs, wls, mjds, plane = _sub_points(data.uv, data.uv_lam, data.uv_dlam, mjd, ns)
    Vf, Jf = _model_jac(jac_cache, model, x, uvs; wl = wls, mjd = mjds, plane_wl = plane)
    p  = length(x)
    # The precision `eval_model_grad` works in: the parameters' and the uv coordinates'.
    C  = Complex{promote_type(float(eltype(x)), eltype(data.uv))}
    V  = reshape(convert(Vector{C}, Vf), N, ns)
    Jd = reshape(convert(Matrix{C}, Jf), N, ns, p)
    if model.smear_coherent
        V  = reshape(vec(sum(V; dims = 2)) ./ ns, N, 1)
        Jd = reshape(dropdims(sum(Jd; dims = 2); dims = 2) ./ ns, N, 1, p)
    end
    r, J = _smeared_obs_jac(V, Jd, data, w)
    if w[6] > 0 && data.nflux > 0
        # The flux is the model's total flux, which smearing does not change. For an
        # uncalibrated flux the scale C is at its optimum, so holding it fixed gives the exact
        # gradient of the χ² minimised over C.
        Fz, gz = _flux_grad_at(model, x, data.uv_lam, mjd)
        fl = _chi2_flux(Fz, data)
        r = vcat(r, sqrt(w[6]) .* (fill(fl.C * Fz, data.nflux) .- data.flux) ./
                    data.flux_err)
        J = vcat(J, sqrt(w[6]) .* fl.C .* repeat(permutedims(gz), data.nflux, 1) ./
                    data.flux_err)
    end
    return r, J
end

# The averaged observables' residuals and Jacobian from concretely typed visibilities: the
# loops run behind this function barrier, since what `eval_model_grad` returns is not inferred.
function _smeared_obs_jac(V::Matrix{Complex{T}}, Jd::Array{Complex{T},3}, data::OIdata,
                          w) where {T}
    N, ns, p = size(Jd)
    ε = eps(T)
    r_parts = Vector{T}[]
    J_parts = Matrix{T}[]

    if w[1] > 0 && data.nv2 > 0
        n = length(data.indx_v2)
        v2 = zeros(T, n); dv2 = zeros(T, n, p)
        @inbounds for (j, q) in enumerate(data.indx_v2), k in 1:ns
            v2[j] += abs2(V[q, k]) / ns
            for l in 1:p
                dv2[j, l] += 2real(conj(V[q, k]) * Jd[q, k, l]) / ns
            end
        end
        push!(r_parts, T(sqrt(w[1])) .* (v2 .- data.v2) ./ data.v2_err)
        push!(J_parts, T(sqrt(w[1])) .* dv2 ./ data.v2_err)
    end

    if (w[2] > 0 && data.nt3amp > 0) || (w[3] > 0 && data.nt3phi > 0)
        n = length(data.indx_t3_1)
        t3 = zeros(Complex{T}, n); dt3 = zeros(Complex{T}, n, p)
        @inbounds for j in 1:n, k in 1:ns
            q1, q2, q3 = data.indx_t3_1[j], data.indx_t3_2[j], data.indx_t3_3[j]
            a, b, c = V[q1, k], V[q2, k], V[q3, k]
            t3[j] += a * b * c / ns
            for l in 1:p
                dt3[j, l] += (b * c * Jd[q1, k, l] + a * c * Jd[q2, k, l] +
                              a * b * Jd[q3, k, l]) / ns
            end
        end
        if w[2] > 0 && data.nt3amp > 0
            A = abs.(t3)
            push!(r_parts, T(sqrt(w[2])) .* (A .- data.t3amp) ./ data.t3amp_err)
            push!(J_parts, T(sqrt(w[2])) .* real.(conj.(t3) ./ max.(A, ε) .* dt3) ./ data.t3amp_err)
        end
        if w[3] > 0 && data.nt3phi > 0
            φ = angle.(t3) .* T(180 / π)
            push!(r_parts, T(sqrt(w[3])) .* _mod360(φ .- data.t3phi) ./ data.t3phi_err)
            push!(J_parts, T(sqrt(w[3])) .* T(180 / π) .*
                           imag.(conj.(t3) ./ max.(abs2.(t3), ε) .* dt3) ./ data.t3phi_err)
        end
    end

    if (w[4] > 0 && data.nvisamp > 0) || (w[5] > 0 && data.nvisphi > 0)
        n = length(data.indx_vis)
        vv = zeros(Complex{T}, n); dvv = zeros(Complex{T}, n, p)
        @inbounds for (j, q) in enumerate(data.indx_vis), k in 1:ns
            vv[j] += V[q, k] / ns
            for l in 1:p
                dvv[j, l] += Jd[q, k, l] / ns
            end
        end
        if w[4] > 0 && data.nvisamp > 0
            A = abs.(vv)
            push!(r_parts, T(sqrt(w[4])) .* (A .- data.visamp) ./ data.visamp_err)
            push!(J_parts, T(sqrt(w[4])) .* real.(conj.(vv) ./ max.(A, ε) .* dvv) ./ data.visamp_err)
        end
        if w[5] > 0 && data.nvisphi > 0
            φ = angle.(vv) .* T(180 / π)
            push!(r_parts, T(sqrt(w[5])) .* _mod360(φ .- data.visphi) ./ data.visphi_err)
            push!(J_parts, T(sqrt(w[5])) .* T(180 / π) .*
                           imag.(conj.(vv) ./ max.(abs2.(vv), ε) .* dvv) ./ data.visphi_err)
        end
    end

    return reduce(vcat, r_parts; init = T[]),
           isempty(J_parts) ? zeros(T, 0, p) : reduce(vcat, J_parts)
end


# ─────────────────────────────────────────────────────────────────────────────
# chi2_flat_reduced  (chi2 / ndof — for display)
# ─────────────────────────────────────────────────────────────────────────────

"""
    chi2_flat_reduced(model, x, data; weights=[1,1,1,0,0,0,0])
        -> T

Return chi2 / ndof (reduced chi-squared).  Useful for assessing goodness of fit.
Works with both monochromatic (OIdata) and polychromatic (Vector{<:OIdata}) data.
"""
function chi2_flat_reduced(model::FlatModel,
                           x::AbstractVector,
                           data::Union{OIdata, Vector{<:OIdata}};
                           weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
    chi2 = chi2_flat(model, x, data; weights)
    return chi2 / _ndof(data, weights)
end


# ─────────────────────────────────────────────────────────────────────────────
# chi2_flat_print  (verbose summary without running a separate forward pass)
# ─────────────────────────────────────────────────────────────────────────────

"""
    chi2_flat_print(model, x, data; weights=[1,1,1,0,0,0,0])

Print a per-observable chi2 breakdown.  Does a single forward pass.
"""
function chi2_flat_print(model::FlatModel,
                         x::AbstractVector,
                         data::Union{OIdata, Vector{<:OIdata}};
                         weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
    chi2_flat(model, x, data; weights, verb=true)
    nothing
end

"""
    model_to_chi2(model, x, data; weights=[1,1,1,0,0,0,0], verb=false)

Alias for `chi2_flat`.  Compute the weighted chi-squared of a parametric
model against data, for consistency with the `model_to_obs` / `model_to_vis` family.
"""
const model_to_chi2 = chi2_flat

"""
    model_to_chi2_fg(model, x, data; weights=[1,1,1,0,0,0,0], verb=false)

Alias for `chi2_flat_fg`.  Compute chi-squared and its gradient w.r.t. `x`,
for consistency with the `model_to_obs` / `model_to_vis` family.
"""
const model_to_chi2_fg = chi2_flat_fg
# ─────────────────────────────────────────────────────────────────────────────
# cvis_to_chi2_noalloc — append this to src/chi2_flat.jl.
#
# WHY THIS EXISTS ALONGSIDE `cvis_to_chi2_f`
#
# `cvis_to_chi2_f` already gives chi2 without the gradient, and is the right
# choice almost everywhere — it also handles the flux term and does proper
# `_chi2_eltype` precision handling.  But it routes through `_chi2_terms`, which
# allocates on nearly every line:
#
#     v2_model    = abs2.(V[data.indx_v2])        # gather + broadcast temporary
#     V1          = V[data.indx_t3_1]             # gather
#     t3amp_model = abs.(t3)                      # temporary
#     dphi        = _mod360(t3phi_model .- ...)   # two more
#
# That is ~10 heap allocations per call, several of them nuv-length arrays.  The
# SQUEEZE sampler calls this ONCE PER PROPOSAL — millions of times per run — so
# the allocation and GC pressure dominate the arithmetic.  Measured on the
# original OITOOLS: 0 bytes and ~9x faster than `cvis_to_chi2_fg`.
#
# So: use `cvis_to_chi2_f` for ordinary work; use this in inner sampler loops.
# The two must agree to round-off, which `test_squeeze.jl` checks.
#
# The chi2 accumulator is Float64 regardless of the data precision `T`: it is a
# sum of ~ndata positive terms, so Float64 costs one scalar and removes
# summation error from the Metropolis acceptance test.
# ─────────────────────────────────────────────────────────────────────────────

@inline _mod360s(x::Float64) = mod(mod(x + 180.0, 360.0) + 360.0, 360.0) - 180.0

const _RAD2DEG = 180.0 / π
const _DEG2RAD = π / 180.0

"""
    _weights_tuple(weights) -> NTuple{7,Float64}

Pad/truncate a weight vector to a stack-allocated 7-tuple, so the sampler's inner
loop never re-allocates (`_pad_weights` returns a heap `Vector`).
"""
function _weights_tuple(weights::AbstractVector{<:Real})
    n = length(weights)
    return ntuple(i -> i <= n ? Float64(weights[i]) : 0.0, Val(7))
end

function _cvis_to_chi2(V::AbstractVector{<:Complex}, data::OIdata,
                       w::NTuple{7,Float64})
    chi2 = 0.0

    # ── V2 ────────────────────────────────────────────────────────────────────
    if w[1] > 0 && data.nv2 > 0
        s = 0.0
        idx = data.indx_v2
        @inbounds for i in eachindex(idx)
            d = Float64(abs2(V[idx[i]])) - Float64(data.v2[i])
            e = Float64(data.v2_err[i])
            s += d * d / (e * e)
        end
        chi2 += w[1] * s
    end

    # ── T3 (one pass computes the triple product for both amp and phi) ────────
    do_t3amp = w[2] > 0 && data.nt3amp > 0
    do_t3phi = w[3] > 0 && data.nt3phi > 0
    if do_t3amp || do_t3phi
        s_amp = 0.0
        s_phi = 0.0
        i1 = data.indx_t3_1; i2 = data.indx_t3_2; i3 = data.indx_t3_3
        @inbounds for i in eachindex(i1)
            t3 = V[i1[i]] * V[i2[i]] * V[i3[i]]
            if do_t3amp
                d = Float64(abs(t3)) - Float64(data.t3amp[i])
                e = Float64(data.t3amp_err[i])
                s_amp += d * d / (e * e)
            end
            if do_t3phi
                phi = Float64(angle(t3)) * _RAD2DEG
                dphi = _mod360s(phi - Float64(data.t3phi[i]))
                e = Float64(data.t3phi_err[i])
                s_phi += dphi * dphi / (e * e)
            end
        end
        chi2 += w[2] * s_amp + w[3] * s_phi
    end

    # ── Absolute VISAMP / VISPHI ──────────────────────────────────────────────
    do_visamp = w[4] > 0 && data.nvisamp > 0
    do_visphi = w[5] > 0 && data.nvisphi > 0
    if do_visamp || do_visphi
        s_amp = 0.0
        s_phi = 0.0
        idx = data.indx_vis
        @inbounds for i in eachindex(idx)
            Vv = V[idx[i]]
            if do_visamp
                d = Float64(abs(Vv)) - Float64(data.visamp[i])
                e = Float64(data.visamp_err[i])
                s_amp += d * d / (e * e)
            end
            if do_visphi
                dphi = _mod360s(Float64(angle(Vv)) * _RAD2DEG - Float64(data.visphi[i]))
                e = Float64(data.visphi_err[i])
                s_phi += dphi * dphi / (e * e)
            end
        end
        chi2 += w[4] * s_amp + w[5] * s_phi
    end

    return chi2
end

"""
    cvis_to_chi2_noalloc(V, data; weights=[1,1,1,0,0,0,0])

Weighted chi2 of model complex visibilities `V` against `data`, computed **without
allocating**.

Returns the same value as [`cvis_to_chi2_f`](@ref) (and as
`cvis_to_chi2_fg(V, data; weights)[1]`) to round-off, but fuses every
observable into a single pass with no gathered temporaries.  Intended for inner
sampler loops that evaluate chi2 millions of times — see `reconstruct_squeeze`.
Prefer [`cvis_to_chi2_f`](@ref) for ordinary use: it is clearer and additionally
handles the flux term (`weights[6]`), which this does not.
"""
function cvis_to_chi2_noalloc(V::AbstractVector{<:Complex}, data::OIdata;
                              weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
    return _cvis_to_chi2(V, data, _weights_tuple(weights))
end
