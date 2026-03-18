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
#       → Float64     (forward only — cheap, for initial checks or gradient-free optimizers)
#
#   chi2_flat_fg(model, x, data; weights, verb, vonmises)
#       → (Float64, Vector{Float64})   (chi2 + gradient — for gradient-based optimizers)
#
#   chi2_flat_reduced(model, x, data; weights, vonmises)
#       → Float64     (chi2 / ndof — for human-readable goodness-of-fit)
#
# Polychromatic interface (data::Vector{OIdata}):
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


# ─────────────────────────────────────────────────────────────────────────────
# mod360: wrap phase residuals to [-180, 180]  (same as oichi2.jl)
# ─────────────────────────────────────────────────────────────────────────────

_mod360(x) = mod.(mod.(x .+ 180.0, 360.0) .+ 360.0, 360.0) .- 180.0


# ─────────────────────────────────────────────────────────────────────────────
# _chi2_terms: compute all chi2 components from model cvis
# Returns a NamedTuple with the values needed by both forward and gradient paths.
# ─────────────────────────────────────────────────────────────────────────────

function _chi2_terms(V::AbstractVector{<:Complex},
                     data::OIdata,
                     weights::AbstractVector{<:Real},
                     vonmises::Bool)

    w = _pad_weights(weights)
    w_v2, w_t3amp, w_t3phi, w_visamp, w_visphi = w[1], w[2], w[3], w[4], w[5]

    # ── V2 ────────────────────────────────────────────────────────────────────
    chi2_v2 = 0.0
    r_v2    = zeros(Float64, 0)     # residual scale: (v2_model - v2) / err²
    v2_model = Float64[]

    if w_v2 > 0 && data.nv2 > 0
        v2_model = abs2.(V[data.indx_v2])
        r_v2     = (v2_model .- data.v2) ./ data.v2_err.^2
        chi2_v2  = sum(r_v2 .* (v2_model .- data.v2))   # = norm(...)^2
    end

    # ── T3 products ───────────────────────────────────────────────────────────
    chi2_t3amp = 0.0
    chi2_t3phi = 0.0
    r_t3amp    = zeros(Float64, 0)
    r_t3phi    = zeros(Float64, 0)
    t3amp_model = Float64[]
    t3phi_model = Float64[]
    V1 = ComplexF64[]; V2c = ComplexF64[]; V3 = ComplexF64[]  # per-leg cvis

    need_t3 = (w_t3amp > 0 && data.nt3amp > 0) || (w_t3phi > 0 && data.nt3phi > 0)
    if need_t3
        V1  = V[data.indx_t3_1]
        V2c = V[data.indx_t3_2]
        V3  = V[data.indx_t3_3]
        t3  = V1 .* V2c .* V3
        t3amp_model = abs.(t3)
        t3phi_model = angle.(t3) .* (180.0 / π)
    end

    if w_t3amp > 0 && data.nt3amp > 0
        r_t3amp   = 2.0 .* (t3amp_model .- data.t3amp) ./ data.t3amp_err.^2
        chi2_t3amp = sum(r_t3amp .* (t3amp_model .- data.t3amp)) / 2   # norm(...)^2
    end

    if w_t3phi > 0 && data.nt3phi > 0
        if !vonmises
            dphi       = _mod360(t3phi_model .- data.t3phi)
            r_t3phi    = (360.0 / π) .* dphi ./ data.t3phi_err.^2
            chi2_t3phi = sum(dphi .* dphi ./ data.t3phi_err.^2)
        else
            dphi_rad   = (t3phi_model .- data.t3phi) .* (π / 180.0)
            r_t3phi    = 2.0 .* data.t3phi_vonmises_err .* sin.(dphi_rad)
            chi2_t3phi = sum(-2.0 .* data.t3phi_vonmises_err .* cos.(dphi_rad)
                             .+ data.t3phi_vonmises_chi2_offset)
        end
    end

    # ── Absolute VISAMP / VISPHI ───────────────────────────────────────────────
    chi2_visamp  = 0.0
    chi2_visphi  = 0.0
    r_visamp     = zeros(Float64, 0)
    r_visphi     = zeros(Float64, 0)
    Vvis         = ComplexF64[]
    visamp_model = Float64[]
    visphi_model = Float64[]

    need_vis = (w_visamp > 0 && data.nvisamp > 0) || (w_visphi > 0 && data.nvisphi > 0)
    if need_vis
        Vvis         = V[data.indx_vis]
        visamp_model = abs.(Vvis)
        visphi_model = angle.(Vvis) .* (180.0 / π)
    end

    if w_visamp > 0 && data.nvisamp > 0
        r_visamp    = 2.0 .* (visamp_model .- data.visamp) ./ data.visamp_err.^2
        chi2_visamp = sum(r_visamp .* (visamp_model .- data.visamp)) / 2
    end

    if w_visphi > 0 && data.nvisphi > 0
        dphi        = _mod360(visphi_model .- data.visphi)
        r_visphi    = (360.0 / π) .* dphi ./ data.visphi_err.^2
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

function _total_chi2(t, weights; chi2_flux=0.0, chi2_diffphase=0.0)
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
# _chi2_flux: compute flux chi2.
# If flux is calibrated (CALSTAT=C), compare model directly to data.
# If uncalibrated (CALSTAT=U), fit optimal scaling C = Σ(fm·fd·w) / Σ(fm²·w).
# ─────────────────────────────────────────────────────────────────────────────

function _chi2_flux(model_flux::Float64, data::OIdata)
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

function _chi2_flux_polychromatic(model_fluxes::Vector{Float64}, data::Vector{<:OIdata})
    fm = Float64[]; fd = Float64[]; fe = Float64[]; chan_idx = Int[]
    nwavs = length(data)
    for i in 1:nwavs
        d = data[i]; d.nflux > 0 || continue
        append!(fm, fill(model_fluxes[i], d.nflux))
        append!(fd, d.flux)
        append!(fe, d.flux_err)
        append!(chan_idx, fill(i, d.nflux))
    end
    isempty(fm) && return (; chi2=0.0, C=1.0, residual=Float64[], chan_idx=Int[])
    calibrated = !isempty(data) && data[1].flux_calibrated
    if calibrated
        residual = (fm .- fd) ./ fe.^2
        chi2 = sum(((fm .- fd) ./ fe).^2)
        C = 1.0
    else
        w = 1.0 ./ fe.^2
        C = sum(fm .* fd .* w) / sum(fm.^2 .* w)
        residual = (C .* fm .- fd) ./ fe.^2
        chi2 = sum(((C .* fm .- fd) ./ fe).^2)
    end
    return (; chi2, C, residual, chan_idx)
end


# ─────────────────────────────────────────────────────────────────────────────
# _verb_print: match the color-coded output format of oichi2.jl
# ─────────────────────────────────────────────────────────────────────────────

function _verb_print(t, data::OIdata, weights::AbstractVector{<:Real}, chi2_total::Float64;
                     chi2_flux::Float64=0.0, C_flux::Float64=NaN,
                     chi2_diffphase::Float64=0.0)
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
        -> Float64

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
                   vonmises::Bool = false)::Float64

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

    t = _chi2_terms(V, data, w, vonmises)

    # Flux chi2
    chi2_fl = 0.0; C_fl = NaN
    if need_flux
        fl = _chi2_flux(model_flux, data)
        chi2_fl = fl.chi2; C_fl = fl.C
    end

    chi2 = _total_chi2(t, w; chi2_flux=chi2_fl)
    verb && _verb_print(t, data, w, chi2; chi2_flux=chi2_fl, C_flux=C_fl)
    return chi2
end


# ─────────────────────────────────────────────────────────────────────────────
# chi2_flat_fg  (chi2 + gradient — monochromatic)
# ─────────────────────────────────────────────────────────────────────────────

"""
    chi2_flat_fg(model, x, data::OIdata; weights=[1,1,1,0,0,0,0], verb=false, vonmises=false)
        -> (Float64, Vector{Float64})

Evaluate chi-squared and its gradient w.r.t. `x`.

The gradient is computed analytically via the Wirtinger chain rule applied to
the complex visibility Jacobian J = ∂V/∂x from `eval_model_grad`.

Returns `(chi2, grad)` where `grad` is a `Vector{Float64}` of length `length(x)`.
"""
function chi2_flat_fg(model::FlatModel,
                      x::AbstractVector,
                      data::OIdata;
                      weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                      verb::Bool = false,
                      vonmises::Bool = false)

    w = _pad_weights(weights)
    w_v2, w_t3amp, w_t3phi, w_visamp, w_visphi, w_flux = w[1], w[2], w[3], w[4], w[5], w[6]
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

    # ── Observable chi2 terms ───────────────────────────────────────────────
    t    = _chi2_terms(V, data, w, vonmises)

    # ── Flux chi2 ───────────────────────────────────────────────────────────
    chi2_fl = 0.0; C_fl = NaN; r_flux = Float64[]
    if need_flux
        fl = _chi2_flux(model_flux, data)
        chi2_fl = fl.chi2; C_fl = fl.C; r_flux = fl.residual
    end

    chi2 = _total_chi2(t, w; chi2_flux=chi2_fl)

    # ── Gradient ────────────────────────────────────────────────────────────
    nparams = length(x)
    g = zeros(Float64, nparams)

    # V2: g += 4 Re( Jᵀ[v2_idx,:] · (r_v2 .* conj(V[v2_idx])) )
    if w_v2 > 0 && data.nv2 > 0
        rhs   = t.r_v2 .* conj.(V[data.indx_v2])
        g .+= w_v2 .* 4 .* real.(transpose(J[data.indx_v2, :]) * rhs)
    end

    # T3amp: g += Re( Σ_k Jᵀ[t3_k,:] · (r_t3amp .* conj(V_k)/|V_k| .* |V_other legs|) )
    if w_t3amp > 0 && data.nt3amp > 0
        a1 = abs.(t.V1); a2 = abs.(t.V2); a3 = abs.(t.V3)
        ε  = eps(Float64)
        safe1 = max.(a1, ε); safe2 = max.(a2, ε); safe3 = max.(a3, ε)
        rhs1  = t.r_t3amp .* conj.(t.V1) ./ safe1 .* a2 .* a3
        rhs2  = t.r_t3amp .* conj.(t.V2) ./ safe2 .* a1 .* a3
        rhs3  = t.r_t3amp .* conj.(t.V3) ./ safe3 .* a1 .* a2
        g .+= w_t3amp .* real.(
            transpose(J[data.indx_t3_1, :]) * rhs1 .+
            transpose(J[data.indx_t3_2, :]) * rhs2 .+
            transpose(J[data.indx_t3_3, :]) * rhs3
        )
    end

    # T3phi: g += Im( Σ_k Jᵀ[t3_k,:] · (r_t3phi .* conj(V_k) / |V_k|²) )
    if w_t3phi > 0 && data.nt3phi > 0
        ε     = eps(Float64)
        safe1 = max.(abs2.(t.V1), ε)
        safe2 = max.(abs2.(t.V2), ε)
        safe3 = max.(abs2.(t.V3), ε)
        rhs1  = t.r_t3phi .* conj.(t.V1) ./ safe1
        rhs2  = t.r_t3phi .* conj.(t.V2) ./ safe2
        rhs3  = t.r_t3phi .* conj.(t.V3) ./ safe3
        g .+= w_t3phi .* imag.(
            transpose(J[data.indx_t3_1, :]) * rhs1 .+
            transpose(J[data.indx_t3_2, :]) * rhs2 .+
            transpose(J[data.indx_t3_3, :]) * rhs3
        )
    end

    # Visamp: g += Re( Jᵀ[vis_idx,:] · (r_visamp .* conj(Vvis) / |Vvis|) )
    if w_visamp > 0 && data.nvisamp > 0
        ε    = eps(Float64)
        safe = max.(abs.(t.Vvis), ε)
        rhs  = t.r_visamp .* conj.(t.Vvis) ./ safe
        g .+= w_visamp .* real.(transpose(J[data.indx_vis, :]) * rhs)
    end

    # Visphi: g += Im( Jᵀ[vis_idx,:] · (r_visphi .* conj(Vvis) / |Vvis|²) )
    if w_visphi > 0 && data.nvisphi > 0
        ε    = eps(Float64)
        safe = max.(abs2.(t.Vvis), ε)
        rhs  = t.r_visphi .* conj.(t.Vvis) ./ safe
        g .+= w_visphi .* imag.(transpose(J[data.indx_vis, :]) * rhs)
    end

    # Flux: g += 2·C·Σ(r_flux)·Re(J_zero)
    # C=1 for calibrated flux; for uncalibrated, C_flux treated as constant
    if need_flux
        g .+= w_flux .* 2.0 .* C_fl .* sum(r_flux) .* real.(vec(J_zero))
    end

    verb && _verb_print(t, data, w, chi2; chi2_flux=chi2_fl, C_flux=C_fl)
    return chi2, g
end


# ─────────────────────────────────────────────────────────────────────────────
# chi2_flat  (forward only — polychromatic)
# ─────────────────────────────────────────────────────────────────────────────

"""
    chi2_flat(model, x, data::Vector{OIdata}; weights=[1,1,1,0,0,0,0], verb=false, vonmises=false)
        -> Float64

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
                   vonmises::Bool = false)::Float64

    w = _pad_weights(weights)
    nwavs = length(data)

    # Evaluate model at each channel's UV points
    need_flux = w[6] > 0 && any(d.nflux > 0 for d in data)
    if need_flux
        cvis = Vector{Vector{ComplexF64}}(undef, nwavs)
        model_fluxes = Vector{Float64}(undef, nwavs)
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
        model_fluxes = Float64[]
    end

    # Per-channel V2, T3, visamp, visphi
    chi2 = 0.0
    for i in 1:nwavs
        t = _chi2_terms(cvis[i], data[i], w, vonmises)
        chi2 += _total_chi2(t, w)  # flux/diffphase weights handled separately below
    end

    # Flux across all channels (single C_flux)
    chi2_fl = 0.0
    if need_flux
        fl = _chi2_flux_polychromatic(model_fluxes, data)
        chi2_fl = fl.chi2
        chi2 += w[6] * chi2_fl
    end

    # Differential phase across channels
    chi2_dp = 0.0
    if w[7] > 0 && nwavs > 1
        chi2_dp = _diffphase_chi2(cvis, data, nwavs)
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
    chi2_flat_fg(model, x, data::Vector{OIdata}; weights=[1,1,1,0,0,0,0], verb=false, vonmises=false)
        -> (Float64, Vector{Float64})

Polychromatic chi-squared and gradient.  Computes per-channel V2/T3/visamp/visphi/flux
gradients plus cross-channel differential phase gradient.
"""
function chi2_flat_fg(model::FlatModel,
                      x::AbstractVector,
                      data::Vector{<:OIdata};
                      weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                      verb::Bool = false,
                      vonmises::Bool = false)

    w = _pad_weights(weights)
    w_v2, w_t3amp, w_t3phi, w_visamp, w_visphi, w_flux, w_dp = w[1], w[2], w[3], w[4], w[5], w[6], w[7]
    nwavs = length(data)
    nparams = length(x)

    need_flux = w_flux > 0 && any(d.nflux > 0 for d in data)
    need_dp   = w_dp > 0 && nwavs > 1

    # ── Evaluate model + Jacobian at all channels ──────────────────────────
    cvis = Vector{Vector{ComplexF64}}(undef, nwavs)
    Jac  = Vector{Matrix{ComplexF64}}(undef, nwavs)
    model_fluxes = Vector{Float64}(undef, nwavs)
    J_zeros = Vector{Vector{ComplexF64}}(undef, nwavs)

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
    chi2 = 0.0
    g = zeros(Float64, nparams)

    for i in 1:nwavs
        V = cvis[i]
        J = Jac[i]
        d = data[i]
        t = _chi2_terms(V, d, w, vonmises)
        chi2 += _total_chi2(t, w)

        # V2 gradient
        if w_v2 > 0 && d.nv2 > 0
            rhs = t.r_v2 .* conj.(V[d.indx_v2])
            g .+= w_v2 .* 4 .* real.(transpose(J[d.indx_v2, :]) * rhs)
        end

        # T3amp gradient
        if w_t3amp > 0 && d.nt3amp > 0
            a1 = abs.(t.V1); a2 = abs.(t.V2); a3 = abs.(t.V3)
            ε = eps(Float64)
            safe1 = max.(a1, ε); safe2 = max.(a2, ε); safe3 = max.(a3, ε)
            rhs1 = t.r_t3amp .* conj.(t.V1) ./ safe1 .* a2 .* a3
            rhs2 = t.r_t3amp .* conj.(t.V2) ./ safe2 .* a1 .* a3
            rhs3 = t.r_t3amp .* conj.(t.V3) ./ safe3 .* a1 .* a2
            g .+= w_t3amp .* real.(
                transpose(J[d.indx_t3_1, :]) * rhs1 .+
                transpose(J[d.indx_t3_2, :]) * rhs2 .+
                transpose(J[d.indx_t3_3, :]) * rhs3
            )
        end

        # T3phi gradient
        if w_t3phi > 0 && d.nt3phi > 0
            ε = eps(Float64)
            safe1 = max.(abs2.(t.V1), ε)
            safe2 = max.(abs2.(t.V2), ε)
            safe3 = max.(abs2.(t.V3), ε)
            rhs1 = t.r_t3phi .* conj.(t.V1) ./ safe1
            rhs2 = t.r_t3phi .* conj.(t.V2) ./ safe2
            rhs3 = t.r_t3phi .* conj.(t.V3) ./ safe3
            g .+= w_t3phi .* imag.(
                transpose(J[d.indx_t3_1, :]) * rhs1 .+
                transpose(J[d.indx_t3_2, :]) * rhs2 .+
                transpose(J[d.indx_t3_3, :]) * rhs3
            )
        end

        # Visamp gradient
        if w_visamp > 0 && d.nvisamp > 0
            ε = eps(Float64)
            safe = max.(abs.(t.Vvis), ε)
            rhs = t.r_visamp .* conj.(t.Vvis) ./ safe
            g .+= w_visamp .* real.(transpose(J[d.indx_vis, :]) * rhs)
        end

        # Visphi gradient (absolute)
        if w_visphi > 0 && d.nvisphi > 0
            ε = eps(Float64)
            safe = max.(abs2.(t.Vvis), ε)
            rhs = t.r_visphi .* conj.(t.Vvis) ./ safe
            g .+= w_visphi .* imag.(transpose(J[d.indx_vis, :]) * rhs)
        end
    end

    # ── Flux chi2 + gradient (across all channels) ─────────────────────────
    chi2_fl = 0.0
    if need_flux
        fl = _chi2_flux_polychromatic(model_fluxes, data)
        chi2_fl = fl.chi2
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
                ε = eps(Float64)
                dT3 = -360.0 / π .* dphi ./ data[i].visphi_err.^2
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
        -> Float64

Return chi2 / ndof (reduced chi-squared).  Useful for assessing goodness of fit.
Works with both monochromatic (OIdata) and polychromatic (Vector{OIdata}) data.
"""
function chi2_flat_reduced(model::FlatModel,
                           x::AbstractVector,
                           data::Union{OIdata, Vector{<:OIdata}};
                           weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                           vonmises::Bool = false)::Float64
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
