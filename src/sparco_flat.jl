# sparco_flat.jl
#
# SPARCO image reconstruction using the flat model infrastructure.
#
# The combined visibility model is:
#
#   V_total(uv, λ) = V_model(uv, λ; params) + W(λ; params) · V_image(uv)
#
# where:
#   V_model  = eval_model(flat_model, params, uv; wl=λ) — parametric sources
#   W(λ)     = chromatic image weight from the resolver (a derived $WL expression)
#   V_image  = NFFT of the grey image (flux-normalized)
#
# The flat model dict defines parametric sources (stars, background, etc.) with
# $WL expressions for chromatic behavior.  The image weight W is a named
# resolver output (e.g. "env,f"), identified by the `w_name` argument.
#
# Gradient paths:
#   Parameters : eval_model_grad Jacobian + ForwardDiff of W through resolver
#   Image      : NFFT adjoint modulated by W (same structure as classic SPARCO)

using Printf, LinearAlgebra, NLopt


# ─────────────────────────────────────────────────────────────────────────────
# Standard SPARCO model templates
# ─────────────────────────────────────────────────────────────────────────────
#
# Star (point source) + resolved background + grey environment image.
# Chromatic flux fractions:
#   star,f = star,f0 · (λ/λ0)^(-star,di) / D
#   bg,f   = bg,f0   · (λ/λ0)^(-4)       / D       (Rayleigh-Jeans)
#   W      = (1 - star,f0 - bg,f0) · (λ/λ0)^(d_env) / D
# where D = sum of all numerators (total flux normalization).

const _SPARCO_EXPRS = Dict{String,Any}(
    "D" => "\$star,f0 * (\$WL/\$wl0)^(-\$star,di) + " *
           "\$bg,f0 * (\$WL/\$wl0)^(-4) + " *
           "(1 - \$star,f0 - \$bg,f0) * (\$WL/\$wl0)^(\$d_env)",
    "star,f" => "\$star,f0 * (\$WL/\$wl0)^(-\$star,di) / \$D",
    "bg,f"   => "\$bg,f0 * (\$WL/\$wl0)^(-4) / \$D",
    "W"      => "(1 - \$star,f0 - \$bg,f0) * (\$WL/\$wl0)^(\$d_env) / \$D",
    "bg,resolved" => true,
)

const _SPARCO_DEFAULTS = Dict{String,Any}(
    "wl0" => 1.6e-6, "star,f0" => 0.5, "star,di" => 4.0,
    "bg,f0" => 0.0, "d_env" => 0.0,
)

const _SPARCO_FREE_PARAMS          = ["star,f0", "bg,f0", "d_env"]
const _SPARCO_FREE_PARAMS_FREE_DI  = ["star,f0", "bg,f0", "d_env", "star,di"]

const _SPARCO_LB = Dict("star,f0" => 0.0, "bg,f0" => 0.0, "d_env" => -20.0, "star,di" => -20.0)
const _SPARCO_UB = Dict("star,f0" => 1.0, "bg,f0" => 1.0, "d_env" =>  20.0, "star,di" =>  20.0)


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    _resolve_w_idx(model, w_name) -> Int

Look up the index of the image-weight output in the resolver.
"""
function _resolve_w_idx(model::FlatModel, w_name::String)
    idx = findfirst(==(w_name), model.all_names)
    isnothing(idx) && error("Image weight '$w_name' not found in model.all_names = $(model.all_names)")
    return idx
end

"""
    _eval_W(model, x, wl, w_idx) -> Vector{T}

Evaluate the image chromatic weight W at each UV-point wavelength.
"""
function _eval_W(model::FlatModel, x::AbstractVector, wl, w_idx::Int)
    pv = model.resolver.fn(x, wl, NaN)
    W = pv[w_idx]
    # Model constants resolve to Float64; bring W to the wavelength grid's precision so a
    # Float32 pipeline is not promoted back to double through the chromatic image weight.
    Tw = float(eltype(wl))
    return W isa Number ? fill(convert(Tw, W), length(wl)) : collect(Tw, W)
end

"""
    _jacobian_W(model, x, wl, w_idx) -> Matrix{T}  (nUV × nparams)

Jacobian ∂W/∂x via ForwardDiff through the resolver.
"""
function _jacobian_W(model::FlatModel, x::AbstractVector, wl, w_idx::Int)
    nUV = length(wl)
    f = x_ -> begin
        T = eltype(x_)
        pv = model.resolver.fn(x_, wl, NaN)
        W = pv[w_idx]
        # Must preserve Dual type for ForwardDiff — no Float64 conversion
        return W isa Number ? fill(T(W), nUV) : Vector{T}(W)
    end
    return ForwardDiff.jacobian(f, collect(float(eltype(x)), x))
end


# ─────────────────────────────────────────────────────────────────────────────
# Forward model (no gradient)
# ─────────────────────────────────────────────────────────────────────────────

function _sparco_flat_forward(model::FlatModel, x, x_img, ftplan, data, w_idx::Int)
    wl = data.uv_lam
    V_model = eval_model(model, x, data.uv; wl=wl)
    W = _eval_W(model, x, wl, w_idx)
    V_image = image_to_vis(x_img, ftplan[1])
    V_total = V_model .+ W .* V_image
    return (; V_total, V_model, V_image, W)
end


# ─────────────────────────────────────────────────────────────────────────────
# chi2  (forward only)
# ─────────────────────────────────────────────────────────────────────────────

"""
    chi2_sparco_flat_f(x_img, params, model, ftplan, data;
                       w_name="env,f", weights=[1,1,1], verb=true, vonmises=false)

Evaluate chi2 for the SPARCO-flat model (no gradient).
"""
function chi2_sparco_flat_f(x_img::AbstractMatrix{<:AbstractFloat},
        params::AbstractVector{<:AbstractFloat},
        model::FlatModel,
        ftplan::AbstractVector{<:NFFT.NFFTPlan},
        data::OIdata;
        w_name::String = "W",
        weights = [1.0, 1.0, 1.0],
        verb::Bool = true,
        vonmises::Bool = false)

    w_idx = _resolve_w_idx(model, w_name)
    m = _sparco_flat_forward(model, params, x_img, ftplan, data, w_idx)

    v2_model = vis_to_v2(m.V_total, data.indx_v2)
    _, t3amp_model, t3phi_model = vis_to_t3(m.V_total, data.indx_t3_1, data.indx_t3_2, data.indx_t3_3)

    Ts = promote_type(float(eltype(x_img)), eltype(data.v2))
    chi2_v2 = zero(Ts); chi2_t3amp = zero(Ts); chi2_t3phi = zero(Ts)
    if (weights[1] > 0) && (data.nv2 > 0)
        chi2_v2 = norm((v2_model - data.v2) ./ data.v2_err)^2
    end
    if (weights[2] > 0) && (data.nt3amp > 0)
        chi2_t3amp = norm((t3amp_model - data.t3amp) ./ data.t3amp_err)^2
    end
    if (weights[3] > 0) && (data.nt3phi > 0)
        if !vonmises
            chi2_t3phi = norm(mod360(t3phi_model - data.t3phi) ./ data.t3phi_err)^2
        else
            chi2_t3phi = sum(-2 * data.t3phi_vonmises_err .* cos.((t3phi_model - data.t3phi) / 180 * pi) .+ data.t3phi_vonmises_chi2_offset)
        end
    end
    if verb
        printstyled(@sprintf("V2: %.2f ", chi2_v2 / data.nv2), color=:red)
        printstyled(@sprintf("T3A: %.2f ", chi2_t3amp / data.nt3amp), color=:blue)
        printstyled(@sprintf("T3P: %.2f ", chi2_t3phi / data.nt3phi), color=:green)
        printstyled(@sprintf("Flux: %.4f ", sum(x_img)), color=:normal)
    end
    return Ts(weights[1]) * chi2_v2 + Ts(weights[2]) * chi2_t3amp + Ts(weights[3]) * chi2_t3phi
end


# ─────────────────────────────────────────────────────────────────────────────
# chi2 + gradient
# ─────────────────────────────────────────────────────────────────────────────

"""
    chi2_sparco_flat_fg(x, g, model, ftplan, data, nparams;
                        w_name="env,f", weights=[1,1,1], verb=true, vonmises=false)

Chi2 and gradient for the combined [params; image_pixels] vector.

The gradient has two parts:
  - Parameter gradient: via J_total = J_model + V_image ⊙ J_W
  - Image gradient:     via NFFT adjoint modulated by W
"""
function chi2_sparco_flat_fg(x::AbstractVector{<:AbstractFloat},
        g::AbstractVector{<:AbstractFloat},
        model::FlatModel,
        ftplan::AbstractVector{<:NFFT.NFFTPlan},
        data::OIdata,
        nparams::Int;
        w_name::String = "W",
        weights = [1.0, 1.0, 1.0],
        verb::Bool = true,
        vonmises::Bool = false)

    params = x[1:nparams]
    nx = Int(sqrt(length(x) - nparams))
    # x holds params and image pixels in one vector, so it is Float64 whenever the model
    # parameters are: bring the image slice back to the operator's precision.
    x_img = to_ft_precision(reshape(x[nparams+1:end], nx, nx), ftplan)
    w_idx = _resolve_w_idx(model, w_name)
    wl = data.uv_lam

    # ── Forward pass ──────────────────────────────────────────────────────────
    V_model, J_model = eval_model_grad(model, params, data.uv; wl=wl)
    W = _eval_W(model, params, wl, w_idx)
    J_W = _jacobian_W(model, params, wl, w_idx)
    V_image = image_to_vis(x_img, ftplan[1])
    V_total = V_model .+ W .* V_image

    # Combined Jacobian for parameters: ∂V_total/∂params = J_model + V_image ⊙ J_W
    J_total = J_model .+ V_image .* J_W   # (nUV, nparams) broadcasting

    # ── Chi2 + adjoint source ─────────────────────────────────────────────────
    w = _pad_weights(weights)
    chi2, g_cvis = cvis_to_chi2_fg(V_total, data; weights=w, vonmises=vonmises)

    # ── Parameter gradient ────────────────────────────────────────────────────
    g[1:nparams] = real.(transpose(J_total) * g_cvis)

    # ── Image gradient via NFFT adjoint ───────────────────────────────────────
    # Identity: real(DFTᵀ · g_cvis) = real(NFFT^H · conj(g_cvis))
    g[nparams+1:end] = vec(real.(adjoint(ftplan[1]) * Complex{ft_eltype(ftplan)}.(conj.(W .* g_cvis))))

    # ── Verbose output ────────────────────────────────────────────────────────
    if verb
        t = _chi2_terms(V_total, data, w, vonmises)
        printstyled(@sprintf("V2: %.2f ", t.chi2_v2 / data.nv2), color=:red)
        printstyled(@sprintf("T3A: %.2f ", t.chi2_t3amp / data.nt3amp), color=:blue)
        printstyled(@sprintf("T3P: %.2f ", t.chi2_t3phi / data.nt3phi), color=:green)
        printstyled(@sprintf("Flux: %.4f ", sum(x_img)), color=:normal)
    end

    return chi2
end


# ─────────────────────────────────────────────────────────────────────────────
# model_and_image_to_vis: public convenience wrapper for the SPARCO forward model.
#
# Returns a NamedTuple (; V_total, V_model, V_image, W).
# ─────────────────────────────────────────────────────────────────────────────

"""
    model_and_image_to_vis(model, params, x_img, ftplan, data; w_name="W")

Compute the combined visibility V_total = V_model + W(λ) · V_image.

Returns a NamedTuple `(; V_total, V_model, V_image, W)`.
"""
function model_and_image_to_vis(model::FlatModel,
                                params::AbstractVector,
                                x_img::AbstractMatrix,
                                ftplan::AbstractVector{<:NFFT.NFFTPlan},
                                data::OIdata;
                                w_name::String = "W")
    w_idx = _resolve_w_idx(model, w_name)
    return _sparco_flat_forward(model, params, x_img, ftplan, data, w_idx)
end


# ─────────────────────────────────────────────────────────────────────────────
# model_and_image_to_chi2_fg: monochromatic variant
#
# When w_name is given, uses _eval_W / _jacobian_W for the image weight.
# When w_name is nothing, falls back to f_env = 1 - V_model(0,0).
# ─────────────────────────────────────────────────────────────────────────────

"""
    model_and_image_to_chi2_fg(model, x_params, image, g_image, ft, data::OIdata;
                                w_name=nothing, weights=[1,1,1], verb=false, vonmises=false)

Hybrid chi2+gradient combining a parametric model with an image (monochromatic).

When `w_name` is provided, the image weight `W(λ)` is evaluated from the named
model output and `V_total = V_model + W(λ) · V_image`.

When `w_name` is `nothing`, uses the SPARCO flux convention:
`V_total = (1 - V_model(0,0)) · V_image + V_model`.

Mutates `g_image` in-place. Returns `(chi2, g_params)`.
"""
function model_and_image_to_chi2_fg(
        model::FlatModel, x_params::AbstractVector,
        image::AbstractMatrix, g_image::AbstractMatrix,
        ft, data::OIdata;
        w_name::Union{String,Nothing} = nothing,
        weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0],
        verb::Bool = false, vonmises::Bool = false)

    w = _pad_weights(weights)

    # 1. Image → unit-flux cvis
    cvis_image = image_to_vis(image, ft)

    # 2. Model → cvis + Jacobian
    wl = hasproperty(data, :uv_lam) ? data.uv_lam : nothing
    cvis_model, J = eval_model_grad(model, x_params, data.uv; wl=wl)

    if !isnothing(w_name)
        # Use named W output from model
        w_idx = _resolve_w_idx(model, w_name)
        W = _eval_W(model, x_params, wl, w_idx)
        J_W = _jacobian_W(model, x_params, wl, w_idx)
        cvis_total = cvis_model .+ W .* cvis_image

        chi2, g_cvis = cvis_to_chi2_fg(cvis_total, data; weights=w, vonmises=vonmises)

        # Image gradient: W * adjoint + flux normalization
        _backprop_image_gradient!(g_image, W .* g_cvis, image, ft)

        # Parameter gradient: J_total = J_model + cvis_image .* J_W
        J_total = J .+ cvis_image .* J_W
        g_params = real.(transpose(J_total) * g_cvis)
    else
        # Fallback: f_env = 1 - V_model(0,0)
        uv_zero = zeros(2, 1)
        wl_zero = isnothing(wl) ? nothing : [wl[1]]
        cvis_zero, J_zero = eval_model_grad(model, x_params, uv_zero; wl=wl_zero)
        f_env = 1.0 - real(cvis_zero[1])

        cvis_total = f_env .* cvis_image .+ cvis_model

        chi2, g_cvis = cvis_to_chi2_fg(cvis_total, data; weights=w, vonmises=vonmises)

        _backprop_image_gradient!(g_image, g_cvis, image, ft, f_env)

        g_params = real.(transpose(J) * g_cvis)

        # f_env gradient chain rule
        g_f_env = real(sum(cvis_image .* g_cvis))
        df_env_dx = -real.(vec(J_zero[1, :]))
        g_params .+= g_f_env .* df_env_dx
    end

    if verb
        t = _chi2_terms(cvis_total, data, w, vonmises)
        _verb_print(t, data, w, chi2)
    end
    return chi2, g_params
end


# ─────────────────────────────────────────────────────────────────────────────
# model_and_image_to_chi2_fg: polychromatic variant (power-law grey image)
#
# Uses _eval_W / _jacobian_W from this file.
# ─────────────────────────────────────────────────────────────────────────────

"""
    model_and_image_to_chi2_fg(model, x_params, image, g_image, ft, data::Vector{<:OIdata};
                                w_name="W", weights=[1,1,1,0,0,0,0], verb=false, vonmises=false)

Polychromatic hybrid chi2+gradient for a parametric model + power-law grey image.

A single grey image is modulated by `W(λ)` at each wavelength channel, following
the SPARCO convention: `V_total(λ) = V_model(λ) + W(λ) · V_image`.

Mutates `g_image` in-place. Returns `(chi2, g_params)`.
"""
function model_and_image_to_chi2_fg(
        model::FlatModel, x_params::AbstractVector,
        image::AbstractMatrix, g_image::AbstractMatrix,
        ft, data::Vector{<:OIdata};
        w_name::String = "W",
        weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        verb::Bool = false, vonmises::Bool = false)

    w = _pad_weights(weights)
    nwavs = length(data)

    # Image → unit-flux cvis (single grey image)
    cvis_image = image_to_vis(image, ft[1])

    w_idx = _resolve_w_idx(model, w_name)
    Ts = promote_type(float(eltype(x_params)), eltype(first(data).uv))
    chi2 = zero(Ts)
    g_params = zeros(Ts, length(x_params))
    g_image .= 0

    for i in 1:nwavs
        wl_i = data[i].uv_lam
        V_model_i, J_i = eval_model_grad(model, x_params, data[i].uv; wl=wl_i)
        W_i = _eval_W(model, x_params, wl_i, w_idx)
        J_W_i = _jacobian_W(model, x_params, wl_i, w_idx)
        V_total_i = V_model_i .+ W_i .* cvis_image

        chi2_i, g_cvis_i = cvis_to_chi2_fg(V_total_i, data[i]; weights=w, vonmises=vonmises)
        chi2 += chi2_i

        # Parameter gradient: J_total = J_model + cvis_image .* J_W
        J_total_i = J_i .+ cvis_image .* J_W_i
        g_params .+= real.(transpose(J_total_i) * g_cvis_i)

        # Image gradient: accumulate across channels
        _backprop_image_gradient!(g_image, W_i .* g_cvis_i, image, ft[1]; accumulate=true)
    end

    if verb
        ndof = _ndof(data, w)
        printstyled(@sprintf("χ²r: %.3f\n", chi2 / ndof), color=:white)
    end
    return chi2, g_params
end


# ─────────────────────────────────────────────────────────────────────────────
# Criterion  (chi2 + regularization)
# ─────────────────────────────────────────────────────────────────────────────

function crit_sparco_flat_fg(x::AbstractVector{<:AbstractFloat},
        g::AbstractVector{<:AbstractFloat},
        model::FlatModel,
        ftplan::AbstractVector{<:NFFT.NFFTPlan},
        data::OIdata,
        nparams::Int;
        w_name::String = "W",
        weights = [1.0, 1.0, 1.0],
        printcolor = :normal,
        regularizers = [],
        verb::Bool = true,
        vonmises::Bool = false)

    chi2 = chi2_sparco_flat_fg(x, g, model, ftplan, data, nparams;
        w_name=w_name, weights=weights, verb=verb, vonmises=vonmises)

    nx = Int(sqrt(length(x) - nparams))
    reg_g = zeros(eltype(x), nx, nx)
    reg_f = regularization(reshape(x[nparams+1:end], nx, nx), reg_g;
        regularizers=regularizers, printcolor=printcolor, verb=verb)
    g[nparams+1:end] += vec(reg_g)

    # Flux-normalization correction for image gradient
    flux = sum(x[nparams+1:end])
    g[nparams+1:end] = (g[nparams+1:end] .- sum(vec(x[nparams+1:end]) .* g[nparams+1:end]) / flux) / flux

    return chi2 + reg_f
end


# ─────────────────────────────────────────────────────────────────────────────
# Parameter optimization  (image fixed)
# ─────────────────────────────────────────────────────────────────────────────

"""
    optimize_sparco_flat_parameters(params_start, model, x_img, ft, data;
                                    w_name="W", weights=[1,1,1], lb, ub,
                                    method=:LN_NELDERMEAD, maxeval=2000)

Optimize SPARCO-flat model parameters with the image held fixed.
Supports gradient-free (`:LN_NELDERMEAD`) and gradient-based (`:LD_LBFGS`, etc.) methods.
For `:LD_*` methods, the parameter gradient is extracted from `chi2_sparco_flat_fg`.
Returns `(minchi2, params_opt, ret)`.
"""
function optimize_sparco_flat_parameters(params_start::AbstractVector{<:AbstractFloat},
        model::FlatModel,
        x_img::AbstractMatrix{<:AbstractFloat},
        ft::AbstractVector{<:NFFT.NFFTPlan},
        data::OIdata;
        w_name::String = "W",
        weights = [1.0, 1.0, 1.0],
        lb::Union{Nothing, Vector{Float64}} = nothing,
        ub::Union{Nothing, Vector{Float64}} = nothing,
        method::Symbol = :LN_NELDERMEAD,
        maxeval::Int = 2000)

    np = length(params_start)
    if lb === nothing
        lb = fill(-Inf, np)
    end
    if ub === nothing
        ub = fill(Inf, np)
    end

    use_grad = startswith(string(method), "LD_")

    function f_obj(p, grad)
        if use_grad && length(grad) > 0
            x_full = [p; vec(x_img)]
            g_full = zeros(length(x_full))
            chi2 = chi2_sparco_flat_fg(x_full, g_full, model, ft, data, np;
                w_name=w_name, weights=weights, verb=false)
            grad .= g_full[1:np]
            return chi2
        else
            return chi2_sparco_flat_f(x_img, p, model, ft, data;
                w_name=w_name, weights=weights, verb=false)
        end
    end

    optimizer = Opt(method, np)
    min_objective!(optimizer, f_obj)
    lower_bounds!(optimizer, lb)
    upper_bounds!(optimizer, ub)
    maxeval!(optimizer, maxeval)
    minchi2, params_opt, ret = optimize(optimizer, collect(Float64, params_start))
    return minchi2, params_opt, ret
end


# ─────────────────────────────────────────────────────────────────────────────
# Reconstruction  (VMLMB image + NelderMead params, alternating)
# ─────────────────────────────────────────────────────────────────────────────

"""
    reconstruct_hybrid(x_start, params_start, model, data, ft;
                              w_name="W", regularizers=[], weights=[1,1,1],
                              maxiter=200, rounds=5, verb=false, vonmises=false,
                              params_lower=nothing, params_upper=nothing, ...)

SPARCO image reconstruction using a flat-model parametric component.

Alternates between:
1. **Image VMLMB**: optimize image pixels with parameters fixed, using
   `model_and_image_to_chi2_fg` + regularization.
2. **Parameter NelderMead**: optimize model parameters with the image fixed,
   using `optimize_sparco_flat_parameters`.

This separation avoids the scale-mismatch problem of jointly optimizing
parameters and image pixels in a single VMLMB run.

Returns `(params_final, x_final)`.
"""
function reconstruct_hybrid(x_start::AbstractMatrix{<:AbstractFloat},
        params_start::AbstractVector{<:AbstractFloat},
        model::FlatModel,
        data::Union{OIdata, Vector{<:OIdata}},
        ft;
        w_name::String = "W",
        printcolor = :normal,
        verb::Bool = false,
        maxiter::Int = 200,
        rounds::Int = 5,
        regularizers = [],
        weights = [1.0, 1.0, 1.0],
        vonmises::Bool = false,
        params_lower::Union{Nothing, Vector{Float64}} = nothing,
        params_upper::Union{Nothing, Vector{Float64}} = nothing,
        ftol = (0, 1e-8),
        xtol = (0, 1e-8),
        gtol = (0, 1e-8))

    nx = size(x_start, 1)
    w = _pad_weights(weights)

    # Unwrap ft for mono data: Matrix{Vector{NFFTPlan}} → Vector{NFFTPlan}
    ft_use = (data isa OIdata && ft isa AbstractMatrix) ? ft[1] : ft

    # The image runs at the precision of the FT operator (Float32 by default); the model
    # parameters stay Float64, being a handful of scalars where the precision is free.
    x_img = to_ft_precision(copy(x_start), ft_use)
    T = eltype(x_img)
    x_params = collect(Float64, params_start)

    for round in 1:rounds
        # ── Phase 1: VMLMB over image pixels (params fixed) ──────────────
        g_image = zeros(T, nx, nx)

        function crit_image(x_flat, g_flat)
            img = reshape(x_flat, nx, nx)
            chi2, _ = model_and_image_to_chi2_fg(model, x_params, img, g_image, ft_use, data;
                w_name=w_name, weights=w, verb=verb, vonmises=vonmises)
            reg_g = zeros(T, nx, nx)
            reg_f = regularization(img, reg_g; regularizers=regularizers,
                printcolor=printcolor, verb=verb)
            g_flat .= vec(g_image .+ reg_g)
            return chi2 + reg_f
        end

        x_img_flat = OptimPackNextGen.vmlmb(crit_image, vec(copy(x_img));
            lower=0, upper=Inf,
            maxiter=maxiter, verb=verb, blmvm=false,
            xtol=xtol, ftol=ftol, gtol=gtol)
        x_img = reshape(x_img_flat, nx, nx)

        # ── Phase 2: NelderMead over parameters (image fixed) ────────────
        # Only for monochromatic data with NFFT plans — use optimize_sparco_flat_parameters
        if data isa OIdata && ft_use isa AbstractVector{<:NFFT.NFFTPlan}
            minchi2, x_params, ret = optimize_sparco_flat_parameters(x_params, model, x_img, ft_use, data;
                w_name=w_name, weights=collect(Float64, weights),
                lb=params_lower, ub=params_upper)
            if verb
                @printf("  Round %d param opt: chi2=%.2f  ret=%s\n", round, minchi2, ret)
            end
        end
    end

    return (x_params, x_img)
end

"""
    reconstruct_sparco_flat(x_start, params_start, model, data, ft; kwargs...)

Deprecated. Use `reconstruct_hybrid` instead.
"""
function reconstruct_sparco_flat(x_start::AbstractMatrix{<:AbstractFloat},
        params_start::AbstractVector{<:AbstractFloat},
        model::FlatModel,
        data::OIdata,
        ft::AbstractVector{<:NFFT.NFFTPlan};
        w_name::String = "W",
        printcolor = :normal,
        verb::Bool = false,
        maxiter::Int = 100,
        regularizers = [],
        weights = [1.0, 1.0, 1.0],
        vonmises::Bool = false,
        params_lower::Union{Nothing, Vector{Float64}} = nothing,
        params_upper::Union{Nothing, Vector{Float64}} = nothing,
        ftol = (0, 1e-8),
        xtol = (0, 1e-8),
        gtol = (0, 1e-8))

    nparams = length(params_start)
    npix = length(x_start)

    crit = (x, g) -> crit_sparco_flat_fg(x, g, model, ft, data, nparams;
        w_name=w_name, regularizers=regularizers, verb=verb,
        weights=weights, vonmises=vonmises)

    # Parameter bounds: default to unconstrained; image pixels are non-negative
    plb = params_lower === nothing ? fill(-Inf, nparams) : params_lower
    pub = params_upper === nothing ? fill(Inf, nparams) : params_upper
    lower = [plb; fill(0.0, npix)]
    upper = [pub; fill(Inf, npix)]

    sol = OptimPackNextGen.vmlmb(crit, [params_start; vec(x_start)]; verb=verb,
        lower=lower, upper=upper, maxiter=maxiter, blmvm=false,
        xtol=xtol, ftol=ftol, gtol=gtol)

    return (sol[1:nparams], reshape(sol[nparams+1:end], size(x_start)))
end


# ─────────────────────────────────────────────────────────────────────────────
# reconstruct_sparco: lean high-level interface
# ─────────────────────────────────────────────────────────────────────────────

"""
    reconstruct_sparco(x_start, data, ft;
                       lambda_ref, star_flux=0.5, bg_flux=0.0,
                       d_env=0.0, star_di=4.0,
                       star_spectral_index_fixed=true,
                       fit_model_first=true,
                       regularizers=[], weights=[1,1,1],
                       maxiter=200, rounds=3, verb=false, ...)

High-level SPARCO image reconstruction with the standard star + background + image model.

Builds the flat-model dict internally from the physical parameters, then calls
`reconstruct_hybrid`.

# Arguments
- `x_start`:  starting image (nx × nx)
- `data`:     `OIdata` or `Vector{<:OIdata}`
- `ft`:       Fourier transform plans from `setup_nfft`
- `lambda_ref`:  reference wavelength in metres (required)
- `star_flux`:   initial stellar flux fraction at λ₀ (default 0.5)
- `bg_flux`:     initial background flux fraction at λ₀ (default 0.0)
- `d_env`:       environment spectral index (default 0.0)
- `star_di`:     stellar spectral index (default 4.0 = Rayleigh-Jeans)
- `star_spectral_index_fixed`:  fix `star_di` during reconstruction (default `true`)
- `fit_model_first`:  run a model-only fit before image reconstruction (default `true`)

Returns a named tuple `(params, image, model, free_params, chi2)`.
"""
function reconstruct_sparco(x_start::AbstractMatrix{<:AbstractFloat},
        data::Union{OIdata, Vector{<:OIdata}},
        ft;
        lambda_ref::Real,
        star_flux::Real = 0.5,
        bg_flux::Real = 0.0,
        d_env::Real = 0.0,
        star_di::Real = 4.0,
        star_spectral_index_fixed::Bool = true,
        fit_model_first::Bool = true,
        weights = [1.0, 1.0, 1.0],
        regularizers = [],
        maxiter::Int = 200,
        rounds::Int = 3,
        verb::Bool = false,
        vonmises::Bool = false,
        printcolor = :normal,
        ftol = (0, 1e-8),
        xtol = (0, 1e-8),
        gtol = (0, 1e-8))

    # Build model dict
    model_dict = merge(_SPARCO_DEFAULTS, _SPARCO_EXPRS)
    model_dict["wl0"]     = Float64(lambda_ref)
    model_dict["star,f0"] = Float64(star_flux)
    model_dict["bg,f0"]   = Float64(bg_flux)
    model_dict["d_env"]   = Float64(d_env)
    model_dict["star,di"] = Float64(star_di)

    free_params = star_spectral_index_fixed ? _SPARCO_FREE_PARAMS : _SPARCO_FREE_PARAMS_FREE_DI
    model = dict_to_model(model_dict, free_params)

    # Bounds
    lb = Float64[_SPARCO_LB[k] for k in free_params]
    ub = Float64[_SPARCO_UB[k] for k in free_params]
    x_params = Float64[model_dict[k] for k in free_params]

    # Optional model-only fit
    if fit_model_first
        data_single = data isa Vector ? data[1] : data
        if verb
            printstyled("SPARCO: fitting model parameters first...\n", color=:cyan)
        end
        result = fit_model(model_dict, free_params, data_single;
            weights=[collect(Float64, weights); 0.0; 0.0; 0.0; 0.0],
            method=:LD_LBFGS,
            lb=Dict(k => _SPARCO_LB[k] for k in free_params),
            ub=Dict(k => _SPARCO_UB[k] for k in free_params),
            maxeval=5000)
        x_params = result.x_opt
        if verb
            println(result)
        end
        # Rebuild model with fitted values
        for (i, k) in enumerate(free_params)
            model_dict[k] = x_params[i]
        end
        model = dict_to_model(model_dict, free_params)
    end

    # Reconstruct
    x_params, x_img = reconstruct_hybrid(x_start, x_params, model, data, ft;
        w_name="W", regularizers=regularizers, weights=weights,
        maxiter=maxiter, rounds=rounds, verb=verb, vonmises=vonmises,
        printcolor=printcolor, params_lower=lb, params_upper=ub,
        ftol=ftol, xtol=xtol, gtol=gtol)

    # Final chi2
    data_single = data isa Vector ? data[1] : data
    ft_vec = ft isa AbstractVector{<:NFFT.NFFTPlan} ? ft :
             ft isa AbstractMatrix ? ft[1] : ft[1]
    chi2 = chi2_sparco_flat_f(x_img, x_params, model, ft_vec, data_single;
        w_name="W", weights=weights, verb=verb)

    return (; params=x_params, image=x_img, model, free_params, chi2,
              param_names=Dict(free_params[i] => x_params[i] for i in eachindex(free_params)))
end
