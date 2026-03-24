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
    _eval_W(model, x, wl, w_idx) -> Vector{Float64}

Evaluate the image chromatic weight W at each UV-point wavelength.
"""
function _eval_W(model::FlatModel, x::AbstractVector, wl, w_idx::Int)
    pv = model.resolver.fn(x, wl, NaN)
    W = pv[w_idx]
    return W isa Number ? fill(Float64(W), length(wl)) : collect(Float64, W)
end

"""
    _jacobian_W(model, x, wl, w_idx) -> Matrix{Float64}  (nUV × nparams)

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
    return ForwardDiff.jacobian(f, collect(Float64, x))
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

    chi2_v2 = 0.0; chi2_t3amp = 0.0; chi2_t3phi = 0.0
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
    return weights[1] * chi2_v2 + weights[2] * chi2_t3amp + weights[3] * chi2_t3phi
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
    x_img = reshape(x[nparams+1:end], nx, nx)
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

    # ── Observables ───────────────────────────────────────────────────────────
    v2_model = vis_to_v2(V_total, data.indx_v2)
    t3, t3amp_model, t3phi_model = vis_to_t3(V_total, data.indx_t3_1, data.indx_t3_2, data.indx_t3_3)

    chi2_v2 = 0.0; chi2_t3amp = 0.0; chi2_t3phi = 0.0

    # Parameter gradient (accumulated per-observable, matching chi2_flat_fg)
    g_params = zeros(Float64, nparams)

    # ── V2 ────────────────────────────────────────────────────────────────────
    g_v2 = 0.0
    if (weights[1] > 0) && (data.nv2 > 0)
        r_v2 = (v2_model .- data.v2) ./ data.v2_err.^2
        chi2_v2 = sum(r_v2 .* (v2_model .- data.v2))

        # Parameter gradient: g += 4 Re(J^T * (r_v2 .* conj(V)))
        rhs = r_v2 .* conj.(V_total[data.indx_v2])
        g_params .+= weights[1] .* 4 .* real.(transpose(J_total[data.indx_v2, :]) * rhs)

        # Image gradient via NFFT adjoint
        g_v2 = real(adjoint(ftplan[3]) * (4 .* r_v2 .* V_total[data.indx_v2] .* W[data.indx_v2]))
    end

    # ── T3amp ─────────────────────────────────────────────────────────────────
    g_t3amp = 0.0
    if (weights[2] > 0) && (data.nt3amp > 0)
        r_t3amp = 2.0 .* (t3amp_model .- data.t3amp) ./ data.t3amp_err.^2
        chi2_t3amp = sum(r_t3amp .* (t3amp_model .- data.t3amp)) / 2

        V1 = V_total[data.indx_t3_1]; V2c = V_total[data.indx_t3_2]; V3 = V_total[data.indx_t3_3]
        a1 = abs.(V1); a2 = abs.(V2c); a3 = abs.(V3)
        ε = eps(Float64)
        safe1 = max.(a1, ε); safe2 = max.(a2, ε); safe3 = max.(a3, ε)

        # Parameter gradient: g += Re(Σ_k J^T[t3_k,:] * rhs_k)
        rhs1 = r_t3amp .* conj.(V1) ./ safe1 .* a2 .* a3
        rhs2 = r_t3amp .* conj.(V2c) ./ safe2 .* a1 .* a3
        rhs3 = r_t3amp .* conj.(V3) ./ safe3 .* a1 .* a2
        g_params .+= weights[2] .* real.(
            transpose(J_total[data.indx_t3_1, :]) * rhs1 .+
            transpose(J_total[data.indx_t3_2, :]) * rhs2 .+
            transpose(J_total[data.indx_t3_3, :]) * rhs3)

        # Image gradient via NFFT adjoint
        g_t3amp = real(adjoint(ftplan[4]) * (r_t3amp .* V1 .* W[data.indx_t3_1] ./ safe1 .* a2 .* a3)) +
                  real(adjoint(ftplan[5]) * (r_t3amp .* V2c .* W[data.indx_t3_2] ./ safe2 .* a1 .* a3)) +
                  real(adjoint(ftplan[6]) * (r_t3amp .* V3 .* W[data.indx_t3_3] ./ safe3 .* a1 .* a2))
    end

    # ── T3phi ─────────────────────────────────────────────────────────────────
    g_t3phi = 0.0
    if (weights[3] > 0) && (data.nt3phi > 0)
        V1 = V_total[data.indx_t3_1]; V2c = V_total[data.indx_t3_2]; V3 = V_total[data.indx_t3_3]

        if !vonmises
            dphi = mod360(t3phi_model .- data.t3phi)
            chi2_t3phi = norm(dphi ./ data.t3phi_err)^2
            r_t3phi = (360.0 / π) .* dphi ./ data.t3phi_err.^2
        else
            dphi_rad = (t3phi_model .- data.t3phi) .* (π / 180.0)
            chi2_t3phi = sum(-2 .* data.t3phi_vonmises_err .* cos.(dphi_rad) .+ data.t3phi_vonmises_chi2_offset)
            r_t3phi = 2.0 .* data.t3phi_vonmises_err .* sin.(dphi_rad)
        end

        ε = eps(Float64)
        safe1 = max.(abs2.(V1), ε)
        safe2 = max.(abs2.(V2c), ε)
        safe3 = max.(abs2.(V3), ε)

        # Parameter gradient: g += Im(Σ_k J^T[t3_k,:] * rhs_k)
        rhs1 = r_t3phi .* conj.(V1) ./ safe1
        rhs2 = r_t3phi .* conj.(V2c) ./ safe2
        rhs3 = r_t3phi .* conj.(V3) ./ safe3
        g_params .+= weights[3] .* imag.(
            transpose(J_total[data.indx_t3_1, :]) * rhs1 .+
            transpose(J_total[data.indx_t3_2, :]) * rhs2 .+
            transpose(J_total[data.indx_t3_3, :]) * rhs3)

        # Image gradient via NFFT adjoint (SPARCO convention)
        dT3 = .-r_t3phi   # sign flip for NFFT adjoint convention
        g_t3phi = imag(adjoint(ftplan[4]) * (dT3 ./ safe1 .* V1 .* W[data.indx_t3_1]) +
                       adjoint(ftplan[5]) * (dT3 ./ safe2 .* V2c .* W[data.indx_t3_2]) +
                       adjoint(ftplan[6]) * (dT3 ./ safe3 .* V3 .* W[data.indx_t3_3]))
    end

    # ── Verbose output ────────────────────────────────────────────────────────
    if verb
        printstyled(@sprintf("V2: %.2f ", chi2_v2 / data.nv2), color=:red)
        printstyled(@sprintf("T3A: %.2f ", chi2_t3amp / data.nt3amp), color=:blue)
        printstyled(@sprintf("T3P: %.2f ", chi2_t3phi / data.nt3phi), color=:green)
        printstyled(@sprintf("Flux: %.4f ", sum(x_img)), color=:normal)
    end

    # ── Assemble gradient ─────────────────────────────────────────────────────
    g[1:nparams] = g_params

    # Image gradient: weighted sum of NFFT adjoint contributions
    g[nparams+1:end] = vec(weights[1] * g_v2 .+ weights[2] * g_t3amp .+ weights[3] * g_t3phi)

    return weights[1] * chi2_v2 + weights[2] * chi2_t3amp + weights[3] * chi2_t3phi
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
# Reconstruction  (VMLMB: joint params + image)
# ─────────────────────────────────────────────────────────────────────────────

"""
    reconstruct_sparco_flat(x_start, params_start, model, data, ft;
                            w_name="env,f", regularizers=[], weights=[1,1,1],
                            maxiter=100, verb=false, vonmises=false, ...)

SPARCO image reconstruction using a flat-model parametric component.

Returns `(params_final, x_final)`.
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
