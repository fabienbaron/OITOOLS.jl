# ═══════════════════════════════════════════════════════════════════════════════
# BSDMM (Block-Structured SDMM) for optical interferometry
#
# Architecture:
#   f(x) = chi2(H(x/sum(x)))          — non-convex data fit (VMLMB, lower=0)
#   g_i(L_i x) = non-smooth regs      — TV, centering (proximal z-updates)
#
# Positivity and flux normalization are handled inside the x-update
# (VMLMB bounds + chain rule), NOT as separate ADMM blocks. Only truly
# non-smooth regularizers (TV, L2smooth, centering) go into proximal blocks.
#
# All L_i = I (identity). Polychromatic mode couples N wavelength
# channels via group sparsity or group TV transspectral regularizers.
#
# Per-block penalty ρ_i adapted via three strategies:
#   :spectral  — AADMM adaptive BB + predictive dual + cross-block balancing
#   :balanced  — residual balancing (primal/dual ratio)
#   :aradmm   — ARADMM spectral + over-relaxation γ
#
# Reference: Moolekamp & Melchior (arXiv:1708.09066), Xu & Taylor (ARADMM),
#            Xu (AMADMM multi-block), Goldstein (FASTA adaptive BB)
# ═══════════════════════════════════════════════════════════════════════════════

# ─── Block data structure ────────────────────────────────────────────────────

mutable struct BSDMMBlock{T<:AbstractFloat}
    name::Symbol
    z::Vector{T}         # auxiliary variable
    u::Vector{T}         # scaled dual variable
    z_prev::Vector{T}    # previous z (for adaptive update)
    u_prev::Vector{T}    # previous u (for adaptive update)
    prox!::Any           # prox!(z, v, rho) modifies z in place
    rho::T               # penalty parameter
    rho_min::T
    rho_max::T
    # Storage for spectral/ARADMM adaptive updates
    z_adp::Vector{T}       # z at previous adaptive-update iteration
    lambda_adp::Vector{T}  # unscaled dual (ρu) at previous adaptive-update iteration
    lhat_adp::Vector{T}    # predictive dual at previous adaptive-update iteration
end

function BSDMMBlock(name::Symbol, n::Int, prox!;
                    rho=10.0, rho_min=1e-3, rho_max=1e8)
    BSDMMBlock{Float64}(name,
        zeros(n), zeros(n), zeros(n), zeros(n),
        prox!, rho, rho_min, rho_max,
        zeros(n), zeros(n), zeros(n))
end

# ─── Adaptive BB selection (Goldstein FASTA paper) ─────────────────────────
# Ported from aradmm_core.jl (Baron/Xu)

function curv_adaptive_BB(αSD, αMG)
    if 2.0 * αMG > αSD
        return αMG
    else
        return αSD - 0.5 * αMG
    end
end

# ─── Spectral penalty update (AADMM/ARADMM) ──────────────────────────────
# Adapted from aradmm_core.jl:aradmm_estimate (two-block) and
# amadmm_core.m (multi-block fallback) for BSDMM consensus form.
#
# For each block i (consensus: A=I, B=-I, b=0):
#   x-block curvature α: from (Δx, Δl_hat)  where l_hat isolates x changes
#   z-block curvature β: from (-Δz, Δλ)
#   Combined: ρ = √(αβ)
#
# Returns (rho_estimate, gamma) where rho_estimate is nothing if curvature
# could not be estimated (caller uses cross-block fallback).

function spectral_update!(blk::BSDMMBlock, x::Vector, x_adp::Vector;
                          ε_cor=0.2, compute_gamma=false,
                          fallback_rho=nothing)
    minval = 1e-30
    rho = blk.rho

    # Current unscaled dual and predictive dual
    # l_hat = what dual would be if z had not changed (isolates x-block)
    # BSDMM loop: z-update, x-update, dual-update
    # u_prev/z_prev are from start of this iteration (before z-update)
    lambda = rho .* blk.u
    lhat = rho .* (blk.u_prev .+ x .- blk.z_prev)

    # Deltas vs previous adaptive-update iteration
    Δlhat = lhat .- blk.lhat_adp
    Δlambda = lambda .- blk.lambda_adp
    Δx = x .- x_adp
    Δz = blk.z .- blk.z_adp

    # Norms
    dl_hat = norm(Δlhat)
    dl = norm(Δlambda)
    du = norm(Δx)
    dv = norm(Δz)

    # x-block curvature (α) — uses predictive dual
    ul_hat = dot(Δx, Δlhat)
    hflag = (du * dl_hat > minval) && (ul_hat > ε_cor * du * dl_hat + minval)
    α = rho  # fallback
    if hflag
        αSD = dl_hat^2 / ul_hat
        αMG = ul_hat / du^2
        α = curv_adaptive_BB(αSD, αMG)
    end

    # z-block curvature (β) — uses actual dual, B=-I so ΔBv = -Δz
    vl = -dot(Δz, Δlambda)   # = dot(-Δz, Δλ)
    gflag = (dv * dl > minval) && (vl > ε_cor * dv * dl + minval)
    β = rho  # fallback
    if gflag
        βSD = dl^2 / vl
        βMG = vl / dv^2
        β = curv_adaptive_BB(βSD, βMG)
    end

    # Combine curvatures
    rho_est = nothing
    gamma = 1.0
    if hflag && gflag
        rho_est = sqrt(α * β)
        if compute_gamma
            gamma = min(1.0 + 2.0 * sqrt(α * β) / (α + β), 1.9)
        end
    elseif hflag
        rho_est = α
        gamma = compute_gamma ? 1.9 : 1.0
    elseif gflag
        rho_est = β
        gamma = compute_gamma ? 1.1 : 1.0
    else
        # Neither estimable — use fallback from other blocks (AMADMM line 129)
        if fallback_rho !== nothing
            rho_est = fallback_rho
        end
        gamma = compute_gamma ? 1.5 : 1.0
    end

    # Apply penalty update
    if rho_est !== nothing
        rho_new = clamp(rho_est, blk.rho_min, blk.rho_max)
        rho_old = blk.rho
        blk.rho = rho_new
        blk.u .*= rho_old / blk.rho
    end

    # Record for next adaptive update
    blk.z_adp .= blk.z
    blk.lambda_adp .= blk.rho .* blk.u  # use current rho (post-update); λ = ρu is invariant
    blk.lhat_adp .= lhat

    return (rho_est, gamma)
end

# ─── Residual balancing ───────────────────────────────────────────────────
# Ported from aradmm_core.jl lines 134-142 (adp==3)

function balanced_update!(blk::BSDMMBlock, x::Vector;
                          beta_scale=2.0, res_scale=0.1)
    r_norm = norm(x .- blk.z)                      # primal residual
    s_norm = blk.rho * norm(blk.z .- blk.z_prev)   # dual residual
    rho_old = blk.rho
    if s_norm < r_norm * res_scale
        blk.rho = clamp(beta_scale * blk.rho, blk.rho_min, blk.rho_max)
    elseif r_norm < s_norm * res_scale
        blk.rho = clamp(blk.rho / beta_scale, blk.rho_min, blk.rho_max)
    end
    if blk.rho != rho_old
        blk.u .*= rho_old / blk.rho
    end
end

# ─── Proximal operator factories ─────────────────────────────────────────────

function make_prox_tv(nx, mu; niter=50)
    return function(z, v, rho)
        v_2d = reshape(v, nx, nx)
        z .= vec(prox_tv(v_2d, mu / rho; niter=niter))
    end
end

function make_prox_l2smooth(nx, mu)
    return function(z, v, rho)
        v_2d = reshape(v, nx, nx)
        z .= vec(prox_l2smooth(v_2d, mu / rho))
    end
end

function make_prox_centering(nx, mu)
    npix = nx * nx
    center = (nx + 1) / 2.0
    p = Float64[mod(i - 1, nx) + 1 - center for i in 1:npix]
    q = Float64[div(i - 1, nx) + 1 - center for i in 1:npix]
    AtA = [dot(p, p) dot(p, q); dot(q, p) dot(q, q)]
    return function(z, v, rho)
        prox_centering!(z, v, mu / rho, p, q, AtA)
    end
end

# ─── Polychromatic proximal operator factories ───────────────────────────────

function make_prox_group_sparsity(nx, nwav, mu)
    return function(z, v, rho)
        v_3d = reshape(v, nx, nx, nwav)
        z .= vec(prox_group_sparsity(v_3d, mu / rho))
    end
end

function make_prox_grouptv(nx, nwav, mu; niter=50)
    return function(z, v, rho)
        v_3d = reshape(v, nx, nx, nwav)
        z .= vec(prox_grouptv(v_3d, mu / rho; niter=niter))
    end
end

function make_prox_tv_poly(nx, nwav, mu; niter=50)
    npix = nx * nx
    return function(z, v, rho)
        for w in 1:nwav
            idx = (w-1)*npix+1 : w*npix
            v_2d = reshape(@view(v[idx]), nx, nx)
            @view(z[idx]) .= vec(prox_tv(v_2d, mu / rho; niter=niter))
        end
    end
end

function make_prox_l2smooth_poly(nx, nwav, mu)
    npix = nx * nx
    return function(z, v, rho)
        for w in 1:nwav
            idx = (w-1)*npix+1 : w*npix
            v_2d = reshape(@view(v[idx]), nx, nx)
            @view(z[idx]) .= vec(prox_l2smooth(v_2d, mu / rho))
        end
    end
end

function make_prox_centering_poly(nx, nwav, mu)
    npix = nx * nx
    center = (nx + 1) / 2.0
    p = Float64[mod(i - 1, nx) + 1 - center for i in 1:npix]
    q = Float64[div(i - 1, nx) + 1 - center for i in 1:npix]
    AtA = [dot(p, p) dot(p, q); dot(q, p) dot(q, q)]
    return function(z, v, rho)
        for w in 1:nwav
            idx = (w-1)*npix+1 : w*npix
            prox_centering!(@view(z[idx]), @view(v[idx]), mu / rho, p, q, AtA)
        end
    end
end

# ─── Main BSDMM reconstruction ──────────────────────────────────────────────

"""
    reconstruct_bsdmm(x_init, data, ft; kwargs...)

Block-Structured SDMM (BSDMM) image reconstruction for optical interferometry.

Solves:  `minimize_x  f(x) + Σ_i g_i(x)`

where `f(x) = χ²(H(x/sum(x)))` is the non-convex data fidelity term (squared
visibilities, closure phases) and each `g_i` is a non-smooth regularizer
(total variation, L2 smoothness, centering) handled via its proximal operator.

The BSDMM splits the problem into alternating steps:
1. **z-updates** (proximal): `z_i = prox_{g_i/ρ_i}(x + u_i)` — one per regularizer
2. **x-update** (VMLMB): minimize `χ² + Σ_i (ρ_i/2)||x - z_i + u_i||²` with `x ≥ 0`
3. **Dual updates**: `u_i ← u_i + x - z_i`

Positivity is enforced via VMLMB lower bounds. Flux normalization uses the
chain-rule gradient `∇_x χ²(H(x/sum(x)))` inside the x-update, plus explicit
re-normalization after each step.

Supports monochromatic (`nx × nx` image, `(1,1)` data) and polychromatic
(`nx × nx × nwav × 1` image, `(nwav,1)` data) modes.

# Arguments
- `x_init`: starting image (`nx × nx` for mono, `nx × nx × nwav × 1` for poly)
- `data`: interferometric data (`Matrix{OIdata}`)
- `ft`: Fourier transform plans from `setup_ft`

# Keyword arguments (monochromatic)
- `mu_reg=0.0`: regularization weight (0 = no regularization)
- `mu_cen=0.0`: centering weight (0 = no centering)
- `reg_type=:tv`: regularization type (`:tv` or `:l2smooth`)

# Keyword arguments (polychromatic)
- `mu_tv=0.0`: per-channel spatial TV/L2smooth weight
- `mu_group=0.0`: transspectral coupling weight (group sparsity or group TV)
- `mu_cen=0.0`: per-channel centering weight
- `reg_type=:tv`: spatial regularizer (`:tv` or `:l2smooth`)
- `group_type=:sparsity`: transspectral regularizer (`:sparsity` or `:tv`)

# Common keyword arguments
- `tv_niter=50`: Chambolle iterations for TV proximal operator
- `rho_init=10.0`: initial ADMM penalty parameter for all blocks
- `rho_min=1e-3`: minimum penalty (for adaptive update)
- `rho_max=1e4`: maximum penalty (for adaptive update)
- `adaptive=false`: adaptive penalty strategy:
  - `false`: fixed ρ throughout
  - `true` or `:spectral`: AADMM spectral (adaptive BB + predictive dual + cross-block balancing)
  - `:balanced`: residual balancing (adjusts ρ based on primal/dual residual ratio)
  - `:aradmm`: ARADMM (spectral penalty + over-relaxation parameter γ)
- `adp_freq=2`: adaptive update frequency (every N outer iterations)
- `adp_start=1`: first iteration for adaptive updates
- `maxit=300`: maximum BSDMM outer iterations
- `x_maxiter=50`: maximum VMLMB iterations per x-update (use small values, e.g. 5)
- `verb=true`: print convergence info every 10 iterations
- `history=false`: when `true`, returns `(image, hist)` where `hist` is a NamedTuple
  with fields `chi2`, `obj`, `r_norm`, `s_norm`, `rho` (Dict), `gamma`

# Returns
- Image with best objective value (same shape as `x_init`)
- If `history=true`: `(image, hist)` tuple

# References
- Moolekamp & Melchior (arXiv:1708.09066) — BSDMM algorithm
- Xu, Taylor & Goldstein — ARADMM / AMADMM adaptive penalty
- Chambolle (2004) — TV proximal operator

See also: [`reconstruct`](@ref), [`prox_tv`](@ref), [`prox_l2smooth`](@ref),
[`prox_group_sparsity`](@ref), [`prox_grouptv`](@ref)
"""
function reconstruct_bsdmm(x_init::AbstractMatrix{<:AbstractFloat},
                           data::AbstractMatrix{<:OIdata},
                           ft::AbstractMatrix; kwargs...)
    @assert size(data) == (1, 1) "reconstruct_bsdmm with 2D image requires (1,1) data"
    return reconstruct_bsdmm(x_init, data[1, 1], ft, data; kwargs...)
end

# ─── Polychromatic dispatch (4D image) ───────────────────────────────────────

function reconstruct_bsdmm(x_init::AbstractArray{<:AbstractFloat,4},
                           data::AbstractMatrix{<:OIdata},
                           ft::AbstractMatrix; kwargs...)
    nwav, nepoch = size(data)
    @assert nepoch == 1 "reconstruct_bsdmm: temporal not supported yet"
    @assert size(x_init, 3) == nwav "x_init 3rd dim ($(size(x_init,3))) must match nwav=$nwav"
    @assert size(x_init, 4) == 1 "x_init 4th dim must be 1"
    return _reconstruct_bsdmm_poly(x_init, data, ft; kwargs...)
end

function reconstruct_bsdmm(x_init, d::OIdata, ft, data;
                           mu_reg=0.0, mu_cen=0.0,
                           reg_type=:tv, tv_niter=50,
                           rho_init=10.0, rho_min=1e-3, rho_max=1e4,
                           adaptive=false,
                           adp_freq=2, adp_start=1,
                           maxit=300, x_maxiter=50, verb=true,
                           history=false)
    nx = size(x_init, 1)
    npix = nx * nx

    # Normalize adaptive keyword
    adp_mode = if adaptive === true || adaptive == :spectral
        :spectral
    elseif adaptive === false
        :none
    elseif adaptive == :balanced
        :balanced
    elseif adaptive == :aradmm
        :aradmm
    else
        error("Unknown adaptive mode: $adaptive. Use false, true, :spectral, :balanced, or :aradmm")
    end

    # Initialize image: positive, unit flux
    x = copy(vec(x_init))
    x .= max.(0, x)
    x ./= sum(x)

    # ── Build blocks (only non-smooth regularizers) ──
    blocks = BSDMMBlock[]

    # Regularization (TV or L2smooth)
    use_reg = mu_reg > 0
    if use_reg
        if reg_type == :tv
            push!(blocks, BSDMMBlock(:reg, npix, make_prox_tv(nx, mu_reg; niter=tv_niter);
                                     rho=rho_init, rho_min=rho_min, rho_max=rho_max))
        elseif reg_type == :l2smooth
            push!(blocks, BSDMMBlock(:reg, npix, make_prox_l2smooth(nx, mu_reg);
                                     rho=rho_init, rho_min=rho_min, rho_max=rho_max))
        end
    end

    # Centering
    use_cen = mu_cen > 0
    if use_cen
        push!(blocks, BSDMMBlock(:cen, npix, make_prox_centering(nx, mu_cen);
                                 rho=rho_init, rho_min=rho_min, rho_max=rho_max))
    end

    # Initialize z and adaptive storage from x
    for blk in blocks
        blk.z .= x
        blk.z_prev .= x
        blk.z_adp .= x
        blk.lambda_adp .= blk.rho .* blk.u
        blk.lhat_adp .= blk.rho .* blk.u
    end

    # ── Adaptive state ──
    gamma = 1.0
    x_adp = copy(x)

    # ── History tracking ──
    if history
        hist = (chi2=Float64[], obj=Float64[], r_norm=Float64[], s_norm=Float64[],
                rho=Dict(blk.name => Float64[] for blk in blocks), gamma=Float64[])
    end

    # ── Tracking ──
    best_obj = Inf
    best_x = copy(x)

    # Built once: the x-update runs a full VMLMB per ADMM iteration, so a workspace allocated
    # per criterion evaluation would be re-allocated thousands of times. Only the 1x1 case
    # reaches the single-cell kernel; anything wider dispatches to the 4-D path instead.
    chi2_cache = size(data) == (1, 1) ? ImageChi2Cache(ft[1], data[1]) : nothing

    for k in 1:maxit
        # ── z-updates (parallelizable) ──
        for blk in blocks
            blk.z_prev .= blk.z
            blk.u_prev .= blk.u
            if adp_mode == :aradmm && gamma > 1.0 + 1e-10
                v = gamma .* x .+ (1.0 - gamma) .* blk.z_prev .+ blk.u
            else
                v = x .+ blk.u
            end
            blk.prox!(blk.z, v, blk.rho)
        end

        # ── x-update: min_x chi2(H(x/s)) + Σ (ρ_i/2)||x - z_i + u_i||² ──
        # Precompute weighted target (valid because all L_i = I)
        rho_sum = sum(blk.rho for blk in blocks)
        target = zeros(npix)
        for blk in blocks
            target .+= blk.rho .* (blk.z .- blk.u)
        end
        if rho_sum > 0
            target ./= rho_sum
        end

        function crit_bsdmm(x_val, g_val)
            x_2d = reshape(x_val, nx, nx)
            g_2d = reshape(g_val, nx, nx)
            # chi2 with chain-rule normalization (flux-invariant)
            fval = chi2_cache === nothing ?
                image_to_chi2_fg(x_2d, g_2d, ft, data; verb=false) :
                image_to_chi2_fg(x_2d, g_2d, ft, data; verb=false, cache=chi2_cache)
            # Augmented Lagrangian penalty toward proximal targets
            if rho_sum > 0
                r = x_val .- target
                fval += (rho_sum / 2) * dot(r, r)
                g_val .+= rho_sum .* r
            end
            return fval
        end

        x .= OptimPackNextGen.vmlmb(crit_bsdmm, x, verb=false,
                                     maxiter=x_maxiter, lower=0, blmvm=false)

        # Re-normalize after x-update to maintain unit flux
        sx = sum(x)
        if sx > 0
            x ./= sx
        end

        # ── Dual updates ──
        for blk in blocks
            if adp_mode == :aradmm && gamma > 1.0 + 1e-10
                blk.u .+= gamma .* x .+ (1.0 - gamma) .* blk.z_prev .- blk.z
            else
                blk.u .+= x .- blk.z
            end
        end

        # ── Adaptive penalty update ──
        if adp_mode != :none && k >= adp_start && mod(k, adp_freq) == 0
            if adp_mode == :spectral || adp_mode == :aradmm
                cg = adp_mode == :aradmm
                # First pass: estimate curvature per block
                estimates = [spectral_update!(blk, x, x_adp; compute_gamma=cg) for blk in blocks]
                # Cross-block fallback (AMADMM): non-estimable blocks use max of estimable
                rho_ests = [e[1] for e in estimates]
                good = filter(!isnothing, rho_ests)
                if !isempty(good) && any(isnothing, rho_ests)
                    fb = maximum(good)
                    for (i, blk) in enumerate(blocks)
                        if isnothing(rho_ests[i])
                            spectral_update!(blk, x, x_adp; compute_gamma=cg, fallback_rho=fb)
                        end
                    end
                end
                if adp_mode == :aradmm
                    gammas = [e[2] for e in estimates]
                    gamma = clamp(sum(gammas) / length(gammas), 1.0, 1.9)
                end
            elseif adp_mode == :balanced
                for blk in blocks
                    balanced_update!(blk, x)
                end
            end
            x_adp .= x
        end

        # ── Convergence monitoring ──
        r_norm_sq = 0.0
        for blk in blocks
            r_norm_sq += sum((x .- blk.z).^2)
        end
        s_norm_sq = 0.0
        for blk in blocks
            s_norm_sq += blk.rho^2 * sum((blk.z .- blk.z_prev).^2)
        end
        r_norm = sqrt(r_norm_sq)
        s_norm = sqrt(s_norm_sq)

        # Track best objective
        x_2d = reshape(x, nx, nx)
        chi2_x = image_to_chi2(x_2d, ft, data; verb=false)
        reg_val = use_reg ? (reg_type == :tv ? tv_norm(x_2d) : l2smooth_norm(x_2d)) : 0.0
        obj_x = chi2_x + (use_reg ? mu_reg * reg_val : 0.0)
        if obj_x < best_obj
            best_obj = obj_x
            best_x .= vec(x)
        end

        # ── History tracking ──
        if history
            push!(hist.chi2, chi2_x)
            push!(hist.obj, obj_x)
            push!(hist.r_norm, r_norm)
            push!(hist.s_norm, s_norm)
            for blk in blocks
                push!(hist.rho[blk.name], blk.rho)
            end
            push!(hist.gamma, gamma)
        end

        # ── Display ──
        if verb && (mod(k, 10) == 0 || k <= 5)
            @printf("Iter %3d | ", k)
            image_to_chi2(x_2d, ft, data; verb=true)
            rho_str = join([@sprintf("%s=%.1e", blk.name, blk.rho) for blk in blocks], " ")
            gamma_str = adp_mode == :aradmm ? @sprintf("  γ=%.3f", gamma) : ""
            @printf("  obj=%.2f  |r|=%.2e |s|=%.2e  reg=%.4f  flux=%.4f  ρ: %s%s\n",
                    obj_x, r_norm, s_norm, reg_val, sum(x), rho_str, gamma_str)
        end
    end

    # Under `verb`, like every other line this prints. It used to be unconditional, which made
    # `verb = false` a half-promise: a caller that asked for silence still got one line, and a
    # precompile workload or a batch script had to redirect stdout to get the silence it asked
    # for.
    verb && @printf("Best objective: %.2f\n", best_obj)
    result = reshape(best_x, nx, nx)
    return history ? (result, hist) : result
end

# ═══════════════════════════════════════════════════════════════════════════════
# Polychromatic BSDMM
# ═══════════════════════════════════════════════════════════════════════════════

function _reconstruct_bsdmm_poly(x_init, data, ft;
                                  mu_tv=0.0, mu_group=0.0, mu_cen=0.0,
                                  reg_type=:tv, group_type=:sparsity,
                                  tv_niter=50,
                                  rho_init=10.0, rho_min=1e-3, rho_max=1e4,
                                  adaptive=false, adp_freq=2, adp_start=1,
                                  maxit=300, x_maxiter=50, verb=true,
                                  history=false)
    nx = size(x_init, 1)
    nwav = size(x_init, 3)
    npix = nx * nx
    npix_total = npix * nwav

    # Normalize adaptive keyword
    adp_mode = if adaptive === true || adaptive == :spectral
        :spectral
    elseif adaptive === false
        :none
    elseif adaptive == :balanced
        :balanced
    elseif adaptive == :aradmm
        :aradmm
    else
        error("Unknown adaptive mode: $adaptive")
    end

    # Initialize: positive, per-channel unit flux
    x = copy(vec(x_init[:,:,:,1]))
    x .= max.(0, x)
    for w in 1:nwav
        idx = (w-1)*npix+1 : w*npix
        s = sum(@view(x[idx]))
        if s > 0; @view(x[idx]) ./= s; end
    end

    # ── Build blocks ──
    blocks = BSDMMBlock[]

    use_tv = mu_tv > 0
    if use_tv
        if reg_type == :tv
            push!(blocks, BSDMMBlock(:tv, npix_total,
                  make_prox_tv_poly(nx, nwav, mu_tv; niter=tv_niter);
                  rho=rho_init, rho_min=rho_min, rho_max=rho_max))
        elseif reg_type == :l2smooth
            push!(blocks, BSDMMBlock(:tv, npix_total,
                  make_prox_l2smooth_poly(nx, nwav, mu_tv);
                  rho=rho_init, rho_min=rho_min, rho_max=rho_max))
        end
    end

    use_group = mu_group > 0
    if use_group
        if group_type == :sparsity
            push!(blocks, BSDMMBlock(:group, npix_total,
                  make_prox_group_sparsity(nx, nwav, mu_group);
                  rho=rho_init, rho_min=rho_min, rho_max=rho_max))
        elseif group_type == :tv
            push!(blocks, BSDMMBlock(:group, npix_total,
                  make_prox_grouptv(nx, nwav, mu_group; niter=tv_niter);
                  rho=rho_init, rho_min=rho_min, rho_max=rho_max))
        end
    end

    use_cen = mu_cen > 0
    if use_cen
        push!(blocks, BSDMMBlock(:cen, npix_total,
              make_prox_centering_poly(nx, nwav, mu_cen);
              rho=rho_init, rho_min=rho_min, rho_max=rho_max))
    end

    # Initialize z and adaptive storage from x
    for blk in blocks
        blk.z .= x
        blk.z_prev .= x
        blk.z_adp .= x
        blk.lambda_adp .= blk.rho .* blk.u
        blk.lhat_adp .= blk.rho .* blk.u
    end

    # ── Adaptive state ──
    gamma = 1.0
    x_adp = copy(x)

    # ── History tracking ──
    if history
        hist = (chi2=Float64[], obj=Float64[], r_norm=Float64[], s_norm=Float64[],
                rho=Dict(blk.name => Float64[] for blk in blocks), gamma=Float64[])
    end

    # ── Tracking ──
    best_obj = Inf
    best_x = copy(x)

    caches = [ImageChi2Cache(ft[w,t], data[w,t]) for w in axes(data,1), t in axes(data,2)]

    for k in 1:maxit
        # ── z-updates ──
        for blk in blocks
            blk.z_prev .= blk.z
            blk.u_prev .= blk.u
            if adp_mode == :aradmm && gamma > 1.0 + 1e-10
                v = gamma .* x .+ (1.0 - gamma) .* blk.z_prev .+ blk.u
            else
                v = x .+ blk.u
            end
            blk.prox!(blk.z, v, blk.rho)
        end

        # ── x-update ──
        rho_sum = sum(blk.rho for blk in blocks)
        target = zeros(npix_total)
        for blk in blocks
            target .+= blk.rho .* (blk.z .- blk.u)
        end
        if rho_sum > 0
            target ./= rho_sum
        end

        function crit_bsdmm_poly(x_val, g_val)
            x_4d = reshape(x_val, nx, nx, nwav, 1)
            g_4d = reshape(g_val, nx, nx, nwav, 1)
            fval = image_to_chi2_fg(x_4d, g_4d, ft, data; verb=false, caches=caches)
            if rho_sum > 0
                r = x_val .- target
                fval += (rho_sum / 2) * dot(r, r)
                g_val .+= rho_sum .* r
            end
            return fval
        end

        x .= OptimPackNextGen.vmlmb(crit_bsdmm_poly, x, verb=false,
                                     maxiter=x_maxiter, lower=0, blmvm=false)

        # Per-channel flux renormalization
        for w in 1:nwav
            idx = (w-1)*npix+1 : w*npix
            sx = sum(@view(x[idx]))
            if sx > 0; @view(x[idx]) ./= sx; end
        end

        # ── Dual updates ──
        for blk in blocks
            if adp_mode == :aradmm && gamma > 1.0 + 1e-10
                blk.u .+= gamma .* x .+ (1.0 - gamma) .* blk.z_prev .- blk.z
            else
                blk.u .+= x .- blk.z
            end
        end

        # ── Adaptive penalty update ──
        if adp_mode != :none && k >= adp_start && mod(k, adp_freq) == 0
            if adp_mode == :spectral || adp_mode == :aradmm
                cg = adp_mode == :aradmm
                estimates = [spectral_update!(blk, x, x_adp; compute_gamma=cg) for blk in blocks]
                rho_ests = [e[1] for e in estimates]
                good = filter(!isnothing, rho_ests)
                if !isempty(good) && any(isnothing, rho_ests)
                    fb = maximum(good)
                    for (i, blk) in enumerate(blocks)
                        if isnothing(rho_ests[i])
                            spectral_update!(blk, x, x_adp; compute_gamma=cg, fallback_rho=fb)
                        end
                    end
                end
                if adp_mode == :aradmm
                    gammas = [e[2] for e in estimates]
                    gamma = clamp(sum(gammas) / length(gammas), 1.0, 1.9)
                end
            elseif adp_mode == :balanced
                for blk in blocks
                    balanced_update!(blk, x)
                end
            end
            x_adp .= x
        end

        # ── Convergence monitoring ──
        r_norm_sq = 0.0
        for blk in blocks
            r_norm_sq += sum((x .- blk.z).^2)
        end
        s_norm_sq = 0.0
        for blk in blocks
            s_norm_sq += blk.rho^2 * sum((blk.z .- blk.z_prev).^2)
        end
        r_norm = sqrt(r_norm_sq)
        s_norm = sqrt(s_norm_sq)

        # Track best objective
        x_3d = reshape(x, nx, nx, nwav)
        chi2_x = image_to_chi2(reshape(x, nx, nx, nwav, 1), ft, data; verb=false)
        tv_val = use_tv ? sum(tv_norm(@view(x_3d[:,:,w])) for w in 1:nwav) : 0.0
        group_val = use_group ? (group_type == :sparsity ? group_sparsity_norm(x_3d) : grouptv_norm(x_3d)) : 0.0
        obj_x = chi2_x + mu_tv * tv_val + mu_group * group_val
        if obj_x < best_obj
            best_obj = obj_x
            best_x .= vec(x)
        end

        # ── History tracking ──
        if history
            push!(hist.chi2, chi2_x)
            push!(hist.obj, obj_x)
            push!(hist.r_norm, r_norm)
            push!(hist.s_norm, s_norm)
            for blk in blocks
                push!(hist.rho[blk.name], blk.rho)
            end
            push!(hist.gamma, gamma)
        end

        # ── Display ──
        if verb && (mod(k, 10) == 0 || k <= 5)
            @printf("Iter %3d | ", k)
            image_to_chi2(reshape(x, nx, nx, nwav, 1), ft, data; verb=true)
            rho_str = join([@sprintf("%s=%.1e", blk.name, blk.rho) for blk in blocks], " ")
            gamma_str = adp_mode == :aradmm ? @sprintf("  γ=%.3f", gamma) : ""
            @printf("  obj=%.2f  |r|=%.2e |s|=%.2e  tv=%.4f  grp=%.4f  ρ: %s%s\n",
                    obj_x, r_norm, s_norm, tv_val, group_val, rho_str, gamma_str)
        end
    end

    verb && @printf("Best objective: %.2f\n", best_obj)
    result = reshape(best_x, nx, nx, nwav, 1)
    return history ? (result, hist) : result
end
