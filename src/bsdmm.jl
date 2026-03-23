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
# Currently all L_i = I (identity). Future: transspectral operators for
# polychromatic (N images coupled across wavelength channels).
#
# Per-block penalty ρ_i adapted via Barzilai-Borwein spectral estimates
# with correlation safeguards (ARADMM-style).
#
# Reference: Moolekamp & Melchior (arXiv:1708.09066), Xu & Taylor (ARADMM)
# ═══════════════════════════════════════════════════════════════════════════════

# ─── Block data structure ────────────────────────────────────────────────────

mutable struct BSDMMBlock{T<:AbstractFloat}
    name::Symbol
    z::Vector{T}         # auxiliary variable
    u::Vector{T}         # scaled dual variable
    z_prev::Vector{T}    # previous z (for spectral update)
    u_prev::Vector{T}    # previous u (for spectral update)
    prox!::Any           # prox!(z, v, rho) modifies z in place
    rho::T               # penalty parameter
    rho_min::T
    rho_max::T
end

function BSDMMBlock(name::Symbol, n::Int, prox!;
                    rho=10.0, rho_min=1e-3, rho_max=1e8)
    BSDMMBlock{Float64}(name,
        zeros(n), zeros(n), zeros(n), zeros(n),
        prox!, rho, rho_min, rho_max)
end

# ─── Spectral rho update (Barzilai-Borwein with ARADMM safeguards) ──────────

function spectral_rho_update!(blk::BSDMMBlock; ε_cor=0.2, smooth=0.9)
    s = blk.z .- blk.z_prev     # Δz
    y = blk.u .- blk.u_prev     # Δu
    ss = dot(s, s)
    sy = dot(s, y)

    if ss < 1e-30
        return
    end

    # Correlation safeguard
    yy = dot(y, y)
    if sy < ε_cor * sqrt(ss * yy) + 1e-30
        return
    end

    # BB spectral estimate
    rho_bb = clamp(sy / ss, blk.rho_min, blk.rho_max)

    # Smooth update
    rho_old = blk.rho
    blk.rho = smooth * rho_old + (1 - smooth) * rho_bb

    # Rescale dual variable to maintain λ = ρ u
    blk.u .*= rho_old / blk.rho
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

Currently supports monochromatic (1×1) data only.

# Arguments
- `x_init`: starting image (`nx × nx` matrix)
- `data`: interferometric data (`Matrix{OIdata}`, size `(1,1)`)
- `ft`: Fourier transform plans from `setup_ft`

# Keyword arguments
- `mu_reg=0.0`: regularization weight (0 = no regularization)
- `mu_cen=0.0`: centering weight (0 = no centering)
- `reg_type=:tv`: regularization type (`:tv` or `:l2smooth`)
- `tv_niter=50`: Chambolle iterations for TV proximal operator
- `rho_init=10.0`: initial ADMM penalty parameter for all blocks
- `rho_min=1e-3`: minimum penalty (for adaptive update)
- `rho_max=1e4`: maximum penalty (for adaptive update)
- `adaptive=true`: enable per-block Barzilai-Borwein spectral ρ adaptation
- `maxit=300`: maximum BSDMM outer iterations
- `x_maxiter=50`: maximum VMLMB iterations per x-update (use small values, e.g. 5)
- `verb=true`: print convergence info every 10 iterations

# Returns
- `nx × nx` image (positive, unit flux) with best objective value

# References
- Moolekamp & Melchior (arXiv:1708.09066) — BSDMM algorithm
- Xu & Taylor — ARADMM adaptive penalty
- Chambolle (2004) — TV proximal operator

See also: [`reconstruct`](@ref), [`prox_tv`](@ref), [`prox_l2smooth`](@ref)
"""
function reconstruct_bsdmm(x_init::AbstractMatrix{<:AbstractFloat},
                           data::AbstractMatrix{<:OIdata},
                           ft::AbstractMatrix; kwargs...)
    @assert size(data) == (1, 1) "reconstruct_bsdmm currently supports mono (1,1) data only"
    return reconstruct_bsdmm(x_init, data[1, 1], ft, data; kwargs...)
end

function reconstruct_bsdmm(x_init, d::OIdata, ft, data;
                           mu_reg=0.0, mu_cen=0.0,
                           reg_type=:tv, tv_niter=50,
                           rho_init=10.0, rho_min=1e-3, rho_max=1e4,
                           adaptive=true,
                           maxit=300, x_maxiter=50, verb=true)
    nx = size(x_init, 1)
    npix = nx * nx

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

    # Initialize z from x
    for blk in blocks
        blk.z .= x
        blk.z_prev .= x
    end

    # ── Tracking ──
    best_obj = Inf
    best_x = copy(x)

    for k in 1:maxit
        # ── z-updates (parallelizable) ──
        for blk in blocks
            blk.z_prev .= blk.z
            blk.u_prev .= blk.u
            v = x .+ blk.u
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
            fval = image_to_chi2_fg(x_2d, g_2d, ft, data; verb=false)
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
            blk.u .+= x .- blk.z
        end

        # ── Adaptive rho (spectral BB) ──
        if adaptive && k > 1
            for blk in blocks
                spectral_rho_update!(blk)
            end
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

        # ── Display ──
        if verb && (mod(k, 10) == 0 || k <= 5)
            @printf("Iter %3d | ", k)
            image_to_chi2(x_2d, ft, data; verb=true)
            rho_str = join([@sprintf("%s=%.1e", blk.name, blk.rho) for blk in blocks], " ")
            @printf("  obj=%.2f  |r|=%.2e |s|=%.2e  reg=%.4f  flux=%.4f  ρ: %s\n",
                    obj_x, r_norm, s_norm, reg_val, sum(x), rho_str)
        end
    end

    @printf("Best objective: %.2f\n", best_obj)
    return reshape(best_x, nx, nx)
end
