# ═══════════════════════════════════════════════════════════════════════════════
# ADMM image reconstruction for optical interferometry
# ═══════════════════════════════════════════════════════════════════════════════

# Centering regularizer for ADMM (divides gradient by flux, unlike reg_centering in oichi2.jl)
function admm_reg_centering(x_2d, g_vec)
    nx = size(x_2d, 1)
    flux = sum(x_2d)
    if flux <= 0; g_vec .= 0; return 0.0; end
    xvals = collect(1:nx)
    cx = sum(xvals' * x_2d) / flux
    cy = sum(x_2d * xvals) / flux
    center = (nx + 1) / 2.0
    f = (cx - center)^2 + (cy - center)^2
    xx = [mod(i - 1, nx) + 1 for i in 1:nx*nx]
    yy = [div(i - 1, nx) + 1 for i in 1:nx*nx]
    g_vec .= (2 * (cx - center) .* xx .+ 2 * (cy - center) .* yy) ./ flux
    return f
end

# NFFT helpers (internal)
function nfft_forward(plan, x_vec, nx)
    plan * ComplexF64.(reshape(x_vec, nx, nx))
end

function nfft_adjoint(plan, v, nx)
    real(vec(adjoint(plan) * v))
end

"""
    reconstruct_admm(x_init, d::OIdata, nfft_plan, ft, data, nx; kwargs...)

ADMM image reconstruction with V2, T3phi, and regularization proximal operators.

Uses NFFT for the forward model and VMLMB for the x-update with flux normalization.
The penalty parameter `rho` must be large enough (~10000) for convergence.

# Arguments
- `x_init`: initial image (nx × nx matrix)
- `d::OIdata`: interferometric data for the single (wavelength, epoch) cell
- `nfft_plan`: NFFT plan (first element of `setup_ft` output cell)
- `ft`: full Fourier transform setup from `setup_ft`
- `data`: full data matrix from `readoifits`
- `nx::Int`: image size in pixels

# Keyword arguments
- `rho_v2=10000.0`: ADMM penalty for V2 block
- `rho_t3=10000.0`: ADMM penalty for T3phi block
- `rho_reg=10.0`: ADMM penalty for regularization block
- `mu_reg=0.0`: regularization weight (0 = no regularization)
- `mu_cen=0.0`: centering regularization weight
- `reg_type=:tv`: regularization type (`:tv` or `:l2smooth`)
- `tv_niter=50`: number of Chambolle iterations for TV proximal
- `use_t3=true`: include closure phases
- `maxit=300`: maximum ADMM iterations
- `x_maxiter=50`: maximum VMLMB iterations per x-update
- `verb=true`: print convergence info
"""
function reconstruct_admm(x_init, d::OIdata, nfft_pl, ft, data, nx::Int;
                          rho_v2=10000.0, rho_t3=10000.0, rho_reg=10.0,
                          mu_reg=0.0, mu_cen=0.0, reg_type=:tv, tv_niter=50,
                          use_t3=true,
                          maxit=300, x_maxiter=50, verb=true)
    npix = nx * nx
    nuv = d.nuv
    use_reg = mu_reg > 0

    # Compute T3 UV indices
    indx_t3 = sort(unique(vcat(d.indx_t3_1, d.indx_t3_2, d.indx_t3_3)))

    x = copy(vec(x_init))
    x .= max.(0, x); x ./= sum(x)

    Hx = nfft_forward(nfft_pl, x, nx)

    # V2 block
    v2_z = copy(Hx)
    u_v2 = zeros(ComplexF64, nuv)

    # T3phi block
    t3_z = copy(Hx)
    u_t3 = zeros(ComplexF64, nuv)

    # Regularization block (image domain)
    reg_z = reshape(copy(x), nx, nx)
    u_reg = zeros(nx, nx)

    best_L = Inf
    best_x = copy(x)

    for k in 1:maxit
        # ── V2 proximal ──
        v2_z = prox_v2(Hx .+ u_v2, d.v2, d.v2_err, d.indx_v2, rho_v2)

        # ── T3phi proximal ──
        if use_t3
            t3_z = prox_t3phi(Hx .+ u_t3, d.t3phi, d.t3phi_err,
                              indx_t3, d.indx_t3_1, d.indx_t3_2, d.indx_t3_3,
                              rho_t3)
        end

        # ── Regularization proximal ──
        if use_reg
            x_img = reshape(x, nx, nx)
            if reg_type == :tv
                reg_z .= prox_tv(x_img .+ u_reg, mu_reg / rho_reg; niter=tv_niter)
            elseif reg_type == :l2smooth
                reg_z .= prox_l2smooth(x_img .+ u_reg, mu_reg / rho_reg)
            end
        end

        # ── x-update: VMLMB with normalization ──
        target_v2 = v2_z .- u_v2
        target_t3 = t3_z .- u_t3
        target_reg = use_reg ? vec(reg_z .- u_reg) : zeros(0)

        function crit_x(x_val, g)
            flux = sum(x_val)
            if flux <= 0; g .= 0; return 0.0; end
            y = x_val ./ flux
            Hx_val = nfft_forward(nfft_pl, y, nx)

            # Data terms
            res_v2 = Hx_val .- target_v2
            f = (rho_v2 / 2.0) * real(dot(res_v2, res_v2))
            g_vis = rho_v2 .* res_v2

            if use_t3
                res_t3 = Hx_val .- target_t3
                f += (rho_t3 / 2.0) * real(dot(res_t3, res_t3))
                g_vis .+= rho_t3 .* res_t3
            end

            g_raw = nfft_adjoint(nfft_pl, g_vis, nx)

            # Regularization penalty term
            if use_reg
                res_reg = y .- target_reg
                f += (rho_reg / 2.0) * dot(res_reg, res_reg)
                g_raw .+= rho_reg .* res_reg
            end

            # Centering regularizer on normalized image (before chain rule)
            if mu_cen > 0
                g_cen = zeros(npix)
                f_cen = admm_reg_centering(reshape(y, nx, nx), g_cen)
                f += mu_cen * f_cen
                g_raw .+= mu_cen .* g_cen
            end

            # Chain rule for x/sum(x) normalization
            g .= (g_raw .- dot(g_raw, x_val) / flux) ./ flux
            return f
        end

        x .= OptimPackNextGen.vmlmb(crit_x, x, verb=false, lower=0,
                                     maxiter=x_maxiter, blmvm=false)

        flux = sum(x); if flux > 0; x ./= flux; end
        Hx .= nfft_forward(nfft_pl, x, nx)

        # ── Dual updates ──
        u_v2 .+= Hx .- v2_z
        if use_t3
            u_t3 .+= Hx .- t3_z
        end
        if use_reg
            u_reg .+= reshape(x, nx, nx) .- reg_z
        end

        # ── Track best objective = chi2(x) + mu_reg * reg(x) ──
        x_2d = reshape(x, nx, nx)
        chi2_x = image_to_chi2(x_2d, ft, data; verb=false)
        reg_val = use_reg ? (reg_type == :tv ? tv_norm(x_2d) : l2smooth_norm(x_2d)) : 0.0
        obj_x = chi2_x + (use_reg ? mu_reg * reg_val : 0.0)
        if obj_x < best_L
            best_L = obj_x
            best_x .= x
        end

        # ── Display ──
        if verb && (mod(k, 10) == 0 || k <= 5)
            @printf("Iter %3d | ", k)
            image_to_chi2(x_2d, ft, data; verb=true)
            r_v2 = norm(Hx .- v2_z)
            r_t3 = use_t3 ? norm(Hx .- t3_z) : 0.0
            r_reg = use_reg ? norm(x_2d .- reg_z) : 0.0
            @printf("  obj=%.2f  |r_v2|=%.2e |r_t3|=%.2e |r_reg|=%.2e  reg(x)=%.4f\n",
                    obj_x, r_v2, r_t3, r_reg, reg_val)
        end
    end
    @printf("Best objective: %.2f\n", best_L)
    return reshape(best_x, nx, nx)
end
