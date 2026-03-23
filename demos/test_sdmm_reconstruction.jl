#
# ADMM image reconstruction for optical interferometry
#
# ADMM with V2, T3phi, and TV proximal operators:
#   v2    = prox_v2(Hx + u_v2, ρ_v2)              — V2 fit
#   t3    = prox_t3phi(Hx + u_t3, ρ_t3)           — T3phi fit
#   s     = prox_{μ/ρ_tv · TV}(x + u_tv)          — isotropic TV
#   x     = VMLMB min Σ ρ_i||H(x/s)-t_i||² + ρ_tv||x/s - t_tv||²  s.t. x≥0
#   dual updates
#
# Uses NFFT. ρ must be large enough (~10000) for convergence.

using OITOOLS
using OptimPackNextGen
using LinearAlgebra, Printf, FFTW

# ═══════════════════════════════════════════════════════════════════════════════
# Proximal operators
# ═══════════════════════════════════════════════════════════════════════════════
 function realcuberoot(xx::Real)
     sign(xx) * abs(xx)^(1//3)
 end


function paintercubicroots(c::Real, d::Real)
    q = c / 3.0
    r = d / 2.0
    discriminant = q^3 + r^2
    if discriminant >= 0
        s = realcuberoot(r + sqrt(discriminant))
        t = realcuberoot(r - sqrt(discriminant))
        root = [s + t]
    else
        rho = sqrt(r^2 - discriminant)
        cubeRootrho = realcuberoot(rho)
        thetaOn3 = acos(r / rho) / 3
        crRhoCosThetaOn3 = cubeRootrho * cos(thetaOn3)
        crRhoSinThetaOn3 = cubeRootrho * sin(thetaOn3)
        root = zeros(3)
        root[1] = 2 * crRhoCosThetaOn3
        root[2] = -crRhoCosThetaOn3 - sqrt(3) * crRhoSinThetaOn3
        root[3] = -crRhoCosThetaOn3 + sqrt(3) * crRhoSinThetaOn3
    end
    return root
end

function prox_v2(z0, v2_data, v2_err, indx_v2, rho_v2)
    z = copy(z0)
    for n in eachindex(v2_data)
        idx = indx_v2[n]
        v0 = abs(z[idx])
        ang = angle(z[idx])
        sig2 = v2_err[n]^2
        tmp1 = sig2 * rho_v2 / 4.0
        c = tmp1 - v2_data[n]
        d = tmp1 * v0
        sol = max.(0.0, paintercubicroots(c, d))
        cst = (1.0 / sig2) .* (sol.^2 .- v2_data[n]).^2 .+
              (rho_v2 / 2.0) .* (sol .- v0).^2
        _, b = findmin(cst)
        z[idx] = sol[b] * cis(ang)
    end
    return z
end

# T3phi proximal (Haniff least-squares via VMLMB on phases)
function mod2pi_wrap(x)
    mod.(mod.(x .+ pi, 2pi) .+ 2pi, 2pi) .- pi
end

function crit_t3phi_haniff(g, phi, phi0, V0, rho, t3phi_rad, t3phi_err_rad,
                            indx_t3, indx_t3_1, indx_t3_2, indx_t3_3)
    Vopt = max.(0, V0 .* cos.(phi .- phi0))
    res = mod2pi_wrap((phi[indx_t3_1] .+ phi[indx_t3_2] .+ phi[indx_t3_3]) .- t3phi_rad) ./ t3phi_err_rad
    f = sum(res.^2) + 0.5 * rho * norm(Vopt .* cis.(phi .- phi0) .- V0)^2
    g .= 0
    g[indx_t3] .= rho .* Vopt[indx_t3] .* V0[indx_t3] .* sin.(phi[indx_t3] .- phi0[indx_t3])
    dres = 2 .* res ./ t3phi_err_rad
    g[indx_t3_1] .+= dres
    g[indx_t3_2] .+= dres
    g[indx_t3_3] .+= dres
    return f
end

function prox_t3phi(z0, t3phi_deg, t3phi_err_deg,
                     indx_t3, indx_t3_1, indx_t3_2, indx_t3_3, rho_t3phi)
    z = copy(z0)
    V0 = abs.(z)
    phi0 = angle.(z)
    t3phi_rad = t3phi_deg .* (pi / 180.0)
    t3phi_err_rad = t3phi_err_deg .* (pi / 180.0)
    crit = (phi, g) -> crit_t3phi_haniff(g, phi, phi0, V0, rho_t3phi,
                            t3phi_rad, t3phi_err_rad,
                            indx_t3, indx_t3_1, indx_t3_2, indx_t3_3)
    phi_opt = OptimPackNextGen.vmlmb(crit, copy(phi0), verb=false, maxiter=200, blmvm=false)
    z[indx_t3] .= max.(0, V0[indx_t3] .* cos.(phi_opt[indx_t3] .- phi0[indx_t3])) .*
                  cis.(phi_opt[indx_t3])
    return z
end

# ═══════════════════════════════════════════════════════════════════════════════
# Isotropic TV proximal (Chambolle 2004 dual algorithm)
# ═══════════════════════════════════════════════════════════════════════════════

# Forward differences with Neumann BC
function grad_op!(g_h, g_v, x)
    nx = size(x, 1)
    @views g_h[1:nx-1, :] .= x[2:nx, :] .- x[1:nx-1, :]
    g_h[nx, :] .= 0
    @views g_v[:, 1:nx-1] .= x[:, 2:nx] .- x[:, 1:nx-1]
    g_v[:, nx] .= 0
end

# Divergence = negative adjoint of gradient (backward differences)
function div_op!(d, p_h, p_v)
    nx = size(d, 1)
    # horizontal component
    d[1, :] .= p_h[1, :]
    @views d[2:nx-1, :] .= p_h[2:nx-1, :] .- p_h[1:nx-2, :]
    @views d[nx, :] .= .-p_h[nx-1, :]
    # vertical component
    @views d[:, 1] .+= p_v[:, 1]
    @views d[:, 2:nx-1] .+= p_v[:, 2:nx-1] .- p_v[:, 1:nx-2]
    @views d[:, nx] .+= .-p_v[:, nx-1]
end

# prox_{λ·TV}(y) = argmin_x (1/2)||x - y||² + λ TV(x)
# Chambolle's projected gradient ascent on the dual variable p
function prox_tv(y::AbstractMatrix, lambda; niter=100, tau=1.0/8.0)
    nx = size(y, 1)
    p_h = zeros(nx, nx)
    p_v = zeros(nx, nx)
    d = zeros(nx, nx)
    g_h = zeros(nx, nx)
    g_v = zeros(nx, nx)

    nrm = zeros(nx, nx)

    for n in 1:niter
        # d = div(p) - y/λ
        div_op!(d, p_h, p_v)
        d .-= y ./ lambda

        # gradient of d
        grad_op!(g_h, g_v, d)

        # update p + project onto unit ball
        p_h .+= tau .* g_h
        p_v .+= tau .* g_v
        @. nrm = max(1.0, sqrt(p_h^2 + p_v^2))
        p_h ./= nrm
        p_v ./= nrm
    end

    # x = y - λ div(p)
    div_op!(d, p_h, p_v)
    return y .- lambda .* d
end

# ═══════════════════════════════════════════════════════════════════════════════
# Centering regularizer: penalize COG deviation from image center
# ═══════════════════════════════════════════════════════════════════════════════
function local_reg_centering(x_2d, g_vec)
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

# ═══════════════════════════════════════════════════════════════════════════════
# Quadratic smoothness proximal via FFT
# prox_{λ·(1/2)||∇x||²}(y) = (I + λ·L)^{-1} y
# where L is the discrete Laplacian (periodic BC), diagonalized by FFT
# ═══════════════════════════════════════════════════════════════════════════════
function prox_l2smooth(y::AbstractMatrix, lambda)
    nx = size(y, 1)
    # Eigenvalues of discrete Laplacian with periodic BC
    freq = [0:nx-1;]
    eig_h = 2 .- 2 .* cos.(2pi .* freq ./ nx)  # 1D eigenvalues
    eig_2d = eig_h .+ eig_h'                     # 2D eigenvalues (separable)
    # Solve in Fourier domain: x̂ = ŷ / (1 + λ·eigenvalue)
    Y = fft(y)
    X = Y ./ (1 .+ lambda .* eig_2d)
    return real(ifft(X))
end

# L2 smoothness norm: (1/2)||∇x||²
function l2smooth_norm(x_img)
    nx = size(x_img, 1)
    g_h = zeros(nx, nx)
    g_v = zeros(nx, nx)
    grad_op!(g_h, g_v, x_img)
    return 0.5 * sum(g_h.^2 .+ g_v.^2)
end

# Isotropic TV value of image
function tv_norm(x_img)
    nx = size(x_img, 1)
    g_h = zeros(nx, nx)
    g_v = zeros(nx, nx)
    grad_op!(g_h, g_v, x_img)
    return sum(sqrt.(g_h.^2 .+ g_v.^2))
end

# chi2_v2 evaluated at visibility vector z (not at image)
function chi2_v2_at(z, v2_data, v2_err, indx_v2)
    v2_model = abs2.(z[indx_v2])
    return sum(((v2_model .- v2_data) ./ v2_err).^2)
end

# chi2_t3phi evaluated at visibility vector z
function chi2_t3phi_at(z, t3phi_data, t3phi_err, indx_t3_1, indx_t3_2, indx_t3_3)
    t3 = z[indx_t3_1] .* z[indx_t3_2] .* z[indx_t3_3]
    t3phi_model = angle.(t3) .* (180.0 / pi)
    res = mod.(mod.(t3phi_model .- t3phi_data .+ 180.0, 360.0) .+ 360.0, 360.0) .- 180.0
    return sum((res ./ t3phi_err).^2)
end

# ═══════════════════════════════════════════════════════════════════════════════
# NFFT helpers
# ═══════════════════════════════════════════════════════════════════════════════
function nfft_forward(plan, x_vec, nx)
    plan * ComplexF64.(reshape(x_vec, nx, nx))
end

function nfft_adjoint(plan, v, nx)
    real(vec(adjoint(plan) * v))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Load data
# ═══════════════════════════════════════════════════════════════════════════════
oifitsfile = joinpath(@__DIR__, "data", "2004-data1.oifits")
data = readoifits(oifitsfile)
d = data[1,1]

nx = 64
pixsize = 0.2

ft_nfft = setup_ft(data, nx, pixsize)
nfft_plans = ft_nfft[1,1]
nfft_plan = nfft_plans[1]
nuv = d.nuv

println("Data: nv2=$(d.nv2), nt3phi=$(d.nt3phi), nt3amp=$(d.nt3amp), nuv=$(nuv)")

indx_t3 = sort(unique(vcat(d.indx_t3_1, d.indx_t3_2, d.indx_t3_3)))

# ═══════════════════════════════════════════════════════════════════════════════
# ADMM reconstruction
# ═══════════════════════════════════════════════════════════════════════════════
function admm_reconstruct(x_init, d, nfft_pl, ft_nfft, data, nx, indx_t3;
                          rho_v2=10000.0, rho_t3=10000.0, rho_reg=10.0,
                          mu_reg=0.0, mu_cen=0.0, reg_type=:tv, tv_niter=50,
                          use_t3=true,
                          maxit=300, x_maxiter=50, verb=true)
    npix = nx * nx
    nuv = d.nuv
    use_reg = mu_reg > 0

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
                f_cen = local_reg_centering(reshape(y, nx, nx), g_cen)
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

        # ── Compute augmented Lagrangian ──
        L_obj_v2 = chi2_v2_at(v2_z, d.v2, d.v2_err, d.indx_v2)
        L_obj_t3 = use_t3 ? chi2_t3phi_at(t3_z, d.t3phi, d.t3phi_err,
                              d.indx_t3_1, d.indx_t3_2, d.indx_t3_3) : 0.0
        x_2d = reshape(x, nx, nx)
        if use_reg
            if reg_type == :tv
                L_obj_reg = mu_reg * tv_norm(reg_z)
            else
                L_obj_reg = mu_reg * l2smooth_norm(reg_z)
            end
        else
            L_obj_reg = 0.0
        end

        L_aug_v2 = (rho_v2/2) * real(dot(Hx .- v2_z .+ u_v2, Hx .- v2_z .+ u_v2)) -
                   (rho_v2/2) * real(dot(u_v2, u_v2))
        L_aug_t3 = use_t3 ? ((rho_t3/2) * real(dot(Hx .- t3_z .+ u_t3, Hx .- t3_z .+ u_t3)) -
                             (rho_t3/2) * real(dot(u_t3, u_t3))) : 0.0
        L_aug_reg = use_reg ? ((rho_reg/2) * norm(x_2d .- reg_z .+ u_reg)^2 -
                               (rho_reg/2) * norm(u_reg)^2) : 0.0

        L_total = L_obj_v2 + L_obj_t3 + L_obj_reg + L_aug_v2 + L_aug_t3 + L_aug_reg

        # Track best objective = chi2(x) + mu_reg * reg(x)
        chi2_x = image_to_chi2(x_2d, ft_nfft, data; verb=false)
        reg_val = use_reg ? (reg_type == :tv ? tv_norm(x_2d) : l2smooth_norm(x_2d)) : 0.0
        obj_x = chi2_x + (use_reg ? mu_reg * reg_val : 0.0)
        if obj_x < best_L
            best_L = obj_x
            best_x .= x
        end

        # ── Display ──
        if verb && (mod(k, 10) == 0 || k <= 5)
            @printf("Iter %3d | ", k)
            image_to_chi2(x_2d, ft_nfft, data; verb=true)
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

# ═══════════════════════════════════════════════════════════════════════════════
# Run
# ═══════════════════════════════════════════════════════════════════════════════
x_start = gaussian2d(nx, nx, nx / 6)

println("\n=== Starting chi2 ===")
image_to_chi2(x_start, ft_nfft, data; verb=true)

# V2 + T3phi + TV regularization
println("\n=== V2 + T3phi + TV ADMM ===")
x_tv = admm_reconstruct(x_start, d, nfft_plan, ft_nfft, data, nx, indx_t3;
                          rho_v2=10000.0, rho_t3=10000.0,
                          rho_reg=10000.0, mu_reg=1e5, mu_cen=1e2,
                          reg_type=:tv, tv_niter=50,
                          use_t3=true,
                          maxit=300, x_maxiter=50)
println("\nV2+T3phi+TV final:")
image_to_chi2(x_tv, ft_nfft, data; verb=true)

# Save ADMM result as PNG
using PyPlot
imdisp(x_tv, pixsize=pixsize, figtitle="ADMM reconstruction")
savefig("admm_result.png", dpi=150, bbox_inches="tight")
println("Saved admm_result.png")

# # Reference: OITOOLS reconstruct with TV + centering
# println("\n=== OITOOLS reconstruct (baseline) ===")
# x_ref = reconstruct(copy(x_start), data, ft_nfft;
#                      regularizers=[["centering", 1e5], ["tv", 3e2]],
#                      verb=true, maxiter=500)
# println("Reference:")
# image_to_chi2(x_ref, ft_nfft, data; verb=true)
# imdisp(x_ref, pixsize=pixsize, figtitle="OITOOLS reference")
# savefig("oitools_ref.png", dpi=150, bbox_inches="tight")
# println("Saved oitools_ref.png")
