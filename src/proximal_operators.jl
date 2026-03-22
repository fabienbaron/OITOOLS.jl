# ═══════════════════════════════════════════════════════════════════════════════
# Proximal operators for ADMM image reconstruction
# ═══════════════════════════════════════════════════════════════════════════════

using FFTW

# ── V2 proximal (PAINTER cubic roots) ────────────────────────────────────────

function realcuberoot(xx::Real)
    sign(xx) * abs(xx)^(1/3)
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

"""
    prox_v2(z0, v2_data, v2_err, indx_v2, rho)

Proximal operator for squared visibility using PAINTER cubic roots.
Modifies visibility amplitudes at `indx_v2` to balance V2 data fit
against proximity to input `z0`, with penalty weight `rho`.
"""
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

# ── T3phi proximal (Haniff least-squares via VMLMB on phases) ────────────────

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

"""
    prox_t3phi(z0, t3phi_deg, t3phi_err_deg, indx_t3, indx_t3_1, indx_t3_2, indx_t3_3, rho)

Proximal operator for closure phase using Haniff least-squares (VMLMB).
Optimizes phases of complex visibilities `z0` to fit closure phase data
while staying close to input, with penalty weight `rho`.
"""
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

# ── TV operators and proximal (Chambolle 2004) ───────────────────────────────

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
    d[1, :] .= p_h[1, :]
    @views d[2:nx-1, :] .= p_h[2:nx-1, :] .- p_h[1:nx-2, :]
    @views d[nx, :] .= .-p_h[nx-1, :]
    @views d[:, 1] .+= p_v[:, 1]
    @views d[:, 2:nx-1] .+= p_v[:, 2:nx-1] .- p_v[:, 1:nx-2]
    @views d[:, nx] .+= .-p_v[:, nx-1]
end

"""
    prox_tv(y, lambda; niter=100, tau=1/8)

Proximal operator for isotropic total variation (Chambolle 2004).
Solves `argmin_x (1/2)||x - y||² + λ TV(x)` via dual projected gradient ascent.
"""
function prox_tv(y::AbstractMatrix, lambda; niter=100, tau=1.0/8.0)
    nx = size(y, 1)
    p_h = zeros(nx, nx)
    p_v = zeros(nx, nx)
    d = zeros(nx, nx)
    g_h = zeros(nx, nx)
    g_v = zeros(nx, nx)
    nrm = zeros(nx, nx)

    for n in 1:niter
        div_op!(d, p_h, p_v)
        d .-= y ./ lambda
        grad_op!(g_h, g_v, d)
        p_h .+= tau .* g_h
        p_v .+= tau .* g_v
        @. nrm = max(1.0, sqrt(p_h^2 + p_v^2))
        p_h ./= nrm
        p_v ./= nrm
    end

    div_op!(d, p_h, p_v)
    return y .- lambda .* d
end

# ── L2 smoothness proximal via FFT ───────────────────────────────────────────

"""
    prox_l2smooth(y, lambda)

Proximal operator for L2 smoothness via FFT.
Solves `argmin_x (1/2)||x - y||² + (λ/2)||∇x||²` using the discrete Laplacian
eigendecomposition with periodic boundary conditions.
"""
function prox_l2smooth(y::AbstractMatrix, lambda)
    nx = size(y, 1)
    freq = [0:nx-1;]
    eig_h = 2 .- 2 .* cos.(2pi .* freq ./ nx)
    eig_2d = eig_h .+ eig_h'
    Y = fft(y)
    X = Y ./ (1 .+ lambda .* eig_2d)
    return real(ifft(X))
end

# ── Regularization norms ─────────────────────────────────────────────────────

"""
    tv_norm(x_img)

Isotropic total variation of a 2D image.
"""
function tv_norm(x_img)
    nx = size(x_img, 1)
    g_h = zeros(nx, nx)
    g_v = zeros(nx, nx)
    grad_op!(g_h, g_v, x_img)
    return sum(sqrt.(g_h.^2 .+ g_v.^2))
end

"""
    l2smooth_norm(x_img)

L2 smoothness norm (half squared gradient norm) of a 2D image.
"""
function l2smooth_norm(x_img)
    nx = size(x_img, 1)
    g_h = zeros(nx, nx)
    g_v = zeros(nx, nx)
    grad_op!(g_h, g_v, x_img)
    return 0.5 * sum(g_h.^2 .+ g_v.^2)
end

# ── Chi-squared at visibility vector ─────────────────────────────────────────

"""
    chi2_v2_at(z, v2_data, v2_err, indx_v2)

V2 chi-squared evaluated at a complex visibility vector.
"""
function chi2_v2_at(z, v2_data, v2_err, indx_v2)
    v2_model = abs2.(z[indx_v2])
    return sum(((v2_model .- v2_data) ./ v2_err).^2)
end

"""
    chi2_t3phi_at(z, t3phi_data, t3phi_err, indx_t3_1, indx_t3_2, indx_t3_3)

Closure phase chi-squared evaluated at a complex visibility vector.
"""
function chi2_t3phi_at(z, t3phi_data, t3phi_err, indx_t3_1, indx_t3_2, indx_t3_3)
    t3 = z[indx_t3_1] .* z[indx_t3_2] .* z[indx_t3_3]
    t3phi_model = angle.(t3) .* (180.0 / pi)
    res = mod.(mod.(t3phi_model .- t3phi_data .+ 180.0, 360.0) .+ 360.0, 360.0) .- 180.0
    return sum((res ./ t3phi_err).^2)
end
