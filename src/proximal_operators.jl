# ═══════════════════════════════════════════════════════════════════════════════
# Proximal operators for BSDMM image reconstruction
# ═══════════════════════════════════════════════════════════════════════════════

import FFTW: fft, ifft

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

# ── Centering proximal (Woodbury) ─────────────────────────────────────────────

"""
    prox_centering!(x, y, lambda, p, q, AtA)

Proximal operator for centering: `argmin_x (1/2)||x-y||² + λ[(p'x)² + (q'x)²]`.
Uses Woodbury formula with `p` = centered column indices, `q` = centered row indices,
`AtA` = pre-computed `[p'p p'q; q'p q'q]` (2×2 matrix).
Result stored in `x`.
"""
function prox_centering!(x::AbstractVector, y::AbstractVector, lambda,
                         p::AbstractVector, q::AbstractVector, AtA::AbstractMatrix)
    c1 = dot(p, y)
    c2 = dot(q, y)
    M = [1 + 2lambda * AtA[1,1]  2lambda * AtA[1,2];
         2lambda * AtA[2,1]      1 + 2lambda * AtA[2,2]]
    s = M \ [c1, c2]
    @. x = y - 2lambda * (s[1] * p + s[2] * q)
    return x
end
