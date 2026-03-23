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

# ── Group sparsity proximal (L2,1 on pixel values across channels) ────────────

"""
    prox_group_sparsity(y, lambda)

Proximal operator for group sparsity (L2,1 mixed norm) on a 3D cube `(nx, ny, nwav)`.
Block soft-thresholding: at each spatial pixel, the across-channel vector is
shrunk toward zero, promoting shared support across wavelength channels.

Solves `argmin_x (1/2)||x - y||² + λ Σ_{i,j} ||x[i,j,:]||₂`.
"""
function prox_group_sparsity(y::AbstractArray{<:AbstractFloat,3}, lambda)
    nx, ny, nwav = size(y)
    x = copy(y)
    for i in 1:nx, j in 1:ny
        nrm = 0.0
        for w in 1:nwav
            nrm += y[i,j,w]^2
        end
        nrm = sqrt(nrm)
        if nrm > 0
            scale = max(0.0, 1.0 - lambda / nrm)
            for w in 1:nwav
                x[i,j,w] = scale * y[i,j,w]
            end
        end
    end
    return x
end

"""
    group_sparsity_norm(x_cube)

L2,1 mixed norm of a 3D image cube: `Σ_{i,j} ||x[i,j,:]||₂`.
"""
function group_sparsity_norm(x_cube::AbstractArray{<:AbstractFloat,3})
    nx, ny, nwav = size(x_cube)
    f = 0.0
    for i in 1:nx, j in 1:ny
        s = 0.0
        for w in 1:nwav
            s += x_cube[i,j,w]^2
        end
        f += sqrt(s)
    end
    return f
end

# ── Group TV proximal (L2,1 on spatial gradients across channels) ─────────────

"""
    prox_grouptv(y, lambda; niter=100, tau=1/8)

Proximal operator for group (vectorial) total variation on a 3D cube `(nx, ny, nwav)`.
Encourages shared edge locations across all wavelength channels via Chambolle dual
projection with coupled L2,1 norm across the spectral dimension.

Solves `argmin_x (1/2)||x - y||² + λ GroupTV(x)`.
"""
function prox_grouptv(y::AbstractArray{<:AbstractFloat,3}, lambda;
                      niter=100, tau=1.0/8.0)
    nx, ny, nwav = size(y)
    p_h = zeros(nx, ny, nwav)
    p_v = zeros(nx, ny, nwav)
    d   = zeros(nx, ny, nwav)
    g_h = zeros(nx, ny, nwav)
    g_v = zeros(nx, ny, nwav)
    nrm = zeros(nx, ny)

    for n in 1:niter
        for w in 1:nwav
            div_op!(@view(d[:,:,w]), @view(p_h[:,:,w]), @view(p_v[:,:,w]))
        end
        d .-= y ./ lambda
        for w in 1:nwav
            grad_op!(@view(g_h[:,:,w]), @view(g_v[:,:,w]), @view(d[:,:,w]))
        end
        p_h .+= tau .* g_h
        p_v .+= tau .* g_v
        # Coupled projection: L2,1 norm across wavelengths
        nrm .= 0.0
        for w in 1:nwav
            @. nrm += @view(p_h[:,:,w])^2 + @view(p_v[:,:,w])^2
        end
        @. nrm = max(1.0, sqrt(nrm))
        for w in 1:nwav
            @view(p_h[:,:,w]) ./= nrm
            @view(p_v[:,:,w]) ./= nrm
        end
    end

    for w in 1:nwav
        div_op!(@view(d[:,:,w]), @view(p_h[:,:,w]), @view(p_v[:,:,w]))
    end
    return y .- lambda .* d
end

"""
    grouptv_norm(x_cube)

Group (vectorial) total variation of a 3D image cube:
`Σ_{i,j} √( Σ_λ [|∇_h x_{i,j,λ}|² + |∇_v x_{i,j,λ}|²] )`.
"""
function grouptv_norm(x_cube::AbstractArray{<:AbstractFloat,3})
    nx, ny, nwav = size(x_cube)
    g_h = zeros(nx, ny)
    g_v = zeros(nx, ny)
    S = zeros(nx, ny)
    for w in 1:nwav
        grad_op!(g_h, g_v, @view(x_cube[:,:,w]))
        S .+= g_h.^2 .+ g_v.^2
    end
    return sum(sqrt.(S))
end
