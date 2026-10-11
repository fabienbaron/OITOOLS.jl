# model_jacobian.jl
#
# The model's Jacobian ∂V/∂x in closed form (PLAN_ModelJacobian.md, J3):
#
#     V_i = K(ρ_i) Σ_c f_c S_c W_c(u'_i, v'_i; θ_c),     ∂V_i/∂x_k = Σ_q ∂V_i/∂p_q · ∂p_q/∂x_k
#
# Layer 1, ∂V/∂p, comes from the component kernels of model_kernels.jl: for each component a
# matrix G whose columns are ∂(f S W)/∂p for its flux, offsets, inclination, position angle and
# own parameters. Layer 2, ∂p/∂x, is the resolver run once on dual numbers: free parameters are
# a selection, fixed ones zero, derived ones the derivative of the user's expression, and a
# chromatic parameter carries one row per uv point. The two are contracted per component, and
# the spatial kernel K is applied last.
#
# Everything is generic in the element type T of the parameters and baselines, and allocates
# nothing for analytic components once a `ModelJacWorkspace` exists. Radial profiles and
# visibility functions still allocate inside their kernels (the transform products, the dual
# numbers through the user's expression).

struct _ModelJacTag end

# Typed by T alone. The dual numbers the resolver runs on carry the number of free parameters
# in their type; held untyped and used behind function barriers (`_seed!`, `_extract!`), they
# leave every other function here compiled once per precision, not once per count.
#
# Mutable so that it passes by reference through the per-component dynamic call: an immutable
# one this size would be copied to the heap at every component.
"""
    ModelJacWorkspace(model, x, uv)

Buffers for [`eval_model_jac!`](@ref): sized for `model`, for `length(x)` free parameters and
`size(uv, 2)` uv points, in the precision `eval_model` would use for `x` and `uv`. One per
thread; reusable for any parameter values with the same model and the same number of points.
"""
mutable struct ModelJacWorkspace{T<:AbstractFloat}
    n      :: Int                         # uv points
    nf     :: Int                         # free parameters
    xd     :: Any                         # x, the resolver's scalars and its per-point columns,
    sd     :: Any                         # as dual numbers (Vector{D}, Vector{D},
    cd     :: Any                         # Vector{Vector{D}})
    rp     :: ResolvedParams{T}           # the resolved values,
    Rs     :: Matrix{T}                   # ∂(scalar)/∂x, nscalar × nf,
    Rc     :: Vector{Matrix{T}}           # ∂(column)/∂x, one n × nf matrix per chromatic parameter
    ρ      :: Vector{T}                   # |uv|, cycles/rad
    u2     :: Vector{T}                   # the component's own (u', v', ρ') after incl/PA
    v2     :: Vector{T}
    ρ2     :: Vector{T}
    du_di  :: Vector{T}                   # ∂u'/∂incl, ∂u'/∂PA, ∂v'/∂PA (radians); ∂v'/∂incl = 0
    du_dpa :: Vector{T}
    dv_dpa :: Vector{T}
    Bh     :: Vector{Float64}             # ρ' in cycles/mas, or Mλ: as `eval_model` forms them
    ψ      :: Vector{T}                   # the uv position angle, for azimuthal modes
    Wr     :: Vector{T}                   # a real kernel's visibility, ∂/∂params, ∂/∂ρ
    Dr     :: Matrix{T}
    Dρ     :: Vector{T}
    W      :: Vector{Complex{T}}          # the component's visibility, complex
    G      :: Matrix{Complex{T}}          # ∂(f S W)/∂p: f, x, y, incl, pa, then its own
    Du     :: Vector{Complex{T}}          # ∂W/∂u', ∂W/∂v' (crescent); ∂W/∂B, ∂W/∂ψ (profiles)
    Dv     :: Vector{Complex{T}}
    pbuf   :: Vector{T}                   # a component's parameters: scalars,
    pcols  :: Vector{Vector{T}}           # or one value per point
    pvec   :: Vector{Vector{T}}
    qidx   :: Vector{Int}                 # resolver index of each column of G (0 = absent)
    F      :: Vector{T}                   # the total flux Σ f at each point, and ∂F/∂x (n × nf):
    JF     :: Matrix{T}                   # V and J are divided by it at the end
    fidx   :: Vector{Int}                 # each component's flux in the resolver (0 = f is 1)
end

const _G_FIXED = 5      # columns of G before a component's own parameters: f, x, y, incl, pa

_ncomp_params(c::AnalyticSpec) = length(c.param_idx)
_ncomp_params(c::VisFuncSpec)  = length(c.param_idx)
_ncomp_params(c::HankelSpec)   = 2 + length(c.profile_param_idx) + 2 * length(c.az_modes)
_ncomp_params(c::ImageSpec)    = 0

function ModelJacWorkspace(model::FlatModel, x::AbstractVector, uv::AbstractMatrix)
    T  = promote_type(float(eltype(x)), float(eltype(uv)))
    n  = size(uv, 2); nf = length(x)
    D  = ForwardDiff.Dual{_ModelJacTag, T, nf}
    r  = model.resolver
    np = maximum((_ncomp_params(c) for c in model.components); init = 0)
    vt() = Vector{T}(undef, n)
    vc() = Vector{Complex{T}}(undef, n)
    rp = ResolvedParams{T}(Vector{T}(undef, r.nscalar), [vt() for _ in 1:r.nchromatic], r.slot, n)
    return ModelJacWorkspace{T}(n, nf, Vector{D}(undef, nf), Vector{D}(undef, r.nscalar),
        [Vector{D}(undef, n) for _ in 1:r.nchromatic], rp, Matrix{T}(undef, r.nscalar, nf),
        [Matrix{T}(undef, n, nf) for _ in 1:r.nchromatic],
        vt(), vt(), vt(), vt(), vt(), vt(), vt(), Vector{Float64}(undef, n), vt(),
        vt(), Matrix{T}(undef, n, np), vt(),
        vc(), Matrix{Complex{T}}(undef, n, _G_FIXED + np), vc(), vc(),
        Vector{T}(undef, np), [vt() for _ in 1:np], Vector{Vector{T}}(undef, np),
        zeros(Int, _G_FIXED + np), vt(), Matrix{T}(undef, n, nf),
        Int[c.f_idx for c in model.components])
end

"""
    eval_model_jac!(V, J, ws, model, x, uv; wl = nothing, mjd = nothing, plane_wl = wl) -> (V, J)

The model's visibilities `V` (length N) and their Jacobian `J` (N × length(x), complex) with
respect to the free parameters, in place, from closed-form derivatives. `ws` is a
[`ModelJacWorkspace`](@ref) for this model, `x` and `uv`. Same conventions as
[`eval_model`](@ref); [`eval_model_grad_ad`](@ref) computes the same by automatic
differentiation.
"""
function eval_model_jac!(V::AbstractVector{Complex{T}}, J::AbstractMatrix{Complex{T}},
                         ws::ModelJacWorkspace{T}, model::FlatModel, x::AbstractVector,
                         uv::AbstractMatrix; wl = nothing, mjd = nothing, n::Int = 0,
                         plane_wl = wl, normalize::Bool = true) where {T}
    n == 0 || throw(ArgumentError("eval_model_jac! evaluates radial profiles at order 0 only"))
    N = size(uv, 2)
    (N == ws.n && length(V) == N && size(J) == (N, length(x)) && ws.nf == length(x)) ||
        throw(DimensionMismatch("eval_model_jac!: V, J or the workspace do not match x and uv"))
    _wl  = isnothing(wl)  ? NaN : wl
    _mjd = isnothing(mjd) ? NaN : mjd
    _resolve_jac!(ws, model, x, _wl, _mjd)
    _eval_model_jac!(V, J, ws, model, uv, _wl, _mjd, plane_wl, normalize)
    return V, J
end

"""
    eval_model_jac(model, x, uv; wl = nothing, mjd = nothing, plane_wl = wl) -> (V, J)

[`eval_model_jac!`](@ref) with its own buffers.
"""
function eval_model_jac(model::FlatModel, x::AbstractVector, uv::AbstractMatrix; kw...)
    ws = ModelJacWorkspace(model, x, uv)
    T  = eltype(ws.ρ)
    V  = Vector{Complex{T}}(undef, size(uv, 2))
    J  = Matrix{Complex{T}}(undef, size(uv, 2), length(x))
    return eval_model_jac!(V, J, ws, model, x, uv; kw...)
end

"""
    ModelJacCache()

Buffers for a caller that evaluates the Jacobian of one model many times: a
[`ModelJacWorkspace`](@ref) and `V`, `J`, rebuilt only when the model, the number of uv
points, the number of free parameters or the precision changes. Pass it to
[`eval_model_jac!`](@ref) or as `jac_cache` to [`model_to_chi2_fg`](@ref) and
`residuals_flat_jac`. The `V` and `J` returned belong to the cache and the next call
overwrites them. A smeared χ² keeps its sub-wavelength points here too. A cache is not shared
between tasks: make one per task.
"""
mutable struct ModelJacCache
    model :: Any
    ws    :: Any
    V     :: Any
    J     :: Any
    sub   :: Any
end
ModelJacCache() = ModelJacCache(nothing, nothing, nothing, nothing, nothing)

"""
    eval_model_jac!(cache::ModelJacCache, model, x, uv; wl, mjd, plane_wl) -> (V, J)

[`eval_model_jac!`](@ref) into the cache's buffers, which are rebuilt when they do not fit.
"""
function eval_model_jac!(c::ModelJacCache, model::FlatModel, x::AbstractVector,
                         uv::AbstractMatrix; kw...)
    T = promote_type(float(eltype(x)), float(eltype(uv)))
    N = size(uv, 2)
    ws = c.ws
    if !(c.model === model && ws isa ModelJacWorkspace{T} && ws.n == N && ws.nf == length(x))
        c.model = model
        c.ws = ws = ModelJacWorkspace(model, x, uv)
        c.V = Vector{Complex{T}}(undef, N)
        c.J = Matrix{Complex{T}}(undef, N, length(x))
    end
    return eval_model_jac!(c.V, c.J, ws, model, x, uv; kw...)
end

# `eval_model_grad`, through a cache when the caller holds one and the parameters are plain
# numbers (a Dual goes to the AD path, as `eval_model_grad` decides).
_model_jac(::Nothing, model, x, uv; kw...) = eval_model_grad(model, x, uv; kw...)
function _model_jac(c::ModelJacCache, model, x, uv; kw...)
    (eltype(x) <: Union{AbstractFloat, Integer} && eltype(uv) <: Union{AbstractFloat, Integer}) ||
        return eval_model_grad(model, x, uv; kw...)
    return eval_model_jac!(c, model, x, uv; kw...)
end

# Layer 2: the resolver on dual numbers, one partial per free parameter, then their values and
# partials copied into `rp`, `Rs` and `Rc`. The values are those `eval_model` resolves, bit for
# bit: dual arithmetic computes the same values.
function _resolve_jac!(ws::ModelJacWorkspace, model, x, wl, mjd)
    _seed!(ws.xd, x)
    # A new generated function for every model, called as `eval_model` calls it: dynamically,
    # on arguments that are all heap objects, so nothing is boxed and nothing recompiles.
    model.resolver.fill!(ws.sd, ws.cd, ws.xd, wl, mjd)
    _extract!(ws, ws.sd, ws.cd)
    return ws
end

function _seed!(xd::Vector{D}, x) where {T, NF, D<:ForwardDiff.Dual{_ModelJacTag, T, NF}}
    @inbounds for k in 1:NF
        xd[k] = D(T(x[k]), ForwardDiff.Partials(ntuple(j -> T(j == k), Val(NF))))
    end
    return nothing
end

function _extract!(ws::ModelJacWorkspace{T}, sd::Vector{D}, cd::Vector{Vector{D}}) where
                  {T, NF, D<:ForwardDiff.Dual{_ModelJacTag, T, NF}}
    @inbounds for j in eachindex(sd)
        ws.rp.s[j] = ForwardDiff.value(sd[j])
        for k in 1:NF
            ws.Rs[j, k] = ForwardDiff.partials(sd[j], k)
        end
    end
    @inbounds for c in eachindex(cd)
        col = cd[c]; vals = ws.rp.cols[c]; R = ws.Rc[c]
        for i in 1:ws.n
            vals[i] = ForwardDiff.value(col[i])
        end
        for k in 1:NF, i in 1:ws.n
            R[i, k] = ForwardDiff.partials(col[i], k)
        end
    end
    return nothing
end

# A resolved parameter at point i, and its derivative with respect to free parameter k.
@inline function _pv(rp::ResolvedParams, idx::Int, i::Int)
    s = rp.slot[idx]
    return s > 0 ? rp.s[s] : rp.cols[-s][i]
end

function _eval_model_jac!(V, J, ws::ModelJacWorkspace{T}, model, uv, _wl, _mjd, plane_wl,
                          normalize::Bool = true) where {T}
    N = ws.n
    @inbounds for i in 1:N
        u, v = uv[1, i], uv[2, i]
        ws.ρ[i] = sqrt(u^2 + v^2)        # in the baselines' precision, as `eval_model`
    end
    fill!(V, zero(Complex{T})); fill!(J, zero(Complex{T}))
    for (c, comp) in enumerate(model.components)
        _component_jac!(V, J, ws, comp, model.chromatic[c], uv, _wl, _mjd, plane_wl)
    end
    normalize && _normalise_jac!(V, J, ws, model)
    model.kernel_idx != 0 && _kernel_jac!(V, J, ws, model.kernel_idx)
    return V, J
end

# V/F and its Jacobian, F = Σ f: V ← V/F, then J ← (J − V ∂F/∂xᵀ)/F with that normalised V,
# as `eval_model` divides.
function _normalise_jac!(V, J, ws::ModelJacWorkspace{T}, model) where {T}
    F, JF, rp = ws.F, ws.JF, ws.rp
    fill!(F, zero(T)); fill!(JF, zero(T))
    # the components' flux indices, read from the workspace: `model.components` is abstractly
    # typed, and a field read through it would allocate
    @inbounds for fi in ws.fidx
        if fi == 0
            F .+= one(T)
            continue
        end
        s = rp.slot[fi]
        if s > 0
            v = rp.s[s]
            for i in 1:ws.n
                F[i] += v
            end
            for k in 1:ws.nf
                r = ws.Rs[s, k]
                iszero(r) && continue
                for i in 1:ws.n
                    JF[i, k] += r
                end
            end
        else
            col = rp.cols[-s]; R = ws.Rc[-s]
            for i in 1:ws.n
                F[i] += col[i]
            end
            for k in 1:ws.nf, i in 1:ws.n
                JF[i, k] += R[i, k]
            end
        end
    end
    @inbounds for i in 1:ws.n
        V[i] /= F[i]
    end
    @inbounds for k in 1:ws.nf, i in 1:ws.n
        J[i, k] = (J[i, k] - V[i] * JF[i, k]) / F[i]
    end
    return V, J
end

# One component: geometry, its kernel into G, then V += f S W and J += G ∂p/∂x.
function _component_jac!(V, J, ws::ModelJacWorkspace{T}, comp, chrom::Bool, uv, _wl, _mjd,
                         plane_wl) where {T}
    np = _ncomp_params(comp)
    q  = ws.qidx
    q[1] = comp.f_idx; q[2] = comp.x_idx; q[3] = comp.y_idx
    q[4] = comp.incl_idx; q[5] = comp.pa_idx
    _own_param_idx!(q, comp)
    oriented = comp.incl_idx != 0 || comp.pa_idx != 0
    _geometry!(ws, comp, uv, oriented)
    _kernel_into_G!(ws, comp, chrom, np, oriented, uv, _wl, _mjd, plane_wl)
    _flux_shift!(V, ws, comp, uv, np)
    _contract!(J, ws, _G_FIXED + np)
    return nothing
end

function _own_param_idx!(q, comp::Union{AnalyticSpec, VisFuncSpec})
    for (j, k) in enumerate(comp.param_idx)
        q[_G_FIXED + j] = k
    end
end
# A profile's columns: d_min, d_max (0 when fixed), its parameters, then the azimuthal modes.
function _own_param_idx!(q, comp::HankelSpec)
    np = 2 + length(comp.profile_param_idx)
    q[_G_FIXED + 1] = comp.dmin_idx
    q[_G_FIXED + 2] = comp.dmax_idx
    for (j, k) in enumerate(comp.profile_param_idx)
        q[_G_FIXED + 2 + j] = k
    end
    for (m, az) in enumerate(comp.az_modes)
        q[_G_FIXED + np + 2m - 1] = az.amp_idx
        q[_G_FIXED + np + 2m]     = az.phi_idx
    end
end
_own_param_idx!(q, ::ImageSpec) = nothing

# (u', v', ρ') and their derivatives in incl and PA (radians), PMOIRED's convention as in
# `_stretch`: u' = cos i (cos PA u − sin PA v), v' = sin PA u + cos PA v.
function _geometry!(ws::ModelJacWorkspace{T}, comp, uv, oriented::Bool) where {T}
    rp = ws.rp
    @inbounds for i in 1:ws.n
        u, v = uv[1, i], uv[2, i]
        if oriented
            ir = (comp.incl_idx != 0 ? _pv(rp, comp.incl_idx, i) : zero(T)) * π / 180
            pr = (comp.pa_idx   != 0 ? _pv(rp, comp.pa_idx,   i) : zero(T)) * π / 180
            si, ci = sin(ir), cos(ir); sp, cp = sin(pr), cos(pr)
            a  = cp * u - sp * v
            u2 = ci * a
            v2 = sp * u + cp * v
            ws.u2[i] = u2; ws.v2[i] = v2; ws.ρ2[i] = sqrt(u2^2 + v2^2)
            ws.du_di[i]  = -si * a
            ws.du_dpa[i] = ci * (-sp * u - cp * v)
            ws.dv_dpa[i] = a
        else
            ws.u2[i] = u; ws.v2[i] = v; ws.ρ2[i] = ws.ρ[i]
        end
    end
    return ws
end

# ∂ρ'/∂incl, ∂ρ'/∂PA (per radian); zero at ρ' = 0, where W is stationary.
@inline function _dρ(ws, i)
    r = ws.ρ2[i]
    r == 0 && return (zero(r), zero(r))
    return (ws.u2[i] * ws.du_di[i] / r, (ws.u2[i] * ws.du_dpa[i] + ws.v2[i] * ws.dv_dpa[i]) / r)
end

# The component's parameters as the kernels take them, in the first `length(idx)` entries: a
# vector of scalars, or (chromatic) of per-point columns, a scalar parameter spread into a
# buffer. The kernels index their parameters by position, so the buffers pass whole.
function _kernel_params(ws::ModelJacWorkspace{T}, idx, chrom::Bool) where {T}
    rp = ws.rp
    if !chrom
        p = ws.pbuf
        @inbounds for (j, k) in enumerate(idx)
            p[j] = rp.s[rp.slot[k]]
        end
        return p
    end
    p = ws.pvec
    @inbounds for (j, k) in enumerate(idx)
        s = rp.slot[k]
        if s > 0
            fill!(ws.pcols[j], rp.s[s]); p[j] = ws.pcols[j]
        else
            p[j] = rp.cols[-s]
        end
    end
    return p
end

const _DEG = π / 180

# ── analytic, radial: W(ρ'), ∂W/∂params, ∂W/∂ρ' ──────────────────────────────────────────
function _kernel_into_G!(ws::ModelJacWorkspace{T}, comp::AnalyticSpec{K}, chrom, np, oriented,
                         uv, _wl, _mjd, plane_wl) where {T, K}
    p = _kernel_params(ws, comp.param_idx, chrom)
    _vis_jac!(ws.Wr, view(ws.Dr, :, 1:np), ws.Dρ, Val(K), p, ws.ρ2)
    G = ws.G
    @inbounds for i in 1:ws.n
        ws.W[i] = ws.Wr[i]
        for j in 1:np
            G[i, _G_FIXED + j] = ws.Dr[i, j]
        end
        if oriented
            dri, drp = _dρ(ws, i)
            G[i, 4] = ws.Dρ[i] * dri * _DEG
            G[i, 5] = ws.Dρ[i] * drp * _DEG
        end
    end
    return ws
end

# ── crescent: W(u', v') ──────────────────────────────────────────────────────────────────
function _kernel_into_G!(ws::ModelJacWorkspace{T}, comp::AnalyticSpec{:crescent}, chrom, np,
                         oriented, uv, _wl, _mjd, plane_wl) where {T}
    p = _kernel_params(ws, comp.param_idx, chrom)
    G = ws.G
    _crescent_jac!(ws.W, view(G, :, _G_FIXED+1:_G_FIXED+np), ws.Du, ws.Dv, p, ws.u2, ws.v2)
    if oriented
        @inbounds for i in 1:ws.n
            G[i, 4] = ws.Du[i] * ws.du_di[i] * _DEG
            G[i, 5] = (ws.Du[i] * ws.du_dpa[i] + ws.Dv[i] * ws.dv_dpa[i]) * _DEG
        end
    end
    return ws
end

# ── a written visibility function of B (Mλ) ──────────────────────────────────────────────
function _kernel_into_G!(ws::ModelJacWorkspace{T}, comp::VisFuncSpec, chrom, np, oriented,
                         uv, _wl, _mjd, plane_wl) where {T}
    p = view(_kernel_params(ws, comp.param_idx, chrom), 1:np)    # the expression takes np
    @inbounds for i in 1:ws.n
        ws.Bh[i] = ws.ρ2[i] / 1e6
    end
    G = ws.G
    _visfunc_jac!(ws.W, view(G, :, _G_FIXED+1:_G_FIXED+np), ws.Du, comp.vis_fn, p, ws.Bh,
                  _wl, _mjd)
    if oriented
        @inbounds for i in 1:ws.n
            dri, drp = _dρ(ws, i)
            G[i, 4] = ws.Du[i] * dri / 1e6 * _DEG
            G[i, 5] = ws.Du[i] * drp / 1e6 * _DEG
        end
    end
    return ws
end

# ── a radial profile, with azimuthal modes ──────────────────────────────────────────────
function _kernel_into_G!(ws::ModelJacWorkspace{T}, comp::HankelSpec, chrom, np, oriented,
                         uv, _wl, _mjd, plane_wl) where {T}
    rp = ws.rp
    # [d_min, d_max, profile params...], scalars always
    dmin = comp.dmin_idx != 0 ? T(rp.s[rp.slot[comp.dmin_idx]]) : T(comp.d_min)
    dmax = comp.dmax_idx != 0 ? T(rp.s[rp.slot[comp.dmax_idx]]) : T(comp.d_max)
    prof = vcat(dmin, dmax, T[rp.s[rp.slot[k]] for k in comp.profile_param_idx])
    npp = length(prof)
    @inbounds for i in 1:ws.n
        ws.Bh[i] = ws.ρ2[i] / _MAS2RAD_PM
        ws.ψ[i]  = atan(ws.v2[i], ws.u2[i])
    end
    # A size no parameter moves has zero columns, and is left out of what follows: neither
    # ∂I/∂d, nor the size's motion of the jumps, nor the kernels' part is computed for it.
    sized = _depends(ws, comp.dmin_idx) || _depends(ws, comp.dmax_idx)
    c0 = sized ? 1 : 3
    st = _hankel_setup(comp, ws.Bh, prof; cached = !oriented,      # the nodes follow the data
                       jumpcols = c0:npp)
    s, r, μ, g = st.s, st.r, st.μ, st.g
    Ti = promote_type(T, eltype(r))
    Ifun(pp) = (v = comp.profile_fn(s, μ, pp...); v isa Number ? fill(v, length(s)) : v)
    I  = convert(Vector{Ti}, Ifun(prof))
    nm = length(comp.az_modes)
    if sized
        dI = Matrix{Ti}(ForwardDiff.jacobian(Ifun, prof))
    else
        ws.G[:, _G_FIXED+1:_G_FIXED+2] .= zero(eltype(ws.G))
        dI = npp == 2 ? zeros(Ti, length(s), 0) :
             Matrix{Ti}(ForwardDiff.jacobian(p -> (v = comp.profile_fn(s, μ, dmin, dmax, p...);
                                                   v isa Number ? fill(v, length(s)) : v),
                                             prof[3:end]))
    end
    jumps = st.jumps
    modes = Int[az.order for az in comp.az_modes]
    cache = oriented || !st.cacheable ? nothing : comp.kcache
    D = view(ws.G, :, _G_FIXED+c0:_G_FIXED+npp+2nm)
    if chrom
        amps = Vector{T}[_pcol(ws, az.amp_idx) for az in comp.az_modes]
        phis = Vector{T}[_pcol(ws, az.phi_idx) for az in comp.az_modes]
        _hankel_jac!(ws.W, D, ws.Du, ws.Dv, r, g, I, dI, ws.Bh, ws.ψ, modes, amps, phis, cache,
                     jumps)
    else
        amps = T[_pv(rp, az.amp_idx, 1) for az in comp.az_modes]
        phis = T[_pv(rp, az.phi_idx, 1) for az in comp.az_modes]
        _hankel_jac!(ws.W, D, ws.Du, ws.Dv, r, g, I, dI, ws.Bh, ws.ψ, modes, amps, phis, cache,
                     jumps)
    end
    # r = (d_max/2)·s: the kernels' part of ∂/∂d_max is (B/d_max)·∂/∂B, every order included
    if sized
        @inbounds for i in 1:ws.n
            D[i, 2] += ws.Bh[i] / dmax * ws.Du[i]
        end
    end
    G = ws.G
    if oriented
        @inbounds for i in 1:ws.n
            dri, drp = _dρ(ws, i)
            r2 = ws.ρ2[i]^2
            dψi, dψp = r2 == 0 ? (zero(T), zero(T)) :
                (-ws.v2[i] * ws.du_di[i] / r2,
                 (ws.u2[i] * ws.dv_dpa[i] - ws.v2[i] * ws.du_dpa[i]) / r2)
            G[i, 4] = (ws.Du[i] * dri / _MAS2RAD_PM + ws.Dv[i] * dψi) * _DEG
            G[i, 5] = (ws.Du[i] * drp / _MAS2RAD_PM + ws.Dv[i] * dψp) * _DEG
        end
    end
    return ws
end

_pcol(ws, idx) = (s = ws.rp.slot[idx]; s > 0 ? fill(ws.rp.s[s], ws.n) : ws.rp.cols[-s])

# Whether the resolved parameter `idx` moves with the fit vector (0: not a parameter).
function _depends(ws::ModelJacWorkspace, idx::Int)
    idx == 0 && return false
    s = ws.rp.slot[idx]
    return s > 0 ? any(!iszero, view(ws.Rs, s, :)) : true
end

# ── an image: cached visibilities, no parameters of its own ────────────────────────────────
function _kernel_into_G!(ws::ModelJacWorkspace{T}, comp::ImageSpec, chrom, np, oriented, uv,
                         _wl, _mjd, plane_wl) where {T}
    Vi = model_image_vis(comp.image, uv; wl = plane_wl, mjd = _mjd, cache = comp.cache,
                         name = comp.comp_name)
    ws.W .= Vi
    return ws
end

# Flux and offset: V += f S W; G's columns become ∂(f S W)/∂p.
function _flux_shift!(V, ws::ModelJacWorkspace{T}, comp, uv, np) where {T}
    rp = ws.rp; G = ws.G
    shifted  = comp.x_idx != 0 || comp.y_idx != 0
    oriented = comp.incl_idx != 0 || comp.pa_idx != 0
    @inbounds for i in 1:ws.n
        f = comp.f_idx != 0 ? _pv(rp, comp.f_idx, i) : one(T)
        S = one(Complex{T})
        if shifted
            dx = comp.x_idx != 0 ? _pv(rp, comp.x_idx, i) : zero(T)
            dy = comp.y_idx != 0 ? _pv(rp, comp.y_idx, i) : zero(T)
            S  = cis(-2π * (uv[1, i] * dx + uv[2, i] * dy) / _MAS2RAD_PM)
        end
        fW  = f * ws.W[i]                   # (f W) S, `eval_model`'s order
        fSW = shifted ? fW * S : fW
        V[i] += fSW
        G[i, 1] = S * ws.W[i]
        if shifted
            G[i, 2] = fSW * (-2π * im * uv[1, i] / _MAS2RAD_PM)
            G[i, 3] = fSW * (-2π * im * uv[2, i] / _MAS2RAD_PM)
        end
        fS = f * S
        if oriented
            G[i, 4] *= fS; G[i, 5] *= fS
        end
        for j in 1:np
            G[i, _G_FIXED + j] *= fS
        end
    end
    return V
end

# J[:, k] += Σ_q G[:, q] ∂p_q/∂x_k, over the columns whose parameter exists.
function _contract!(J, ws::ModelJacWorkspace{T}, nq) where {T}
    G = ws.G; slot = ws.rp.slot
    @inbounds for qq in 1:nq
        idx = ws.qidx[qq]
        idx == 0 && continue
        s = slot[idx]
        if s > 0
            for k in axes(J, 2)
                r = ws.Rs[s, k]
                iszero(r) && continue
                for i in 1:ws.n
                    J[i, k] += G[i, qq] * r
                end
            end
        else
            R = ws.Rc[-s]
            for k in axes(J, 2), i in 1:ws.n
                J[i, k] += G[i, qq] * R[i, k]
            end
        end
    end
    return J
end

# The spatial kernel, K = exp(−t²/ln 2), t = π (FWHM/2) ρ / MAS2RAD, on the raw baselines:
# J ← K J + V ∂K/∂FWHM ∂FWHM/∂x, then V ← K V.
function _kernel_jac!(V, J, ws::ModelJacWorkspace{T}, kidx) where {T}
    s = ws.rp.slot[kidx]
    @inbounds for i in 1:ws.n
        fw = _pv(ws.rp, kidx, i)
        t  = π * (fw / 2) * ws.ρ[i] / _MAS2RAD_PM
        K  = exp(-t^2 / log(2))
        dK = K * (-2t / log(2)) * π * ws.ρ[i] / (2 * _MAS2RAD_PM)
        for k in axes(J, 2)
            P = s > 0 ? ws.Rs[s, k] : ws.Rc[-s][i, k]
            J[i, k] = K * J[i, k] + V[i] * dK * P
        end
        V[i] *= K
    end
    return V, J
end
