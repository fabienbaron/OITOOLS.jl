# The companion search: a grid over a companion's position, the companion fitted at every node,
# the map of χ² that results, and how significant its best node is against no companion at all.
#
# The primary is whatever the model holds -- an image component, a uniform disc, anything --
# and the search only moves the companion. Two paths:
#
#   * fast, when the companion's flux is the only free parameter (CANDID's `chi2Map`): the
#     model is linear in that flux, so the primary's visibilities and the companion's at the
#     origin are computed once, and at each node only the companion's phase changes. V² is then
#     a quadratic and a triple product a cubic in the flux, sub-wavelength by sub-wavelength, and
#     a node is a 1-D fit of a polynomial-built χ².
#   * general, for anything else free at a node (a companion diameter, the primary's own
#     parameters): VMLMB with bounds, the companion's position held at the node.
#
# Threads run over the nodes; nothing in the loop takes a lock.

using OptimPackNextGen: BraDi

"""
    CompanionSearch

What [`companion_search`](@ref) returns.

| field | meaning |
|---|---|
| `x`, `y` | the grid axes, mas (`x` East, `y` North) |
| `chi2r` | `(length(x), length(y))`, reduced χ² at each node; `NaN` where not searched |
| `params` | `(length(free), length(x), length(y))`, the free parameters fitted at each node |
| `free` | their names, `companion,f` among them |
| `companion` | the companion component's name |
| `best` | the best node: `(x, y, sep, pa, f, dmag, chi2r, params)`; `pa` in degrees East of North, `dmag` = −2.5 log₁₀(f/(1−f)) |
| `chi2r_null` | reduced χ² of the model without the companion (`f = 0`) |
| `nsigma` | the best node's significance against it (CANDID's convention) |
| `ndof` | the degrees of freedom the significance uses: data points − 1 |
| `fit` | a [`FitResult`](@ref) for the best node, the companion's position fixed |
| `path` | `:fast` or `:general` |
| `solver` | the 1-D solver of the fast path (`:bradi` or `:halley`), `:vmlmb` for the general path |
| `fallbacks` | nodes where Halley did not converge and BraDi was used instead |
| `nnodes`, `seconds` | how many nodes were fitted, and how long it took |
"""
struct CompanionSearch
    x          :: Vector{Float64}
    y          :: Vector{Float64}
    chi2r      :: Matrix{Float64}
    params     :: Array{Float64,3}
    free       :: Vector{String}
    companion  :: String
    best       :: NamedTuple
    chi2r_null :: Float64
    nsigma     :: Float64
    ndof       :: Int
    fit        :: FitResult
    path       :: Symbol
    solver     :: Symbol
    fallbacks  :: Int
    nnodes     :: Int
    seconds    :: Float64
end

"""The fitted companion flux at every node, `(length(x), length(y))`."""
companion_flux_map(s::CompanionSearch) =
    s.params[findfirst(==(s.companion * ",f"), s.free), :, :]

function Base.show(io::IO, s::CompanionSearch)
    b = s.best
    @printf(io, "CompanionSearch (%s path, %s): %d nodes in %.1f s\n", s.path, s.solver,
            s.nnodes, s.seconds)
    @printf(io, "  best: x = %+.3f, y = %+.3f mas (sep %.3f mas, PA %.2f°)\n",
            b.x, b.y, b.sep, b.pa)
    @printf(io, "        %s,f = %.5f (Δm = %.3f), χ²ᵣ = %.4f\n", s.companion, b.f, b.dmag,
            b.chi2r)
    @printf(io, "  without companion: χ²ᵣ = %.4f;  significance %.1fσ (ndof %d)\n",
            s.chi2r_null, s.nsigma, s.ndof)
end

"""
    detection_sigma(chi2r_null, chi2r_best, ndof) -> Float64

How many σ the better model is preferred by, as CANDID (Gallenne et al. 2015) computes it: the
best model is taken to describe the data (its χ²ᵣ rescaled to 1), the χ² of the null model is
then `ndof · χ²ᵣ,null / χ²ᵣ,best`, and its probability under a χ² distribution with `ndof`
degrees of freedom is converted to the equivalent number of Gaussian σ. Not exact — it trusts
the error bars up to a common factor — which is why it is a convention to state.
"""
function detection_sigma(chi2r_null::Real, chi2r_best::Real, ndof::Integer)
    x = ndof * chi2r_null / chi2r_best
    p = gamma_inc(ndof / 2, x / 2)[2]           # survival function of χ²(ndof)
    p > 0 && return sqrt(2) * erfcinv(p)
    # Beyond Float64: the tail's leading term, log p ≈ (k/2−1)·log(x/2) − x/2 − logΓ(k/2),
    # and n_σ ≈ √(−2 log p) for a Gaussian tail.
    lp = (ndof / 2 - 1) * log(x / 2) - x / 2 - loggamma(ndof / 2)
    return sqrt(-2lp)
end

"""
    companion_search(model_dict, data; companion = "comp", free = [companion * ",f"],
                     lb = Dict(), ub = Dict(), fmax = 0.5,
                     step = nothing, rmin = nothing, rmax = nothing, mask = 0.01,
                     weights = [1, 1, 1], solver = :bradi, threads = true,
                     on_grid = nothing) -> CompanionSearch

Search for a companion around the model's primary: on a grid of positions, fit the companion
at each node with its position held there, map χ², and report the best node and its
significance against the model without a companion.

**The model** is the primary plus a component `companion` with keys `companion,f`,
`companion,x` and `companion,y` (their values are only starting points; x and y are set by the
grid). The primary's flux should be written as `1 - \$companion,f`, so that the total is 1 at
every node. A point companion is a component with no shape key; give it `ud` (and free it) for
a resolved one. Bandwidth smearing follows the model's `"smear"` setting.

**Free at each node**: `free`, which must include `companion,f`, bounded to `[0, fmax]`.
- Only `companion,f` free → the fast path, everything else fixed at the dict's values. Fit the
  primary first (without the companion) and pass its fitted values.
- Anything else → the general path, VMLMB with `lb`/`ub` (defaults from `default_bounds`).

**The grid** (defaults from the data):
- `step` λ/(2 B_max) — half of PMOIRED's λ/B_max, since the position is not refined at a node;
- `rmax` R λ/(2 B_max), the bandwidth-smearing radius, when the data have channel widths;
- `rmin` the extent of the primary's analytic components (offset plus radius), at least `step`;
- `mask`: nodes where the primary's image components, convolved with a beam of FWHM
  λ/B_max, exceed `mask` times their peak, are not searched. `0` turns it off.

`solver` (fast path): `:bradi`, a global 1-D minimisation on [0, fmax] (OptimPackNextGen's
BraDi), or `:halley`, Halley's method from the best of a coarse sampling, falling back to BraDi
at any node where it does not converge in the bounds.

`weights` are those of `model_to_chi2`; the fast path takes V², T3amp, T3φ, VISAMP and VISPHI.

`path = :general` uses VMLMB even when only the flux is free, which is how the two paths are
compared.

`on_grid`, when given, is called once the grid is laid out and before any node is fitted, with
`(; x, y, chi2r, params, nnodes)`: the axes, and the arrays the search then fills node by node
(`NaN` until a node is done). It is how a window draws the map while the search runs; the
arrays are written by the search's threads while the caller reads them.
"""
function companion_search(model_dict::AbstractDict, data::OIdata;
                          companion::AbstractString = "comp",
                          free::AbstractVector{<:AbstractString} = [companion * ",f"],
                          lb = Dict{String,Float64}(), ub = Dict{String,Float64}(),
                          fmax::Real = 0.5,
                          step = nothing, rmin = nothing, rmax = nothing, mask::Real = 0.01,
                          weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0],
                          solver::Symbol = :bradi, path::Symbol = :auto,
                          threads::Bool = true, on_grid = nothing)
    t0 = time()
    md = Dict{String,Any}(String(k) => v for (k, v) in model_dict)
    cn = String(companion)
    kf, kx, ky = cn * ",f", cn * ",x", cn * ",y"
    for k in (kf, kx, ky)
        haskey(md, k) || throw(ArgumentError(
            "the model has no $k: the companion is a component `$cn` with f, x and y"))
    end
    free = String.(free)
    kf in free || throw(ArgumentError("$kf must be free: it is what each node fits"))
    (kx in free || ky in free) && throw(ArgumentError(
        "$kx and $ky are set by the grid, not fitted at a node"))
    solver in (:bradi, :halley) || throw(ArgumentError("solver is :bradi or :halley"))
    path in (:auto, :general) || throw(ArgumentError("path is :auto or :general"))
    w = _pad_weights(weights)
    any(>(0), w[6:7]) && throw(ArgumentError(
        "the search fits V², T3 and VISAMP/VISPHI; OI_FLUX and differential phase are not taken"))
    md[kx] = 0.0; md[ky] = 0.0

    # ── the grid ──────────────────────────────────────────────────────────────────────
    g = companion_grid(md, data; companion = cn, step, rmin, rmax, mask)
    xs, ys, nodes = g.x, g.y, g.nodes

    nfree = length(free)
    chi2r  = fill(NaN, length(xs), length(ys))
    params = fill(NaN, nfree, length(xs), length(ys))
    ndata  = _ndof(data, w)
    nfit   = nfree + 2                                 # the free parameters and the position
    fallbacks = Threads.Atomic{Int}(0)
    on_grid === nothing || on_grid((; x = xs, y = ys, chi2r, params, nnodes = length(nodes)))

    if free == [kf] && path === :auto
        path = :fast
        fp = _fast_problem(md, data, cn, w)
        work = () -> FastWork(fp)
        run = function (ws, i, j)
            f, c2, fb = _fast_node(ws, fp, xs[i], ys[j], fmax, solver)
            fb && Threads.atomic_add!(fallbacks, 1)
            chi2r[i, j] = c2 / (ndata - nfit)
            params[1, i, j] = f
        end
        # The primary alone: any node at flux 0.
        chi2_null = _fast_chi2(0.0, _fast_coefficients!(FastWork(fp), fp, 0.0, 0.0), fp)
        nnull = nfree - 1
    else
        path = :general
        gp = _general_problem(md, data, cn, free, lb, ub, fmax, w)
        work = ModelJacCache
        run = function (cache, i, j)
            p, c2 = _general_node(gp, xs[i], ys[j]; cache)
            chi2r[i, j] = c2 / (ndata - nfit)
            params[:, i, j] .= p
        end
        chi2_null = _general_null(gp)
        nnull = nfree - 1
    end

    # One task per share of the nodes, each with its own scratch, so nothing is shared but the
    # output arrays, each node writing only its own entries. The nodes are visited in a
    # golden-ratio stride and the shares interleave, so at any moment the fitted nodes are
    # spread over the whole grid: a window drawing the map sees it fill evenly instead of row by
    # row, and the nodes that cost more are spread across the tasks. No node's result depends
    # on the order.
    nt = threads ? Threads.nthreads() : 1
    order = _spread_order(length(nodes))
    @sync for chunk in (view(order, k:nt:length(order)) for k in 1:min(nt, length(order)))
        Threads.@spawn begin
            ws = work()
            for q in chunk
                i, j = nodes[q]
                run(ws, i, j)
            end
        end
    end

    # ── the best node ─────────────────────────────────────────────────────────────────
    ib = argmin(c -> isnan(chi2r[c]) ? Inf : chi2r[c], CartesianIndices(chi2r))
    bx, by = xs[ib[1]], ys[ib[2]]
    pb = params[:, ib[1], ib[2]]
    fb = pb[findfirst(==(kf), free)]
    best = (x = bx, y = by, sep = hypot(bx, by), pa = mod(atand(bx, by), 360), f = fb,
            dmag = -2.5 * log10(fb / (1 - fb)), chi2r = chi2r[ib], params = pb)
    chi2r_null = chi2_null / (ndata - nnull)
    ndof = ndata - 1

    bd = merge(md, Dict{String,Any}(kx => bx, ky => by), Dict{String,Any}(zip(free, pb)))
    fm = dict_to_model(bd, free)
    c2b = model_to_chi2(fm, pb, data; weights = w)
    fit = FitResult(pb, free, c2b, c2b / (ndata - nfit), ndata, length(nodes), :GRID, fm)

    return CompanionSearch(xs, ys, chi2r, params, free, cn, best, chi2r_null,
                           detection_sigma(chi2r_null, best.chi2r, ndof), ndof, fit, path,
                           path === :fast ? solver : :vmlmb, fallbacks[], length(nodes),
                           time() - t0)
end

"""
    companion_grid(model_dict, data; companion = "comp", step = nothing, rmin = nothing,
                   rmax = nothing, mask = 0.01) -> NamedTuple

The grid [`companion_search`](@ref) searches, with the same defaults, without searching it:
`(; x, y, nodes, step, rmin, rmax)`, where `nodes` are the `(i, j)` indices into `x` and `y`
of the nodes that will be fitted (inside the radius range, outside the mask). For a panel that
shows how large a search is before it runs.
"""
function companion_grid(model_dict::AbstractDict, data::OIdata; companion::AbstractString = "comp",
                        step = nothing, rmin = nothing, rmax = nothing, mask::Real = 0.01)
    md = Dict{String,Any}(String(k) => v for (k, v) in model_dict)
    cn = String(companion)
    md[cn * ",x"] = 0.0; md[cn * ",y"] = 0.0
    ρ = hypot.(data.uv[1, :], data.uv[2, :])
    bmax = maximum(ρ)                                  # cycles/rad: B_max/λ at its shortest λ
    st = step === nothing ? _MAS2RAD_PM / (2bmax) : Float64(step)
    R = median(data.uv_lam ./ data.uv_dlam)
    rmx = if rmax !== nothing
        Float64(rmax)
    elseif isfinite(R) && R > 0
        R * _MAS2RAD_PM / (2bmax)
    else
        throw(ArgumentError("the data have no channel widths to set rmax from; give rmax"))
    end
    m_probe = dict_to_model(md, String[])
    prim_extent = _model_extent(Dict(k => v for (k, v) in md
                                     if !startswith(k, cn * ",") &&
                                        !haskey(md, first(split(k, ',')) * ",image")),
                                [c for c in _component_names(md) if c != cn])
    rmn = rmin === nothing ? max(st, prim_extent) : Float64(rmin)
    # Symmetric about the origin, so a node's mirror image is a node too.
    nh = floor(Int, rmx / st)
    xs = collect((-nh:nh) .* st); ys = copy(xs)
    masked = _image_mask(m_probe, md, cn, xs, ys, mask, bmax)
    nodes = [(i, j) for i in eachindex(xs), j in eachindex(ys)
             if rmn <= hypot(xs[i], ys[j]) <= rmx && !masked[i, j]]
    isempty(nodes) && throw(ArgumentError(
        "no node between rmin = $(round(rmn; digits = 3)) and rmax = $(round(rmx; digits = 3)) mas " *
        "outside the mask"))
    return (; x = xs, y = ys, nodes, step = st, rmin = rmn, rmax = rmx)
end

# 1:n visited as k·p mod n with p coprime to n near n/φ: consecutive entries land far apart,
# and any prefix is spread evenly over 1:n.
function _spread_order(n::Integer)
    n <= 2 && return collect(1:n)
    p = max(1, round(Int, n * 0.6180339887498949))
    while gcd(p, n) != 1
        p += 1
    end
    return [mod(k * p, n) + 1 for k in 0:n-1]
end

# ─────────────────────────────────────────────────────────────────────────────────────
# the mask: nodes on a bright part of an image primary
# ─────────────────────────────────────────────────────────────────────────────────────

function _image_mask(m::FlatModel, md, cn, xs, ys, thr, bmax)
    masked = falses(length(xs), length(ys))
    thr > 0 || return masked
    for comp in m.components
        comp isa ImageSpec && comp.comp_name != cn || continue
        mi = comp.image
        img = dropdims(sum(mi.image; dims = (3, 4)); dims = (3, 4))
        nx, ny = size(img)
        # A Gaussian beam of FWHM λ/B_max, in pixels, applied by FFT.
        σpix = _MAS2RAD_PM / bmax / mi.pixsize / (2 * sqrt(2 * log(2)))
        fx = fftfreq(nx); fy = fftfreq(ny)
        G = [exp(-2 * (π * σpix)^2 * (a^2 + b^2)) for a in fx, b in fy]
        conv = real.(FFTW.ifft(FFTW.fft(img) .* G))
        peak = maximum(conv)
        ox = mi.offset[1] + (haskey(md, comp.comp_name * ",x") ? _num(md, comp.comp_name * ",x") : 0.0)
        oy = mi.offset[2] + (haskey(md, comp.comp_name * ",y") ? _num(md, comp.comp_name * ",y") : 0.0)
        cx, cy = nx ÷ 2 + 1, ny ÷ 2 + 1
        for (i, x) in enumerate(xs), (j, y) in enumerate(ys)
            px = round(Int, (x - ox) / mi.pixsize) + cx
            py = round(Int, (y - oy) / mi.pixsize) + cy
            (1 <= px <= nx && 1 <= py <= ny) || continue
            conv[px, py] > thr * peak && (masked[i, j] = true)
        end
    end
    return masked
end

_num(md, k) = (try _resolve_numeric(k, md) catch; 0.0 end)

# ─────────────────────────────────────────────────────────────────────────────────────
# the fast path
# ─────────────────────────────────────────────────────────────────────────────────────

# Everything a node needs that does not depend on the node: the primary's visibilities A0
# (companion flux 0) and the companion's at the origin C0 (flux 1, the primary's 0), N × nsub;
# each point's uv scale at each sub-wavelength; the data the χ² compares with.
# In the data's precision `T`: a Float32 dataset is searched in Float32.
struct FastProblem{T<:AbstractFloat}
    A0       :: Matrix{Complex{T}}
    C0       :: Matrix{Complex{T}}
    S        :: Matrix{T}
    u        :: Vector{T}
    v        :: Vector{T}
    coherent :: Bool
    w        :: Vector{T}
    # V²
    i2 :: Vector{Int};  v2 :: Vector{T};  v2e :: Vector{T}
    # triple products
    t1 :: Vector{Int};  t2 :: Vector{Int};  t3 :: Vector{Int}
    t3a :: Vector{T};  t3ae :: Vector{T};  t3p :: Vector{T};  t3pe :: Vector{T}
    # VISAMP / VISPHI
    iv :: Vector{Int};  va :: Vector{T};  vae :: Vector{T}
    vp :: Vector{T};  vpe :: Vector{T}
end

# One thread's scratch: a node's companion visibilities and the polynomial coefficients of its
# observables in f. V² = a + 2b f + c f²; T3 = c0 + c1 f + c2 f² + c3 f³; VIS = α + β f.
struct FastWork{T<:AbstractFloat}
    C  :: Matrix{Complex{T}}
    a  :: Vector{T};  b :: Vector{T};  c :: Vector{T}
    c0 :: Vector{Complex{T}}; c1 :: Vector{Complex{T}}
    c2 :: Vector{Complex{T}}; c3 :: Vector{Complex{T}}
    α  :: Vector{Complex{T}}; β :: Vector{Complex{T}}
end
FastWork(fp::FastProblem{T}) where {T} = FastWork{T}(similar(fp.C0),
    zeros(T, length(fp.i2)), zeros(T, length(fp.i2)), zeros(T, length(fp.i2)),
    (zeros(Complex{T}, length(fp.t1)) for _ in 1:4)...,
    zeros(Complex{T}, length(fp.iv)), zeros(Complex{T}, length(fp.iv)))

function _fast_problem(md, data, cn, w)
    kf = cn * ",f"
    T = float(eltype(data.uv))
    m = dict_to_model(md, String[])
    smeared = _smeared(m, data)
    ns = smeared ? smear_samples(m, data) : 1
    S = smeared ? T.(sub_scales(data.uv_lam, data.uv_dlam, ns)) : ones(T, size(data.uv, 2), 1)
    function vis(f)
        mm = dict_to_model(merge(md, Dict{String,Any}(kf => f)), String[])
        V = smeared ? eval_model_sub(mm, T[], data.uv; wl = data.uv_lam,
                                     mjd = data.uv_mjd, dlam = data.uv_dlam, nsub = ns) :
                      reshape(eval_model(mm, T[], data.uv; wl = data.uv_lam,
                                         mjd = data.uv_mjd), :, 1)
        return Matrix{Complex{T}}(V)
    end
    A0, C0, H = vis(0.0), vis(1.0), vis(0.5)
    maximum(abs.(H .- (A0 .+ C0) ./ 2)) < 1e-6 * max(1.0, maximum(abs.(C0))) ||
        throw(ArgumentError(
            "the model is not linear in $kf, which the fast path needs: write the primary's " *
            "flux as 1 - \$$kf, or free another parameter to use the general path"))
    on(k, n) = w[k] > 0 && n > 0
    F(x) = T.(x)
    return FastProblem{T}(A0, C0, S, F(data.uv[1, :]), F(data.uv[2, :]),
        smeared && m.smear_coherent, T.(w),
        on(1, data.nv2) ? data.indx_v2 : Int[], F(data.v2), F(data.v2_err),
        on(2, data.nt3amp) || on(3, data.nt3phi) ? data.indx_t3_1 : Int[],
        data.indx_t3_2, data.indx_t3_3,
        F(data.t3amp), F(data.t3amp_err), F(data.t3phi), F(data.t3phi_err),
        on(4, data.nvisamp) || on(5, data.nvisphi) ? data.indx_vis : Int[],
        F(data.visamp), F(data.visamp_err), F(data.visphi), F(data.visphi_err))
end

# The coefficients at the node (x, y): the companion is its origin visibilities times its phase
# at each sub-wavelength; B = C − A. Averages over the sub-wavelengths follow the smearing mode:
# of the observables (incoherent), or of the visibilities first (coherent).
function _fast_coefficients!(ws::FastWork{T}, fp::FastProblem{T}, x, y) where {T}
    A, C, S = fp.A0, ws.C, fp.S
    N, ns = size(A)
    k2π = T(-2π / _MAS2RAD_PM); xT, yT = T(x), T(y)
    @inbounds for k in 1:ns, i in 1:N
        C[i, k] = fp.C0[i, k] * cis(k2π * (fp.u[i] * xT + fp.v[i] * yT) * S[i, k])
    end
    nk = fp.coherent ? 1 : ns
    Am(i, k) = fp.coherent ? _colmean(A, i) : A[i, k]
    Bm(i, k) = fp.coherent ? _colmean(C, i) - _colmean(A, i) : C[i, k] - A[i, k]
    @inbounds for (j, p) in enumerate(fp.i2)
        a = b = c = zero(T)
        for k in 1:nk
            Ak, Bk = Am(p, k), Bm(p, k)
            a += abs2(Ak); b += real(conj(Ak) * Bk); c += abs2(Bk)
        end
        ws.a[j], ws.b[j], ws.c[j] = a / nk, b / nk, c / nk
    end
    @inbounds for j in eachindex(fp.t1)
        p1, p2, p3 = fp.t1[j], fp.t2[j], fp.t3[j]
        c0 = c1 = c2 = c3 = zero(Complex{T})
        for k in 1:nk
            A1, A2, A3 = Am(p1, k), Am(p2, k), Am(p3, k)
            B1, B2, B3 = Bm(p1, k), Bm(p2, k), Bm(p3, k)
            c0 += A1 * A2 * A3
            c1 += B1 * A2 * A3 + A1 * B2 * A3 + A1 * A2 * B3
            c2 += A1 * B2 * B3 + B1 * A2 * B3 + B1 * B2 * A3
            c3 += B1 * B2 * B3
        end
        ws.c0[j], ws.c1[j], ws.c2[j], ws.c3[j] = c0 / nk, c1 / nk, c2 / nk, c3 / nk
    end
    @inbounds for (j, p) in enumerate(fp.iv)
        α = β = zero(Complex{T})
        for k in 1:ns
            α += A[p, k]; β += C[p, k] - A[p, k]
        end
        ws.α[j], ws.β[j] = α / ns, β / ns
    end
    return ws
end

_colmean(M, i) = (s = zero(eltype(M)); @inbounds for k in axes(M, 2); s += M[i, k]; end; s / size(M, 2))

_wrap180(d) = mod(d + 180, 360) - 180

# χ² at flux f from the coefficients, without allocating. The terms are `_chi2_terms`'.
function _fast_chi2(f, ws::FastWork{T}, fp::FastProblem{T}) where {T}
    w = fp.w
    f = T(f)
    χ = zero(T)
    if w[1] > 0
        @inbounds for j in eachindex(fp.i2)
            r = (ws.a[j] + f * (2ws.b[j] + f * ws.c[j]) - fp.v2[j]) / fp.v2e[j]
            χ += w[1] * r^2
        end
    end
    if (w[2] > 0 || w[3] > 0) && !isempty(fp.t1)
        @inbounds for j in eachindex(fp.t1)
            P = ws.c0[j] + f * (ws.c1[j] + f * (ws.c2[j] + f * ws.c3[j]))
            if w[2] > 0
                χ += w[2] * ((abs(P) - fp.t3a[j]) / fp.t3ae[j])^2
            end
            if w[3] > 0
                χ += w[3] * (_wrap180(angle(P) * T(180 / π) - fp.t3p[j]) / fp.t3pe[j])^2
            end
        end
    end
    if (w[4] > 0 || w[5] > 0) && !isempty(fp.iv)
        @inbounds for j in eachindex(fp.iv)
            V = ws.α[j] + f * ws.β[j]
            w[4] > 0 && (χ += w[4] * ((abs(V) - fp.va[j]) / fp.vae[j])^2)
            w[5] > 0 && (χ += w[5] * (_wrap180(angle(V) * T(180 / π) - fp.vp[j]) / fp.vpe[j])^2)
        end
    end
    return χ
end

function _fast_node(ws::FastWork{T}, fp::FastProblem{T}, x, y, fmax, solver) where {T}
    _fast_coefficients!(ws, fp, x, y)
    χ2 = f -> _fast_chi2(f, ws, fp)
    if solver === :halley
        f, ok = _halley(ws, fp, fmax)
        ok && return (f, χ2(f), false)
        f, c = _bradi(χ2, fmax, T)
        return (f, c, true)
    end
    f, c = _bradi(χ2, fmax, T)
    return (f, c, false)
end

function _bradi(χ2::F, fmax, ::Type{T} = Float64) where {F,T}
    xm, fm, _, _, _ = BraDi.minimize(T, χ2, range(zero(T), T(fmax); length = 9);
                                     atol = 4eps(T) * T(fmax))
    # BraDi samples the interior; the bound itself is a candidate (no companion at all).
    c0 = χ2(zero(T))
    return c0 <= fm ? (zero(T), c0) : (T(xm), T(fm))
end

# χ² and its first three derivatives in f, in closed form. Every observable is a polynomial in
# f, so each residual r and its derivatives follow directly:
#   V² = a + 2bf + cf²: dr/df = (2b + 2cf)/σ, d²r/df² = 2c/σ, d³r/df³ = 0;
#   a triple product T (cubic) or a visibility (linear), through L = log T, whose derivatives
#   are L₁ = T₁/T, L₂ = T₂/T − L₁², L₃ = T₃/T − 3L₁L₂ − L₁³: the phase is Im L, and the
#   amplitude |T| = exp(Re L), whose derivatives are |T|ρ₁, |T|(ρ₂ + ρ₁²) and
#   |T|(ρ₃ + 3ρ₁ρ₂ + ρ₁³) with ρₙ = Re Lₙ.
# Each term w·r² then contributes 2w·r r₁, 2w(r₁² + r r₂) and 2w(3r₁r₂ + r r₃).
# One term w·r² of χ², added with its derivatives to the running (χ, χ₁, χ₂, χ₃).
@inline _acc3(a, wk, r, r1, r2, r3) =
    (a[1] + wk * r^2, a[2] + 2wk * r * r1, a[3] + 2wk * (r1^2 + r * r2),
     a[4] + 2wk * (3r1 * r2 + r * r3))

# The amplitude and phase terms of a complex polynomial value P with derivatives P1, P2, P3.
@inline function _acc3_complex(a, wa, wp, P, P1, P2, P3, da, σa, dp, σp)
    L1 = P1 / P
    L2 = P2 / P - L1^2
    L3 = P3 / P - 3L1 * L2 - L1^3
    A = abs(P)
    ρ1, ρ2, ρ3 = real(L1), real(L2), real(L3)
    wa > 0 && (a = _acc3(a, wa, (A - da) / σa, A * ρ1 / σa, A * (ρ2 + ρ1^2) / σa,
                         A * (ρ3 + 3ρ1 * ρ2 + ρ1^3) / σa))
    k = real(typeof(P))(180 / π)
    wp > 0 && (a = _acc3(a, wp, _wrap180(angle(P) * k - dp) / σp, k * imag(L1) / σp,
                         k * imag(L2) / σp, k * imag(L3) / σp))
    return a
end

function _fast_chi2_d3(f, ws::FastWork{T}, fp::FastProblem{T}) where {T}
    w = fp.w
    f = T(f)
    z = zero(T)
    a = (z, z, z, z)
    if w[1] > 0
        @inbounds for j in eachindex(fp.i2)
            σ = fp.v2e[j]
            a = _acc3(a, w[1], (ws.a[j] + f * (2ws.b[j] + f * ws.c[j]) - fp.v2[j]) / σ,
                      (2ws.b[j] + 2ws.c[j] * f) / σ, 2ws.c[j] / σ, z)
        end
    end
    if w[2] > 0 || w[3] > 0
        @inbounds for j in eachindex(fp.t1)
            c0, c1, c2, c3 = ws.c0[j], ws.c1[j], ws.c2[j], ws.c3[j]
            P  = c0 + f * (c1 + f * (c2 + f * c3))
            P1 = c1 + f * (2c2 + 3c3 * f)
            a = _acc3_complex(a, w[2], w[3], P, P1, 2c2 + 6c3 * f, 6c3,
                              fp.t3a[j], fp.t3ae[j], fp.t3p[j], fp.t3pe[j])
        end
    end
    if w[4] > 0 || w[5] > 0
        @inbounds for j in eachindex(fp.iv)
            P = ws.α[j] + f * ws.β[j]
            a = _acc3_complex(a, w[4], w[5], P, ws.β[j], zero(P), zero(P),
                              fp.va[j], fp.vae[j], fp.vp[j], fp.vpe[j])
        end
    end
    return a
end

# Halley's method on χ²'(f) = 0, from the best of a coarse sampling of [0, fmax], with the
# derivatives above. Converged means: a step below tolerance, at a minimum (χ²'' > 0), inside
# the bounds.
# The tolerance is relative to fmax and never finer than the precision can resolve.
function _halley(ws::FastWork{T}, fp::FastProblem{T}, fmax; maxit = 30,
                 tol = max(1e-10, 8eps(T))) where {T}
    f = zero(T); best = _fast_chi2(f, ws, fp)
    fmax = T(fmax)
    for k in 1:8
        g = T(k * fmax / 8); c = _fast_chi2(g, ws, fp)
        c < best && (f = g; best = c)
    end
    for _ in 1:maxit
        _, g, g1, g2 = _fast_chi2_d3(f, ws, fp)
        den = 2g1^2 - g * g2
        (den == 0 || !isfinite(den)) && return (f, false)
        Δ = 2g * g1 / den
        f -= Δ
        (0 <= f <= fmax) || return (clamp(f, zero(T), fmax), false)
        abs(Δ) < tol * max(fmax, 1e-12) && return (f, _fast_chi2_d3(f, ws, fp)[3] > 0)
    end
    return (f, false)
end

# ─────────────────────────────────────────────────────────────────────────────────────
# the general path
# ─────────────────────────────────────────────────────────────────────────────────────

function _general_problem(md, data, cn, free, lb, ub, fmax, w)
    kx, ky, kf = cn * ",x", cn * ",y", cn * ",f"
    all = vcat([kx, ky], free)
    m = dict_to_model(md, all)
    dlb, dub = default_bounds(md, free; data)
    lo = [k == kf ? 0.0 : Float64(get(lb, k, get(dlb, k, -Inf))) for k in free]
    hi = [k == kf ? Float64(fmax) : Float64(get(ub, k, get(dub, k, Inf))) for k in free]
    p0 = [k == kf ? 0.01 * fmax : Float64(md[k]) for k in free]
    p0 = clamp.(p0, lo, hi)
    return (; m, data, w, lo, hi, p0, ifree = findfirst(==(kf), free))
end

function _general_node(gp, x, y; pstart = gp.p0, fixf = false, cache = ModelJacCache())
    function crit(p, g)
        c, gfull = model_to_chi2_fg(gp.m, vcat(x, y, p), gp.data; weights = gp.w,
                                    jac_cache = cache)
        g .= gfull[3:end]
        fixf && (g[gp.ifree] = 0)
        return c
    end
    p = OptimPackNextGen.vmlmb(crit, copy(pstart); lower = gp.lo, upper = gp.hi,
                               maxiter = 500, verb = false, blmvm = false)
    return p, model_to_chi2(gp.m, vcat(x, y, p), gp.data; weights = gp.w)
end

# Without a companion: f = 0, the other free parameters fitted.
function _general_null(gp)
    p0 = copy(gp.p0); p0[gp.ifree] = 0.0
    lo = copy(gp.lo); hi = copy(gp.hi); lo[gp.ifree] = 0.0; hi[gp.ifree] = 0.0
    gq = merge(gp, (; lo, hi))
    _, c = _general_node(gq, 0.0, 0.0; pstart = p0, fixf = true)
    return c
end
