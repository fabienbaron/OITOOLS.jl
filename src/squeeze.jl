# ─────────────────────────────────────────────────────────────────────────────
# squeeze.jl — SQUEEZE-style MCMC image reconstruction.
#
# The image is not a pixel-intensity vector: it is a bag of `nelements` discrete
# flux quanta ("elements"), each carrying flux 1/nelements, living on an nx×nx
# integer pixel grid.  The dense image is the histogram of element positions.
#
# This buys, for free:
#   * positivity and fixed total flux (a histogram is >= 0 and sums to nelements)
#   * O(nuv) rank-1 dchi2 per move — moving one quantum is a single subtract-add
#     on a maintained model-visibility vector
#   * non-convex regularizers (L0) — no gradient is ever taken
#   * a genuine posterior (mean image), not just a MAP point
#
# Reference: SQUEEZE (F. Baron), ~/SOFTWARE/squeeze/src/squeeze.c
# ─────────────────────────────────────────────────────────────────────────────

using Random, Printf
using SpecialFunctions: besselj1, logabsgamma

# ── Pixel indexing convention ────────────────────────────────────────────────
#
# `setup_dft` (oichi2.jl) builds a (nuv, nx, nx) array whose dim 2 is the x span
# and dim 3 the y span, then reshapes to (nuv, nx^2).  Column-major reshape means
# column p corresponds to (ix, iy) with
#
#     p = (iy - 1) * nx + ix
#
# which is exactly Julia's linear index into an nx×nx matrix `histogram[ix, iy]`.
# `image_to_vis` uses `vec(x)`, so this matches the rest of OITOOLS bit for bit.

@inline _pixel_index(ix::Integer, iy::Integer, nx::Integer) = (iy - 1) * nx + ix

# Torus wrap: 1-based modular arithmetic.  This is what makes the SMALL / MEDIUM /
# ANYWHERE proposals symmetric, hence Metropolis-valid without a Hastings ratio.
@inline _wrap(i::Integer, nx::Integer) = mod1(i, nx)

# ─────────────────────────────────────────────────────────────────────────────
# Constants — ported verbatim from squeeze.h so the annealing behaviour matches.
# ─────────────────────────────────────────────────────────────────────────────

const SQ_DAMPING_TIME             = 25.0    # acceptance EWMA timescale
const SQ_TEMP_CHANGE_TIME         = 1000.0
const SQ_FLAT_CHI2_MULT           = 3.0
const SQ_BURN_IN_FRAC             = 5
const SQ_MPROB_LOW                = 0.2
const SQ_MPROB_HIGH               = 0.45
const SQ_TARGET_MPROB             = 0.3
const SQ_INITIAL_SCALED_CHI2_MULT = 2.0
const SQ_INITIAL_TMIN             = 1.0
const SQ_PARAM_DAMPING_TIME       = 10.0   # squeeze.h:85
const SQ_STEPSIZE_ADJUST_TIME     = 20.0   # squeeze.h:88

# Step types.  Ordering matters: the adaptation increments/decrements this.
const STEP_SMALL    = 0
const STEP_MEDIUM   = 1
const STEP_ANYWHERE = 2
const STEP_COPYCAT  = 3
const STEP_MAX      = 3

"""
    default_resync(T) -> Int

Number of accepted moves after which the maintained `raw_im_vis` is recomputed
from scratch, to bound accumulated round-off in the rank-1 subtract-add.

The relative drift after `N` accepted moves is `~ eps(T) * sqrt(2N)`; holding it
under `δ = 1e-5` (about 1% of a tight V² error bar) gives `N = δ²/(2 eps(T)²)`.
That is ~1e21 for Float64 (unreachable — never re-syncs in practice) and ~3.5e3
for Float32.
"""
function default_resync(::Type{T}) where {T<:AbstractFloat}
    n = 1e-10 / (2 * Float64(eps(T))^2)
    return n >= typemax(Int64) ? typemax(Int64) : max(round(Int64, n), 100)
end

"""
    SqueezeState{T}

Per-chain sampler state.  Everything is thread-private and preallocated; the
proposal loop allocates nothing.  The DFT matrix `A` is *not* stored here — it is
read-only and shared across all chains.
"""
mutable struct SqueezeState{T<:AbstractFloat}
    nx::Int
    nelements::Int
    nuv::Int

    # Element positions.  Int32 rather than C's UInt16: narrow ints buy nothing
    # here and cost widening conversions on every index computation.
    element_x::Vector{Int32}
    element_y::Vector{Int32}
    element_p::Vector{Int32}          # cached linear column index into A

    histogram::Matrix{Int32}          # exact integer counts, [ix, iy]

    # Maintained visibilities and their proposal twins.  Accept swaps references.
    raw_im_vis::Vector{Complex{T}}    # = A * vec(histogram), unnormalised
    new_raw_im_vis::Vector{Complex{T}}
    mod_vis::Vector{Complex{T}}
    new_mod_vis::Vector{Complex{T}}

    # Parametric model contribution (zero / all-ones without a model)
    param_vis::Vector{Complex{T}}
    fluxratio_image::Vector{T}
    new_param_vis::Vector{Complex{T}}
    new_fluxratio_image::Vector{T}

    # Scalars
    lLikelihood::Float64
    lPrior::Float64
    temperature::Float64
    prob_movement::Float64            # acceptance EWMA, adapted toward ~0.3
    naccepted::Int64
    nproposed::Int64
    since_resync::Int64
end

"""
    SqueezeState(T, nx, nelements, nuv)

Allocate an empty chain state.  Element positions are left at zero; use
[`init_elements!`](@ref) to populate them.
"""
function SqueezeState(::Type{T}, nx::Integer, nelements::Integer, nuv::Integer) where {T<:AbstractFloat}
    return SqueezeState{T}(
        nx, nelements, nuv,
        zeros(Int32, nelements), zeros(Int32, nelements), zeros(Int32, nelements),
        zeros(Int32, nx, nx),
        zeros(Complex{T}, nuv), zeros(Complex{T}, nuv),
        zeros(Complex{T}, nuv), zeros(Complex{T}, nuv),
        zeros(Complex{T}, nuv), ones(T, nuv),
        zeros(Complex{T}, nuv), ones(T, nuv),
        0.0, 0.0, 1.0, 0.3, 0, 0, 0)
end

const SQUEEZE_STARTS = (:point_source, :random)

"""
    init_elements!(s, rng; x_start=:point_source)

Place the `nelements` quanta on the grid and rebuild `histogram`, `element_p`.

`x_start` selects the starting configuration:

| value | meaning |
|---|---|
| `:point_source` | all quanta on the **centre pixel** — SQUEEZE's default when `-i` is not given (a Dirac, `squeeze.c:496-501`) |
| `:random` | quanta scattered uniformly over the grid (an OITOOLS addition; useful for giving multiple chains genuinely different starts) |
| an `nx × nx` image | digitised into quanta, the equivalent of C's `-i` |

Strings are accepted for the symbols (`"point_source"`, `"random"`).

Image digitisation follows C's deterministic flux-accumulation scan
(`squeeze.c:450-474`): walk the pixels in column-major order accumulating
`pixel · nelements`, emitting one quantum each time the accumulator passes 1.
Being deterministic, it gives every chain the same start — as in C, which copies
chain 0's elements to all the others.  `:random` is the option to reach for when
you want the chains to explore independently.
"""
function init_elements!(s::SqueezeState, rng::AbstractRNG; x_start = :point_source)
    nx = s.nx
    fill!(s.histogram, 0)

    start = x_start isa AbstractString ? Symbol(x_start) : x_start

    if start isa Symbol
        if start === :point_source
            # C uses axis_len/2 in 0-based indexing (squeeze.c:499); centpos in
            # cent_change is 0.5*axis_len to match, "so the starting dirac falls
            # on a pixel for an even number of pixels".
            c = nx ÷ 2 + 1
            for i in 1:s.nelements
                _place!(s, i, c, c)
            end
        elseif start === :random
            for i in 1:s.nelements
                _place!(s, i, rand(rng, 1:nx), rand(rng, 1:nx))
            end
        else
            throw(ArgumentError("unknown x_start $(repr(start)); use one of " *
                                "$(SQUEEZE_STARTS), or an $(nx)x$(nx) image"))
        end
        return s
    end

    start isa AbstractMatrix || throw(ArgumentError(
        "x_start must be one of $(SQUEEZE_STARTS) or an $(nx)x$(nx) image, " *
        "got $(typeof(x_start))"))
    size(start) == (nx, nx) ||
        throw(DimensionMismatch("x_start is $(size(start)), expected ($nx, $nx)"))
    any(<(0), start) && throw(ArgumentError("x_start image must be non-negative"))
    tot = sum(start)
    tot > 0 || throw(ArgumentError("x_start image must have positive total flux"))

    # Normalised so one full scan emits exactly `nelements` quanta.  (C does not
    # normalise and instead wraps around the image until it has enough, which is
    # the same thing for a unit-flux image and merely slower otherwise.)
    scale = s.nelements / Float64(tot)
    acc = 0.0
    n = 0
    p = 0
    npix = nx * nx
    while n < s.nelements
        if acc >= 1.0
            acc -= 1.0
            n += 1
            iy, ix = divrem(p - 1, nx)
            _place!(s, n, ix + 1, iy + 1)
        else
            p += 1
            p > npix && (p = 1)          # wrap, as C does
            acc += Float64(start[p]) * scale
        end
    end
    return s
end

@inline function _place!(s::SqueezeState, i::Integer, ix::Integer, iy::Integer)
    s.element_x[i] = ix
    s.element_y[i] = iy
    s.element_p[i] = _pixel_index(ix, iy, s.nx)
    s.histogram[ix, iy] += Int32(1)
    return nothing
end

"""
    resync_raw_im_vis!(s, A)

Recompute `raw_im_vis` from scratch as `A * vec(histogram)`, discarding any
round-off accumulated by the rank-1 updates.  O(nuv * nx²).
"""
function resync_raw_im_vis!(s::SqueezeState{T}, A::AbstractMatrix{Complex{T}}) where {T}
    fill!(s.raw_im_vis, zero(Complex{T}))
    h = s.histogram
    @inbounds for p in eachindex(h)
        c = h[p]
        c == 0 && continue
        f = T(c)
        col = view(A, :, p)
        @simd for k in eachindex(s.raw_im_vis)
            s.raw_im_vis[k] += f * col[k]
        end
    end
    s.since_resync = 0
    return s.raw_im_vis
end

"""
    rank1_update!(dest, src, A, p_new, p_old)

`dest .= src .+ A[:, p_new] .- A[:, p_old]` — the O(nuv) core of the sampler.

`A` is `nuv × nx²` so its columns are contiguous in Julia's column-major layout;
this is a unit-stride SIMD loop with no gather and no allocation.
"""
@inline function rank1_update!(dest::Vector{Complex{T}}, src::Vector{Complex{T}},
                               A::AbstractMatrix{Complex{T}},
                               p_new::Integer, p_old::Integer) where {T}
    @inbounds begin
        anew = view(A, :, p_new)
        aold = view(A, :, p_old)
        @simd for k in eachindex(dest)
            dest[k] = src[k] + anew[k] - aold[k]
        end
    end
    return dest
end

"""
    form_mod_vis!(dest, raw, param_vis, fluxratio, nelements)

`dest .= param_vis .+ raw .* fluxratio_image ./ nelements` — the model complex
visibility handed to the likelihood.  Mirrors SQUEEZE's
`compute_model_visibilities_fromelements`.
"""
@inline function form_mod_vis!(dest::Vector{Complex{T}}, raw::Vector{Complex{T}},
                               pv::Vector{Complex{T}}, fr::Vector{T},
                               nelements::Integer) where {T}
    invn = one(T) / T(nelements)
    @inbounds @simd for k in eachindex(dest)
        dest[k] = pv[k] + raw[k] * (fr[k] * invn)
    end
    return dest
end

@inline form_mod_vis!(dest::Vector{Complex{T}}, raw::Vector{Complex{T}},
                      s::SqueezeState{T}) where {T} =
    form_mod_vis!(dest, raw, s.param_vis, s.fluxratio_image, s.nelements)

# ─────────────────────────────────────────────────────────────────────────────
# Regularizers.
#
# Ported directly from regularizations.c.  These act on the INTEGER histogram
# with a `/flux` normalisation (flux = nelements), which is a different
# convention from OITOOLS' differentiable regularizers in oichi2.jl — those are
# defined on the normalised pixel image and are gradient-oriented.  Do not
# substitute one for the other.
# ─────────────────────────────────────────────────────────────────────────────

"""
    SqueezeRegs

Regularizer weights, their running values, and the mutable centroid state that
the centering penalty needs.  `lPrior = Σ λ_i · value_i`.
"""
mutable struct SqueezeRegs
    λ_l0::Float64
    λ_tv::Float64
    λ_entropy::Float64
    λ_compactness::Float64
    λ_centering::Float64

    v_l0::Float64
    v_tv::Float64
    v_entropy::Float64
    v_compactness::Float64
    v_centering::Float64
    λ_modelparam::Float64
    v_modelparam::Float64
    λ_priorimage::Float64
    v_priorimage::Float64
    prior_penalty::Vector{Float64}    # per-pixel -log(prior); empty when unused

    centroid_x::Float64               # running displacement from the start
    centroid_y::Float64
    fov::Float64
    cent_mult::Float64
    tv_eps::Float64
end

SqueezeRegs(; l0=0.0, tv=0.0, entropy=0.0, compactness=0.0, centering=0.0,
              priorimage=0.0, prior_penalty=Float64[],
              fov=1.0, cent_mult=1.0, tv_eps=0.0) =
    SqueezeRegs(l0, tv, entropy, compactness, centering,
                0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, priorimage, 0.0, prior_penalty,
                0.0, 0.0, fov, cent_mult, tv_eps)

@inline _any_spatial(r::SqueezeRegs) =
    r.λ_l0 > 0 || r.λ_tv > 0 || r.λ_entropy > 0 || r.λ_compactness > 0

@inline lprior(r::SqueezeRegs) =
    r.λ_l0 * r.v_l0 + r.λ_tv * r.v_tv + r.λ_entropy * r.v_entropy +
    r.λ_compactness * r.v_compactness + r.λ_centering * r.v_centering +
    r.λ_modelparam * r.v_modelparam + r.λ_priorimage * r.v_priorimage

# --- L0: count of nonzero pixels (regularizations.c:256).  Not flux-normalised.
reg_l0(h::AbstractMatrix{<:Integer}) = Float64(count(>(0), h))

# --- Entropy: SQUEEZE's `entropy_full` (regularizations.c:28), a quantised /
#     Poisson prior Σ_{x>0} lgamma(x).  This is NOT OITOOLS' `entropy`
#     (Σ x log x, oichi2.jl:574) — different prior, different scaling.
@inline _ent(s::Real) = s < 1e-7 ? 0.0 : logabsgamma(Float64(s))[1]

function reg_entropy(h::AbstractMatrix{<:Integer})
    r = 0.0
    @inbounds for i in eachindex(h)
        h[i] > 0 && (r += _ent(h[i]))
    end
    return r
end

# --- TV: regularizations.c:149.  Backward differences, image borders ignored
#     (the C loops start at 1), offset by -sqrt(eps)*nx*ny, divided by flux.
function reg_tv(h::AbstractMatrix{<:Integer}, flux::Real, eps::Float64)
    nx, ny = size(h)
    L1g = -sqrt(eps) * (nx * ny)
    @inbounds for j in 2:ny, i in 2:nx
        dx = Float64(h[i, j]) - Float64(h[i-1, j])
        dy = Float64(h[i, j]) - Float64(h[i, j-1])
        L1g += sqrt(dx * dx + dy * dy + eps * eps)
    end
    return L1g / flux
end

# --- Compactness: regularizations.c:199.  Σ radial²·x², not flux-normalised.
@inline function _compact_w(i::Integer, j::Integer, nx::Integer, ny::Integer)
    di = i - 0.5 * (nx - 1)
    dj = j - 0.5 * (ny - 1)
    return (di * di + dj * dj) / (nx * ny)
end

function reg_compactness(h::AbstractMatrix{<:Integer})
    nx, ny = size(h)
    r = 0.0
    @inbounds for j in 1:ny, i in 1:nx
        c = Float64(h[i, j])
        c == 0 && continue
        r += _compact_w(i - 1, j - 1, nx, ny) * c * c
    end
    return r
end

"""
    cent_change!(r, new_x, new_y, old_x, old_y, nx)

Centroid / second-moment penalty *delta* (regularizations.c:6).  Mutates the
running centroid, so calling it with `(old, new)` swapped exactly reverses it —
which is how a rejected move is undone.
"""
@inline function cent_change!(r::SqueezeRegs, new_x::Integer, new_y::Integer,
                              old_x::Integer, old_y::Integer, nx::Integer)
    old_sumsqr = r.centroid_x^2 + r.centroid_y^2
    r.centroid_x += Float64(new_x - old_x)
    r.centroid_y += Float64(new_y - old_y)
    centpos = 0.5 * nx
    new_sumsqr = r.centroid_x^2 + r.centroid_y^2
    fov = r.fov
    return ((new_sumsqr - old_sumsqr) * r.cent_mult +
            ((new_x - centpos)^2 - (old_x - centpos)^2 +
             (new_y - centpos)^2 - (old_y - centpos)^2) * fov * fov) *
           4.0 / (nx * nx)
end

"""
    recompute_regs!(r, s)

Recompute every enabled regularizer's running value from the histogram, and
reset the running centroid.  Called at chain start and by the drift checks.
"""
function recompute_regs!(r::SqueezeRegs, s::SqueezeState)
    h = s.histogram
    nx = s.nx
    r.λ_l0           > 0 && (r.v_l0           = reg_l0(h))
    r.λ_entropy      > 0 && (r.v_entropy      = reg_entropy(h))
    r.λ_tv           > 0 && (r.v_tv           = reg_tv(h, s.nelements, r.tv_eps))
    r.λ_compactness  > 0 && (r.v_compactness  = reg_compactness(h))
    if r.λ_priorimage > 0
        acc = 0.0
        @inbounds for i in 1:s.nelements
            acc += r.prior_penalty[s.element_p[i]]
        end
        r.v_priorimage = acc
    end
    if r.λ_centering > 0
        # C initialises centroid_image_x/y to 0 and reg_value[REG_CENTERING] to 0
        # (squeeze.c:747, :752), and `compute_regularizers` never touches CENTERING.
        # So the centroid is a running *displacement from the starting configuration*,
        # not an absolute position — starting it from Σx would make the penalty
        # O(nelements²·nx²) instead of O(1) and swamp the likelihood entirely.
        r.centroid_x = 0.0
        r.centroid_y = 0.0
        r.v_centering = 0.0
    end
    return r
end

# ─────────────────────────────────────────────────────────────────────────────
# SPARCO parametric model (Kluska et al. 2012) + over-resolved background.
#
# Ported from ~/SOFTWARE/squeeze/src/models/modelcode_sparco.c.  The image (the
# "environment") is reconstructed grey; all the chromatism lives in this model.
# Per uv point:
#
#            f_star (λ/λ0)^-4 V_star  +  f_bg (λ/λ0)^bg_ind · 0  +  f_env V_env
#   V_tot = ────────────────────────────────────────────────────────────────────
#            f_star (λ/λ0)^-4         +  f_bg (λ/λ0)^bg_ind      +  f_env
#
# The background is over-resolved (V = 0): it contributes flux but no visibility,
# so it only dilutes.  `param_vis` gets the star term; `fluxratio_image` gets the
# environment's flux fraction, which is what scales the reconstructed image.
# ─────────────────────────────────────────────────────────────────────────────

const MAS_RAD = 206264806.2          # squeeze.h:57
const NPARAMS_SPARCO = 6
const SPARCO_PARAM_NAMES = (:f_star, :ud, :env_indx, :lambda0, :f_bg, :bg_indx)

"""
    SqueezeSparco(; f_star, ud, env_indx, lambda0, f_bg, bg_indx, free, stepsize)

SPARCO model state for [`reconstruct_squeeze`](@ref).

Parameters, in C's order:

| # | name | meaning |
|---|---|---|
| 1 | `f_star`   | stellar flux fraction at `lambda0` |
| 2 | `ud`       | uniform-disc diameter of the star, mas (`0` = point source) |
| 3 | `env_indx` | flux power-law index of the environment (the image) |
| 4 | `lambda0`  | reference wavelength, **metres**, must be > 0 |
| 5 | `f_bg`     | over-resolved background flux fraction at `lambda0` |
| 6 | `bg_indx`  | flux power-law index of the background |

`free` lists the parameters the sampler may vary, by name or index; everything
else is held fixed.  `stepsize` is the two-point lattice half-width per parameter
and is adapted toward a 0.3 acceptance rate during the run.

The model is only meaningful on **polychromatic** data: with a single wavelength
every `(λ/λ0)^x` is a constant and the three indices are exactly degenerate.
"""
mutable struct SqueezeSparco
    params::Vector{Float64}
    stepsize::Vector{Float64}
    free::Vector{Int}
    prob_pmovement::Vector{Float64}
end

function SqueezeSparco(; f_star::Real = 0.5, ud::Real = 0.0, env_indx::Real = 0.0,
                         lambda0::Real, f_bg::Real = 0.0, bg_indx::Real = 0.0,
                         free = (:f_star,),
                         stepsize = nothing)
    lambda0 > 0 || throw(ArgumentError("SPARCO lambda0 must be > 0 (metres); got $lambda0"))
    p = Float64[f_star, ud, env_indx, lambda0, f_bg, bg_indx]

    idx = Int[]
    for f in free
        if f isa Integer
            1 <= f <= NPARAMS_SPARCO ||
                throw(ArgumentError("SPARCO parameter index $f out of range 1:$NPARAMS_SPARCO"))
            push!(idx, Int(f))
        else
            j = findfirst(==(Symbol(f)), SPARCO_PARAM_NAMES)
            j === nothing && throw(ArgumentError(
                "unknown SPARCO parameter $(f); valid names: $(SPARCO_PARAM_NAMES)"))
            push!(idx, j)
        end
    end
    allunique(idx) || throw(ArgumentError("duplicate entries in `free`: $free"))

    # Defaults chosen so each parameter moves on its own natural scale.
    ss = stepsize === nothing ?
         Float64[0.01, 0.05, 0.05, 0.0, 0.01, 0.05] :
         collect(Float64, stepsize)
    length(ss) == NPARAMS_SPARCO ||
        throw(ArgumentError("stepsize must have $NPARAMS_SPARCO entries, got $(length(ss))"))
    for j in idx
        ss[j] > 0 || throw(ArgumentError(
            "parameter $(SPARCO_PARAM_NAMES[j]) is free but its stepsize is $(ss[j]); " *
            "a zero stepsize proposes no move"))
    end

    return SqueezeSparco(p, ss, idx, fill(SQ_TARGET_MPROB, NPARAMS_SPARCO))
end

"""
    sparco_vis!(param_vis, fluxratio_image, params, data) -> lPriorModel

Fill the star term and the image's flux fraction for every uv point, and return
the box-prior penalty (`1e99` per violated bound, exactly as C does — an
unacceptable proposal rather than a hard error).
"""
function sparco_vis!(param_vis::Vector{Complex{T}}, fluxratio_image::Vector{T},
                     params::Vector{Float64}, data::OIdata{T}) where {T}
    f_star_ref, ud, env_ind, lambda_ref, f_bg_ref, bg_ind = params

    @inbounds for i in 1:data.nuv
        # Uniform-disc star; ud <= 0 means a point source (V = 1).
        vis_star = 1.0
        if ud > 0
            b = sqrt(Float64(data.uv[1,i])^2 + Float64(data.uv[2,i])^2)
            t = pi * ud / MAS_RAD * b + 1e-15
            vis_star = 2.0 * besselj1(t) / t
        end
        r = Float64(data.uv_lam[i]) / lambda_ref
        f_star = f_star_ref * r^(-4.0)
        f_bg   = f_bg_ref * r^bg_ind
        # NOTE: C uses the *wavelength-scaled* f_bg here, not f_bg_ref
        # (modelcode_sparco.c).  The two agree at λ = λ0 and diverge away from it.
        # Ported as-is; see the handoff note.
        f_env  = (1.0 - f_star_ref - f_bg) * r^env_ind
        tot = f_star + f_env + f_bg
        fluxratio_image[i] = T(f_env / tot)
        param_vis[i] = Complex{T}(f_star * vis_star / tot)
    end

    # Box priors (modelcode_sparco.c): fluxes >= 0, UD >= 0, indices in (-5, 5].
    lp = 0.0
    params[1] >= 0 || (lp += 1e99)
    params[2] >= 0 || (lp += 1e99)
    params[5] >= 0 || (lp += 1e99)
    (-5 < params[3] <= 5) || (lp += 1e99)
    (-5 < params[6] <= 5) || (lp += 1e99)
    # No prior on lambda0: it is validated (> 0) at construction, since we divide by it.
    return lp
end

# ─────────────────────────────────────────────────────────────────────────────
# Live monitoring (optional).
#
# Off by default and free when off: the only cost on the sampling path is one
# integer comparison per sweep.  When on, the plotting itself is what costs
# (~40 ms per update), so keep `monitor` large enough that updates are rare
# compared with a sweep.
#
# Anything that draws must run on the main thread: a Python call from a spawned
# task segfaults inside PythonCall. `reconstruct_squeeze` therefore pins the one
# monitored chain to the calling thread and spawns the rest;
# `reconstruct_squeeze_tempered` has no such split and runs single-threaded.
# ─────────────────────────────────────────────────────────────────────────────

"""
    SqueezeMonitor

Live-display state for a SQUEEZE run.  Holds the update interval, the figure
titles (so successive updates reuse the same windows rather than opening new
ones), and the history needed for the trace plots.
"""
mutable struct SqueezeMonitor
    every::Int
    pixsize::Float64
    title::String
    colormap::String
    iters::Vector{Int}
    chi2r::Vector{Float64}
    temperature::Vector{Float64}
    params::Vector{Vector{Float64}}
    free::Vector{Int}
    ndf::Float64               # so the trace panel can show reduced chi2
    trace2_label::String       # what the second panel is: "T" or "round"
end

function SqueezeMonitor(every::Integer; pixsize::Real = -1.0,
                        title::AbstractString = "SQUEEZE",
                        colormap::AbstractString = "gist_heat",
                        free::Vector{Int} = Int[],
                        trace2_label::AbstractString = "T")
    return SqueezeMonitor(Int(every), Float64(pixsize), String(title), String(colormap),
                          Int[], Float64[], Float64[], Vector{Float64}[], free, 1.0,
                          String(trace2_label))
end

"""
    monitor_update!(mon, image, iter; chi2r, temperature, params)

Redraw the live displays: the current image and a trace panel with χ²r, temperature (or
Pigeons round) and each free model parameter against iteration.

The drawing lives in `oiplot.jl`, i.e. in the PythonPlot extension — `SqueezeMonitor` itself
holds only plain data, so the sampler needs no plotting stack unless monitoring is on. Load
PythonPlot to enable this method.
"""
function monitor_update! end


# ─────────────────────────────────────────────────────────────────────────────
# The sampler.
# ─────────────────────────────────────────────────────────────────────────────

"""
    default_nelements(data, nx)

C's heuristic (squeeze.c:298): `2·ceil(nx · ndata^(1/3))`, floored at 500.
"""
function default_nelements(data::OIdata, nx::Integer)
    ndata = data.nv2 + data.nvisamp + data.nvisphi + data.nt3amp + data.nt3phi
    n = 2 * ceil(Int, nx * cbrt(max(ndata, 1)))
    return max(n, 500)
end

"""
    flat_chi2(data, weights; seed=1)

χ² of a near-zero-visibility model — SQUEEZE's `get_flat_chi2` (squeeze.c:1895),
which uses `exp(-10 + i·φ)` with random phase.  Used only to cap the initial
annealing temperature.
"""
function flat_chi2(data::OIdata{T}, w::NTuple{7,Float64}; seed::Integer = 1) where {T}
    rng = Xoshiro(seed)
    V = [Complex{T}(exp(-10.0) * cis(rand(rng, 0:999) * 0.006283)) for _ in 1:data.nuv]
    return _cvis_to_chi2(V, data, w, false)
end

"""
    squeeze_ndof(data, weights; nparams_free=0, centering=false)

`ndf` for the annealing schedule: the weighted data count (via OITOOLS' `_ndof`)
plus free model parameters, plus 2 when centering is active — matching C.
"""
function squeeze_ndof(data::OIdata, weights; nparams_free::Integer = 0, centering::Bool = false)
    return _ndof(data, weights) + nparams_free + (centering ? 2 : 0)
end

"""
    _draw_step(rng, steptype, s, i, nx) -> (xstep, ystep)

Draw a proposal displacement of the given type, redrawing until it is non-zero
(squeeze.c:1010-1055).  SMALL/MEDIUM/ANYWHERE are symmetric on the torus, so no
Hastings ratio is needed; COPYCAT is not (see `reconstruct_squeeze`).
"""
@inline function _draw_step(rng::AbstractRNG, steptype::Int, s::SqueezeState,
                            i::Integer, nx::Integer)
    while true
        xstep = 0; ystep = 0
        if steptype == STEP_SMALL
            dir = rand(rng, 0:3)
            if dir > 1
                ystep = (dir - 2) * 2 - 1
            else
                xstep = dir * 2 - 1
            end
        elseif steptype == STEP_MEDIUM
            xstep = rand(rng, -2:2)
            ystep = rand(rng, -2:2)
        elseif steptype == STEP_ANYWHERE
            xstep = rand(rng, 0:nx-1) - nx ÷ 2
            ystep = rand(rng, 0:nx-1) - nx ÷ 2
        else # STEP_COPYCAT — land on a random donor element's pixel
            d = rand(rng, 1:s.nelements)
            xstep = s.element_x[d] - s.element_x[i]
            ystep = s.element_y[d] - s.element_y[i]
        end
        (xstep != 0 || ystep != 0) && return (xstep, ystep)
    end
end

"""
    _proposal_regs!(regs, s, ix_old, iy_old, ix_new, iy_new) -> (newvals, new_lPrior)

Regularizer values for the proposed state.  Entropy, L0, compactness, the prior
image and centering have exact O(1) deltas; TV is recomputed the way C does, with
the histogram temporarily moved.  Centering mutates the running centroid and MUST
be reversed by the caller on rejection.
"""
function _proposal_regs!(regs::SqueezeRegs, s::SqueezeState,
                         ix_old::Integer, iy_old::Integer,
                         ix_new::Integer, iy_new::Integer)
    h = s.histogram
    c_old = h[ix_old, iy_old]
    c_new = h[ix_new, iy_new]

    v_l0 = regs.v_l0
    v_ent = regs.v_entropy
    v_tv = regs.v_tv
    v_cp = regs.v_compactness
    v_ct = regs.v_centering
    v_pi = regs.v_priorimage

    if regs.λ_priorimage > 0
        # Exact O(1) delta: only this element's pixel changed (squeeze.c:1073).
        v_pi = regs.v_priorimage -
               regs.prior_penalty[_pixel_index(ix_old, iy_old, s.nx)] +
               regs.prior_penalty[_pixel_index(ix_new, iy_new, s.nx)]
    end
    if regs.λ_entropy > 0
        v_ent = regs.v_entropy - _ent(c_old) + _ent(c_old - 1) +
                                 _ent(c_new + 1) - _ent(c_new)
    end
    if regs.λ_l0 > 0
        v_l0 = regs.v_l0 + (c_new < 1 ? 1.0 : 0.0) - (c_old == 1 ? 1.0 : 0.0)
    end
    if regs.λ_compactness > 0
        nx, ny = size(h)
        w_old = _compact_w(ix_old - 1, iy_old - 1, nx, ny)
        w_new = _compact_w(ix_new - 1, iy_new - 1, nx, ny)
        v_cp = regs.v_compactness +
               w_old * ((c_old - 1)^2 - c_old^2) + w_new * ((c_new + 1)^2 - c_new^2)
    end
    if regs.λ_tv > 0
        # No cheap exact delta in C either: move, recompute, move back.
        h[ix_old, iy_old] -= Int32(1)
        h[ix_new, iy_new] += Int32(1)
        v_tv = reg_tv(h, s.nelements, regs.tv_eps)
        h[ix_old, iy_old] += Int32(1)
        h[ix_new, iy_new] -= Int32(1)
    end
    if regs.λ_centering > 0
        v_ct = regs.v_centering +
               regs.fov * cent_change!(regs, ix_new, iy_new, ix_old, iy_old, s.nx)
    end

    return (v_l0, v_tv, v_ent, v_cp, v_ct, v_pi),
           regs.λ_l0 * v_l0 + regs.λ_tv * v_tv + regs.λ_entropy * v_ent +
           regs.λ_compactness * v_cp + regs.λ_centering * v_ct +
           regs.λ_priorimage * v_pi + regs.λ_modelparam * regs.v_modelparam
end

"""
    run_chain!(s, regs, A, data, w, rng; kwargs...) -> NamedTuple

One independent annealing chain.  Returns its diagnostics plus the accumulated
post-burn-in mean image.
"""
function run_chain!(s::SqueezeState{T}, regs::SqueezeRegs,
                    A::AbstractMatrix{Complex{T}}, data::OIdata{T},
                    w::NTuple{7,Float64}, rng::AbstractRNG;
                    niter::Int, ndf::Float64, chi2_temp::Float64,
                    chi2_target::Float64, tmin::Float64,
                    f_anywhere::Float64, f_copycat::Float64,
                    resync_every::Int64, vonmises::Bool,
                    flatchi2::Float64, model::Union{Nothing,SqueezeSparco} = nothing,
                    verb::Bool = false, print_every::Int = 0, chain::Int = 1,
                    monitor::Union{Nothing,SqueezeMonitor} = nothing) where {T}

    nx = s.nx
    nelements = s.nelements
    nparams_free = model === nothing ? 0 : length(model.free)
    # C draws a target from [0, nelements + nparams_free*PARAMS_PER_ELT), so each
    # free parameter is proposed twice as often as a single element (squeeze.c:1005).
    ntargets = nelements + 2 * nparams_free

    # ── Initial state ─────────────────────────────────────────────────────────
    if model !== nothing
        regs.v_modelparam = sparco_vis!(s.param_vis, s.fluxratio_image, model.params, data)
        regs.λ_modelparam = 1.0
        copyto!(s.new_param_vis, s.param_vis)
        copyto!(s.new_fluxratio_image, s.fluxratio_image)
    end
    resync_raw_im_vis!(s, A)
    form_mod_vis!(s.mod_vis, s.raw_im_vis, s)
    s.lLikelihood = 0.5 * _cvis_to_chi2(s.mod_vis, data, w, vonmises)
    recompute_regs!(regs, s)
    s.lPrior = lprior(regs)

    # Initial temperature (squeeze.c:797-808)
    tcap = SQ_FLAT_CHI2_MULT * flatchi2 / ndf
    s.temperature = 2.0 * s.lLikelihood / ndf / SQ_INITIAL_SCALED_CHI2_MULT / chi2_temp
    s.temperature > tcap && (s.temperature = tcap)
    converged_early = false
    if s.temperature < SQ_INITIAL_TMIN
        s.temperature = SQ_INITIAL_TMIN
        converged_early = true
    end

    # Burn-in: C initialises this to niter and lowers it on convergence, which
    # can leave a non-converging chain with zero samples.  We instead start at
    # "the last 1/BURN_IN_FRAC of sweeps" so the mean image is always defined,
    # and let the convergence rule lower it further.
    burn_in = converged_early ? niter ÷ SQ_BURN_IN_FRAC :
                                niter - max(1, niter ÷ SQ_BURN_IN_FRAC)

    mean_image = zeros(Int64, nx, nx)     # Int64: a Float32 accumulator would
    nsamples = 0                          # silently stall past 2^24

    steptype = STEP_MEDIUM
    s.prob_movement = SQ_TARGET_MPROB
    s.naccepted = 0; s.nproposed = 0; s.since_resync = 0

    verb && print_every > 0 &&
        _print_squeeze_diagnostics(stdout, chain, 0, niter, s, regs, model,
                                   data, w, vonmises, burn_in)

    for sweep in 1:niter
        for _ in 1:nelements
            # ── Parametric-model move ────────────────────────────────────────
            if nparams_free > 0 && (target = rand(rng, 1:ntargets)) > nelements
                _param_move!(s, regs, model, data, w, rng;
                             vonmises = vonmises, ndf = ndf, target = target,
                             nelements = nelements)
                continue
            end

            i = rand(rng, 1:nelements)

            xstep, ystep = _draw_step(rng, steptype, s, i, nx)
            ix_old = Int(s.element_x[i]); iy_old = Int(s.element_y[i])
            p_old  = Int(s.element_p[i])
            ix_new = _wrap(ix_old + xstep, nx)
            iy_new = _wrap(iy_old + ystep, nx)
            p_new  = _pixel_index(ix_new, iy_new, nx)

            newvals, new_lPrior = _proposal_regs!(regs, s, ix_old, iy_old, ix_new, iy_new)

            rank1_update!(s.new_raw_im_vis, s.raw_im_vis, A, p_new, p_old)
            form_mod_vis!(s.new_mod_vis, s.new_raw_im_vis, s)
            new_lLikelihood = 0.5 * _cvis_to_chi2(s.new_mod_vis, data, w, vonmises)

            s.prob_movement *= (1.0 - 1.0 / SQ_DAMPING_TIME)
            s.nproposed += 1

            # Metropolis: temperature tempers the LIKELIHOOD ONLY, never the prior.
            tt = (new_lLikelihood - s.lLikelihood) / s.temperature + (new_lPrior - s.lPrior)

            if log(rand(rng)) < -tt
                s.lLikelihood = new_lLikelihood
                s.lPrior = new_lPrior
                s.raw_im_vis, s.new_raw_im_vis = s.new_raw_im_vis, s.raw_im_vis
                s.mod_vis, s.new_mod_vis = s.new_mod_vis, s.mod_vis
                s.histogram[ix_old, iy_old] -= Int32(1)
                s.histogram[ix_new, iy_new] += Int32(1)
                s.element_x[i] = ix_new; s.element_y[i] = iy_new; s.element_p[i] = p_new
                regs.v_l0, regs.v_tv, regs.v_entropy, regs.v_compactness,
                    regs.v_centering, regs.v_priorimage = newvals
                s.prob_movement += 1.0 / SQ_DAMPING_TIME
                s.naccepted += 1
                s.since_resync += 1

                # Temperature relaxes toward chi2r/chi2_temp — on ACCEPTED moves
                # only, so the effective timescale is TEMP_CHANGE_TIME/acceptance.
                # Do not "fix" this: updating every proposal overcools.
                chi2r = 2.0 * s.lLikelihood / ndf
                s.temperature *= 1.0 + (1.0 / s.temperature / SQ_TEMP_CHANGE_TIME) *
                                 (chi2r - chi2_temp * s.temperature) *
                                 (chi2r - chi2_target) * ndf / (2.0 * s.lLikelihood)
                if s.temperature < tmin || chi2r < chi2_target
                    burn_in = min(burn_in, sweep + niter ÷ SQ_BURN_IN_FRAC)
                end
                s.temperature = clamp(s.temperature, tmin, tcap)

                if s.since_resync >= resync_every
                    resync_raw_im_vis!(s, A)
                    form_mod_vis!(s.mod_vis, s.raw_im_vis, s)
                    s.lLikelihood = 0.5 * _cvis_to_chi2(s.mod_vis, data, w, vonmises)
                end
            else
                # Reverse the centroid mutation made by cent_change!
                regs.λ_centering > 0 &&
                    cent_change!(regs, ix_old, iy_old, ix_new, iy_new, nx)
            end

            # ── Adapt the step type (squeeze.c:1310-1360) ─────────────────────
            if 2.0 * s.lLikelihood < SQ_FLAT_CHI2_MULT * flatchi2
                if steptype == STEP_COPYCAT || steptype == STEP_ANYWHERE
                    steptype = STEP_SMALL
                elseif s.prob_movement > SQ_MPROB_HIGH
                    steptype = min(steptype + 1, STEP_MAX)
                elseif s.prob_movement < SQ_MPROB_LOW
                    steptype = max(steptype - 1, 0)
                end
            else
                steptype = STEP_MEDIUM
            end
            # Fixed-fraction overrides; ANYWHERE applied last so it wins.
            rand(rng) < f_copycat  && (steptype = STEP_COPYCAT)
            rand(rng) < f_anywhere && (steptype = STEP_ANYWHERE)
        end

        if sweep >= burn_in
            @inbounds for p in eachindex(mean_image)
                mean_image[p] += s.histogram[p]
            end
            nsamples += 1
        end
        if verb && print_every > 0 && (sweep % print_every == 0 || sweep == niter)
            _print_squeeze_diagnostics(stdout, chain, sweep, niter, s, regs, model,
                                       data, w, vonmises, burn_in)
        end
        if monitor !== nothing && sweep % monitor.every == 0
            monitor_update!(monitor, Float64.(s.histogram) ./ nelements, sweep;
                            chi2r = 2 * s.lLikelihood / ndf,
                            temperature = s.temperature,
                            params = model === nothing ? nothing : model.params)
        end
    end

    # Guard against maintained-state drift going unnoticed.
    resync_raw_im_vis!(s, A)
    form_mod_vis!(s.mod_vis, s.raw_im_vis, s)
    final_chi2 = _cvis_to_chi2(s.mod_vis, data, w, vonmises)

    img = nsamples > 0 ? Float64.(mean_image) ./ (nsamples * nelements) :
                         Float64.(s.histogram) ./ nelements

    # Regularizer values of the final histogram, recomputed from scratch: these
    # double as a check that the running deltas have not drifted.
    regs_final = (l0 = reg_l0(s.histogram),
                  tv = reg_tv(s.histogram, nelements, regs.tv_eps),
                  entropy = reg_entropy(s.histogram),
                  compactness = reg_compactness(s.histogram),
                  centering = regs.v_centering,
                  priorimage = regs.v_priorimage)

    return (image = img,
            params = model === nothing ? Float64[] : copy(model.params),
            stepsize = model === nothing ? Float64[] : copy(model.stepsize),
            last_image = Float64.(s.histogram) ./ nelements,
            chi2r_last = final_chi2 / ndf,
            temperature = s.temperature,
            acceptance = s.naccepted / max(s.nproposed, 1),
            nsamples = nsamples, burn_in = burn_in,
            regs_final = regs_final)
end

"""
    _param_move!(s, regs, model, data, w, rng; …)

One parametric-model proposal (squeeze.c:1161-1192).  Unlike an element move
this recomputes `param_vis`/`fluxratio_image` from scratch — there is no rank-1
shortcut — and rescales the whole image contribution:

    new_mod_vis = raw_im_vis .* new_fluxratio ./ nelements .+ new_param_vis

⚠️ The `params[j] ± stepsize[j]` two-point lattice is C's annealing proposal.  The
lattice spacing is rescaled by the adaptation on every proposal, so the chain can
never return exactly to a previous parameter value: the move is **not reversible**
and the resulting spread of `params` is not a posterior width.  That is acceptable
for an optimizer; under tempering the Gaussian variant is used instead
(`_param_step_tempered!`).
"""
function _param_move!(s::SqueezeState{T}, regs::SqueezeRegs, model::SqueezeSparco,
                      data::OIdata{T}, w::NTuple{7,Float64}, rng::AbstractRNG;
                      vonmises::Bool, ndf::Float64, target::Int,
                      nelements::Int) where {T}

    j = model.free[((target - nelements - 1) >> 1) + 1]   # PARAMS_PER_ELT = 2
    old = model.params[j]
    model.params[j] = old + model.stepsize[j] * (rand(rng, Bool) ? 1.0 : -1.0)

    new_v_modelparam = sparco_vis!(s.new_param_vis, s.new_fluxratio_image,
                                   model.params, data)
    form_mod_vis!(s.new_mod_vis, s.raw_im_vis,
                  s.new_param_vis, s.new_fluxratio_image, nelements)
    new_lLikelihood = 0.5 * _cvis_to_chi2(s.new_mod_vis, data, w, vonmises)

    new_lPrior = s.lPrior + regs.λ_modelparam * (new_v_modelparam - regs.v_modelparam)

    # C decays the EWMA and rescales the lattice BEFORE the accept test
    # (squeeze.c:1184-1187), so the step in force is the one adapted on the
    # previous proposal's outcome.  Kept in that order.
    model.prob_pmovement[j] *= (1.0 - 1.0 / SQ_PARAM_DAMPING_TIME)
    model.stepsize[j] *= 1.0 +
        (model.prob_pmovement[j] - SQ_TARGET_MPROB) / SQ_STEPSIZE_ADJUST_TIME
    s.nproposed += 1

    tt = (new_lLikelihood - s.lLikelihood) / s.temperature + (new_lPrior - s.lPrior)
    if log(rand(rng)) < -tt
        s.lLikelihood = new_lLikelihood
        s.lPrior = new_lPrior
        regs.v_modelparam = new_v_modelparam
        s.param_vis, s.new_param_vis = s.new_param_vis, s.param_vis
        s.fluxratio_image, s.new_fluxratio_image = s.new_fluxratio_image, s.fluxratio_image
        s.mod_vis, s.new_mod_vis = s.new_mod_vis, s.mod_vis
        model.prob_pmovement[j] += 1.0 / SQ_PARAM_DAMPING_TIME
        s.naccepted += 1
    else
        model.params[j] = old
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# Parallel tempering (Pigeons.jl).
#
# A deliberate *behavioural change* from the C code, not a port.  C runs an
# ad-hoc schedule that tracks chi2r (squeeze.c:684-965); here Pigeons owns the
# ladder, the non-reversible swaps and the schedule adaptation (via the global
# communication barrier Λ), and returns logZ.
#
# The path is
#
#     π_β(x)  ∝  exp( β · ( -lLikelihood(x) - lPrior(x) ) )
#
# i.e. the reference (β = 0) is the UNIFORM distribution over element
# configurations, which can be sampled i.i.d. exactly.  This tempers the
# likelihood *and* the prior together, whereas C's annealing tempers only the
# likelihood.  What matters is that β = 1 is the correct posterior; the path
# between is a design choice, and an exactly-i.i.d.-samplable reference is what
# makes PT's guarantees hold.
#
# The types live here so they exist without Pigeons loaded; all the Pigeons glue
# is in ext/OITOOLSPigeonsExt.jl.
# ─────────────────────────────────────────────────────────────────────────────

"""
    SqueezeTarget{T}

Pigeons target wrapping a SQUEEZE reconstruction problem.  Immutable and shared
across replicas; everything mutable lives in [`SqueezePTState`](@ref).
"""
struct SqueezeTarget{T<:AbstractFloat}
    data::OIdata{T}
    A::Matrix{Complex{T}}
    nx::Int
    nelements::Int
    weights::NTuple{7,Float64}
    vonmises::Bool
    regularizers::Vector{Any}
    prior_penalty::Vector{Float64}
    fov::Float64
    cent_mult::Float64
    tv_eps::Float64
    x_start::Any
    model::Union{Nothing,SqueezeSparco}   # template; each replica gets a copy
    param_lo::Vector{Float64}             # reference box for the free parameters
    param_hi::Vector{Float64}
end

"""
    SqueezePTState{T}

Per-replica state: the element configuration, its regularizer bookkeeping, a
private model copy, and running accumulators of the β = 1 samples this replica
has visited.

Pigeons swaps *chain indices* between replicas rather than moving states, so a
state stays with its replica for the whole run and its preallocated buffers are
never copied.  Every replica therefore spends some time at β = 1; summing
`mean_accum` over replicas gives exactly the set of target-chain samples.
"""
mutable struct SqueezePTState{T<:AbstractFloat}
    s::SqueezeState{T}
    regs::SqueezeRegs
    model::Union{Nothing,SqueezeSparco}   # thread-private copy, or nothing
    mean_accum::Matrix{Int64}     # Int64: a Float32 accumulator stalls past 2^24
    param_accum::Vector{Float64}  # running sum of SPARCO params over beta=1 visits
    nsamples::Int
    last_round::Int               # accumulators are reset when the round changes
end

"""
    SqueezeUniformReference

The β = 0 end of the path: the uniform distribution over element configurations
(times a uniform box over any free SPARCO parameters).  Its log density is
constant on its support, and it is sampled exactly.
"""
struct SqueezeUniformReference{T<:AbstractFloat}
    target::SqueezeTarget{T}
end

"""
    default_sparco_bounds(model) -> (lo, hi)

Box for the β = 0 reference over SPARCO parameters.

Parallel tempering needs a reference that can be sampled i.i.d., so the free
parameters need a *proper* (bounded) uniform reference — SPARCO's own box priors
leave `f_star` and `ud` unbounded above, which is improper.  These defaults are
deliberately generous; override with `param_bounds` when you know better, since a
box that is too wide wastes chains and one that excludes the truth biases the
answer.
"""
function default_sparco_bounds(model::SqueezeSparco)
    lo = [0.0,  0.0, -5.0, model.params[4], 0.0, -5.0]
    hi = [1.0, 20.0,  5.0, model.params[4], 1.0,  5.0]
    return lo, hi
end

@inline function _in_box(p::Vector{Float64}, lo::Vector{Float64}, hi::Vector{Float64},
                         free::Vector{Int})
    @inbounds for j in free
        (lo[j] <= p[j] <= hi[j]) || return false
    end
    return true
end

"""
    SqueezeExplorer

Element-move kernel at a fixed β, used as a custom Pigeons explorer.

The mixture is **fixed** (not adapted), which keeps the kernel a valid
β-invariant MH transition without any vanishing-adaptation argument.
`STEP_COPYCAT` is deliberately absent: it lands on a random donor element's
pixel, so `P(propose o→p) ∝ count(p)`, which is asymmetric and irreversible when
`count(o) == 1`.  That is an effective optimizer heuristic for annealing and an
invalid sampler move — see `reconstruct_squeeze`'s `f_copycat`.
"""
Base.@kwdef struct SqueezeExplorer
    n_passes::Int = 1          # sweeps of `nelements` proposals per Pigeons step
    p_small::Float64 = 0.7
    p_medium::Float64 = 0.25
    p_anywhere::Float64 = 0.05
    monitor::Union{Nothing,SqueezeMonitor} = nothing   # live display, or nothing
end

"""
    squeeze_logdensity(st) -> Float64

Un-tempered log density of a state: `-(lLikelihood + lPrior)`.  The maintained
values are kept current by the explorer, so this is O(1).
"""
squeeze_logdensity(st::SqueezePTState) = -(st.s.lLikelihood + st.s.lPrior)

"""
    squeeze_refresh!(st, target)

Recompute the maintained visibilities, likelihood and regularizers from the
current element configuration and model parameters.  Used after an i.i.d.
reference draw and to bound rank-1 drift.
"""
function squeeze_refresh!(st::SqueezePTState{T}, target::SqueezeTarget{T}) where {T}
    if st.model !== nothing
        st.regs.v_modelparam = sparco_vis!(st.s.param_vis, st.s.fluxratio_image,
                                           st.model.params, target.data)
        st.regs.λ_modelparam = 1.0
        copyto!(st.s.new_param_vis, st.s.param_vis)
        copyto!(st.s.new_fluxratio_image, st.s.fluxratio_image)
    end
    resync_raw_im_vis!(st.s, target.A)
    form_mod_vis!(st.s.mod_vis, st.s.raw_im_vis, st.s)
    st.s.lLikelihood = 0.5 * _cvis_to_chi2(st.s.mod_vis, target.data,
                                           target.weights, target.vonmises)
    recompute_regs!(st.regs, st.s)
    st.s.lPrior = lprior(st.regs)
    return st
end

"""
    squeeze_step!(st, target, explorer, beta, rng)

One sweep of `nelements` proposals at inverse temperature `beta`.

The Metropolis test is `log(u) < -beta * (Δ lLikelihood + Δ lPrior)`, the
acceptance for `π_β ∝ exp(β · logdensity)`.  At `beta = 1` this is the posterior;
at `beta = 0` every move is accepted, recovering the uniform reference.
"""
function squeeze_step!(st::SqueezePTState{T}, target::SqueezeTarget{T},
                       explorer::SqueezeExplorer, beta::Float64,
                       rng::AbstractRNG) where {T}
    s, regs = st.s, st.regs
    nx, nelements, A = s.nx, s.nelements, target.A
    p_sm, p_md = explorer.p_small, explorer.p_small + explorer.p_medium
    model = st.model
    nparams_free = model === nothing ? 0 : length(model.free)
    ntargets = nelements + 2 * nparams_free

    for _ in 1:(explorer.n_passes * nelements)
        # ── Parametric-model move ────────────────────────────────────────────
        if nparams_free > 0 && (tgt = rand(rng, 1:ntargets)) > nelements
            _param_step_tempered!(st, target, beta, tgt, nelements, rng)
            continue
        end

        i = rand(rng, 1:nelements)

        u = rand(rng)
        steptype = u < p_sm ? STEP_SMALL : (u < p_md ? STEP_MEDIUM : STEP_ANYWHERE)
        xstep, ystep = _draw_step(rng, steptype, s, i, nx)

        ix_old = Int(s.element_x[i]); iy_old = Int(s.element_y[i])
        p_old  = Int(s.element_p[i])
        ix_new = _wrap(ix_old + xstep, nx)
        iy_new = _wrap(iy_old + ystep, nx)
        p_new  = _pixel_index(ix_new, iy_new, nx)

        newvals, new_lPrior = _proposal_regs!(regs, s, ix_old, iy_old, ix_new, iy_new)
        rank1_update!(s.new_raw_im_vis, s.raw_im_vis, A, p_new, p_old)
        form_mod_vis!(s.new_mod_vis, s.new_raw_im_vis, s)
        new_lLikelihood = 0.5 * _cvis_to_chi2(s.new_mod_vis, target.data,
                                              target.weights, target.vonmises)

        s.nproposed += 1
        Δ = (new_lLikelihood - s.lLikelihood) + (new_lPrior - s.lPrior)
        if log(rand(rng)) < -beta * Δ
            s.lLikelihood = new_lLikelihood
            s.lPrior = new_lPrior
            s.raw_im_vis, s.new_raw_im_vis = s.new_raw_im_vis, s.raw_im_vis
            s.mod_vis, s.new_mod_vis = s.new_mod_vis, s.mod_vis
            s.histogram[ix_old, iy_old] -= Int32(1)
            s.histogram[ix_new, iy_new] += Int32(1)
            s.element_x[i] = ix_new; s.element_y[i] = iy_new; s.element_p[i] = p_new
            regs.v_l0, regs.v_tv, regs.v_entropy, regs.v_compactness,
                regs.v_centering, regs.v_priorimage = newvals
            s.naccepted += 1
        else
            regs.λ_centering > 0 && cent_change!(regs, ix_old, iy_old, ix_new, iy_new, nx)
        end
    end
    return st
end

"""
    _param_step_tempered!(st, target, beta, tgt, nelements, rng)

One SPARCO parameter proposal inside the tempered kernel.

Two differences from the annealing move in [`_param_move!`](@ref), both required
for the kernel to be a valid β-invariant MH transition — and both prescribed by
the C code, which switches proposals when tempering is enabled:

  * a **symmetric Gaussian** random walk `p + σ·randn()` rather than the
    two-point `p ± σ` lattice.  The lattice's spacing is rescaled on every
    proposal, so the chain can never return exactly to a previous value: the move
    is not reversible, and the spread of visited values is not a posterior width.
  * the stepsize is **frozen** — no adaptation, so no vanishing-adaptation
    argument is needed.

Proposals leaving the reference box are rejected outright, which is what makes
the β = 0 reference a proper uniform-on-box distribution (SPARCO's own priors
leave `f_star` and `ud` unbounded above, which cannot be sampled i.i.d.).
"""
function _param_step_tempered!(st::SqueezePTState{T}, target::SqueezeTarget{T},
                               beta::Float64, tgt::Int, nelements::Int,
                               rng::AbstractRNG) where {T}
    model = st.model
    s, regs = st.s, st.regs
    j = model.free[((tgt - nelements - 1) >> 1) + 1]      # PARAMS_PER_ELT = 2
    old = model.params[j]

    model.params[j] = old + model.stepsize[j] * randn(rng)
    s.nproposed += 1

    # Outside the reference box the proposal has zero density at every β.
    if !_in_box(model.params, target.param_lo, target.param_hi, model.free)
        model.params[j] = old
        return st
    end

    new_v_modelparam = sparco_vis!(s.new_param_vis, s.new_fluxratio_image,
                                   model.params, target.data)
    form_mod_vis!(s.new_mod_vis, s.raw_im_vis,
                  s.new_param_vis, s.new_fluxratio_image, nelements)
    new_lLikelihood = 0.5 * _cvis_to_chi2(s.new_mod_vis, target.data,
                                          target.weights, target.vonmises)
    new_lPrior = s.lPrior + regs.λ_modelparam * (new_v_modelparam - regs.v_modelparam)

    Δ = (new_lLikelihood - s.lLikelihood) + (new_lPrior - s.lPrior)
    if log(rand(rng)) < -beta * Δ
        s.lLikelihood = new_lLikelihood
        s.lPrior = new_lPrior
        regs.v_modelparam = new_v_modelparam
        s.param_vis, s.new_param_vis = s.new_param_vis, s.param_vis
        s.fluxratio_image, s.new_fluxratio_image = s.new_fluxratio_image, s.fluxratio_image
        s.mod_vis, s.new_mod_vis = s.new_mod_vis, s.mod_vis
        s.naccepted += 1
    else
        model.params[j] = old
    end
    return st
end

"""
    reconstruct_squeeze_tempered(data, ft; kwargs...)

Parallel-tempered SQUEEZE reconstruction, built on Pigeons.jl.

Requires Pigeons to be loaded:

```julia
using OITOOLS, Pigeons
img, diag = reconstruct_squeeze_tempered(data, ft; n_rounds = 10, n_chains = 10)
```

Unlike [`reconstruct_squeeze`](@ref), which anneals (an optimizer that floors at
T = 1), this samples the posterior at β = 1 with a full temperature ladder, and
returns a log normalising constant `logZ` alongside the mean image.
"""
function reconstruct_squeeze_tempered end

"""
    reconstruct_squeeze(data, ft; kwargs...)           -> (image, diagnostics)
    reconstruct_squeeze(x_start, data, ft; kwargs...)  -> (image, diagnostics)

SQUEEZE-style MCMC image reconstruction.

The image is represented as a bag of `nelements` discrete flux quanta on the
pixel grid, sampled by Metropolis-Hastings with simulated annealing.  This gives
positivity and fixed total flux for free, admits non-convex regularizers such as
L0, and returns a posterior mean image rather than a MAP point.

`x_start` selects the starting configuration; omit it for SQUEEZE's own default.

| value | meaning |
|---|---|
| `:point_source` (default) | all quanta on the centre pixel — a Dirac, as when C is run without `-i` |
| `:random` | quanta scattered uniformly (an OITOOLS addition; gives multiple chains genuinely different starts) |
| an `nx × nx` image | digitised into quanta — C's `-i` |

Strings work too (`"point_source"`, `"random"`).  Note the default is a point
source, **not** a random or flat image: the annealing is expected to spread the
flux out from it.

`ft` follows the same convention as `reconstruct` / `reconstruct_bsmem`: the image
geometry is taken *from the Fourier operator*, not from keywords.

- Pass a **DFT matrix** (`setup_ft(data, nx, pixsize; mode="dft")` or `setup_dft`)
  and it is used directly — this is what the sampler wants, and nothing is rebuilt.
- Pass an **NFFT plan** and `nx`/`pixsize` are recovered from it exactly as
  `reconstruct_bsmem` does, then a DFT matrix is built internally.  The sampler
  cannot use an NFFT plan: it needs O(nuv) access to a *single column* per move
  (the rank-1 update), which a transform-based operator cannot provide.

# Keyword arguments
- `nelements`: number of flux quanta (default [`default_nelements`](@ref)).
- `niter`: sweeps; each sweep is `nelements` proposals.
- `nchains`: independent annealing restarts, run on separate threads.
- `regularizers`: OITOOLS-style `[["l0", λ], ["tv", λ], ...]`.  Supported names
  are `l0`, `tv`, `entropy`, `compactness`, `centering`, `priorimage`.  **These
  are SQUEEZE's forms on the integer histogram, not OITOOLS' differentiable
  versions.**
- `weights`: observable weights, default `[1.0,1.0,1.0]` = V², T3amp, T3phi, as in
  `reconstruct` (padded to 7 internally).
- `prior_image`: an `nx × nx` map of per-pixel prior *probabilities* (C's `-p`).
  It becomes an additive penalty `-log(p)` per quantum, with `p <= 0` mapped to
  `1e12` — which is what makes it a **mask**: no quantum can ever sit on a
  zero-prior pixel.  Its weight defaults to 1 and can be set explicitly with
  `["priorimage", λ]` in `regularizers`.  Supplying one disables auto-centering,
  since the prior already fixes the position (C warns about this combination).
- `pixsize`: only needed to override the value recovered from `ft`, or to supply
  one when `ft` is a bare DFT matrix and `outfile` is set (the FITS WCS header).
- `tmin`, `chi2_temp`, `chi2_target`: annealing schedule (C's `-tm`, `-ct`, `-fc`).
- `f_anywhere`, `f_copycat`: fixed proposal fractions.
- `model`: an optional [`SqueezeSparco`](@ref) chromatic star + background model.
- `print_every`: C-format text diagnostic line every n sweeps.
- `monitor`: if > 0, redraw a live image and trace display every `monitor`
  sweeps (chain 1 only).  Off by default, and free when off — one integer
  comparison per sweep.  Enabling it pins chain 1 to the calling thread, because
  drawing must happen on the main thread; the other chains still run in parallel.
  `monitor_colormap` selects the image colour map.
- `seed`: per-chain streams are a pure function of `(seed, chain)`, so results
  are reproducible regardless of thread scheduling.

# Notes
`COPYCAT` proposals are asymmetric and irreversible; they are an effective
optimizer heuristic for annealing but are not a valid sampler move.  Set
`f_copycat=0` if you want the post-burn-in samples to be a defensible posterior.
"""
function reconstruct_squeeze(x_start, data::OIdata{T}, ft;
        pixsize::Real = -1.0,
        nelements::Integer = 0,
        niter::Integer = 200,
        nchains::Integer = 1,
        regularizers = [],
        weights = [1.0, 1.0, 1.0],
        vonmises::Bool = false,
        tmin::Real = 1.0,
        chi2_temp::Real = 4.0,
        chi2_target::Real = 0.0,
        f_anywhere::Real = 0.05,
        f_copycat::Real = 0.1,
        fov::Real = 1.0,
        cent_mult::Real = 1.0,
        tv_eps::Real = 0.0,
        resync_every::Integer = default_resync(T),
        auto_centering::Bool = true,
        model::Union{Nothing,SqueezeSparco} = nothing,
        prior_image = nothing,
        outfile::AbstractString = "",
        print_every::Integer = 0,
        monitor::Integer = 0,
        monitor_colormap::AbstractString = "gist_heat",
        seed::Integer = 12345,
        verb::Bool = true) where {T<:AbstractFloat}

    A, nx, pixsize = _squeeze_operator(ft, data, pixsize)

    if x_start isa AbstractMatrix && size(x_start, 1) != nx
        throw(DimensionMismatch("x_start is $(size(x_start,1))×$(size(x_start,2)) but " *
                                "`ft` describes an $(nx)×$(nx) image"))
    end

    # Prior image (C's -p): a per-pixel prior probability map, turned into an
    # additive penalty -log(p).  Pixels with p <= 0 get 1e12, which is what makes
    # this a *mask*: an element can never sit there.
    prior_penalty = Float64[]
    if prior_image !== nothing
        size(prior_image) == (nx, nx) || throw(DimensionMismatch(
            "prior_image is $(size(prior_image)) but `ft` describes an $(nx)x$(nx) image"))
        prior_penalty = [p > 0 ? -log(Float64(p)) : 1e12 for p in vec(prior_image)]
    end

    nelem = nelements > 0 ? Int(nelements) : default_nelements(data, nx)
    w = _weights_tuple(weights)

    # C enables centering by default when there is no parametric model to fit
    # (squeeze.c:567-570, DEFAULT_CENT_MULT=1).  Without it nothing pins the image
    # against the translation degeneracy of the bispectrum.  Set
    # `auto_centering=false` to opt out.  A parametric model or a prior image
    # already pins the position, so C turns it off in those cases too.
    regs_have_centering = any(r -> lowercase(String(r[1])) == "centering", regularizers)
    auto_centering = auto_centering && model === nothing && prior_image === nothing
    local effective_regs = (auto_centering && !regs_have_centering) ?
                     vcat(collect(Any, regularizers), Any[Any["centering", 1.0]]) :
                     collect(Any, regularizers)
    # C defaults reg_param[REG_PRIORIMAGE] to 1 when -p is given without a weight
    # (squeeze.c:552-553).
    if prior_image !== nothing &&
       !any(r -> lowercase(String(r[1])) == "priorimage", effective_regs)
        effective_regs = vcat(collect(Any, effective_regs), Any[Any["priorimage", 1.0]])
    end
    regs_of = () -> _parse_squeeze_regularizers(effective_regs; fov, cent_mult, tv_eps,
                                                prior_penalty = prior_penalty)

    centering_on = regs_of().λ_centering > 0
    nparams_free = model === nothing ? 0 : length(model.free)
    ndf = Float64(squeeze_ndof(data, weights; nparams_free = nparams_free,
                               centering = centering_on))
    fchi2 = flat_chi2(data, w)

    if verb
        @printf("SQUEEZE: nx=%d  pixsize=%.3f mas  nuv=%d  nelements=%d  ndf=%.0f\n",
                nx, pixsize, data.nuv, nelem, ndf)
        @printf("         T=%s  A=%.1f MB  niter=%d  nchains=%d  chi2r(flat)=%.2f\n",
                T, sizeof(A) / 1e6, niter, nchains, fchi2 / ndf)
    end

    results = Vector{Any}(undef, nchains)
    nw = min(Threads.nthreads(), nchains)

    # Only chain 1 ever builds a monitor (see `run_one` below), so drawing is never
    # concurrent. What it must not be is *off the main thread*: a Python call from a
    # `Threads.@spawn`ed task segfaults inside PythonCall's `PyTuple_New`, with or without
    # an interactive backend. This is a property of calling CPython from a non-main thread,
    # not of any one bridge — PythonCall behaves exactly as PyCall did here.
    #
    # So a monitored run keeps chain 1 on the calling thread and spawns the others, instead
    # of dropping every chain to serial as it used to.
    monitored_inline = monitor > 0 && nw > 1
    if monitored_inline && Threads.threadid() != 1
        @info "reconstruct_squeeze: monitoring is enabled but this call is not on the " *
              "main thread, so the live display cannot be pinned there; running chains " *
              "serially. Set monitor=0 for the fully threaded path."
        monitored_inline = false
        nw = 1
    end

    function run_one(c)
        s = SqueezeState(T, nx, nelem, data.nuv)
        rng = Xoshiro(hash((seed, c)))
        init_elements!(s, rng; x_start = x_start)
        # Each chain mutates params/stepsize/prob_pmovement, so give it its own copy.
        m = model === nothing ? nothing :
            SqueezeSparco(copy(model.params), copy(model.stepsize),
                          copy(model.free), copy(model.prob_pmovement))
        mon = (monitor > 0 && c == 1) ?
              SqueezeMonitor(monitor; pixsize = pixsize,
                             title = @sprintf("SQUEEZE chain %d", c),
                             colormap = monitor_colormap,
                             free = m === nothing ? Int[] : copy(m.free)) : nothing
        mon === nothing || (mon.ndf = ndf)
        results[c] = run_chain!(s, regs_of(), A, data, w, rng;
                                niter = Int(niter), ndf = ndf,
                                chi2_temp = Float64(chi2_temp),
                                chi2_target = Float64(chi2_target),
                                tmin = Float64(tmin),
                                f_anywhere = Float64(f_anywhere),
                                f_copycat = Float64(f_copycat),
                                resync_every = Int64(resync_every),
                                vonmises = vonmises, flatchi2 = fchi2, model = m,
                                verb = verb, chain = c, monitor = mon,
                                print_every = Int(print_every) > 0 ? Int(print_every) :
                                              (nchains == 1 ? max(1, Int(niter) ÷ 10) : 0))
    end

    if nw <= 1
        for c in 1:nchains; run_one(c); end
    elseif monitored_inline
        # Chain 1 owns the display, so it runs here on the main thread while the rest are
        # spawned. That recovers nchains-1 chains of parallelism versus going serial.
        rest   = collect(2:nchains)
        nwk    = min(nw, length(rest))
        chunks = [rest[k:nwk:end] for k in 1:nwk]
        tasks  = [Threads.@spawn (for c in chunks[k]; run_one(c); end) for k in 1:nwk]
        run_one(1)
        foreach(wait, tasks)
    else
        chunks = [c:nw:nchains for c in 1:nw]
        tasks = [Threads.@spawn (for c in chunks[k]; run_one(c); end) for k in 1:nw]
        foreach(wait, tasks)
    end

    # Rank chains by the chi2r of their MEAN image (not the last sample).
    chi2r_mean = Vector{Float64}(undef, nchains)
    pv = zeros(Complex{T}, data.nuv); fr = ones(T, data.nuv)
    for c in 1:nchains
        V = image_to_vis(T.(results[c].image), A)     # normalised to unit flux
        if model !== nothing
            # Re-apply this chain's fitted model, or the mean image is scored
            # against the data without the star/background it was fitted with.
            sparco_vis!(pv, fr, results[c].params, data)
            @inbounds @simd for k in eachindex(V)
                V[k] = pv[k] + V[k] * fr[k]
            end
        end
        chi2r_mean[c] = _cvis_to_chi2(V, data, w, vonmises) / ndf
    end
    best = argmin(chi2r_mean)

    if verb
        for c in 1:nchains
            @printf("  chain %2d: chi2r(mean)=%9.4f  chi2r(last)=%9.4f  acc=%.3f  nsamples=%d%s\n",
                    c, chi2r_mean[c], results[c].chi2r_last, results[c].acceptance,
                    results[c].nsamples, c == best ? "  <-- best" : "")
        end
    end

    diagnostics = (chi2r_mean = chi2r_mean,
                   chi2r_last = [r.chi2r_last for r in results],
                   acceptance = [r.acceptance for r in results],
                   temperature = [r.temperature for r in results],
                   nsamples = [r.nsamples for r in results],
                   burn_in = [r.burn_in for r in results],
                   best_chain = best,
                   images = [r.image for r in results],
                   regs_final = [r.regs_final for r in results],
                   params = [r.params for r in results],
                   stepsize = [r.stepsize for r in results],
                   param_names = SPARCO_PARAM_NAMES,
                   ndf = ndf, nelements = nelem, flat_chi2r = fchi2 / ndf)

    # `pixsize` gives writefits the CDELT/CRPIX WCS header; without it the image
    # is written headerless (writefits' own `pixsize=-1` sentinel).
    isempty(outfile) || writefits(results[best].image, outfile; pixsize = pixsize)

    return results[best].image, diagnostics
end

"""
    reconstruct_squeeze(data, ft; kwargs...)

Start from SQUEEZE's default configuration — a point source at the centre of the
grid, the equivalent of running the C code without `-i`.  Identical to passing
`:point_source` as `x_start`.
"""
reconstruct_squeeze(data::OIdata, ft; kwargs...) =
    reconstruct_squeeze(:point_source, data, ft; kwargs...)

reconstruct_squeeze(data::AbstractArray{<:OIdata}, ft; kwargs...) =
    reconstruct_squeeze(:point_source, data, ft; kwargs...)

"""
    reconstruct_squeeze(x_start, data::AbstractArray{<:OIdata}, ft; kwargs...)

Convenience dispatch matching the rest of the `reconstruct*` family, so the usual
`data = readoifits(file); ft = setup_ft(data, nx, pixsize; mode="dft")` pattern
works without indexing.  Only the monochromatic, single-epoch case is supported.
"""
function reconstruct_squeeze(x_start, data::AbstractArray{<:OIdata}, ft; kwargs...)
    length(data) == 1 || throw(ArgumentError(
        "reconstruct_squeeze: polychromatic / multi-epoch reconstruction is not " *
        "implemented (got a $(size(data)) OIdata array). Pass a single OIdata, " *
        "e.g. data[1,1], to reconstruct one channel."))
    return reconstruct_squeeze(x_start, first(data), ft; kwargs...)
end

function _parse_squeeze_regularizers(regularizers; fov = 1.0, cent_mult = 1.0, tv_eps = 0.0,
                                     prior_penalty = Float64[])
    r = SqueezeRegs(; fov = Float64(fov), cent_mult = Float64(cent_mult),
                      tv_eps = Float64(tv_eps), prior_penalty = prior_penalty)
    for spec in regularizers
        length(spec) == 2 ||
            throw(ArgumentError("each regularizer must be [name, λ], got $spec"))
        name = lowercase(String(spec[1])); λ = Float64(spec[2])
        if     name == "l0";           r.λ_l0 = λ
        elseif name == "tv";           r.λ_tv = λ
        elseif name == "entropy";      r.λ_entropy = λ
        elseif name == "compactness";  r.λ_compactness = λ
        elseif name == "centering";    r.λ_centering = λ
        elseif name == "priorimage";   r.λ_priorimage = λ
        else
            throw(ArgumentError("unknown SQUEEZE regularizer \"$name\"; supported: " *
                                "l0, tv, entropy, compactness, centering, priorimage"))
        end
    end
    return r
end

"""
    _squeeze_operator(ft, data, pixsize) -> (A, nx, pixsize)

Resolve the image geometry and the dense DFT matrix from whatever Fourier
operator the caller passed, following `reconstruct_bsmem`'s convention of taking
`nx`/`pixsize` from `ft` rather than from keywords.

A DFT matrix is used as-is.  An NFFT plan only yields the geometry — the sampler
needs single-column access to the transform (`rank1_update!`), which no
transform-based operator provides, so a DFT matrix is built and the memory cost
is reported.
"""
function _squeeze_operator(ft, data::OIdata{T}, pixsize::Real) where {T}
    # Unwrap the 1×1 containers setup_ft returns for the monochromatic case.
    while (ft isa AbstractArray) && !(ft isa AbstractMatrix{<:Complex}) &&
          !(ft isa AbstractVector{<:NFFT.NFFTPlan}) && length(ft) == 1
        ft = first(ft)
    end
    # setup_nfft returns SIX plans: [all-uv, vis, v2, t3_1, t3_2, t3_3].  The
    # first is the full-uv one, which is the only one whose geometry we need —
    # this is the same `ft[1]` convention reconstruct_bsmem uses (oimem.jl).
    ft isa AbstractVector{<:NFFT.NFFTPlan} && (ft = first(ft))

    if ft isa AbstractMatrix{<:Complex}
        size(ft, 1) == data.nuv ||
            throw(DimensionMismatch("DFT matrix has $(size(ft,1)) rows but data.nuv = $(data.nuv)"))
        nx = isqrt(size(ft, 2))
        nx * nx == size(ft, 2) ||
            throw(ArgumentError("DFT matrix has $(size(ft,2)) columns, not a perfect square"))
        eltype(ft) === Complex{T} ||
            throw(ArgumentError("DFT matrix is $(eltype(ft)) but data is OIdata{$T}; " *
                                "build it from the same data (see `setup_dft`)"))
        return ft, nx, Float64(pixsize)

    elseif ft isa NFFT.NFFTPlan
        # Same recovery as reconstruct_bsmem (oimem.jl): setup_nfft stores
        # ft.k[:, j] = pixsize_rad * [u_j; v_j].
        nx = ft.N[1]
        j = argmax(vec(sqrt.(sum(data.uv .^ 2, dims = 1))))
        ps = pixsize > 0 ? Float64(pixsize) :
             (norm(ft.k[:, j]) / norm(data.uv[:, j])) / (pi / 180.0 / 3600000.0)
        A = setup_dft(data, nx, ps)
        @info @sprintf("reconstruct_squeeze: NFFT plan given; building the DFT matrix the \
                        sampler needs (nx=%d, pixsize=%.4f mas, %.1f MB)", nx, ps, sizeof(A)/1e6)
        return A, nx, ps
    end

    throw(ArgumentError("reconstruct_squeeze: `ft` must be a DFT matrix " *
                        "(setup_ft(...; mode=\"dft\")) or an NFFT plan, got $(typeof(ft))"))
end

"""
    _print_squeeze_diagnostics(io, chain, iter, niter, s, regs, model, data, w, vonmises, burn_in)

Per-sweep diagnostic line, matching C's `print_diagnostics` (squeeze.c:3234) field
for field.  `lPost`/`lPrior`/`lLike` are the un-reduced log-posterior terms; the
per-observable figures are **reduced** χ² (divided by their own counts, not by
ndf), as in C.  Regularizer entries show `λ·value`.  Colours follow OITOOLS' own
convention (V2 red, T3A blue, T3P green), which happens to be C's too.
"""
function _print_squeeze_diagnostics(io::IO, chain::Int, iter::Int, niter::Int,
                                    s::SqueezeState{T}, regs::SqueezeRegs,
                                    model::Union{Nothing,SqueezeSparco},
                                    data::OIdata{T}, w::NTuple{7,Float64},
                                    vonmises::Bool, burn_in::Int) where {T}
    lPost = s.lLikelihood + s.lPrior
    @printf(io, "Chain: %d lPost:%8.1f lPrior:%8.1f lLike:%9.1f ",
            chain, lPost, s.lPrior, s.lLikelihood)

    # Per-observable reduced chi2.  This is the allocating `_chi2_terms` path, but
    # it runs once per `print_every` sweeps, not once per proposal.
    t = _chi2_terms(s.mod_vis, data, collect(w), vonmises)
    w[1] > 0 && data.nv2    > 0 && printstyled(io, @sprintf("V2:%5.2f ",  t.chi2_v2    / data.nv2),    color=:red)
    w[2] > 0 && data.nt3amp > 0 && printstyled(io, @sprintf("T3A:%5.2f ", t.chi2_t3amp / data.nt3amp), color=:blue)
    w[3] > 0 && data.nt3phi > 0 && printstyled(io, @sprintf("T3P:%5.2f ", t.chi2_t3phi / data.nt3phi), color=:green)
    w[4] > 0 && data.nvisamp > 0 && printstyled(io, @sprintf("VA:%5.2f ",  t.chi2_visamp / data.nvisamp), color=:cyan)
    w[5] > 0 && data.nvisphi > 0 && printstyled(io, @sprintf("VP:%5.2f ",  t.chi2_visphi / data.nvisphi), color=:magenta)

    regs.λ_l0          > 0 && @printf(io, "l0:%5.2f ",   regs.λ_l0 * regs.v_l0)
    regs.λ_tv          > 0 && @printf(io, "tv:%5.2f ",   regs.λ_tv * regs.v_tv)
    regs.λ_entropy     > 0 && @printf(io, "en:%5.2f ",   regs.λ_entropy * regs.v_entropy)
    regs.λ_compactness > 0 && @printf(io, "cp:%5.2f ",   regs.λ_compactness * regs.v_compactness)
    regs.λ_priorimage  > 0 && @printf(io, "prior:%7.2f ", regs.λ_priorimage * regs.v_priorimage)
    if regs.λ_centering > 0
        @printf(io, "cent:%5.2f ", regs.λ_centering * regs.v_centering)
        @printf(io, "XY:(%5.2f,%5.2f) ", regs.centroid_x / s.nelements,
                                          regs.centroid_y / s.nelements)
    end

    if iter > 0
        @printf(io, "E: %5d MPr: %4.2f T: %5.2f B:%d Iter: %4d of %4d\n",
                s.nelements, s.prob_movement, s.temperature, burn_in, iter, niter)
    else
        @printf(io, "E: %5d MPr: %4.2f T: %5.2f -- INITIAL\n",
                s.nelements, s.prob_movement, s.temperature)
    end

    if model !== nothing
        @printf(io, "Chain: %d Model Parameters:\n", chain)
        for j in 1:NPARAMS_SPARCO
            # `MPr` per parameter is the acceptance EWMA that drives the stepsize.
            # If it sits near 0 the parameter is stuck: the image has adapted to the
            # current value and no step of any usable size is acceptable.  Without
            # this figure that failure is invisible.
            if j in model.free
                @printf(io, "%d - %-9s %10.5g +/- %9.5g  MPr: %4.2f\n", j - 1,
                        SPARCO_PARAM_NAMES[j], model.params[j], model.stepsize[j],
                        model.prob_pmovement[j])
            else
                @printf(io, "%d - %-9s %10.5g +/- %9.5g   (fixed)\n", j - 1,
                        SPARCO_PARAM_NAMES[j], model.params[j], model.stepsize[j])
            end
        end
        println(io)
    end
    return nothing
end
