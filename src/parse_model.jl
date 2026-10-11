# parse_model.jl
#
# Compile a flat parameter dict into a FlatModel struct that can be evaluated
# and differentiated without any Dict allocation in the hot path.
#
# Load order (within OITOOLS module):
#   include("resolvers.jl")          # SharedUtils, RGF, HandRolled, Symbolic
#   include("hankel.jl")             # hankel_vis, compile_profile, trapz_weights
#                                    #   (includes hankel_chainrules.jl internally)
#   include("model_chainrules.jl")   # vis_ud, vis_ldlin, vis_ldquad, vis_ldpow
#   include("parse_model.jl")        # this file
#
# Concepts
# --------
# A *model* is a flat Dict{String,Any} whose keys follow the convention
# "component_name,parameter_name".  Example:
#
#   Dict(
#     "star,ud"       => 0.8,               # UD diameter in mas
#     "star,f"        => 0.6,               # flux fraction
#     "ring,profile"  => "exp(-(\$R/\$scale)^2/2)",  # Hankel profile
#     "ring,scale"    => 1.5,               # profile param (mas)
#     "ring,udout"    => 10.0,              # r grid outer edge = udout/2
#     "ring,f"        => "1 - \$star,f",   # derived flux fraction
#   )
#
# Component type is inferred from which geometry key is present:
#   "ud"      -> analytic uniform disk        (vis_ud from model_chainrules.jl)
#   "fwhm"    -> analytic Gaussian            (vis_gaussian below)
#   "ldlin"   -> linear limb-darkened disk    (vis_ldlin)
#   "ldquad"  -> quadratic limb-darkened disk (vis_ldquad)
#   "ldpow"   -> power-law limb-darkened disk (vis_ldpow)
#   "profile" -> arbitrary radial profile via Hankel transform
#   (none of the above, only "f") -> unresolved point source
#
# All scalar params (numeric or derived-expression) go through the resolver.
# The "profile" string is compiled separately via compile_profile.
#
# Units
# -----
# uv      : 2×N matrix of spatial frequencies in cycles/rad  (OITOOLS convention)
# ρ       : sqrt(u²+v²) in cycles/rad,  passed to analytic vis functions
# B_mas   : ρ / MAS2RAD, passed to Hankel functions (cycles/mas)
# diameter: mas (angular)
# r grid  : mas

using LinearAlgebra: dot, mul!
using ForwardDiff

const _MAS2RAD_PM = 2.0626480624709636e8    # 1 radian in milli-arcseconds


# ─────────────────────────────────────────────────────────────────────────────
# Geometry key → component kind
# ─────────────────────────────────────────────────────────────────────────────

const _GEOMETRY_KEYS = ["image", "visfunc", "profile", "diamin", "fwhmin", "crin",
                       "ud", "fwhm", "ldlin", "ldquad", "ldsqrt", "ldclaret4",
                       "ldpow", "resolved"]

const _KIND_GEOMETRY = Dict{Symbol,String}(
    :image         => "image",        # a FITS image or cube (image_component.jl)
    :visfunc       => "visfunc",      # V(B) written out, in $BL (Mlambda)
    :hankel        => "profile",
    :ring          => "diamin",       # uniform ring (analytic)
    :gaussian_ring => "fwhmin",       # bi-Gaussian ring (analytic)
    :crescent      => "crin",         # crescent (analytic, two offset UDs)
    :ud            => "ud",
    :gaussian      => "fwhm",
    :ldlin         => "ldlin",
    :ldquad        => "ldquad",
    :ldsqrt        => "ldsqrt",       # square-root law
    :ldclaret4     => "ldclaret4",    # Claret four-parameter law
    :ldpow         => "ldpow",
    :resolved      => "resolved",     # fully resolved (V=1 at B=0, 0 at B>0), e.g. halo
    :point         => "",             # no geometry key — unresolved point
)

# Geometry params consumed by each analytic kind (in vis function order).
# These are the *suffix* parts of the keys (prefix is the component name).
const _ANALYTIC_PARAM_SUFFIXES = Dict{Symbol,Vector{String}}(
    :ring          => ["diamin", "diamout"],                   # inner, outer diameter
    :gaussian_ring => ["fwhmin", "fwhmout"],                   # inner, outer FWHM
    :crescent      => ["crin", "crout", "croff", "crprojang"], # diameters, offset, PA
    :ud            => ["ud"],
    :gaussian      => ["fwhm"],
    :ldlin         => ["ldlin", "u"],       # diameter, u
    :ldquad        => ["ldquad", "u", "w"], # diameter, u, w
    # Square root: I = 1 - c(1-mu) - d(1-sqrt(mu)). Named u/w for the same reason ldquad does --
    # u is the linear coefficient in every law that has one.
    :ldsqrt        => ["ldsqrt", "u", "w"],
    :ldclaret4     => ["ldclaret4", "c1", "c2", "c3", "c4"],
    :ldpow         => ["ldpow", "alpha"],   # diameter, alpha
    :resolved      => [],                   # fully resolved: no geometry params
    :point         => [],
)


# Default bounds by parameter suffix, for the dict interface.
#
# Keyed by SUFFIX, which is what lets it answer "what are sensible bounds for the key
# `star,ud`?" -- the question the removed `init_bounds` could not, being keyed by function
# object. The values are deliberately generous starting points, not physics: their
# job is to give a fit somewhere sane to begin and to satisfy samplers that require finite
# bounds (`fit_model_ultranest` errors without them), not to encode prior knowledge.
"""
Default upper bound on an angular size, in mas, when no data is supplied.

Long-baseline interferometry senses scales between roughly λ/2B_max and λ/B_min. Measured on
three real datasets, that is 0.5–1.2 mas at the fine end and 18–41 mas at the coarse end: a
source much larger than λ/B_min is fully resolved out and its size is not constrained at all.
100 mas therefore leaves generous headroom over anything an array can measure, while staying a
real bound. Pass `data` to `default_bounds` to size this from the actual uv coverage instead.
"""
const DEFAULT_MAX_SIZE_MAS = 100.0

const _SIZE_KEYS = ("ud", "fwhm", "ldlin", "ldquad", "ldsqrt", "ldclaret4", "ldpow",
                    "diamin", "diamout", "diam",
                    "fwhmin", "fwhmout", "crin", "crout", "udout", "r_max", "spatial_kernel")
const _OFFSET_KEYS = ("x", "y", "croff")

const _DEFAULT_PARAM_BOUNDS = Dict{String,Tuple{Float64,Float64}}(
    # angular sizes, mas
    (k => (0.0, DEFAULT_MAX_SIZE_MAS) for k in _SIZE_KEYS)...,
    # positions, mas: a companion can sit outside the largest scale the array resolves, so
    # these are looser than the sizes
    (k => (-DEFAULT_MAX_SIZE_MAS, DEFAULT_MAX_SIZE_MAS) for k in _OFFSET_KEYS)...,
    # limb darkening
    "u"       => (-1.0, 1.0), "w"       => (-1.0, 1.0), "alpha"   => (0.0, 3.0),
    # Claret's four coefficients are individually unbounded in sign and routinely exceed 1 --
    # they trade against each other, and only the profile they sum to has to stay positive.
    "c1" => (-3.0, 3.0), "c2" => (-3.0, 3.0), "c3" => (-3.0, 3.0), "c4" => (-3.0, 3.0),
    # orientation
    "pa"      => (-180.0, 180.0), "projang" => (-180.0, 180.0),
    "crprojang" => (-180.0, 180.0),
    "incl"    => (0.0, 90.0),
    # A flux is a relative weight -- the model divides by the total -- so it has no upper
    # bound; a fraction of a shape does.
    "f"       => (0.0, Inf),  "thick"   => (0.0, 1.0),
    # The crescent's offset is a FRACTION of the largest offset that keeps the subtracted disk
    # inside the outer one: `off = croff*(crout - crin)/2`, so 0 is a concentric ring and 1 is
    # internally tangent. Outside [0, 1] the hole leaves the disk and the model stops being a
    # crescent, which is a bound rather than a preference.
    "croff"   => (0.0, 1.0),
)

"""
    max_angular_scale(data) -> Float64

Largest angular scale the uv coverage is sensitive to, in mas: λ/B_min, from the shortest
baseline present. A source larger than this is fully resolved and its size is unconstrained,
which makes it the natural ceiling for a size bound.
"""
function max_angular_scale(data)
    B = Float64[]
    for d in (data isa AbstractArray ? vec(data) : [data])
        u = Float64.(d.uv[1, :]); v = Float64.(d.uv[2, :])
        append!(B, sqrt.(u .^ 2 .+ v .^ 2))
    end
    filter!(b -> b > 0 && isfinite(b), B)
    isempty(B) && return DEFAULT_MAX_SIZE_MAS
    return (1 / minimum(B)) / MAS      # uv is in cycles/rad
end

"""
    default_bounds(model_dict, list_free_params; data = nothing, max_size = nothing) -> (lb, ub)

Suggest lower and upper bounds for each free parameter, as two `Dict{String,Float64}` ready
to pass to any fitter's `lb`/`ub` keyword.

Bounds are chosen from the parameter's suffix — the part after the first comma, so
`"star,ud"` is looked up as `"ud"`. Azimuthal-mode keys (`"c,az amp1"`, `"c,az projang1"`)
are recognised by prefix. Anything unrecognised, including bare global names, gets
`(-Inf, Inf)`; callers that need finite bounds (nested sampling) should check for that.

Angular sizes and offsets are the only bounds worth tuning, and they are the ones the data can
tune. Pass `data` and the ceiling becomes `2 λ/B_min` from the actual uv coverage — the array
cannot constrain a source larger than the scale its shortest baseline senses. Without `data`
the ceiling is [`DEFAULT_MAX_SIZE_MAS`](@ref); `max_size` overrides both.

```julia
lb, ub = default_bounds(model_dict, free)              # 100 mas ceiling
lb, ub = default_bounds(model_dict, free; data = data) # sized to this uv coverage
```

Everything else is physical rather than instrumental — fluxes in `[0,Inf)`, since a flux is a
weight relative to the others, `thick` in `[0,1]`, inclination in `[0,90]`, position angles in
`[-180,180]` — and does not depend on the data. A nested sampler needs a finite upper bound on
every free flux, since the bounds are its prior.

See also [`display_model`](@ref), which validates values against whatever bounds you end up
with, and [`max_angular_scale`](@ref).
"""
function default_bounds(model_dict::Dict{String}, list_free_params::AbstractVector{<:AbstractString};
                        data = nothing, max_size = nothing)
    # Size the angular bounds from the data when it is available: the coverage says what can
    # be constrained far better than any constant can. Twice λ/B_min leaves room for a
    # marginally-resolved source without admitting sizes the array cannot measure.
    smax = if max_size !== nothing
        Float64(max_size)
    elseif data !== nothing
        max(2 * max_angular_scale(data), 1.0)
    else
        DEFAULT_MAX_SIZE_MAS
    end

    lb = Dict{String,Float64}()
    ub = Dict{String,Float64}()
    for p in list_free_params
        i = findfirst(==(','), p)
        suffix = i === nothing ? p : p[nextind(p, i):end]
        b = if suffix in _SIZE_KEYS
                (0.0, smax)
            elseif suffix in _OFFSET_KEYS
                (-smax, smax)
            elseif haskey(_DEFAULT_PARAM_BOUNDS, suffix)
                _DEFAULT_PARAM_BOUNDS[suffix]
            elseif startswith(suffix, "az amp")
                (0.0, 1.0)
            elseif startswith(suffix, "az projang")
                (-180.0, 180.0)
            else
                (-Inf, Inf)
            end
        lb[p], ub[p] = b
    end
    return lb, ub
end


# ─────────────────────────────────────────────────────────────────────────────
# Analytic visibility functions
# ─────────────────────────────────────────────────────────────────────────────
# All take (params::Vector{Float64}, ρ_rad::Vector{Float64}) -> Vector{Float64}
# where ρ_rad is in cycles/rad and params are in mas / dimensionless.
#
# These call the ChainRules-compatible implementations from model_chainrules.jl.

function _vis_ud(params::AbstractVector, ρ::AbstractVector)
    vis_ud(params[1], ρ)
end

function _vis_gaussian(params::AbstractVector, ρ::AbstractVector)
    # V = exp(-π² FWHM² ρ² / (4 ln2)), the standard Gaussian transform: mfit computes
    # `exp(-((pi*a*rho)**2)/(4*log(2)))` with a = FWHM, and a numerical Hankel transform of
    # exp(-4 ln2 r²/FWHM²) agrees with it to 4e-8.
    #
    # The `/ln2` form takes a HALF-width, not σ. Passing σ = FWHM/(2√(2ln2)) made the component
    # narrower than its stated FWHM by √(2ln2) = 1.177 -- measured in the image domain, a
    # Gaussian asked for FWHM 3 mas came out 2.52 mas -- so a fitted FWHM read ~18% high.
    fwhm = params[1]
    hwhm = fwhm / 2
    t    = @. π * hwhm * ρ / _MAS2RAD_PM
    return @. exp(-t^2 / log(2))   # = exp(-π² FWHM² ρ²/MAS2RAD²/(4 ln2))
end

function _vis_ldlin(params::AbstractVector, ρ::AbstractVector)
    vis_ldlin(params[1], params[2], ρ)
end

function _vis_ldquad(params::AbstractVector, ρ::AbstractVector)
    vis_ldquad(params[1], params[2], params[3], ρ)
end

function _vis_ldsqrt(params::AbstractVector, ρ::AbstractVector)
    vis_ldsqrt(params[1], params[2], params[3], ρ)
end

function _vis_ldclaret4(params::AbstractVector, ρ::AbstractVector)
    vis_ldclaret4(params[1], params[2], params[3], params[4], params[5], ρ)
end

function _vis_ldpow(params::AbstractVector, ρ::AbstractVector)
    vis_ldpow(params[1], params[2], ρ)
end

# ForwardDiff specialization: the order ν = α/2 + 1 moves with α, so the derivatives come from
# `_ldpow_d`'s closed forms and are packed into Duals. ρ carries partials too when an
# inclination or position angle is differentiated.
function _vis_ldpow(params::AbstractVector{D}, ρ::AbstractVector) where {T, V, N, D<:ForwardDiff.Dual{T,V,N}}
    θ  = ForwardDiff.value(params[1])
    α  = ForwardDiff.value(params[2])
    ρv = ForwardDiff.value.(ρ)
    W, dV_dζ, dV_dα = _ldpow_vjac(θ, α, ρv)
    θ_p = ForwardDiff.partials(params[1])
    α_p = ForwardDiff.partials(params[2])
    return [ForwardDiff.Dual{T}(W[j], dV_dζ[j] * π * ρv[j] / MAS2RAD * θ_p + dV_dα[j] * α_p +
                                      dV_dζ[j] * π * θ / MAS2RAD * _partials_of(ρ[j], D))
            for j in eachindex(ρ)]
end
_partials_of(r::ForwardDiff.Dual, ::Type{D}) where {D} = ForwardDiff.partials(r)
_partials_of(r::Real, ::Type{D}) where {D} = zero(ForwardDiff.partials(zero(D)))

function _vis_point(params::AbstractVector, ρ::AbstractVector)
    ones(eltype(ρ), length(ρ))
end

# Fully resolved (a background): its own flux at zero baseline, 1, and 0 at every baseline the
# array measures. The spike at zero is what puts its flux in an image, as a uniform level.
function _vis_resolved(params::AbstractVector, ρ::AbstractVector)
    @. ifelse(iszero(ρ), one(eltype(ρ)), zero(eltype(ρ)))
end

# Uniform ring: area-weighted difference of two uniform disks
function _vis_ring(params::AbstractVector, ρ::AbstractVector)
    din  = params[1]   # inner diameter (mas)
    dout = params[2]   # outer diameter (mas)
    V_out = vis_ud(dout, ρ)
    V_in  = vis_ud(din, ρ)
    return @. (dout^2 * V_out - din^2 * V_in) / (dout^2 - din^2)
end

# Bi-Gaussian ring: flux-weighted difference of two Gaussians
function _vis_gaussian_ring(params::AbstractVector, ρ::AbstractVector)
    fwhmin  = params[1]
    fwhmout = params[2]
    # Half-widths, for the reason given in `_vis_gaussian`. The flux weighting below is a ratio
    # of these squared, so the correction cancels there and only the widths move.
    σ_in  = fwhmin  / 2
    σ_out = fwhmout / 2
    t_in  = @. π * σ_in  * ρ / _MAS2RAD_PM
    t_out = @. π * σ_out * ρ / _MAS2RAD_PM
    V_in  = @. exp(-t_in^2  / log(2))
    V_out = @. exp(-t_out^2 / log(2))
    return @. (σ_out^2 * V_out - σ_in^2 * V_in) / (σ_out^2 - σ_in^2)
end

# Crescent: two offset uniform disks (returns Complex visibility).
# Takes (u,v) separately because the internal offset creates direction-dependent
# phase shifts.  Called from eval_model with a special branch.
function _vis_crescent(params::AbstractVector, u::AbstractVector, v::AbstractVector)
    crin      = params[1]  # inner diameter (mas)
    crout     = params[2]  # outer diameter (mas)
    croff     = params[3]  # offset factor 0..1 (0 = concentric ring)
    crprojang = params[4]  # PA of thinnest part (deg, 0=N, 90=E)

    # broadcast throughout: any of these may be a per-point vector when the
    # parameter is given as a $WL / $MJD / $BL expression
    crpa = @. crprojang * π / 180
    off  = @. croff * (crout - crin) / 2   # offset distance (mas)

    # Outer disk shifted away from thin side, inner toward thin side
    offXo = @. -off / 2 * sin(crpa)
    offYo = @. -off / 2 * cos(crpa)
    offXi = @.  off / 2 * sin(crpa)
    offYi = @.  off / 2 * cos(crpa)

    ρ = @. sqrt(u^2 + v^2)
    V_out = vis_ud(crout, ρ)
    V_in  = vis_ud(crin,  ρ)

    phase_out = @. -2π * (u * offXo + v * offYo) / _MAS2RAD_PM
    phase_in  = @. -2π * (u * offXi + v * offYi) / _MAS2RAD_PM

    return @. (crout^2 * V_out * cis(phase_out) -
               crin^2  * V_in  * cis(phase_in)) / (crout^2 - crin^2)
end

const _ANALYTIC_VIS = Dict{Symbol, Function}(
    :ring          => _vis_ring,
    :gaussian_ring => _vis_gaussian_ring,
    :ldsqrt        => _vis_ldsqrt,
    :ldclaret4     => _vis_ldclaret4,
    # :crescent handled specially in eval_model (needs u,v not just ρ)
    :ud            => _vis_ud,
    :gaussian      => _vis_gaussian,
    :ldlin         => _vis_ldlin,
    :ldquad        => _vis_ldquad,
    :ldpow         => _vis_ldpow,
    :resolved      => _vis_resolved,
    :point         => _vis_point,
)

# The same table as methods, so `AnalyticSpec{K}` reaches its kernel by static dispatch.
for (k, fn) in _ANALYTIC_VIS
    @eval @inline _vis_kernel(::Val{$(QuoteNode(k))}, params, ρ) = $fn(params, ρ)
end


# ─────────────────────────────────────────────────────────────────────────────
# Component specification structs
# ─────────────────────────────────────────────────────────────────────────────

abstract type AbstractComponentSpec end

"""
Analytic component: uniform disk, Gaussian, limb-darkened disk, point, resolved. `K` is the
kind as a type parameter, so evaluation dispatches to the kind's visibility function
statically once the component's own type is known.
"""
struct AnalyticSpec{K} <: AbstractComponentSpec
    comp_name ::String
    kind      ::Symbol
    vis_fn    ::Function          # (params::Vector, ρ_rad::Vector) -> Vector
    f_idx     ::Int               # index into resolver output (flux fraction)
    param_idx ::Vector{Int}       # indices of geometry params in resolver output
    incl_idx  ::Int               # inclination in deg (0 = not present)
    pa_idx    ::Int               # position angle in deg (0 = not present)
    x_idx     ::Int               # x offset in mas (0 = not present)
    y_idx     ::Int               # y offset in mas (0 = not present)
end

struct AzMode
    order    ::Int     # harmonic order (1, 2, 3, …)
    amp_idx  ::Int     # index into resolver output for amplitude
    phi_idx  ::Int     # index into resolver output for phase (degrees)
end

# A component's profile as the quadrature sees it: a function of s = r/r_max, whose first two
# parameters are the inner and outer diameters, `(s, μ, d_min, d_max, p...) -> I` at
# R = (d_max/2)·s. `f` is the profile compiled with `size_args`. The size is then a parameter
# like any other: the probe sees breakpoints fixed in R move in s, and ForwardDiff takes
# ∂I/∂d at fixed s.
struct _SizedProfile{F} <: Function
    f::F
end
(p::_SizedProfile)(s, μ, dmin, dmax, pp...) = p.f((dmax / 2) .* s, μ, dmin / 2, dmax / 2, pp...)

"""
Hankel component: arbitrary radial profile.

The profile is compiled to an AD-transparent closure at dict_to_model time. It is evaluated
on the nodes of `_hankel_nodes` over [r_min, r_max]: panels between the breakpoints the
profile shows, and a count that follows the baselines the component sees and the profile's
own convergence (`_hankel_setup`), unless `c,nr` fixes it. The diameters are read each call
(`dmin_idx`, `dmax_idx`, or the numbers `d_min`, `d_max` when no key gives them), so the size
can be fitted; the quadrature is in s = r/r_max, and the profile's parameter vector is
`[d_min, d_max, profile params...]`.
"""
struct HankelSpec{T<:AbstractFloat} <: AbstractComponentSpec
    comp_name        ::String
    f_idx            ::Int
    profile_fn       ::_SizedProfile    # (s, mu, d_min, d_max, profile_params...) -> I
    profile_param_idx::Vector{Int}      # indices in resolver output
    dmin_idx         ::Int              # `c,diamin` in the resolver (0 = d_min is fixed)
    dmax_idx         ::Int              # `c,udout`/`diamout`/`diam` (0 = `c,r_max`, fixed)
    d_min            ::T                # the diameters the component was built with, mas
    d_max            ::T
    r_max            ::T                # d_max / 2
    s_min            ::T                # r_min / r_max
    nr               ::Int              # `c,nr`: a fixed node count, or 0 = from the baselines
    workspace        ::HankelWorkspace{T}
    az_modes         ::Vector{AzMode}   # azimuthal variation modes (empty = none)
    kcache           ::HankelKernelCache{T}  # the nodes, per count; the matrices, per uv set
    incl_idx         ::Int
    pa_idx           ::Int
    x_idx            ::Int
    y_idx            ::Int
end

"""
A component whose VISIBILITY is written out, rather than its brightness.

The counterpart of `HankelSpec`: that one compiles I(r) and transforms it, this one is handed
V(B) and evaluates it. `\$BL` is the baseline in Mλ, which is the unit the observable plots use
and the one a formula from a paper is normally written in.
"""
struct VisFuncSpec <: AbstractComponentSpec
    comp_name    ::String
    f_idx        ::Int
    vis_fn       ::Function       # (B_Mlambda, WL_m, MJD, params...) -> V
    param_idx    ::Vector{Int}
    incl_idx     ::Int
    pa_idx       ::Int
    x_idx        ::Int
    y_idx        ::Int
end

"""
A FITS image or cube as a component: `c,image` names the file, `c,pixsize` (mas) overrides
its header and `c,norm` (`plane` or `relative`) its `NORMALIZ`.

The image is read once, when the model is compiled, and its visibilities are cached per uv set
(`ImageVisCache`), since nothing a fit varies changes them: `f`, `x` and `y` are applied after.
Inclination and position angle would, and are refused. `incl_idx`/`pa_idx` are always 0 and
exist so that the evaluation loop reads every component the same way.
"""
struct ImageSpec <: AbstractComponentSpec
    comp_name ::String
    f_idx     ::Int
    image     ::ModelImage
    cache     ::ImageVisCache
    incl_idx  ::Int
    pa_idx    ::Int
    x_idx     ::Int
    y_idx     ::Int
end

"""Top-level compiled model."""
struct FlatModel
    resolver   ::RGF.Resolver                  # concrete; its generated functions are the dynamic calls
    components ::Vector{AbstractComponentSpec}
    all_names  ::Vector{String}                # resolver output names (for debug)
    list_free_params ::Vector{String}
    kernel_idx ::Int                           # spatial_kernel FWHM index (0 = none)
    smear      ::Int                           # sub-wavelengths per channel (< 2 = none, -1 = auto)
    smear_radius   ::Float64                   # mas: the model's extent, for smear = "auto"
    smear_coherent ::Bool                      # average V over the channel, not the observables
    chromatic  ::Vector{Bool}                  # per component: any of its parameters per uv point
end

"""
    ResolvedParams{T}

A model's parameter values for one evaluation, concretely typed whatever the model: scalars
in `s`, per-uv-point values in `cols` (one `Vector{T}` per chromatic parameter), and `slot`,
the resolver's map from a parameter's index in `all_names` to `s[k]` (`k > 0`) or `cols[-k]`.
"""
struct ResolvedParams{T}
    s    :: Vector{T}
    cols :: Vector{Vector{T}}
    slot :: Vector{Int}
    n    :: Int
end

# The resolver is the model's one dynamic call: a RuntimeGeneratedFunction, a new one for
# every compiled model. It writes into storage typed by `T` alone, so nothing downstream
# depends on the model instance and a rebuilt model reuses every compiled method.
function _resolve(model::FlatModel, x, wl, mjd, n::Int, ::Type{T}) where {T}
    r = model.resolver
    s = Vector{T}(undef, r.nscalar)
    cols = [Vector{T}(undef, n) for _ in 1:r.nchromatic]
    r.fill!(s, cols, x, wl, mjd)
    return ResolvedParams{T}(s, cols, r.slot, n)
end

# A parameter's value: a scalar for a component evaluated on the scalar path (`Val(false)`),
# one value per uv point on the chromatic path (`Val(true)`), where a scalar is spread.
@inline _par(rp::ResolvedParams, i::Int, ::Val{false}) = rp.s[rp.slot[i]]
@inline function _par(rp::ResolvedParams, i::Int, ::Val{true})
    k = rp.slot[i]
    return k > 0 ? fill(rp.s[k], rp.n) : rp.cols[-k]
end


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

"""Extract unique component name prefixes from a flat param dict."""
function _component_names(model_dict::Dict{String})
    names = String[]
    for k in keys(model_dict)
        idx = findfirst(',', k)
        isnothing(idx) && continue
        push!(names, k[1:idx-1])
    end
    # Sorted, because `keys(::Dict)` has no defined order and changes it when keys are
    # inserted: adding `incl` and `pa` to one component reshuffled the whole list, so the
    # panel's component order moved when nothing about the model had.
    return sort!(unique(names))
end

"""Identify the component kind from which geometry key is present."""
function _identify_kind(comp_name::String, model_dict::Dict{String})::Symbol
    for gk in _GEOMETRY_KEYS
        haskey(model_dict, "$comp_name,$gk") && return get(
            Dict("image"=>:image, "visfunc"=>:visfunc,
                 "profile"=>:hankel, "diamin"=>:ring, "fwhmin"=>:gaussian_ring,
                 "ldsqrt"=>:ldsqrt, "ldclaret4"=>:ldclaret4,
                 "crin"=>:crescent, "ud"=>:ud, "fwhm"=>:gaussian,
                 "ldlin"=>:ldlin, "ldquad"=>:ldquad, "ldpow"=>:ldpow,
                 "resolved"=>:resolved), gk, :unknown)
    end
    # No geometry key — must have f, treat as unresolved point
    haskey(model_dict, "$comp_name,f") && return :point
    error("Cannot identify component type for '$comp_name': " *
          "no recognized geometry key (ud, fwhm, profile, ldlin, ldquad, ldpow)")
end

"""
Build a name→index lookup from the resolver's all_names vector.
Used at parse time; result is stored as Int indices in ComponentSpec.
"""
function _name_to_idx(all_names::Vector{String})::Dict{String,Int}
    Dict(n => i for (i, n) in enumerate(all_names))
end

"""
Recursively resolve a model_dict value to a numeric Float64.
If it's already a Number, return it. If it's a string expression containing
\$-references, substitute them recursively and evaluate.
Used at parse time for grid parameters (diamin, diamout, etc.) that may be expressions.
"""
function _resolve_numeric(key::AbstractString, model_dict::Dict{String}, visited::Set{String}=Set{String}())::Float64
    k = String(key)
    haskey(model_dict, k) || error("Key '$k' not found in model_dict")
    v = model_dict[k]
    v isa Number && return Float64(v)
    v isa AbstractString || error("Key '$k' has unsupported type $(typeof(v))")
    k in visited && error("Circular reference detected for '$k'")
    push!(visited, k)
    expr = v
    # Replace $ref with resolved numeric values
    resolved = replace(expr, _REF_PATTERN => function(m)
        ref = m[2:end]  # strip leading $
        string(_resolve_numeric(ref, model_dict, visited))
    end)
    result = try
        eval(Meta.parse(resolved))
    catch e
        error("Cannot evaluate expression '$resolved' (from key '$key'): $e")
    end
    return Float64(result)
end

"""
Determine the r grid outer boundary for a Hankel component.
Priority: "udout" key (outer diameter/2), then "r_max" key.
Errors if neither is present.
"""
function _hankel_r_max(comp_name::String, model_dict::Dict{String})::Float64
    udout_key   = "$comp_name,udout"
    diam_key    = "$comp_name,diam"
    diamout_key = "$comp_name,diamout"
    r_max_key   = "$comp_name,r_max"
    if haskey(model_dict, udout_key)
        return _resolve_numeric(udout_key, model_dict) / 2.0
    elseif haskey(model_dict, diamout_key)
        return _resolve_numeric(diamout_key, model_dict) / 2.0
    elseif haskey(model_dict, diam_key)
        return _resolve_numeric(diam_key, model_dict) / 2.0
    elseif haskey(model_dict, r_max_key)
        return _resolve_numeric(r_max_key, model_dict)
    else
        error("Hankel component '$comp_name' requires '$udout_key', '$diam_key', " *
              "'$diamout_key', or '$r_max_key' to define the r grid.")
    end
end

"""
The node count `c,nr` fixes for a Hankel component, or 0 when it is chosen from the baselines
the component sees (`_hankel_n`).
"""
function _hankel_Nr(comp_name::String, model_dict::Dict{String})::Int
    k = "$comp_name,nr"
    haskey(model_dict, k) || return 0
    n = Int(model_dict[k])
    n >= 2 || error("$k must be at least 2 quadrature nodes, got $n")
    return n
end

"""Determine r_min for a Hankel component's r grid (0 by default, diamin/2 for rings)."""
function _hankel_r_min(comp_name::String, model_dict::Dict{String})::Float64
    diamin_key = "$comp_name,diamin"
    if haskey(model_dict, diamin_key)
        return _resolve_numeric(diamin_key, model_dict) / 2.0
    end
    return 0.0
end

"""
Keys to exclude from the scalar resolver dict.
Profile strings are compiled separately; udout/r_max/nr are grid config,
not runtime params.
"""
const _HANKEL_META_SUFFIXES = ["profile", "r_max", "nr"]

function _resolver_dict(model_dict::Dict{String}, comp_names::Vector{String})
    # Exclude profile strings, Hankel grid-config keys, and resolved markers
    # from the resolver.  udout IS kept — it may be a fit param.
    # `smear` is a setting read when the model is compiled, not a parameter.
    exclude = Set{String}(["smear", "smear_mode"])
    for cn in comp_names
        push!(exclude, "$cn,profile")
        # A visibility expression is compiled by `compile_visfunc` against the baseline, not
        # resolved against `(x, WL, MJD)`. Left in, the resolver would try to evaluate it as a
        # derived parameter and fail on the first function the expression names.
        push!(exclude, "$cn,visfunc")
        push!(exclude, "$cn,r_max")
        push!(exclude, "$cn,nr")
        push!(exclude, "$cn,resolved")
        # An image component's file and how to read it are consumed when the model is
        # compiled; a path is not an expression.
        if haskey(model_dict, "$cn,image")
            push!(exclude, "$cn,image", "$cn,pixsize", "$cn,norm")
        end
    end
    return Dict{String,Any}(k => v for (k, v) in model_dict if k ∉ exclude)
end


# ─────────────────────────────────────────────────────────────────────────────
# Sugar preprocessing
# ─────────────────────────────────────────────────────────────────────────────

"""
    _preprocess_sugar(model_dict, comp_names) -> Dict{String,Any}

Expand shorthand notations to canonical keys:
- `diam + thick` → `diamin + diamout`  (when no `profile` present)
"""
function _preprocess_sugar(model_dict::Dict{String}, comp_names::Vector{String})
    d = copy(model_dict)
    for cn in comp_names
        diam_key  = "$cn,diam"
        thick_key = "$cn,thick"
        # Only expand if both diam and thick are present, no profile, and
        # diamin/diamout are not already defined
        if haskey(d, diam_key) && haskey(d, thick_key) &&
           !haskey(d, "$cn,profile") &&
           !haskey(d, "$cn,diamin") && !haskey(d, "$cn,diamout")
            dv = d[diam_key]
            tv = d[thick_key]
            if dv isa Number && tv isa Number
                d["$cn,diamout"] = Float64(dv)
                d["$cn,diamin"]  = Float64(dv * (1 - tv))
            else
                # Expression-based: create derived expressions
                d["$cn,diamout"] = "\$$cn,diam"
                d["$cn,diamin"]  = "\$$cn,diam * (1 - \$$cn,thick)"
            end
            # Keep diam and thick — they may be fit params referenced by the expressions
        end

        # ── spectrum → f  (wavelength-dependent flux via $WL) ─────────────
        spec_key = "$cn,spectrum"
        f_key    = "$cn,f"
        if haskey(d, spec_key) && !haskey(d, f_key)
            d[f_key] = d[spec_key]      # promote spectrum expression to f
            delete!(d, spec_key)
        end
    end
    return d
end


# ─────────────────────────────────────────────────────────────────────────────
# dict_to_model
# ─────────────────────────────────────────────────────────────────────────────


# Suffixes every component may carry, the sugar `dict_to_model` expands, and the keys a radial
# profile reads on top of its expression.
const _UNIVERSAL_SUFFIXES = Set(["f", "x", "y", "incl", "pa", "projang"])
const _SUGAR_SUFFIXES     = Set(["diam", "thick", "spectrum"])
const _HANKEL_SUFFIXES    = Set(["profile", "r_max", "nr", "udout", "diamin", "diamout"])
const _IMAGE_SUFFIXES     = Set(["image", "pixsize", "norm"])

"""
    reserved_suffix(suffix) -> Bool

Whether the parser gives `component,suffix` a meaning of its own, for ANY kind of component.

A parameter added under such a name would not be a new parameter. `ud` on a Gaussian turns it
into a uniform disc, because `_identify_kind` reads the shape off the keys, and `u` is a
limb-darkening coefficient. Anything else is free for a component to carry as a parameter of
its own, used from an expression.
"""
reserved_suffix(suffix::AbstractString) =
    suffix in _UNIVERSAL_SUFFIXES || suffix in _SUGAR_SUFFIXES || suffix in _HANKEL_SUFFIXES ||
    suffix in _IMAGE_SUFFIXES || suffix in _GEOMETRY_KEYS ||
    any(s -> suffix in s, values(_ANALYTIC_PARAM_SUFFIXES)) ||
    startswith(suffix, "az amp") || startswith(suffix, "az projang")

"""
Warn about `component,key` pairs the parser never consumes.

A misspelt geometry key does not error — it changes what the component *is*. Writing
`"ring,gaussian_ring"` alongside `"ring,fwhm"` yields a plain Gaussian, silently, because
`_identify_kind` matches on `fwhm` and the unknown key is ignored. Two shipped demos did
exactly this and fitted the wrong model for as long as they existed.

Keys referenced from a `\$`-expression anywhere in the dict are always legitimate (that is
how radial-profile parameters and shared globals work), so they are never reported. A reference
that is one of the expression's implicit variables (`\$BL` in a visibility function, `\$R` in a
profile) is that variable, so it does not count as a use of a component's `B` or `R`.
"""
function _warn_unrecognised_keys(pd::Dict{String})
    referenced = Set{String}()
    for (k, v) in pd
        v isa AbstractString || continue
        union!(referenced, setdiff(extract_refs(v), implicit_vars(k)))
    end

    unknown = String[]
    for cn in _component_names(pd)
        kind = _identify_kind(cn, pd)
        ok = union(_UNIVERSAL_SUFFIXES, _SUGAR_SUFFIXES,
                   Set{String}(get(_ANALYTIC_PARAM_SUFFIXES, kind, String[])))
        kind === :hankel && union!(ok, _HANKEL_SUFFIXES)
        kind === :image  && union!(ok, _IMAGE_SUFFIXES)
        gk = get(_KIND_GEOMETRY, kind, "")
        isempty(gk) || push!(ok, gk)

        prefix = cn * ","
        for k in keys(pd)
            startswith(k, prefix) || continue
            suffix = k[(length(prefix) + 1):end]
            suffix in ok && continue
            startswith(suffix, "az amp") && continue
            startswith(suffix, "az projang") && continue
            (suffix in referenced || k in referenced) && continue
            push!(unknown, k)
        end
    end

    if !isempty(unknown)
        @warn "Model keys not used by the parser — check for typos. An unrecognised " *
              "geometry key does not error, it silently changes the component type." *
              "" keys = sort(unknown)
    end
    return unknown
end

"""
    dict_to_model(model_dict, list_free_params; nB_workspace=100) -> FlatModel

Compile a flat parameter dict into a `FlatModel` ready for evaluation.

Arguments
---------
model_dict   : Dict{String,Any} with "component,param" keys.
               Values may be Float64 (fixed or free) or String (expression).
               Special string key "component,profile" triggers the Hankel pathway.
list_free_params  : Vector{String} of keys to optimize (must be numeric in model_dict).
nB_workspace : Pre-allocated baseline grid size for HankelWorkspace buffers.
               Should be ≥ the number of baselines you will evaluate at.

Returns
-------
A FlatModel whose eval_model / eval_model_grad methods can be called in a loop
with zero Dict allocation.

Example
-------
    model_dict = Dict(
        "star,ud"      => 0.8,
        "star,f"       => 0.6,
        "ring,profile" => "exp(-(\\\$R / \\\$scale)^2 / 2)",
        "ring,scale"   => 1.5,
        "ring,udout"   => 12.0,
        "ring,f"       => "1 - \\\$star,f",
    )
    list_free_params = ["star,ud", "star,f", "ring,scale"]
    model = dict_to_model(model_dict, list_free_params)
"""
function dict_to_model(model_dict::Dict{String},
                      list_free_params::Vector{String};
                      nB_workspace::Int = 100,
                      T::Type{<:AbstractFloat} = Float64)::FlatModel

    comp_names = _component_names(model_dict)

    # ── Expand sugar (diam+thick → diamin+diamout) ───────────────────────────
    pd = _preprocess_sugar(model_dict, comp_names)

    # ── Report keys the parser will ignore ───────────────────────────────────
    _warn_unrecognised_keys(pd)

    # ── Build resolver dict (numeric params only) ────────────────────────────
    res_dict = _resolver_dict(pd, comp_names)

    # ── Build RGF resolver ───────────────────────────────────────────────────
    resolver = RGF.build_resolver(res_dict, list_free_params)
    n2i      = _name_to_idx(resolver.all_names)

    # ── Spatial kernel ───────────────────────────────────────────────────────
    kernel_idx = get(n2i, "spatial_kernel", 0)

    # ── Bandwidth smearing ───────────────────────────────────────────────────
    "smear" in list_free_params &&
        error("smear is the number of sub-wavelengths per channel, a setting, not a parameter " *
              "to fit")
    sm = get(pd, "smear", 0)
    sm = sm == "auto" ? -1 :
         (sm isa Real && isinteger(sm) && sm >= 0) ? Int(sm) :
         error("smear must be a whole number of sub-wavelengths per channel (0 or 1 for none) " *
               "or \"auto\", not $(repr(sm))")
    smode = get(pd, "smear_mode", "incoherent")
    smode in ("incoherent", "coherent") ||
        error("smear_mode must be \"incoherent\" (the observables averaged over a channel) or " *
              "\"coherent\" (the visibility averaged, then the observables), not $(repr(smode))")
    smear_radius = sm == -1 ? _model_extent(pd, comp_names) : 0.0

    # ── Build component specs ────────────────────────────────────────────────
    components = AbstractComponentSpec[]

    for cn in comp_names
        kind = _identify_kind(cn, pd)

        # Flux fraction index (default to 1.0 if not in resolver)
        f_key = "$cn,f"
        f_idx = get(n2i, f_key, 0)   # 0 means "use 1.0"

        # ── Geometry modifiers (shared by all component types) ───────────
        incl_idx = get(n2i, "$cn,incl", 0)
        pa_idx   = get(n2i, "$cn,pa", 0)
        if pa_idx == 0
            pa_idx = get(n2i, "$cn,projang", 0)   # PMOIRED alias
        end
        x_idx = get(n2i, "$cn,x", 0)
        y_idx = get(n2i, "$cn,y", 0)

        if kind == :image
            push!(components, _image_spec(cn, pd, f_idx, incl_idx, pa_idx, x_idx, y_idx))

        elseif kind == :visfunc
            vis_expr = pd["$cn,visfunc"]
            vis_expr isa String ||
                error("$cn,visfunc must be a String expression")

            # Same discovery rule as a profile: anything that is not an implicit variable is a
            # parameter, resolved component-qualified first and then global.
            # `\$BL` is implicit HERE and nowhere else. It is not in `IMPLICIT_VARS` because the
            # resolver's signature is `(x, WL, MJD)` and has no baseline to give a derived
            # parameter; a visibility expression is the one place a baseline exists, and it is
            # supplied by `eval_model` rather than resolved.
            raw_refs = extract_refs(vis_expr)
            vf_refs  = filter(r -> r ∉ IMPLICIT_VARS && r != "BL", raw_refs)
            vf_names = map(vf_refs) do r
                if occursin(",", r)
                    r
                else
                    q = "$cn,$r"
                    haskey(pd, q) ? q : r
                end
            end
            vis_fn = compile_visfunc(vis_expr, vf_refs)
            vf_idx = Int[]
            for nm in vf_names
                haskey(n2i, nm) ||
                    error("Visibility parameter '$nm' (from '$cn,visfunc') not found in " *
                          "resolver. Add it to the model dict.")
                push!(vf_idx, n2i[nm])
            end

            push!(components, VisFuncSpec(cn, f_idx, vis_fn, vf_idx,
                                          incl_idx, pa_idx, x_idx, y_idx))

        elseif kind == :hankel

            # Profile expression and parameter names
            profile_expr = pd["$cn,profile"]
            profile_expr isa String ||
                error("$cn,profile must be a String expression")

            # Extract param names referenced in the profile (excluding implicit vars)
            raw_refs  = extract_refs(profile_expr)
            prof_refs = filter(r -> r ∉ IMPLICIT_VARS, raw_refs)
            # Qualify unqualified refs with component name, but fall back to
            # global name if the component-qualified name doesn't exist
            prof_param_names = map(prof_refs) do r
                if occursin(",", r)
                    r
                else
                    qualified = "$cn,$r"
                    haskey(pd, qualified) ? qualified : r
                end
            end

            # Compile the profile closure (AD-transparent), the size among its arguments:
            # `$RMIN`, `$RMAX`, `$DMIN` and `$DMAX` follow the diameters each call. The size the
            # component is built with sets its nodes until the diameters change.
            r_min = _hankel_r_min(cn, pd)
            r_max = _hankel_r_max(cn, pd)
            r_max > r_min >= 0 ||
                error("Hankel component '$cn': the profile needs 0 ≤ r_min < r_max, " *
                      "got r_min = $r_min, r_max = $r_max")
            Nr    = _hankel_Nr(cn, pd)

            profile_fn = _SizedProfile(compile_profile(profile_expr, prof_refs; size_args = true))
            # The diameters' keys in the resolver, read each call; `c,r_max` is a fixed setting.
            dmin_idx = get(n2i, "$cn,diamin", 0)
            kmax = findfirst(k -> haskey(pd, "$cn,$k"), ("udout", "diamout", "diam"))
            dmax_idx = kmax === nothing ? 0 :
                       get(n2i, "$cn," * ("udout", "diamout", "diam")[kmax], 0)

            # Indices of profile params in resolver output
            prof_idx = map(prof_param_names) do k
                haskey(n2i, k) ||
                    error("Profile param '$k' not found in resolver. " *
                          "It must be a numeric param in model_dict.")
                n2i[k]
            end

            ws = HankelWorkspace(Nr > 0 ? Nr : _hankel_n(0), nB_workspace, length(prof_idx); T=T)

            # Detect azimuthal variation modes: "cn,az amp1", "cn,az projang1", etc.
            az_modes = AzMode[]
            for k in keys(pd)
                m = match(Regex("^$(cn),az amp(\\d+)\$"), k)
                if m !== nothing
                    order = parse(Int, m.captures[1])
                    amp_key = "$cn,az amp$order"
                    phi_key = "$cn,az projang$order"
                    haskey(n2i, amp_key) ||
                        error("Azimuthal param '$amp_key' not found in resolver.")
                    haskey(n2i, phi_key) ||
                        error("Azimuthal param '$phi_key' not found in resolver " *
                              "(required when '$amp_key' is present).")
                    push!(az_modes, AzMode(order, n2i[amp_key], n2i[phi_key]))
                end
            end
            sort!(az_modes; by = m -> m.order)

            push!(components, HankelSpec(
                cn, f_idx, profile_fn, prof_idx, dmin_idx, dmax_idx, T(2r_min), T(2r_max),
                T(r_max), T(r_min / r_max), Nr, ws, az_modes,
                HankelKernelCache{T}(),
                incl_idx, pa_idx, x_idx, y_idx
            ))

        elseif kind == :crescent
            # Crescent: handled by AnalyticSpec with kind=:crescent
            # Uses _vis_crescent in eval_model (needs u,v, not just ρ)
            suffixes  = _ANALYTIC_PARAM_SUFFIXES[:crescent]
            param_idx = map(suffixes) do s
                k = "$cn,$s"
                haskey(n2i, k) ||
                    error("Crescent param '$k' not found in resolver.")
                n2i[k]
            end

            push!(components, AnalyticSpec{:crescent}(
                cn, :crescent, _vis_point, f_idx, param_idx,
                incl_idx, pa_idx, x_idx, y_idx
            ))

        else
            # Analytic component (ud, gaussian, ring, gaussian_ring, ldlin, etc.)
            suffixes  = _ANALYTIC_PARAM_SUFFIXES[kind]
            param_idx = map(suffixes) do s
                k = "$cn,$s"
                haskey(n2i, k) ||
                    error("Analytic param '$k' not found in resolver.")
                n2i[k]
            end

            vis_fn = _ANALYTIC_VIS[kind]

            push!(components, AnalyticSpec{kind}(
                cn, kind, vis_fn, f_idx, param_idx,
                incl_idx, pa_idx, x_idx, y_idx
            ))
        end
    end

    # A component is evaluated on the chromatic path when any parameter it reads is per uv
    # point. A radial profile's own parameters must be scalars: its I(r) is on the r grid.
    ischrom(i) = i != 0 && resolver.slot[i] < 0
    chromatic = map(components) do c
        idx = [c.f_idx, c.incl_idx, c.pa_idx, c.x_idx, c.y_idx]
        if c isa AnalyticSpec || c isa VisFuncSpec
            append!(idx, c.param_idx)
        elseif c isa HankelSpec
            any(ischrom, c.profile_param_idx) && error(
                "$(c.comp_name): a radial profile's parameters cannot depend on \$WL or \$MJD " *
                "(its I(r) is computed once, on the r grid); give the component a chromatic " *
                "flux instead")
            (ischrom(c.dmin_idx) || ischrom(c.dmax_idx)) && error(
                "$(c.comp_name): a radial profile's size cannot depend on \$WL or \$MJD " *
                "(its nodes are laid out once per call, for one size)")
            for az in c.az_modes
                push!(idx, az.amp_idx, az.phi_idx)
            end
        end
        any(ischrom, idx)
    end
    return FlatModel(resolver, components, resolver.all_names, list_free_params, kernel_idx,
                     sm, smear_radius, smode == "coherent", chromatic)
end

# The image component's spec: the file read with the overrides its keys give.
function _image_spec(cn, pd, f_idx, incl_idx, pa_idx, x_idx, y_idx)
    path = pd["$cn,image"]
    path isa AbstractString ||
        error("$cn,image must be the path of a FITS file, not $(repr(path))")
    (incl_idx != 0 || pa_idx != 0) &&
        error("$cn is an image: it is used at the orientation it was made with, so it takes " *
              "no incl or pa. Rotate or tilt the image itself and save it again.")
    ps = get(pd, "$cn,pixsize", nothing)
    ps === nothing || ps isa Real ||
        error("$cn,pixsize must be a number of mas per pixel, not $(repr(ps))")
    nm = get(pd, "$cn,norm", nothing)
    norm = nm === nothing ? nothing :
           nm in ("plane", "relative") ? Symbol(nm) :
           error("$cn,norm must be \"plane\" or \"relative\", not $(repr(nm))")
    mi = read_model_image(path; pixsize = ps, normalization = norm)
    return ImageSpec(cn, f_idx, mi, ImageVisCache(), incl_idx, pa_idx, x_idx, y_idx)
end

# How far the model reaches from the phase centre, in mas, at the values the dict holds: each
# component's offset plus half its largest size. A value that cannot be resolved without data
# (an expression in `\$WL`, say) counts as 0. `smear = "auto"` sizes its sampling from this.
function _model_extent(pd, comp_names)
    num(k) = haskey(pd, k) ? (try _resolve_numeric(k, pd) catch; 0.0 end) : 0.0
    r = 0.0
    for cn in comp_names
        off  = hypot(num("$cn,x"), num("$cn,y"))
        sz  = maximum((abs(num("$cn,$s")) for s in _SIZE_KEYS); init = 0.0)
        if haskey(pd, "$cn,image")
            mi = try read_model_image(pd["$cn,image"]; pixsize = get(pd, "$cn,pixsize", nothing))
                 catch; nothing end
            mi === nothing || (sz = max(sz, size(mi.image, 1) * mi.pixsize))
        end
        r = max(r, off + sz / 2)
    end
    return r
end

"""Backward-compatible alias for [`dict_to_model`](@ref)."""
const parse_model = dict_to_model


# ─────────────────────────────────────────────────────────────────────────────
# eval_model  (forward pass, no gradient)
# ─────────────────────────────────────────────────────────────────────────────

"""
    eval_model(model, x, uv; wl=nothing, mjd=nothing, n=0, normalize=true) -> Vector{Complex{T}}

Evaluate the model visibility at the given parameter vector `x` and baselines `uv`:
Σ fᵢ Vᵢ / Σ fᵢ, normalised by the total flux at each point ([`model_flux`](@ref)), so a
component's `f` is its flux in any unit and V(0) = 1.

Arguments
---------
model : FlatModel from dict_to_model
x     : current free-parameter values, length = length(model.list_free_params)
uv    : 2×N matrix of (u,v) spatial frequencies in cycles/rad
wl    : wavelength per UV point (metres), length N  — enables `\$WL` (microns) in expressions
mjd   : MJD per UV point, length N  — enables `\$MJD` in expressions
n     : Bessel order for Hankel components (default 0)
plane_wl : the wavelength an image component chooses its plane by (default `wl`)
normalize : divide by the total flux (default). `false` returns Σ fᵢ Vᵢ, for a model that is
        one part of a larger sum, beside an image (SPARCO, VI).

`eval_model` returns one visibility per uv point and does not smear. A model with
`"smear" => n` is averaged over its channels by [`eval_model_sub`](@ref) and the χ² functions.

Returns
-------
V : complex visibility vector, length N
"""
function eval_model(model::FlatModel,
                    x::AbstractVector,
                    uv::AbstractMatrix;
                    wl  = nothing,
                    mjd = nothing,
                    n::Int = 0,
                    plane_wl = wl,
                    normalize::Bool = true)

    # `$WL`/`$MJD` are NaN when the caller passes none, as in a derived parameter.
    _wl  = isnothing(wl)  ? NaN : wl
    _mjd = isnothing(mjd) ? NaN : mjd

    # Precision follows the caller: the parameters and the uv coordinates. Float32 against
    # Float32 gives Float32, Float64 gives Float64, and a ForwardDiff Dual parameter vector
    # promotes to Dual so gradients flow. The resolved parameters are stored in it too.
    T  = eltype(x) <: Real ? promote_type(eltype(x), float(eltype(uv))) : Float64
    rp = _resolve(model, x, _wl, _mjd, size(uv, 2), T)
    return _eval_model(model, rp, uv, _wl, _mjd, n, plane_wl, normalize)
end

# Everything after the resolver, behind a function barrier: concretely typed by the
# parameters' element type and the uv matrix, never by the model instance.
function _eval_model(model::FlatModel, rp::ResolvedParams{T}, uv::AbstractMatrix,
                     _wl, _mjd, n::Int, plane_wl, normalize::Bool = true) where {T}
    nB = size(uv, 2)
    u_raw = @view uv[1,:]
    v_raw = @view uv[2,:]
    ρ  = @. sqrt(u_raw^2 + v_raw^2)   # cycles/rad
    B  = ρ ./ _MAS2RAD_PM             # cycles/mas (for Hankel)
    V  = zeros(Complex{T}, nB)

    # One dynamic dispatch per component, on its type and on whether any of its parameters is
    # per uv point (the chromatic path, where every one of its parameters is a vector).
    for (c, comp) in enumerate(model.components)
        _add_component!(V, comp, rp, Val(model.chromatic[c]), uv, u_raw, v_raw, ρ, B,
                        _wl, _mjd, n, plane_wl)
    end
    # Normalised by the total flux at each point: a component's `f` is its flux in any unit,
    # and fractions summing to one are the case where this divides by 1. Not when the model is
    # a part of a larger sum (`normalize = false`: an image beside it, SPARCO and VI).
    normalize && (V ./= _total_flux(model, rp, nB))

    # ── Spatial kernel (Gaussian smoothing in visibility space) ─────────
    if model.kernel_idx != 0
        k = rp.slot[model.kernel_idx]
        k > 0 ? _apply_kernel!(V, rp.s[k], ρ) : _apply_kernel!(V, rp.cols[-k], ρ)
    end
    return V
end

# Σ f over the components, one value per point (a flux may be chromatic), in the
# parameters' type.
function _total_flux(model::FlatModel, rp::ResolvedParams{T}, n::Int) where {T}
    F = zeros(T, n)
    for (c, comp) in enumerate(model.components)
        if comp.f_idx == 0
            F .+= one(T)
        elseif model.chromatic[c]
            F .+= _par(rp, comp.f_idx, Val(true))
        else
            F .+= _par(rp, comp.f_idx, Val(false))
        end
    end
    return F
end

"""
    model_flux(model, x; wl = nothing, mjd = nothing) -> Vector

The model's total flux, Σ f over its components, at each wavelength (and epoch) given, or
once without them. `eval_model` divides the visibilities by it; it is what an OI_FLUX
measurement is compared with.
"""
function model_flux(model::FlatModel, x::AbstractVector; wl = nothing, mjd = nothing)
    n   = wl !== nothing ? length(wl) : mjd !== nothing ? length(mjd) : 1
    _wl = isnothing(wl) ? NaN : wl
    _mjd = isnothing(mjd) ? NaN : mjd
    T  = eltype(x) <: Real ? float(eltype(x)) : Float64
    rp = _resolve(model, x, _wl, _mjd, n, T)
    return _total_flux(model, rp, n)
end

function _apply_kernel!(V, kern_fwhm, ρ)
    σ_kern = kern_fwhm / 2          # half-width; see `_vis_gaussian`
    t_kern = @. π * σ_kern * ρ / _MAS2RAD_PM
    @. V *= exp(-t_kern^2 / log(2))
    return V
end

@inline _one(::ResolvedParams{T}, ::Val{false}) where {T} = one(T)
@inline _one(rp::ResolvedParams{T}, ::Val{true}) where {T} = fill(one(T), rp.n)
@inline _zero(::ResolvedParams{T}, ::Val{false}) where {T} = zero(T)
@inline _zero(rp::ResolvedParams{T}, ::Val{true}) where {T} = fill(zero(T), rp.n)

function _add_component!(V, comp, rp, ch, uv, u_raw, v_raw, ρ, B, _wl, _mjd, n, plane_wl)
    # ── UV stretching for inclination / position angle ──────────────
    geom = (comp.incl_idx != 0 || comp.pa_idx != 0) ? _stretch(comp, rp, ch, u_raw, v_raw) :
                                                        (u_raw, v_raw, ρ, B)
    # ── Compute visibility for this component ──────────────────────
    Vi = _component_vis(comp, rp, ch, geom, uv, _wl, _mjd, n, plane_wl)
    f  = comp.f_idx == 0 ? _one(rp, ch) : _par(rp, comp.f_idx, ch)
    _shift_add!(V, f, Vi, comp, rp, ch, u_raw, v_raw)
    return V
end

function _stretch(comp, rp, ch, u_raw, v_raw)
    incl_rad = (comp.incl_idx != 0 ? _par(rp, comp.incl_idx, ch) : _zero(rp, ch)) * π / 180
    pa_rad   = (comp.pa_idx   != 0 ? _par(rp, comp.pa_idx,   ch) : _zero(rp, ch)) * π / 180
    # PMOIRED convention: rot = -PA, compress u' by cos(incl)
    # u' = cos(incl) * ( cos(PA)*u - sin(PA)*v )
    # v' = sin(PA)*u + cos(PA)*v
    u_comp = @. cos(incl_rad) * (cos(pa_rad) * u_raw - sin(pa_rad) * v_raw)
    v_comp = @. sin(pa_rad) * u_raw + cos(pa_rad) * v_raw
    ρ_comp = @. sqrt(u_comp^2 + v_comp^2)
    B_comp = ρ_comp ./ _MAS2RAD_PM
    return (u_comp, v_comp, ρ_comp, B_comp)
end

# ── Position offset → phase shift, and the flux-weighted sum ──────────
function _shift_add!(V, f, Vi, comp, rp, ch, u_raw, v_raw)
    if comp.x_idx != 0 || comp.y_idx != 0
        dx = comp.x_idx != 0 ? _par(rp, comp.x_idx, ch) : _zero(rp, ch)
        dy = comp.y_idx != 0 ? _par(rp, comp.y_idx, ch) : _zero(rp, ch)
        phase = @. -2π * (u_raw * dx + v_raw * dy) / _MAS2RAD_PM
        @. V += f * Vi * cis(phase)
    else
        @. V += f * Vi
    end
    return V
end

_params(rp, idx, ch) = [_par(rp, i, ch) for i in idx]

function _component_vis(comp::ImageSpec, rp, ch, geom, uv, _wl, _mjd, n, plane_wl)
    # Computed once per uv set and cached: only f and the shift change in a fit. The plane by
    # `plane_wl`, which is `wl` except within a smeared channel: there every sub-wavelength
    # belongs to the channel's own plane.
    return model_image_vis(comp.image, uv; wl = plane_wl, mjd = _mjd, cache = comp.cache,
                           name = comp.comp_name)
end

function _component_vis(comp::VisFuncSpec, rp, ch, geom, uv, _wl, _mjd, n, plane_wl)
    # `$BL` is in Mλ: ρ is cycles/rad, which IS baseline/λ, and the observable plots divide it
    # by 1e6 and label the axis Mλ. A formula copied from a paper is written in those units,
    # so the expression is evaluated in them. `$WL` and `$MJD` are the same per-uv-point
    # values the resolver was given, or NaN when the caller passed none.
    ρ_comp = geom[3]
    return comp.vis_fn(ρ_comp ./ 1e6, _wl, _mjd, _params(rp, comp.param_idx, ch)...)
end

function _component_vis(comp::AnalyticSpec{:crescent}, rp, ch, geom, uv, _wl, _mjd, n, plane_wl)
    # Crescent needs (u,v) for internal offset phase shifts
    return _vis_crescent(_params(rp, comp.param_idx, ch), geom[1], geom[2])
end

function _component_vis(comp::AnalyticSpec{K}, rp, ch, geom, uv, _wl, _mjd, n,
                        plane_wl) where {K}
    return _vis_kernel(Val(K), _params(rp, comp.param_idx, ch), geom[3])
end

# A jump that moves with the parameters, for `_hankel_jac!`: where (in s and in mas), its
# size I(b⁺) − I(b⁻), and ∂b/∂p in s.
const _HankelJump = @NamedTuple{bs::Float64, br::Float64, ΔI::Float64, dbdp::Vector{Float64}}

"""
    _hankel_setup(comp::HankelSpec{T}, B, q; cached=true) -> (s, μ, g, r, cacheable, jumps)

The nodes a profile component is integrated on, in s = r/r_max (`s`, `μ`, weights `g`), and
their radii `r` (mas) at the current size, for the baselines `B` (cycles/mas) and the
parameters `q = [d_min, d_max, profile params...]`:

- what the profile is like (`_hankel_probe`: breakpoints, a singular centre) and how many
  nodes it needs (`_hankel_converge`, from `_hankel_n` for B_max·2r_max lobes) are found
  the first time, and again whenever a uv set is first seen (`cached`); `c,nr` fixes the
  count;
- the nodes and matrices are kept for one size, and dropped when the diameters change;
- breakpoints that move with the parameters (the size among them, for a breakpoint fixed in
  R) are located again when the parameters change, carried as Duals when `q` is, so that
  ForwardDiff differentiates their motion; the nodes are then not cached (`cacheable` is
  false). So is the inner edge s_min = d_min/d_max when the size carries partials; and `r`
  follows the size then too;
- `jumps` lists what moves with the parameters `q[jumpcols]`, for the closed-form Jacobian
  (`_hankel_jac!`), with ∂b/∂q for those columns: the jumps that move with them, and the inner
  edge when they include the diameters and a key gives the size.
"""
function _hankel_setup(comp::HankelSpec{T}, B::AbstractVector, q; cached::Bool = true,
                       jumpcols::UnitRange{Int} = 1:length(q)) where {T}
    cache = comp.kcache
    qv = T[ForwardDiff.value(x) for x in q]
    rm, smin = _hankel_size(comp, qv)
    n0 = comp.nr > 0 ? comp.nr :
         _hankel_n(isempty(B) ? 0.0 : 2 * rm * maximum(b -> abs(ForwardDiff.value(b)), B))
    lay, m = lock(cache.lock) do
        _hankel_resize!(cache, qv)
        l = cache.layout
        if l === nothing || (cached && !_hankel_moving(l) && !_hankel_has_uv(cache, B) &&
                             !_hankel_seen(cache, B))
            _hankel_check!(comp, B, qv, n0, rm, smin)
        end
        cache.layout::_HankelLayout, cache.mult
    end
    n = comp.nr > 0 ? comp.nr : m * n0
    br = lay.breaks
    if _hankel_moving(lay)
        lay, n, br = _hankel_relocate!(comp, B, qv, lay, n, n0, rm, smin)
    end
    jumps = _hankel_jumps(comp, qv, lay, br, rm, smin, jumpcols)
    sized = _hankel_hasgrad(q[1]) || _hankel_hasgrad(q[2])
    if (_hankel_moving(lay) && !(eltype(q) <: AbstractFloat)) || (sized && smin > 0)
        return _hankel_dual_nodes(comp, q, qv, lay, br, n, n0, smin, jumps)
    end
    s, μ, g = lock(cache.lock) do
        get!(() -> _hankel_nodes(T, smin, n; layout = lay, breaks = br, n0),
             cache.nodes, (n, n0))
    end::NTuple{3,Vector{T}}
    r = sized ? (q[2] / 2) .* s : rm .* s
    return (; s, μ, g, r, cacheable = !sized, jumps)
end

# r_max and s_min at the diameters qv[1:2]: the component's own numbers at the size it was
# built with, so that a fixed size gives the nodes it always did.
function _hankel_size(comp::HankelSpec{T}, qv) where {T}
    qv[1] == comp.d_min && qv[2] == comp.d_max && return comp.r_max, comp.s_min
    return T(qv[2] / 2), T(qv[1] / qv[2])
end

# Whether x carries a derivative: a Dual with a partial that is not zero.
_hankel_hasgrad(x::ForwardDiff.Dual) = !iszero(ForwardDiff.partials(x))
_hankel_hasgrad(::Real) = false

# The nodes and matrices are for one size; drop them when the diameters change. Under the lock.
function _hankel_resize!(cache::HankelKernelCache, qv)
    length(cache.dsize) == 2 && cache.dsize[1] == qv[1] && cache.dsize[2] == qv[2] && return
    empty!(cache.nodes); empty!(cache.entries)
    cache.dsize = Float64[qv[1], qv[2]]
    return
end

# Whether the count was checked on the baselines `B`: the uv sets the cache holds matrices
# for, or, after a change of size dropped them, those it was checked on. Under the lock.
_hankel_seen(cache::HankelKernelCache, B::AbstractVector) = hash(B) in cache.uvs

function _hankel_has_uv(cache::HankelKernelCache{T}, B::AbstractVector{TB}) where {T,TB}
    TK = promote_type(TB, T)
    for e in cache.entries
        e isa _HankelUV{TB,TK} || continue
        length(e.B) == length(B) && e.B == B && return true
    end
    return false
end

# The breakpoints at the parameters qv, for a layout whose breakpoints move: located again
# when the parameters have changed since the last call, and the nodes and matrices dropped
# when they have moved. Probed again if one is not found near where it was.
function _hankel_relocate!(comp::HankelSpec{T}, B, qv, lay, n, n0, rm, smin) where {T}
    cache = comp.kcache
    known = lock(cache.lock) do
        cache.last_pp == qv ? copy(cache.breaks) : nothing
    end
    br = known !== nothing ? known : _hankel_breaks_now(lay, comp.profile_fn, qv, 1, smin, T)
    if br === nothing                       # a breakpoint moved too far: probe again
        lay, m = lock(cache.lock) do
            _hankel_check!(comp, B, qv, n0, rm, smin)
            cache.layout::_HankelLayout, cache.mult
        end
        n = comp.nr > 0 ? comp.nr : m * n0
        br = lay.breaks
    end
    lock(cache.lock) do
        if cache.breaks != br               # they moved: the nodes and matrices go with them
            empty!(cache.nodes); empty!(cache.entries)
            cache.breaks = copy(br)
        end
        cache.last_pp = Float64.(qv)
    end
    return lay, n, br
end

# What moves with the parameters qv[cols], for the closed-form Jacobian: each jump that moves
# with one of them, its height from 100 eps(T) either side; and the inner edge
# s_min = d_min/d_max when `cols` holds the diameters and a key gives the size, a jump from 0
# to I(s_min). ∂b/∂q for the columns `cols`.
function _hankel_jumps(comp::HankelSpec{T}, qv, lay::_HankelLayout, br, rm, smin,
                       cols::UnitRange{Int}) where {T}
    jumps = _HankelJump[]
    isempty(cols) && return jumps
    pf = comp.profile_fn
    for i in eachindex(br)
        lay.kinds[i] === :jump && any(x -> x != 0, view(lay.dbdp, i, cols)) || continue
        δ = 100 * eps(T)
        ΔI = _hankel_profile1(pf, qv, 1, T(br[i] + δ)) - _hankel_profile1(pf, qv, 1, T(br[i] - δ))
        push!(jumps, (; bs = br[i], br = Float64(rm) * br[i], ΔI, dbdp = lay.dbdp[i, cols]))
    end
    if first(cols) == 1 && smin > 0 && (comp.dmin_idx != 0 || comp.dmax_idx != 0)
        d = zeros(length(cols)); d[1] = 1 / qv[2]; d[2] = -smin / qv[2]
        push!(jumps, (; bs = Float64(smin), br = Float64(rm * smin),
                       ΔI = _hankel_profile1(pf, qv, 1, T(smin)), dbdp = d))
    end
    return jumps
end

# The nodes as Duals, moving with the parameters: each breakpoint b + Σ ∂b/∂qⱼ (qⱼ − value(qⱼ)),
# and the inner edge d_min/d_max when the size carries partials. Not cached.
function _hankel_dual_nodes(comp::HankelSpec{T}, q, qv, lay::_HankelLayout, br, n, n0, smin,
                            jumps) where {T}
    qc = collect(q)
    bd = [T(br[i]) + sum((T(lay.dbdp[i, j]) * (qc[j] - qv[j]) for j in eachindex(qc));
                         init = zero(eltype(qc))) for i in eachindex(br)]
    sm = smin > 0 ? smin + (qc[1] / qc[2] - ForwardDiff.value(qc[1] / qc[2])) : smin
    s, μ, g = _hankel_nodes(T, sm, n; layout = lay, breaks = bd, n0)
    return (; s, μ, g, r = (qc[2] / 2) .* s, cacheable = false, jumps)
end

# The probe and the count, for these baselines and parameter values. Called under the lock.
function _hankel_check!(comp::HankelSpec{T}, B, qv, n0, rm, smin) where {T}
    cache = comp.kcache
    lay = _hankel_probe(comp.profile_fn, qv, 1, smin, T)
    old = cache.layout
    if old === nothing || old.breaks != lay.breaks || old.sing != lay.sing
        empty!(cache.nodes); empty!(cache.entries)
    end
    cache.layout = lay
    cache.breaks = copy(lay.breaks); cache.last_pp = Float64.(qv)
    if !isempty(B)
        push!(cache.uvs, hash(B))
        length(cache.uvs) > 16 && popfirst!(cache.uvs)
    end
    comp.nr > 0 && return cache
    nodes(n) = _hankel_nodes(T, smin, n; layout = lay, n0)
    m, agree, ok = _hankel_converge(nodes, comp.profile_fn, qv, B .* rm, n0; tol = _hankel_tol(T))
    cache.mult = m
    if !ok && !cache.warned
        cache.warned = true
        @warn "Profile component '$(comp.comp_name)' did not converge: V agrees only to " *
              "$(round(agree; sigdigits = 2)) between $(m * n0 ÷ 2) and $(m * n0) nodes. " *
              "The profile has $(_hankel_describe(lay, rm)). Set '$(comp.comp_name),nr' to " *
              "fix the node count."
    end
    return cache
end

# The profile's parameter vector [d_min, d_max, profile params...] from the resolved values:
# the diameters from their keys, or the component's numbers.
function _hankel_q(comp::HankelSpec, rp, pp)
    dmin = comp.dmin_idx != 0 ? _par(rp, comp.dmin_idx, Val(false)) : comp.d_min
    dmax = comp.dmax_idx != 0 ? _par(rp, comp.dmax_idx, Val(false)) : comp.d_max
    return vcat(dmin, dmax, pp)
end

function _component_vis(comp::HankelSpec, rp, ch, geom, uv, _wl, _mjd, n, plane_wl)
    # A profile's parameters are scalars (`dict_to_model` refuses chromatic ones).
    q = _hankel_q(comp, rp, _params(rp, comp.profile_param_idx, Val(false)))
    oriented = comp.incl_idx != 0 || comp.pa_idx != 0
    st = _hankel_setup(comp, geom[4], q; cached = !oriented, jumpcols = 1:0)
    I_raw = comp.profile_fn(st.s, st.μ, q...)
    # A vector of a known type, whatever the generated profile returned: a constant profile
    # like "1" returns a scalar, and the call itself is not inferred.
    Ti = promote_type(eltype(q), eltype(st.s))
    I = I_raw isa Number ? fill(Ti(I_raw), length(st.s)) : convert(Vector{Ti}, I_raw)
    return _hankel_component_vis(comp, I::Vector{Ti}, st.r, st.g, rp, ch, geom, n,
                                 !oriented && st.cacheable)
end

function _hankel_component_vis(comp::HankelSpec, I, r, g, rp, ch, geom, n, cached::Bool)
    u_comp, v_comp, _, B_comp = geom
    # The data's own baselines: the transform matrices are cached per uv set. With an
    # inclination or position angle B moves with the parameters, and with breakpoints that
    # move with them the nodes do: the transform is then direct.
    if cached
        cache = comp.kcache
        N_h = dot(g, I)
        V0  = (_hkernel(cache, r, g, B_comp, n, false) * I) ./ N_h
        isempty(comp.az_modes) && return V0
        Vi  = complex.(_hkernel(cache, r, g, B_comp, 0, false) * I ./ N_h)
        ψ   = @. atan(v_comp, u_comp)
        for az in comp.az_modes
            amp_val = _par(rp, az.amp_idx, ch)
            phi_rad = _par(rp, az.phi_idx, ch) * π / 180
            Hn = (_hkernel(cache, r, g, B_comp, az.order, false) * I) ./ N_h
            @. Vi += amp_val * (-1im)^az.order * Hn * cos(az.order * (ψ + phi_rad - π/2))
        end
        return Vi
    end
    isempty(comp.az_modes) && return _hankel_wvis(I, r, g, B_comp; n)
    # Azimuthal variations: V = H₀/N + Σ amp·(-j)^n·Hₙ/N·cos(n·(ψ+φ-π/2))
    #
    # PMOIRED's convention, term for term: `_Vazvar` in its oimodels.py computes
    # `amp * (-1j)**n * Hankel(n) * cos(n*(_PA + phi*pi/180 - pi/2))` with
    # `_PA = arctan2(_v, _u)`, and its image-domain brightness is the matching
    # `1 + amp*cos(n*(PA + phi*pi/180 - pi/2))`. A `projang` therefore means the
    # same angle in both packages and a model can move between them unchanged.
    N_h = dot(g, I)
    Vi  = complex.(_hankel_wtransform(I, r, g, B_comp; n=0) ./ N_h)
    ψ   = @. atan(v_comp, u_comp)    # UV position angle
    for az in comp.az_modes
        amp_val = _par(rp, az.amp_idx, ch)
        phi_rad = _par(rp, az.phi_idx, ch) * π / 180
        Hn = _hankel_wtransform(I, r, g, B_comp; n=az.order) ./ N_h
        @. Vi += amp_val * (-1im)^az.order * Hn *
                 cos(az.order * (ψ + phi_rad - π/2))
    end
    return Vi
end

"""
    model_to_vis(model, x, uv; wl=nothing, mjd=nothing, n=0)

Alias for [`eval_model`](@ref).  Compute complex visibilities from a parametric
model, for consistency with the `model_to_obs` / `model_to_chi2` family.
"""
const model_to_vis = eval_model


# ─────────────────────────────────────────────────────────────────────────────
# eval_model_sub  (bandwidth smearing: the model across each channel)
# ─────────────────────────────────────────────────────────────────────────────

"""
    sub_scales(wl, dlam, nsub) -> Matrix{Float64}

How a channel is sampled for bandwidth smearing: `N × nsub` factors by which each point's uv
(its baseline over the channel's central wavelength `wl`) is multiplied at each
sub-wavelength, which is `wl / factor`. The channel `wl ± dlam/2` is cut into `nsub` equal
steps in wavenumber and sampled at their midpoints. A point with no width keeps factor 1.
"""
function sub_scales(λ::AbstractVector, Δ::AbstractVector, ns::Integer)
    T = float(promote_type(eltype(λ), eltype(Δ)))
    S = ones(T, length(λ), ns)
    for k in 1:ns, i in eachindex(λ)
        λi, Δi = T(λ[i]), T(Δ[i])
        isfinite(Δi) && 0 < Δi < 2λi || continue
        σlo, σhi = 1 / (λi + Δi / 2), 1 / (λi - Δi / 2)
        S[i, k] = λi * (σlo + (k - 0.5) / ns * (σhi - σlo))
    end
    return S
end

"""
    eval_model_sub(model, x, uv; wl, dlam, mjd = nothing, nsub = model.smear, n = 0)
        -> Matrix{Complex{T}}

The model's visibilities across each spectral channel: an `N × nsub` matrix whose column `k`
holds every uv point at its channel's `k`-th sub-wavelength.

A channel is `wl ± dlam/2` (metres, per point or one value for all), sampled at the midpoints
of `nsub` equal steps in WAVENUMBER, a flat spectrum across the channel. A point's uv,
which is its baseline over the channel's central wavelength, is rescaled to each
sub-wavelength, and `\$WL` takes the sub-wavelength. An image component keeps the channel's
own plane at every sub-wavelength. A point with no width (`dlam` 0 or not finite) repeats its
central value.

The observables are averaged over the columns — V² as the mean of |V|², a triple product as
the mean of the products — which is what a power-spectrum estimator measures across a channel
and how PMOIRED and CANDID model bandwidth smearing. The χ² functions, `residuals_flat` and
`model_to_obs` do that for a model with `"smear" => nsub` (2 or more) when the data carry
channel widths.
"""
function eval_model_sub(model::FlatModel, x::AbstractVector, uv::AbstractMatrix;
                        wl, dlam, mjd = nothing, nsub::Integer = model.smear, n::Int = 0)
    N  = size(uv, 2)
    ns = max(Int(nsub), 1)
    λ  = wl isa AbstractVector ? wl : fill(wl, N)
    Δ  = dlam isa AbstractVector ? dlam : fill(dlam, N)
    length(λ) == N && length(Δ) == N ||
        throw(DimensionMismatch("eval_model_sub: $N uv points, $(length(λ)) wavelengths, " *
                                "$(length(Δ)) channel widths"))
    uvs, wls, mjds, plane = _sub_points(uv, λ, Δ, mjd, ns)
    V = eval_model(model, x, uvs; wl = wls, mjd = mjds, n, plane_wl = plane)
    return reshape(V, N, ns)
end

# The sub-wavelength points of `eval_model_sub`, all `nsub` of them in one array: point `i` at
# sub-wavelength `k` is column `(k-1)N + i`. With the wavelength each is evaluated at, its
# epoch, and the channel wavelength its image plane is chosen by.
function _sub_points(uv, λ, Δ, mjd, ns)
    N = size(uv, 2)
    S = sub_scales(λ, Δ, ns)
    uvs = similar(uv, size(uv, 1), N * ns)
    wls = Vector{eltype(S)}(undef, N * ns)
    for k in 1:ns, i in 1:N
        j = (k - 1) * N + i
        @views uvs[:, j] .= uv[:, i] .* S[i, k]
        wls[j] = λ[i] / S[i, k]
    end
    mjds = mjd === nothing ? nothing : repeat(mjd isa AbstractVector ? mjd : fill(mjd, N), ns)
    return uvs, wls, mjds, repeat(λ, ns)
end


# ─────────────────────────────────────────────────────────────────────────────
# eval_model_grad  (Jacobian ∂V/∂x)
# ─────────────────────────────────────────────────────────────────────────────

"""
    eval_model_grad(model, x, uv; wl=nothing, mjd=nothing, n=0) -> (V, J)

Evaluate model visibility and its Jacobian w.r.t. `x`.

Returns
-------
V : complex visibility vector (length nB)
J : Jacobian, Matrix{Complex{T}} of shape (nB, length(x)); T follows x and uv

From closed-form derivatives, by [`eval_model_jac`](@ref). Parameters that are not
floating-point numbers (a ForwardDiff Dual, when a Hessian is taken through this function)
and Bessel orders `n ≠ 0` go through [`eval_model_grad_ad`](@ref).
"""
function eval_model_grad(model::FlatModel, x::AbstractVector, uv::AbstractMatrix;
                         n::Int = 0, kw...)
    (n == 0 && eltype(x) <: Union{AbstractFloat, Integer} &&
     eltype(uv) <: Union{AbstractFloat, Integer}) ||
        return eval_model_grad_ad(model, x, uv; n, kw...)
    return eval_model_jac(model, x, uv; kw...)
end

"""
    eval_model_grad_ad(model, x, uv; wl=nothing, mjd=nothing, n=0) -> (V, J)

[`eval_model_grad`](@ref) by ForwardDiff through [`eval_model`](@ref): the reference the
closed-form [`eval_model_jac!`](@ref) is tested against.
"""
function eval_model_grad_ad(model::FlatModel,
                            x::AbstractVector,
                            uv::AbstractMatrix;
                            wl  = nothing,
                            mjd = nothing,
                            n::Int = 0,
                            plane_wl = wl,
                            normalize::Bool = true)

    nB    = size(uv, 2)

    # ForwardDiff works on real arrays, so the real and imaginary parts are stacked.
    # wl and mjd are captured as constants — not differentiated.
    function eval_real(x_)
        V = eval_model(model, x_, uv; wl, mjd, n, plane_wl, normalize)
        return [real.(V); imag.(V)]   # 2nB-vector, eltype follows x and uv
    end

    # Seed the dual numbers in the caller's precision so the Jacobian comes back in it too;
    # Float64 here would silently force the whole gradient path to double.
    J_ri = ForwardDiff.jacobian(eval_real, collect(float(eltype(x)), x))  # (2nB, n_par)
    J_re = J_ri[1:nB, :]
    J_im = J_ri[nB+1:end, :]
    J    = complex.(J_re, J_im)

    V = eval_model(model, x, uv; wl, mjd, n, plane_wl)
    return V, J
end


# ─────────────────────────────────────────────────────────────────────────────
# Convenience: show model summary
# ─────────────────────────────────────────────────────────────────────────────

function Base.show(io::IO, model::FlatModel)
    println(io, "FlatModel with $(length(model.components)) component(s):")
    for comp in model.components
        if comp isa AnalyticSpec
            println(io, "  $(comp.comp_name) [$(comp.kind)]  " *
                    "param_idx=$(comp.param_idx)  f_idx=$(comp.f_idx)")
        elseif comp isa ImageSpec
            mi = comp.image
            println(io, "  $(comp.comp_name) [image]  $(basename(mi.path))  " *
                    "$(join(size(mi.image), "×"))  $(mi.pixsize) mas/px  f_idx=$(comp.f_idx)")
        elseif comp isa HankelSpec
            println(io, "  $(comp.comp_name) [hankel]  " *
                    "nodes=$(comp.nr > 0 ? comp.nr : "auto")  r_max=$(comp.r_max)  " *
                    "profile_param_idx=$(comp.profile_param_idx)  f_idx=$(comp.f_idx)")
        end
    end
    if model.kernel_idx != 0
        println(io, "  spatial_kernel: idx=$(model.kernel_idx)")
    end
    model.smear >= 2 && println(io, "  smear: $(model.smear) sub-wavelengths per channel" *
                                    (model.smear_coherent ? ", coherent" : ""))
    model.smear == -1 && println(io, "  smear: auto, for an extent of " *
                                     "$(round(model.smear_radius; sigdigits = 3)) mas" *
                                     (model.smear_coherent ? ", coherent" : ""))
    println(io, "  list_free_params: $(model.list_free_params)")
    println(io, "  all_names:  $(model.all_names)")
end
