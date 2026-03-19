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

const _GEOMETRY_KEYS = ["profile", "diamin", "fwhmin", "crin",
                       "ud", "fwhm", "ldlin", "ldquad", "ldpow", "resolved"]

const _KIND_GEOMETRY = Dict{Symbol,String}(
    :hankel        => "profile",
    :ring          => "diamin",       # uniform ring (analytic)
    :gaussian_ring => "fwhmin",       # bi-Gaussian ring (analytic)
    :crescent      => "crin",         # crescent (analytic, two offset UDs)
    :ud            => "ud",
    :gaussian      => "fwhm",
    :ldlin         => "ldlin",
    :ldquad        => "ldquad",
    :ldpow         => "ldpow",
    :resolved      => "resolved",     # fully resolved (V=0 at B>0), e.g. halo
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
    :ldpow         => ["ldpow", "alpha"],   # diameter, alpha
    :resolved      => [],                   # fully resolved: no geometry params
    :point         => [],
)


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
    # Gaussian: V = exp(-π²σ²ρ²/ln2 * (1/MAS2RAD)²)
    # FWHM → σ = FWHM / (2√(2ln2))
    fwhm = params[1]
    σ    = fwhm / (2 * sqrt(2 * log(2)))
    t    = @. π * σ * ρ / _MAS2RAD_PM
    return @. exp(-t^2 / log(2))   # = exp(-π²σ²ρ²/MAS2RAD²/ln2)
end

function _vis_ldlin(params::AbstractVector, ρ::AbstractVector)
    vis_ldlin(params[1], params[2], ρ)
end

function _vis_ldquad(params::AbstractVector, ρ::AbstractVector)
    vis_ldquad(params[1], params[2], params[3], ρ)
end

function _vis_ldpow(params::AbstractVector, ρ::AbstractVector)
    vis_ldpow(params[1], params[2], ρ)
end

# ForwardDiff specialization: besselj(ν, x) returns NaN when ν is a Dual,
# so we compute the analytic derivatives manually and construct Dual output.
function _vis_ldpow(params::AbstractVector{D}, ρ::AbstractVector) where {T, V, N, D<:ForwardDiff.Dual{T,V,N}}
    θ = ForwardDiff.value(params[1])
    α = ForwardDiff.value(params[2])

    ν  = α/2 + 1
    ζ  = @. π * θ * ρ / MAS2RAD
    safe    = ζ .> 1e-8
    ζ_safe  = @. ifelse(safe, ζ, one(ζ))

    Jν  = @. ifelse(safe, besselj(ν, ζ),   zero(ζ))
    Jν1 = @. ifelse(safe, besselj(ν+1, ζ), zero(ζ))
    Gν1 = gamma(ν+1)

    V_val = @. ifelse(safe, Gν1 * Jν * (ζ_safe/2)^(-ν), one(ζ))

    # ∂V/∂θ per baseline (chain rule through ζ = πθρ/C)
    dV_dζ  = @. ifelse(safe, -Gν1 * Jν1 * 2^ν / ζ_safe^ν, zero(ζ))
    dV_dθ  = @. dV_dζ * π * ρ / MAS2RAD

    # ∂V/∂α per baseline (chain rule through ν = α/2 + 1)
    dJν_dν = map((z, s) -> s ? _dbesselj_dnu(ν, z) : 0.0, ζ, safe)
    ψν1    = digamma(ν + 1)
    dV_dν  = @. ifelse(safe,
        ψν1 * V_val + Gν1 * 2^ν * (dJν_dν - log(ζ_safe) * Jν) / ζ_safe^ν,
        zero(ζ))
    dV_dα  = dV_dν ./ 2

    # Propagate partials: dV[j] = dV_dθ[j] * ∂θ + dV_dα[j] * ∂α
    θ_p = ForwardDiff.partials(params[1])
    α_p = ForwardDiff.partials(params[2])
    return [ForwardDiff.Dual{T}(V_val[j], dV_dθ[j] * θ_p + dV_dα[j] * α_p)
            for j in eachindex(ρ)]
end

function _vis_point(params::AbstractVector, ρ::AbstractVector)
    ones(eltype(ρ), length(ρ))
end

function _vis_resolved(params::AbstractVector, ρ::AbstractVector)
    zeros(eltype(ρ), length(ρ))
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
    σ_in  = fwhmin  / (2 * sqrt(2 * log(2)))
    σ_out = fwhmout / (2 * sqrt(2 * log(2)))
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

    crpa = crprojang * π / 180
    off  = croff * (crout - crin) / 2   # offset distance (mas)

    # Outer disk shifted away from thin side, inner toward thin side
    offXo = -off / 2 * sin(crpa)
    offYo = -off / 2 * cos(crpa)
    offXi =  off / 2 * sin(crpa)
    offYi =  off / 2 * cos(crpa)

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
    # :crescent handled specially in eval_model (needs u,v not just ρ)
    :ud            => _vis_ud,
    :gaussian      => _vis_gaussian,
    :ldlin         => _vis_ldlin,
    :ldquad        => _vis_ldquad,
    :ldpow         => _vis_ldpow,
    :resolved      => _vis_resolved,
    :point         => _vis_point,
)


# ─────────────────────────────────────────────────────────────────────────────
# Component specification structs
# ─────────────────────────────────────────────────────────────────────────────

abstract type AbstractComponentSpec end

"""Analytic component: uniform disk, Gaussian, limb-darkened disk, point, resolved."""
struct AnalyticSpec <: AbstractComponentSpec
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

"""
Hankel component: arbitrary radial profile.

The profile is compiled to an AD-transparent closure at parse_model time.
The r grid and HankelWorkspace are fixed; only the profile param values
change between chi2 evaluations.
"""
struct AzMode
    order    ::Int     # harmonic order (1, 2, 3, …)
    amp_idx  ::Int     # index into resolver output for amplitude
    phi_idx  ::Int     # index into resolver output for phase (degrees)
end

struct HankelSpec <: AbstractComponentSpec
    comp_name        ::String
    f_idx            ::Int
    profile_fn       ::Function         # (r, mu, profile_params...) -> I
    profile_param_idx::Vector{Int}      # indices in resolver output
    r_grid           ::Vector{Float64}
    mu_grid          ::Vector{Float64}  # sqrt(1-(r/r_max)^2)
    workspace        ::HankelWorkspace
    az_modes         ::Vector{AzMode}   # azimuthal variation modes (empty = none)
    incl_idx         ::Int
    pa_idx           ::Int
    x_idx            ::Int
    y_idx            ::Int
end

"""Top-level compiled model."""
struct FlatModel
    resolver   ::Any                           # RGF.Resolver
    components ::Vector{AbstractComponentSpec}
    all_names  ::Vector{String}                # resolver output names (for debug)
    fit_params ::Vector{String}
    kernel_idx ::Int                           # spatial_kernel FWHM index (0 = none)
end


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

"""Extract unique component name prefixes from a flat param dict."""
function _component_names(param_dict::Dict{String})
    names = String[]
    for k in keys(param_dict)
        idx = findfirst(',', k)
        isnothing(idx) && continue
        push!(names, k[1:idx-1])
    end
    return unique(names)
end

"""Identify the component kind from which geometry key is present."""
function _identify_kind(comp_name::String, param_dict::Dict{String})::Symbol
    for gk in _GEOMETRY_KEYS
        haskey(param_dict, "$comp_name,$gk") && return get(
            Dict("profile"=>:hankel, "diamin"=>:ring, "fwhmin"=>:gaussian_ring,
                 "crin"=>:crescent, "ud"=>:ud, "fwhm"=>:gaussian,
                 "ldlin"=>:ldlin, "ldquad"=>:ldquad, "ldpow"=>:ldpow,
                 "resolved"=>:resolved), gk, :unknown)
    end
    # No geometry key — must have f, treat as unresolved point
    haskey(param_dict, "$comp_name,f") && return :point
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
Recursively resolve a param_dict value to a numeric Float64.
If it's already a Number, return it. If it's a string expression containing
\$-references, substitute them recursively and evaluate.
Used at parse time for grid parameters (diamin, diamout, etc.) that may be expressions.
"""
function _resolve_numeric(key::AbstractString, param_dict::Dict{String}, visited::Set{String}=Set{String}())::Float64
    k = String(key)
    haskey(param_dict, k) || error("Key '$k' not found in param_dict")
    v = param_dict[k]
    v isa Number && return Float64(v)
    v isa AbstractString || error("Key '$k' has unsupported type $(typeof(v))")
    k in visited && error("Circular reference detected for '$k'")
    push!(visited, k)
    expr = v
    # Replace $ref with resolved numeric values
    resolved = replace(expr, r"\$([A-Za-z_][A-Za-z0-9_,]*)" => function(m)
        ref = m[2:end]  # strip leading $
        string(_resolve_numeric(ref, param_dict, visited))
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
function _hankel_r_max(comp_name::String, param_dict::Dict{String})::Float64
    udout_key   = "$comp_name,udout"
    diam_key    = "$comp_name,diam"
    diamout_key = "$comp_name,diamout"
    r_max_key   = "$comp_name,r_max"
    if haskey(param_dict, udout_key)
        return _resolve_numeric(udout_key, param_dict) / 2.0
    elseif haskey(param_dict, diamout_key)
        return _resolve_numeric(diamout_key, param_dict) / 2.0
    elseif haskey(param_dict, diam_key)
        return _resolve_numeric(diam_key, param_dict) / 2.0
    elseif haskey(param_dict, r_max_key)
        return _resolve_numeric(r_max_key, param_dict)
    else
        error("Hankel component '$comp_name' requires '$udout_key', '$diam_key', " *
              "'$diamout_key', or '$r_max_key' to define the r grid.")
    end
end

"""Determine Nr for a Hankel component's r grid (default 100)."""
function _hankel_Nr(comp_name::String, param_dict::Dict{String})::Int
    k = "$comp_name,nr"
    haskey(param_dict, k) ? Int(param_dict[k]) : 100
end

"""Determine r_min for a Hankel component's r grid (0 by default, diamin/2 for rings)."""
function _hankel_r_min(comp_name::String, param_dict::Dict{String})::Float64
    diamin_key = "$comp_name,diamin"
    if haskey(param_dict, diamin_key)
        return _resolve_numeric(diamin_key, param_dict) / 2.0
    end
    return 0.0
end

"""
Keys to exclude from the scalar resolver dict.
Profile strings are compiled separately; udout/r_max/nr are grid config,
not runtime params.
"""
const _HANKEL_META_SUFFIXES = ["profile", "r_max", "nr"]

function _resolver_dict(param_dict::Dict{String}, comp_names::Vector{String})
    # Exclude profile strings, Hankel grid-config keys, and resolved markers
    # from the resolver.  udout IS kept — it may be a fit param.
    exclude = Set{String}()
    for cn in comp_names
        push!(exclude, "$cn,profile")
        push!(exclude, "$cn,r_max")
        push!(exclude, "$cn,nr")
        push!(exclude, "$cn,resolved")
    end
    return Dict{String,Any}(k => v for (k, v) in param_dict if k ∉ exclude)
end


# ─────────────────────────────────────────────────────────────────────────────
# Sugar preprocessing
# ─────────────────────────────────────────────────────────────────────────────

"""
    _preprocess_sugar(param_dict, comp_names) -> Dict{String,Any}

Expand shorthand notations to canonical keys:
- `diam + thick` → `diamin + diamout`  (when no `profile` present)
"""
function _preprocess_sugar(param_dict::Dict{String}, comp_names::Vector{String})
    d = copy(param_dict)
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
# parse_model
# ─────────────────────────────────────────────────────────────────────────────

"""
    parse_model(param_dict, fit_params; nB_workspace=100) -> FlatModel

Compile a flat parameter dict into a `FlatModel` ready for evaluation.

Arguments
---------
param_dict   : Dict{String,Any} with "component,param" keys.
               Values may be Float64 (fixed or free) or String (expression).
               Special string key "component,profile" triggers the Hankel pathway.
fit_params   : Vector{String} of keys to optimize (must be numeric in param_dict).
nB_workspace : Pre-allocated baseline grid size for HankelWorkspace buffers.
               Should be ≥ the number of baselines you will evaluate at.

Returns
-------
A FlatModel whose eval_model / eval_model_grad methods can be called in a loop
with zero Dict allocation.

Example
-------
    param = Dict(
        "star,ud"      => 0.8,
        "star,f"       => 0.6,
        "ring,profile" => "exp(-(\\\$R / \\\$scale)^2 / 2)",
        "ring,scale"   => 1.5,
        "ring,udout"   => 12.0,
        "ring,f"       => "1 - \\\$star,f",
    )
    fit_params = ["star,ud", "star,f", "ring,scale"]
    model = parse_model(param, fit_params)
"""
function parse_model(param_dict::Dict{String},
                     fit_params::Vector{String};
                     nB_workspace::Int = 100)::FlatModel

    comp_names = _component_names(param_dict)

    # ── Expand sugar (diam+thick → diamin+diamout) ───────────────────────────
    pd = _preprocess_sugar(param_dict, comp_names)

    # ── Build resolver dict (numeric params only) ────────────────────────────
    res_dict = _resolver_dict(pd, comp_names)

    # ── Build RGF resolver ───────────────────────────────────────────────────
    resolver = RGF.build_resolver(res_dict, fit_params)
    n2i      = _name_to_idx(resolver.all_names)

    # ── Spatial kernel ───────────────────────────────────────────────────────
    kernel_idx = get(n2i, "spatial_kernel", 0)

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

        if kind == :hankel

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

            # Compile the profile closure (AD-transparent)
            profile_fn = compile_profile(profile_expr, prof_refs)

            # Indices of profile params in resolver output
            prof_idx = map(prof_param_names) do k
                haskey(n2i, k) ||
                    error("Profile param '$k' not found in resolver. " *
                          "It must be a numeric param in param_dict.")
                n2i[k]
            end

            # Build r grid (starts at diamin/2 for rings with profile)
            r_min = _hankel_r_min(cn, pd)
            r_max = _hankel_r_max(cn, pd)
            Nr    = _hankel_Nr(cn, pd)
            r     = collect(range(r_min, r_max; length=Nr))
            mu    = @. sqrt(max(0.0, 1.0 - (r / r_max)^2))

            ws = HankelWorkspace(Nr, nB_workspace, length(prof_idx))

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
                cn, f_idx, profile_fn, prof_idx, r, mu, ws, az_modes,
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

            push!(components, AnalyticSpec(
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

            push!(components, AnalyticSpec(
                cn, kind, vis_fn, f_idx, param_idx,
                incl_idx, pa_idx, x_idx, y_idx
            ))
        end
    end

    return FlatModel(resolver, components, resolver.all_names, fit_params, kernel_idx)
end


# ─────────────────────────────────────────────────────────────────────────────
# eval_model  (forward pass, no gradient)
# ─────────────────────────────────────────────────────────────────────────────

"""
    eval_model(model, x, uv; wl=nothing, mjd=nothing, n=0) -> Vector{ComplexF64}

Evaluate the model visibility at the given parameter vector `x` and baselines `uv`.

Arguments
---------
model : FlatModel from parse_model
x     : current free-parameter values, length = length(model.fit_params)
uv    : 2×N matrix of (u,v) spatial frequencies in cycles/rad
wl    : wavelength per UV point (metres), length N  — enables \$WL in expressions
mjd   : MJD per UV point, length N  — enables \$MJD in expressions
n     : Bessel order for Hankel components (default 0)

Returns
-------
V : complex visibility vector, length N
"""
function eval_model(model::FlatModel,
                    x::AbstractVector,
                    uv::AbstractMatrix;
                    wl  = nothing,
                    mjd = nothing,
                    n::Int = 0)

    # Resolve all params → raw vector (no Dict allocation)
    # Elements may be scalars or vectors (when expressions reference $WL/$MJD)
    _wl  = isnothing(wl)  ? NaN : wl
    _mjd = isnothing(mjd) ? NaN : mjd
    pv = model.resolver.fn(x, _wl, _mjd)

    nB = size(uv, 2)
    u_raw = @view uv[1,:]
    v_raw = @view uv[2,:]
    ρ  = @. sqrt(u_raw^2 + v_raw^2)   # cycles/rad
    B  = ρ ./ _MAS2RAD_PM             # cycles/mas (for Hankel)

    # eltype(pv) may be Any for chromatic models (mix of scalars and vectors);
    # always use Float64 for the visibility accumulator — broadcasting handles
    # scalar/vector f and Vi correctly.
    T  = eltype(x) <: Real ? promote_type(eltype(x), Float64) : Float64
    V  = zeros(Complex{T}, nB)

    for comp in model.components
        f = comp.f_idx == 0 ? 1.0 : pv[comp.f_idx]

        # ── UV stretching for inclination / position angle ──────────────
        if comp.incl_idx != 0 || comp.pa_idx != 0
            incl_rad = (comp.incl_idx != 0 ? pv[comp.incl_idx] : 0.0) * π / 180
            pa_rad   = (comp.pa_idx   != 0 ? pv[comp.pa_idx]   : 0.0) * π / 180
            # PMOIRED convention: rot = -PA, compress u' by cos(incl)
            # u' = cos(incl) * ( cos(PA)*u - sin(PA)*v )
            # v' = sin(PA)*u + cos(PA)*v
            u_comp = @. cos(incl_rad) * (cos(pa_rad) * u_raw - sin(pa_rad) * v_raw)
            v_comp = @. sin(pa_rad) * u_raw + cos(pa_rad) * v_raw
            ρ_comp = @. sqrt(u_comp^2 + v_comp^2)
            B_comp = ρ_comp ./ _MAS2RAD_PM
        else
            u_comp = u_raw
            v_comp = v_raw
            ρ_comp = ρ
            B_comp = B
        end

        # ── Compute visibility for this component ──────────────────────
        if comp isa AnalyticSpec && comp.kind == :crescent
            # Crescent needs (u,v) for internal offset phase shifts
            params = pv[comp.param_idx]
            Vi = _vis_crescent(params, u_comp, v_comp)

        elseif comp isa AnalyticSpec
            params = pv[comp.param_idx]          # small Vector, on stack
            Vi = comp.vis_fn(params, ρ_comp)

        elseif comp isa HankelSpec
            prof_params = pv[comp.profile_param_idx]
            I_raw = comp.profile_fn(comp.r_grid, comp.mu_grid, prof_params...)
            # Ensure I is a vector (constant profiles like "1" return a scalar)
            I = I_raw isa Number ? fill(Float64(I_raw), length(comp.r_grid)) : I_raw

            if isempty(comp.az_modes)
                Vi = hankel_vis(I, comp.r_grid, B_comp; n)
            else
                # Azimuthal variations: V = H₀/N + Σ amp·(-j)^n·Hₙ/N·cos(n·(ψ+φ+π/2))
                # The +π/2 (vs -π/2 in PMOIRED's image-domain formula) accounts
                # for the sign flip in the Fourier image↔visibility relationship.
                N_h = hankel_norm(I, comp.r_grid)
                Vi  = complex.(hankel_transform(I, comp.r_grid, B_comp; n=0) ./ N_h)
                ψ   = @. atan(v_comp, u_comp)    # UV position angle
                for az in comp.az_modes
                    amp_val = pv[az.amp_idx]
                    phi_rad = pv[az.phi_idx] * π / 180
                    Hn = hankel_transform(I, comp.r_grid, B_comp; n=az.order) ./ N_h
                    @. Vi += amp_val * (-1im)^az.order * Hn *
                             cos(az.order * (ψ + phi_rad + π/2))
                end
            end
        end

        # ── Position offset → phase shift ──────────────────────────────
        if comp.x_idx != 0 || comp.y_idx != 0
            dx = comp.x_idx != 0 ? pv[comp.x_idx] : 0.0
            dy = comp.y_idx != 0 ? pv[comp.y_idx] : 0.0
            phase = @. -2π * (u_raw * dx + v_raw * dy) / _MAS2RAD_PM
            @. V += f * Vi * cis(phase)
        else
            @. V += f * Vi
        end
    end

    # ── Spatial kernel (Gaussian smoothing in visibility space) ─────────
    if model.kernel_idx != 0
        kern_fwhm = pv[model.kernel_idx]
        σ_kern = kern_fwhm / (2 * sqrt(2 * log(2)))
        t_kern = @. π * σ_kern * ρ / _MAS2RAD_PM
        @. V *= exp(-t_kern^2 / log(2))
    end

    return V
end

"""
    model_to_vis(model, x, uv; wl=nothing, mjd=nothing, n=0)

Alias for [`eval_model`](@ref).  Compute complex visibilities from a parametric
model, for consistency with the `model_to_obs` / `model_to_chi2` family.
"""
const model_to_vis = eval_model


# ─────────────────────────────────────────────────────────────────────────────
# eval_model_grad  (Jacobian ∂V/∂x via ForwardDiff)
# ─────────────────────────────────────────────────────────────────────────────

"""
    eval_model_grad(model, x, uv; wl=nothing, mjd=nothing, n=0) -> (V, J)

Evaluate model visibility and its Jacobian w.r.t. `x`.

Returns
-------
V : complex visibility vector (length nB)
J : Jacobian, Matrix{ComplexF64} of shape (nB, length(x))

Uses ForwardDiff for the resolver and analytic components; for Hankel
components uses the cached-K chain rule via hankel_vis_fwd! + pullback.
This is more efficient than naive ForwardDiff through the full Hankel loop
when nB is large, since Bessel evaluations are shared across parameters.
"""
function eval_model_grad(model::FlatModel,
                          x::AbstractVector,
                          uv::AbstractMatrix;
                          wl  = nothing,
                          mjd = nothing,
                          n::Int = 0)

    nB    = size(uv, 2)
    n_par = length(x)

    # ── Strategy: ForwardDiff through the full eval_model call ──────────────
    # This is correct and simple.  For Hankel components with many params,
    # the hankel_vis_fwd! + pullback strategy (Strategy B in the benchmark)
    # would be faster; that upgrade path is left for a later iteration.
    #
    # We split real and imaginary parts because ForwardDiff works on real
    # arrays; ComplexF64 dual numbers are not directly supported.
    # wl and mjd are captured as constants — not differentiated.

    function eval_real(x_)
        V = eval_model(model, x_, uv; wl, mjd, n)
        return [real.(V); imag.(V)]   # 2nB-vector of Float64
    end

    J_ri = ForwardDiff.jacobian(eval_real, collect(Float64, x))  # (2nB, n_par)
    J_re = J_ri[1:nB, :]
    J_im = J_ri[nB+1:end, :]
    J    = complex.(J_re, J_im)

    V = eval_model(model, x, uv; wl, mjd, n)
    return V, J
end


# ─────────────────────────────────────────────────────────────────────────────
# eval_model_grad_hankel!  (zero-allocation Hankel gradient, for later)
# ─────────────────────────────────────────────────────────────────────────────
# Stub — to be filled when the full chi2 loop is implemented.
# The pattern is:
#   for each HankelSpec:
#     N = hankel_vis_fwd!(ws.V, ws.K, I, r, B)
#     ForwardDiff.jacobian!(ws.dI_dp, p -> profile_fn(r, mu, p...), prof_params)
#     hankel_grad_params!(dL_dp_comp, ȳ_f, ws.V, ws.K, N, ws.w, r, ws.dI_dp)
#     accumulate into full ∂L/∂x via resolver Jacobian


# ─────────────────────────────────────────────────────────────────────────────
# Convenience: show model summary
# ─────────────────────────────────────────────────────────────────────────────

function Base.show(io::IO, model::FlatModel)
    println(io, "FlatModel with $(length(model.components)) component(s):")
    for comp in model.components
        if comp isa AnalyticSpec
            println(io, "  $(comp.comp_name) [$(comp.kind)]  " *
                    "param_idx=$(comp.param_idx)  f_idx=$(comp.f_idx)")
        elseif comp isa HankelSpec
            Nr = length(comp.r_grid)
            println(io, "  $(comp.comp_name) [hankel]  " *
                    "Nr=$Nr  r_max=$(comp.r_grid[end])  " *
                    "profile_param_idx=$(comp.profile_param_idx)  f_idx=$(comp.f_idx)")
        end
    end
    if model.kernel_idx != 0
        println(io, "  spatial_kernel: idx=$(model.kernel_idx)")
    end
    println(io, "  fit_params: $(model.fit_params)")
    println(io, "  all_names:  $(model.all_names)")
end
