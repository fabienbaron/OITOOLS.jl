# chi2_map.jl
#
# Grid search over two free parameters, and the χ² map it produces.
#
# Load order (within OITOOLS module):
#   include("parse_model.jl")        # FlatModel, parse_model
#   include("chi2_flat.jl")          # chi2_flat, _ndof
#   include("fit_model.jl")          # FitResult
#   include("chi2_map.jl")           # this file
#
# Why a grid at all, when there are four optimisers
# -------------------------------------------------
# A gradient fitter reports one point and a covariance around it. That is the right answer only
# when the χ² surface is a single well. For two correlated angular parameters it frequently is
# not: `demos/example_model_fitting_limb_darkening.jl` fits α Cen A from a seed of 8 mas and
# lands on 8.31, while seeds of 3 and 5 mas walk out to the bound instead — the surface has
# more than one basin and no error bar says so.
#
# A map says so. It is also the only one of the fitters that cannot miss the minimum inside the
# region it covers, because it looks everywhere in that region, and it needs neither a gradient
# nor a starting guess.
#
# The cost is exponential in the number of parameters, so this is deliberately limited to two.
# That is not a placeholder for a general N-D grid: three parameters at 60 points a side is
# 216 000 evaluations for a picture nobody can read.

"""
    Chi2Map

A χ² surface over two free parameters, with everything needed to draw it and to report the
best point on it.

| field | meaning |
|---|---|
| `p1`, `p2` | the two parameter names, as they appear in `list_free_params` |
| `v1`, `v2` | the grid values along each axis |
| `chi2` | `(length(v1), length(v2))` matrix of raw weighted χ² |
| `ndof` | data points, so `chi2 ./ ndof` is the reduced map |
| `x_opt` | the FULL free-parameter vector at the grid minimum |
| `list_free_params` | names for `x_opt`, in fit-vector order |
| `model` | the compiled model the map was computed with |

Parameters not on the two axes are held at their starting values throughout: a map is a slice
through the χ² surface, not a profile likelihood, and reading it as one overstates how well the
two plotted parameters are determined.

See also [`chi2_map`](@ref) and [`delta_chi2_levels`](@ref).
"""
struct Chi2Map
    p1     ::String
    p2     ::String
    v1     ::Vector{Float64}
    v2     ::Vector{Float64}
    chi2   ::Matrix{Float64}
    ndof   ::Int
    x_opt  ::Vector{Float64}
    list_free_params ::Vector{String}
    model  ::FlatModel
end

"""
    delta_chi2_levels

Δχ² for the 68.3, 95.4 and 99.7% confidence regions of **two** jointly estimated parameters:
`2.30`, `6.17`, `11.8`.

These are the two-parameter values, not the one-parameter `1, 4, 9`. Contouring a
two-parameter map at 1σ = 1 draws a region about a third too small, which is the usual way a
χ² map is misread.
"""
const delta_chi2_levels = (2.30, 6.17, 11.83)

function Base.show(io::IO, m::Chi2Map)
    lo = minimum(m.chi2)
    @printf(io, "Chi2Map: %s × %s, %d × %d grid\n", m.p1, m.p2, length(m.v1), length(m.v2))
    @printf(io, "  %s ∈ [%.4g, %.4g],  %s ∈ [%.4g, %.4g]\n",
            m.p1, first(m.v1), last(m.v1), m.p2, first(m.v2), last(m.v2))
    @printf(io, "  best χ² = %.4f  (χ²r = %.4f, ndof = %d)\n", lo, lo / m.ndof, m.ndof)
    for (i, p) in enumerate(m.list_free_params)
        @printf(io, "  %-20s = %10.4f%s\n", p, m.x_opt[i],
                (p == m.p1 || p == m.p2) ? "" : "   (held)")
    end
end

"""
    FitResult(m::Chi2Map)

The best grid point, in the same shape every other fitter returns.

`n_evals` is the grid size, because that is exactly how many χ² evaluations were spent, and
`ret` is `:GRID` — the search cannot fail to converge, so a convergence code would be a
fiction. The value is only ever as good as the grid is fine: the true minimum lies within half
a cell of the point reported.
"""
FitResult(m::Chi2Map) = FitResult(copy(m.x_opt), copy(m.list_free_params),
                                  minimum(m.chi2), minimum(m.chi2) / m.ndof, m.ndof,
                                  length(m.chi2), :GRID, m.model)

"""
    chi2_map(model_dict, list_free_params, data, p1, p2; kwargs...) -> Chi2Map

Evaluate χ² on a grid over the two free parameters `p1` and `p2`, holding every other free
parameter at its starting value.

```julia
m = read_model_file("LD_power.toml")
map = chi2_map(m.model, m.free, data, "star,ldpow", "star,alpha"; m.lb, m.ub, weights = [1.0,0,0])
res = FitResult(map)          # the best grid point, as any other fitter would report it
```

| keyword | default | meaning |
|---|---|---|
| `n1`, `n2` | `60` | grid points along each axis |
| `range1`, `range2` | from `lb`/`ub` | `(lo, hi)` for each axis, overriding the bounds |
| `lb`, `ub` | `Dict()` | bounds; used as the default axis extent |
| `weights` | `[1,1,1,0,0,0,0]` | as everywhere else in model fitting — seven long |
| `vonmises` | `false` | passed through to `chi2_flat` |
| `nB_workspace` | `100` | passed through to `parse_model` |

An axis with no finite bound and no explicit range is an error rather than a guess: a grid
needs an extent, and inventing one would quietly decide the answer.

`p1` and `p2` must both be free. Mapping a fixed parameter would be a perfectly reasonable
thing to want, but `x[i] ↔ list_free_params[i]` is what makes the result a `FitResult` at all,
so free it first.
"""
function chi2_map(model_dict::Dict{String}, list_free_params::AbstractVector{<:AbstractString},
                  data::Union{OIdata, Vector{<:OIdata}},
                  p1::AbstractString, p2::AbstractString;
                  n1::Integer = 60, n2::Integer = 60,
                  range1 = nothing, range2 = nothing,
                  lb = Dict{String,Float64}(), ub = Dict{String,Float64}(),
                  weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                  vonmises::Bool = false,
                  nB_workspace::Int = 100)

    free = String.(list_free_params)
    k1, k2 = String(p1), String(p2)
    k1 == k2 && error("chi2_map needs two different parameters; both are '$k1'")
    i1 = findfirst(==(k1), free)
    i2 = findfirst(==(k2), free)
    i1 === nothing && error("'$k1' is not a free parameter; free it before mapping it")
    i2 === nothing && error("'$k2' is not a free parameter; free it before mapping it")
    (n1 > 1 && n2 > 1) || error("chi2_map needs at least two points per axis, got $n1 × $n2")

    axis(k, r) = begin
        r === nothing || return (Float64(r[1]), Float64(r[2]))
        lo, hi = Float64(get(lb, k, -Inf)), Float64(get(ub, k, Inf))
        (isfinite(lo) && isfinite(hi)) ||
            error("No range for '$k': it has no finite bounds, so pass range1/range2 to say " *
                  "where the grid should run")
        lo < hi || error("Range for '$k' is empty: $lo is not below $hi")
        return (lo, hi)
    end
    lo1, hi1 = axis(k1, range1)
    lo2, hi2 = axis(k2, range2)

    fm = parse_model(model_dict, free; nB_workspace)
    x0 = Float64[model_dict[p] for p in free]

    v1 = collect(range(lo1, hi1; length = Int(n1)))
    v2 = collect(range(lo2, hi2; length = Int(n2)))
    chi2 = Matrix{Float64}(undef, length(v1), length(v2))

    # One scratch vector, reused. Single-threaded on purpose: a Hankel component's workspace
    # buffers live on the compiled model and are shared by every evaluation, so parallel
    # evaluation of one `FlatModel` is a data race rather than a speed-up.
    x = copy(x0)
    for j in eachindex(v2)
        x[i2] = v2[j]
        for i in eachindex(v1)
            x[i1] = v1[i]
            chi2[i, j] = chi2_flat(fm, x, data; weights, verb = false, vonmises)
        end
    end

    best = argmin(chi2)
    x_opt = copy(x0)
    x_opt[i1] = v1[best[1]]
    x_opt[i2] = v2[best[2]]

    return Chi2Map(k1, k2, v1, v2, chi2, _ndof(data, weights), x_opt, free, fm)
end
