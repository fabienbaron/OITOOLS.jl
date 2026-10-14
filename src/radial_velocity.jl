# Radial velocities fitted with the interferometry: a table of measurements, and their residuals
# against the orbit of the model's component on an orbit (`c,x = c,y = "orbit"`), with the
# velocity elements of PMOIRED's `_orbit` (`c,orb Ka`, `Kb`, `gamma`, or `M` + `q` with `incl`).

"""
    RVData

Radial velocities: `mjd`, `rv` and `sigma` (km/s), and `star` — `:a` (the primary, at the
origin) or `:b` (the secondary, the component on an orbit) — one entry per measurement.
"""
struct RVData
    mjd   :: Vector{Float64}
    rv    :: Vector{Float64}
    sigma :: Vector{Float64}
    star  :: Vector{Symbol}
    path  :: String
end

Base.length(r::RVData) = length(r.mjd)

"""
    read_rv(path) -> RVData

A CSV table of radial velocities with the header `MJD,RV,sigma,star`: time (MJD), velocity and
its error (km/s, positive receding), and the star, `a` (primary) or `b` (secondary). Lines
starting with `#` are comments.
"""
function read_rv(path::AbstractString)
    lines = [strip(l) for l in eachline(path) if !isempty(strip(l)) && !startswith(strip(l), "#")]
    isempty(lines) && throw(ArgumentError("$(basename(path)): no radial velocities"))
    head = lowercase.(strip.(split(lines[1], ',')))
    function col(n)
        i = findfirst(==(n), head)
        i === nothing && throw(ArgumentError("$(basename(path)): no column \"$n\" in the " *
                                             "header (MJD,RV,sigma,star)"))
        return i
    end
    im, ir, is, it = col("mjd"), col("rv"), col("sigma"), col("star")
    mjd, rv, σ, star = Float64[], Float64[], Float64[], Symbol[]
    for (k, l) in enumerate(lines[2:end])
        f = strip.(split(l, ','))
        s = Symbol(lowercase(f[it]))
        s in (:a, :b) || throw(ArgumentError("$(basename(path)), row $k: star is a or b, not $(f[it])"))
        push!(mjd, parse(Float64, f[im])); push!(rv, parse(Float64, f[ir]))
        push!(σ, parse(Float64, f[is])); push!(star, s)
    end
    all(>(0), σ) || throw(ArgumentError("$(basename(path)): every sigma must be positive"))
    return RVData(mjd, rv, σ, star, String(path))
end

# The component on an orbit (its x is "orbit", or the expression `_expand_orbit!` writes).
function _orbit_component(model_dict::AbstractDict)
    cs = [cn for cn in _component_names(Dict{String,Any}(model_dict))
          if (v = get(model_dict, "$cn,x", nothing); v isa AbstractString &&
              (v == "orbit" || startswith(v, "orbit_x(")))]
    isempty(cs) && throw(ArgumentError("radial velocities need a component on an orbit " *
                                       "(c,x = c,y = \"orbit\")"))
    length(cs) == 1 || throw(ArgumentError("radial velocities with several components on " *
                                           "orbits ($(join(cs, ", "))) are not supported"))
    return cs[1]
end

# Where the RV model's inputs sit in the resolver's output, for component `cn`.
struct RVIndex
    P::Int; MJD0::Int; e::Int; omega::Int; gamma::Int
    Ka::Int; Kb::Int                      # 0 when given by M and q
    M::Int; q::Int; incl::Int
end

function _rv_index(fm, cn::AbstractString)
    n2i = Dict(n => i for (i, n) in enumerate(fm.all_names))
    g(k) = get(n2i, "$cn,orb $k", 0)
    need(k) = (i = g(k); i == 0 ? throw(ArgumentError(
        "radial velocities need $cn,orb $k")) : i)
    ix = RVIndex(need("P"), need("MJD0"), need("e"), need("omega"),
                 g("gamma"), g("Ka"), g("Kb"), g("M"), g("q"), g("incl"))
    (ix.Ka != 0 && ix.Kb != 0) || (ix.M != 0 && ix.q != 0 && ix.incl != 0) ||
        throw(ArgumentError("radial velocities need $cn,orb Ka and Kb, or $cn,orb M and q"))
    return ix
end

"""
    rv_residuals(fm, x, rv, index) -> Vector

(RV − model)/σ for each measurement, the model at the free parameters `x`. Generic in the
number type, so ForwardDiff differentiates it through the resolver.
"""
function rv_residuals(fm, x, rv::RVData, ix::RVIndex)
    v = fm.resolver.fn(x)
    P, MJD0, e, ω = v[ix.P], v[ix.MJD0], v[ix.e], v[ix.omega]
    γ = ix.gamma == 0 ? zero(P) : v[ix.gamma]
    Ka, Kb = ix.Ka != 0 ? (v[ix.Ka], v[ix.Kb]) :
             orbit_semiamplitudes(P, e, v[ix.incl], v[ix.M], v[ix.q])
    return [(rv.rv[k] - orbit_rv(rv.mjd[k], rv.star[k], P, MJD0, e, ω,
                                 rv.star[k] === :a ? Ka : Kb, γ)) / rv.sigma[k]
            for k in eachindex(rv.mjd)]
end

rv_jacobian(fm, x, rv::RVData, ix::RVIndex) =
    ForwardDiff.jacobian(p -> rv_residuals(fm, p, rv, ix), x)

"A bootstrap draw of the measurements: as many, with replacement."
resample_rv(rv::RVData, rng) =
    (k = rand(rng, 1:length(rv), length(rv));
     RVData(rv.mjd[k], rv.rv[k], rv.sigma[k], rv.star[k], rv.path))
