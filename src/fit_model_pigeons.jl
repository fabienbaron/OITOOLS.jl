# `fit_model_pigeons`: non-reversible parallel tempering over a model's free parameters, through
# Pigeons.jl. Included by ext/OITOOLSPigeonsExt.jl, which supplies Pigeons; the function is
# declared and documented in src/fit_model_samplers.jl. From ROTIR's `src/fit_pigeons.jl`.

const _LDP   = Pigeons.LogDensityProblems
const _LDPAD = Pigeons.LogDensityProblemsAD
const _Dists = Pigeons.Distributions

"""
    ModelPigeonsTarget{F,G,R,ORD}

The tempered target: the log-posterior on the unconstrained coordinates `z` (`logπ`), its
value and gradient (`logπ_grad`), the reference to anneal from and the starting point. The
reference, the start and the differentiation order travel with the value (fields and a type
parameter) rather than as methods defined during the call: a method defined inside the fit is
in a newer world age than the `pigeons` call that follows it, which would not see it.
"""
struct ModelPigeonsTarget{F,G,R,ORD}
    logπ      :: F
    logπ_grad :: G
    D         :: Int
    reference :: R
    z0        :: Vector{Float64}
end

ModelPigeonsTarget(logπ, logπ_grad, D::Int, reference, z0::Vector{Float64}, order::Int) =
    ModelPigeonsTarget{typeof(logπ),typeof(logπ_grad),typeof(reference),order}(
        logπ, logπ_grad, D, reference, z0)

# Pigeons calls a log potential as a function as well as through LogDensityProblems.
(m::ModelPigeonsTarget)(z) = _LDP.logdensity(m, z)

_LDP.dimension(m::ModelPigeonsTarget) = m.D
# Order 0 under the slice explorer, 1 under MALA.
_LDP.capabilities(::Type{<:ModelPigeonsTarget{F,G,R,ORD}}) where {F,G,R,ORD} =
    _LDP.LogDensityOrder{ORD}()
Pigeons.default_reference(m::ModelPigeonsTarget) = m.reference
Pigeons.initialization(m::ModelPigeonsTarget, ::Random.AbstractRNG, ::Int) = copy(m.z0)
function _LDP.logdensity(m::ModelPigeonsTarget, z)
    v = Float64(m.logπ(z))
    return isfinite(v) ? v : -1e30
end
function _LDP.logdensity_and_gradient(m::ModelPigeonsTarget, z)
    v, g = m.logπ_grad(z)
    (isfinite(v) && all(isfinite, g)) || return -1e30, zeros(length(z))
    return Float64(v), Vector{Float64}(g)
end

# The closed-form gradient reaches Pigeons only through this method: Pigeons routes every
# target through `LogDensityProblemsAD.ADgradient` and would otherwise differentiate
# `logdensity` itself. Both the kind and the target are typed, or the call is ambiguous with
# Pigeons' own `ADgradient(kind, log_potential, replica)` (ROTIR's finding).
_LDPAD.ADgradient(kind::ADTypes.AbstractADType, lp::ModelPigeonsTarget, r::Pigeons.Replica;
                  kwargs...) = Pigeons.BufferedAD(lp, r.recorders.buffers)

function _fit_model_pigeons(model_dict::Dict{String}, list_free_params::Vector{String},
                            data::OIdata;
                            lb = Dict{String,Float64}(), ub = Dict{String,Float64}(),
                            weights = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0], constraints = [],
                            nB_workspace = nothing,
                            n_chains::Integer = 10, n_rounds::Integer = 10,
                            explorer::Symbol = :slice, variational::Bool = false,
                            n_chains_variational::Integer = 0,
                            variational_first_round::Integer = 2,
                            reference_sigma::Real = 3.0, seed = nothing,
                            multithreaded::Bool = Threads.nthreads() > 1, verb::Bool = true)
    explorer in (:slice, :mala) ||
        throw(ArgumentError("explorer is :slice or :mala, not :$explorer"))
    fm, cons, idx, lower, upper, x0 = _sampler_setup(model_dict, list_free_params, data;
                                                     lb, ub, constraints, nB_workspace)
    D = length(x0)
    to_θ, to_z, ljac, dθdz, dljac = _box_transform(lower, upper)
    loglike = _sampler_loglike(fm, data, weights, cons, idx)
    loglike_grad = _sampler_loglike_grad(fm, data, weights, cons, idx)
    # `lo + exp(z)` overflows for a bound on one side only when the explorer steps far out:
    # zero density there, before the model sees an infinite parameter.
    logπ = function (z)
        θ = to_θ(z)
        all(isfinite, θ) || return -Inf
        return Float64(loglike(θ)) + ljac(z)
    end
    logπ_grad = function (z)
        θ = to_θ(z)
        all(isfinite, θ) || return -Inf, zeros(length(z))
        v, g = loglike_grad(θ)
        return v + ljac(z), g .* dθdz(z) .+ dljac(z)
    end
    # Started at the model's values, strictly inside the box (a margin on a finite scale: the
    # box's width, or the value's size when a side is open).
    scale = [isfinite(upper[i] - lower[i]) ? upper[i] - lower[i] : max(abs(x0[i]), 1.0) for i in 1:D]
    θ0 = [clamp(x0[i], isfinite(lower[i]) ? lower[i] + 1e-9 * scale[i] : -Inf,
                isfinite(upper[i]) ? upper[i] - 1e-9 * scale[i] : Inf) for i in 1:D]
    z0 = Vector{Float64}(to_z(θ0))
    reference = Pigeons.DistributionLogPotential(
        _Dists.MvNormal(zeros(D), LinearAlgebra.Diagonal(fill(Float64(reference_sigma)^2, D))))
    target = ModelPigeonsTarget(logπ, logπ_grad, D, reference, z0, explorer === :mala ? 1 : 0)
    isfinite(_LDP.logdensity(target, z0)) ||
        error("fit_model_pigeons: the starting point has zero posterior density")
    expl = explorer === :mala ? Pigeons.AutoMALA(default_autodiff_backend = ADTypes.AutoForwardDiff()) :
                                Pigeons.SliceSampler()
    # The variational reference does nothing until `first_tuning_round`; Pigeons' default (6)
    # would leave a run of 6 rounds bit-identical to a fixed reference.
    gref = Pigeons.GaussianReference(first_tuning_round = Int(variational_first_round))
    vkw = n_chains_variational > 0 ? (variational = gref, n_chains_variational = Int(n_chains_variational)) :
          variational              ? (variational = gref,) : NamedTuple()
    kw = (target = target, n_chains = Int(n_chains), n_rounds = Int(n_rounds), multithreaded,
          explorer = expl, record = [Pigeons.traces, Pigeons.round_trip], show_report = verb,
          vkw...)
    pt = seed === nothing ? Pigeons.pigeons(; kw...) : Pigeons.pigeons(; kw..., seed = Int(seed))

    # Each sample carries the D coordinates and then the log potential.
    Z = collect(Pigeons.get_sample(pt))
    Θ = reduce(vcat, (transpose(to_θ(z[1:D])) for z in Z))
    logls = [Float64(loglike(view(Θ, k, :))) for k in axes(Θ, 1)]
    ss = Pigeons.stepping_stone(pt)
    # The stepping stone is log ∫ L(θ) dθ (the reference is normalised and the Jacobian is in
    # the target); the nested samplers' evidence is log ∫ L(θ) π(θ) dθ with π uniform on the
    # bounds, so the bounds' log volume comes off. An infinite bound has no such volume.
    logv = all(isfinite, lower) && all(isfinite, upper) ? sum(log.(upper .- lower)) : NaN
    logz = ss - logv
    rt = try
        Pigeons.n_round_trips(pt)
    catch
        -1
    end
    if verb
        @printf("Pigeons: %d chains, 2^%d scans, %s explorer, %d samples, log(Z) = %.4f, %d round trips\n",
                n_chains, n_rounds, explorer, size(Θ, 1), logz, rt)
    end
    0 <= rt < 3 && @warn "Pigeons: $(rt) round trip(s) — the chains barely traversed the " *
                         "temperature ladder, so other modes may not have been visited. " *
                         "Raise n_rounds or n_chains."
    return _sampler_result(Θ, logls, fm, list_free_params, data, weights;
                           logz, logzerr = rt > 0 ? 1 / sqrt(rt) : NaN,
                           diagnostics = (round_trips = rt, n_chains = Int(n_chains),
                                          n_rounds = Int(n_rounds), explorer = explorer,
                                          stepping_stone = ss),
                           result = pt, sampler = :pigeons)
end
