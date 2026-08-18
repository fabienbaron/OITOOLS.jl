"""
    OITOOLSPigeonsExt

Parallel tempering for SQUEEZE reconstructions, on Pigeons.jl.
Install at ext/OITOOLSPigeonsExt.jl. Loaded automatically by `using OITOOLS, Pigeons`.

Two facts about Pigeons shape the design:
  * A swap exchanges `replica.chain`, NOT the state (`swap.jl:_swap!`), so a
    SqueezePTState stays with its replica all run.  Its preallocated buffers are
    never copied, and each replica accumulates the β=1 samples it visits —
    summing those over replicas gives exactly the target-chain sample set, with
    no custom recorder needed.
  * The reference must be i.i.d.-samplable for PT's guarantees.  Ours is the
    uniform distribution over element configurations (times a uniform box over
    free SPARCO parameters), sampled exactly.
"""
module OITOOLSPigeonsExt

using OITOOLS
using Pigeons
using Random
using Printf

using OITOOLS: SqueezeTarget, SqueezePTState, SqueezeUniformReference,
               SqueezeExplorer, SqueezeState, SqueezeRegs,
               squeeze_logdensity, squeeze_refresh!, squeeze_step!,
               init_elements!, _parse_squeeze_regularizers, _squeeze_operator,
               _weights_tuple, _cvis_to_chi2, squeeze_ndof, default_nelements,
               resync_raw_im_vis!, form_mod_vis!, lprior, sparco_vis!,
               SqueezeSparco, default_sparco_bounds, _in_box,
               SqueezeMonitor, monitor_update!

# ── Log potentials ───────────────────────────────────────────────────────────
(target::SqueezeTarget)(st::SqueezePTState) = squeeze_logdensity(st)

# Constant inside the support, -Inf outside.  With a SPARCO model the support is
# the parameter box, which is what makes the reference a *proper* uniform
# distribution and hence exactly i.i.d.-samplable.
function (ref::SqueezeUniformReference)(st::SqueezePTState)
    st.model === nothing && return 0.0
    return _in_box(st.model.params, ref.target.param_lo, ref.target.param_hi,
                   st.model.free) ? 0.0 : -Inf
end

# ── Target interface ─────────────────────────────────────────────────────────

function Pigeons.initialization(target::SqueezeTarget{T}, rng::AbstractRNG,
                                replica_index::Int) where {T}
    s = SqueezeState(T, target.nx, target.nelements, target.data.nuv)
    init_elements!(s, rng; x_start = target.x_start)
    regs = _parse_squeeze_regularizers(target.regularizers;
                                       fov = target.fov, cent_mult = target.cent_mult,
                                       tv_eps = target.tv_eps,
                                       prior_penalty = target.prior_penalty)
    # Each replica gets its own model copy: params/stepsize are mutated in place.
    m = target.model === nothing ? nothing :
        SqueezeSparco(copy(target.model.params), copy(target.model.stepsize),
                      copy(target.model.free), copy(target.model.prob_pmovement))
    st = SqueezePTState{T}(s, regs, m, zeros(Int64, target.nx, target.nx),
                           zeros(Float64, OITOOLS.NPARAMS_SPARCO), 0, 0)
    squeeze_refresh!(st, target)
    return st
end

Pigeons.default_reference(target::SqueezeTarget) = SqueezeUniformReference(target)
Pigeons.default_explorer(::SqueezeTarget) = SqueezeExplorer()

function Pigeons.sample_iid!(ref::SqueezeUniformReference, replica, shared)
    st = replica.state
    init_elements!(st.s, replica.rng; x_start = :random)
    if st.model !== nothing
        for j in st.model.free
            lo = ref.target.param_lo[j]; hi = ref.target.param_hi[j]
            st.model.params[j] = lo + (hi - lo) * rand(replica.rng)
        end
    end
    squeeze_refresh!(st, ref.target)
    return st
end

# ── Explorer interface ───────────────────────────────────────────────────────

function Pigeons.step!(explorer::SqueezeExplorer, replica, shared)
    lp = Pigeons.find_log_potential(replica, shared.tempering, shared)
    beta = _beta_of(lp)
    st = replica.state
    target = _target_of(lp)

    squeeze_step!(st, target, explorer, beta, replica.rng)

    # Pigeons' early rounds are TUNING rounds: the schedule is still adapting, so
    # their samples are burn-in.  Reset the accumulators when the round changes,
    # leaving only the final (fully tuned) round.  Averaging over all rounds gave
    # chi2r = 151 instead of 0.78.
    round = shared.iterators.round
    if round != st.last_round
        fill!(st.mean_accum, 0)
        fill!(st.param_accum, 0.0)
        st.nsamples = 0
        st.last_round = round
    end

    if beta == 1.0
        h = st.s.histogram
        @inbounds for p in eachindex(h)
            st.mean_accum[p] += h[p]
        end
        st.model === nothing || (st.param_accum .+= st.model.params)
        st.nsamples += 1

        mon = explorer.monitor
        if mon !== nothing
            # Pigeons restarts `scan` at 1 each round (round r does 2^r scans),
            # so plotting against it makes traces double back.  Rounds 1..r-1
            # contribute 2^1 + ... + 2^(r-1) = 2^r - 2 scans.
            round_i = shared.iterators.round
            cumulative = (1 << round_i) - 2 + shared.iterators.scan
            if cumulative % mon.every == 0
                monitor_update!(mon, Float64.(st.s.histogram) ./ length(st.s.element_x),
                                cumulative;
                                chi2r = 2 * st.s.lLikelihood / mon.ndf,
                                temperature = round_i,
                                params = st.model === nothing ? nothing : st.model.params)
            end
        end
    end
    return st
end

Pigeons.explorer_recorder_builders(::SqueezeExplorer) = []

# The path is InterpolatingPath(reference, target); pull β and the target out.
_beta_of(lp) = lp.beta
_target_of(lp) = lp.path.target

# ── Driver ───────────────────────────────────────────────────────────────────

function OITOOLS.reconstruct_squeeze_tempered(
        x_start, data::OIdata{T}, ft;
        pixsize::Real = -1.0,
        nelements::Integer = 0,
        n_rounds::Integer = 10,
        n_chains::Integer = 10,
        regularizers = [],
        weights = [1.0, 1.0, 1.0],
        vonmises::Bool = false,
        prior_image = nothing,
        fov::Real = 1.0, cent_mult::Real = 1.0, tv_eps::Real = 0.0,
        auto_centering::Bool = true,
        model::Union{Nothing,SqueezeSparco} = nothing,
        param_bounds = nothing,
        n_passes::Integer = 1,
        multithreaded::Bool = true,
        monitor::Integer = 0,
        monitor_colormap::AbstractString = "gist_heat",
        outfile::AbstractString = "",
        seed::Integer = 1,
        verb::Bool = true,
        pigeons_kwargs...) where {T<:AbstractFloat}

    A, nx, pixsize = _squeeze_operator(ft, data, pixsize)
    nelem = nelements > 0 ? Int(nelements) : default_nelements(data, nx)
    w = _weights_tuple(weights)

    prior_penalty = Float64[]
    if prior_image !== nothing
        size(prior_image) == (nx, nx) || throw(DimensionMismatch(
            "prior_image is $(size(prior_image)) but `ft` describes an $(nx)x$(nx) image"))
        prior_penalty = [p > 0 ? -log(Float64(p)) : 1e12 for p in vec(prior_image)]
    end

    regs_have_centering = any(r -> lowercase(String(r[1])) == "centering", regularizers)
    effective = collect(Any, regularizers)
    # As in the annealed path: a parametric model already pins the position.
    auto_centering = auto_centering && model === nothing
    if auto_centering && !regs_have_centering && prior_image === nothing
        push!(effective, Any["centering", 1.0])
    end
    if prior_image !== nothing &&
       !any(r -> lowercase(String(r[1])) == "priorimage", effective)
        push!(effective, Any["priorimage", 1.0])
    end

    centering_on = any(r -> lowercase(String(r[1])) == "centering" && r[2] > 0, effective)
    nparams_free = model === nothing ? 0 : length(model.free)
    ndf = Float64(squeeze_ndof(data, weights; nparams_free = nparams_free,
                               centering = centering_on))

    lo, hi = if model === nothing
        Float64[], Float64[]
    elseif param_bounds === nothing
        default_sparco_bounds(model)
    else
        length(param_bounds) == 2 || throw(ArgumentError(
            "param_bounds must be (lo, hi), each with 6 entries"))
        collect(Float64, param_bounds[1]), collect(Float64, param_bounds[2])
    end
    if model !== nothing
        for j in model.free
            lo[j] < hi[j] || throw(ArgumentError(
                "parameter $(j) is free but its reference box is empty " *
                "(lo=$(lo[j]), hi=$(hi[j])); widen `param_bounds`"))
            lo[j] <= model.params[j] <= hi[j] || throw(ArgumentError(
                "initial value $(model.params[j]) for parameter $(j) is outside " *
                "its reference box [$(lo[j]), $(hi[j])]"))
        end
    end

    target = SqueezeTarget{T}(data, A, nx, nelem, w, vonmises, effective,
                              prior_penalty, Float64(fov), Float64(cent_mult),
                              Float64(tv_eps), x_start, model, lo, hi)

    if verb
        @printf("SQUEEZE (tempered): nx=%d  nuv=%d  nelements=%d  ndf=%.0f\n",
                nx, data.nuv, nelem, ndf)
        @printf("                    %d chains x %d rounds (2^%d scans), T=%s, A=%.1f MB\n",
                n_chains, n_rounds, n_rounds, T, sizeof(A) / 1e6)
    end

    mon = nothing
    if monitor > 0
        mon = SqueezeMonitor(Int(monitor); pixsize = pixsize,
                             title = "SQUEEZE (tempered)",
                             colormap = monitor_colormap,
                             free = model === nothing ? Int[] : copy(model.free),
                             trace2_label = "round")
        mon.ndf = ndf
        if multithreaded
            @info "reconstruct_squeeze_tempered: monitoring enabled, disabling " *
                  "multithreading — the live display has to draw from the main thread, " *
                  "and a Python call from a worker thread segfaults. Set monitor=0 to " *
                  "restore the threaded path."
            multithreaded = false
        end
    end

    pt = pigeons(; target = target,
                   reference = SqueezeUniformReference(target),
                   explorer = SqueezeExplorer(; n_passes = Int(n_passes), monitor = mon),
                   n_rounds = Int(n_rounds),
                   n_chains = Int(n_chains),
                   multithreaded = multithreaded,
                   seed = Int(seed),
                   show_report = verb,
                   pigeons_kwargs...)

    states = [r.state for r in pt.replicas]
    total = zeros(Int64, nx, nx)
    ptot = zeros(Float64, OITOOLS.NPARAMS_SPARCO)
    nsamples = 0
    for st in states
        total .+= st.mean_accum
        ptot .+= st.param_accum
        nsamples += st.nsamples
    end
    nsamples > 0 || error("reconstruct_squeeze_tempered: no beta=1 samples were " *
                          "recorded; increase n_rounds")

    img = Float64.(total) ./ (nsamples * nelem)

    # Posterior-mean parameters: the average over beta = 1 visits ONLY.  Averaging
    # over replicas would mix in the cold chains, whose parameters are essentially
    # draws from the reference box.
    params_mean = model === nothing ? Float64[] : ptot ./ nsamples

    V = image_to_vis(T.(img), A)
    if model !== nothing
        # Score the mean image WITH its fitted model, or the star and background
        # are simply missing from the comparison (this gave chi2r 1351 vs 0.11).
        pv = zeros(Complex{T}, data.nuv); fr = ones(T, data.nuv)
        sparco_vis!(pv, fr, params_mean, data)
        @inbounds @simd for k in eachindex(V)
            V[k] = pv[k] + V[k] * fr[k]
        end
    end
    chi2r_mean = _cvis_to_chi2(V, data, w, vonmises) / ndf
    lz = try
        Pigeons.stepping_stone(pt)
    catch
        NaN
    end

    if verb
        @printf("  chi2r(mean image) = %.4f   logZ = %.4f   beta=1 samples = %d\n",
                chi2r_mean, lz, nsamples)
        @printf("  global barrier    = %.4f  (round trips are ~ n_scans / barrier^2)\n",
                Pigeons.global_barrier(pt))
    end

    isempty(outfile) || writefits(img, outfile; pixsize = pixsize)

    diagnostics = (chi2r_mean = chi2r_mean,
                   logZ = lz,
                   global_barrier = Pigeons.global_barrier(pt),
                   nsamples = nsamples,
                   ndf = ndf,
                   nelements = nelem,
                   schedule = pt.shared.tempering.schedule.grids,
                   acceptance = [st.s.naccepted / max(st.s.nproposed, 1) for st in states],
                   params = params_mean,
                   params_by_replica = model === nothing ? Vector{Float64}[] :
                                       [copy(st.model.params) for st in states],
                   param_names = OITOOLS.SPARCO_PARAM_NAMES,
                   pt = pt)

    return img, diagnostics
end

OITOOLS.reconstruct_squeeze_tempered(data::OIdata, ft; kwargs...) =
    OITOOLS.reconstruct_squeeze_tempered(:random, data, ft; kwargs...)

OITOOLS.reconstruct_squeeze_tempered(data::AbstractArray{<:OIdata}, ft; kwargs...) =
    OITOOLS.reconstruct_squeeze_tempered(:random, first(data), ft; kwargs...)

function OITOOLS.reconstruct_squeeze_tempered(x_start, data::AbstractArray{<:OIdata}, ft;
                                              kwargs...)
    length(data) == 1 || throw(ArgumentError(
        "reconstruct_squeeze_tempered: polychromatic reconstruction is not implemented"))
    return OITOOLS.reconstruct_squeeze_tempered(x_start, first(data), ft; kwargs...)
end

end # module
