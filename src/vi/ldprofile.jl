# ============================================================================
# Error bars on a free-form limb-darkening profile, by MGVI / geoVI.
#
# Latents z = [ξ_θ; ξ_f]. The diameter is log-normal about θ₀, θ = θ₀ exp(s_θ ξ_θ). The profile
# lives on M nodes even in μ, joined by straight lines (the `ldpiecewise` disc), with
#     log I(μ_k) = m_k + f_k,     I(μ_k)/I(1) = exp(m_k − m_M + f_k − f_M),
# m the log of the profile the disc had and f a Gaussian process: f = L ξ_f with K = LLᵀ a
# Matérn-3/2 kernel (`:fixed`), or VarInf's 1-D correlated field on a grid padded to 2M so
# that the limb and the centre are not neighbours through the periodic boundary (`:learned`).
# The transformation is the whitened residual vector of the observables, whose Jacobian is
# the model's closed-form one times ∂(θ, nodes)/∂z.
# ============================================================================

mutable struct LDProfileProblem <: AbstractInferenceProblem
    fm      :: FlatModel
    data    :: OIdata
    weights :: Vector{Float64}
    jkw     :: NamedTuple
    M       :: Int
    θ0      :: Float64
    sθ      :: Float64
    m       :: Vector{Float64}          # log I₀ at the nodes
    L       :: Matrix{Float64}          # :fixed, the kernel's Cholesky factor (M × M)
    mcf     :: Any                      # :learned, the CorrelatedField; else nothing
    nfield  :: Int                      # latents of f
    ndata   :: Int
    zc      :: Vector{Float64}          # the z the cache below was computed at
    r       :: Vector{Float64}
    Jr      :: Matrix{Float64}          # ∂r/∂(θ, nodes)
    JT      :: Matrix{Float64}          # ∂r/∂z
end

# f at the M nodes and ∂f/∂ξ_f.
function _ld_field(p::LDProfileProblem, zf::AbstractVector)
    p.mcf === nothing && return p.L * zf, p.L
    f = vec(p.mcf(zf))[1:p.M]
    Jf = Matrix{Float64}(undef, p.M, p.nfield)
    e = zeros(p.nfield)
    for j in 1:p.nfield
        e[j] = 1.0
        Jf[:, j] = vec(mcf_jvp(p.mcf, zf, e))[1:p.M]
        e[j] = 0.0
    end
    return f, Jf
end

# The model's parameters [θ; I_1 … I_{M-1}] (I_M = 1) at z, and their Jacobian in z.
function _ld_params(p::LDProfileProblem, z::AbstractVector)
    M = p.M
    zf = view(z, 2:length(z))
    f, Jf = _ld_field(p, zf)
    # Held within ±8σ of the prior (θ) and e^±30 (the nodes), where the derivative is zero and
    # the prior pulls back: a Newton step of the samplers can try an absurd diameter, whose ζ
    # lies past the kernel's table and costs a quadrature per uv point.
    ξθ = clamp(z[1], -8.0, 8.0)
    θ = p.θ0 * exp(p.sθ * ξθ)
    x = Vector{Float64}(undef, M)
    Jx = zeros(M, length(z))
    x[1] = θ; Jx[1, 1] = abs(z[1]) < 8 ? θ * p.sθ : 0.0
    for k in 1:(M - 1)
        a = p.m[k] - p.m[M] + f[k] - f[M]
        c = exp(clamp(a, -30.0, 30.0))
        x[k+1] = c
        abs(a) < 30 && @views Jx[k+1, 2:end] .= c .* (Jf[k, :] .- Jf[M, :])
    end
    return x, Jx
end

function _ld_update!(p::LDProfileProblem, z::AbstractVector)
    z == p.zc && return p
    x, Jx = _ld_params(p, z)
    residuals_flat_jac!(p.r, p.Jr, p.fm, x, p.data; weights = p.weights, p.jkw...)
    mul!(p.JT, p.Jr, Jx)
    copyto!(p.zc, z)
    return p
end

VarInf.latent_size(p::LDProfileProblem) = 1 + p.nfield
VarInf.data_size(p::LDProfileProblem) = p.ndata
VarInf.transformation(p::LDProfileProblem, z::AbstractVector) = copy(_ld_update!(p, z).r)
VarInf.right_sqrt_metric(p::LDProfileProblem, z::AbstractVector, v::AbstractVector) =
    _ld_update!(p, z).JT * v
VarInf.left_sqrt_metric(p::LDProfileProblem, z::AbstractVector, w::AbstractVector) =
    _ld_update!(p, z).JT' * w
function VarInf.energy_and_gradient(p::LDProfileProblem, z::AbstractVector)
    _ld_update!(p, z)
    return 0.5 * (sum(abs2, p.r) + sum(abs2, z)), p.JT' * p.r .+ z
end
VarInf.latent_blocks(p::LDProfileProblem) =
    [("diameter", 1:1), ("profile", 2:(1 + p.nfield))]

_matern32(d, ℓ, a) = (s = sqrt(3) * d / ℓ; a^2 * (1 + s) * exp(-s))

# The problem for `component` of `model_dict` on `data`, with what the driver reports from it:
# the model with the disc on the grid, its free parameters, the grid.
function _ld_vi_problem(model_dict::AbstractDict, component::AbstractString, data::OIdata;
                        ngrid::Integer = 24, prior::Symbol = :fixed,
                        length_scale::Real = 0.3, amplitude::Real = 0.3,
                        fluct = (0.3, 0.2), slope = (-4.0, 1.0), diameter_std::Real = 0.1,
                        weights = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    prior in (:fixed, :learned) || throw(ArgumentError("prior is :fixed or :learned, not :$prior"))
    cn = String(component)
    M = Int(ngrid)
    M >= 3 || throw(ArgumentError("ngrid must be at least 3"))
    md = Dict{String,Any}(String(k) => v for (k, v) in model_dict)
    ld_free!(md, String[], cn; n = M, basis = :piecewise)
    dkey = "$cn,ldpiecewise"
    free = vcat([dkey], ["$cn,ld$k" for k in 1:(M - 1)])
    I0 = max.(ld_profile(md, cn, ld_nodes(M)), 1e-3)
    m = log.(I0)
    θ0 = Float64(OITOOLS._resolve_numeric(dkey, md))
    fm = parse_model(md, free; nB_workspace = size(data.uv, 2))
    w = Float64.(OITOOLS._pad_weights(weights))
    x0 = vcat(θ0, I0[1:end-1] ./ I0[end])
    ndata = length(residuals_flat(fm, x0, data; weights = w))
    μ = ld_nodes(M)

    L, mcf, nfield = if prior === :fixed
        K = [_matern32(abs(a - b), length_scale, amplitude) for a in μ, b in μ]
        Matrix(cholesky(Symmetric(K + 1e-10 * amplitude^2 * I)).L), nothing, M
    else
        cfg = CorrFieldConfig(; slope_prior = slope, fluct_prior = fluct)
        g = Axis1DInfo(2M, 1.0 / (M - 1))
        c = CorrelatedField([g], [cfg])
        zeros(0, 0), c, VarInf.latent_size(c)
    end
    nz = 1 + nfield
    p = LDProfileProblem(fm, data, w, OITOOLS._jac_cache_kw(data), M, θ0, Float64(diameter_std),
                         m, L, mcf, nfield, ndata, fill(NaN, nz), zeros(ndata),
                         zeros(ndata, M), zeros(ndata, nz))
    return p, md, free
end

# The posterior's maximum. Visibilities past the first null leave several minima in the
# diameter, and the one the disc's law sits in need not be the profile's: the diameter is
# scanned over ±22 % (±2σ of its prior) with the profile at its prior mean, and the three lowest
# minima taken to the maximum; the lowest of those is returned.
function _ld_vi_map(p::LDProfileProblem, verb::Bool)
    nz = VarInf.latent_size(p)
    zs(ξ) = (z = zeros(nz); z[1] = ξ; z)
    ξs = -2.0:0.05:2.0
    Es = [VarInf.energy_and_gradient(p, zs(ξ))[1] for ξ in ξs]
    minima = [i for i in eachindex(Es) if (i == 1 || Es[i] <= Es[i-1]) &&
                                          (i == length(Es) || Es[i] <= Es[i+1])]
    sort!(minima; by = i -> Es[i])
    zmap, emap = zs(0.0), Inf
    for i in first(minima, 3)
        zi, info = VarInf.reconstruct_map_with_info(p; z0 = zs(ξs[i]), maxiter = 500, verb)
        info.energy < emap && ((zmap, emap) = (zi, info.energy))
    end
    return zmap
end

_ld_chi2r(p::LDProfileProblem, z) = sum(abs2, VarInf.transformation(p, z)) / p.ndata

# VI from `z0`, the posterior's maximum, and `n_draws` antithetic pairs about its centre
# (refined for geoVI): the centre, `:mean` or `:map`, and the samples as offsets from it.
function _ld_vi_samples(p::LDProfileProblem, z0, chi2r_map, method::Symbol, n_iterations,
                        n_samples, n_draws, verb::Bool; warn::Bool = true)
    nz = VarInf.latent_size(p)
    # A few dozen latents: the metric solves are run to convergence, which the drivers'
    # defaults (sized for images) do not do; an inexact sample lands off the posterior.
    z, _ = method === :geovi ?
        reconstruct_geovi(p; z0, n_iterations, n_samples, verb,
                          cg_maxiter = 4nz, cg_tol = 1e-6, kl_absdelta = 1e-3) :
        reconstruct_mgvi(p; z0, n_iterations, n_samples, verb,
                         cg_maxiter = 4nz, cg_tol = 1e-6, kl_absdelta = 1e-3)

    # A variational mean far worse than the posterior's maximum is a run that diverged: the
    # samples are then drawn about the maximum instead.
    centre = :mean
    c = _ld_chi2r(p, z)
    if c > 2 * chi2r_map
        warn && @warn "fit_ld_profile_vi: the variational mean diverged (χ²r " *
                      "$(round(c; digits = 3)) against $(round(chi2r_map; digits = 3)) at " *
                      "the posterior's maximum); the samples are drawn about the maximum"
        z = z0
        centre = :map
    end

    lin, ms, _ = VarInf._draw_samples_geovi(p, z, n_draws; cg_maxiter = 4nz, cg_tol = 1e-6)
    samples = Vector{Vector{Float64}}()
    for k in 1:n_draws
        if method === :geovi
            push!(samples, VarInf._refine_sample_geovi(p, z, lin[k]; metric_sample = ms[k],
                                                       verb = false))
            push!(samples, VarInf._refine_sample_geovi(p, z, -lin[k]; metric_sample = -ms[k],
                                                       verb = false))
        else
            push!(samples, lin[k]); push!(samples, -lin[k])
        end
    end
    return z, centre, samples
end

# The samples as profiles: I/I(1) and I/⟨I⟩ (columns) and diameters, the non-finite left out.
function _ld_vi_profiles(p::LDProfileProblem, z, samples)
    M = p.M
    g0 = OITOOLS._ld_table(:piecewise, M, 0.0).g0       # ∫ b_k μ dμ
    Is = Matrix{Float64}(undef, M, length(samples)); Fs = similar(Is)
    θs = Vector{Float64}(undef, length(samples))
    for (j, s) in enumerate(samples)
        x, _ = _ld_params(p, z .+ s)
        Ij = vcat(x[2:end], 1.0)
        Is[:, j] = Ij
        Fs[:, j] = Ij ./ (2 * dot(g0, Ij))
        θs[j] = x[1]
    end
    keep = [isfinite(θs[j]) && all(isfinite, view(Is, :, j)) for j in eachindex(θs)]
    return Is[:, keep], Fs[:, keep], θs[keep], count(!, keep)
end

# Bagged VI: the same run on `nboot` block-bootstrap replicates of the data (configurations
# drawn with replacement), each from the full data's maximum with its error scale, the
# replicates' samples pooled with equal weight. Replicate `b` is seeded `seed0 + b`, so the
# pool does not depend on the threads.
function _ld_vi_bagged(p::LDProfileProblem, zmap, model_dict, cn, data; nboot, seed0,
                       method, n_iterations, n_samples, n_draws, kw, verb)
    blocks = OITOOLS.data_blocks(data; granularity = :config)
    length(blocks) >= 2 ||
        throw(ArgumentError("bootstrap needs at least 2 configurations, the data have $(length(blocks))"))
    out = Vector{Any}(nothing, nboot)
    next = Threads.Atomic{Int}(1); done = Threads.Atomic{Int}(0)
    function work()
        while (b = Threads.atomic_add!(next, 1)) <= nboot
            try
                Random.seed!(seed0 + b)
                db = OITOOLS.resample_blocks(data, blocks; rng = Xoshiro(seed0 + b))
                pb, _, _ = _ld_vi_problem(model_dict, cn, db; kw...)
                pb.weights .= p.weights
                zb = VarInf.reconstruct_map(pb; z0 = zmap, maxiter = 500, verb = false)
                z, centre, s = _ld_vi_samples(pb, zb, _ld_chi2r(pb, zb), method, n_iterations,
                                              n_samples, n_draws, false; warn = false)
                Is, Fs, θs, _ = _ld_vi_profiles(pb, z, s)
                out[b] = (Is, Fs, θs, centre)
            catch err
                err isa InterruptException && rethrow()
            end
            n = Threads.atomic_add!(done, 1) + 1
            verb && @printf("\rbootstrap: %d/%d replicates", n, nboot)
        end
    end
    nw = min(Threads.nthreads(), nboot)
    nw == 1 ? work() : foreach(wait, [Threads.@spawn(work()) for _ in 1:nw])
    verb && println()
    ok = [o for o in out if o !== nothing && length(o[3]) >= 2]
    length(ok) >= 2 || error("fit_ld_profile_vi: fewer than 2 bootstrap replicates ran")
    length(ok) < nboot && @warn "fit_ld_profile_vi: $(nboot - length(ok)) of $nboot bootstrap replicates failed and are left out"
    nmap = count(o -> o[4] === :map, ok)
    nmap > 0 && @warn "fit_ld_profile_vi: $nmap of $(length(ok)) bootstrap replicates diverged and are drawn about their maximum"
    # equal weight: each replicate brings as many samples as the smallest
    n = minimum(length(o[3]) for o in ok)
    return reduce(hcat, [o[1][:, 1:n] for o in ok]), reduce(hcat, [o[2][:, 1:n] for o in ok]),
           reduce(vcat, [o[3][1:n] for o in ok]), length(ok)
end

function OITOOLS.fit_ld_profile_vi(model_dict::AbstractDict, component::AbstractString,
                                   data::OIdata;
                                   ngrid::Integer = 24, prior::Symbol = :fixed,
                                   length_scale::Real = 0.3, amplitude::Real = 0.3,
                                   fluct = (0.3, 0.2), slope = (-4.0, 1.0),
                                   diameter_std::Real = 0.1, method::Symbol = :geovi,
                                   n_iterations::Integer = 6, n_samples::Integer = 4,
                                   n_draws::Integer = 50,
                                   weights = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                                   scale_errors::Bool = false, bootstrap::Integer = 0,
                                   seed = nothing, verb::Bool = false)
    method in (:geovi, :mgvi) || throw(ArgumentError("method is :geovi or :mgvi, not :$method"))
    (bootstrap == 0 || bootstrap >= 2) ||
        throw(ArgumentError("bootstrap is 0 (off) or at least 2 replicates, not $bootstrap"))
    t0 = time()
    seed === nothing || Random.seed!(seed)
    seed0 = seed === nothing ? UInt64(rand(UInt32)) : UInt64(seed)
    cn = String(component)
    kw = (; ngrid, prior, length_scale, amplitude, fluct, slope, diameter_std, weights)
    p, md, free = _ld_vi_problem(model_dict, cn, data; kw...)
    M = p.M; μ = ld_nodes(M); dkey = free[1]; w = p.weights; fm = p.fm; ndata = p.ndata

    zmap = _ld_vi_map(p, verb)
    chi2r_map = _ld_chi2r(p, zmap)
    # Error bars that are too small (χ²r > 1 at the free-form profile's best fit) make the
    # profile follow the calibration's systematics with a narrow band: scaled by √χ²r there.
    error_scale = 1.0
    if scale_errors && chi2r_map > 1
        error_scale = sqrt(chi2r_map)
        p.weights ./= chi2r_map
        fill!(p.zc, NaN)
        zmap = VarInf.reconstruct_map(p; z0 = zmap, maxiter = 500, verb)
        chi2r_map = _ld_chi2r(p, zmap)
    end

    z, centre, samples = _ld_vi_samples(p, zmap, chi2r_map, method, n_iterations, n_samples,
                                        n_draws, verb)
    Is, Fs, θs, nbad = _ld_vi_profiles(p, z, samples)
    length(θs) >= 2 || error("fit_ld_profile_vi: the posterior samples are not finite; " *
                             "try prior = :fixed or method = :geovi")
    nbad > 0 && @warn "fit_ld_profile_vi: $nbad of $(length(samples)) samples are not finite and are left out"
    # The model at the variational centre, a point of the posterior. The samples' mean is not
    # one: θ and the limb trade along a curved valley, and its average lies off it.
    xc, _ = _ld_params(p, z)
    md[dkey] = xc[1]
    for k in 1:(M - 1); md["$cn,ld$k"] = xc[k+1]; end
    chi2r = sum(abs2, residuals_flat(fm, xc, data; weights = w)) / ndata

    stats(A) = (vec(mean(A; dims = 2)), vec(std(A; dims = 2)),
                [quantile(A[i, :], 0.16) for i in axes(A, 1)],
                [quantile(A[i, :], 0.84) for i in axes(A, 1)])
    result(Is, Fs, θs, nboot, single) =
        LDProfileVIResult(cn, method, prior, μ, stats(Is)..., Is, stats(Fs)..., Fs,
                          mean(θs), std(θs), θs, chi2r, chi2r_map, error_scale, centre, ndata,
                          md, free, nboot, single, time() - t0)
    bootstrap == 0 && return result(Is, Fs, θs, 0, nothing)

    single = result(Is, Fs, θs, 0, nothing)
    Ib, Fb, θb, nok = _ld_vi_bagged(p, zmap, model_dict, cn, data; nboot = Int(bootstrap),
                                    seed0, method, n_iterations, n_samples,
                                    n_draws = max(4, cld(n_draws, bootstrap)), kw, verb)
    return result(Ib, Fb, θb, nok, single)
end
