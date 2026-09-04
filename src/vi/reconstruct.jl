# ============================================================================
# reconstruct_*: MAP, MGVI, geoVI and the hybrid schedule.
# ============================================================================

# The generic reconstruct_map / reconstruct_mgvi entry points live in
# VarInf.jl. The SkyModelParams-typed methods below are thin domain wrappers:
# they build an InterferometricProblem, delegate the algorithm to VarInf, and
# add interferometry-specific reporting + posterior image computation.

"""
    reconstruct_map(p, ft, data; kwargs...) -> (xi_opt, image_opt)

MAP reconstruction: minimize posterior energy via L-BFGS.
Returns optimized latent vector and corresponding image cube.
"""
function reconstruct_map(p::SkyModelParams, ft, data;
                         weights=[1.0, 1.0, 1.0],
                         maxiter=200, verb=true,
                         xi0=nothing)
    prob = _build_interferometric_problem(p, ft, data; weights=weights, verb=verb)
    z_opt = reconstruct_map(prob; z0=xi0, maxiter=maxiter, verb=verb)
    image_opt = sky_forward(z_opt[1:_latent_size(p)], p)
    if verb
        report_chi2(image_opt, p, ft, data)
        report_latents(z_opt, p)
    end
    return z_opt, image_opt
end


"""
    reconstruct_mgvi(p, ft, data; kwargs...) -> (xi_center, image_mean, image_std, samples)

Full MGVI reconstruction with posterior sampling.

Algorithm per iteration:
1. Warm-start from MAP (first iteration only)
2. Build implicit Hessian H at current center via finite-diff gradient
3. Draw residual samples δ_k by solving H·δ_k = n_k (n_k ~ N(0,I)) via CG
4. Minimize sample-averaged energy (KL divergence) w.r.t. center using
   antithetic pairs: E_KL(m) = (1/2K) Σ_k [E(m+δ_k) + E(m-δ_k)]
5. Repeat with updated Hessian

Returns: optimized center, posterior mean image, posterior std image, residual samples.
"""
function reconstruct_mgvi(p::SkyModelParams, ft, data;
                          weights=[1.0, 1.0, 1.0],
                          n_iterations=6,
                          n_samples=3,
                          map_maxiter=100,
                          kl_maxiter=80,
                          cg_maxiter=30,
                          cg_tol=0.1,
                          damping=1.0,
                          verb=true)
    prob = _build_interferometric_problem(p, ft, data; weights=weights, verb=verb)

    # Per-iteration interferometric reporting.
    cb = (pr, z, samples, iter) -> begin
        verb || return nothing
        xi = z[1:_latent_size(p)]
        report_chi2(sky_forward(xi, p), p, ft, data)
        report_latents(xi, p)
        return nothing
    end

    z, samples = reconstruct_mgvi(prob;
                                   n_iterations=n_iterations,
                                   n_samples=n_samples,
                                   map_maxiter=map_maxiter,
                                   kl_maxiter=kl_maxiter,
                                   cg_maxiter=cg_maxiter,
                                   cg_tol=cg_tol,
                                   damping=damping,
                                   iter_callback=cb,
                                   verb=verb)

    image_mean, image_std = _posterior_image_stats(z, samples, p; antithetic=true)

    if verb
        @printf("  Mean image range: [%.4e, %.4e]\n",
                minimum(image_mean), maximum(image_mean))
        @printf("  Std image range:  [%.4e, %.4e]\n",
                minimum(image_std), maximum(image_std))
    end

    return z, image_mean, image_std, samples
end


# Domain glue for the interferometric inference problem.
#
# Three pieces:
#   • geovi_transformation / _right_sqrt_metric / _left_sqrt_metric — the
#     `InterferometricProblem` protocol implementations (`T(z) = model/σ`,
#     plus its JVP/VJP). These are called by VarInf's generic algorithm
#     helpers via the protocol dispatch in src/inference_problem.jl.
#   • reconstruct_geovi(p, ft, data; ...)  and
#     reconstruct_hybrid(p, ft, data; ...) — SkyModelParams-typed thin
#     wrappers that build an InterferometricProblem, run a domain MAP
#     warm-start, hand the algorithm to VarInf's generic
#     `reconstruct_geovi(prob; …)` / `reconstruct_hybrid(prob; …)` via an
#     `iter_callback` for per-iteration chi²/minisanity reporting, and
#     finally collapse samples into a posterior image mean/std.

# ============================================================================
# GeoVI metric operations  (InterferometricProblem protocol implementations)
# ============================================================================

"""
    geovi_transformation(z, p, ctx) -> normalized_residuals

Coordinate transformation T(z) = model_observables(z) / σ.
`z = [xi; params]` is the extended latent vector (or just `xi` if no model).
Loops over wavelength channels using per-channel `ObservationConfig`.
"""
function geovi_transformation(z::Vector{Float64}, p::SkyModelParams,
                              ctx::ObsContext)
    xi, params = _split_latent(z, p, ctx)
    image = sky_forward(xi, p)  # npix × npix × nf
    nf = length(p.freq)
    parts = Vector{Float64}[]
    for c in 1:nf
        obs_model_c, _ = observe(image[:, :, c], ctx.obs_vec[c]; params=params)
        push!(parts, obs_model_c)
    end
    if ctx.dp !== nothing
        push!(parts, vec(diffphi_forward(image, ctx.dp)))
    end
    return vcat(parts...) ./ ctx.sigma_vec
end

"""
    geovi_right_sqrt_metric(z, v, p, ctx) -> J_T(z) · v

JVP of the transformation w.r.t. extended latent `z = [xi; params]`.
Loops over wavelength channels using per-channel `ObservationConfig`.
"""
function geovi_right_sqrt_metric(z::Vector{Float64}, v::Vector{Float64},
                                 p::SkyModelParams, ctx::ObsContext)
    xi, params = _split_latent(z, p, ctx)
    n_xi = _latent_size(p)
    v_xi = v[1:n_xi]
    np = n_model_params(ctx)
    v_params = np > 0 ? v[n_xi+1:end] : nothing

    image, d_image = sky_forward_jvp(xi, v_xi, p)
    nf = length(p.freq)
    parts = Vector{Float64}[]
    for c in 1:nf
        _, _, d_obs_c = observe_jvp(image[:, :, c], d_image[:, :, c], ctx.obs_vec[c];
                                     params=params, d_params=v_params)
        push!(parts, d_obs_c)
    end
    if ctx.dp !== nothing
        push!(parts, vec(diffphi_jvp(image, d_image, ctx.dp)))
    end
    return vcat(parts...) ./ ctx.sigma_vec
end

"""
    geovi_left_sqrt_metric(z, v, p, ctx) -> J_T(z)' · v

VJP of the transformation w.r.t. extended latent `z = [xi; params]`.
Loops over wavelength channels, splitting v by per-channel data sizes.
"""
function geovi_left_sqrt_metric(z::Vector{Float64}, v::Vector{Float64},
                                p::SkyModelParams, ctx::ObsContext)
    xi, params = _split_latent(z, p, ctx)
    nf = length(p.freq)
    np = n_model_params(ctx)

    image = sky_forward(xi, p)  # npix × npix × nf
    v_scaled = v ./ ctx.sigma_vec

    g_image = zeros(p.npix, p.npix, nf)
    g_params_total = np > 0 ? zeros(np) : nothing
    offset = 0
    for c in 1:nf
        n_c = ctx.obs_vec[c].nv2 + ctx.obs_vec[c].nt3amp + ctx.obs_vec[c].nt3phi
        v_c = v_scaled[offset+1:offset+n_c]
        offset += n_c

        result_c = observe_adjoint(v_c, image[:, :, c], ctx.obs_vec[c]; params=params)
        if np > 0
            g_image_c, g_params_c = result_c
            g_image[:, :, c] .= g_image_c
            g_params_total .+= g_params_c
        else
            g_image[:, :, c] .= result_c
        end
    end

    # Differential phase adjoint (cross-channel)
    if ctx.dp !== nothing
        n_dp = ctx.dp.ncommon * nf
        v_dp = reshape(v_scaled[offset+1:offset+n_dp], ctx.dp.ncommon, nf)
        g_image .+= diffphi_adjoint(v_dp, image, ctx.dp)
    end

    g_xi = sky_adjoint(g_image, xi, p)
    return np > 0 ? vcat(g_xi, g_params_total) : g_xi
end


# ============================================================================
# Domain reconstruct_* wrappers
# ============================================================================

# Build the per-iteration reporting callback used by both reconstruct_geovi
# and reconstruct_hybrid. The callback closes over (p, ft, data) so VarInf
# never has to touch them.
function _iter_report_callback(p::SkyModelParams, ft, data; verb::Bool=true)
    n_xi = _latent_size(p)
    return function (_prob, z, samples, _iter)
        verb || return nothing
        xi = z[1:n_xi]
        report_chi2(sky_forward(xi, p), p, ft, data)
        # `samples` from reconstruct_geovi/hybrid stores antithetic pairs as
        # adjacent entries; the +δ half is a fair set of base samples for
        # minisanity.
        base = [s[1:n_xi] for s in samples[1:2:end]]
        minisanity(xi, base, p)
        return nothing
    end
end


"""
    reconstruct_geovi(p, ft, data; kwargs...)
        -> (xi_center, image_mean, image_std, samples)

Full GeoVI reconstruction with nonlinearly refined posterior sampling.
"""
function reconstruct_geovi(p::SkyModelParams, ft, data;
                           model::Union{Nothing, FlatModel}=nothing,
                           params_start::Union{Nothing, AbstractVector{Float64}}=nothing,
                           weights=[1.0, 1.0, 1.0],
                           n_iterations=6,
                           n_samples=iter -> iter <= 2 ? 2 : 4,
                           map_maxiter=200,
                           kl_maxiter=35,
                           kl_absdelta=0.5,
                           cg_maxiter=100,
                           cg_tol=0.01,
                           geo_newton_maxiter=10,
                           geo_cg_maxiter=50,
                           geo_tol=1e-5,
                           verb=true)
    prob = _build_interferometric_problem(p, ft, data;
                                           model=model, weights=weights, verb=verb)

    verb && println("=== GeoVI: MAP warm-start ===")
    xi, _ = reconstruct_map(p, ft, data;
                             weights=weights, maxiter=map_maxiter, verb=verb)
    z0 = params_start !== nothing ? vcat(xi[1:_latent_size(p)], params_start) : xi

    cb = _iter_report_callback(p, ft, data; verb=verb)

    z, samples = reconstruct_geovi(prob;
                                    z0=z0,
                                    n_iterations=n_iterations,
                                    n_samples=n_samples,
                                    kl_maxiter=kl_maxiter,
                                    kl_absdelta=kl_absdelta,
                                    cg_maxiter=cg_maxiter,
                                    cg_tol=cg_tol,
                                    geo_newton_maxiter=geo_newton_maxiter,
                                    geo_cg_maxiter=geo_cg_maxiter,
                                    geo_tol=geo_tol,
                                    iter_callback=cb, verb=verb)

    verb && println("\n=== Computing posterior statistics ===")
    image_mean, image_std = _posterior_image_stats(z, samples, p; antithetic=false)
    if verb
        @printf("  Mean image range: [%.4e, %.4e]\n", minimum(image_mean), maximum(image_mean))
        @printf("  Std image range:  [%.4e, %.4e]\n", minimum(image_std), maximum(image_std))
        if maximum(image_std) > 0
            @printf("  Peak SNR: %.1f\n", maximum(image_mean) / maximum(image_std))
        end
    end
    return z, image_mean, image_std, samples
end


"""
    reconstruct_hybrid(p, ft, data; kwargs...)
        -> (xi_center, image_mean, image_std, samples)

Hybrid reconstruction matching NIFTy's `optimize_kl` strategy.
KL minimization uses Newton-CG with sample-averaged posterior metric.
"""
function reconstruct_hybrid(p::SkyModelParams, ft, data;
                            model::Union{Nothing, FlatModel}=nothing,
                            params_start::Union{Nothing, AbstractVector{Float64}}=nothing,
                            weights=[1.0, 1.0, 1.0],
                            n_mgvi=10,
                            n_geovi=6,
                            n_samples=iter -> iter <= 2 ? 2 : 4,
                            map_maxiter=200,
                            kl_maxiter=35,
                            kl_absdelta=0.5,
                            cg_maxiter=100,
                            cg_tol=0.01,
                            geo_newton_maxiter=10,
                            geo_cg_maxiter=50,
                            geo_tol=1e-5,
                            frozen_ranges::Union{Nothing, Vector{UnitRange{Int}}}=nothing,
                            sample_mode=nothing,
                            verb=true)
    prob = _build_interferometric_problem(p, ft, data;
                                           model=model, weights=weights, verb=verb)

    verb && println("=== Hybrid: MAP warm-start ===")
    xi, _ = reconstruct_map(p, ft, data;
                             weights=weights, maxiter=map_maxiter, verb=verb)
    z0 = params_start !== nothing ? vcat(xi[1:_latent_size(p)], params_start) : xi

    cb = _iter_report_callback(p, ft, data; verb=verb)

    z, samples = reconstruct_hybrid(prob;
                                     z0=z0,
                                     n_mgvi=n_mgvi, n_geovi=n_geovi,
                                     n_samples=n_samples,
                                     kl_maxiter=kl_maxiter,
                                     kl_absdelta=kl_absdelta,
                                     cg_maxiter=cg_maxiter,
                                     cg_tol=cg_tol,
                                     geo_newton_maxiter=geo_newton_maxiter,
                                     geo_cg_maxiter=geo_cg_maxiter,
                                     geo_tol=geo_tol,
                                     frozen_ranges=frozen_ranges,
                                     sample_mode=sample_mode,
                                     iter_callback=cb, verb=verb)

    verb && println("\n=== Computing posterior statistics ===")
    image_mean, image_std = _posterior_image_stats(z, samples, p; antithetic=false)
    if verb
        @printf("  Mean image range: [%.4e, %.4e]\n", minimum(image_mean), maximum(image_mean))
        @printf("  Std image range:  [%.4e, %.4e]\n", minimum(image_std), maximum(image_std))
        if maximum(image_std) > 0
            @printf("  Peak SNR: %.1f\n", maximum(image_mean) / maximum(image_std))
        end
    end
    return z, image_mean, image_std, samples
end


"""
    _require_float64(ft)

Refuse Float32 plans up front, with the fix, instead of a `MethodError` ten frames down.

`readoifits` and `setup_ft` default to **Float32**, and that is the right default for image
reconstruction. Variational inference is not image reconstruction: VarInf's correlated field is
built on `ComplexF64` FFT plans and its CG / Newton-CG solvers are conditioned for double
precision, so the whole latent space here is Float64. Promoting silently would mean rebuilding
every NFFT plan behind the caller's back — expensive, and surprising when it shows up as a
pause rather than an error.

This was the single thing that broke when the package was last run against a current OITOOLS,
and it surfaced as `_obs_to_g_cvis(::Vector{Float64}, ::Vector{ComplexF32}, …)` — a message
that says nothing about what to do.
"""
function _require_float64(ft)
    T = try
        ft_eltype(ft)
    catch
        return nothing        # not a shape we can inspect; let the caller proceed
    end
    T === Float64 && return nothing
    error("""
        Variational inference runs at Float64, and these Fourier plans are $(T).

        Re-read the data and rebuild the plans at double precision:

            data = readoifits(file; T = Float64)
            ft   = setup_ft(data, nx, pixsize)

        `readoifits` and `setup_ft` default to Float32, which is right for image
        reconstruction and not for this.
        """)
end

# Shared helpers for the SkyModelParams-typed reconstruct_* domain wrappers:
# build an InterferometricProblem from (p, ft, data) (including optional
# diff-phase config) and compute posterior image statistics from a list of
# samples in either of the two conventions used by the MGVI / GeoVI loops.

"""
    _build_interferometric_problem(p, ft, data; model=nothing,
                                    weights=[1, 1, 1], verb=true)
        -> InterferometricProblem

Build the per-channel `ObservationConfig`s, the concatenated σ vector,
the optional cross-channel differential-phase config, the `ObsContext`,
and finally wrap everything in an `InterferometricProblem` ready to feed
into VarInf's generic reconstruct_* methods.
"""
function _build_interferometric_problem(p::SkyModelParams, ft, data;
                                         model::Union{Nothing, FlatModel}=nothing,
                                         weights::AbstractVector=[1.0, 1.0, 1.0],
                                         verb::Bool=true)
    _require_float64(ft)
    nf = length(p.freq)
    obs_vec = ObservationConfig[]
    for c in 1:nf
        ft_ch   = ft   isa AbstractMatrix ? ft[c, 1]   : (nf == 1 ? ft : ft[c])
        data_ch = data isa AbstractMatrix ? data[c, 1] : (nf == 1 ? data : data[c])
        push!(obs_vec, ObservationConfig(ft_ch, data_ch; model=model))
    end
    sigma_vec = vcat([vcat(o.v2_err, o.t3amp_err, o.t3phi_err) for o in obs_vec]...)

    dp = (ft isa AbstractMatrix && data isa AbstractMatrix) ?
         build_diffphase_config(ft, data, nf) : nothing
    if dp !== nothing
        sigma_vec = vcat(sigma_vec, vec(dp.diffphi_err))
        verb && @printf("  DiffPhase enabled: %d common baselines × %d channels\n",
                        dp.ncommon, nf)
    end

    ctx = ObsContext(obs_vec, sigma_vec, dp)
    return InterferometricProblem(p, ft, data, ctx; weights=weights)
end


"""
    _posterior_image_stats(z, samples, p; antithetic) -> (image_mean, image_std)

Aggregate posterior image statistics from `samples` by forward-modelling
`z + δ_k` through `sky_forward`. Two sample conventions:

- `antithetic = true`  (MGVI):
    Each entry of `samples` represents one side of an antithetic pair; we
    expand to z ± samples[k] before forward-modelling.
- `antithetic = false` (GeoVI / Hybrid):
    Entries of `samples` are already explicit pairs (the +δ and −δ refinements
    after Newton-CG sample refinement), so we use each one as-is.
"""
function _posterior_image_stats(z::AbstractVector{Float64},
                                 samples::Vector{Vector{Float64}},
                                 p::SkyModelParams;
                                 antithetic::Bool)
    n_xi = _latent_size(p)
    nf = length(p.freq)
    image_sum  = zeros(p.npix, p.npix, nf)
    image_sum2 = zeros(p.npix, p.npix, nf)
    count = 0
    function _accumulate!(z_full)
        xi = z_full[1:n_xi]
        img = sky_forward(xi, p)
        image_sum  .+= img
        image_sum2 .+= img .^ 2
        count += 1
        return nothing
    end
    if antithetic
        for s in samples
            _accumulate!(z .+ s)
            _accumulate!(z .- s)
        end
    else
        for s in samples
            _accumulate!(z .+ s)
        end
    end
    image_mean = image_sum ./ count
    image_std  = sqrt.(max.(image_sum2 ./ count .- image_mean .^ 2, 0.0))
    return image_mean, image_std
end
