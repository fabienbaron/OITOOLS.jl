# ============================================================================
# The inference problem: likelihood energy, and the VarInf protocol wrappers.
# ============================================================================

"""
    _split_latent(z, p, obs) -> (xi, params_or_nothing)

Split the extended latent vector into image latents and parametric model params.
If no model is present, `z` is just `xi` and params is `nothing`.
"""
function _split_latent(z::AbstractVector{Float64}, p::SkyModelParams,
                       obs::ObservationConfig)
    n_xi = _latent_size(p)
    np = n_model_params(obs)
    if np > 0
        @assert length(z) == n_xi + np "Expected latent of size $(n_xi + np), got $(length(z))"
        return z[1:n_xi], z[n_xi+1:end]
    else
        return z, nothing
    end
end

"""
    total_latent_size(p, obs) -> Int

Total latent vector size: image latents + parametric model params.
"""
total_latent_size(p::SkyModelParams, obs::ObservationConfig) =
    _latent_size(p) + n_model_params(obs)
total_latent_size(p::SkyModelParams, obs_vec::Vector{ObservationConfig}) =
    total_latent_size(p, obs_vec[1])
total_latent_size(p::SkyModelParams, ctx::ObsContext) =
    total_latent_size(p, ctx.obs_vec[1])

function _split_latent(z::AbstractVector{Float64}, p::SkyModelParams,
                       obs_vec::Vector{ObservationConfig})
    _split_latent(z, p, obs_vec[1])
end
_split_latent(z::AbstractVector{Float64}, p::SkyModelParams, ctx::ObsContext) =
    _split_latent(z, p, ctx.obs_vec[1])

"""
    energy_fg(z, p, ft, data, obs; weights) -> (energy, g_z)

Compute posterior energy and gradient w.r.t. the extended latent vector
`z = [xi; params]` (or just `xi` if no parametric model).

Energy = chi²(image(xi), params) + 0.5*||xi||² + 0.5*||params||²

For hybrid models, uses OITOOLS's `model_and_image_to_chi2_fg` with SPARCO flux.
The parametric prior is a standard Gaussian (centered at 0, unit variance) on the
raw parameter latents — same as the image prior.
"""
function energy_fg(z::AbstractVector{Float64}, p::SkyModelParams,
                   ft, data, obs::ObservationConfig;
                   weights=[1.0, 1.0, 1.0])
    xi, params = _split_latent(z, p, obs)
    nf = length(p.freq)
    image = sky_forward(xi, p)  # npix × npix × nf

    g_image = zeros(p.npix, p.npix, nf)
    chi2 = 0.0
    g_params = params !== nothing ? zeros(length(params)) : nothing

    if params !== nothing
        # Hybrid mode: use OITOOLS's model_and_image_to_chi2_fg
        for ch in 1:nf
            g_ch = zeros(p.npix, p.npix)
            ft_ch   = ft   isa AbstractMatrix ? ft[ch,1]   : (nf == 1 ? ft : ft[ch])
            data_ch = data isa AbstractMatrix ? data[ch,1] : (nf == 1 ? data : data[ch])
            chi2_ch, g_params_ch = model_and_image_to_chi2_fg(
                obs.model, params, image[:, :, ch], g_ch, ft_ch, data_ch;
                weights=weights, verb=false)
            chi2 += chi2_ch
            g_image[:, :, ch] .= g_ch
            g_params .+= g_params_ch
        end
    else
        # Image-only mode
        for ch in 1:nf
            g_ch = zeros(p.npix, p.npix)
            ft_ch   = ft   isa AbstractMatrix ? ft[ch,1]   : (nf == 1 ? ft : ft[ch])
            data_ch = data isa AbstractMatrix ? data[ch,1] : (nf == 1 ? data : data[ch])
            chi2 += image_to_chi2_fg(image[:, :, ch], g_ch, ft_ch, data_ch;
                                      weights=weights, verb=false)
            g_image[:, :, ch] .= g_ch
        end
    end

    # Backprop through sky model: g_xi = J' * g_image
    g_xi = sky_adjoint(g_image, xi, p)

    # Gaussian prior on image latents
    prior = 0.5 * dot(xi, xi)
    g_xi .+= xi

    if g_params !== nothing
        # Gaussian prior on parametric params
        prior += 0.5 * dot(params, params)
        g_params .+= params
        return chi2 + prior, vcat(g_xi, g_params)
    else
        return chi2 + prior, g_xi
    end
end

# Backward-compatible signature without obs (image-only)
function energy_fg(xi::AbstractVector{Float64}, p::SkyModelParams,
                   ft, data;
                   weights=[1.0, 1.0, 1.0])
    nf = length(p.freq)
    image = sky_forward(xi, p)

    g_image = zeros(p.npix, p.npix, nf)
    chi2 = 0.0
    for ch in 1:nf
        g_ch = zeros(p.npix, p.npix)
        ft_ch   = ft   isa AbstractMatrix ? ft[ch,1]   : (nf == 1 ? ft : ft[ch])
        data_ch = data isa AbstractMatrix ? data[ch,1] : (nf == 1 ? data : data[ch])
        chi2 += image_to_chi2_fg(image[:, :, ch], g_ch, ft_ch, data_ch;
                                  weights=weights, verb=false)
        g_image[:, :, ch] .= g_ch
    end

    g_xi = sky_adjoint(g_image, xi, p)
    prior = 0.5 * dot(xi, xi)
    g_xi .+= xi
    return chi2 + prior, g_xi
end

"""
    energy_fg(z, p, ft, data, ctx; weights) -> (energy, g_z)

ObsContext dispatch: delegates to the per-channel `energy_fg` for V²/T3,
then adds diffphi chi² and gradient if `ctx.dp !== nothing`.
"""
function energy_fg(z::AbstractVector{Float64}, p::SkyModelParams,
                   ft, data, ctx::ObsContext;
                   weights=[1.0, 1.0, 1.0])
    # Per-channel chi2 via existing dispatch
    obs1 = ctx.obs_vec[1]
    chi2_base, g_base = energy_fg(z, p, ft, data, obs1; weights=weights)

    if ctx.dp === nothing
        return chi2_base, g_base
    end

    # Add diffphi chi²: Σ [mod360(φ_model - φ_data) / σ]²
    xi, params = _split_latent(z, p, ctx)
    image = sky_forward(xi, p)
    dp = ctx.dp
    nf = length(dp.ft_vis)

    dphi_model = diffphi_forward(image, dp)
    resid = mod360(dphi_model .- dp.diffphi_data)
    chi2_dp = sum((resid ./ dp.diffphi_err) .^ 2)

    # Gradient: ∂chi2_dp/∂φ = 2·mod360(φ_model - φ_data)/σ², then chain through adjoint
    g_dphi = 2.0 .* resid ./ (dp.diffphi_err .^ 2)
    g_image_dp = diffphi_adjoint(g_dphi, image, dp)
    g_xi_dp = sky_adjoint(g_image_dp, xi, p)

    g_total = copy(g_base)
    n_xi = _latent_size(p)
    g_total[1:n_xi] .+= g_xi_dp
    return chi2_base + chi2_dp, g_total
end

"""
    energy_fg!(xi, g_xi, p, ft, data; weights) -> energy

In-place version compatible with OptimPackNextGen.vmlmb interface.
Fills g_xi with the gradient and returns the scalar energy.
"""
function energy_fg!(xi::AbstractVector{Float64}, g_xi::AbstractVector{Float64},
                    p::SkyModelParams, ft, data;
                    weights=[1.0, 1.0, 1.0])
    e, g = energy_fg(xi, p, ft, data; weights=weights)
    g_xi .= g
    return e
end


# ============================================================================
# Concrete AbstractInferenceProblem subtypes — InterferometricProblem and
# PointSourceProblem — with method implementations that delegate to existing
# free functions (`energy_fg`, `geovi_*`, `ps_energy_fg`, `ps_geovi_*`).
#
# The abstract type + function stubs live in src/abstract_inference_problem.jl
# (included earlier so geovi.jl can type-annotate prob::AbstractInferenceProblem).
# ============================================================================


# ============================================================================
# InterferometricProblem: sky image + NFFT observation + diffphase
# ============================================================================

"""
    InterferometricProblem(p, ft, data, ctx; weights=[1.0,1.0,1.0])

Wraps a polychromatic sky-image inference problem (SkyModelParams + per-channel
OITOOLS plans + ObsContext) as an `AbstractInferenceProblem`. The extended
latent z is `[xi; params]` when a parametric model is attached, else just `xi`.
"""
struct InterferometricProblem <: AbstractInferenceProblem
    p::SkyModelParams
    ft::Any
    data::Any
    ctx::ObsContext
    weights::Vector{Float64}
end

InterferometricProblem(p::SkyModelParams, ft, data, ctx::ObsContext;
                       weights::AbstractVector=[1.0, 1.0, 1.0]) =
    InterferometricProblem(p, ft, data, ctx, Vector{Float64}(weights))

energy_and_gradient(prob::InterferometricProblem, z::AbstractVector{Float64}) =
    energy_fg(z, prob.p, prob.ft, prob.data, prob.ctx; weights=prob.weights)

latent_size(prob::InterferometricProblem) =
    total_latent_size(prob.p, prob.ctx)

transformation(prob::InterferometricProblem, z::AbstractVector{Float64}) =
    geovi_transformation(Vector{Float64}(z), prob.p, prob.ctx)

right_sqrt_metric(prob::InterferometricProblem,
                  z::AbstractVector{Float64}, v::AbstractVector{Float64}) =
    geovi_right_sqrt_metric(Vector{Float64}(z), Vector{Float64}(v),
                            prob.p, prob.ctx)

left_sqrt_metric(prob::InterferometricProblem,
                 z::AbstractVector{Float64}, v::AbstractVector{Float64}) =
    geovi_left_sqrt_metric(Vector{Float64}(z), Vector{Float64}(v),
                           prob.p, prob.ctx)

data_size(prob::InterferometricProblem) = length(prob.ctx.sigma_vec)


# ============================================================================
# PointSourceProblem: analytic multi-point-source visibility model
# ============================================================================

"""
    PointSourceProblem(ps, data; weights=[1.0,1.0,1.0])

Wraps a polychromatic point-source inference problem (PointSourceParams + data)
as an `AbstractInferenceProblem`. The latent z has length `2N + N·nf`
(positions + per-channel log-fluxes).
"""
struct PointSourceProblem <: AbstractInferenceProblem
    ps::PointSourceParams
    data::Any
    weights::Vector{Float64}
end

PointSourceProblem(ps::PointSourceParams, data;
                   weights::AbstractVector=[1.0, 1.0, 1.0]) =
    PointSourceProblem(ps, data, Vector{Float64}(weights))

energy_and_gradient(prob::PointSourceProblem, z::AbstractVector{Float64}) =
    ps_energy_fg(z, prob.ps, prob.data; weights=prob.weights)

latent_size(prob::PointSourceProblem) = ps_latent_size(prob.ps)

transformation(prob::PointSourceProblem, z::AbstractVector{Float64}) =
    ps_geovi_transformation(Vector{Float64}(z), prob.ps)

right_sqrt_metric(prob::PointSourceProblem,
                  z::AbstractVector{Float64}, v::AbstractVector{Float64}) =
    ps_geovi_right_sqrt_metric(Vector{Float64}(z), Vector{Float64}(v), prob.ps)

left_sqrt_metric(prob::PointSourceProblem,
                 z::AbstractVector{Float64}, v::AbstractVector{Float64}) =
    ps_geovi_left_sqrt_metric(Vector{Float64}(z), Vector{Float64}(v), prob.ps)

data_size(prob::PointSourceProblem) = length(prob.ps.sigma_vec)
