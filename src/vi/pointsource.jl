# ============================================================================
# Point sources: analytic visibilities, energy, and their own geoVI problem.
# ============================================================================

# ============================================================================
# Polychromatic multi-point-source model
# ============================================================================
#
# Represents a target as N point sources, each with a flux vector across
# spectral channels.  Complex visibilities computed analytically from the
# point source positions and fluxes (no image grid, no NFFT).
# Flux-normalized so V(0)=1, matching OITOOLS V² calibration.
# ============================================================================

# ── Struct ──────────────────────────────────────────────────────────────────

"""
    PointSourceParams

Parameters for a polychromatic multi-point-source model.
The latent vector `z ~ N(0, I)` maps to physical parameters via:
  - positions: `x_i = pos_mean[1,i] + pos_std * z[i]` (mas → radians)
  - log-fluxes: `logf[i,c] = logf_mean[i,c] + logf_std * z[2N+(c-1)*N+i]`
  - fluxes: `f[i,c] = exp(logf[i,c])`
"""
struct PointSourceParams
    N::Int                             # number of point sources
    nf::Int                            # number of spectral channels
    freq::Vector{Float64}              # channel frequencies (Hz)
    pos_mean::Matrix{Float64}          # (2, N) prior mean positions in mas
    pos_std::Float64                   # position prior scale (mas)
    logf_mean::Matrix{Float64}         # (N, nf) log-flux prior means
    logf_std::Float64                  # log-flux prior scale
    obs_vec::Vector{ObservationConfig} # per-channel observation configs
    sigma_vec::Vector{Float64}         # concatenated errors for GeoVI
end

"""
    PointSourceParams(N, data, ft; pos_init, pos_std, logf_mean, logf_std)

Construct from OITOOLS data and ft plans.

`pos_init` is a `(2, N)` matrix of initial positions in mas.
`logf_mean` defaults to zeros (equal prior flux for all sources/channels).
"""
function PointSourceParams(N::Int, data, ft;
                           pos_init::Matrix{Float64},
                           pos_std::Float64=0.5,
                           logf_mean::Union{Nothing, Matrix{Float64}}=nothing,
                           logf_std::Float64=1.0)
    nf = data isa AbstractMatrix ? size(data, 1) : 1

    freq = Float64[]
    obs_vec = ObservationConfig[]
    for c in 1:nf
        ft_ch   = ft   isa AbstractMatrix ? ft[c,1]   : (nf == 1 ? ft : ft[c])
        data_ch = data isa AbstractMatrix ? data[c,1] : (nf == 1 ? data : data[c])
        push!(obs_vec, ObservationConfig(ft_ch, data_ch))
        push!(freq, 3e8 / mean(data_ch.uv_lam))
    end

    sigma_vec = vcat([vcat(o.v2_err, o.t3amp_err, o.t3phi_err) for o in obs_vec]...)

    lf = logf_mean === nothing ? zeros(N, nf) : logf_mean
    @assert size(pos_init) == (2, N) "pos_init must be (2, $N)"
    @assert size(lf) == (N, nf) "logf_mean must be ($N, $nf)"

    return PointSourceParams(N, nf, freq, copy(pos_init), pos_std,
                             copy(lf), logf_std, obs_vec, sigma_vec)
end

# ── Latent layout ───────────────────────────────────────────────────────────

"""
    ps_latent_size(ps) -> Int

Total latent vector length: `2N + N*nf`.
"""
ps_latent_size(ps::PointSourceParams) = 2 * ps.N + ps.N * ps.nf

"""
    ps_unpack(z, ps) -> (x_rad, y_rad, flux)

Map latent `z` to physical point-source parameters.
Returns positions in **radians** and flux matrix `(N, nf)`.
"""
function ps_unpack(z::AbstractVector{Float64}, ps::PointSourceParams)
    N, nf = ps.N, ps.nf
    x_rad = (ps.pos_mean[1,:] .+ ps.pos_std .* z[1:N]) .* MAS2RAD
    y_rad = (ps.pos_mean[2,:] .+ ps.pos_std .* z[N+1:2N]) .* MAS2RAD
    flux = Matrix{Float64}(undef, N, nf)
    for c in 1:nf
        off = 2N + (c - 1) * N
        @inbounds for i in 1:N
            flux[i, c] = exp(ps.logf_mean[i, c] + ps.logf_std * z[off + i])
        end
    end
    return x_rad, y_rad, flux
end

"""
    ps_unpack_tangent(v, ps, flux) -> (dx_rad, dy_rad, dflux)

Tangent of `ps_unpack`: map tangent vector `v` to physical tangents.
`flux` is the primal flux matrix (needed for the exp chain rule).
"""
function ps_unpack_tangent(v::AbstractVector{Float64}, ps::PointSourceParams,
                           flux::Matrix{Float64})
    N, nf = ps.N, ps.nf
    scale_pos = ps.pos_std * MAS2RAD
    dx_rad = scale_pos .* v[1:N]
    dy_rad = scale_pos .* v[N+1:2N]
    dflux = Matrix{Float64}(undef, N, nf)
    for c in 1:nf
        off = 2N + (c - 1) * N
        @inbounds for i in 1:N
            dflux[i, c] = flux[i, c] * ps.logf_std * v[off + i]
        end
    end
    return dx_rad, dy_rad, dflux
end

"""
    ps_adjoint_chain!(g_z, g_x, g_y, g_flux, flux, ps, channel)

Accumulate the VJP of `ps_unpack` for one channel into `g_z`.
Chains g_x, g_y (position gradients in radians) and g_flux (flux gradient)
back through the prior scales and exp transform.
"""
function ps_adjoint_chain!(g_z::AbstractVector{Float64},
                           g_x::Vector{Float64}, g_y::Vector{Float64},
                           g_flux::Vector{Float64},
                           flux::Matrix{Float64},
                           ps::PointSourceParams, channel::Int)
    N = ps.N
    scale_pos = ps.pos_std * MAS2RAD
    @inbounds for i in 1:N
        g_z[i]     += g_x[i] * scale_pos
        g_z[N + i] += g_y[i] * scale_pos
    end
    off = 2N + (channel - 1) * N
    @inbounds for i in 1:N
        g_z[off + i] += g_flux[i] * flux[i, channel] * ps.logf_std
    end
    return nothing
end


# ── Complex visibilities (flux-normalized) ──────────────────────────────────

"""
    ps_cvis(x_rad, y_rad, flux_c, obs) -> cvis

Compute flux-normalized complex visibilities for one channel:
  `V[k] = (1/F) * Σ_i f_i * exp(-2πi(u[k]*x_i + v[k]*y_i))`
where `F = Σ_i f_i`.  Direct analytic evaluation — no transform grid.
"""
function ps_cvis(x_rad::AbstractVector{Float64}, y_rad::AbstractVector{Float64},
                 flux_c::AbstractVector{Float64}, obs::ObservationConfig)
    uv = obs.uv
    nuv = size(uv, 2)
    N = length(x_rad)
    inv_flux = 1.0 / sum(flux_c)
    cvis = zeros(ComplexF64, nuv)
    @inbounds for i in 1:N
        fi = flux_c[i] * inv_flux
        xi, yi = x_rad[i], y_rad[i]
        for k in 1:nuv
            phase = -2.0 * π * (uv[1, k] * xi + uv[2, k] * yi)
            cvis[k] += fi * cis(phase)
        end
    end
    return cvis
end


# ── Observable extraction from cvis ─────────────────────────────────────────

"""
    ps_observables(cvis, obs) -> (obs_model, obs_err)

Extract V², T3amp, T3phi from complex visibilities (same as `observe` but
without an image).
"""
function ps_observables(cvis::AbstractVector{ComplexF64}, obs::ObservationConfig)
    v2_model = abs2.(cvis[obs.indx_v2])
    t3 = cvis[obs.indx_t3_1] .* cvis[obs.indx_t3_2] .* cvis[obs.indx_t3_3]
    t3amp_model = abs.(t3)
    t3phi_model = angle.(t3) .* (180.0 / π)
    obs_model = vcat(v2_model, t3amp_model, t3phi_model)
    obs_err   = vcat(obs.v2_err, obs.t3amp_err, obs.t3phi_err)
    return obs_model, obs_err
end


# ── JVP of visibility computation ────────────────────────────────────────────

"""
    ps_cvis_jvp(x_rad, y_rad, flux_c, dx_rad, dy_rad, dflux_c, obs)
        -> (cvis, d_cvis)

Forward-mode tangent of `ps_cvis`.
Returns the primal cvis and its tangent d_cvis (both flux-normalized).
"""
function ps_cvis_jvp(x_rad::AbstractVector{Float64}, y_rad::AbstractVector{Float64},
                     flux_c::AbstractVector{Float64},
                     dx_rad::AbstractVector{Float64}, dy_rad::AbstractVector{Float64},
                     dflux_c::AbstractVector{Float64},
                     obs::ObservationConfig)
    uv = obs.uv
    nuv = size(uv, 2)
    N = length(x_rad)
    F = sum(flux_c)
    dF = sum(dflux_c)
    inv_F = 1.0 / F

    cvis_raw  = zeros(ComplexF64, nuv)
    d_cvis_raw = zeros(ComplexF64, nuv)

    @inbounds for i in 1:N
        fi, dfi = flux_c[i], dflux_c[i]
        xi, yi = x_rad[i], y_rad[i]
        dxi, dyi = dx_rad[i], dy_rad[i]
        for k in 1:nuv
            phase  = -2.0 * π * (uv[1, k] * xi  + uv[2, k] * yi)
            dphase = -2.0 * π * (uv[1, k] * dxi + uv[2, k] * dyi)
            e = cis(phase)
            cvis_raw[k]  += fi * e
            d_cvis_raw[k] += dfi * e + fi * e * (im * dphase)
        end
    end

    # Quotient rule: cvis = cvis_raw / F
    inv_F2 = inv_F * inv_F
    cvis   = cvis_raw .* inv_F
    d_cvis = (d_cvis_raw .* F .- cvis_raw .* dF) .* inv_F2

    return cvis, d_cvis
end


# ── VJP of visibility computation ────────────────────────────────────────────

"""
    ps_cvis_adjoint(g_cvis, x_rad, y_rad, flux_c, obs)
        -> (g_x, g_y, g_flux)

Adjoint (VJP) of `ps_cvis`.
Returns gradients w.r.t. physical parameters (radians for positions,
raw flux for g_flux).  Caller must chain through `ps_adjoint_chain!`
to map back to latent z.
"""
function ps_cvis_adjoint(g_cvis::AbstractVector{ComplexF64},
                         x_rad::AbstractVector{Float64},
                         y_rad::AbstractVector{Float64},
                         flux_c::AbstractVector{Float64},
                         obs::ObservationConfig)
    uv = obs.uv
    nuv = size(uv, 2)
    N = length(x_rad)
    F = sum(flux_c)
    inv_F = 1.0 / F

    # Adjoint of normalization: g_cvis_raw = g_cvis / F
    # Position adjoint uses f_norm = f/F (normalization doesn't depend on pos)
    g_f_raw = zeros(N)
    g_x = zeros(N)
    g_y = zeros(N)

    @inbounds for i in 1:N
        f_norm_i = flux_c[i] * inv_F
        xi, yi = x_rad[i], y_rad[i]
        sf = 0.0 + 0.0im
        sx = 0.0 + 0.0im
        sy = 0.0 + 0.0im
        for k in 1:nuv
            phase = -2.0 * π * (uv[1, k] * xi + uv[2, k] * yi)
            e = cis(phase)
            gc = conj(g_cvis[k])
            sf += gc * e
            sx += gc * uv[1, k] * e
            sy += gc * uv[2, k] * e
        end
        g_f_raw[i] = real(sf)                   # unnormalized DFT adjoint
        g_x[i] = f_norm_i * real((-2.0π * im) * sx)
        g_y[i] = f_norm_i * real((-2.0π * im) * sy)
    end

    # Flux normalization correction (same pattern as _g_cvis_to_g_image):
    # g_flux[i] = (g_f_raw[i] - correction) / F
    f_norm = flux_c .* inv_F
    correction = dot(g_f_raw, f_norm)
    g_flux = (g_f_raw .- correction) .* inv_F

    return g_x, g_y, g_flux
end


# ── Energy (posterior = chi² + Gaussian prior) ──────────────────────────────

"""
    ps_energy_fg(z, ps, data; weights) -> (energy, g_z)

Compute posterior energy and gradient for the point-source model.

Energy = Σ_c chi²_c(z) + 0.5 ||z||²

where `chi²_c` includes weighted V², T3amp, T3phi residuals.
"""
function ps_energy_fg(z::AbstractVector{Float64}, ps::PointSourceParams,
                      data; weights=[1.0, 1.0, 1.0])
    N, nf = ps.N, ps.nf
    n = ps_latent_size(ps)
    g_z = zeros(n)
    chi2 = 0.0

    x_rad, y_rad, flux = ps_unpack(z, ps)

    for c in 1:nf
        obs = ps.obs_vec[c]
        data_c = data isa AbstractMatrix ? data[c, 1] : data

        # Forward
        cvis_c = ps_cvis(x_rad, y_rad, flux[:, c], obs)
        obs_model, obs_err = ps_observables(cvis_c, obs)
        obs_data = observe_data(data_c)

        # Residuals with phase wrapping
        nv2, nt3 = obs.nv2, obs.nt3amp
        resid = obs_model .- obs_data
        @views resid[nv2+nt3+1:end] .= mod360.(resid[nv2+nt3+1:end])

        # Weighted chi²
        w_v2   = weights[1]
        w_t3a  = weights[2]
        w_t3p  = weights[3]
        for j in 1:nv2
            r = resid[j] / obs_err[j]
            chi2 += w_v2 * r * r
        end
        for j in 1:nt3
            r = resid[nv2 + j] / obs_err[nv2 + j]
            chi2 += w_t3a * r * r
        end
        for j in 1:nt3
            r = resid[nv2 + nt3 + j] / obs_err[nv2 + nt3 + j]
            chi2 += w_t3p * r * r
        end

        # Observable gradient
        g_obs = zeros(nv2 + 2 * nt3)
        @inbounds for j in 1:nv2
            g_obs[j] = 2.0 * w_v2 * resid[j] / (obs_err[j]^2)
        end
        @inbounds for j in 1:nt3
            g_obs[nv2 + j] = 2.0 * w_t3a * resid[nv2 + j] / (obs_err[nv2 + j]^2)
        end
        @inbounds for j in 1:nt3
            g_obs[nv2 + nt3 + j] = 2.0 * w_t3p * resid[nv2 + nt3 + j] / (obs_err[nv2 + nt3 + j]^2)
        end

        # VJP: observables → cvis → physical params → z
        g_cvis = _obs_to_g_cvis(g_obs, cvis_c, obs)
        g_x, g_y, g_flux = ps_cvis_adjoint(g_cvis, x_rad, y_rad, flux[:, c], obs)
        ps_adjoint_chain!(g_z, g_x, g_y, g_flux, flux, ps, c)
    end

    # Gaussian prior: 0.5 ||z||²
    prior = 0.5 * dot(z, z)
    g_z .+= z

    return chi2 + prior, g_z
end

"""
    ps_energy_fg!(z, g_z, ps, data; weights) -> energy

In-place version for OptimPackNextGen.vmlmb interface.
"""
function ps_energy_fg!(z::AbstractVector{Float64}, g_z::AbstractVector{Float64},
                       ps::PointSourceParams, data;
                       weights=[1.0, 1.0, 1.0])
    e, g = ps_energy_fg(z, ps, data; weights=weights)
    g_z .= g
    return e
end


# ============================================================================
# Point-source domain glue
# ============================================================================
#
# Three GeoVI metric operations specific to the point-source forward model
# (the protocol implementations), the PointSourcePosterior struct, chi²
# reporting, and the `reconstruct_pointsource` entry point.
#
# The VI algorithm bodies live in VarInf.jl; reconstruct_pointsource is a
# thin wrapper that builds a PointSourceProblem, warm-starts via the
# generic reconstruct_map(prob; ...), then delegates the hybrid MGVI →
# GeoVI loop to reconstruct_hybrid(prob; iter_callback=...) for the
# point-source-specific per-iteration chi²/parameter reporting. The
# PointSourcePosterior aggregation stays here because it requires the
# physical (mas / flux fraction) unpack.

# ── GeoVI metric operations (PointSourceProblem protocol implementations) ──

"""
    ps_geovi_transformation(z, ps) -> normalized_residuals

Coordinate transformation T(z) = obs_model(z) / σ.
"""
function ps_geovi_transformation(z::Vector{Float64}, ps::PointSourceParams)
    x_rad, y_rad, flux = ps_unpack(z, ps)
    parts = Vector{Float64}[]
    for c in 1:ps.nf
        cvis = ps_cvis(x_rad, y_rad, flux[:, c], ps.obs_vec[c])
        obs_model, _ = ps_observables(cvis, ps.obs_vec[c])
        push!(parts, obs_model)
    end
    return vcat(parts...) ./ ps.sigma_vec
end

"""
    ps_geovi_right_sqrt_metric(z, v, ps) -> J_T(z) · v

JVP of the transformation w.r.t. latent z.
"""
function ps_geovi_right_sqrt_metric(z::Vector{Float64}, v::Vector{Float64},
                                     ps::PointSourceParams)
    x_rad, y_rad, flux = ps_unpack(z, ps)
    dx_rad, dy_rad, dflux = ps_unpack_tangent(v, ps, flux)

    parts = Vector{Float64}[]
    for c in 1:ps.nf
        cvis, d_cvis = ps_cvis_jvp(x_rad, y_rad, flux[:, c],
                                     dx_rad, dy_rad, dflux[:, c],
                                     ps.obs_vec[c])
        _, _, d_obs = _cvis_to_obs_jvp(cvis, d_cvis, ps.obs_vec[c])
        push!(parts, d_obs)
    end
    return vcat(parts...) ./ ps.sigma_vec
end

"""
    ps_geovi_left_sqrt_metric(z, v, ps) -> J_T(z)' · v

VJP of the transformation w.r.t. latent z.
"""
function ps_geovi_left_sqrt_metric(z::Vector{Float64}, v::Vector{Float64},
                                    ps::PointSourceParams)
    x_rad, y_rad, flux = ps_unpack(z, ps)
    v_scaled = v ./ ps.sigma_vec

    n = ps_latent_size(ps)
    g_z = zeros(n)
    offset = 0
    for c in 1:ps.nf
        obs = ps.obs_vec[c]
        n_c = obs.nv2 + obs.nt3amp + obs.nt3phi
        v_c = v_scaled[offset+1:offset+n_c]
        offset += n_c

        cvis = ps_cvis(x_rad, y_rad, flux[:, c], obs)
        g_cvis = _obs_to_g_cvis(v_c, cvis, obs)
        g_x, g_y, g_flux = ps_cvis_adjoint(g_cvis, x_rad, y_rad, flux[:, c], obs)
        ps_adjoint_chain!(g_z, g_x, g_y, g_flux, flux, ps, c)
    end
    return g_z
end


# ── Posterior output ─────────────────────────────────────────────────────────

"""
    PointSourcePosterior

Posterior summary: per-source positions and per-channel fluxes with
uncertainties.
"""
struct PointSourcePosterior
    x_mas::Vector{Float64}        # mean x positions (N,)
    y_mas::Vector{Float64}        # mean y positions (N,)
    x_std::Vector{Float64}        # x uncertainties (N,)
    y_std::Vector{Float64}        # y uncertainties (N,)
    flux::Matrix{Float64}         # mean flux fractions (N, nf)
    flux_std::Matrix{Float64}     # flux uncertainties (N, nf)
end


# ── Chi² reporting ───────────────────────────────────────────────────────────

"""
    ps_report_chi2(z, ps, data; weights, verb) -> chi2_total

Report per-channel and per-observable chi² (reduced).
"""
function ps_report_chi2(z::Vector{Float64}, ps::PointSourceParams, data;
                        weights=[1.0, 1.0, 1.0], verb=true)
    x_rad, y_rad, flux = ps_unpack(z, ps)
    chi2_total = 0.0
    n_total = 0

    for c in 1:ps.nf
        obs = ps.obs_vec[c]
        data_c = data isa AbstractMatrix ? data[c, 1] : data

        cvis = ps_cvis(x_rad, y_rad, flux[:, c], obs)
        obs_model, obs_err = ps_observables(cvis, obs)
        obs_data = observe_data(data_c)

        nv2, nt3 = obs.nv2, obs.nt3amp
        resid = obs_model .- obs_data
        @views resid[nv2+nt3+1:end] .= mod360.(resid[nv2+nt3+1:end])

        chi2_v2  = nv2 > 0 ? sum((resid[1:nv2] ./ obs_err[1:nv2]) .^ 2) : 0.0
        chi2_t3a = nt3 > 0 ? sum((resid[nv2+1:nv2+nt3] ./ obs_err[nv2+1:nv2+nt3]) .^ 2) : 0.0
        chi2_t3p = nt3 > 0 ? sum((resid[nv2+nt3+1:end] ./ obs_err[nv2+nt3+1:end]) .^ 2) : 0.0

        n_c = nv2 + 2 * nt3
        chi2_c = weights[1] * chi2_v2 + weights[2] * chi2_t3a + weights[3] * chi2_t3p
        chi2_total += chi2_c
        n_total += n_c

        if verb
            @printf("  ch%d: χ²_v2=%.2f/%d  χ²_t3a=%.2f/%d  χ²_t3p=%.2f/%d\n",
                    c, chi2_v2, nv2, chi2_t3a, nt3, chi2_t3p, nt3)
        end
    end

    if verb
        @printf("  Total weighted χ² = %.2f / %d  (χ²_r = %.3f)\n",
                chi2_total, n_total, chi2_total / max(n_total, 1))
    end
    return chi2_total
end


# ── Main reconstruction ─────────────────────────────────────────────────────

"""
    reconstruct_pointsource(ps, data; kwargs...)
        -> (z, posterior, samples)

Full point-source reconstruction using hybrid MGVI + GeoVI. Builds a
`PointSourceProblem` and dispatches to the generic VI helpers defined in
geovi.jl. The MAP warm-start, per-iteration chi² + parameter reporting,
and the final aggregation into a `PointSourcePosterior` stay here because
they are point-source specific.
"""
function reconstruct_pointsource(ps::PointSourceParams, data;
                                  weights=[1.0, 1.0, 1.0],
                                  n_mgvi=6,
                                  n_geovi=4,
                                  n_samples=iter -> iter <= 2 ? 2 : 4,
                                  map_maxiter=200,
                                  kl_maxiter=35,
                                  kl_absdelta=0.5,
                                  cg_maxiter=100,
                                  cg_tol=0.01,
                                  geo_newton_maxiter=10,
                                  geo_cg_maxiter=50,
                                  geo_tol=1e-5,
                                  z0::Union{Nothing, Vector{Float64}}=nothing,
                                  verb=true)
    n = ps_latent_size(ps)
    prob = PointSourceProblem(ps, data; weights=weights)

    # ── MAP warm-start ──
    verb && println("=== Point-source MAP warm-start ===")
    z_init = z0 !== nothing ? Vector{Float64}(z0) : 0.1 .* randn(n)
    z_map = reconstruct_map(prob; z0=z_init, maxiter=map_maxiter, verb=verb)
    if verb
        ps_report_chi2(z_map, ps, data; weights=weights)
        x_r, y_r, fl = ps_unpack(z_map, ps)
        _ps_print_params(x_r, y_r, fl, ps)
    end

    # Per-iteration point-source-specific reporting (called by the generic
    # reconstruct_hybrid after each KL step).
    cb = (_prob, z, _samples, _iter) -> begin
        verb || return nothing
        ps_report_chi2(z, ps, data; weights=weights)
        x_r, y_r, fl = ps_unpack(z, ps)
        _ps_print_params(x_r, y_r, fl, ps)
        return nothing
    end

    # The default sample_mode in reconstruct_hybrid is exactly this:
    #   :linear_resample for iter ≤ n_mgvi, :nonlinear_update afterwards.
    z, samples = reconstruct_hybrid(prob;
                                     z0=z_map,
                                     n_mgvi=n_mgvi, n_geovi=n_geovi,
                                     n_samples=n_samples,
                                     kl_maxiter=kl_maxiter,
                                     kl_absdelta=kl_absdelta,
                                     cg_maxiter=cg_maxiter,
                                     cg_tol=cg_tol,
                                     geo_newton_maxiter=geo_newton_maxiter,
                                     geo_cg_maxiter=geo_cg_maxiter,
                                     geo_tol=geo_tol,
                                     iter_callback=cb, verb=verb)

    # ── Posterior aggregation (point-source physical-space) ──
    verb && println("\n=== Computing posterior statistics ===")
    posterior = _ps_posterior_stats(z, samples, ps)

    if verb
        @printf("\n=== Point Source Posterior ===\n")
        for i in 1:ps.N
            @printf("Source %d: x=%+.4f±%.4f  y=%+.4f±%.4f mas\n",
                    i, posterior.x_mas[i], posterior.x_std[i],
                    posterior.y_mas[i], posterior.y_std[i])
            for c in 1:ps.nf
                @printf("  ch%d (%.3f μm): f=%.4f±%.4f\n",
                        c, 3e8 / ps.freq[c] * 1e6,
                        posterior.flux[i, c], posterior.flux_std[i, c])
            end
        end
    end

    return z, posterior, samples
end

# Aggregate posterior x/y (in mas) and per-channel flux fractions from a
# sample set whose entries are already explicit antithetic pairs (the
# convention used by reconstruct_hybrid).
function _ps_posterior_stats(z::AbstractVector{Float64},
                              samples::Vector{Vector{Float64}},
                              ps::PointSourceParams)
    N, nf = ps.N, ps.nf
    n_samp = length(samples)
    x_sum  = zeros(N);  x_sum2 = zeros(N)
    y_sum  = zeros(N);  y_sum2 = zeros(N)
    f_sum  = zeros(N, nf); f_sum2 = zeros(N, nf)

    for k in 1:n_samp
        z_k = z .+ samples[k]
        x_r, y_r, fl = ps_unpack(z_k, ps)
        x_mas = x_r ./ MAS2RAD
        y_mas = y_r ./ MAS2RAD
        x_sum  .+= x_mas;  x_sum2 .+= x_mas .^ 2
        y_sum  .+= y_mas;  y_sum2 .+= y_mas .^ 2
        for c in 1:nf
            ftot = sum(fl[:, c])
            f_sum[:, c]  .+= fl[:, c] ./ ftot
            f_sum2[:, c] .+= (fl[:, c] ./ ftot) .^ 2
        end
    end

    x_mean = x_sum ./ n_samp
    y_mean = y_sum ./ n_samp
    x_std  = sqrt.(max.(x_sum2 ./ n_samp .- x_mean .^ 2, 0.0))
    y_std  = sqrt.(max.(y_sum2 ./ n_samp .- y_mean .^ 2, 0.0))
    f_mean = f_sum ./ n_samp
    f_std  = sqrt.(max.(f_sum2 ./ n_samp .- f_mean .^ 2, 0.0))

    return PointSourcePosterior(x_mean, y_mean, x_std, y_std, f_mean, f_std)
end


# ── Internal helpers ─────────────────────────────────────────────────────────

function _ps_print_params(x_rad, y_rad, flux, ps)
    N, nf = ps.N, ps.nf
    x_mas = x_rad ./ MAS2RAD
    y_mas = y_rad ./ MAS2RAD
    for i in 1:N
        ftot_str = join([@sprintf("%.3f", flux[i, c] / sum(flux[:, c]))
                         for c in 1:nf], ", ")
        @printf("  src%d: (%+.3f, %+.3f) mas  f=[%s]\n",
                i, x_mas[i], y_mas[i], ftot_str)
    end
end
