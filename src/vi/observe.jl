# ============================================================================
# Observation operator: image -> visibilities -> observables, and its adjoint and JVP.
# ============================================================================

using NFFT

"""
    ObservationConfig

Holds the NFFT plans and data indices needed to compute observables
from an image, in a form suitable for the differentiable forward model.

For hybrid models (image + parametric components), also stores UV coordinates
and an optional `FlatModel` for parametric visibility evaluation.
"""
struct ObservationConfig
    ft_uv::NFFTPlan       # all UV points
    ft_v2::NFFTPlan       # V2 baselines
    ft_t3_1::NFFTPlan     # T3 triangle leg 1
    ft_t3_2::NFFTPlan     # T3 triangle leg 2
    ft_t3_3::NFFTPlan     # T3 triangle leg 3
    indx_v2::Vector{Int64}
    indx_t3_1::Vector{Int64}
    indx_t3_2::Vector{Int64}
    indx_t3_3::Vector{Int64}
    v2_err::Vector{Float64}
    t3amp_err::Vector{Float64}
    t3phi_err::Vector{Float64}
    nv2::Int
    nt3amp::Int
    nt3phi::Int
    uv::Matrix{Float64}               # (2, nuv) — UV coordinates for eval_model
    model::Union{Nothing, FlatModel}   # optional parametric component
end

"""
    ObservationConfig(ft, data; model=nothing) -> ObservationConfig

Construct from OITOOLS ft plans and OIdata.

Accepts either:
- `ft::Vector{NFFTPlan}` + single OIdata  (legacy `setup_nfft` output)
- `ft::Matrix` + `data::Matrix{OIdata}`   (`setup_ft`/`readoifits` output — extracts [1,1])

Pass `model=FlatModel(...)` to enable hybrid image + parametric reconstruction.
"""
# Named fields, not `ft[1], ft[3], ft[4]…`. `NFFTCell` is an `AbstractVector` whose order
# happens to be (uv, vis, v2, t3_1, t3_2, t3_3), so the positional form worked — and hid the
# fact that the type had changed underneath it entirely. A reordering of the struct would have
# silently swapped V2 for a T3 leg.
function ObservationConfig(ft::NFFTCell, data;
                           model::Union{Nothing, FlatModel}=nothing)
    ObservationConfig(
        ft.uv, ft.v2, ft.t3_1, ft.t3_2, ft.t3_3,
        data.indx_v2, data.indx_t3_1, data.indx_t3_2, data.indx_t3_3,
        data.v2_err, data.t3amp_err, data.t3phi_err,
        data.nv2, data.nt3amp, data.nt3phi,
        Matrix{Float64}(data.uv),
        model
    )
end

# Matrix wrapper: extract [1,1] element from setup_ft / readoifits output
function ObservationConfig(ft::AbstractMatrix, data::AbstractMatrix;
                           model::Union{Nothing, FlatModel}=nothing)
    ObservationConfig(ft[1,1], data[1,1]; model=model)
end


# ============================================================================
# Internal: compute combined cvis from image (+ optional parametric model)
# ============================================================================

"""
    plan_precision(obs) -> Type

Element type of the NFFT plans this configuration was built on.

PRECISION HERE IS TWO NUMBERS, NOT ONE, and conflating them is a 20-60% error rather than a
rounding difference. The plans may be Float32 — `readoifits` and `setup_ft` default to it, and
at Float32 the transform is twice as fast — while the ADJOINT SOURCE must be accumulated at
Float64. A single uv point collects contributions from several observables (it appears in V²
and in all three legs of a T3) which are large and of opposite sign; summing them at Float32
loses the cancellation. Measured against a Float64 reference, gradient error at nx = 512:
0.57 with a Float32 accumulator, 8.0e-04 with a Float64 one.

So there are exactly three boundaries, and each converts once:

  * `_combined_cvis` and `observe_jvp` — image enters the plan at the plan's precision;
  * `_obs_to_g_cvis` — the cotangent is accumulated at Float64, whatever the plan is;
  * `_g_cvis_to_g_image` — that accumulator is rounded to the plan's precision once, at the
    adjoint.

This is the split `ImageChi2Cache{T,W}` makes in the image kernel, with T the plan and W the
accumulation. `scatter_obs_cotangent!` already accepts the mixed pair and floors its divisions
at `eps` of the VISIBILITIES, which is the right one of the two.
"""
plan_precision(obs::ObservationConfig) = ft_eltype(obs.ft_uv)

"""
    _combined_cvis(image, obs; params=nothing) -> cvis

Compute complex visibilities from image, optionally adding parametric model
contribution with SPARCO flux weighting.
"""
function _combined_cvis(image::AbstractMatrix{<:Real}, obs::ObservationConfig;
                        params::Union{Nothing, AbstractVector}=nothing)
    flux = sum(image)
    img_norm = Complex{plan_precision(obs)}.(image ./ flux)
    cvis_image = obs.ft_uv * img_norm

    if obs.model !== nothing && params !== nothing
        cvis_model = eval_model(obs.model, params, obs.uv)
        # SPARCO: f_env = 1 - total_param_flux (V(0)=1 for all analytic shapes)
        f_env = 1.0 - real(eval_model(obs.model, params, zeros(2, 1))[1])
        return f_env .* cvis_image .+ cvis_model
    else
        return cvis_image
    end
end


# ============================================================================
# Forward: image (+ params) → observables
# ============================================================================

"""
    observe(image, obs; params=nothing) -> (obs_model, obs_err)

Compute model observables (V², T3amp, T3phi) from an image.

For hybrid models, pass `params` — the parametric model parameters. The combined
visibility is `V_total = f_env · V_image + V_model` (SPARCO convention).
"""
function observe(image::AbstractMatrix{<:Real}, obs::ObservationConfig;
                 params::Union{Nothing, AbstractVector}=nothing)
    cvis = _combined_cvis(image, obs; params=params)

    # V²
    v2_model = abs2.(cvis[obs.indx_v2])

    # T3
    t3 = cvis[obs.indx_t3_1] .* cvis[obs.indx_t3_2] .* cvis[obs.indx_t3_3]
    t3amp_model = abs.(t3)
    t3phi_model = angle.(t3) .* (180.0 / π)  # degrees

    obs_model = vcat(v2_model, t3amp_model, t3phi_model)
    obs_err = vcat(obs.v2_err, obs.t3amp_err, obs.t3phi_err)
    return obs_model, obs_err
end

"""
    observe_data(data) -> Vector{Float64}

Extract the observed data vector (V², T3amp, T3phi) from OIdata,
in the same order as `observe` returns model observables.
"""
function observe_data(data)
    return vcat(data.v2, data.t3amp, data.t3phi)
end

"""
    n_model_params(obs) -> Int

Number of free parametric model parameters (0 if no model).
"""
n_model_params(obs::ObservationConfig) =
    obs.model === nothing ? 0 : length(obs.model.list_free_params)
n_model_params(obs_vec::Vector{ObservationConfig}) = n_model_params(obs_vec[1])


# ============================================================================
# ObsContext: opaque bundle of all observation data for GeoVI
# ============================================================================

"""
    ObsContext

Bundles per-channel observation configs, the concatenated error vector,
and optional cross-channel differential phase config into a single
opaque struct that GeoVI algorithm functions pass through without inspecting.

Only the 3 core GeoVI metric functions (transformation, right/left sqrt metric)
unpack this struct. All other GeoVI functions treat it as opaque.
"""
struct ObsContext
    obs_vec::Vector{ObservationConfig}
    sigma_vec::Vector{Float64}
    dp::Union{Nothing, DiffPhaseConfig}
end

n_model_params(ctx::ObsContext) = n_model_params(ctx.obs_vec[1])


# ============================================================================
# Internal: observable cotangent → cvis cotangent (shared by all adjoint paths)
# ============================================================================

"""
    _obs_to_g_cvis(g_obs, cvis, obs) -> g_cvis

Backpropagate observable-space cotangent to cvis-space cotangent.
This is the adjoint of the cvis → (V², T3amp, T3phi) extraction.
"""
function _obs_to_g_cvis(g_obs::AbstractVector{<:Real},
                        cvis::AbstractVector{<:Complex},
                        obs::ObservationConfig)
    nv2 = obs.nv2
    nt3 = obs.nt3amp
    # ComplexF64 whatever `cvis` is: this is the accumulator, and following the plan's
    # precision here is exactly the mistake `plan_precision` documents. It costs 16 bytes per
    # uv point -- under a megabyte on the largest set here -- and no arithmetic.
    g_cvis = zeros(ComplexF64, length(cvis))

    # OITOOLS' scatter, not a second copy of it. Two conventions have to be bridged here, and
    # both are real differences rather than taste:
    #
    #   CONJUGATION. OITOOLS defines g_cvis so that dL/dz = real(transpose(A) * g_cvis) for
    #   V = A z, which puts conj(V) in the V² term and -im on the phases. This file grew the
    #   conjugate of that — V and +im — so the result is conjugated on the way out. Neither is
    #   wrong; what was wrong was having both spellings of the same derivative in two packages,
    #   free to drift apart by exactly this sign.
    #
    #   DEGREES vs RADIANS. `g_obs`'s phase block is a cotangent per DEGREE, which is what the
    #   old inline code's stray 180/pi was doing. The scatter takes radians, so the conversion
    #   becomes an explicit scale instead of a bare constant in the middle of an expression.
    scatter_obs_cotangent!(g_cvis, cvis, obs;
        g_v2     = view(g_obs, 1:nv2),
        g_t3amp  = view(g_obs, nv2+1:nv2+nt3),
        g_t3phi  = view(g_obs, nv2+nt3+1:nv2+2*nt3),
        scale_t3phi = 180.0 / π)
    return conj!(g_cvis)
end


# ============================================================================
# Internal: cvis cotangent → image gradient (NFFT adjoint + flux normalization)
# ============================================================================

"""
    _g_cvis_to_g_image(g_cvis, image, obs; scale=1.0) -> g_image

Backpropagate cvis cotangent to image gradient via NFFT adjoint,
including flux normalization correction. Optionally scale by `scale`
(used for SPARCO f_env weighting).
"""
function _g_cvis_to_g_image(g_cvis::AbstractVector{<:Complex},
                            image::AbstractMatrix{<:Real},
                            obs::ObservationConfig;
                            scale::Float64=1.0)
    flux = sum(image)
    inv_flux = 1.0 / flux
    # The one rounding of the accumulator, here at the plan boundary. Everything after it is
    # Float64 again, since `inv_flux` and `scale` are.
    src = convert(Vector{Complex{plan_precision(obs)}}, g_cvis)
    g_img_norm = scale .* real.(adjoint(obs.ft_uv) * src)
    img_norm_real = image .* inv_flux
    correction = sum(g_img_norm .* img_norm_real)
    return (g_img_norm .- correction) .* inv_flux
end


# ============================================================================
# Adjoint (VJP): observable cotangent → image (+ params) gradient
# ============================================================================

"""
    observe_adjoint(g_obs, image, obs; params=nothing) -> g_image  OR  (g_image, g_params)

Adjoint (VJP) of the observe operator: maps a cotangent vector in observable-space
back to cotangent vectors in image-space (and optionally parameter-space).

When `params=nothing` (image-only mode), returns `g_image::Matrix{Float64}`.
When `params` is provided (hybrid mode), returns `(g_image, g_params)`.

Used by GeoVI's `left_sqrt_metric` to compute J_T' · v.
"""
function observe_adjoint(g_obs::AbstractVector{Float64},
                         image::AbstractMatrix{Float64},
                         obs::ObservationConfig;
                         params::Union{Nothing, AbstractVector{Float64}}=nothing)
    # Forward pass to get combined cvis
    cvis = _combined_cvis(image, obs; params=params)

    # Observable → cvis adjoint (shared logic)
    g_cvis = _obs_to_g_cvis(g_obs, cvis, obs)

    if obs.model !== nothing && params !== nothing
        # Hybrid mode: split g_cvis back to image and params

        # Compute quantities needed for param gradient
        wl = nothing  # monochromatic for now
        _, J = eval_model_grad(obs.model, params, obs.uv; wl=wl)
        uv_zero = zeros(2, 1)
        cvis_zero, J_zero = eval_model_grad(obs.model, params, uv_zero; wl=wl)
        f_env = 1.0 - real(cvis_zero[1])

        # Image gradient: f_env * NFFT adjoint with flux normalization
        g_image = _g_cvis_to_g_image(g_cvis, image, obs; scale=f_env)

        # Param gradient: direct model contribution
        # Contract: dL = Re(dot(g_cvis, d_cvis)) = Re(g_cvis^H · d_cvis)
        # For d_cvis = J·dp: g_params = Re(J^H · g_cvis)
        g_params = real.(adjoint(J) * g_cvis)

        # f_env chain rule: ∂loss/∂f_env · ∂f_env/∂params
        # d_cvis from f_env = cvis_image · df_env, so g_f_env = Re(g_cvis^H · cvis_image)
        flux = sum(image)
        cvis_image = obs.ft_uv * Complex{Float64}.(image ./ flux)
        g_f_env = real(dot(g_cvis, cvis_image))
        df_env_dx = -real.(vec(J_zero[1, :]))
        g_params .+= g_f_env .* df_env_dx

        return g_image, g_params
    else
        # Image-only mode (original behavior)
        return _g_cvis_to_g_image(g_cvis, image, obs)
    end
end


# ── Jacobian-vector products ────────────────────────────────────────────────
#
# Kept beside the forward and adjoint they differentiate, as `diffphase.jl` already does.
# Split into a file of their own they were easy to forget when an operator changed.

"""
    sky_forward_jvp(xi, v, p) -> (image, d_image)

Compute sky_forward(xi, p) and its JVP d(sky_forward)/d(xi) * v simultaneously.
"""
function sky_forward_jvp(xi::AbstractVector{<:Real}, v::AbstractVector{<:Real},
                         p::SkyModelParams)
    xi_sp, xi_sc, sp, sc, x_alpha, x_logF0 = _unpack(xi, p)
    v_sp, v_sc, vsp, vsc, v_alpha, v_logF0 = _unpack(v, p)
    sp_slope, sp_fluct, sp_flex, sp_asp, sp_spec = sp
    sc_slope, sc_fluct, sc_flex, sc_asp, sc_spec = sc
    vsp_slope, vsp_fluct, vsp_flex, vsp_asp, vsp_spec = vsp
    vsc_slope, vsc_fluct, vsc_flex, vsc_asp, vsc_spec = vsc
    nf = length(p.freq)

    # Spatial: amplitude spectrum + tangent
    amp_sp = amplitude_spectrum(sp_slope, sp_fluct, sp_flex, sp_asp, sp_spec,
                                p.spatial, p)
    d_amp_sp = amplitude_spectrum_jvp(vsp_slope, vsp_fluct, vsp_flex, vsp_asp, vsp_spec,
                                      sp_slope, sp_fluct, sp_flex, sp_asp, sp_spec,
                                      p.spatial, p)

    kernel_sp   = amp_sp[p.bin_index]
    d_kernel_sp = d_amp_sp[p.bin_index]

    # Spatial correlated field with product rule — transient FFTs via scratch.
    # Hold F_xi (cs2) and F_v (cs3); compute d_tau0_raw first (needs both), then
    # tau0_raw (may then overwrite F_xi).
    p.cs1 .= xi_sp;  mul!(p.cs2, p.P, p.cs1)     # F_xi → cs2
    p.cs1 .= v_sp;   mul!(p.cs3, p.P, p.cs1)     # F_v  → cs3
    @. p.cs3 = p.cs3 * kernel_sp + p.cs2 * d_kernel_sp
    mul!(p.cs1, p.iP, p.cs3)
    d_tau0_raw = real.(p.cs1)
    p.cs2 .*= kernel_sp
    mul!(p.cs1, p.iP, p.cs2)
    tau0_raw = real.(p.cs1)

    # Zero-mean projection (linear)
    wmean_val   = sum(p.D .* tau0_raw) / p.D_sum
    d_wmean_val = sum(p.D .* d_tau0_raw) / p.D_sum
    tau0   = p.D .* (tau0_raw .- wmean_val)
    d_tau0 = p.D .* (d_tau0_raw .- d_wmean_val)

    # Spectral: amplitude spectrum + tangent
    # DC mode is zero — overall flux is handled by logF0 + alpha
    amp_sc = amplitude_spectrum(sc_slope, sc_fluct, sc_flex, sc_asp, sc_spec,
                                p.spectral, p)
    d_amp_sc = amplitude_spectrum_jvp(vsc_slope, vsc_fluct, vsc_flex, vsc_asp, vsc_spec,
                                      sc_slope, sc_fluct, sc_flex, sc_asp, sc_spec,
                                      p.spectral, p)
    amp_sc[1] = 0.0
    d_amp_sc[1] = 0.0

    kernel_sc   = amp_sc[p.bin_index]
    d_kernel_sc = d_amp_sc[p.bin_index]

    # Per-channel spectral field with product rule
    delta_tau   = Array{Float64}(undef, p.npix, p.npix, nf)
    d_delta_tau = Array{Float64}(undef, p.npix, p.npix, nf)
    for c in 1:nf
        p.cs1 .= view(xi_sc, :, :, c);  mul!(p.cs2, p.P, p.cs1)   # F_c  → cs2
        p.cs1 .= view(v_sc, :, :, c);   mul!(p.cs3, p.P, p.cs1)   # dF_c → cs3
        @. p.cs3 = p.cs3 * kernel_sc + p.cs2 * d_kernel_sc
        mul!(p.cs1, p.iP, p.cs3)
        @. d_delta_tau[:, :, c] = real(p.cs1)
        p.cs2 .*= kernel_sc
        mul!(p.cs1, p.iP, p.cs2)
        @. delta_tau[:, :, c] = real(p.cs1)
    end

    # Build log-intensity and tangent
    Y   = Array{Float64}(undef, p.npix, p.npix, nf)
    d_Y = Array{Float64}(undef, p.npix, p.npix, nf)
    for c in 1:nf
        @. Y[:, :, c]  = x_logF0 + tau0 + delta_tau[:, :, c] +
                          x_alpha * p.log_freq_ratio[c]
        @. d_Y[:, :, c] = v_logF0 + d_tau0 + d_delta_tau[:, :, c] +
                           v_alpha * p.log_freq_ratio[c]
    end

    clamp_mask = (Y .> -30.0) .& (Y .< 30.0)
    @. Y = clamp(Y, -30.0, 30.0)
    d_Y .*= clamp_mask

    image   = Array{Float64}(undef, p.npix, p.npix, nf)
    d_image = Array{Float64}(undef, p.npix, p.npix, nf)
    for c in 1:nf
        @. image[:, :, c]   = p.D * exp(Y[:, :, c])
        @. d_image[:, :, c] = image[:, :, c] * d_Y[:, :, c]
    end

    return image, d_image
end


# ============================================================================
# Internal: cvis + d_cvis → observable tangents
# ============================================================================

"""
    _cvis_to_obs_jvp(cvis, d_cvis, obs) -> (obs_model, obs_err, d_obs_model)

Compute observables and their tangents from combined complex visibilities.
Shared by both image-only and hybrid observe_jvp.
"""
function _cvis_to_obs_jvp(cvis::AbstractVector{<:Complex},
                          d_cvis::AbstractVector{<:Complex},
                          obs::ObservationConfig)
    cv2 = cvis[obs.indx_v2]
    dv2 = d_cvis[obs.indx_v2]
    v2_model = abs2.(cv2)
    d_v2 = 2.0 .* real.(conj.(cv2) .* dv2)

    z1 = cvis[obs.indx_t3_1]
    z2 = cvis[obs.indx_t3_2]
    z3 = cvis[obs.indx_t3_3]
    dz1 = d_cvis[obs.indx_t3_1]
    dz2 = d_cvis[obs.indx_t3_2]
    dz3 = d_cvis[obs.indx_t3_3]

    t3 = z1 .* z2 .* z3
    d_t3 = dz1 .* z2 .* z3 .+
           z1 .* dz2 .* z3 .+
           z1 .* z2 .* dz3

    t3amp_model = abs.(t3)
    t3amp_safe = t3amp_model .+ 1e-30
    d_t3amp = real.(conj.(t3) .* d_t3) ./ t3amp_safe

    t3abs2 = abs2.(t3) .+ 1e-60
    t3phi_model = angle.(t3) .* (180.0 / π)
    d_t3phi = imag.(conj.(t3) .* d_t3) ./ t3abs2 .* (180.0 / π)

    obs_model = vcat(v2_model, t3amp_model, t3phi_model)
    obs_err = vcat(obs.v2_err, obs.t3amp_err, obs.t3phi_err)
    d_obs_model = vcat(d_v2, d_t3amp, d_t3phi)

    return obs_model, obs_err, d_obs_model
end


# ============================================================================
# observe JVP: image (+ params) → observables with tangents
# ============================================================================

"""
    observe_jvp(image, d_image, obs; params=nothing, d_params=nothing)
        -> (obs_model, obs_err, d_obs_model)

Compute observe(image, obs) and its JVP simultaneously.

For hybrid models, pass `params` and `d_params` to include parametric model
tangent contributions. The combined visibility tangent follows the SPARCO formula:
`d_V_total = f_env · d_V_image + d_f_env · V_image + d_V_model`.
"""
function observe_jvp(image::AbstractMatrix{Float64}, d_image::AbstractMatrix{Float64},
                     obs::ObservationConfig;
                     params::Union{Nothing, AbstractVector{Float64}}=nothing,
                     d_params::Union{Nothing, AbstractVector{Float64}}=nothing)
    # Image cvis + tangent (flux normalization)
    flux = sum(image)
    d_flux = sum(d_image)
    inv_flux = 1.0 / flux
    img_norm = image .* inv_flux
    d_img_norm = (d_image .* flux .- image .* d_flux) .* (inv_flux * inv_flux)

    cvis_image   = obs.ft_uv * Complex{plan_precision(obs)}.(img_norm)
    d_cvis_image = obs.ft_uv * Complex{plan_precision(obs)}.(d_img_norm)

    if obs.model !== nothing && params !== nothing && d_params !== nothing
        # Hybrid: add parametric model contribution
        wl = nothing  # monochromatic for now
        cvis_model, J = eval_model_grad(obs.model, params, obs.uv; wl=wl)
        d_cvis_model = J * d_params

        uv_zero = zeros(2, 1)
        cvis_zero, J_zero = eval_model_grad(obs.model, params, uv_zero; wl=wl)
        f_env = 1.0 - real(cvis_zero[1])
        d_f_env = -real(dot(J_zero[1, :], d_params))

        cvis = f_env .* cvis_image .+ cvis_model
        d_cvis = f_env .* d_cvis_image .+ d_f_env .* cvis_image .+ d_cvis_model
    else
        cvis = cvis_image
        d_cvis = d_cvis_image
    end

    return _cvis_to_obs_jvp(cvis, d_cvis, obs)
end
