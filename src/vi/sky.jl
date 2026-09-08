# ============================================================================
# Sky model: correlated-field image with a limb-darkened disc weight.
# ============================================================================

const MAS2RAD = π / (180.0 * 3600.0 * 1000.0)

"""
    _check_weight(w, npix) -> Matrix{Float64}

A caller-supplied support map, checked and rescaled the way `limb_weight` returns one.

Every requirement here is one the model cannot survive without. A negative weight would make
`D .* exp(field)` a negative intensity; an all-zero map divides the zero-mean projection by
`D_sum = 1e-12`; a wrong shape reaches the FFTs as a `DimensionMismatch` several calls later,
where it says nothing about the file that caused it.
"""
function _check_weight(w::AbstractMatrix{<:Real}, npix::Int)
    size(w) == (npix, npix) || throw(ArgumentError(
        "the sky-prior weight is $(size(w,1))×$(size(w,2)) but npix is $npix"))
    all(isfinite, w) || throw(ArgumentError("the sky-prior weight has non-finite pixels"))
    any(<(0), w) && throw(ArgumentError(
        "the sky-prior weight has negative pixels; it multiplies an intensity"))
    m = maximum(w)
    m > 0 || throw(ArgumentError("the sky-prior weight is everywhere zero"))
    return Float64.(w) ./ m
end

"""
    limb_weight(nx, ny, R_pix, u) -> Matrix{Float64}

Limb-darkened disk weight in pixel space.
`R_pix` = radius in pixels, `u` = limb-darkening coefficient.
"""
function limb_weight(nx::Int, ny::Int, R_pix::Float64; u::Float64=0.5)
    cx, cy = (nx - 1) / 2.0, (ny - 1) / 2.0
    w = zeros(nx, ny)
    for j in 1:ny, i in 1:nx
        xc = (j - 1 - cy) / R_pix
        yc = (i - 1 - cx) / R_pix
        r2 = xc^2 + yc^2
        if r2 <= 1.0
            mu = sqrt(max(1.0 - r2, 0.0))
            w[i, j] = (1.0 - u) + u * mu
        end
    end
    m = maximum(w)
    return m > 0 ? w ./ m : w
end

"""
    SkyModelParams

Fixed hyperparameters for the sky model (not inferred).
Both spatial and spectral fields use NIFTy-style correlated fields
with learned power spectra.
"""
struct SkyModelParams{TP,TIP}
    grid::FourierGridInfo  # Fourier-grid data (npix, dx, bin_index, mode_*, n_bins, ...)
    freq::Vector{Float64}  # frequencies in Hz
    D::Matrix{Float64}     # limb-darkened disk weight (npix x npix)
    D_sum::Float64         # sum of D weights

    # Correlated field configurations
    spatial::CorrFieldConfig
    spectral::CorrFieldConfig

    log_freq_ratio::Vector{Float64}

    # Cached out-of-place FFT plans over the npix×npix image grid + three complex
    # scratch buffers. The spatial/spectral correlated-field transforms route
    # their transient FFTs through mul!(scratch, plan, ·) so those outputs
    # allocate nothing (three buffers: the JVP/adjoint hold two spectra at once).
    # Safe because inference runs single-threaded. Parameterised on the plan
    # types so the cached plans stay type-stable.
    P::TP
    iP::TIP
    cs1::Matrix{ComplexF64}
    cs2::Matrix{ComplexF64}
    cs3::Matrix{ComplexF64}
end

# Forward old top-level grid field names to the embedded FourierGridInfo so
# existing call sites (`p.npix`, `p.bin_index`, etc.) keep working. The
# amplitude_spectrum* family is duck-typed on these names and works on either
# a SkyModelParams (forwarded) or a FourierGridInfo (direct).
const _GRID_FORWARD_FIELDS = (:npix, :dx, :n_bins, :bin_index,
                              :mode_multiplicity, :mode_lengths,
                              :rel_log_mode_lengths, :log_volume,
                              :total_volume)

@inline function Base.getproperty(p::SkyModelParams, sym::Symbol)
    if sym in _GRID_FORWARD_FIELDS
        return getfield(getfield(p, :grid), sym)
    end
    return getfield(p, sym)
end

@inline function Base.propertynames(p::SkyModelParams, private::Bool=false)
    return (_GRID_FORWARD_FIELDS..., fieldnames(typeof(p))...)
end

"""
    SkyModelParams(npix, pixsize, freq; R_mas, u, weight,
                   spatial_slope_prior, spatial_fluct_prior, ...,
                   spectral_slope_prior, spectral_fluct_prior, ...)

Construct sky model parameters from physical hyperparameters.
`pixsize` is the pixel size in milliarcseconds (converted to radians internally).
The spatial and spectral fields each get their own correlated field config.
Set `spatial_flex_prior=nothing` to disable IWP.

The weight `D` is the image's SUPPORT and envelope, not a display choice: the model is
`D .* exp(field)`, so `D = 0` forces a pixel to zero flux and the zero-mean projection is taken
against `D` as well. `R_mas` and `u` build a limb-darkened disc for it, which is the default;
pass `weight` instead to supply any non-negative `npix × npix` map — a mask, a photosphere, a
previous reconstruction — and `R_mas` is then unused. Whatever is given is rescaled to a
maximum of 1, as `limb_weight` returns, so the flux parameters keep their meaning.
"""
function SkyModelParams(npix::Int, pixsize::Float64, freq::Vector{Float64};
                        R_mas::Float64=NaN, u::Float64=0.0,
                        weight::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
                        spatial_slope_prior::Tuple{Float64,Float64}=(-4.0, 1.0),
                        spatial_fluct_prior::Tuple{Float64,Float64}=(1.6487, 2.1612),  # value-space lognormal (≡ old log μ=0, σ=1)
                        spatial_flex_prior::Union{Nothing,Tuple{Float64,Float64}}=nothing,
                        spatial_asp_prior::Union{Nothing,Tuple{Float64,Float64}}=nothing,
                        spectral_slope_prior::Tuple{Float64,Float64}=(-4.0, 1.0),
                        spectral_fluct_prior::Tuple{Float64,Float64}=(0.03423, 0.01824),  # value-space lognormal (≡ old log μ=-3.5, σ=0.5)
                        spectral_flex_prior::Union{Nothing,Tuple{Float64,Float64}}=nothing,
                        spectral_asp_prior::Union{Nothing,Tuple{Float64,Float64}}=nothing,
                        # Deprecated aliases (will be removed)
                        slope_prior=nothing, fluct_prior=nothing,
                        flex_prior=nothing, asp_prior=nothing,
                        ch_slope_prior=nothing, ch_fluct_prior=nothing,
                        ch_flex_prior=nothing, ch_asp_prior=nothing)
    dx = pixsize * MAS2RAD

    # Support old kwarg names as fallback
    sp_slope = something(slope_prior, spatial_slope_prior)
    sp_fluct = something(fluct_prior, spatial_fluct_prior)
    sp_flex  = flex_prior !== nothing ? flex_prior : spatial_flex_prior
    sp_asp   = asp_prior  !== nothing ? asp_prior  : spatial_asp_prior
    sc_slope = something(ch_slope_prior, spectral_slope_prior)
    sc_fluct = something(ch_fluct_prior, spectral_fluct_prior)
    sc_flex  = ch_flex_prior !== nothing ? ch_flex_prior : spectral_flex_prior
    sc_asp   = ch_asp_prior  !== nothing ? ch_asp_prior  : spectral_asp_prior

    D = if weight === nothing
        isfinite(R_mas) || throw(ArgumentError(
            "SkyModelParams needs either R_mas (a limb-darkened disc) or weight (a map)"))
        limb_weight(npix, npix, R_mas / pixsize; u=u)
    else
        _check_weight(weight, npix)
    end
    D_sum = sum(D) + 1e-12

    grid = FourierGridInfo(npix, dx)

    sp_cfg = CorrFieldConfig(; slope_prior=sp_slope, fluct_prior=sp_fluct,
                               flex_prior=sp_flex, asp_prior=sp_asp)
    sc_cfg = CorrFieldConfig(; slope_prior=sc_slope, fluct_prior=sc_fluct,
                               flex_prior=sc_flex, asp_prior=sc_asp)

    lf = log.(freq)
    nu0 = exp(sum(lf) / length(lf))
    log_freq_ratio = log.(freq ./ nu0)

    buf = zeros(ComplexF64, npix, npix)        # ESTIMATE plans ⇒ array-agnostic
    P   = plan_fft(buf)
    iP  = plan_ifft(buf)
    return SkyModelParams(grid, freq, D, D_sum, sp_cfg, sc_cfg, log_freq_ratio,
                          P, iP, zeros(ComplexF64, npix, npix),
                          zeros(ComplexF64, npix, npix), zeros(ComplexF64, npix, npix))
end

# --- Latent vector layout ---
# xi[1 : n²]                    = xi_sp         (spatial white noise)
# xi[n²+1 : n²+n²*nf]          = xi_sc         (per-channel spectral white noise)
# [spatial hypers: xi_slope, xi_fluct, [xi_flex, xi_asp, xi_spectrum]]
# [spectral hypers: xi_slope, xi_fluct, [xi_flex, xi_asp, xi_spectrum]]
# xi[end-1]                     = x_alpha       (spectral index)
# xi[end]                       = x_logF0       (log flux)

function _latent_size(p::SkyModelParams)
    n2 = p.npix^2
    nf = length(p.freq)
    n_sp = _spectral_latent_size(p.spatial, p.n_bins)
    n_sc = _spectral_latent_size(p.spectral, p.n_bins)
    return n2 + n2 * nf + n_sp + n_sc + 2
end

"""
    frozen_spectral_range(p::SkyModelParams) -> Vector{UnitRange{Int}}

Return the latent index ranges that should be frozen for monochromatic data:
  1. spectral xi (degenerate with spatial xi when nf=1)
  2. spectral hyperparameters (slope, fluct — inert when spectral fluct≈0)

Spatial hyperparameters are NOT frozen — they remain free for inference.
"""
function frozen_spectral_range(p::SkyModelParams)
    n2 = p.npix^2
    nf = length(p.freq)
    n_sp = _spectral_latent_size(p.spatial, p.n_bins)
    n_sc = _spectral_latent_size(p.spectral, p.n_bins)
    # Layout: [xi_sp(n2) | xi_sc(n2*nf) | sp_hypers(n_sp) | sc_hypers(n_sc) | alpha | flux]
    spectral_xi_range   = (n2 + 1):(n2 + n2 * nf)
    spectral_hyp_range  = (n2 + n2 * nf + n_sp + 1):(n2 + n2 * nf + n_sp + n_sc)
    return [spectral_xi_range, spectral_hyp_range]
end

# Keep old name as deprecated alias
frozen_chromatic_range(p::SkyModelParams) = frozen_spectral_range(p)

function _unpack_spectral(xi::AbstractVector, offset::Int,
                          cfg::CorrFieldConfig, n_bins::Int)
    xi_slope = xi[offset + 1]
    xi_fluct = xi[offset + 2]
    if cfg.use_iwp
        xi_flex = xi[offset + 3]
        xi_asp  = xi[offset + 4]
        n_iwp_flat = 2 * (n_bins - 2)
        xi_spectrum = reshape(view(xi, offset+5:offset+4+n_iwp_flat), n_bins - 2, 2)
    else
        xi_flex = 0.0
        xi_asp  = 0.0
        xi_spectrum = zeros(0, 2)
    end
    return xi_slope, xi_fluct, xi_flex, xi_asp, xi_spectrum
end

function _unpack(xi::AbstractVector, p::SkyModelParams)
    n2 = p.npix^2
    nf = length(p.freq)
    xi_sp = reshape(view(xi, 1:n2), p.npix, p.npix)
    xi_sc = reshape(view(xi, n2+1:n2+n2*nf), p.npix, p.npix, nf)

    offset = n2 + n2 * nf
    sp = _unpack_spectral(xi, offset, p.spatial, p.n_bins)
    offset += _spectral_latent_size(p.spatial, p.n_bins)
    sc = _unpack_spectral(xi, offset, p.spectral, p.n_bins)

    x_alpha = xi[end-1]
    x_logF0 = xi[end]

    return xi_sp, xi_sc, sp, sc, x_alpha, x_logF0
end

"""
    report_latents(xi, p) -> nothing

Print inferred physical parameters at the center point.
"""
function report_latents(xi::AbstractVector, p::SkyModelParams)
    _, _, sp, sc, x_alpha, x_logF0 = _unpack(xi, p)
    sp_slope, sp_fluct, sp_flex, sp_asp, _ = sp
    sc_slope, sc_fluct, sc_flex, sc_asp, _ = sc

    slope_sp = p.spatial.slope_mean + p.spatial.slope_std * sp_slope
    fluct_sp = exp(p.spatial.fluct_mean + p.spatial.fluct_std * sp_fluct)
    slope_sc = p.spectral.slope_mean + p.spectral.slope_std * sc_slope
    fluct_sc = exp(p.spectral.fluct_mean + p.spectral.fluct_std * sc_fluct)

    @printf("    spatial:   slope=%.2f  fluct=%.3e", slope_sp, fluct_sp)
    if p.spatial.use_iwp
        flex_sp = exp(p.spatial.flex_mean + p.spatial.flex_std * sp_flex)
        asp_sp  = exp(p.spatial.asp_mean  + p.spatial.asp_std  * sp_asp)
        @printf("  flex=%.3e  asp=%.3e", flex_sp, asp_sp)
    end
    println()
    @printf("    spectral:  slope=%.2f  fluct=%.3e", slope_sc, fluct_sc)
    if p.spectral.use_iwp
        flex_sc = exp(p.spectral.flex_mean + p.spectral.flex_std * sc_flex)
        asp_sc  = exp(p.spectral.asp_mean  + p.spectral.asp_std  * sc_asp)
        @printf("  flex=%.3e  asp=%.3e", flex_sc, asp_sc)
    end
    println()
    @printf("    spectral index=%.3f  log flux=%.3f\n", x_alpha, x_logF0)
end

"""
    _prior_stats(group) -> (rchi2, mean)

Compute reduced chi-squared and mean for a latent variable group.
"""
function _prior_stats(group::AbstractArray)
    v = vec(group)
    ndof = length(v)
    if ndof == 0
        return 0.0, 0.0, 0
    end
    rchi2 = dot(v, v) / ndof
    avg = sum(v) / ndof
    return rchi2, avg, ndof
end

"""
    _latent_groups(xi, p) -> Vector{Tuple{String, AbstractArray}}

Return named latent variable groups for minisanity reporting.
Each group corresponds to a set of N(0,1) latent variables.
"""
function _latent_groups(xi::AbstractVector, p::SkyModelParams)
    groups = Tuple{String, AbstractArray}[]
    n2 = p.npix^2
    nf = length(p.freq)

    # Spatial Fourier excitation field (npix²)
    push!(groups, ("spatial xi", reshape(view(xi, 1:n2), p.npix, p.npix)))
    # Spectral Fourier excitation field (npix² × nfreq)
    push!(groups, ("spectral xi", reshape(view(xi, n2+1:n2+n2*nf), p.npix, p.npix, nf)))

    # Spatial power spectrum parameters
    offset = n2 + n2 * nf
    push!(groups, ("spatial slope", view(xi, offset+1:offset+1)))
    push!(groups, ("spatial fluct", view(xi, offset+2:offset+2)))
    if p.spatial.use_iwp
        push!(groups, ("spatial flex", view(xi, offset+3:offset+3)))
        push!(groups, ("spatial asperity", view(xi, offset+4:offset+4)))
        n_iwp = 2 * (p.n_bins - 2)
        push!(groups, ("spatial spectrum", view(xi, offset+5:offset+4+n_iwp)))
    end
    offset += _spectral_latent_size(p.spatial, p.n_bins)

    # Spectral power spectrum parameters
    push!(groups, ("spectral slope", view(xi, offset+1:offset+1)))
    push!(groups, ("spectral fluct", view(xi, offset+2:offset+2)))
    if p.spectral.use_iwp
        push!(groups, ("spectral flex", view(xi, offset+3:offset+3)))
        push!(groups, ("spectral asperity", view(xi, offset+4:offset+4)))
        n_iwp = 2 * (p.n_bins - 2)
        push!(groups, ("spectral spectrum", view(xi, offset+5:offset+4+n_iwp)))
    end

    # Global parameters
    push!(groups, ("spectral index", view(xi, length(xi)-1:length(xi)-1)))
    push!(groups, ("log flux", view(xi, length(xi):length(xi))))
    return groups
end

"""
    minisanity(center, samples, p) -> nothing

NIFTy-style prior residual report: reduced chi² and mean per latent group,
averaged over antithetic sample pairs (center ± δ_k).
Reports mean ± std across all sample evaluations.
"""
function minisanity(center::AbstractVector, samples::AbstractVector{<:AbstractVector},
                    p::SkyModelParams)
    # Determine which groups are frozen (prior std = 0)
    frozen_names = _frozen_group_names(p)
    DIM = "\033[90m"   # dark grey ANSI
    RST = "\033[0m"

    # Print inferred physical parameters at center
    _, _, sp, sc, x_alpha, x_logF0 = _unpack(center, p)
    sp_slope, sp_fluct, sp_flex, sp_asp, _ = sp
    sc_slope, sc_fluct, sc_flex, sc_asp, _ = sc

    slope_sp = p.spatial.slope_mean + p.spatial.slope_std * sp_slope
    fluct_sp = exp(p.spatial.fluct_mean + p.spatial.fluct_std * sp_fluct)
    slope_sc = p.spectral.slope_mean + p.spectral.slope_std * sc_slope
    fluct_sc = exp(p.spectral.fluct_mean + p.spectral.fluct_std * sc_fluct)

    sp_frozen = p.spatial.slope_std == 0.0 && p.spatial.fluct_std == 0.0
    sc_frozen = p.spectral.slope_std == 0.0 && p.spectral.fluct_std == 0.0

    @printf("    Inferred parameters (center):\n")
    c1 = sp_frozen ? DIM : ""
    c1r = sp_frozen ? RST : ""
    @printf("      %sspatial:   slope=%.2f  fluct=%.3e", c1, slope_sp, fluct_sp)
    if p.spatial.use_iwp
        flex_sp = exp(p.spatial.flex_mean + p.spatial.flex_std * sp_flex)
        asp_sp  = exp(p.spatial.asp_mean  + p.spatial.asp_std  * sp_asp)
        @printf("  flex=%.3e  asp=%.3e", flex_sp, asp_sp)
    end
    print(c1r)
    println()
    c2 = sc_frozen ? DIM : ""
    c2r = sc_frozen ? RST : ""
    @printf("      %sspectral:  slope=%.2f  fluct=%.3e", c2, slope_sc, fluct_sc)
    if p.spectral.use_iwp
        flex_sc = exp(p.spectral.flex_mean + p.spectral.flex_std * sc_flex)
        asp_sc  = exp(p.spectral.asp_mean  + p.spectral.asp_std  * sc_asp)
        @printf("  flex=%.3e  asp=%.3e", flex_sc, asp_sc)
    end
    print(c2r)
    println()
    @printf("      spectral index=%.3f  log flux=%.3f\n", x_alpha, x_logF0)

    # Collect stats per group across all sample evaluations
    groups_ref = _latent_groups(center, p)
    n_groups = length(groups_ref)
    group_names = [g[1] for g in groups_ref]
    group_ndofs = [length(vec(g[2])) for g in groups_ref]

    rchi2_all = [Float64[] for _ in 1:n_groups]
    mean_all  = [Float64[] for _ in 1:n_groups]

    for k in eachindex(samples)
        for s in (1.0, -1.0)
            xi_s = center .+ s .* samples[k]
            groups = _latent_groups(xi_s, p)
            for (ig, (_, grp)) in enumerate(groups)
                rchi2, avg, ndof = _prior_stats(grp)
                if ndof > 0
                    push!(rchi2_all[ig], rchi2)
                    push!(mean_all[ig], avg)
                end
            end
        end
    end

    println("    Prior residual(s):")
    for ig in 1:n_groups
        ndof = group_ndofs[ig]
        if ndof == 0
            continue
        end
        rc = rchi2_all[ig]
        mn = mean_all[ig]
        rc_mean = sum(rc) / length(rc)
        rc_std  = length(rc) > 1 ? std(rc; corrected=true) : 0.0
        mn_mean = sum(mn) / length(mn)
        mn_std  = length(mn) > 1 ? std(mn; corrected=true) : 0.0
        is_frozen = group_names[ig] in frozen_names
        pre  = is_frozen ? DIM : ""
        post = is_frozen ? RST : ""
        @printf("      %s%-16s  red.χ²=%7.2f±%5.2f  avg=%+8.3f±%5.3f  #dof=%d%s\n",
                pre, group_names[ig], rc_mean, rc_std, mn_mean, mn_std, ndof, post)
    end
end

"""
    _frozen_group_names(p) -> Set{String}

Return the set of latent group names whose prior std is zero (frozen).
"""
function _frozen_group_names(p::SkyModelParams)
    frozen = Set{String}()
    # Spatial xi is frozen when fluct ≈ 0 (no field amplitude)
    if p.spatial.fluct_std == 0.0
        push!(frozen, "spatial xi")
    end
    if p.spatial.slope_std == 0.0
        push!(frozen, "spatial slope")
    end
    if p.spatial.fluct_std == 0.0
        push!(frozen, "spatial fluct")
    end
    if p.spatial.use_iwp
        if p.spatial.flex_std == 0.0; push!(frozen, "spatial flex"); end
        if p.spatial.asp_std  == 0.0; push!(frozen, "spatial asperity"); end
    end
    # Spectral
    if p.spectral.fluct_std == 0.0
        push!(frozen, "spectral xi")
    end
    if p.spectral.slope_std == 0.0
        push!(frozen, "spectral slope")
    end
    if p.spectral.fluct_std == 0.0
        push!(frozen, "spectral fluct")
    end
    if p.spectral.use_iwp
        if p.spectral.flex_std == 0.0; push!(frozen, "spectral flex"); end
        if p.spectral.asp_std  == 0.0; push!(frozen, "spectral asperity"); end
    end
    return frozen
end

"""
    report_chi2(image, p, ft, data) -> nothing

Print per-observable reduced chi-squared breakdown.
Accepts both matrix types (`setup_ft`/`readoifits`) and legacy single-element types.
"""
function report_chi2(image::AbstractArray{Float64,3}, p::SkyModelParams, ft, data)
    nf = length(p.freq)
    # Unwrap matrix types: linear indexing into Matrix{OIdata} / Matrix{Vector{NFFTPlan}}
    _ft(c)   = ft   isa AbstractMatrix ? ft[c]   : (nf == 1 ? ft : ft[c])
    _data(c) = data isa AbstractMatrix ? data[c] : (nf == 1 ? data : data[c])

    if nf == 1
        d = _data(1)
        chi2_v2, chi2_t3a, chi2_t3p = _chi2_components(image[:,:,1], _ft(1), d)
        @printf("    Likelihood: V²=%.2f  T3amp=%.2f  T3phi=%.2f\n",
                chi2_v2/d.nv2, chi2_t3a/d.nt3amp, chi2_t3p/d.nt3phi)
    else
        chi2_total = 0.0
        ndof_total = 0
        for c in 1:nf
            d = _data(c)
            chi2_v2, chi2_t3a, chi2_t3p = _chi2_components(image[:,:,c], _ft(c), d)
            chi2_total += chi2_v2 + chi2_t3a + chi2_t3p
            ndof_total += d.nv2 + d.nt3amp + d.nt3phi
        end
        @printf("    Likelihood: total red.χ²=%.2f\n", chi2_total/ndof_total)
    end
end

"""
    _chi2_components(image, ft_plans, data) -> (chi2_v2, chi2_t3a, chi2_t3p)

Compute per-observable chi² for a single channel.
`ft_plans` is a Vector{NFFTPlan} (6 plans from setup_nfft); `ft_plans[1]` is the full UV plan.
"""
function _chi2_components(image::AbstractMatrix{Float64}, ft_plans, data)
    flux = sum(image)
    img_norm = image ./ flux
    cvis_all = ft_plans[1] * Complex{Float64}.(img_norm)
    v2_model = abs2.(cvis_all[data.indx_v2])
    chi2_v2 = sum(((v2_model .- data.v2) ./ data.v2_err) .^ 2)
    t3 = cvis_all[data.indx_t3_1] .* cvis_all[data.indx_t3_2] .* cvis_all[data.indx_t3_3]
    chi2_t3a = sum(((abs.(t3) .- data.t3amp) ./ data.t3amp_err) .^ 2)
    t3phi_model = angle.(t3) .* (180.0 / π)
    dphi = t3phi_model .- data.t3phi
    dphi = mod.(dphi .+ 180.0, 360.0) .- 180.0
    chi2_t3p = sum((dphi ./ data.t3phi_err) .^ 2)
    return chi2_v2, chi2_t3a, chi2_t3p
end

"""
    sky_forward(xi, p::SkyModelParams) -> Array{Float64,3}

Map latent vector to image cube (npix, npix, nfreq).
Both spatial and spectral fields use learned amplitude spectra.
"""
function sky_forward(xi::AbstractVector{<:Real}, p::SkyModelParams)
    xi_sp, xi_sc, sp, sc, x_alpha, x_logF0 = _unpack(xi, p)
    sp_slope, sp_fluct, sp_flex, sp_asp, sp_spec = sp
    sc_slope, sc_fluct, sc_flex, sc_asp, sc_spec = sc
    nf = length(p.freq)

    # Spatial correlated field
    amp_sp = amplitude_spectrum(sp_slope, sp_fluct, sp_flex, sp_asp, sp_spec,
                                p.spatial, p)
    kernel_sp = amp_sp[p.bin_index]
    # tau0_raw = real(ifft(fft(xi_sp) .* kernel_sp)) — transient transforms via scratch
    p.cs1 .= xi_sp
    mul!(p.cs2, p.P, p.cs1)
    p.cs2 .*= kernel_sp
    mul!(p.cs1, p.iP, p.cs2)
    tau0_raw = real.(p.cs1)

    # Zero-mean inside disk
    wmean_val = sum(p.D .* tau0_raw) / p.D_sum
    tau0 = p.D .* (tau0_raw .- wmean_val)

    # Spectral correlated field (shared spectrum, per-channel noise)
    # DC mode is zero — overall flux is handled by logF0 + alpha
    amp_sc = amplitude_spectrum(sc_slope, sc_fluct, sc_flex, sc_asp, sc_spec,
                                p.spectral, p)
    amp_sc[1] = 0.0
    kernel_sc = amp_sc[p.bin_index]
    delta_tau = Array{Float64}(undef, p.npix, p.npix, nf)
    for c in 1:nf
        p.cs1 .= view(xi_sc, :, :, c)
        mul!(p.cs2, p.P, p.cs1)
        p.cs2 .*= kernel_sc
        mul!(p.cs1, p.iP, p.cs2)
        @. delta_tau[:, :, c] = real(p.cs1)
    end

    # Build log-intensity cube
    Y = Array{Float64}(undef, p.npix, p.npix, nf)
    for c in 1:nf
        @. Y[:, :, c] = x_logF0 + tau0 + delta_tau[:, :, c] +
                         x_alpha * p.log_freq_ratio[c]
    end
    @. Y = clamp(Y, -30.0, 30.0)

    image = Array{Float64}(undef, p.npix, p.npix, nf)
    for c in 1:nf
        @. image[:, :, c] = p.D * exp(Y[:, :, c])
    end

    return image
end

"""
    sky_adjoint(g_image, xi, p::SkyModelParams) -> Vector{Float64}

Compute J'*g: adjoint of sky model Jacobian applied to image-space gradient.
"""
function sky_adjoint(g_image::AbstractArray{Float64,3},
                     xi::AbstractVector{Float64}, p::SkyModelParams)
    xi_sp, xi_sc, sp, sc, x_alpha, x_logF0 = _unpack(xi, p)
    sp_slope, sp_fluct, sp_flex, sp_asp, sp_spec = sp
    sc_slope, sc_fluct, sc_flex, sc_asp, sc_spec = sc
    nf = length(p.freq)
    n2 = p.npix^2

    # --- Recompute forward intermediates ---
    amp_sp = amplitude_spectrum(sp_slope, sp_fluct, sp_flex, sp_asp, sp_spec,
                                p.spatial, p)
    kernel_sp = amp_sp[p.bin_index]
    # F_xi = fft(xi_sp) — persists (g_kernel_sp below); ifft transient via scratch
    p.cs1 .= xi_sp
    mul!(p.cs2, p.P, p.cs1)
    F_xi = copy(p.cs2)
    p.cs2 .*= kernel_sp
    mul!(p.cs1, p.iP, p.cs2)
    tau0_raw = real.(p.cs1)
    wmean_val = sum(p.D .* tau0_raw) / p.D_sum
    tau0 = p.D .* (tau0_raw .- wmean_val)

    amp_sc = amplitude_spectrum(sc_slope, sc_fluct, sc_flex, sc_asp, sc_spec,
                                p.spectral, p)
    amp_sc[1] = 0.0  # DC mode handled by logF0 + alpha
    kernel_sc = amp_sc[p.bin_index]

    delta_tau = Array{Float64}(undef, p.npix, p.npix, nf)
    F_sc = Array{ComplexF64}(undef, p.npix, p.npix, nf)
    for c in 1:nf
        p.cs1 .= view(xi_sc, :, :, c)
        mul!(p.cs2, p.P, p.cs1)
        F_sc[:, :, c] .= p.cs2                 # persists (g_kernel_sc below)
        p.cs2 .*= kernel_sc
        mul!(p.cs1, p.iP, p.cs2)
        @. delta_tau[:, :, c] = real(p.cs1)
    end

    Y = Array{Float64}(undef, p.npix, p.npix, nf)
    for c in 1:nf
        @. Y[:, :, c] = x_logF0 + tau0 + delta_tau[:, :, c] +
                         x_alpha * p.log_freq_ratio[c]
    end
    Y_clamped = clamp.(Y, -30.0, 30.0)

    # --- Backward pass ---
    g_Y = Array{Float64}(undef, p.npix, p.npix, nf)
    for c in 1:nf
        @. g_Y[:, :, c] = p.D * exp(Y_clamped[:, :, c]) * g_image[:, :, c]
    end
    g_Y .*= (Y .> -30.0) .& (Y .< 30.0)

    g_logF0 = sum(g_Y)

    g_alpha = 0.0
    for c in 1:nf
        g_alpha += p.log_freq_ratio[c] * sum(view(g_Y, :, :, c))
    end

    # --- Spectral field adjoint ---
    g_xi = zeros(length(xi))
    g_kernel_sc_2d = zeros(p.npix, p.npix)
    for c in 1:nf
        p.cs1 .= view(g_Y, :, :, c)
        mul!(p.cs2, p.P, p.cs1)                  # G_sc
        # g_kernel_sc accumulates from all channels — uses G_sc + F_sc[c] (before
        # the buffer is reused for the ifft below)
        Fc = view(F_sc, :, :, c)
        @. g_kernel_sc_2d += real(conj(Fc) * p.cs2) / n2
        # g_xi_sc = real(ifft(kernel_sc .* G_sc))
        p.cs2 .*= kernel_sc
        mul!(p.cs1, p.iP, p.cs2)
        g_xi[n2 + (c-1)*n2 + 1 : n2 + c*n2] .= vec(real.(p.cs1))
    end

    # Scatter-add spectral kernel gradient to radial bins
    g_amp_sc = zeros(p.n_bins)
    for j in 1:p.npix, i in 1:p.npix
        g_amp_sc[p.bin_index[i,j]] += g_kernel_sc_2d[i,j]
    end
    g_amp_sc[1] = 0.0  # DC mode is fixed at zero

    g_sc_slope, g_sc_fluct, g_sc_flex, g_sc_asp, g_sc_spectrum, _ =
        amplitude_spectrum_adjoint(g_amp_sc, sc_slope, sc_fluct, sc_flex, sc_asp,
                                   sc_spec, p.spectral, p)

    # --- Spatial field adjoint ---
    g_tau0 = dropdims(sum(g_Y, dims=3), dims=3)

    # Zero-mean projection adjoint (self-adjoint)
    g_tau0_raw = p.D .* g_tau0 .- p.D .* (sum(p.D .* g_tau0) / p.D_sum)

    p.cs1 .= g_tau0_raw
    mul!(p.cs2, p.P, p.cs1)                       # G_raw
    # g_kernel_sp uses G_raw + F_xi — compute before reusing the buffer for ifft
    g_kernel_sp_2d = real.(conj.(F_xi) .* p.cs2) ./ n2
    # g_xi_sp = real(ifft(kernel_sp .* G_raw))
    p.cs2 .*= kernel_sp
    mul!(p.cs1, p.iP, p.cs2)
    g_xi[1:n2] .= vec(real.(p.cs1))

    g_amp_sp = zeros(p.n_bins)
    for j in 1:p.npix, i in 1:p.npix
        g_amp_sp[p.bin_index[i,j]] += g_kernel_sp_2d[i,j]
    end

    g_sp_slope, g_sp_fluct, g_sp_flex, g_sp_asp, g_sp_spectrum, _ =
        amplitude_spectrum_adjoint(g_amp_sp, sp_slope, sp_fluct, sp_flex, sp_asp,
                                   sp_spec, p.spatial, p)

    # --- Pack spectral parameter gradients ---
    offset = n2 + n2 * nf
    _pack_spectral_grad!(g_xi, offset, g_sp_slope, g_sp_fluct, g_sp_flex, g_sp_asp,
                         g_sp_spectrum, p.spatial, p.n_bins)
    offset += _spectral_latent_size(p.spatial, p.n_bins)
    _pack_spectral_grad!(g_xi, offset, g_sc_slope, g_sc_fluct, g_sc_flex, g_sc_asp,
                         g_sc_spectrum, p.spectral, p.n_bins)

    g_xi[end-1] = g_alpha
    g_xi[end]   = g_logF0

    return g_xi
end

function _pack_spectral_grad!(g_xi, offset, g_slope, g_fluct, g_flex, g_asp,
                              g_spectrum, cfg::CorrFieldConfig, n_bins::Int)
    g_xi[offset + 1] = g_slope
    g_xi[offset + 2] = g_fluct
    if cfg.use_iwp
        g_xi[offset + 3] = g_flex
        g_xi[offset + 4] = g_asp
        n_iwp_flat = 2 * (n_bins - 2)
        g_xi[offset + 5 : offset + 4 + n_iwp_flat] .= vec(g_spectrum)
    end
end
