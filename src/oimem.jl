"""
oimem.jl – Optical-interferometry operators for MaximENT, with OITOOLS integration.

This module provides the application-specific operators connecting MaximENT to
OIFITS data (power spectra + bispectra).  It can be used either with raw arrays
(low-level) or directly with an `OIdata` struct from OITOOLS (high-level).

Operators:
- `compute_data!`, `linearised_fwd!`, `linearised_adj!` — hot-path, allocation-free.
- `set_dataspace!` — populate data vector with full elliptic / classic bispectrum
  error models and triple-amplitude extrapolation from power spectra.
- ICF = identity by construction (hidden-space = visible-space).

Data vector layout  (length  npow + 2·nbis)
--------------------------------------------
  s.data[1 .. npow]              power-spectrum values   |V_i|²
  s.data[npow+2i-1 .. npow+2i]  bispectrum triple i:    Re / Im
                                 of  V_ab · V_bc · V_ca · φ_i
                                 where φ_i = exp(−i·ϕᵢ), ϕᵢ = measured closure phase

  Note: npow = nuv (all unique UV slots); flag_pow controls which are active.

High-level usage (OITOOLS)
--------------------------
    ctx, s, p = maxent_setup(data, nx, xyint, prior_image)
    result     = maxent_reconstruct!(ctx, s, p; maxiter=200, verbose=true)
    image      = reshape(result, nx, nx)

Low-level usage (raw arrays)
-----------------------------
    ctx = imaging_context(u, v, bs_pnt, bs_sign, closure_phases_deg,
                        flag_pow, flag_bs, nx, xyint)
    linfwd!, linadj!, fwd! = make_operators(ctx)

    s = MaximENTState(nx*nx, npow + 2*nbis)
    p = MaximENTParams(methd=[1,1,1], aim=1.0)
    s.model   .= prior_image
    s.h       .= prior_image
    s.data    .= packed_data_vector
    s.acc_vec .= packed_accuracy_vector
    reconstruct!(s, p, maxiter, linfwd!, linadj!, fwd!; verbose=true)
"""

using LinearAlgebra
using NFFT
using Printf



# ─────────────────────────────────────────────────────────────────────────────
# Physical constants
# ─────────────────────────────────────────────────────────────────────────────

"""Radians per milliarcsecond.  MAS = π / (180 · 3 600 000)."""
const MAS = (π / 180.0) / 3_600_000.0

# ─────────────────────────────────────────────────────────────────────────────
# Utility helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    auto_pixsize(data; oversampling=3.0) → Float64

Return a pixel size in mas giving `oversampling` pixels per fringe at the
longest baseline.
"""
function auto_pixsize(data; oversampling=3.0)
    uvmax = maximum(sqrt.(data.uv[1,:].^2 .+ data.uv[2,:].^2))
    return 1.0 / (2.0 * MAS * oversampling * uvmax)
end

"""
    gaussian_prior(nx, pixsize_mas; fwhm_mas, flux=1.0) → Vector{Float64}

Centred Gaussian prior matching bsmem.c modeltype=3.
`fwhm_mas` defaults to 1/5 of the image FoV.
Pixels are clipped at 1e-8 to avoid log(0) in the entropy.
"""
function gaussian_prior(nx::Int, pixsize_mas::Float64;
                        fwhm_mas::Float64 = nx * pixsize_mas / 5.0,
                        flux::Float64 = 1.0)
    sigma_pix = fwhm_mas / (2.0 * sqrt(2.0 * log(2.0))) / pixsize_mas
    cx = nx / 2.0
    model = [exp(-((i - cx)^2 + (j - cx)^2) / (2 * sigma_pix^2))
             for j in 0:nx-1, i in 0:nx-1] |> vec
    model .= max.(model, 1e-8)
    model .*= flux / sum(model)
    return model
end

# ─────────────────────────────────────────────────────────────────────────────
# ImagingContext
# ─────────────────────────────────────────────────────────────────────────────

"""
    ImagingContext

Per-reconstruction workspace shared by compute_data!, linearised_fwd!, and linearised_adj!.
Construct via `imaging_context(...)`.

`T` is the floating-point precision of the Fourier operator and its complex scratch —
`Float32` roughly halves the plan's memory and speeds up the transform. The MaximENT
solver itself remains `Float64` regardless (see `maxent_reconstruct!`); the widening at
`compute_data!`/`linearised_fwd!` and the narrowing at `linearised_adj!` convert at the
context boundary, which costs O(nx²) against an O(nuv·m² + nx² log nx) transform.

Fields
------
- `npow`, `nbis`, `nuv`, `nx` : problem dimensions
- `xyint`    : pixel size [mas] — always Float64; the uv→k scaling chains three
  multiplies off u ~ 1e8 λ and must not be computed in reduced precision
- `plan`     : 2-D NFFTPlan mapping nx×nx image ↔ nuv complex visibilities
- `visi`     : current complex visibilities (set by compute_data!, read by linearised_fwd!/linearised_adj!)
- `bs_pnt`   : (3 × nbis) 1-based UV-table indices; rows = [ab; bc; ca]
- `bs_sign`  : (3 × nbis) Int8 conjugation flags ±1 per bispectrum leg
- `phasor`   : exp(−i·ϕ_data_rad) per triple
- `flag_pow` : Bool vector length npow  (false = flagged / no measurement)
- `flag_bs`  : Bool vector length nbis  (false = flagged / excluded)
- `_f_hat`   : nx×nx Complex{T} scratch for NFFT input/output
- `_dvisi`   : npow Complex{T} scratch for adjoint accumulation
"""
mutable struct ImagingContext{T<:AbstractFloat}
    npow     :: Int
    nbis     :: Int
    nuv      :: Int
    nx       :: Int
    xyint    :: Float64
    plan     :: NFFTPlan{T,2,1}   # concrete: a bare `NFFTPlan` makes every mul! a dynamic dispatch
    visi     :: Vector{Complex{T}}
    bs_pnt   :: Matrix{Int}    # 3 × nbis
    bs_sign  :: Matrix{Int8}   # 3 × nbis
    phasor   :: Vector{Complex{T}}
    flag_pow :: Vector{Bool}
    flag_bs  :: Vector{Bool}
    _f_hat   :: Matrix{Complex{T}}
    _dvisi   :: Vector{Complex{T}}
end

# ─────────────────────────────────────────────────────────────────────────────
# Low-level constructor (raw arrays)
# ─────────────────────────────────────────────────────────────────────────────

"""
    imaging_context(u, v, bs_pnt, bs_sign, closure_phases_deg,
                  flag_pow, flag_bs, nx, xyint) → ImagingContext

Construct a `ImagingContext` from raw interferometry arrays and precompute the
NFFT plan.

Arguments
---------
- `u`, `v`             : UV coordinates in **wavelengths**, length nuv.
- `bs_pnt`             : (nbis × 3) integer matrix; columns [ab_idx, bc_idx, ca_idx],
                         **1-based** indices into the UV table.
- `bs_sign`            : (nbis × 3) Int8 matrix; ±1 conjugation flags.
                         −1 means that baseline is from the conjugated half-plane.
- `closure_phases_deg` : measured closure phases in degrees, length nbis.
- `flag_pow`           : Bool vector length npow; `true` = use this power-spectrum slot.
- `flag_bs`            : Bool vector length nbis; `true` = use this closure triple.
- `nx`                 : image side in pixels (image is nx × nx).
- `xyint`              : pixel size in milliarcseconds.

NFFT convention
---------------
Normalised frequencies follow the C bsmem convention:

    k[1, i] = v[i] · xyint · MAS      (first NFFT dimension  ↔ v, N–S baseline)
    k[2, i] = u[i] · xyint · MAS      (second NFFT dimension ↔ u, E–W baseline)

Kaiser–Bessel window: half-width m = 6, oversampling σ = 2 (matches
`nfft_init_guru` with PRE_FULL_PSI and nn = 2·NN in the original C code).
"""
function imaging_context(u         :: AbstractVector{<:Real},
                       v         :: AbstractVector{<:Real},
                       bs_pnt    :: AbstractMatrix{<:Integer},  # nbis × 3
                       bs_sign   :: AbstractMatrix{<:Integer},  # nbis × 3
                       closure_phases_deg :: AbstractVector{<:Real},
                       flag_pow  :: AbstractVector{Bool},
                       flag_bs   :: AbstractVector{Bool},
                       nx        :: Int,
                       xyint     :: Float64;
                       T         :: Type{<:AbstractFloat} = Float64,
                       fftflags  = FFT_FLAGS)
    nuv  = length(u)
    nbis = size(bs_pnt, 1)
    @assert length(v) == nuv
    @assert size(bs_pnt)  == (nbis, 3)
    @assert size(bs_sign) == (nbis, 3)
    @assert length(closure_phases_deg) == nbis
    @assert length(flag_bs) == nbis
    @assert length(flag_pow) == nuv

    scale = xyint * MAS

    # Prepend a (u=0, v=0) point as UV slot 1 (matching the C convention in
    # oifits.c: uv[0]=(0,0), pow[0]=1, powerr[0]=fluxerr, then real data from
    # index 1 onward). V(0,0) = Σ image = total flux.
    #
    # eltype(k) alone fixes the plan's precision. Build the node positions in Float64 and
    # narrow only on store: u·xyint·MAS chains three multiplies off u ~ 1e8 λ, and doing
    # that chain in Float32 would cost ~7 significant digits on the node positions
    # themselves. Same idiom as setup_nfft (oichi2.jl).
    # Row 1 is u and row 2 is v, matching `data.uv` and the plans `setup_nfft` builds. NFFT node
    # row 1 indexes the FIRST array dimension, so this is what fixes the image plane: with the
    # rows the other way round the reconstruction came back transposed against every other
    # engine and against `imdisp`, which reads dimension 1 as α. It also decides how `x_start`
    # and a prior image are READ, so those were being consumed transposed too.
    #
    # Measured on 2004-simulated: these nodes now match `setup_nfft`'s to 5e-9 where they
    # differed by 0.12, and scoring the result through the shared criterion in all eight
    # orientations inverts as it should -- the transpose used to beat the image as returned
    # (719 against 817) and now loses to it (856 against 772).
    #
    # This does NOT bring BSMEM's own χ² into line with `image_to_chi2` (500 against 772 on
    # that run), and it was never going to: BSMEM scores the bispectrum through its elliptic
    # noise approximation rather than the shared criterion's residuals, so the two numbers
    # measure different things on the same image. Orientation was one bug, not both.
    k64 = Matrix{Float64}(undef, 2, nuv + 1)
    k64[:, 1]       .= 0.0          # zero-baseline flux point (slot 1)
    k64[1, 2:end]   .= Float64.(u) .* scale   # u ↔ E–W baseline
    k64[2, 2:end]   .= Float64.(v) .* scale   # v ↔ N–S baseline
    k = T === Float64 ? k64 : T.(k64)

    npow = nuv + 1
    flag_pow_ext = vcat(true, Vector{Bool}(flag_pow))  # flux slot always active

    # FFTW planning strategy: see `FFT_FLAGS` in oichi2.jl for the measurements motivating
    # MEASURE over NFFT's unset-flags default (which makes FFTW use ESTIMATE).
    plan = NFFTPlan(k, (nx, nx); m=6, σ=2.0,
                    precompute=NFFT.POLYNOMIAL,
                    blocking=true, sortNodes=false, fftflags)

    # Transpose to 3×nbis for fast column access in the hot-path loops
    bs_pnt_T  = Matrix{Int}( permutedims(bs_pnt)  )
    bs_sign_T = Matrix{Int8}(permutedims(bs_sign) )

    # Phase conversion in Float64, narrowed on store (same reasoning as k).
    phasor = Complex{T}.(exp.(-im .* (Float64.(closure_phases_deg) .* (π / 180.0))))

    ImagingContext{T}(
        npow, nbis, nuv, nx, Float64(xyint),
        plan,
        zeros(Complex{T}, npow),
        bs_pnt_T, bs_sign_T,
        phasor,
        flag_pow_ext, Vector{Bool}(flag_bs),
        zeros(Complex{T}, nx, nx),
        zeros(Complex{T}, npow)
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# High-level constructor from OIdata (OITOOLS)
# ─────────────────────────────────────────────────────────────────────────────

"""
    imaging_context(data::OIdata, nx, xyint) → ImagingContext

Construct a `ImagingContext` directly from an OITOOLS `OIdata` struct.

The UV table (length nuv) is taken from `data.uv_lam`.  Power-spectrum slots
1..nuv are flagged active where a V² measurement exists and is not marked bad
in `data.v2_flag`.  Bispectrum triples use `data.indx_t3_1/2/3` for UV
indices; in OITOOLS convention **negative** index values indicate conjugated
baselines (sign = −1).

OITOOLS is assumed to have already applied Nyquist filtering and bad-baseline
masking before calling this function (via its own `select_data` / `filter_data`
pipeline).

Data layout: npow = nuv → MaximENT4 data vector has length nuv + 2·nt3phi.
"""
function imaging_context(data,            # OIdata from OITOOLS
                       nx    :: Int,
                       xyint :: Float64;
                       T     :: Type{<:AbstractFloat} = Float64,
                       fftflags = FFT_FLAGS)

    u = data.uv[1, :]   # E–W, wavelengths
    v = data.uv[2, :]   # N–S, wavelengths
    nuv = length(u)

    # Power-spectrum flags for real UV points (slots 2..nuv+1 after prepending flux)
    # OITOOLS indx_v2 is 1-based into data.uv; shift by +1 for our extended table.
    flag_pow = zeros(Bool, nuv)   # for real UV slots only; flux slot added in constructor
    for j in eachindex(data.v2_flag)
        !data.v2_flag[j] && (flag_pow[data.indx_v2[j]] = true)
    end

    # Bispectrum tables: shift indices by +1; negative = conjugated (OITOOLS convention)
    nbis    = length(data.t3phi)
    bs_pnt  = Matrix{Int}( undef, nbis, 3)
    bs_sign = Matrix{Int8}(undef, nbis, 3)
    flag_bs = trues(nbis)

    for i in 1:nbis
        for (col, idx) in enumerate((data.indx_t3_1[i], data.indx_t3_2[i], data.indx_t3_3[i]))
            bs_pnt[i, col]  = abs(idx) + 1   # +1 for prepended flux slot
            bs_sign[i, col] = Int8(sign(idx))
        end
        data.t3_flag[i] && (flag_bs[i] = false)
    end

    imaging_context(u, v, bs_pnt, bs_sign,
                  data.t3phi,      # closure phases in degrees
                  flag_pow, flag_bs, nx, xyint; T, fftflags)
end

# ─────────────────────────────────────────────────────────────────────────────
# set_dataspace!  –  port of set_maximent_dataspace() in bsmem.c
# ─────────────────────────────────────────────────────────────────────────────

"""
    set_dataspace!(s, ctx, data;
                   biserrtype=:full,
                   force_extrapolate=false,
                   flux_err=0.01,
                   verbose=true)

Fill `s.data` and `s.acc_vec` from OIFITS observables.  The total data vector
length must equal `ctx.npow + 2·ctx.nbis`.

Power spectrum (positions 1..npow)
-----------------------------------
    s.data[indx_v2[j]]    = v2[j]
    s.acc_vec[indx_v2[j]] = 1 / v2_err[j]

Bad measurements (err ≤ 0 or NaN) → both slots → 0 (excluded from χ²).

Bispectrum (positions npow+2i-1 and npow+2i for triple i)
----------------------------------------------------------
Two error models (argument `biserrtype`):

  `:full`    — Full elliptic approximation (Meimon 2005 eq. 14, 2009 appendix E):
                   err_rad = √{½(|B|²+σ_A²)(1+e^{−2σ_φ²}) − |B|²e^{−σ_φ²}}
                   err_tan = √{½(|B|²+σ_A²)(1−e^{−2σ_φ²})}
               Data amplitude is bias-corrected:
                   s.data[Re] = |B|·exp(−½σ_φ²),   s.data[Im] = 0.

  `:classic` — First-order approximation:
                   err_rad = t3amp_err
                   err_tan = |B|·t3phi_err·π/180

    s.acc_vec[Re] = 1/err_rad,  s.acc_vec[Im] = 1/err_tan.
    Bad measurements → both → 0.

Triple-amplitude extrapolation
-------------------------------
When `t3amp_err[i] ≤ 0` or `Inf` (missing), or always when
`force_extrapolate=true`, the triple amplitude is estimated from the three
matching power-spectrum measurements using an unbiased estimator
(Meimon 2009 appendix).  If any of the three power-spectrum slots is also
missing, the triple is zeroed (flagged).
"""
function set_dataspace!(s   :: MaximENTState,
                        ctx :: ImagingContext,
                        data;           # OIdata from OITOOLS
                        biserrtype        :: Symbol = :full,
                        force_extrapolate :: Bool   = false,
                        flux_err          :: Float64 = 0.01,
                        verbose           :: Bool   = true)

    npow = ctx.npow   # = nuv + 1  (last slot = flux constraint)
    nbis = ctx.nbis
    fill!(s.data,    0.0)
    fill!(s.acc_vec, 0.0)

    # ── Power spectrum ────────────────────────────────────────────────────────
    # Slot 1 = flux constraint (filled below). Real V² slots are indx_v2[j]+1.
    nv2 = length(data.v2)
    for j in 1:nv2
        i = data.indx_v2[j] + 1   # +1: slot 1 reserved for flux constraint
        e = data.v2_err[j]
        if e > 0.0 && isfinite(e) && !data.v2_flag[j]
            s.data[i]    = data.v2[j]
            s.acc_vec[i] = 1.0 / e
        else
            e <= 0.0 && @printf("Warning: V² err[%d] ≤ 0 — slot zeroed\n", j)
        end
    end

    # Scratch arrays for triple-amplitude extrapolation lookup
    pow_val = copy(s.data[1:npow])
    pow_err = zeros(Float64, npow)
    for j in 1:nv2
        i = data.indx_v2[j] + 1   # +1 for flux slot
        pow_err[i] = (data.v2_err[j] > 0.0 && !data.v2_flag[j]) ? data.v2_err[j] : 0.0
    end

    # ── Bispectrum ────────────────────────────────────────────────────────────
    # Under `verbose`, like everything else BSMEM prints. Unconditional, this made
    # `verbose = false` a half-promise: a caller asking for silence still got a line.
    verbose && biserrtype == :full &&
        println("Bispectrum noise:\tFull elliptic approximation")
    verbose && biserrtype == :classic &&
        println("Bispectrum noise:\tClassic elliptic approximation")
    warn_extrap = false

    for i in 1:nbis
        ab = ctx.bs_pnt[1, i]
        bc = ctx.bs_pnt[2, i]
        ca = ctx.bs_pnt[3, i]

        bisamp_i    = data.t3amp[i]
        bisamperr_i = data.t3amp_err[i]
        bisphserr_i = data.t3phi_err[i]

        need_extrap = force_extrapolate ||
                      !(bisamperr_i > 0.0 && isfinite(bisamperr_i))

        if need_extrap
            p1 = pow_val[ab]; e1 = pow_err[ab]
            p2 = pow_val[bc]; e2 = pow_err[bc]
            p3 = pow_val[ca]; e3 = pow_err[ca]

            if e1 > 0.0 && e2 > 0.0 && e3 > 0.0
                # Unbiased triple-amplitude from V² (Meimon 2009)
                sq1  = (p1 + sqrt(p1^2 + 2e1^2)) / 2;  sqe1 = 1.0 / (1/sq1 + 2*(3sq1-p1)/e1^2)
                sq2  = (p2 + sqrt(p2^2 + 2e2^2)) / 2;  sqe2 = 1.0 / (1/sq2 + 2*(3sq2-p2)/e2^2)
                sq3  = (p3 + sqrt(p3^2 + 2e3^2)) / 2;  sqe3 = 1.0 / (1/sq3 + 2*(3sq3-p3)/e3^2)
                bisamp_i    = sqrt(sq1 * sq2 * sq3)
                bisamperr_i = abs(bisamp_i) * sqrt(sqe1/sq1 + sqe2/sq2 + sqe3/sq3)
                warn_extrap || (warn_extrap = true;
                    println("WARNING: ≥1 triple amplitude extrapolated from V²"))
            else
                println("WARNING: triple amplitude extrapolation failed for triple $i — flagged")
                # leave acc_vec at 0 (already initialised), skip
                continue
            end
        end

        # ── Error radii ───────────────────────────────────────────────────────
        σ_A  = bisamperr_i
        σ_φ² = (bisphserr_i * π / 180.0)^2

        if biserrtype === :full
            B²   = bisamp_i^2
            eA²  = σ_A^2
            eσ2  = exp(-2σ_φ²)
            eσ   = exp(-σ_φ²)
            err_rad = sqrt(max(0.0, 0.5*(B² + eA²)*(1 + eσ2) - B²*eσ))
            err_tan = sqrt(max(0.0, 0.5*(B² + eA²)*(1 - eσ2)))
            s.data[npow + 2i - 1] = bisamp_i * exp(-0.5σ_φ²)
        else  # :classic
            err_rad = σ_A
            err_tan = abs(bisamp_i * bisphserr_i * π / 180.0)
            s.data[npow + 2i - 1] = bisamp_i
        end
        s.data[npow + 2i] = 0.0   # Im part of phasor-rotated bispectrum

        # ── Final bad-data check ──────────────────────────────────────────────
        bad = !(bisamperr_i > 0.0 && isfinite(bisamperr_i)) ||
              !(bisphserr_i > 0.0 && isfinite(bisphserr_i)) ||
              data.t3_flag[i]
        if bad
            s.data[npow + 2i - 1]    = 0.0
            s.data[npow + 2i]        = 0.0
            s.acc_vec[npow + 2i - 1] = 0.0
            s.acc_vec[npow + 2i]     = 0.0
        else
            s.acc_vec[npow + 2i - 1] = err_rad > 0.0 ? 1.0 / err_rad : 0.0
            s.acc_vec[npow + 2i]     = err_tan > 0.0 ? 1.0 / err_tan : 0.0
        end
    end

    # ── Flux constraint: slot 1 = zero-baseline, V(0,0) = total flux = 1.0 ──
    s.data[1]    = 1.0
    s.acc_vec[1] = flux_err > 0.0 ? 1.0 / (2.0 * flux_err) : 0.0

    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# maxent_setup  –  one-stop convenience constructor
# ─────────────────────────────────────────────────────────────────────────────

"""
    maxent_setup(data, nx, xyint, prior_image=nothing;
                methd=[1,1,1], nrand=1, iseed=0,
                aim=1.0, rate=1.0, utol=0.1, alpha=1.0,
                biserrtype=:full,
                force_extrapolate=false,
                mackay_alpha=false, ritz_alpha=false)
    → (ctx::ImagingContext, s::MaximENTState, p::MaximENTParams)

Build everything needed for `maxent_reconstruct!` in a single call.

`prior_image` may be a length-nx² vector or an nx×nx matrix; it is normalised
to unit flux.  Pass `nothing` (default) for a flat uniform prior.

MaximENT4 parameters map directly to `MaximENTParams` keyword arguments:
  `methd`  = [regularisation, entropy, noise_model, nonlinearity]
  `aim`    = target α (1.0 for classic Bayesian, or fixed α value)
  `utol`   = convergence tolerance on χ²
"""
function maxent_setup(data,
                     nx          :: Int,
                     xyint       :: Float64,
                     prior_image = nothing;
                     methd       :: Vector{Int} = [1, 1, 1],
                     nrand       :: Int         = 1,
                     iseed       :: Int         = 0,
                     aim         :: Float64     = 1.0,
                     rate        :: Float64     = 1.0,
                     utol        :: Float64     = 0.1,
                     alpha       :: Float64     = 1.0,
                     biserrtype        :: Symbol = :full,
                     force_extrapolate :: Bool   = false,
                     flux_err          :: Float64 = 0.01,
                     mackay_alpha      :: Bool    = false,
                     ritz_alpha        :: Bool    = false,
                     verbose     :: Bool        = true,
                     T           :: Type{<:AbstractFloat} = Float64,
                     fftflags    = FFT_FLAGS)

    ctx  = imaging_context(data, nx, xyint; T, fftflags)
    npix = nx * nx
    ndat = ctx.npow + 2 * ctx.nbis

    s = MaximENTState{T}(npix, ndat)
    p = MaximENTParams(; methd, nrand, iseed, aim, rate, utol, alpha, mackay_alpha, ritz_alpha)

    model = if isnothing(prior_image)
        fill(1.0 / npix, npix)
    else
        vf = vec(Float64.(prior_image))
        vf ./ max(sum(vf), eps())
    end
    s.model .= model
    s.h     .= model

    set_dataspace!(s, ctx, data;
                   biserrtype, force_extrapolate, flux_err, verbose)

    return ctx, s, p
end

# ─────────────────────────────────────────────────────────────────────────────
# maxent_reconstruct!
# ─────────────────────────────────────────────────────────────────────────────

"""
    maxent_reconstruct!(ctx, s, p; maxiter=200, verbose=true) → Matrix{Float64}

Run the MaximENT iterative reconstruction loop and return the result as an
`nx × nx` `Matrix{Float64}`, normalised to unit flux.

`ctx`, `s`, and `p` are all modified in-place.  The same result is accessible
afterwards as `reshape(s.h, ctx.nx, ctx.nx)`.  Call `maxent_setup`
to build all three from an `OIdata` struct.

# Arguments

- `ctx::ImagingContext` — operator context: NFFT plan, bispectrum tables,
  conjugation signs, phasors, and visibility scratch arrays.
- `s::MaximENTState` — algorithm state: current image `s.h`, prior `s.model`,
  packed data vector `s.data`, accuracy vector `s.acc_vec`, and internal
  workspace arrays.
- `p::MaximENTParams` — hyper-parameters: entropy mode, noise model, stopping
  criterion (`methd`), regularisation strength `alpha`, convergence tolerance
  `utol`, and α-update strategy.

# Keyword arguments

- `maxiter` — maximum number of outer MaximENT iterations. Default: `200`.
- `verbose` — print per-iteration χ², entropy, and α diagnostics to stdout.
  Default: `true`.
- `history` — if a `Vector` is passed, the per-iteration diagnostic trace is appended to
  it (entropy `S`, `chisq`, `alpha`, `omega`, `ntrans`, `istat`, evidence bounds, …).
  A far more precise convergence fingerprint than the `verbose` printout, which is
  rounded to 4-5 significant digits. Default: `nothing` (no trace, no cost).

# Example

```julia
prior = gaussian_prior(nx, pixsize; fwhm_mas = nx * pixsize / 5)
ctx, s, p = maxent_setup(data, nx, pixsize, prior; nrand=10)
image = maxent_reconstruct!(ctx, s, p; maxiter=200, verbose=true)
imdisp(image; pixsize)
```

See also: `maxent_setup`, [`reconstruct_bsmem`](@ref).
"""
function maxent_reconstruct!(ctx :: ImagingContext,
                             s   :: MaximENTState,
                             p   :: MaximENTParams;
                             maxiter :: Int  = 200,
                             verbose :: Bool = true,
                             history :: Union{Nothing,Vector} = nothing)
    linfwd! = (δim, δdat) -> linearised_fwd!(ctx, δim, δdat)
    linadj! = (δdat, δim) -> linearised_adj!(ctx, δdat, δim)
    fwd!    = (im,  dat)  -> compute_data!(ctx, im, dat)

    # Warm-up: populate ctx.visi from the starting model
    # before the first MaxEnt iteration
    compute_data!(ctx, s.h, s.d_w2)

    h = reconstruct!(s, p, maxiter, linfwd!, linadj!, fwd!; verbose)
    history === nothing || append!(history, h)

    return reshape(s.h, ctx.nx, ctx.nx)
end

# ─────────────────────────────────────────────────────────────────────────────
# reconstruct_bsmem  –  high-level interface matching reconstruct()
# ─────────────────────────────────────────────────────────────────────────────

"""
    reconstruct_bsmem(x_start, data, ft; kwargs...) → Matrix

Maximum Entropy image reconstruction using power spectra (V²) and closure
phases (T3φ), with the same calling convention as `reconstruct`.

The MaximENT solver itself always runs in `Float64` — it builds its own NFFT plan,
and its entropy, Ritz-value and trust-region steps are precision-sensitive — but the
image is returned at the precision of `ft`, so it composes with the rest of the
pipeline (`Float32` by default).

This is the recommended high-level entry point.  Pixel size and image size
are inferred from the OITOOLS NFFT plan `ft`; you do not need to call
`maxent_setup` or `maxent_reconstruct!` directly.

# Arguments

- `x_start` — starting image: an nx×nx matrix or length-nx² vector.  Used as
  the initial iterate and, unless overridden via `regularizers`, as the entropy
  prior.
- `data::OIdata` — interferometric data loaded by `readoifits`.
- `ft` — NFFT plan produced by `setup_nfft(data, nx, pixsize)`.

# Keyword arguments

## Reconstruction control
- `method` — MaximENT stopping criterion.  Either a 4-element `Int` vector
  passed directly as `methd`, or one of these shorthands:
  - `1` → `[1,1,1,2]` classic known noise
  - `4` → `[4,1,1,2]` χ² = N  (default)
- `maxiter` — maximum outer iterations. Default: `200`.
- `verbose` — print iteration diagnostics. Default: `true`.

## Prior / entropy
- `regularizers` — list of regularization entries.  The first entry of the
  form `["mem", model]` supplies the entropy prior image; all other entries
  are ignored by MaxEnt.  If no `"mem"` entry is found, `x_start` is used as
  the prior.

## Data-space options
- `flux_err` — relative flux uncertainty for the zero-baseline constraint.
  Default: `1e-5`.
- `biserrtype` — bispectrum noise model: `:full` (Meimon 2005 elliptic,
  default) or `:classic` (first-order approximation).
- `force_extrapolate` — always estimate triple amplitudes from V² rather than
  using measured values. Default: `false`.

## MaximENT hyper-parameters
- `nrand` — number of random probe vectors for evidence estimation. Default: `10`.
- `aim`, `rate`, `utol`, `alpha` — MaximENT convergence parameters.
- `mackay_alpha` — use MacKay fixed-point α update. Default: `false`.
- `ritz_alpha` — use Ritz-value bisection α update. Default: `false`.

# Example

```julia
data   = readoifits("data/2004-data1.oifits")[1,1]
nx     = 128
pixsize = 0.05   # mas/pixel
ft      = setup_nfft(data, nx, pixsize)
prior   = gaussian_prior(nx, pixsize; fwhm_mas = nx * pixsize / 5)
x = reconstruct_bsmem(prior, data, ft;
                      regularizers = [["mem", prior]],
                      method       = [1, 1, 1, 2],
                      maxiter      = 100,
                      flux_err     = 1e-5)
imdisp(x; pixsize)
```

See also: `maxent_setup`, `maxent_reconstruct!`.
"""
function reconstruct_bsmem(x_start, data::OIdata, ft;
                            regularizers      = [],
                            method            = [4, 1, 1, 2],
                            maxiter  :: Int   = 200,
                            verbose  :: Bool  = true,
                            flux_err :: Float64 = 1e-5,
                            biserrtype        :: Symbol = :full,
                            force_extrapolate :: Bool   = false,
                            nrand  :: Int     = 10,
                            iseed  :: Int     = 0,
                            aim    :: Float64 = 1.0,
                            rate   :: Float64 = 1.0,
                            utol   :: Float64 = 0.1,
                            alpha        :: Float64 = 1.0,
                            mackay_alpha :: Bool    = false,
                            ritz_alpha   :: Bool    = false,
                            fftflags               = FFT_FLAGS,
                            history :: Union{Nothing,Vector} = nothing)

    # ── Extract nx and pixsize from the OITOOLS NFFT plan ─────────────────────
    # setup_nfft stores ft[1].k[:, j] = pixsize_rad * [u_j; v_j],
    # so norm(ft[1].k[:, j]) = pixsize_rad * norm(data.uv[:, j]).
    nx  = ft[1].N[1]
    j   = argmax(vec(sqrt.(sum(data.uv .^ 2, dims=1))))
    pixsize = (norm(ft[1].k[:, j]) / norm(data.uv[:, j])) / MAS

    # ── Parse method (scalar shorthand or full vector) ─────────────────────────
    methd_vec = method isa Integer ? [Int(method), 1, 1, 2] : collect(Int, method)

    # ── Extract prior from regularizers ───────────────────────────────────────
    # Accepts ["mem", prior_image]; falls back to x_start if absent.
    prior_image = vec(Float64.(x_start))
    for reg in regularizers
        if !isempty(reg) && reg[1] == "mem" && length(reg) >= 2 && !isnothing(reg[2])
            prior_image = vec(Float64.(reg[2]))
            break
        end
    end

    # ── Build context and run ─────────────────────────────────────────────────
    # The Fourier operator follows the precision of `ft` (Float32 by default since
    # readoifits does); the MaximENT solver itself stays Float64 either way.
    ctx, s, p = maxent_setup(data, nx, pixsize, prior_image;
                             methd=methd_vec, nrand, iseed, aim, rate, utol, alpha,
                             biserrtype, force_extrapolate, flux_err, mackay_alpha, ritz_alpha,
                             verbose, T = ft_eltype(ft), fftflags)

    # The solver is Float64 throughout; hand the image back at the caller's precision.
    return to_ft_precision(maxent_reconstruct!(ctx, s, p; maxiter, verbose, history), ft)
end

# ─────────────────────────────────────────────────────────────────────────────
# Polychromatic MaximENT
# ─────────────────────────────────────────────────────────────────────────────

"""
    PolyImagingContext

Multi-wavelength wrapper: one `ImagingContext` per spectral channel, with
bookkeeping for the concatenated data vector used by MaximENT.

Fields
------
- `nwav`    : number of wavelength channels
- `nx`      : image side length (pixels), same for all channels
- `xyint`   : pixel size [mas], same for all channels
- `ctxs`    : `Vector{ImagingContext{T}}`, one per channel
- `ndat_ch` : number of data slots per channel (`npow + 2*nbis`)
- `offset`  : cumulative data offset; `offset[w]` = sum of ndat_ch[1:w-1]

`T` is the Fourier-operator precision, as for [`ImagingContext`](@ref). This is where
`Float32` pays off most: the plan and complex scratch are replicated per channel.
"""
struct PolyImagingContext{T<:AbstractFloat}
    nwav    :: Int
    nx      :: Int
    xyint   :: Float64
    ctxs    :: Vector{ImagingContext{T}}
    ndat_ch :: Vector{Int}
    offset  :: Vector{Int}
end

"""
    PolyImagingContext(data_channels, nx, xyint; T=Float64) → PolyImagingContext

Construct from a vector of `OIdata` (one per wavelength channel).
"""
function PolyImagingContext(data_channels::AbstractVector, nx::Int, xyint::Float64;
                            T::Type{<:AbstractFloat} = Float64, fftflags = FFT_FLAGS)
    nwav = length(data_channels)
    ctxs = ImagingContext{T}[imaging_context(data_channels[w], nx, xyint; T, fftflags) for w in 1:nwav]
    ndat_ch = [ctx.npow + 2 * ctx.nbis for ctx in ctxs]
    offset  = [sum(ndat_ch[1:w-1]) for w in 1:nwav]
    return PolyImagingContext{T}(nwav, nx, xyint, ctxs, ndat_ch, offset)
end

# -- Polychromatic operators -------------------------------------------------

function compute_data_poly!(pctx::PolyImagingContext,
                            image_cube::AbstractVector{<:AbstractFloat},
                            data_vec::AbstractVector{<:AbstractFloat})
    npix = pctx.nx * pctx.nx
    for w in 1:pctx.nwav
        img_w = @view image_cube[(w-1)*npix+1 : w*npix]
        dat_w = @view data_vec[pctx.offset[w]+1 : pctx.offset[w]+pctx.ndat_ch[w]]
        compute_data!(pctx.ctxs[w], img_w, dat_w)
    end
    return nothing
end

function linearised_fwd_poly!(pctx::PolyImagingContext,
                              δimage::AbstractVector{<:AbstractFloat},
                              δdata::AbstractVector{<:AbstractFloat})
    npix = pctx.nx * pctx.nx
    for w in 1:pctx.nwav
        δimg_w = @view δimage[(w-1)*npix+1 : w*npix]
        δdat_w = @view δdata[pctx.offset[w]+1 : pctx.offset[w]+pctx.ndat_ch[w]]
        linearised_fwd!(pctx.ctxs[w], δimg_w, δdat_w)
    end
    return nothing
end

function linearised_adj_poly!(pctx::PolyImagingContext,
                              δdata::AbstractVector{<:AbstractFloat},
                              δimage::AbstractVector{<:AbstractFloat})
    npix = pctx.nx * pctx.nx
    for w in 1:pctx.nwav
        δdat_w = @view δdata[pctx.offset[w]+1 : pctx.offset[w]+pctx.ndat_ch[w]]
        δimg_w = @view δimage[(w-1)*npix+1 : w*npix]
        linearised_adj!(pctx.ctxs[w], δdat_w, δimg_w)
    end
    return nothing
end

# -- Polychromatic data packing ----------------------------------------------

"""
    set_dataspace_poly!(s, pctx, data_channels; kwargs...)

Fill the concatenated data/accuracy vectors for all wavelength channels.
Delegates to the monochromatic `set_dataspace!` for each channel.
"""
function set_dataspace_poly!(s::MaximENTState,
                             pctx::PolyImagingContext,
                             data_channels::AbstractVector;
                             biserrtype::Symbol = :full,
                             force_extrapolate::Bool = false,
                             flux_err::Float64 = 0.01)
    npix = pctx.nx * pctx.nx
    for w in 1:pctx.nwav
        ndat_w = pctx.ndat_ch[w]
        s_tmp = MaximENTState(npix, ndat_w)
        set_dataspace!(s_tmp, pctx.ctxs[w], data_channels[w];
                       biserrtype=biserrtype,
                       force_extrapolate=force_extrapolate,
                       flux_err=flux_err)
        off = pctx.offset[w]
        s.data[off+1 : off+ndat_w]    .= s_tmp.data
        s.acc_vec[off+1 : off+ndat_w] .= s_tmp.acc_vec
    end
    return nothing
end

# -- Polychromatic setup -----------------------------------------------------

"""
    maxent_setup_poly(data_channels, nx, xyint, prior_cube=nothing; kwargs...)
        → (pctx::PolyImagingContext, s::MaximENTState, p::MaximENTParams)

Build everything needed for polychromatic MaximENT reconstruction.

`data_channels` is a `Vector{<:OIdata}` with one entry per wavelength channel.

`prior_cube` may be:
- `nothing`           → flat uniform prior (replicated per channel)
- an `nx×nx` matrix   → replicated identically across all channels
- an `nx×nx×nwav` array → used directly (normalised per-channel to unit flux)
"""
function maxent_setup_poly(data_channels::AbstractVector,
                           nx::Int, xyint::Float64,
                           prior_cube = nothing;
                           methd       :: Vector{Int} = [1, 1, 1],
                           nrand       :: Int         = 1,
                           iseed       :: Int         = 0,
                           aim         :: Float64     = 1.0,
                           rate        :: Float64     = 1.0,
                           utol        :: Float64     = 0.1,
                           alpha       :: Float64     = 1.0,
                           biserrtype        :: Symbol = :full,
                           force_extrapolate :: Bool   = false,
                           flux_err          :: Float64 = 0.01,
                           mackay_alpha      :: Bool    = false,
                           ritz_alpha        :: Bool    = false,
                           T           :: Type{<:AbstractFloat} = Float64,
                           fftflags    = FFT_FLAGS)
    nwav = length(data_channels)
    pctx = PolyImagingContext(data_channels, nx, xyint; T, fftflags)

    npix = nx * nx
    nhid = npix * nwav
    ndat = sum(pctx.ndat_ch)

    s = MaximENTState{T}(nhid, ndat)
    p = MaximENTParams(; methd, nrand, iseed, aim, rate, utol, alpha,
                         mackay_alpha, ritz_alpha)

    # Build prior model: per-channel unit-flux normalisation
    model = if isnothing(prior_cube)
        repeat(fill(1.0 / npix, npix), nwav)
    elseif ndims(prior_cube) == 2 && length(prior_cube) == npix
        m = vec(Float64.(prior_cube))
        m ./= max(sum(m), eps())
        repeat(m, nwav)
    else
        m = vec(Float64.(prior_cube))
        @assert length(m) == nhid "prior_cube must have $nhid elements (nx²×nwav), got $(length(m))"
        for w in 1:nwav
            idx = (w-1)*npix+1 : w*npix
            mw = @view m[idx]
            mw ./= max(sum(mw), eps())
        end
        m
    end
    s.model .= model
    s.h     .= model

    set_dataspace_poly!(s, pctx, data_channels;
                        biserrtype, force_extrapolate, flux_err)
    return pctx, s, p
end

# -- Polychromatic reconstruct -----------------------------------------------

"""
    maxent_reconstruct_poly!(pctx, s, p; maxiter=200, verbose=true) → Array{Float64,3}

Run the MaximENT iterative reconstruction loop for polychromatic data and
return the result as an `nx × nx × nwav` array.

Pass a `Vector` as `history` to collect the per-iteration diagnostic trace; see
[`maxent_reconstruct!`](@ref).
"""
function maxent_reconstruct_poly!(pctx :: PolyImagingContext,
                                  s    :: MaximENTState,
                                  p    :: MaximENTParams;
                                  maxiter :: Int  = 200,
                                  verbose :: Bool = true,
                                  history :: Union{Nothing,Vector} = nothing)
    linfwd! = (δim, δdat) -> linearised_fwd_poly!(pctx, δim, δdat)
    linadj! = (δdat, δim) -> linearised_adj_poly!(pctx, δdat, δim)
    fwd!    = (im,  dat)  -> compute_data_poly!(pctx, im, dat)

    # Warm-up: populate each ctx.visi from the starting model
    compute_data_poly!(pctx, s.h, s.d_w2)

    h = reconstruct!(s, p, maxiter, linfwd!, linadj!, fwd!; verbose)
    history === nothing || append!(history, h)

    return reshape(s.h, pctx.nx, pctx.nx, pctx.nwav)
end

# -- Polychromatic reconstruct_bsmem entry points ----------------------------

"""
    reconstruct_bsmem(x_start_3d, data_channels, ft_channels; kwargs...) → Array{Float64,3}

Polychromatic Maximum Entropy reconstruction.  Same keyword arguments as the
monochromatic `reconstruct_bsmem`, but accepts:

- `x_start_3d` — `nx × nx × nwav` starting image cube
- `data_channels` — `Vector{<:OIdata}` (one per wavelength channel)
- `ft_channels` — `Vector` of NFFT plans (one per channel, from a column of `setup_ft`)

Returns an `nx × nx × nwav` image cube normalised to unit flux per channel.
"""
function reconstruct_bsmem(x_start::AbstractArray{<:AbstractFloat,3},
                           data_channels::AbstractVector{<:OIdata},
                           ft_channels::AbstractVector;
                           regularizers      = [],
                           method            = [4, 1, 1, 2],
                           maxiter  :: Int   = 200,
                           verbose  :: Bool  = true,
                           flux_err :: Float64 = 1e-5,
                           biserrtype        :: Symbol = :full,
                           force_extrapolate :: Bool   = false,
                           nrand  :: Int     = 10,
                           iseed  :: Int     = 0,
                           aim    :: Float64 = 1.0,
                           rate   :: Float64 = 1.0,
                           utol   :: Float64 = 0.1,
                           alpha        :: Float64 = 1.0,
                           mackay_alpha :: Bool    = false,
                           ritz_alpha   :: Bool    = false,
                           fftflags               = FFT_FLAGS,
                           history :: Union{Nothing,Vector} = nothing)

    nx   = size(x_start, 1)
    nwav = size(x_start, 3)
    @assert length(data_channels) == nwav "x_start has $nwav channels but $(length(data_channels)) data entries"

    # Extract pixsize from first channel FT plan
    ft1 = ft_channels[1]
    plan1 = ft1 isa AbstractVector ? ft1[1] : ft1
    j = argmax(vec(sqrt.(sum(data_channels[1].uv .^ 2, dims=1))))
    pixsize = (norm(plan1.k[:, j]) / norm(data_channels[1].uv[:, j])) / MAS

    # Parse method
    methd_vec = method isa Integer ? [Int(method), 1, 1, 2] : collect(Int, method)

    # Extract prior from regularizers
    prior_cube = vec(Float64.(x_start))
    for reg in regularizers
        if !isempty(reg) && reg[1] == "mem" && length(reg) >= 2 && !isnothing(reg[2])
            prior_cube = vec(Float64.(reg[2]))
            break
        end
    end

    pctx, s, p = maxent_setup_poly(data_channels, nx, pixsize, prior_cube;
                                   methd=methd_vec, nrand, iseed, aim, rate, utol, alpha,
                                   biserrtype, force_extrapolate, flux_err,
                                   mackay_alpha, ritz_alpha,
                                   T = ft_eltype(ft_channels), fftflags)

    # The solver is Float64 throughout; hand the cube back at the caller's precision.
    return to_ft_precision(maxent_reconstruct_poly!(pctx, s, p; maxiter, verbose, history),
                           ft_channels)
end

"""
    reconstruct_bsmem(x_start_4d, data_matrix, ft_matrix; kwargs...) → Array{<:AbstractFloat,4}

Polychromatic convenience wrapper accepting 4D image and `Matrix{OIdata}`.
Requires single epoch (`size(data,2) == 1`).
"""
function reconstruct_bsmem(x_start::AbstractArray{<:AbstractFloat,4},
                           data::AbstractMatrix{<:OIdata},
                           ft::AbstractMatrix;
                           kwargs...)
    nwav, nepoch = size(data)
    @assert nepoch == 1 "Polychromatic BSMEM currently supports single-epoch only (got nepoch=$nepoch)"
    x3d = x_start[:, :, :, 1]
    data_channels = [data[w, 1] for w in 1:nwav]
    ft_channels   = [ft[w, 1]   for w in 1:nwav]
    result_3d = reconstruct_bsmem(x3d, data_channels, ft_channels; kwargs...)
    return reshape(result_3d, size(result_3d,1), size(result_3d,2), nwav, 1)
end

# -- Matrix convenience methods (auto-dispatch mono vs poly) ------------------

"""
    reconstruct_bsmem(x_start_2d, data::AbstractMatrix{OIdata}, ft; kwargs...)

Convenience method accepting the raw `Matrix{OIdata}` from `readoifits` and
a 2D starting image.  If `size(data,1) == 1` (monochromatic), delegates to the
monochromatic path; otherwise wraps into 3D/Vector form for polychromatic.
"""
function reconstruct_bsmem(x_start::AbstractMatrix{<:AbstractFloat},
                           data::AbstractMatrix{<:OIdata},
                           ft::AbstractMatrix;
                           kwargs...)
    nwav, nepoch = size(data)
    @assert nepoch == 1 "Polychromatic BSMEM currently supports single-epoch only (got nepoch=$nepoch)"
    if nwav == 1
        # Monochromatic — unwrap and delegate
        return reconstruct_bsmem(x_start, data[1,1], ft[1,1] isa AbstractVector ? ft[1,1] : ft[1]; kwargs...)
    else
        # Polychromatic — expand x_start to cube (replicate across channels)
        nx = size(x_start, 1)
        x3d = zeros(nx, nx, nwav)
        m = Float64.(x_start); m ./= max(sum(m), eps())
        for w in 1:nwav
            x3d[:,:,w] .= m
        end
        data_channels = [data[w,1] for w in 1:nwav]
        ft_channels   = [ft[w,1]   for w in 1:nwav]
        return reconstruct_bsmem(x3d, data_channels, ft_channels; kwargs...)
    end
end

# Mono: accept Matrix data/ft with 2D image, ft as Vector (from ft[1])
function reconstruct_bsmem(x_start::AbstractMatrix{<:AbstractFloat},
                           data::AbstractMatrix{<:OIdata},
                           ft::AbstractVector;
                           kwargs...)
    nwav, nepoch = size(data)
    @assert nepoch == 1 "BSMEM currently supports single-epoch only (got nepoch=$nepoch)"
    if nwav == 1
        return reconstruct_bsmem(x_start, data[1,1], ft; kwargs...)
    else
        error("Polychromatic data requires a 3D starting image (nx, nx, nwav)")
    end
end

# -- Prior cube helper -------------------------------------------------------

"""
    gaussian_prior_cube(nx, pixsize_mas, nwav; fwhm_mas=nx*pixsize_mas/5)

Build an `nx × nx × nwav` Gaussian prior image cube with identical spatial
profile in each channel, each normalised to unit flux.
"""
function gaussian_prior_cube(nx::Int, pixsize_mas::Float64, nwav::Int;
                             fwhm_mas::Float64 = Float64(nx) * pixsize_mas / 5.0)
    mono = reshape(gaussian_prior(nx, pixsize_mas; fwhm_mas=fwhm_mas), nx, nx)
    cube = zeros(nx, nx, nwav)
    for w in 1:nwav
        cube[:, :, w] .= mono
    end
    return cube
end

# ─────────────────────────────────────────────────────────────────────────────
# make_operators
# ─────────────────────────────────────────────────────────────────────────────

"""
    make_operators(ctx::ImagingContext) → (linfwd!, linadj!, fwd!)

Return the three operator closures for use with `reconstruct!`:

    linfwd!(δimage, δdata)   – linearised forward operator
    linadj!(δdata, δimage)   – adjoint operator
    fwd!(image,  data)       – nonlinear forward model; updates ctx.visi
"""
function make_operators(ctx::ImagingContext)
    linfwd! = (δim, δdat) -> linearised_fwd!(ctx, δim, δdat)
    linadj! = (δdat, δim) -> linearised_adj!(ctx, δdat, δim)
    fwd!    = (im,  dat)  -> compute_data!(ctx, im, dat)
    return linfwd!, linadj!, fwd!
end

# ─────────────────────────────────────────────────────────────────────────────
# Internal: sign-aware visibility lookup helpers
# ─────────────────────────────────────────────────────────────────────────────

# Forward pass (linearised_fwd): conjugate if sign < 0
@inline _visleg(visi, pnt::Int, sgn::Int8) =
    sgn < 0 ? conj(visi[pnt]) : visi[pnt]

# Adjoint pass (linearised_adj): conjugate if sign ≥ 0  (conj first, then maybe double-conj)
@inline _visleg_adj(visi, pnt::Int, sgn::Int8) =
    sgn < 0 ? visi[pnt] : conj(visi[pnt])

# ─────────────────────────────────────────────────────────────────────────────
# compute_data!  –  nonlinear forward model
# image[1..nx²]  →  data[1..npow+2·nbis]
# ─────────────────────────────────────────────────────────────────────────────

function compute_data!(ctx   :: ImagingContext,
                 image :: AbstractVector{<:AbstractFloat},
                 data  :: AbstractVector{<:AbstractFloat})
    (; npow, nbis, nuv, nx, plan, visi, _f_hat, bs_pnt, bs_sign, phasor, flag_pow, flag_bs) = ctx

    _f_hat .= reshape(image, nx, nx)
    mul!(visi, plan, _f_hat)

    @views data[1:npow] .= flag_pow .* abs2.(visi[1:npow])

    @inbounds for i in 1:nbis
        if flag_bs[i]
            Vab = _visleg(visi, bs_pnt[1,i], bs_sign[1,i])
            Vbc = _visleg(visi, bs_pnt[2,i], bs_sign[2,i])
            Vca = _visleg(visi, bs_pnt[3,i], bs_sign[3,i])
            t   = Vab * Vbc * Vca * phasor[i]
            data[npow + 2i - 1] = real(t)
            data[npow + 2i    ] = imag(t)
        else
            data[npow + 2i - 1] = 0.0
            data[npow + 2i    ] = 0.0
        end
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# linearised_fwd!  –  linearised forward operator
# δimage[1..nx²]  →  δdata[1..npow+2·nbis]
# ─────────────────────────────────────────────────────────────────────────────

function linearised_fwd!(ctx    :: ImagingContext,
                δimage :: AbstractVector{<:AbstractFloat},
                δdata  :: AbstractVector{<:AbstractFloat})
    (; npow, nbis, nx, plan, visi, _f_hat, _dvisi,
       bs_pnt, bs_sign, phasor, flag_pow, flag_bs) = ctx

    _f_hat .= reshape(δimage, nx, nx)
    mul!(_dvisi, plan, _f_hat)

    @views δdata[1:npow] .= flag_pow .* 2 .* real.(_dvisi[1:npow] .* conj.(visi[1:npow]))

    @inbounds for i in 1:nbis
        if flag_bs[i]
            V0ab = _visleg(visi,   bs_pnt[1,i], bs_sign[1,i])
            V0bc = _visleg(visi,   bs_pnt[2,i], bs_sign[2,i])
            V0ca = _visleg(visi,   bs_pnt[3,i], bs_sign[3,i])
            dVab = _visleg(_dvisi, bs_pnt[1,i], bs_sign[1,i])
            dVbc = _visleg(_dvisi, bs_pnt[2,i], bs_sign[2,i])
            dVca = _visleg(_dvisi, bs_pnt[3,i], bs_sign[3,i])
            t = phasor[i] * (dVab*V0bc*V0ca + V0ab*dVbc*V0ca + V0ab*V0bc*dVca)
            δdata[npow + 2i - 1] = real(t)
            δdata[npow + 2i    ] = imag(t)
        else
            δdata[npow + 2i - 1] = 0.0
            δdata[npow + 2i    ] = 0.0
        end
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# linearised_adj!  –  adjoint of linearised_fwd!
# δdata[1..npow+2·nbis]  →  δimage[1..nx²]
#
# Adjoint derivation
# ------------------
# Power-spectrum: linearised_fwd returns 2 Re(δV·V*).
#   Adjoint w.r.t. δV:  2·δdata·V  → accumulated into _dvisi.
#
# Bispectrum: linearised_fwd packs [Re(t), Im(t)], t = φ·(δVab·Vbc·Vca + …).
#   Adjoint undoes the phasor (×conj(φ)) and applies the *reversed*
#   conjugation convention (_visleg_adj vs _visleg):
#   conj(V[pnt]) first and then optionally double-conjugates.
# ─────────────────────────────────────────────────────────────────────────────

function linearised_adj!(ctx    :: ImagingContext,
                δdata  :: AbstractVector{<:AbstractFloat},
                δimage :: AbstractVector{<:AbstractFloat})
    (; npow, nbis, plan, visi, _f_hat, _dvisi,
       bs_pnt, bs_sign, phasor, flag_pow, flag_bs) = ctx

    fill!(_dvisi, zero(eltype(_dvisi)))

    @views _dvisi[1:npow] .+= flag_pow .* 2 .* δdata[1:npow] .* visi[1:npow]

    @inbounds for i in 1:nbis
        if flag_bs[i]
            ab = bs_pnt[1,i];  bc = bs_pnt[2,i];  ca = bs_pnt[3,i]
            sab = bs_sign[1,i]; sbc = bs_sign[2,i]; sca = bs_sign[3,i]

            V0ab = _visleg_adj(visi, ab, sab)
            V0bc = _visleg_adj(visi, bc, sbc)
            V0ca = _visleg_adj(visi, ca, sca)

            bs  = complex(δdata[npow + 2i - 1], δdata[npow + 2i]) * conj(phasor[i])

            dVab = bs * V0bc * V0ca;  sab < 0 && (dVab = conj(dVab))
            dVbc = bs * V0ca * V0ab;  sbc < 0 && (dVbc = conj(dVbc))
            dVca = bs * V0ab * V0bc;  sca < 0 && (dVca = conj(dVca))

            _dvisi[ab] += dVab
            _dvisi[bc] += dVbc
            _dvisi[ca] += dVca
        end
    end

    mul!(_f_hat, adjoint(plan), _dvisi)
    # `vec(real.(_f_hat))` would materialise a whole nx×nx temporary: vec() is a plain
    # call, so it breaks broadcast fusion and forces real.() to allocate. Taking vec()
    # first (a cheap reshape view) lets real.() fuse straight into the store.
    δimage .= real.(vec(_f_hat))
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# Utility: set_mask!  (port of SETMSK in bsmemf.f)
# ─────────────────────────────────────────────────────────────────────────────

"""
    set_mask!(f, j, k)

Set vector `f` to a box support: 1 for `j ≤ i ≤ k`, 0 elsewhere (1-based).
"""
function set_mask!(f::AbstractVector, j::Int, k::Int)
    fill!(f, zero(eltype(f)))
    f[j:k] .= one(eltype(f))
    return f
end

# ─────────────────────────────────────────────────────────────────────────────
# Utility: log_message  (port of UINFO in bsmemf.f)
# ─────────────────────────────────────────────────────────────────────────────

"""
    log_message(msg, mlevel; level=12)

Print `msg` to stdout if `mlevel` is within the verbosity gate `level`.

Two-digit encoding:
  tens digit  : threshold for numerical diagnostics   (mlevel ≥ 10)
  units digit : threshold for progress diagnostics    (mlevel  < 10)
"""
function log_message(msg::AbstractString, mlevel::Int; level::Int=12)
    if mlevel >= 10
        mlevel > level  && return
    else
        (mlevel % 10) > (level % 10) && return
    end
    println(msg)
    return nothing
end

