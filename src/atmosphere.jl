# atmosphere.jl — photometric bands, atmospheric seeing, and the AO Strehl model.
#
# The Strehl model is a port of JMMC's `Band.strehl` from jmal
# (https://github.com/JMMC-OpenDev/jmal, src/main/java/fr/jmmc/jmal/Band.java), which is what
# ASPRO 2's NoiseService uses to set the flux entering an interferometric combiner. Following
# it exactly is deliberate: it lets us reproduce ASPRO's published sensitivity curves, which
# is the only independent check available for CHARA.
#
# jmal cites le Louarn et al. 1998, MNRAS 295, 756 and AMB-IGR-011 p.5.

# ─── Photometric bands ────────────────────────────────────────────────────────

"""
    PhotometricBand

One photometric band: its name, central wavelength and width in metres, the photon flux of a
zero-magnitude star (`f0`, photons/s/m²/m), and the default Strehl ceiling jmal applies in
that band. See [`PHOTOMETRIC_BANDS`](@ref).
"""
struct PhotometricBand
    name::String
    lambda::Float64      # central wavelength [m]
    width::Float64       # bandwidth [m]
    f0::Float64          # zero-magnitude flux [photons / s / m^2 / m]
    strehl_max::Float64  # default Strehl ceiling for this band (jmal per-band table)
end

# Zero points derived from the Bessell et al. (1998) / Cohen et al. (2003) f_lambda values for
# a 0-magnitude star, converted to photon flux as n = f_lambda / (hc/lambda) and expressed in
# photons/s/m^2/m:
#
#     band    f_lambda [erg/cm^2/s/A]    n [ph/s/cm^2/A]    F0 [ph/s/m^2/um]
#     V       3.631e-9                   1005               1.005e11
#     H       1.138e-10                    93.4             9.338e9
#     K       3.961e-11                    43.7             4.367e9
#
# The V figure is the textbook "a zeroth-magnitude star delivers ~1000 photons/s/cm^2/A".
#
# NOTE: this supersedes an earlier table in simulate.jl whose values were low by 2.8x at V,
# growing to 52x at N. That error propagated straight into every simulated error bar, and it
# is independently confirmed here by the AO knee position: reproducing ASPRO's Strehl-vs-V-mag
# curve requires ~17 photons per subaperture per WFS frame at V=10.5, which these values give
# (17.0) and the old ones did not (6.2).
"""
    PHOTOMETRIC_BANDS

The photometric bands known to OITOOLS, from U to Q, as [`PhotometricBand`](@ref) entries.

Zero points are derived from the Bessell et al. (1998) / Cohen et al. (2003) `f_lambda`
values for a 0-magnitude star, converted to photon flux. The V entry is the textbook
"a zeroth-magnitude star delivers ~1000 photons/s/cm²/Å".
"""
const PHOTOMETRIC_BANDS = [
    PhotometricBand("U",    0.360e-6, 0.066e-6, 7.576e10 * 1e6, 0.40),
    PhotometricBand("B",    0.440e-6, 0.094e-6, 1.400e11 * 1e6, 0.45),
    PhotometricBand("V",    0.550e-6, 0.088e-6, 1.005e11 * 1e6, 0.50),
    PhotometricBand("R",    0.640e-6, 0.138e-6, 7.014e10 * 1e6, 0.65),
    PhotometricBand("G_rp", 0.750e-6, 0.280e-6, 4.833e10 * 1e6, 0.70),
    PhotometricBand("I",    0.790e-6, 0.149e-6, 4.478e10 * 1e6, 0.75),
    PhotometricBand("J",    1.220e-6, 0.300e-6, 1.933e10 * 1e6, 0.77),
    PhotometricBand("H",    1.630e-6, 0.350e-6, 9.338e9  * 1e6, 0.84),
    PhotometricBand("K",    2.190e-6, 0.410e-6, 4.367e9  * 1e6, 0.93),
    PhotometricBand("L",    3.450e-6, 0.550e-6, 1.230e9  * 1e6, 0.972),
    PhotometricBand("M",    4.750e-6, 0.300e-6, 4.902e8  * 1e6, 0.985),
    PhotometricBand("N",   10.400e-6, 5.000e-6, 5.026e7  * 1e6, 0.996),
    PhotometricBand("Q",   20.100e-6, 8.000e-6, 6.577e6  * 1e6, 0.999),
]

"""
    band_for_wavelength(λ)

The `PhotometricBand` whose central wavelength is closest to `λ` (in metres).
"""
function band_for_wavelength(λ)
    best = PHOTOMETRIC_BANDS[1]
    bestd = abs(log(λ / best.lambda))
    for b in PHOTOMETRIC_BANDS
        d = abs(log(λ / b.lambda))
        if d < bestd
            bestd = d
            best = b
        end
    end
    return best
end

"""
    band_by_name(name)

Look up a `PhotometricBand` by name (case-insensitive), or `nothing`.
"""
function band_by_name(name::AbstractString)
    key = uppercase(String(name))
    for b in PHOTOMETRIC_BANDS
        uppercase(b.name) == key && return b
    end
    return nothing
end

# Zero-point flux for a mag=0 star in photons/s/m^2/um, log-linear in log(λ) between the
# tabulated band centres. Kept in these units for backwards compatibility.
const _F0_λ_μm = [b.lambda * 1e6 for b in PHOTOMETRIC_BANDS[3:end]]     # V .. Q
const _F0_phot = [b.f0 * 1e-6    for b in PHOTOMETRIC_BANDS[3:end]]

"""
    zero_point_flux(λ)

Photon flux of a zero-magnitude star at wavelength `λ` (metres), in photons/s/m²/µm.

Interpolated log-linearly in `log(λ)` between the band centres of [`PHOTOMETRIC_BANDS`](@ref),
and clamped to the V and Q endpoints outside that range.
"""
function zero_point_flux(λ_meters)
    λ_μm = clamp(λ_meters * 1e6, _F0_λ_μm[1], _F0_λ_μm[end])
    logλ = log.(_F0_λ_μm)
    logF = log.(_F0_phot)
    for i in 1:length(_F0_λ_μm)-1
        if λ_μm <= _F0_λ_μm[i+1]
            t = (log(λ_μm) - logλ[i]) / (logλ[i+1] - logλ[i])
            return exp(logF[i] + t * (logF[i+1] - logF[i]))
        end
    end
    return _F0_phot[end]
end

# ─── Atmospheric transmission hook ────────────────────────────────────────────

"""
    atm_transmission(λ)

Fraction of light transmitted by the atmosphere at wavelength `λ` (metres). Defaults to 1.0
at every wavelength.

ASPRO convolves an ESO SkyCalc table (`skytable_mean_atm.fits`, computed for Paranal) over
each spectral channel. That table is not redistributable and hand-rolling telluric bands
would be worse than saying nothing, so this is a hook: pass a `transmission` column in a
wavelength config, or redefine this, to model tellurics. Bands to care about at CHARA are
O2 A at 760 nm and H2O at 820-950 nm (SPICA), and H2O at 1350-1420 nm (MIRC-X J).
"""
atm_transmission(λ) = one(float(typeof(λ)))

# ─── Seeing / turbulence ──────────────────────────────────────────────────────

# jmal uses the Airy-FWHM constant 1.029 (not 0.98) in the seeing -> r0 conversion, and works
# in arcseconds; 206264.806 is arcsec per radian.
const _R0_FACTOR = 1.0289939700094716812e-6 * 206264.806   # = 0.212245 m x arcsec / um

"""
    fried_parameter(seeing, λ; elevation_deg=90.0)

Fried parameter r0 in metres at wavelength `λ` (metres) for a given `seeing` (arcsec, defined
at 500 nm), scaled to the line of sight at `elevation_deg`.

`r0 ∝ λ^(6/5) · airmass^(-3/5)`.
"""
function fried_parameter(seeing, λ; elevation_deg=90.0)
    r0_V = _R0_FACTOR * (0.5 / seeing) * airmass(elevation_deg)^(-3/5)
    return r0_V * (λ * 1e6 / 0.5)^(6/5)
end

"""
    seeing_from_r0(r0)

Seeing in arcsec (at 500 nm) corresponding to a Fried parameter `r0` in metres, using the
same constant as [`fried_parameter`](@ref).
"""
seeing_from_r0(r0) = _R0_FACTOR * 0.5 / r0

# ─── Adaptive optics ──────────────────────────────────────────────────────────

"""
    AOConfig

Adaptive-optics parameters for a facility's telescopes, consumed by [`strehl_ratio`](@ref).

Defaults are ASPRO 2's `AO_CHARA` setup as of the CHARA 2026A configuration.

- `band`: wavefront-sensor band. "V" is special-cased to a 0.5 um reference wavelength.
- `n_subpupils`: subapertures, sets the WFS photon collecting area.
- `n_actuators`: corrected modes, sets the fitting error. Falls back to `n_subpupils` if 0.
- `dit`: WFS integration time [s].
- `ron`: WFS detector read noise [e-].
- `qe`: WFS quantum efficiency.
- `transmission`: fraction of starlight diverted to the WFS, i.e. lost to the science
  channel. 0 means the WFS is fed by a separate path.
- `mag_offset`: added to the guide-star magnitude before computing WFS flux.
- `strehl_max`: Strehl ceiling. **0 does not mean "uncapped"** — it selects the per-band
  default from [`PHOTOMETRIC_BANDS`](@ref), matching ASPRO.
"""
Base.@kwdef struct AOConfig
    band::String          = "V"
    n_subpupils::Int      = 36
    n_actuators::Int      = 24
    dit::Float64          = 2.0e-3
    ron::Float64          = 5.0
    qe::Float64           = 0.7
    transmission::Float64 = 0.0
    mag_offset::Float64   = 0.0
    strehl_max::Float64   = 0.0
end

"""
    strehl_ratio(λ, D, seeing, t0, mag_ao, ao::AOConfig; elevation_deg=90.0)

Strehl ratio at wavelength `λ` (metres) for a telescope of diameter `D` (metres), given
`seeing` (arcsec at 500 nm), atmospheric coherence time `t0` (seconds), and adaptive optics
guiding on a star of magnitude `mag_ao` in the WFS band.

Port of jmal's `Band.strehl`. Five variance terms are summed and the result is the sum of a
diffraction-limited core and a seeing-limited halo:

    S = exp(-σ²) + (1 - exp(-σ²)) / (1 + (D/r0_λ)²)

The halo is why the Strehl does not collapse to zero on faint guide stars: as σ² → ∞ the
Strehl tends to `1/(1 + (D/r0_λ)²)`, the relative peak intensity of an uncorrected
long-exposure PSF. Reproduces ASPRO's CHARA curves to a median 0.7% (worst case 3.1% in H/K).

Two details that are easy to get wrong and are deliberate here: the servo term uses `t0`
**unscaled** — neither shifted to the science wavelength nor corrected for airmass — and the
subaperture pitch is `sqrt(π(D/2)²/n_actuators)`, not `D/sqrt(n_actuators)`.
"""
function strehl_ratio(λ, D, seeing, t0, mag_ao, ao::AOConfig; elevation_deg=90.0)
    n_act = ao.n_actuators > 0 ? ao.n_actuators : ao.n_subpupils
    telarea = π * (D/2)^2

    # WFS reference wavelength: V is special-cased to 500 nm, as in jmal.
    ao_band = band_by_name(ao.band)
    λ_ao = (ao_band === nothing || uppercase(ao.band) == "V") ? 0.5e-6 : ao_band.lambda
    wfs_band = ao_band === nothing ? band_by_name("V") : ao_band

    r0_λ = fried_parameter(seeing, λ; elevation_deg=elevation_deg)

    ds2 = telarea / ao.n_subpupils          # subaperture area [m^2]
    ds  = sqrt(telarea / n_act)             # actuator pitch  [m]
    ratio = λ / λ_ao

    # Photons per subaperture per WFS frame from the guide star.
    nphot = ao.qe * wfs_band.f0 * wfs_band.width * ao.dit * 10.0^(-0.4 * (mag_ao + ao.mag_offset)) * ds2

    # 1. fitting + aliasing, at the science wavelength
    σ2 = 0.54 * (ds / r0_λ)^(5/3)

    # 2. servo lag. t0 is used as given: not scaled to λ, not airmass-corrected.
    σ2 += 0.962 * (ao.dit / t0)^(5/3)

    if nphot > 0
        # 3. WFS photon noise
        σ2 += (4π^2 / 3) * ratio^(-2) / nphot
        # 4. WFS read noise
        σ2 += (8π^2 / 9) * ratio^(-2) * (ao.ron / nphot)^2 *
              (1 + ratio^(12/5) * (ds/r0_λ)^2)^2 * ratio^(-12/5) * (ds/r0_λ)^(-2)
    else
        σ2 = Inf
    end

    # 5. fixed term: a Strehl ceiling. strehl_max == 0 selects the per-band default.
    smax = ao.strehl_max > 0 ? ao.strehl_max : band_for_wavelength(λ).strehl_max
    σ2 += -log(smax)

    core = exp(-σ2)
    halo = 1.0 / (1.0 + (D / r0_λ)^2)
    return clamp(core + (1 - core) * halo, 0.0, 1.0)
end

"""
    coupling_efficiency(λ, D, seeing, t0, mag_ao, ao; elevation_deg=90.0)

Fraction of the light collected by one telescope that reaches the combiner, i.e. the Strehl
ratio when an [`AOConfig`](@ref) is available and the seeing-limited estimate `min(1,(r0/D)²)`
when `ao === nothing`.

The fallback is what this package used unconditionally before the AO model existed. It is
appropriate only for an uncorrected aperture and underestimates a modern AO-equipped array by
roughly 5x in H and 20x in R at CHARA, so facilities that have AO should carry an `[ao]`
block in their config.
"""
function coupling_efficiency(λ, D, seeing, t0, mag_ao, ao::Union{AOConfig,Nothing}; elevation_deg=90.0)
    if ao === nothing
        @info("No [ao] block in the facility config: modelling the telescopes as " *
              "uncorrected apertures with coupling min(1,(r0/D)^2). For an AO-equipped " *
              "array this underestimates the flux by roughly 5x in H and 20x in R.",
              maxlog=1)
        r0_λ = fried_parameter(seeing, λ; elevation_deg=elevation_deg)
        return min(1.0, (r0_λ / D)^2)
    end
    return strehl_ratio(λ, D, seeing, t0, mag_ao, ao; elevation_deg=elevation_deg)
end
