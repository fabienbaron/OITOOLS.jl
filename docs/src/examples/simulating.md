# Simulating observations

OITOOLS writes synthetic OIFITS files from an image or a parametric model, either on the uv
coverage of an existing file or from an array, a combiner and a list of observation times.

## From an existing OIFITS file

[`simulate_from_oifits`](@ref) keeps the uv points, wavelengths, flags and headers of an input
file and replaces the observables:

```julia
using OITOOLS
simulate_from_oifits("data/Beauty_Contest_2004/2004-data1.oifits", "data/sim.oifits";
                     image="data/images/2004true.fits", pixsize=0.101)
```

A flat-dict parametric model can be used instead of an image:

```julia
params = Dict("star,ud" => 3.0, "star,f" => 1.0)
model = dict_to_model(params, String[])

simulate_from_oifits("data/Beauty_Contest_2004/2004-data1.oifits", "data/sim.oifits";
                     flat_model=model, flat_params=Float64[])
```

`mode` sets the error bars:

| `mode` | error bars |
|---|---|
| `"copy_errors"` (default) | the input's |
| `"copy_snr"` | the input's signal-to-noise, applied to the new observables |
| `"noise_model"` | from `errors` (`v2_multit`, `v2_addit` and the T3 equivalents) |

`noise = false` writes the noiseless observables with the error bars filled in. The output path
must differ from the input. See `demos/example_simulate_observations_from_OIFITS.jl`.

## From observation times and array geometry

[`simulate`](@ref) takes four configurations: facility, target, combiner and wavelength setup.
Built-in ones are TOML files in `src/configs/`.

### Configuration files

Built-in configurations are read by name, without the `.toml` extension.

**Facility**: array layout, telescope positions, atmospheric conditions.

```julia
facility = read_facility_file("CHARA")
```

| Config name | Interferometer | Telescopes |
|---|---|---|
| `CHARA` | CHARA array | 6×1 m |
| `VLTI_UT` | VLTI Unit Telescopes | 4×8.2 m |
| `VLTI_AT_small` | VLTI ATs, small configuration | 4×1.8 m |
| `VLTI_AT_medium` | VLTI ATs, medium configuration | 4×1.8 m |
| `VLTI_AT_large` | VLTI ATs, large configuration | 4×1.8 m |

**Target**: coordinates and proper motion.

```julia
target = read_obs_file("default_obs")
```

From SIMBAD:

```julia
ra, dec = ra_dec_from_simbad("Vega")   # decimal degrees, matching TargetConfig.raep0
target = TargetConfig(target="Vega", raep0=ra, decep0=dec)
```

[`simbad_target`](@ref) returns the other `TargetConfig` fields and the magnitudes used by the
noise model:

```julia
t = simbad_target("Vega")
target = TargetConfig(target = t.main_id, raep0 = t.ra, decep0 = t.dec,
                      pmra = t.pmra, pmdec = t.pmdec,
                      parallax = t.plx, spectyp = t.sptype)

simulate(facility, target, combiner, wavelength, dates, "vega.oifits";
         mag = t.mags["H"], mag_ao = t.mags["V"])
```

Bands without a SIMBAD measurement are `NaN`.

**Combiner**: transmission, detector, calibration errors.

```julia
combiner = read_comb_file("MIRCX")
```

| Config name | Instrument | Array | Band |
|---|---|---|---|
| `GRAVITY` | GRAVITY | VLTI | K |
| `MATISSE_LM` | MATISSE | VLTI | L+M |
| `MATISSE_N` | MATISSE | VLTI | N |
| `MIRCX` | MIRC-X | CHARA | H |
| `MYSTIC` | MYSTIC | CHARA | K |
| `SPICA` | SPICA | CHARA | V |

**Wavelength**: the spectral setup of a combiner mode.

```julia
wave = read_wave_file("MIRCX_LOWH")
```

| Config name | Combiner | Mode | Band |
|---|---|---|---|
| `GRAVITY_LOWK` | GRAVITY | Low spectral resolution | K |
| `MATISSE_LOWL` | MATISSE_LM | Low spectral resolution | L |
| `MATISSE_LOWN` | MATISSE_N | Low spectral resolution | N |
| `MIRCX_HIGHH` | MIRCX | High spectral resolution (R ≈ 1035, 230 channels) | H |
| `MIRCX_LOWH` | MIRCX | Low spectral resolution | H |
| `MIRCX_LOWJ` | MIRCX | Low spectral resolution | J |
| `MIRCX_MEDIUMH` | MIRCX | Medium spectral resolution (R ≈ 190, 42 channels) | H |
| `MYSTIC_LOWK` | MYSTIC | Low spectral resolution | K |
| `SPICA_LR` | SPICA | Low resolution | V |

### Simulating from an image

```julia
using Dates

# Observation times: every 15 minutes over a 5.5-hour window
dates = collect(DateTime(2024,8,13,3,0,0):Minute(15):DateTime(2024,8,13,8,30,0))

facility = read_facility_file("CHARA")
target   = read_obs_file("default_obs")
combiner = read_comb_file("MIRCX")
wave     = read_wave_file("MIRCX_LOWH")

simulate(facility, target, combiner, wave, dates, "sim_image.oifits";
         image="data/images/2004true.fits", pixsize=0.101)
```

### Simulating from a parametric model

```julia
params = Dict(
    "star,ud"    => 3.0,
    "star,f"     => 0.7,
    "disk,f"     => "1 - \$star,f",
    "disk,diamout" => 10.0,
    "disk,profile" => "exp(-(\$R/3.0)^2)",
)
model = dict_to_model(params, String[])

simulate(facility, target, combiner, wave, dates, "sim_model.oifits";
         flat_model=model, flat_params=Float64[])
```

See `demos/example_simulate_observations_from_model.jl` and
`demos/example_simulate_observations_from_image.jl`.

### Polychromatic simulation

`demos/example_simulate_polychromatic_disk.jl` simulates a chromatic, time-variable disc with
an off-axis companion, and writes OI_VIS2, OI_VIS (differential phases), OI_T3 and OI_FLUX.

A FITS file is read by its keywords. Each uv point is observed with the plane of its channel and
its date, so the cube's wavelengths and epochs must cover the spectral setup's channels and the
dates; a channel with no plane is an error that lists the cube's wavelengths. A file whose
keywords do not give its axes is refused; set them in the GUI (Tools › FITS Image Handling, Set
keywords, Save FITS) or with [`fits_set_keywords`](@ref). Each plane of a cube is normalised to
unit flux for the visibilities; the plane totals weight the photon count per channel and fill
`OI_FLUX`. A model may contain an image component (`"env,image" => path`); a SPARCO star plus
image is passed as `flat_model`.

`wavelength = :image` takes the spectral setup from the cube's channels.

### A model cube from elsewhere: AMHRA's Be star

AMHRA's "Be model star + Keplerian disk" (Delaa et al. 2011) returns a cube with the pixel scale
in radians, the wavelengths in metres (`CUNIT3 = 'm'`) and no `CTYPE3`. Its keywords are set,
keeping its reference pixel, and the cube is simulated on its own channels
(`demos/example_simulate_amhra_cube.jl`, PMOIRED's example 4 observed with CHARA):

```julia
fits_set_keywords("output_Be.fits", "be_cube.fits"; pixsize = 2.65845e-10 * 180 / π * 3.6e6,
                  wavelengths = [2.1661e-6 + (k - 25) * 2.7e-10 for k in 1:51])
simulate(read_facility_file("CHARA"), TargetConfig(raep0 = 0.0, decep0 = -24.0),
         read_comb_file("MYSTIC"), :image, dates, "be.oifits"; image = "be_cube.fits",
         noise = FractionalNoise(v2 = 1e-4, visamp = 1e-4, t3amp = 1e-4, visphi = 0.01,
                                 t3phi = 0.01, flux = 1e-4))
```

Without `CTYPE1` the cube is read as stored, first axis East, as PMOIRED's `readCube` does; the
two agree on |V| to 4e-5, closure phases to 0.005° and the flux spectrum to 4e-5. This puts the
continuum disc at PA 120° (−60°). AMHRA's input PA of 60° corresponds to `east_left = true`.

[`FractionalNoise`](@ref) sets the error bars as PMOIRED's `oifake`: fractions of V², |V|, T3
amplitude and flux, phase errors in degrees, each observable drawn independently.

## Keyword arguments to `simulate`

| keyword | default | meaning |
|---|---|---|
| `image` / `pixsize` | `""` / from the file | truth image: a FITS path read by its keywords (2-D, wavelength cube, time series or 4-D; [`read_model_image`](@ref)), or an array (2-D, or 3-D one plane per channel); `pixsize` in mas overrides the header (an array's default is 0.1) |
| `flat_model` / `flat_params` | `nothing` | parametric model instead of an image |
| `mag` | `2.0` | target magnitude: a number, a `Dict("H"=>1.8, …)` of band magnitudes, or one value per spectral channel |
| `mag_ao` | from `mag` | guide-star magnitude in the AO wavefront-sensor band |
| `noise` | `true` | the photon-noise model; `false` writes the model with its computed error bars; a [`FractionalNoise`](@ref) fixed relative errors (PMOIRED's `oifake`) |
| `debias` | `true` | subtract the `2σ²` bias from `V²` |
| `systematics` | `true` | draw the calibration error, correlated per baseline per night (see below) |
| `fringe_tracker` | `true` | use the instrument's fringe tracker, if it has one |
| `n_samples` | `100` | Monte-Carlo samples for the T3 error bars; `0` uses the analytic form |
| `seed` / `rng` | `nothing` | reproducible noise realisation |
| `dft` | `false` | DFT instead of NFFT for an array `image` |
| `observability` | `nothing` | observability filter, see below |
| `nonoise` | — | deprecated spelling of `noise=false` |

Coordinates follow the OIFITS standard: `target.raep0` and `target.decep0` are in **degrees**.

## Noise model

Photons per telescope, per spectral channel, per frame:

```
N = F0(λ) · 10^(-m(λ)/2.5) · A_tel · δλ · DIT
    · T_atm(λ)             atmospheric transmission (`atm_transmission`, 1.0 by default)
    · facility.throughput  telescopes and beam train only
    · combiner.transmission   array + instrument, excluding QE, Strehl and atmosphere
    · flux_frac            split between interferometric and photometric channels
    · QE
    · S(λ, elevation, m_ao)   Strehl ratio, or min(1,(r₀/D)²) without AO
```

### Matching ASPRO

The CHARA combiner configs (`MIRCX.toml`, `MYSTIC.toml`, `SPICA.toml`) are transcribed from
ASPRO 2's `aspro-conf/…/CHARA.xml`. `CombinerConfig` fields correspond to ASPRO's
`FocalInstrumentSetup`:

| ASPRO XML | `CombinerConfig` | note |
|---|---|---|
| `transmission` | `transmission` | array and instrument; divided by the atmosphere, see below |
| `instrumentVisibility` | `instrument_visibility` | |
| `dit` | `dit` | fixed; shortened only to avoid saturation |
| `defaultTotalIntegrationTime` | `total_int_time` | |
| `detectorSaturation` | `detector_saturation` | |
| `ron`, `quantumEfficiency` | `read_noise`, `quantum_efficiency` | QE defaults to 1.0 in both |
| `nbPixInterferometry` / `nbPixPhotometry` | `n_pix_fringe` / `n_pix_photometry` | |
| `fracFluxInInterferometry` / `…Photometry` | `flux_frac_fringes` / `flux_frac_photometry` | |
| `instrumentVisibilityBias` | `vis_cal_err` | ASPRO stores percent: `1` → `0.01`; the V² systematic is twice this |
| `instrumentPhaseBias` | `phase_cal_err` | degrees |

`FacilityConfig.throughput` is 1.0 for CHARA, since `transmission` covers the array.

!!! note "Transmission values"
    ASPRO's `ConfigurationManager` divides each setup's `<transmission>` by the band-mean
    atmospheric transmission (`includeAtmosphereCorrection`, true by default, not overridden for
    CHARA). The configs store the divided value: MIRC-X's `0.01` becomes
    `0.01 / 0.8219 = 0.0121671`. The divisor is the mean over the MIRC-X/MYSTIC J-to-K range, so
    MYSTIC uses the same one; SPICA's is 0.9333. The XML value without the division gives
    σ(V²) 1.44× larger at H = 8.

Differences from ASPRO:

- **SPICA's Strehl.** ASPRO's `NoiseService` uses fixed Strehl ratios for SPICA (0.25 / 0.15 /
  0.10 by seeing), while its Strehl plots use the AO model. `SPICA.toml` sets
  `strehl_model = "fixed_spica"` to match ASPRO's noise; `strehl_model = "ao"` uses the AO model
  (about 0.28 at 1″).
- **Telluric absorption.** `atm_transmission` returns 1.0 at all wavelengths; ASPRO uses an ESO
  SkyCalc table (computed for Paranal). Per-channel σ is therefore smaller inside telluric bands:
  at H = 8 MIRC-X agrees to 0.99 on the median and 0.44 in the H₂O channel at 1.4152 µm.
- **Photometric zero point.** `zero_point_flux` interpolates between band centres; ASPRO uses
  one zero point per band. Channel-by-channel ratios therefore differ from 1.00 above R ≈ 150.

**Fringe trackers.** A combiner config with a `[fringe_tracker]` table integrates up to the
tracker's `max_integration` instead of `dit`, multiplied by the tracker's visibility loss:

```toml
[fringe_tracker]      # SPICA.toml; last in the file, as a TOML table takes every key after it
band = "H"            # the tracker's band
mag_limit = 10.0
max_integration = 0.2 # seconds, against SPICA's 20 ms `dit`
instrument_visibility = 0.8
mode = "FringeTrack"  # "GroupTrack": visibility loss without the longer integration
```

As in ASPRO, the tracker is used when the target is within `mag_limit` in the tracker's band
and at least one frame fits before saturation. `fringe_tracker = false` turns it off. SPICA is
the only shipped combiner with a tracker; ASPRO marks it as required.

`S` comes from [`strehl_ratio`](@ref), a port of JMMC's `Band.strehl`, which reproduces ASPRO
2's CHARA Strehl curves to a median 0.7%. It needs an `[ao]` table in the facility config;
without one `S = min(1,(r₀/D)²)`, about 5× lower than the AO value in H and 20× in R.

The photon noise is drawn on the complex visibility, and all observables are computed from it,
so `VISAMP² == VIS2`.

### Calibration errors

`systematics = true` draws `vis_cal_err` and `phase_cal_err` as correlated errors, in three
terms:

| term | drawn once per |
|---|---|
| amplitude gain, on the complex visibility (V², VISAMP, T3AMP) | baseline × night |
| closure-phase bias | triangle × night |
| differential-phase bias | baseline × night × channel |

A night runs from local noon to local noon. The closure phase includes its own bias, so it is
not the sum of the three baseline phases. When the systematics dominate, χ²ᵣ of the true model
scatters by `√(2/n_baselines)`. `systematics = false` puts the calibration errors in the error
bars only, without drawing them.

T3AMP and T3PHI error bars are estimated from `n_samples` Monte-Carlo draws (default 100);
`n_samples = 0` uses the first-order analytic form. Reduced χ² of the true model on a resolved
disc:

| median SNR(V²) | 6.3 | 1.0 | 0.09 | 0.02 |
|---|---|---|---|---|
| T3AMP, analytic | 1.04 | 2.16 | 27.6 | 159.6 |
| T3AMP, sampled  | 1.01 | 1.02 | 1.03 | 1.03 |
| T3PHI, analytic | 1.02 | 0.96 | 0.58 | 0.57 |
| T3PHI, sampled  | 1.01 | 1.01 | 1.00 | 1.00 |

### Checking it

Run `demos/validate_noise_model.jl` to print the Strehl comparison against ASPRO and the
predicted σ(V²)/σ(CP) against magnitude for MIRC-X, MYSTIC and SPICA.

`test/test_aspro_noise.jl` pins the model against ASPRO 2's own error bars, stored per channel
in `test/references/aspro_noise.csv` and regenerated with `tools/AsproNoise.java` (ASPRO run
under Xvfb) and `tools/compare_errors.jl`. Ratio of σ(V²) to ASPRO's, at the best elevation of
the night:

| photons/tel/DIT | MIRC-X Low_H | MYSTIC Low_K | SPICA LR |
|---|---|---|---|
| ≥ 60 (mag 0–4) | 1.00 | 1.00 | 1.00 |
| ~10 (mag 6) | 1.00 | 1.00 | 1.02 |
| ~1.5 (mag 8) | 0.99 | 1.02 | 1.06 |
| ≲ 1 (mag 10) | — | 1.04 | 1.16 |

On bright targets σ is set by the calibration floor; the faint rows test the photon budget.

## Observability filtering

By default `simulate` uses every epoch passed to it, whatever the target's elevation or the
delay lines. To filter the epochs, call [`observable_epochs`](@ref):

```julia
dates_ok, mask, report = observable_epochs(facility, target, dates;
                                           min_elevation = 30.0,
                                           pops = [1,3,5,2,4,1])   # from best_pop
simulate(facility, target, combiner, wave, dates_ok, "sim.oifits"; flat_model=model)
```

or pass the same options to `simulate`:

```julia
simulate(facility, target, combiner, wave, dates, "sim.oifits";
         flat_model=model, observability=(min_elevation=30.0,))
```

Without `pops` there is no delay-line check. [`best_pop`](@ref) proposes POP configurations.

## Observation planning

`obs_plan` (needs PythonPlot) computes the night and draws the Gantt chart:

```julia
facility = read_facility_file("CHARA")
ra, dec  = ra_dec_from_simbad("Vega")        # decimal degrees
config   = [1, 1, 1, 1, 1, 2]                # 0=unused, 1=use, 2=reference cart
pop      = [1, 1, 1, 1, 1, 1]                # one POP per telescope, 1:5

obs_plan("Vega", facility, ra, dec, DateTime(2026, 6, 3), pop, config;
         alt_limit = 30.0, savefile = "vega.png")
```

The underlying functions:

```julia
# Dark window, in decimal UT hours. `zenith` is in degrees: 102 is nautical twilight.
dusk_rise, dusk_set = sunrise_sunset(DateTime(2026, 6, 3), facility.lat, facility.lon)

# Everything for one night: LST, hour angle, altitude, azimuth, Moon separation
obs = night_observability(facility, ra, dec, DateTime(2026, 6, 3); alt_limit = 30.0)

# Altitude and azimuth in degrees; arguments (dec, lat, ha), `ha` in hours
altitude, azimuth = alt_az(dec, facility.lat, obs.ha)

# Delay-line feasibility for one POP choice, and the best POP choices
d = in_delay(facility, dec, obs.ha, config, pop)
results = best_pop(facility, dec, obs.ha, config; n_best = 5)
print_pop_results(facility, config, results)
```

See `demos/example_chara_plan.jl`.
