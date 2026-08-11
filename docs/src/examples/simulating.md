# Simulating observations

OITOOLS can generate synthetic OIFITS datasets from a parametric model or
image, either by copying the UV coverage of an existing file or by building
observations from scratch using array geometry and observation times.

## From an existing OIFITS file

The simplest approach reuses the UV coverage and noise properties of a real
dataset:

```julia
using OITOOLS
simulate_from_oifits("data/2004-data1.oifits", "data/sim.oifits";
                     image="data/2004true.fits", pixsize=0.101)
```

A flat-dict parametric model can be used instead of an image:

```julia
params = Dict("star,ud" => 3.0, "star,f" => 1.0)
model = dict_to_model(params, String[])

simulate_from_oifits("data/2004-data1.oifits", "data/sim.oifits";
                     flat_model=model, flat_params=Float64[])
```

See `example_simulate_observations_from_OIFITS.jl`.

## From observation times and array geometry

To simulate a full night of observations at a specific interferometer, you
need four configuration objects: facility, target, combiner, and wavelength
setup. These are read from TOML files shipped with OITOOLS in `src/configs/`.

### Configuration files

The `.toml` extension is optional — OITOOLS resolves built-in config names
automatically.

**Facility** — array layout, telescope positions, atmospheric conditions:

```julia
facility = read_facility_file("CHARA")
```

| Config name | Interferometer | Telescopes |
|---|---|---|
| `CHARA` | CHARA array | 6×1 m |
| `VLTI_UT` | VLTI Unit Telescopes | 4×8.2 m |
| `VLTI_AT_small` | VLTI ATs — small config | 4×1.8 m |
| `VLTI_AT_medium` | VLTI ATs — medium config | 4×1.8 m |
| `VLTI_AT_large` | VLTI ATs — large config | 4×1.8 m |

**Target** — celestial coordinates and proper motion:

```julia
target = read_obs_file("default_obs")
```

You can also query SIMBAD directly:

```julia
ra, dec = ra_dec_from_simbad("Vega")
target = TargetConfig(target="Vega", raep0=ra, decep0=dec)
```

**Combiner** — beam combiner properties (throughput, read noise, calibration errors):

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

**Wavelength** — spectral setup for a given combiner mode:

```julia
wave = read_wave_file("MIRCX_LOWH")
```

| Config name | Combiner | Mode | Band |
|---|---|---|---|
| `GRAVITY_LOWK` | GRAVITY | Low spectral resolution | K |
| `MATISSE_LOWL` | MATISSE_LM | Low spectral resolution | L |
| `MATISSE_LOWN` | MATISSE_N | Low spectral resolution | N |
| `MIRCX_LOWH` | MIRCX | Low spectral resolution | H |
| `MIRCX_LOWJ` | MIRCX | Low spectral resolution | J |
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
         image="data/2004true.fits", pixsize=0.101)
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

See `example_simulate_observations_from_model.jl` and
`example_simulate_observations_from_image.jl`.

### Polychromatic simulation

`example_simulate_polychromatic_disk.jl` demonstrates simulating a chromatic,
time-variable disk with an off-axis companion. The companion introduces
wavelength-dependent differential phases, producing non-zero OI_VIS signals.
The example writes OI_VIS2, OI_VIS (with differential phases), OI_T3, and
OI_FLUX tables.

For an image **cube**, the spectrum matters: each plane is normalised to unit total flux (as
it must be, for the visibilities to be correct), but the plane-to-plane totals are captured
first and used to weight the photon count per channel and to fill `OI_FLUX`.

## Keyword arguments to `simulate`

| keyword | default | meaning |
|---|---|---|
| `image` / `pixsize` | `""` / `0.1` | truth image (2-D) or cube (3-D), and its pixel scale in mas |
| `flat_model` / `flat_params` | `nothing` | parametric model instead of an image |
| `mag` | `2.0` | target magnitude: a number, a `Dict("H"=>1.8, …)` of band magnitudes, or one value per spectral channel |
| `mag_ao` | from `mag` | guide-star magnitude in the AO wavefront-sensor band |
| `noise` | `true` | add noise; `false` writes the model with its computed error bars |
| `debias` | `true` | subtract the `2σ²` bias from `V²`, as a real pipeline does |
| `n_samples` | `100` | Monte-Carlo samples for the T3 error bars; `0` uses the analytic form |
| `seed` / `rng` | `nothing` | make the realisation reproducible |
| `observability` | `nothing` | opt-in observability filtering, see below |
| `nonoise` | — | deprecated spelling of `noise=false` |

Coordinates follow the OIFITS standard: `target.raep0` and `target.decep0` are in **degrees**.

## Noise model

Photons per telescope, per spectral channel, per frame:

```
N = F0(λ) · 10^(-m(λ)/2.5) · A_tel · δλ · DIT
    · T_atm(λ)             atmospheric transmission (`atm_transmission`, 1.0 by default)
    · facility.throughput  telescopes and beam train only
    · T_ins(λ)             combiner (`throughput_fringes` / `throughput_photometry`)
    · flux_frac            split between interferometric and photometric channels
    · QE
    · S(λ, elevation, m_ao)   Strehl ratio
```

`S` comes from [`strehl_ratio`](@ref), a port of JMMC's `Band.strehl`, and reproduces ASPRO 2's
published CHARA Strehl curves to a median 0.7%. It needs an `[ao]` block in the facility
config; without one the code falls back to the seeing-limited coupling `min(1,(r₀/D)²)`,
which underestimates an AO-equipped array by roughly 5× in H and 20× in R.

Noise is drawn **once**, on the complex visibility, and every observable is derived from that
one perturbation. `VISAMP² == VIS2` exactly and the closure phase is exactly the sum of the
three baseline phases — which is not true if each observable is noised independently.

T3AMP and T3PHI error bars are estimated by sampling (`n_samples`, default 100). The analytic
closure errors are a small-error expansion and are only adequate while every baseline is well
detected; measured reduced chi² against the true model, on a resolved disc:

| median SNR(V²) | 6.3 | 1.0 | 0.09 | 0.02 |
|---|---|---|---|---|
| T3AMP, analytic | 1.04 | 2.16 | 27.6 | 159.6 |
| T3AMP, sampled  | 1.01 | 1.02 | 1.03 | 1.03 |
| T3PHI, analytic | 1.02 | 0.96 | 0.58 | 0.57 |
| T3PHI, sampled  | 1.01 | 1.01 | 1.00 | 1.00 |

Sampling costs no measurable time, so it is the default; set `n_samples=0` for the analytic
form.

Run `demos/validate_noise_model.jl` to print the Strehl comparison against ASPRO and the
predicted σ(V²)/σ(CP) against magnitude for MIRC-X, MYSTIC and SPICA.

## Observability filtering (opt-in)

`simulate` is a pure uv-coverage simulator: by default every epoch you pass is used, whether
or not the target was above the horizon or within the delay lines. That is usually what you
want when generating data to test image reconstruction.

To apply real constraints, either filter first:

```julia
dates_ok, mask, report = observable_epochs(facility, target, dates;
                                           min_elevation = 30.0,
                                           pops = [1,3,5,2,4,1])   # from best_pop
simulate(facility, target, combiner, wave, dates_ok, "sim.oifits"; flat_model=model)
```

or pass the same options through:

```julia
simulate(facility, target, combiner, wave, dates, "sim.oifits";
         flat_model=model, observability=(min_elevation=30.0,))
```

POP configurations are never chosen for you — omit `pops` and no delay-line check is done at
all; run [`best_pop`](@ref) yourself if you want a recommendation.

## Observation planning

OITOOLS provides tools for checking delay-line feasibility and producing
Gantt charts for a given target and night:

```julia
gantt_onenight(facility, target, date)
```

Additional planning utilities:

```julia
# Sunrise/sunset times
rise, set = sunrise_sunset(facility, date)

# Hour angle and altitude for observability
ha = jd_to_hour_angle(jd, target.raep0, facility.lon)
altitude, azimuth = alt_az(ha, target.decep0, facility.lat)

# Delay-line limits
opd_min, opd_max = opd_limits(facility, target, ha)
```

See `example_chara_plan.jl`.
