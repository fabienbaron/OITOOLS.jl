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
facility = read_facility_file("CHARA_new")
```

| Config name | Interferometer | Telescopes |
|---|---|---|
| `CHARA` | CHARA array | 6×1 m |
| `CHARA_new` | CHARA array (updated) | 6×1 m |
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
| `MIRC` | MIRC (legacy) | CHARA | H |
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
| `MIRC_LOWH` | MIRC | H-band prism | H |
| `SPICA_LR` | SPICA | Low resolution | V |

### Simulating from an image

```julia
using Dates

# Observation times: every 15 minutes over a 5.5-hour window
dates = collect(DateTime(2024,8,13,3,0,0):Minute(15):DateTime(2024,8,13,8,30,0))

facility = read_facility_file("CHARA_new")
target   = read_obs_file("default_obs")
combiner = read_comb_file("MIRCX")
wave     = read_wave_file("MIRCX")

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
