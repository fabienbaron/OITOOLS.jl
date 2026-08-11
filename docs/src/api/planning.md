# Observation Planning

| Function | Description |
|----------|-------------|
| `night_observability(facility, ra, dec, date)` | Compute observability windows, moon separation, etc. |
| `best_pop(facility, dec, ha, config)` | Find best POP configurations for a target |
| `print_pop_results(facility, config, results)` | Pretty-print `best_pop` results |
| `obs_plan(target, facility, ra, dec, date, pop, config)` | Gantt-style observability plot |
| `chara_plan(target, facility, ra, dec, date, pop, config)` | CHARA-plan style delay plot |
| `compute_delays(facility, dec, ha, pop, config)` | Compute OPD delays for given POP/config |
| `in_delay(facility, dec, ha, pop, config)` | Check which baselines are within delay limits |
| `observable_epochs(facility, target, dates)` | Filter epochs on elevation and (optionally) delay lines — opt-in, used by `simulate` |
| `moon_illumination(jd)` | Fractional lunar illumination for a Julian Date |
| `moon_radec(jd)` | Approximate lunar RA/Dec for a Julian Date |
| `angular_separation(ra1, dec1, ra2, dec2)` | Angular separation in degrees |

```@docs
night_observability
best_pop
print_pop_results
obs_plan
chara_plan
compute_delays
in_delay
moon_illumination
moon_radec
angular_separation
```

## Simulation

| Function | Description |
|----------|-------------|
| `simulate(facility, target, combiner, wave, dates, outfile)` | Simulate observations from array geometry |
| `simulate_from_oifits(infile, outfile)` | Simulate using UV coverage of an existing OIFITS file |
| `get_uv(facility, target, combiner, wave, dates)` | Compute UV coordinates for a given setup |
| `read_facility_file(file)` | Read facility configuration (TOML) |
| `read_obs_file(file)` | Read target configuration (TOML) |
| `read_comb_file(file)` | Read combiner configuration (TOML) |
| `read_wave_file(file)` | Read wavelength configuration (TOML) |
| `observable_epochs(facility, target, dates; min_elevation, pops)` | Opt-in observability filter for `simulate` |
| `strehl_ratio(λ, D, seeing, t0, mag_ao, ao)` | AO Strehl ratio (port of JMMC `Band.strehl`) |
| `coupling_efficiency(λ, D, seeing, t0, mag_ao, ao)` | Strehl, or seeing-limited coupling when `ao === nothing` |
| `fried_parameter(seeing, λ; elevation_deg)` | Fried parameter r₀ at a wavelength and airmass |
| `seeing_from_r0(r0)` | Seeing in arcsec from a Fried parameter |
| `atm_transmission(λ)` | Atmospheric transmission hook (1.0 by default) |
| `band_for_wavelength(λ)` / `band_by_name(name)` | Photometric band lookup |
| `zero_point_flux(λ)` | Zero-magnitude photon flux, ph/s/m²/µm |
| `airmass(alt_deg)` | Plane-parallel airmass from altitude |
| `datetime_to_jd(dt)` / `datetime_to_mjd(dt)` | Julian and Modified Julian Date |
| `gantt_onenight(facility, target, date)` | Observation planning Gantt chart |
| `sunrise_sunset(facility, date)` | Sunrise/sunset times for a facility |
| `alt_az(ha, dec, lat)` | Hour angle to altitude/azimuth |
| `opd_limits(facility, target, ha)` | Delay-line OPD limits |

| Type | Description |
|------|-------------|
| `FacilityConfig` | Array facility configuration (telescopes, positions, atmosphere) |
| `TargetConfig` | Target configuration (coordinates, proper motion) |
| `CombinerConfig` | Beam combiner configuration (throughput, noise, calibration) |
| `WaveConfig` | Wavelength/spectral configuration |

```@docs
observable_epochs
AOConfig
CombinerConfig
strehl_ratio
coupling_efficiency
fried_parameter
seeing_from_r0
atm_transmission
PhotometricBand
PHOTOMETRIC_BANDS
band_for_wavelength
band_by_name
zero_point_flux
airmass
datetime_to_jd
datetime_to_mjd
```
