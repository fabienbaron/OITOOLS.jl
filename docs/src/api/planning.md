# Observation Planning

| Function | Description |
|----------|-------------|
| `night_observability(facility, ra, dec, date)` | Observability windows, Moon separation, etc. `ra` and `dec` in **decimal degrees** |
| `best_pop(facility, dec, ha, config)` | Find best POP configurations for a target |
| `print_pop_results(facility, config, results)` | Pretty-print `best_pop` results |
| `obs_plan(target, facility, ra, dec, date, pop, config)` | Gantt-style observability plot |
| `chara_plan(target, facility, ra, dec, date, pop, config)` | CHARA-plan style delay plot |
| `compute_delays(facility, dec, ha, config, pop)` | Compute OPD delays for a given config/POP. `config`: 0=unused, 1=use, 2=reference cart; `pop`: one entry per telescope in 1:5 |
| `in_delay(facility, dec, ha, config, pop)` | Check which baselines are within delay limits (same argument order as `compute_delays`) |
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
| `get_uv(l, h, λ, δ, baselines)` | UV coordinates from latitude `l` (rad), hour angles `h` (rad, row vector), wavelengths `λ` (m), declination `δ` (rad) and a 3×n baseline matrix (m). Returns `(nuv, uv, u_M, v_M, w_M)`; `uv` is in cycles/rad |
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
| `gantt_onenight(targetname, obsdate, lst, lst_midnight, az, alt, good_alt, good_delay; good_twilight, figsize, savefile, show_alt)` | Gantt chart from precomputed `night_observability` output. Usually reached through `obs_plan` |
| `sunrise_sunset(obsdate, latitude, longitude; zenith=102.0)` | Rise/set times in decimal UT hours. `zenith` in degrees: 90°50′ sunrise/sunset, 96° civil, 102° nautical (default), 108° astronomical twilight |
| `alt_az(dec_deg, lat_deg, ha_hours)` | Altitude/azimuth in degrees. Note the argument order: declination first, then latitude, then hour angle **in hours** |
| `opd_limits(base, alt, az)` | Projected optical path for a baseline vector. `alt`/`az` in **radians** — note `alt_az` returns degrees, so they do not compose directly |

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
