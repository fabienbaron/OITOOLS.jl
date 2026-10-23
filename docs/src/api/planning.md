# Observation Planning

| Function | Description |
|----------|-------------|
| `night_observability(facility, ra, dec, date)` | One night's LST, hour angle, altitude, azimuth, delay and Moon constraints; `ra`, `dec` in decimal degrees |
| `best_pop(facility, dec, ha, config)` | The POP configurations with the most usable time |
| `print_pop_results(facility, config, results)` | Print `best_pop` results |
| `pop_plan(facility, dec, ha, config; observable, reference, max_changes, max_drop, change_minutes)` | A night's sequence of POP settings and dropped telescopes, ranked by time with every telescope, then with fewer, then fewest changes |
| `pop_plans(facility, dec, ha, config; …)` | The best `pop_plan` for each number of changes, from one search |
| `print_pop_plan(facility, plan)` | Print a `pop_plan` result |
| `obs_plan(target, facility, ra, dec, date, pop, config)` | Gantt chart of the night (needs PythonPlot) |
| `chara_plan(target, facility, ra, dec, date, pop, config)` | Delay plot in the style of CHARA's planner (needs PythonPlot) |
| `compute_delays(facility, dec, ha, config, pop)` | Compute OPD delays for a given config/POP. `config`: 0=unused, 1=use, 2=reference cart; `pop`: one entry per telescope in 1:5 |
| `in_delay(facility, dec, ha, config, pop)` | Check which baselines are within delay limits (same argument order as `compute_delays`) |
| `delay_ha_intervals(facility, dec, config, pop)` | The hour-angle intervals where every baseline is within its delay limits, solved analytically |
| `horizon_limit(facility, az; config)` | Lowest elevation the telescopes in use can see at that azimuth, from the facility's horizon profiles |
| `night_window(date, lat, lon; zenith)` | The dark window bracketing local midnight. `zenith` picks the twilight: 96 civil, 102 nautical, 108 astronomical |
| `observable_epochs(facility, target, dates; min_elevation, pops)` | The epochs within elevation and, with `pops`, delay limits; also through `simulate(…; observability)` |
| `index_runs(idx)` | Consecutive runs in a sorted index vector, as `(first, last)` pairs |
| `moon_illumination(jd)` | Fractional lunar illumination for a Julian Date |
| `moon_radec(jd)` | Approximate lunar RA/Dec for a Julian Date |
| `angular_separation(ra1, dec1, ra2, dec2)` | Angular separation in degrees |

```@docs
night_observability
best_pop
print_pop_results
pop_plan
pop_plans
PopPlan
PopSegment
print_pop_plan
obs_plan
chara_plan
compute_delays
in_delay
delay_ha_intervals
horizon_limit
night_window
moon_illumination
moon_radec
angular_separation
index_runs
```

## Simulation

| Function | Description |
|----------|-------------|
| `simulate(facility, target, combiner, wave, dates, outfile)` | Simulate observations from array geometry |
| `simulate_from_oifits(infile, outfile)` | Simulate using UV coverage of an existing OIFITS file |
| `wave_from_image(file_or_model)` | The spectral setup a cube defines (one channel per plane), for `simulate(…, :image, …)` |
| `get_uv(l, h, λ, δ, baselines)` | UV coordinates from latitude `l` (rad), hour angles `h` (rad, row vector), wavelengths `λ` (m), declination `δ` (rad) and a 3×n baseline matrix (m). Returns `(nuv, uv, u_M, v_M, w_M)`; `uv` is in cycles/rad |
| `read_facility_file(file)` | Read facility configuration (TOML) |
| `read_obs_file(file)` | Read target configuration (TOML) |
| `read_comb_file(file)` | Read combiner configuration (TOML) |
| `read_wave_file(file)` | Read wavelength configuration (TOML) |
| `list_configs([dir])` | Enumerate available facility / combiner / wavelength / target configs, classified by content, plus the wavelength→combiner map |
| `predict_errors(facility, combiner, wavelength; mag, visamp, elevation_deg)` | Predicted per-channel σ(V²), σ(CP), σ(visamp), σ(visphi), Strehl, DIT and photon counts, without simulating |
| `strehl_ratio(λ, D, seeing, t0, mag_ao, ao)` | AO Strehl ratio (port of JMMC `Band.strehl`) |
| `coupling_efficiency(λ, D, seeing, t0, mag_ao, ao)` | Strehl, or seeing-limited coupling when `ao === nothing` |
| `fried_parameter(seeing, λ; elevation_deg)` | Fried parameter r₀ at a wavelength and airmass |
| `seeing_from_r0(r0)` | Seeing in arcsec from a Fried parameter |
| `atm_transmission(λ)` | Atmospheric transmission hook (1.0 by default) |
| `band_for_wavelength(λ)` / `band_by_name(name)` | Photometric band lookup |
| `zero_point_flux(λ)` | Zero-magnitude photon flux, ph/s/m²/µm |
| `airmass(alt_deg)` | Plane-parallel airmass from altitude |
| `datetime_to_jd(dt)` / `datetime_to_mjd(dt)` | Julian and Modified Julian Date |
| `gantt_onenight(targetname, obsdate, lst, lst_midnight, az, alt, good_alt, good_delay; good_twilight, figsize, savefile, show_alt)` | Gantt chart from `night_observability` output, as drawn by `obs_plan` (needs PythonPlot) |
| `sunrise_sunset(obsdate, latitude, longitude; zenith=102.0)` | Rise/set times in decimal UT hours. `zenith` in degrees: 90°50′ sunrise/sunset, 96° civil, 102° nautical (default), 108° astronomical twilight |
| `alt_az(dec_deg, lat_deg, ha_hours)` | Altitude and azimuth in degrees; arguments declination, latitude (degrees), hour angle (hours) |
| `opd_limits(base, alt, az)` | Projected optical path of a baseline vector; `alt`, `az` in radians (`alt_az` returns degrees) |

| Type | Description |
|------|-------------|
| `FacilityConfig` | Array facility configuration (telescopes, positions, atmosphere) |
| `TargetConfig` | Target configuration (coordinates, proper motion) |
| `CombinerConfig` | Beam combiner configuration (throughput, noise, calibration) |
| `WaveConfig` | Wavelength/spectral configuration |

```@docs
simulate
FractionalNoise
simulate_from_oifits
wave_from_image
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
read_facility_file
read_obs_file
read_comb_file
read_wave_file
list_configs
predict_errors
OITOOLS.fringe_tracker_state
```
