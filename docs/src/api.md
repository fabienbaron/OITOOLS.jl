# API Reference

This page lists all user-facing functions, types, and constants in OITOOLS.

## Reading OIFITS data

```@docs
OIdata
readoifits
readoifits_multiepochs
readoifits_multicolors
filter_data
set_data_filter
```

| Function | Description |
|----------|-------------|
| `list_oifits_targets(file)` | List all target names in an OIFITS file |
| `readfits(file)` | Read a FITS image into a matrix |
| `writefits(data, file)` | Write a matrix to a FITS image |
| `oifits_prep(data; kwargs...)` | Inflate error bars (additive/relative floors, multiplicative scaling) |
| `updatefits_aspro(in, out, pixsize)` | Add ASPRO-compatible WCS headers to a FITS image |
| `remove_redundant_uv!(data; uvtol)` | Merge redundant UV points in-place |

## Plotting

```@docs
uvplot
plot_v2
plot_t3phi
plot_t3amp
plot_visphi
plot_flux
imdisp
imdisp_polychromatic
```

| Function | Description |
|----------|-------------|
| `set_oiplot_defaults()` | Apply consistent matplotlib style settings |
| `plot_diffphi(data)` | Plot differential phases |
| `plot_v2_residuals(data, v2_model)` | Plot V\u00B2 residuals against a model |
| `plot_t3phi_residuals(data, t3phi_model)` | Plot closure-phase residuals |
| `plot_t3amp_residuals(data, t3amp_model)` | Plot triple-amplitude residuals |
| `plot_v2_multifile(data_vec)` | Overlay V\u00B2 from multiple datasets |
| `imdisp_temporal(cube, nepochs)` | Display a time-variable image cube |
| `plot_facility(facility)` | Plot telescope positions from a `FacilityConfig` |

## Model fitting

```@docs
FlatModel
parse_model
eval_model
eval_model_grad
display_model
fit_model
fit_model_lsqfit
fit_model_ultranest
model_to_obs
model_to_image
model_to_sed
resample_data
```

| Type | Description |
|------|-------------|
| `FitResult` | Result from `fit_model` (fields: `x_opt`, `chi2r`, `model`, ...) |
| `LsqFitResult` | Result from `fit_model_lsqfit` (adds `stderror`, `covar`, `converged`) |
| `UltraNestResult` | Result from `fit_model_ultranest` (adds `logz`, `logzerr`, `posterior`) |

## Visibility functions

Analytic visibility functions for standard source geometries.
All take baseline spatial frequency arguments and return complex visibilities.

| Function | Model |
|----------|-------|
| `visibility_ud(b, diam)` | Uniform disk |
| `visibility_ldlin(b, diam, u)` | Linear limb-darkened disk |
| `visibility_ldquad(b, diam, u, w)` | Quadratic limb-darkened disk |
| `visibility_ldquad_alt(b, diam, u, w)` | Quadratic LD disk (alternative convention) |
| `visibility_ldpow(b, diam, alpha)` | Power-law limb-darkened disk |
| `visibility_ldsquareroot(b, diam, u, w)` | Square-root limb-darkened disk |
| `visibility_annulus(b, din, dout)` | Uniform annulus |
| `visibility_ellipse_uniform(u, v, a, b, pa)` | Uniform ellipse |
| `visibility_ellipse_quad(u, v, a, b, pa, c1, c2)` | Quadratic LD ellipse |
| `visibility_thin_ring(b, diam)` | Infinitely thin ring |
| `visibility_Gaussian_ring(b, diam, fwhm)` | Gaussian ring |
| `visibility_Gaussian_ring_az(...)` | Gaussian ring with azimuthal modulation |
| `visibility_Lorentzian_ring(b, diam, fwhm)` | Lorentzian ring |
| `visibility_GaussianLorentzian_ring_az(...)` | Gaussian-Lorentzian ring with azimuthal modulation |
| `bb(T, wl)` | Planck function (black-body spectral radiance) |

## Image reconstruction

### Setup

| Function | Description |
|----------|-------------|
| `setup_dft(data, nx, pixsize)` | Build a DFT matrix for exact Fourier imaging |
| `setup_dft_polychromatic(data, nx, pixsize)` | DFT setup for polychromatic data |
| `setup_nfft(data, nx, pixsize)` | Build an NFFT plan for fast Fourier imaging |
| `setup_nfft_multiepochs(data, nx, pixsize)` | NFFT setup for multi-epoch data |
| `setup_nfft_polychromatic(data, nx, pixsize)` | NFFT setup for polychromatic data |

### Forward model

| Function | Description |
|----------|-------------|
| `image_to_vis_dft(x, dft)` | Image to complex visibilities via DFT |
| `image_to_vis(x, plan)` | Image to complex visibilities via NFFT |
| `image_to_v2(x, ft, data)` | Image to squared visibilities |
| `image_to_t3phi(x, ft, data)` | Image to closure phases |
| `image_to_t3amp(x, ft, data)` | Image to triple amplitudes |
| `image_to_obs(x, ft, data)` | Image to all observables (V², T3phi, T3amp) |
| `vis_to_v2(cvis, indx)` | Complex visibilities to squared visibilities |
| `vis_to_t3(cvis, i1, i2, i3)` | Complex visibilities to triple product, T3amp, T3phi |
| `observables(cvis, data)` | Complex visibilities to (V², T3phi, T3amp) |

### Chi-squared and criterion

| Function | Description |
|----------|-------------|
| `chi2_f(x, ft, data)` | Image chi-squared (value only) |
| `chi2_fg(x, g, ft, data)` | Image chi-squared with gradient |
| `chi2_vis_nfft_f(x, plan, data)` | Chi-squared via NFFT (value only) |
| `chi2_vis_dft_fg(x, g, dft, data)` | Chi-squared + gradient via DFT |
| `chi2_vis_nfft_fg(x, g, plan, data)` | Chi-squared + gradient via NFFT |
| `chi2_polychromatic_f(x, ft, data)` | Polychromatic chi-squared (value only) |
| `crit_f(x, ft, data)` | Criterion = chi-squared + regularization (value only) |
| `crit_fg(x, g, ft, data)` | Criterion + gradient (includes normalization correction) |
| `crit_polychromatic_fg(x, g, ft, data)` | Polychromatic criterion + gradient |
| `crit_multitemporal_fg(x, g, ft, data)` | Multi-temporal criterion + gradient |

### Reconstruction

| Function | Description |
|----------|-------------|
| `reconstruct(x0, data, ft)` | Monochromatic image reconstruction (VMLMB) |
| `reconstruct_multitemporal(x0, data, ft)` | Multi-epoch reconstruction |
| `reconstruct_polychromatic(x0, data, ft)` | Polychromatic reconstruction |
| `reconstruct_sparco_gray(x0, data, ft)` | SPARCO grey image reconstruction |
| `gaussian2d(nx, ny, sigma)` | Generate a 2D Gaussian starting image |

### Maximum entropy (BSMEM)

```@docs
reconstruct_bsmem
maxent_setup
maxent_reconstruct!
auto_pixsize
gaussian_prior
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
| `gantt_onenight(facility, target, date)` | Observation planning Gantt chart |
| `sunrise_sunset(facility, date)` | Sunrise/sunset times for a facility |
| `alt_az(ha, dec, lat)` | Hour angle to altitude/azimuth |
| `opd_limits(facility, target, ha)` | Delay-line OPD limits |
| `query_target_from_simbad(name)` | Query SIMBAD for target information |
| `ra_dec_from_simbad(name)` | Get RA/Dec from SIMBAD |

| Type | Description |
|------|-------------|
| `FacilityConfig` | Array facility configuration (telescopes, positions, atmosphere) |
| `TargetConfig` | Target configuration (coordinates, proper motion) |
| `CombinerConfig` | Beam combiner configuration (throughput, noise, calibration) |
| `WaveConfig` | Wavelength/spectral configuration |

## Writing OIFITS

| Function | Description |
|----------|-------------|
| `write_oi_header(file)` | Write OIFITS primary header |
| `write_oi_array(file, facility)` | Write OI\_ARRAY table |
| `write_oi_target(file, target)` | Write OI\_TARGET table |
| `write_oi_wavelength(file, wave)` | Write OI\_WAVELENGTH table |
| `write_oi_vis2(file, data)` | Write OI\_VIS2 table |
| `write_oi_t3(file, data)` | Write OI\_T3 table |
| `oifits_check(file)` | Validate an OIFITS file |
| `oifits_merge(files, outfile)` | Merge multiple OIFITS files |
| `oifits_filter(infile, outfile)` | Filter an OIFITS file |

## PMOIRED compatibility

```@docs
pmoired_to_julia
pmoired_to_julia_file
```
