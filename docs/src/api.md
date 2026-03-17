# API Reference

This page lists all user-facing functions, types, and constants in OITOOLS.
Symbols with full docstrings are shown in detail; the rest are listed with
brief descriptions.

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
| `plot_v2_residuals(data, v2_model)` | Plot V² residuals against a model |
| `plot_t3phi_residuals(data, t3phi_model)` | Plot closure-phase residuals |
| `plot_t3amp_residuals(data, t3amp_model)` | Plot triple-amplitude residuals |
| `plot_v2_multifile(data_vec)` | Overlay V² from multiple datasets |
| `imdisp_temporal(cube, nepochs)` | Display a time-variable image cube |
| `plot_facility(facility)` | Plot telescope positions from a `FacilityConfig` |

## Model definition and parsing

```@docs
FlatModel
AnalyticSpec
parse_model
eval_model
eval_model_grad
```

| Type / Function | Description |
|-----------------|-------------|
| `HankelSpec` | Component specification for Hankel-transform profiles |
| `AbstractComponentSpec` | Abstract supertype for `AnalyticSpec` and `HankelSpec` |

## Model fitting

```@docs
fit_model
fit_model_lsqfit
fit_model_ultranest
display_model
model_to_obs
model_to_image
model_to_sed
resample_data
chi2_flat
chi2_flat_fg
residuals_flat
residuals_flat_jac
```

| Type | Description |
|------|-------------|
| `FitResult` | Result struct from `fit_model` (fields: `x_opt`, `chi2r`, `model`, …) |
| `LsqFitResult` | Result struct from `fit_model_lsqfit` (adds `stderror`, `covar`, `converged`) |
| `UltraNestResult` | Result struct from `fit_model_ultranest` (adds `logz`, `logzerr`, `posterior`) |

## Visibility functions

Analytic visibility functions for standard source geometries.
All take baseline spatial frequency arguments and return complex visibilities.

| Function | Model |
|----------|-------|
| `visibility_ud(b, diam)` | Uniform disk |
| `visibility_ldlin(b, diam, u)` | Linear limb-darkened disk |
| `visibility_ldquad(b, diam, u, w)` | Quadratic limb-darkened disk |
| `visibility_ldquad_alt(b, diam, u, w)` | Quadratic LD disk (alternative convention) |
| `visibility_ldpow(b, diam, α)` | Power-law limb-darkened disk |
| `visibility_ldsquareroot(b, diam, u, w)` | Square-root limb-darkened disk |
| `visibility_annulus(b, din, dout)` | Uniform annulus |
| `visibility_ellipse_uniform(u, v, a, b, pa)` | Uniform ellipse |
| `visibility_ellipse_quad(u, v, a, b, pa, c1, c2)` | Quadratic LD ellipse |
| `visibility_thin_ring(b, diam)` | Infinitely thin ring |
| `visibility_Gaussian_ring(b, diam, fwhm)` | Gaussian ring |
| `visibility_Gaussian_ring_az(b, diam, fwhm, az_amp, az_pa)` | Gaussian ring with azimuthal modulation |
| `visibility_Lorentzian_ring(b, diam, fwhm)` | Lorentzian ring |
| `visibility_GaussianLorentzian_ring_az(...)` | Gaussian–Lorentzian ring with azimuthal modulation |
| `bb(T, λ)` | Planck function (black-body spectral radiance) |

### AD-compatible visibility functions

These use ChainRules rrules for automatic differentiation:

| Function | Model |
|----------|-------|
| `vis_ud(b, diam)` | Uniform disk (AD-compatible) |
| `vis_ldlin(b, diam, u)` | Linear LD disk (AD-compatible) |
| `vis_ldquad(b, diam, u, w)` | Quadratic LD disk (AD-compatible) |
| `vis_ldpow(b, diam, α)` | Power-law LD disk (AD-compatible) |
| `visibility_ud_d(b, diam)` | UD visibility + analytic derivative |
| `visibility_ldlin_d(b, diam, u)` | Linear LD + analytic derivative |
| `visibility_ldquad_d(b, diam, u, w)` | Quadratic LD + analytic derivative |
| `visibility_ldpow_d(b, diam, α)` | Power-law LD + analytic derivative |

## Hankel transforms

```@docs
hankel_transform
hankel_norm
hankel_vis
hankel_vis_interp
hankel_vis_full
eval_profile
compile_profile
trapz
trapz_weights
```

## Image reconstruction

### Forward model and chi-squared

| Function | Description |
|----------|-------------|
| `setup_dft(data, nx, pixsize)` | Build a DFT matrix for exact Fourier imaging |
| `setup_dft_polychromatic(data, nx, pixsize)` | DFT setup for polychromatic data |
| `setup_nfft(data, nx, pixsize)` | Build an NFFT plan for fast Fourier imaging |
| `setup_nfft_multiepochs(data, nx, pixsize)` | NFFT setup for multi-epoch data |
| `setup_nfft_polychromatic(data, nx, pixsize)` | NFFT setup for polychromatic data |
| `image_to_vis_dft(x, dft)` | Image → complex visibilities via DFT |
| `image_to_vis(x, plan)` | Image → complex visibilities via NFFT |
| `image_to_v2(x, ft, data)` | Image → squared visibilities |
| `image_to_t3phi(x, ft, data)` | Image → closure phases |
| `image_to_t3amp(x, ft, data)` | Image → triple amplitudes |
| `image_to_obs(x, ft, data)` | Image → all observables (V², T3φ, T3amp) |
| `observables(cvis, data)` | Complex visibilities → (V², T3φ, T3amp) |
| `chi2_f(x, ft, data)` | Image → χ² (value only) |
| `chi2_fg(x, g, ft, data)` | Image → χ² with gradient |
| `chi2_vis_nfft_f(x, plan, data)` | χ² via NFFT (value only) |
| `chi2_vis_dft_fg(x, g, dft, data)` | χ² + gradient via DFT |
| `chi2_vis_nfft_fg(x, g, plan, data)` | χ² + gradient via NFFT |
| `chi2_polychromatic_f(x, ft, data)` | Polychromatic χ² (value only) |
| `crit_f(x, ft, data)` | Criterion = χ² + regularization (value only) |
| `crit_fg(x, g, ft, data)` | Criterion + gradient (includes normalization correction) |
| `crit_polychromatic_fg(x, g, ft, data)` | Polychromatic criterion + gradient |
| `crit_multitemporal_fg(x, g, ft, data)` | Multi-temporal criterion + gradient |
| `reconstruct(x0, data, ft)` | Monochromatic image reconstruction (VMLMB) |
| `reconstruct_multitemporal(x0, data, ft)` | Multi-epoch reconstruction |
| `reconstruct_polychromatic(x0, data, ft)` | Polychromatic reconstruction |
| `gaussian2d(nx, ny, σ)` | Generate a 2D Gaussian starting image |
| `cdg(x)` | Centre-of-gravity of an image |
| `mod360(ϕ)` | Wrap angle to [−180, 180] |

### Regularization

| Function | Description |
|----------|-------------|
| `regularization(x, g; regularizers)` | Evaluate compound regularization (value + gradient) |
| `reg_centering(x, g)` | Soft centroid constraint |
| `tvsq(x, g)` | Total variation squared |
| `tv(x, g)` | Total variation |
| `l1l2(x, g; α)` | Edge-preserving L1–L2 smoothness |
| `l1l2w(x, g)` | Weighted L1–L2 smoothness |
| `l1hyp(x, g)` | L1 hyperbolic smoothness |
| `l2sq(x, g)` | Quadratic (L2²) smoothness |
| `compactness(x, g)` | Compactness penalty |
| `radial_variance(x, g; H, G)` | Radial variance (requires `setup_radial_reg`) |
| `entropy(x, g)` | Maximum entropy |
| `reg_support(x, g; prior)` | Support constraint (penalise flux outside mask) |
| `setup_radial_reg(nx, pixsize)` | Precompute H, G matrices for `radial_variance` |

### SPARCO

| Function | Description |
|----------|-------------|
| `chi2_sparco_f(x, ft, data)` | SPARCO χ² (value only) |
| `chi2_sparco_fg(x, g, ft, data)` | SPARCO χ² + gradient |
| `reconstruct_sparco_gray(x0, data, ft)` | SPARCO grey image reconstruction |
| `optimize_sparco_parameters(data, ft)` | Optimize SPARCO stellar parameters |
| `disk(nx, r)` | Generate a uniform-disk image of radius `r` pixels |

## Maximum entropy (BSMEM)

```@docs
reconstruct_bsmem
maxent_reconstruct!
maxent_setup
MaximENTParams
auto_pixsize
gaussian_prior
```

| Function | Description |
|----------|-------------|
| `reconstruct!(ctx, state, params)` | Low-level MaxEnt reconstruction loop |
| `maxent_step!(state, params)` | Single MaxEnt iteration step |

## Simulation

| Function | Description |
|----------|-------------|
| `simulate(facility, target, combiner, wave, dates, outfile)` | Simulate observations from array geometry |
| `simulate_from_oifits(infile, outfile)` | Simulate using UV coverage from an existing OIFITS file |
| `get_uv(facility, target, combiner, wave, dates)` | Compute UV coordinates for given setup |
| `get_uv_indxes(data)` | Get UV index arrays for V² and T3 |
| `prep_arrays(facility, combiner, dates)` | Prepare array geometry for simulation |
| `read_obs_file(file)` | Read a target configuration TOML file → `TargetConfig` |
| `read_comb_file(file)` | Read a combiner configuration TOML file → `CombinerConfig` |
| `read_wave_file(file)` | Read a wavelength configuration TOML file → `WaveConfig` |
| `read_facility_file(file)` | Read a facility configuration TOML file → `FacilityConfig` |
| `get_v2_baselines(facility, combiner)` | Get V² baseline pairs |
| `get_t3_baselines(facility, combiner)` | Get T3 baseline triplets |
| `v2mapt3(v2_baselines, t3_baselines)` | Map V² baselines to T3 legs |
| `hour_angle_calc(facility, target, dates)` | Compute hour angles |
| `hours_to_date(hours, date)` | Convert decimal hours to DateTime |
| `sunrise_sunset(facility, date)` | Sunrise/sunset times for a facility |
| `mjd_to_utdate(mjd)` | Convert MJD to UT date string |
| `dates_to_jd(dates)` | Convert DateTime array to Julian dates |
| `jd_to_hour_angle(jd, ra, lon)` | Julian date → hour angle |
| `opd_limits(facility, target, ha)` | Delay-line OPD limits |
| `alt_az(ha, dec, lat)` | Hour angle → altitude/azimuth |
| `geometric_delay(facility, target, ha)` | Geometric delay per baseline |
| `cart_delay(facility, target, ha)` | Cart/delay-line positions |
| `query_target_from_simbad(name)` | Query SIMBAD for target information |
| `ra_dec_from_simbad(name)` | Get RA/Dec from SIMBAD |
| `get_baselines(facility)` | Get all baseline vectors |
| `recenter(x)` | Recentre an image by shifting the peak to centre |
| `gantt_onenight(facility, target, date)` | Observation planning Gantt chart |

| Type | Description |
|------|-------------|
| `FacilityConfig` | Array facility configuration (telescopes, positions, atmosphere) |
| `TargetConfig` | Target configuration (coordinates, proper motion) |
| `CombinerConfig` | Beam combiner configuration (throughput, noise, calibration) |
| `WaveConfig` | Wavelength/spectral configuration |

| Function | Description |
|----------|-------------|
| `facility_info(facility)` | Print facility configuration summary |
| `obsv_info(target)` | Print target configuration summary |
| `combiner_info(combiner)` | Print combiner configuration summary |
| `wave_info(wave)` | Print wavelength configuration summary |

## Writing OIFITS

| Function | Description |
|----------|-------------|
| `write_oi_header(file)` | Write OIFITS primary header |
| `write_oi_array(file, facility)` | Write OI_ARRAY table |
| `write_oi_target(file, target)` | Write OI_TARGET table |
| `write_oi_wavelength(file, wave)` | Write OI_WAVELENGTH table |
| `write_oi_vis2(file, data)` | Write OI_VIS2 table |
| `write_oi_t3(file, data)` | Write OI_T3 table |
| `oifits_check(file)` | Validate an OIFITS file (via oifitslib) |
| `oifits_merge(files, outfile)` | Merge multiple OIFITS files |
| `oifits_filter(infile, outfile)` | Filter an OIFITS file |

## PMOIRED compatibility

```@docs
pmoired_to_julia
pmoired_to_julia_file
```

## Utilities

| Function | Description |
|----------|-------------|
| `gaussianwrapped_to_vonmises_fast(μ, σ)` | Convert wrapped Gaussian to von Mises parameters |
| `logbesselI0(x)` | Log of modified Bessel function I₀ |
| `vis_to_v2(cvis, indx)` | Complex visibilities → squared visibilities |
| `vis_to_t3(cvis, i1, i2, i3)` | Complex visibilities → triple product, T3amp, T3φ |
