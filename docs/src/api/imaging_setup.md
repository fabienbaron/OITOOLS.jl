# Imaging (Setup)

## Fourier transform plans

| Function | Description |
|----------|-------------|
| `setup_ft(data, nx, pixsize; mode="nfft")` | Set up Fourier transform plans (NFFT or DFT) for any `Matrix{OIdata}` |
| `setup_dft(data, nx, pixsize)` | Build a DFT matrix for a single `OIdata` |
| `setup_nfft(data, nx, pixsize)` | Build an NFFT plan for a single `OIdata` |
| `gaussian2d(nx, ny, sigma)` | Generate a 2D Gaussian starting image |

```@docs
setup_ft
gaussian2d
```

## Forward model

| Function | Description |
|----------|-------------|
| `image_to_vis(x, ft)` | Image to complex visibilities (DFT or NFFT) |
| `image_to_v2(x, ft, data)` | Image to squared visibilities |
| `image_to_t3phi(x, ft, data)` | Image to closure phases |
| `image_to_t3amp(x, ft, data)` | Image to triple amplitudes |
| `image_to_obs(x, ft, data)` | Image to all observables (V², T3phi, T3amp) |
| `image_to_residuals(x, ft, data)` | Image to normalised residuals (model - data) / error |
| `vis_to_v2(cvis, indx)` | Complex visibilities to squared visibilities |
| `vis_to_t3(cvis, i1, i2, i3)` | Complex visibilities to triple product, T3amp, T3phi |
| `observables(cvis, data)` | Complex visibilities to (V², T3phi, T3amp) |

## Chi-squared and criterion

| Function | Description |
|----------|-------------|
| `image_to_chi2(x, ft, data)` | Image chi-squared (value only). Accepts 2D or 4D images, `OIdata` or `Matrix{OIdata}` |
| `image_to_chi2_fg(x, g, ft, data)` | Image chi-squared with gradient. Same dispatch flexibility |
| `crit_f(x, ft, data)` / `image_to_crit(x, ft, data)` | Criterion (χ² + regularization)/ndof (forward only, no gradient) |
| `crit_fg(x, g, ft, data)` | Criterion (χ² + regularization)/ndof + gradient |

```@docs
image_to_residuals
image_to_chi2
image_to_chi2_fg
crit_fg(::AbstractArray{<:AbstractFloat,4}, ::AbstractArray{<:AbstractFloat,4}, ::AbstractMatrix, ::AbstractMatrix{<:OITOOLS.OIdata})
crit_f(::AbstractArray{<:AbstractFloat,4}, ::AbstractMatrix, ::AbstractMatrix{<:OITOOLS.OIdata})
```

## L-curve and regularizer scaling

| Function | Description |
|----------|-------------|
| `lcurve(x_start, data, ft, reg_name, mu_range)` | Sweep regularizer weight and return L-curve data |
| `lcurve_elbow(chi2_vals, reg_vals)` | Find L-curve corner (maximum discrete curvature) |
| `lcurve_normalize(chi2_vals, reg_vals)` | Normalize L-curve axes by elbow-point values |
| `scale_regularizers(regularizers, pixsize)` | Rescale regularizer weights for pixel-size independence |
| `scale_tv_epsilon(pixsize)` | Scale TV relaxation parameter ε for pixel-size independence |

```@docs
lcurve
lcurve_elbow
lcurve_normalize
scale_regularizers
scale_tv_epsilon
```

## Deprecated

These functions still work but emit deprecation warnings. Use the unified API instead.

| Deprecated | Replacement |
|------------|-------------|
| `setup_nfft_polychromatic` | `setup_ft` |
| `setup_dft_polychromatic` | `setup_ft(data, nx, pixsize; mode="dft")` |
| `setup_nfft_multiepochs` | `setup_ft` |
| `chi2_f` | `image_to_chi2` |
| `chi2_fg` | `image_to_chi2_fg` |
| `chi2_polychromatic_f` | `image_to_chi2` with 4D image + `Matrix{OIdata}` |
| `crit_polychromatic_fg` | `crit_fg` with 4D image + `Matrix{OIdata}` |
| `crit_multitemporal_fg` | `crit_fg` with 4D image + `Matrix{OIdata}` |
| `reconstruct_polychromatic` | `reconstruct` with 4D image + `Matrix{OIdata}` |
| `reconstruct_multitemporal` | `reconstruct` with 4D image + `Matrix{OIdata}` |
| `reconstruct_sparco_gray` | `reconstruct_hybrid` with flat model dict |
| `reconstruct_sparco_multi` | `reconstruct_hybrid` with flat model dict |
| `chi2_sparco_multi_f` | `chi2_sparco_flat_f` with flat model dict |
| `optimize_sparco_multi_parameters` | `optimize_sparco_flat_parameters` with flat model dict |
