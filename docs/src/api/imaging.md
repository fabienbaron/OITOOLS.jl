# Imaging

## Setup

| Function | Description |
|----------|-------------|
| `setup_dft(data, nx, pixsize)` | Build a DFT matrix for exact Fourier imaging |
| `setup_dft_polychromatic(data, nx, pixsize)` | DFT setup for polychromatic data |
| `setup_nfft(data, nx, pixsize)` | Build an NFFT plan for fast Fourier imaging |
| `setup_nfft_multiepochs(data, nx, pixsize)` | NFFT setup for multi-epoch data |
| `setup_nfft_polychromatic(data, nx, pixsize)` | NFFT setup for polychromatic data |
| `gaussian2d(nx, ny, sigma)` | Generate a 2D Gaussian starting image |

## Forward model

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

## Chi-squared and criterion

| Function | Description |
|----------|-------------|
| `chi2_f(x, ft, data)` | Image chi-squared (value only); `ft` can be a DFT matrix or NFFT plan |
| `chi2_fg(x, g, ft, data)` | Image chi-squared with gradient; dispatches on `ft` type |
| `chi2_polychromatic_f(x, ft, data)` | Polychromatic chi-squared (value only) |
| `crit_f(x, ft, data)` | Criterion = chi-squared + regularization (value only) |
| `crit_fg(x, g, ft, data)` | Criterion + gradient (includes normalization correction) |
| `crit_polychromatic_fg(x, g, ft, data)` | Polychromatic criterion + gradient |
| `crit_multitemporal_fg(x, g, ft, data)` | Multi-temporal criterion + gradient |

## Reconstruction

| Function | Description |
|----------|-------------|
| `reconstruct(x0, data, ft)` | Monochromatic image reconstruction (VMLMB) |
| `reconstruct_multitemporal(x0, data, ft)` | Multi-epoch reconstruction |
| `reconstruct_polychromatic(x0, data, ft)` | Polychromatic reconstruction |

## SPARCO

| Function | Description |
|----------|-------------|
| `reconstruct_sparco_gray(x0, params, data, ft)` | SPARCO grey image reconstruction |
| `reconstruct_sparco_multi(x0, params, nsources, data, ft)` | Multi-source SPARCO reconstruction |
| `reconstruct_sparco_flat(x0, params, model, data, ft)` | SPARCO reconstruction using flat model |
| `chi2_sparco_multi_f(x, params, nsources, ft, data)` | Multi-source SPARCO chi-squared (forward only) |
| `chi2_sparco_flat_f(x, params, model, ft, data)` | Flat-model SPARCO chi-squared (forward only) |
| `optimize_sparco_multi_parameters(params, nsources, x, ft, data)` | Optimize multi-source SPARCO parameters (image fixed) |
| `optimize_sparco_flat_parameters(params, model, x, ft, data)` | Optimize flat-model SPARCO parameters (image fixed) |

```@docs
reconstruct_sparco_gray
reconstruct_sparco_multi
reconstruct_sparco_flat
chi2_sparco_multi_f
chi2_sparco_flat_f
optimize_sparco_multi_parameters
optimize_sparco_flat_parameters
```

## Maximum entropy (BSMEM)

| Function | Description |
|----------|-------------|
| `reconstruct_bsmem(x0, data, ft)` | BSMEM image reconstruction |
| `auto_pixsize(data)` | Estimate pixel size from data |
| `gaussian_prior(nx, fwhm)` | Generate a Gaussian prior image |

```@docs
reconstruct_bsmem
auto_pixsize
gaussian_prior
```
