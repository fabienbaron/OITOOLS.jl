# Imaging

## Setup

| Function | Description |
|----------|-------------|
| `setup_ft(data, nx, pixsize; mode="nfft")` | Set up Fourier transform plans (NFFT or DFT) for any `Matrix{OIdata}` |
| `setup_dft(data, nx, pixsize)` | Build a DFT matrix for a single `OIdata` |
| `setup_nfft(data, nx, pixsize)` | Build an NFFT plan for a single `OIdata` |
| `gaussian2d(nx, ny, sigma)` | Generate a 2D Gaussian starting image |

```@docs
setup_ft
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

## Reconstruction

| Function | Description |
|----------|-------------|
| `reconstruct(x0_4d, data, ft)` | Unified image reconstruction (mono/poly/temporal, VMLMB) |
| `reconstruct(x0_2d, data, ft)` | Monochromatic convenience wrapper |

```@docs
reconstruct(::AbstractArray{<:AbstractFloat,4}, ::AbstractMatrix{<:OITOOLS.OIdata}, ::AbstractMatrix)
reconstruct(::AbstractMatrix{<:AbstractFloat}, ::OITOOLS.OIdata, ::Any)
```

## ADMM reconstruction

| Function | Description |
|----------|-------------|
| `reconstruct_admm(x0, d, nfft_plan, ft, data, nx)` | ADMM image reconstruction with proximal operators |
| `prox_v2(z0, v2_data, v2_err, indx_v2, rho)` | Proximal operator for squared visibilities |
| `prox_t3phi(z0, t3phi_deg, t3phi_err_deg, ...)` | Proximal operator for closure phases |
| `prox_tv(y, lambda)` | Proximal operator for isotropic total variation |
| `prox_l2smooth(y, lambda)` | Proximal operator for L2 smoothness |
| `tv_norm(x_img)` | Isotropic total variation of a 2D image |
| `l2smooth_norm(x_img)` | L2 smoothness norm of a 2D image |
| `chi2_v2_at(z, v2_data, v2_err, indx_v2)` | V2 chi-squared at a visibility vector |
| `chi2_t3phi_at(z, t3phi_data, t3phi_err, ...)` | Closure phase chi-squared at a visibility vector |

```@docs
reconstruct_admm
prox_v2
prox_t3phi
prox_tv
prox_l2smooth
tv_norm
l2smooth_norm
chi2_v2_at
chi2_t3phi_at
```

## SPARCO

SPARCO functions work with a single `OIdata` — use `data[1]` and `ft[1]`.

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

BSMEM works with a single `OIdata` — use `data[1]` and `ft[1]`.

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
