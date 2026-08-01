# Model Fitting

## Parametric models

| Function | Description |
|----------|-------------|
| `dict_to_model(model_dict, list_free_params)` | Compile a flat parameter dict into a `FlatModel` |
| `model_to_vis(model, x, uv)` | Evaluate complex visibilities for a model (alias: `eval_model`) |
| `eval_model_grad(model, x, uv)` | Evaluate visibilities + Jacobian |
| `display_model(model_dict, list_free_params)` | Pretty-print model parameters |
| `fit_model(model, x0, data)` | Fit a `FlatModel` via NLopt |
| `fit_model_lsqfit(model, x0, data)` | Fit a `FlatModel` via Levenberg-Marquardt |
| `fit_model_ultranest(model, data; lb, ub)` | Fit a `FlatModel` via nested sampling (UltraNest) |
| `model_to_obs(model, x, data)` | Compute observables (V², T3amp, T3phi) from a model |
| `model_to_residuals(model, x, data)` | Compute normalised residuals (model - data) / error |
| `model_to_chi2(model, x, data)` | Compute weighted chi² (alias: `chi2_flat`) |
| `model_to_chi2_fg(model, x, data)` | Compute chi² + gradient (alias: `chi2_flat_fg`) |
| `model_to_image(model, x; nx, pixsize)` | Synthesize a model image via inverse FFT |
| `model_to_sed(model, x, wl_grid)` | Compute spectral energy distribution |
| `model_to_flux(model, x; wl)` | Total flux at zero baseline: `real(V(0,0))` |

## Uncertainty estimation by resampling

| Function | Description |
|----------|-------------|
| `bootstrap_fit(model_dict, list_free_params, data)` | Nonparametric block bootstrap: refit replicates in which blocks of data are resampled |
| `bootstrap_driver(fitfun, x_opt, list_free_params)` | The model-agnostic replicate loop and statistics behind `bootstrap_fit`, for callers with their own fitter or resampling unit |
| `data_blocks(data; granularity)` | Partition data into resampling blocks (`:config`, `:epoch`, `:point`) |
| `resample_blocks(data, blocks; mode)` | One bootstrap replicate (`:replacement`, `:halfsample`, `:weights`; `:pmoired` reproduces PMOIRED's scheme and is biased low by √2) |
| `block_counts(nblocks, mode)` | Block multiplicities drawn by a resampling scheme |
| `block_weights(nblocks)` | Continuous block weights for the multiplier (Bayesian) bootstrap |
| `apply_block_weights(data, blocks, w)` | Build the weighted replicate (error bars scaled by 1/√w) |
| `apply_block_counts(data, blocks, counts)` | Build the replicate for given block multiplicities |
| `perturb_data(data)` | Add Gaussian noise drawn from the error bars — a simulation utility, not an uncertainty estimator (was `resample_data`) |

`bootstrap_fit` resamples *which observations are used*, in blocks of (MJD,
telescope configuration), and therefore responds to correlated calibration
errors and to mis-stated error bars — neither of which the analytic covariance
of `fit_model_lsqfit` can see.  See the model-fitting examples page, and
`demos/bootstrap_validation` for the calibration test behind that statement.

```@docs
FlatModel
dict_to_model
parse_model
model_to_vis
eval_model
eval_model_grad
display_model
fit_model
fit_model_lsqfit
fit_model_ultranest
model_to_obs
model_to_residuals
model_to_chi2
model_to_chi2_fg
model_to_image
model_to_sed
model_to_flux
bootstrap_fit
bootstrap_driver
data_blocks
resample_blocks
block_counts
apply_block_counts
block_weights
apply_block_weights
perturb_data
resample_data
BootstrapResult
DataBlocks
```

| Type | Description |
|------|-------------|
| `FitResult` | Result from `fit_model` (fields: `x_opt`, `chi2r`, `model`, ...) |
| `LsqFitResult` | Result from `fit_model_lsqfit` (adds `stderror`, `covar`, `converged`) |
| `UltraNestResult` | Result from `fit_model_ultranest` (adds `logz`, `logzerr`, `posterior`) |
| `BootstrapResult` | Result from `bootstrap_fit` (adds `samples`, `median`, `sigma_minus`, `sigma_plus`, `covar`) |
| `DataBlocks` | Partition of an `OIdata` into resampling blocks |

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
