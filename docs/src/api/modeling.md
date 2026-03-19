# Model Fitting

## Parametric models

| Function | Description |
|----------|-------------|
| `parse_model(dict, fit_params)` | Compile a flat parameter dict into a `FlatModel` |
| `model_to_vis(model, x, uv)` | Evaluate complex visibilities for a model (alias: `eval_model`) |
| `eval_model_grad(model, x, uv)` | Evaluate visibilities + Jacobian |
| `display_model(dict, fit_params)` | Pretty-print model parameters |
| `fit_model(dict, fit_params, data)` | Fit model via gradient descent (NLopt) |
| `fit_model_lsqfit(dict, fit_params, data)` | Fit model via Levenberg-Marquardt (LsqFit) |
| `fit_model_ultranest(dict, fit_params, data)` | Fit model via nested sampling (UltraNest) |
| `model_to_obs(model, x, data)` | Compute observables (V², T3amp, T3phi) from a model |
| `model_to_chi2(model, x, data)` | Compute weighted chi² (alias: `chi2_flat`) |
| `model_to_chi2_fg(model, x, data)` | Compute chi² + gradient (alias: `chi2_flat_fg`) |
| `model_to_image(model, x; nx, pixsize)` | Synthesize a model image via inverse FFT |
| `model_to_sed(model, x, wl_grid)` | Compute spectral energy distribution |
| `resample_data(data)` | Bootstrap resample data (add Gaussian noise from error bars) |

```@docs
FlatModel
parse_model
model_to_vis
eval_model
eval_model_grad
display_model
fit_model
fit_model_lsqfit
fit_model_ultranest
model_to_obs
model_to_chi2
model_to_chi2_fg
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
