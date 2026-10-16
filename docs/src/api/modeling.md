# Model Fitting

## Fitters

| Function | Description |
|----------|-------------|
| `fit_model(model, x0, data)` | Fit a `FlatModel` via NLopt |
| `fit_model_lsqfit(model, x0, data)` | Fit a `FlatModel` via Levenberg-Marquardt |
| `fit_model_nested(model, data; lb, ub)` | Nested sampling: posterior and log-evidence; needs finite bounds |
| `fit_model_ultranest(model, data; lb, ub)` | `fit_model_nested` pinned to the UltraNest backend |
| `nested_backend()`, `set_nested_backend!(b)` | The default nested sampler: `:nautilus` or `:ultranest` |
| `fit_model_ensemble(model, free, data; lb, ub)` | Posterior samples by Ensemble Slice Sampling, with R̂ and τ |
| `fit_model_pigeons(model, free, data; lb, ub)` | Posterior samples and log-evidence by parallel tempering; needs `using Pigeons` |
| `chi2_map(model_dict, free, data, p1, p2)` | χ² surface over a grid of two free parameters |

The packages each fitter needs are listed under [What each optimiser needs](@ref).

| Type | Description |
|------|-------------|
| `FitResult` | Result from `fit_model` (fields: `x_opt`, `chi2r`, `model`, ...) |
| `LsqFitResult` | Result from `fit_model_lsqfit` (adds `stderror`, `covar`, `converged`) |
| `NestedResult` | Result from `fit_model_nested` (adds `logz`, `logzerr`, `posterior`, `result`, `backend`) |
| `UltraNestResult` | Alias for `NestedResult` |
| `BootstrapResult` | Result from `bootstrap_fit` (adds `samples`, `median`, `sigma_minus`, `sigma_plus`, `covar`) |
| `Chi2Map` | χ² surface over two free parameters from `chi2_map`; `FitResult(map)` gives its best grid point |
| `DataBlocks` | Partition of an `OIdata` into resampling blocks |
| `ModelConstraint` | One relation between model parameters, enforced during fitting |

```@docs
fit_model
fit_model_lsqfit
fit_model_nested
fit_model_ultranest
nested_backend
set_nested_backend!
chi2_map
delta_chi2_levels
FitResult
NestedResult
UltraNestResult
Chi2Map
```

## Bounds, constraints and model files

`lb`/`ub` bound each parameter. Relations between parameters (a ring's outer diameter larger
than its inner one, two fluxes summing to one) are [`ModelConstraint`](@ref)s.

| Function | Description |
|----------|-------------|
| `default_bounds(model_dict, free; data, max_size)` | Suggested `lb`/`ub` per free parameter; with `data`, angular sizes are capped at `2 λ/B_min` |
| `max_angular_scale(data)` | The largest angular scale the shortest baseline senses, in mas |
| `ModelConstraint(lhs, op, rhs; tol)` | A relation between parameters, with `op` one of `<`, `<=`, `>`, `>=`, `=` |
| `parse_constraints(specs)` | Constraints from `ModelConstraint`s, `(lhs, op, rhs[, tol])` tuples (PMOIRED's layout) or dicts |
| `check_constraints(constraints, model_dict)` | Which constraints the starting model already satisfies |
| `read_model_file(path)` | Read a TOML model file into `(; model, free, lb, ub, constraints, priors, name)` |
| `write_model_file(path, model_dict; free, lb, ub, constraints, priors)` | Write a TOML model file |

`fit_model` passes constraints to NLopt as nonlinear constraints, which hold at the optimum; an
algorithm that does not support them (including the default `:LD_LBFGS`) is wrapped in
`:AUGLAG`. `fit_model_lsqfit`, `fit_model_nested` and the samplers add a one-sided quadratic
penalty on the normalised violation, as PMOIRED's `prior` list does; a steep χ² can outweigh it.

A TOML model file holds the model dict, the free parameters, the bounds, the constraints and the
priors.

```toml
free = ["star,ud", "disk,pa"]

[model]
"star,ud"   = 6.5
"disk,fwhm" = "$star,ud * 3"

[bounds]
"star,ud" = [0.0, 20.0]

[[constraints]]
param = "disk,diamout"
op    = ">"
value = "disk,diamin"
tol   = 0.001

[[priors]]
expr   = "star,ud"
target = 6.0
sigma  = 0.5
```

```julia
m = read_model_file("binary.toml")
res = fit_model(m.model, m.free, data; m.lb, m.ub, m.constraints, m.priors)
```

`free` is a top-level key; `[model]` holds the model dict as it is, global keys included.

```@docs
default_bounds
DEFAULT_MAX_SIZE_MAS
max_angular_scale
ModelConstraint
parse_constraints
check_constraints
DEFAULT_CONSTRAINT_TOL
read_model_file
write_model_file
```

## Measurement model

`fit_model_lsqfit(…; measurement)` sets how the model is compared with spectrally resolved
data: differential phases referenced to the continuum of each spectral vector
(`visphi = :continuum`), and PMOIRED's spectral kernel (`spectral_kernel`). A spectral vector is
one table row (a baseline, triangle or telescope at one time) across the wavelengths.

```@docs
measurement_model
MeasurementModel
spectral_groups
SpectralGroups
```

## Many starts

`multistart_fit` is PMOIRED's `gridFit`: Levenberg-Marquardt from every node of a grid
(`grid = ["p" => (min, max, step)]`) and/or from uniform or normal draws (`rand`, `randn`,
`nstarts`), without the starts outside `constrain`. The fits are sorted as in PMOIRED's
`analyseGrid`: failed fits, fits without uncertainties and fits that did not move are set aside,
and the rest are grouped into distinct minima, lowest χ²ᵣ first. An explored parameter that is
fixed in the model is held at each start's value (PMOIRED's `doNotFit`), so a grid over held
parameters gives a χ² profile. The fits run on threads; start i draws from `Xoshiro(seed + i)`,
so the result does not depend on the number of threads. For a component on an orbit, minima
whose MJD0 differ by whole periods are merged.

```julia
r = multistart_fit(m, ["B,f", "B,x", "B,y"], data; weights = [0, 0, 1, 1, 0, 0, 0],
                   grid = ["B,x" => (-12, 12, 1.5), "B,y" => (-12, 12, 1.5)],
                   constrain = [("sqrt(\$B,x^2 + \$B,y^2)", "<=", 12)])
r.minima[1].fit          # the global minimum, an LsqFitResult
```

```@docs
multistart_fit
MultiStartResult
MultiStartMinimum
```

## Posterior samplers

Two samplers besides nested sampling return a posterior for the free parameters, with the
bounds as a uniform prior and −(χ² + constraint penalty)/2 as the log-likelihood: Ensemble Slice
Sampling (the `zeus` algorithm, gradient-free) and non-reversible parallel tempering with
Pigeons.jl, for multimodal posteriors. Both return a [`SamplerResult`](@ref); its `posterior` is
drawn by `plot_corner_makie` and the Modeling tab's corner view.

```@docs
fit_model_ensemble
fit_model_pigeons
SamplerResult
ensemble_slice
OITOOLS.rhat
OITOOLS.iat
```

## Uncertainty estimation by resampling

| Function | Description |
|----------|-------------|
| `bootstrap_fit(model_dict, list_free_params, data)` | Nonparametric block bootstrap: refit replicates in which blocks of data are resampled |
| `bootstrap_driver(fitfun, x_opt, list_free_params)` | The replicate loop and statistics of `bootstrap_fit`, for another fitter or resampling unit |
| `data_blocks(data; granularity)` | Partition data into resampling blocks (`:config`, `:epoch`, `:point`) |
| `resample_blocks(data, blocks; mode)` | One bootstrap replicate (`:replacement`, `:halfsample`, `:weights`; `:pmoired` is PMOIRED's scheme, √2 too small) |
| `block_counts(nblocks, mode)` | Block multiplicities drawn by a resampling scheme |
| `block_weights(nblocks)` | Continuous block weights for the multiplier (Bayesian) bootstrap |
| `apply_block_weights(data, blocks, w)` | Build the weighted replicate (error bars scaled by 1/√w) |
| `apply_block_counts(data, blocks, counts)` | Build the replicate for given block multiplicities |
| `perturb_data(data)` | Add Gaussian noise drawn from the error bars (deprecated alias `resample_data`) |

`bootstrap_fit` resamples blocks of data, by default one per (MJD, telescope configuration),
so its error bars include correlated calibration errors and wrong error bars, which the
analytic covariance of `fit_model_lsqfit` does not. See [Uncertainties](@ref) for the
calibration study.

```@docs
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
