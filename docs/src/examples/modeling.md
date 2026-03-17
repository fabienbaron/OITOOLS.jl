# Model fitting

OITOOLS supports parametric model fitting using a flat-dictionary interface
compatible with [PMOIRED](https://github.com/amerand/PMOIRED).
Models are assembled from named components, each described by key-value pairs
in a `Dict{String,Any}`. Expression strings with `$`-references let parameters
depend on each other (e.g. `"1 - \$star,f"`).

## Defining a model

Parameters use flat keys of the form `"component,parameter"`:

```julia
using OITOOLS
data = readoifits("data/AlphaCenA.oifits")[1,1]

params = Dict(
    "star,ud"  => 8.0,    # uniform disk diameter (mas)
    "star,f"   => 1.0,    # flux fraction
)
fit_params = ["star,ud"]  # parameters to optimize

display_model(params, fit_params)
```

### Component types

The component type is inferred from which geometry key is present:

| Key(s) | Type | Description |
|--------|------|-------------|
| `ud` | Uniform disk | Diameter in mas |
| `fwhm` | Gaussian | FWHM in mas |
| `ldlin` + `u` | Linear limb-darkened disk | |
| `ldquad` + `u`, `w` | Quadratic limb-darkened disk | |
| `ldpow` + `alpha` | Power-law limb-darkened disk | |
| `diamin`, `diamout` | Uniform ring | Inner/outer diameters |
| `crin`, `crout`, `croff`, `crprojang` | Crescent | Two offset disks |
| `profile` + `diamout` | Hankel profile | Arbitrary radial profile via Hankel transform |
| *(only `f`)* | Point source | Unresolved |
| `resolved` | Resolved background | Fully resolved (V=0) |

### Common parameters (all components)

| Key | Description |
|-----|-------------|
| `f` | Flux fraction (can be an expression) |
| `x`, `y` | Position offset in mas |
| `incl` | Inclination (degrees) |
| `pa` or `projang` | Position angle (degrees) |
| `spectrum` | Spectral law string (e.g. `"($WL/2.2e-6)^-4"`) |

### Expression strings

Values can be strings referencing other parameters with `$`:

```julia
params = Dict(
    "star,f"     => 0.7,
    "disk,f"     => "\$star,f * 0.2",       # 20% of star flux
    "disk,diamout" => "\$star,ud * 8",      # 8× the stellar diameter
)
```

Available implicit variables in expressions: `\$WL` (wavelength in m),
`\$MJD` (modified Julian date), `\$B` (baseline in m).

## PMOIRED compatibility

The parameter dictionary format is designed to be compatible with PMOIRED.
Use `pmoired_to_julia()` to convert Python dict literals directly:

```julia
julia_str = pmoired_to_julia("{'star,ud': 3.2, 'ring,f': '1 - \$star,f'}")
params = eval(Meta.parse(julia_str))
```

See `example_model_fitting_pmoired_models.jl` for a gallery of PMOIRED model
types (rings, crescents, Gaussians, azimuthal modulations) reproduced in OITOOLS,
and `example_model_fitting_pmoired_conversion.jl` for conversion examples.

## Hankel models (radial profiles)

For components that cannot be described by simple analytic visibility functions,
OITOOLS uses a numerical Hankel transform of an arbitrary radial profile.
This is ideal for circumstellar disks, envelopes, and limb-darkened models
with custom intensity laws.

Define a profile expression string using `\$R` (radius in mas) and `\$MU`
(cosine of the zenith angle), plus any scalar parameters:

```julia
# Parabolic ring profile (YSO disk)
params = Dict(
    "disk,profile"  => "max(0, 1 - (2*(\$R - \$Rmid)/\$width)^2)",
    "disk,diamout"  => 10.0,     # outer radius (mas) — sets the r-grid
    "disk,nr"       => 200,      # radial grid points (default: 100)
    "disk,Rmid"     => 3.0,      # peak radius
    "disk,width"    => 2.0,      # ring width
    "disk,f"        => 0.5,
    "disk,incl"     => 45.0,     # inclination
    "disk,pa"       => 120.0,    # position angle
    "star,f"        => "1 - \$disk,f",
)
```

Azimuthal modulations can be added with harmonic coefficients:

```julia
params["disk,az amp1"]     = 0.3    # first azimuthal harmonic amplitude
params["disk,az projang1"] = 45.0   # harmonic orientation (degrees)
```

### Demo examples

- **`example_model_fitting_v1295aql.jl`** — Fits three YSO disk models
  (parabolic ring, double sigmoid, α-sigmoid) to V1295 Aql (HD 190073)
  H-band data, reproducing Ibrahim et al. 2023 (ApJ, 947, 68).
  Demonstrates chromatic models with wavelength-dependent spectra.

- **`example_model_fitting_pmoired_models.jl`** — Gallery of all supported
  component types at various inclinations with azimuthal variations.

## Fitting with NLopt (gradient-based)

`fit_model` uses NLopt for optimisation. Gradient-based methods (default
`:LD_LBFGS`) use the analytic Wirtinger-derivative chain rule; gradient-free
methods like `:LN_NELDERMEAD` are available for non-smooth profiles.

```julia
result = fit_model(params, fit_params, data;
    weights = [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],  # V² + T3φ
    method  = :LD_LBFGS,
    maxeval = 500,
    verb    = true)

println("Best χ²/ν = ", result.chi2r)
println("Parameters: ", result.x_opt)
```

The 7-element `weights` vector controls per-observable weighting:
`[V², T3amp, T3φ, visamp, visφ, flux, diffφ]`.

### Bounds and priors

```julia
lb = Dict("star,ud" => 0.1)   # lower bounds
ub = Dict("star,ud" => 20.0)  # upper bounds

# Gaussian priors: Dict of (mean, sigma)
priors = Dict("star,ud" => (8.5, 0.5))

result = fit_model(params, fit_params, data; lb, ub, priors)
```

## Fitting with LsqFit (Levenberg-Marquardt)

`fit_model_lsqfit` uses the Levenberg-Marquardt algorithm from LsqFit.jl.
It provides **parameter covariance** and **1σ error bars** from the Jacobian:

```julia
result = fit_model_lsqfit(params, fit_params, data;
    weights = [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    lb = Dict("star,ud" => 0.1),
    ub = Dict("star,ud" => 20.0),
    maxIter = 200)

println("Best χ²/ν = ", result.chi2r)
println("Parameters: ", result.x_opt)
println("1σ errors:  ", result.stderror)
println("Covariance:\n", result.covar)
println("Converged:  ", result.converged)
```

The analytic Jacobian is computed via the Wirtinger chain rule through the
full observable pipeline (V² → T3 → residuals), so no finite differences
are needed.

## Bayesian inference (UltraNest)

`fit_model_ultranest` performs nested sampling via
[UltraNest](https://johannesbuchner.github.io/UltraNest/) (Python, called
through PyCall). It returns the Bayesian log-evidence for model comparison:

```julia
result = fit_model_ultranest(params, fit_params, data;
    lb = Dict("star,ud" => 0.1),
    ub = Dict("star,ud" => 20.0),    # bounds required for all fit params
    min_num_live_points = 100,
    cornerplot = true)                # produce corner plot

println("log(Z) = ", result.logz, " ± ", result.logzerr)
println("Best χ²/ν = ", result.chi2r)
```

The posterior samples are available as `result.posterior` (matrix of
`nsamples × nparams`).

## Inspecting results

### Model observables

```julia
obs = model_to_obs(result.model, result.x_opt, data)
# obs.v2, obs.t3amp, obs.t3phi, obs.visamp, obs.visphi
```

### Synthetic images

```julia
img = model_to_image(result.model, result.x_opt;
    nx=256, pixsize=0.1, wl=1.65e-6)
imdisp(img)
```

### Spectral energy distribution

```julia
wl_grid = range(1.5e-6, 2.5e-6, length=100)
total_flux, component_fluxes = model_to_sed(result.model, result.x_opt, wl_grid)
```

## Bootstrap error estimation

Generate resampled datasets with Gaussian noise from the error bars and
re-fit to estimate parameter uncertainties:

```julia
nboot = 200
results = [fit_model(params, fit_params, resample_data(data);
    lb=lb, ub=ub, maxeval=200).x_opt for _ in 1:nboot]
```

See `example_bootstrap_fit.jl`.
