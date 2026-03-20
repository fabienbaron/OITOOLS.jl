# Model fitting

OITOOLS supports parametric model fitting using a flat-dictionary interface
compatible with [PMOIRED](https://github.com/amerand/PMOIRED).

## Workflow overview

Model fitting in OITOOLS follows three steps:

1. **Define a model dictionary** (`model_dict`) — a `Dict{String,Any}` listing
   all parameters for every component. Values are either numbers (free or fixed)
   or expression strings with `\$`-references for derived quantities.
2. **Choose which parameters to fit** (`list_free_params`) — a `Vector{String}`
   naming the subset of keys in `model_dict` that the optimizer is allowed to vary.
   Optionally define lower/upper bounds (`lb`, `ub`) for each.
3. **Fit** — calling `fit_model` (or `fit_model_lsqfit`, `fit_model_ultranest`)
   compiles the dictionary into an efficient `FlatModel` internally, then
   optimizes the free parameters against the data.

The result contains the best-fit values (`result.x_opt`) and the compiled
`FlatModel` (`result.model`), which you can pass to `model_to_obs`,
`model_to_image`, etc.

| Object | Type | Description | Key functions |
|--------|------|-------------|---------------|
| `model_dict` | `Dict{String,Any}` | Full model specification: all parameters, formulas, flags | `display_model`, `fit_model`, `fit_model_lsqfit`, `fit_model_ultranest` |
| `list_free_params` | `Vector{String}` | Names of parameters to optimize (must be numeric in `model_dict`) | same as above |
| `model` | `FlatModel` | Compiled model returned by `parse_model` or `result.model` | `model_to_obs`, `model_to_image`, `model_to_chi2`, `model_to_sed`, `eval_model` |

!!! tip
    You rarely need to call `parse_model` yourself — the fitting functions
    do it internally. Use it only when you want to evaluate a model without
    fitting (e.g. to compute observables or render an image at known parameter
    values).

## Defining a model dictionary

Parameters use flat keys of the form `"component,parameter"`:

```julia
using OITOOLS
data = readoifits("data/AlphaCenA.oifits")[1,1]

model_dict = Dict{String,Any}(
    "star,ud"  => 8.0,    # uniform disk diameter (mas)
    "star,f"   => 1.0,    # flux fraction
)
list_free_params = ["star,ud"]  # parameters to optimize

display_model(model_dict, list_free_params)
```

The component name (e.g. `"star"`, `"disk"`, `"ring"`) is arbitrary — you
choose it. The component *type* is determined automatically from which
geometry key is present.

## Component types

### Uniform disk

A single key `ud` sets the diameter in mas:

```julia
model_dict = Dict{String,Any}("star,ud" => 1.0, "star,f" => 1.0)
```

### Gaussian

A single key `fwhm` sets the full-width at half-maximum in mas:

```julia
model_dict = Dict{String,Any}("g,fwhm" => 0.5, "g,f" => 1.0)
```

### Limb-darkened disks

Three limb-darkening laws are available. The geometry key gives the
diameter in mas:

```julia
# Linear: I(μ) = 1 - u(1-μ)
model_dict = Dict{String,Any}("star,ldlin" => 1.0, "star,u" => 0.3, "star,f" => 1.0)

# Quadratic: I(μ) = 1 - u(1-μ) - w(1-μ)²
model_dict = Dict{String,Any}("star,ldquad" => 1.0, "star,u" => 0.2, "star,w" => 0.1, "star,f" => 1.0)

# Power-law: I(μ) = μ^α
model_dict = Dict{String,Any}("star,ldpow" => 1.0, "star,alpha" => 0.5, "star,f" => 1.0)
```

### Uniform ring

Defined by inner and outer diameters in mas:

```julia
model_dict = Dict{String,Any}("ring,diamin" => 0.6, "ring,diamout" => 1.0, "ring,f" => 1.0)
```

A convenience sugar is available using `diam` (outer diameter) and `thick`
(fractional thickness 0–1):

```julia
# Equivalent to diamin = 2.0*(1-0.3) = 1.4, diamout = 2.0
model_dict = Dict{String,Any}("ring,diam" => 2.0, "ring,thick" => 0.3, "ring,f" => 1.0)
```

### Gaussian ring

A bi-Gaussian ring defined by inner and outer FWHM in mas:

```julia
model_dict = Dict{String,Any}("gr,fwhmin" => 0.2, "gr,fwhmout" => 0.5, "gr,f" => 1.0)
```

### Crescent

Two offset uniform disks producing a crescent shape:

```julia
model_dict = Dict{String,Any}(
    "cr,crin"      => 0.8,    # inner disk diameter (mas)
    "cr,crout"     => 1.0,    # outer disk diameter (mas)
    "cr,croff"     => 0.8,    # offset factor (0 = concentric ring, 1 = max offset)
    "cr,crprojang" => 120.0,  # PA of thinnest part (degrees)
    "cr,f"         => 1.0,
)
```

### Point source

A component with only `f` (no geometry key) is an unresolved point source:

```julia
model_dict = Dict{String,Any}("companion,f" => 0.05, "companion,x" => 2.0, "companion,y" => -1.0)
```

### Fully resolved background

Sets V = 0 at all non-zero baselines (fully resolved flux):

```julia
model_dict = Dict{String,Any}("bg,resolved" => true, "bg,f" => 0.1)
```

### Summary table

| Key(s) | Type | Description |
|--------|------|-------------|
| `ud` | Uniform disk | Diameter in mas |
| `fwhm` | Gaussian | FWHM in mas |
| `ldlin` + `u` | Linear limb-darkened disk | |
| `ldquad` + `u`, `w` | Quadratic limb-darkened disk | |
| `ldpow` + `alpha` | Power-law limb-darkened disk | |
| `diamin` + `diamout` | Uniform ring | Inner/outer diameters in mas |
| `diam` + `thick` | Uniform ring (sugar) | Outer diameter + fractional thickness |
| `fwhmin` + `fwhmout` | Gaussian ring | Inner/outer FWHM in mas |
| `crin`, `crout`, `croff`, `crprojang` | Crescent | Two offset disks |
| `profile` + `diamout` | Hankel profile | Arbitrary radial profile (see below) |
| *(only `f`)* | Point source | Unresolved |
| `resolved` | Resolved background | Fully resolved (V=0) |

## Common parameters

Every component supports these optional keys:

| Key | Description |
|-----|-------------|
| `f` | Flux fraction (can be a number or an expression string) |
| `x`, `y` | Position offset in mas (east, north) |
| `incl` | Inclination in degrees (0 = face-on) |
| `pa` or `projang` | Position angle in degrees (north through east) |
| `spectrum` | Spectral law string; promoted to `f` internally (e.g. `"(\$WL/2.2e-6)^-4"`) |
| `spatial_kernel` | Gaussian smoothing FWHM in visibility space (global, not per-component) |

### Geometric transformations

Inclination and position angle project the component on sky:

```julia
# Inclined uniform disk
model_dict = Dict{String,Any}(
    "star,ud"      => 1.0,
    "star,incl"    => 60.0,    # degrees from face-on
    "star,projang" => 30.0,    # PA in degrees
    "star,f"       => 1.0,
)

# Offset companion
model_dict = Dict{String,Any}(
    "star,ud" => 1.0, "star,f" => 0.9,
    "comp,f"  => 0.1, "comp,x" => 3.0, "comp,y" => -1.5,
)
```

## Expression strings

Parameter values can be strings referencing other parameters with `\$`.
References are resolved in topological order, so forward and backward
references both work:

```julia
model_dict = Dict{String,Any}(
    "star,ud"      => 3.0,
    "star,f"       => 0.7,
    "disk,f"       => "1 - \$star,f",            # complement of star flux
    "disk,diamout" => "\$star,ud * 8",           # 8× the stellar diameter
    "disk,diamin"  => "\$disk,diamout * 0.5",    # half the outer diameter
)
```

### Implicit variables

Three implicit variables are available in expressions:

| Variable | Description |
|----------|-------------|
| `\$WL` | Wavelength in metres (set per data point during fitting) |
| `\$MJD` | Modified Julian date (set per data point during fitting) |
| `\$B` | Baseline length in metres |

These enable chromatic and time-variable models:

```julia
model_dict = Dict{String,Any}(
    "star,f"        => 1.0,
    "star,spectrum"  => "(\$WL/2.2e-6)^(-4)",    # Rayleigh–Jeans spectrum
    "disk,f"        => 0.3,
    "disk,spectrum"  => "(\$WL/2.2e-6)^(-2)",    # greyer spectrum
    "disk,diamout"  => 10.0,
    "disk,profile"  => "exp(-(\$R/3.0)^2)",
)
```

### Shared parameters across components

Bare keys (without a `component,` prefix) act as global variables that
multiple components can reference:

```julia
model_dict = Dict{String,Any}(
    "PA"              => 60.0,        # global position angle
    "INC"             => 45.0,        # global inclination
    "inner,fwhm"      => 1.0,
    "inner,projang"   => "\$PA",
    "inner,incl"      => "\$INC",
    "outer,diamin"    => "2 * \$inner,fwhm",
    "outer,diamout"   => "3 * \$outer,diamin",
    "outer,projang"   => "\$PA",
    "outer,incl"      => "\$INC",
    "outer,f"         => 0.5,
)
```

## Combining multiple components

Models with multiple components are built by giving each component a
different name prefix. The total model visibility is the flux-weighted
sum of all component visibilities:

```julia
model_dict = Dict{String,Any}(
    # Compact star
    "star,fwhm"      => 0.1,
    "star,spectrum"   => "\$WL^(-3)",
    # Inclined disk with profile and azimuthal modulation
    "disk,diamin"    => 0.5,
    "disk,diamout"   => 1.0,
    "disk,profile"   => "\$R^(-2)",
    "disk,az amp1"   => 1.0,
    "disk,az projang1" => 60.0,
    "disk,projang"   => 45.0,
    "disk,incl"      => -30.0,
    "disk,x"         => -0.05,
    "disk,y"         => 0.05,
    "disk,spectrum"   => "5 * \$WL^(-2)",
)
```

Visualise the model image and SED:

```julia
params = parse_model(model_dict, String[])
img = model_to_image(params, Float64[]; nx=128, pixsize=0.02, wl=1.65e-6)
imdisp(img; pixsize=0.02)

wl_grid = collect(range(1.0e-6, 2.5e-6; length=200))
f_total, f_comps = model_to_sed(params, Float64[], wl_grid)
```

See `example_model_fitting_pmoired_models.jl` for a full gallery of all
component types and combinations.

## Hankel models (radial profiles)

For components that cannot be described by simple analytic visibility functions,
OITOOLS uses a numerical Hankel transform of an arbitrary radial profile.
This is ideal for circumstellar disks, envelopes, and limb-darkened models
with custom intensity laws.

### Profile expressions

Define a profile expression string using `\$R` (radius in mas) and `\$MU`
(cosine of the zenith angle), plus any scalar parameters:

```julia
# Disk with intensity ∝ μ^0.5 (limb-darkening)
model_dict = Dict{String,Any}(
    "star,diam"    => 1.0,
    "star,profile" => "\$MU^0.5",
    "star,f"       => 1.0,
)

# Ring with power-law profile
model_dict = Dict{String,Any}(
    "ring,udout"   => 1.0,
    "ring,profile" => "(\$R > 0.25) * \$R^(-0.5)",
    "ring,f"       => 1.0,
)

# Parabolic ring profile (YSO disk)
model_dict = Dict{String,Any}(
    "disk,profile"  => "max(0, 1 - (2*(\$R - \$Rmid)/\$width)^2)",
    "disk,diamout"  => 10.0,     # outer diameter (mas) — sets the r-grid extent
    "disk,nr"       => 200,      # radial grid points (default: 100)
    "disk,Rmid"     => 3.0,      # peak radius (custom parameter)
    "disk,width"    => 2.0,      # ring width (custom parameter)
    "disk,f"        => 0.5,
    "disk,incl"     => 45.0,
    "disk,pa"       => 120.0,
    "star,f"        => "1 - \$disk,f",
)
```

### Hankel-specific parameters

| Key | Description |
|-----|-------------|
| `profile` | Radial intensity expression (uses `\$R`, `\$MU`, and any custom params) |
| `diamout` or `udout` | Outer diameter in mas (sets r-grid extent) |
| `diamin` | Inner diameter in mas (optional, sets r-grid inner boundary for rings) |
| `nr` | Number of radial grid points (default 100; increase for sharp features) |
| `r_max` | Alternative to `diamout` for setting the outer boundary |

### Azimuthal modulations

Azimuthal variations are added with harmonic coefficients `az ampN` and
`az projangN` for harmonic order N = 1, 2, 3, ...:

```julia
model_dict = Dict{String,Any}(
    "disk,diamin"       => 1.0,
    "disk,diamout"      => 3.0,
    "disk,profile"      => "1",
    "disk,projang"      => -20.0,
    "disk,incl"         => 60.0,
    "disk,f"            => 1.0,
    "disk,az amp1"      => 0.3,     # m=1 harmonic amplitude
    "disk,az projang1"  => 45.0,    # m=1 orientation (degrees)
    "disk,az amp2"      => 0.1,     # m=2 harmonic amplitude
    "disk,az projang2"  => 90.0,    # m=2 orientation (degrees)
)
```

Multiple harmonics combine additively to produce asymmetric features
such as one-armed spirals or brightness asymmetries in disks.

### Spiral patterns

By building N concentric rings whose `az projang1` rotates with radius,
you can create spiral patterns. Each ring is a thin annulus with an m=1
azimuthal modulation; the progressive phase offset produces a spiral arm.
Add a `spatial_kernel` to smooth the result.

See `example_model_fitting_pmoired_models.jl` (section 5b) for a complete
spiral example with 5 rings using shared global parameters (`Din`, `Dout`,
`INCL`, `PROJANG`, `pitch`, `PAin`).

### YSO disk fitting example

`example_model_fitting_v1295aql.jl` fits three YSO disk models
(parabolic ring, double sigmoid, α-sigmoid) to V1295 Aql (HD 190073)
H-band data, reproducing Ibrahim et al. 2023 (ApJ, 947, 68).
It demonstrates chromatic models with wavelength-dependent spectra.

## Spatial smoothing

The global parameter `spatial_kernel` applies Gaussian smoothing (in
visibility space) to the model. This is useful for softening sharp
edges:

```julia
model_dict = Dict{String,Any}(
    "cr,crin" => 0.5, "cr,crout" => 1.0,
    "cr,croff" => 0.8, "cr,crprojang" => 120.0,
    "cr,incl" => 45.0, "cr,projang" => 30.0,
    "cr,f" => 1.0,
    "spatial_kernel" => 0.2,    # Gaussian FWHM in mas
)
```

## PMOIRED compatibility

The parameter dictionary format is designed to be directly compatible with
[PMOIRED](https://github.com/amerand/PMOIRED). Models written for PMOIRED
can be ported to OITOOLS with minimal changes.

### Converting Python dicts to Julia

`pmoired_to_dict()` converts a PMOIRED Python dict literal string directly
to a Julia `Dict`:

```julia
model_dict = pmoired_to_dict("{'star,ud': 3.2, 'ring,f': '1 - \$star,f'}")
```

The lower-level `pmoired_to_julia()` returns the Julia source string instead,
if you need to inspect or edit it before evaluation.

It handles:

| Python | Julia |
|--------|-------|
| `{ ... }` | `Dict( ... )` |
| `'key': value` | `"key" => value` |
| `'expr with \$ref'` | `raw"expr with \$ref"` |
| `True` / `False` / `None` | `true` / `false` / `nothing` |
| `[...]` (lists) | `[...]` (arrays) |
| `# comment` | `# comment` |

### Converting entire scripts

`pmoired_to_julia_file()` converts a Python/PMOIRED script file
line-by-line:

```julia
pmoired_to_julia_file("pmoired_model.py", "julia_model.jl")
```

### Full conversion example

Given a PMOIRED model:

```python
# PMOIRED (Python)
param = {
    'star,ud':    3.2,
    'star,f':     0.7,
    'ring,udout': '\$star,ud * 8',
    'ring,f':     '1 - \$star,f',
    'ring,incl':  30.0,
}
fitOnly = ['star,ud', 'star,f', 'ring,udout', 'ring,f']
```

Convert and use in Julia:

```julia
using OITOOLS

# Option 1: inline conversion
model_dict = pmoired_to_dict("{'star,ud': 3.2, 'star,f': 0.7, 'ring,udout': '\$star,ud * 8', 'ring,f': '1 - \$star,f', 'ring,incl': 30.0}")

# Option 2: file conversion
pmoired_to_julia_file("pmoired_model.py", "julia_model.jl")
include("julia_model.jl")

# Option 3: write it directly in Julia (recommended)
model_dict = Dict{String,Any}(
    "star,ud"    => 3.2,
    "star,f"     => 0.7,
    "ring,udout" => raw"$star,ud * 8",
    "ring,f"     => raw"1 - $star,f",
    "ring,incl"  => 30.0,
)
list_free_params = ["star,ud", "star,f", "ring,incl"]
```

!!! note
    Expression strings containing `$` must use `raw"..."` or `\$` in Julia
    to prevent string interpolation. `pmoired_to_julia()` handles this
    automatically.

See `example_model_fitting_pmoired_conversion.jl` for more conversion
examples, and `example_model_fitting_pmoired_models.jl` for a gallery
of all PMOIRED model types reproduced in OITOOLS.

## Fitting with NLopt (gradient-based)

`fit_model` uses NLopt for optimisation. Gradient-based methods (default
`:LD_LBFGS`) use the analytic Wirtinger-derivative chain rule; gradient-free
methods like `:LN_NELDERMEAD` are available for non-smooth profiles.

```julia
result = fit_model(model_dict, list_free_params, data;
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

result = fit_model(model_dict, list_free_params, data; lb, ub, priors)
```

## Fitting with LsqFit (Levenberg-Marquardt)

`fit_model_lsqfit` uses the Levenberg-Marquardt algorithm from LsqFit.jl.
It provides **parameter covariance** and **1σ error bars** from the Jacobian:

```julia
result = fit_model_lsqfit(model_dict, list_free_params, data;
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
result = fit_model_ultranest(model_dict, list_free_params, data;
    lb = Dict("star,ud" => 0.1),
    ub = Dict("star,ud" => 20.0),    # bounds required for all list_free_params
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
imdisp(img; pixsize=0.1)
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
results = [fit_model(model_dict, list_free_params, resample_data(data);
    lb=lb, ub=ub, maxeval=200).x_opt for _ in 1:nboot]
```

See `example_bootstrap_fit.jl`.
