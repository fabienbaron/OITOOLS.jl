# Model fitting

Models are flat dictionaries in [PMOIRED](https://github.com/amerand/PMOIRED)'s format.

## Workflow

1. **Model dictionary** (`model_dict`): a `Dict{String,Any}` of every component's parameters.
   Values are numbers or expression strings with `\$` references.
2. **Free parameters** (`list_free_params`): the keys of `model_dict` the fit varies, with
   optional bounds `lb`, `ub`.
3. **Fit**: `fit_model` (or `fit_model_lsqfit`, `fit_model_nested`, …) compiles the dictionary
   into a `FlatModel` and fits the free parameters.

The result holds the best-fit values (`result.x_opt`) and the compiled model (`result.model`),
for `model_to_obs`, `model_to_image`, etc.

| Object | Type | Description | Key functions |
|--------|------|-------------|---------------|
| `model_dict` | `Dict{String,Any}` | Full model specification: all parameters, formulas, flags | `display_model` |
| `list_free_params` | `Vector{String}` | Names of parameters to optimize (must be numeric in `model_dict`) | `display_model` |
| `model` | `FlatModel` | Compiled model from `dict_to_model(model_dict, list_free_params)` or `result.model` | `fit_model`, `fit_model_lsqfit`, `fit_model_nested`, `model_to_obs`, `model_to_image`, `model_to_chi2`, `model_to_sed`, `eval_model` |

To compile a model dictionary into a `FlatModel`:

```julia
model = dict_to_model(model_dict, list_free_params)
x0 = Float64[model_dict[p] for p in list_free_params]
```

`model` and `x0` are then passed to the fitting and evaluation functions.

## Defining a model dictionary

Parameters use flat keys of the form `"component,parameter"`:

```julia
using OITOOLS, PythonPlot
data = readoifits("data/AlphaCenA.oifits")[1]

model_dict = Dict{String,Any}(
    "star,ud"  => 8.0,    # uniform disk diameter (mas)
    "star,f"   => 1.0,    # flux fraction
)
list_free_params = ["star,ud"]  # parameters to optimize

display_model(model_dict, list_free_params)
```

The component name (`"star"`, `"disk"`, …) is free; the component type is set by its geometry
key.

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

Five laws; the geometry key names the law and gives the diameter in mas:

```julia
# Linear: I(μ) = 1 - u(1-μ)
model_dict = Dict{String,Any}("star,ldlin" => 1.0, "star,u" => 0.3, "star,f" => 1.0)

# Quadratic: I(μ) = 1 - u(1-μ) - w(1-μ)²
model_dict = Dict{String,Any}("star,ldquad" => 1.0, "star,u" => 0.2, "star,w" => 0.1, "star,f" => 1.0)

# Square root: I(μ) = 1 - u(1-μ) - w(1-√μ)
model_dict = Dict{String,Any}("star,ldsqrt" => 1.0, "star,u" => 0.3, "star,w" => 0.2, "star,f" => 1.0)

# Power-law: I(μ) = μ^α   (Hestroffer)
model_dict = Dict{String,Any}("star,ldpow" => 1.0, "star,alpha" => 0.5, "star,f" => 1.0)

# Four-parameter: I(μ) = 1 - Σₖ cₖ(1 - μ^{k/2})   (Claret)
model_dict = Dict{String,Any}("star,ldclaret4" => 1.0, "star,c1" => 0.4, "star,c2" => -0.2,
                              "star,c3" => 0.5, "star,c4" => -0.15, "star,f" => 1.0)
```

!!! warning "Coefficients belong to a law"
    `u` is the linear coefficient of both the quadratic and the square-root law, but a value
    fitted under one law does not apply to another. Claret tabulates each law's coefficients
    separately, and a power-law α is a different quantity from the linear `u`.

### Uniform ring

Defined by inner and outer diameters in mas:

```julia
model_dict = Dict{String,Any}("ring,diamin" => 0.6, "ring,diamout" => 1.0, "ring,f" => 1.0)
```

or by `diam` (outer diameter) and `thick` (fractional thickness, 0–1):

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

Its flux is part of the total that the visibilities are divided by, so it lowers the other
components' visibilities by (F − f_bg)/F, and it is part of the flux compared with OI_FLUX. In
`model_to_image` it is a uniform level.

### Summary table

| Key(s) | Type | Description |
|--------|------|-------------|
| `ud` | Uniform disk | Diameter in mas |
| `fwhm` | Gaussian | FWHM in mas |
| `ldlin` + `u` | Linear limb-darkened disk | `I(μ) = 1 - u(1-μ)` |
| `ldquad` + `u`, `w` | Quadratic limb-darkened disk | `1 - u(1-μ) - w(1-μ)²` |
| `ldsqrt` + `u`, `w` | Square-root limb-darkened disk | `1 - u(1-μ) - w(1-√μ)` |
| `ldpow` + `alpha` | Power-law limb-darkened disk | `μ^α` (Hestroffer) |
| `ldclaret4` + `c1`…`c4` | Four-parameter limb-darkened disk | `1 - Σₖ cₖ(1 - μ^{k/2})` (Claret) |
| `diamin` + `diamout` | Uniform ring | Inner/outer diameters in mas |
| `diam` + `thick` | Uniform ring (sugar) | Outer diameter + fractional thickness |
| `fwhmin` + `fwhmout` | Gaussian ring | Inner/outer FWHM in mas |
| `crin`, `crout`, `croff`, `crprojang` | Crescent | Two offset disks |
| `profile` + `diamout` | Hankel profile | Arbitrary radial profile (see below) |
| *(only `f`)* | Point source | Unresolved |
| `resolved` | Resolved background | Fully resolved (V=0) |

!!! note "FWHM"
    `fwhm`, `fwhmin`/`fwhmout` and `spatial_kernel` are full widths at half maximum in mas. The
    visibility is

    ```
    V(B) = exp(-π² · FWHM² · B² / (4 ln 2))
    ```

    as in [mfit](https://github.com/jsy1001/mfit) (`exp(-((pi*a*rho)**2)/(4*log(2)))` for
    `a = FWHM`).

## Common parameters

Every component supports these optional keys:

| Key | Description |
|-----|-------------|
| `f` | Flux fraction (can be a number or an expression string) |
| `x`, `y` | Position offset in mas (east, north) |
| `incl` | Inclination in degrees (0 = face-on) |
| `pa` or `projang` | Position angle in degrees (north through east) |
| `spectrum` | Spectral law string; promoted to `f` internally (e.g. `"(\$WL/2.2)^-4"`) |
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

A value can be an expression referring to other parameters with `\$`, in any order:

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

These implicit variables are available in expressions:

| Variable | Where | Description |
|----------|-------|-------------|
| `\$WL` | any expression | Wavelength in microns, one value per uv point during fitting (every `wl` argument and the data stay in metres) |
| `\$MJD` | any expression | Modified Julian date, one value per uv point |
| `\$R` | `profile` expressions only | Radius in mas, on the component's radial grid |
| `\$MU` | `profile` expressions only | `sqrt(1 - (R/r_max)^2)`, for limb-darkening-style profiles |
| `\$BL` | `visfunc` expressions only | Baseline length in Mλ, one value per uv point |

`B` is not a variable, so `\$disk,B` is an ordinary parameter, as PMOIRED's notebooks write a
spectral slope.

`\$MJD` is read from `data.uv_mjd`, which is `Float64` for any element type `T`.

A parameter whose expression uses `\$WL` or `\$MJD`, directly or through another parameter, is
evaluated per data point; the others are evaluated once. Any parameter can depend on them:
diameters and position angles as well as fluxes.

A chromatic model:

```julia
model_dict = Dict{String,Any}(
    "star,f"        => 1.0,
    "star,spectrum"  => "(\$WL/2.2)^(-4)",    # Rayleigh–Jeans spectrum
    "disk,f"        => 0.3,
    "disk,spectrum"  => "(\$WL/2.2)^(-2)",    # greyer spectrum
    "disk,diamout"  => 10.0,
    "disk,profile"  => "exp(-(\$R/3.0)^2)",
)
```

A diameter varying with wavelength:

```julia
model_dict = Dict{String,Any}(
    # a diameter that grows linearly with wavelength
    "star,d0"    => 3.0,
    "star,slope" => 0.05,
    "star,ud"    => "\$star,d0 * (1 + \$star,slope * (\$WL/1.6 - 1))",
    "star,f"     => 1.0,
)
list_free_params = ["star,d0", "star,slope"]   # both are fitted as usual
```

This applies to every parameter of every analytic component, for example a limb-darkening
coefficient varying with time (`"star,u" => "0.3 + 0.01*(\$MJD - 60000)"`), and the fitters
compute the gradients.

### Shared parameters across components

Keys without a `component,` prefix are global parameters that components can refer to:

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

Each component has its own name prefix. The model visibility is the flux-weighted sum of the
components' visibilities:

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

The model image and SED:

```julia
params = dict_to_model(model_dict, String[])
img = model_to_image(params, Float64[]; nx=128, pixsize=0.02, wl=1.65e-6)
imdisp(img; pixsize=0.02)

wl_grid = collect(range(1.0e-6, 2.5e-6; length=200))
f_total, f_comps = model_to_sed(params, Float64[], wl_grid)
```

`demos/example_model_fitting_pmoired_models.jl` shows every component type.

## Hankel models (radial profiles)

A component with a `profile` key has a radial intensity profile, transformed numerically
(Hankel transform): disks, envelopes, other limb-darkening laws.

### Profile expressions

The profile is an expression in `\$R` (radius in mas), `\$MU` (μ) and other parameters:

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
| `diam`, `diamout` or `udout` | Outer diameter in mas, where the profile ends. It can be fitted |
| `diamin` | Inner diameter in mas (optional: the profile starts there, for a ring). It can be fitted |
| `nr` | Optional: fixes the number of quadrature nodes (at least 2). Without it the count is chosen from the data and the profile, see below |
| `r_max` | Alternative to `diamout`: the outer radius, a fixed number |

The diameters are parameters like any other: free them, or write them as expressions, and the
profile follows. `\$RMIN`, `\$RMAX`, `\$DMIN` and `\$DMAX` in the profile are the current
radii and diameters, and `\$MU` = √(1 − (R/R_max)²) the current μ, so `\$MU^\$alpha` with
`diam` free fits as `ldpow` does. With a free size (or inclination) the transform is recomputed
at every evaluation: about 0.3 ms at 600 baselines instead of a few µs.

### How a profile is integrated

A profile is integrated with Gauss–Legendre quadrature in panels:

- **Near the edge** μ = t² is used as the variable, so a limb-darkened edge (μ^α, which has an
  infinite slope at the rim) becomes smooth.
- **Before integrating, the profile is probed** at its current parameter values. Jumps
  (`(\$R >= \$rin)`), kinks (`max(0, …)`) and steep transitions become panel boundaries, and a
  centre where the profile is infinite (`\$R^(-1.5)`) gets a change of variable.
- **A jump whose position is a fitted parameter** is followed as it moves, and its derivative
  includes the term from the moving boundary. So does a jump at a fixed radius
  (`\$R >= 2.0`) when the size is fitted, and the inner edge when `diamin` is.
- **The node count** starts from the resolution of the data (about 32 nodes up to 7 lobes,
  then 1.1 per lobe), and is doubled until the visibilities and their parameter derivatives
  agree to 1e-8 (1e-6 in Float32).

The count is set when a dataset is first seen; if a fit then makes a feature much sharper, set
`nr`. For a ring, `diamin` starts the grid at the inner edge.

### Azimuthal modulations

Azimuthal harmonics have amplitudes `az ampN` and orientations `az projangN`, N = 1, 2, 3, …:

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

The harmonics add: the brightness is multiplied by `1 + Σ ampN cos(…)`.

### Spiral patterns

Concentric thin rings with an m = 1 modulation whose `az projang1` rotates with radius make a
spiral arm; `spatial_kernel` smooths it. Section 5b of
`demos/example_model_fitting_pmoired_models.jl` builds one from 5 rings and global parameters
(`Din`, `Dout`, `INCL`, `PROJANG`, `pitch`, `PAin`).

### YSO disk fitting example

`demos/example_model_fitting_v1295aql.jl` fits three chromatic YSO disk models (parabolic
ring, double sigmoid, α-sigmoid) to V1295 Aql (HD 190073) H-band data, as Ibrahim et al. 2023
(ApJ, 947, 68).

## Spatial smoothing

The global parameter `spatial_kernel` convolves the model with a Gaussian of that FWHM (mas):

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

The dictionary format is [PMOIRED](https://github.com/amerand/PMOIRED)'s.

### Converting Python dicts to Julia

`pmoired_to_dict()` converts a PMOIRED dict literal to a Julia `Dict`:

```julia
model_dict = pmoired_to_dict("{'star,ud': 3.2, 'ring,f': '1 - \$star,f'}")
```

`pmoired_to_julia()` returns the Julia source instead. The conversions:

| Python | Julia |
|--------|-------|
| `{ ... }` | `Dict( ... )` |
| `'key': value` | `"key" => value` |
| `'expr with \$ref'` | `raw"expr with \$ref"` |
| `True` / `False` / `None` | `true` / `false` / `nothing` |
| `[...]` (lists) | `[...]` (arrays) |
| `# comment` | `# comment` |

### Converting entire scripts

`pmoired_to_julia_file()` converts a script line by line:

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
using OITOOLS, PythonPlot

# Option 1: inline conversion
model_dict = pmoired_to_dict("{'star,ud': 3.2, 'star,f': 0.7, 'ring,udout': '\$star,ud * 8', 'ring,f': '1 - \$star,f', 'ring,incl': 30.0}")

# Option 2: file conversion
pmoired_to_julia_file("pmoired_model.py", "julia_model.jl")
include("julia_model.jl")

# Option 3: written in Julia
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
    In Julia, expressions containing `$` are written `raw"..."` or with `\$`;
    `pmoired_to_julia()` does this.

See `demos/example_model_fitting_pmoired_conversion.jl` and
`demos/example_model_fitting_pmoired_models.jl`.

## Fitting with NLopt (gradient-based)

`fit_model` uses NLopt. Gradient methods (default `:LD_LBFGS`) use the analytic gradient;
gradient-free methods such as `:LN_NELDERMEAD` are also available.

```julia
result = fit_model(model_dict, list_free_params, data;
    weights = [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],  # V² + T3φ
    method  = :LD_LBFGS,
    maxeval = 500,
    verb    = true)

println("Best χ²/ν = ", result.chi2r)
println("Parameters: ", result.x_opt)
```

`weights` has seven entries:
`[V², T3amp, T3φ, visamp, visφ, flux, diffφ]`.

### Bounds and priors

```julia
lb = Dict("star,ud" => 0.1)   # lower bounds
ub = Dict("star,ud" => 20.0)  # upper bounds

# Gaussian priors: (expression, target, sigma); the expression may combine parameters
priors = [("star,ud",           8.5, 0.5),      # prior on one parameter
          ("star,f + disk,f",   1.0, 0.01)]     # prior on their sum

result = fit_model(model_dict, list_free_params, data; lb, ub, priors)
```

Each prior adds `(value - target)^2 / sigma^2` to the objective. Priors need the `Dict` form of
`fit_model` (the `FlatModel` form raises an error).

`default_bounds(model_dict, list_free_params)` suggests bounds for every free parameter, and
`display_model` checks values against them. With `data`, angular sizes are bounded by
`2 λ/B_min`:

```julia
lb, ub = default_bounds(model_dict, list_free_params; data)
```

### Constraints

Relations between parameters are `ModelConstraint`s:

```julia
constraints = [ModelConstraint("ring,fwhmout", ">", "ring,fwhmin"),
               ModelConstraint("star,f + disk,f", "=", 1.0)]

result = fit_model(model_dict, list_free_params, data; lb, ub, constraints)
```

`op` is one of `<`, `<=`, `>`, `>=`, `=`. Either side may be an expression, and the right may
be a number. `tol` defaults to `1e-3`.

`fit_model` passes constraints to NLopt as nonlinear constraints, which hold at the optimum; an
algorithm that does not support them (including the default `:LD_LBFGS`) is wrapped in
`:AUGLAG`. `fit_model_lsqfit`, `fit_model_nested` and the samplers use a one-sided quadratic
penalty, which a steep χ² can outweigh.

`check_constraints(constraints, model_dict)` reports which constraints the starting model
satisfies.

### Saving a model and its fit settings

A TOML model file holds the model dict, the free parameters, the bounds, the constraints and
the priors.

```julia
write_model_file("binary.toml", model_dict; free = list_free_params, lb, ub, constraints, priors)

m = read_model_file("binary.toml")
result = fit_model(m.model, m.free, data; m.lb, m.ub, m.constraints, m.priors)
```

```toml
free = ["star,ud", "disk,pa"]

[model]
"star,ud"   = 6.5
"disk,fwhm" = "$star,ud * 3"

[bounds]
"star,ud" = [0.0, 20.0]

[[constraints]]
param = "disk,fwhm"
op    = ">"
value = "star,ud"
tol   = 0.001

[[priors]]
expr   = "star,ud"
target = 6.0
sigma  = 0.5
```

Keys are written sorted. To save a fitted model, merge the result into the dict:

```julia
fitted = merge(model_dict, Dict(zip(result.list_free_params, result.x_opt)))
write_model_file("binary_fitted.toml", fitted; free = result.list_free_params, lb, ub)
```

## Fitting with LsqFit (Levenberg-Marquardt)

`fit_model_lsqfit` (Levenberg-Marquardt) also returns the covariance and 1σ uncertainties,
scaled to χ²ᵣ = 1:

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

The Jacobian is analytic.

## Bayesian inference (nested sampling)

`fit_model_nested` returns posterior samples and the log-evidence, for comparing models. Every
free parameter needs bounds, which define the uniform prior. Two samplers implement it, each
enabled by loading its package:

| `using` | backend | notes |
|---|---|---|
| `Nautilus` | `:nautilus` | [Nautilus.jl](https://github.com/fabienbaron/Nautilus.jl), importance nested sampling in Julia; used by the compiled application |
| `PythonCall` | `:ultranest` | [UltraNest](https://johannesbuchner.github.io/UltraNest/), through PythonCall |

```julia
using OITOOLS, Nautilus

result = fit_model_nested(model_dict, list_free_params, data;
    lb = Dict("star,ud" => 0.1),
    ub = Dict("star,ud" => 20.0),    # bounds required for all list_free_params
    n_live = 500)

println("log(Z) = ", result.logz, " ± ", result.logzerr)
println("Best χ²/ν = ", result.chi2r)
println("sampler  = ", result.backend)
```

`result.posterior` is an `nsamples × nparams` matrix of equally weighted samples.

`nested_backend()` reports the default sampler and `set_nested_backend!` changes it. With both
loaded, the two `logz` should agree within `logzerr`:

```julia
using OITOOLS, Nautilus, PythonCall           # both backends

rj = fit_model_nested(model_dict, list_free_params, data; lb, ub, backend = :nautilus)
ru = fit_model_nested(model_dict, list_free_params, data; lb, ub, backend = :ultranest)
abs(rj.logz - ru.logz) <= 3 * sqrt(rj.logzerr^2 + ru.logzerr^2)
```

`fit_model_ultranest` is `fit_model_nested` with `backend = :ultranest`.

## Posterior samplers

`fit_model_ensemble` (Ensemble Slice Sampling, in OITOOLS) and `fit_model_pigeons` (parallel
tempering, needs `using Pigeons`) return posterior samples with the bounds as a uniform prior.
`fit_model_pigeons` also returns the log-evidence on `fit_model_nested`'s scale.

```julia
r = fit_model_ensemble(model_dict, list_free_params, data; lb, ub, nsweep = 2000, seed = 1)
r.median, r.q16, r.q84        # per parameter
r.diagnostics.rhat            # R̂ per parameter

using Pigeons
r = fit_model_pigeons(model_dict, list_free_params, data; lb, ub, n_rounds = 12,
                      variational = true, n_chains_variational = 10)
r.logz, r.diagnostics
```

`plot_corner_makie(r.posterior, r.list_free_params)` (with PairPlots) draws any of these
posteriors.

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

## Uncertainty estimation

| Method | Assumes | Accounts for |
|--------|---------|--------------|
| `fit_model_lsqfit` covariance | error bars correct and uncorrelated; model locally linear | the quoted errors, rescaled to χ²ᵣ = 1 |
| `bootstrap_fit` | many exchangeable blocks of data | correlated calibration errors, wrong error bars, a bad night or baseline |

With correlated calibration errors the analytic covariance is 3–5× too small; the
[Uncertainties](@ref) page has the study.

### Block bootstrap

```julia
boot = bootstrap_fit(model_dict, list_free_params, data;
    lb=lb, ub=ub, weights=[1.0, 0.0, 0.0], nboot=500, seed=42)

boot.median        # median of the bootstrap distribution
boot.sigma         # (84th − 16th percentile) / 2
boot.sigma_plus    # asymmetric error bars
boot.samples       # nboot × npar matrix, for corner plots or derived quantities
boot.covar         # parameter covariance
```

The resampling unit is set by `granularity`:

- `:config` (default) — one block per (MJD, baseline / triangle / telescope), all channels
  together (PMOIRED's spectral vector).
- `:epoch` — one block per MJD; needs many epochs.
- `:point` — one block per data point; underestimates the uncertainties on real data.

`mode` selects how block multiplicities are drawn:

| Mode | Scheme | Measured ratio, many blocks | few blocks | Cost |
|------|--------|------------------|------------|------|
| `:replacement` (default) | multinomial | 0.87–1.01 | 0.96–1.14 | 1 dataset |
| `:halfsample` | balanced repeated replication: a random half | **0.94–1.06** | 1.24–1.56 | ½ dataset |
| `:weights` | multiplier (Bayesian) bootstrap: continuous block weights on the error bars | 0.78–0.99 | 0.71–0.85 | 1 dataset |
| `:pmoired` | PMOIRED's two independent half-samples, fitted jointly | 0.62–0.74 ⚠ | 0.67–0.79 ⚠ | 1 dataset |

"Measured ratio" is the quoted σ divided by the true parameter scatter over 120
simulated realisations (1.00 = calibrated; see [Uncertainties](@ref)). `:pmoired` reproduces
PMOIRED's scheme for comparison; its multiplicities have variance ½, so its error bars are about
√2 too small, and `bootstrap_fit` warns on its first use.

Derived quantities come from the samples:

```julia
sep = sqrt.(boot.samples[boot.mask, 1].^2 .+ boot.samples[boot.mask, 2].^2)
println("separation = ", median(sep), " ± ", 0.5*(quantile(sep, 0.84) - quantile(sep, 0.16)))
```

### Adding noise to the data

`perturb_data(data)` adds a Gaussian draw from each error bar to every observable. Refitting
such replicates reproduces the analytic covariance; in the study it was 11–25× too small with
calibration systematics, and too small by the underestimation factor when the error bars were
underestimated. `resample_data` is a deprecated alias.

`demos/example_bootstrap_fit.jl` compares the schemes.
