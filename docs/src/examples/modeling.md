# Model fitting

OITOOLS supports analytic model fitting via several optimisation backends.
Models are assembled from one or more components, each with a type and a set of
free parameters.

## Creating a model

```julia
using OITOOLS
data = readoifits("data/AlphaCenA.oifits")[1,1]

model = create_model(create_component(type="ldpow", name="Star"))
display(model)   # lists free parameters and current values
```

Component types include: `"ud"` (uniform disc), `"ldpow"` (power-law limb darkening),
`"ldlin"` (linear limb darkening), plus point source, binary, and Gaussian components.

## Optimisation

```julia
# Nelder–Mead via NLopt (default)
minf, minx, cvis_model, result = fit_model_nlopt(data, model)

# Global search (ISRES genetic algorithm)
minf, minx, cvis_model, result = fit_model_nlopt(data, model, fitter=:GN_ISRES)

# Levenberg–Marquardt
minf, minx, cvis_model, result = fit_model_leastsq(data, model)
```

The `weights=[w_v2, w_t3amp, w_t3phi]` keyword controls per-observable weighting.

See `example_model_fitting_limb_darkening.jl`.

## Bootstrap error estimation

```julia
boot_results = bootstrap_fit(data, model, nboot=200, resample=:baseline)
```

See `example_bootstrap_fit.jl`.

## Bayesian model selection (UltraNest)

```julia
logz, result = fit_model_ultranest(data, model)
```

Nested sampling returns the log-evidence for model comparison (e.g. comparing
limb-darkening laws). See `example_model_fitting_functions.jl`.
