# Model fitting data

## Simple model fitting

`example_limb_darkening_fit.jl` : given OIFITS data, do model-fitting (uniform disc, limb-darkened disc)
`example_satlas_fit.jl`         : model OIFITS data using a SATLAS model (open-source stellar atmosphere model code)

## Bootstrapping errors in fits

`example_bootstrap_fit.jl`      : use the boostrap method to estimate error bars

## Bayesian model selection

`example_nested_sampling_fit.jl`: use Bayesian model selection via nested sampling to compare limb-darkening laws
