# OITOOLS

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/fabienbaron/OITOOLS.jl)
[![Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build](https://github.com/fabienbaron/OITOOLS.jl/workflows/CI/badge.svg)](https://github.com/fabienbaron/OITOOLS.jl/actions)

**OITOOLS** is a [Julia](http://julialang.org/) package for optical interferometry,
developed by Prof. Fabien Baron (Georgia State University) and collaborators.
It covers:

- **Reading** OIFITS v1/v2 files: quality filtering, spectral and temporal binning, OI_CORR,
  checking and repairing files, writing OIFITS v2
- **Plotting** uv coverage and observables, with matplotlib or Makie
- **Model fitting** with NLopt, Levenberg–Marquardt, nested sampling (Nautilus, UltraNest),
  ensemble slice sampling and parallel tempering (Pigeons); companion search, orbits, block
  bootstrap
- **Image reconstruction** by gradient descent (VMLMB), maximum entropy (BSMEM), BSDMM and
  MCMC (SQUEEZE), including polychromatic, time-variable and SPARCO modes
- **Simulating** OIFITS datasets from images or models, with a noise model matched to ASPRO
- **Observation planning**: observability, delay lines, POP selection
- A **graphical interface** covering all of the above ([The GUI](@ref gui-install))

Source code: [github.com/fabienbaron/OITOOLS.jl](https://github.com/fabienbaron/OITOOLS.jl)

## Forward model pipeline

Images and parametric models have the same functions:

| Operation | Image domain | Model domain |
|-----------|-------------|-------------|
| → complex visibilities | `image_to_vis(x, ft)` | `model_to_vis(model, x, uv)` |
| → observables | `image_to_obs(x, ft, data)` | `model_to_obs(model, x, data)` |
| → residuals | `image_to_residuals(x, ft, data)` | `model_to_residuals(model, x, data)` |
| → chi² | `image_to_chi2(x, ft, data)` | `model_to_chi2(model, x, data)` |
| → chi² + gradient | `image_to_chi2_fg(x, g, ft, data)` | `model_to_chi2_fg(model, x, data)` |

!!! note

    OITOOLS is under active development. Contributions to the code, documentation,
    and demos are welcome.
