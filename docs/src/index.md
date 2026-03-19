# OITOOLS

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/fabienbaron/OITOOLS.jl)
[![Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build](https://github.com/fabienbaron/OITOOLS.jl/workflows/CI/badge.svg)](https://github.com/fabienbaron/OITOOLS.jl/actions)

**OITOOLS** is a [Julia](http://julialang.org/) package for optical interferometry,
developed by Prof. Fabien Baron (Georgia State University) and collaborators.
It covers the full data analysis workflow:

- **Reading** OIFITS v1/v2 files with automatic quality filtering, spectral and
  temporal binning, and OI_CORR support
- **Plotting** UV coverage, V², closure phases, triple amplitudes, and flux spectra
- **Model fitting** with NLopt, Levenberg–Marquardt, and UltraNest (Bayesian evidence)
- **Image reconstruction** via gradient descent (VMLMB) with NFFT, including
  polychromatic, time-variable, and SPARCO modes
- **Simulating** synthetic OIFITS datasets from images or analytic models

Source code: [github.com/fabienbaron/OITOOLS.jl](https://github.com/fabienbaron/OITOOLS.jl)

## Forward model pipeline

Both image-based and parametric-model-based workflows follow the same pattern:

| Operation | Image domain | Model domain |
|-----------|-------------|-------------|
| → complex visibilities | `image_to_vis(x, ft)` | `model_to_vis(model, x, uv)` |
| → observables | `image_to_obs(x, ft, data)` | `model_to_obs(model, x, data)` |
| → chi² | `image_to_chi2(x, ft, data)` | `model_to_chi2(model, x, data)` |
| → chi² + gradient | `image_to_chi2_fg(x, g, ft, data)` | `model_to_chi2_fg(model, x, data)` |

!!! note

    OITOOLS is under active development. Contributions to the code, documentation,
    and demos are welcome.
