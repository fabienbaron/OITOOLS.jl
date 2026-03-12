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

!!! note

    OITOOLS is under active development. Contributions to the code, documentation,
    and demos are welcome.

## API index

```@index
Pages = ["examples/reading.md", "examples/plotting.md", "examples/imaging.md"]
Modules = [OITOOLS]
```
