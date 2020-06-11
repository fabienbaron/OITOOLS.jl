# OITOOLS: the All-in-One Tool Package for Optical Interferometry

| **Documentation**               | **License**                     |**Build**                      |
|:--------------------------------|:--------------------------------|:------------------------------|
| [![][doc-dev-img]][doc-dev-url] | [![][license-img]][license-url] | [![][travis-img]][travis-url] |

OITOOLS is a Julia package to read, plot, model-fit and image optical interferometric data coming from astronomical arrays such as CHARA, VLTI, and NPOI.

[doc-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[doc-dev-url]: https://fabienbaron.github.io/OITOOLS.jl/dev

[license-url]: ./LICENSE.md
[license-img]: http://img.shields.io/badge/license-GPL3-brightgreen.svg?style=flat

[travis-img]: https://travis-ci.com/fabienbaron/OITOOLS.jl.svg?branch=master
[travis-url]: https://travis-ci.com/fabienbaron/OITOOLS.jl

## Coming soon / TODO list

- Observation planning: validate the entire chain, check daylight saving days, handle multiple targets,

- Model fitting: multiple component models, polychromatism

- Simulations: better noise models, examples for polychromatic and time-variable simulations.

- Imaging: example with SPARCO model, port of SQUEEZE to Julia, image comparison tools.
