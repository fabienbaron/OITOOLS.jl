# Installation

OITOOLS needs Julia 1.13 or later. It is not in the General registry; some of its dependencies
(OIFITS, OptimPackNextGen) are in [EmmtRegistry](https://github.com/emmt/EmmtRegistry):

```julia
using Pkg
pkg"registry add General"
pkg"registry add https://github.com/emmt/EmmtRegistry"
Pkg.add(url = "https://github.com/fabienbaron/OITOOLS.jl.git")
```

`using OITOOLS` loads no plotting package, no Python and no GUI. These come from optional
packages, each enabling part of the package when it is loaded:

| package | enables | install |
|:--|:--|:--|
| PythonPlot | the matplotlib figures (`uvplot`, `plot_v2`, `imdisp`, …) | `Pkg.add("PythonPlot")` |
| PythonCall | UltraNest (`fit_model_nested(…; backend = :ultranest)`) | `Pkg.add("PythonCall")` |
| Nautilus | nested sampling in Julia | `Pkg.add(url = "https://github.com/fabienbaron/Nautilus.jl.git")` |
| Pigeons | parallel tempering (`fit_model_pigeons`, `reconstruct_squeeze_tempered`) | `Pkg.add("Pigeons")` |
| PairPlots | corner plots (`plot_corner_makie`) | `Pkg.add("PairPlots")` |
| VarInf | variational inference for imaging | `Pkg.add(url = "https://github.com/fabienbaron/VarInf.jl.git")` |
| GLMakie, QMLMakie, QML | the GUI and the `_makie` figures | `Pkg.add(["GLMakie", "QMLMakie", "QML"])` |
| GLFW_jll | native Wayland for the GUI | `Pkg.add("GLFW_jll")` |

Main dependencies: OIFITS.jl and AstroFITS (files), NFFT and FFTW (Fourier transforms),
OptimPackNextGen (VMLMB, image reconstruction), NLopt and LsqFit (model fitting),
NearestNeighbors (uv deduplication), SpecialFunctions, ForwardDiff and
RuntimeGeneratedFunctions (models).

## Python

Python is used only through PythonCall (UltraNest) and PythonPlot (matplotlib). When either is
first loaded, CondaPkg installs `CondaPkg.toml`'s packages (python, matplotlib-base, numpy,
scipy, ultranest) into a project-local environment. `matplotlib-base` has no GUI backend;
PythonPlot installs full `matplotlib`.

To use an existing Python instead, set `JULIA_PYTHONCALL_EXE` to its interpreter; matplotlib and
ultranest must then be installed there.

## What each optimiser needs

A function whose package is not loaded raises a `MethodError`, or, for `fit_model_nested` and
`fit_model_pigeons`, an error naming the package.

### Model fitting

| method | function | needs |
|---|---|---|
| NLopt (any `LD_`/`LN_` method) | `fit_model(…; method = :LD_LBFGS)` | — |
| Levenberg–Marquardt | `fit_model_lsqfit` | — |
| Grid search, two parameters | `chi2_map` | — |
| Bootstrap | `bootstrap_fit` | — |
| Ensemble slice sampling | `fit_model_ensemble` | — |
| Nested sampling, Julia | `fit_model_nested` | Nautilus |
| Nested sampling, UltraNest | `fit_model_nested(…; backend = :ultranest)`, `fit_model_ultranest` | PythonCall |
| Parallel tempering | `fit_model_pigeons` | Pigeons |

With both nested samplers loaded, `nested_backend()` reports the default and
`set_nested_backend!` changes it.

### Imaging

| method | function | needs |
|---|---|---|
| VMLMB with regularisers | `reconstruct` | — |
| Maximum entropy | `reconstruct_bsmem` | — |
| ADMM | `reconstruct_bsdmm` | — |
| SPARCO, chromatic hybrid | `reconstruct_sparco`, `reconstruct_hybrid` | — |
| SQUEEZE (simulated annealing) | `reconstruct_squeeze` | — |
| SQUEEZE, parallel tempering | `reconstruct_squeeze_tempered` | Pigeons |
| Variational inference | `reconstruct_map`, `reconstruct_mgvi` | VarInf |

### Figures

| figure | function | needs |
|---|---|---|
| matplotlib figures | see [Plotting](@ref) | PythonPlot |
| Makie figures (`*_makie`) | see [Plotting](@ref) | a Makie backend |
| Corner plot | `plot_corner_makie` | PairPlots and a Makie backend |
| UltraNest corner plot | `plot_ultranest_corner` | PythonPlot |
| SQUEEZE monitor | `reconstruct_squeeze(…; monitor = n)` | PythonPlot |

## [Plotting](@id plotting-install)

```julia
Pkg.add("PythonPlot")
using OITOOLS, PythonPlot
uvplot(data)
```

Loading PythonPlot activates the extension that provides `uvplot`, `plot_v2`, `plot_t3phi`,
`plot_t3amp`, `plot_visamp`, `plot_visphi`, `plot_diffphi`, `plot_flux`, `plot_obs`,
`plot_multi`, the `*_residuals` functions, `plot_residuals`, `plot_v2_multifile`,
`plot_facility`, `imdisp`, `imdisp_multi`, `plot_ultranest_corner`, `set_oiplot_defaults`,
`gantt_onenight`, `obs_plan`, `chara_plan` and `empty_night`. Without it these functions have no
methods.

For scripts without a display, set the backend before loading PythonPlot:

```julia
ENV["MPLBACKEND"] = "Agg"
using OITOOLS, PythonPlot
```

`OBS_PLOT_SPECS`, `oiplot_colors` and `canonical_color` (axis labels, groupings, palette) are in
the core package, for other plotting front-ends.

## [The GUI](@id gui-install)

From a clone of the repository, in its root directory:

```bash
julia -e 'using Pkg; Pkg.Registry.add("General"); Pkg.Registry.add(url = "https://github.com/emmt/EmmtRegistry")'   # once
julia --project=bin -e 'using Pkg; Pkg.instantiate()'   # once
julia --project=bin bin/oitoolsgui.jl [file.oifits]
```

From Julia, with GLMakie, QMLMakie and QML installed:

```julia
using OITOOLS
oitoolsgui()                       # or oitoolsgui("mydata.oifits")
```

Call `oitoolsgui()` before `using GLMakie` in that session.

A sysimage makes the start faster:

```bash
julia --project=bin bin/compile.jl                  # once; tens of minutes
julia --project=bin --sysimage bin/oitoolsgui.so bin/oitoolsgui.jl [file.oifits]
```

Rebuild it after updating OITOOLS or Julia. `app/build.jl` builds a standalone application, and
`app/appimage.sh` packs it as a Linux AppImage.
