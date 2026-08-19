# Installation

## Dependencies

| Package | Purpose |
|---------|---------|
| [OIFITS.jl](https://github.com/emmt/OIFITS.jl) | Read OIFITS v1/v2 files |
| NFFT | Non-uniform fast Fourier transform for imaging |
| OptimPackNextGen | VMLMB gradient optimizer for image reconstruction (Éric Thiébaut) |
| NearestNeighbors | KD-tree UV deduplication |
| SparseArrays | OI_CORR correlation matrix storage |
| NLopt | Local and global optimizers (Nelder–Mead, ISRES, …) |
| LsqFit | Levenberg–Marquardt least-squares fitting |
| UltraNest | Bayesian model selection via nested sampling |
| SpecialFunctions | Bessel functions for visibility calculations |
| PythonCall | Julia↔Python bridge (SIMBAD queries, UltraNest) |
| PythonPlot | **optional** — matplotlib plotting, see [Plotting](@ref plotting-install) |
| Crayons | Coloured terminal output |

## Step 1: Python packages

Nothing to do. OITOOLS declares its Python dependencies in `CondaPkg.toml`, and
[PythonCall](https://github.com/JuliaPy/PythonCall.jl) provisions them automatically into a
project-local environment the first time you `using OITOOLS`:

| Package | Used by |
|---|---|
| matplotlib-base | `fit_model_ultranest`'s corner plot, and the plotting extension below |
| astroquery | `magnitudes_from_simbad`, `ra_dec_from_simbad` |
| ultranest | `fit_model_ultranest` |
| scipy | imported at module level by `ultranest.plot` |

`matplotlib-base`, not `matplotlib`: on conda-forge the latter is a metapackage adding
PySide6, which pulls in Qt, libclang and libLLVM — about 510 MB unpacked for a GUI binding
nothing here imports. PythonPlot declares full `matplotlib` itself, so an interactive backend
arrives with it if you want one.

!!! tip "Reusing an existing Python"

    Set `JULIA_PYTHONCALL_EXE` to your interpreter to reuse an existing installation. Note
    that doing so opts out of `CondaPkg`, so you become responsible for keeping matplotlib,
    astroquery and ultranest installed and mutually consistent.

## Step 2: Julia packages

Add the EmmtRegistry (required for OptimPackNextGen, OIFITS, and related packages),
then install all dependencies:

```julia
using Pkg
pkg"registry add General"
pkg"registry add https://github.com/emmt/EmmtRegistry"
Pkg.add([
    "AstroTime", "CFITSIO", "Crayons", "Dates", "DelimitedFiles",
    "FFTW", "FITSIO", "Glob", "LaTeXStrings", "LinearAlgebra",
    "LsqFit", "Match", "NFFT", "NLopt", "NearestNeighbors",
    "Parameters", "Random", "SparseArrays",
    "SpecialFunctions", "Statistics", "UltraNest",
    "ArrayTools", "LazyAlgebra", "OptimPackNextGen", "OIFITS",
])
Pkg.add(url="https://github.com/fabienbaron/OITOOLS.jl.git")
```

Installation typically takes 2–10 minutes depending on your system.

Verify the installation:

```julia
using OITOOLS
```

## [Step 3: Plotting](@id plotting-install)

**The plotting functions live in a package extension and require PythonPlot.** `using OITOOLS`
alone gives you `uvplot`, `plot_v2`, `imdisp` and the rest as functions with no methods, and
calling one raises a `MethodError`.

```julia
using Pkg
Pkg.add("PythonPlot")

using OITOOLS, PythonPlot     # both: PythonPlot is what activates the extension
uvplot(data)
```

The extension (`OITOOLSPythonPlotExt`) loads automatically once PythonPlot is present; there
is nothing else to enable. It supplies every plotting entry point:

`uvplot`, `plot_v2`, `plot_t3phi`, `plot_t3amp`, `plot_visamp`, `plot_visphi`, `plot_diffphi`,
`plot_flux`, `plot_obs`, `plot_multi`, the `*_residuals` family, `plot_residuals`,
`plot_v2_multifile`, `plot_facility`, `imdisp`, `imdisp_multi`, `plot_ultranest_corner`,
`set_oiplot_defaults`, and the observation planners `gantt_onenight`, `obs_plan`,
`chara_plan`, `empty_night`.

Keeping them out of the core package means `using OITOOLS` starts no Python plotting stack:
it does not import matplotlib, and it maps no Qt into the process — which matters on a
headless machine, in CI, and beside any other Qt application.

To pick a non-interactive backend, for scripts or CI, set `MPLBACKEND` before loading
PythonPlot:

```julia
ENV["MPLBACKEND"] = "Agg"
using OITOOLS, PythonPlot
```

Observable metadata — `OBS_PLOT_SPECS`, `oiplot_colors`, `canonical_color` — stays in the core
package, so a front-end that draws with something other than matplotlib can share the axis
labels, groupings and palette without loading Python at all.

