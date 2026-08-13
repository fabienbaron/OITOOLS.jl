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
| PythonCall / PythonPlot | Matplotlib-based plotting |
| Crayons | Coloured terminal output |

## Step 1: Python packages

Nothing to do. OITOOLS declares its Python dependencies in `CondaPkg.toml`, and
[PythonCall](https://github.com/JuliaPy/PythonCall.jl) provisions them automatically into a
project-local environment the first time you `using OITOOLS`:

| Package | Used by |
|---|---|
| matplotlib | all plotting (`oiplot.jl`) |
| astroquery | `magnitudes_from_simbad`, `ra_dec_from_simbad` |
| ultranest | `fit_model_ultranest` |

To pick a non-interactive plotting backend (for scripts or CI), set `MPLBACKEND` before
loading OITOOLS:

```julia
ENV["MPLBACKEND"] = "Agg"     # or "qt5agg" for an interactive window
using OITOOLS
```

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

