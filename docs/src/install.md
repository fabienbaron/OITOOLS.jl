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
| PyCall / PyPlot | Matplotlib-based plotting |
| Crayons | Coloured terminal output |

## Step 1: Python packages

OITOOLS uses a Julia-managed Conda environment for UltraNest and Astroquery.
Paste the following into the Julia REPL:

```julia
ENV["PYTHON"] = ""
ENV["MPLBACKEND"] = "qt5agg"
using Pkg
Pkg.add("Conda")
using Conda
Conda.add("ultranest",  channel="conda-forge")
Conda.add("astroquery", channel="astropy")
```

!!! tip

    If you already have a Python installation you want to use, omit the two
    `ENV[...]` lines so Julia reuses your existing Python.

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
    "Parameters", "PyCall", "PyPlot", "Random", "SparseArrays",
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

