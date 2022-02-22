# OITOOLS Installation

## Overview of Dependencies

| Package       | Usage     | Algorithm/Functions |
| ------------- |:-------------:|:-------------:|
| NFFT      | Compute Fourier transform | Non Equispaced Fourier Transform |
| OptimPackNextGen | Image reconstruction | VLMBM (Éric Thiébaut)
| LsqFit      | Model fitting | Levenberg-Marquardt
| UltraNest | Model fitting and model comparison  |  Nested sampling    |
| NLopt     | Model fitting | Several local (Nelder-Mead) and global (Genetic) optimizers |
| SpecialFunctions | Complex visibility calculations | Bessel Functions
| NearestNeighbors | simplify uv sampling | KD trees |
| OIFITS | data import | read OIFITS files|

## Step 1: installation of Python Packages (UltraNest, Astroquery)

OITOOLS downloads then uses a Conda installation with the UltraNest ad Astroquery packages. For this you should copy/paste the following line into the REPL:
```julia
using Pkg; Pkg.add("Conda"); using Conda; Conda.add("ultranest", channel="conda-forge"); Conda.add("astroquery", channel="astropy");
```
These operations should take a couple of minutes to complete.

### Step 2: installation of Julia Packages

Because some of OITOOLS dependencies are not registered packages, we elect to go through ```Pkg()``` again rather than the activate/instantiate mechanism of Julia. Here again, copy/paste the following line into the REPL:
```julia
using Pkg; Pkg.add(["CFITSIO","AstroTime","DelimitedFiles","Documenter","DocumenterTools","FITSIO","Glob","LaTeXStrings","LinearAlgebra","NFFT","NLopt","UltraNest","LsqFit","NearestNeighbors","PyCall","PyPlot","Random","SparseArrays","SpecialFunctions","Statistics","Parameters"]); Pkg.add(url="https://github.com/fabienbaron/OIFITS.jl", rev="t4"); Pkg.add(url="https://github.com/emmt/ArrayTools.jl.git");Pkg.add(url="https://github.com/emmt/LazyAlgebra.jl.git"); Pkg.add(url="https://github.com/emmt/OptimPackNextGen.jl.git");Pkg.add(url="https://github.com/fabienbaron/OITOOLS.jl.git")
```
Installation may take between 2-10 minutes depending on OS and computer performance.

# Installation FAQ

### How do I install/use Python?

OITOOLS uses matplotlib through the Julia library [PyPlot](https://github.com/JuliaPy/PyPlot.jl), which will use your system Python. If you do not want to deal with python package management, using Julia's Conda package instead of the system-wide python is recommended. To do so, just type:

```julia
ENV["PYTHON"]=""  # this will tell Julia there is currently no usable python
```

before installing the Conda package via:
```julia
]add Conda
```

### Which Matplotlib backend should I use?

You will have to select a Matplotlib backend for rendering plots. Qt5 works well in our experience. To set it up you will need to call

```julia
ENV["MPLBACKEND"]="qt5agg"
```

before using Matplotlib (or even better, before installing PyPlot).
