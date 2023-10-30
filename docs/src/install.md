# OITOOLS Installation

## Overview of Dependencies

| Package       | Usage     | Algorithm/Functions |
| ------------- |:-------------:|:-------------:|
| NFFT      | Compute Fourier transform | Non Equispaced Fourier Transform |
| OptimPackNextGen | Image reconstruction | VMLMB (Éric Thiébaut)
| LsqFit      | Model fitting | Levenberg-Marquardt
| UltraNest | Model fitting and model comparison  |  Nested sampling    |
| NLopt     | Model fitting | Several local (Nelder-Mead) and global (Genetic) optimizers |
| SpecialFunctions | Complex visibility calculations | Bessel Functions
| NearestNeighbors | simplify uv sampling | KD trees |
| OIFITS | data import | read OIFITS files|


## Step 1: installation of Python Packages (UltraNest, Astroquery)

OITOOLS downloads then uses a Conda installation with the UltraNest ad Astroquery packages. For this you should copy/paste the following line into the REPL:
```julia
ENV["PYTHON"]=""; ENV["MPLBACKEND"]="qt5agg"; using Pkg; Pkg.add("Conda"); using Conda; Conda.add("ultranest", channel="conda-forge"); Conda.add("astroquery", channel="astropy");
```
These operations should take a couple of minutes to complete.
Feel free to remove the two ENV settings at the beginning of the line if you have already Python setup on your machine and do not want to use the Julia Conda python.


## Step 2: installation of Julia Packages

Because some of OITOOLS dependencies are not registered packages, we elect to go through ```Pkg()``` again rather than the activate/instantiate mechanism of Julia. Here again, copy/paste the following line into the REPL:
```julia
using Pkg; Pkg.add(["CFITSIO","AstroTime","Dates","DelimitedFiles","Documenter","DocumenterTools","FITSIO","Glob","LaTeXStrings","LinearAlgebra","FFTW", "NFFT","NLopt","UltraNest","LsqFit","NearestNeighbors","PyCall","PyPlot","Random","SparseArrays","SpecialFunctions","Statistics","Parameters"]); Pkg.add(url="https://github.com/fabienbaron/OIFITS.jl", rev="t4"); Pkg.add(url="https://github.com/emmt/ArrayTools.jl.git");Pkg.add(url="https://github.com/emmt/LazyAlgebra.jl.git"); Pkg.add(url="https://github.com/emmt/OptimPackNextGen.jl.git");Pkg.add(url="https://github.com/fabienbaron/OITOOLS.jl.git")
```
Installation may take between 2-10 minutes depending on OS and computer performance.

