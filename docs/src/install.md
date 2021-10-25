# OITOOLS Installation


## Python Prerequisites

### UltraNest

As of version 0.6, OITOOLS uses UltraNest for model-fitting. You may install it as part of your local python install, or you may use Julia's Conda (see FAQ below), in which case you can type:

```julia    
using Conda
Conda.add("ultranest", channel="conda-forge")
```

### Astroquery (optional)

If you're going to use OITOOLS to simulate observations, you may need to look up targets with SIMBAD. OITOOLS already includes some code to call Astroquery. To use it you need to add astroquery to the python installation used by Julia. If you use Julia's Conda:

```julia    
using Conda
Conda.add("astroquery", channel="astropy")
```

## Installation of Dependencies

Overview of dependency usage:

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

### How to install dependencies: options 1, in a Julia environment

Simply clone OITOOLS using ```git clone```, cd to the project directory and use the package manager (```]``` key) then:
```
(v1.6) pkg> activate .

(OITOOLS) pkg> instantiate
```
This will install the packages in the same state that is given by the OITOOLS Manifest.toml. Otherwise, it will resolve the latest versions of the dependencies compatible with the project.

### How to install dependencies: options 2, global installation

In case the previous instantiation does not work, or if you want to install all dependencies globally, you can use the package manager (```]``` key) then do:

```julia
add https://github.com/fabienbaron/OIFITS.jl#t4 https://github.com/emmt/ArrayTools.jl.git https://github.com/emmt/LazyAlgebra.jl.git https://github.com/emmt/OptimPackNextGen.jl.git CFITSIO AstroTime DelimitedFiles Documenter DocumenterTools FITSIO LaTeXStrings LinearAlgebra NFFT NLopt UltraNest LsqFit NearestNeighbors PyCall PyPlot Random SparseArrays SpecialFunctions Statistics Parameters https://github.com/fabienbaron/OITOOLS.jl.git
```

# Installation FAQ

### How do I install/use Python?
    OITOOLS uses matplotlib through the Julia library [PyPlot](https://github.com/JuliaPy/PyPlot.jl), which will use your system Python. If you do not want to deal with python package management, using Julia's Conda package instead of the system-wide python is recommended. To do so, just type:

    ```julia
    ENV["PYTHON"]=""  # this will tell Julia there is currently no usable python
    ```

    before installing the Conda package ```]add Conda```.

### Which Matplotlib backend should I use?
    You will have to select a Matplotlib backend for rendering plots. Qt5 works well in our experience. To set it up you will need to call

    ```julia
    ENV["MPLBACKEND"]="qt5agg"
    ```

    before using Matplotlib (or even better, before installing PyPlot).
