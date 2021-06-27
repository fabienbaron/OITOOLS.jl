# OITOOLS Installation


## Python Prerequisites

### UltraNest

As of version 0.6, OITOOLS uses UltraNest for model-fitting. You may install it as part of your local python install, or you may use Julia's Conda (see FAQ below), in which case you can type:
```    
using Conda
Conda.add("ultranest", channel="ultranest")
```

### Astroquery (optional)

If you're going to use OITOOLS to simulate observations, you may need to look up targets with SIMBAD. OITOOLS already includes some code to call Astroquery. To use it you need to add astroquery to the python installation used by Julia. If you use Julia's Conda:

```    
using Conda
Conda.add("astropy", channel="astropy")
Conda.add("astroquery", channel="astropy")
```

## Installation of Dependencies

### Options 1: in a Julia environment

Simply clone OITOOLS using ```git clone```, cd to the project directory and use the package manager (```]``` key) then:
```
(v1.6) pkg> activate .

(OITOOLS) pkg> instantiate
```
This will install the packages in the same state that is given by the OITOOLS Manifest.toml. Otherwise, it will resolve the latest versions of the dependencies compatible with the project.

### Options 2: global Installation of dependencies

In case the previous instantiation does not work, or if you want to install all dependencies globally, you can use the package manager (```]``` key) then do:

``` add https://github.com/fabienbaron/OIFITS.jl#t4 https://github.com/emmt/ArrayTools.jl.git https://github.com/emmt/LazyAlgebra.jl.git https://github.com/emmt/OptimPackNextGen.jl.git CFITSIO Dates DelimitedFiles Documenter DocumenterTools FITSIO LaTeXStrings LinearAlgebra NFFT NLopt UltraNest LsqFit NearestNeighbors PyCall PyPlot Random SparseArrays SpecialFunctions Statistics Parameters https://github.com/fabienbaron/OITOOLS.jl.git```

!!! info
    OITOOLS uses the ```OIFITS```, ```NFFT```, ```SpecialFunctions``` and ```NearestNeighbors``` packages. For model fitting, ```NLopt``` (derivative-free local and global optimizers), ```LsqFit``` (Levenberg-Marquardt with local error estimation) ```UltraNest``` (MCMC, model selection and best error estimates) are used. For image reconstruction, we use ```OptimPackNextGen``` written by Éric Thiébaut.


# FAQ for Installation

### Python
    OITOOLS uses matplotlib through the library ```PyPlot```, which will call whatever version of Python you have installed.
    If you are starting from scratch, we recommend using Julia's Conda package for Python instead of your OS-wide python. To do so, just type:
    ```
    ENV["PYTHON"]=""  # this will tell Julia there is currently no usable python
    ```
    before installing the Conda package ```]add Conda```.

### Matplotlib backend
    You will have to select a Matplotlib backend for rendering plots. Qt5 works well in our experience. To set it up you will need to call
    ```
    ENV["MPLBACKEND"]="qt5agg"
    ```
    before using Matplotlib (or even better, before installing PyPlot).
