# Installation

OITOOLS uses matplotlib through the library ```PyPlot```.
If you are starting from a fresh julia install, we recommend using the internal julia library instead of the OS-wide python matplotlib. To do so, just type:
```
ENV["PYTHON"]=""
ENV["MPLBACKEND"]="qt5agg"
```
before proceeding with the following package installation.
Use the package manager (```]``` key) then do:

``` add PyCall PyPlot LaTeXStrings CFITSIO FITSIO Libdl NLopt NFFT SpecialFunctions NearestNeighbors https://github.com/fabienbaron/OIFITS.jl#t4 https://github.com/emmt/ArrayTools.jl.git https://github.com/emmt/LazyAlgebra.jl.git https://github.com/emmt/OptimPackNextGen.jl.git https://github.com/fabienbaron/OITOOLS.jl.git```

!!! info

    OITOOLS uses the ```OIFITS```, ```NFFT```, ```SpecialFunctions``` and ```NearestNeighbors``` packages. For model fitting, ```NLopt``` (derivative-free local and global optimizers) and ```Multinest``` (model selection) are used. ```DNest4``` is likely to replace ```Multinest``` soon. For image reconstruction, we use ```OptimPackNextGen``` written by Éric Thiébaut.

To install Astroquery:
```    
using Conda
Conda.add("astropy", channel="astropy")
Conda.add("astroquery", channel="astropy")
```
