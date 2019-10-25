# Requirements

As of March 2018, Julia 1.1 is required.

From a fresh Julia 1.1 installation, use the package manager (```]``` key) then do:

``` add PyCall PyPlot LaTeXStrings FITSIO Libdl NLopt NFFT SpecialFunctions NearestNeighbors https://github.com/fabienbaron/OIFITS.jl#t4 https://github.com/emmt/LazyAlgebra.jl.git https://github.com/emmt/OptimPackNextGen.jl.git```

OITOOLS uses the ```OIFITS```, ```NFFT```, ```SpecialFunctions``` and ```NearestNeighbors``` packages. For model fitting, ```NLopt``` (derivative-free local and global optimizers) and ```Multinest``` (model selection) are used. ```DNest4``` is likely to replace ```Multinest``` soon. For image reconstruction, we use ```OptimPackNextGen```.

# Installation

Add the package:
``` add https://github.com/fabienbaron/OITOOLS.jl.git```

Check that everything works by doing:
``` using OITOOLS```
