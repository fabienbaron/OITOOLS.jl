# OITOOLS: Modeling Tools for Optical Interferometry

# Requirements

The ```OIFITS```, ```NFFT```, ```SpecialFunctions``` and ```NearestNeighbors``` packages are required. For model fitting, ```NLopt``` (derivative-free local and global optimizers) and ```Multinest``` (model selection). ```DNest4``` is likely to replace ```Multinest``` soon.
For image reconstruction, ```OptimPackNextGen```.

From a fresh Julia 0.7 installation, use the package manager (```]``` key) then do:

``` add FITSIO#master``` (as of September 2018)

``` add https://github.com/fabienbaron/OIFITS.jl.git```

``` add https://github.com/emmt/LazyAlgebra.jl.git```

``` add https://github.com/fabienbaron/OptimPackNextGen.jl.git```

# Installation

Add the package:
``` add https://github.com/fabienbaron/OITOOLS.jl.git```

# Demos

We provide several demo files:
* example1.jl: given OIFITS data and a model image, compute the chi2, and plot the interferometric observables
* example2.jl: model OIFITS data using model-fitting (uniform disc, limb-darkened disc)
* example3.jl: model OIFITS data using a SATLAS model (open-source stellar atmosphere model code)
* example4.jl: use Bayesian model selection to pick the best limb-darkening law
* example5.jl: simple image reconstruction framework (mostly for teaching & experimenting purpose)
* example6.jl: writing OIFITS data to file
* example_image_reconstruction_*.jl: various examples of image reconstruction
