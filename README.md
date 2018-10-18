# OITOOLS: Modeling Tools for Optical Interferometry

# Requirements

As of September 2018, Julia 1.0 is required.

The ```OIFITS```, ```NFFT```, ```SpecialFunctions``` and ```NearestNeighbors``` packages are required. For model fitting, ```NLopt``` (derivative-free local and global optimizers) and ```Multinest``` (model selection). ```DNest4``` is likely to replace ```Multinest``` soon.
For image reconstruction, ```OptimPackNextGen```.

From a fresh Julia 1.0 installation, use the package manager (```]``` key) then do:

``` add PyCall```

``` add PyPlot```

``` add FITSIO```

``` add Libdl```

``` add NLopt``` (optional, for model fitting)

``` add https://github.com/fabienbaron/OIFITS.jl.git```

``` add https://github.com/emmt/LazyAlgebra.jl.git``` (optional, for image reconstruction)

``` add https://github.com/emmt/OptimPackNextGen.jl.git``` (optional, for image reconstruction)

Check that everything works by doing:
``` using OITOOLS```

# Installation

Add the package:
``` add https://github.com/fabienbaron/OITOOLS.jl.git```

# Demos
We provide several demo files in the demos/ subfolder
* example1.jl: given OIFITS data and a model image, compute the chi2, and plot the interferometric observables
* example2.jl: model OIFITS data using model-fitting (uniform disc, limb-darkened disc)
* example3.jl: model OIFITS data using a SATLAS model (open-source stellar atmosphere model code)
* example4.jl: use Bayesian model selection to pick the best limb-darkening law
* example5.jl: simple image reconstruction framework (mostly for teaching & experimenting purpose)
* example6.jl: writing OIFITS data to file
* example_image_reconstruction_*.jl: various examples of image reconstruction
