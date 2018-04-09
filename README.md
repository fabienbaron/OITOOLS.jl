# OITOOLS: Modeling Tools for Optical Interferometry

# Requirements

The ```OIFITS```, ```NFFT```, ```SpecialFunctions``` and ```NearestNeighbors``` packages are required and will be automatically installed the first time running OITOOLS if they are missing.

We use two optional packages for model fitting are ```NLopt``` (derivative-free local and global optimizers) and ```Multinest``` (model selection). ```DNest4``` is likely to replace ```Multinest``` soon.

For image reconstruction, we use ```OptimPack```.

# Installation

Just clone the directory:
``` git clone https://github.com/fabienbaron/OITOOLS.jl.git```

You can use ```include("oitools.jl")``` to start using the OITOOLS functions.

# Usage

We provide several example files:
* example1.jl: given OIFITS data and a model image, compute the chi2, and plot the interferometric observables
* example2.jl: model OIFITS data using model-fitting (uniform disc, limb-darkened disc)
* example3.jl: model OIFITS data using a SATLAS model (open-source stellar atmosphere model code)
* example4.jl: use Bayesian model selection to pick the best limb-darkening law
* example5.jl: simple image reconstruction framework (mostly for teaching & experimenting purpose)
* example6.jl: writing OIFITS data to file
