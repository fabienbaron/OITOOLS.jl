# OITOOLS: Modeling Tools for Optical Interferometry

# Installation

Just clone the directory:
``` git clone https://github.com/fabienbaron/OITOOLS.jl.git```

You can use ```include("oitools.jl")``` to start using the OITOOLS functions.

# Usage

We provide several example files:
* example1.jl: given OIFITS data and a model image, compute the chi2, and plot the interferometric observables
* example2.jl: model OIFITS data using model-fitting (uniform disc, limb-darkened disc)
* example3.jl: model OIFITS data using a SATLAS model (open-source stellar atmosphere model code)
* example4.jl: coming soon (use Bayesian model selection to pick the best limb-dareking law)
* example5.jl: simple image reconstruction framework (mostly for teaching & experimenting purpose)
* example6.jl: coming soon (writing OIFITS data to file)
