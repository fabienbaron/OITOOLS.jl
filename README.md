# OITOOLS: Modeling Tools for Optical Interferometry

# Requirements

As of September 2018, Julia 1.0 is required.

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


OITOOLS uses the ```OIFITS```, ```NFFT```, ```SpecialFunctions``` and ```NearestNeighbors``` packages. For model fitting, ```NLopt``` (derivative-free local and global optimizers) and ```Multinest``` (model selection) are used. ```DNest4``` is likely to replace ```Multinest``` soon. For image reconstruction, we use ```OptimPackNextGen```.

# Installation

Add the package:
``` add https://github.com/fabienbaron/OITOOLS.jl.git```

# Demos
We provide several demo files in the demos/ subfolder:
* example_image_and_model.jl    : given OIFITS data and a model image, compute the chi2, and plot the interferometric observables
* example_limb_darkening_fit.jl : given OIFITS data, do model-fitting (uniform disc, limb-darkened disc)
* example_satlas_fit.jl         : model OIFITS data using a SATLAS model (open-source stellar atmosphere model code)
* example_nested_sampling_fit.jl: use Bayesian model selection via nested sampling to compare limb-darkening laws
* example_bootstrap_fit.jl      : use the boostrap method to estimate error bars
* example_npoi_target_filter.jl : how to select only a given target within an OIFITS, and filter bad SNR data
* example_fakedata_hourangle.jl : simulate observations from target image and Hour Angle, write OIFITS data to file
* example_fakedata_databased.jl : simulate observations from target image and already existing OIFITS
* example_image_reconstruction_dft.jl  : gradient-based image reconstruction using the exact DFT
* example_image_reconstruction_nfft.jl : gradient-based image reconstruction using fast yet accurate NFFT
* example_image_reconstruction_lcurve.jl : l-curve method to determine the regularization factor
* example_image_reconstruction_epll.jl  : gradient-based image reconstruction using machine-learned priors (work in progress)
* example_image_reconstruction_multitemporal.jl : gradient-based image reconstruction for time-variable images, with temporal regularization
* example_image_reconstruction_multiwavelength.jl : (upcoming) gradient-based image reconstruction for spectrally dependent images, with transpectral regularization
* example_oifitslib.jl                  : (upcoming) an interface to John Young's OIFITSLIB utilities (oimerge, oifilter, oicheck)

