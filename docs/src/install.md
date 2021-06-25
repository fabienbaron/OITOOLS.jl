# Requirement: a working matplotlib

OITOOLS uses matplotlib through the library ```PyPlot```.
If you know you have a working PyPlot installation, just skip this section.

Otherwise, or if you are starting from a fresh julia install, we recommend using the internal julia library instead of the OS-wide python matplotlib. To do so, just type:
```
ENV["PYTHON"]=""
ENV["MPLBACKEND"]="qt5agg"
```
before proceeding with the rest of the package installation.

# Installation of Dependencies in an Environment

Simply clone OITOOLS using ```git clone```, cd to the project directory and call

(v1.6) pkg> activate .

(OITOOLS) pkg> instantiate

This will install the packages in the same state that is given by the OITOOLS Manifest.toml. Otherwise, it will resolve the latest versions of the dependencies compatible with the project.

# Global Installation of Dependencies + OITOOLS

In case the previous instantiation does not work, or if you want to install all dependencies globally, you can use the package manager (```]``` key) then do:

``` add https://github.com/fabienbaron/OIFITS.jl#t4 https://github.com/emmt/ArrayTools.jl.git https://github.com/emmt/LazyAlgebra.jl.git https://github.com/emmt/OptimPackNextGen.jl.git CFITSIO Dates DelimitedFiles Documenter DocumenterTools FITSIO LaTeXStrings LinearAlgebra NFFT NLopt UltraNest LsqFit NearestNeighbors PyCall PyPlot Random SparseArrays SpecialFunctions Statistics https://github.com/fabienbaron/OITOOLS.jl.git```

!!! info

    OITOOLS uses the ```OIFITS```, ```NFFT```, ```SpecialFunctions``` and ```NearestNeighbors``` packages. For model fitting, ```NLopt``` (derivative-free local and global optimizers), ```LsqFit``` (Levenberg-Marquardt with local error estimation) ```UltraNest``` (MCMC, model selection and best error estimates) are used. For image reconstruction, we use ```OptimPackNextGen``` written by Éric Thiébaut.


# Extra steps

If you're going to use OITOOLS to simulate observations, you may need to look up targets with SIMBAD. OITOOLS already includes some code to call Astroquery. To use it you need to add astroquery to the python installation used by Julia.

```    
using Conda
Conda.add("astropy", channel="astropy")
Conda.add("astroquery", channel="astropy")
```
