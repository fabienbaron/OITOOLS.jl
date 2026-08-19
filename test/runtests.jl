using OITOOLS, Test, LinearAlgebra

# Loads the OITOOLSPythonPlotExt extension. Plotting is no longer in the core package, so
# without this every figure test would call a stub function with no methods. MPLBACKEND is
# set first because PythonPlot picks its backend at import time, and the default probe
# loads Qt.
ENV["MPLBACKEND"] = "Agg"
using PythonPlot

# _bidiag_svd calls LAPACK and reductions can reassociate with thread count, while the
# BSMEM gate is bit-for-bit. Pin threads so the comparison is reproducible.
BLAS.set_num_threads(1)

@testset "OITOOLS" begin
    include("test_bsmem_precision.jl")   # cheap structural checks first
    include("test_chromatic_params.jl")  # $WL / $MJD / $B model parameters
    include("test_model_gradients.jl")   # analytic model gradients vs finite differences
    include("test_vis_functions.jl")     # legacy (param, uv) visibility function API
    include("test_bootstrap.jl")         # resampling / uncertainty estimation
    include("test_oifits_tdim.jl")       # OIFITS writing interoperability
    include("test_simulate.jl")          # simulate(): geometry, noise model, observability
    include("test_planning.jl")          # astrometry, twilight, observability conventions
    include("test_python_boundary.jl")   # Julia<->Python crossings (plots, UltraNest, SIMBAD)
    include("test_plotting.jl")          # every figure: renders, plotted values, options
    include("test_squeeze.jl")           # SQUEEZE MCMC sampler
    include("test_squeeze_tempering.jl") # SQUEEZE + Pigeons (skipped if absent)
    include("test_bsmem_regression.jl")  # the bit-for-bit numeric gate
end
