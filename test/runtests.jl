using OITOOLS, Test, LinearAlgebra

# _bidiag_svd calls LAPACK and reductions can reassociate with thread count, while the
# BSMEM gate is bit-for-bit. Pin threads so the comparison is reproducible.
BLAS.set_num_threads(1)

@testset "OITOOLS" begin
    include("test_bsmem_precision.jl")   # cheap structural checks first
    include("test_chromatic_params.jl")  # $WL / $MJD / $B model parameters
    include("test_bootstrap.jl")         # resampling / uncertainty estimation
    include("test_oifits_tdim.jl")       # OIFITS writing interoperability
    include("test_bsmem_regression.jl")  # the bit-for-bit numeric gate
end
