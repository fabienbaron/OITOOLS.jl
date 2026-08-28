# The Nautilus.jl backend of `fit_model_nested`.
#
# Pure Julia, like `:nestedsamplers`, so this is a nested sampler a PackageCompiler build can
# contain — and unlike that one it is an importance sampler, which is where the reduction in
# χ² evaluations comes from.
#
# Structured exactly like ext/OITOOLSUltraNestExt.jl: the implementation lives in
# src/fit_model_nautilus.jl and knows nothing about backend dispatch, and the adapter below is
# the whole connection.
module OITOOLSNautilusExt

using OITOOLS
using Nautilus
using Printf, Random

# Internals src/fit_model_nautilus.jl reaches for. Listed rather than glob-imported so that
# moving one shows up here as a load error rather than an `UndefVarError` mid-fit.
using OITOOLS: OIdata, FlatModel, NestedResult, ModelConstraint,
               _bounds_vec, _ndof, chi2_flat, parse_model,
               parse_constraints, augment_constraints!, constraint_penalty,
               plot_corner_makie

include(joinpath(pkgdir(OITOOLS), "src", "fit_model_nautilus.jl"))

OITOOLS._fit_nested(::Val{:nautilus}, args...; kwargs...) =
    _nautilus_fit(args...; kwargs...)

__init__() = OITOOLS._register_nested_backend!(:nautilus)

end # module
