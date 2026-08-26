# The NestedSamplers.jl backend of `fit_model_nested`.
#
# Pure Julia, so this is the nested sampler a PackageCompiler build can contain — PythonCall
# and its conda environment cannot go into one.
#
# Structured exactly like ext/OITOOLSUltraNestExt.jl: the implementation lives in
# src/fit_model_nestedsamplers.jl and knows nothing about backend dispatch, and the adapter
# below is the whole connection.
module OITOOLSNestedSamplersExt

using OITOOLS
using NestedSamplers
using NestedSamplers: Bounds, Proposals, Nested, NestedModel,
                      prior_transform_and_loglikelihood
using Printf, Random
using Random: AbstractRNG

# Internals src/fit_model_nestedsamplers.jl reaches for. Listed rather than glob-imported so
# that moving one shows up here as a load error rather than an `UndefVarError` mid-fit.
using OITOOLS: OIdata, FlatModel, NestedResult, ModelConstraint,
               _bounds_vec, _ndof, chi2_flat, parse_model,
               parse_constraints, augment_constraints!, constraint_penalty,
               plot_corner_makie

include(joinpath(pkgdir(OITOOLS), "src", "fit_model_nestedsamplers.jl"))

OITOOLS._fit_nested(::Val{:nestedsamplers}, args...; kwargs...) =
    _nestedsamplers_fit(args...; kwargs...)

__init__() = OITOOLS._register_nested_backend!(:nestedsamplers)

end # module
