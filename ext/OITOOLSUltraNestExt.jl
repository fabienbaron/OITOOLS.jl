# The UltraNest backend of `fit_model_nested`.
#
# Loaded automatically when the caller has PythonCall. It is PythonCall, not `ultranest`, that
# triggers this: extensions can only be keyed on Julia packages, and a Python package is not
# one. The Python side is checked when the sampler actually runs, which is also where a missing
# conda environment produces a message worth reading.
#
# `src/fit_model_ultranest.jl` holds the implementation and knows nothing about backend
# dispatch; the adapter below is the whole connection. That split is deliberate — the
# implementation file reads as "how UltraNest is driven" and nothing else.
module OITOOLSUltraNestExt

using OITOOLS
using PythonCall, Printf

# Internals `src/fit_model_ultranest.jl` reaches for. Listed rather than glob-imported so that
# moving one of them shows up here as a load error instead of a runtime `UndefVarError` in the
# middle of a fit.
using OITOOLS: OIdata, FlatModel, NestedResult,
               _batch_loglike, _bounds_vec, _ndof, chi2_flat, parse_model,
               parse_constraints, augment_constraints!, constraint_penalty,
               plot_ultranest_corner

include(joinpath(pkgdir(OITOOLS), "src", "fit_model_ultranest.jl"))

OITOOLS._fit_nested(::Val{:ultranest}, args...; kwargs...) =
    _ultranest_fit(args...; kwargs...)

__init__() = OITOOLS._register_nested_backend!(:ultranest)

end # module
