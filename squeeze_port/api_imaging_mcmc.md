# Imaging (MCMC)

SQUEEZE-style image reconstruction by Markov-chain Monte Carlo, in two flavours:
simulated annealing ([`reconstruct_squeeze`](@ref)) and parallel tempering
(`reconstruct_squeeze_tempered`, which also returns a log evidence). See the
[SQUEEZE (MCMC)](../examples/squeeze.md) guide for the method and worked
examples.

Unlike the gradient reconstructions, the image is represented as a bag of
`nelements` discrete flux quanta on the pixel grid. Positivity and total flux are
exact by construction, non-convex regularizers such as L0 are available, and the
result is a posterior mean image rather than a MAP point.

| Function | Description |
|----------|-------------|
| [`reconstruct_squeeze`](@ref)`(x_start, data, ft)` | MCMC image reconstruction by simulated annealing |
| [`SqueezeSparco`](@ref)`(; f_star, ud, …)` | SPARCO chromatic star + background model |
| `reconstruct_squeeze_tempered(data, ft)` | Parallel-tempered version (needs `using Pigeons`) |
| [`default_nelements`](@ref)`(data, nx)` | Default number of flux quanta for a dataset |

Both drivers accept `monitor = n` to redraw the image and a χ²r / parameter trace
panel every `n` iterations, and `print_every = n` for the C-format text
diagnostic line. Both are off by default.

## Regularizers

Passed as `regularizers = [["l0", λ], ["tv", λ], …]`, the same convention as
[`reconstruct`](@ref).

| Name | Form on the integer histogram | Typical λ |
|------|-------------------------------|-----------|
| `"l0"` | count of nonzero pixels | 1 – 20 |
| `"tv"` | `Σ √(dx²+dy²) / nelements`, backward differences | 200 – 2000 |
| `"entropy"` | `Σ_{x>0} lgamma(x)` (quantised / Poisson prior) | 1 – 10 |
| `"compactness"` | `Σ r²·x²` | 20 – 500 |
| `"centering"` | centroid / second-moment penalty | 1 (auto) |
| `"priorimage"` | `Σ -log(prior[pixel])` over quanta; see `prior_image` | 1 (auto) |

!!! warning "These are not OITOOLS' differentiable regularizers"
    These are SQUEEZE's forms, evaluated on the **integer histogram** with a
    `/nelements` normalisation. They are not the gradient-oriented versions in
    `oichi2.jl` (`tv`, `entropy`, `compactness`), which are defined on the
    normalised pixel image. The definitions differ, so the useful λ ranges differ
    too — by up to three orders of magnitude between regularizers, as the table
    shows. Do not carry a λ over from [`reconstruct`](@ref).

## API

```@docs
reconstruct_squeeze
reconstruct_squeeze_tempered
SqueezeSparco
default_nelements
```

!!! note "Pigeons is a weak dependency"
    `reconstruct_squeeze_tempered` lives in a package extension and only becomes
    available once you `using Pigeons`. It is not installed with OITOOLS, and CI
    does not exercise it — Pigeons pulls in MPI.jl, which is not worth building
    across every CI configuration for an opt-in feature.
