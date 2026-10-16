# Imaging (MCMC)

SQUEEZE image reconstruction by Markov-chain Monte Carlo: simulated annealing
([`reconstruct_squeeze`](@ref)) or parallel tempering (`reconstruct_squeeze_tempered`, which
also returns the log evidence). See [SQUEEZE (MCMC)](../examples/squeeze.md).

The image is `nelements` flux elements on the pixel grid, so it is non-negative with a fixed
total flux, non-convex regularizers such as L0 are available, and the result is a mean image
over samples.

| Function | Description |
|----------|-------------|
| [`reconstruct_squeeze`](@ref)`(x_start, data, ft)` | MCMC image reconstruction by simulated annealing |
| [`SqueezeSparco`](@ref)`(; f_star, ud, …)` | SPARCO chromatic star + background model |
| `reconstruct_squeeze_tempered(data, ft)` | Parallel-tempered version (needs `using Pigeons`) |
| [`default_nelements`](@ref)`(data, nx)` | Default number of flux elements for a dataset |

Both accept `monitor = n` (draw the image and χ²r / parameter traces every `n` iterations) and
`print_every = n` (SQUEEZE's text line); both are off by default.

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
| `"priorimage"` | `Σ -log(prior[pixel])` over elements; see `prior_image` | 1 (auto) |

!!! warning "Not the same as `reconstruct`'s regularizers"
    These are SQUEEZE's forms, on the integer histogram with a `/nelements` normalisation.
    `reconstruct`'s `tv`, `entropy` and `compactness` are defined on the normalised image, so
    their λ values do not carry over.

## API

```@docs
reconstruct_squeeze
reconstruct_squeeze_tempered
SqueezeSparco
default_nelements
```

`reconstruct_squeeze_tempered` needs `using Pigeons`, which is not installed with OITOOLS.
