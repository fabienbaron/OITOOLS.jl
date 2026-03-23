# BSDMM reconstruction

OITOOLS includes a **Block-Structured SDMM (BSDMM)** image reconstruction
algorithm that handles non-smooth regularizers (total variation, L2
smoothness) via proximal operators, while keeping the non-convex data
fidelity term (chi-squared) in a VMLMB-based x-update.

## Algorithm overview

Standard gradient-based reconstruction (`reconstruct`) minimises a
composite objective directly with VMLMB, which requires differentiable
regularizers.  BSDMM instead splits the problem using the alternating
direction method of multipliers (ADMM):

```
minimize_x  f(x) + sum_i  g_i(x)
```

where:

- **f(x) = chi2(H(x / sum(x)))** is the data fidelity (squared visibilities,
  closure phases, triple amplitudes).  It is non-convex and handled by
  VMLMB with positivity bounds and chain-rule flux normalisation.

- **g_i(x)** are non-smooth regularizers (TV, L2 smoothness, centering),
  each handled by its own proximal operator.

Each BSDMM iteration performs three steps:

1. **z-updates** (proximal, parallelisable):
   `z_i = prox_{g_i / rho_i}(x + u_i)` for each regularizer block.

2. **x-update** (VMLMB):
   minimise `chi2(H(x/sum(x))) + sum_i (rho_i/2) ||x - z_i + u_i||^2`
   subject to `x >= 0`.

3. **Dual updates**: `u_i <- u_i + x - z_i`.

The penalty parameters `rho_i` can optionally be adapted per block using
Barzilai-Borwein spectral estimates (ARADMM-style).

### Why BSDMM?

| Feature | `reconstruct` (VMLMB) | `reconstruct_bsdmm` |
|---------|----------------------|---------------------|
| TV regularization | smooth approximation | exact (Chambolle proximal) |
| Positivity | VMLMB bounds | VMLMB bounds |
| Flux normalisation | chain-rule gradient | chain-rule + re-normalisation |
| Multiple regularizers | sum of smooth terms | independent proximal blocks |
| Extensibility | requires differentiable reg | any proximal operator |

BSDMM is particularly useful when you want exact (non-smooth) total
variation or plan to add custom non-differentiable constraints.

## Quick start

```julia
using OITOOLS

data = readoifits("data/2004-data1.oifits")
nx      = 64
pixsize = 0.2
ft      = setup_ft(data, nx, pixsize)

x0 = gaussian2d(nx, nx, nx / 6)

x = reconstruct_bsdmm(x0, data, ft;
                       mu_reg   = 1e5,        # TV regularization weight
                       mu_cen   = 1e2,        # centering weight
                       reg_type = :tv,         # :tv or :l2smooth
                       rho_init = 1e4,         # ADMM penalty
                       adaptive = false,       # fixed rho
                       maxit    = 500,         # outer BSDMM iterations
                       x_maxiter = 5)          # VMLMB steps per x-update
imdisp(x; pixsize)
```

See `example_image_reconstruction_bsdmm.jl` for a complete script that
compares BSDMM with the standard `reconstruct`.

## Regularization types

| `reg_type` | Proximal | Description |
|------------|----------|-------------|
| `:tv` | `prox_tv` (Chambolle 2004) | Isotropic total variation — edge-preserving, non-smooth |
| `:l2smooth` | `prox_l2smooth` (FFT) | Quadratic smoothness — `(1/2)||nabla x||^2`, smooth |

Both can be combined with the centering constraint (`mu_cen > 0`).

## Tuning parameters

### Regularization weight (`mu_reg`)

Controls the strength of spatial regularization.  Larger values produce
smoother images; smaller values allow more detail but risk overfitting.
Typical range: `1e3` to `1e6` for TV.

### ADMM penalty (`rho_init`)

Controls the coupling between the x-update and the proximal blocks.
A good starting point is `rho_init = 1e4`.

### Inner iterations (`x_maxiter`)

The number of VMLMB iterations per x-update.  **Use small values**
(3-10).  With too many inner iterations, the x-update overfits the data
at each step, and the regularization proximal cannot keep up.
Recommended: `x_maxiter = 5`.

### Adaptive penalty (`adaptive`)

When `adaptive = true`, each block's penalty `rho_i` is updated using
Barzilai-Borwein spectral estimates with correlation safeguards.  This
can improve convergence but may also cause instability for some
problems.  Start with `adaptive = false`.

## API reference

See the [Imaging](@ref) API reference for full docstrings of
`reconstruct_bsdmm`, proximal operators, and regularization norms.
