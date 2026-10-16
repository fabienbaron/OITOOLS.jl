# BSDMM reconstruction

`reconstruct_bsdmm` uses BSDMM (Block-Simultaneous Direction Method of Multipliers, Moolekamp &
Melchior, arXiv:1708.09066): non-smooth regularizers (total variation, L2 smoothness, centering)
are applied through their proximal operators, and the χ² is minimised by VMLMB.

## Algorithm

The problem is

```
minimize_x  f(x) + sum_i  g_i(x)
```

where:

- `f(x) = chi2(H(x / sum(x)))` is the data term (V², closure phases, triple amplitudes),
  minimised by VMLMB with positivity bounds;
- `g_i(x)` are the regularizers, each with its proximal operator.

Each iteration has three steps:

1. **z-updates** (proximal):
   `z_i = prox_{g_i / rho_i}(x + u_i)` for each regularizer block.

2. **x-update** (VMLMB):
   minimise `chi2(H(x/sum(x))) + sum_i (rho_i/2) ||x - z_i + u_i||^2`
   subject to `x >= 0`.

3. **Dual updates**: `u_i <- u_i + x - z_i`.

The penalties `rho_i` can be adapted per block (see `adaptive` below).

Compared with `reconstruct`:

| | `reconstruct` (VMLMB) | `reconstruct_bsdmm` |
|---------|----------------------|---------------------|
| TV regularization | smooth approximation | exact (Chambolle proximal) |
| Positivity | VMLMB bounds | VMLMB bounds |
| Flux normalisation | chain-rule gradient | chain-rule + re-normalisation |
| Multiple regularizers | sum of smooth terms | independent proximal blocks |
| Regularizers | differentiable | any with a proximal operator |

## Quick start

```julia
using OITOOLS, PythonPlot

data = readoifits("data/Beauty_Contest_2004/2004-data1.oifits")
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

`demos/example_image_reconstruction_bsdmm.jl` compares BSDMM with `reconstruct`.

## Regularization types

| `reg_type` | Proximal | Description |
|------------|----------|-------------|
| `:tv` | `prox_tv` (Chambolle 2004) | Isotropic total variation |
| `:l2smooth` | `prox_l2smooth` (FFT) | Quadratic smoothness, `(1/2)||∇x||²` |

Both can be combined with the centering constraint (`mu_cen > 0`).

## Tuning parameters

### Regularization weight (`mu_reg`)

The weight of the spatial regularizer; for TV, typically `1e3` to `1e6`.

### ADMM penalty (`rho_init`)

The coupling between the x-update and the proximal blocks; `1e4` is a starting value.

### Inner iterations (`x_maxiter`)

VMLMB iterations per x-update. Use small values (3–10, for example 5): with many inner
iterations each x-update fits the data before the proximal steps act.

### Adaptive penalty (`adaptive`)

The `adaptive` keyword selects the penalty adaptation strategy:

| `adaptive` | Strategy | Description |
|---|---|---|
| `false` | Fixed | Constant ρ (default) |
| `true` / `:spectral` | Spectral (AADMM) | Adaptive BB selection + predictive dual + cross-block balancing |
| `:balanced` | Residual balancing | Adjusts ρ based on primal/dual residual ratio |
| `:aradmm` | ARADMM | Spectral penalty + over-relaxation parameter γ |

`adp_freq` (default 2) sets how many outer iterations pass between updates, and `adp_start`
(default 1) the first iteration updated.

### Convergence history (`history`)

When `history = true`, `reconstruct_bsdmm` returns `(image, hist)` where
`hist` is a NamedTuple with per-iteration traces: `chi2`, `obj`, `r_norm`,
`s_norm`, `rho` (Dict of per-block vectors) and `gamma`. `demos/test_adaptive_strategies.jl`
compares the strategies.

## Polychromatic reconstruction

For polychromatic data, `group_type` selects the transspectral regularizer:

| `group_type` | Regularizer | Effect |
|---|---|---|
| `:sparsity` (default) | Group sparsity (L2,1 on pixels) | Shared support across channels |
| `:tv` | Group TV (vectorial TV) | Shared edges across channels |

```julia
data = readoifits("data/V1295_Aql.oifits"; polychromatic=true)
nwav = size(data, 1)
ft   = setup_ft(data, nx, pixsize)
x0   = reshape(repeat(x_gray, 1, 1, nwav), nx, nx, nwav, 1)

x = reconstruct_bsdmm(x0, data, ft;
                       mu_tv=1e5, mu_group=1e4, mu_cen=1e2,
                       group_type=:sparsity,  # or :tv
                       rho_init=1e4, maxit=500, x_maxiter=5)
```

`data` is a `(nwav, 1)` `Matrix{OIdata}` and `x0` an `(nx, nx, nwav, 1)` cube. `mu_group = 0`
reconstructs the channels independently. See
`demos/example_image_reconstruction_bsdmm_polychromatic_v1295aql.jl`.

The docstrings are in [Imaging (BSDMM)](@ref).
