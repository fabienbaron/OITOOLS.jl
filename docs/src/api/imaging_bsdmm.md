# Imaging (BSDMM)

BSDMM (Block-Sparse Direction Method of Multipliers) reconstruction with
proximal regularizers.

| Function | Description |
|----------|-------------|
| `reconstruct_bsdmm(x0, data, ft)` | BSDMM image reconstruction with proximal regularizers |
| `prox_tv(y, lambda)` | Proximal operator for isotropic total variation (Chambolle 2004) |
| `prox_l2smooth(y, lambda)` | Proximal operator for L2 smoothness (FFT-based) |
| `prox_centering!(x, y, lambda, p, q, AtA)` | Proximal operator for centering (Woodbury formula) |
| `tv_norm(x_img)` | Isotropic total variation of a 2D image |
| `l2smooth_norm(x_img)` | L2 smoothness norm of a 2D image |
| [`prox_group_sparsity`](@ref)`(y, lambda)` | Proximal operator for L2,1 group sparsity (shared support across channels) |
| [`prox_grouptv`](@ref)`(y, lambda)` | Proximal operator for group TV (shared edges across channels) |
| `group_sparsity_norm(x_cube)` | L2,1 mixed norm of a 3D image cube |
| `grouptv_norm(x_cube)` | Group total variation of a 3D image cube |

```@docs
reconstruct_bsdmm
prox_tv
prox_l2smooth
prox_centering!
tv_norm
l2smooth_norm
prox_group_sparsity
prox_grouptv
group_sparsity_norm
grouptv_norm
```
