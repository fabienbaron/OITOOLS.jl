# Imaging (SPARCO)

SPARCO (SPectral Analysis of Resolved Compact Objects) decomposes the scene into
a parametric model (stars, background) + a reconstructed image (environment).

| Function | Description |
|----------|-------------|
| `reconstruct_sparco(x0, data, ft; lambda_ref, ...)` | High-level SPARCO: set physical params, get image back |
| `chi2_sparco_flat_f(x, params, model, ft, data)` | Flat-model SPARCO chi-squared (forward only) |
| `optimize_sparco_flat_parameters(params, model, x, ft, data)` | Optimize flat-model SPARCO parameters (image fixed) |
| `model_and_image_to_vis(model, params, x_img, ft, data)` | Combined visibility: `V_model + W(λ) · V_image` |

```@docs
reconstruct_sparco
chi2_sparco_flat_f
optimize_sparco_flat_parameters
model_and_image_to_vis
```
