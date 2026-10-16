# Imaging (+ Modeling)

A parametric model for the compact components plus a reconstructed image for the extended
emission.

| Function | Description |
|----------|-------------|
| `reconstruct_hybrid(x0, params, model, data, ft)` | Alternates image VMLMB and parameter Nelder-Mead |
| `reconstruct_sparco_flat(x0, params, model, data, ft)` | Deprecated; use `reconstruct_hybrid` |
| `model_and_image_to_chi2_fg(model, params, image, g_image, ft, data)` | Hybrid model+image chi² + gradient (mono or poly) |
| `cvis_to_chi2_f(V, data)` | Chi² from complex visibilities (value only) |
| `cvis_to_chi2_noalloc(V, data)` | Chi² from complex visibilities, without allocation |
| `cvis_to_chi2_fg(V, data)` | Chi² and its gradient with respect to the complex visibilities |

```@docs
reconstruct_hybrid
reconstruct_sparco_flat
model_and_image_to_chi2_fg
cvis_to_chi2_f
cvis_to_chi2_fg
cvis_to_chi2_noalloc
```
