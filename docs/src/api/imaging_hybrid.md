# Imaging (+ Modeling)

Hybrid model + image reconstruction: a parametric model handles compact/unresolved
components while an image captures the extended emission.

| Function | Description |
|----------|-------------|
| `reconstruct_hybrid(x0, params, model, data, ft)` | Alternating image VMLMB + parameter NelderMead |
| `reconstruct_sparco_flat(x0, params, model, data, ft)` | Joint VMLMB over params + image (legacy) |
| `model_and_image_to_chi2_fg(model, params, image, g_image, ft, data)` | Hybrid model+image chi² + gradient (mono or poly) |
| `cvis_to_chi2_f(V, data)` | Chi² from complex visibilities (value only) |
| `cvis_to_chi2_fg(V, data)` | Chi² + complex adjoint source from complex visibilities (building block) |

```@docs
reconstruct_hybrid
reconstruct_sparco_flat
model_and_image_to_chi2_fg
cvis_to_chi2_f
cvis_to_chi2_fg
```
