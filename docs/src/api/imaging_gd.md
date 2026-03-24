# Imaging (Gradient Descent)

VMLMB-based image reconstruction with differentiable regularizers.

| Function | Description |
|----------|-------------|
| `reconstruct(x0_4d, data, ft)` | Unified image reconstruction (mono/poly/temporal, VMLMB) |
| `reconstruct(x0_2d, data, ft)` | Monochromatic convenience wrapper |

```@docs
reconstruct(::AbstractArray{<:AbstractFloat,4}, ::AbstractMatrix{<:OITOOLS.OIdata}, ::AbstractMatrix)
reconstruct(::AbstractMatrix{<:AbstractFloat}, ::OITOOLS.OIdata, ::Any)
```
