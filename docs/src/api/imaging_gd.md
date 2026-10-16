# Imaging (Gradient Descent)

Image reconstruction by VMLMB with differentiable regularizers.

| Function | Description |
|----------|-------------|
| `reconstruct(x0_4d, data, ft)` | Image reconstruction: monochromatic, polychromatic or time series |
| `reconstruct(x0_2d, data, ft)` | Monochromatic, with a 2-D image and one `OIdata` |

```@docs
reconstruct(::AbstractArray{<:AbstractFloat,4}, ::AbstractMatrix{<:OITOOLS.OIdata}, ::AbstractMatrix)
reconstruct(::AbstractMatrix{<:AbstractFloat}, ::OITOOLS.OIdata, ::Any)
```
