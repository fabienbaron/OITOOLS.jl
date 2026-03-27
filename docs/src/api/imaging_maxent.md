# Imaging (Maximum Entropy)

BSMEM (BiSpectrum Maximum Entropy Method) reconstruction.

You can pass `data` and `ft` directly from `readoifits` / `setup_ft` —
`reconstruct_bsmem` auto-dispatches based on the input dimensions.
For explicit control, pass a single `OIdata` (mono) or a 3D starting image
with `Vector{OIdata}` (polychromatic).

| Function | Description |
|----------|-------------|
| `reconstruct_bsmem(x0, data, ft)` | BSMEM image reconstruction (mono or polychromatic) |
| `auto_pixsize(data)` | Estimate pixel size from data |
| `gaussian_prior(nx, fwhm)` | Generate a 2D Gaussian prior image |
| `gaussian_prior_cube(nx, pixsize, nwav)` | Generate a 3D Gaussian prior cube for polychromatic BSMEM |

```@docs
reconstruct_bsmem
auto_pixsize
gaussian_prior
gaussian_prior_cube
```
