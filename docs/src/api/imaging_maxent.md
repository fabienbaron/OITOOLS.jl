# Imaging (Maximum Entropy)

BSMEM (BiSpectrum Maximum Entropy Method) reconstruction.

`reconstruct_bsmem` takes `data` and `ft` from `readoifits` / `setup_ft` and chooses the mode
from their dimensions; it also takes a single `OIdata` (monochromatic) or a 3-D starting image
with a `Vector{OIdata}` (polychromatic).

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
