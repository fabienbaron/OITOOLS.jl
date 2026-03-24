# Imaging (Maximum Entropy)

BSMEM (BiSpectrum Maximum Entropy Method) reconstruction.
Works with a single `OIdata` -- use `data[1]` and `ft[1]`.

| Function | Description |
|----------|-------------|
| `reconstruct_bsmem(x0, data, ft)` | BSMEM image reconstruction |
| `auto_pixsize(data)` | Estimate pixel size from data |
| `gaussian_prior(nx, fwhm)` | Generate a Gaussian prior image |

```@docs
reconstruct_bsmem
auto_pixsize
gaussian_prior
```
