# Image reconstruction

OITOOLS implements gradient-based image reconstruction using the VMLMB optimizer
from [OptimPackNextGen.jl](https://github.com/emmt/OptimPackNextGen.jl). The
forward model uses the NFFT for speed or the exact DFT for accuracy.

## Basic monochromatic reconstruction

```julia
using OITOOLS
data = readoifits("data/2004-data1.oifits")[1,1]

nx      = 64       # image size (pixels per side)
pixsize = 0.1      # pixel scale (mas)
ft      = setup_nfft(data, nx, pixsize)

x0 = gaussian2d(nx, nx, nx/6)   # Gaussian starting image
regularizers = [["centering", 1e4], ["l1l2", 7e6, 1e-3]]
x  = reconstruct(x0, data, ft; regularizers, maxiter=500)
imdisp(x; pixsize)
```

See `example_image_reconstruction_basic.jl`.

## Regularizations

Several regularizations can be combined in the `regularizers` list:

| Name | Description |
|------|-------------|
| `"centering"` | Soft centroid constraint |
| `"l1l2"` | Edge-preserving smoothness (L1–L2 prior) |
| `"tvsq"` | Total variation squared |
| `"tv"` | Total variation |
| `"l2sq"` | Quadratic smoothness |

Use `example_image_reconstruction_lcurve.jl` to find a good regularization weight
via the L-curve method, and `example_image_reconstruction_regularization_types.jl`
to compare regularization strategies side by side.

## Polychromatic reconstruction

When the data span multiple spectral channels, pass `polychromatic=true` to
`readoifits` to get one bin per channel, then use the polychromatic setup:

```julia
data = vec(readoifits("data/MWC480.oifits"; polychromatic=true))
ft   = setup_nfft_polychromatic(data, nx, pixsize)
x    = reconstruct_polychromatic(x0_cube, data, ft; regularizers)
imdisp_polychromatic(x; wavs, pixsize)
```

See `example_image_reconstruction_polychromatic_MWC480.jl`.

## Time-variable reconstruction

For a sequence of epochs loaded with `readoifits_multiepochs`, use temporal
regularization to link successive frames:

```julia
nepochs, tepochs, data = readoifits_multiepochs(files)
ft = setup_nfft_multiepochs(data, nx, pixsize)
x  = reconstruct_multitemporal(x0_cube, data, ft; regularizers)
imdisp_temporal(x, nepochs; pixsize)
```

## SPARCO grey reconstruction

For targets with a spectrally-variable environment (e.g. T Tauri discs), the
SPARCO framework separates a grey circumstellar component from an unresolved
stellar point source:

```julia
x = reconstruct_sparco_gray(x0, data, ft; regularizers)
```

See `example_image_reconstruction_sparco_grey.jl`.
