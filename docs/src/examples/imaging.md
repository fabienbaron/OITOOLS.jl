# Image reconstruction

OITOOLS implements gradient-based image reconstruction using the VMLMB optimizer
from [OptimPackNextGen.jl](https://github.com/emmt/OptimPackNextGen.jl). The
forward model uses the NFFT for speed or the exact DFT for accuracy.

## Basic monochromatic reconstruction

```julia
using OITOOLS, PythonPlot
data = readoifits("data/2004-data1.oifits")

nx      = 64       # image size (pixels per side)
pixsize = 0.1      # pixel scale (mas)
ft      = setup_ft(data, nx, pixsize)

x0 = gaussian2d(nx, nx, nx/6)   # Gaussian starting image
regularizers = [["centering", 1e4], ["l1l2", 7e6, 1e-3]]
x  = reconstruct(x0, data, ft; regularizers, maxiter=500)
imdisp(x; pixsize)
```

See `example_image_reconstruction_basic.jl`.

## Regularizations

Several regularizations can be combined in the `regularizers` list:

| Name | Syntax | Description |
|------|--------|-------------|
| `"centering"` | `["centering", μ]` | Soft centroid constraint |
| `"l1l2"` | `["l1l2", μ, α]` | Edge-preserving smoothness (L1–L2 prior with threshold α) |
| `"l1l2w"` | `["l1l2w", μ]` | Weighted L1–L2 smoothness |
| `"l1hyp"` | `["l1hyp", μ]` | L1 hyperbolic smoothness |
| `"tvsq"` | `["tvsq", μ]` | Total variation squared |
| `"tv"` | `["tv", μ]` | Total variation |
| `"l2sq"` | `["l2sq", μ]` | Quadratic smoothness |
| `"compactness"` | `["compactness", μ]` or `["compactness", μ, w]` | Compactness with optional weight map |
| `"radialvar"` | `["radialvar", μ, H, G]` | Radial variance (requires precomputed H, G from `setup_radial_reg`) |
| `"entropy"` | `["entropy", μ]` | Maximum entropy |
| `"support"` | `["support", μ, mask]` | Support constraint (penalise flux outside mask) |

Use `example_image_reconstruction_lcurve.jl` to find a good regularization weight
via the L-curve method, and `example_image_reconstruction_regularization_types.jl`
to compare regularization strategies side by side.

## Polychromatic reconstruction

When the data span multiple spectral channels, pass `polychromatic=true` to
`readoifits` to get one bin per channel:

```julia
data = readoifits("data/MWC480.oifits"; polychromatic=true)  # N×1 Matrix{OIdata}
ft   = setup_ft(data, nx, pixsize)
nwav = size(data, 1)
x0   = reshape(repeat(x0_gray, 1, 1, nwav), nx, nx, nwav, 1)  # 4D image cube
x    = reconstruct(x0, data, ft; regularizers, transspectral_regularizers)
imdisp_multi(x[:,:,:,1]; pixsize)
```

See `example_image_reconstruction_polychromatic_MWC480.jl`.

## Time-variable reconstruction

For a sequence of epochs loaded with `readoifits_multiepochs`, use temporal
regularization to link successive frames:

```julia
data = readoifits_multiepochs(files)   # 1×M Matrix{OIdata}
ft   = setup_ft(data, nx, pixsize)
nepochs = size(data, 2)
x0   = reshape(repeat(x0_gray, 1, 1, nepochs), nx, nx, 1, nepochs)  # 4D
x    = reconstruct(x0, data, ft; regularizers, temporal_regularizers)
imdisp_multi(x[:,:,1,:]; labels=["Epoch $i" for i in 1:nepochs], pixsize)
```

See the [Imaging (Gradient Descent)](@ref) and [Imaging (Setup)](@ref) API
references for full docstrings.
