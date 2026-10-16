# Image reconstruction

`reconstruct` minimises χ² plus regularization with the VMLMB optimizer of
[OptimPackNextGen.jl](https://github.com/emmt/OptimPackNextGen.jl). The Fourier transform is an
NFFT (default) or a DFT (`setup_ft(…; mode = "dft")`).

## Basic monochromatic reconstruction

```julia
using OITOOLS, PythonPlot
data = readoifits("data/Beauty_Contest_2004/2004-data1.oifits")

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

Regularizers are combined in the `regularizers` list:

| Name | Syntax | Description |
|------|--------|-------------|
| `"centering"` | `["centering", μ]` | Soft centroid constraint |
| `"l1l2"` | `["l1l2", μ, α]` | Edge-preserving smoothness (L1–L2 prior with threshold α) |
| `"l1l2w"` | `["l1l2w", μ]` | `Σ (x − log(1 + x))` on the pixels |
| `"l1hyp"` | `["l1hyp", μ]` | Hyperbolic ℓ1 norm of the pixels (sparsity) |
| `"tvsq"` | `["tvsq", μ]` | Squared gradient norm (quadratic smoothness) |
| `"tv"` | `["tv", μ]` or `["tv", μ, ε]` | Total variation, with relaxation ε (default 1e-8) |
| `"l2sq"` | `["l2sq", μ]` | `Σ x²` (Tikhonov) |
| `"compactness"` | `["compactness", μ]` or `["compactness", μ, w]` | Compactness with optional weight map |
| `"radialvar"` | `["radialvar", μ, H, G]` | Radial variance (requires precomputed H, G from `setup_radial_reg`) |
| `"entropy"` | `["entropy", μ]` | `Σ (x log x − x)` |
| `"support"` | `["support", μ, mask]` | `Σ x²` over the pixels where `mask ≤ 0` |

`example_image_reconstruction_lcurve.jl` chooses a weight by the L-curve
([`lcurve`](@ref)), and `example_image_reconstruction_regularization_types.jl` compares
regularizers.

## Polychromatic reconstruction

`readoifits(…; polychromatic=true)` gives one bin per spectral channel:

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

For epochs read with `readoifits_multiepochs`, temporal regularizers link successive frames:

```julia
data = readoifits_multiepochs(files)   # 1×M Matrix{OIdata}
ft   = setup_ft(data, nx, pixsize)
nepochs = size(data, 2)
x0   = reshape(repeat(x0_gray, 1, 1, nepochs), nx, nx, 1, nepochs)  # 4D
x    = reconstruct(x0, data, ft; regularizers, temporal_regularizers)
imdisp_multi(x[:,:,1,:]; labels=["Epoch $i" for i in 1:nepochs], pixsize)
```

The docstrings are in [Imaging (Gradient Descent)](@ref) and [Imaging (Setup)](@ref).
