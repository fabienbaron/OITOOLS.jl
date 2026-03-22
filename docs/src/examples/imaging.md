# Image reconstruction

OITOOLS implements gradient-based image reconstruction using the VMLMB optimizer
from [OptimPackNextGen.jl](https://github.com/emmt/OptimPackNextGen.jl). The
forward model uses the NFFT for speed or the exact DFT for accuracy.

## Basic monochromatic reconstruction

```julia
using OITOOLS
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

## SPARCO grey reconstruction

For targets with a spectrally-variable environment (e.g. T Tauri discs), the
SPARCO framework separates a grey circumstellar component from an unresolved
stellar point source. SPARCO functions work with a single `OIdata`:

```julia
x = reconstruct_sparco_gray(x0, data[1], ft[1]; regularizers)
```

See `example_image_reconstruction_sparco_grey.jl`.

## Maximum Entropy reconstruction (OIMEM)

OITOOLS includes a pure-Julia port of the classic **BSMEM** Maximum Entropy
algorithm (Gull, Skilling & Buscher 1991), implemented in `src/maximent.jl`
(core MaximENT engine) and `src/oimem.jl` (interferometry operators and
OITOOLS integration).

### Quick start

```julia
using OITOOLS
data    = readoifits("data/2004-data1.oifits")
nx      = 128
pixsize = 0.05   # mas/pixel
ft      = setup_ft(data, nx, pixsize)

prior = gaussian_prior(nx, pixsize; fwhm_mas = nx * pixsize / 5)
x = reconstruct_bsmem(prior, data, ft;
                      method       = [1, 1, 1, 2],
                      maxiter      = 100,
                      flux_err     = 1e-5)
imdisp(x; pixsize)
```

`reconstruct_bsmem` has the same interface as `reconstruct` and is the
recommended entry point.  It infers pixel size and image dimensions from the
NFFT plan `ft`, builds the imaging context and MaximENT state internally, and
returns an `nx × nx` image normalised to unit flux.

### Stopping criteria (`method`)

The `method` keyword maps to the MaximENT `methd` parameter:

| `method` | Behaviour |
|----------|-----------|
| `[1,1,1,2]` | Classic known noise — stop when χ² ≈ `aim × N_data` |
| `[4,1,1,2]` | χ² = N regardless of noise estimate (default) |
| `[2,1,1,2]` | Classic unknown noise — estimate σ from data |
| `[3,1,1,2]` | Fixed α supplied externally |

### α-update strategies

By default the Skilling evidence maximisation sets α automatically.  Two
alternatives are available:

```julia
reconstruct_bsmem(...; mackay_alpha=true)   # MacKay fixed-point update
reconstruct_bsmem(...; ritz_alpha=true)     # Ritz-value bisection
```

### Pixel size helper

When the optimal pixel size is unknown, `auto_pixsize` computes it from the
longest baseline with a given oversampling factor:

```julia
pixsize = auto_pixsize(data; oversampling=3.0)
```

See `example_image_reconstruction_oimem.jl` for a complete worked example.

See the [Imaging](@ref) API reference for full docstrings.
