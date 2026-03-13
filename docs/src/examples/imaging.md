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

## Maximum Entropy reconstruction (OIMEM)

OITOOLS includes a pure-Julia port of the classic **BSMEM** Maximum Entropy
algorithm (Gull, Skilling & Buscher 1991), implemented in `src/maximent.jl`
(core MaximENT engine) and `src/oimem.jl` (interferometry operators and
OITOOLS integration).

### Quick start

```julia
using OITOOLS
data    = readoifits("data/2004-data1.oifits")[1,1]
nx      = 128
pixsize = 0.05   # mas/pixel
ft      = setup_nfft(data, nx, pixsize)

prior = gaussian_prior(nx, pixsize; fwhm_mas = nx * pixsize / 5)
x = reconstruct_bsmem(prior, data, ft;
                      regularizers = [["mem", prior]],
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

### Low-level interface

For more control, call `maxent_setup` and `maxent_reconstruct!` directly:

```julia
ctx, s, p = maxent_setup(data, nx, pixsize, prior; nrand=10)
image = maxent_reconstruct!(ctx, s, p; maxiter=200, verbose=true)
```

See `example_image_reconstruction_oimem.jl` for a complete worked example.

See the [API Reference](@ref) for full docstrings.
