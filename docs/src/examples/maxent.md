# Maximum Entropy reconstruction

`reconstruct_bsmem` is a Julia port of BSMEM, maximum entropy reconstruction with the MaximENT
algorithm (Gull, Charter & Skilling). The MaximENT engine is in `src/maximent.jl`, the
interferometry operators in `src/oimem.jl`.

## Quick start

```julia
using OITOOLS, PythonPlot
data    = readoifits("data/Beauty_Contest_2004/2004-data1.oifits")
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

`reconstruct_bsmem` has the same arguments as `reconstruct`. The pixel size and image size come
from `ft`; the result is an `nx × nx` image normalised to unit flux.

## Stopping criteria (`method`)

The `method` keyword maps to the MaximENT `methd` parameter:

| `method` | Behaviour |
|----------|-----------|
| `[1,1,1,2]` | Classic known noise — stop when χ² ≈ `aim × N_data` |
| `[4,1,1,2]` | χ² = N regardless of noise estimate (default) |
| `[2,1,1,2]` | Classic unknown noise — estimate σ from data |
| `[3,1,1,2]` | Fixed α supplied externally |

## α-update strategies

By default α is set by Skilling's evidence maximisation. The alternative is MacKay's
fixed-point update:

```julia
reconstruct_bsmem(...; mackay_alpha=true)   # MacKay fixed-point update
```

It uses ω directly, so the `Good` estimate (about 500× smaller than the exact trace in these
runs) cancels.

## Pixel size helper

`auto_pixsize` gives a pixel size from the longest baseline and an oversampling factor:

```julia
pixsize = auto_pixsize(data; oversampling=3.0)
```

See `example_image_reconstruction_oimem.jl` for a complete worked example.

## Polychromatic BSMEM

For several channels, `reconstruct_bsmem` takes a 3-D starting image `(nx, nx, nwav)` and
per-channel vectors of data and plans. Each channel has its own flux constraint and operators;
the prior cube links the channels.

```julia
using OITOOLS, PythonPlot
data = readoifits("data/polychromatic.oifits"; polychromatic=true, filter_bad_data=true)
nwav = size(data, 1)
nx   = 64
pixsize = 0.1

ft = setup_ft(data, nx, pixsize)
data_channels = data[:, 1]   # Vector{OIdata}, one per channel
ft_channels   = ft[:, 1]     # one plan per channel

prior_cube = gaussian_prior_cube(nx, pixsize, nwav; fwhm_mas=nx*pixsize/5)

x = reconstruct_bsmem(prior_cube, data_channels, ft_channels;
                      method   = [1, 1, 1, 2],
                      maxiter  = 100,
                      flux_err = 1e-5)
# x is an nx x nx x nwav cube, each channel normalised to unit flux
imdisp(x[:,:,1]; pixsize)
```

The channels' data are concatenated into one vector and the image has `nx*nx*nwav` pixels; the
entropy is per pixel and channel, relative to the prior cube.

The docstrings are in [Imaging (Maximum Entropy)](@ref).
