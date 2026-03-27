# Maximum Entropy reconstruction

OITOOLS includes a pure-Julia port of the classic **BSMEM** Maximum Entropy
algorithm (Gull, Skilling & Buscher 1991), implemented in `src/maximent.jl`
(core MaximENT engine) and `src/oimem.jl` (interferometry operators and
OITOOLS integration).

## Quick start

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

## Stopping criteria (`method`)

The `method` keyword maps to the MaximENT `methd` parameter:

| `method` | Behaviour |
|----------|-----------|
| `[1,1,1,2]` | Classic known noise — stop when χ² ≈ `aim × N_data` |
| `[4,1,1,2]` | χ² = N regardless of noise estimate (default) |
| `[2,1,1,2]` | Classic unknown noise — estimate σ from data |
| `[3,1,1,2]` | Fixed α supplied externally |

## α-update strategies

By default the Skilling evidence maximisation sets α automatically.  Two
alternatives are available:

```julia
reconstruct_bsmem(...; mackay_alpha=true)   # MacKay fixed-point update
reconstruct_bsmem(...; ritz_alpha=true)     # Ritz-value bisection
```

## Pixel size helper

When the optimal pixel size is unknown, `auto_pixsize` computes it from the
longest baseline with a given oversampling factor:

```julia
pixsize = auto_pixsize(data; oversampling=3.0)
```

See `example_image_reconstruction_oimem.jl` for a complete worked example.

## Polychromatic BSMEM

For multi-wavelength data, `reconstruct_bsmem` accepts a 3D starting image
(`nx, nx, nwav`) and per-channel data/FT vectors. Each wavelength channel gets
its own flux constraint and imaging operators; spectral smoothness is enforced
through the prior cube.

```julia
using OITOOLS
data = readoifits("data/polychromatic.oifits", filter_bad_data=true)
nwav = size(data, 1)
nx   = 64
pixsize = 0.1

# Extract per-channel data and FT plans
ft = setup_nfft(data, nx, pixsize)
data_channels = data[:, 1]   # Vector{OIdata}, one per wavelength
ft_channels   = ft[:, 1]     # Vector of FT plans

# Build a spectrally-smooth prior cube
prior_cube = gaussian_prior_cube(nx, pixsize, nwav; fwhm_mas=nx*pixsize/5)

# Reconstruct
x = reconstruct_bsmem(prior_cube, data_channels, ft_channels;
                      method   = [1, 1, 1, 2],
                      maxiter  = 100,
                      flux_err = 1e-5)
# x is an nx x nx x nwav cube, each channel normalised to unit flux
imdisp(x[:,:,1]; pixsize)
```

The polychromatic path concatenates per-channel data into one flat vector and
uses `nhid = nx*nx*nwav`, so the MaximENT engine in `maximent.jl` runs
unmodified. Entropy is computed per pixel per channel, and the spectrally-smooth
prior naturally penalizes spectral departures.

See the [Imaging (Maximum Entropy)](@ref) API reference for full docstrings.
