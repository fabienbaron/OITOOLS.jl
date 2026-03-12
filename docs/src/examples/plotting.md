# Plotting OIFITS data

All plot functions accept either a single `OIdata` or an array `AbstractArray{<:OIdata}`
(e.g., the 2D output of `readoifits`, a polychromatic array, or a multi-epoch vector).
Call `set_oiplot_defaults()` to apply consistent matplotlib style settings.

## UV coverage

```julia
uvplot(data)                       # colour by baseline length
uvplot(data; color="wav")          # colour by wavelength
uvplot(data; color="mjd")          # colour by MJD
uvplot(data; color="baseline", flipx=true, square=true, figsize=(10,10))
```

Key keyword arguments: `color`, `figsize`, `minuv`, `maxuv`, `square`, `legend_below`,
`figtitle`, `cmap`, `flipx`.

## Squared visibilities

```julia
plot_v2(data)                      # linear scale
plot_v2(data; logplot=true)        # log scale
plot_v2(data; color="wav")         # colour by wavelength
plot_v2(data; color="mjd")         # colour by timestamp
```

Key arguments: `figsize`, `logplot`, `color`, `markopt`, `legend_below`, `figtitle`.

## Closure phases

```julia
plot_t3phi(data)                   # x-axis = shortest baseline in triangle
plot_t3phi(data; t3base="max")     # x-axis = longest baseline
plot_t3phi(data; t3base="geom")    # x-axis = geometric mean baseline
```

## Triple amplitudes

```julia
plot_t3amp(data)
plot_t3amp(data; t3base="max")
```

## Flux spectra

```julia
plot_flux(data)                    # colour by station; "Calibrated" for CALSTAT=C data
plot_flux(data; color="mjd")       # scatter coloured by MJD with colorbar
plot_flux(data; color="red")       # fixed colour for all points
```

`plot_flux` maps `flux_sta_index == 0` (calibrated OI_FLUX) to the label
`"Calibrated"` and uses the station name otherwise.

Key arguments: `figsize`, `color`, `markopt`, `legend_below`.

## Differential / complex visibilities

```julia
plot_visphi(data)                  # visibility phases
plot_diffphi(data)                 # differential phases (vector of OIdata)
```

## Combined V² + closure-phase panel

```julia
plot_v2_and_t3phi_wav(data)
plot_v2_and_t3phi_wav(data; logplot=true, figsize=(14,14))
```

## Residuals vs model

Pass a precomputed model vector alongside the data:

```julia
plot_v2_residuals(data, v2_model)
plot_t3phi_residuals(data, t3phi_model)
plot_t3amp_residuals(data, t3amp_model)
```

## Multi-file V² comparison

```julia
files_data = [data1, data2, data3]          # Vector{OIdata}
plot_v2_multifile(files_data)
```

## Image display

```julia
imdisp(image; pixsize=0.5, colormap="gist_heat", use_colorbar=true)
imdisp_polychromatic(image_cube; wavs=wavs, pixsize=0.5)
imdisp_temporal(image_cube, nepochs; pixsize=0.5)
```

`pixsize` is in mas/pixel. `beamsize` (mas) draws a PSF indicator.

## Facility / array layout

```julia
plot_facility(facility)
```
