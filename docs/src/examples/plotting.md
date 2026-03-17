# Plotting OIFITS data

All plot functions accept either a single `OIdata` or an array `AbstractArray{<:OIdata}`
(e.g., the 2D output of `readoifits`, a polychromatic array, or a multi-epoch vector).
Call `set_oiplot_defaults()` to apply consistent matplotlib style settings.
Use `set_oiplot_defaults(compact=true)` for smaller plots suitable for stacking
in a multi-panel figure (smaller fonts, markers, error bars, and no in-plot title).

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
plot_t3phi(data)                   # x-axis = geometric mean baseline (default)
plot_t3phi(data; t3base="max")     # x-axis = longest baseline
plot_t3phi(data; color="mjd")      # colour by MJD
```

Key arguments: `figsize`, `color`, `markopt`, `legend_below`, `t3base`, `figtitle`.

## Triple amplitudes

```julia
plot_t3amp(data)
plot_t3amp(data; logplot=true)     # log scale
plot_t3amp(data; color="mjd")      # colour by MJD
```

Key arguments: `figsize`, `logplot`, `color`, `markopt`, `legend_below`, `t3base`, `figtitle`.

## Visibility amplitudes and phases

```julia
plot_visamp(data)                  # adapts title to AMPTYP (absolute, differential, correlated flux)
plot_visamp(data; logplot=true)    # log scale
plot_visamp(data; color="wav")     # colour by wavelength
plot_visphi(data)                  # auto-detects PHITYP: absolute → single panel, differential → per-baseline subplots
plot_visphi(data; color="mjd")     # colour by MJD (absolute phases)
```

Key arguments: `figsize`, `logplot` (visamp only), `color`, `markopt`, `legend_below`, `figtitle`.

## Flux spectra

```julia
plot_flux(data)                    # coloured by wavelength (default)
plot_flux(data; color="station")   # colour by station; "Calibrated" for CALSTAT=C data
plot_flux(data; color="mjd")       # scatter coloured by MJD with colorbar
```

`plot_flux` maps `flux_sta_index == 0` (calibrated OI_FLUX) to the label
`"Calibrated"` and uses the station name otherwise.

Key arguments: `figsize`, `color`, `markopt`, `legend_below`, `figtitle`.

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
imdisp_multi(image_cube; labels=string.(wavs), pixsize=0.5)           # polychromatic
imdisp_multi(image_cube; labels=["Epoch $i" for i in 1:n], pixsize=0.5) # multi-temporal
```

`pixsize` is in mas/pixel. `beamsize` (mas) draws a PSF indicator.

## Facility / array layout

Plot the telescope positions from a `FacilityConfig`:

```julia
facility = read_facility_file("CHARA_new.toml")
plot_facility(facility)
```

See the [API Reference](@ref) for full docstrings.
