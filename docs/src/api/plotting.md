# Plotting

The matplotlib figures need `using PythonPlot`; the `_makie` figures need a Makie backend.

| Function | Description |
|----------|-------------|
| `uvplot(data)` | Plot UV coverage |
| `plot_v2(data)` | Plot squared visibilities vs baseline |
| `plot_t3phi(data)` | Plot closure phases vs baseline |
| `plot_t3amp(data)` | Plot triple amplitudes vs baseline |
| `plot_visamp(data)` | Plot visibility amplitudes |
| `plot_visphi(data)` | Plot visibility phases |
| `plot_diffphi(data)` | Plot differential phases |
| `plot_flux(data)` | Plot flux vs wavelength |
| `plot_obs(data; obs, ...)` | Multi-panel figure with shared x-axis and legend (alias `plot_multi`) |
| `plot_v2_residuals(data, v2_model)` | Plot V² residuals against a model |
| `plot_t3phi_residuals(data, t3phi_model)` | Plot closure-phase residuals |
| `plot_t3amp_residuals(data, t3amp_model)` | Plot triple-amplitude residuals |
| `plot_visamp_residuals(data, visamp_model)` | Plot visibility-amplitude residuals |
| `plot_visphi_residuals(data, visphi_model)` | Plot visibility-phase residuals |
| `plot_residuals(data, obs)` | Combined residual plot for all observables |
| `plot_v2_multifile(data_vec)` | Overlay V² from multiple datasets |
| `plot_facility(facility)` | Plot telescope positions from a `FacilityConfig` |
| `imdisp(image; pixsize, ...)` | Display a single image |
| `imdisp_multi(cube; labels, ...)` | Display multiple images side by side |
| `set_oiplot_defaults()` | Apply matplotlib style settings |
| `plot_ultranest_corner(result)` | Corner plot of an UltraNest result |

```@docs
set_oiplot_defaults
uvplot
plot_v2
plot_t3phi
plot_t3amp
plot_visamp
plot_visphi
plot_diffphi
plot_flux
plot_obs
plot_residuals
imdisp
imdisp_multi
```

## Makie plotting

The same figures drawn with Makie, without Python: `CairoMakie` for files, `GLMakie` for a
window. The GUI and the compiled application use these.

| Function | Description |
|----------|-------------|
| `uvplot_makie(data)` | uv coverage, every sampled spatial frequency |
| `plot_v2_makie(data)` | V² against baseline; also `plot_t3phi_makie`, `plot_t3amp_makie`, `plot_visamp_makie`, `plot_visphi_makie`, `plot_flux_makie`, `plot_diffphi_makie` |
| `plot_observable_makie(data, kind)` | the same, the observable given by `kind`, a key of `OBS_SPECS` |
| `imdisp_makie(image)` | one image, East left and North up |
| `imdisp_multi_makie(cube)` | a stack of images as a grid |
| `plot_residuals_makie(model, x, data)` | data over model, normalised residuals below |
| `plot_obs_makie(data)` | several observables stacked on a shared x axis |
| `plot_v2_multifile_makie(datasets)` | V² from several files, one legend entry per baseline |
| `plot_facility_makie(facility)` | telescope positions to scale |
| `plot_corner_makie(posterior, names)` | corner plot of posterior samples (any sampler); needs `using PairPlots` |

Each name is the matplotlib one with `_makie` appended. Each returns a `PlotData`, whose
`figure` field can be passed to `Makie.save`:

```julia
using OITOOLS, CairoMakie
data = readoifits("data.oifits")[1, 1]
Makie.save("uv.pdf", uvplot_makie(data).figure)
Makie.save("v2.pdf", plot_v2_makie(data; logscale = true).figure)
```

Colours, groupings and axis labels match the matplotlib figures: `OBS_SPECS` holds the values
of `OBS_PLOT_SPECS`, both use the `oiplot_colors` palette, and `test/gui/plotport.jl` compares
the two backends' figures.

No Makie version: `set_oiplot_defaults` (matplotlib rcParams), the per-observable
`plot_*_residuals` functions (`plot_residuals_makie` draws all observables) and
`plot_ultranest_corner` (use `plot_corner_makie`).

```@docs
uvplot_makie
plot_observable_makie
plot_v2_makie
plot_t3phi_makie
plot_t3amp_makie
plot_visamp_makie
plot_visphi_makie
plot_flux_makie
plot_diffphi_makie
imdisp_makie
imdisp_multi_makie
plot_residuals_makie
plot_obs_makie
plot_v2_multifile_makie
plot_facility_makie
plot_corner_makie
plot_ultranest_corner
PlotData
```

## Graphical interface

[`oitoolsgui`](@ref) launches the GUI (see [The GUI](@ref gui-install)); the other functions are
the steps it runs.

```@docs
oitoolsgui
gui
configure_graphics!
configure_qt_platform!
is_wsl
resource_dir
```
