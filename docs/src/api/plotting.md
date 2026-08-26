# Plotting

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
| `plot_obs(data; obs, ...)` | Multi-panel figure with shared x-axis and legend |
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

The same figures, drawn in Makie instead of matplotlib. Load any Makie backend to enable them —
`CairoMakie` for vector output with no GPU and no display, `GLMakie` for a window — and no
Python is involved at any point. This is the plotting a PackageCompiler build can contain, and
it is what the GUI draws with.

| Function | Description |
|----------|-------------|
| `uvplot_makie(data)` | uv coverage, every sampled spatial frequency |
| `plot_v2_makie(data)` | V² against baseline — and `plot_t3phi_makie`, `plot_t3amp_makie`, `plot_visamp_makie`, `plot_visphi_makie`, `plot_flux_makie`, `plot_diffphi_makie` |
| `plot_observable_makie(data, kind)` | the same, with the observable chosen at run time; `kind` is any key of `OBS_SPECS` |
| `imdisp_makie(image)` | one image, East left and North up — the counterpart of `imdisp` |
| `imdisp_multi_makie(cube)` | a stack of images as a grid — the counterpart of `imdisp_multi` |
| `plot_residuals_makie(model, x, data)` | data over model with normalised residuals beneath — the counterpart of `plot_residuals` |
| `plot_corner_makie(posterior, names)` | corner plot of a posterior; needs `using PairPlots` |

`plot_ultranest_corner` is the matplotlib counterpart of `plot_corner_makie`, and takes UltraNest's
own result object rather than a sample matrix — so it serves that one sampler where
`plot_corner_makie` serves any of them.

| `plot_obs_makie(data)` | several observables stacked on a shared x axis — the counterpart of `plot_obs`/`plot_multi` |
| `plot_v2_multifile_makie(datasets)` | V² from several files, one legend entry per baseline — the counterpart of `plot_v2_multifile` |
| `plot_facility_makie(facility)` | telescope positions to scale — the counterpart of `plot_facility` |

Every name is the matplotlib one with `_makie` appended, so a function you already know leads
straight to its counterpart. They are separate names rather than extra methods on the
matplotlib functions because the two backends return different objects — a matplotlib figure
and a Makie one — and one name over both would have to lie about which. Each returns a
`PlotData`, whose `figure` field is ready for `Makie.save`:

```julia
using OITOOLS, CairoMakie
data = readoifits("data.oifits")[1, 1]
Makie.save("uv.pdf", uvplot_makie(data).figure)
Makie.save("v2.pdf", plot_v2_makie(data; logscale = true).figure)
```

Colours, groupings, axis labels and the palette come from the same tables the matplotlib
functions use — `OBS_SPECS` and `oiplot_colors` in the core package — so a baseline is the same
colour in both, and `test/gui/plotport.jl` asserts it stays that way.

Matplotlib-only for now: `set_oiplot_defaults` (it sets rcParams, which have no Makie
equivalent), the per-observable `plot_*_residuals` functions — `plot_residuals_makie` draws all of
them at once — and `plot_ultranest_corner`, which `plot_corner_makie` replaces.

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
