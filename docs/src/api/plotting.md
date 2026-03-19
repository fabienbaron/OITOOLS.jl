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
