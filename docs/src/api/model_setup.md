# Model Setup

## Parametric models

| Function | Description |
|----------|-------------|
| `dict_to_model(model_dict, list_free_params)` | Compile a flat parameter dict into a `FlatModel` |
| `model_to_vis(model, x, uv)` | Evaluate complex visibilities for a model (alias: `eval_model`) |
| `eval_model_grad(model, x, uv)` | Evaluate visibilities + Jacobian |
| `eval_model_jac!(V, J, ws, model, x, uv)` | Visibilities + Jacobian from closed-form derivatives, in place (`ws = ModelJacWorkspace(model, x, uv)`) |
| `ModelJacCache()` | The Jacobian's buffers kept between calls; `jac_cache` of `chi2_flat_fg` and `residuals_flat_jac` |
| `eval_model_sub(model, x, uv; wl, dlam)` | Visibilities across each spectral channel, for bandwidth smearing |
| `display_model(model_dict, list_free_params)` | Pretty-print model parameters |
| `model_to_obs(model, x, data)` | Compute observables (V², T3amp, T3phi) from a model |
| `model_to_residuals(model, x, data)` | Compute normalised residuals (model - data) / error |
| `model_to_chi2(model, x, data)` | Compute weighted chi² (alias: `chi2_flat`) |
| `model_to_chi2_fg(model, x, data)` | Compute chi² + gradient (alias: `chi2_flat_fg`) |
| `model_to_image(model, x; nx, pixsize)` | Synthesize a model image via inverse FFT |
| `model_to_sed(model, x, wl_grid)` | Compute spectral energy distribution |
| `model_to_flux(model, x; wl)` | Total flux at zero baseline: `real(V(0,0))` |

The packages each fitter needs are listed under [What each optimiser needs](@ref).

```@docs
FlatModel
dict_to_model
parse_model
model_to_vis
eval_model
model_flux
eval_model_grad
eval_model_grad_ad
eval_model_jac!
eval_model_jac
ModelJacWorkspace
ModelJacCache
eval_model_sub
smear_samples
display_model
model_warnings
model_to_obs
model_to_residuals
model_to_chi2
model_to_chi2_fg
model_to_image
model_to_sed
model_to_flux
```

## Bandwidth smearing

A model with the global key `"smear" => n` (2 or more) is averaged over each spectral channel
wherever it meets data: `model_to_chi2`, its gradient, `residuals_flat`, `model_to_obs`, every
fitter built on them, and `simulate_from_oifits`. Each channel, `EFF_WAVE ± EFF_BAND/2`, is
sampled at `n` sub-wavelengths, evenly in wavenumber; V² is the channel mean of |V|², a triple
product the mean of the products, a visibility the mean visibility, as in PMOIRED and CANDID.
Smearing affects a companion near or beyond the smearing radius `R λ / (2 B_max)` and fills in
visibility nulls. Without the key, or with 0 or 1, there is no smearing. SPARCO's image + model
criterion refuses a smeared model.

`"smear" => "auto"` chooses `n` per dataset by CANDID's rule (four samples per cycle of the
largest phase change across a channel, at least 3) for a point at the model's extent: each
component's offset plus half its size, at the values the dict holds when the model is compiled.
[`smear_samples`](@ref) reports `n`.

`"smear_mode" => "coherent"` averages the visibility over the channel and forms the observables
from that average, as for an ABCD or coherent-flux measurement (GRAVITY) and in CANDID. The
default, `"incoherent"`, corresponds to power-spectrum estimators (fringe-scanning combiners).
Smearing of V² and closure phase depends on the fringe encoding and the estimator (Lachaume &
Berger 2013, MNRAS 435, 2501). For a binary the two modes differ by f²(1 − sinc²).

## Spectral lines

A component's flux can carry spectral lines, with PMOIRED's keys and formulas: `c,line_<n>_wl0`
(µm), `c,line_<n>_f` (the peak added to the continuum; < 0 for absorption) and one profile with
its width in nm: `c,line_<n>_lorentzian` (FWHM), `c,line_<n>_gaussian` (FWHM, with an optional
`c,line_<n>_power`), `c,line_<n>_truncexp` (with `c,line_<n>_pow`) or `c,line_<n>_trunclor`
(one-sided: redward for a positive width, blueward for a negative one). The component's flux at each
wavelength is its continuum (`c,f`, or 1) plus its lines, and weights its visibility, so the
lines appear in V, the phases, the closure phases and the flux. Each key is an ordinary
parameter; a line centre written as a velocity:
`"A,line_brg_wl0" => "2.16612 * (1 + \$vA / 2.998e5)"`.

```@docs
line_lorentzian
line_gaussian
line_truncexp
line_trunclor
```

## Visibility functions

Visibilities through the model dictionary:

```julia
model = dict_to_model(Dict{String,Any}("star,ud" => 6.5, "star,f" => 1.0), String[])
V     = eval_model(model, Float64[], data.uv)
```

The component types and their parameters are listed in [Component types](@ref).

`eval_model` returns `Complex{T}` (`T` from `x` and `uv`); the imaginary part is zero for a
centred, symmetric component.

Build the model once and pass the varying values as free parameters: the result is the same as
with the values in the dict, and `dict_to_model` takes about 15 µs against 3 µs for evaluating a
uniform disc at 300 uv points:

```julia
model = dict_to_model(Dict{String,Any}("star,ud" => 6.5, "star,f" => 1.0), ["star,ud"])
V     = eval_model(model, [diameter], data.uv)     # vary the value, keep the model
```

### One component at a time

The analytic visibility of one centred component, without a model dictionary; the dictionary
interface uses the same functions. They are `public`, not exported: `OITOOLS.vis_ud`.

| Function | Law |
|----------|-----|
| `vis_ud(θ, ρ)` | uniform disc |
| `vis_ldlin(θ, u, ρ)` | linear, `I(μ) = 1 - u(1-μ)` |
| `vis_ldquad(θ, u, w, ρ)` | quadratic, `I(μ) = 1 - u(1-μ) - w(1-μ)²` |
| `vis_ldsqrt(θ, c, d, ρ)` | square root, `I(μ) = 1 - c(1-μ) - d(1-√μ)` |
| `vis_ldclaret4(θ, c1, c2, c3, c4, ρ)` | Claret four-parameter |
| `vis_ldpow(θ, α, ρ)` | power law, `I(μ) = μ^α` |

`θ` is in mas and `ρ = √(u²+v²)` in cycles/rad, the unit of `OIdata.uv`:
`ρ = sqrt.(uv[1,:].^2 .+ uv[2,:].^2)`. The result is real. Any parameter may be a vector (one
value per wavelength). Gaussians and other elliptical components go through the dictionary.

```@docs
OITOOLS.vis_ud
OITOOLS.vis_ldlin
OITOOLS.vis_ldquad
OITOOLS.vis_ldsqrt
OITOOLS.vis_ldclaret4
OITOOLS.vis_ldpow
```
