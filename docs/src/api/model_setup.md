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

### Free-form limb darkening

A disc whose I(μ) is set by its values at N nodes μ_k = (k-1)/(N-1):

```julia
md = Dict{String,Any}("star,ldspline" => 8.4,         # diameter, mas
                      "star,ld1" => 0.60, "star,ld2" => 0.70,
                      "star,ld3" => 0.80, "star,ld4" => 0.90,   # I(μ_k)/I(1); I(1) = 1
                      "star,f" => 1.0)
```

`ldspline` interpolates a natural cubic spline in μ through the nodes, `ldpiecewise` straight
lines. N is the number of `ld<k>` keys plus one. The visibility is exact for any node values:
the Hankel transform of each basis function is tabulated once per (basis, N), with its
derivative, so the cost and the Jacobian are those of a closed-form law. A linear law is
represented exactly by either basis.

| Function | Description |
|----------|-------------|
| `ld_free!(model_dict, free, comp; n, basis, monotonic, lb, ub)` | Make a disc free-form, its nodes taken from its current profile; `monotonic` fits the ratios `ldq<k>` ∈ [0, 1] |
| `ld_profile(model_dict, comp, μ)` | I(μ)/I(1) of any disc law or free-form disc |
| `ld_nodes(n)` | The node positions in μ |

```@docs
ld_free!
ld_profile
ld_nodes
```

### Error bars on the profile

`fit_ld_profile_vi` (needs `using VarInf`) returns I(μ) with a band. The disc's profile is put
on `ngrid` nodes (24) joined by straight lines, log I = log I₀ + f with I₀ the disc's current
profile and f a Gaussian process (Matérn-3/2, length 0.3 in μ, amplitude 0.3 in log I); the
diameter is log-normal about its current value with a 10 % width. The run starts at the
posterior's maximum, found from a scan of the diameter over ±22 %, and draws its samples with
geoVI (or MGVI).

```julia
using OITOOLS, VarInf
m = Dict{String,Any}("star,ldlin" => 8.4, "star,u" => 0.3, "star,f" => 1.0)
r0 = fit_model_lsqfit(m, ["star,ldlin", "star,u"], data)
m["star,ldlin"], m["star,u"] = r0.x_opt
r = fit_ld_profile_vi(m, "star", data; scale_errors = true)
r.μ, r.I_mean, r.I_lo, r.I_hi          # I/I(1) and its 16–84 % band
r.diameter, r.diameter_std
```

- Start from a law fitted to the same data. The profile's prior is centred on that law.
- Visibilities past the first null leave several minima in the diameter. The scan covers
  ±22 % of the starting value. `chi2r` far above `chi2r_map` in the result means the run did
  not converge.
- The band assumes the error bars. With `scale_errors = true`, they are multiplied by √χ²r when
  the free-form profile's best fit has χ²r > 1 (`error_scale` in the result).
- `bootstrap = n` (n ≥ 2) is bagged VI: the run is repeated on n block-bootstrap replicates
  of the data (configurations drawn with replacement) and their samples pooled. The band then
  holds the spread of the profile between replicates as well, which error bars or a model
  that do not describe the data widen. `single` in the result is the run on the full data.
  The replicates run on all threads; 20 replicates take about 3 s on 4 threads for α Cen A.
- Closure phases are not used by default: those of a symmetric disc are 0 or 180°, which V²
  already gives.
- Calibration, on data simulated with α Cen A's uv coverage and noise (σ(V²) = 2.5 %·V² +
  3.5e-4): with profiles drawn from the prior, the truth lies in the 16–84 % band at 66–68 % of
  the nodes for a disc resolved into the third lobe, and at 57 % with data in the first lobe
  only; against Hamiltonian Monte Carlo on the same posterior the band is 1.0–1.13 times as
  wide and the mean within 0.25σ. For profiles from the usual laws the band is conservative
  (90–98 %). Bagged over 20 replicates, the band holds the truth at 69–73 % of the nodes in
  both cases (⟨z²⟩ 0.88–0.93), and at 96–100 % for the usual laws. On α Cen A, whose error
  bars are too small by 2.2, the bagged band without `scale_errors` is 1.25 times the band of
  a single run with it.
- MGVI is faster; geoVI follows the curved posterior that the diameter and the limb
  intensity make together.
- `prior = :learned` infers the Gaussian process's spectrum (VarInf's correlated field)
  instead of fixing it. Its bands were not calibrated in the same tests (36–39 % coverage).

```@docs
fit_ld_profile_vi
LDProfileVIResult
```
