# Model Fitting

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
| `fit_model(model, x0, data)` | Fit a `FlatModel` via NLopt |
| `fit_model_lsqfit(model, x0, data)` | Fit a `FlatModel` via Levenberg-Marquardt |
| `fit_model_nested(model, data; lb, ub)` | Fit by nested sampling — posterior **and** log-evidence; needs finite bounds |
| `fit_model_ultranest(model, data; lb, ub)` | `fit_model_nested` pinned to the UltraNest backend |
| `nested_backend()`, `set_nested_backend!(b)` | Which nested sampler is in force: `:nautilus` or `:ultranest` |

Which package each fitter needs is tabulated under
[What each optimiser needs](@ref) — `fit_model`, `fit_model_lsqfit`, `chi2_map` and
`bootstrap_fit` need nothing beyond OITOOLS; `fit_model_nested` needs a sampler.
| `chi2_map(model_dict, free, data, p1, p2)` | Grid-search two free parameters and return the χ² surface; needs no gradient and no starting guess |
| `model_to_obs(model, x, data)` | Compute observables (V², T3amp, T3phi) from a model |
| `model_to_residuals(model, x, data)` | Compute normalised residuals (model - data) / error |
| `model_to_chi2(model, x, data)` | Compute weighted chi² (alias: `chi2_flat`) |
| `model_to_chi2_fg(model, x, data)` | Compute chi² + gradient (alias: `chi2_flat_fg`) |
| `model_to_image(model, x; nx, pixsize)` | Synthesize a model image via inverse FFT |
| `model_to_sed(model, x, wl_grid)` | Compute spectral energy distribution |
| `model_to_flux(model, x; wl)` | Total flux at zero baseline: `real(V(0,0))` |

### Bandwidth smearing

A model with the global key `"smear" => n` (2 or more) is averaged over each spectral channel
wherever it meets data: `model_to_chi2`, its gradient, `residuals_flat`, `model_to_obs`, every
fitter built on them, and `simulate_from_oifits`. Each channel, `EFF_WAVE ± EFF_BAND/2`, is
sampled at `n` sub-wavelengths, evenly in wavenumber; V² is the channel mean of |V|², a triple
product the mean of the products, a visibility the mean visibility. That is what a
power-spectrum estimator measures across a channel, and how PMOIRED and CANDID model it. It
matters for a companion near or beyond the smearing radius `R λ / (2 B_max)` and for a
visibility null, which it fills in. Absent, 0 or 1, nothing changes. SPARCO's image + model
criterion refuses a smeared model.

`"smear" => "auto"` chooses `n` per dataset by CANDID's rule — four samples per cycle of the
largest phase change across a channel, at least 3 — for a point at the model's extent (each
component's offset plus half its size, as the dict holds them when the model is compiled);
[`smear_samples`](@ref) reports the number.

`"smear_mode" => "coherent"` averages the visibility over the channel instead and forms the
observables from that average: what an ABCD / coherent-flux measurement gives (GRAVITY), and
how CANDID models smearing. The default, `"incoherent"`, suits power-spectrum estimators
(fringe-scanning combiners). Lachaume & Berger (2013, MNRAS 435, 2501) show the smearing of
V² and closure phase depends on the fringe encoding and the estimator; each is an
approximation to a given instrument. For a binary the two differ by f²(1 − sinc²).

## Companion search

| Function | Description |
|----------|-------------|
| `companion_search(model_dict, data; companion, free, ...)` | Grid search for a companion around any primary (a disc, an image, anything): the companion fitted at every node, a χ² map, the best node and its significance |
| `companion_grid(model_dict, data; companion, step, rmin, rmax, mask)` | The grid `companion_search` would search, with its defaults, and how many nodes it holds |
| `companion_flux_map(search)` | The companion's fitted flux at every node |
| `detection_sigma(chi2r_null, chi2r_best, ndof)` | CANDID's significance convention: the null model's χ² against a χ² distribution, as Gaussian σ |

Only the companion's flux free is the fast path (CANDID's `chi2Map`): the model is linear in
it, so a node is a 1-D fit (BraDi, or Halley with a BraDi fallback). Anything else free at a
node — the primary's diameter, a companion diameter — runs VMLMB with bounds. The grid
defaults come from the data: step λ/(2B_max), out to the smearing radius Rλ/(2B_max); an image
primary's bright pixels are masked. Bandwidth smearing follows the model's `"smear"`.

```@docs
companion_search
companion_grid
CompanionSearch
companion_flux_map
detection_sigma
```

## Spectral lines

A component's flux can carry spectral lines, with PMOIRED's keys and formulas: `c,line_<n>_wl0`
(µm), `c,line_<n>_f` (the peak added to the continuum; < 0 for absorption) and one profile with
its width in nm — `c,line_<n>_lorentzian` (FWHM), `c,line_<n>_gaussian` (FWHM, with an optional
`c,line_<n>_power`), `c,line_<n>_truncexp` (with `c,line_<n>_pow`) or `c,line_<n>_trunclor`
(one-sided: redward for a positive width, blueward for a negative one). The component's flux is
then its continuum (`c,f`, or 1) plus its lines at each wavelength, which weights its
visibility, so the lines show in V, the phases, the closure phases and the flux. Each key is a
parameter like any other: EX5 writes a centre as a velocity,
`"A,line_brg_wl0" => "2.16612 * (1 + \$vA / 2.998e5)"`.

```@docs
line_lorentzian
line_gaussian
line_truncexp
line_trunclor
```

## Many starts

`multistart_fit` is PMOIRED's `gridFit`: Levenberg-Marquardt from every node of a grid
(`grid = ["p" => (min, max, step)]`) and/or from uniform or normal draws (`rand`, `randn`,
`nstarts`), starts outside `constrain` dropped, then the fits sorted as PMOIRED's `analyseGrid`
does — set aside when failed, without uncertainties or unmoved, the rest grouped into distinct
minima, lowest χ²ᵣ first. An explored parameter that is fixed in the model is held at each
start's value (PMOIRED's grid parameter in `doNotFit`): a grid over held parameters is a χ²
profile, every node fitting the rest. The fits run on threads, one compiled model each; start i draws from
`Xoshiro(seed + i)`, so the result does not depend on the thread count. For a component on an
orbit, minima whose MJD0 differ by whole periods are one.

```julia
r = multistart_fit(m, ["B,f", "B,x", "B,y"], data; weights = [0, 0, 1, 1, 0, 0, 0],
                   grid = ["B,x" => (-12, 12, 1.5), "B,y" => (-12, 12, 1.5)],
                   constrain = [("sqrt(\$B,x^2 + \$B,y^2)", "<=", 12)])
r.minima[1].fit          # the global minimum, an LsqFitResult
```

```@docs
multistart_fit
MultiStartResult
MultiStartMinimum
```

## Orbits and radial velocities

A component on an orbit takes PMOIRED's keys and PMOIRED's definition (`_orbit`):
`"c,x" = "c,y" = "orbit"`, with the elements `c,orb P` (days), `MJD0` (periastron, MJD), `e`,
`omega`, `OMEGA`, `incl` (degrees) and `a` (mas) — or `M` (M☉) and `plx` (mas) in place of `a`.
Its position is computed at each data point's own MJD, once per epoch (an expression of `$MJD`
alone is evaluated once per run of equal MJDs), in Float64 whatever the data's precision. ω is
the secondary's, relative to the
primary (orbitize!'s convention, 180° from Gallenne et al. 2023's). Each element is a parameter
like any other: fixed, free, or an expression.

Radial velocities ([`read_rv`](@ref), a CSV `MJD,RV,sigma,star` with star `a` or `b`) are fitted
with the interferometry by `fit_model_lsqfit(…; rv)`, `fit_model(…; rv)` and
`bootstrap_fit(…; rv)`, against the velocity elements `c,orb Ka`, `Kb`, `gamma` — or `M` and
`q` (= Mb/Ma) with `incl`. `chi2` and `ndof` count them; `chi2_rv` is their part.

```julia
m = Dict{String,Any}("1,ud" => 0.2, "2,f" => 0.37, "2,x" => "orbit", "2,y" => "orbit",
                     "2,orb P" => 42.38, "2,orb MJD0" => 59785.45, "2,orb e" => 0.623,
                     "2,orb omega" => 114.1, "2,orb OMEGA" => 77.4, "2,orb incl" => 71.0,
                     "2,orb a" => 3.75, "2,orb Ka" => 50.3, "2,orb Kb" => 58.8, "2,orb gamma" => 14.9)
rv = read_rv("demos/data/HD210763/rv_gallenne2023.csv")
r  = fit_model_lsqfit(m, ["2,orb e", "2,orb a", "2,orb Ka", "2,orb Kb", "2,orb gamma"], data; rv)
```

| Function | Description |
|----------|-------------|
| `orbit_x(t, P, MJD0, e, omega, OMEGA, incl, a)`, `orbit_y(…)` | The secondary's offset East and North (mas) at MJD `t` |
| `orbit_rv(t, star, P, MJD0, e, omega, K, gamma)` | A star's radial velocity (km/s) |
| `orbit_a(P, M, plx)`, `orbit_semiamplitudes(P, e, incl, M, q)` | Kepler's third law: the semi-major axis, and Ka, Kb |
| `kepler_E(M, e)` | The eccentric anomaly, with exact derivatives for dual numbers |
| `read_rv(path)`, `RVData` | Radial velocities from a CSV table |

```@docs
orbit_x
orbit_y
orbit_rv
orbit_a
orbit_semiamplitudes
kepler_E
read_rv
RVData
```

## Bounds, constraints and model files

`lb`/`ub` describe a box, one parameter at a time. A relation *between* parameters — a ring's
outer diameter exceeding its inner one, two flux fractions summing to one — is not a box, and
needs [`ModelConstraint`](@ref).

| Function | Description |
|----------|-------------|
| `default_bounds(model_dict, free; data, max_size)` | Suggested `lb`/`ub` per free parameter; pass `data` and angular sizes are capped at `2 λ/B_min` from the actual uv coverage |
| `max_angular_scale(data)` | The largest angular scale the shortest baseline senses, in mas |
| `ModelConstraint(lhs, op, rhs; tol)` | A relation between parameters, with `op` one of `<`, `<=`, `>`, `>=`, `=` |
| `parse_constraints(specs)` | Accept constraints as `ModelConstraint`s, `(lhs, op, rhs[, tol])` tuples (the PMOIRED layout) or dicts |
| `check_constraints(constraints, model_dict)` | Which constraints the starting model already satisfies |
| `read_model_file(path)` | Read a TOML model file into `(; model, free, lb, ub, constraints, priors, name)` |
| `write_model_file(path, model_dict; free, lb, ub, constraints, priors)` | Write all of that back out |

`fit_model` hands constraints to NLopt as real nonlinear constraints, so they hold at the
optimum rather than being encouraged there; an algorithm that cannot take them (including the
default `:LD_LBFGS`) is wrapped in `:AUGLAG` rather than replaced. `fit_model_lsqfit` and
`fit_model_nested` have no such machinery and use a one-sided quadratic penalty on the
normalised violation, matching PMOIRED's `prior` list — soft, and so able to lose to a steep
χ². The distinction is worth knowing before choosing a fitter for a constrained model.

A model dict alone does not describe a fit: the free list, the bounds, the constraints and the
priors all change the answer and none of them lived in a file before. A TOML model file carries
all five.

```toml
free = ["star,ud", "disk,pa"]

[model]
"star,ud"   = 6.5
"disk,fwhm" = "$star,ud * 3"

[bounds]
"star,ud" = [0.0, 20.0]

[[constraints]]
param = "disk,diamout"
op    = ">"
value = "disk,diamin"
tol   = 0.001

[[priors]]
expr   = "star,ud"
target = 6.0
sigma  = 0.5
```

```julia
m = read_model_file("binary.toml")
res = fit_model(m.model, m.free, data; m.lb, m.ub, m.constraints, m.priors)
```

`free` is a top-level key rather than a member of `[model]` so that `[model]` mirrors the model
dict exactly — a model may hold a bare global key, and one named `free` would otherwise be
eaten by the free-parameter list.

## Uncertainty estimation by resampling

| Function | Description |
|----------|-------------|
| `bootstrap_fit(model_dict, list_free_params, data)` | Nonparametric block bootstrap: refit replicates in which blocks of data are resampled |
| `bootstrap_driver(fitfun, x_opt, list_free_params)` | The model-agnostic replicate loop and statistics behind `bootstrap_fit`, for callers with their own fitter or resampling unit |
| `data_blocks(data; granularity)` | Partition data into resampling blocks (`:config`, `:epoch`, `:point`) |
| `resample_blocks(data, blocks; mode)` | One bootstrap replicate (`:replacement`, `:halfsample`, `:weights`; `:pmoired` reproduces PMOIRED's scheme and is biased low by √2) |
| `block_counts(nblocks, mode)` | Block multiplicities drawn by a resampling scheme |
| `block_weights(nblocks)` | Continuous block weights for the multiplier (Bayesian) bootstrap |
| `apply_block_weights(data, blocks, w)` | Build the weighted replicate (error bars scaled by 1/√w) |
| `apply_block_counts(data, blocks, counts)` | Build the replicate for given block multiplicities |
| `perturb_data(data)` | Add Gaussian noise drawn from the error bars — a simulation utility, not an uncertainty estimator (was `resample_data`) |

`bootstrap_fit` resamples *which observations are used*, in blocks of (MJD,
telescope configuration), and therefore responds to correlated calibration
errors and to mis-stated error bars — neither of which the analytic covariance
of `fit_model_lsqfit` can see.  See the model-fitting examples page, and
`demos/bootstrap_validation` for the calibration test behind that statement.

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
default_bounds
DEFAULT_MAX_SIZE_MAS
max_angular_scale
ModelConstraint
parse_constraints
check_constraints
DEFAULT_CONSTRAINT_TOL
read_model_file
write_model_file
fit_model
fit_model_lsqfit
fit_model_nested
fit_model_ultranest
nested_backend
set_nested_backend!
chi2_map
delta_chi2_levels
model_warnings
model_to_obs
model_to_residuals
model_to_chi2
model_to_chi2_fg
model_to_image
model_to_sed
model_to_flux
bootstrap_fit
bootstrap_driver
data_blocks
resample_blocks
block_counts
apply_block_counts
block_weights
apply_block_weights
perturb_data
resample_data
BootstrapResult
NestedResult
UltraNestResult
FitResult
Chi2Map
DataBlocks
```

| Type | Description |
|------|-------------|
| `FitResult` | Result from `fit_model` (fields: `x_opt`, `chi2r`, `model`, ...) |
| `LsqFitResult` | Result from `fit_model_lsqfit` (adds `stderror`, `covar`, `converged`) |
| `NestedResult` | Result from `fit_model_nested` (adds `logz`, `logzerr`, `posterior`, `result`, `backend`) |
| `UltraNestResult` | Alias for `NestedResult` |
| `BootstrapResult` | Result from `bootstrap_fit` (adds `samples`, `median`, `sigma_minus`, `sigma_plus`, `covar`) |
| `Chi2Map` | χ² surface over two free parameters from `chi2_map`; `FitResult(map)` gives its best grid point |
| `DataBlocks` | Partition of an `OIdata` into resampling blocks |
| `ModelConstraint` | One relation between model parameters, enforced during fitting |

## Visibility functions

Evaluate a model's visibilities through the model dictionary:

```julia
model = dict_to_model(Dict{String,Any}("star,ud" => 6.5, "star,f" => 1.0), String[])
V     = eval_model(model, Float64[], data.uv)
```

The geometry `kind` names are listed with their parameters under **Parametric models** above.

`eval_model` returns `ComplexF64`, with the imaginary part exactly zero for a component at the
centre — a component's phase comes from its offset. Code that only wants the modulus should take
`real.()` or `abs.()` at the boundary rather than carrying complex arithmetic downstream.

**Build the model once, outside the loop.** Naming the parameters free and passing their values
to `eval_model` gives bit-identical results to baking them into the dict, and `dict_to_model` is
~10 µs against ~3 µs for the evaluation itself — so rebuilding it per call dominates at the few
hundred uv points of a single epoch:

```julia
model = dict_to_model(Dict{String,Any}("star,ud" => 6.5, "star,f" => 1.0), ["star,ud"])
V     = eval_model(model, [diameter], data.uv)     # vary the value, keep the model
```

### One component at a time

For a single centred component these evaluate the analytic visibility directly, without a model
dictionary. They are what the dictionary interface calls underneath, so there is no second
implementation to disagree with it. `public`, not exported: reach them as `OITOOLS.vis_ud`.

| Function | Law |
|----------|-----|
| `vis_ud(θ, ρ)` | uniform disc |
| `vis_ldlin(θ, u, ρ)` | linear, `I(μ) = 1 - u(1-μ)` |
| `vis_ldquad(θ, u, w, ρ)` | quadratic, `I(μ) = 1 - u(1-μ) - w(1-μ)²` |
| `vis_ldsqrt(θ, c, d, ρ)` | square root, `I(μ) = 1 - c(1-μ) - d(1-√μ)` |
| `vis_ldclaret4(θ, c1, c2, c3, c4, ρ)` | Claret four-parameter |
| `vis_ldpow(θ, α, ρ)` | power law, `I(μ) = μ^α` |

`θ` is in **mas** and `ρ = √(u²+v²)` in **cycles/rad**, the unit `OIdata.uv` is stored in, so
`ρ = sqrt.(uv[1,:].^2 .+ uv[2,:].^2)` needs no conversion. The result is **real**. Any parameter
may be a vector, for a quantity that varies with wavelength.

There is deliberately no Gaussian here: an ellipse needs an inclination and a position angle, and
applying those belongs to the model geometry. Go through the dictionary for it.

```@docs
OITOOLS.vis_ud
OITOOLS.vis_ldlin
OITOOLS.vis_ldquad
OITOOLS.vis_ldsqrt
OITOOLS.vis_ldclaret4
OITOOLS.vis_ldpow
```

!!! note
    The standalone `visibility_*(param, uv)` functions were removed in 0.13.2 — the defect was
    the calling convention, a packed parameter vector read by index, not the per-component
    functions themselves. `vis_ud` and the rest above are that same physics with named
    arguments. Porting an old call: `visibility_ud([D], uv)` → `vis_ud(D, ρ)`, and
    `visibility_Gaussian([FWHM, i, ϕ], uv)` → a dictionary with `"c,fwhm"`, `"c,incl"` and
    **`"c,pa" => 90 - ϕ`** — the old angle was measured from the u axis, `pa` is a position
    angle, and mapping it straight across misorients the ellipse without raising anything.
