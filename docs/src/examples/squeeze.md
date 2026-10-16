# SQUEEZE (MCMC)

`reconstruct_squeeze` is a Julia port of SQUEEZE: image reconstruction by Markov-chain Monte
Carlo with simulated annealing. `reconstruct_squeeze_tempered` samples the same posterior by
parallel tempering (Pigeons.jl).

The image is `nelements` flux elements of `1/nelements` each, on an `nx × nx` pixel grid; the
result is the normalised histogram of their positions. Consequences:

- the image is non-negative and its total flux fixed, without constraints;
- moving one element changes the model visibilities by one column difference of the Fourier
  matrix, an `O(nuv)` update;
- no gradient is needed, so non-convex regularizers such as L0 are available;
- the result is a mean over samples, with a spread across chains.

## Basic use

The image size and pixel size come from the Fourier operator, as for
[`reconstruct_bsmem`](@ref):

```julia
using OITOOLS, PythonPlot

data = readoifits("2004-data1.oifits")
ft   = setup_ft(data, 64, 0.2; mode = "dft")

img, diag = reconstruct_squeeze(data, ft;
                                niter   = 400,
                                nchains = 4,
                                seed    = 1)

println("chi2r per chain: ", diag.chi2r_mean)
imdisp(img, pixsize = 0.2)
```

Use `mode = "dft"`: each move needs one column of the transform. An NFFT plan is accepted, but a
DFT matrix is then built from it.

The DFT matrix takes `nuv · nx² · 2 · sizeof(T)` bytes:

| | nx=64 | nx=128 | nx=256 |
|---|---|---|---|
| nuv=500, `Float64` | 33 MB | 131 MB | 524 MB |
| nuv=2000, `Float64` | 131 MB | 524 MB | 2.1 GB |
| `Float32` | ÷2 | ÷2 | ÷2 |

Data read with `T = Float32` halve the memory and run about twice as fast, to the same χ²r.

## Starting image

| first argument | start |
|---|---|
| none, or `:point_source` | all elements on the central pixel (SQUEEZE's default) |
| `:random` | elements placed uniformly at random |
| an image | the image digitised into elements (SQUEEZE's `-i`) |

Strings (`"point_source"`, `"random"`) are accepted. An image start is deterministic, so all
chains start from the same configuration; use `:random` for independent starts.

```julia
reconstruct_squeeze(:random, data, ft; niter = 400, nchains = 4)
```

## Prior images and masks

`prior_image` (SQUEEZE's `-p`) is an `nx × nx` map of per-pixel prior probabilities. It adds
`-log(p)` per element; pixels with `p <= 0` cannot hold an element.

```julia
c = (nx + 1) / 2
mask = [sqrt((i-c)^2 + (j-c)^2) <= 20 ? 1.0 : 0.0 for i in 1:nx, j in 1:nx]

img, diag = reconstruct_squeeze(data, ft; niter = 300, prior_image = mask)
@assert sum(img[mask .== 0]) == 0
```

Its weight is 1, or λ with `["priorimage", λ]` in `regularizers`. A prior image turns
auto-centering off.

On the 2004 data (source about 5 mas across), masks of decreasing radius:

| mask radius | χ²r |
|---|---|
| 12 px (2.4 mas) | 739 |
| 18 px (3.6 mas) | 1.11 |
| 25 px (5.0 mas) | 0.750 |
| none | 0.759 |

## Output

`reconstruct_squeeze` returns `(image, diag)`. The image is the best chain's mean over the
samples after burn-in, normalised to unit flux. `diag`, per chain:

| field | |
|-------|---------|
| `chi2r_mean` | reduced χ² of the chain's mean image; the best chain is chosen on it |
| `chi2r_last` | reduced χ² of the last sample |
| `acceptance` | acceptance rate |
| `temperature` | final temperature |
| `nsamples`, `burn_in` | number of samples, and the sweep where sampling started |
| `regs_final` | each regularizer's value on the final image |
| `params`, `stepsize` | SPARCO parameters (empty without a model) |
| `best_chain`, `images`, `ndf` | best chain, all chain images, degrees of freedom |

`outfile = "result.fits"` writes the best chain's image, with its WCS when `pixsize` is given.

## Number of sweeps

A sweep is `nelements` proposals. The annealing has a fixed point at `chi2r = chi2_temp · T`; a
chain that has not cooled past it stays at a large χ²r. On the 2004 Beauty Contest data (64 px,
0.2 mas):

| `niter` | chains stuck at χ²r > 2 |
|---------|---------------------------|
| 100 | 1 in 8 |
| 250 | 0 in 8 |
| 500 | 0 in 8 |

Use a few hundred sweeps or more, and several chains.

## Monitoring

`print_every = n` prints a line every `n` sweeps, in SQUEEZE's format:

```
Chain: 1 lPost:  1132.6 lPrior:   481.9 lLike:   650.7 V2: 0.83 T3P: 0.61 \
l0:417.00 tv:64.89 E:  1230 MPr: 0.61 T:  1.00 B:19 Iter:   60 of   60
```

`lPost`, `lPrior`, `lLike`: log-posterior terms; `V2`, `T3A`, `T3P`: reduced χ² per observable
(each divided by its own count); regularizers: `λ·value`; `E`: element count; `MPr`: acceptance
average; `T`: temperature; `B`: burn-in sweep.

`monitor = n` (with PythonPlot) draws the image and traces of χ²r, temperature (or Pigeons
round) and the free SPARCO parameters every `n` sweeps or scans:

```julia
reconstruct_squeeze(data, ft; niter = 400, pixsize = 0.2, monitor = 50)
reconstruct_squeeze_tempered(data, ft; n_rounds = 11, pixsize = 0.2, monitor = 50)
```

An update takes about 40 ms. The display has to be drawn from the main thread: with `monitor > 0`,
`reconstruct_squeeze` runs chain 1 on the main thread and the others on threads (all serially if
called from another thread), and `reconstruct_squeeze_tempered` runs without multithreading.

## Threads and reproducibility

Chains run in parallel on Julia's threads (`julia -t N`). Each chain's random numbers depend only
on `(seed, chain)`, so a seed gives identical results whatever the number of threads.

## SPARCO

[`SqueezeSparco`](@ref) adds a parametric star and an over-resolved background to the image
(Kluska et al. 2012):

```
        f_star (λ/λ0)^-4 V_star  +  f_bg (λ/λ0)^bg_indx · 0  +  f_env V_env
V_tot = ───────────────────────────────────────────────────────────────────
        f_star (λ/λ0)^-4         +  f_bg (λ/λ0)^bg_indx      +  f_env
```

The star is a uniform disc (`ud = 0` for a point source); the background adds flux without
visibility. The image is the environment, grey; the chromatism is in the model, so one `OIdata`
with several wavelengths is enough.

```julia
data = readoifits("V1295_Aql.oifits"; filter_bad_data = true)
ft   = setup_ft(data, 64, 0.125; mode = "dft")

model = SqueezeSparco(f_star  = 0.48,
                      ud      = 0.0,
                      env_indx = -0.1,
                      lambda0 = 1.6e-6,
                      free    = (:f_star, :env_indx),
                      stepsize = [0.01, 0.0, 0.05, 0.0, 0.0, 0.0])

img, diag = reconstruct_squeeze(data, ft;
                                niter   = 400, nchains = 4,
                                weights = [1.0, 0.0, 1.0],   # this dataset's T3amp is noisy
                                regularizers = [["tv", 100.0]],
                                model = model, seed = 1)

println("f_star = ", diag.params[diag.best_chain][1])
```

`free` names the parameters that vary (by name or index). On V1295 Aql, χ²r goes from about 69
(image alone) to about 0.9.

Centering is on without a model and off with one; `auto_centering = false` turns it off.

!!! warning "SPARCO needs several wavelengths"
    On monochromatic data a central point star cannot be told apart from a central component of
    the image, and `f_star` is not constrained (in a test, `f_star = 0.2` fitted better than the
    true 0.4). With one channel, fix the model parameters.

!!! warning "Annealed parameters are not posterior samples"
    Under annealing the parameter moves are not reversible, so the spread of parameters across
    chains is not an uncertainty. Use the tempered version for a posterior, or
    [`bootstrap_fit`](@ref) / [`fit_model`](@ref).

Shape parameters, `ud` above all, often do not move under annealing: once the image has adapted
to a value, a change of `ud` with the image fixed is rejected (on V1295 Aql, 0.05 mas costs
ΔlogL ≈ 2.8 × 10⁴), and the step size shrinks to zero. Alternatives:

1. fix `ud` at a grid of values and compare `chi2r_mean`;
2. start from a model fit ([`fit_model`](@ref), [`optimize_sparco_flat_parameters`](@ref));
3. check `MPr` in the `print_every` output: near zero, the parameter is stuck;
4. use the tempered version (below).

## Parallel tempering

`reconstruct_squeeze_tempered` (needs `using Pigeons`) runs a temperature ladder with Pigeons.jl
and returns the evidence:

```julia
using OITOOLS, Pigeons

img, diag = reconstruct_squeeze_tempered(data, ft; n_rounds = 11, n_chains = 30)

diag.logZ            # log evidence (stepping stone)
diag.global_barrier  # Λ
diag.schedule        # the β ladder
```

It takes the same model options as `reconstruct_squeeze` (regularizers, `prior_image`,
`weights`, `nelements`, `model`). Round *i* makes 2ⁱ scans; only the last round's samples are
kept.

The tempered family is `π_β(x) ∝ exp(β · (−lLikelihood(x) − lPrior(x)))`, so at β = 0 the
elements are uniformly distributed, a reference sampled exactly. Unlike SQUEEZE's annealing, the
prior is tempered with the likelihood, and the asymmetric `COPYCAT` move is not used; the
explorer is a fixed mixture of the SMALL, MEDIUM and ANYWHERE moves.

**Λ** (`diag.global_barrier`): round trips between β = 0 and 1 take about Λ² scans; use
`n_chains` of order Λ. Λ is estimated from swap rates and comes out low with too few chains. On
the 2004 data at 64 px:

| `n_chains` | Λ | χ²r |
|---|---|---|
| 8 | 6.96 | 0.79 |
| 16 | 14.9 | 1.01 |
| 30 | 27.1 | 0.77 |

If `logZ` changes when chains are added, the run has not converged.

**SPARCO under tempering**: `model` is accepted, and shape parameters move: at low β the
likelihood is weak, and swaps carry the parameters to β = 1. On simulated V1295 Aql data with a
true `ud` of 2.0 mas:

| start | annealed | tempered |
|---|---|---|
| `ud = 4.0` | 3.73 | 2.15 |
| `ud = 0.5` | 1.93 | 2.17 |

`f_star` comes out at 0.4697 and 0.4698 from starts of 0.15 and 0.80. Parameter moves are a
symmetric Gaussian random walk with a fixed step. `diag.params` is the posterior mean over the
β = 1 samples; `diag.params_by_replica` holds each replica's current values, mostly at low β.

Free parameters need a bounded box, the reference being uniform in it: `default_sparco_bounds`,
or `param_bounds = (lo, hi)`.

## Comparison

| | method | result |
|---|---|---|
| [`reconstruct`](@ref) | gradient descent (VMLMB) | MAP image, differentiable regularizers |
| [`reconstruct_bsmem`](@ref) | maximum entropy | MAP image, entropy prior |
| [`reconstruct_bsdmm`](@ref) | ADMM | MAP image, non-smooth regularizers |
| `reconstruct_squeeze` | MCMC, annealing | mean image, non-convex regularizers |
| `reconstruct_squeeze_tempered` | MCMC, parallel tempering | mean image and `logZ` |

## Not implemented

Wavelet regularizers, polychromatic cubes, median and mode images, the full chain history, and
bandwidth smearing.
