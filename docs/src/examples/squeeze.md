# SQUEEZE (MCMC)

`reconstruct_squeeze` is a Julia port of the SQUEEZE algorithm: optical
interferometric image reconstruction by Markov-chain Monte Carlo rather than by
gradient descent.

## The idea

The image is **not** a vector of pixel intensities. It is a bag of `nelements`
discrete flux quanta ("elements"), each carrying flux `1/nelements`, living on an
`nx × nx` integer pixel grid. The image you get back is the *histogram* of element
positions, normalised.

That representation buys four things for free:

- **Positivity and fixed total flux.** A histogram is non-negative and sums to
  `nelements` by construction. No positivity constraint, no flux normalisation
  term, no way for either to be violated.
- **An `O(nuv)` update per move.** Moving one quantum from pixel `p_old` to
  `p_new` changes the model visibility by one column difference of the Fourier
  operator. No transform is recomputed.
- **Non-convex regularizers.** No gradient is ever taken, so L0 (a literal count
  of nonzero pixels) is available — true compressed sensing rather than its
  convex L1 surrogate.
- **A posterior**, not a single MAP point: a mean image over post-burn-in
  samples, and a spread across chains.

The cost is that it is a sampler, so it needs sweeps. See
[How many sweeps?](@ref) below — it is the parameter that matters most.

## Basic use

The geometry comes from the Fourier operator, exactly as in
[`reconstruct_bsmem`](@ref) — there are no `nx`/`pixsize` keywords:

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

!!! note "Use `mode=\"dft\"`"
    The sampler needs `O(nuv)` access to a **single column** of the transform per
    move, which a transform-based operator cannot provide. An NFFT plan is
    accepted — `nx` and `pixsize` are recovered from it as `reconstruct_bsmem`
    does — but a DFT matrix is then built internally anyway, so you may as well
    build it once yourself.

### Memory

The DFT matrix is `nuv · nx² · 2 · sizeof(T)` bytes:

| | nx=64 | nx=128 | nx=256 |
|---|---|---|---|
| nuv=500, `Float64` | 33 MB | 131 MB | 524 MB |
| nuv=2000, `Float64` | 131 MB | 524 MB | 2.1 GB |
| same, `Float32` | ÷2 | ÷2 | ÷2 |

Reading with `T=Float32` halves this and is roughly 2× faster on the inner update
(it is memory-bandwidth bound), converging to the same χ²r. Use `Float64` while
developing, then measure.

## Where the chain starts

Omitting `x_start`, as above, gives SQUEEZE's own default: **a point source at the
centre of the grid**, every quantum on one pixel. The annealing is expected to
spread the flux out from there. Pass it explicitly to choose:

```julia
reconstruct_squeeze(:point_source, data, ft; …)   # the default: central Dirac (C without -i)
reconstruct_squeeze(:random,       data, ft; …)   # scattered uniformly over the grid
reconstruct_squeeze(previous_img,  data, ft; …)   # digitise an image (C's -i)
```

Strings work too (`"point_source"`, `"random"`).

An image start is **deterministic** — it is digitised by accumulating
`pixel · nelements` and emitting a quantum each time the total passes 1 — so all
chains begin identically, as they do in C. `:random` is the one to reach for when
you want several chains to explore independently.

## Prior images and masks

`prior_image` (C's `-p`) takes an `nx × nx` map of per-pixel prior
*probabilities*. It becomes an additive penalty of `-log(p)` per quantum, and any
pixel with `p <= 0` gets `1e12` — which makes it a **hard mask**: no quantum can
ever sit there.

```julia
c = (nx + 1) / 2
mask = [sqrt((i-c)^2 + (j-c)^2) <= 20 ? 1.0 : 0.0 for i in 1:nx, j in 1:nx]

img, diag = reconstruct_squeeze(data, ft; niter = 300, prior_image = mask)
@assert sum(img[mask .== 0]) == 0    # exactly zero, by construction
```

The weight defaults to 1 and can be set with `["priorimage", λ]` in
`regularizers`. Supplying a prior image disables auto-centering, since the prior
already fixes the position — the C code prints a warning for that combination.

Tightening a mask past the source degrades the fit, which is a useful way to
bound the emission's extent. On the 2004 data (source ~5 mas across):

| mask radius | χ²r |
|---|---|
| 12 px (2.4 mas) | 739 |
| 18 px (3.6 mas) | 1.11 |
| 25 px (5.0 mas) | 0.750 |
| none | 0.759 |

## Output

`reconstruct_squeeze` returns `(image, diagnostics)`. The image is the mean over
post-burn-in samples of the best chain, normalised to unit total flux. The
diagnostics named tuple carries, per chain:

| Field | Meaning |
|-------|---------|
| `chi2r_mean` | reduced χ² of that chain's **mean image** — chains are ranked on this |
| `chi2r_last` | reduced χ² of the last sample |
| `acceptance` | overall acceptance rate |
| `temperature` | final annealing temperature |
| `nsamples`, `burn_in` | number of accumulated samples, and the sweep they started |
| `regs_final` | each regularizer's value on the final histogram |
| `params`, `stepsize` | fitted SPARCO parameters (empty without a model) |
| `best_chain`, `images`, `ndf` | index of the best chain, all chain images, degrees of freedom |

Pass `outfile = "result.fits"` to write the best chain's mean image, and
`pixsize` alongside it so the CDELT/CRPIX WCS header is written.

## How many sweeps?

Each sweep is `nelements` proposals. **`niter` matters more than any other
setting.** The annealing schedule has a fixed point at `chi2r = chi2_temp · T`,
and a chain that has not yet cooled past it looks converged while sitting at a
large χ²r. On the 2004 Beauty Contest data (64 px @ 0.2 mas):

| `niter` | chains stalled at χ²r > 2 |
|---------|---------------------------|
| 100 | 1 in 8 |
| 250 | 0 in 8 |
| 500 | 0 in 8 |

Use at least a few hundred sweeps, and run several chains — they are independent
restarts and cost nothing but threads. If a regularized run looks badly fitted,
raise `niter` before suspecting the regularizer.

## Monitoring

### Text diagnostics

`print_every = n` prints a per-sweep diagnostic line in the same format as the C
code:

```
Chain: 1 lPost:  1132.6 lPrior:   481.9 lLike:   650.7 V2: 0.83 T3P: 0.61 \
l0:417.00 tv:64.89 E:  1230 MPr: 0.61 T:  1.00 B:19 Iter:   60 of   60
```

`lPost`/`lPrior`/`lLike` are the un-reduced log-posterior terms; `V2`, `T3A`,
`T3P` are **per-observable** reduced χ² (divided by their own counts, not by
`ndf`); regularizer entries show `λ·value`, their contribution to `lPrior`; `E`
is the element count, `MPr` the acceptance EWMA, `T` the temperature, `B` the
burn-in sweep.

### Live displays

Both drivers can redraw the image and a trace panel while they run:

```julia
reconstruct_squeeze(data, ft; niter = 400, pixsize = 0.2, monitor = 50)
reconstruct_squeeze_tempered(data, ft; n_rounds = 11, pixsize = 0.2, monitor = 50)
```

`monitor = n` updates every `n` sweeps (annealing) or every `n` cumulative scans
(tempering); `0`, the default, disables it. Two windows appear: the current image
via [`imdisp`](@ref), and traces of χ²r, the temperature (or Pigeons round) and
every free SPARCO parameter against iteration. Passing `pixsize` gives the image
its mas axes; `monitor_colormap` picks the colour map.

Cost: nothing when off (one integer comparison per sweep) and about **40 ms per
update** when on — so a 400-sweep run at `monitor = 50` costs roughly 33% extra,
and `monitor = 100` about half that. Choose the interval so updates stay rare
compared with the run.

!!! warning "Monitoring forces serial execution"
    The Python plotting bridge is not thread-safe, so anything that draws must
    run on the main thread. Rather than risk a segfault, `reconstruct_squeeze`
    runs its chains serially when `monitor > 0`, and
    `reconstruct_squeeze_tempered` sets `multithreaded = false`. Both log a
    message saying so. Use `monitor = 0` for the threaded path.

## Threading and reproducibility

Chains run one per task over `Threads.@spawn`; start Julia with `-t 4` or more.
Each chain's random stream is a pure function of `(seed, chain index)`, so a
given `seed` reproduces byte-identical results regardless of how the threads are
scheduled.

## SPARCO

[`SqueezeSparco`](@ref) adds a chromatic parametric model in front of the image
(Kluska et al. 2012), plus an over-resolved background:

```
        f_star (λ/λ0)^-4 V_star  +  f_bg (λ/λ0)^bg_indx · 0  +  f_env V_env
V_tot = ───────────────────────────────────────────────────────────────────
        f_star (λ/λ0)^-4         +  f_bg (λ/λ0)^bg_indx      +  f_env
```

The star is a uniform disc (`ud = 0` for a point source); the background is
over-resolved, so it contributes flux but no visibility and only dilutes. The
reconstructed image is the *environment*, and it is **grey** — all the chromatism
lives in the model. That is why one `OIdata` spanning several wavelengths is
enough: the chromatism enters per uv point through `uv_lam`, not through a
polychromatic image cube.

```julia
data = readoifits("2019_v1295Aql.WL_SMOOTH.A.oifits"; filter_bad_data = true)
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

`free` names the parameters the sampler may vary (by name or index); everything
else is held fixed. On V1295 Aql this takes χ²r from ~69 (image alone) to ~0.9.

Centering is enabled automatically when there is no model, and disabled
automatically when there is one — the model already pins the image against the
translation degeneracy. This mirrors the C code; `auto_centering = false` opts
out.

!!! warning "SPARCO needs chromatic leverage"
    On **monochromatic** data a central point star is degenerate with a central
    component of the image, and `f_star` is not identifiable — in a controlled
    test, fixing `f_star = 0.2` gave a *lower* χ²r than the true 0.4. With several
    wavelengths the χ² profile minimises at the truth. If your data is a single
    channel, fix the model parameters rather than fitting them.

!!! warning "Annealed parameters are not posterior samples"
    Under annealing the parametric proposal is the C code's two-point lattice,
    whose spacing is rescaled on every proposal. The move is therefore
    irreversible, and the spread of fitted parameters across chains is a spread
    of optimizer endpoints, not an uncertainty. Use tempering (below) for a real
    posterior, or quote errors from [`bootstrap_fit`](@ref) / [`fit_model`](@ref).

### Parameters that will not move

Flux-like parameters (`f_star`, `f_bg`, the spectral indices) travel well, because
changing them rescales the image's flux fraction analytically within the same
proposal. **Shape parameters — `ud` above all — often will not move at all.**

The reason is coupling, not step size. Once the image has adapted itself to the
current `ud`, changing `ud` with the image held fixed is catastrophic: on V1295
Aql, a 0.05 mas nudge costs ΔlogL ≈ 2.8 × 10⁴, an acceptance probability of
exactly zero. The adaptation then shrinks the step (measured: 0.5 → 8.7 × 10⁻³)
without ever recovering, and the parameter sits where it started.

So do not expect a global search over shape parameters under annealing. Instead:

1. **Fix `ud` and scan it.** Run chains at a grid of fixed values and compare
   `chi2r_mean` — chains are independent, so this parallelises for free.
2. **Seed from a model fit.** Use [`fit_model`](@ref) or
   [`optimize_sparco_flat_parameters`](@ref) to get near the optimum first, then
   let `reconstruct_squeeze` refine locally.
3. Watch the per-parameter `MPr` in the `print_every` output. If it sits near
   zero, that parameter is stuck and its value is just its starting value.
4. **Or temper instead** — see [SPARCO under tempering](@ref) below. That is the
   principled fix rather than a workaround.

## Parallel tempering

`reconstruct_squeeze` **anneals**: it is an optimizer that floors at T = 1 and
then samples. `reconstruct_squeeze_tempered` instead runs a full temperature
ladder with [Pigeons.jl](https://pigeons.run), which owns the schedule, the
non-reversible swaps and their adaptation, and returns a log normalising constant.

Pigeons is a **weak dependency** — install it yourself and load it:

```julia
using OITOOLS, Pigeons          # loading Pigeons enables the method

img, diag = reconstruct_squeeze_tempered(data, ft; n_rounds = 11, n_chains = 30)

diag.logZ            # log normalising constant, by stepping stone
diag.global_barrier  # Λ: see below
diag.schedule        # the adapted β ladder
```

Everything `reconstruct_squeeze` accepts for the model — regularizers,
`prior_image`, `weights`, `nelements` — works here too. Round count is
`n_rounds`, and round *i* performs 2ⁱ scans, so cost grows geometrically; only
the **final** round's samples are kept, the earlier ones being tuning.

### The path, and how it differs from the C code

The tempered family is

```
π_β(x) ∝ exp( β · ( −lLikelihood(x) − lPrior(x) ) )
```

so the reference at β = 0 is the **uniform distribution over element
configurations**, which can be sampled exactly — every quantum placed uniformly
at random. That exactness matters: parallel tempering's guarantees rest on the
reference being i.i.d.-samplable.

Two deliberate differences from the C annealing:

- It tempers the likelihood **and** the prior together, where C's temperature
  scales only the likelihood. Only β = 1 needs to be the correct posterior; the
  path between is a design choice.
- `COPYCAT` moves are **not used**. They land on a random donor element's pixel,
  so `P(propose o→p) ∝ count(p)` — asymmetric, and irreversible once the last
  quantum leaves a pixel. That is an effective optimizer heuristic and an invalid
  sampler move. The tempered explorer is a fixed mixture of the symmetric
  SMALL/MEDIUM/ANYWHERE moves, with no adaptation, so it is a valid β-invariant
  kernel without any vanishing-adaptation argument.

### Reading Λ, the global communication barrier

`diag.global_barrier` is the single number that tells you whether the run was
adequate: round trips between β = 0 and β = 1 take roughly Λ² scans, and you want
`n_chains` of order Λ.

Λ is a property of the *path*, not of the run — but it is **estimated from the
swap rates, so too few chains bias it low**. On the 2004 data at 64 px:

| `n_chains` | Λ reported | χ²r |
|---|---|---|
| 8 | 6.96 | 0.79 |
| 16 | 14.9 | 1.01 |
| 30 | 27.1 | 0.77 |

The barrier is genuinely ≈ 27 here; the smaller runs simply could not resolve it.
Going from a uniform random image to a sharp posterior over a 4096-dimensional
discrete space is a long way, so expect Λ of a few tens and budget chains
accordingly. If `logZ` shifts materially when you add chains, it has not
converged.

### SPARCO under tempering

`reconstruct_squeeze_tempered` accepts the same `model` keyword, and it solves the
problem above: at low β the likelihood is weak, so the image is free and the
parameter can move; swaps then carry that information up to β = 1. Measured on
synthetic V1295 Aql data with a true `ud` of 2.0 mas:

| started at | annealed | tempered |
|---|---|---|
| `ud = 4.0` | **3.73** (stuck) | **2.15** |
| `ud = 0.5` | 1.93 | 2.17 |

Annealing gives two answers differing by a factor of two depending on where it
started; tempering gives the same answer from either side. `f_star` likewise
recovers to 0.4697 and 0.4698 from starts of 0.15 and 0.80.

Two things change under tempering, both required for the kernel to be valid and
both prescribed by the C code:

- the parameter proposal is a **symmetric Gaussian** random walk, not the
  two-point `p ± σ` lattice (which is irreversible);
- the stepsize is **frozen** — no adaptation, so no vanishing-adaptation argument
  is needed.

`diag.params` is then a genuine **posterior mean** over the β = 1 visits, not an
optimizer endpoint. (`diag.params_by_replica` gives each replica's current values,
but most replicas sit at low β where the parameters are near-draws from the
reference box — do not average those.)

!!! note "Free parameters need a bounded reference box"
    PT needs a reference it can sample i.i.d., so the free parameters need a
    *proper* uniform prior. SPARCO's own box priors leave `f_star` and `ud`
    unbounded above, so `reconstruct_squeeze_tempered` uses a box —
    `default_sparco_bounds`, generous by default — and rejects proposals leaving
    it. Override with `param_bounds = (lo, hi)`. A box that is too wide wastes
    chains; one that excludes the truth biases the answer.

### When to use which

`reconstruct_squeeze` is much cheaper and is the right tool for getting an image.
Reach for the tempered version when you want a defensible posterior mean, an
evidence estimate for model comparison, or when you suspect the annealing is
getting trapped in a mode. On the 2004 data the two agree in morphology, and the
tempered mean is visibly smoother — as a genuine posterior mean over thousands of
samples should be.

## Relation to the other reconstructors

| | Method | Gives you |
|---|---|---|
| [`reconstruct`](@ref) | gradient descent (VMLMB) | MAP point, differentiable regularizers |
| [`reconstruct_bsmem`](@ref) | maximum entropy | MAP point, entropy prior |
| [`reconstruct_bsdmm`](@ref) | ADMM / proximal | MAP point, non-smooth regularizers |
| `reconstruct_squeeze` | MCMC + annealing | posterior mean, non-convex regularizers |
| `reconstruct_squeeze_tempered` | MCMC + parallel tempering (Pigeons) | posterior mean **and** `logZ` |

Use `reconstruct_squeeze` when you want L0 or another non-convex prior, when
positivity and flux conservation must be exact, or when you want a posterior
rather than a point estimate. Use the gradient reconstructors when you want speed
and a differentiable regularizer is adequate.

## Not yet implemented

Wavelet regularizers, polychromatic image cubes, MEDIAN/MODE images, the
full-chain history cube, and bandwidth smearing in the Fourier kernel.
