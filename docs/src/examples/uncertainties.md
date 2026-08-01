# Uncertainties

How much can you trust a fitted parameter?  OITOOLS offers two answers, and
they do not measure the same thing:

| Method | What it assumes | What it can see |
|--------|-----------------|-----------------|
| `fit_model_lsqfit` covariance | error bars correct and uncorrelated; model locally linear | nothing beyond the quoted errors, rescaled so that χ²ᵣ = 1 |
| [`bootstrap_fit`](@ref) | blocks of data are numerous and exchangeable | correlated calibration errors, mis-stated error bars, a bad night or a bad baseline |

This page reports a simulation study that measures how well each of them
actually performs, and compares them with
[PMOIRED](https://github.com/amerand/PMOIRED) on identical data.  The code is
in `demos/bootstrap_validation`; `RESULTS.md` there carries the same tables.

## How the calibration is measured

An uncertainty is a claim about a distribution you never see: the spread of
best-fit values you *would* get if you repeated the observation.  So the study
repeats it.  120 independent noise realisations ("universes") of the same
observation are generated, each is fitted, and what each method claims about a
single universe is compared with what the ensemble actually does:

```
sigma_true    = std over universes of the best-fit value      <- the truth
sigma_report  = mean over universes of the quoted 1 sigma     <- the claim
ratio         = sigma_report / sigma_true                     <- 1.00 = calibrated
coverage      = fraction of universes with |best - truth| < sigma_report
                                                              <- 0.683 = calibrated
```

A ratio below 1 means overconfident error bars.  `sigma_true` is itself
estimated from 120 universes, so it carries about 6.5% noise — but that noise
is *common* to every method within a regime, so the differences between methods
in a column are far better determined than 6.5%.

The model is a binary at CHARA/MIRC H-band resolution (a partially resolved
primary plus an unresolved companion at 1.44 mas with 5% of the flux), fitted
for four parameters against V² and T3PHI.  The uv coverage, MJDs, telescope
configurations and error bars come from a real dataset; only the observable
values are replaced.

Five regimes stress different assumptions:

| Regime | What is simulated | What it tests |
|--------|-------------------|---------------|
| `ideal` | independent Gaussian noise at exactly the quoted error bars | control: everything should be calibrated |
| `systematics` | same noise, plus a calibration offset per (epoch, baseline/triangle) shared by all wavelengths: 5% on V², 1.5° on T3PHI | the case block bootstrapping exists for |
| `underestimated` | noise generated at 2× the quoted error bars | whether χ²ᵣ rescaling and resampling recover the true scatter |
| `fewblocks` | ideal noise, but a single snapshot: 18 blocks instead of 180 | where a block bootstrap runs out of blocks |
| `mismatch` | the truth contains a third component that nobody fits | no method can be right; how do they fail? |

## Choosing a resampling scheme

`bootstrap_fit(...; mode=...)` selects how block multiplicities are drawn.
Wall times are per replicate, single-threaded, on the benchmark dataset (1440
data points, 180 blocks, 4 free parameters); reproduce with
`demos/bootstrap_validation/benchmark.jl`.

| Mode | Scheme | 180 blocks | 18 blocks | ms / replicate |
|------|--------|-----------|-----------|----------------|
| `:replacement` (default) | multinomial — the textbook bootstrap | 0.87–1.01 | **0.96–1.14** | 3.4 |
| `:halfsample` | balanced repeated replication: one random half | **0.94–1.06** | 1.24–1.56 | **2.8** |
| `:weights` | multiplier (Bayesian) bootstrap: continuous block weights applied by scaling the error bars | 0.78–0.99 | 0.71–0.85 | 3.1 |
| `:pmoired` | PMOIRED's two independent half-samples, fitted jointly | 0.62–0.74 ⚠ | 0.67–0.79 ⚠ | 3.2 |

The block columns give the calibration ratio; bold marks the best in each.

!!! warning "`:pmoired` is biased low by √2 — comparison only"
    `:pmoired` is not an alternative estimator, it is a reproduction of
    PMOIRED's construction, kept so that a PMOIRED result can be cross-checked
    inside OITOOLS.  Its block multiplicities have variance ½ instead of 1, so
    it returns error bars that are about √2 too small: 0.62–0.79 of the true
    scatter in every regime tested, with a 1σ coverage near 0.50 instead of
    0.683.  Never quote an uncertainty obtained with it.  `bootstrap_fit` emits
    a warning the first time it is used in a session.

- **`:halfsample`** is the best-calibrated scheme wherever blocks are numerous,
  and the cheapest, because each replicate fits half the data.  Its failure mode
  with few blocks is over-conservatism, which is the safe direction.
- **`:replacement`** is the default: the textbook bootstrap, the most uniform
  across regimes, never off by more than 14% including the single-snapshot case.
- **`:weights`** is consistently the most optimistic of the three, and clearly
  so when blocks are scarce.  It is included for completeness.
- **`:pmoired`** exists only to reproduce PMOIRED results, and is biased low by
  about √2 (see the warning above and the tables below).

A replicate costs about the same as the original fit: 3.4 ms against 3.8 ms for
a single Levenberg-Marquardt fit, of which only 0.85 ms is the resampling —
the rest is the fit itself, so the scheme you choose barely matters for
runtime.  `bootstrap_fit` spreads replicates over all available Julia threads;
500 replicates of this problem take about 0.6 s on 16 threads.

Granularity matters as much as the scheme.  `:config` (the default, one block
per MJD and telescope configuration) works from 18 blocks upward; `:epoch`
becomes erratic when there are few epochs (0.92–1.83 on a single snapshot);
`:point` destroys the correlation structure and underestimates by design.

## A note on MJD precision

Blocking is keyed on the MJD, and `Float32` cannot resolve it: near MJD 55000 its
representable values are 3.9 × 10⁻³ d apart — **5.6 minutes** — so exposures taken
minutes apart would collapse onto one value and their blocks would silently merge.

`readoifits` therefore stores `v2_mjd`, `t3_mjd`, `vis_mjd` and `flux_mjd` as
`Float64` regardless of `T`, so the block structure is identical whether you read
with the default `T=Float32` or with `T=Float64`. (`uv_mjd` follows `T` — it is
the array handed to the model evaluator as `$MJD`, where the extra precision buys
nothing.) Nothing to do on your side; the numbers below are unaffected by the
precision you choose.

## Results

`ratio = quoted sigma / true scatter`, one column per parameter; the last
column is the coverage averaged over the four parameters (0.683 when
calibrated).


### Regime: `ideal`

| method | A,ud | B,f | B,x | B,y | coverage |
|---|---|---|---|---|---|
| OITOOLS  analytic (LM) | 0.96 | 0.98 | 0.98 | 1.04 | 0.66 |
| OITOOLS  parametric MC | 0.91 | 0.96 | 0.96 | 1.03 | 0.63 |
| OITOOLS  block bootstrap | 0.89 | 0.95 | 0.95 | 1.01 | 0.64 |
| OITOOLS  epoch bootstrap | 0.90 | 0.97 | 0.92 | 1.01 | 0.62 |
| OITOOLS  weighted bootstrap | 0.88 | 0.95 | 0.87 | 0.99 | 0.61 |
| OITOOLS  half-sample (BRR) | 0.94 | 0.98 | 1.01 | 1.06 | 0.66 |
| OITOOLS  PMOIRED scheme | 0.65 | 0.70 | 0.67 | 0.74 | 0.49 |
| PMOIRED  analytic (dpfit) | 0.96 | 0.98 | 0.98 | 1.04 | 0.65 |
| PMOIRED  bootstrapFit | 0.64 | 0.69 | 0.66 | 0.74 | 0.47 |
| PMOIRED  bootstrapFit (no start jitter) | 0.65 | 0.69 | 0.67 | 0.74 | 0.48 |

### Regime: `systematics`

| method | A,ud | B,f | B,x | B,y | coverage |
|---|---|---|---|---|---|
| OITOOLS  analytic (LM) | 0.18 | 0.34 | 0.27 | 0.31 | 0.22 |
| OITOOLS  parametric MC | 0.04 | 0.09 | 0.07 | 0.08 | 0.06 |
| OITOOLS  block bootstrap | 0.88 | 0.92 | 0.87 | 0.96 | 0.61 |
| OITOOLS  epoch bootstrap | 0.89 | 0.90 | 0.85 | 0.94 | 0.61 |
| OITOOLS  weighted bootstrap | 0.82 | 0.90 | 0.78 | 0.92 | 0.58 |
| OITOOLS  half-sample (BRR) | 0.99 | 0.98 | 0.94 | 1.00 | 0.65 |
| OITOOLS  PMOIRED scheme | 0.68 | 0.69 | 0.62 | 0.70 | 0.50 |
| PMOIRED  analytic (dpfit) | 0.18 | 0.34 | 0.27 | 0.31 | 0.21 |
| PMOIRED  bootstrapFit | 0.66 | 0.67 | 0.63 | 0.70 | 0.50 |
| PMOIRED  bootstrapFit (no start jitter) | 0.67 | 0.68 | 0.62 | 0.70 | 0.49 |

### Regime: `underestimated`

| method | A,ud | B,f | B,x | B,y | coverage |
|---|---|---|---|---|---|
| OITOOLS  analytic (LM) | 1.01 | 0.96 | 1.00 | 0.95 | 0.68 |
| OITOOLS  parametric MC | 0.48 | 0.47 | 0.49 | 0.47 | 0.37 |
| OITOOLS  block bootstrap | 0.96 | 0.91 | 0.97 | 0.90 | 0.64 |
| OITOOLS  epoch bootstrap | 0.94 | 0.93 | 0.93 | 0.90 | 0.63 |
| OITOOLS  weighted bootstrap | 0.92 | 0.90 | 0.86 | 0.90 | 0.63 |
| OITOOLS  half-sample (BRR) | 1.00 | 0.94 | 1.03 | 0.95 | 0.65 |
| OITOOLS  PMOIRED scheme | 0.68 | 0.68 | 0.67 | 0.66 | 0.51 |
| PMOIRED  analytic (dpfit) | 1.01 | 0.96 | 0.99 | 0.95 | 0.67 |
| PMOIRED  bootstrapFit | 0.67 | 0.65 | 0.67 | 0.65 | 0.49 |
| PMOIRED  bootstrapFit (no start jitter) | 0.67 | 0.67 | 0.68 | 0.65 | 0.50 |

### Regime: `fewblocks`

| method | A,ud | B,f | B,x | B,y | coverage |
|---|---|---|---|---|---|
| OITOOLS  analytic (LM) | 1.05 | 1.02 | 0.88 | 1.02 | 0.70 |
| OITOOLS  parametric MC | 1.01 | 0.99 | 0.86 | 1.00 | 0.69 |
| OITOOLS  block bootstrap | 0.96 | 1.06 | 0.97 | 1.14 | 0.69 |
| OITOOLS  epoch bootstrap | 0.92 | 1.73 | 1.83 | 1.62 | 0.70 |
| OITOOLS  weighted bootstrap | 0.81 | 0.85 | 0.71 | 0.75 | 0.57 |
| OITOOLS  half-sample (BRR) | 1.33 | 1.26 | 1.24 | 1.56 | 0.79 |
| OITOOLS  PMOIRED scheme | 0.70 | 0.76 | 0.67 | 0.79 | 0.52 |
| PMOIRED  analytic (dpfit) | 1.03 | 1.00 | 0.87 | 1.00 | 0.70 |
| PMOIRED  bootstrapFit | 0.70 | 0.72 | 0.64 | 0.71 | 0.54 |
| PMOIRED  bootstrapFit (no start jitter) | 0.69 | 0.71 | 0.64 | 0.70 | 0.52 |

### Regime: `mismatch`

| method | A,ud | B,f | B,x | B,y | coverage |
|---|---|---|---|---|---|
| OITOOLS  analytic (LM) | 2.76 | 2.30 | 2.77 | 2.65 | 0.05 |
| OITOOLS  parametric MC | 1.02 | 0.87 | 1.02 | 1.03 | 0.01 |
| OITOOLS  block bootstrap | 11.11 | 5.80 | 4.79 | 6.01 | 0.25 |
| OITOOLS  epoch bootstrap | 11.61 | 2.01 | 2.49 | 2.44 | 0.04 |
| OITOOLS  weighted bootstrap | 10.22 | 5.57 | 5.44 | 6.83 | 0.25 |
| OITOOLS  half-sample (BRR) | 10.11 | 5.38 | 6.22 | 6.63 | 0.25 |
| OITOOLS  PMOIRED scheme | 6.27 | 3.96 | 4.15 | 4.63 | 0.23 |
| PMOIRED  analytic (dpfit) | 2.76 | 2.30 | 2.71 | 2.64 | 0.05 |
| PMOIRED  bootstrapFit | 6.76 | 4.10 | 3.62 | 4.60 | 0.22 |
| PMOIRED  bootstrapFit (no start jitter) | 6.94 | 4.08 | 3.66 | 4.55 | 0.22 |

## What the numbers say

1. **The two packages' analytic errors are the same number.** `OITOOLS analytic
   (LM)` and `PMOIRED analytic (dpfit)` agree to two decimals in every regime.
   Both rescale the covariance to χ²ᵣ = 1, so the two implementations are
   cross-validated against each other.

2. **PMOIRED's bootstrap is overconfident by a factor 1/√2, in every regime.**
   `PMOIRED bootstrapFit` lands at 0.64–0.74 everywhere, and OITOOLS' emulation
   of the same scheme (`mode=:pmoired`) reproduces it to ±0.02 in all 20
   (regime, parameter) cells.  The deficit is the *scheme* — two independent
   half-samples give block multiplicities with variance ½ instead of 1 — not an
   implementation detail.  Coverage is ~0.50 instead of 0.683.

3. **Only resampling sees correlated systematics.** With a 5% per-(epoch,
   baseline) calibration error the analytic covariance is 3–5× too small *even
   after* its χ²ᵣ rescaling, which is itself a factor 3.8 here.

4. **χ²ᵣ rescaling handles pure underestimation.** With noise at 2× the quoted
   errors, both analytic estimates come back calibrated (0.95–1.01).

5. **Adding noise to the data is not a bootstrap.** Refitting
   [`perturb_data`](@ref) replicates ("parametric Monte Carlo") is dominated
   everywhere: identical to the analytic covariance in the well-behaved regimes,
   exactly half in `underestimated` — the factor it was lied to about — and
   11–25× too small under systematics, because by construction it believes the
   error bars.  This is why it is not offered as an uncertainty API.

6. **No uncertainty is meaningful when the model is wrong.** With an unmodelled
   third component the best fit is biased by 10–36 true σ; the bootstraps
   inflate their error bars by 5–11× (they do sense the unmodelled structure)
   and coverage still collapses to 0.22–0.25, the analytic errors to 0.05.
   Resampling measures precision, never accuracy — a small error bar on a bad
   model is still a bad answer.

## Cross-package performance

Same file, same model, same number of replicates:

| | single fit | bootstrap, 1 process/thread | bootstrap, all cores |
|---|---|---|---|
| OITOOLS `bootstrap_fit` | 3.8 ms | 3.4 ms / replicate | 1.2 ms / replicate (16 threads) |
| PMOIRED `bootstrapFit` | 44 ms | 78 ms / replicate | 30 ms / replicate (24 processes) |

PMOIRED's replicate costs more than its own single fit because each replicate
fits a dataset built from two copies; its multi-process scaling is limited by
process start-up and by pickling the data to the workers.

## Reproducing

```bash
cd demos/bootstrap_validation
OUT=/tmp/bootstrap_validation

julia --project=../.. make_universes.jl \$OUT 120   # simulate
python3 sanitize.py \$OUT                           # make files readable by PMOIRED
julia --project=../.. -t auto run_oitools.jl \$OUT 100
PMOIRED_MAX_THREADS=12 python3 run_pmoired.py \$OUT 100
python3 analyze.py \$OUT                            # the tables above
julia --project=../.. benchmark.jl \$OUT            # the wall times above
```

See `demos/bootstrap_validation/README.md` for the details, including the two
interoperability fixes needed to make PMOIRED read OITOOLS-written OIFITS files.
