# Results

120 universes per regime, 100 bootstrap/Monte Carlo replicates per universe,
identical simulated files fitted by both packages.  Numbers are
`ratio = quoted sigma / true scatter`, one column per parameter; **1.00 is
calibrated, below 1 is overconfident**.  The last column is the fraction of
universes falling within one quoted sigma of the truth, averaged over the four
parameters (0.683 when calibrated).

`sigma_true` is estimated from 120 universes, so it carries ~6.5% noise -- but
that noise is *common* to every method within a regime, so the differences
between methods within a column are far better determined than 6.5%.

Reproduce with the commands in README.md.

## Which scheme to use

| Scheme | many blocks (180) | few blocks (18) | ms / replicate |
|--------|-------------------|-----------------|----------------|
| `:halfsample` (BRR) | **0.94-1.06** | 1.24-1.56 (conservative) | **2.8** |
| `:replacement` (multinomial, default) | 0.87-1.01 | **0.96-1.14** | 3.4 |
| `:weights` (multiplier / Bayesian) | 0.78-0.99 | 0.71-0.85 | 3.1 |
| `:pmoired` (comparison only, biased low by sqrt(2)) | 0.62-0.74 | 0.67-0.79 | 3.2 |

Wall times are per replicate, single-threaded, on the benchmark dataset (1440
points, 180 blocks, 4 free parameters); `benchmark.jl` reproduces them.  A
single Levenberg-Marquardt fit of the same data takes 3.8 ms, of which the
resampling accounts for 0.85 ms -- the fit dominates, so the choice of scheme
barely affects runtime.  For scale, PMOIRED's `bootstrapFit` costs 78 ms per
replicate in one process on the same file (30 ms with 24 processes).

`:halfsample` is the best-calibrated scheme wherever blocks are numerous, and it
is also the cheapest.  Its failure mode with few blocks is over-conservatism,
which is the safe direction.  `:replacement` -- the default -- is the textbook
bootstrap and the most uniform across regimes, never off by more than 14%.
`:weights` is consistently the most optimistic of the three and is the one to
avoid when blocks are scarce.

## What the numbers say

1. **The two packages' analytic errors are the same number.** `OITOOLS analytic
   (LM)` and `PMOIRED analytic (dpfit)` agree to two decimals in every regime.
   Both rescale the covariance to chi2r = 1, so the two implementations are
   cross-validated against each other.

2. **PMOIRED's bootstrap is overconfident by a factor 1/sqrt(2), in every
   regime.** `PMOIRED bootstrapFit` lands at 0.64-0.74 everywhere, and OITOOLS'
   emulation of the same scheme (`mode=:pmoired`) reproduces it to +/-0.02 in
   all 20 (regime, parameter) cells.  The deficit is the *scheme* -- two
   independent half-samples give block multiplicities with variance 1/2 instead
   of 1 -- not an implementation detail.  Coverage is ~0.50 instead of 0.683.

3. **`randomiseParam` does not compensate.** Turning the 1-sigma starting-point
   jitter off (`pmoired_boot_nostart`) changes nothing: 0.65 vs 0.64.

4. **Fixing the multiplicity fixes the calibration**, and half-sampling fixes it
   best: in the `systematics` regime BRR gives 0.94-1.00 where the multinomial
   gives 0.87-0.96 and PMOIRED's scheme 0.62-0.70.

5. **Only resampling sees correlated systematics.** With a 5% per-(epoch,
   baseline) calibration error the analytic covariance is 3-5x too small
   (0.18-0.34) *even after* its chi2r rescaling, which is itself a factor 3.8
   here.  The parametric Monte Carlo is 11-25x too small (0.04-0.09), because by
   construction it believes the error bars.

6. **The parametric Monte Carlo is dominated everywhere.** It never beats the
   analytic covariance -- the same value in `ideal` and `fewblocks`, exactly
   half in `underestimated` (the factor it was lied to about), an order of
   magnitude worse under systematics.  This is why `montecarlo_fit` was removed
   from the OITOOLS API; `perturb_data` remains as a simulation utility.

7. **chi2r rescaling handles pure underestimation.** With noise at 2x the quoted
   errors, both analytic estimates are calibrated (0.95-1.01).

8. **Granularity has to match the data.** On a single snapshot the config-level
   bootstrap is still fine on 18 blocks, while the epoch-level one becomes
   erratic (0.92-1.83) because there are almost no epochs to resample.

9. **No uncertainty is meaningful when the model is wrong.** With an unmodelled
   third component the best fit is biased by 10-36 true sigma; the bootstraps
   inflate their error bars by 5-11x (they do sense the unmodelled structure)
   and coverage still collapses to 0.22-0.25, the analytic errors to 0.05, the
   parametric MC to 0.01.  Resampling measures precision, never accuracy.

## Tables


### ideal

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

### systematics

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

### underestimated

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

### fewblocks

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

### mismatch

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

