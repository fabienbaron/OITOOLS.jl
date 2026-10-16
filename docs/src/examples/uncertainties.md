# Uncertainties

| Method | Assumes | Accounts for |
|--------|---------|--------------|
| `fit_model_lsqfit` covariance | error bars correct and uncorrelated; model locally linear | the quoted errors, rescaled so that χ²ᵣ = 1 |
| [`bootstrap_fit`](@ref) | many exchangeable blocks of data | correlated calibration errors, wrong error bars, a bad night or baseline |

The posterior samplers ([`fit_model_nested`](@ref), [`fit_model_ensemble`](@ref),
[`fit_model_pigeons`](@ref)) use the error bars as given; they are not part of the study below.

This page reports a simulation study of the two methods, compared with
[PMOIRED](https://github.com/amerand/PMOIRED) on the same data. The code is in
`demos/bootstrap_validation`, where `RESULTS.md` has the same tables.

## Method

120 independent noise realisations ("universes") of one observation are generated and fitted.
The error bar each method gives on one universe is compared with the scatter of the best fits
over all universes:

```
sigma_true    = std over universes of the best-fit value      <- the truth
sigma_report  = mean over universes of the quoted 1 sigma     <- the claim
ratio         = sigma_report / sigma_true                     <- 1.00 = calibrated
coverage      = fraction of universes with |best - truth| < sigma_report
                                                              <- 0.683 = calibrated
```

A ratio below 1 means the error bars are too small. `sigma_true` from 120 universes has about
6.5% noise; it is the same for all methods in a regime, so differences between methods in a
column are smaller than that.

The model is a binary at CHARA/MIRC H-band resolution (a partially resolved
primary plus an unresolved companion at 1.44 mas with 5% of the flux), fitted
for four parameters to V² and T3PHI. The uv coverage, MJDs, telescope configurations and error
bars come from a real dataset; the observables are replaced.

| Regime | Simulated noise |
|--------|-----------------|
| `ideal` | independent Gaussian noise at the quoted error bars |
| `systematics` | as `ideal`, plus a calibration offset per (epoch, baseline or triangle), common to all wavelengths: 5% on V², 1.5° on T3PHI |
| `underestimated` | noise at 2× the quoted error bars |
| `fewblocks` | as `ideal`, one snapshot: 18 blocks instead of 180 |
| `mismatch` | the truth has a third component that is not fitted |

## Choosing a resampling scheme

`bootstrap_fit(...; mode=...)` sets how the block multiplicities are drawn. Times are per
replicate, single-threaded, on the benchmark dataset (1440 data points, 180 blocks, 4 free
parameters; `demos/bootstrap_validation/benchmark.jl`).

| Mode | Scheme | 180 blocks | 18 blocks | ms / replicate |
|------|--------|-----------|-----------|----------------|
| `:replacement` (default) | multinomial | 0.87–1.01 | **0.96–1.14** | 3.4 |
| `:halfsample` | balanced repeated replication: a random half | **0.94–1.06** | 1.24–1.56 | **2.8** |
| `:weights` | multiplier (Bayesian) bootstrap: continuous block weights, applied to the error bars | 0.78–0.99 | 0.71–0.85 | 3.1 |
| `:pmoired` | PMOIRED's two independent half-samples, fitted jointly | 0.62–0.74 ⚠ | 0.67–0.79 ⚠ | 3.2 |

The block columns give the range of the ratio over the regimes and parameters with that number
of blocks; bold is the best per column.

!!! warning "`:pmoired` gives error bars √2 too small"
    `:pmoired` reproduces PMOIRED's scheme, for comparing with PMOIRED results. Its block
    multiplicities have variance ½ instead of 1, so its error bars are 0.62–0.79 of the true
    scatter in every regime, with 1σ coverage near 0.50. `bootstrap_fit` warns on its first use
    in a session.

- `:halfsample` is the best calibrated with many blocks, and the fastest, since each replicate
  fits half the data. With few blocks its error bars are too large.
- `:replacement` is within 14% in all regimes, including the single snapshot.
- `:weights` gives the smallest error bars of the three, most so with few blocks.

A replicate takes 3.4 ms, against 3.8 ms for one Levenberg-Marquardt fit; 0.85 ms of it is the
resampling. With `threaded = true` (default) replicates run on all Julia threads: 500
replicates take about 0.6 s on 16 threads.

`granularity` sets the blocks. `:config` (default) is one block per MJD and telescope
configuration, and is calibrated from 18 blocks up. `:epoch` (one block per MJD) gives
0.92–1.83 on a single snapshot. `:point` (one block per data point) ignores correlations
between points. `bootstrap_fit` warns when there are fewer than 10 blocks.

## MJD precision

Blocks are keyed on the MJD. `readoifits` stores `v2_mjd`, `t3_mjd`, `vis_mjd`, `flux_mjd` and
`uv_mjd` as `Float64` for any `T` (Float32 values near MJD 55000 are 5.6 minutes apart), so the
blocks are the same with `T = Float32` and `T = Float64`.

## Results

`ratio = quoted sigma / true scatter`, one column per parameter; the last column is the coverage
averaged over the four parameters (0.683 when calibrated).

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

## Summary

1. `OITOOLS analytic (LM)` and `PMOIRED analytic (dpfit)` agree to two decimals in every regime;
   both rescale the covariance to χ²ᵣ = 1.
2. `PMOIRED bootstrapFit` gives 0.64–0.74 in every regime, with coverage near 0.50.
   `mode = :pmoired` reproduces it to ±0.02 in all 20 (regime, parameter) cells.
3. With correlated calibration errors (`systematics`), the analytic covariance is 3–5× too small
   after χ²ᵣ rescaling (a factor 3.8 here). The block bootstraps give 0.78–1.00.
4. With noise at 2× the quoted errors, the analytic errors give 0.95–1.01.
5. Refitting [`perturb_data`](@ref) replicates (parametric Monte Carlo) equals the analytic
   covariance in `ideal` and `fewblocks`, gives half in `underestimated`, and is 11–25× too small
   with systematics. It is not offered as an uncertainty method.
6. With an unmodelled third component the best fit is off by 10–36 true σ. The bootstraps'
   error bars grow 5–11×, and coverage is 0.22–0.25 for the bootstraps and 0.05 for the
   analytic errors.

## Cross-package performance

Same file, same model, same number of replicates:

| | single fit | bootstrap, 1 process/thread | bootstrap, all cores |
|---|---|---|---|
| OITOOLS `bootstrap_fit` | 3.8 ms | 3.4 ms / replicate | 1.2 ms / replicate (16 threads) |
| PMOIRED `bootstrapFit` | 44 ms | 78 ms / replicate | 30 ms / replicate (24 processes) |

A PMOIRED replicate fits a dataset built from two copies of the data.

## Reproducing

```bash
cd demos/bootstrap_validation
OUT=/tmp/bootstrap_validation

julia --project=../.. make_universes.jl $OUT 120   # simulate
python3 sanitize.py $OUT                           # make files readable by PMOIRED
julia --project=../.. -t auto run_oitools.jl $OUT 100
PMOIRED_MAX_THREADS=12 python3 run_pmoired.py $OUT 100
python3 analyze.py $OUT                            # the tables above
julia --project=../.. benchmark.jl $OUT            # the wall times above
```

`demos/bootstrap_validation/README.md` describes the steps, including what `sanitize.py`
changes for PMOIRED to read the files.
