# Are the quoted uncertainties right?

A head-to-head calibration test of the uncertainty estimators in **OITOOLS** and
**PMOIRED**, on identical simulated data.

The idea is the only one that actually settles the question: generate many
independent realisations of the same observation ("universes"), fit each of
them, and compare what each method *claims* about a single universe with what
the ensemble of universes actually *does*.

```
sigma_true    = std over universes of the best-fit value      <- the truth
sigma_report  = mean over universes of the quoted 1 sigma     <- the claim
ratio         = sigma_report / sigma_true                     <- 1.00 = calibrated
coverage      = fraction of universes with |best - truth| < sigma_report
                                                              <- 0.683 = calibrated
bias          = (mean best - truth) / sigma_true
```

## The model

A binary at CHARA/MIRC H-band resolution (λ/2B ≈ 0.47 mas): a partially
resolved primary and an unresolved companion 1.44 mas away with 5% of the flux.
Four free parameters (`A,ud`, `B,f`, `B,x`, `B,y`), fitted against V² and T3PHI.
The uv coverage, MJDs, telescope configurations and error bars are taken from a
real dataset (`demos/data/iota_peg6t.oifits`); only the observable values are
replaced.

**The total flux is pinned to 1** via `A,f = 1 - B,f`. This is not cosmetic:
PMOIRED normalises the model visibility by the total flux and OITOOLS does not,
so the two packages fit the same function *only* when the fluxes sum to 1. With
this constraint, PMOIRED recovers an OITOOLS-generated noiseless dataset with
χ²ᵣ = 0.003 and parameters correct to five digits — the two model
implementations agree, and any difference in the results is a difference
between the *uncertainty* methods, not between the models.

## The regimes

| Regime | What is simulated | What it tests |
|--------|-------------------|---------------|
| `ideal` | independent Gaussian noise at exactly the quoted error bars | control: every method should be calibrated |
| `systematics` | same noise, plus a calibration offset per (epoch, baseline/triangle), shared by all wavelengths: 5% on V², 1.5° on T3PHI | the case block bootstrapping exists for |
| `underestimated` | noise generated at 2× the quoted error bars | whether χ²ᵣ rescaling and resampling recover the true scatter |
| `fewblocks` | ideal noise, but a single snapshot (18 blocks instead of 180) | where a block bootstrap runs out of blocks |
| `mismatch` | the truth contains a third component that nobody fits | no method can be right about the truth; how do they fail? |

## Methods compared

| Key | Method |
|-----|--------|
| `lsqfit` | OITOOLS `fit_model_lsqfit` covariance (LM, rescaled by χ²ᵣ) |
| `montecarlo` | parametric Monte Carlo — refits of `perturb_data` replicates. No longer part of the OITOOLS API (see RESULTS.md); computed inline by `run_oitools.jl` so the evidence stays reproducible |
| `boot_config` | OITOOLS `bootstrap_fit` — block bootstrap over (MJD, config), with replacement |
| `boot_epoch` | OITOOLS `bootstrap_fit` — block bootstrap over whole epochs |
| `boot_weights` | OITOOLS `bootstrap_fit` with `mode=:weights` — multiplier (Bayesian) bootstrap, continuous block weights |
| `boot_halfsample` | OITOOLS `bootstrap_fit` with `mode=:halfsample` — balanced repeated replication, one random half per replicate |
| `boot_pmoired` | OITOOLS `bootstrap_fit` with `mode=:pmoired` — PMOIRED's two-half-samples scheme on the same blocks. Biased low by ~sqrt(2); present to isolate the scheme from the implementation, never to quote from |
| `pmoired_fit` | PMOIRED `doFit` covariance (`dpfit`, `normalizedUncer=True`) |
| `pmoired_boot` | PMOIRED `bootstrapFit` with its defaults |
| `pmoired_boot_nostart` | same with `randomiseParam=False`, isolating the resampling scatter from the starting-point jitter |

## Running it

```bash
cd demos/bootstrap_validation
OUT=/tmp/bootstrap_validation

# 1. simulate (a few minutes; ~400 kB per universe)
julia --project=../.. make_universes.jl $OUT 120

# 2. make the files readable by PMOIRED (see "Interoperability" below)
python3 sanitize.py $OUT

# 3. fit everything — these two are independent and can run concurrently
julia --project=../.. -t auto run_oitools.jl $OUT 100
PMOIRED_MAX_THREADS=12 python3 run_pmoired.py $OUT 100

# 4. report
python3 analyze.py $OUT            # add --csv or --md for machine-readable output

# 5. wall time per replicate for each scheme
julia --project=../.. benchmark.jl $OUT
```

The second argument of the runners is the number of bootstrap/Monte Carlo
replicates per universe (100 is plenty: the σ of a bootstrap with B replicates
is itself uncertain by ~1/√(2B) = 7%, and that averages away over the
universes).

`PMOIRED_DIR` overrides where PMOIRED is imported from (default
`~/SOFTWARE/PMOIRED`).

## Interoperability notes

Two things had to be worked around to get PMOIRED to read OITOOLS output; both
are worth knowing about independently of this study.

1. **`TDIM` on scalar columns.** OIFITS.jl (through FITSIO/cfitsio) writes a
   `TDIMn = (1)` card for every column, including scalar ones. astropy honours
   it and returns `TARGET_ID` with shape (nrow, 1) instead of (nrow,), after
   which PMOIRED dies in `loadOI` with
   `TypeError: cannot use 'numpy.ndarray' as a set element`. Every file written
   by `simulate` or `simulate_from_oifits` is affected. `sanitize.py` strips the
   redundant cards; `OITOOLS.oifits_fix_tdim(file)` does the same from Julia and
   is the better fix for files you intend to share.

2. **Module-level imports in PMOIRED.** `pmoired/oimodels.py` imports `serial`
   (only its optional Arduino progress-bar helper uses it) and `oifake.py`
   imports `astroquery` (only for online catalogue queries). `stubs/` provides
   inert stand-ins so the harness runs without installing either.
