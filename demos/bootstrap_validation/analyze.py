#!/usr/bin/env python3
"""Summarise the validation study.

    python3 analyze.py <outdir> [--csv]

Reads <outdir>/results_oitools.tsv and <outdir>/results_pmoired.tsv and, for
every (regime, method, parameter), compares what the method *claims* with what
the ensemble of universes actually *does*:

    sigma_true    std over universes of the best-fit value
                  (computed per package, from that package's own fits)
    sigma_report  mean over universes of the 1 sigma the method quotes
    ratio         sigma_report / sigma_true   -> 1.00 = calibrated
    coverage      fraction of universes with |best - truth| < sigma_report
                  -> 0.683 for a calibrated 1 sigma Gaussian
    bias          (mean best - truth) / sigma_true

A method can only be blamed for `ratio`.  `coverage` also folds in the bias of
the estimator itself, which is why the `mismatch` regime destroys it for every
method at once.
"""
import os
import sys
from collections import defaultdict

import numpy as np

TRUTH = {"A,ud": 0.80, "B,f": 0.05, "B,x": 1.20, "B,y": -0.80}
PARAMS = ["A,ud", "B,f", "B,x", "B,y"]
REGIMES = ["ideal", "systematics", "underestimated", "fewblocks", "mismatch"]
METHOD_ORDER = [
    "lsqfit", "montecarlo", "boot_config", "boot_epoch", "boot_weights",
    "boot_halfsample", "boot_pmoired",
    "pmoired_fit", "pmoired_boot", "pmoired_boot_nostart",
]
PACKAGE = {m: ("PMOIRED" if m.startswith("pmoired") else "OITOOLS") for m in METHOD_ORDER}
LABEL = {
    "lsqfit":               "OITOOLS  analytic (LM)",
    "montecarlo":           "OITOOLS  parametric MC",
    "boot_config":          "OITOOLS  block bootstrap",
    "boot_epoch":           "OITOOLS  epoch bootstrap",
    "boot_weights":         "OITOOLS  weighted bootstrap",
    "boot_halfsample":      "OITOOLS  half-sample (BRR)",
    "boot_pmoired":         "OITOOLS  PMOIRED scheme",
    "pmoired_fit":          "PMOIRED  analytic (dpfit)",
    "pmoired_boot":         "PMOIRED  bootstrapFit",
    "pmoired_boot_nostart": "PMOIRED  bootstrapFit (no start jitter)",
}


def load(path, rows):
    if not os.path.exists(path):
        print("missing:", path)
        return
    with open(path) as fh:
        next(fh)
        for line in fh:
            regime, universe, method, param, best, sigma, chi2r = line.rstrip("\n").split("\t")
            rows[(regime, method, param)].append(
                (int(universe), float(best), float(sigma), float(chi2r)))


def main():
    outdir = sys.argv[1]
    rows = defaultdict(list)
    load(os.path.join(outdir, "results_oitools.tsv"), rows)
    load(os.path.join(outdir, "results_pmoired.tsv"), rows)
    if not rows:
        sys.exit("no results found in " + outdir)

    # per-package truth scatter: the estimator is the same for every method of a
    # package, so sigma_true is a property of (regime, package, param)
    sigma_true = {}
    nuni = {}
    chi2 = {}
    for (regime, method, param), v in rows.items():
        key = (regime, PACKAGE.get(method, "?"), param)
        best = np.array([x[1] for x in v])
        sigma_true.setdefault(key, np.std(best, ddof=1))
        nuni.setdefault(key, len(best))
        chi2.setdefault((regime, PACKAGE.get(method, "?")), np.median([x[3] for x in v]))

    csv = "--csv" in sys.argv
    md = "--md" in sys.argv
    if csv:
        print("regime,method,param,n,sigma_true,sigma_report,ratio,coverage,bias")

    if md:
        # compact summary: ratio per parameter, plus mean coverage, per regime
        for regime in REGIMES:
            present = [m for m in METHOD_ORDER if (regime, m, PARAMS[0]) in rows]
            if not present:
                continue
            print("\n### %s\n" % regime)
            print("| method | " + " | ".join(PARAMS) + " | coverage |")
            print("|" + "---|" * (len(PARAMS) + 2))
            for method in present:
                cells, covs = [], []
                for param in PARAMS:
                    v = rows[(regime, method, param)]
                    best = np.array([x[1] for x in v])
                    sig = np.array([x[2] for x in v])
                    st = sigma_true[(regime, PACKAGE[method], param)]
                    good = np.isfinite(sig) & (sig > 0)
                    cells.append("%.2f" % (np.mean(sig[good]) / st))
                    covs.append(np.mean(np.abs(best - TRUTH[param]) < sig))
                print("| %s | %s | %.2f |"
                      % (LABEL.get(method, method), " | ".join(cells), np.mean(covs)))
        return

    for regime in REGIMES:
        present = [m for m in METHOD_ORDER if (regime, m, PARAMS[0]) in rows]
        if not present:
            continue
        if not csv:
            n = nuni.get((regime, "OITOOLS", PARAMS[0]), 0)
            m = nuni.get((regime, "PMOIRED", PARAMS[0]), 0)
            print("\n" + "=" * 100)
            print("REGIME: %s   (%d OITOOLS universes, %d PMOIRED universes,"
                  "  median chi2r: OITOOLS %.2f / PMOIRED %.2f)"
                  % (regime, n, m,
                     chi2.get((regime, "OITOOLS"), float("nan")),
                     chi2.get((regime, "PMOIRED"), float("nan"))))
            print("=" * 100)
            print("%-42s %-8s %10s %10s %7s %9s %7s"
                  % ("method", "param", "sig_true", "sig_report", "ratio", "coverage", "bias"))
            print("-" * 100)

        for method in present:
            for param in PARAMS:
                v = rows[(regime, method, param)]
                if not v:
                    continue
                best = np.array([x[1] for x in v])
                sig = np.array([x[2] for x in v])
                st = sigma_true[(regime, PACKAGE[method], param)]
                good = np.isfinite(sig) & (sig > 0)
                sr = np.mean(sig[good]) if good.any() else float("nan")
                cov = np.mean(np.abs(best - TRUTH[param]) < sig)
                bias = (np.mean(best) - TRUTH[param]) / st if st > 0 else float("nan")
                if csv:
                    print("%s,%s,%s,%d,%.6g,%.6g,%.3f,%.3f,%.3f"
                          % (regime, method, param, len(v), st, sr, sr / st, cov, bias))
                else:
                    print("%-42s %-8s %10.4g %10.4g %7.2f %9.3f %7.2f"
                          % (LABEL.get(method, method) if param == PARAMS[0] else "",
                             param, st, sr, sr / st, cov, bias))
            if not csv:
                print("-" * 100)

    if not csv:
        print("\nratio    = quoted sigma / true scatter    (1.00 = calibrated,"
              " <1 = overconfident, >1 = conservative)")
        print("coverage = fraction of universes within 1 quoted sigma of the truth"
              "  (0.683 = calibrated)")
        print("bias     = (mean best fit - truth) / true scatter")


if __name__ == "__main__":
    main()
