#!/usr/bin/env python3
"""Fit every simulated universe with PMOIRED and record the 1σ it reports.

    python3 run_pmoired.py <outdir> [nboot] [regime ...]

Output: <outdir>/results_pmoired.tsv (long format)

    regime  universe  method  param  best  sigma  chi2r

Methods
    pmoired_fit         analytic covariance from dpfit.leastsqFit, rescaled to
                        chi2r = 1 (normalizedUncer=True, PMOIRED's default)
    pmoired_boot        oi.bootstrapFit(...) with PMOIRED's defaults
                        (two half-samples per replicate, starting point
                        randomised by 1 sigma)
    pmoired_boot_nostart  same, with randomiseParam=False, which isolates the
                        resampling scatter from the starting-point scatter

The model is identical to the one OITOOLS fits: the primary flux is pinned by
'A,f' = '1-$B,f' so that the total flux is 1, which is the only way the two
packages evaluate the same function (PMOIRED normalises the visibility by the
total flux, OITOOLS does not).
"""
import os
import sys
import time

# PMOIRED imports `serial` and `astroquery` at module level; neither is needed
# here.  demos/bootstrap_validation/stubs provides inert stand-ins.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "stubs"))
PMOIRED_DIR = os.environ.get("PMOIRED_DIR", os.path.expanduser("~/SOFTWARE/PMOIRED"))
sys.path.insert(0, PMOIRED_DIR)

import numpy as np  # noqa: E402
import pmoired  # noqa: E402

# leave cores free when the OITOOLS runner is going at the same time
if "PMOIRED_MAX_THREADS" in os.environ:
    pmoired.MAX_THREADS = int(os.environ["PMOIRED_MAX_THREADS"])

REGIMES = ["ideal", "systematics", "underestimated", "fewblocks", "mismatch"]
FREE = ["A,ud", "B,f", "B,x", "B,y"]
MODEL = {
    "A,ud": 0.80,
    "A,f": "1-$B,f",
    "B,ud": 0.0,
    "B,f": 0.05,
    "B,x": 1.20,
    "B,y": -0.80,
}


def fit_one(path, nboot):
    """Return {method: (best, sigma)} plus the reduced chi2 of the single fit."""
    oi = pmoired.OI(path, verbose=False)
    oi.setupFit({"obs": ["V2", "T3PHI"]})
    oi.doFit(MODEL, doNotFit=["B,ud"], verbose=0)
    best = [float(oi.bestfit["best"][p]) for p in FREE]
    chi2r = float(oi.bestfit["chi2"])
    out = {"pmoired_fit": (best, [float(oi.bestfit["uncer"][p]) for p in FREE])}

    for method, kw in (("pmoired_boot", {}),
                       ("pmoired_boot_nostart", {"randomiseParam": False})):
        oi.bootstrapFit(nboot, verbose=0, **kw)
        out[method] = (best, [float(oi.boot["uncer"][p]) for p in FREE])
    return out, chi2r


def main():
    outdir = sys.argv[1]
    nboot = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    which = sys.argv[3:] if len(sys.argv) > 3 else REGIMES

    rows = []
    t0 = time.time()
    for regime in which:
        d = os.path.join(outdir, regime)
        if not os.path.isdir(d):
            continue
        files = sorted(f for f in os.listdir(d) if f.endswith(".oifits"))
        print("%-15s : %d universes" % (regime, len(files)), flush=True)
        for u, f in enumerate(files, start=1):
            try:
                res, chi2r = fit_one(os.path.join(d, f), nboot)
            except Exception as exc:                      # noqa: BLE001
                print("  universe %d failed: %r" % (u, exc), flush=True)
                continue
            for method, (best, sigma) in res.items():
                for p, b, s in zip(FREE, best, sigma):
                    rows.append((regime, u, method, p, b, s, chi2r))
            if u % 10 == 0:
                print("  %4d/%d  (%.1f s elapsed)" % (u, len(files), time.time() - t0),
                      flush=True)

    out = os.path.join(outdir, "results_pmoired.tsv")
    with open(out, "w") as fh:
        fh.write("regime\tuniverse\tmethod\tparam\tbest\tsigma\tchi2r\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print("\nWrote %d rows to %s  (%.1f s)" % (len(rows), out, time.time() - t0))


if __name__ == "__main__":
    main()
