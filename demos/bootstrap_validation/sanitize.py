#!/usr/bin/env python3
"""Make OIFITS files written by OITOOLS readable by PMOIRED (and astropy).

OIFITS.jl (through FITSIO/cfitsio) emits a TDIMn card for *every* column,
including scalar ones:

    TTYPE1 = TARGET_ID / TFORM1 = 1I / TDIM1 = (1)

astropy honours TDIM and returns a column of shape (nrow, 1) instead of
(nrow,).  PMOIRED then does `set(hdu.data['TARGET_ID'])` and dies with

    TypeError: cannot use 'numpy.ndarray' as a set element

TDIM is redundant for one-dimensional columns, so we simply drop those cards.
Files are rewritten in place; anything already clean is left untouched.

    python3 sanitize.py <file-or-directory> [...]
"""
import os
import sys

from astropy.io import fits


def sanitize(path):
    changed = False
    with fits.open(path, mode="update") as hdul:
        for hdu in hdul:
            hdr = getattr(hdu, "header", None)
            if hdr is None or "TFIELDS" not in hdr:
                continue
            for i in range(1, hdr["TFIELDS"] + 1):
                key = "TDIM%d" % i
                if key in hdr and hdr[key].count(",") == 0:
                    del hdr[key]
                    changed = True
        if changed:
            hdul.flush()
    return changed


def walk(target):
    if os.path.isdir(target):
        for root, _, files in os.walk(target):
            for f in sorted(files):
                if f.endswith((".oifits", ".fits")):
                    yield os.path.join(root, f)
    else:
        yield target


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    n = k = 0
    for target in sys.argv[1:]:
        for path in walk(target):
            n += 1
            k += bool(sanitize(path))
    print("sanitize: %d file(s) scanned, %d modified" % (n, k))
