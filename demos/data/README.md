# Demo data

Sample OIFITS files for the demos, the documentation examples and the test suite.

## Files the test suite requires

`test/` reads its fixtures from here rather than keeping a second copy — the two directories
held 42 MB of byte-identical files, `OBJECT2_K.oifits` alone twice at 33 MB. That makes the
following load-bearing: **deleting one breaks `Pkg.test()`, not just a demo.**

| file | used by |
|---|---|
| `BC2004/2004-data1.oifits` | bsmem precision + regression, oifits TDIM, SQUEEZE, tempering, Makie ext, `bin/trace.jl` |
| `AlphaCenA.oifits` | nested sampling, bootstrap, constraints |
| `2019_v1295Aql.WL_SMOOTH.A.oifits` | SQUEEZE polychromatic, tempering |
| `polaris.oifits` | GUI drawing-path tests |
| `BC2026/OBJECT1_N.oifits` | bsmem polychromatic cases |

`test/gui/runtests.jl` checks these exist before running and names any that do not, so a
missing fixture reports itself rather than surfacing as a confusing read error.

Everything else here is demo material and can be pruned freely.
