# SQUEEZE port — recovery bundle

The working tree was overwritten by an rsync from voyager, which destroyed the
in-place integration. These files are the port, staged out of the way so the
tree stays clean. **They are not wired in.** See "Integration" below — this
voyager version of OITOOLS has diverged and a straight copy will not work.

## What is here

| file | install to | lines | state |
|---|---|---|---|
| `squeeze.jl` | `src/squeeze.jl` | 1771 | complete, final version (all later fixes folded in) |
| `OITOOLSPigeonsExt.jl` | `ext/OITOOLSPigeonsExt.jl` | 321 | complete |
| `cvis_to_chi2_noalloc.jl` | append to `src/chi2_flat.jl` | 134 | complete |
| `test_squeeze.jl` | `test/test_squeeze.jl` | 452 | complete |
| `test_squeeze_tempering.jl` | `test/test_squeeze_tempering.jl` | 217 | complete, skips without Pigeons |
| `api_imaging_mcmc.md` | `docs/src/api/imaging_mcmc.md` | 60 | replaces the "*Coming soon.*" stub |
| `examples_squeeze.md` | `docs/src/examples/squeeze.md` | 427 | complete |
| `demo_squeeze.jl` | `demos/example_image_reconstruction_squeeze.jl` | | complete |
| `demo_squeeze_sparco_v1295.jl` | `demos/..._squeeze_sparco_v1295.jl` | | complete |
| `demo_squeeze_tempered.jl` | `demos/..._squeeze_tempered.jl` | | complete |

Everything is recovered. Nothing is wired in.

## Integration (do NOT just copy — this tree has diverged)

1. **Append `cvis_to_chi2_noalloc.jl` to `src/chi2_flat.jl`** and export
   `cvis_to_chi2_noalloc`.
   This version of OITOOLS already has `cvis_to_chi2_f` (a forward-only chi2 that
   also handles the flux term) — but it routes through `_chi2_terms`, which STILL
   allocates: `abs2.(V[data.indx_v2])`, `V[data.indx_t3_1]`, `abs.(t3)`,
   `_mod360(...)` are all gathers or broadcast temporaries, ~10 heap allocations
   per call, several of them nuv-length. The sampler calls this once per proposal,
   millions of times, so that dominates. Hence a separate allocation-free copy.
   Verified free of name clashes against this tree: `_mod360s`, `_RAD2DEG`,
   `_DEG2RAD`, `_weights_tuple`, `_cvis_to_chi2`, `cvis_to_chi2_noalloc` are all
   unused here. `squeeze.jl` calls `_cvis_to_chi2` / `_weights_tuple`, so it drops
   in with no edits.
   (Note `_pad_weights` in this tree is now **two-argument**: `_pad_weights(w, T)`.
   `_weights_tuple` is separate and does not collide.)

2. `src/OITOOLS.jl`: add `include("squeeze.jl")` after `include("bsdmm.jl")`, and
   `export reconstruct_squeeze, reconstruct_squeeze_tempered, SqueezeSparco, default_nelements`

3. `Project.toml`: add

       [weakdeps]
       Pigeons = "0eb8d820-af6a-4919-95ae-11206f830c31"

       [extensions]
       OITOOLSPigeonsExt = "Pigeons"

   and `Pigeons = "0.4"` under `[compat]`.

4. **This tree uses PythonCall/PythonPlot, not PyCall/PyPlot** (commit 56753f8,
   "Big jump PyPlot->PythonPlot and PyCall->PythonCall"). The live-monitor code in
   `squeeze.jl` calls `PyPlot.isinteractive()` and `PyPlot.pause(0.001)` — those
   two lines must become `PythonPlot.`-qualified (or be dropped). Everything else
   the monitor uses (`figure`, `clf`, `subplot`, `plot`, `ylabel`, `yscale`,
   `grid`, `subplots_adjust`, `xlabel`, `imdisp`) has the same name under
   PythonPlot. `L"..."` still comes from LaTeXStrings.

5. **The environment does not currently load.** `Pkg.resolve()` succeeds, but
   `using OITOOLS` dies in `fit_model.jl:20` initialising PythonCall — CondaPkg is
   trying to change `python`/`numpy`/`matplotlib`/`scipy` and add `corner`, and
   failing. This is rsync fallout, not related to the port. Nothing can be tested
   until it is fixed.

## Two one-line doc edits still to make by hand

- `docs/make.jl`: add `"SQUEEZE (MCMC)" => "examples/squeeze.md",` to the Guides
  list, after the BSDMM entry.
- `docs/src/api/imaging_hybrid.md`: add `cvis_to_chi2_noalloc` to the function
  table and the `@docs` block, beside `cvis_to_chi2_fg`.

Note `docs/make.jl` uses `checkdocs = :exports`, so **every exported symbol must
appear in a `@docs` block or the docs build fails**. The exports to add are
`reconstruct_squeeze`, `reconstruct_squeeze_tempered`, `SqueezeSparco`,
`default_nelements` (all covered by `api_imaging_mcmc.md`) and
`cvis_to_chi2_noalloc` (covered by the `imaging_hybrid.md` edit above).

## Verified results before the loss (2004-data1.oifits, 64px @ 0.2 mas)

- annealed chi2r 0.75 (4 chains x 400 sweeps); flat-image chi2r 3983
- tempered chi2r 0.771, logZ ~ -1187, global barrier Lambda ~ 27 (30 chains x 2^11)
- rank-1 incremental vs full recompute: 1.3e-14 relative after 200k moves (Float64)
- cvis_to_chi2: 0 alloc, ~9x faster than cvis_to_chi2_fg
- SPARCO on real V1295 Aql: chi2r 68.9 (image only) -> 0.94; free fit f_star = 0.446 vs published 0.48
- SPARCO ud recovery, true 2.0 mas: annealed from 4.0 sticks at 3.73; tempered reaches 2.15
- full suite 655/655 without Pigeons, 40/40 tempering tests with it

## Known caveats carried in the code

- `f_env` uses the wavelength-scaled `f_bg`, not `f_bg_ref`, matching
  `modelcode_sparco.c`. Owner to confirm which is intended.
- Annealing's parametric proposal cannot travel for shape parameters (the image
  adapts to the current value; a 0.05 mas ud nudge costs dlogL ~ 2.8e4). Use
  tempering, or scan the parameter on a grid.
- Only the final Pigeons round is kept; earlier rounds are tuning.
