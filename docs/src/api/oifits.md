# OIFITS Handling

## Reading

| Function | Description |
|----------|-------------|
| `readoifits(file)` | Read an OIFITS file into an `OIdata` struct |
| `readoifits_multiepochs(file)` | Read multi-epoch OIFITS data |
| `list_oifits_targets(file)` | List all target names in an OIFITS file |
| `filter_data(data; kwargs...)` | Filter data by baseline, wavelength, etc. |
| `set_data_filter(data; kwargs...)` | Set persistent data filters |
| `readfits(file)` | Read a FITS image into a matrix |
| `writefits(data, file)` | Write a matrix to a FITS image |
| `oifits_prep(data; kwargs...)` | Inflate error bars (additive/relative floors, multiplicative scaling) |
| `updatefits_aspro(in, out, pixsize)` | Add ASPRO-compatible WCS headers to a FITS image |

```@docs
OIdata
readoifits
readoifits_multiepochs
filter_data
set_data_filter
```

## Writing

| Function | Description |
|----------|-------------|
| `oifits_check(file)` | Validate an OIFITS file |
| `oifits_merge(files, outfile)` | Merge multiple OIFITS files |
| `oifits_filter(infile, outfile)` | Filter an OIFITS file |
| `oifits_fix_tdim(file)` | Drop the redundant `TDIM` cards cfitsio writes on scalar columns (needed for astropy/PMOIRED to read the file) |

```@docs
oifits_fix_tdim
```

## PMOIRED compatibility

| Function | Description |
|----------|-------------|
| `pmoired_to_dict(s)` | Convert a PMOIRED model string directly to a Julia `Dict` |
| `pmoired_to_julia(s)` | Convert a PMOIRED model string to Julia source code (returns `String`) |
| `pmoired_to_julia_file(infile, outfile)` | Convert a PMOIRED notebook snippet file |
| `dict_to_pmoired(d)` | Render an OITOOLS model dict as PMOIRED source (returns `String`) |
| `dict_to_pmoired_file(d, outfile)` | Write an OITOOLS model dict out as PMOIRED source |

```@docs
pmoired_to_dict
pmoired_to_julia
pmoired_to_julia_file
dict_to_pmoired
dict_to_pmoired_file
```

## Utilities

| Function | Description |
|----------|-------------|
| `recenter(x; mask, max)` | Recenter an image by circular shift to centroid or peak |
| `sexagesimal_to_degrees(s)` | Parse `"-46 28 00.5"` or `"12:34:56"` to decimal degrees (sign applies to the whole value) |

```@docs
recenter
sexagesimal_to_degrees
```

## SIMBAD

| Function | Description |
|----------|-------------|
| `simbad_target(name)` | Full SIMBAD record in one request: coordinates, proper motion, parallax, radial velocity, spectral type, object type, magnitudes |
| `ra_dec_from_simbad(name)` | RA and Dec alone, in **decimal degrees** |
| `magnitudes_from_simbad(name)` | Photometric magnitudes alone, B through N |
| `simbad_tap(adql)` | Any ADQL query against SIMBAD's TAP service; returns column names and rows |
| `SIMBAD_BANDS` | The ten bands `simbad_target` reports, in panel order |
| `SIMBAD_TAP_URL` | The TAP endpoint queries go to |

`simbad_target` is the one to reach for when planning. Proper motion places the target at the
epoch actually being observed, parallax turns an angular diameter into a physical one, and the
spectral type is what a surface-brightness relation needs to predict a diameter before anything
is measured:

```julia
t = simbad_target("Vega")
t.main_id           # "* alf Lyr"
t.ra, t.dec         # 279.2347, 38.7837   — degrees
t.plx               # 130.23 mas
t.sptype            # "A0V"
t.mags["K"]         # 0.129
t.mags["N"]         # NaN — SIMBAD has no N magnitude for this target
```

Coordinates come back in **decimal degrees**, and a value SIMBAD does not hold is `NaN` (or
`""` for a string) rather than an error: not knowing a parallax is a fact about the target, and
has to be distinguishable from the query having failed. `mags` always carries all ten
`SIMBAD_BANDS` as keys, so a panel iterating them never has to test for presence as well as for
`NaN`.

It is also the cheapest call. Astrometry and photometry arrive together from a single ADQL
query, which is what SIMBAD asks clients to do; `ra_dec_from_simbad` and
`magnitudes_from_simbad` are wrappers around it, so asking for both costs two requests where
`simbad_target` costs one. The transport is `Downloads`, a standard library — no Python is
involved.

Failures are reported apart from one another, because they call for different actions: the
service being unreachable, the service refusing the query, and the query succeeding with no
rows. Only the last is a problem with the name you asked for.

`simbad_tap` is the escape hatch for anything else in the
[SIMBAD schema](https://simbad.cds.unistra.fr/simbad/tap/tapsearch.html). It returns the header
and the rows as `String`s, exactly as the service sent them:

```julia
cols, rows = simbad_tap("SELECT TOP 5 main_id, plx_value FROM basic WHERE plx_value > 500")
```

```@docs
simbad_target
ra_dec_from_simbad
magnitudes_from_simbad
simbad_tap
SIMBAD_BANDS
SIMBAD_TAP_URL
```
