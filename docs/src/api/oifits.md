# OIFITS Handling

## Reading

| Function | Description |
|----------|-------------|
| `readoifits(file)` | Read an OIFITS file into an `OIdata` struct; `binning` averages channels, `insname` / `targetname` read one instrument / target |
| `readoifits_multiepochs(file)` | Read multi-epoch OIFITS data |
| `list_oifits_targets(file)` | List all target names in an OIFITS file |
| `filter_data(data, idx; drop)` | Remove the points `set_data_filter` selects; `drop` removes whole observables (`:t3amp`, `:visphi`, …) |
| `set_data_filter(data; kwargs...)` | Set persistent data filters |
| `readfits(file)` | Read a FITS image into a matrix |
| `writefits(data, file; pixsize, wavelengths, bandwidths, mjds, header)` | Write an image or a cube to FITS, with its WCS, in the convention shared with ROTIR |
| `read_model_image(file)` | Read a FITS image or cube for use as a model component: East-first, with its pixel scale, wavelengths and epochs |
| `model_image_plane(image, λ, mjd)` | Which plane of a model image a data point uses |
| `model_image_vis(image, uv; wl, mjd)` | Its complex visibilities at uv points, plane by plane |
| `fits_pixsize(file)` | Pixel scale of a FITS image in mas, from its WCS keywords, or `nothing` |
| `oifits_prep(data; kwargs...)` | Inflate error bars (additive/relative floors, multiplicative scaling) |
| `write_oifits(path, data; overwrite)` | Write a dataset, filtered or binned, to an OIFITS v2 file that reads back to the same points |
| `merge_oidata(datas)`, `readoifits(files::Vector)` | Several files as one dataset: stations matched by name, every point kept |
| `updatefits_aspro(in, out, pixsize)` | Add ASPRO-compatible WCS headers to a FITS image |
| `OITOOLS.load_oifits(file; notes)` | The `OIDataSet` every reader starts from, read table by table (below) |

```@docs
OIdata
writefits
fits_pixsize
ModelImage
read_model_image
model_image_plane
model_image_vis
readoifits
readoifits_multiepochs
filter_data
set_data_filter
oifits_prep
write_oifits
merge_oidata
OITOOLS.load_oifits
```

### Files that break the standard

Every reader opens a file through `load_oifits`, which reads it table by table. A table that
does not read as the standard defines it is read the closest way that keeps its data, and the
file is otherwise read as declared:

| In the file | Read as |
|---|---|
| a table without a column its declared revision requires (`FOV` in OI_ARRAY) | the revision whose columns it has |
| a revision OIFITS does not define (OI_ARRAY at 5) | the highest defined revision that reads |
| OI_FLUX with its data in `FLUX` (GRAVITY) | `FLUX` as `FLUXDATA` |
| OI_VIS with amplitudes only, or phases only | the missing half "not measured" (see below) |
| rows whose TARGET_ID is not in OI_TARGET | left out; a table left without rows is skipped |
| a ZIP archive holding one FITS file | that FITS file |

`readoifits` lists what it worked around in one warning per file (`warn = true`), and
`validate_oifits` reports each item.

**A missing half.** A T3 amplitude and its closure phase are one point, as are a VIS amplitude
and its phase. When a point carries one half only, the other is held as `0 ± Inf`: it is not
counted (`nt3amp`, `nvisphi`, …), contributes nothing to χ², is not drawn, and is written back
as NaN. `filter_data(…; drop = [:t3amp])` makes the same of a half the file does carry.

## Inspecting, checking and splitting files

| Function | Description |
|----------|-------------|
| `oifits_headers(file)` | Every HDU of a FITS file: extension, instrument, revision, rows, columns, image dimensions, and every header card in order |
| `validate_oifits(file)` | A list of `(severity, table, message)`: what does not follow the standard, what was read around, and suspect values |
| `repair_oifits(file, out; targetname)` | A copy of the file as clean OIFITS v2, one target at a time |
| `split_by_time(data; gap, binwidth, edges)` | One `OIdata` per epoch, the epochs taken from every observable's times |
| `split_by_instrument(data)` | One `OIdata` per instrument, from the table each point was read from |
| `split_by_target(file; kwargs...)` | One dataset per target, each read with `targetname` |
| `fits_set_keywords(file, out; pixsize, east_left, axis3, wavelengths, mjds, refpix)` | A FITS image or cube written with the keywords `read_model_image` reads; the file's reference pixel and other cards kept |
| `fits_refpix(file)` | Its reference pixel `(CRPIX1, CRPIX2)`, or `nothing` |
| `fits_keyword_cube(file; east_left, axis3)` | Its array as an `(East, North, wavelength, epoch)` cube |

```@docs
oifits_headers
validate_oifits
repair_oifits
split_by_time
split_by_instrument
split_by_target
fits_set_keywords
fits_keyword_cube
fits_refpix
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
