# OIFITS Handling

## Reading

| Function | Description |
|----------|-------------|
| `readoifits(file)` | Read an OIFITS file into an `OIdata` struct; `binning` averages channels, `insname` / `targetname` read one instrument / target |
| `readoifits_multiepochs(file)` | Read multi-epoch OIFITS data |
| `list_oifits_targets(file)` | List all target names in an OIFITS file |
| `filter_data(data, idx; drop)` | Remove the points `set_data_filter` selects; `drop` removes whole observables (`:t3amp`, `:visphi`, …) |
| `set_data_filter(data; kwargs...)` | The indices of the points outside the given ranges (wavelength, MJD, baseline, …), for `filter_data` |
| `readfits(file)` | Read a FITS image into a matrix |
| `writefits(data, file; pixsize, wavelengths, bandwidths, mjds, header)` | Write an image or a cube to FITS, with its WCS, in the convention shared with ROTIR |
| `read_model_image(file)` | Read a FITS image or cube for use as a model component: East-first, with its pixel scale, wavelengths and epochs |
| `model_image_plane(image, λ, mjd)` | Which plane of a model image a data point uses |
| `model_image_vis(image, uv; wl, mjd)` | Its complex visibilities at uv points, plane by plane |
| `fits_pixsize(file)` | Pixel scale of a FITS image in mas, from its WCS keywords, or `nothing` |
| `oifits_prep(data; kwargs...)` | Inflate error bars (additive/relative floors, multiplicative scaling) |
| `write_oifits(path, data; targetname, overwrite)` | Write a dataset, filtered or binned, to an OIFITS v2 file |
| `merge_oidata(datas)`, `readoifits(files::Vector)` | Several files as one dataset: stations matched by name, every point kept |
| `updatefits_aspro(in, out, pixsize)` | Add ASPRO-compatible WCS headers to a FITS image |
| `OITOOLS.load_oifits(file; notes)` | The file as an `OIDataSet`, read table by table (below) |

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

All readers open a file through `load_oifits`, which reads it table by table. These deviations
from the standard are read as follows:

| In the file | Read as |
|---|---|
| a table without a column its declared revision requires (`FOV` in OI_ARRAY) | the revision whose columns it has |
| a revision OIFITS does not define (OI_ARRAY at 5) | the highest defined revision that reads |
| OI_FLUX with its data in `FLUX` (GRAVITY) | `FLUX` as `FLUXDATA` |
| OI_VIS with amplitudes only, or phases only | the missing half "not measured" (see below) |
| rows whose TARGET_ID is not in OI_TARGET | left out; a table left without rows is skipped |
| a ZIP archive holding one FITS file | that FITS file |

`readoifits` lists these in one warning per file (`warn = true`); `validate_oifits` reports each.

**Missing half.** A T3 amplitude and its closure phase are one point, as are a VIS amplitude and
its phase. When a point has only one half, the other is stored as `0 ± Inf`: it is not counted
(`nt3amp`, `nvisphi`, …), not in χ², not plotted, and written back as NaN.
`filter_data(…; drop = [:t3amp])` does the same to a half the file has.

## Inspecting, checking and splitting files

| Function | Description |
|----------|-------------|
| `oifits_headers(file)` | Every HDU of a FITS file: extension, instrument, revision, rows, columns, image dimensions, and every header card in order |
| `oifits_table(file, hdu; per_channel)` | A table HDU's contents as the file stores them; an OI data table per channel, with λ, stations, B and B/λ joined |
| `validate_oifits(file)` | A list of `(severity, table, message)`: what does not follow the standard, what was read around, and suspect values |
| `repair_oifits(file, out; targetname)` | A copy of the file as clean OIFITS v2, one target at a time |
| `split_by_time(data; gap, binwidth, edges)` | One `OIdata` per epoch, the epochs taken from every observable's times |
| `split_by_instrument(data)` | One `OIdata` per instrument, from the table each point was read from |
| `split_by_target(file; kwargs...)` | One dataset per target, each read with `targetname` |
| `join_targets(files, out)` | Several OIFITS files written as one, a target per file listed in one OI_TARGET |
| `fits_set_keywords(file, out; pixsize, east_left, axis3, wavelengths, mjds, refpix)` | A FITS image or cube written with the keywords `read_model_image` reads; the file's reference pixel and other cards kept |
| `fits_refpix(file)` | Its reference pixel `(CRPIX1, CRPIX2)`, or `nothing` |
| `fits_keyword_cube(file; east_left, axis3)` | Its array as an `(East, North, wavelength, epoch)` cube |

```@docs
oifits_headers
oifits_table
validate_oifits
repair_oifits
split_by_time
split_by_instrument
split_by_target
join_targets
fits_set_keywords
fits_keyword_cube
fits_refpix
```

## File utilities

| Function | Description |
|----------|-------------|
| `oifits_fix_tdim(file)` | Remove the `TDIM` cards cfitsio writes on scalar columns, in place (astropy and PMOIRED read those columns as 2-D) |
| `oifits_check(; infile)` | Run oifitslib's `oifits-check` on a file |
| `oifits_merge(; outfile, infiles)` | Run oifitslib's `oifits-merge`; `infiles` a vector or a glob pattern |
| `oifits_filter(; outfile, infile, …)` | Run oifitslib's `oifits-filter` with its selection options |

The last three call the command-line tools of John Young's
[oifitslib](https://github.com/jsy1001/oifitslib), which must be on the `PATH`.

```@docs
oifits_fix_tdim
```

## PMOIRED compatibility

| Function | Description |
|----------|-------------|
| `pmoired_to_dict(s)` | Convert a PMOIRED model string to a Julia `Dict` |
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
| `simbad_target(name)` | SIMBAD record in one request: coordinates, proper motion, parallax, radial velocity, spectral type, object type, magnitudes |
| `ra_dec_from_simbad(name)` | RA and Dec in decimal degrees |
| `magnitudes_from_simbad(name)` | Magnitudes, B through N |
| `simbad_tap(adql)` | Any ADQL query against SIMBAD's TAP service; returns column names and rows |
| `SIMBAD_BANDS` | The ten bands `simbad_target` reports, in panel order |
| `SIMBAD_TAP_URL` | The TAP endpoint queries go to |


```julia
t = simbad_target("Vega")
t.main_id           # "* alf Lyr"
t.ra, t.dec         # 279.2347, 38.7837 (degrees)
t.plx               # 130.23 mas
t.sptype            # "A0V"
t.mags["K"]         # 0.129
t.mags["N"]         # NaN: no N magnitude in SIMBAD
```

Coordinates are in decimal degrees. A value SIMBAD does not hold is `NaN` (`""` for a string).
`mags` has all ten `SIMBAD_BANDS` as keys. `ra_dec_from_simbad` and `magnitudes_from_simbad`
each make the same single request as `simbad_target`. Requests go through `Downloads` (no
Python). An unreachable service, a refused query and a query with no rows raise errors with different
messages.

`simbad_tap` runs any query on the
[SIMBAD schema](https://simbad.cds.unistra.fr/simbad/tap/tapsearch.html) and returns the header
and the rows as `String`s:

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
