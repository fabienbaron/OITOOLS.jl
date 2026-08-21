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
| `query_target_from_simbad(name)` | Query SIMBAD for target information |
| `ra_dec_from_simbad(name)` | Resolve a target through SIMBAD; returns RA and Dec in **decimal degrees** |
| `sexagesimal_to_degrees(s)` | Parse `"-46 28 00.5"` or `"12:34:56"` to decimal degrees (sign applies to the whole value) |
| `magnitudes_from_simbad(name)` | Query SIMBAD for photometric magnitudes (B through N) |
| `simbad_target(name)` | Full SIMBAD record: coordinates, proper motion, parallax, spectral type, magnitudes |
| `SIMBAD_BANDS` | Photometric bands `simbad_target` asks for, in panel order |

`simbad_target` is the one to reach for when planning: proper motion places the target at the
epoch actually being observed, parallax turns an angular diameter into a physical one, and the
spectral type is what a surface-brightness relation needs to predict a diameter before anything
is measured. Coordinates come back in **decimal degrees**, and a value SIMBAD does not hold is
`NaN` rather than an error — not knowing a parallax is a fact about the target, and has to be
distinguishable from the query having failed.

```@docs
recenter
magnitudes_from_simbad
ra_dec_from_simbad
sexagesimal_to_degrees
simbad_target
SIMBAD_BANDS
```
