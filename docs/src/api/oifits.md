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

```@docs
pmoired_to_dict
pmoired_to_julia
pmoired_to_julia_file
```

## Utilities

| Function | Description |
|----------|-------------|
| `recenter(x; mask, max)` | Recenter an image by circular shift to centroid or peak |
| `query_target_from_simbad(name)` | Query SIMBAD for target information |
| `ra_dec_from_simbad(name)` | Get RA/Dec from SIMBAD |
| `magnitudes_from_simbad(name)` | Query SIMBAD for photometric magnitudes (V, J, H, K, L, M, N) |

```@docs
recenter
magnitudes_from_simbad
```
