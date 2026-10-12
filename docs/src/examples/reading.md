# Reading OIFITS files

## Basic usage

`readoifits` loads an OIFITS file and returns a 2D array of `OIdata{T}` indexed as
`[nwavbin, ntimebin]`. The simplest call reads everything into a single bin:

```julia
using OITOOLS
data = readoifits("mystar.oifits")        # returns Array{OIdata{Float32}, 2}
data[1]                                   # the single OIdata bin
```

The observables are `Float32` by default, which halves memory and speeds up fitting and
imaging. Use `T=Float64` for double precision:

```julia
data = readoifits("mystar.oifits"; T=Float64)
```

MJDs (`v2_mjd`, `t3_mjd`, `vis_mjd`, `flux_mjd`, `mean_mjd`) are `Float64` whatever `T` is;
`uv_mjd`, the time a `$MJD` model is evaluated at, is stored in `T` and resolves ~5.6 min at
`Float32`.

## Keyword reference

| Keyword | Default | Description |
|---------|---------|-------------|
| **Target selection** | | |
| `targetname` | `""` | Select a single target by name (all targets if empty) |
| `insname` | `""` | Select one instrument (an OI\_WAVELENGTH table's INSNAME; all if empty) |
| **Spectral / temporal binning** | | |
| `spectralbin` | `[[]]` | Vector of `[λ_min, λ_max]` windows (metres) |
| `temporalbin` | `[[]]` | Vector of `[mjd_min, mjd_max]` windows |
| `polychromatic` | `false` | One spectral bin per channel (overrides `spectralbin`) |
| `merge_oi_wavelength` | `false` | Deduplicate compatible OI\_WAVELENGTH tables (polychromatic mode) |
| `splitting` | `false` | Force multi-bin mode even with default bin vectors |
| `get_specbin_file` | `true` | Auto-derive spectral bins from the file |
| `get_timebin_file` | `true` | Auto-derive temporal bins from the file |
| **Averaging channels** | | |
| `binning` | `0` | Average every N channels into one (an integer), or the channels in given λ windows; `0` keeps every channel |
| `binmode` | `:block` | `:block` (1/σ² weights, circular-mean phases, scatter check) or `:pmoired` (PMOIRED's `binning`) |
| `bincorr` | `1` | Correlation between adjacent channels: averaged errors are multiplied by √`bincorr` |
| `binstats` | `nothing` | A `Dict` that receives, per observable, the χ²/dof of the averaged channels about a line in λ |
| **Observable selection** | | |
| `use_vis` | `true` | Load complex visibilities |
| `use_v2` | `true` | Load squared visibilities |
| `use_t3` | `true` | Load triple products |
| `use_flux` | `true` | Load flux spectra |
| **Quality filtering** | | |
| `filter_bad_data` | `true` | Apply quality cuts (flags, NaN, SNR, range) on load |
| `ignore_flags` | `false` | Keep points the file flags |
| `force_full_vis` | `false` | Require both visamp and visphi to be valid |
| `force_full_t3` | `false` | Require both t3amp and t3phi to be valid |
| `cutoff_minv2` | `-1` | Minimum V² (points below are discarded) |
| `cutoff_maxv2` | `2.0` | Maximum V² |
| `cutoff_mint3amp` | `-1.0` | Minimum T3 amplitude |
| `cutoff_maxt3amp` | `1.5` | Maximum T3 amplitude |
| `filter_v2_snr_threshold` | `0.01` | Minimum \|V²/σ\| to keep |
| `special_filter_diffvis` | `false` | Keep only vis points common across all spectral bins |
| **UV deduplication** | | |
| `redundance_remove` | `true` | Merge UV points closer than `uvtol` |
| `uvtol` | `200.0` | Merge radius in cycles/rad (B/λ) |
| **Output** | | |
| `warn` | `true` | Print warnings about non-standard files |
| `verbose` | `true` | Print summary of loaded tables |
| `T` | `Float32` | Element type of the observables |

## The OIdata struct

Each `OIdata{T}` bin holds flat vectors for every observable:

| Field | Description |
|-------|-------------|
| `v2`, `v2_err`, `v2_baseline`, `v2_lam`, `v2_mjd` | Squared visibilities |
| `t3phi`, `t3phi_err`, `t3amp`, `t3amp_err` | Closure phases and triple amplitudes |
| `visamp`, `visamp_err`, `visphi`, `visphi_err` | Complex visibility amplitude/phase |
| `flux`, `flux_err`, `flux_lam`, `flux_sta_index` | Flux spectra (0 = calibrated) |
| `uv`, `uv_lam`, `uv_baseline` | UV coverage (2×nuv matrix, λ in metres) |
| `indx_v2`, `indx_t3_1/2/3`, `indx_vis` | Indices mapping observables → UV points |
| `sta_name`, `tel_name`, `sta_index` | Station / telescope names |
| `mean_mjd` | Mean MJD of the bin (always `Float64`) |
| `nv2`, `nt3phi`, `nt3amp`, `nvisamp`, `nvisphi`, `nflux`, `nuv` | Counts of the points holding a measurement |
| `tables` | One entry per table read: `(extname, insname, arrname, amptyp, phityp)` |
| `v2_tab`, `t3_tab`, `vis_tab`, `flux_tab` | For each point, the index of its table in `tables` |

Call `display(data)` or `display(data[1])` to print a compact summary.

## Correlated data (OI_CORR)

OITOOLS reads OIFITSv2 `OI_CORR` tables and stores the correlation matrix
and per-point indices in the `OIdata` struct:

| Field | Description |
|-------|-------------|
| `v2_corr`, `v2_corr_idx` | V² correlation matrix and per-point index |
| `t3phi_corr`, `t3phi_corr_idx` | Closure-phase correlation |
| `t3amp_corr`, `t3amp_corr_idx` | Triple-amplitude correlation |
| `visamp_corr`, `visamp_corr_idx` | Visibility amplitude correlation |
| `visphi_corr`, `visphi_corr_idx` | Visibility phase correlation |
| `flux_corr`, `flux_corr_idx` | Flux correlation |

`*_corr_idx[i] == 0` means point `i` has no associated correlation row.

Multi-baseline correlation blocks (e.g. GRAVITY, where `NDATA > nwave`) are
automatically detected and expanded so that each baseline within a timestamp
maps to the correct slice of the correlation matrix.

See `example_chi2_correlated.jl` for a worked example computing whitened
(correlated) chi-squared from GRAVITY data.

## Selecting a target

Files with multiple targets require a `targetname` filter:

```julia
data = readoifits("multi_target.oifits"; targetname="Betelgeuse")
```

Use `list_oifits_targets("multi_target.oifits")` to see all available targets.

## Instruments

Each point records the table it was read from (`tables`, `v2_tab`, …), so a dataset knows the
instrument of every point. `insname` reads one instrument; `split_by_instrument` divides a
dataset already read, merged or filtered:

```julia
sc = readoifits("gravity.fits"; insname = "GRAVITY_SC")
parts = split_by_instrument(readoifits("gravity.fits")[1])   # ["GRAVITY_SC" => …, "GRAVITY_FT" => …]
```

`write_oifits` writes one OI\_WAVELENGTH table per instrument, under its name and with its own
channels, so a written file reads back by instrument as the original does.

## Polychromatic (per-channel) mode

To keep each spectral channel as a separate bin, pass `polychromatic=true`. The
bin boundaries are derived automatically from midpoints between adjacent channel
centres, so each stored λ value falls in exactly one bin:

```julia
data = readoifits("my_spectrum.oifits"; polychromatic=true)
# data is nchannels×1; data[k,1] is the k-th channel bin
```

A file with several OI\_WAVELENGTH tables gets one bin per distinct channel: tables sharing a
grid share their bins. `merge_oi_wavelength=true` checks that the tables' channels overlap and
uses the first table's grid.

To specify bins manually:

```julia
spectralbin = [[1.6e-6, 1.8e-6], [2.0e-6, 2.4e-6]]   # two K-band bins
data = readoifits("my_spectrum.oifits"; spectralbin)
```

Similarly for time bins:

```julia
temporalbin = [[58000.0, 58100.0], [58200.0, 58300.0]]
data = readoifits("mystar.oifits"; temporalbin)
```

## Averaging channels

`binning` averages adjacent channels into one when the file is read; `spectralbin` and
`polychromatic` group channels into bins and keep each one.

```julia
data = readoifits("gravity.fits"; binning = 20)                   # every 20 channels → 1
data = readoifits("gravity.fits"; binning = 20, bincorr = 2)      # adjacent channels correlated
data = readoifits("gravity.fits"; binning = 20, binmode = :pmoired)
stats = Dict{Symbol,Float64}()
readoifits("gravity.fits"; binning = 20, binstats = stats)        # χ²/dof of each average
```

In `:block` mode the average is 1/σ²-weighted, a phase is averaged as a weighted circular mean,
and the error is the error of the weighted mean times √`bincorr`. `binstats` reports the scatter
of the averaged channels about a straight line in λ within each group: near 1 when the errors
describe the scatter, above 1 when they are underestimated, below 1 when adjacent channels are
correlated. OI\_CORR describes the file's own channels and is not used when averaging.

## Quality filtering

`filter_bad_data=true` (the default) removes flagged points and applies sanity
cuts. Key thresholds can be adjusted:

```julia
data = readoifits("mystar.oifits";
    filter_bad_data        = true,
    cutoff_minv2           = -0.05,
    cutoff_maxv2           = 1.05,
    cutoff_mint3amp        = -0.1,
    cutoff_maxt3amp        = 1.5,
    filter_v2_snr_threshold = 1.0)   # require SNR > 1 for V²
```

To selectively disable loading of specific observable types:

```julia
data = readoifits("mystar.oifits"; use_flux=false, use_vis=false)
```

### Verbose output

By default `readoifits` prints a data summary. Pass `verb=false` to suppress
all output, or `verbose=true` to additionally print OI_VIS amplitude/phase
type information (AMPTYP/PHITYP):

```julia
data = readoifits("mystar.oifits"; verbose=true)   # print table summary
data = readoifits("mystar.oifits"; warn=false)      # suppress warnings
```

## UV deduplication

Redundant UV points (within `uvtol` cycles/rad, i.e. B/λ; default 200) are merged:

```julia
data = readoifits("mystar.oifits"; redundance_remove=true, uvtol=100.0)
```

Disable with `redundance_remove=false`.

## Post-load filtering

After loading, use `set_data_filter` + `filter_data` to apply additional cuts
without re-reading the file:

```julia
idx = set_data_filter(data[1];
    filter_bad_data      = true,
    baseline_range       = [10.0e6, 300.0e6],
    wav_range            = [1.6e-6, 1.8e-6])
data_filtered = filter_data(data[1], idx)
```

`drop` removes observables. A T3 amplitude and its closure phase are one point, as are a VIS
amplitude and its phase: dropping one half keeps the point for the other, the dropped half "not
measured" (`0 ± Inf`: not counted, ignored by χ², not drawn, written as NaN).

```julia
no_t3amp = filter_data(data[1]; drop = [:t3amp])
v2_only  = filter_data(data[1]; drop = [:t3amp, :t3phi, :visamp, :visphi, :flux])
```

`split_by_time` divides a dataset into epochs, from the times of every observable:

```julia
parts, mjd = split_by_time(data[1])                    # a new epoch after a 0.5-day gap
parts, mjd = split_by_time(data[1]; binwidth = 1/24)   # hourly windows
```

## Reading multiple files

**Several files as one dataset** — stations matched by name, every point kept:

```julia
data = readoifits(["night1.oifits", "night2.oifits"])
```

**Multi-epoch** — each file is one epoch; returns a `Matrix{OIdata}`:

```julia
files = ["epoch1.oifits", "epoch2.oifits", "epoch3.oifits"]
data = readoifits_multiepochs(files)   # 1×3 Matrix{OIdata}
```

**Multi-colour** — each file is a different waveband; use `readoifits_multiepochs` with one file per band:

```julia
files = ["H_band.oifits", "K_band.oifits"]
data = readoifits_multiepochs(files)   # 1×2 Matrix — one column per file
```

## Error inflation

`oifits_prep` inflates error bars in-place, which is useful for preventing
underestimated uncertainties from dominating the fit:

```julia
oifits_prep(data[1];
    min_v2_err_add   = 0.01,    # additive floor for V² errors
    min_v2_err_rel   = 0.02,    # relative floor (2% of |V²|)
    v2_err_mult      = 1.0,     # multiplicative scaling
    min_t3phi_err_add = 1.0,    # additive floor for T3φ errors (degrees)
    t3phi_err_mult   = 1.0,
    quad             = false)   # true = combine in quadrature
```

By default errors are replaced when the floor exceeds the original; set
`quad=true` to add floors in quadrature instead.

## ASPRO-compatible FITS headers

`updatefits_aspro` adds WCS metadata (pixel scale, reference pixel) to a FITS
image so that ASPRO can read it:

```julia
updatefits_aspro("input.fits", "output.fits", 0.1)   # pixsize in mas
```

## Non-standard files

Files are read table by table (`load_oifits`). A table that does not follow the standard is read
the closest way that keeps its data — a declared revision without the columns it requires, or a
revision that does not exist, is read at the revision that fits; GRAVITY's `FLUX` column is read
as `FLUXDATA`; an OI\_VIS with amplitudes only has its phases "not measured"; rows whose
TARGET\_ID is not in OI\_TARGET are left out — and a ZIP archive holding one FITS file is read
through it. `readoifits` lists what it did in one warning per file:

```julia
data = readoifits("2022-08-03-calibrated_cleaned.fits")   # a ZIP archive
```

`validate_oifits(file)` reports the same, along with references to missing stations or targets
and values with unusable errors; `repair_oifits(file, out)` writes a standard OIFITS v2 copy;
`oifits_headers(file)` lists every HDU and its header cards.

## OIFITSv1 compatibility

Files following the older OIFITS v1 standard (0-based station indices) are loaded
automatically. A warning is printed indicating that station indices have been
reindexed from 1 internally.
