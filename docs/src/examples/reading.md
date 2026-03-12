# Reading OIFITS files

## Basic usage

`readoifits` loads an OIFITS file and returns a 2D array of `OIdata{T}` indexed as
`[nwavbin, ntimebin]`. The simplest call reads everything into a single bin:

```julia
using OITOOLS
data = readoifits("mystar.oifits")        # returns Array{OIdata{Float64}, 2}
data[1,1]                                 # the single bin
```

Use `T=Float32` to halve memory usage for large datasets:

```julia
data = readoifits("mystar.oifits"; T=Float32)
```

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
| `nv2`, `nt3phi`, `nt3amp`, `nvisamp`, `nvisphi`, `nflux`, `nuv` | Counts |

Correlation matrices from `OI_CORR` tables (empty sparse when absent):

| Field | Description |
|-------|-------------|
| `v2_corr`, `v2_corr_idx` | V² correlation matrix and per-point index |
| `t3phi_corr`, `t3phi_corr_idx` | Closure-phase correlation |
| `t3amp_corr`, `t3amp_corr_idx` | Triple-amplitude correlation |
| `visamp_corr`, `visamp_corr_idx` | Visibility amplitude correlation |
| `visphi_corr`, `visphi_corr_idx` | Visibility phase correlation |
| `flux_corr`, `flux_corr_idx` | Flux correlation |

`*_corr_idx[i] == 0` means point `i` has no associated correlation row.

Call `display(data)` or `display(data[1,1])` to print a compact summary.

## Selecting a target

Files with multiple targets require a `targetname` filter:

```julia
data = readoifits("multi_target.oifits"; targetname="Betelgeuse")
```

Use `list_oifits_targets("multi_target.oifits")` to see all available targets.

## Polychromatic (per-channel) mode

To keep each spectral channel as a separate bin, pass `polychromatic=true`. The
bin boundaries are derived automatically from midpoints between adjacent channel
centres, so each stored λ value falls in exactly one bin:

```julia
data = readoifits("my_spectrum.oifits"; polychromatic=true)
# data is nchannels×1; data[k,1] is the k-th channel bin
```

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

## UV deduplication

Redundant UV points (within `uvtol` metres by default 200 m) are merged:

```julia
data = readoifits("mystar.oifits"; redundance_remove=true, uvtol=100.0)
```

Disable with `redundance_remove=false`.

## Post-load filtering

After loading, use `set_data_filter` + `filter_data` to apply additional cuts
without re-reading the file:

```julia
idx = set_data_filter(data[1,1];
    filter_bad_data      = true,
    baseline_range       = [10.0e6, 300.0e6],
    wav_range            = [1.6e-6, 1.8e-6])
data_filtered = filter_data(data[1,1], idx)
```

## Reading multiple files

**Multi-epoch** — each file is one epoch; returns `(nepochs, tepochs, data)`:

```julia
files = ["epoch1.oifits", "epoch2.oifits", "epoch3.oifits"]
nepochs, tepochs, data = readoifits_multiepochs(files)
```

**Multi-colour** — each file is a different waveband; returns a vector of `OIdata`:

```julia
files = ["H_band.oifits", "K_band.oifits"]
data = readoifits_multicolors(files)
```

## OIFITSv1 compatibility

Files following the older OIFITS v1 standard (0-based station indices) are loaded
automatically. A warning is printed indicating that station indices have been
reindexed from 1 internally.

## Internal pipeline

For reference, `readoifits` calls the following internal functions in order:

1. `load_oifits` — opens the file via OIFITS.jl (falls back to `hack_revn=1` for v1 files)
2. `collect_station_info` — builds a sparse station-index conversion table from all `OI_ARRAY` tables
3. `read_flux/vis/v2/t3_tables` — flattens all relevant HDUs into merged NamedTuples
4. `slice_to_bin` — selects the spectral/temporal window and assembles the UV plane
5. `filter_bad_observables!` — quality cuts, in-place
6. `remove_redundant_uv!` — KDTree-based UV deduplication, in-place
7. `make_oidata` — packages everything into the final `OIdata{T}` struct
