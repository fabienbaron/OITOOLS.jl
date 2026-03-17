# Running examples

Example scripts are provided in the `demos/` folder and are meant to be run
from inside that directory. Data files live in `demos/data/`.

## Key demos

### Data I/O and plotting

| Script | Description |
|--------|-------------|
| `example_plot_observables.jl` | Load data and produce UV coverage, V², T3φ plots |
| `example_filter_data.jl` | Post-load filtering by baseline and wavelength range |
| `example_identify_oifits_points.jl` | Interactive identification of data points |
| `example_chi2_correlated.jl` | Correlated χ² using OI_CORR from GRAVITY data |

### Model fitting

| Script | Description |
|--------|-------------|
| `example_model_fitting_limb_darkening.jl` | Fit limb-darkened disk models to V² data |
| `example_model_fitting_v1295aql.jl` | YSO disk fitting (Hankel profiles) — Ibrahim et al. 2023 |
| `example_model_fitting_pmoired_models.jl` | Gallery of PMOIRED-compatible model types |
| `example_model_fitting_pmoired_conversion.jl` | Convert PMOIRED Python dicts to Julia |
| `example_model_fitting_2004.jl` | Model fitting on Beauty Contest 2004 data |
| `example_model_fitting_functions.jl` | Comparing visibility functions and fitting backends |
| `example_grid_search_limb_darkening.jl` | Grid search over limb-darkening parameters |
| `example_companion_gridsearch.jl` | Binary companion grid search |
| `example_companion_fastgridsearch.jl` | Fast binary companion detection |
| `example_binary_gridsearch.jl` | Binary system grid search |
| `example_bootstrap_fit.jl` | Bootstrap error estimation |

### Image reconstruction

| Script | Description |
|--------|-------------|
| `example_image_reconstruction_basic.jl` | Monochromatic reconstruction with NFFT |
| `example_image_reconstruction_polychromatic_MWC480.jl` | Polychromatic reconstruction with spectral regularization |
| `example_image_reconstruction_betaLyrae.jl` | CHARA/MIRC β Lyrae reconstruction |
| `example_image_reconstruction_regularization_types.jl` | Compare regularization types |
| `example_image_reconstruction_regularization_weight.jl` | Regularization weight selection |
| `example_image_reconstruction_lcurve.jl` | L-curve for regularization weight |
| `example_image_reconstruction_support.jl` | Support-constrained reconstruction |
| `example_image_reconstruction_sparco_grey.jl` | SPARCO grey reconstruction |
| `example_image_reconstruction_sparco_grey_v1295.jl` | SPARCO grey for V1295 Aql |
| `example_image_reconstruction_oimem.jl` | Maximum entropy (BSMEM) reconstruction |
| `example_image_reconstruction_pigru.jl` | Pi Gru reconstruction |
| `example_dirty_beam_and_image.jl` | Dirty beam and dirty image computation |

### Simulation

| Script | Description |
|--------|-------------|
| `example_simulate_observations_from_OIFITS.jl` | Simulate data matching an existing OIFITS file |
| `example_simulate_observations_from_model.jl` | Simulate from a parametric model with array geometry |
| `example_simulate_observations_from_image.jl` | Simulate from a FITS image |
| `example_simulate_polychromatic_disk.jl` | Polychromatic disk with companion (OI_VIS + differential phases) |
| `example_chara_plan.jl` | CHARA observation planning and Gantt chart |
| `example_npoi_target_selection.jl` | NPOI target selection |

## Demo data

| File | Description |
|------|-------------|
| `2004-data1.oifits` | SPIE Beauty Contest 2004 synthetic dataset — used in most reconstruction examples |
| `2004true.fits` | Ground-truth image for the above |
| `betlyr6t.oifits` | CHARA/MIRC data — β Lyrae |
| `AlphaCenA.oifits` | VLTI/PIONIER data — used in model fitting examples |
| `MWC480.oifits` | Polychromatic dataset for multi-channel reconstruction |
| `rho_Cas_example.oifits` | Example file for plotting demos |
| `BC2026/OBJECT1_LM.oifits` | Beauty Contest 2026 — object 1, LM band |
| `BC2026/OBJECT1_N.oifits` | Beauty Contest 2026 — object 1, N band |
| `BC2026/OBJECT2_K.oifits` | Beauty Contest 2026 — object 2, K band (GRAVITY, with OI_CORR) |
