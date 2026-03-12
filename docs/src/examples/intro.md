# Running examples

Example scripts are provided in the `demos/` folder and are meant to be run
from inside that directory. Data files live in `demos/data/`.

## Key demos

| Script | Description |
|--------|-------------|
| `example_plot_observables.jl` | Load data and produce UV coverage, V², T3φ plots |
| `example_filter_data.jl` | Post-load filtering by baseline and wavelength range |
| `example_model_fitting_limb_darkening.jl` | Fit a limb-darkened disc to V² data |
| `example_image_reconstruction_basic.jl` | Monochromatic image reconstruction with NFFT |
| `example_image_reconstruction_polychromatic_MWC480.jl` | Polychromatic reconstruction with spectral regularization |
| `example_simulate_observations_from_OIFITS.jl` | Simulate data matching an existing OIFITS file |
| `example_chara_plan.jl` | CHARA observation planning and Gantt chart |

## Demo data

| File | Description |
|------|-------------|
| `2004-data1.oifits` | CHARA/MIRC dataset — used in most reconstruction examples |
| `2004true.fits` | Ground-truth image for the above |
| `AlphaCenA.oifits` | VLTI/PIONIER data — used in model fitting examples |
| `MWC480.oifits` | Polychromatic dataset for multi-channel reconstruction |
| `rho_Cas_example.oifits` | Example file for plotting demos |
