# API Reference

This page lists all user-facing functions, types, and constants in OITOOLS.

## Reading OIFITS data

```@docs
OIdata
readoifits
readoifits_multiepochs
readoifits_multicolors
list_oifits_targets
readfits
writefits
oifits_prep
updatefits_aspro
remove_redundant_uv!
filter_data
set_data_filter
```

## Plotting

```@docs
set_oiplot_defaults
uvplot
onclickidentify
plot_v2
plot_t3phi
plot_t3amp
plot_v2_and_t3phi_wav
plot_visphi
plot_diffphi
plot_flux
plot_v2_residuals
plot_t3phi_residuals
plot_t3amp_residuals
plot_v2_multifile
imdisp
imdisp_temporal
imdisp_polychromatic
```

## Model definition and parsing

```@docs
FlatModel
parse_model
eval_model
eval_model_grad
```

## Model fitting

```@docs
fit_model
fit_model_lsqfit
fit_model_ultranest
FitResult
LsqFitResult
UltraNestResult
display_model
model_to_obs
model_to_image
model_to_sed
resample_data
chi2_flat
chi2_flat_fg
residuals_flat
residuals_flat_jac
```

## Visibility functions

```@docs
visibility_ud
visibility_ldlin
visibility_ldquad
visibility_ldquad_alt
visibility_ldpow
visibility_ldsquareroot
visibility_annulus
visibility_ellipse_uniform
visibility_ellipse_quad
visibility_thin_ring
visibility_Gaussian_ring
visibility_Gaussian_ring_az
visibility_Lorentzian_ring
visibility_GaussianLorentzian_ring_az
bb
```

### AD-compatible visibility functions

```@docs
vis_ud
vis_ldlin
vis_ldquad
vis_ldpow
visibility_ud_d
visibility_ldlin_d
visibility_ldquad_d
visibility_ldpow_d
```

## Hankel transforms

```@docs
hankel_transform
hankel_norm
hankel_vis
hankel_vis_interp
hankel_vis_full
eval_profile
compile_profile
trapz
trapz_weights
```

## Image reconstruction (chi2 / regularization)

```@docs
setup_dft
setup_dft_polychromatic
setup_nfft
setup_nfft_multiepochs
setup_nfft_polychromatic
image_to_vis_dft
image_to_vis
image_to_v2
image_to_t3phi
image_to_t3amp
image_to_obs
observables
chi2_f
chi2_fg
chi2_vis_nfft_f
chi2_vis_dft_fg
chi2_vis_nfft_fg
crit_f
crit_fg
crit_polychromatic_fg
crit_multitemporal_fg
reconstruct
reconstruct_multitemporal
reconstruct_polychromatic
gaussian2d
cdg
mod360
```

### Regularization

```@docs
regularization
reg_centering
tvsq
tv
l1l2
l1l2w
l1hyp
l2sq
compactness
radial_variance
entropy
reg_support
setup_radial_reg
```

### SPARCO

```@docs
chi2_sparco_f
chi2_sparco_fg
reconstruct_sparco_gray
optimize_sparco_parameters
disk
```

## Maximum entropy (BSMEM)

```@docs
reconstruct_bsmem
maxent_reconstruct!
maxent_setup
reconstruct!
maxent_step!
MaximENTParams
auto_pixsize
gaussian_prior
```

## Simulation

```@docs
simulate
simulate_from_oifits
get_uv
get_uv_indxes
prep_arrays
read_array_file
read_obs_file
read_comb_file
read_wave_file
get_v2_baselines
get_t3_baselines
v2mapt3
hour_angle_calc
hours_to_date
sunrise_sunset
mjd_to_utdate
dates_to_jd
jd_to_hour_angle
opd_limits
alt_az
geometric_delay
cart_delay
query_target_from_simbad
ra_dec_from_simbad
get_baselines
recenter
gantt_onenight
FacilityConfig
TargetConfig
CombinerConfig
WaveConfig
facility_info
obsv_info
combiner_info
wave_info
read_facility_file
```

## Writing OIFITS

```@docs
write_oi_header
write_oi_array
write_oi_target
write_oi_wavelength
write_oi_vis2
write_oi_t3
oifits_check
oifits_merge
oifits_filter
```

## PMOIRED compatibility

```@docs
pmoired_to_julia
pmoired_to_julia_file
```

## Utilities

```@docs
gaussianwrapped_to_vonmises_fast
logbesselI0
vis_to_v2
vis_to_t3
```
