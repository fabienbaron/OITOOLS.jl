# SPARCO reconstruction

SPARCO models the target as a star (a point source with a power-law spectrum), an optional
over-resolved background, and a grey image for the environment with its own spectral index.
The wavelengths are kept in one `OIdata`; the chromatism is in the model. `reconstruct_sparco`
takes the physical parameters:

```julia
using OITOOLS, PythonPlot
data    = readoifits("data/V1295_Aql.oifits"; filter_bad_data=true)[1,1]
nx      = 64
pixsize = 0.125
ft      = setup_nfft(data, nx, pixsize)
x0      = gaussian2d(nx, nx, nx / 6)

result = reconstruct_sparco(x0, data, ft;
    lambda_ref  = 1.6e-6,
    star_flux   = 0.48,
    bg_flux     = 0.0,
    d_env       = -0.1,
    weights     = [1.0, 0.0, 1.0],
    regularizers = [["tv", 1e2]],
    maxiter     = 200,
    rounds      = 3,
    verb        = true)

imdisp(result.image; pixsize)
println("Final params: ", result.params)
println("χ² = ", result.chi2)
```

`result` is a named tuple with fields `params`, `image`, `model`, `free_params`, `chi2` and
`param_names` (a `Dict` of the fitted values, `star,f0`, `bg,f0`, `bg,frac`, `d_env`).

The background is fitted as `bg,frac`, its fraction of the flux the star leaves
(`bg,f0 = (1 − star,f0) · bg,frac`). With `star,f0` and `bg,frac` in [0, 1], the star,
background and image fluxes at `lambda_ref` are each in [0, 1] and sum to 1.

## Keyword reference

| Keyword | Default | Description |
|---------|---------|-------------|
| `lambda_ref` | required | Reference wavelength (m) |
| `star_flux` | `0.5` | Stellar flux fraction at `lambda_ref` |
| `bg_flux` | `0.0` | Background flux fraction at `lambda_ref` (`star_flux + bg_flux ≤ 1`) |
| `d_env` | `0.0` | Environment spectral index |
| `star_di` | `4.0` | Stellar spectral index (Rayleigh–Jeans default) |
| `star_spectral_index_fixed` | `true` | If `false`, fit the stellar spectral index too |
| `fit_model_first` | `false` | Fit the parameters with the starting image before the first image step |
| `rounds` | `3` | Number of alternating image/parameter rounds |
| `regularizers` | `[]` | Image regularizers |
| `weights` | `[1, 1, 1]` | Observable weights `[V², T3amp, T3φ]` |
| `maxiter` | `200` | VMLMB iterations per image step |

See `demos/example_image_reconstruction_sparco_grey_v1295.jl`, and the
[Imaging (SPARCO)](@ref) API reference.
