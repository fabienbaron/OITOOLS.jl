# SPARCO reconstruction

For targets with a spectrally-variable environment (e.g. T Tauri discs), the
SPARCO framework separates a grey circumstellar component from an unresolved
stellar point source. The high-level `reconstruct_sparco` takes physical
parameters directly:

```julia
using OITOOLS
data = readoifits("data/2019_v1295Aql.WL_SMOOTH.A.oifits";
                  polychromatic=true, filter_bad_data=true)
nx      = 64
pixsize = 0.25
ft      = setup_ft(data, nx, pixsize)
x0      = gaussian2d(nx, nx, nx / 6)

result = reconstruct_sparco(x0, data[1,1], ft;
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
println("χ²/ν = ", result.chi2)
```

The `result` is a named tuple with fields `params`, `image`, `model`,
`free_params`, `chi2`, and `param_names`.

## Keyword reference

| Keyword | Default | Description |
|---------|---------|-------------|
| `lambda_ref` | `1.6e-6` | Reference wavelength (m) |
| `star_flux` | `0.5` | Stellar flux fraction at `lambda_ref` |
| `bg_flux` | `0.0` | Background flux fraction at `lambda_ref` |
| `d_env` | `0.0` | Environment spectral index |
| `star_di` | `4.0` | Stellar spectral index (Rayleigh–Jeans default) |
| `star_spectral_index_fixed` | `true` | If `false`, fit the stellar spectral index too |
| `fit_model_first` | `true` | Run a model-only fit before image reconstruction |
| `rounds` | `3` | Number of alternating image/parameter rounds |
| `regularizers` | `[["centering", 1e4]]` | Image regularizers |
| `weights` | `[1,1,1,0,0,0,0]` | Observable weights `[V², T3amp, T3φ, ...]` |

See `example_image_reconstruction_sparco_grey_v1295_newAPI.jl`.

See the [Imaging (SPARCO)](@ref) API reference for full docstrings.
