# Hybrid reconstruction

A parametric model for the compact components plus an image for the extended emission. The
image (VMLMB) and the model parameters (Nelder-Mead) are optimised in alternation.

## Custom model + image

`reconstruct_hybrid` takes any model dict:

```julia
using OITOOLS
data    = readoifits("data/V1295_Aql.oifits"; filter_bad_data=true)[1,1]
nx      = 64
pixsize = 0.125
ft      = setup_nfft(data, nx, pixsize)
x0      = gaussian2d(nx, nx, nx / 6)

model_dict = Dict{String,Any}(
    "wl0"     => 1.6,                # µm, as $WL
    "star,f0" => 0.5,
    "star,di" => 4.0,
    "bg,f0"   => 0.0,
    "d_env"   => 0.0,
    "D"       => raw"$star,f0 * ($WL/$wl0)^(-$star,di) + $bg,f0 * ($WL/$wl0)^(-4) + (1 - $star,f0 - $bg,f0) * ($WL/$wl0)^($d_env)",
    "star,f"  => raw"$star,f0 * ($WL/$wl0)^(-$star,di) / $D",
    "bg,f"    => raw"$bg,f0 * ($WL/$wl0)^(-4) / $D",
    "W"       => raw"(1 - $star,f0 - $bg,f0) * ($WL/$wl0)^($d_env) / $D",
    "bg,resolved" => true,
)
free_params = ["star,f0", "bg,f0", "d_env"]

model   = parse_model(model_dict, free_params)
params0 = Float64[model_dict[p] for p in free_params]

params, x_img = reconstruct_hybrid(x0, params0, model, data, ft;
    w_name       = "W",
    regularizers = [["tv", 1e2]],
    rounds       = 3,
    maxiter      = 200,
    verb         = true)
```

`w_name = "W"` names the model output that gives the image's weight W(λ), so that
`V = V_model + W(λ) · V_image`.

`reconstruct_sparco` builds the standard star + background + environment dict and calls
`reconstruct_hybrid`. `reconstruct_sparco_flat` (joint VMLMB over pixels and parameters) is
deprecated.

The docstrings are in [Imaging (+ Modeling)](@ref).
