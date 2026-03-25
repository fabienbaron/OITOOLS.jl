# Hybrid reconstruction

Hybrid reconstruction combines a parametric model (handling compact/unresolved
components) with an image (capturing extended emission). The optimizer
alternates between VMLMB for image pixels and NelderMead for model parameters.

## Custom model + image

Use `reconstruct_hybrid` with a flat model dict for arbitrary parametric
model + image combinations:

```julia
using OITOOLS
data = readoifits("data/2019_v1295Aql.WL_SMOOTH.A.oifits";
                  polychromatic=true, filter_bad_data=true)
nx      = 64
pixsize = 0.25
ft      = setup_ft(data, nx, pixsize)
x0      = gaussian2d(nx, nx, nx / 6)

# Define a custom model dict
model_dict = Dict{String,Any}(
    "wl0"     => 1.6e-6,
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

params, x_img = reconstruct_hybrid(x0, params0, model, data[1,1], ft;
    w_name       = "W",
    regularizers = [["tv", 1e2]],
    rounds       = 3,
    maxiter      = 200,
    verb         = true)
```

The `w_name="W"` keyword tells the optimizer which model output gives the
image weight function W(λ). This is used to compute the chromatic flux
normalisation for the image gradient.

## When to use hybrid vs SPARCO

`reconstruct_sparco` is a convenience wrapper around `reconstruct_hybrid` that
builds the standard SPARCO model dict internally. Use `reconstruct_hybrid`
directly when:

- You need a non-standard chromatic model (e.g. multiple point sources with
  different spectral indices, or a custom flux decomposition)
- You want full control over the model dict and free parameters
- Your model does not fit the star + background + environment template

## Legacy joint optimisation

`reconstruct_sparco_flat` performs joint VMLMB optimisation over both image
pixels and model parameters simultaneously. This is the older approach --
`reconstruct_hybrid` (alternating VMLMB + NelderMead) is generally more
robust because it avoids scale mismatch between pixel and parameter gradients.

```julia
params, x_img = reconstruct_sparco_flat(x0, params0, model, data[1,1], ft;
    w_name       = "W",
    regularizers = [["tv", 1e2]],
    maxiter      = 500)
```

See the [Imaging (+ Modeling)](@ref) API reference for full docstrings.
