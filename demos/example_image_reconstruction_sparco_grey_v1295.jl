# SPARCO image reconstruction of V1295 Aql (HD 190073) using flat model API
#
# Polychromatic SPARCO: parametric star + resolved background + grey environment image
# with chromatic flux normalization. The image weight W(λ) accounts for spectral
# differences between star and environment.
#
# Physical model:
#   - Star:        point source, flux fraction f_star, spectrum ∝ (λ/λ0)^(-4)
#   - Background:  over-resolved (V=0), flux fraction f_bg, spectrum ∝ (λ/λ0)^(-4)
#   - Environment: grey NFFT image, flux fraction 1 - f_star - f_bg,
#                  spectrum ∝ (λ/λ0)^(d_env)
#
# Normalization: V_total = (f_star·V_star·α + f_env·V_image·β) / (f_star·α + f_bg·α + f_env·β)
# where α = (λ/λ0)^(-4) and β = (λ/λ0)^(d_env).
#
# Note: T3amp is disabled (weights=[1, 0, 1]) — this dataset has noisy T3amp.

using OITOOLS
using Statistics, Printf

# ═════════════════════════════════════════════════════════════════════════════
# Load data and setup
# ═════════════════════════════════════════════════════════════════════════════

oifitsfile = joinpath(@__DIR__, "data", "2019_v1295Aql.WL_SMOOTH.A.oifits")
data = readoifits(oifitsfile, filter_bad_data=true)[1,1]

pixsize = 0.125  # mas/pixel
nx = 64
ft = setup_nfft(data, nx, pixsize)

λ0 = 1.6e-6  # H-band reference wavelength (m)
println("Data: nv2=$(data.nv2)  nt3=$(data.nt3amp)  nuv=$(data.nuv)")
println("λ range: ", extrema(data.uv_lam))

# ═════════════════════════════════════════════════════════════════════════════
# Define flat model
# ═════════════════════════════════════════════════════════════════════════════

model_dict = Dict{String,Any}(
    # Reference wavelength (fixed)
    "wl0"     => λ0,

    # Star: point source at origin (no x, y, ud keys → point source)
    "star,f0" => 0.5,                # flux fraction at λ0
    "star,di" => 4.0,                # spectral index (fixed, standard star)

    # Background: over-resolved (V=0), same spectrum as star
    "bg,resolved" => true,
    "bg,f0"       => 0.0,            # flux fraction at λ0

    # Environment spectral index (d_env=0 → grey)
    "d_env"   => 0.0,

    # Normalization denominator: sum of all chromatic flux contributions
    "D" => "\$star,f0 * (\$WL/\$wl0)^(-\$star,di) + " *
           "\$bg,f0 * (\$WL/\$wl0)^(-4) + " *
           "(1 - \$star,f0 - \$bg,f0) * (\$WL/\$wl0)^(\$d_env)",

    # Normalized chromatic flux fractions
    "star,f" => "\$star,f0 * (\$WL/\$wl0)^(-\$star,di) / \$D",
    "bg,f"   => "\$bg,f0 * (\$WL/\$wl0)^(-4) / \$D",

    # Image weight W(λ): chromatic weight for the grey environment image
    "W" => "(1 - \$star,f0 - \$bg,f0) * (\$WL/\$wl0)^(\$d_env) / \$D",
)

# Free parameters to optimize
free_params = ["star,f0", "bg,f0", "d_env"]

# Weights: V² and T3phi only (T3amp disabled — noisy for this dataset)
weights = [1.0, 0.0, 1.0]

# ═════════════════════════════════════════════════════════════════════════════
# Phase 0: Fit model first (star dominates the data)
# ═════════════════════════════════════════════════════════════════════════════

println("\n=== Phase 0: model-only fit (star dominates) ===")
result = fit_model(model_dict, free_params, data;
    weights=[weights; 0.0; 0.0; 0.0; 0.0],
    method=:LD_LBFGS,
    lb=Dict("star,f0" => 0.01, "bg,f0" => 0.0, "d_env" => -10.0),
    ub=Dict("star,f0" => 0.99, "bg,f0" => 0.5, "d_env" =>  10.0),
    maxeval=5000)
println(result)

# Use model-fit parameters as starting point
x0 = result.x_opt
model = parse_model(model_dict, free_params)
println("Starting params from model fit: ", collect(zip(free_params, x0)))

# ═════════════════════════════════════════════════════════════════════════════
# Initial chi2 with random image
# ═════════════════════════════════════════════════════════════════════════════

x_start = rand(nx, nx)
x_start ./= sum(x_start)
println("\n=== Initial chi2 ===")
chi2 = chi2_sparco_flat_f(x_start, x0, model, ft, data; w_name="W", weights=weights, verb=true)
@printf("\nTotal chi2 = %.2f\n", chi2)

# ═════════════════════════════════════════════════════════════════════════════
# Reconstruction: alternating image + parameter optimization
# ═════════════════════════════════════════════════════════════════════════════

# Bounds for parameters: [star,f0, bg,f0, d_env]
lb_params = [0.01, 0.0, -10.0]
ub_params = [0.99, 0.5,  10.0]

# ── Phase 1: no regularization ────────────────────────────────────────────
println("\n=== Phase 1: no regularization (200 iters) ===")
params, x = reconstruct_sparco_flat(x_start, x0, model, data, ft;
    w_name="W", regularizers=[], weights=weights, verb=true, maxiter=200,
    params_lower=lb_params, params_upper=ub_params)

println("\n  Optimizing parameters...")
minchi2, params, ret = optimize_sparco_flat_parameters(params, model, x, ft, data;
    w_name="W", weights=weights, lb=lb_params, ub=ub_params)
@printf("  NelderMead: chi2=%.2f  ret=%s\n", minchi2, ret)
for (i, k) in enumerate(free_params)
    @printf("    %10s = %.6f\n", k, params[i])
end

params, x = reconstruct_sparco_flat(x, params, model, data, ft;
    w_name="W", regularizers=[], weights=weights, verb=true, maxiter=200,
    params_lower=lb_params, params_upper=ub_params)

# ── Phase 2: decreasing TV regularization ─────────────────────────────────
for (phase, mu_tv) in enumerate([1e7, 1e5, 1e3, 1e1])
    global params, x, minchi2, ret
    @printf("\n=== Phase %d: tvsq μ=%.0e (200 iters) ===\n", phase + 1, mu_tv)
    regularizers = [["tvsq", mu_tv]]

    params, x = reconstruct_sparco_flat(x, params, model, data, ft;
        w_name="W", regularizers=regularizers, weights=weights, verb=true, maxiter=200,
        params_lower=lb_params, params_upper=ub_params)

    println("\n  Optimizing parameters...")
    minchi2, params, ret = optimize_sparco_flat_parameters(params, model, x, ft, data;
        w_name="W", weights=weights, lb=lb_params, ub=ub_params)
    @printf("  NelderMead: chi2=%.2f  ret=%s\n", minchi2, ret)
    for (i, k) in enumerate(free_params)
        @printf("    %10s = %.6f\n", k, params[i])
    end
end

# ═════════════════════════════════════════════════════════════════════════════
# Final results
# ═════════════════════════════════════════════════════════════════════════════

println("\n=== Final results ===")
chi2_final = chi2_sparco_flat_f(x, params, model, ft, data;
    w_name="W", weights=weights, verb=true)
@printf("\nFinal chi2 = %.2f\n", chi2_final)
println("Final params:")
for (i, k) in enumerate(free_params)
    @printf("  %10s = %.6f\n", k, params[i])
end
@printf("  %10s = %.6f  (derived: 1 - star,f0 - bg,f0)\n", "f_env", 1 - params[1] - params[2])
println("Image flux = ", sum(x))

imdisp(x, pixsize=pixsize, figtitle="V1295 Aql — SPARCO flat model")
