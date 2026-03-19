##############################################################################
#  example_model_fitting_limb_darkening_new.jl
#
#  Showcases the new flat-dict model fitting interface on Alpha Centauri A
#  PIONIER data (Kervella et al. 2017, arXiv:1610.06185).
#
#  Demonstrates:
#   1. display_model — inspect model before fitting
#   2. fit_model — NLopt with gradient (LD_LBFGS)
#   3. fit_model — gradient-free (LN_NELDERMEAD)
#   4. model_to_chi2 — compute chi2 at arbitrary parameter values
#   5. model_to_obs + plot_v2_residuals — compare model vs data
#   6. model_to_image + imdisp — visualise the fitted model
#   7. fit_model_ultranest — Bayesian nested sampling (optional, slow)
#
#  Reference results (Kervella+ 2017, H band):
#   Alpha Cen A: θ_LD = 8.502 ± 0.038 mas, α = 0.1404 ± 0.0050
#
##############################################################################

using OITOOLS
using Printf

# ── Load data ─────────────────────────────────────────────────────────────────
oifitsfile = joinpath(@__DIR__, "data", "AlphaCenA.oifits")
data = readoifits(oifitsfile)[1,1]

# ===========================================================================
# 1. Uniform disc fit (single parameter)
# ===========================================================================

println("\n" * "="^70)
println("  1. Uniform disc fit")
println("="^70)

param_ud = Dict{String,Any}(
    "star,ud" => 8.0,    # starting diameter (mas)
    "star,f"  => 1.0,    # flux fraction (fixed)
)
fit_params_ud = ["star,ud"]
lb_ud = Dict("star,ud" => 5.0)
ub_ud = Dict("star,ud" => 12.0)

# Inspect the model setup
display_model(param_ud, fit_params_ud; lb=lb_ud, ub=ub_ud)

# Fit with gradient-based optimizer (V² only)
result_ud = fit_model(param_ud, fit_params_ud, data;
    lb=lb_ud, ub=ub_ud, weights=[1.0, 0, 0])

println(result_ud)
@printf("  UD diameter = %.4f mas,  χ²r = %.3f\n",
        result_ud.x_opt[1], result_ud.chi2r)


# ===========================================================================
# 2. Power-law limb-darkened disc (2 parameters, gradient-based)
# ===========================================================================

println("\n" * "="^70)
println("  2. Power-law LD fit (LD_LBFGS)")
println("="^70)

param_ld = Dict{String,Any}(
    "star,ldpow" => result_ud.x_opt[1],  # seed from UD fit
    "star,alpha" => 0.15,                 # LD exponent
    "star,f"     => 1.0,
)
fit_params_ld = ["star,ldpow", "star,alpha"]
lb_ld = Dict("star,ldpow" => 5.0,  "star,alpha" => 0.0)
ub_ld = Dict("star,ldpow" => 12.0, "star,alpha" => 0.5)

display_model(param_ld, fit_params_ld; lb=lb_ld, ub=ub_ld)

result_ld = fit_model(param_ld, fit_params_ld, data;
    lb=lb_ld, ub=ub_ld, weights=[1.0, 0, 0])

println(result_ld)
@printf("  Fitted:    θ_LD = %.4f mas,  α = %.4f,  χ²r = %.3f\n",
        result_ld.x_opt[1], result_ld.x_opt[2], result_ld.chi2r)
@printf("  Published: θ_LD = 8.502 ± 0.038 mas,  α = 0.1404 ± 0.0050\n")


# ===========================================================================
# 3. Same model, gradient-free optimizer (Nelder-Mead)
# ===========================================================================

println("\n" * "="^70)
println("  3. Power-law LD fit (LN_NELDERMEAD)")
println("="^70)

result_nm = fit_model(param_ld, fit_params_ld, data;
    lb=lb_ld, ub=ub_ld, weights=[1.0, 0, 0],
    method=:LN_NELDERMEAD, maxeval=5000)

println(result_nm)
@printf("  Nelder-Mead: θ_LD = %.4f,  α = %.4f,  χ²r = %.3f\n",
        result_nm.x_opt[1], result_nm.x_opt[2], result_nm.chi2r)


# ===========================================================================
# 4. Manual chi2 evaluation
# ===========================================================================

println("\n" * "="^70)
println("  4. Manual χ² evaluation")
println("="^70)

# Compile model and evaluate chi2 at published values
model_ld = parse_model(param_ld, fit_params_ld)
chi2_pub = model_to_chi2(model_ld, [8.502, 0.1404], data; weights=[1.0, 0, 0])
@printf("  χ² at published values (8.502, 0.1404): %.2f\n", chi2_pub)

chi2_fit = model_to_chi2(model_ld, result_ld.x_opt, data; weights=[1.0, 0, 0])
@printf("  χ² at best-fit values:                  %.2f\n", chi2_fit)


# ===========================================================================
# 5. Model observables and residual plot
# ===========================================================================

println("\n" * "="^70)
println("  5. Model observables & residual plot")
println("="^70)

obs = model_to_obs(result_ld.model, result_ld.x_opt, data)
@printf("  V² model: %d points, range [%.4f, %.4f]\n",
        length(obs.v2), minimum(obs.v2), maximum(obs.v2))

# Plot V² residuals (data - model)
plot_v2_residuals(data, obs.v2; logplot=true)


# ===========================================================================
# 6. Model image
# ===========================================================================

println("\n" * "="^70)
println("  6. Model image synthesis")
println("="^70)

img = model_to_image(result_ld.model, result_ld.x_opt;
                     nx=128, pixsize=0.1, oversample=2)
@printf("  Image size: %d × %d,  peak = %.4e,  sum = %.4f\n",
        size(img)..., maximum(img), sum(img))

imdisp(img, pixsize=0.1)


# ===========================================================================
# 7. Linear LD law (ldlin) — alternative model
# ===========================================================================

println("\n" * "="^70)
println("  7. Linear LD fit")
println("="^70)

param_ldlin = Dict{String,Any}(
    "star,ldlin" => result_ud.x_opt[1],
    "star,u"     => 0.3,    # linear LD coefficient
    "star,f"     => 1.0,
)
fit_params_ldlin = ["star,ldlin", "star,u"]
lb_ldlin = Dict("star,ldlin" => 5.0, "star,u" => 0.0)
ub_ldlin = Dict("star,ldlin" => 12.0, "star,u" => 1.0)

result_ldlin = fit_model(param_ldlin, fit_params_ldlin, data;
    lb=lb_ldlin, ub=ub_ldlin, weights=[1.0, 0, 0])

println(result_ldlin)
@printf("  Linear LD: θ = %.4f mas,  u = %.4f,  χ²r = %.3f\n",
        result_ldlin.x_opt[1], result_ldlin.x_opt[2], result_ldlin.chi2r)


# ===========================================================================
# 8. UltraNest nested sampling (commented out — slow but powerful)
# ===========================================================================

println("\n" * "="^70)
println("  8. UltraNest nested sampling")
println("="^70)

result_un = fit_model_ultranest(param_ld, fit_params_ld, data;
    lb=lb_ld, ub=ub_ld, weights=[1.0, 0, 0],
    min_num_live_points=200, verb=true, cornerplot=true)

println(result_un)
@printf("  UltraNest: θ_LD = %.4f,  α = %.4f,  log(Z) = %.2f ± %.2f\n",
        result_un.x_opt[1], result_un.x_opt[2],
        result_un.logz, result_un.logzerr)

println("\nDone.")
