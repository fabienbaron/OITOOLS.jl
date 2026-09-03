##############################################################################
#  example_model_fitting_limb_darkening.jl
#
#  Limb-darkened disc fitting on the VLTI/PIONIER observations of α Centauri
#  A and B published by Kervella, Bigot, Gallenne & Thévenin (2017),
#  A&A 597, A137 — the paper shipped alongside this repository as
#  `aa29505-16.pdf`.
#
#  Part I  — tutorial: the flat-dict model interface on α Cen A.
#  Part II — validation: every limb-darkening law OITOOLS implements, fitted
#            with every model-fitting algorithm OITOOLS offers, on both stars,
#            against Table 3 of the paper.
#
#  ── The wavelength scale, which you must apply before comparing ────────────
#
#  The OIFITS files here carry the raw PIONIER wavelength scale. Kervella et al.
#  calibrate it against the binary HD 123999 and multiply it by
#
#       γ = 1.00481 ± 0.00412                              (their Eq. 1)
#
#  before quoting any angular diameter (§2.2.4). An angular diameter scales
#  linearly with the assumed wavelength, so a diameter fitted from these files
#  must be multiplied by γ to be compared with Table 3. Limb-darkening
#  coefficients are dimensionless and are unaffected.
#
#  Skipping this makes OITOOLS look 0.04 mas low on every law at once. It is a
#  property of the data, not of the fit.
#
#  Run headless with:  MPLBACKEND=Agg julia --project=.. <this file>
##############################################################################

using OITOOLS
using Printf, Statistics
# The pure-Julia nested sampler. Add `using PythonCall` as well to have the
# comparison below run UltraNest too and check the two against each other.
using Nautilus

const DATADIR = joinpath(@__DIR__, "data")
const GAMMA   = 1.00481          # PIONIER wavelength scaling, Kervella+2017 Eq. 1
const V2ONLY  = [1.0, 0, 0]      # the paper fits squared visibilities only

data = readoifits(joinpath(DATADIR, "AlphaCenA.oifits"); warn=false, verbose=false)[1,1]

##############################################################################
#  PART I — tutorial
##############################################################################

# ===========================================================================
# 1. Uniform disc fit (single parameter)
# ===========================================================================

println("\n" * "="^78)
println("  1. Uniform disc fit")
println("="^78)

model_dict_ud = Dict{String,Any}(
    "star,ud" => 8.0,    # starting diameter (mas)
    "star,f"  => 1.0,    # flux fraction (fixed)
)
list_free_params_ud = ["star,ud"]
lb_ud = Dict("star,ud" => 5.0)
ub_ud = Dict("star,ud" => 12.0)

display_model(model_dict_ud, list_free_params_ud; lb=lb_ud, ub=ub_ud)

result_ud = fit_model(model_dict_ud, list_free_params_ud, data;
    lb=lb_ud, ub=ub_ud, weights=V2ONLY)

println(result_ud)
@printf("  UD diameter = %.4f mas  (×γ = %.4f)   χ²r = %.3f\n",
        result_ud.x_opt[1], GAMMA*result_ud.x_opt[1], result_ud.chi2r)
@printf("  Kervella+2017 Table 3:  θ_UD = 8.347 ± 0.004 ± 0.035 mas,  χ²r = 15.23\n")


# ===========================================================================
# 2. Power-law limb-darkened disc (2 parameters, gradient-based)
# ===========================================================================

println("\n" * "="^78)
println("  2. Power-law LD fit (LD_LBFGS)")
println("="^78)

model_dict_ld = Dict{String,Any}(
    "star,ldpow" => result_ud.x_opt[1],  # seed from UD fit
    "star,alpha" => 0.15,                # LD exponent
    "star,f"     => 1.0,
)
list_free_params_ld = ["star,ldpow", "star,alpha"]
lb_ld = Dict("star,ldpow" => 5.0,  "star,alpha" => 0.0)
ub_ld = Dict("star,ldpow" => 12.0, "star,alpha" => 0.5)

display_model(model_dict_ld, list_free_params_ld; lb=lb_ld, ub=ub_ld)

result_ld = fit_model(model_dict_ld, list_free_params_ld, data;
    lb=lb_ld, ub=ub_ld, weights=V2ONLY)

println(result_ld)
@printf("  Fitted:    θ_LD = %.4f mas (×γ = %.4f),  α = %.4f,  χ²r = %.3f\n",
        result_ld.x_opt[1], GAMMA*result_ld.x_opt[1], result_ld.x_opt[2], result_ld.chi2r)
@printf("  Published: θ_LD = 8.502 ± 0.006 ± 0.036 mas,  α = 0.1404 ± 0.0050\n")


# ===========================================================================
# 3. Same model, gradient-free optimizer (Nelder-Mead)
# ===========================================================================

println("\n" * "="^78)
println("  3. Power-law LD fit (LN_NELDERMEAD)")
println("="^78)

result_nm = fit_model(model_dict_ld, list_free_params_ld, data;
    lb=lb_ld, ub=ub_ld, weights=V2ONLY,
    method=:LN_NELDERMEAD, maxeval=5000)

println(result_nm)
@printf("  Nelder-Mead: θ_LD = %.4f,  α = %.4f,  χ²r = %.3f\n",
        result_nm.x_opt[1], result_nm.x_opt[2], result_nm.chi2r)
@printf("  vs LD_LBFGS: Δθ = %+.2e mas,  Δα = %+.2e,  Δχ² = %+.2e\n",
        result_nm.x_opt[1]-result_ld.x_opt[1], result_nm.x_opt[2]-result_ld.x_opt[2],
        result_nm.chi2-result_ld.chi2)


# ===========================================================================
# 4. Manual chi2 evaluation
# ===========================================================================

println("\n" * "="^78)
println("  4. Manual χ² evaluation")
println("="^78)

model_ld = dict_to_model(model_dict_ld, list_free_params_ld)
@printf("  χ² at the published values, uncorrected (8.5020, 0.1404): %8.2f\n",
        model_to_chi2(model_ld, [8.502, 0.1404], data; weights=V2ONLY))
@printf("  χ² at the published values, θ/γ        (%.4f, 0.1404): %8.2f\n",
        8.502/GAMMA, model_to_chi2(model_ld, [8.502/GAMMA, 0.1404], data; weights=V2ONLY))
@printf("  χ² at the OITOOLS best fit             (%.4f, %.4f): %8.2f\n",
        result_ld.x_opt..., model_to_chi2(model_ld, result_ld.x_opt, data; weights=V2ONLY))
println("  (the first line is why γ matters: forgetting it costs ~870 in χ²)")


# ===========================================================================
# 5. Model observables and residual plot
# ===========================================================================

println("\n" * "="^78)
println("  5. Model observables & residual plot")
println("="^78)

obs = model_to_obs(result_ld.model, result_ld.x_opt, data)
@printf("  V² model: %d points, range [%.4f, %.4f]\n",
        length(obs.v2), minimum(obs.v2), maximum(obs.v2))

plot_v2_residuals(data, obs.v2; logplot=true)


# ===========================================================================
# 6. Model image
# ===========================================================================

println("\n" * "="^78)
println("  6. Model image synthesis")
println("="^78)

img = model_to_image(result_ld.model, result_ld.x_opt; nx=128, pixsize=0.1, oversample=2)
@printf("  Image size: %d × %d,  peak = %.4e,  sum = %.4f\n",
        size(img)..., maximum(img), sum(img))

imdisp(img, pixsize=0.1)


##############################################################################
#  PART II — validation against Kervella et al. 2017, Table 3
##############################################################################

# ── What the paper fits, and what OITOOLS can reproduce ─────────────────────
#
# Table 3 lists nine LD models. OITOOLS' flat-dict interface implements four of
# them (`_GEOMETRY_KEYS` in src/parse_model.jl):
#
#     uniform disc        -> "ud"
#     linear              -> "ldlin"   + "u"          I(μ)/I(1) = 1 - u(1-μ)
#     quadratic           -> "ldquad"  + "u", "w"     I(μ)/I(1) = 1 - u(1-μ) - w(1-μ)²
#     power law           -> "ldpow"   + "alpha"      I(μ)/I(1) = μ^α  (Hestroffer 1997)
#
# The other five are out of reach here:
#
#     square root   — available as the `ldsqrt` geometry key in the model dictionary;
#                     the old standalone `visibility_ldsquareroot` was never wired into
#                     parse_model, so no dict-based fit can select it.
#     4-parameter   — not implemented (the paper holds its coefficients fixed at
#                     Claret & Bloemen 2011 values rather than fitting them).
#     3D / scaled solar — tabulated intensity profiles, not parametric laws.
#
# The linear law appears twice in Table 3: once with u fixed to the Claret & Bloemen
# (2011) value, once with u free. Both are reproduced below, since fixing a parameter
# is just leaving it out of `list_free_params`.

# Kervella+2017 Table 3.  θ in mas as (value, σ_stat, σ_syst); NaN = not fitted there.
const PUBLISHED = Dict(
  "AlphaCenA" => Dict(
     "uniform"       => (θ=(8.347, 0.004, 0.035), chi2r=15.23, pars=NamedTuple()),
     "linear (u fixed)"=>(θ=(8.505, 0.003, 0.036), chi2r= 5.50, pars=(u=(0.2392, NaN),)),
     "linear (u free)" => (θ=(8.458, 0.005, 0.035), chi2r= 4.24, pars=(u=(0.1761, 0.0062),)),
     "quadratic"     => (θ=(8.451, 0.013, 0.035), chi2r= 4.25, pars=(u=(0.191, 0.026), w=(-0.031, 0.054))),
     "power law"     => (θ=(8.502, 0.006, 0.036), chi2r= 3.90, pars=(alpha=(0.1404, 0.0050),)),
  ),
  "AlphaCenB" => Dict(
     "uniform"       => (θ=(5.883, 0.003, 0.025), chi2r=14.26, pars=NamedTuple()),
     "linear (u fixed)"=>(θ=(6.003, 0.002, 0.025), chi2r= 4.69, pars=(u=(0.2698, NaN),)),
     "linear (u free)" => (θ=(5.962, 0.003, 0.025), chi2r= 2.89, pars=(u=(0.1907, 0.0048),)),
     # the paper fits no two-LD-parameter model to B: its coverage stops in the
     # second lobe, so it cannot constrain more than one LD parameter (§3.1.1).
     "quadratic"     => (θ=(NaN, NaN, NaN),        chi2r=NaN,   pars=(u=(NaN,NaN), w=(NaN,NaN))),
     "power law"     => (θ=(5.999, 0.004, 0.025), chi2r= 3.33, pars=(alpha=(0.1545, 0.0044),)),
  ),
)

# u fixed to the Claret & Bloemen (2011) H-band value, per star (Table 3, note b)
const U_CLARET = Dict("AlphaCenA" => 0.2392, "AlphaCenB" => 0.2698)

"""
Model dict / free-parameter list / bounds for one LD law, seeded near `θ0`.

`θ0` is per-star on purpose. The χ² of a well-resolved disc is strongly multimodal —
α Cen B is resolved into its third lobe — and starting a B fit at A's 8 mas drops the
uniform-disc fit into a spurious minimum at 10.7 mas with χ²r ≈ 803. That is a property
of the χ² surface, not a bug: seed the diameter within ~20% of the truth.
"""
function ld_case(law, θ0, star)
    if law == "uniform"
        Dict{String,Any}("s,ud"=>θ0, "s,f"=>1.0), ["s,ud"],
            Dict("s,ud"=>0.6θ0), Dict("s,ud"=>1.5θ0)
    elseif law == "linear (u fixed)"
        Dict{String,Any}("s,ldlin"=>θ0, "s,u"=>U_CLARET[star], "s,f"=>1.0), ["s,ldlin"],
            Dict("s,ldlin"=>0.6θ0), Dict("s,ldlin"=>1.5θ0)
    elseif law == "linear (u free)"
        Dict{String,Any}("s,ldlin"=>θ0, "s,u"=>0.2, "s,f"=>1.0), ["s,ldlin","s,u"],
            Dict("s,ldlin"=>0.6θ0,"s,u"=>0.0), Dict("s,ldlin"=>1.5θ0,"s,u"=>1.0)
    elseif law == "quadratic"
        Dict{String,Any}("s,ldquad"=>θ0, "s,u"=>0.2, "s,w"=>0.0, "s,f"=>1.0),
            ["s,ldquad","s,u","s,w"],
            Dict("s,ldquad"=>0.6θ0,"s,u"=>-1.0,"s,w"=>-1.0),
            Dict("s,ldquad"=>1.5θ0,"s,u"=> 1.0,"s,w"=> 1.0)
    elseif law == "power law"
        Dict{String,Any}("s,ldpow"=>θ0, "s,alpha"=>0.15, "s,f"=>1.0), ["s,ldpow","s,alpha"],
            Dict("s,ldpow"=>0.6θ0,"s,alpha"=>0.0), Dict("s,ldpow"=>1.5θ0,"s,alpha"=>0.6)
    else
        error("unknown law $law")
    end
end

const LAWS = ["uniform", "linear (u fixed)", "linear (u free)", "quadratic", "power law"]
const STARS = [("AlphaCenA", 8.4), ("AlphaCenB", 5.9)]

# Every fitting algorithm OITOOLS exposes. `x_opt` and `chi2` are common to all;
# `sigma` is whatever independent uncertainty each one can produce (NLopt: none).
function run_all_fitters(md, free, dat, lb, ub; run_nested=true)
    out = Pair{String,Any}[]

    r = fit_model(md, free, dat; lb, ub, weights=V2ONLY)
    push!(out, "fit_model :LD_LBFGS" => (x=r.x_opt, σ=fill(NaN,length(free)), chi2=r.chi2, chi2r=r.chi2r))

    r = fit_model(md, free, dat; lb, ub, weights=V2ONLY, method=:LN_NELDERMEAD, maxeval=20000)
    push!(out, "fit_model :NELDERMEAD" => (x=r.x_opt, σ=fill(NaN,length(free)), chi2=r.chi2, chi2r=r.chi2r))

    r = fit_model(md, free, dat; lb, ub, weights=V2ONLY, method=:LN_COBYLA, maxeval=20000)
    push!(out, "fit_model :LN_COBYLA" => (x=r.x_opt, σ=fill(NaN,length(free)), chi2=r.chi2, chi2r=r.chi2r))

    r = fit_model_lsqfit(md, free, dat; lb, ub, weights=V2ONLY)
    push!(out, "fit_model_lsqfit" => (x=r.x_opt, σ=r.stderror, chi2=r.chi2, chi2r=r.chi2r))

    # Both nested samplers, whichever are loaded. Neither is a dependency of OITOOLS, so this
    # asks rather than assumes -- and running both is the point: two independent samplers on
    # one χ² surface either agree, or one of them has not converged.
    if run_nested
        for backend in (:nestedsamplers, :ultranest)
            backend in OITOOLS.NESTED_BACKENDS_LOADED || continue
            r = fit_model_nested(md, free, dat; backend, lb, ub, weights=V2ONLY,
                                 verb=false, cornerplot=false)
            push!(out, "fit_model_nested :$backend" =>
                       (x=r.x_opt, σ=vec(std(r.posterior, dims=1)),
                        chi2=r.chi2, chi2r=r.chi2r))
        end
    end
    return out
end

println("\n\n" * "#"^78)
println("#  PART II — validation against Kervella et al. 2017, A&A 597, A137, Table 3")
println("#"^78)
@printf("""
#
#  Diameters below are printed BOTH as fitted (raw PIONIER wavelength scale, as
#  stored in the OIFITS) and multiplied by γ = %.5f, which is the form comparable
#  with the paper. "pull" is (θ·γ − θ_published) / σ_stat.
#
""", GAMMA)

for (star, θ0) in STARS
    dat = readoifits(joinpath(DATADIR, "$star.oifits"); warn=false, verbose=false)[1,1]
    @printf("\n%s\n  %s   (N_V² = %d)\n%s\n", "="^78, star, dat.nv2, "="^78)

    for law in LAWS
        md, free, lb, ub = ld_case(law, θ0, star)
        ref = PUBLISHED[star][law]
        # A ±NaN in a results table reads as a broken fit; here it only ever means "this
        # number does not exist" — the optimiser reports no uncertainty, or the paper held
        # the parameter fixed, or it never fitted this law to this star.
        fmt_pm(v, σ) = isnan(σ) ? @sprintf("%.4f", v) : @sprintf("%.4f±%.4f", v, σ)

        @printf("\n  ── %s ──  free: %s\n", law, join(free, ", "))
        @printf("    %-22s %10s %10s %8s   %-34s %7s\n",
                "algorithm", "θ_fit", "θ×γ", "pull", "LD parameters", "χ²r")

        # LD coefficients that were held fixed rather than fitted, so the "linear (u fixed)"
        # rows say what u they were held at instead of showing an empty column.
        heldstr = join([@sprintf("%s=%.4f (fixed)", split(k,",")[2], md[k])
                        for k in ("s,u","s,w","s,alpha") if haskey(md,k) && !(k in free)], "  ")

        for (name, r) in run_all_fitters(md, free, dat, lb, ub)
            pullstr = isnan(ref.θ[1]) ? "     —  " :
                      @sprintf("%+7.1fσ", (GAMMA*r.x[1] - ref.θ[1]) / ref.θ[2])
            parstr = join(vcat([string(split(free[i],",")[2], "=", fmt_pm(r.x[i], r.σ[i]))
                                for i in 2:length(free)],
                               isempty(heldstr) ? String[] : [heldstr]), "  ")
            @printf("    %-22s %10.4f %10.4f %8s   %-34s %7.3f\n",
                    name, r.x[1], GAMMA*r.x[1], pullstr, parstr, r.chi2r)
        end

        println("    " * "-"^84)
        if isnan(ref.θ[1])
            @printf("    %-22s %10s %10s %8s   %-34s %7s\n",
                    "Kervella+2017", "—", "—", "", "not fitted for this star", "—")
        else
            refpar = join([string(k, "=", fmt_pm(v[1], v[2])) *
                           (isnan(v[2]) ? " (fixed)" : "") for (k,v) in pairs(ref.pars)], "  ")
            @printf("    %-22s %10s %10.4f %8s   %-34s %7.2f\n",
                    "Kervella+2017", "—", ref.θ[1], "", refpar, ref.chi2r)
            @printf("    %-22s %10s %10s %8s   σ_stat = ±%.4f,  σ_syst = ±%.3f\n",
                    "", "", "", "", ref.θ[2], ref.θ[3])
        end
    end
end

# ===========================================================================
#  Independent uncertainties: bootstrap resampling on the headline fit
# ===========================================================================
#
# The paper quotes σ_stat = ±0.006 mas on θ and ±0.0050 on α for α Cen A.
# `fit_model_lsqfit` gets its σ from the Jacobian and `fit_model_nested` from the
# posterior; `bootstrap_fit` resamples the data instead, so it is the one estimate that
# does not assume the χ² surface is quadratic near the minimum.

println("\n" * "="^78)
println("  Bootstrap uncertainties — α Cen A, power law")
println("="^78)

md_b, free_b, lb_b, ub_b = ld_case("power law", 8.4, "AlphaCenA")
boot = bootstrap_fit(md_b, free_b, data; nboot=500, weights=V2ONLY,
                     lb=lb_b, ub=ub_b, seed=1234, verb=false)
for (i, p) in enumerate(free_b)
    @printf("  %-12s = %.4f ± %.4f  (bootstrap: -%.4f / +%.4f, %d resamples)\n",
            split(p,",")[2], boot.x_opt[i], boot.sigma[i],
            boot.sigma_minus[i], boot.sigma_plus[i], 500)
end
@printf("  Kervella+2017 σ_stat:  θ ± 0.006 mas,  α ± 0.0050\n")

# ===========================================================================
#  Corner plot for the headline fit
# ===========================================================================

println("\n" * "="^78)
println("  Nested sampling posterior — α Cen A, power law  ($(nested_backend()))")
println("="^78)

result_un = fit_model_nested(md_b, free_b, data;
    lb=lb_b, ub=ub_b, weights=V2ONLY,
    verb=false, cornerplot=true)

println(result_un)
@printf("  θ_LD = %.4f (×γ = %.4f),  α = %.4f,  log(Z) = %.2f ± %.2f\n",
        result_un.x_opt[1], GAMMA*result_un.x_opt[1], result_un.x_opt[2],
        result_un.logz, result_un.logzerr)

println("\nDone.")
