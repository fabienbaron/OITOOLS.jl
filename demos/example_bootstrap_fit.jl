#
# UNCERTAINTY ESTIMATION BY RESAMPLING
#
# Compares, on the same data and the same model:
#   - the analytic (Levenberg-Marquardt) covariance,
#   - a parametric Monte Carlo (noise drawn from the error bars),
#   - the nonparametric block bootstrap, at three granularities and under
#     three resampling schemes (including PMOIRED's).
#
# Run with several threads to speed up the replicates:  julia -t auto
#
using OITOOLS
using Statistics, Printf

# VLTI/PIONIER example
oifitsfile = "./data/AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1];

# Model: linear limb-darkened disk
model = Dict{String,Any}(
    "star,ldlin" => 8.0,
    "star,u"     => 0.3,
    "star,f"     => 1.0,
)
free_params = ["star,ldlin", "star,u"]
lb = Dict("star,ldlin" => 5.0, "star,u" => 0.0)
ub = Dict("star,ldlin" => 12.0, "star,u" => 1.0)
weights = [1.0, 0.0, 0.0]      # V² only

# ── Analytic errors ─────────────────────────────────────────────────────────
# LsqFit rescales the covariance by the reduced chi2, i.e. it inflates the
# error bars until chi2r = 1.
fit = fit_model_lsqfit(model, free_params, data; lb=lb, ub=ub, weights=weights)
println(fit)

# ── How is the data split into blocks? ──────────────────────────────────────
println(data_blocks(data; granularity=:config))
println(data_blocks(data; granularity=:epoch))

# ── Block bootstrap ─────────────────────────────────────────────────────────
# The values are left untouched; what is resampled is which (MJD, baseline)
# spectral vectors enter the fit.  This responds to correlated calibration
# errors and to underestimated error bars.
boot = bootstrap_fit(model, free_params, data;
    lb=lb, ub=ub, weights=weights, nboot=500, seed=42)
println(boot)

# ── Everything side by side ─────────────────────────────────────────────────
nboot = 500
@printf("\n%-32s  %-24s %-24s\n", "method", free_params[1], free_params[2])
@printf("%-32s  %8.5f +/- %-11.5f %8.5f +/- %-11.5f\n", "analytic (LM, chi2r-scaled)",
        fit.x_opt[1], fit.stderror[1], fit.x_opt[2], fit.stderror[2])

for (label, kw) in [
        ("block bootstrap (config)",    (mode=:replacement, granularity=:config)),
        ("half-sampling  (config)",     (mode=:halfsample,  granularity=:config)),
        ("weighted       (config)",     (mode=:weights,     granularity=:config)),
        ("PMOIRED scheme [biased!]",    (mode=:pmoired,     granularity=:config)),
        ("block bootstrap (epoch)",     (mode=:replacement, granularity=:epoch)),
        ("point bootstrap (no blocks)", (mode=:replacement, granularity=:point)),
    ]
    b = bootstrap_fit(model, free_params, data;
        lb=lb, ub=ub, weights=weights, nboot=nboot, seed=1234, verb=false, kw...)
    s = [std(b.samples[b.mask, j]) for j in eachindex(free_params)]
    @printf("%-32s  %8.5f +/- %-11.5f %8.5f +/- %-11.5f\n",
            label, b.median[1], s[1], b.median[2], s[2])
end

# Expected pattern on real (correlated) data:
#   - the analytic errors are the most optimistic: they only know the quoted
#     error bars, rescaled to chi2r = 1;
#   - the per-point bootstrap is barely better: dropping single channels does
#     not probe the correlated part of the error budget;
#   - the block bootstrap is 2-4x larger, and the epoch bootstrap larger still;
#   - PMOIRED's scheme sits ~1/sqrt(2) below the block bootstrap by
#     construction (block multiplicity variance 1/2 instead of 1).  It is here
#     to cross-check PMOIRED results, not to quote error bars from.

# ── Derived quantities: use the samples, not the error bars ─────────────────
# (propagates the full covariance; here, the LD-corrected "Rosseland" size)
θ_ld = boot.samples[boot.mask, 1]
u    = boot.samples[boot.mask, 2]
θ_ud = θ_ld .* sqrt.((1 .- 7 .* u ./ 15) ./ (1 .- u ./ 3))
@printf("\nequivalent UD diameter = %.4f +%.4f -%.4f mas\n",
        median(θ_ud), quantile(θ_ud, 0.84) - median(θ_ud),
        median(θ_ud) - quantile(θ_ud, 0.16))

# ── Corner plot of the bootstrap samples ────────────────────────────────────
using PythonPlot
figure("Bootstrap", figsize=(6,6)); clf()
subplot(2,2,1); PythonPlot.hist(boot.samples[boot.mask,1], 40); xlabel(free_params[1])
subplot(2,2,4); PythonPlot.hist(boot.samples[boot.mask,2], 40); xlabel(free_params[2])
subplot(2,2,3); plot(boot.samples[boot.mask,1], boot.samples[boot.mask,2], ".", ms=2)
xlabel(free_params[1]); ylabel(free_params[2])
tight_layout()
