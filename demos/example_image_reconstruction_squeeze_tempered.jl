#
# SQUEEZE with parallel tempering, on Pigeons.jl.
# Install at demos/example_image_reconstruction_squeeze_tempered.jl
#
# `reconstruct_squeeze` anneals: it is an optimizer that floors at T = 1 and then
# samples.  This runs a full temperature ladder instead, so you get a defensible
# posterior mean plus a log normalising constant (logZ) for model comparison.
#
# Pigeons is a WEAK dependency of OITOOLS — install it yourself:
#     julia> using Pkg; Pkg.add("Pigeons")
# Loading it enables `reconstruct_squeeze_tempered`.
#
# Run with several threads:  julia -t 4 <this file>
#
using OITOOLS
using Pigeons              # <- this is what enables the method

oifitsfile = joinpath(@__DIR__, "data", "2004-data1.oifits")
nx, pixsize = 64, 0.2

data = readoifits(oifitsfile; T=Float64)
ft   = setup_ft(data, nx, pixsize; mode="dft")

# ── Annealing, for comparison ────────────────────────────────────────────────
ann, da = reconstruct_squeeze(data, ft; niter=400, nchains=4, seed=1, verb=false)
println("annealed : chi2r = ", round(minimum(da.chi2r_mean), digits=3))

# ── Parallel tempering ───────────────────────────────────────────────────────
# Budget chains by the global communication barrier Λ: round trips between β=0
# and β=1 take about Λ² scans, and you want n_chains of order Λ.  Λ is a property
# of the path, but it is ESTIMATED from swap rates, so too few chains report it
# too low.  On this dataset the true value is ~27, so 8 chains would mislead you.
img, diag = reconstruct_squeeze_tempered(data, ft;
                                         n_rounds = 11,     # round i does 2^i scans
                                         n_chains = 30,
                                         seed     = 1,
                                         pixsize  = pixsize,
                                         outfile  = "squeeze_2004_tempered.fits",
                                         verb     = true)

println("\ntempered : chi2r = ", round(diag.chi2r_mean, digits=3))
println("logZ             = ", round(diag.logZ, digits=2))
println("global barrier Λ = ", round(diag.global_barrier, digits=2),
        "   (round trips ~ n_scans / Λ²)")
println("beta=1 samples   = ", diag.nsamples, "   (final round only; earlier rounds are tuning)")
println("beta schedule    = ", round.(diag.schedule, digits=4))

# If logZ moves materially when you add chains or rounds, it has not converged.

# ── Regularizers and masks work exactly as in the annealed version ───────────
c = (nx + 1) / 2
mask = [sqrt((i-c)^2 + (j-c)^2) <= 25 ? 1.0 : 0.0 for i in 1:nx, j in 1:nx]
masked, dm = reconstruct_squeeze_tempered(data, ft; n_rounds=9, n_chains=16,
                                          prior_image=mask, seed=1, verb=false)
println("\ntempered + mask  : chi2r = ", round(dm.chi2r_mean, digits=3),
        ", flux outside mask = ", sum(masked[mask .== 0]))

# ── SPARCO under tempering ───────────────────────────────────────────────────
# This is the principled fix for a shape parameter that annealing cannot move:
# at low beta the likelihood is weak, the image is free, and the parameter can
# travel; swaps then carry that up to beta = 1.  The proposal here is a symmetric
# Gaussian with a frozen stepsize (reversible), not the annealing lattice, so
# `diag.params` is a genuine posterior mean.
#
#   v1295 = readoifits(joinpath(@__DIR__, "data", "2019_v1295Aql.WL_SMOOTH.A.oifits");
#                      T=Float64, filter_bad_data=true)
#   ftv   = setup_ft(v1295, 32, 0.25; mode="dft")
#   model = SqueezeSparco(f_star=0.30, env_indx=0.0, lambda0=1.6e-6,
#                         free=(:f_star, :env_indx), stepsize=[0.02,0,0.05,0,0,0])
#   _, dv = reconstruct_squeeze_tempered(v1295, ftv; n_rounds=9, n_chains=16,
#                                        model=model, weights=[1.0,0.0,1.0])
#   println("posterior mean f_star = ", dv.params[1])

# ── Live monitoring (optional) ───────────────────────────────────────────────
#   reconstruct_squeeze_tempered(data, ft; n_rounds=10, n_chains=16,
#                                pixsize=pixsize, monitor=50)

imdisp_multi(cat(ann, img, dims=3);
             labels = ["annealed (optimizer)", "tempered (posterior mean)"],
             pixsize = pixsize, use_colorbar = true,
             figtitle = "SQUEEZE — annealing vs parallel tempering")
