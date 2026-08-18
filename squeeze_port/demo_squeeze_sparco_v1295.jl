#
# SQUEEZE + SPARCO on V1295 Aql (HD 190073) — SQUEEZE's own SPARCO example.
# Install at demos/example_image_reconstruction_squeeze_sparco_v1295.jl
#
# Physical model (Kluska et al. 2012, plus an over-resolved background):
#   - Star:        uniform disc (0 = point source), flux fraction f_star ∝ (λ/λ0)^-4
#   - Background:  over-resolved (V = 0), flux fraction f_bg ∝ (λ/λ0)^bg_indx
#   - Environment: the reconstructed GREY image, flux 1 - f_star - f_bg,
#                  scaled by (λ/λ0)^env_indx
#
# All the chromatism lives in the model; the image is grey.  That is why a single
# OIdata spanning several wavelengths is enough — the chromatism enters per uv
# point through `uv_lam`, not through a polychromatic image cube.
#
# It is also what makes f_star identifiable: on MONOCHROMATIC data a central point
# star is degenerate with a central component of the image, and no amount of
# sampling will separate them.
#
# T3amp is disabled (weights = [1, 0, 1]) — this dataset's T3amp is noisy.
#
# Run with several threads:  julia -t 4 <this file>
#
using OITOOLS

oifitsfile = joinpath(@__DIR__, "data", "2019_v1295Aql.WL_SMOOTH.A.oifits")
nx, pixsize = 64, 0.125
lambda0 = 1.6e-6
weights = [1.0, 0.0, 1.0]

data = readoifits(oifitsfile; T=Float64, filter_bad_data=true)
ft   = setup_ft(data, nx, pixsize; mode="dft")

# Heads-up on memory: the sampler needs the dense DFT matrix, which is
# nuv * nx^2 * 2 * sizeof(T) bytes.  This dataset has nuv ~ 4900, so nx=64 in
# Float64 is ~320 MB.  Read with T=Float32 to halve that.

nchains = max(Threads.nthreads(), 2)

# ── Image only, no model: this is what SPARCO is for ─────────────────────────
_, d0 = reconstruct_squeeze(data, ft; niter = 300, nchains = nchains,
                            weights = weights, regularizers = [["tv", 100.0]],
                            seed = 1, verb = false)
println("image only, no SPARCO : chi2r = ", round(minimum(d0.chi2r_mean), digits=3))

# ── With SPARCO, f_star fixed at the published value ─────────────────────────
model_fixed = SqueezeSparco(f_star = 0.48, ud = 0.0, env_indx = -0.1,
                            lambda0 = lambda0, free = ())
_, d1 = reconstruct_squeeze(data, ft; niter = 300, nchains = nchains,
                            weights = weights, regularizers = [["tv", 100.0]],
                            model = model_fixed, seed = 1, verb = false)
println("SPARCO, f_star = 0.48 : chi2r = ", round(minimum(d1.chi2r_mean), digits=3))

# ── With SPARCO, fitting f_star and the environment index ────────────────────
# `free` names the parameters the sampler may vary; `stepsize` is the two-point
# lattice half-width, adapted toward a 0.3 acceptance rate during the run.
model_free = SqueezeSparco(f_star = 0.30, ud = 0.0, env_indx = 0.0,
                           lambda0 = lambda0,
                           free = (:f_star, :env_indx),
                           stepsize = [0.01, 0.0, 0.05, 0.0, 0.0, 0.0])

img, diag = reconstruct_squeeze(data, ft;
                                niter        = 400,
                                nchains      = nchains,
                                weights      = weights,
                                regularizers = [["tv", 100.0]],
                                model        = model_free,
                                seed         = 1,
                                pixsize      = pixsize,
                                outfile      = "squeeze_v1295_sparco.fits",
                                print_every  = 50,   # C prints every sweep
                                verb         = true)

b = diag.best_chain
println("\nSPARCO, f_star free   : chi2r = ", round(minimum(diag.chi2r_mean), digits=3))
for (j, nm) in enumerate(diag.param_names)
    println("  $nm = ", round(diag.params[b][j], sigdigits=5),
            j in model_free.free ? "" : "   (fixed)")
end
println("spread across chains  : f_star = ",
        round.([p[1] for p in diag.params], digits=4))

# ⚠️ The parametric proposal is C's two-point lattice, whose spacing is rescaled
# on every proposal.  That makes the move irreversible, so the spread above is a
# spread of optimizer endpoints, NOT a posterior uncertainty.  It is also why a
# parameter started far from its optimum may never travel: the adaptation shrinks
# the step toward zero.  Start near a plausible value (e.g. from `fit_model`),
# widen `stepsize`, or use reconstruct_squeeze_tempered — which fixes this.

imdisp(img, pixsize = pixsize, use_colorbar = true,
       figtitle = "V1295 Aql — SQUEEZE + SPARCO")
