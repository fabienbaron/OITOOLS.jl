#
# SQUEEZE-style MCMC image reconstruction on the 2004 Beauty Contest data.
# Install at demos/example_image_reconstruction_squeeze.jl
#
# Unlike the gradient reconstructions (`reconstruct`, `reconstruct_bsmem`,
# `reconstruct_bsdmm`), the image here is a bag of discrete flux quanta sampled
# by Metropolis-Hastings with simulated annealing.  Positivity and total flux are
# exact by construction, non-convex regularizers such as L0 are available, and
# the result is a posterior mean image rather than a MAP point.
#
# Run with several threads to use several chains:  julia -t 4 <this file>
#
using OITOOLS

oifitsfile = joinpath(@__DIR__, "data", "2004-data1.oifits")
nx, pixsize = 64, 0.2

# Float64 while exploring; Float32 halves the DFT matrix and is ~2x faster on the
# rank-1 update, at the cost of periodic re-syncs (see `default_resync`).
data = readoifits(oifitsfile; T=Float64)

# The sampler needs O(nuv) access to a SINGLE COLUMN of the transform per move,
# so build the operator in "dft" mode.  Passing an NFFT plan also works — nx and
# pixsize are recovered from it, as in reconstruct_bsmem — but the DFT matrix is
# then built internally anyway.
ft = setup_ft(data, nx, pixsize; mode="dft")

nchains = max(Threads.nthreads(), 2)

# ── Unregularized (centering is enabled automatically, as in the C code) ──────
# Omitting x_start gives SQUEEZE's default: a point source at the grid centre.
img, diag = reconstruct_squeeze(data, ft;
                                niter   = 400,
                                nchains = nchains,
                                seed    = 1,
                                pixsize = pixsize,          # for the FITS WCS header
                                outfile = "squeeze_2004.fits")

println("\nbest chain      : ", diag.best_chain)
println("chi2r per chain : ", round.(diag.chi2r_mean, digits=3))
println("acceptance      : ", round.(diag.acceptance, digits=3))
println("total flux      : ", sum(img))

# ── With regularizers ────────────────────────────────────────────────────────
# NOTE these are SQUEEZE's forms evaluated on the INTEGER histogram, not
# OITOOLS' differentiable regularizers from oichi2.jl — the definitions and
# hence the useful lambda ranges differ.  Order-of-magnitude starting points on
# this dataset (nx=64, nelements~1000):
#
#   l0           ~ 1-20      (counts nonzero pixels; delta per move is +-1)
#   tv           ~ 200-2000  (divided by nelements, so deltas are ~1e-3)
#   entropy      ~ 1-10      (sum lgamma(counts))
#   compactness  ~ 20-500
#   centering    ~ 1         (auto-enabled; C's DEFAULT_CENT_MULT)
#
# Push lambda too far and the prior dominates the likelihood: chi2r climbs, the
# anneal needs more sweeps and chains, and L0 in particular will shatter a smooth
# ring into isolated hot pixels long before chi2r warns you.
img_reg, diag_reg = reconstruct_squeeze(data, ft;
                                        niter   = 400,
                                        nchains = nchains,
                                        seed    = 1,
                                        regularizers = [["l0", 1.0], ["tv", 200.0]],
                                        pixsize = pixsize,
                                        outfile = "squeeze_2004_reg.fits")

println("\nregularized chi2r : ", round(minimum(diag_reg.chi2r_mean), digits=3))
println("nonzero pixels    : ", count(>(0), img_reg), " (unregularized: ",
        count(>(0), img), ")")

# ── Other starting configurations ────────────────────────────────────────────
#   reconstruct_squeeze(:random,      data, ft; ...)   # uniform scatter
#   reconstruct_squeeze(previous_img, data, ft; ...)   # digitise an image (C's -i)

# ── Prior image / mask (C's -p) ──────────────────────────────────────────────
#   c = (nx + 1) / 2
#   mask = [sqrt((i-c)^2 + (j-c)^2) <= 25 ? 1.0 : 0.0 for i in 1:nx, j in 1:nx]
#   reconstruct_squeeze(data, ft; niter=400, prior_image=mask)
# Pixels with p <= 0 get a 1e12 penalty, so no quantum can ever sit there.

# ── Live monitoring (optional) ───────────────────────────────────────────────
# Redraw the image and a chi2r/temperature trace every 50 sweeps while running:
#
#   reconstruct_squeeze(data, ft; niter=400, pixsize=pixsize, monitor=50)
#
# Costs ~40 ms per update (so ~33% on a 400-sweep run at monitor=50) and forces
# serial chains, because the Python plotting bridge is not thread-safe.
# monitor=0 (default) is free.

# ── Display ──────────────────────────────────────────────────────────────────
# imdisp/imdisp_multi use OITOOLS' standard orientation (Monnier convention,
# East left) — use them rather than calling imshow directly, or the orientation
# will not match the rest of the package.
imdisp(img, pixsize=pixsize, figtitle="SQUEEZE — unregularized", use_colorbar=true)
imdisp_multi(cat(img, img_reg, dims=3);
             labels = ["unregularized", "l0=1 + tv=200"],
             pixsize = pixsize, use_colorbar = true,
             figtitle = "SQUEEZE — 2004 data")
