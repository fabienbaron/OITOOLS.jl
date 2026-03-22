#
# Very Basic Image reconstruction code
#
using OITOOLS
oifitsfile = "./data/2004-data1.oifits"
data = readoifits(oifitsfile)

# ── Starting image ───────────────────────────────────────────────────────────
nx = 128
pixsize = 0.101
ft = setup_ft(data, nx, pixsize)  # nx=256, pixsize=0.101 mas

x_start = gaussian2d(nx, nx, nx/6)
image_to_chi2(x_start, ft, data; verb=true)

# ── Reconstruction ───────────────────────────────────────────────────────────
regularizers = [["centering", 1e3], ["l1l2", 100.0, 1e-3]]
x = reconstruct(x_start, data, ft; regularizers, verb=true, maxiter=500)
imdisp(x; pixsize)

# ── Check the fit ────────────────────────────────────────────────────────────
image_to_chi2(x, ft, data; verb=true)
plot_residuals(x, ft, data)

# Uncomment if you want to write the result
# writefits(x, "reconstruction.fits")
