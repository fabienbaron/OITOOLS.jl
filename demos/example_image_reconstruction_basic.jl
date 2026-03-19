#
# Very Basic Image reconstruction code
#
using OITOOLS
oifitsfile = "./data/2004-data1.oifits"
data = readoifits(oifitsfile, warn=false)[1,1]

display(data)
plot_obs(data)
plot_obs(data, color="wav")

# ── Fourier transform setup ──────────────────────────────────────────────────
pixsize = 0.2  # size of a pixel in milliarcseconds
nx = 64        # width of image (number of pixels)
ft = setup_nfft(data, nx, pixsize)
# One could also use the DFT instead of NFFT, but NFFT is much faster for larger images
# dft = setup_dft(data, nx, pixsize)

# ── Starting image ───────────────────────────────────────────────────────────
x_start = gaussian2d(nx, nx, nx/6)
image_to_chi2(x_start, ft, data; verb=true)

# ── Reconstruction ───────────────────────────────────────────────────────────
regularizers = [["centering", 1e4], ["l1l2", 7e7, 1e-3]]
x = reconstruct(x_start, data, ft; regularizers, verb=true, maxiter=500)
imdisp(x; pixsize)

# ── Check the fit ────────────────────────────────────────────────────────────
image_to_chi2(x, ft, data; verb=true)
plot_residuals(x, ft, data)

# Uncomment if you want to write the result
# writefits(x, "reconstruction.fits")
