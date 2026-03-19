using OITOOLS

# ── Load image and data ──────────────────────────────────────────────────────
fitsfile   = "./data/2004true.fits"
pixsize    = 0.101   # mas/pixel
x_true     = readfits(fitsfile)
nx         = size(x_true, 1)

imdisp(x_true; pixsize, tickinterval=1.0, beamsize=1.0,
       beamlocation=[0.85, 0.85], use_colorbar=true)

oifitsfile = "./data/2004-data1.oifits"
data       = readoifits(oifitsfile, warn=)[1,1]

# ── Inspect the data ─────────────────────────────────────────────────────────
display(data)
uvplot(data)
plot_obs(data)                                  # multi-panel overview
plot_obs(data; obs=["V2", "T3PHI", "T3AMP"])    # custom selection
plot_obs(data; color="wav")                     # colour by wavelength

# ── Forward model: image → observables ───────────────────────────────────────

# Setup Fourier transform (NFFT — fast)
ft = setup_nfft(data, nx, pixsize)

# Chi-squared
chi2 = image_to_chi2(x_true, ft, data; verb=true)

# All observables as a NamedTuple: (v2, t3amp, t3phi, visamp, visphi)
obs = image_to_obs(x_true, ft, data)
println("Observable fields: ", keys(obs))
println("V² points:         ", length(obs.v2))
println("T3φ points:        ", length(obs.t3phi))

# ── Residuals ────────────────────────────────────────────────────────────────

# Normalised residuals: (model - data) / error
res = image_to_residuals(x_true, ft, data)
println("V² residual std:  ", round(std(res.v2), digits=3))
println("T3φ residual std: ", round(std(res.t3phi), digits=3))

# ── Residual plots ───────────────────────────────────────────────────────────

# Individual residual plots
plot_v2_residuals(data, obs.v2)
plot_t3phi_residuals(data, obs.t3phi)
plot_t3amp_residuals(data, obs.t3amp)

# Combined residual plot — all observables in one figure
plot_residuals(data, obs)

# Convenience: pass image directly (computes obs internally)
plot_residuals(x_true, ft, data)

# ── DFT comparison ───────────────────────────────────────────────────────────

# Setup exact DFT (slower but useful for small problems)
dft = setup_dft(data, nx, pixsize)
chi2_dft = image_to_chi2(x_true, dft, data; verb=true)
obs_dft  = image_to_obs(x_true, dft, data)
plot_residuals(data, obs_dft)
