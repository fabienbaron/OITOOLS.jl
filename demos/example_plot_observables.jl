using OITOOLS

# Read BC2026 = Fake MATISSE LM-band data
oifitsfile = "./data/BC2026/OBJECT1_LM.oifits"
data = readoifits(oifitsfile)

# ── UV coverage ──────────────────────────────────────────────────────────────
uvplot(data, color="baseline")
uvplot(data, color="wav")
uvplot(data, color="mjd")

# ── Squared visibilities ─────────────────────────────────────────────────────
plot_v2(data, color="baseline")
plot_v2(data, color="wav")
plot_v2(data, color="mjd")
plot_v2(data, logplot=true)

# ── Closure phases ───────────────────────────────────────────────────────────
plot_t3phi(data)                        # geometric mean baseline (default)
plot_t3phi(data, t3base="max")          # longest baseline
plot_t3phi(data, color="wav")
plot_t3phi(data, color="mjd")

# ── Triple amplitudes ────────────────────────────────────────────────────────
plot_t3amp(data)
plot_t3amp(data, logplot=true, color="wav")
plot_t3amp(data, color="mjd")

# ── Visibility amplitudes ────────────────────────────────────────────────────
plot_visamp(data)
plot_visamp(data, logplot=true, color="wav")
plot_visamp(data, color="mjd")

# ── Visibility phases ────────────────────────────────────────────────────────
plot_visphi(data)                       # auto-detects absolute vs differential
plot_visphi(data, color="wav")
plot_visphi(data, color="mjd")

# ── Flux (if present) ────────────────────────────────────────────────────────
plot_flux(data)                     # colour by wavelength
plot_flux(data, color="mjd")

# ── Multi-panel overview ───────────────────────────────────────────────────
plot_multi(data)                                             # default: V², T3φ, T3amp
plot_multi(data, obs=["V2", "T3PHI", "VISAMP"])    # custom selection
plot_multi(data, color="wav")                                # colour by wavelength

# ── Compact multi-panel ───────────────────────────────────────────────────
set_oiplot_defaults(compact=true)
plot_multi(data)
plot_multi(data, obs=["V2", "T3PHI", "T3AMP", "VISAMP", "VISPHI"])
set_oiplot_defaults(compact=false)
