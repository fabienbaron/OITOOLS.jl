# compare_errors.jl — put ASPRO 2's error bars beside `predict_errors`, channel by channel.
#
# The other half of `tools/AsproNoise.java`: that one writes an OIFITS carrying ASPRO's own
# uncertainties, this one reads it back and compares against our noise model for the same
# instrument, magnitudes and elevation.
#
#   xvfb-run -a java -cp .:Aspro2-26.09.jar AsproNoise --config "CHARA 2026A" \
#       --instrument MIRCX-MYSTIC --stations "E1 W2 W1 S2 S1 E2" --mode Low_H \
#       --target Vega --ra 18:36:56.336 --dec +38:47:01.28 --magV 0.03 --magH 0.0 \
#       --out vega.oifits
#   julia --project=. tools/compare_errors.jl vega.oifits 0.03 0.0 85 MIRCX MIRCX_LOWH
#
# Arguments: <aspro.oifits> [magV] [mag in the instrument band] [elevation deg] [combiner] [mode]
#
# Four notes on reading the output.
#
# ASPRO reports the same σ at every point of a bright target because both σ(V²) and σ(CP)
# saturate at the calibration bias (`instrumentVisibilityBias`, `instrumentPhaseBias`); only
# faint targets exercise the statistical part. And `predict_errors` takes ONE elevation, while
# the file spans a night, so compare against the BEST (smallest) σ in each channel, which is
# what the table below prints.
#
# The OI_VIS block compares ABSOLUTE σ|V| and σφ: no CHARA instrument asks ASPRO for
# differential quantities (`<oiVisAmpDiff>`/`<oiVisPhiDiff>` are false for SPICA and absent for
# MIRC-X and MYSTIC, so both default to false), and AMPTYP/PHITYP in its output say so. Those
# two columns are a tighter test than V², since σ|V| and σφ come straight from σ_c unsquared.
#
# The ratios are NOT expected to reach 1.00 above R ≈ 150. Two conventions differ, both recorded
# in TODO.md: ASPRO convolves a SkyCalc telluric table per channel where `atm_transmission` is
# flat 1.0, and ASPRO takes one zero point per photometric band where `zero_point_flux`
# interpolates per wavelength. On MIRC-X High_H the first shows up as a handful of channels
# around 1.40 µm and the second as a smooth drift to ~1.5 at the red band edge.

using OITOOLS, Printf, Statistics

if isempty(ARGS)
    println("usage: julia --project=. tools/compare_errors.jl <aspro.oifits> [magV] [mag] [elev] [combiner] [mode]")
    exit(1)
end

path  = ARGS[1]
magV  = length(ARGS) > 1 ? parse(Float64, ARGS[2]) : 0.03
magI  = length(ARGS) > 2 ? parse(Float64, ARGS[3]) : 0.0
elev  = length(ARGS) > 3 ? parse(Float64, ARGS[4]) : 85.0
cname = length(ARGS) > 4 ? ARGS[5] : "MIRCX"
wname = length(ARGS) > 5 ? ARGS[6] : "MIRCX_LOWH"

dm = readoifits(path; warn = false, verbose = false, T = Float64)
d  = dm[1, 1]
@printf("ASPRO file: %s — nv2 %d  nt3amp %d  nt3phi %d\n", basename(path), d.nv2, d.nt3amp, d.nt3phi)

# ── ASPRO's error bars, per spectral channel ──────────────────────────────────
λv2 = Float64.(d.v2_lam);  σv2 = Float64.(d.v2_err);  v2 = Float64.(d.v2)
λt3 = Float64.(d.t3_lam);  σcp = Float64.(d.t3phi_err)
chans = sort(unique(round.(λv2; sigdigits = 7)))
aspro = Dict{Float64,NTuple{2,Float64}}()
@printf("\n%-10s %-11s %-11s %-11s %-11s %-11s\n",
        "λ (µm)", "σV² med", "σV² best", "V² med", "σCP med", "σCP best")
for λ in chans
    m  = abs.(λv2 .- λ) .< 1e-11
    mt = abs.(λt3 .- λ) .< 1e-11
    sv, sc = σv2[m], σcp[mt]
    aspro[λ] = (minimum(sv), isempty(sc) ? NaN : minimum(sc))
    @printf("%-10.4f %-11.4g %-11.4g %-11.4g %-11.4g %-11.4g\n",
            λ * 1e6, median(sv), minimum(sv), median(v2[m]),
            isempty(sc) ? NaN : median(sc), isempty(sc) ? NaN : minimum(sc))
end

# ── our noise model, same instrument, same magnitudes ─────────────────────────
fac  = read_facility_file("CHARA")
comb = read_comb_file(cname)
wav  = read_wave_file(wname)
# Every band from R redward gets the instrument-band magnitude, so `resolve_magnitudes`
# interpolates V→R and then holds flat: a SPICA mode at 0.61-0.90 µm is then driven by the
# magnitude asked for rather than by an extrapolation off V.
mags = Dict("V" => magV, "R" => magI, "I" => magI, "J" => magI, "H" => magI, "K" => magI)
pe   = predict_errors(fac, comb, wav; mag = mags, visamp = 1.0, elevation_deg = elev)

@printf("\nours: %s / %s, V = %.2f, band mag = %.2f, elevation %.0f°\n",
        cname, wname, magV, magI, elev)
@printf("\n%-10s %-11s %-11s %-11s %-11s %-9s %-9s\n",
        "λ (µm)", "ours σV²", "ASPRO σV²", "ours σCP", "ASPRO σCP", "V² ratio", "CP ratio")
for (i, λ) in enumerate(pe.λ)
    aλ = chans[argmin(abs.(chans .- λ))]
    (asv, asc) = aspro[aλ]
    @printf("%-10.4f %-11.4g %-11.4g %-11.4g %-11.4g %-9.2f %-9.2f\n",
            λ * 1e6, pe.sigma_v2[i], asv, pe.sigma_cp[i], asc,
            pe.sigma_v2[i] / asv, pe.sigma_cp[i] / asc)
end
@printf("\nours: strehl %.3f, dit %.4g s, nframes %d, nphot/tel/DIT %.4g\n",
        pe.strehl[1], pe.dit[1], pe.nframes[1], pe.nphot[1])

# ── differential quantities, if the file carries OI_VIS ───────────────────────
if d.nvisamp > 0 || d.nvisphi > 0
    λva = Float64.(d.vis_lam); σva = Float64.(d.visamp_err); σvp = Float64.(d.visphi_err)
    @printf("\nOI_VIS present: amptyp %s, phityp %s  (nvisamp %d, nvisphi %d)\n",
            repr(d.amptyp), repr(d.phityp), d.nvisamp, d.nvisphi)
    @printf("%-10s %-12s %-12s %-12s %-12s %-9s %-9s\n",
            "λ (µm)", "ours σ|V|", "ASPRO σ|V|", "ours σφ(°)", "ASPRO σφ(°)", "|V| ratio", "φ ratio")
    for (i, λ) in enumerate(pe.λ)
        m = abs.(λva .- λ) .< 1e-11
        count(m) == 0 && continue
        av, ap = minimum(σva[m]), minimum(σvp[m])
        @printf("%-10.4f %-12.4g %-12.4g %-12.4g %-12.4g %-9.2f %-9.2f\n",
                λ * 1e6, pe.sigma_visamp[i], av, pe.sigma_visphi[i], ap,
                pe.sigma_visamp[i] / av, pe.sigma_visphi[i] / ap)
    end
end
