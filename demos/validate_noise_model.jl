# validate_noise_model.jl
#
# Prints the two tables you need to judge whether the CHARA noise model is behaving:
#
#   1. Strehl vs seeing and guide magnitude, against ASPRO 2's published CHARA curves.
#   2. sigma(V2) and sigma(CP) vs target magnitude for MIRC-X, MYSTIC and SPICA.
#
# Both are regression checks against ASPRO 2. The instrument parameters in
# src/configs/{MIRCX,MYSTIC,SPICA}.toml are transcribed from ASPRO's aspro-conf CHARA.xml,
# and the Strehl model is a port of jmal's Band.strehl, so these tables should track ASPRO
# rather than needing to be tuned.

using OITOOLS, Printf, Statistics

# ── 1. Strehl vs ASPRO ────────────────────────────────────────────────────────

ao      = AOConfig()                       # ASPRO AO_CHARA defaults
seeings = [0.60, 0.70, 1.00, 1.15, 1.40, 1.80]
t0s     = [5.2, 4.4, 3.2, 2.2, 1.6, 1.0] .* 1e-3

# (bright, faint) plateaus read off ASPRO's reference plots in aspro/doc/strehl/.
# The three R faint values marked (*) sit at the bottom of a plot axis and are unreliable
# readings, not model failures: the model floor there is the seeing-limited halo.
const ASPRO_REF = Dict(
    ("H", 1.57e-6) => [(0.770,0.330),(0.710,0.270),(0.555,0.150),(0.395,0.120),(0.235,0.083),(0.078,0.053)],
    ("K", 2.18e-6) => [(0.878,0.520),(0.828,0.443),(0.683,0.280),(0.518,0.228),(0.337,0.167),(0.140,0.107)],
    ("R", 0.75e-6) => [(0.495,0.080),(0.432,0.060),(0.279,0.030),(0.170,0.012),(0.080,0.008),(0.018,0.005)],
)

println("="^78)
println("1. Strehl ratio vs ASPRO 2 (CHARA, D = 1 m, AO_CHARA)")
println("="^78)
@printf("%-5s %-7s %9s %8s %8s %9s %8s %8s\n",
        "band","seeing","bright","ASPRO","err%","faint","ASPRO","err%")
errs_hk = Float64[]
for ((bn, lam), ref) in sort(collect(ASPRO_REF), by=x->x[1][2])
    for i in eachindex(seeings)
        sb = strehl_ratio(lam, 1.0, seeings[i], t0s[i],  2.0, ao)
        sf = strehl_ratio(lam, 1.0, seeings[i], t0s[i], 16.0, ao)
        rb, rf = ref[i]
        eb, ef = 100*(sb-rb)/rb, 100*(sf-rf)/rf
        bn != "R" && append!(errs_hk, (abs(eb), abs(ef)))
        @printf("%-5s %-7.2f %9.4f %8.3f %8.1f %9.4f %8.3f %8.1f\n",
                bn, seeings[i], sb, rb, eb, sf, rf, ef)
    end
end
@printf("\nH/K agreement: median %.2f%%, worst %.2f%%\n", median(errs_hk), maximum(errs_hk))

print("\nStrehl vs guide magnitude (H, seeing 1.0\"): ")
for m in [2.0, 8.5, 10.5, 11.5, 12.5, 16.0]
    @printf("V=%.1f:%.3f  ", m, strehl_ratio(1.57e-6, 1.0, 1.0, 3.2e-3, m, ao))
end
println("\n   (ASPRO: flat ~0.55 to V=8.5, knee V=10.5-11.5, floor 0.15 by V=12.5)")

print("Strehl vs elevation   (H, seeing 1.0\", V=6): ")
for e in [90.0, 45.0, 30.0, 20.0, 10.0]
    @printf("%.0fdeg:%.3f  ", e, strehl_ratio(1.57e-6, 1.0, 1.0, 3.2e-3, 6.0, ao; elevation_deg=e))
end
println()

# ── 2. Sensitivity per instrument ─────────────────────────────────────────────

println("\n" * "="^78)
println("2. Predicted uncertainties vs magnitude (zenith, seeing 1.0\", per-instrument DIT and total time)")
println("   V2 = 1 (unresolved). Instrument parameters are transcribed from ASPRO CHARA.xml.")
println("="^78)

facility = read_facility_file("CHARA")
cases = [("MIRC-X H", "MIRCX",  "MIRCX_LOWH"),
         ("MIRC-X J", "MIRCX",  "MIRCX_LOWJ"),
         ("MYSTIC K", "MYSTIC", "MYSTIC_LOWK"),
         ("SPICA  R", "SPICA",  "SPICA_LR")]

@printf("%-10s %8s %8s %8s %10s %10s %10s\n",
        "instrument","mag","Strehl","DIT(ms)","N/tel/DIT","sigma(V2)","sigma(CP)")
for (label, combname, wavename) in cases
    comb = read_comb_file(combname)
    wave = read_wave_file(wavename)
    iw   = cld(length(wave.λ), 2)                      # middle channel
    λ, δλ = wave.λ[iw], wave.δλ[iw]
    for mag in [2.0, 5.0, 7.0, 9.0, 11.0]
        mags = fill(mag, length(wave.λ))
        N, DIT, nfr = OITOOLS.photons_per_telescope(mags, wave.λ, wave.δλ, facility, comb,
                                            [90.0], mag; channel=:fringes)
        Np, _, _    = OITOOLS.photons_per_telescope(mags, wave.λ, wave.δλ, facility, comb,
                                            [90.0], mag; channel=:photometry)
        v2st = hcat([[i,j] for i in 1:facility.ntel for j in i+1:facility.ntel]...)
        sq, vc, vk, vp = OITOOLS.correlated_flux_coefficients(N, Np, comb, facility.ntel, v2st)
        σc = OITOOLS.complex_vis_error(1.0, sq[1,1,iw], vc[1,1,iw], vk[1,1,iw], vp[1,1,iw], nfr)
        σ2 = OITOOLS.vis2_error_from_sigma(1.0, σc)
        σ2 = sqrt(σ2^2 + (2*comb.vis_cal_err)^2)
        σcp = sqrt((180/π)^2 * 3 * (σ2/2)^2 + comb.phase_cal_err^2)
        S = comb.strehl_model == "fixed_spica" ? OITOOLS._spica_fixed_strehl(facility.seeing) :
            coupling_efficiency(λ, 1.0, facility.seeing, facility.t0, mag, facility.ao)
        @printf("%-10s %8.1f %8.3f %8.2f %10.1f %10.4f %10.2f\n",
                label, mag, S, DIT*1e3, N[1,1,iw], σ2, σcp)
    end
    println()
end
