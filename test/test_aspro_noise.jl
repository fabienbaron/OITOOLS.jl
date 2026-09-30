# test_aspro_noise.jl — the CHARA noise model against ASPRO 2, channel by channel.
#
# Everything else in this suite checks that the code agrees with itself. This one checks it
# against an INDEPENDENT implementation of the same physics: JMMC's ASPRO 2, which is what
# CHARA proposals are actually written with. It exists because three real defects in the noise
# model were each found by hand, by driving ASPRO under a virtual X server and reading its
# debug log, and none of them was visible to any other test:
#
#   * sigma(CP) carried the calibration bias twice, so it was sqrt(2) too large;
#   * the read-noise term in `complex_vis_error` scaled as N^(-1/4) instead of N^(-1/2);
#   * the shipped transmissions double-counted the mean atmosphere, leaving the throughput
#     ~20% low and sigma(V2) 1.44x too large at H = 8.
#
# The last one sat in three configs while 1637 tests passed. A few seconds here replaces that.
#
# ── The reference ────────────────────────────────────────────────────────────
#
# `references/aspro_noise.csv` holds ASPRO's own error bars: the BEST (smallest) sigma in each
# spectral channel over a full night, which is the highest-elevation point and so what a
# single-elevation `predict_errors` call is comparable to. Regenerate with `tools/AsproNoise.java`:
#
#   xvfb-run -a java -cp .:Aspro2-26.09.jar AsproNoise --config "CHARA 2026A" \
#       --instrument MIRCX-MYSTIC --stations "E1 W2 W1 S2 S1 E2" --mode Low_H --target T \
#       --ra 18:36:56.336 --dec +38:47:01.28 --magV 8.5 --magH 8 --out lad_mircx_8.oifits
#
# and for SPICA add `--ft FringeTrack`, WITHOUT which ASPRO computes it with no fringe tracker —
# a state its own GUI does not offer, since SPICA sets <fringeTrackerRequired>true</...>.
# The `fringe_tracker` column records which was used; getting it wrong moves SPICA by 7x.
#
# ── Why the tolerances are what they are ─────────────────────────────────────
#
# Two modelling differences remain by choice, so the bands are not "as tight as today's code
# happens to be". Both are in TODO.md:
#
#   * tellurics — `atm_transmission` is flat 1.0 where ASPRO convolves a SkyCalc table. That is
#     a per-CHANNEL effect: at H = 8 MIRC-X sits at 0.99 on the median and 0.44 in the H2O
#     channel at 1.4152 um. Hence the assertions are on the MEDIAN over channels, with a loose
#     per-channel sanity band.
#   * the zero point — ours interpolates per wavelength, ASPRO takes one per photometric band.
#
# On the calibration floor neither matters: sigma saturates at 2*vis_cal_err and the photon
# budget cancels, so agreement there is exact and the band is 2%. That also means floor-regime
# agreement proves very little, which is why the faint rows are the ones that earn their keep.

using OITOOLS, Test, Printf, Statistics, DelimitedFiles

@testset "the CHARA noise model against ASPRO 2" begin
    ref = joinpath(@__DIR__, "references", "aspro_noise.csv")
    @test isfile(ref)
    raw, hdr = readdlm(ref, ','; header = true)
    col = Dict(strip(String(h)) => i for (i, h) in enumerate(vec(hdr)))
    fac = read_facility_file("CHARA")

    keys_seen = unique([(String(raw[r, col["combiner"]]), String(raw[r, col["wave_config"]]),
                         Int(raw[r, col["mag"]]), Bool(raw[r, col["fringe_tracker"]]))
                        for r in axes(raw, 1)])
    @test length(keys_seen) == 9          # 3 instruments x 3 magnitudes

    for (cname, wname, mag, ft) in keys_seen
        @testset "$cname $wname mag $mag$(ft ? " +FT" : "")" begin
            rows = [r for r in axes(raw, 1)
                    if String(raw[r, col["combiner"]]) == cname &&
                       Int(raw[r, col["mag"]]) == mag &&
                       Bool(raw[r, col["fringe_tracker"]]) == ft]
            λref = Float64[raw[r, col["lambda_m"]]  for r in rows]
            σv2  = Float64[raw[r, col["sigma_v2"]]  for r in rows]
            σcp  = Float64[raw[r, col["sigma_cp"]]  for r in rows]

            comb = read_comb_file(cname)
            wav  = read_wave_file(wname)
            m    = Dict("V" => mag + 0.5, "R" => float(mag), "I" => float(mag),
                        "J" => float(mag), "H" => float(mag), "K" => float(mag))
            pe = predict_errors(fac, comb, wav; mag = m, visamp = 1.0, elevation_deg = 85.0,
                                fringe_tracker = ft)

            # The reference must describe the same instrument we are about to evaluate.
            @test length(pe.λ) == length(λref)
            @test maximum(abs.(pe.λ .- λref)) < 1e-11

            rv2 = pe.sigma_v2 ./ σv2
            rcp = pe.sigma_cp ./ σcp

            # Which regime each row is in decides the band, and the regime is read off the
            # reference rather than assumed: on the floor ASPRO sits at 2*vis_cal_err.
            floor_v2 = 2 * comb.vis_cal_err
            on_floor = median(σv2) < 1.05 * floor_v2
            tol = on_floor ? 0.02 : 0.25
            @test abs(median(rv2) - 1) < tol
            @test abs(median(rcp) - 1) < tol

            # Per channel, tolerate the telluric outliers but not a broken channel.
            @test all(0.3 .< rv2 .< 3.0)
            @test all(isfinite, rv2)
        end
    end

    @testset "the fringe tracker engages where the instrument has one" begin
        # SPICA is the only shipped combiner with a tracker, and ASPRO requires it. Its effect
        # is the integration time, so that is what this asserts rather than a sigma.
        spica = read_comb_file("SPICA")
        mircx = read_comb_file("MIRCX")
        @test spica.ft_band == "H" && spica.ft_max_dit == 0.2
        @test isempty(mircx.ft_band)

        wav = read_wave_file("SPICA_LR")
        m8  = Dict("V" => 8.5, "R" => 8.0, "H" => 8.0)
        on  = predict_errors(fac, spica, wav; mag = m8, elevation_deg = 85.0, fringe_tracker = true)
        off = predict_errors(fac, spica, wav; mag = m8, elevation_deg = 85.0, fringe_tracker = false)
        @test on.dit[1]  == spica.ft_max_dit        # integrates to the tracker's limit
        @test off.dit[1] == spica.dit               # and to the frame rate without it
        @test median(on.sigma_v2) < median(off.sigma_v2)

        # Beyond the tracker's H limit it must NOT engage, or a faint target gets a sensitivity
        # it cannot have. The limit is in H even though SPICA observes at 0.6-0.9 um.
        faint = predict_errors(fac, spica, wav;
                               mag = Dict("V" => 12.5, "R" => 12.0, "H" => 12.0),
                               elevation_deg = 85.0, fringe_tracker = true)
        @test faint.dit[1] == spica.dit

        # An instrument with no tracker is unaffected by the switch.
        wmx = read_wave_file("MIRCX_LOWH")
        a = predict_errors(fac, mircx, wmx; mag = m8, elevation_deg = 85.0, fringe_tracker = true)
        b = predict_errors(fac, mircx, wmx; mag = m8, elevation_deg = 85.0, fringe_tracker = false)
        @test a.sigma_v2 == b.sigma_v2
    end
end
