# test_simulate.jl — coverage for simulate() and the CHARA noise model.
#
# simulate() previously had no tests at all. These lock down the things that were actually
# broken: it did not run, RA units were inconsistent, closure-phase errors were computed from
# one hour angle, and the observables were noised independently so they contradicted each other.

using OITOOLS, Test, Dates, Statistics, Random
using SpecialFunctions: besselj1
using OIFITS

const _TMP = mktempdir()

function _setup(; comb="MIRCX", wave="MIRCX_LOWH", facility="CHARA", nep=6)
    dates = collect(DateTime(2018,8,13,3,0,0):Minute(30):
                    DateTime(2018,8,13,3,0,0) + Minute(30*(nep-1)))
    t = read_obs_file("default_obs")
    t.target = "TESTSTAR"
    t.raep0  = 15*(20 + 57/60 + 59.4437981/3600)     # degrees
    t.decep0 = 46 + 28/60 + 0.5731825/3600
    return read_facility_file(facility), t, read_comb_file(comb), read_wave_file(wave), dates
end

_ud(d) = dict_to_model(Dict{String,Any}("star,ud"=>d, "star,f"=>1.0), String[])

@testset "simulate" begin

    # ── Step 0: it runs at all, and the time/coordinate plumbing is right ────
    @testset "smoke and epoch plumbing" begin
        f, t, c, w, dates = _setup()
        out = joinpath(_TMP, "smoke.oifits")
        simulate(f, t, c, w, dates, out; flat_model=_ud(2.0), flat_params=Float64[], seed=1)
        @test isfile(out)
        d = readoifits(out; T=Float64, filter_bad_data=false)[1,1]
        @test d.nv2 == 15 * length(dates) * length(w.λ)

        # MJD must be the real MJD, not an AstroTime leftover.
        @test minimum(d.v2_mjd) ≈ 58343.125 atol=1e-6

        # RA is in degrees, and hour_angle_calc consumes hours: a 24h RA shift is a no-op.
        _, ha  = OITOOLS.hour_angle_calc(dates, f.lon, t.raep0/15)
        _, ha2 = OITOOLS.hour_angle_calc(dates, f.lon, t.raep0/15 + 24)
        @test all(abs.(mod.(ha .- ha2 .+ 12, 24) .- 12) .< 1e-9)
        @test -12 <= minimum(ha) && maximum(ha) <= 12
    end

    # ── Geometry ─────────────────────────────────────────────────────────────
    @testset "uv geometry" begin
        f, t, c, w, dates = _setup(nep=4)
        out = joinpath(_TMP, "geom.oifits")
        simulate(f, t, c, w, dates, out; flat_model=_ud(1.0), flat_params=Float64[], noise=false)
        d = readoifits(out; T=Float64, filter_bad_data=false)[1,1]
        # Longest projected baseline cannot exceed the longest physical one.
        Bmax = maximum(sqrt.(sum(abs2, f.sta_xyz[i,:] .- f.sta_xyz[j,:]))
                       for i in 1:f.ntel for j in i+1:f.ntel)
        @test maximum(sqrt.(d.uv[1,:].^2 .+ d.uv[2,:].^2)) <= Bmax / minimum(w.λ) * 1.0001
    end

    # ── Analytic observables, noise off ──────────────────────────────────────
    @testset "analytic observables" begin
        f, t, c, w, dates = _setup(nep=3)
        # Small enough to stay inside the first null of the visibility function, so that
        # V > 0 everywhere and the closure phase of a centro-symmetric source is 0 rather
        # than the equally-correct +-180 you get once V changes sign.
        θmas = 0.3
        out = joinpath(_TMP, "ud.oifits")
        simulate(f, t, c, w, dates, out; flat_model=_ud(θmas), flat_params=Float64[], noise=false)
        d = readoifits(out; T=Float64, filter_bad_data=false)[1,1]

        # A uniform disc has V = 2 J1(x)/x with x = pi * B * theta / lambda.
        θrad = θmas * pi/180/3600/1000
        B    = sqrt.(d.v2_baseline.^2)
        x    = pi .* B .* θrad
        Vth  = 2 .* besselj1.(x) ./ x
        @test maximum(abs.(d.v2 .- Vth.^2)) < 1e-6

        # Centro-symmetric brightness => closure phase identically zero.
        @test maximum(abs.(d.t3phi)) < 1e-4
    end

    # ── Closure-phase convention ────────────────────────────────────────────
    # simulate builds its triangles as legs (i,j), (j,k), (i,k), so the third visibility has
    # to be conjugated. Getting this wrong flips the closure phase but leaves V2 and T3AMP
    # untouched, so only an asymmetric source with a known analytic closure catches it.
    @testset "closure phase convention" begin
        f, t, c, w, dates = _setup(nep=4)
        dx, dy, r = 3.0, 1.0, 0.3                       # mas, mas, flux ratio
        m = dict_to_model(Dict{String,Any}("a,f"=>1/(1+r), "b,f"=>r/(1+r),
                                           "b,x"=>dx, "b,y"=>dy), String[])
        out = joinpath(_TMP, "binary.oifits")
        simulate(f, t, c, w, dates, out; flat_model=m, flat_params=Float64[], noise=false)

        ds  = OIFITS.OIDataSet(out)
        tb  = ds.t3[1]
        lam = ds.instr[1].eff_wave
        mas = pi/180/3600/1000
        Vb(u,v,l) = (1 + r*cis(-2pi*(u*dx*mas + v*dy*mas)/l))/(1+r)

        worst = 0.0
        for k in eachindex(tb.u1coord), j in eachindex(lam)
            u1,v1,u2,v2 = tb.u1coord[k], tb.v1coord[k], tb.u2coord[k], tb.v2coord[k]
            u3,v3 = -(u1+u2), -(v1+v2)                  # OIFITS: the three legs close
            cp = angle(Vb(u1,v1,lam[j])*Vb(u2,v2,lam[j])*Vb(u3,v3,lam[j])) * 180/pi
            d  = tb.t3phi[j,k] - cp;  d -= 360*round(d/360)
            worst = max(worst, abs(d))
        end
        # Tolerance is set well below any convention error (getting the conjugate wrong
        # shifts this by tens of degrees) but above the ~1e-5 deg of floating-point
        # disagreement between the model evaluator and this closed form.
        @test worst < 1e-3
    end

    # ── Step 5: one perturbation, mutually consistent observables ────────────
    @testset "observable self-consistency" begin
        f, t, c, w, dates = _setup(nep=5)
        out = joinpath(_TMP, "cons.oifits")
        simulate(f, t, c, w, dates, out; flat_model=_ud(2.5), flat_params=Float64[],
                 mag=5.0, seed=7)
        d = readoifits(out; T=Float64, filter_bad_data=false)[1,1]
        # VISAMP and VIS2 are the same measurement written two ways. They agree exactly
        # wherever the debiased V2 is non-negative; where debiasing pushes V2 below zero
        # (normal for a resolved baseline at low SNR) VISAMP is clamped at 0 and cannot.
        ok = d.v2 .>= 0
        @test maximum(abs.(d.visamp[ok].^2 .- d.v2[ok])) < 1e-10
        @test all(d.visamp .>= 0)
        # Phases stay in range.
        @test all(-180 .<= d.t3phi .<= 180)
        @test all(-180 .<= d.visphi .<= 180)
    end

    # ── Noise statistics: the error bars describe the actual scatter ─────────
    @testset "noise is calibrated" begin
        f, t, c, w, dates = _setup(nep=5)
        m = _ud(2.5)
        truth = joinpath(_TMP, "truth.oifits");  noisy = joinpath(_TMP, "noisy.oifits")
        simulate(f, t, c, w, dates, truth; flat_model=m, flat_params=Float64[], mag=5.0, noise=false)
        simulate(f, t, c, w, dates, noisy; flat_model=m, flat_params=Float64[], mag=5.0, seed=11)
        d0 = readoifits(truth; T=Float64, filter_bad_data=false)[1,1]
        d1 = readoifits(noisy; T=Float64, filter_bad_data=false)[1,1]

        z = (d1.v2 .- d0.v2) ./ d1.v2_err
        @test abs(mean(z)) < 0.15          # debiased: no systematic offset
        @test 0.85 < std(z) < 1.15         # error bars neither optimistic nor inflated

        zc = d1.t3phi .- d0.t3phi
        zc .= zc .- 360 .* round.(zc ./ 360)
        zc ./= d1.t3phi_err
        @test abs(mean(zc)) < 0.15
        @test 0.75 < std(zc) < 1.25
    end

    @testset "reproducibility" begin
        f, t, c, w, dates = _setup(nep=3)
        m = _ud(2.0)
        o1 = joinpath(_TMP,"r1.oifits"); o2 = joinpath(_TMP,"r2.oifits"); o3 = joinpath(_TMP,"r3.oifits")
        simulate(f,t,c,w,dates,o1; flat_model=m, flat_params=Float64[], seed=99)
        simulate(f,t,c,w,dates,o2; flat_model=m, flat_params=Float64[], seed=99)
        simulate(f,t,c,w,dates,o3; flat_model=m, flat_params=Float64[], seed=100)
        a = readoifits(o1;T=Float64,filter_bad_data=false)[1,1]
        b = readoifits(o2;T=Float64,filter_bad_data=false)[1,1]
        cc= readoifits(o3;T=Float64,filter_bad_data=false)[1,1]
        @test a.v2 == b.v2
        @test a.v2 != cc.v2
    end

    # ── Step 2(a): errors must vary along the track ──────────────────────────
    @testset "errors depend on hour angle" begin
        f, t, c, w, dates = _setup(nep=6)
        out = joinpath(_TMP, "ha.oifits")
        simulate(f, t, c, w, dates, out; flat_model=_ud(2.5), flat_params=Float64[], noise=false)
        d = readoifits(out; T=Float64, filter_bad_data=false)[1,1]
        # Before the fix these were computed at hour angle #1 and tiled, giving very few
        # distinct values. The source is resolved and the elevation changes, so they vary.
        @test length(unique(round.(d.t3phi_err, digits=8))) > length(dates)
        @test length(unique(round.(d.v2_err,    digits=8))) > length(dates)
    end

    # ── Step 2(b): photon background counts every telescope ──────────────────
    @testset "photon background scales with telescope count" begin
        f, t, c, w, dates = _setup(nep=2)
        mags = fill(5.0, length(w.λ))
        N, _, nfr = OITOOLS.photons_per_telescope(mags, w.λ, w.δλ, f, c, [90.0], 5.0)
        Np, _, _  = OITOOLS.photons_per_telescope(mags, w.λ, w.δλ, f, c, [90.0], 5.0;
                                                  channel=:photometry)
        v2st = hcat([[i,j] for i in 1:f.ntel for j in i+1:f.ntel]...)
        _, vc6, _, _ = OITOOLS.correlated_flux_coefficients(N, Np, c, 6, v2st)
        _, vc3, _, _ = OITOOLS.correlated_flux_coefficients(N, Np, c, 3, v2st)
        # More beams overlapping in the interferometric channel => more photon noise.
        @test vc6[1,1,1] > vc3[1,1,1]
    end

    # ── Step 3: Strehl against ASPRO's published CHARA curves ────────────────
    @testset "Strehl vs ASPRO" begin
        ao = AOConfig()
        seeings = [0.60, 0.70, 1.00, 1.15, 1.40, 1.80]
        t0s     = [5.2, 4.4, 3.2, 2.2, 1.6, 1.0] .* 1e-3
        ref = Dict(
          (1.57e-6) => [(0.770,0.330),(0.710,0.270),(0.555,0.150),(0.395,0.120),(0.235,0.083),(0.078,0.053)],
          (2.18e-6) => [(0.878,0.520),(0.828,0.443),(0.683,0.280),(0.518,0.228),(0.337,0.167),(0.140,0.107)])
        worst = 0.0
        for (lam, r) in ref, i in eachindex(seeings)
            sb = strehl_ratio(lam, 1.0, seeings[i], t0s[i],  2.0, ao)
            sf = strehl_ratio(lam, 1.0, seeings[i], t0s[i], 16.0, ao)
            worst = max(worst, abs(sb-r[i][1])/r[i][1], abs(sf-r[i][2])/r[i][2])
        end
        @test worst < 0.05                       # H and K agree to better than 5%

        # Monotonic in seeing, wavelength and airmass; bounded.
        @test strehl_ratio(1.6e-6,1.0,0.6,5.2e-3,5.0,ao) > strehl_ratio(1.6e-6,1.0,1.4,1.6e-3,5.0,ao)
        @test strehl_ratio(2.2e-6,1.0,1.0,3.2e-3,5.0,ao) > strehl_ratio(1.0e-6,1.0,1.0,3.2e-3,5.0,ao)
        @test strehl_ratio(1.6e-6,1.0,1.0,3.2e-3,5.0,ao; elevation_deg=90) >
              strehl_ratio(1.6e-6,1.0,1.0,3.2e-3,5.0,ao; elevation_deg=20)
        for l in [0.6e-6, 1.6e-6, 13.0e-6], D in [1.0, 1.8, 8.0]
            @test 0 <= strehl_ratio(l,D,1.0,3.2e-3,6.0,ao) <= 1
        end

        # No [ao] block must degrade gracefully to the seeing-limited coupling, not error.
        @test coupling_efficiency(1.6e-6,1.0,1.0,3.2e-3,5.0,nothing) ≈
              min(1.0,(fried_parameter(1.0,1.6e-6)/1.0)^2)
    end

    @testset "photometric zero points" begin
        # A zeroth-magnitude star gives ~1000 photons/s/cm^2/Angstrom in V.
        V = band_by_name("V")
        @test V.f0 * 1e-6 ≈ 1.005e11 rtol=0.02          # ph/s/m^2/um
        @test band_for_wavelength(1.57e-6).name == "H"
        @test band_for_wavelength(2.18e-6).name == "K"
        @test zero_point_flux(0.55e-6) ≈ 1.005e11 rtol=0.02
    end

    # ── Step 4: magnitudes per band and per channel ──────────────────────────
    @testset "magnitude specification" begin
        λ = [1.1e-6, 1.6e-6, 2.2e-6]
        @test OITOOLS.resolve_magnitudes(4.0, λ) == fill(4.0, 3)
        @test OITOOLS.resolve_magnitudes([1.0,2.0,3.0], λ) == [1.0,2.0,3.0]
        m = OITOOLS.resolve_magnitudes(Dict("J"=>3.0,"H"=>2.0,"K"=>1.0), λ)
        @test m[1] > m[2] > m[3]                  # brighter towards K, monotonic
        @test m[2] ≈ 2.0 atol=0.1                 # H channel picks up the H magnitude
        @test_throws ArgumentError OITOOLS.resolve_magnitudes([1.0,2.0], λ)

        # A fainter star must produce larger error bars.
        f, t, c, w, dates = _setup(nep=2)
        o5 = joinpath(_TMP,"m5.oifits"); o8 = joinpath(_TMP,"m8.oifits")
        simulate(f,t,c,w,dates,o5; flat_model=_ud(2.0), flat_params=Float64[], mag=5.0, noise=false)
        simulate(f,t,c,w,dates,o8; flat_model=_ud(2.0), flat_params=Float64[], mag=8.0, noise=false)
        d5 = readoifits(o5;T=Float64,filter_bad_data=false)[1,1]
        d8 = readoifits(o8;T=Float64,filter_bad_data=false)[1,1]
        @test mean(d8.v2_err) > mean(d5.v2_err)
    end

    # ── Step 1: observability is opt-in and never on by default ─────────────
    @testset "observability opt-in" begin
        f, t, c, w, dates = _setup(nep=12)
        # Default: every epoch is used, however unobservable.
        out = joinpath(_TMP, "noobs.oifits")
        simulate(f, t, c, w, dates, out; flat_model=_ud(2.0), flat_params=Float64[], noise=false)
        d = readoifits(out; T=Float64, filter_bad_data=false)[1,1]
        @test d.nv2 == 15 * length(dates) * length(w.λ)

        # With a cut, epochs below the limit are dropped, and it matches alt_az directly.
        sel = observable_epochs(f, t, dates; min_elevation=60.0)
        alt, _ = alt_az(t.decep0, f.lat, OITOOLS.hour_angle_calc(dates, f.lon, t.raep0/15)[2])
        @test sel.report.n_out == count(alt .>= 60.0)
        @test length(sel.dates) == sel.report.n_out

        if sel.report.n_out > 0
            out2 = joinpath(_TMP, "obs.oifits")
            simulate(f, t, c, w, dates, out2; flat_model=_ud(2.0), flat_params=Float64[],
                     noise=false, observability=(min_elevation=60.0,))
            d2 = readoifits(out2; T=Float64, filter_bad_data=false)[1,1]
            @test d2.nv2 == 15 * sel.report.n_out * length(w.λ)
        end

        # No constraints at all is a no-op.
        @test observable_epochs(f, t, dates).report.n_out == length(dates)
        @test_throws ArgumentError observable_epochs(f, t, dates; pops=[1,2,3])
    end

    # ── get_baselines: used to throw BoundsError without a reference cart ────
    @testset "get_baselines" begin
        f, _, _, _, _ = _setup()
        nb, xyz, sta, names = OITOOLS.get_baselines(f)
        @test nb == f.ntel*(f.ntel-1)÷2
        @test all(sta[1,:] .< sta[2,:])                 # same i<j convention as simulate
        nb2, _, sta2, _ = OITOOLS.get_baselines(f; config=[1,1,1,1,1,2])
        @test nb2 == f.ntel - 1
        @test all(sta2[1,:] .== 6)
    end

    # ── Step 4: a cube's spectrum drives OI_FLUX and the per-channel noise ───
    @testset "image cube SED" begin
        f, t, c, w, dates = _setup(nep=2)
        nx, nw = 32, length(w.λ)
        cube = zeros(nx, nx, nw)
        for k in 1:nw
            cube[nx÷2+1, nx÷2+1, k] = 1.0 + 2.0*(k-1)/(nw-1)   # rises with wavelength
        end
        out = joinpath(_TMP, "cube.oifits")
        simulate(f, t, c, w, dates, out; image=cube, pixsize=0.2, noise=false, mag=4.0)
        d = readoifits(out; T=Float64, filter_bad_data=false)[1,1]
        fl = reshape(d.flux, nw, :)
        @test fl[end,1] > fl[1,1]                    # OI_FLUX follows the cube spectrum
        @test !all(d.flux .== 1.0)                   # no longer hard-coded to 1
    end

    @testset "deprecated nonoise" begin
        f, t, c, w, dates = _setup(nep=2)
        out = joinpath(_TMP, "dep.oifits")
        simulate(f, t, c, w, dates, out; flat_model=_ud(2.0), flat_params=Float64[], nonoise=true)
        d = readoifits(out; T=Float64, filter_bad_data=false)[1,1]
        # The old spelling still turns noise off, but now writes real error bars rather
        # than the placeholder 1.0 it used to.
        @test all(d.v2_err .< 1.0)
        @test all(isfinite, d.v2)
    end
end
