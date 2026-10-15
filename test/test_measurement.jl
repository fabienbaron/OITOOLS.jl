using OITOOLS, Test

# The measurement model's differential VISPHI: each spectral vector's phase residuals referenced
# to a polynomial fitted over its continuum channels (PMOIRED's DPHI, with the data's errors).
@testset "measurement model: continuum-referenced phases" begin
    src = joinpath(@__DIR__, "..", "demos", "data", "V1295_Aql.oifits")
    truth = Dict{String,Any}("A,ud" => 0.6, "B,f" => 0.3, "B,ud" => 0.2, "B,x" => 2.0, "B,y" => -1.5,
                             "B,line_l_wl0" => 1.65, "B,line_l_f" => 0.8, "B,line_l_gaussian" => 20.0)
    out = joinpath(mktempdir(), "sim.oifits")
    simulate_from_oifits(src, out; mode = "copy_errors", noise = false,
                         flat_model = dict_to_model(truth, String[]), flat_params = Float64[])
    d = readoifits(out; warn = false, verbose = false, T = Float64)[1, 1]
    g = spectral_groups(d)
    @test sort(vcat(g.vis...)) == 1:length(d.vis_mjd) && all(issorted(d.vis_lam[v]) for v in g.vis)

    # the model's own (absolute) phase, each vector moved by its own smooth polynomial, as a
    # pipeline's reference would (the simulated file's VISPHI follows the source's PHITYP)
    dd = deepcopy(d)
    Vt = eval_model(dict_to_model(truth, String[]), Float64[], d.uv; wl = d.uv_lam, mjd = d.uv_mjd)[d.indx_vis]
    dd.visphi .= rad2deg.(angle.(Vt))
    for (j, v) in enumerate(g.vis)
        λ = dd.vis_lam[v] .* 1e6
        dd.visphi[v] .+= (3.0 + j) .+ 2.0 .* (λ .- 1.6) .- 4.0 .* (λ .- 1.6) .^ 2
    end
    w = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    free = ["B,f", "B,x", "B,y", "B,line_l_f"]
    start = merge(copy(truth), Dict{String,Any}("B,f" => 0.25, "B,x" => 2.1, "B,y" => -1.4, "B,line_l_f" => 0.6))
    mm = measurement_model(visphi = :continuum, order = 2)
    ra = fit_model_lsqfit(start, free, dd; weights = w)
    rc = fit_model_lsqfit(start, free, dd; weights = w, measurement = mm)
    @test ra.chi2 > 1e3 * max(rc.chi2, 1e-6)                   # absolute: the polynomials are fitted as signal
    @test rc.x_opt ≈ [0.3, 2.0, -1.5, 0.8] rtol = 1e-3          # referenced: the truth back
    @test rc.ndof == OITOOLS._ndof(dd, w) - 3 * length(g.vis)   # three coefficients per vector

    # the referenced residuals and their Jacobian (a linear map per vector) against finite differences
    fm = dict_to_model(start, free)
    x = Float64[start[k] for k in free]
    ref = OITOOLS._ref_context(mm, fm, dd)
    n = length(OITOOLS.residuals_flat(fm, x, dd; weights = w))
    r = zeros(n); J = zeros(n, length(x))
    OITOOLS.residuals_flat_jac!(r, J, fm, x, dd; weights = w, ref)
    r2 = zeros(n); OITOOLS.residuals_flat!(r2, fm, x, dd; weights = w, ref)
    @test r ≈ r2
    for k in eachindex(x)
        h = 1e-6 * max(1, abs(x[k])); xp = copy(x); xm = copy(x); xp[k] += h; xm[k] -= h
        rp = zeros(n); rm = zeros(n)
        OITOOLS.residuals_flat!(rp, fm, xp, dd; weights = w, ref); OITOOLS.residuals_flat!(rm, fm, xm, dd; weights = w, ref)
        @test J[:, k] ≈ (rp .- rm) ./ (2h) rtol = 1e-5 atol = 1e-6
    end

    # the continuum leaves the line out: the line's channels do not set the reference
    ctx = OITOOLS._ref_context(measurement_model(visphi = :continuum), fm, dd)
    vals = fm.resolver.fn(x)
    v = first(gg for gg in g.vis if length(gg) > 3)
    keep = OITOOLS._continuum_mask(ctx, v, vals)
    λ = dd.vis_lam[v] .* 1e6
    @test all(keep .== (abs.(λ .- 1.65) .>= 1.5 * 20.0 / 1000)) && count(keep) == length(v) - 1
    @test_throws ArgumentError measurement_model(visphi = :nope)

    @testset "PMOIRED's spectral kernel and DPHI (compatibility)" begin
        ker = OITOOLS._pm_kernel(2.0)
        @test length(ker) == 2 * 4 + 3 && sum(ker) ≈ 1 && ker ≈ reverse(ker)
        y = collect(1.0:12.0)
        @test OITOOLS._pm_convolve(fill(3.0, 12), ker) ≈ fill(3.0, 12)            # a constant stays
        c = OITOOLS._pm_convolve(y, ker)
        @test c[6] ≈ sum(ker .* y[1:11])                                            # interior: plain
        p = [reverse(y[1:11]); y; reverse(y[2:12])]                                 # edges mirrored
        @test c[1] ≈ sum(ker .* [reverse(y[1:11]); y][7:17])
        # the smoothed residuals' Jacobian against finite differences
        mk = measurement_model(visphi = :continuum, spectral_kernel = 1.5)
        fm2 = dict_to_model(start, free); x2 = Float64[start[k] for k in free]
        w2 = [1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0]
        refk = OITOOLS._ref_context(mk, fm2, dd)
        n2 = length(OITOOLS.residuals_flat(fm2, x2, dd; weights = w2))
        r2 = zeros(n2); J2 = zeros(n2, length(x2))
        OITOOLS.residuals_flat_jac!(r2, J2, fm2, x2, dd; weights = w2, ref = refk)
        rv = zeros(n2); OITOOLS.residuals_flat!(rv, fm2, x2, dd; weights = w2, ref = refk)
        @test rv ≈ r2
        for k in eachindex(x2)
            h = 1e-6 * max(1, abs(x2[k])); xp = copy(x2); xm = copy(x2); xp[k] += h; xm[k] -= h
            rp = zeros(n2); rm = zeros(n2)
            OITOOLS.residuals_flat!(rp, fm2, xp, dd; weights = w2, ref = refk)
            OITOOLS.residuals_flat!(rm, fm2, xm, dd; weights = w2, ref = refk)
            @test J2[:, k] ≈ (rp .- rm) ./ (2h) rtol = 1e-5 atol = 1e-6
        end
        # pmoired = true: (d − P_w d) − (m − P_u m) over every channel, against a direct computation
        mp = measurement_model(visphi = :continuum, order = 2, pmoired = true)
        wp = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        refp = OITOOLS._ref_context(mp, fm2, dd)
        np_ = length(OITOOLS.residuals_flat(fm2, x2, dd; weights = wp))
        rp = zeros(np_); OITOOLS.residuals_flat!(rp, fm2, x2, dd; weights = wp, ref = refp)
        V = eval_model(fm2, x2, dd.uv; wl = dd.uv_lam, mjd = dd.uv_mjd)[dd.indx_vis]
        want = similar(rp); ok = falses(length(rp))
        for gv in g.vis
            λ = dd.vis_lam[gv] .* 1e6; t = (λ .- sum(λ) / length(λ)) ./ (maximum(λ) - minimum(λ))
            A = hcat(ones(length(gv)), t, t .^ 2); σ = dd.visphi_err[gv]
            m = rad2deg.(angle.(V[gv])); dv = copy(dd.visphi[gv])
            # unwrapping is unique only where neighbouring channels differ by well under 180°
            steps(a) = maximum(abs.(mod.(diff(a) .+ 180, 360) .- 180); init = 0.0)
            (steps(m) < 90 && steps(dv) < 90 && steps(m .- dv) < 90) || continue
            for i in 2:length(gv); m[i] = m[i-1] + mod(m[i] - m[i-1] + 180, 360) - 180; dv[i] = dv[i-1] + mod(dv[i] - dv[i-1] + 180, 360) - 180; end
            m .+= 360 .* round.((dv[1] - m[1]) / 360)                               # the data's branch
            e = (m .- A * (A \ m)) .- (dv .- A * ((A ./ σ) \ (dv ./ σ)))
            want[gv] = e ./ σ; ok[gv] .= true
        end
        @test count(ok) > length(rp) / 2
        @test rp[ok] ≈ want[ok] rtol = 1e-8 atol = 1e-8
        @test_throws ArgumentError measurement_model(spectral_kernel = -1)
    end
end
