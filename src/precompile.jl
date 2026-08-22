# This workload is the bulk of the package's precompilation cost, and that is deliberate.
# Measured on Julia 1.12 by precompiling OITOOLS cold both ways:
#
#                        precompile   cache on disk   `using OITOOLS`   first image_to_chi2
#     with workload         25.1 s        57 MB          2.9-3.3 s            1.33 s
#     workload skipped       4.5 s       5.4 MB           2.87 s             29.7 s
#
# `reconstruct` and `fit_model` were added afterwards: they are what the GUI's Run buttons
# call, and neither is reached by the chi2 entry points above. See the notes beside each.
#
# So it costs ~20 s once per precompile and saves ~28 s in every fresh session that touches
# the image reconstruction path, with no penalty to load time. Do not delete it to make
# `Pkg.precompile` look faster; if it is in the way while iterating on src/, switch it off
# locally instead — see the LocalPreferences.toml note in .gitignore.
#
# Note the `isfile` guard below fails *quietly*: if the data file is ever moved or renamed
# the workload silently becomes a no-op and time-to-first-call regresses ~20x with nothing
# in the build log to say so. demos/data/BC2026/ is tracked, so this holds for a plain
# checkout; it is the reason to be careful about relocating the demo data.
@compile_workload begin
    # BC2026 OBJECT1_N has OI_VIS, OI_VIS2, OI_T3, OI_FLUX across multiple channels
    oifitsfile = joinpath(@__DIR__, "..", "demos", "data", "BC2026", "OBJECT1_N.oifits")
    if isfile(oifitsfile)
        # warn=false as well as verbose=false: OBJECT1_N has station indices absent from its
        # OI_ARRAY, so without it every fresh precompile prints a dozen warnings that have
        # nothing to do with the user's own data.
        data = readoifits(oifitsfile, filter_bad_data=true, polychromatic=true,
                          verbose=false, warn=false, merge_oi_wavelength=true)
        nx = 32
        pixsize = 0.25
        ft = setup_ft(data, nx, pixsize)
        x = gaussian2d(nx, nx, nx ÷ 6)
        nwav = size(data, 1)
        x4d = reshape(repeat(x, 1, 1, nwav), nx, nx, nwav, 1)
        image_to_chi2(x4d, ft, data, verb=false)
        g4d = zero(x4d)
        image_to_chi2_fg(x4d, g4d, ft, data, verb=false)

        # The reconstruction driver itself. `image_to_chi2_fg` above is NOT the same code:
        # `reconstruct` wraps it in `crit_fg` (which folds in the regularisers) and hands that
        # to OptimPackNextGen's VMLMB, and neither is reached by any call above. Measured from
        # the GUI, the first Run spent one to two MINUTES compiling this before ten seconds of
        # actual work. Two iterations are enough: the cost is compilation, not iteration.
        reconstruct(x4d, data, ft; maxiter=2, verb=false)
        # …and once more with a regulariser, since the regularised criterion is a separate
        # specialisation and the panel offers it in the same breath as the plain one.
        reconstruct(x4d, data, ft; maxiter=2, verb=false,
                    regularizers=[["centering", 1e4]])

        # The MONOCHROMATIC read, which is a different branch from the polychromatic one above
        # and is what the GUI takes for an ordinary file. Cheap, and it is the very first thing
        # any session does.
        mono = joinpath(@__DIR__, "..", "demos", "data", "2004-data1.oifits")
        if isfile(mono)
            dmono = readoifits(mono; verbose=false, warn=false)
            ftm = setup_ft(dmono, 32, 0.3)
            image_to_chi2(gaussian2d(32, 32, 32 ÷ 6), ftm, dmono, verb=false)

            # The other reconstructors. The GUI's engine selector reaches all of them, and each
            # is a separate criterion with its own specialisations -- `reconstruct` above shares
            # almost nothing with them. Measured on a fresh session, first call against second:
            #
            #     reconstruct_bsmem    4.55 s of compilation
            #     reconstruct_bsdmm    3.92 s
            #     reconstruct_squeeze  2.11 s
            #
            # Two iterations each; the cost being bought off is compilation, not iteration.
            # Note the argument shapes differ between them and are not interchangeable: BSMEM
            # and SQUEEZE want the bin's cell (`ft[1,1]`) and one `OIdata`, BSDMM the whole
            # plan and the array.
            x32 = Float32.(gaussian2d(32, 32, 32 ÷ 6))
            x64 = Float64.(x32)
            reconstruct_bsmem(x64, dmono[1, 1], ftm[1, 1]; maxiter=2, verbose=false)
            # `mu_cen > 0` is a precondition, not a preference: with every weight at zero ADMM
            # has no proximal block and reduces over an empty collection.
            reconstruct_bsdmm(x32, dmono, ftm; mu_cen=1e2, maxit=2, verb=false)
            # `monitor = 0`: the in-package monitor draws with matplotlib, and precompilation
            # must not reach for Python.
            reconstruct_squeeze(x64, dmono[1, 1], ftm[1, 1];
                                niter=2, nchains=1, verb=false, monitor=0)
        end

        # Planning. Measured, the computation is free and the wait is entirely compilation:
        # `night_observability` costs 1.89 s on first call and 0.1 MILLISECOND thereafter,
        # `in_delay` 1.10 s against 0.0 ms. A coarse grid is used because the step changes how
        # long the workload runs at build time and not which code it compiles.
        #
        # The facility config ships inside the package, so unlike the OIFITS guards above this
        # branch cannot silently become a no-op.
        chara = read_facility_file("CHARA")
        obsv = night_observability(chara, 279.2347, 38.7837, DateTime(2026, 6, 21);
                                   step_minutes = 10)
        cfg4 = [1, 1, 1, 1, 0, 0]
        in_delay(chara, 38.7837, obsv.ha, cfg4, [1, 3, 3, 4, 1, 1])
        best_pop(chara, 38.7837, obsv.ha, cfg4)
        index_runs(obsv.good_alt)

        # Model fitting path
        model_dict = Dict{String,Any}("c1,ud" => 1.0, "c1,f" => 1.0)
        fm = parse_model(model_dict, ["c1,ud"])
        model_to_chi2(fm, [1.0], data[1,1], verb=false)
        model_to_obs(fm, [1.0], data[1,1])

        # The NLopt objective and its gradient, which is what a model fit actually runs — and
        # again not reached by `model_to_chi2` alone.
        fit_model(model_dict, ["c1,ud"], data[1,1];
                  lb=Dict("c1,ud"=>0.1), ub=Dict("c1,ud"=>10.0),
                  maxeval=2, verb=false)
    end
end
