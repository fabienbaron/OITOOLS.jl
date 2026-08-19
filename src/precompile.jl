# This workload is the bulk of the package's precompilation cost, and that is deliberate.
# Measured on Julia 1.12 by precompiling OITOOLS cold both ways:
#
#                        precompile   cache on disk   `using OITOOLS`   first image_to_chi2
#     with workload         25.1 s        57 MB          2.9-3.3 s            1.33 s
#     workload skipped       4.5 s       5.4 MB           2.87 s             29.7 s
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

        # Model fitting path
        model_dict = Dict{String,Any}("c1,ud" => 1.0, "c1,f" => 1.0)
        fm = parse_model(model_dict, ["c1,ud"])
        model_to_chi2(fm, [1.0], data[1,1], verb=false)
        model_to_obs(fm, [1.0], data[1,1])
    end
end
