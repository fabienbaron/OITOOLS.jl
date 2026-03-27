@compile_workload begin
    # BC2026 OBJECT1_N has OI_VIS, OI_VIS2, OI_T3, OI_FLUX across multiple channels
    oifitsfile = joinpath(@__DIR__, "..", "demos", "data", "BC2026", "OBJECT1_N.oifits")
    if isfile(oifitsfile)
        data = readoifits(oifitsfile, filter_bad_data=true, polychromatic=true,
                          verbose=false, merge_oi_wavelength=true)
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
