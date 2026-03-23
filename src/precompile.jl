@compile_workload begin
    # ── Polychromatic path ───────────────────────────────────────────────
    oifitsfile_poly = joinpath(@__DIR__, "..", "demos", "data",
                               "2019_v1295Aql.WL_SMOOTH.A.oifits")
    if isfile(oifitsfile_poly)
        data = readoifits(oifitsfile_poly, filter_bad_data=true, polychromatic=true,
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
    end

    # ── Monochromatic path ───────────────────────────────────────────────
    oifitsfile_gray = joinpath(@__DIR__, "..", "demos", "data", "2004-data1.oifits")
    if isfile(oifitsfile_gray)
        data_gray = readoifits(oifitsfile_gray, filter_bad_data=true, warn=false, verbose=false)
        nx = 32
        pixsize = 0.25
        ft_gray = setup_ft(data_gray, nx, pixsize)
        x = gaussian2d(nx, nx, nx ÷ 6)
        image_to_chi2(x, ft_gray[1,1], data_gray[1,1], verb=false)
        g = zero(x)
        image_to_chi2_fg(x, g, ft_gray[1,1], data_gray[1,1], verb=false)

        # Model fitting path
        model_dict = Dict{String,Any}("c1,ud" => 1.0, "c1,f" => 1.0)
        fm = parse_model(model_dict, ["c1,ud"])
        model_to_chi2(fm, [1.0], data_gray[1,1], verb=false)
    end
end
