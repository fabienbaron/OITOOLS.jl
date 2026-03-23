@compile_workload begin
    oifitsfile = joinpath(@__DIR__, "..", "demos", "data",
                          "2019_v1295Aql.WL_SMOOTH.A.oifits")
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
    end
end
