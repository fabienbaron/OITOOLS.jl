# The six YSO templates of demos/models (Ibrahim et al. 2023, HD 190073): each reads, parses
# without warnings, and keeps the three fluxes summing to 1 at every wavelength; one is fitted to
# V1295 Aql (HD 190073) to the paper's quality.
using OITOOLS, Test

@testset "YSO templates (Ibrahim et al. 2023)" begin
    dir = joinpath(@__DIR__, "..", "demos", "models")
    names = ["YSO_Doughnut", "YSO_DoubleSigmoid", "YSO_AlphaDoubleSigmoid",
             "YSO_SkewedRing", "YSO_SkewedRing2", "YSO_OffCenterStar"]
    uv0 = zeros(2, 3); wl = [1.5e-6, 1.6e-6, 1.7e-6]
    for n in names
        m = read_model_file(joinpath(dir, n * ".toml"))
        @test isempty(OITOOLS.model_warnings(m.model, m.free; lb = m.lb, ub = m.ub))
        @test all(k -> haskey(m.lb, k) && haskey(m.ub, k), m.free)        # every free row bounded
        fm = dict_to_model(m.model, String[])
        @test real.(eval_model(fm, Float64[], uv0; wl)) ≈ ones(3)          # V(0) = 1 at any λ
    end
    d = readoifits(joinpath(@__DIR__, "..", "demos", "data", "V1295_Aql.oifits");
                   warn = false, verbose = false, T = Float64)[1, 1]
    m = read_model_file(joinpath(dir, "YSO_DoubleSigmoid.toml"))
    r = fit_model_lsqfit(m.model, m.free, d; lb = m.lb, ub = m.ub,
                         weights = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    @test r.converged && r.chi2r < 1.2
end
