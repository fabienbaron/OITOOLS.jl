# Nested sampling, both backends.
#
# `fit_model_nested` has two implementations behind one name — Nautilus.jl in pure Julia,
# UltraNest through PythonCall — and the point of these tests is the comparison. A single
# sampler's log-evidence has nothing to check it against; two independent ones agreeing is what
# makes either number usable, and a disagreement here is a real finding rather than a flaky test.
#
# Both are weak dependencies, declared in this package's test target. If one is missing its
# testsets skip rather than fail: the suite must still run in an environment that has only one.

const _NESTED_FILE = joinpath(@__DIR__, "oifits_for_tests", "AlphaCenA.oifits")

# The resampler lives in the extension module, not in OITOOLS itself, so it is reached the
# same way test/gui/ reaches the GUI extension's internals. At file top level, because a
# `const` inside a `@testset begin` block is a local declaration and a syntax error.
const NSEXT = Base.get_extension(OITOOLS, :OITOOLSNautilusExt)

@testset "nested sampling" begin

nested_data = readoifits(_NESTED_FILE; warn = false, verbose = false)[1, 1]
nested_md   = Dict{String,Any}("s,ud" => 6.0, "s,f" => 1.0)
nested_free = ["s,ud"]
nested_lb   = Dict("s,ud" => 4.0)
nested_ub   = Dict("s,ud" => 12.0)

have(b) = b in OITOOLS.NESTED_BACKENDS_LOADED

@testset "backend registry" begin
    # The test target declares both, so both must have registered. If this fails the
    # extensions are not loading, and every skip below would be hiding that.
    @test have(:nautilus)
    @test have(:ultranest)
    @test nested_backend() in (:nautilus, :ultranest)

    # Asking for a backend that is not loaded says so immediately. The alternative — a
    # MethodError from deep inside the dispatch after the arguments have been built — is what
    # this exists to prevent.
    @test_throws ErrorException set_nested_backend!(:no_such_sampler)
    @test_throws ErrorException fit_model_nested(nested_md, nested_free, nested_data;
                                                 backend = :no_such_sampler)

    old = nested_backend()
    @test set_nested_backend!(:nautilus) === :nautilus
    @test nested_backend() === :nautilus
    set_nested_backend!(old)
end

@testset "result type" begin
    # One type, two producers. The old name has to keep resolving or every script and demo
    # written against the UltraNest-only API breaks.
    @test UltraNestResult === NestedResult
    @test :backend in fieldnames(NestedResult)
end

# The comparison. Deliberately a small run: enough live points for the two to agree on a
# one-parameter problem, few enough that the suite does not grow a minute.
if have(:nautilus) && have(:ultranest)
    @testset "the two samplers agree" begin
        rj = fit_model_nested(nested_md, nested_free, nested_data;
                              backend = :nautilus, lb = nested_lb, ub = nested_ub,
                              verb = false, cornerplot = false, seed = 42)
        ru = fit_model_nested(nested_md, nested_free, nested_data;
                              backend = :ultranest, lb = nested_lb, ub = nested_ub,
                              min_num_live_points = 150, verb = false, cornerplot = false)

        @test rj.backend === :nautilus
        @test ru.backend === :ultranest

        for r in (rj, ru)
            @test length(r.x_opt) == 1
            @test r.list_free_params == nested_free
            @test nested_lb["s,ud"] <= r.x_opt[1] <= nested_ub["s,ud"]
            @test size(r.posterior, 2) == 1
            @test size(r.posterior, 1) > 100
            @test isfinite(r.logz) && isfinite(r.logzerr) && r.logzerr > 0
            @test r.chi2r ≈ r.chi2 / r.ndof
        end

        # The maximum-likelihood point. Both samplers are looking at the same χ² surface, so
        # this is a much tighter tolerance than the evidence: a disagreement here means one of
        # them is not fitting what the other is.
        @test rj.x_opt[1] ≈ ru.x_opt[1] rtol = 1e-3
        @test rj.chi2 ≈ ru.chi2 rtol = 1e-3

        # The evidence, at 3σ of the two reported uncertainties combined. This is the check
        # that says the pure-Julia backend is usable at all.
        @test abs(rj.logz - ru.logz) <= 3 * sqrt(rj.logzerr^2 + ru.logzerr^2)
    end
end

if have(:nautilus)
    @testset "FlatModel form and bounds" begin
        fm = dict_to_model(nested_md, nested_free)
        r  = fit_model_nested(fm, nested_data; backend = :nautilus,
                              lb = nested_lb, ub = nested_ub, verb = false, cornerplot = false, seed = 7)
        @test r.backend === :nautilus
        @test nested_lb["s,ud"] <= r.x_opt[1] <= nested_ub["s,ud"]

        # The bounds ARE the prior, so a missing one is an error and not a default. Silently
        # inventing a range would decide the answer.
        @test_throws ErrorException fit_model_nested(nested_md, nested_free, nested_data;
                                                     backend = :nautilus,
                                                     ub = nested_ub, verb = false)
        @test_throws ErrorException fit_model_nested(nested_md, nested_free, nested_data;
                                                     backend = :nautilus,
                                                     lb = nested_lb, verb = false)
    end
end

end   # nested sampling
