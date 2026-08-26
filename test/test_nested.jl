# Nested sampling, both backends.
#
# `fit_model_nested` has two implementations behind one name — NestedSamplers.jl in pure Julia,
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
const NSEXT = Base.get_extension(OITOOLS, :OITOOLSNestedSamplersExt)

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
    @test have(:nestedsamplers)
    @test have(:ultranest)
    @test nested_backend() in (:nestedsamplers, :ultranest)

    # Asking for a backend that is not loaded says so immediately. The alternative — a
    # MethodError from deep inside the dispatch after the arguments have been built — is what
    # this exists to prevent.
    @test_throws ErrorException set_nested_backend!(:no_such_sampler)
    @test_throws ErrorException fit_model_nested(nested_md, nested_free, nested_data;
                                                 backend = :no_such_sampler)

    old = nested_backend()
    @test set_nested_backend!(:nestedsamplers) === :nestedsamplers
    @test nested_backend() === :nestedsamplers
    set_nested_backend!(old)
end

@testset "result type" begin
    # One type, two producers. The old name has to keep resolving or every script and demo
    # written against the UltraNest-only API breaks.
    @test UltraNestResult === NestedResult
    @test :backend in fieldnames(NestedResult)
end

if NSEXT !== nothing
@testset "equal-weight resampling" begin
    R = NSEXT._resample_equal
    rng = Xoshiro(1)

    # Uniform weights must return every sample exactly once: a systematic comb of n teeth over
    # n equal buckets lands one in each, whatever the offset.
    x = reshape(Float64.(1:100), 100, 1)
    w = fill(1 / 100, 100)
    @test sort(vec(R(x, w, rng))) == Float64.(1:100)

    # A weight vector concentrated on one sample yields only that sample
    w2 = zeros(100); w2[7] = 1.0
    @test all(==(7.0), R(x, w2, rng))

    # Shape and column alignment survive
    x3 = hcat(Float64.(1:50), Float64.(51:100))
    r3 = R(x3, fill(1 / 50, 50), rng)
    @test size(r3) == (50, 2)
    @test all(r3[:, 2] .== r3[:, 1] .+ 50)

    # Unnormalised weights are normalised rather than rejected
    @test size(R(x, fill(3.0, 100), rng)) == (100, 1)
    @test_throws ErrorException R(x, zeros(100), rng)
    @test_throws ErrorException R(x, fill(0.1, 99), rng)
end
end

if NSEXT !== nothing
@testset "threaded rejection proposal" begin
    TR   = NSEXT.ThreadedRejection
    pick = NSEXT._pick_proposal

    if Threads.nthreads() > 1
        # Substituted only where `:auto` would have picked `Proposals.Rejection` anyway
        @test pick(:auto, 8, 2) isa TR
        @test pick(:auto, 8, 2).batch == 8
        @test pick(:auto, 8, 12) === :auto      # RWalk territory: not a rejection sampler
        @test pick(:auto, 1, 2)  === :auto      # batch = 1 is the stock serial proposal
    else
        # Nothing to spread the batch over, so it must not pay for batching
        @test pick(:auto, 8, 2) === :auto
    end

    # An explicit choice is never overridden, whatever the thread count
    rw = NestedSamplers.Proposals.RWalk()
    @test pick(rw, 8, 2) === rw

    @test occursin("ThreadedRejection(8)", sprint(show, TR(8)))
end

@testset "batching does not change the answer" begin
    # Batched and serial rejection are statistically equivalent, not bit-identical: the batch
    # consumes draws past the accepted point, so the rng streams diverge. The evidence is
    # therefore compared at its own uncertainty, and the best point — which both find on the
    # same χ² surface — much more tightly.
    #
    # On a single-threaded runner both calls take the same path and this asserts nothing
    # beyond `batch` being accepted, which is itself worth knowing.
    common = (; backend = :nestedsamplers, lb = nested_lb, ub = nested_ub,
                nactive = 150, verb = false, cornerplot = false)
    r1 = fit_model_nested(nested_md, nested_free, nested_data;
                          batch = 1, rng = Xoshiro(7), common...)
    r8 = fit_model_nested(nested_md, nested_free, nested_data;
                          batch = 8, rng = Xoshiro(7), common...)

    @test r1.x_opt ≈ r8.x_opt rtol = 1e-3
    @test r1.chi2  ≈ r8.chi2  rtol = 1e-4
    @test abs(r1.logz - r8.logz) <= 3 * sqrt(r1.logzerr^2 + r8.logzerr^2)
    # The max-likelihood point survives the per-thread slots being reduced correctly: a lost
    # update would show up here as a visibly worse chi2 from the threaded run.
    @test r8.chi2 <= r1.chi2 * (1 + 1e-3)
end
end

# The comparison. Deliberately a small run: enough live points for the two to agree on a
# one-parameter problem, few enough that the suite does not grow a minute.
if have(:nestedsamplers) && have(:ultranest)
    @testset "the two samplers agree" begin
        rj = fit_model_nested(nested_md, nested_free, nested_data;
                              backend = :nestedsamplers, lb = nested_lb, ub = nested_ub,
                              nactive = 150, verb = false, cornerplot = false,
                              rng = Xoshiro(42))
        ru = fit_model_nested(nested_md, nested_free, nested_data;
                              backend = :ultranest, lb = nested_lb, ub = nested_ub,
                              min_num_live_points = 150, verb = false, cornerplot = false)

        @test rj.backend === :nestedsamplers
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

if have(:nestedsamplers)
    @testset "FlatModel form and bounds" begin
        fm = dict_to_model(nested_md, nested_free)
        r  = fit_model_nested(fm, nested_data; backend = :nestedsamplers,
                              lb = nested_lb, ub = nested_ub, nactive = 100,
                              verb = false, cornerplot = false, rng = Xoshiro(7))
        @test r.backend === :nestedsamplers
        @test nested_lb["s,ud"] <= r.x_opt[1] <= nested_ub["s,ud"]

        # The bounds ARE the prior, so a missing one is an error and not a default. Silently
        # inventing a range would decide the answer.
        @test_throws ErrorException fit_model_nested(nested_md, nested_free, nested_data;
                                                     backend = :nestedsamplers,
                                                     ub = nested_ub, verb = false)
        @test_throws ErrorException fit_model_nested(nested_md, nested_free, nested_data;
                                                     backend = :nestedsamplers,
                                                     lb = nested_lb, verb = false)
    end
end

end   # nested sampling
