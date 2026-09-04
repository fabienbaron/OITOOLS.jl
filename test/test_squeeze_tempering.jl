using OITOOLS, Test, Random, Statistics

# ─────────────────────────────────────────────────────────────────────────────
# Parallel tempering for SQUEEZE, via the Pigeons.jl package extension.
#
# Pigeons is a WEAK dependency: it is not installed by CI, because it pulls in
# MPI.jl and that is a poor thing to build across four CI configurations for a
# feature most users will not touch.  These tests therefore skip themselves when
# Pigeons cannot be loaded.  Run them locally with:
#
#     julia --project=. -e 'using Pkg; Pkg.add("Pigeons")' && julia --project=. test/runtests.jl
# ─────────────────────────────────────────────────────────────────────────────

const _HAVE_PIGEONS = try
    @eval using Pigeons
    true
catch
    false
end

@testset "squeeze tempering (Pigeons)" begin

if !_HAVE_PIGEONS
    @info "Pigeons not available — skipping SQUEEZE tempering tests (weak dependency)"
    @test true
else

datafile = joinpath(@__DIR__, "..", "demos", "data", "BC2004", "2004-data1.oifits")
nx, pixsize = 64, 0.2
data = readoifits(datafile; T=Float64, verbose=false, warn=false)[1,1]
ft = setup_dft(data, nx, pixsize)

@testset "extension loads and the entry point is defined" begin
    @test Base.get_extension(OITOOLS, :OITOOLSPigeonsExt) !== nothing
    @test !isempty(methods(reconstruct_squeeze_tempered))
end

@testset "log density and reference" begin
    SQ = OITOOLS
    target = SQ.SqueezeTarget{Float64}(data, ft, nx, 500,
                                       SQ._weights_tuple([1.0,1.0,1.0]), false,
                                       Any[], Float64[], 1.0, 1.0, 0.0, :random,
                                       nothing, Float64[], Float64[])
    st = Pigeons.initialization(target, Xoshiro(1), 1)
    @test st isa SQ.SqueezePTState
    @test sum(st.s.histogram) == 500

    # the target log potential is -(lLikelihood + lPrior), and it must agree with
    # a from-scratch evaluation of the same state
    lp = target(st)
    @test lp ≈ -(st.s.lLikelihood + st.s.lPrior)
    SQ.squeeze_refresh!(st, target)
    @test target(st) ≈ lp rtol=1e-10

    # the reference is the uniform distribution: constant log density
    ref = Pigeons.default_reference(target)
    @test ref isa SQ.SqueezeUniformReference
    @test ref(st) == 0.0
end

@testset "the explorer is beta-invariant in the right limits" begin
    SQ = OITOOLS
    target = SQ.SqueezeTarget{Float64}(data, ft, nx, 400,
                                       SQ._weights_tuple([1.0,1.0,1.0]), false,
                                       Any[], Float64[], 1.0, 1.0, 0.0, :random,
                                       nothing, Float64[], Float64[])
    st = Pigeons.initialization(target, Xoshiro(3), 1)
    ex = SQ.SqueezeExplorer()

    # beta = 0: every proposal is accepted (log(u) < 0 always), so the chain is a
    # pure random walk over configurations -- the uniform reference.
    st.s.naccepted = 0; st.s.nproposed = 0
    SQ.squeeze_step!(st, target, ex, 0.0, Xoshiro(4))
    @test st.s.naccepted == st.s.nproposed

    # beta = 1: acceptance is well below 1 and flux is still conserved
    st.s.naccepted = 0; st.s.nproposed = 0
    SQ.squeeze_step!(st, target, ex, 1.0, Xoshiro(5))
    @test 0 < st.s.naccepted < st.s.nproposed
    @test sum(st.s.histogram) == 400

    # the maintained state must still match a full recompute after stepping
    ll_inc = st.s.lLikelihood
    SQ.squeeze_refresh!(st, target)
    @test st.s.lLikelihood ≈ ll_inc rtol=1e-8
end

@testset "end-to-end tempered reconstruction" begin
    img, d = reconstruct_squeeze_tempered(data, ft; n_rounds=9, n_chains=12,
                                          seed=1, verb=false)
    @test size(img) == (nx, nx)
    @test all(>=(0), img)                     # positivity, by construction
    @test sum(img) ≈ 1.0 rtol=1e-12           # fixed total flux
    @test d.chi2r_mean < 1.5                  # matches the annealed result (~0.75)
    @test isfinite(d.logZ)
    @test d.global_barrier > 0
    @test d.nsamples == 2^9                   # ONLY the final round is kept

    # the tempered posterior mean must agree with the annealed image in structure
    ann, _ = reconstruct_squeeze(data, ft; niter=300, nchains=2, seed=1, verb=false)
    @test cor(vec(img), vec(ann)) > 0.9
end

@testset "regularizers and prior mask under tempering" begin
    _, d1 = reconstruct_squeeze_tempered(data, ft; n_rounds=8, n_chains=10, seed=1,
                                         regularizers=[["l0", 1.0], ["tv", 200.0]],
                                         verb=false)
    @test d1.chi2r_mean < 2.0

    c = (nx + 1) / 2
    mask = [sqrt((i-c)^2 + (j-c)^2) <= 25 ? 1.0 : 0.0 for i in 1:nx, j in 1:nx]
    img, d2 = reconstruct_squeeze_tempered(data, ft; n_rounds=8, n_chains=10, seed=1,
                                           prior_image=mask, verb=false)
    @test sum(img[mask .== 0]) == 0.0         # the mask is hard here too
    @test d2.chi2r_mean < 2.0

    @test_throws ArgumentError reconstruct_squeeze_tempered(data, ft; n_rounds=2,
                                   n_chains=3, regularizers=[["nope", 1.0]], verb=false)
end

@testset "SPARCO under tempering" begin
    polyfile = joinpath(@__DIR__, "..", "demos", "data", "2019_v1295Aql.WL_SMOOTH.A.oifits")
    if !isfile(polyfile)
        @warn "V1295 Aql not found — skipping SPARCO tempering tests"
    else
        SQ = OITOOLS
        nxs, pss, lam0 = 32, 0.25, 1.6e-6
        w1295 = [1.0, 0.0, 1.0]
        pdata = readoifits(polyfile; T=Float64, filter_bad_data=true,
                           verbose=false, warn=false)[1,1]
        Ap = setup_dft(pdata, nxs, pss)
        cc = (nxs + 1) / 2
        img_true = [exp(-((sqrt((i-cc)^2+(j-cc)^2) - 4.0)^2)/4.0) for i in 1:nxs, j in 1:nxs]
        img_true ./= sum(img_true)

        function synth(params)
            pv = zeros(ComplexF64, pdata.nuv); fr = zeros(Float64, pdata.nuv)
            SQ.sparco_vis!(pv, fr, params, pdata)
            V = pv .+ image_to_vis(img_true, Ap) .* fr
            syn = deepcopy(pdata)
            syn.v2 = vis_to_v2(V, pdata.indx_v2)
            _, t3a, t3p = vis_to_t3(V, pdata.indx_t3_1, pdata.indx_t3_2, pdata.indx_t3_3)
            syn.t3amp = t3a; syn.t3phi = t3p
            syn.v2_err = fill(0.02, length(syn.v2))
            syn.t3amp_err = fill(0.02, length(syn.t3amp))
            syn.t3phi_err = fill(2.0, length(syn.t3phi))
            syn
        end

        # the default reference box must be proper (bounded) for every free parameter
        m0 = SqueezeSparco(f_star=0.48, lambda0=lam0, free=(:f_star, :ud))
        lo, hi = SQ.default_sparco_bounds(m0)
        @test all(isfinite, lo) && all(isfinite, hi)
        @test all(j -> lo[j] < hi[j], m0.free)

        syn = synth([0.48, 0.0, -0.1, lam0, 0.0, 0.0])

        # f_star must converge to the truth from either side, start-independently
        res = map((0.15, 0.80)) do start
            m = SqueezeSparco(f_star=start, env_indx=-0.1, lambda0=lam0,
                              free=(:f_star,), stepsize=[0.02,0,0,0,0,0])
            reconstruct_squeeze_tempered(syn, Ap; n_rounds=8, n_chains=12, seed=1,
                                         model=m, weights=w1295, verb=false)[2]
        end
        for d in res
            @test abs(d.params[1] - 0.48) < 0.06
            @test d.chi2r_mean < 2.0
            @test length(d.params) == 6
        end
        @test abs(res[1].params[1] - res[2].params[1]) < 0.03   # start-independent

        # A shape parameter that annealing cannot move (the image adapts to it, so
        # no step is acceptable) DOES travel under tempering, because at low beta
        # the likelihood is weak and swaps carry the information up to beta = 1.
        synud = synth([0.48, 2.0, -0.1, lam0, 0.0, 0.0])
        m = SqueezeSparco(f_star=0.48, ud=4.0, env_indx=-0.1, lambda0=lam0,
                          free=(:ud,), stepsize=[0,0.1,0,0,0,0])
        _, dt = reconstruct_squeeze_tempered(synud, Ap; n_rounds=8, n_chains=12,
                                             seed=1, model=m, weights=w1295, verb=false)
        _, da = reconstruct_squeeze(synud, Ap; niter=400, nchains=2, seed=1,
                                    model=SqueezeSparco(f_star=0.48, ud=4.0, env_indx=-0.1,
                                                        lambda0=lam0, free=(:ud,),
                                                        stepsize=[0,0.1,0,0,0,0]),
                                    weights=w1295, verb=false)
        @test abs(dt.params[2] - 2.0) < 0.6              # tempering reaches the truth
        @test abs(dt.params[2] - 2.0) < abs(da.params[da.best_chain][2] - 2.0)

        # a model disables auto-centering here too
        mf = SqueezeSparco(f_star=0.48, lambda0=lam0, free=(), stepsize=[0.02,0,0,0,0,0])
        _, dm = reconstruct_squeeze_tempered(syn, Ap; n_rounds=6, n_chains=6, seed=1,
                                             model=mf, weights=w1295, verb=false)
        # `params` is a posterior MEAN, so even fixed values pick up summation
        # round-off (0.48 summed over 2^n samples then divided).
        @test dm.params ≈ [0.48, 0.0, 0.0, lam0, 0.0, 0.0] rtol=1e-12

        # the reference box must be validated
        bad = SqueezeSparco(f_star=0.48, lambda0=lam0, free=(:f_star,), stepsize=[0.02,0,0,0,0,0])
        @test_throws ArgumentError reconstruct_squeeze_tempered(syn, Ap; n_rounds=3,
                          n_chains=4, model=bad, weights=w1295, verb=false,
                          param_bounds=([0.6,0,-5,lam0,0,-5], [1.0,20,5,lam0,1,5]))
        @test_throws ArgumentError reconstruct_squeeze_tempered(syn, Ap; n_rounds=3,
                          n_chains=4, model=bad, weights=w1295, verb=false,
                          param_bounds=([0.5,0,-5,lam0,0,-5], [0.5,20,5,lam0,1,5]))
    end
end

@testset "reproducibility" begin
    a, da = reconstruct_squeeze_tempered(data, ft; n_rounds=7, n_chains=8, seed=3,
                                         multithreaded=false, verb=false)
    b, db = reconstruct_squeeze_tempered(data, ft; n_rounds=7, n_chains=8, seed=3,
                                         multithreaded=false, verb=false)
    @test a == b
    @test da.logZ == db.logZ
end

end # _HAVE_PIGEONS
end # testset
