using OITOOLS, Test, Random
using PythonPlot           # live-monitor testset: pyplot.switch_backend, plotclose

# ─────────────────────────────────────────────────────────────────────────────
# SQUEEZE MCMC image reconstruction.
#
# The critical test is the rank-1 one: the entire performance premise of the
# sampler is that a maintained visibility vector updated by subtract-add stays
# equal to a full recomputation.  If that fails, nothing else about the method
# holds.
# ─────────────────────────────────────────────────────────────────────────────

const SQ = OITOOLS
const SQ_DATAFILE = joinpath(@__DIR__, "oifits_for_tests", "2004-data1.oifits")

@testset "squeeze" begin

if !isfile(SQ_DATAFILE)
    @warn "2004-data1.oifits not found — skipping squeeze tests"
else

nx, pixsize = 64, 0.2
weights = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]

@testset "cvis_to_chi2_noalloc matches cvis_to_chi2_fg" begin
    for T in (Float64, Float32)
        data = readoifits(SQ_DATAFILE; T=T, verbose=false, warn=false)[1,1]
        rng = Xoshiro(42)
        V = Complex{T}.(randn(rng, T, data.nuv) .* T(0.3) .+ 1,
                        randn(rng, T, data.nuv) .* T(0.3))
        for w in ([1.0,1.0,1.0,0,0,0,0], [1.0,0,0,0,0,0,0], [0,1.0,1.0,0,0,0,0])
            ref = cvis_to_chi2_fg(V, data; weights=w)[1]
            @test cvis_to_chi2_noalloc(V, data; weights=w) ≈ ref rtol=(T === Float64 ? 1e-12 : 1e-5)
        end
        # The forward-only path must not allocate: it runs once per proposal, and a
        # single nuv-length temporary would cost thousands of bytes.  Bounded rather
        # than == 0 so a different inlining decision on another platform/Julia
        # version is not a failure -- anything real is orders of magnitude above 64.
        wt = SQ._weights_tuple(weights)
        SQ._cvis_to_chi2(V, data, wt, false)
        @test (@allocated SQ._cvis_to_chi2(V, data, wt, false)) < 64
    end
end

@testset "pixel convention matches image_to_vis" begin
    # histogram[ix,iy] must transform exactly as OITOOLS' own image path, or every
    # reconstruction is silently transposed/flipped.
    for T in (Float64, Float32)
        data = readoifits(SQ_DATAFILE; T=T, verbose=false, warn=false)[1,1]
        A = setup_dft(data, nx, pixsize)
        s = SQ.SqueezeState(T, nx, 500, data.nuv)
        SQ.init_elements!(s, Xoshiro(7); x_start=:random)
        SQ.resync_raw_im_vis!(s, A)
        ref = image_to_vis(T.(s.histogram), A)          # normalises by sum(x)
        @test maximum(abs, s.raw_im_vis ./ T(500) .- ref) <
              (T === Float64 ? 1e-12 : 1e-5) * maximum(abs, ref)
    end
end

@testset "rank-1 update stays equal to a full recompute" begin
    for T in (Float64, Float32)
        data = readoifits(SQ_DATAFILE; T=T, verbose=false, warn=false)[1,1]
        A = setup_dft(data, nx, pixsize)
        nelements = 2000
        for nmoves in (100, 10_000, 50_000)
            s = SQ.SqueezeState(T, nx, nelements, data.nuv)
            rng = Xoshiro(7)
            SQ.init_elements!(s, rng; x_start=:random)
            SQ.resync_raw_im_vis!(s, A)
            for _ in 1:nmoves
                i = rand(rng, 1:nelements)
                ixo, iyo, po = Int(s.element_x[i]), Int(s.element_y[i]), Int(s.element_p[i])
                ixn = SQ._wrap(ixo + rand(rng, -2:2), nx)
                iyn = SQ._wrap(iyo + rand(rng, -2:2), nx)
                pn = SQ._pixel_index(ixn, iyn, nx)
                pn == po && continue
                SQ.rank1_update!(s.new_raw_im_vis, s.raw_im_vis, A, pn, po)
                s.raw_im_vis, s.new_raw_im_vis = s.new_raw_im_vis, s.raw_im_vis
                s.histogram[ixo, iyo] -= Int32(1); s.histogram[ixn, iyn] += Int32(1)
                s.element_x[i] = ixn; s.element_y[i] = iyn; s.element_p[i] = pn
            end
            ref = A * Complex{T}.(vec(s.histogram))
            # drift after N accepted moves is bounded by ~eps(T)*sqrt(2N)*nelements
            tol = 10 * eps(T) * sqrt(2 * nmoves) * nelements
            @test maximum(abs, s.raw_im_vis .- ref) < tol
            @test sum(s.histogram) == nelements
        end
        # the hot loop must not allocate
        s = SQ.SqueezeState(T, nx, 500, data.nuv)
        SQ.init_elements!(s, Xoshiro(1); x_start=:random)
        SQ.rank1_update!(s.new_raw_im_vis, s.raw_im_vis, A, 5, 6)
        @test (@allocated SQ.rank1_update!(s.new_raw_im_vis, s.raw_im_vis, A, 5, 6)) < 64
        SQ.form_mod_vis!(s.mod_vis, s.raw_im_vis, s)
        @test (@allocated SQ.form_mod_vis!(s.mod_vis, s.raw_im_vis, s)) < 64
    end
end

@testset "starting configurations" begin
    data = readoifits(SQ_DATAFILE; T=Float64, verbose=false, warn=false)[1,1]
    ft = setup_dft(data, nx, pixsize)
    ne = 500

    # :point_source — every quantum on the centre pixel, matching C's default
    # (squeeze.c:496-501).  This is NOT a random or flat start.
    s = SQ.SqueezeState(Float64, nx, ne, data.nuv)
    SQ.init_elements!(s, Xoshiro(1); x_start=:point_source)
    c = nx ÷ 2 + 1
    @test s.histogram[c, c] == ne
    @test sum(s.histogram) == ne
    @test count(>(0), s.histogram) == 1

    # omitting x_start entirely must give exactly that
    s2 = SQ.SqueezeState(Float64, nx, ne, data.nuv)
    SQ.init_elements!(s2, Xoshiro(1))
    @test s2.histogram == s.histogram
    @test reconstruct_squeeze(data, ft; niter=25, nchains=1, seed=4, verb=false)[1] ==
          reconstruct_squeeze(:point_source, data, ft; niter=25, nchains=1, seed=4, verb=false)[1]

    # :random — spread out, and strings are accepted for the symbols
    s3 = SQ.SqueezeState(Float64, nx, ne, data.nuv)
    SQ.init_elements!(s3, Xoshiro(1); x_start=:random)
    @test sum(s3.histogram) == ne
    @test count(>(0), s3.histogram) > 100
    s4 = SQ.SqueezeState(Float64, nx, ne, data.nuv)
    SQ.init_elements!(s4, Xoshiro(1); x_start="random")
    @test s4.histogram == s3.histogram

    # an image is digitised (C's -i): deterministic, flux-proportional
    img = zeros(Float64, nx, nx); img[10, 20] = 3.0; img[40, 50] = 1.0
    s5 = SQ.SqueezeState(Float64, nx, ne, data.nuv)
    SQ.init_elements!(s5, Xoshiro(1); x_start=img)
    @test sum(s5.histogram) == ne
    @test count(>(0), s5.histogram) == 2
    @test s5.histogram[10, 20] > s5.histogram[40, 50]        # 3:1 flux ratio
    @test isapprox(s5.histogram[10, 20] / ne, 0.75; atol=0.01)
    # deterministic: the rng must not matter for an image start
    s6 = SQ.SqueezeState(Float64, nx, ne, data.nuv)
    SQ.init_elements!(s6, Xoshiro(999); x_start=img)
    @test s6.histogram == s5.histogram

    @test_throws ArgumentError SQ.init_elements!(s, Xoshiro(1); x_start=:nosuchstart)
    @test_throws ArgumentError SQ.init_elements!(s, Xoshiro(1); x_start=fill(-1.0, nx, nx))
    @test_throws ArgumentError SQ.init_elements!(s, Xoshiro(1); x_start=zeros(nx, nx))
    @test_throws DimensionMismatch SQ.init_elements!(s, Xoshiro(1); x_start=zeros(8, 8))
end

@testset "prior image (mask)" begin
    data = readoifits(SQ_DATAFILE; T=Float64, verbose=false, warn=false)[1,1]
    ft = setup_dft(data, nx, pixsize)
    c = (nx + 1) / 2

    # A zero-prior pixel gets a 1e12 penalty, so no quantum can ever sit there:
    # this is what makes the prior act as a hard mask (C's -p).
    mask = [sqrt((i-c)^2 + (j-c)^2) <= 20 ? 1.0 : 0.0 for i in 1:nx, j in 1:nx]
    img, d = reconstruct_squeeze(data, ft; niter=200, nchains=1, seed=1,
                                 prior_image=mask, verb=false)
    @test sum(img[mask .== 0]) == 0.0            # exactly zero, not just small
    @test sum(img) ≈ 1.0 rtol=1e-12
    @test d.regs_final[1].priorimage == 0.0      # inside the mask log(1)=0
    @test d.regs_final[1].centering == 0.0       # a prior disables auto-centering

    # Tightening the mask past the source must degrade the fit monotonically.
    chi2 = map((12, 20, 30)) do R
        m = [sqrt((i-c)^2 + (j-c)^2) <= R ? 1.0 : 0.0 for i in 1:nx, j in 1:nx]
        minimum(reconstruct_squeeze(data, ft; niter=250, nchains=1, seed=1,
                                    prior_image=m, verb=false)[2].chi2r_mean)
    end
    @test chi2[1] > chi2[2] > chi2[3]

    # penalty map: -log(p) inside, 1e12 where p <= 0
    soft = fill(0.5, nx, nx); soft[1, 1] = 0.0
    _, ds = reconstruct_squeeze(data, ft; niter=20, nchains=1, seed=1,
                                prior_image=soft, verb=false)
    # every quantum sits on a p=0.5 pixel => total penalty = nelements * log(2)
    @test ds.regs_final[1].priorimage ≈ ds.nelements * log(2) rtol=1e-10

    @test_throws DimensionMismatch reconstruct_squeeze(data, ft; niter=2,
                                       prior_image=zeros(8, 8), verb=false)
    _, dw = reconstruct_squeeze(data, ft; niter=20, nchains=1, seed=1,
                                prior_image=fill(0.5, nx, nx),
                                regularizers=[["priorimage", 2.0]], verb=false)
    @test dw.regs_final[1].priorimage ≈ dw.nelements * log(2) rtol=1e-10
    @test_throws ArgumentError reconstruct_squeeze(data, ft;
                                   regularizers=[["priorimage", 1.0], ["nope", 1.0]], verb=false)
end

@testset "regularizers" begin
    data = readoifits(SQ_DATAFILE; T=Float64, verbose=false, warn=false)[1,1]
    ft = setup_dft(data, nx, pixsize)
    h = zeros(Int32, 8, 8)
    @test SQ.reg_l0(h) == 0
    h[3,3] = 2; h[4,4] = 1
    @test SQ.reg_l0(h) == 2
    @test SQ.reg_entropy(h) ≈ SQ._ent(2) + SQ._ent(1)     # lgamma(2)+lgamma(1) = 0
    @test SQ.reg_compactness(h) > 0
    @test SQ.reg_tv(h, 3, 0.0) > 0

    # a delta-function at the centre must be more compact than one at a corner
    a = zeros(Int32, 8, 8); a[4,4] = 10
    b = zeros(Int32, 8, 8); b[1,1] = 10
    @test SQ.reg_compactness(a) < SQ.reg_compactness(b)

    # λ = 0 must be bit-identical to omitting the regularizer entirely
    base = reconstruct_squeeze(data, ft; niter=30, nchains=1, seed=3,
                               auto_centering=false, verb=false)[1]
    for nm in ("l0", "tv", "entropy", "compactness")
        z = reconstruct_squeeze(data, ft; niter=30, nchains=1, seed=3,
                                regularizers=[[nm, 0.0]], auto_centering=false, verb=false)[1]
        @test z == base
    end

    @test_throws ArgumentError reconstruct_squeeze(data, ft;
                                   regularizers=[["nosuchreg", 1.0]], verb=false)
end

@testset "live monitoring" begin
    # matplotlib is already imported by the time this runs (OITOOLS loads oiplot.jl
    # at package load), so setting ENV["MPLBACKEND"] here would do nothing — the
    # backend has to be switched at runtime, and restored afterwards so an
    # interactive session is not left headless.
    _oldbackend = string(pyplot.get_backend())
    pyplot.switch_backend("Agg")
    data = readoifits(SQ_DATAFILE; T=Float64, verbose=false, warn=false)[1,1]
    ft = setup_dft(data, nx, pixsize)

    # monitoring must not change the answer, only draw it
    a, _ = reconstruct_squeeze(data, ft; niter=60, nchains=1, seed=5, verb=false)
    b, _ = reconstruct_squeeze(data, ft; niter=60, nchains=1, seed=5, pixsize=pixsize,
                               monitor=20, verb=false)
    @test a == b

    m = SqueezeSparco(f_star=0.3, lambda0=Float64(data.uv_lam[1]),
                      free=(:f_star,), stepsize=[0.01,0,0,0,0,0])
    img, _ = reconstruct_squeeze(data, ft; niter=40, nchains=1, seed=5, pixsize=pixsize,
                                 monitor=20, model=m, verb=false)
    @test sum(img) ≈ 1.0 rtol=1e-12

    mon = SQ.SqueezeMonitor(10; pixsize=pixsize, title="unit test", free=[1])
    mon.ndf = 100.0
    SQ.monitor_update!(mon, fill(1/(nx*nx), nx, nx), 10;
                       chi2r=1.5, temperature=2.0, params=[0.5,0,0,1.6e-6,0,0])
    @test mon.iters == [10]
    @test mon.chi2r == [1.5]
    @test length(mon.params) == 1

    plotclose("all")
    pyplot.switch_backend(_oldbackend)
end

@testset "geometry comes from ft, as in reconstruct_bsmem" begin
    datam = readoifits(SQ_DATAFILE; T=Float64, verbose=false, warn=false)  # Matrix{OIdata}
    data  = datam[1,1]
    dft   = setup_dft(data, nx, pixsize)

    A, nxr, psr = SQ._squeeze_operator(dft, data, -1.0)
    @test A === dft && nxr == nx

    A2, nx2, _ = SQ._squeeze_operator(setup_ft(datam, nx, pixsize; mode="dft"), data, -1.0)
    @test nx2 == nx && size(A2) == size(dft)

    # the whole house pattern must work end to end
    img, _ = reconstruct_squeeze(datam, setup_ft(datam, nx, pixsize; mode="dft");
                                 niter=20, nchains=1, seed=1, verb=false)
    @test sum(img) ≈ 1.0 rtol=1e-12

    # setup_nfft returns SIX plans; ft[1] is the full-uv one
    nfft = setup_ft(datam, nx, pixsize; mode="nfft")
    A3, nx3, ps3 = SQ._squeeze_operator(nfft, data, -1.0)
    @test nx3 == nx
    @test ps3 ≈ pixsize rtol=1e-6
    @test size(A3) == (data.nuv, nx*nx)
    @test maximum(abs, A3 .- dft) < 1e-10

    d32 = readoifits(SQ_DATAFILE; T=Float32, verbose=false, warn=false)[1,1]
    @test_throws ArgumentError SQ._squeeze_operator(dft, d32, -1.0)            # wrong eltype
    @test_throws DimensionMismatch SQ._squeeze_operator(dft[1:end-1, :], data, -1.0)
    @test_throws ArgumentError SQ._squeeze_operator(dft[:, 1:end-1], data, -1.0)
    @test_throws ArgumentError SQ._squeeze_operator("not an operator", data, -1.0)
end

@testset "end-to-end reconstruction" begin
    data = readoifits(SQ_DATAFILE; T=Float64, verbose=false, warn=false)[1,1]
    ft = setup_dft(data, nx, pixsize)
    img, d = reconstruct_squeeze(data, ft; niter=250, nchains=2, seed=1, verb=false)
    @test size(img) == (nx, nx)
    @test all(>=(0), img)                 # positivity, by construction
    @test sum(img) ≈ 1.0 rtol=1e-12       # fixed total flux, by construction
    @test minimum(d.chi2r_mean) < 1.5     # the 2004 data reaches chi2r ~ 0.75
    @test d.best_chain == argmin(d.chi2r_mean)
    @test all(0 .< d.acceptance .< 1)
    @test all(d.nsamples .> 0)
end

@testset "reproducibility" begin
    data = readoifits(SQ_DATAFILE; T=Float64, verbose=false, warn=false)[1,1]
    ft = setup_dft(data, nx, pixsize)
    a, ga = reconstruct_squeeze(data, ft; niter=40, nchains=3, seed=7, verb=false)
    b, gb = reconstruct_squeeze(data, ft; niter=40, nchains=3, seed=7, verb=false)
    c, _  = reconstruct_squeeze(data, ft; niter=40, nchains=3, seed=8, verb=false)
    @test a == b                          # identical regardless of thread scheduling
    @test ga.chi2r_mean == gb.chi2r_mean
    @test a != c
end

@testset "precision is a type parameter, not a constant" begin
    for T in (Float64, Float32)
        data = readoifits(SQ_DATAFILE; T=T, verbose=false, warn=false)[1,1]
        img, d = reconstruct_squeeze(data, setup_dft(data, nx, pixsize); niter=250,
                                     nchains=1, seed=1, verb=false)
        @test eltype(setup_dft(data, nx, pixsize)) === Complex{T}
        @test sum(img) ≈ 1.0 rtol=1e-10
        @test d.chi2r_mean[1] < 1.5
    end
    @test SQ.default_resync(Float64) > 1e18
    @test 1e3 < SQ.default_resync(Float32) < 1e5
end

@testset "sparco model" begin
    data = readoifits(SQ_DATAFILE; T=Float64, verbose=false, warn=false)[1,1]
    pv = zeros(ComplexF64, data.nuv); fr = zeros(Float64, data.nuv)
    lam0 = Float64(data.uv_lam[1])

    # point source, no background: the star keeps f, the image keeps 1-f
    lp = SQ.sparco_vis!(pv, fr, [0.4, 0.0, 0.0, lam0, 0.0, 0.0], data)
    @test lp == 0.0
    @test all(v -> real(v) ≈ 0.4, pv)
    @test all(v -> imag(v) == 0, pv)
    @test all(v -> v ≈ 0.6, fr)
    @test all(i -> real(pv[i]) + fr[i] ≈ 1.0, 1:data.nuv)

    # a uniform-disc star must equal OITOOLS' own closed form
    SQ.sparco_vis!(pv, fr, [1.0, 1.5, 0.0, lam0, 0.0, 0.0], data)
    @test maximum(abs, pv .- visibility_ud([1.5], data.uv)) < 1e-9

    SQ.sparco_vis!(pv, fr, [0.0, 0.0, 0.0, lam0, 0.0, 0.0], data)
    @test all(iszero, pv) && all(≈(1.0), fr)

    # an over-resolved background dilutes without contributing visibility
    SQ.sparco_vis!(pv, fr, [0.0, 0.0, 0.0, lam0, 0.25, 0.0], data)
    @test all(iszero, pv)
    @test all(v -> v ≈ 0.75, fr)

    for bad in ([-0.1,0,0,lam0,0,0], [0.5,-1,0,lam0,0,0], [0.5,0,6.0,lam0,0,0],
                [0.5,0,0,lam0,-0.1,0], [0.5,0,0,lam0,0,-6.0])
        @test SQ.sparco_vis!(pv, fr, bad, data) >= 1e99
    end
    @test SQ.sparco_vis!(pv, fr, [0.5,1.0,2.0,lam0,0.1,1.0], data) == 0.0

    @test_throws ArgumentError SqueezeSparco(lambda0 = 0.0)
    @test_throws ArgumentError SqueezeSparco(lambda0 = lam0, free = (:nosuch,))
    @test_throws ArgumentError SqueezeSparco(lambda0 = lam0, free = (:f_star, :f_star))
    @test_throws ArgumentError SqueezeSparco(lambda0 = lam0, free = (99,))
    @test_throws ArgumentError SqueezeSparco(lambda0 = lam0, free = (:ud,),
                                             stepsize = [0.01,0,0,0,0,0])
    m = SqueezeSparco(lambda0 = lam0, free = (:f_star, :ud))
    @test m.free == [1, 2]
    @test SqueezeSparco(lambda0 = lam0, free = (2, 1)).free == [2, 1]
end

@testset "sparco recovers known parameters (chromatic)" begin
    # V1295 Aql (HD 190073) — SQUEEZE's own SPARCO example.  Reference setup:
    #   lambda0 = 1.6 um, f_star = 0.48, f_bg = 0, env_indx = -0.1,
    #   weights = [1, 0, 1] (this dataset's T3amp is noisy).
    # SPARCO's chromatism is per-uv-point via uv_lam and the image is grey, so ONE
    # OIdata spanning several wavelengths is all that is needed.  On MONOCHROMATIC
    # data a central point star is degenerate with a central image component and
    # f_star is NOT identifiable.
    polyfile   = joinpath(@__DIR__, "oifits_for_tests", "2019_v1295Aql.WL_SMOOTH.A.oifits")
    v1295_lam0 = 1.6e-6
    v1295_w    = [1.0, 0.0, 1.0]
    if !isfile(polyfile)
        @warn "2019_v1295Aql.WL_SMOOTH.A.oifits not found — skipping SPARCO recovery tests"
    else
        nxs, pss = 32, 0.25          # coarser than the demo's 64 @ 0.125, for speed
        data = readoifits(polyfile; T=Float64, filter_bad_data=true,
                          verbose=false, warn=false)[1,1]
        @test length(unique(data.uv_lam)) > 1      # the whole point
        A = setup_dft(data, nxs, pss)
        lam0 = v1295_lam0

        c = (nxs + 1) / 2
        img_true = [exp(-((sqrt((i-c)^2+(j-c)^2) - 4.0)^2)/4.0) for i in 1:nxs, j in 1:nxs]
        img_true ./= sum(img_true)

        function synth(params)
            pv = zeros(ComplexF64, data.nuv); fr = zeros(Float64, data.nuv)
            SQ.sparco_vis!(pv, fr, params, data)
            Vtot = pv .+ image_to_vis(img_true, A) .* fr
            syn = deepcopy(data)
            syn.v2 = vis_to_v2(Vtot, data.indx_v2)
            _, t3a, t3p = vis_to_t3(Vtot, data.indx_t3_1, data.indx_t3_2, data.indx_t3_3)
            syn.t3amp = t3a; syn.t3phi = t3p
            syn.v2_err = fill(0.01, length(syn.v2))
            syn.t3amp_err = fill(0.01, length(syn.t3amp))
            syn.t3phi_err = fill(1.0, length(syn.t3phi))
            return syn, Vtot
        end

        syn, Vtot = synth([0.48, 0.0, -0.1, lam0, 0.0, 0.0])
        @test cvis_to_chi2_noalloc(Vtot, syn; weights=v1295_w) < 1e-16

        for start in (0.25, 0.70)
            m = SqueezeSparco(f_star=start, env_indx=-0.1, lambda0=lam0,
                              free=(:f_star,), stepsize=[0.01,0,0,0,0,0])
            _, d = reconstruct_squeeze(syn, A; niter=400, nchains=2, seed=2,
                                       model=m, weights=v1295_w, verb=false)
            @test all(p -> abs(p[1] - 0.48) < 0.06, d.params)
        end

        m = SqueezeSparco(f_star=0.25, env_indx=1.0, lambda0=lam0,
                          free=(:f_star, :env_indx), stepsize=[0.01,0,0.05,0,0,0])
        _, d = reconstruct_squeeze(syn, A; niter=500, nchains=2, seed=5,
                                   model=m, weights=v1295_w, verb=false)
        b = d.best_chain
        @test abs(d.params[b][1] - 0.48) < 0.06
        @test abs(d.params[b][3] - (-0.1)) < 0.50

        # with chromatic leverage the chi2 profile minimises AT the truth --
        # on monochromatic data it does not (the star/image degeneracy)
        prof = map((0.25, 0.48, 0.75)) do fs
            m = SqueezeSparco(f_star=fs, env_indx=-0.1, lambda0=lam0, free=(),
                              stepsize=[0.01,0,0,0,0,0])
            minimum(reconstruct_squeeze(syn, A; niter=300, nchains=2, seed=3,
                                        model=m, weights=v1295_w, verb=false)[2].chi2r_mean)
        end
        @test prof[2] < prof[1] && prof[2] < prof[3]

        # a model disables auto-centering, as C does (squeeze.c:567-573)
        m = SqueezeSparco(f_star=0.48, lambda0=lam0, free=(), stepsize=[0.01,0,0,0,0,0])
        _, dm = reconstruct_squeeze(syn, A; niter=20, nchains=1, seed=1,
                                    model=m, weights=v1295_w, verb=false)
        @test dm.regs_final[1].centering == 0.0
        m2 = SqueezeSparco(f_star=0.48, lambda0=lam0, free=(:f_star,), stepsize=[0.01,0,0,0,0,0])
        _, dm2 = reconstruct_squeeze(syn, A; niter=20, nchains=1, seed=1,
                                     model=m2, weights=v1295_w, verb=false)
        @test dm2.ndf == dm.ndf + 1                 # free params enter ndf
        @test dm.params[1] == [0.48, 0.0, 0.0, lam0, 0.0, 0.0]   # fixed params untouched

        # chains must not share model state
        m3 = SqueezeSparco(f_star=0.30, lambda0=lam0, free=(:f_star,), stepsize=[0.01,0,0,0,0,0])
        _, d3 = reconstruct_squeeze(syn, A; niter=60, nchains=3, seed=9,
                                    model=m3, weights=v1295_w, verb=false)
        @test m3.params[1] == 0.30                 # the caller's model is not mutated
        @test length(unique(p[1] for p in d3.params)) > 1
    end
end

@testset "warm start and argument validation" begin
    data = readoifits(SQ_DATAFILE; T=Float64, verbose=false, warn=false)[1,1]
    ft = setup_dft(data, nx, pixsize)
    img, _ = reconstruct_squeeze(fill(1.0, nx, nx), data, ft; niter=30, nchains=1,
                                 seed=2, verb=false)
    @test sum(img) ≈ 1.0 rtol=1e-12
    @test_throws ArgumentError reconstruct_squeeze(fill(-1.0, nx, nx), data, ft;
                                                   niter=2, verb=false)
    @test_throws DimensionMismatch reconstruct_squeeze(fill(1.0, 32, 32), data, ft;
                                                       niter=2, verb=false)
end

end # isfile
end # testset
