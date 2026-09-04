# Structural precision tests for the T-parameterised BSMEM/MaximENT.
#
# These are cheap and catch what the numeric gate cannot:
#
#  * NFFT does NOT throw on mismatched element types — `mul!(::Vector{ComplexF32},
#    ::NFFTPlan{Float32}, ::Matrix{ComplexF64})` silently converts. So a stray Float64
#    buffer is invisible to the type system and must be asserted explicitly.
#  * A Float64 scalar inside a Float32 broadcast widens every element without allocating,
#    which no tool reports. summarysize ratios catch the array-level version of that.

using Test, OITOOLS, NFFT, LinearAlgebra
using Random          # seeded RNG: see the reductions testset below

const _MEM = OITOOLS

@testset "MaximENTState constructor shape" begin
    # Distinct primes: any nhid/ndat mis-ordering in the positional constructor is
    # impossible to miss.
    for T in (Float64, Float32)
        s = _MEM.MaximENTState{T}(7, 11)
        for f in _MEM._MES_HID_FIELDS
            @test length(getfield(s, f)) == 7
            @test eltype(getfield(s, f)) === T
        end
        for f in _MEM._MES_DAT_FIELDS
            @test length(getfield(s, f)) == 11
            @test eltype(getfield(s, f)) === T
        end
        # ritz_* stay Float64: they belong to the evidence/Ritz island.
        @test eltype(s.ritz_λ) === Float64
        @test eltype(s.ritz_w) === Float64
        # Fails if a field is added without updating the groups above or the constructor.
        @test length(_MEM._MES_HID_FIELDS) + length(_MEM._MES_DAT_FIELDS) + 2 ==
              fieldcount(_MEM.MaximENTState)
    end
    @test eltype(_MEM.MaximENTState(7, 11).h) === Float64  # back-compat default
end

@testset "_tol reproduces the Skilling formula" begin
    # The historical inline expression, which _tol replaced. If a compiler ever contracts
    # `4/3 - 1` this fails loudly instead of the gate drifting silently.
    @test _MEM._tol(Float64) === abs((4.0 / 3.0 - 1.0) - 1.0 / 3.0)
    @test _MEM._tol(Float32) === Float64(abs((4.0f0 / 3.0f0 - 1.0f0) - 1.0f0 / 3.0f0))
    @test _MEM._tol(Float64) === 5.551115123125783e-17   # frozen at refactor time
    @test _MEM._tol() === _MEM._tol(Float64)
end

@testset "constants are identity at Float64" begin
    # The Float64 path must be unchanged *by construction*, not by testing.
    @test _MEM._eps(Float64)   === 1e-20
    @test _MEM._tiny(Float64)  === 1e-30
    @test _MEM._floor(Float64) === 1e-12
    # ...and representable (not subnormal) at Float32, so no recalibration is needed.
    for f in (_MEM._eps, _MEM._tiny, _MEM._floor)
        @test !issubnormal(f(Float32))
        @test f(Float32) != 0
    end
end

@testset "Float64-accumulating reductions" begin
    # Seeded. This used to call randn() unseeded and assert that naive Float32 accumulation
    # was *inaccurate enough* (relative error > 1e-8) to make the point. That threshold sits
    # inside the lower tail of the distribution: measured over 40 draws of 200_000 elements,
    # the error ranged 6.7e-9 to 3.4e-7 with a median of 1.1e-7, so the assertion failed on
    # about 8% of runs. A file feeding a bit-for-bit gate has no business being random.
    rng = Xoshiro(20260818)

    # Float64 methods must forward to the original calls, or the bit-gate breaks.
    x = randn(rng, 1000); y = randn(rng, 1000)
    @test _MEM._dot64(x, y) === dot(x, y)
    @test _MEM._sum64(x) === sum(x)

    # ...and the Float32 methods must beat naive Float32 accumulation by orders.
    #
    # The contrast is a hand-written sequential Float32 loop, not `dot(x32, x32)`. BLAS is
    # free to accumulate a Float32 dot product however it likes, and implementations differ:
    # measured on the same input, OpenBLAS on x86 gives a relative error of 1.07e-7 while
    # Apple's BLAS on aarch64 gives 2.90e-8. The claim being tested is about accumulation
    # WIDTH, so the comparison has to be against an accumulation this file controls.
    naive32(v) = (s = 0.0f0; for a in v; s += a * a; end; s)

    x64 = randn(rng, 200_000); x32 = Float32.(x64)
    exact  = dot(x64, x64)
    err_64 = abs(_MEM._dot64(x32, x32) - exact) / exact      # Float64 accumulation
    err_32 = abs(Float64(naive32(x32)) - exact) / exact      # Float32 accumulation
    @test err_64 < 1e-9
    # Compare the two against each other rather than against an absolute floor: "Float64
    # accumulation wins by orders of magnitude" is the actual claim, and stating it as a
    # ratio holds for any draw or vector length. The measured ratio is ~1.4e5, so the
    # threshold has three orders of headroom.
    @test err_32 > 100 * err_64
end

@testset "precision flows to plan, scratch and state" begin
    data = readoifits(joinpath(@__DIR__, "..", "demos", "data", "BC2004", "2004-data1.oifits");
                      verbose=false, warn=false, T=Float32)
    ft = setup_ft(data, 32, 0.4)
    @test OITOOLS.ft_eltype(ft) === Float32
    ctx, s, p = _MEM.maxent_setup(data[1,1], 32, 0.4, vec(zeros(Float64, 32*32) .+ 1/1024);
                                  methd=[4,1,1,2], T=Float32)
    # NFFT would not have thrown on a mismatch — assert instead.
    @test ctx.plan isa NFFT.NFFTPlan{Float32,2,1}
    for f in (:visi, :phasor, :_f_hat, :_dvisi)
        @test eltype(getfield(ctx, f)) === ComplexF32
    end
    for f in (_MEM._MES_HID_FIELDS..., _MEM._MES_DAT_FIELDS...)
        @test eltype(getfield(s, f)) === Float32
    end
    # Deliberate Float64 islands.
    @test p.tol === _MEM._tol(Float64)
    @test ctx.xyint isa Float64
end

# ─────────────────────────────────────────────────────────────────────────────
# Precision flows through the chi2 / gradient stack too.
#
# These guard the class of bug where a hardcoded Float64 buffer, literal or signature
# silently promotes a Float32 pipeline back to double — or, worse, fails to dispatch at all.
# The polychromatic-plus-flux combination is called out explicitly because it is the one
# path that needs *both* a Float32 model-flux vector and a Float32 data vector.
# ─────────────────────────────────────────────────────────────────────────────

@testset "precision flows through chi2 and gradient" begin
    f = joinpath(@__DIR__, "..", "demos", "data", "BC2004", "2004-data1.oifits")
    d32 = readoifits(f; T=Float32, filter_bad_data=false, verbose=false, warn=false)[1,1]
    d64 = readoifits(f; T=Float64, filter_bad_data=false, verbose=false, warn=false)[1,1]

    V32 = rand(ComplexF32, d32.nuv)
    V64 = ComplexF64.(V32)

    # cvis building blocks
    @test cvis_to_chi2_f(V32, d32) isa Float32
    @test cvis_to_chi2_f(V64, d64) isa Float64
    c32, g32 = cvis_to_chi2_fg(V32, d32)
    c64, g64 = cvis_to_chi2_fg(V64, d64)
    @test c32 isa Float32 && eltype(g32) === ComplexF32
    @test c64 isa Float64 && eltype(g64) === ComplexF64
    # f and fg must agree on the shared (flux-free) part, exactly
    @test cvis_to_chi2_f(V32, d32) == c32
    @test cvis_to_chi2_f(V64, d64) == c64
    # and the two precisions must agree to Float32 accuracy
    @test abs(c32 - c64) / abs(c64) < 1e-5

    # model path: eval_model, chi2_flat, chi2_flat_fg
    m = dict_to_model(Dict{String,Any}("s,ud"=>3.0, "s,f"=>1.0), ["s,ud"])
    @test eltype(eval_model(m, Float32[3.0], d32.uv)) === ComplexF32
    @test eltype(eval_model(m, Float64[3.0], d64.uv)) === ComplexF64
    # mixed inputs promote rather than silently truncating
    @test eltype(eval_model(m, Float64[3.0], d32.uv)) === ComplexF64

    # `wl` and `mjd` are NOT in the accumulator's type, which is promote_type(eltype(x),
    # float(eltype(uv))) and nothing else. `uv_mjd` is stored in T and is documented as coarse
    # (~5.6 min at Float32) for that reason: keeping it in T is a memory choice, not a
    # precision one, and a Float64 MJD handed to a Float32 call must not promote the result.
    @test eltype(eval_model(m, Float32[3.0], d32.uv;
                            wl = d32.uv_lam, mjd = d32.uv_mjd)) === ComplexF32
    @test eltype(eval_model(m, Float32[3.0], d32.uv;
                            wl = d32.uv_lam, mjd = Float64.(d32.uv_mjd))) === ComplexF32
    @test eltype(eval_model(m, Float32[3.0], d32.uv;
                            wl = Float64.(d32.uv_lam), mjd = d32.uv_mjd)) === ComplexF32

    @test OITOOLS.chi2_flat(m, Float32[3.0], d32) isa Float32
    @test OITOOLS.chi2_flat(m, Float64[3.0], d64) isa Float64
    fc32, fg32 = OITOOLS.chi2_flat_fg(m, Float32[3.0], d32)
    @test fc32 isa Float32 && eltype(fg32) === Float32

    # Hankel components must honour the requested precision, not force Float64.
    hm = dict_to_model(Dict{String,Any}("r,profile"=>"exp(-(\$R/2.0)^2)",
                                        "r,udout"=>10.0, "r,f"=>1.0), String[]; T=Float32)
    @test eltype(hm.components[1].workspace.V) === Float32
    @test eltype(hm.components[1].r_grid) === Float32
    @test eltype(eval_model(hm, Float32[], d32.uv)) === ComplexF32

    # Plot entry points must accept Float32 arrays (they took Array{Float64,N} before).
    @test hasmethod(uvplot, Tuple{Matrix{Float32}})
    @test hasmethod(imdisp_multi, Tuple{Array{Float32,3}})
end

@testset "polychromatic chi2 precision, including the flux term" begin
    f = joinpath(@__DIR__, "..", "demos", "data", "simulation-2004image.oifits")
    if isfile(f)
        m = dict_to_model(Dict{String,Any}("s,ud"=>3.0, "s,f"=>1.0), ["s,ud"])
        wf = [1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0]     # flux term ON
        p32 = vec(readoifits(f; T=Float32, filter_bad_data=false, polychromatic=true,
                             verbose=false, warn=false))
        p64 = vec(readoifits(f; T=Float64, filter_bad_data=false, polychromatic=true,
                             verbose=false, warn=false))
        # _chi2_flux_polychromatic used to demand Vector{Float64} and threw a MethodError here
        c32 = OITOOLS.chi2_flat(m, Float32[3.0], p32; weights=wf)
        c64 = OITOOLS.chi2_flat(m, Float64[3.0], p64; weights=wf)
        @test c32 isa Float32
        @test c64 isa Float64
        @test abs(c32 - c64) / abs(c64) < 1e-4

        gc32, gg32 = OITOOLS.chi2_flat_fg(m, Float32[3.0], p32; weights=wf)
        @test gc32 isa Float32 && eltype(gg32) === Float32
    end
end

@testset "sparco_flat precision" begin
    f = joinpath(@__DIR__, "..", "demos", "data", "BC2004", "2004-data1.oifits")
    spec = Dict{String,Any}("s,ud"=>3.0, "s,f"=>0.5, "W"=>1.0)
    chi2s = Dict{DataType,Any}()
    for T in (Float32, Float64)
        d  = readoifits(f; T=T, filter_bad_data=false, verbose=false, warn=false)[1,1]
        ft = setup_nfft(d, 32, 0.4)
        x_img = ones(T, 32, 32) ./ (32*32)
        m  = dict_to_model(spec, String[]; T=T)

        # The chromatic image weight W used to be forced to Float64, which promoted
        # everything downstream of it.
        W = OITOOLS._eval_W(m, T[], d.uv_lam, OITOOLS._resolve_w_idx(m, "W"))
        @test eltype(W) === T

        c = chi2_sparco_flat_f(x_img, T[], m, ft, d; verb=false)
        @test c isa T
        chi2s[T] = c

        xv = vcat(T[], vec(x_img)); g = similar(xv)
        c2 = OITOOLS.chi2_sparco_flat_fg(xv, g, m, ft, d, 0; verb=false)
        @test c2 isa T
        @test eltype(g) === T
    end
    @test abs(chi2s[Float32] - chi2s[Float64]) / abs(chi2s[Float64]) < 1e-5
end
