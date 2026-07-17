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
    # Float64 methods must forward to the original calls, or the bit-gate breaks.
    x = randn(1000); y = randn(1000)
    @test _MEM._dot64(x, y) === dot(x, y)
    @test _MEM._sum64(x) === sum(x)
    # ...and the Float32 methods must beat naive Float32 accumulation by orders.
    x64 = randn(200_000); x32 = Float32.(x64)
    exact = dot(x64, x64)
    @test abs(_MEM._dot64(x32, x32) - exact) / exact < 1e-9   # Float64 accumulation
    @test abs(Float64(dot(x32, x32)) - exact) / exact > 1e-8  # naive Float32, for contrast
end

@testset "precision flows to plan, scratch and state" begin
    data = readoifits(joinpath(@__DIR__, "oifits_for_tests", "2004-data1.oifits");
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
