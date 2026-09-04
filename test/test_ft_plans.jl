# The Fourier transform plan types: OIft, NFFTCell, DFTCell.
#
# The assertion that matters is §"named plans match their index sets": each NFFT plan is built
# from one observable's uv indices, so the number of points it holds is an independent witness
# that a field's NAME matches the transform it holds. Checking the field order against
# `setup_nfft` would only prove the code agrees with itself.
#
# This is the guard for the criterion code reading `ftplan.v2` instead of `ftplan[3]`. The
# bit-for-bit BSMEM gate would also catch a mis-mapping, but it self-skips whenever the Julia
# version, thread count or CPU differs from the captured baseline, so it cannot be relied on
# to be running.

using OITOOLS, Test, LinearAlgebra

@testset "Fourier transform plans" begin
    data = readoifits(joinpath(@__DIR__, "..", "demos", "data", "BC2004", "2004-data1.oifits");
                      warn = false, verbose = false)
    d = data[1, 1]

    @testset "named plans match their index sets" begin
        c = setup_nfft(d, 64, 0.2)
        @test size(c.uv.k,   2) == d.nuv
        @test size(c.vis.k,  2) == length(d.indx_vis)
        @test size(c.v2.k,   2) == length(d.indx_v2)
        @test size(c.t3_1.k, 2) == length(d.indx_t3_1)
        @test size(c.t3_2.k, 2) == length(d.indx_t3_2)
        @test size(c.t3_3.k, 2) == length(d.indx_t3_3)
    end

    @testset "positional access is the same object as named" begin
        c = setup_nfft(d, 64, 0.2)
        for (i, f) in enumerate((:uv, :vis, :v2, :t3_1, :t3_2, :t3_3))
            @test c[i] === getfield(c, f)
        end
        @test length(c) == 6
        @test_throws BoundsError c[7]
    end

    @testset "OIft records what was built" begin
        for (mode, sym) in (("nfft", :nfft), ("dft", :dft))
            ft = setup_ft(data, 32, 0.25; mode = mode)
            @testset "$mode" begin
                @test ft isa OIft
                @test ft.mode    == sym
                @test ft.nx      == 32
                @test ft.pixsize == 0.25
                @test size(ft)   == size(data)
                # indexes like the matrix it wraps, so every `ft::AbstractMatrix` signature
                # keeps matching
                @test ft isa AbstractMatrix
                @test ft[1] === ft[1, 1]
            end
        end
        @test_throws ErrorException setup_ft(data, 32, 0.25; mode = "wavelet")
    end

    @testset "cells subtype what the call signatures dispatch on" begin
        nf = setup_ft(data, 32, 0.25)[1, 1]
        df = setup_ft(data, 32, 0.25; mode = "dft")[1, 1]
        @test nf isa AbstractVector{<:OITOOLS.NFFT.NFFTPlan}   # image_to_vis dispatches on this
        @test df isa AbstractMatrix{<:Complex}                 # and on this
        @test size(df, 1) == d.nuv
        @test size(df, 2) == 32^2
    end

    @testset "a DFT cell refuses per-observable access" begin
        df = setup_ft(data, 32, 0.25; mode = "dft")[1, 1]
        # A DFT holds one kernel over every uv point; answering `.v2` with something would be
        # answering a question it cannot answer.
        for f in (:uv, :vis, :v2, :t3_1, :t3_2, :t3_3)
            @test_throws ErrorException getproperty(df, f)
        end
    end

    @testset "precision flows through" begin
        d64 = readoifits(joinpath(@__DIR__, "..", "demos", "data", "BC2004", "2004-data1.oifits");
                         warn = false, verbose = false, T = Float64)
        @test setup_ft(data, 16, 0.5) isa OIft{Float32}
        @test setup_ft(d64,  16, 0.5) isa OIft{Float64}
    end

    @testset "show says mode, geometry and channel count" begin
        ft  = setup_ft(data, 64, 0.3)
        out = sprint(show, MIME"text/plain"(), ft)
        @test occursin("OIft{Float32}", out)
        @test occursin("NFFT", out)
        @test occursin("64×64", out)
        @test occursin("0.3", out)
        @test occursin("uv points", out)
        @test occursin("$(d.nuv)", out)
        @test occursin("DFT", sprint(show, MIME"text/plain"(), setup_ft(data, 64, 0.3; mode="dft")))
    end

    @testset "displaying a foreign type does not dispatch into OITOOLS" begin
        # Defining `display` for a type this package does not own — `AbstractMatrix`, say —
        # would change how every matrix in the session prints, here and in any other package.
        # `OIdata` methods are fine and expected: that type IS ours.
        for x in (rand(2, 2), Matrix{Float64}(undef, 0, 0), [1, 2, 3])
            @test !occursin("OITOOLS", string(Base.which(display, (typeof(x),)).module))
        end
        @test (display(Matrix{Float64}(undef, 0, 0)); true)   # and must not throw
        # the summary belongs to `show` on our own type, which is the legitimate hook
        @test !isempty(sprint(show, MIME"text/plain"(), setup_ft(data, 16, 0.5)))
    end
end
