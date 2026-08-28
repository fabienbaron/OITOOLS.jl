# BSMEM/MaximENT regression check.
#
# Every case still runs, and produces an image that is finite, non-negative and carries
# positive flux, with a non-empty convergence trace. Those hold on any machine, and they catch
# the breakage that actually happens — a method that throws, a solver that returns NaN, a plan
# wired to the wrong observable.
#
# They are deliberately not assertions about the NUMBERS. MaximENT bisects on alpha
# (maximent.jl:409, :1149), so a 1-ulp shift flips a branch and diverges the whole trajectory:
# a loose tolerance means nothing here, and an exact one is only defined on the precise
# machine, Julia patch, thread count and `--check-bounds` setting it was cut under — which is
# never the machine running the tests. `BSMEM_CASES` earns its keep by driving every distinct
# path through the solver instead; see `bsmem_cases.jl` for which case covers what.

using Test

include(joinpath(@__DIR__, "bsmem_cases.jl"))

@testset "BSMEM regression" begin
    for case in BSMEM_CASES
        @testset "$(case.name)" begin
            got = run_bsmem_case(case; T = Float64)

            # Runs at all. `run_bsmem_case` swallows the exception and reports it here, so a
            # failure names what went wrong instead of aborting the file.
            @test got.ok || (@error("$(case.name) threw: $(get(got, :err, "?"))"); false)
            got.ok || return

            img = got.image
            @test all(isfinite, img)
            @test all(>=(0), img)      # MaxEnt works in log-space; a negative pixel is a bug
            @test sum(img) > 0
            @test !isempty(got.history)
        end
    end
end
