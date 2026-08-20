# BSMEM/MaximENT regression check.
#
# What this asserts: every case still runs, and produces an image that is finite, non-negative
# and carries positive flux. Those hold on any machine, and they catch the breakage that
# actually happens — a method that throws, a solver that returns NaN, a plan wired to the wrong
# observable.
#
# What it no longer asserts by default: bit-for-bit equality against a frozen baseline.
#
# The reasoning for the bit gate was sound — the solver bisects on alpha (maximent.jl:409,
# :1149), so a 1-ulp shift flips a branch and diverges the whole trajectory, which makes a
# loose tolerance meaningless. But the same chaos makes the comparison valid only on the exact
# machine, Julia patch, thread count and `--check-bounds` setting the baseline was captured
# under, and in practice that is never the machine running the tests. Measured here, all eight
# cases reported:
#
#     bit-for-bit gate skipped — codegen/threading mismatch, baseline not comparable —
#     check_bounds: 1 -> 0, nthreads: 1 -> 16, julia: 1.12.6 -> 1.12.7, cpu: <absent> -> znver5
#
# So it compared nothing, while carrying baselines in the repository and failing CI whenever
# the environment drifted. A check that is skipped everywhere it runs is not a safety net.
#
# It is kept behind a switch rather than deleted, because it is the right tool for one specific
# job: verifying that a pure-typing refactor changed no numbers. Re-capture and run it on ONE
# machine, deliberately:
#
#     julia --project=. test/capture_bsmem_baseline.jl
#     OITOOLS_BSMEM_BITGATE=1 julia --project=. test/runtests.jl

using Test, Serialization

include(joinpath(@__DIR__, "bsmem_cases.jl"))

"Opt in with `OITOOLS_BSMEM_BITGATE=1`, having captured a baseline on this machine first."
const BITGATE = get(ENV, "OITOOLS_BSMEM_BITGATE", "0") in ("1", "true", "yes")

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

    if BITGATE
        @testset "bit-for-bit against the captured baseline" begin
            compared = false
            for case in BSMEM_CASES
                path = baseline_path(case)
                isfile(path) || continue
                ref = try
                    deserialize(path)
                catch e
                    @info "$(case.name): baseline unreadable — $(sprint(showerror, e))" maxlog = 8
                    continue
                end
                compared = true
                got = run_bsmem_case(case; T = Float64)
                cmp = compare_bsmem(ref, got; rtol = 0.0)
                @test cmp.pass || (@error("$(case.name): $(cmp.msg)"); false)
            end
            # Asking for the gate and getting nothing compared is the failure the old version
            # could not report: it passed silently with every case skipped.
            @test compared
        end
    end
end
