# BSMEM/MaximENT regression gate.
#
# Compares the current code against the frozen Float64 baseline in test/references/.
# Capture the baseline first:  julia --project=. test/capture_bsmem_baseline.jl
#
# The gate is bit-for-bit by default. That is deliberate and not excessive: the solver
# contains discrete bisections on alpha (maximent.jl:409, :1149) where a 1-ulp shift
# flips a branch, changes the CG trust region, and diverges the whole trajectory — so a
# loose tolerance would be neither meaningful nor achievable. Refactor phases that are
# pure typing MUST reproduce Float64 exactly; if they don't, something changed.

using Test, Serialization

include(joinpath(@__DIR__, "bsmem_cases.jl"))

@testset "BSMEM regression (T=Float64, bit-for-bit)" begin
    any_baseline = false
    for case in BSMEM_CASES
        path = baseline_path(case)
        if !isfile(path)
            @info "no baseline for $(case.name) — run test/capture_bsmem_baseline.jl" maxlog = 1
            continue
        end
        any_baseline = true
        @testset "$(case.name)" begin
            ref = deserialize(path)
            got = run_bsmem_case(case; T = Float64)
            cmp = compare_bsmem(ref, got; rtol = 0.0)
            @test cmp.pass || (@error("$(case.name): $(cmp.msg)"); false)
        end
    end
    @test any_baseline  # the gate must not silently pass by having nothing to compare
end
