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
#
# But bit-for-bit equality is ONLY defined on the platform the baseline was captured on
# (linux x86_64) with a matching codegen env. This solver's chaos amplifies the ~1e-16
# reduction differences between architectures/BLAS and between Julia versions into whole-
# trajectory divergence — so on any other platform, or any other Julia patch, a bit compare
# is guaranteed to differ for reasons that are NOT a regression. We therefore assert the
# bit gate only where it is meaningful and smoke-test (does it still run?) everywhere else.
# The baselines live at linux x86_64 (see test/references/*.jls env); capture there.

using Test, Serialization

include(joinpath(@__DIR__, "bsmem_cases.jl"))

# The one platform the baselines are valid on. Off it, bit-for-bit is undefined (see above).
const BITGATE_PLATFORM = Sys.islinux() && Sys.ARCH === :x86_64

@testset "BSMEM regression (T=Float64, bit-for-bit)" begin
    any_baseline = false
    for case in BSMEM_CASES
        path = baseline_path(case)
        if !isfile(path)
            @info "no baseline for $(case.name) — run test/capture_bsmem_baseline.jl" maxlog = 1
            continue
        end
        ref = try
            deserialize(path)
        catch e
            # Almost always a Serialization format-version mismatch after a Julia upgrade:
            # these baselines are not forward-readable (see capture_bsmem_baseline.jl — the
            # permanent fix is to switch to JLD2). Skip this case rather than error the whole
            # suite; if EVERY baseline is unreadable, the `any_baseline` assertion below fails
            # loudly and asks for a re-cut.
            @info "$(case.name): baseline unreadable, skipping — $(sprint(showerror, e))" maxlog = 8
            continue
        end
        any_baseline = true
        @testset "$(case.name)" begin
            got = run_bsmem_case(case; T = Float64)
            cmp = compare_bsmem(ref, got; rtol = 0.0)
            env_ok = get(cmp, :comparable, true)  # false ⟺ codegen env differs from baseline
            if BITGATE_PLATFORM && env_ok
                @test cmp.pass || (@error("$(case.name): $(cmp.msg)"); false)
            else
                # Not the capture platform, or a different codegen env: bit-for-bit is
                # undefined here. Verify the case still runs the same way (success/failure),
                # and skip the exact comparison rather than report a bogus regression.
                reason = env_ok ? "platform is not linux x86_64" : cmp.msg
                @info "$(case.name): bit-for-bit gate skipped — $reason" maxlog = 8
                @test got.ok == ref.ok
                @test_skip cmp.pass
            end
        end
    end
    @test any_baseline  # the gate must not silently pass by having nothing to compare
end
