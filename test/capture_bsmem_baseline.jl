# Freeze the BSMEM/MaximENT baseline.
#
# Run this ONCE against the pre-refactor code, then never again until the reference is
# deliberately re-cut:
#
#     julia --project=. test/capture_bsmem_baseline.jl
#
# The refactor gate (test/runtests.jl) compares against what this writes. Re-running it
# after a change would silently bless that change — which is exactly what the gate exists
# to prevent. Delete a baseline file explicitly if you mean to re-cut it.
#
# Baseline history (each re-cut must be a deliberate, recorded decision):
#   1. Initial capture, pre-refactor.
#   2. Re-cut when FFT_FLAGS became FFTW.MEASURE (was NFFT's unset default = FFTW.ESTIMATE).
#      FFTW selects a different, equally valid algorithm, differing at ~1e-15 from nx=128 up;
#      only mono_default_m4_nx128 moved (16 vs 17 iterations) — every nx<=64 case stayed
#      bit-identical, exactly as predicted. Motivation: ~1.3-1.6x on Float64 alone, and it
#      turns Float32 from a regression into a win at nx>=128. See FFT_FLAGS in oichi2.jl.
#   3. Re-cut on a developer workstation, then REVERTED. The new baselines matched that
#      machine but not the CI runner, so CI failed all 8 cases (max|Δ| 1e-16 to 8e-5, and
#      mono_default_m4 took 17 iterations instead of 16). The lesson is the warning below:
#      these files encode the capture hardware, so they must be cut where CI runs.
#
# BLAS threading is pinned because _bidiag_svd calls LAPACK and reductions can
# reassociate with thread count, which would make the baseline machine-dependent.
#
# WARNING — the baseline is HARDWARE-specific, and `env` does not record that.
# `codegen_env()` captures check_bounds/opt_level/nthreads/blas_threads/fft_flags/julia, but
# NOT the CPU. Two machines with identical Julia settings therefore look "comparable" while
# producing images that differ at 1e-16 to 1e-4 (this solver is chaotic; a one-ULP difference
# in an early reduction changes the iterate trajectory and even the iteration count).
#
# Consequence: capture on hardware matching CI (a GitHub linux x86_64 runner), NOT on a
# developer workstation. Re-cutting these on a different machine makes the gate pass locally
# and fail on CI for all 8 cases — which is exactly what happened once; see history entry 3.
#
# Storage format: the Serialization stdlib (.jls), chosen only because it adds no
# dependency. Its format is NOT guaranteed stable across Julia versions, so a Julia
# upgrade may make existing baselines unreadable and force a re-cut. That fails loudly
# (`deserialize` throws) rather than silently passing, and `env` pins the Julia version,
# so it cannot go unnoticed — but the pre-upgrade reference is lost when it happens.
# If this becomes annoying, switch to JLD2: version-stable, no text bulk, drop-in for the
# serialize/deserialize calls here and in test_bsmem_regression.jl.

using LinearAlgebra
BLAS.set_num_threads(1)

include(joinpath(@__DIR__, "bsmem_cases.jl"))

mkpath(REF_DIR)
println("Capturing BSMEM baseline (BLAS threads = $(BLAS.get_num_threads()))")
println()

for case in BSMEM_CASES
    path = baseline_path(case)
    if isfile(path)
        println("  SKIP  $(case.name)  (exists — delete it to re-cut)")
        continue
    end
    t = @elapsed res = run_bsmem_case(case; T = Float64)
    serialize(path, res)
    if res.ok
        println("  OK    $(rpad(case.name, 22)) $(length(res.history)) iters, ",
                "final chisq=", round(res.history[end].chisq, sigdigits = 8),
                ", ntrans=", res.history[end].ntrans, "  [$(round(t, digits=1))s]")
    else
        println("  ERR   $(rpad(case.name, 22)) recorded failure: $(res.err)  [$(round(t, digits=1))s]")
    end
end

println()
println("Baseline written to $REF_DIR")
