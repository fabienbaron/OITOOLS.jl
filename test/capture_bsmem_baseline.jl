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
#
# BLAS threading is pinned because _bidiag_svd calls LAPACK and reductions can
# reassociate with thread count, which would make the baseline machine-dependent.

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
