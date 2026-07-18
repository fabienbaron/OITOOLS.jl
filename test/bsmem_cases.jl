# Shared BSMEM/MaximENT regression cases.
#
# These pin the solver's behaviour so that refactors (notably the ongoing
# parameterization on a float type T) can be shown to be behaviour-preserving.
#
# Each case drives a distinct path through the solver:
#   - method 4 is the default: no RNG, and it bypasses evidence_estimate!/Ritz entirely,
#     so these runs are bit-for-bit deterministic.
#   - methods 1-3 drive trace_estimate_init!/evidence_estimate!/lanczos_bidiag!, and
#     consume the Xoshiro(iseed) stream.
#   - mackay_alpha / ritz_alpha select alternative alpha-update strategies. ritz_alpha
#     needs Ritz pairs from evidence_estimate!, so it only makes sense with method 1-3.
#
# A case that errors is recorded as such rather than skipped: preserving an existing
# failure is still a regression signal.
#
# ── Why this gate is bit-for-bit, and why it pins codegen ────────────────────────────
# MaximENT is chaotic: perturbations at the last bit amplify into whole-trajectory
# divergence. Measured on 2004-data1 at nx=64, method=[4,1,1,2] — a case with NO RNG,
# fully deterministic within a process — merely switching `--check-bounds`:
#
#     --check-bounds=auto :  16 iters, chisq=465.126, ntrans=385
#     --check-bounds=yes  :  16 iters, chisq=496.672, ntrans=353
#
# Bounds checking disables @inbounds SIMD, which reassociates reductions at ~1e-16 —
# and that becomes a 7% chisq difference. Consequences:
#   1. A tolerance-based gate is meaningless here; only bit-for-bit is well-defined.
#   2. A baseline is valid ONLY under the codegen flags that produced it, so we record
#      them and refuse to compare across a mismatch rather than report a false failure.
# Note Pkg.test() defaults to --check-bounds=yes, so baselines are captured that way.

using OITOOLS, Serialization, LinearAlgebra, FFTW

const REF_DIR  = joinpath(@__DIR__, "references")
const DATA_DIR = joinpath(@__DIR__, "oifits_for_tests")

# The gate pins FFTW to ESTIMATE even though the package default is MEASURE (see FFT_FLAGS
# in oichi2.jl). MEASURE *times* candidate algorithms at plan time, so under varying system
# load it can pick a different — equally valid — one, changing results at ~1e-15. Measured:
# two identical gate runs at nx=128 differed by 6.4e-7 and 2.8e-7 in the final image.
# ESTIMATE is a fixed heuristic and therefore bit-reproducible.
#
# This is deliberate: the gate exists to detect *our* changes to the solver, so FFTW's
# algorithm choice must be held fixed rather than left as an uncontrolled variable. It costs
# nothing in coverage — every line of solver code is exercised identically either way.
const GATE_FFTFLAGS = FFTW.ESTIMATE

"""
Environment a bit-for-bit baseline is only valid under.

`fft_flags` is included because FFTW may pick a different (equally valid) algorithm per
planning strategy — MEASURE vs ESTIMATE differ at ~1e-15 from nx=128 up — and this solver
amplifies that into a different trajectory. Recording it means switching strategy reports
itself as an environment change rather than masquerading as a broken refactor.
"""
codegen_env() = (; check_bounds = Int(Base.JLOptions().check_bounds),
                   opt_level    = Int(Base.JLOptions().opt_level),
                   nthreads     = Threads.nthreads(),
                   blas_threads = BLAS.get_num_threads(),
                   fft_flags    = string(GATE_FFTFLAGS),
                   julia        = string(VERSION))

const BSMEM_CASES = (
    (name = "mono_default_m4", file = "2004-data1.oifits", nx = 64, pixsize = 0.2, nwav = 1,
     kw = (; method = [4, 1, 1, 2], maxiter = 40)),
    (name = "mono_default_m4_nx128", file = "2004-data1.oifits", nx = 128, pixsize = 0.1, nwav = 1,
     kw = (; method = [4, 1, 1, 2], maxiter = 20)),
    (name = "mono_method1", file = "2004-data1.oifits", nx = 64, pixsize = 0.2, nwav = 1,
     kw = (; method = [1, 1, 1, 2], maxiter = 20, nrand = 10, iseed = 0)),
    (name = "mono_method2", file = "2004-data1.oifits", nx = 64, pixsize = 0.2, nwav = 1,
     kw = (; method = [2, 1, 1, 2], maxiter = 20, nrand = 10, iseed = 0)),
    (name = "mono_method3", file = "2004-data1.oifits", nx = 64, pixsize = 0.2, nwav = 1,
     kw = (; method = [3, 1, 1, 2], maxiter = 20, nrand = 10, iseed = 0)),
    (name = "mono_mackay", file = "2004-data1.oifits", nx = 64, pixsize = 0.2, nwav = 1,
     kw = (; method = [4, 1, 1, 2], maxiter = 30, mackay_alpha = true)),
    (name = "mono_ritz", file = "2004-data1.oifits", nx = 64, pixsize = 0.2, nwav = 1,
     kw = (; method = [2, 1, 1, 2], maxiter = 20, nrand = 10, iseed = 0, ritz_alpha = true)),
    (name = "poly_object1n", file = "OBJECT1_N.oifits", nx = 32, pixsize = 0.5, nwav = 4,
     kw = (; method = [4, 1, 1, 2], maxiter = 15)),
)

"Load data + FT plans for a case at precision `T`."
function bsmem_case_inputs(case; T = Float64)
    poly = case.nwav > 1
    data = readoifits(joinpath(DATA_DIR, case.file); verbose = false, warn = false, T,
                      polychromatic = poly, merge_oi_wavelength = poly)
    if poly
        @assert size(data, 1) >= case.nwav "$(case.file) has $(size(data,1)) channels, need $(case.nwav)"
        data = data[1:case.nwav, :]
    end
    ft = setup_ft(data, case.nx, case.pixsize; fftflags = GATE_FFTFLAGS)
    x0 = gaussian2d(case.nx, case.nx, case.nx / 6; T)
    if poly
        return (repeat(x0, 1, 1, case.nwav),
                [data[w, 1] for w in 1:case.nwav],
                [ft[w, 1] for w in 1:case.nwav])
    else
        return x0, data, ft
    end
end

"Run one case, returning its image and convergence trace (or the error it raised)."
function run_bsmem_case(case; T = Float64)
    x0, data, ft = bsmem_case_inputs(case; T)
    hist = NamedTuple[]
    try
        img = reconstruct_bsmem(x0, data, ft; verbose = false, history = hist,
                                fftflags = GATE_FFTFLAGS, case.kw...)
        return (; ok = true, image = Array(img), history = hist, env = codegen_env())
    catch e
        return (; ok = false, err = first(split(sprint(showerror, e), "\n")), env = codegen_env())
    end
end

baseline_path(case) = joinpath(REF_DIR, "bsmem_baseline_$(case.name).jls")

"""
    compare_bsmem(ref, got; rtol)

Compare a case result against its baseline. Returns `(; pass, msg)`, plus `comparable=false`
when the recorded codegen env differs (in which case a bit-for-bit verdict is meaningless and
the caller should skip, not fail). `rtol=0` demands bit-for-bit equality (the gate for
pure-typing refactor phases).

Note the env does NOT record OS/arch — bit-for-bit equality is only defined on the platform the
baseline was captured on (linux x86_64), because this solver is chaotic and reductions differ
across architecture and BLAS. The caller must additionally gate on platform; see
`test_bsmem_regression.jl`.
"""
function compare_bsmem(ref, got; rtol = 0.0)
    # A bit-for-bit baseline is only meaningful under the codegen it was captured with
    # (see the header): report the mismatch rather than a misleading numeric failure.
    if rtol == 0 && hasproperty(ref, :env) && ref.env != got.env
        diff = [string(k, ": ", getfield(ref.env, k), " -> ", getfield(got.env, k))
                for k in keys(ref.env) if getfield(ref.env, k) != getfield(got.env, k)]
        return (; pass = false, comparable = false,
            msg = "codegen/threading mismatch, baseline not comparable — " * join(diff, ", ") *
                  ". Re-run under the baseline's settings (Pkg.test() uses --check-bounds=yes), " *
                  "or delete test/references/bsmem_baseline_*.jls and re-capture.")
    end
    ref.ok == got.ok || return (; pass = false,
        msg = "ok mismatch: baseline ok=$(ref.ok) now ok=$(got.ok)" *
              (got.ok ? "" : " ($(got.err))"))
    ref.ok || return (; pass = ref.err == got.err,
        msg = "error text: baseline «$(ref.err)» now «$(got.err)»")

    size(ref.image) == size(got.image) ||
        return (; pass = false, msg = "image size $(size(ref.image)) vs $(size(got.image))")
    length(ref.history) == length(got.history) ||
        return (; pass = false, msg = "iteration count $(length(ref.history)) vs $(length(got.history))")

    if rtol == 0
        ref.image == got.image ||
            return (; pass = false, msg = "image not bit-identical (max|Δ|=$(maximum(abs.(Float64.(got.image) .- ref.image))))")
    else
        d = norm(Float64.(got.image) .- ref.image) / max(norm(ref.image), eps())
        d <= rtol || return (; pass = false, msg = "image rel L2 $d > rtol $rtol")
    end

    # The trace is the sharper fingerprint: it catches a divergence long before it
    # shows up in the final image.
    for (k, (r, g)) in enumerate(zip(ref.history, got.history))
        for f in keys(r)
            rv, gv = getfield(r, f), getfield(g, f)
            rv isa Real && gv isa Real || continue
            if rtol == 0
                rv === gv || return (; pass = false, msg = "iter $k field $f: $rv vs $gv (not bit-identical)")
            else
                d = abs(gv - rv) / max(abs(rv), eps())
                d <= rtol || return (; pass = false, msg = "iter $k field $f: $rv vs $gv (rel $d > $rtol)")
            end
        end
    end
    return (; pass = true, msg = "ok")
end
