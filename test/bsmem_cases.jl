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
# ── Why the cases pin FFTW's planning strategy ──────────────────────────────────────
# MaximENT is chaotic: perturbations at the last bit amplify into whole-trajectory
# divergence. Measured on 2004-data1 at nx=64, method=[4,1,1,2] -- a case with NO RNG,
# fully deterministic within a process -- merely switching `--check-bounds`:
#
#     --check-bounds=auto :  16 iters, chisq=465.126, ntrans=385
#     --check-bounds=yes  :  16 iters, chisq=496.672, ntrans=353
#
# Bounds checking disables @inbounds SIMD, which reassociates reductions at ~1e-16 -- and that
# becomes a 7% chisq difference. So a case that fails here is worth being able to re-run and
# see fail again, which is what pinning the one input we control buys.

using OITOOLS, LinearAlgebra, FFTW

const DATA_DIR = joinpath(@__DIR__, "oifits_for_tests")

# The cases pin FFTW to ESTIMATE even though the package default is MEASURE (see FFT_FLAGS
# in oichi2.jl). MEASURE *times* candidate algorithms at plan time, so under varying system
# load it can pick a different — equally valid — one, changing results at ~1e-15. Measured:
# two identical runs at nx=128 differed by 6.4e-7 and 2.8e-7 in the final image. ESTIMATE is a
# fixed heuristic and therefore reproducible.
#
# The cases exist to detect *our* changes to the solver, so FFTW's algorithm choice is held
# fixed rather than left as an uncontrolled variable. It costs nothing in coverage — every line
# of solver code is exercised identically either way.
const CASE_FFTFLAGS = FFTW.ESTIMATE


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
    ft = setup_ft(data, case.nx, case.pixsize; fftflags = CASE_FFTFLAGS)
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
                                fftflags = CASE_FFTFLAGS, case.kw...)
        return (; ok = true, image = Array(img), history = hist)
    catch e
        return (; ok = false, err = _describe(e))
    end
end

"""
One line naming what actually went wrong.

`TaskFailedException` is unwrapped first: `showerror` puts the words "TaskFailedException" on
line one and the real cause below, so summarising by first line alone reports only that a
thread failed and never why.
"""
function _describe(e)
    while e isa TaskFailedException
        e = e.task.result
    end
    return first(split(sprint(showerror, e), "\n"))
end

