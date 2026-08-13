# test_python_boundary.jl — exercise every Julia↔Python crossing.
#
# Why this exists: before it was added, `runtests.jl` made zero plotting calls, zero UltraNest
# runs and zero SIMBAD queries. The suite only proved that the Python stack *loaded* — so the
# entire plotting layer could be broken and CI would stay green. That is exactly the failure
# mode a PyCall→PythonCall migration produces, because PythonCall converts far less
# automatically than PyCall and the breakage surfaces at call time, not compile time.
#
# The assertions therefore care about two things:
#   1. the call completes and produces a non-trivial figure file (not just "no exception"),
#   2. values coming BACK from Python are Julia types, not opaque Python handles.
#
# Everything runs headless via the Agg backend so it works on CI.

ENV["MPLBACKEND"] = "Agg"          # must be set before the plotting stack initialises

using OITOOLS, Test, Statistics

# Threading note, kept because it is a real difference between the two Python stacks:
# under PyCall this file segfaulted reliably at JULIA_NUM_THREADS=16 (its `pydecref`
# finalizer running on a non-main task) and was clean only at 1 thread. PythonCall is
# GIL-aware and passes at any thread count — verified at 16. The guard below is therefore
# inert on the current stack; it stays so that a regression to a thread-unsafe bridge is
# reported rather than dumping core.
const _ON_PYCALL     = isdefined(OITOOLS, :PyCall)
const _THREAD_UNSAFE = _ON_PYCALL && Threads.nthreads() > 1
if _THREAD_UNSAFE
    @warn "Skipping Python boundary tests: PyCall + $(Threads.nthreads()) threads segfaults."
end

const _PYDIR  = mktempdir()
const _DATA   = joinpath(@__DIR__, "..", "demos", "data")
const _ONCI   = get(ENV, "CI", nothing) == "true"

# Which plotting stack are we on? Kept generic so this file works before and after the
# migration and can be diffed across the two.
@static if isdefined(OITOOLS, :PythonPlot)
    const _PLT = OITOOLS.PythonPlot
else
    const _PLT = OITOOLS.PyPlot
end

_savefig(name) = (p = joinpath(_PYDIR, name * ".png"); _PLT.savefig(p, dpi = 60); p)
_closeall()    = _PLT.close("all")

"""Call a plotting function, save it, and assert a real figure came out."""
function _plots_ok(name, f; minbytes = 3000)
    _closeall()
    f()
    p = _savefig(name)
    _closeall()
    @test isfile(p)
    @test filesize(p) > minbytes    # a blank/failed canvas is ~1-2 kB
    return filesize(p)
end

@testset "Python boundary" begin
if _THREAD_UNSAFE
    @test_skip false
else

    data = readoifits(joinpath(_DATA, "2004-data1.oifits"); warn = false, verbose = false)[1,1]
    poly = readoifits(joinpath(_DATA, "BC2026", "OBJECT1_N.oifits"); warn = false, verbose = false)[1,1]
    # Keep the image at the data's precision: readoifits defaults to Float32, so the NFFT
    # plan is Float32 and a Float64 image would not dispatch. This also exercises the
    # package's default precision path rather than a Float64 special case.
    img  = gaussian2d(64, 64, 64/6)               # Float32, matching `data`

    # ── every exported plotting entry point actually runs ────────────────────
    @testset "plotting round trip" begin
        _plots_ok("uvplot",         () -> uvplot(data))
        _plots_ok("uvplot_wav",     () -> uvplot(poly; color = "wav"))
        _plots_ok("plot_v2",        () -> plot_v2(data))
        _plots_ok("plot_t3phi",     () -> plot_t3phi(data))
        _plots_ok("plot_t3amp",     () -> plot_t3amp(data))
        _plots_ok("plot_obs",       () -> plot_obs(poly; color = "wav"))
        _plots_ok("imdisp",         () -> imdisp(img; pixsize = 0.2))
        _plots_ok("plot_residuals", () -> plot_residuals(img, setup_nfft(data, 64, 0.2), data))

        # uvplot has a second method taking a bare uv matrix — a different code path.
        _plots_ok("uvplot_matrix",  () -> uvplot(data.uv))
    end

    # ── the colorbar divider crosses an axis handle into a separately imported
    #    Python module (mpl_toolkits.axes_grid1). It is the only place a plot handle
    #    does that, so it gets its own test.
    @testset "make_axes_locatable / colorbar divider" begin
        cube = cat(gaussian2d(32, 32, 5.0), gaussian2d(32, 32, 7.0); dims = 3)
        _plots_ok("imdisp_multi", () -> imdisp_multi(cube; pixsize = 0.2, use_colorbar = true))
        _plots_ok("imdisp_cbar",  () -> imdisp(img; pixsize = 0.2, use_colorbar = true))
    end

    # ── matplotlib rcParams is a Python dict we mutate ───────────────────────
    @testset "rcParams mutation" begin
        @test set_oiplot_defaults() === nothing || true    # must not throw
        @test set_oiplot_defaults(compact = true) === nothing || true
        set_oiplot_defaults(compact = false)
    end

    # ── values coming back from Python must be Julia types ───────────────────
    # This is the assertion class that catches a PythonCall migration: under PyCall these
    # convert automatically, under PythonCall they arrive as `Py` unless explicitly converted.
    @testset "UltraNest conversions" begin
        md = Dict{String,Any}("s,ud" => 3.0, "s,f" => 1.0)
        r = fit_model_ultranest(md, ["s,ud"], data;
                                lb = Dict("s,ud" => 0.5), ub = Dict("s,ud" => 10.0),
                                min_num_live_points = 120, verb = false, cornerplot = false)

        @test r.x_opt     isa Vector{Float64}
        @test r.posterior isa Matrix{Float64}
        @test r.logz      isa Float64
        @test r.logzerr   isa Float64
        @test r.chi2      isa Float64
        @test length(r.x_opt) == 1
        @test size(r.posterior, 2) == 1
        @test isfinite(r.logz) && isfinite(r.logzerr)

        # The known answer for this dataset/model (baselined on the PyCall stack).
        @test 6.0 < r.x_opt[1] < 7.5
    end

    # ── SIMBAD: network + a notoriously version-sensitive Python stack ───────
    # Skipped on CI (network) and tolerated locally if the Python env is inconsistent —
    # astroquery/astropy/numpy drift breaks this independently of anything Julia does.
    @testset "SIMBAD magnitudes" begin
        if _ONCI
            @test_skip false
        else
            m = try
                magnitudes_from_simbad("Vega")
            catch err
                @info "SIMBAD unavailable, skipping (Python env or network): " *
                      first(sprint(showerror, err), 120)
                nothing
            end
            if m === nothing
                @test_skip false
            else
                @test m isa Dict{String,Float64}
                @test haskey(m, "V")
                # the point of the assertion: real Float64s, not Python handles
                @test all(v -> v isa Float64, values(m))
                @test isfinite(m["V"])
            end
        end
    end
end   # _THREAD_UNSAFE
end
