# benchmark.jl — wall time per replicate for each uncertainty method
#
#   julia --project=../.. benchmark.jl <outdir> [nboot]
#
# Single-threaded on purpose: the point is the cost of one replicate, not how
# well it parallelises (all the resampling schemes parallelise identically).
# Run `make_universes.jl` first — the benchmark uses one of its files.

include("common.jl")

const OUTDIR = length(ARGS) >= 1 ? ARGS[1] :
               joinpath(tempdir(), "oitools_bootstrap_validation")
const NBOOT  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 100

file = joinpath(OUTDIR, "ideal", "universe_0001.oifits")
isfile(file) || error("run make_universes.jl first (missing $file)")
data = readoifits(file; T=Float64)[1,1]
blk  = data_blocks(data)

println("Data:    ", data.nv2, " V2 + ", data.nt3phi, " T3PHI points, ",
        length(blk), " blocks, ", length(FREE), " free parameters")
println("Threads: ", Threads.nthreads(), " (benchmark runs single-threaded)")
println("Replicates per method: ", NBOOT, "\n")

common = (lb=LB, ub=UB, weights=WEIGHTS, verb=false, threaded=false)

# warm up the compiler
fit_model_lsqfit(FIT_MODEL, FREE, data; lb=LB, ub=UB, weights=WEIGHTS)
bootstrap_fit(FIT_MODEL, FREE, data; nboot=2, seed=1, common...)

t_fit = @elapsed for _ in 1:20
    fit_model_lsqfit(FIT_MODEL, FREE, data; lb=LB, ub=UB, weights=WEIGHTS)
end
t_fit /= 20

@printf("%-28s %10s %12s %10s\n", "method", "total (s)", "per repl (ms)", "vs 1 fit")
@printf("%-28s %10.3f %12.1f %10.1f\n", "single fit (LM + covar)", t_fit, 1000*t_fit, 1.0)

for (label, kw) in [
        ("bootstrap :replacement", (mode=:replacement, granularity=:config)),
        ("bootstrap :halfsample",  (mode=:halfsample,  granularity=:config)),
        ("bootstrap :weights",     (mode=:weights,     granularity=:config)),
        ("bootstrap :pmoired",     (mode=:pmoired,     granularity=:config)),
        ("bootstrap :epoch blocks",(mode=:replacement, granularity=:epoch)),
    ]
    t = @elapsed bootstrap_fit(FIT_MODEL, FREE, data; nboot=NBOOT, seed=7, common..., kw...)
    @printf("%-28s %10.3f %12.1f %10.1f\n", label, t, 1000*t/NBOOT, (t/NBOOT)/t_fit)
end

# and the same work spread over all available threads, for reference
if Threads.nthreads() > 1
    t = @elapsed bootstrap_fit(FIT_MODEL, FREE, data; nboot=NBOOT, seed=7,
                               lb=LB, ub=UB, weights=WEIGHTS, verb=false)
    @printf("\n%-28s %10.3f %12.1f  (on %d threads)\n",
            "bootstrap :replacement", t, 1000*t/NBOOT, Threads.nthreads())
end
