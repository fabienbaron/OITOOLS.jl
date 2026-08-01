# run_oitools.jl — fit every simulated universe with OITOOLS and record, for
# each uncertainty method, the 1σ it reports.
#
#   julia --project=../.. -t auto run_oitools.jl <outdir> [nboot] [regime ...]
#
# Output: <outdir>/results_oitools.tsv  (long format)
#   regime  universe  method  param  best  sigma  chi2r
#
# Methods
#   lsqfit      analytic covariance from Levenberg-Marquardt, rescaled by chi2r
#   montecarlo  refit of perturb_data replicates (parametric, trusts the error bars)
#   boot_config block bootstrap over (MJD, baseline/triangle), with replacement
#   boot_epoch  block bootstrap over whole epochs
#   boot_weights multiplier (Bayesian) bootstrap: continuous block weights
#   boot_halfsample balanced repeated replication: one random half per replicate
#   boot_pmoired PMOIRED's two-half-samples scheme, same blocks as boot_config

include("common.jl")

"""
    montecarlo_sigma(data, x0; nboot, seed) -> Vector{Float64}

Parametric Monte Carlo: refit `nboot` copies of the data perturbed by Gaussian
draws from the quoted error bars, and return the 16/84 half-width of the
resulting parameter distribution.
"""
function montecarlo_sigma(data::OIdata, x0::Vector{Float64}; nboot::Int, seed::Int)
    fm = parse_model(FIT_MODEL, FREE; nB_workspace=size(data.uv, 2))
    S  = fill(NaN, nboot, length(FREE))
    Threads.@threads for i in 1:nboot
        rng = Xoshiro(hash((seed, i)))
        try
            r = fit_model_lsqfit(deepcopy(fm), x0, perturb_data(data; rng=rng);
                                 lb=LB, ub=UB, weights=WEIGHTS)
            S[i, :] = r.x_opt
        catch
        end
    end
    ok = [all(isfinite, @view S[i, :]) for i in 1:nboot]
    return [0.5 * (quantile(S[ok, j], 0.84) - quantile(S[ok, j], 0.16))
            for j in eachindex(FREE)]
end

const OUTDIR = ARGS[1]
const NBOOT  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 200
const WHICH  = length(ARGS) >= 3 ? ARGS[3:end] : REGIMES

rows = Vector{Any}[]
t_start = time()

for regime in WHICH
    dir   = joinpath(OUTDIR, regime)
    files = sort(filter(f -> endswith(f, ".oifits"), readdir(dir; join=true)))
    isempty(files) && continue
    @printf("%-15s : %d universes\n", regime, length(files))

    for (u, file) in enumerate(files)
        data = readoifits(file; T=Float64)[1,1]

        common = (lb=LB, ub=UB, weights=WEIGHTS, verb=false)

        f = fit_model_lsqfit(FIT_MODEL, FREE, data; lb=LB, ub=UB, weights=WEIGHTS)
        best  = f.x_opt
        chi2r = f.chi2r

        methods = Dict{String,Vector{Float64}}()
        methods["lsqfit"] = f.stderror

        # Parametric Monte Carlo.  No longer part of the OITOOLS API (it was
        # never better than the analytic covariance and collapses under
        # systematics) but kept here inline so the evidence stays reproducible.
        methods["montecarlo"] = montecarlo_sigma(data, best; nboot=NBOOT, seed=u)

        for (name, kw) in (("boot_config",     (mode=:replacement, granularity=:config)),
                           ("boot_epoch",      (mode=:replacement, granularity=:epoch)),
                           ("boot_weights",    (mode=:weights,     granularity=:config)),
                           ("boot_halfsample", (mode=:halfsample,  granularity=:config)),
                           ("boot_pmoired",    (mode=:pmoired,     granularity=:config)))
            b = bootstrap_fit(FIT_MODEL, FREE, data; nboot=NBOOT, seed=u, common..., kw...)
            methods[name] = b.sigma
        end

        for (name, sig) in methods, (j, p) in enumerate(FREE)
            push!(rows, Any[regime, u, name, p, best[j], sig[j], chi2r])
        end

        if u % 10 == 0
            @printf("  %4d/%d  (%.1f s elapsed)\n", u, length(files), time() - t_start)
        end
    end
end

out = joinpath(OUTDIR, "results_oitools.tsv")
open(out, "w") do io
    println(io, "regime\tuniverse\tmethod\tparam\tbest\tsigma\tchi2r")
    for r in rows
        println(io, join(r, "\t"))
    end
end
@printf("\nWrote %d rows to %s  (%.1f s)\n", length(rows), out, time() - t_start)
