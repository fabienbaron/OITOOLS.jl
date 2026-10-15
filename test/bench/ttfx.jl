# Time to first X: in this fresh process, loading OITOOLS, then the first call of each typical
# task (each includes its compilation). Appends one row per step to the CSV given:
#
#     for i in 1 2 3; do julia --project=. -t 1 test/bench/ttfx.jl ttfx.csv; done
#
# `test/bench/compare.jl` reads these rows with the suite's (median over the runs).
const T0 = time()
using OITOOLS
const T_USING = time() - T0
using Printf
const OUT = isempty(ARGS) ? "ttfx_$(VERSION).csv" : ARGS[1]
const DATA = joinpath(@__DIR__, "..", "..", "demos", "data")
rows = Tuple{String,Float64}[("using OITOOLS", T_USING)]
step(name, f) = (t = @elapsed(f()); push!(rows, (name, t)); t)

step("first readoifits (AX Cir)", () -> readoifits(joinpath(DATA, "AX_Cir.oifits"); T = Float64,
                                                    warn = false, verbose = false))
ax = readoifits(joinpath(DATA, "AX_Cir.oifits"); T = Float64, warn = false, verbose = false)[1, 1]
step("first fit_model_lsqfit (binary + smearing)", () -> fit_model_lsqfit(
    Dict{String,Any}("star,ud" => 0.8, "star,f" => "1 - \$c,f", "c,f" => 0.01, "c,x" => 6.0,
                     "c,y" => -28.0, "smear" => "auto"), ["star,ud", "c,f", "c,x", "c,y"], ax;
    weights = [1.0, 0.0, 1.0, 0, 0, 0, 0]))
bc = readoifits(joinpath(DATA, "Beauty_Contest_2004", "2004-data1.oifits"); T = Float64,
                warn = false, verbose = false)
step("first setup_ft + reconstruct (64², 50 iterations)", () -> begin
    ft = setup_ft(bc, 64, 0.2)
    x0 = gaussian2d(64, 64, 64 / 6)
    reconstruct(x0 ./ sum(x0), bc, ft; regularizers = [["centering", 1e5], ["l1l2", 1e7, 1e-3]],
                verb = false, maxiter = 50)
end)
push!(rows, ("total (process)", time() - T0))
mkpath(dirname(abspath(OUT)))
open(OUT, "a") do io
    for (n, t) in rows
        println(io, "ttfx,\"", n, "\",", @sprintf("%.3f", 1e3 * t), ",0,0,NaN,1,", VERSION)
    end
end
for (n, t) in rows
    @printf("%-52s %8.2f s\n", n, t)
end
