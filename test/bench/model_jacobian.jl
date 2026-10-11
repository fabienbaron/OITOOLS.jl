# Baseline and progress benchmarks for the model-Jacobian work (PLAN_ModelJacobian.md):
# BenchmarkTools timings and JET report counts over a zoo of models, in Float64 and Float32.
# Run by hand, not part of runtests.jl:
#
#     julia --project=. -t 8 test/bench/model_jacobian.jl            # everything
#     julia --project=. -t 8 test/bench/model_jacobian.jl quick      # the zoo only
#
# BenchmarkTools and JET come from the default environment, which stacks on the project.
using OITOOLS, BenchmarkTools, JET, Printf
const O = OITOOLS

include(joinpath(@__DIR__, "zoo.jl"))
const QUICK = "quick" in ARGS

t(b) = median(b).time; m(b) = median(b).memory; a(b) = median(b).allocs
cell(b) = @sprintf("%9s %9s %6d", BenchmarkTools.prettytime(t(b)), BenchmarkTools.prettymemory(m(b)), a(b))
jet(f) = try length(JET.get_reports(f())) catch e; -1 end

println("model jacobian benchmarks — ", Threads.nthreads(), " threads, ", VERSION)
for T in (Float64, Float32)
    d = readoifits(AXCIR; T, warn = false, verbose = false)[1, 1]
    uv, wl, mjd = d.uv, d.uv_lam, d.uv_mjd
    println("\n── $T, ", size(uv, 2), " uv points")
    @printf("%-24s %-28s %-28s %-28s %-28s %-28s %s\n", "model", "eval_model", "eval_model_grad",
            "eval_model_jac!", "model_to_chi2_fg", "residuals_flat_jac", "JET eval/grad/jac!")
    for (name, md, free) in ZOO
        fm = dict_to_model(md, free)
        x = T.([md[k] for k in free])
        b1 = @benchmark eval_model($fm, $x, $uv; wl = $wl, mjd = $mjd)
        b2 = @benchmark eval_model_grad($fm, $x, $uv; wl = $wl, mjd = $mjd)
        ws = ModelJacWorkspace(fm, x, uv)
        V = Vector{Complex{T}}(undef, size(uv, 2)); J = Matrix{Complex{T}}(undef, size(uv, 2), length(x))
        b5 = @benchmark eval_model_jac!($V, $J, $ws, $fm, $x, $uv; wl = $wl, mjd = $mjd)
        b3 = @benchmark model_to_chi2_fg($fm, $x, $d)
        b4 = @benchmark O.residuals_flat_jac($fm, $x, $d)
        j1 = jet(() -> @report_opt target_modules = (O,) eval_model(fm, x, uv; wl, mjd))
        j2 = jet(() -> @report_opt target_modules = (O,) eval_model_grad(fm, x, uv; wl, mjd))
        j3 = jet(() -> @report_opt target_modules = (O,) eval_model_jac!(V, J, ws, fm, x, uv; wl, mjd))
        @printf("%-24s %s   %s   %s   %s   %s   %d/%d/%d\n", name, cell(b1), cell(b2), cell(b5),
                cell(b3), cell(b4), j1, j2, j3)
    end
end

QUICK && exit()

println("\n── end-to-end fits (Float64, AX Cir, V² + T3φ)")
d = readoifits(AXCIR; T = Float64, warn = false, verbose = false)[1, 1]
w = [1.0, 0.0, 1.0]
for (name, md, free) in ZOO[[1, 4, 8]]
    start = merge(md, Dict{String,Any}(k => md[k] * 1.05 for k in free))
    for (fitter, f) in (("fit_model LD_MMA", () -> fit_model(start, free, d; weights = w, method = :LD_MMA)),
                        ("fit_model_lsqfit", () -> fit_model_lsqfit(start, free, d; weights = w)))
        r = f()
        b = @benchmark $f() samples = 10 evals = 1
        @printf("%-10s %-18s %s   χ²ᵣ %.5f   x = %s\n", name, fitter, cell(b), r.chi2r,
                join((@sprintf("%.6g", v) for v in r.x_opt), ", "))
    end
end

println("\n── companion search, general path (AX Cir, ud + flux free per node)")
for sm in (0, "auto")
    md = Dict{String,Any}("star,ud" => 0.83, "star,f" => "1 - \$comp,f", "comp,f" => 0.009,
                          "comp,x" => 0.0, "comp,y" => 0.0, "smear" => sm)
    for th in (false, true)
        b = @benchmark companion_search($md, $d; weights = $w, free = ["comp,f", "star,ud"],
                                        threads = $th) samples = 3 evals = 1
        @printf("smear=%-5s threads=%-5s %s\n", sm, th, cell(b))
    end
end

println("\n── time to first evaluation of a new model (fresh process)")
script = """
using OITOOLS
d = readoifits("$AXCIR"; T = Float64, warn = false, verbose = false)[1, 1]
t1 = @elapsed (m = dict_to_model(Dict{String,Any}("s,ud" => 0.9, "c,f" => 0.1, "c,x" => 3.0, "c,y" => 1.0), ["s,ud", "c,f"]);
               eval_model(m, [0.9, 0.1], d.uv; wl = d.uv_lam))
t2 = @elapsed eval_model_grad(m, [0.9, 0.1], d.uv; wl = d.uv_lam)
t3 = @elapsed model_to_chi2_fg(m, [0.9, 0.1], d)
t4 = @elapsed eval_model_jac(m, [0.9, 0.1], d.uv; wl = d.uv_lam)
# The same structure rebuilt with other values: a new resolver, which must not recompile.
t5 = @elapsed (m2 = dict_to_model(Dict{String,Any}("s,ud" => 0.8, "c,f" => 0.2, "c,x" => -3.0, "c,y" => 2.0), ["s,ud", "c,f"]);
               eval_model_jac(m2, [0.8, 0.2], d.uv; wl = d.uv_lam))
println("first eval_model ", round(t1; digits = 3), " s, first eval_model_grad ", round(t2; digits = 3),
        " s, first model_to_chi2_fg ", round(t3; digits = 3), " s, first eval_model_jac ",
        round(t4; digits = 3), " s, rebuilt model eval_model_jac ", round(t5; digits = 3), " s")
"""
run(`$(Base.julia_cmd()) --project=$(pkgdir(OITOOLS)) -e $script`)
