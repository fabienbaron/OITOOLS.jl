# Whole-problem benchmarks: reading files, fitting models, searching, reconstructing images,
# simulating — the calls and data the test suites already use (demos/data). Each row: the
# median time, memory and allocations (BenchmarkTools), and a result value, so that two runs
# (two Julia versions, two commits) are seen to compute the same thing. Run by hand, one thread:
#
#     JULIA_LOAD_PATH="@:<env with BenchmarkTools>:@stdlib" julia --project=. -t 1 \
#         test/bench/suite.jl [quick] [out.csv]
#
# `test/bench/compare.jl a.csv b.csv` puts two runs side by side.
using OITOOLS, BenchmarkTools, Printf, LinearAlgebra
const O = OITOOLS
BLAS.set_num_threads(1)

const QUICK = "quick" in ARGS
const OUT = something(findfirst(a -> endswith(a, ".csv"), ARGS) |> i -> i === nothing ? nothing : ARGS[i],
                      "bench_$(VERSION).csv")
const DATA = joinpath(@__DIR__, "..", "..", "demos", "data")
const ROWS = Any[]

# `f` is benchmarked; `result(f())` is the row's value, computed once outside the timing
function bench(group, name, f; result = x -> NaN, seconds = QUICK ? 3 : 10, samples = 100)
    f()                                                    # compiled before timing
    b = @benchmark $f() seconds = seconds samples = samples evals = 1
    r = result(f())
    t = median(b)
    push!(ROWS, (group, name, t.time / 1e6, t.memory, t.allocs, Float64(r), length(b.times)))
    @printf("%-9s %-46s %11s %11s %9d   %-14.10g (%d samples)\n", group, name,
            BenchmarkTools.prettytime(t.time), BenchmarkTools.prettymemory(t.memory), t.allocs,
            r, length(b.times))
    flush(stdout)
end

rd(f; T = Float64, kw...) = readoifits(f; T, warn = false, verbose = false, kw...)

println("OITOOLS benchmarks — Julia ", VERSION, ", ", Threads.nthreads(), " thread(s), ",
        QUICK ? "quick" : "full")

# ── reading ─────────────────────────────────────────────────────────────────────────────
axf = joinpath(DATA, "AX_Cir.oifits")
grav = joinpath(DATA, "HD210763", "2022-08-03-calibrated_cleaned.fits")
bench("read", "readoifits AX Cir, Float64", () -> rd(axf); result = d -> d[1, 1].nv2)
bench("read", "readoifits AX Cir, Float32", () -> rd(axf; T = Float32); result = d -> d[1, 1].nv2)
bench("read", "readoifits GRAVITY (HD 210763), Float64", () -> rd(grav); result = d -> d[1, 1].nv2)
bench("read", "readoifits GRAVITY (HD 210763), Float32", () -> rd(grav; T = Float32);
      result = d -> d[1, 1].nv2)

# ── models ──────────────────────────────────────────────────────────────────────────────
include(joinpath(@__DIR__, "zoo.jl"))
ax = rd(axf)[1, 1]
for (name, md, free) in (QUICK ? ZOO[1:3] : ZOO)
    fm = dict_to_model(md, free)
    x = Float64[md[k] for k in free]
    bench("model", "eval_model $(name)", () -> eval_model(fm, x, ax.uv; wl = ax.uv_lam, mjd = ax.uv_mjd);
          result = v -> sum(abs, v))
    ws = ModelJacWorkspace(fm, x, ax.uv)
    V = Vector{ComplexF64}(undef, size(ax.uv, 2)); J = Matrix{ComplexF64}(undef, size(ax.uv, 2), length(x))
    bench("model", "eval_model_jac! $(name)",
          () -> (eval_model_jac!(V, J, ws, fm, x, ax.uv; wl = ax.uv_lam, mjd = ax.uv_mjd); J);
          result = J -> sum(abs, J))
end

w3 = [1.0, 0.0, 1.0, 0, 0, 0, 0]
bench("fit", "lsqfit uniform disc (AX Cir)",
      () -> fit_model_lsqfit(Dict{String,Any}("star,ud" => 0.8, "smear" => "auto"), ["star,ud"], ax; weights = w3);
      result = r -> r.chi2)
binary = Dict{String,Any}("star,ud" => 0.8, "star,f" => "1 - \$companion,f", "companion,f" => 0.01,
                          "companion,x" => 6.0, "companion,y" => -28.0, "smear" => "auto")
bfree = ["star,ud", "companion,f", "companion,x", "companion,y"]
bench("fit", "lsqfit binary + smearing (AX Cir, EX3)", () -> fit_model_lsqfit(binary, bfree, ax; weights = w3);
      result = r -> r.chi2)
aca = rd(joinpath(DATA, "AlphaCenA.oifits"))[1, 1]
ldpow = Dict{String,Any}("star,ldpow" => 8.3, "star,alpha" => 0.2, "star,f" => 1.0)
bench("fit", "lsqfit power-law limb darkening (α Cen A)",
      () -> fit_model_lsqfit(ldpow, ["star,ldpow", "star,alpha"], aca; weights = [1.0, 0, 0, 0, 0, 0, 0]);
      result = r -> r.chi2)
cs = Dict{String,Any}("star,ud" => 0.93, "star,f" => "1 - \$comp,f", "comp,f" => 0.01,
                      "comp,x" => 0.0, "comp,y" => 0.0, "smear" => "auto")
bench("search", "companion_search, default grid (AX Cir)",
      () -> companion_search(cs, ax; weights = [1.0, 0.0, 1.0], free = ["comp,f", "star,ud"]);
      result = g -> minimum(filter(isfinite, g.chi2r)), seconds = QUICK ? 5 : 30, samples = 10)
lcfree = ["star,ldlin", "star,u"]
lin = Dict{String,Any}("star,ldlin" => 8.0, "star,u" => 0.3, "star,f" => 1.0)
bench("search", "bootstrap_fit 40 replicates (α Cen A)",
      () -> bootstrap_fit(lin, lcfree, aca; lb = Dict("star,ldlin" => 5.0, "star,u" => 0.0),
                          ub = Dict("star,ldlin" => 12.0, "star,u" => 1.0), weights = [1.0, 0.0, 0.0],
                          nboot = 40, seed = 99, verb = false, threaded = false);
      result = b -> b.median[1], seconds = QUICK ? 5 : 30, samples = 10)
# a binary through V1295 Aql's uv points, noiseless; the multistart test's grid
src = joinpath(DATA, "V1295_Aql.oifits")
truth = Dict{String,Any}("A,ud" => 0.4, "B,f" => 0.3, "B,ud" => 0.1, "B,x" => 3.0, "B,y" => -2.0)
simfile = joinpath(mktempdir(), "binary.oifits")
redirect_stdout(devnull) do
    simulate_from_oifits(src, simfile; mode = "copy_errors", noise = false,
                         flat_model = dict_to_model(truth, String[]), flat_params = Float64[])
end
vb = rd(simfile)[1, 1]
ms = merge(copy(truth), Dict{String,Any}("B,f" => 0.5, "B,x" => 0.0, "B,y" => 0.0))
bench("search", "multistart_fit 7×7 grid (binary, V1295 Aql uv)",
      () -> multistart_fit(ms, ["B,f", "B,x", "B,y"], vb; grid = ["B,x" => (-6, 6, 2), "B,y" => (-6, 6, 2)],
                           lb = Dict("B,f" => 0.0), ub = Dict("B,f" => 1.0), verbose = false, threaded = false);
      result = r -> r.minima[1].fit.chi2, seconds = QUICK ? 5 : 30, samples = 10)

# ── imaging ─────────────────────────────────────────────────────────────────────────────
bc = rd(joinpath(DATA, "Beauty_Contest_2004", "2004-data1.oifits"))
nx, ps = 64, 0.2
bench("imaging", "setup_ft NFFT 64² (Beauty Contest 2004)", () -> setup_ft(bc, nx, ps); result = _ -> nx)
ft = setup_ft(bc, nx, ps)
x0 = gaussian2d(nx, nx, nx / 6); x0 ./= sum(x0)
bench("imaging", "image_to_chi2 64²", () -> image_to_chi2(x0, ft, bc; verb = false); result = c -> c)
x4 = reshape(copy(x0), nx, nx, 1, 1); g4 = similar(x4)
bench("imaging", "image_to_chi2_fg 64²", () -> image_to_chi2_fg(x4, g4, ft, bc); result = c -> c)
regs = [["centering", 1e5], ["l1l2", 1e7, 1e-3]]
bench("imaging", "reconstruct 200 iterations 64², centering + l1l2",
      () -> reconstruct(x0, bc, ft; regularizers = regs, verb = false, maxiter = 200);
      result = x -> image_to_chi2(x, ft, bc; verb = false), seconds = QUICK ? 5 : 30, samples = 10)
bench("imaging", "reconstruct_bsdmm 15 iterations 64²",
      () -> reconstruct_bsdmm(x0, bc, ft; mu_reg = 1e-3, mu_cen = 1e-3, maxit = 15, x_maxiter = 5,
                              verb = false);
      result = x -> image_to_chi2(x isa Tuple ? x[1] : x, ft, bc; verb = false),
      seconds = QUICK ? 5 : 30, samples = 10)
mwc = rd(joinpath(DATA, "MWC480.oifits"); filter_bad_data = true, polychromatic = true)
nw = size(mwc, 1)
ftp = setup_ft(mwc, nx, 0.1)
xp0 = Float64.(repeat(gaussian2d(nx, nx, nx / 8), 1, 1, nw, 1)); xp0 ./= sum(xp0[:, :, 1, 1])
bench("imaging", "reconstruct polychromatic 50 iterations ($(nw) λ, MWC 480)",
      () -> reconstruct(xp0, mwc, ftp; regularizers = [["centering", 1e7]],
                        transspectral_regularizers = [["transspectral_structnorm", 1e2],
                                                      ["transspectral_tvsq", 1e3]],
                        maxiter = 50, verb = false);
      result = x -> sum(x), seconds = QUICK ? 5 : 30, samples = 10)

# ── simulating ──────────────────────────────────────────────────────────────────────────
cube = joinpath(mktempdir(), "cube.fits")
λc = [1.52e-6, 1.6e-6, 1.68e-6]
g(σ) = (c = 17; [exp(-((i - c)^2 + (j - c)^2) / (2σ^2)) for i in 1:33, j in 1:33])
writefits(cat((g(2 + k) for k in 1:3)...; dims = 3), cube; pixsize = 0.1, wavelengths = λc)
fac = read_facility_file("CHARA"); tgt = read_obs_file("default_obs"); comb = read_comb_file("MIRCX")
using Dates
dates = collect(DateTime(2024, 8, 13, 3):Minute(30):DateTime(2024, 8, 13, 8))
simout = joinpath(mktempdir(), "sim.oifits")
bench("simulate", "simulate 3-channel cube, CHARA/MIRC-X, 11 epochs",
      () -> redirect_stdout(devnull) do
          simulate(fac, tgt, comb, :image, dates, simout; image = cube, seed = 1, n_samples = 20)
      end; result = _ -> rd(simout)[1, 1].nv2, seconds = QUICK ? 5 : 20, samples = 20)

mkpath(dirname(abspath(OUT)))
open(OUT, "w") do io
    println(io, "group,name,median_ms,memory_bytes,allocs,result,samples,julia")
    for r in ROWS
        println(io, join((r[1], "\"" * r[2] * "\"", @sprintf("%.6f", r[3]), r[4], r[5],
                          @sprintf("%.15g", r[6]), r[7], VERSION), ","))
    end
end
println("wrote ", OUT)
