# Two benchmark runs side by side (CSV files from suite.jl and ttfx.jl): the median time of
# each row in both, their ratio (b / a), and the result values where they differ.
#
#     julia test/bench/compare.jl a.csv b.csv
using Printf

function load(f)
    rows = Dict{Tuple{String,String},Vector{Tuple{Float64,Float64,Float64,String}}}()
    order = Tuple{String,String}[]
    for l in eachline(f)
        startswith(l, "group,") && continue
        m = match(r"^([^,]*),\"(.*)\",([^,]*),([^,]*),([^,]*),([^,]*),([^,]*),(.*)$", l)
        m === nothing && continue
        k = (String(m[1]), String(m[2]))
        haskey(rows, k) || (rows[k] = []; push!(order, k))
        push!(rows[k], (parse(Float64, m[3]), parse(Float64, m[4]), parse(Float64, m[6]), String(m[8])))
    end
    return rows, order
end
med(v) = (s = sort(v); n = length(s); isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2)
fmt(ms) = ms >= 1000 ? @sprintf("%.2f s", ms / 1000) : ms >= 1 ? @sprintf("%.2f ms", ms) : @sprintf("%.1f μs", 1000ms)

a, oa = load(ARGS[1]); b, _ = load(ARGS[2])
va = first(a[oa[1]])[4]; vb = first(b[first(keys(b))])[4]
@printf("%-9s %-52s %12s %12s %7s %9s  %s\n", "group", "task", va, vb, "ratio", "mem ratio", "result")
for k in oa
    haskey(b, k) || continue
    ta = med(first.(a[k])); tb = med(first.(b[k]))
    ma = a[k][1][2]; mb = b[k][1][2]
    ra = a[k][1][3]; rb = b[k][1][3]
    same = (isnan(ra) && isnan(rb)) || isapprox(ra, rb; rtol = 1e-10)
    flag = tb / ta > 1.10 ? "  slower" : tb / ta < 0.90 ? "  faster" : ""
    @printf("%-9s %-52s %12s %12s %7.3f %9s  %s%s\n", k[1], k[2], fmt(ta), fmt(tb), tb / ta,
            ma > 0 ? @sprintf("%.3f", mb / ma) : "", same ? "same" : @sprintf("%.10g → %.10g", ra, rb), flag)
end
