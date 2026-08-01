# make_universes.jl — generate the simulated datasets for the validation study
#
#   julia --project=../.. make_universes.jl [outdir] [nuniverses]
#
# Writes <outdir>/<regime>/universe_XXXX.oifits, one independent realisation of
# the same observation per file.

include("common.jl")

const OUTDIR = length(ARGS) >= 1 ? ARGS[1] :
               joinpath(tempdir(), "oitools_bootstrap_validation")
const NUNI   = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 100

mkpath(OUTDIR)
println("Output directory:     ", OUTDIR)
println("Universes per regime: ", NUNI)

templates = Dict(
    :many => joinpath(OUTDIR, "template_many.oifits"),
    :few  => joinpath(OUTDIR, "template_few.oifits"),
)
isfile(templates[:many]) || cp(SRC_TEMPLATE, templates[:many])
isfile(templates[:few])  || make_template_few(SRC_TEMPLATE, templates[:few])

# The flat-array ordering assumed by write_realization must match readoifits
for (k, f) in templates
    c = check_ordering(f, joinpath(OUTDIR, "roundtrip_check.oifits"))
    @printf("template %-4s : nv2=%4d nt3=%4d  round-trip error: V2 %.2e  T3PHI %.2e\n",
            string(k), c.nv2, c.nt3, c.v2, c.t3phi)
    (c.v2 < 1e-6 && c.t3phi < 1e-4) ||
        error("round-trip check failed for $f — flat array ordering mismatch")
end
rm(joinpath(OUTDIR, "roundtrip_check.oifits"), force=true)

for regime in REGIMES
    tmpl = regime == "fewblocks" ? templates[:few] : templates[:many]
    data = read_template(tmpl)
    blk  = data_blocks(data; granularity=:config)
    dir  = joinpath(OUTDIR, regime)
    mkpath(dir)
    @printf("%-15s : %3d blocks, %4d V2 + %4d T3PHI points\n",
            regime, length(blk), data.nv2, data.nt3phi)
    for u in 1:NUNI
        out = joinpath(dir, @sprintf("universe_%04d.oifits", u))
        isfile(out) && continue
        rng = Xoshiro(hash((regime, u)))
        v2, v2e, t3a, t3ae, t3p, t3pe = make_universe(data, blk, regime, rng)
        write_realization(tmpl, out, v2, v2e, t3a, t3ae, t3p, t3pe)
    end
end

write_table(joinpath(OUTDIR, "truth.tsv"), FREE, [truth_values()])
println("\nTruth: ", join([@sprintf("%s=%.4f", p, v) for (p, v) in zip(FREE, truth_values())], "  "))
println("Done.")
