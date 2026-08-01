# common.jl — shared definitions for the bootstrap validation study
#
# The study asks a simple question: when a package quotes a 1σ uncertainty on a
# fitted parameter, is that number right?  The only way to know is to generate
# many independent realisations of the same observation ("universes"), fit each
# one, and compare the *actual* scatter of the best-fit parameters across
# universes with the uncertainty each method reports from a single universe.
#
#   sigma_true(p)   = std over universes of the best-fit value of p
#   sigma_method(p) = mean over universes of what the method reports
#   ratio           = sigma_method / sigma_true      (1 = calibrated)
#   coverage        = fraction of universes with |p_fit - p_true| < sigma_method
#                     (0.683 for a calibrated 1σ)
#
# Regimes, each stressing a different assumption:
#
#   ideal          error bars correct, noise independent.  Control: every
#                  method should be calibrated here.
#   systematics    same statistical noise, plus a calibration offset shared by
#                  all wavelengths of a given (epoch, baseline/triangle).  The
#                  quoted error bars know nothing about it.
#   underestimated noise generated at 2x the quoted error bars.
#   fewblocks      ideal noise, but a single snapshot: ~13 resampling blocks,
#                  where a block bootstrap is expected to struggle.
#   mismatch       the truth contains a third component that neither package
#                  fits.  No error estimate can be right about the truth here;
#                  the question is how each behaves.
#
# The uv coverage, MJDs and error bars come from a real CHARA/MIRC dataset
# (demos/data/iota_peg6t.oifits); only the observable values are replaced.

using OITOOLS
using OIFITS
using Random, Printf, Statistics, DelimitedFiles

# ── Model ────────────────────────────────────────────────────────────────────
# A binary at CHARA/H-band resolution (λ/2B ≈ 0.47 mas): partially resolved
# primary + unresolved companion at ~1.4 mas.

# The fluxes are constrained to sum to 1 ("A,f" = 1 - "B,f").  This matters for
# the head-to-head: PMOIRED normalises the model visibility by the total flux
# and OITOOLS does not, so the two packages fit the *same* function only when
# the total flux is pinned to 1.  With this constraint the two agree to
# chi2r = 0.003 on noiseless data (see README).

const TRUTH = Dict{String,Any}(
    "A,ud" => 0.80,
    "A,f"  => "1 - \$B,f",
    "B,ud" => 0.0,
    "B,f"  => 0.05,
    "B,x"  => 1.20,
    "B,y"  => -0.80,
)

# Extra component present in the `mismatch` truth only (undetected third star)
const TRUTH_MISMATCH = Dict{String,Any}(
    "A,ud" => 0.80,
    "A,f"  => "1 - \$B,f - \$C,f",
    "B,ud" => 0.0,
    "B,f"  => 0.05,
    "B,x"  => 1.20,
    "B,y"  => -0.80,
    "C,ud" => 0.0,
    "C,f"  => 0.012,
    "C,x"  => -2.0,
    "C,y"  => 1.50,
)

const FIT_MODEL = copy(TRUTH)

const FREE = ["A,ud", "B,f", "B,x", "B,y"]
# generous bounds, but B,x > 0 to exclude the (x,y) -> (-x,-y) mirror solution
const LB   = Dict("A,ud" => 0.05, "B,f" => 0.002, "B,x" => 0.20, "B,y" => -3.0)
const UB   = Dict("A,ud" => 3.00, "B,f" => 0.300, "B,x" => 3.00, "B,y" =>  0.2)
const WEIGHTS = [1.0, 0.0, 1.0]      # V2 and T3PHI, no T3AMP

const REGIMES = ["ideal", "systematics", "underestimated", "fewblocks", "mismatch"]

const SRC_TEMPLATE = joinpath(@__DIR__, "..", "data", "iota_peg6t.oifits")

# Correlated-calibration amplitudes for the `systematics` regime
const V2_SYS_RMS    = 0.05     # 5% multiplicative on V2, per (epoch, baseline)
const T3PHI_SYS_RMS = 1.5      # 1.5 degrees on T3PHI, per (epoch, triangle)

# ── Templates ────────────────────────────────────────────────────────────────

"""
    make_template_few(src, out) -> String

Single-snapshot version of `src`: keep only the rows of the earliest MJD.
"""
function make_template_few(src::String, out::String)
    ds = OIFITS.OIDataSet(src)
    mjd0 = minimum(minimum(db.mjd) for db in ds.vis2)
    keep(mjd) = abs(mjd - mjd0) < 1e-4

    for db in ds.vis2
        k = findall(keep, db.mjd)
        db.target_id = db.target_id[k]; db.time = db.time[k]; db.mjd = db.mjd[k]
        db.int_time  = db.int_time[k]
        db.vis2data  = db.vis2data[:, k]; db.vis2err = db.vis2err[:, k]
        db.ucoord    = db.ucoord[k]; db.vcoord = db.vcoord[k]
        db.sta_index = db.sta_index[:, k]; db.flag = db.flag[:, k]
    end
    for db in ds.t3
        k = findall(keep, db.mjd)
        db.target_id = db.target_id[k]; db.time = db.time[k]; db.mjd = db.mjd[k]
        db.int_time  = db.int_time[k]
        db.t3amp     = db.t3amp[:, k]; db.t3amperr = db.t3amperr[:, k]
        db.t3phi     = db.t3phi[:, k]; db.t3phierr = db.t3phierr[:, k]
        db.u1coord   = db.u1coord[k]; db.v1coord = db.v1coord[k]
        db.u2coord   = db.u2coord[k]; db.v2coord = db.v2coord[k]
        db.sta_index = db.sta_index[:, k]; db.flag = db.flag[:, k]
    end

    OIFITS.write(out, ds; overwrite=true)
    return out
end

# ── Writing a realisation ────────────────────────────────────────────────────

"""
    write_realization(template, outfile, v2, v2err, t3amp, t3amperr, t3phi, t3phierr)

Copy `template` and replace the OI_VIS2 / OI_T3 payload.  The flat vectors must
follow the ordering used by `readoifits`, i.e. table by table, each reshaped as
(nwavelength, nrow) — the same walk as `simulate_from_oifits`.
`check_ordering` below verifies that this assumption holds.
"""
function write_realization(template::String, outfile::String,
                           v2, v2err, t3amp, t3amperr, t3phi, t3phierr)
    ds = OIFITS.OIDataSet(template)

    off = 0
    for db in ds.vis2
        nrows = length(db.ucoord); nwavs = size(db.vis2data, 1); n = nwavs * nrows
        db.vis2data = reshape(Float64.(v2[off+1:off+n]),    nwavs, nrows)
        db.vis2err  = reshape(Float64.(v2err[off+1:off+n]), nwavs, nrows)
        db.flag     = falses(nwavs, nrows)
        off += n
    end
    off == length(v2) || error("write_realization: V2 length mismatch ($off vs $(length(v2)))")

    off = 0
    for db in ds.t3
        nrows = length(db.u1coord); nwavs = size(db.t3amp, 1); n = nwavs * nrows
        db.t3amp    = reshape(Float64.(t3amp[off+1:off+n]),    nwavs, nrows)
        db.t3amperr = reshape(Float64.(t3amperr[off+1:off+n]), nwavs, nrows)
        db.t3phi    = reshape(Float64.(t3phi[off+1:off+n]),    nwavs, nrows)
        db.t3phierr = reshape(Float64.(t3phierr[off+1:off+n]), nwavs, nrows)
        db.flag     = falses(nwavs, nrows)
        off += n
    end
    off == length(t3phi) || error("write_realization: T3 length mismatch")

    OIFITS.write(outfile, ds; overwrite=true)
    return outfile
end

"""
    read_template(file) -> OIdata

Read a template without dropping anything, so that the flat arrays map 1:1 onto
the file layout, and repair error bars that are missing or non-positive (the
simulated universes overwrite every value anyway, and carry no flags).
"""
function read_template(file::String)
    d = readoifits(file; T=Float64, filter_bad_data=false)[1,1]
    fix!(e) = begin
        good = filter(x -> isfinite(x) && x > 0, e)
        med  = isempty(good) ? 1.0 : median(good)
        for i in eachindex(e)
            (isfinite(e[i]) && e[i] > 0) || (e[i] = med)
        end
    end
    fix!(d.v2_err); fix!(d.t3amp_err); fix!(d.t3phi_err)
    # values themselves are replaced by the model, but NaNs would poison the
    # systematics multiplication
    for v in (d.v2, d.t3amp, d.t3phi)
        for i in eachindex(v); isfinite(v[i]) || (v[i] = 0.0); end
    end
    d.nv2 = length(d.v2); d.nt3amp = length(d.t3amp); d.nt3phi = length(d.t3phi)
    return d
end

"""
    check_ordering(template, workfile) -> NamedTuple

Round-trip test: write the noiseless model through `write_realization`, read it
back with `readoifits`, and compare with the model evaluated on the template.
Any mismatch in the flat-array ordering shows up here.
"""
function check_ordering(template::String, workfile::String)
    data = read_template(template)
    m    = dict_to_model(TRUTH, String[])
    obs  = model_to_obs(m, Float64[], data)
    write_realization(template, workfile,
                      obs.v2, data.v2_err, obs.t3amp, data.t3amp_err,
                      obs.t3phi, data.t3phi_err)
    back = read_template(workfile)
    return (v2 = maximum(abs.(back.v2 .- obs.v2)),
            t3phi = maximum(abs.(back.t3phi .- obs.t3phi)),
            nv2 = back.nv2, nt3 = back.nt3phi)
end

# ── Noise recipes ────────────────────────────────────────────────────────────

"""
    make_universe(data, blocks, regime, rng)

Model observables for `regime` plus noise according to that regime's recipe.
The returned error bars are always the *quoted* ones from the template: what
changes between regimes is how the actual noise relates to them.
"""
function make_universe(data::OIdata, blocks::DataBlocks, regime::String,
                       rng::AbstractRNG)
    truth_dict = regime == "mismatch" ? TRUTH_MISMATCH : TRUTH
    m   = dict_to_model(truth_dict, String[])
    obs = model_to_obs(m, Float64[], data)

    v2       = Float64.(collect(obs.v2))
    t3amp    = Float64.(collect(obs.t3amp))
    t3phi    = Float64.(collect(obs.t3phi))
    v2err    = Float64.(collect(data.v2_err))
    t3amperr = Float64.(collect(data.t3amp_err))
    t3phierr = Float64.(collect(data.t3phi_err))

    # correlated calibration errors: one draw per (epoch, baseline/triangle),
    # shared by every wavelength channel of that block — applied to the noiseless
    # observables, as a real transfer-function error would be
    if regime == "systematics"
        for b in 1:length(blocks)
            iv2 = blocks.idx_v2[b]
            isempty(iv2) || (v2[iv2] .*= (1 + V2_SYS_RMS * randn(rng)))
            it3 = blocks.idx_t3[b]
            isempty(it3) || (t3phi[it3] .+= T3PHI_SYS_RMS * randn(rng))
        end
    end

    # statistical noise; the `underestimated` regime draws it twice too large
    scale = regime == "underestimated" ? 2.0 : 1.0
    v2    .+= scale .* v2err    .* randn(rng, length(v2))
    t3amp .+= scale .* t3amperr .* randn(rng, length(t3amp))
    t3phi .+= scale .* t3phierr .* randn(rng, length(t3phi))

    return v2, v2err, t3amp, t3amperr, t3phi, t3phierr
end

# ── Bookkeeping ──────────────────────────────────────────────────────────────

truth_values() = [Float64(TRUTH[p]) for p in FREE]

function write_table(path::String, header::Vector{String}, rows)
    open(path, "w") do io
        println(io, join(header, "\t"))
        for r in rows
            println(io, join([x isa Real ? (isfinite(x) ? (@sprintf "%.10g" x) : "NaN") : string(x)
                              for x in r], "\t"))
        end
    end
    return path
end
