#!/usr/bin/env julia
# test_readoifits.jl
#
# Usage:
#   julia test_readoifits.jl /path/to/oifits/dir          # recurse into dir
#   julia test_readoifits.jl file1.oifits file2.fits ...  # explicit files
#
# Output: a table to stdout + a CSV log at test_readoifits_results.csv
#
# For each file it attempts:
#   1. list_oifits_targets()          — pure I/O, tests OIFITS loading
#   2. readoifits()                   — full pipeline, first target, default options
#   3. Sanity checks on OIdata:
#        UV: nuv>0, no NaN, index bounds
#        V2: no NaN in values/errors, index bounds
#        T3: no NaN in amp/phi, index bounds
#        VIS: no NaN in visamp/visphi, no non-positive errors, index bounds
#        FLUX: no NaN in flux/flux_err, flux_sta_index length consistent
#   4. readoifits(T=Float32)          — counts + value closeness vs Float64
#   5. readoifits(special_filter_diffvis=true) — if file has OI_VIS data

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))   # adjust if running from a different dir

# ---- Load OITOOLS ----------------------------------------------------------
include(joinpath(@__DIR__, "../src/OITOOLS.jl"))
using Main.OITOOLS

using Printf, Dates

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

struct FileResult
    path          :: String
    targets       :: Vector{String}
    nv2           :: Int
    nt3amp        :: Int
    nt3phi        :: Int
    nvis          :: Int
    nflux         :: Int
    nuv           :: Int
    mean_mjd      :: Float64
    f32_ok        :: Bool
    diffvis_ok    :: Union{Bool, Nothing}   # nothing = no OI_VIS in file
    warnings      :: Vector{String}
    error         :: Union{String, Nothing}
    elapsed_s     :: Float64
end

# Float32 casting tolerance: values should agree to ~1e-5 relative
const F32_RTOL = 1e-5

function sanity_check(data::OITOOLS.OIdata{T}) where T
    issues = String[]

    # ---- UV plane -----------------------------------------------------------
    data.nuv == 0 && push!(issues, "nuv=0 (no UV points)")
    data.nuv > 0  && any(isnan.(data.uv))          && push!(issues, "NaN in uv matrix")

    # ---- V2 -----------------------------------------------------------------
    if data.nv2 > 0
        any(isnan.(data.v2))                        && push!(issues, "NaN in v2 (unflagged)")
        any(isnan.(data.v2_err))                    && push!(issues, "NaN in v2_err")
        any(data.v2_err .<= 0)                      && push!(issues, "non-positive v2_err")
        maximum(data.indx_v2) > data.nuv            && push!(issues, "indx_v2 out of range")
    end

    # ---- T3 -----------------------------------------------------------------
    if data.nt3amp > 0 || data.nt3phi > 0
        data.nt3amp > 0 && any(isnan.(data.t3amp))  && push!(issues, "NaN in t3amp (unflagged)")
        data.nt3amp > 0 && any(isnan.(data.t3amp_err)) && push!(issues, "NaN in t3amp_err")
        data.nt3phi > 0 && any(isnan.(data.t3phi))  && push!(issues, "NaN in t3phi (unflagged)")
        data.nt3phi > 0 && any(isnan.(data.t3phi_err)) && push!(issues, "NaN in t3phi_err")
        maximum(data.indx_t3_1) > data.nuv          && push!(issues, "indx_t3_1 out of range")
        maximum(data.indx_t3_2) > data.nuv          && push!(issues, "indx_t3_2 out of range")
        maximum(data.indx_t3_3) > data.nuv          && push!(issues, "indx_t3_3 out of range")
    end

    # ---- OI_VIS -------------------------------------------------------------
    if data.nvisamp > 0
        any(isnan.(data.visamp))                    && push!(issues, "NaN in visamp (unflagged)")
        any(isnan.(data.visamp_err))                && push!(issues, "NaN in visamp_err")
        any(data.visamp_err .<= 0)                  && push!(issues, "non-positive visamp_err")
    end
    if data.nvisphi > 0
        any(isnan.(data.visphi))                    && push!(issues, "NaN in visphi (unflagged)")
        any(isnan.(data.visphi_err))                && push!(issues, "NaN in visphi_err")
        any(data.visphi_err .<= 0)                  && push!(issues, "non-positive visphi_err")
    end
    if data.nvisamp > 0 || data.nvisphi > 0
        !isempty(data.indx_vis) && maximum(data.indx_vis) > data.nuv &&
            push!(issues, "indx_vis out of range")
    end

    # ---- OI_FLUX ------------------------------------------------------------
    if data.nflux > 0
        any(isnan.(data.flux))                      && push!(issues, "NaN in flux (unflagged)")
        any(isnan.(data.flux_err))                  && push!(issues, "NaN in flux_err")
        any(data.flux_err .<= 0)                    && push!(issues, "non-positive flux_err")
        length(data.flux_sta_index) != data.nflux   &&
            push!(issues, "flux_sta_index length $(length(data.flux_sta_index)) ≠ nflux $(data.nflux)")
    end

    return issues
end

# Check Float32 values are close to Float64 reference (within casting tolerance)
function check_f32_values(d64::OITOOLS.OIdata{Float64}, d32::OITOOLS.OIdata{Float32})
    issues = String[]
    function relerr(a64, a32, name)
        isempty(a64) && return
        maxerr = maximum(abs.(a64 .- Float64.(a32)) ./ (abs.(a64) .+ 1e-30))
        maxerr > F32_RTOL && push!(issues, "Float32 $name max rel error = $(@sprintf("%.2e", maxerr))")
    end
    d64.nv2    > 0 && (relerr(d64.v2,    d32.v2,    "v2");
                       relerr(d64.v2_err, d32.v2_err, "v2_err"))
    d64.nt3amp > 0 && (relerr(d64.t3amp, d32.t3amp, "t3amp");
                       relerr(d64.t3amp_err, d32.t3amp_err, "t3amp_err"))
    d64.nt3phi > 0 &&  relerr(d64.t3phi, d32.t3phi, "t3phi")
    d64.nvisamp > 0 && relerr(d64.visamp, d32.visamp, "visamp")
    d64.nvisphi > 0 && relerr(d64.visphi, d32.visphi, "visphi")
    d64.nflux  > 0 &&  relerr(d64.flux,  d32.flux,  "flux")
    return issues
end

function test_file(path::String) :: FileResult
    warnings = String[]
    t0 = time()

    # Step 1: list targets
    targets = String[]
    try
        targets = OITOOLS.list_oifits_targets(path)
    catch e
        return FileResult(path, targets, 0, 0, 0, 0, 0, 0, 0.0, false, nothing,
                          warnings, "list_targets failed: $(sprint(showerror, e))",
                          time() - t0)
    end

    isempty(targets) && return FileResult(path, targets, 0, 0, 0, 0, 0, 0, 0.0, false, nothing,
                                          warnings, "no targets found", time() - t0)

    targetname = targets[1]

    # Step 2: readoifits Float64
    data = nothing
    try
        result = readoifits(path; targetname, warn=false)
        data = result[1, 1]
    catch e
        return FileResult(path, targets, 0, 0, 0, 0, 0, 0, 0.0, false, nothing,
                          warnings, "readoifits failed: $(sprint(showerror, e))",
                          time() - t0)
    end

    # Step 3: sanity checks
    append!(warnings, sanity_check(data))

    # Step 4: Float32 — counts + value closeness
    f32_ok = false
    try
        r32  = readoifits(path; targetname, warn=false, T=Float32)
        d32  = r32[1, 1]
        count_ok = (d32.nv2 == data.nv2) && (d32.nt3amp == data.nt3amp) &&
                   (d32.nvisphi == data.nvisphi) && (d32.nflux == data.nflux)
        if !count_ok
            push!(warnings, "Float32 count mismatch: " *
                  "nv2 $(d32.nv2)/$(data.nv2), nt3amp $(d32.nt3amp)/$(data.nt3amp), " *
                  "nvisphi $(d32.nvisphi)/$(data.nvisphi), nflux $(d32.nflux)/$(data.nflux)")
        end
        val_issues = check_f32_values(data, d32)
        append!(warnings, val_issues)
        f32_ok = count_ok && isempty(val_issues)
    catch e
        push!(warnings, "Float32 path failed: $(sprint(showerror, e))")
    end

    # Step 5: differential visibility filter (only if file has OI_VIS)
    diffvis_ok = nothing
    if data.nvisamp > 0 || data.nvisphi > 0
        diffvis_ok = false
        try
            rd = readoifits(path; targetname, warn=false, special_filter_diffvis=true)
            dd = rd[1, 1]
            # filtered count must be ≤ original (cross-bin intersection can only shrink it)
            if dd.nvisphi > data.nvisphi
                push!(warnings, "diffvis filter increased nvisphi $(data.nvisphi)→$(dd.nvisphi)")
            else
                diffvis_ok = true
            end
        catch e
            push!(warnings, "special_filter_diffvis failed: $(sprint(showerror, e))")
        end
    end

    return FileResult(path, targets,
                      data.nv2, data.nt3amp, data.nt3phi,
                      data.nvisamp, data.nflux, data.nuv, data.mean_mjd,
                      f32_ok, diffvis_ok, warnings, nothing, time() - t0)
end

# ---------------------------------------------------------------------------
# Collect files
# ---------------------------------------------------------------------------

args = isempty(ARGS) ? ["."] : ARGS
oifits_files = String[]
for arg in args
    if isdir(arg)
        for (root, _, files) in walkdir(arg)
            for f in files
                if any(endswith(f, ext) for ext in (".oifits", ".fits", ".fit", ".OIFITS", ".FITS"))
                    push!(oifits_files, joinpath(root, f))
                end
            end
        end
    elseif isfile(arg)
        push!(oifits_files, arg)
    else
        @warn "Argument '$arg' is neither a file nor a directory, skipping."
    end
end

isempty(oifits_files) && (println("No OIFITS files found."); exit(0))
sort!(oifits_files)
println("Found $(length(oifits_files)) file(s) to test.\n")

# ---------------------------------------------------------------------------
# Run tests
# ---------------------------------------------------------------------------

results = FileResult[]
for (i, path) in enumerate(oifits_files)
    print(@sprintf("[%3d/%3d] %-60s ", i, length(oifits_files), basename(path)))
    flush(stdout)
    r = test_file(path)
    push!(results, r)
    if isnothing(r.error)
        status = isempty(r.warnings) ? "✓ OK" : "⚠ WARN"
        diffvis_str = isnothing(r.diffvis_ok) ? "" : (r.diffvis_ok ? " diffvis✓" : " diffvis✗")
        println(@sprintf("%s  nV2=%-5d nT3=%-5d nVIS=%-5d nFlux=%-5d nuv=%-5d%s  %.2fs",
                         status, r.nv2, r.nt3amp, r.nvis, r.nflux, r.nuv,
                         diffvis_str, r.elapsed_s))
        for w in r.warnings
            println(@sprintf("          ↳ %s", w))
        end
    else
        println(@sprintf("✗ FAIL  %.2fs", r.elapsed_s))
        println(@sprintf("          ↳ %s", r.error))
    end
end

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

n_ok         = count(r -> isnothing(r.error) && isempty(r.warnings), results)
n_warn       = count(r -> isnothing(r.error) && !isempty(r.warnings), results)
n_fail       = count(r -> !isnothing(r.error), results)
n_f32_fail   = count(r -> isnothing(r.error) && !r.f32_ok, results)
n_diffvis_fail = count(r -> isnothing(r.error) && r.diffvis_ok === false, results)

println("\n" * "="^70)
println(@sprintf("RESULTS:  %d passed  |  %d warned  |  %d failed  |  %d Float32 issues  |  %d diffvis issues",
                 n_ok, n_warn, n_fail, n_f32_fail, n_diffvis_fail))
println("="^70)

if n_fail > 0
    println("\nFailed files:")
    for r in results
        isnothing(r.error) && continue
        println("  $(r.path)")
        println("    $(r.error)")
    end
end

if n_warn > 0
    println("\nFiles with warnings:")
    for r in results
        (isnothing(r.error) && !isempty(r.warnings)) || continue
        println("  $(r.path)  [targets: $(join(r.targets, ", "))]")
        for w in r.warnings
            println("    ↳ $w")
        end
    end
end

# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------

csv_path = joinpath(dirname(first(oifits_files)), "test_readoifits_results.csv")
open(csv_path, "w") do io
    println(io, "path,targets,nv2,nt3amp,nt3phi,nvisamp,nflux,nuv,mean_mjd,f32_ok,diffvis_ok,warnings,error,elapsed_s")
    for r in results
        targets_str  = join(r.targets, ";")
        warn_str     = join(r.warnings, "; ")
        err_str      = isnothing(r.error) ? "" : r.error
        diffvis_str  = isnothing(r.diffvis_ok) ? "n/a" : string(r.diffvis_ok)
        println(io, join([r.path, targets_str, r.nv2, r.nt3amp, r.nt3phi,
                          r.nvis, r.nflux, r.nuv,
                          @sprintf("%.5f", r.mean_mjd),
                          r.f32_ok, diffvis_str, warn_str, err_str,
                          @sprintf("%.3f", r.elapsed_s)], ","))
    end
end
println("\nResults written to: $csv_path")
