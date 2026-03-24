#
# Architecture:
#   load_oifits(file)            → OIDataSet   (I/O only, try/catch for v1 files)
#   collect_station_info(...)    → NamedTuple  (station name/index mapping)
#   read_flux_tables(...)        → NamedTuple  (flat arrays, all obs merged)
#   read_vis_tables(...)         → NamedTuple
#   read_v2_tables(...)          → NamedTuple
#   read_t3_tables(...)          → NamedTuple
#   BinData{T}                   → mutable struct for per-bin workspace
#   slice_to_bin(...)            → BinData     (select + assemble UV plane)
#   filter_bad_observables!(bd)  → BinData     (filter + UV prune, in-place)
#   remove_redundant_uv!(bd)     → BinData     (UV deduplication, in-place)
#   make_oidata(bd, ...)         → OIdata{T}   (final packaging)
#   readoifits(file; ...)        → Array{OIdata{T}}
#
# emmt OIFITS.jl v1 API:
#   ds = read(OIDataSet, file)   or   ds = OIDataSet(file)
#   ds.instr  → Vector{OI_WAVELENGTH};  ds.array  → Vector{OI_ARRAY}
#   ds.vis2   → Vector{OI_VIS2};        ds.t3     → Vector{OI_T3}
#   ds.vis    → Vector{OI_VIS};         ds.flux   → Vector{OI_FLUX}
#   ds.target → single OI_TARGET (indexable as vector of OITargetEntry)
#     tgt[i].target_id,  tgt[i].target  (not array-valued columns)
#
# OIdata{T}: numeric arrays stored as type T (default Float64).
# Use T=Float32 for ~50 % memory reduction.

using OIFITS, FITSIO, Statistics, SparseArrays, Crayons
using NearestNeighbors

# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

function rm_redundance_kdtree(uv, uvtol)
    indx_redundance = collect(1:size(uv,2))
    kdtree = KDTree(uv)
    for value in indx_redundance
        redundance = inrange(kdtree, uv[:,value], uvtol)
        @inbounds indx_redundance[redundance] .= minimum(redundance)
    end
    tokeep = unique(indx_redundance)
    indx_red_conv = indexin(indx_redundance, tokeep)
    return indx_red_conv, tokeep
end

# Find the OI_WAVELENGTH block whose insname matches db.insname
function _find_wav(db, wavtables, wavtableref)
    idx = findfirst(==(db.insname), wavtableref)
    isnothing(idx) && error("OI_WAVELENGTH table not found for insname=$(db.insname)")
    return wavtables[idx]
end

"""
    merge_wav_tables(wavtables; T=Float64, verbose=true)

Check whether multiple OI_WAVELENGTH tables share compatible spectral channels
(i.e. their bandpass intervals overlap channel-by-channel) and return a single
deduplicated wavelength vector.

Returns `eff_wave_merged::Vector{T}` — the sorted effective wavelengths from the
first table (used as reference).

Throws an error if the tables have different channel counts or if any channel
pair has no bandpass overlap.
"""
function merge_wav_tables(wavtables; T::Type{<:AbstractFloat}=Float64, verbose=true)
    ref = wavtables[1]
    ref_lo = T.(ref.eff_wave) .- T.(ref.eff_band) ./ 2
    ref_hi = T.(ref.eff_wave) .+ T.(ref.eff_band) ./ 2
    for k in 2:length(wavtables)
        wt = wavtables[k]
        if length(wt.eff_wave) != length(ref.eff_wave)
            error("merge_oi_wavelength: OI_WAVELENGTH table $k has $(length(wt.eff_wave)) channels vs $(length(ref.eff_wave)) in table 1 — tables are incompatible.")
        end
        wt_lo = T.(wt.eff_wave) .- T.(wt.eff_band) ./ 2
        wt_hi = T.(wt.eff_wave) .+ T.(wt.eff_band) ./ 2
        for j in eachindex(ref_lo)
            if wt_lo[j] >= ref_hi[j] || wt_hi[j] <= ref_lo[j]
                error("merge_oi_wavelength: channel $j has no bandpass overlap between table 1 ($(ref_lo[j])–$(ref_hi[j]) m) and table $k ($(wt_lo[j])–$(wt_hi[j]) m) — tables are incompatible.")
            end
        end
    end
    verbose && printstyled("  merge_oi_wavelength: $(length(wavtables)) tables share $(length(ref.eff_wave)) compatible channels — merged.\n", color=:green)
    return sort(T.(ref.eff_wave))
end

# Convert per-baseline (u, v) coordinates to spatial frequencies u/λ, v/λ.
# coords: vector of length nbaselines, each element a scalar coordinate in metres
# lam:    wavelength vector of length nwave (metres)
# Returns flat vectors of length nbaselines×nwave, with baselines varying slowest.
function _uv_lambda(ucoords::AbstractVector{T}, vcoords::AbstractVector{T},
                    lam::AbstractVector{T}) where T
    u = vcat([ucoords[i] ./ lam for i in eachindex(ucoords)]...)
    v = vcat([vcoords[i] ./ lam for i in eachindex(vcoords)]...)
    return u, v
end

# ---------------------------------------------------------------------------
# OI_CORR helpers
# ---------------------------------------------------------------------------

# Build a symmetric sparse correlation matrix from an OI_CORR table entry.
# OI_CORR stores only the upper-triangle off-diagonal (j > i per the standard);
# diagonal is implicitly 1.  We explicitly store both triangles so the sparse
# matrix is symmetric.  (Symmetric wrapper deferred to Cholesky decomposition.)
function _build_corr_sparse(ct; T::Type{<:AbstractFloat}=Float64)
    n    = ct.ndata
    iis  = vcat(Int64.(ct.iindx), Int64.(ct.jindx), Int64.(1:n))
    jjs  = vcat(Int64.(ct.jindx), Int64.(ct.iindx), Int64.(1:n))
    vals = vcat(T.(ct.corr),      T.(ct.corr),      ones(T, n))
    return sparse(iis, jjs, vals, n, n)
end

# Safely access a named property from an OIFITS table object.
# Returns nothing if the field is missing, uninitialized, or otherwise errors.
_safe_field(obj, field::Symbol) = try; getproperty(obj, field); catch; nothing; end

# Expand per-observation corrindx to a flat per-(wave,obs) index vector.
# Per OIFITSv2 §7.2: index for channel j = CORRINDX + j − 1.
# corrindx_obs[r] = 1-based start in OI_CORR matrix for obs r (0 = no corr).
# offset: added to all non-zero indices for block-diagonal assembly.
# corr_ndata: size of the OI_CORR matrix — used for bounds validation.
# mjd_obs: optional MJD vector (one per obs) to detect multi-baseline timestamp
#          groups (GRAVITY-style files where CORRINDX=1 for all rows but
#          corr_ndata = nbase × nwave).
function _expand_corrindx(corrindx_obs, nwave::Int, offset::Int, corr_ndata::Int;
                          mjd_obs::Union{Nothing,AbstractVector}=nothing)
    nobs   = length(corrindx_obs)
    result = zeros(Int64, nwave * nobs)

    # Detect multi-baseline correlation blocks:
    # When corr_ndata > nwave and is an exact multiple, multiple rows within
    # one timestamp collectively span the full block.  GRAVITY writes
    # CORRINDX=1 for every row, so we must assign per-row offsets within
    # each timestamp group based on row order.
    nrows_per_group = corr_ndata ÷ nwave
    multi_baseline  = (nrows_per_group > 1 && corr_ndata == nrows_per_group * nwave
                       && !isnothing(mjd_obs))

    # Track position within each timestamp group (keyed by CORRINDX + MJD)
    group_counter = Dict{Tuple{Int64,Float64}, Int}()

    for obs in 1:nobs
        base = Int64(corrindx_obs[obs])
        base == 0 && continue

        # In multi-baseline mode, compute the intra-group row index (0-based)
        row_offset = 0
        if multi_baseline
            mjd_val = Float64(mjd_obs[obs])
            key = (base, mjd_val)
            row_offset = get(group_counter, key, 0)
            group_counter[key] = row_offset + 1
            if row_offset >= nrows_per_group
                @warn("CORRINDX expansion: timestamp group at MJD=$mjd_val has " *
                      "more than $nrows_per_group rows — skipping excess.")
                continue
            end
        end

        # Effective base index for this row
        effective_base = base + row_offset * nwave

        # Bounds check
        max_raw = effective_base + nwave - 1
        if max_raw > corr_ndata
            @warn("CORRINDX expansion out of bounds: obs $obs has " *
                  "effective_base=$effective_base, nwave=$nwave → " *
                  "max index $max_raw > NDATA=$corr_ndata. Skipping.")
            continue
        end
        effective_base += offset
        for w in 1:nwave
            result[(obs-1)*nwave + w] = effective_base + (w - 1)
        end
    end
    return result
end

# For one data table: produce the flat corr-index vector for nwave*nobs points.
# Uses a Dict{String,Int} (corrname_offsets) to map each CORRNAME to its stable
# block-diagonal offset, independent of the order data tables are encountered.
# Returns zeros(Int64, nwave*nobs) if correlation is unavailable for this table.
# mjd_col: optional MJD vector (one per obs in tid_ok) for multi-baseline detection.
function _process_corrindx(corrname_val, corrindx_col, tid_ok, nwave::Int,
                            correlt::Dict,
                            corrname_offsets::Dict{String,Int},
                            seen_corrnames::Vector{String};
                            mjd_obs::Union{Nothing,AbstractVector}=nothing)
    nobs = length(tid_ok)
    cn   = isnothing(corrname_val) ? "" : strip(string(corrname_val))
    (isempty(cn) || !haskey(correlt, cn) || isnothing(corrindx_col)) &&
        return zeros(Int64, nwave * nobs)
    cidx_obs = Int64.(corrindx_col[tid_ok])
    # Register this CORRNAME if first encounter; offset is determined by
    # the order CORRNAMEs appear in the OI_CORR tables (via correlt iteration
    # order), NOT by the order data tables happen to reference them.
    if !haskey(corrname_offsets, cn)
        push!(seen_corrnames, cn)
        corrname_offsets[cn] = sum(size(correlt[k], 1) for k in seen_corrnames[1:end-1]; init=0)
    end
    offset     = corrname_offsets[cn]
    corr_ndata = size(correlt[cn], 1)
    return _expand_corrindx(cidx_obs, nwave, offset, corr_ndata; mjd_obs)
end

# Build the final block-diagonal corr sparse matrix from the accumulated list.
function _assemble_corr(seen_corrnames::Vector{String}, correlt::Dict,
                        T::Type{<:AbstractFloat})
    isempty(seen_corrnames) && return spzeros(T, 0, 0)
    parts = [correlt[cn] for cn in seen_corrnames]
    length(parts) == 1 && return parts[1]
    return blockdiag(parts...)
end

# ---------------------------------------------------------------------------
# OIdata{T}: central data container
# ---------------------------------------------------------------------------
"""
    OIdata{T<:AbstractFloat}

Central data container produced by `readoifits`. All observable arrays are flat
vectors of length `nobs`; UV coordinates are stored in a `2×nuv` matrix.

# Observable arrays
- `v2`, `v2_err` — squared visibilities and errors
- `v2_baseline`, `v2_lam`, `v2_dlam`, `v2_mjd`, `v2_flag` — baseline (m), λ (m), Δλ (m), MJD, flag
- `t3phi`, `t3phi_err`, `t3amp`, `t3amp_err` — closure phases (deg) and triple amplitudes
- `t3_baseline` — geometric-mean baseline (m); `t3_maxbaseline` — longest side (m)
- `t3_lam`, `t3_dlam`, `t3_mjd`, `t3_flag`
- `visamp`, `visamp_err`, `visphi`, `visphi_err` — complex visibility amplitude and phase (deg)
- `vis_baseline`, `vis_lam`, `vis_dlam`, `vis_mjd`, `vis_flag`
- `flux`, `flux_err`, `flux_lam`, `flux_dlam`, `flux_mjd`, `flux_flag`
- `flux_sta_index` — station index for each flux point; `0` means calibrated (OI_FLUX CALSTAT=C)
- `flux_calibrated` — `true` if OI_FLUX has CALSTAT="C" (calibrated source spectrum / SED)

# UV plane
- `uv` — `2×nuv` matrix of (u, v) spatial frequencies in cycles/m (i.e. baseline/λ)
- `uv_lam`, `uv_dlam`, `uv_mjd`, `uv_baseline`
- `indx_v2`, `indx_vis` — index of each V²/vis point into the UV array
- `indx_t3_1`, `indx_t3_2`, `indx_t3_3` — UV indices for the three legs of each triangle

# Station / telescope metadata
- `sta_name`, `tel_name`, `sta_index` — station names, telescope names, station indices
- `v2_sta_index` — `2×nv2` matrix; `t3_sta_index` — `3×nt3` matrix; `vis_sta_index` — `2×nvis`

# Sizes
- `nv2`, `nt3amp`, `nt3phi`, `nvisamp`, `nvisphi`, `nflux`, `nuv`

# Correlation matrices (from OI_CORR; empty sparse matrices when absent)
- `v2_corr`, `v2_corr_idx` — V² correlation matrix and per-point 1-based index
- `t3amp_corr`, `t3amp_corr_idx`, `t3phi_corr`, `t3phi_corr_idx`
- `visamp_corr`, `visamp_corr_idx`, `visphi_corr`, `visphi_corr_idx`
- `flux_corr`, `flux_corr_idx`
- `*_corr_idx[i] == 0` means point `i` has no associated correlation row.

# Other
- `mean_mjd::Float64` — mean MJD of the bin (always `Float64` regardless of `T`)
- `filename` — path to the originating OIFITS file
"""
mutable struct OIdata{T<:AbstractFloat}
    # Complex visibilities
    visamp::Vector{T};             visamp_err::Vector{T}
    visphi::Vector{T};             visphi_err::Vector{T}
    vis_baseline::Vector{T};       vis_mjd::Vector{T}
    vis_lam::Vector{T};            vis_dlam::Vector{T}
    vis_flag::Vector{Bool}
    # V2
    v2::Vector{T};                 v2_err::Vector{T}
    v2_baseline::Vector{T};        v2_mjd::Vector{T}
    mean_mjd::Float64              # always Float64 — MJD precision matters
    v2_lam::Vector{T};             v2_dlam::Vector{T}
    v2_flag::Vector{Bool}
    # T3
    t3amp::Vector{T};              t3amp_err::Vector{T}
    t3phi::Vector{T};              t3phi_err::Vector{T}
    t3phi_vonmises_err::Vector{T}; t3phi_vonmises_chi2_offset::Vector{T}
    t3_baseline::Vector{T};        t3_maxbaseline::Vector{T}
    t3_mjd::Vector{T};             t3_lam::Vector{T};    t3_dlam::Vector{T}
    t3_flag::Vector{Bool}
    # OIFlux
    flux::Vector{T};               flux_err::Vector{T}
    flux_mjd::Vector{T};           flux_lam::Vector{T};  flux_dlam::Vector{T}
    flux_flag::Vector{Bool};       flux_sta_index::Vector{Int64}
    flux_calibrated::Bool          # true if CALSTAT="C" (calibrated source spectrum)
    # UV coverage (columns = uv points, rows = [u; v])
    uv::Matrix{T};                 uv_lam::Vector{T};    uv_dlam::Vector{T}
    uv_mjd::Vector{T};             uv_baseline::Vector{T}
    # Sizes
    nflux::Int64;  nvisamp::Int64; nvisphi::Int64
    nv2::Int64;    nt3amp::Int64;  nt3phi::Int64;  nuv::Int64
    # UV indices for each observable
    indx_vis::Vector{Int64};       indx_v2::Vector{Int64}
    indx_t3_1::Vector{Int64};      indx_t3_2::Vector{Int64}; indx_t3_3::Vector{Int64}
    # Station / telescope info (shared across all bins)
    sta_name::Vector{String};      tel_name::Vector{String}
    sta_index::Vector{Int64}
    vis_sta_index::Matrix{Int64};  v2_sta_index::Matrix{Int64}
    t3_sta_index::Matrix{Int64}
    # Correlation matrices from OI_CORR (empty sparse if absent).
    # *_corr: the raw OI_CORR sparse matrix (size ndata×ndata, per the standard).
    # *_corr_idx: per-data-point 1-based index into *_corr (0 = no correlation).
    v2_corr::SparseMatrixCSC{T,Int64};       v2_corr_idx::Vector{Int64}
    t3amp_corr::SparseMatrixCSC{T,Int64};    t3amp_corr_idx::Vector{Int64}
    t3phi_corr::SparseMatrixCSC{T,Int64};    t3phi_corr_idx::Vector{Int64}
    visamp_corr::SparseMatrixCSC{T,Int64};   visamp_corr_idx::Vector{Int64}
    visphi_corr::SparseMatrixCSC{T,Int64};   visphi_corr_idx::Vector{Int64}
    flux_corr::SparseMatrixCSC{T,Int64};     flux_corr_idx::Vector{Int64}
    # OI_VIS metadata (OIFITSv2 §6)
    amptyp::String                           # "absolute", "differential", "correlated flux", or ""
    phityp::String                           # "absolute", "differential", or ""
    filename::String
end

function Base.display(data::OIdata{T}) where T
    println("Mean MJD: $(data.mean_mjd)  [eltype: $T]")
    println("Wavelength range: $(minimum(data.uv_lam)) - $(maximum(data.uv_lam))")
    println("nflux: $(data.nflux)$(data.nflux > 0 ? (data.flux_calibrated ? " (calibrated)" : " (uncalibrated)") : "") | nuv: $(data.nuv) | nvisamp: $(data.nvisamp) | " *
            "nvisphi: $(data.nvisphi) | nv2: $(data.nv2) | nt3amp: $(data.nt3amp) | nt3phi: $(data.nt3phi)")
end

function Base.display(data::Array{<:OIdata})
    println("Original data file: $(data[1].filename)")
    nwav = size(data, 1)
    nepo = size(data, 2)
    println("Number of wavelength bins: $nwav")
    println("Number of time/epoch bins: $nepo")
    # Blue-to-red palette (short λ → long λ), interpolated for any nwav
    _palette = [20, 28, 112, 148, 220, 208, 166, 196]
    if nwav == 1
        printcolors = [_palette[end]]
    else
        printcolors = [_palette[1 + round(Int, (i-1)/(nwav-1) * (length(_palette)-1))] for i in 1:nwav]
    end
    if size(data) == (1,1)
        display(data[1,1])
    elseif nepo == 1
        for i in 1:nwav
            print(Crayon(foreground=printcolors[i], bold=true), "Wavelength: $i/$nwav\n")
            display(data[i,1])
        end
    else
        # Aggregate summary across all bins
        tot_nv2 = sum(d.nv2 for d in data)
        tot_nt3amp = sum(d.nt3amp for d in data)
        tot_nt3phi = sum(d.nt3phi for d in data)
        tot_nvisamp = sum(d.nvisamp for d in data)
        tot_nvisphi = sum(d.nvisphi for d in data)
        tot_nflux = sum(d.nflux for d in data)
        println("Totals: nv2=$tot_nv2  nt3amp=$tot_nt3amp  nt3phi=$tot_nt3phi" *
                "  nvisamp=$tot_nvisamp  nvisphi=$tot_nvisphi  nflux=$tot_nflux")
    end
end

# ---------------------------------------------------------------------------
# set_data_filter / filter_data — post-load filtering (unchanged logic)
# ---------------------------------------------------------------------------

"""
    set_data_filter(data; kwargs...) -> [uv_bad, vis_bad, v2_bad, t3_bad, flux_bad]

Compute lists of indices to discard from a loaded `OIdata` bin without modifying it.
Pass the result to `filter_data` to obtain a filtered copy.

# Keyword arguments
- `wav_range` — wavelength window(s) in metres, e.g. `[1.6e-6, 1.8e-6]` or a vector of
  windows `[[1.6e-6,1.8e-6],[2.0e-6,2.4e-6]]`. Default: keep all.
- `mjd_range` — MJD window(s), same format. Default: keep all.
- `baseline_range` — `[min, max]` baseline in cycles/m. Default: keep all.
- `filter_bad_data` — apply quality cuts (flags, NaN, SNR, amplitude range). Default: `false`.
- `filter_vis`, `filter_v2`, `filter_t3amp`, `filter_t3phi`, `filter_flux` — enable cuts per observable type.
- `cutoff_minv2`, `cutoff_maxv2` — V² range cut. Default: `(-1, 2.0)`.
- `cutoff_mint3amp`, `cutoff_maxt3amp` — T3 amplitude range cut. Default: `(-1.0, 1.5)`.
- `filter_v2_snr_threshold` — minimum |V²/σ| to keep. Default: `0.01`.
- `force_full_vis` — require both amplitude and phase to be valid (default: either).
- `force_full_t3` — require both T3amp and T3phi to be valid (default: either).
- `special_filter_diffvis` — enable differential visibility filtering mode.
- `uv_bad` — pre-supplied list of UV indices to remove.
- `filter_visphi`, `filter_visamp` — enable visibility phase/amplitude filtering.

Returns `[uv_bad, vis_bad, v2_bad, t3_bad, flux_bad]` — five `Vector{Int64}` of indices to discard.
"""
function set_data_filter(data::OIdata{T};
        wav_range::Union{Vector{Float64}, Vector{Vector{Float64}}} = [-1.0, 1e99],
        mjd_range::Union{Vector{Float64}, Vector{Vector{Float64}}} = [-1.0, 1e99],
        baseline_range::Vector{Float64} = [0.0, 1e99],
        filter_bad_data = false, filter_vis = true, filter_v2 = true,
        filter_t3amp = true, filter_t3phi = true, filter_flux = true,
        cutoff_minv2 = -1, cutoff_maxv2 = 2.0, cutoff_mint3amp = -1.0, cutoff_maxt3amp = 1.5,
        special_filter_diffvis = false, force_full_vis = false, force_full_t3 = false,
        filter_v2_snr_threshold = 0.01, uv_bad = Int64[],
        filter_visphi = false, filter_visamp = false) where T

    !isempty(wav_range) && typeof(wav_range) == Vector{Float64} && (wav_range = [wav_range])
    !isempty(mjd_range) && typeof(mjd_range) == Vector{Float64} && (mjd_range = [mjd_range])

    use_vis  = ((data.nvisphi > 0) && filter_visphi) || ((data.nvisamp > 0) && filter_visamp)
    use_v2   = (data.nv2 > 0) && filter_v2
    use_t3   = ((data.nt3amp > 0) && filter_t3amp) || ((data.nt3phi > 0) && filter_t3phi)
    use_flux = (data.nflux > 0) && filter_flux

    vis_bad = Int64[]; v2_bad = Int64[]; t3_bad = Int64[]; flux_bad = Int64[]

    if filter_bad_data
        if use_vis
            va_ok = (.!isnan.(data.visamp)) .& (.!isnan.(data.visamp_err)) .& (data.visamp_err .> 0)
            vp_ok = (.!isnan.(data.visphi)) .& (.!isnan.(data.visphi_err)) .& (data.visphi_err .> 0)
            vis_good = force_full_vis ?
                findall(.!data.vis_flag .& (va_ok .& vp_ok)) :
                findall(.!data.vis_flag .& (va_ok .| vp_ok))
            special_filter_diffvis && (vis_good = findall(data.vis_flag .!= 2))
            vis_bad = setdiff(1:length(data.vis_flag), vis_good)
        end
        if use_v2
            v2_good = findall(
                (.!data.v2_flag) .& (data.v2_err .> 0) .& (data.v2_err .< 1) .&
                (data.v2 .> cutoff_minv2) .& (data.v2 .< cutoff_maxv2) .&
                .!isnan.(data.v2) .& .!isnan.(data.v2_err) .&
                (abs.(data.v2 ./ data.v2_err) .> filter_v2_snr_threshold))
            v2_bad = setdiff(1:length(data.v2_flag), v2_good)
        end
        if use_t3
            ta_ok = (.!isnan.(data.t3amp)) .& (.!isnan.(data.t3amp_err)) .& (data.t3amp_err .> 0)
            tp_ok = (.!isnan.(data.t3phi)) .& (.!isnan.(data.t3phi_err)) .& (data.t3phi_err .> 0)
            t3_good = force_full_t3 ?
                findall(.!data.t3_flag .& (ta_ok .& tp_ok) .& (data.t3amp .> cutoff_mint3amp) .& (data.t3amp .< cutoff_maxt3amp)) :
                findall(.!data.t3_flag .& (ta_ok .| tp_ok))
            t3_bad = setdiff(1:length(data.t3_flag), t3_good)
        end
        if use_flux
            flux_good = findall(
                (.!data.flux_flag) .&
                .!isnan.(data.flux) .& .!isnan.(data.flux_err) .&
                (data.flux_err .> 0))
            flux_bad = setdiff(1:length(data.flux_flag), flux_good)
        end
    end

    baseline_range != [0.0, 1e99] && (uv_bad = union(uv_bad,
        findall(data.uv_baseline .< baseline_range[1]),
        findall(data.uv_baseline .> baseline_range[2])))
    uv_good = setdiff(1:data.nuv, uv_bad)
    !isempty(mjd_range) && (uv_good = intersect(uv_good,
        vcat([findall(mjd_range[i][1] .<= data.uv_mjd .<= mjd_range[i][2]) for i in eachindex(mjd_range)]...)))
    !isempty(wav_range) && (uv_good = intersect(uv_good,
        vcat([findall(wav_range[i][1] .<= data.uv_lam .<= wav_range[i][2]) for i in eachindex(wav_range)]...)))
    uv_bad = setdiff(1:data.nuv, uv_good)

    (data.nvisamp > 0 || data.nvisphi > 0) &&
        (vis_bad = union(vis_bad, findall([data.indx_vis[i] ∉ uv_good for i in eachindex(data.indx_vis)])))
    data.nv2 > 0 &&
        (v2_bad = union(v2_bad, findall([data.indx_v2[i] ∉ uv_good for i in eachindex(data.indx_v2)])))
    (data.nt3amp > 0 || data.nt3phi > 0) &&
        (t3_bad = union(t3_bad, findall([
            data.indx_t3_1[i] ∉ uv_good || data.indx_t3_2[i] ∉ uv_good || data.indx_t3_3[i] ∉ uv_good
            for i in eachindex(data.indx_t3_1)])))

    # Flux is not tied to the UV grid, but filter by wavelength and MJD ranges
    if data.nflux > 0
        flux_good_range = collect(1:data.nflux)
        if !isempty(wav_range) && wav_range != [[-1.0, 1e99]]
            flux_good_range = intersect(flux_good_range,
                vcat([findall(wav_range[i][1] .<= data.flux_lam .<= wav_range[i][2]) for i in eachindex(wav_range)]...))
        end
        if !isempty(mjd_range) && mjd_range != [[-1.0, 1e99]]
            flux_good_range = intersect(flux_good_range,
                vcat([findall(mjd_range[i][1] .<= data.flux_mjd .<= mjd_range[i][2]) for i in eachindex(mjd_range)]...))
        end
        flux_bad = union(flux_bad, setdiff(1:data.nflux, flux_good_range))
    end

    return [uv_bad, vis_bad, v2_bad, t3_bad, flux_bad]
end

"""
    filter_data(data, indexes_to_discard) -> OIdata

Return a deep copy of `data` with the specified points removed. `indexes_to_discard`
must be the five-element vector `[uv_bad, vis_bad, v2_bad, t3_bad, flux_bad]` returned by
`set_data_filter`. UV points that become unreferenced after removing observables are
pruned automatically and all index arrays are remapped.

# Example
```julia
idx = set_data_filter(data[1,1]; filter_bad_data=true, baseline_range=[5e6, 300e6])
clean = filter_data(data[1,1], idx)
```
"""
function filter_data(data_in::OIdata{T}, indexes_to_discard = Int64[]) where T
    data = deepcopy(data_in)
    good_uv_vis = Int64[]; good_uv_v2 = Int64[]
    good_uv_t3_1 = Int64[]; good_uv_t3_2 = Int64[]; good_uv_t3_3 = Int64[]

    if data.nvisamp > 0 || data.nvisphi > 0
        vis_good          = setdiff(1:length(data.indx_vis), indexes_to_discard[2])
        good_uv_vis       = data.indx_vis[vis_good]
        data.visamp       = data.visamp[vis_good];    data.visamp_err   = data.visamp_err[vis_good]
        data.visphi       = data.visphi[vis_good];    data.visphi_err   = data.visphi_err[vis_good]
        data.vis_baseline = data.vis_baseline[vis_good]
        data.vis_mjd      = data.vis_mjd[vis_good];   data.vis_lam      = data.vis_lam[vis_good]
        data.vis_dlam     = data.vis_dlam[vis_good];  data.vis_flag     = data.vis_flag[vis_good]
        data.vis_sta_index= data.vis_sta_index[:, vis_good]
        !isempty(data.visamp_corr_idx) && (data.visamp_corr_idx = data.visamp_corr_idx[vis_good])
        !isempty(data.visphi_corr_idx) && (data.visphi_corr_idx = data.visphi_corr_idx[vis_good])
        data.nvisamp      = count(i -> !isnan(data.visamp[i]) && !isnan(data.visamp_err[i]) && data.visamp_err[i] > 0, eachindex(data.visamp))
        data.nvisphi      = count(i -> !isnan(data.visphi[i]) && !isnan(data.visphi_err[i]) && data.visphi_err[i] > 0, eachindex(data.visphi))
    end
    if data.nv2 > 0
        v2_good          = setdiff(1:length(data.indx_v2), indexes_to_discard[3])
        good_uv_v2       = data.indx_v2[v2_good]
        data.v2          = data.v2[v2_good];          data.v2_err      = data.v2_err[v2_good]
        data.v2_baseline = data.v2_baseline[v2_good]; data.v2_mjd      = data.v2_mjd[v2_good]
        data.v2_lam      = data.v2_lam[v2_good];      data.v2_dlam     = data.v2_dlam[v2_good]
        data.v2_flag     = data.v2_flag[v2_good];     data.v2_sta_index= data.v2_sta_index[:, v2_good]
        !isempty(data.v2_corr_idx) && (data.v2_corr_idx = data.v2_corr_idx[v2_good])
        data.nv2         = length(data.v2)
    end
    if data.nt3amp > 0 || data.nt3phi > 0
        t3_good             = setdiff(1:length(data.indx_t3_1), indexes_to_discard[4])
        good_uv_t3_1 = data.indx_t3_1[t3_good]; good_uv_t3_2 = data.indx_t3_2[t3_good]
        good_uv_t3_3 = data.indx_t3_3[t3_good]
        data.t3amp          = data.t3amp[t3_good];    data.t3amp_err      = data.t3amp_err[t3_good]
        data.t3phi          = data.t3phi[t3_good];    data.t3phi_err      = data.t3phi_err[t3_good]
        data.t3_baseline    = data.t3_baseline[t3_good]; data.t3_maxbaseline = data.t3_maxbaseline[t3_good]
        data.t3_mjd         = data.t3_mjd[t3_good];  data.t3_lam         = data.t3_lam[t3_good]
        data.t3_dlam        = data.t3_dlam[t3_good];  data.t3_flag        = data.t3_flag[t3_good]
        data.t3_sta_index   = data.t3_sta_index[:, t3_good]
        !isempty(data.t3amp_corr_idx) && (data.t3amp_corr_idx = data.t3amp_corr_idx[t3_good])
        !isempty(data.t3phi_corr_idx) && (data.t3phi_corr_idx = data.t3phi_corr_idx[t3_good])
        data.nt3amp         = count(i -> !isnan(data.t3amp[i]) && !isnan(data.t3amp_err[i]) && data.t3amp_err[i] > 0, eachindex(data.t3amp))
        data.nt3phi         = count(i -> !isnan(data.t3phi[i]) && !isnan(data.t3phi_err[i]) && data.t3phi_err[i] > 0, eachindex(data.t3phi))
    end
    if data.nflux > 0 && length(indexes_to_discard) >= 5
        flux_good            = setdiff(1:data.nflux, indexes_to_discard[5])
        data.flux            = data.flux[flux_good];        data.flux_err      = data.flux_err[flux_good]
        data.flux_mjd        = data.flux_mjd[flux_good];    data.flux_lam      = data.flux_lam[flux_good]
        data.flux_dlam       = data.flux_dlam[flux_good];   data.flux_flag     = data.flux_flag[flux_good]
        data.flux_sta_index  = data.flux_sta_index[flux_good]
        !isempty(data.flux_corr_idx) && (data.flux_corr_idx = data.flux_corr_idx[flux_good])
        data.nflux           = length(data.flux)
    end
    uv_sel = falses(data.nuv)
    isempty(good_uv_vis)  || (uv_sel[good_uv_vis]  .= true)
    isempty(good_uv_v2)   || (uv_sel[good_uv_v2]   .= true)
    isempty(good_uv_t3_1) || (uv_sel[good_uv_t3_1] .= true)
    isempty(good_uv_t3_2) || (uv_sel[good_uv_t3_2] .= true)
    isempty(good_uv_t3_3) || (uv_sel[good_uv_t3_3] .= true)
    uv_sel[indexes_to_discard[1]] .= false
    iconv = cumsum(uv_sel); sel = findall(uv_sel)
    data.uv          = data.uv[:, sel];      data.uv_lam      = data.uv_lam[sel]
    data.uv_dlam     = data.uv_dlam[sel];    data.uv_mjd      = data.uv_mjd[sel]
    data.uv_baseline = data.uv_baseline[sel]; data.nuv        = size(data.uv, 2)
    isempty(good_uv_vis)  || (data.indx_vis   = iconv[good_uv_vis])
    isempty(good_uv_v2)   || (data.indx_v2    = iconv[good_uv_v2])
    isempty(good_uv_t3_1) || (data.indx_t3_1  = iconv[good_uv_t3_1])
    isempty(good_uv_t3_2) || (data.indx_t3_2  = iconv[good_uv_t3_2])
    isempty(good_uv_t3_3) || (data.indx_t3_3  = iconv[good_uv_t3_3])
    return data
end

"""
    remove_redundant_uv!(data; uvtol=2e2)

Merge redundant UV points in-place using a KD-tree with tolerance `uvtol` (cycles/rad).
Remaps all observable indices accordingly.
"""
function remove_redundant_uv!(data::OIdata{T}; uvtol = 2e2) where T
    iconv, tokeep    = rm_redundance_kdtree(data.uv, uvtol)
    data.uv          = data.uv[:, tokeep];  data.uv_lam  = data.uv_lam[tokeep]
    data.uv_dlam     = data.uv_dlam[tokeep]; data.uv_mjd = data.uv_mjd[tokeep]
    data.uv_baseline = data.uv_baseline[tokeep]; data.nuv = length(tokeep)
    data.nvisphi > 0 || data.nvisamp > 0 && (data.indx_vis = iconv[data.indx_vis])
    data.nv2 > 0    && (data.indx_v2  = iconv[data.indx_v2])
    if data.nt3amp > 0 || data.nt3phi > 0
        data.indx_t3_1 = iconv[data.indx_t3_1]
        data.indx_t3_2 = iconv[data.indx_t3_2]
        data.indx_t3_3 = iconv[data.indx_t3_3]
    end
    return data
end

# ===========================================================================
# LOADING + PARSING HELPERS
# ===========================================================================

# ---------------------------------------------------------------------------
# load_oifits: open an OIFITS file, falling back to hack_revn=1 for v1 files
# ---------------------------------------------------------------------------
function load_oifits(filename)
    try
        read(OIDataSet, filename)
    catch
        read(OIDataSet, filename, hack_revn=1)
    end
end

# ---------------------------------------------------------------------------
# _read_vis_header_keywords: read AMPTYP / PHITYP from OI_VIS FITS headers.
# OIFITS.jl may not populate these fields (e.g. when OI_REVN = 1), so we
# fall back to reading the raw FITS headers with FITSIO.
# Returns (amptyp, phityp) strings — first non-empty value found wins.
# ---------------------------------------------------------------------------
function _read_vis_header_keywords(filename::String)
    amptyp = ""; phityp = ""
    try
        f = FITS(filename)
        for i in 1:length(f)
            hdr = read_header(f[i])
            get(hdr, "EXTNAME", "") == "OI_VIS" || continue
            at = get(hdr, "AMPTYP", "")
            pt = get(hdr, "PHITYP", "")
            isempty(amptyp) && !isempty(at) && (amptyp = lowercase(strip(at)))
            isempty(phityp) && !isempty(pt) && (phityp = lowercase(strip(pt)))
            !isempty(amptyp) && !isempty(phityp) && break
        end
        close(f)
    catch
    end
    return amptyp, phityp
end

# ---------------------------------------------------------------------------
# collect_station_info: parse OI_ARRAY tables and build station index mapping
#
# Returns a NamedTuple:
#   new_station_name, new_telescope_name, new_station_index  — merged station lists
#   conversion_index[itable, old_idx] → new_idx             — sparse mapping
#   station_index_offset                                     — 0 or 1 (OIFITSv1 compat)
# ---------------------------------------------------------------------------
function collect_station_info(arraytables, v2tables, t3tables, arraytableref; warn=true)
    array_ntables  = length(arraytables)
    station_name   = [arraytables[i].sta_name  for i in 1:array_ntables]
    telescope_name = [arraytables[i].tel_name  for i in 1:array_ntables]
    station_index  = [arraytables[i].sta_index for i in 1:array_ntables]

    station_index_offset = 0
    if minimum(vcat(station_index...)) == 0
        revn_max = maximum(db.revn for db in arraytables)
        if warn
            if revn_max == 2
                @warn("This file does not follow OIFITSv2 standard — station indexing should start at 1, not 0.")
            else
                @warn("OIFITSv1 detected. OITOOLS will internally reindex stations from 1.")
            end
        end
        station_index_offset = 1
    end

    # Detect station indices referenced in V2/T3 tables but absent from OI_ARRAY
    unknown_station_names   = String[]
    unknown_tel_names       = String[]
    unknown_station_indexes = Int64[]

    for (itable, db) in enumerate(v2tables)
        iarray = findfirst(==(db.arrname), arraytableref)
        if !isnothing(iarray)
            corresp = arraytables[iarray].sta_index
            for jj in unique(db.sta_index)
                if jj ∉ corresp
                    @warn("V2 table $itable refers to station index $jj, not found in OI_ARRAY=$(db.arrname)")
                    push!(unknown_station_names, "UNKN$jj")
                    push!(unknown_tel_names,     "UNKN$jj")
                    push!(unknown_station_indexes, jj)
                end
            end
        else
            @warn("V2 table $itable is missing its corresponding OI_ARRAY $(db.arrname)")
        end
    end

    for (itable, db) in enumerate(t3tables)
        iarray = findfirst(==(db.arrname), arraytableref)
        if !isnothing(iarray)
            corresp = arraytables[iarray].sta_index
            for jj in unique(db.sta_index)
                if jj ∉ corresp
                    @warn("T3 table $itable refers to station index $jj, not found in OI_ARRAY=$(db.arrname)")
                    push!(unknown_station_names, "UNKN$jj")
                    push!(unknown_tel_names,     "UNKN$jj")
                    push!(unknown_station_indexes, jj)
                end
            end
        else
            @warn("T3 table $itable is missing its corresponding OI_ARRAY $(db.arrname)")
        end
    end

    if !isempty(unknown_station_names)
        nr = indexin(unique(unknown_station_indexes), unknown_station_indexes)
        unknown_station_indexes = unknown_station_indexes[nr]
        unknown_station_names   = unknown_station_names[nr]
        unknown_tel_names       = unknown_tel_names[nr]
        sp = sortperm(unknown_station_indexes)
        unknown_station_indexes = unknown_station_indexes[sp]
        unknown_station_names   = unknown_station_names[sp]
        unknown_tel_names       = unknown_tel_names[sp]
        @warn("Unknown stations: $unknown_station_names with indexes = $unknown_station_indexes")
    end

    station_names_all   = vec(vcat(station_name...,   unknown_station_names))
    telescope_names_all = vec(vcat(telescope_name..., unknown_tel_names))
    station_indexes_all = vec(vcat(station_index...,  unknown_station_indexes))

    list_stations = unique(station_names_all)
    nstations     = length(list_stations)
    new_station_name   = Vector{String}(undef, nstations)
    new_telescope_name = Vector{String}(undef, nstations)
    new_station_index  = zeros(Int64, nstations)

    for istation in 1:nstations
        name = list_stations[istation]
        loc  = findall(station_names_all .== name)
        tel  = unique(telescope_names_all[loc])
        length(tel) > 1 && @warn("Multiple telescopes at station $name — plots may be misleading")
        new_station_name[istation]   = name
        new_telescope_name[istation] = tel[1]
        new_station_index[istation]  = istation
    end

    maxidx = maximum(station_indexes_all) + station_index_offset
    conversion_index = spzeros(Int64, array_ntables, maxidx)
    for itable in 1:array_ntables
        for istation in 1:length(station_name[itable])
            name    = station_name[itable][istation]
            oldindx = station_index[itable][istation]
            newindx = findfirst(==(name), new_station_name)
            conversion_index[itable, station_index_offset + oldindx] = newindx
        end
        if !isempty(unknown_station_names)
            nmax = sum(conversion_index[itable, :] .!= 0)
            conversion_index[itable, unknown_station_indexes] .= nmax+1:maxidx
        end
    end

    return (; new_station_name, new_telescope_name, new_station_index,
              conversion_index, station_index_offset)
end

# ===========================================================================
# TABLE READERS — each returns a NamedTuple of flat (all-obs merged) arrays
# ===========================================================================

# ---------------------------------------------------------------------------
# read_flux_tables
#
# OIFITS2 §7.1 defines two cases distinguished by the CALSTAT keyword:
#
#   CALSTAT = "C"  (Calibrated): source spectrum, no ARRNAME or STA_INDEX.
#                  flux_sta_index will be zeros (no telescope association).
#
#   CALSTAT = "U"  (Uncalibrated): per-telescope flux.  ARRNAME and STA_INDEX
#                  are mandatory.  We convert station indices the same way as
#                  vis/v2/t3.  flux_sta_index stores the converted station for
#                  every (wavelength × observation) point.
#
# OI_FLUX also has TARGET_ID, so we apply the same target filter as other tables.
# ---------------------------------------------------------------------------
function read_flux_tables(fluxtables, targetid_filter,
                          wavtables, wavtableref,
                          arraytableref, conversion_index, station_index_offset;
                          correlt::Dict=Dict(), T::Type{<:AbstractFloat}=Float64)
    flux_all = T[]; flux_err_all = T[]; flux_mjd_all = T[]
    flux_lam_all = T[]; flux_dlam_all = T[]; flux_flag_all = Bool[]
    flux_sta_index_all = Int64[]
    flux_calibrated = false   # will be set true if any table has CALSTAT="C"
    flux_corr_idx_parts = Vector{Int64}[]
    seen_flux_names = String[]; flux_corr_offsets = Dict{String,Int}()

    for db in fluxtables
        tid_ok = findall(sum([db.target_id .== id for id in targetid_filter], dims=1)[1] .> 0)
        isempty(tid_ok) && continue

        wav   = _find_wav(db, wavtables, wavtableref)
        nwave = length(wav.eff_wave)
        nobs  = length(tid_ok)

        fdata = T.(db.fluxdata[:, tid_ok])
        ferr  = T.(db.fluxerr[:,  tid_ok])
        fflag = Bool.(db.flag[:, tid_ok])
        fmjd  = repeat(T.(db.mjd[tid_ok])', nwave)
        flam  = repeat(T.(wav.eff_wave), 1, nobs)
        fdlam = repeat(T.(wav.eff_band), 1, nobs)

        # Station index: zero for calibrated spectra, converted index for uncalibrated
        calstat = uppercase(strip(db.calstat))
        if calstat == "C"
            flux_calibrated = true
        end
        if calstat == "U"
            iarray = findfirst(==(db.arrname), arraytableref)
            if !isnothing(iarray)
                raw_si = Int.(Matrix(conversion_index[iarray,
                    station_index_offset .+ db.sta_index[tid_ok]']))
            else
                raw_si = 1000 .+ station_index_offset .+ db.sta_index[tid_ok]'
            end
            # raw_si is 1×nobs; repeat over wavelengths → nwave×nobs, then flatten
            fsi = vec(repeat(raw_si, nwave, 1))
        else  # "C": calibrated source spectrum, no telescope association
            fsi = zeros(Int64, nwave * nobs)
        end

        push!(flux_corr_idx_parts, _process_corrindx(
            _safe_field(db, :corrname), _safe_field(db, :corrindx_fluxdata),
            tid_ok, nwave, correlt, flux_corr_offsets, seen_flux_names;
            mjd_obs=db.mjd[tid_ok]))

        append!(flux_all,          vec(fdata));  append!(flux_err_all,      vec(ferr))
        append!(flux_mjd_all,      vec(fmjd));   append!(flux_lam_all,      vec(flam))
        append!(flux_dlam_all,     vec(fdlam));  append!(flux_flag_all,     vec(fflag))
        append!(flux_sta_index_all, fsi)
    end
    flux_corr_idx = isempty(seen_flux_names) ? Int64[] : vcat(flux_corr_idx_parts...)
    flux_corr     = _assemble_corr(seen_flux_names, correlt, T)
    return (; flux=flux_all, flux_err=flux_err_all, flux_mjd=flux_mjd_all,
              flux_lam=flux_lam_all, flux_dlam=flux_dlam_all, flux_flag=flux_flag_all,
              flux_sta_index=flux_sta_index_all, flux_calibrated, flux_corr, flux_corr_idx)
end

# ---------------------------------------------------------------------------
# read_vis_tables
# Returns vis_uv as N×2 matrix (column 1 = u/λ, column 2 = v/λ)
# ---------------------------------------------------------------------------
function read_vis_tables(vistables, targetid_filter, wavtables, wavtableref,
                         arraytableref, conversion_index, station_index_offset;
                         correlt::Dict=Dict(), T::Type{<:AbstractFloat}=Float64)
    visamp_all = T[]; visamp_err_all = T[]
    visphi_all = T[]; visphi_err_all = T[]
    vis_mjd_all = T[]; vis_lam_all = T[]; vis_dlam_all = T[]
    vis_flag_all = Bool[]
    vis_uv_all   = Matrix{T}(undef, 0, 2)   # N×2
    vis_baseline_all  = T[]
    vis_sta_index_all = Matrix{Int64}(undef, 2, 0)  # 2×N
    visamp_corr_idx_parts = Vector{Int64}[]; visphi_corr_idx_parts = Vector{Int64}[]
    seen_va_names = String[]; va_corr_offsets = Dict{String,Int}()
    seen_vp_names = String[]; vp_corr_offsets = Dict{String,Int}()
    amptyp_all = ""; phityp_all = ""

    for (itable, db) in enumerate(vistables)
        # Read AMPTYP / PHITYP (OIFITSv2 §6).  OIFITS.jl may leave these
        # undefined for REVN=1 tables, so fall back to empty string.
        _at = isdefined(db, :amptyp) ? string(db.amptyp) : ""
        _pt = isdefined(db, :phityp) ? string(db.phityp) : ""
        if !isempty(_at) && isempty(amptyp_all); amptyp_all = _at; end
        if !isempty(_pt) && isempty(phityp_all); phityp_all = _pt; end
        tid_ok = findall(sum([db.target_id .== id for id in targetid_filter], dims=1)[1] .> 0)
        wav    = _find_wav(db, wavtables, wavtableref)
        nwave  = length(wav.eff_wave)
        vamp   = T.(db.visamp[:,    tid_ok])
        verr   = T.(db.visamperr[:, tid_ok])
        vphi   = T.(db.visphi[:,    tid_ok])
        vperr  = T.(db.visphierr[:, tid_ok])
        vmjd   = repeat(T.(db.mjd[tid_ok])', nwave)
        vflag  = Bool.(db.flag[:, tid_ok])
        uc     = T.(db.ucoord[tid_ok]);  vc = T.(db.vcoord[tid_ok])
        vlam   = repeat(T.(wav.eff_wave), 1, length(tid_ok))
        vdlam  = repeat(T.(wav.eff_band), 1, length(tid_ok))
        vu, vv = _uv_lambda(uc, vc, T.(wav.eff_wave))

        iarray = findfirst(==(db.arrname), arraytableref)
        si = if !isnothing(iarray)
            Matrix(conversion_index[iarray,
                station_index_offset .+ repeat(db.sta_index[:, tid_ok], outer=[nwave, 1])])
        else
            1000 .+ station_index_offset .+
                repeat(db.sta_index[:, tid_ok], outer=[nwave, 1])
        end

        cn = _safe_field(db, :corrname)
        _mjd_tid = db.mjd[tid_ok]
        push!(visamp_corr_idx_parts, _process_corrindx(
            cn, _safe_field(db, :corrindx_visamp), tid_ok, nwave, correlt, va_corr_offsets, seen_va_names;
            mjd_obs=_mjd_tid))
        push!(visphi_corr_idx_parts, _process_corrindx(
            cn, _safe_field(db, :corrindx_visphi), tid_ok, nwave, correlt, vp_corr_offsets, seen_vp_names;
            mjd_obs=_mjd_tid))

        append!(visamp_all, vec(vamp));   append!(visamp_err_all, vec(verr))
        append!(visphi_all, vec(vphi));   append!(visphi_err_all, vec(vperr))
        append!(vis_mjd_all,  vec(vmjd)); append!(vis_lam_all,  vec(vlam))
        append!(vis_dlam_all, vec(vdlam)); append!(vis_flag_all, vec(vflag))
        vis_uv_all = vcat(vis_uv_all, hcat(vu, vv))
        append!(vis_baseline_all, sqrt.(vu.^2 .+ vv.^2))
        vis_sta_index_all = hcat(vis_sta_index_all, reshape(si, 2, div(length(si), 2)))
    end
    visamp_corr_idx = isempty(seen_va_names) ? Int64[] : vcat(visamp_corr_idx_parts...)
    visphi_corr_idx = isempty(seen_vp_names) ? Int64[] : vcat(visphi_corr_idx_parts...)
    visamp_corr = _assemble_corr(seen_va_names, correlt, T)
    visphi_corr = _assemble_corr(seen_vp_names, correlt, T)
    return (; visamp=visamp_all, visamp_err=visamp_err_all,
              visphi=visphi_all, visphi_err=visphi_err_all,
              vis_mjd=vis_mjd_all, vis_lam=vis_lam_all, vis_dlam=vis_dlam_all,
              vis_flag=vis_flag_all, vis_uv=vis_uv_all,
              vis_baseline=vis_baseline_all, vis_sta_index=vis_sta_index_all,
              visamp_corr, visamp_corr_idx, visphi_corr, visphi_corr_idx,
              amptyp=amptyp_all, phityp=phityp_all)
end

# ---------------------------------------------------------------------------
# read_v2_tables
# ---------------------------------------------------------------------------
function read_v2_tables(v2tables, targetid_filter, wavtables, wavtableref,
                        arraytableref, conversion_index, station_index_offset;
                        correlt::Dict=Dict(), T::Type{<:AbstractFloat}=Float64)
    v2_all = T[]; v2_err_all = T[]
    v2_mjd_all = T[]; v2_lam_all = T[]; v2_dlam_all = T[]
    v2_flag_all = Bool[]
    v2_uv_all   = Matrix{T}(undef, 0, 2)   # N×2
    v2_baseline_all  = T[]
    v2_sta_index_all = Matrix{Int64}(undef, 2, 0)  # 2×N
    v2_corr_idx_parts = Vector{Int64}[]
    seen_v2_names = String[]; v2_corr_offsets = Dict{String,Int}()

    for (itable, db) in enumerate(v2tables)
        tid_ok = findall(sum([db.target_id .== id for id in targetid_filter], dims=1)[1] .> 0)
        wav    = _find_wav(db, wavtables, wavtableref)
        nwave  = length(wav.eff_wave)
        v2d    = T.(db.vis2data[:, tid_ok])
        v2e    = T.(db.vis2err[:,  tid_ok])
        v2mjd  = repeat(T.(db.mjd[tid_ok])', nwave)
        v2flag = Bool.(db.flag[:, tid_ok])
        uc     = T.(db.ucoord[tid_ok]);  vc = T.(db.vcoord[tid_ok])
        vlam   = repeat(T.(wav.eff_wave), 1, length(tid_ok))
        vdlam  = repeat(T.(wav.eff_band), 1, length(tid_ok))
        vu, vv = _uv_lambda(uc, vc, T.(wav.eff_wave))

        iarray = findfirst(==(db.arrname), arraytableref)
        si = if !isnothing(iarray)
            Matrix(conversion_index[iarray,
                station_index_offset .+ repeat(db.sta_index[:, tid_ok], outer=[nwave, 1])])
        else
            1000 .+ station_index_offset .+
                repeat(db.sta_index[:, tid_ok], outer=[nwave, 1])
        end

        push!(v2_corr_idx_parts, _process_corrindx(
            _safe_field(db, :corrname), _safe_field(db, :corrindx_vis2data),
            tid_ok, nwave, correlt, v2_corr_offsets, seen_v2_names;
            mjd_obs=db.mjd[tid_ok]))

        append!(v2_all, vec(v2d));    append!(v2_err_all, vec(v2e))
        append!(v2_mjd_all, vec(v2mjd))
        append!(v2_lam_all, vec(vlam)); append!(v2_dlam_all, vec(vdlam))
        append!(v2_flag_all, vec(v2flag))
        v2_uv_all = vcat(v2_uv_all, hcat(vu, vv))
        append!(v2_baseline_all, sqrt.(vu.^2 .+ vv.^2))
        v2_sta_index_all = hcat(v2_sta_index_all, reshape(si, 2, div(length(si), 2)))
    end
    v2_corr_idx = isempty(seen_v2_names) ? Int64[] : vcat(v2_corr_idx_parts...)
    v2_corr     = _assemble_corr(seen_v2_names, correlt, T)
    return (; v2=v2_all, v2_err=v2_err_all, v2_mjd=v2_mjd_all,
              v2_lam=v2_lam_all, v2_dlam=v2_dlam_all, v2_flag=v2_flag_all,
              v2_uv=v2_uv_all, v2_baseline=v2_baseline_all, v2_sta_index=v2_sta_index_all,
              v2_corr, v2_corr_idx)
end

# ---------------------------------------------------------------------------
# read_t3_tables
# UV legs stored as separate vectors (t3_u1/v1, t3_u2/v2, t3_u3/v3) since
# each T3 has 3 legs; they are NOT combined into a single uv matrix here.
# ---------------------------------------------------------------------------
function read_t3_tables(t3tables, targetid_filter, wavtables, wavtableref,
                        arraytableref, conversion_index, station_index_offset;
                        correlt::Dict=Dict(), T::Type{<:AbstractFloat}=Float64)
    t3amp_all = T[]; t3amp_err_all = T[]
    t3phi_all = T[]; t3phi_err_all = T[]
    t3_mjd_all = T[]; t3_lam_all = T[]; t3_dlam_all = T[]
    t3_flag_all = Bool[]
    t3_u1_all = T[]; t3_v1_all = T[]
    t3_u2_all = T[]; t3_v2_all = T[]
    t3_u3_all = T[]; t3_v3_all = T[]
    t3_baseline_all    = T[]
    t3_maxbaseline_all = T[]
    t3_sta_index_all   = Matrix{Int64}(undef, 3, 0)  # 3×N
    t3amp_corr_idx_parts = Vector{Int64}[]; t3phi_corr_idx_parts = Vector{Int64}[]
    seen_ta_names = String[]; ta_corr_offsets = Dict{String,Int}()
    seen_tp_names = String[]; tp_corr_offsets = Dict{String,Int}()

    for (itable, db) in enumerate(t3tables)
        tid_ok = findall(sum([db.target_id .== id for id in targetid_filter], dims=1)[1] .> 0)
        wav    = _find_wav(db, wavtables, wavtableref)
        nwave  = length(wav.eff_wave)
        t3a    = T.(db.t3amp[:,    tid_ok]); t3ae = T.(db.t3amperr[:, tid_ok])
        t3p    = T.(db.t3phi[:,    tid_ok]); t3pe = T.(db.t3phierr[:,  tid_ok])
        t3mjd  = repeat(T.(db.mjd[tid_ok])', nwave)
        t3flag = Bool.(db.flag[:, tid_ok])
        u1c = T.(db.u1coord[tid_ok]); v1c = T.(db.v1coord[tid_ok])
        u2c = T.(db.u2coord[tid_ok]); v2c = T.(db.v2coord[tid_ok])
        u3c = -(u1c .+ u2c);          v3c = -(v1c .+ v2c)
        tlam  = repeat(T.(wav.eff_wave), 1, length(tid_ok))
        tdlam = repeat(T.(wav.eff_band), 1, length(tid_ok))
        u1, v1 = _uv_lambda(u1c, v1c, T.(wav.eff_wave))
        u2, v2 = _uv_lambda(u2c, v2c, T.(wav.eff_wave))
        u3, v3 = _uv_lambda(u3c, v3c, T.(wav.eff_wave))
        b1 = sqrt.(u1.^2 .+ v1.^2); b2 = sqrt.(u2.^2 .+ v2.^2); b3 = sqrt.(u3.^2 .+ v3.^2)

        iarray = findfirst(==(db.arrname), arraytableref)
        si = if !isnothing(iarray)
            Matrix(conversion_index[iarray,
                station_index_offset .+ repeat(db.sta_index[:, tid_ok], outer=[nwave, 1])])
        else
            1000 .+ station_index_offset .+
                repeat(db.sta_index[:, tid_ok], outer=[nwave, 1])
        end

        cn = _safe_field(db, :corrname)
        _mjd_tid = db.mjd[tid_ok]
        push!(t3amp_corr_idx_parts, _process_corrindx(
            cn, _safe_field(db, :corrindx_t3amp), tid_ok, nwave, correlt, ta_corr_offsets, seen_ta_names;
            mjd_obs=_mjd_tid))
        push!(t3phi_corr_idx_parts, _process_corrindx(
            cn, _safe_field(db, :corrindx_t3phi), tid_ok, nwave, correlt, tp_corr_offsets, seen_tp_names;
            mjd_obs=_mjd_tid))

        append!(t3amp_all, vec(t3a));   append!(t3amp_err_all, vec(t3ae))
        append!(t3phi_all, vec(t3p));   append!(t3phi_err_all, vec(t3pe))
        append!(t3_mjd_all, vec(t3mjd)); append!(t3_lam_all, vec(tlam))
        append!(t3_dlam_all, vec(tdlam)); append!(t3_flag_all, vec(t3flag))
        append!(t3_u1_all, u1); append!(t3_v1_all, v1)
        append!(t3_u2_all, u2); append!(t3_v2_all, v2)
        append!(t3_u3_all, u3); append!(t3_v3_all, v3)
        append!(t3_baseline_all,    (b1 .* b2 .* b3) .^ (one(T)/3))
        append!(t3_maxbaseline_all, max.(b1, b2, b3))
        t3_sta_index_all = hcat(t3_sta_index_all, reshape(si, 3, div(length(si), 3)))
    end
    t3amp_corr_idx = isempty(seen_ta_names) ? Int64[] : vcat(t3amp_corr_idx_parts...)
    t3phi_corr_idx = isempty(seen_tp_names) ? Int64[] : vcat(t3phi_corr_idx_parts...)
    t3amp_corr = _assemble_corr(seen_ta_names, correlt, T)
    t3phi_corr = _assemble_corr(seen_tp_names, correlt, T)
    return (; t3amp=t3amp_all, t3amp_err=t3amp_err_all,
              t3phi=t3phi_all, t3phi_err=t3phi_err_all,
              t3_mjd=t3_mjd_all, t3_lam=t3_lam_all, t3_dlam=t3_dlam_all,
              t3_flag=t3_flag_all,
              t3_u1=t3_u1_all, t3_v1=t3_v1_all,
              t3_u2=t3_u2_all, t3_v2=t3_v2_all,
              t3_u3=t3_u3_all, t3_v3=t3_v3_all,
              t3_baseline=t3_baseline_all, t3_maxbaseline=t3_maxbaseline_all,
              t3_sta_index=t3_sta_index_all,
              t3amp_corr, t3amp_corr_idx, t3phi_corr, t3phi_corr_idx)
end

# ===========================================================================
# BIN-LEVEL DATA HANDLING
# ===========================================================================

# ---------------------------------------------------------------------------
# BinData{T}: mutable workspace for one (spectral × temporal) bin.
# Holds all per-bin arrays. Passed into filter_bad_observables! and
# remove_redundant_uv! which modify it in place, then into make_oidata.
# ---------------------------------------------------------------------------
mutable struct BinData{T<:AbstractFloat}
    # VIS
    visamp::Vector{T};         visamp_err::Vector{T}
    visphi::Vector{T};         visphi_err::Vector{T}
    vis_mjd::Vector{T};        vis_lam::Vector{T};    vis_dlam::Vector{T}
    vis_flag::Vector{Bool};    vis_baseline::Vector{T}
    vis_sta_index::Matrix{Int64}
    indx_vis::Vector{Int64};   nvisamp::Int64;         nvisphi::Int64
    # V2
    v2::Vector{T};             v2_err::Vector{T}
    v2_mjd::Vector{T};         v2_lam::Vector{T};     v2_dlam::Vector{T}
    v2_flag::Vector{Bool};     v2_baseline::Vector{T}
    v2_sta_index::Matrix{Int64}
    indx_v2::Vector{Int64};    nv2::Int64
    # T3
    t3amp::Vector{T};          t3amp_err::Vector{T}
    t3phi::Vector{T};          t3phi_err::Vector{T}
    t3_mjd::Vector{T};         t3_lam::Vector{T};     t3_dlam::Vector{T}
    t3_flag::Vector{Bool};     t3_baseline::Vector{T}; t3_maxbaseline::Vector{T}
    t3_sta_index::Matrix{Int64}
    indx_t3_1::Vector{Int64};  indx_t3_2::Vector{Int64}; indx_t3_3::Vector{Int64}
    nt3amp::Int64;             nt3phi::Int64
    # Flux
    flux::Vector{T};           flux_err::Vector{T}
    flux_mjd::Vector{T};       flux_lam::Vector{T};   flux_dlam::Vector{T}
    flux_flag::Vector{Bool};   flux_sta_index::Vector{Int64};  nflux::Int64
    flux_calibrated::Bool      # true if CALSTAT="C"
    # UV plane (2×nuv: row 1 = u/λ, row 2 = v/λ)
    uv::Matrix{T};             uv_lam::Vector{T};     uv_dlam::Vector{T}
    uv_mjd::Vector{T};         uv_baseline::Vector{T}
    nuv::Int64;                mean_mjd::Float64
    # Correlation matrices (empty sparse if no OI_CORR in file)
    v2_corr::SparseMatrixCSC{T,Int64};       v2_corr_idx::Vector{Int64}
    t3amp_corr::SparseMatrixCSC{T,Int64};    t3amp_corr_idx::Vector{Int64}
    t3phi_corr::SparseMatrixCSC{T,Int64};    t3phi_corr_idx::Vector{Int64}
    visamp_corr::SparseMatrixCSC{T,Int64};   visamp_corr_idx::Vector{Int64}
    visphi_corr::SparseMatrixCSC{T,Int64};   visphi_corr_idx::Vector{Int64}
    flux_corr::SparseMatrixCSC{T,Int64};     flux_corr_idx::Vector{Int64}
    # OI_VIS metadata
    amptyp::String;                          phityp::String
end

# ---------------------------------------------------------------------------
# slice_to_bin: create a BinData from the flat raw arrays + boolean bin masks.
# Also assembles the UV plane (vis columns first, then v2, then 3×t3 legs).
# ---------------------------------------------------------------------------
function slice_to_bin(raw_vis, raw_v2, raw_t3, raw_flux,
                      bin_vis, bin_v2, bin_t3, bin_flux;
                      use_vis::Bool, use_v2::Bool, use_t3::Bool, use_flux::Bool,
                      T::Type{<:AbstractFloat}=Float64)

    # -- VIS --
    if use_vis
        visamp     = raw_vis.visamp[bin_vis];     visamp_err = raw_vis.visamp_err[bin_vis]
        visphi     = raw_vis.visphi[bin_vis];     visphi_err = raw_vis.visphi_err[bin_vis]
        vis_mjd    = raw_vis.vis_mjd[bin_vis];    vis_lam    = raw_vis.vis_lam[bin_vis]
        vis_dlam   = raw_vis.vis_dlam[bin_vis];   vis_flag   = raw_vis.vis_flag[bin_vis]
        vis_baseline = raw_vis.vis_baseline[bin_vis]
        vis_uv     = hcat(raw_vis.vis_uv[bin_vis, 1], raw_vis.vis_uv[bin_vis, 2])'  # 2×nvis
        vis_sta_index = raw_vis.vis_sta_index[:, bin_vis]
        nvisamp    = count(i -> !isnan(visamp[i]) && !isnan(visamp_err[i]) && visamp_err[i] > 0, eachindex(visamp))
        nvisphi    = count(i -> !isnan(visphi[i]) && !isnan(visphi_err[i]) && visphi_err[i] > 0, eachindex(visphi))
        indx_vis   = collect(1:length(visamp))
        visamp_corr_idx = isempty(raw_vis.visamp_corr_idx) ? Int64[] : raw_vis.visamp_corr_idx[bin_vis]
        visphi_corr_idx = isempty(raw_vis.visphi_corr_idx) ? Int64[] : raw_vis.visphi_corr_idx[bin_vis]
    else
        visamp = T[]; visamp_err = T[]; visphi = T[]; visphi_err = T[]
        vis_mjd = T[]; vis_lam = T[]; vis_dlam = T[]
        vis_flag = Bool[]; vis_baseline = T[]; vis_uv = Matrix{T}(undef, 2, 0)
        vis_sta_index = Matrix{Int64}(undef, 2, 0)
        nvisamp = 0; nvisphi = 0; indx_vis = Int64[]
        visamp_corr_idx = Int64[]; visphi_corr_idx = Int64[]
    end
    visamp_corr = raw_vis.visamp_corr;  visphi_corr = raw_vis.visphi_corr

    # -- V2 --
    if use_v2
        v2      = raw_v2.v2[bin_v2];       v2_err   = raw_v2.v2_err[bin_v2]
        v2_mjd  = raw_v2.v2_mjd[bin_v2];   v2_lam   = raw_v2.v2_lam[bin_v2]
        v2_dlam = raw_v2.v2_dlam[bin_v2];  v2_flag  = raw_v2.v2_flag[bin_v2]
        v2_baseline = raw_v2.v2_baseline[bin_v2]
        v2_uv   = hcat(raw_v2.v2_uv[bin_v2, 1], raw_v2.v2_uv[bin_v2, 2])'  # 2×nv2
        v2_sta_index = raw_v2.v2_sta_index[:, bin_v2]
        nv2     = length(v2)
        indx_v2 = collect(length(visamp) .+ (1:nv2))
        v2_corr_idx = isempty(raw_v2.v2_corr_idx) ? Int64[] : raw_v2.v2_corr_idx[bin_v2]
    else
        v2 = T[]; v2_err = T[]; v2_mjd = T[]; v2_lam = T[]; v2_dlam = T[]
        v2_flag = Bool[]; v2_baseline = T[]; v2_uv = Matrix{T}(undef, 2, 0)
        v2_sta_index = Matrix{Int64}(undef, 2, 0)
        nv2 = 0; indx_v2 = Int64[]
        v2_corr_idx = Int64[]
    end
    v2_corr = raw_v2.v2_corr

    # -- T3 --
    if use_t3
        t3amp    = raw_t3.t3amp[bin_t3];    t3amp_err  = raw_t3.t3amp_err[bin_t3]
        t3phi    = raw_t3.t3phi[bin_t3];    t3phi_err  = raw_t3.t3phi_err[bin_t3]
        t3_mjd   = raw_t3.t3_mjd[bin_t3];   t3_lam   = raw_t3.t3_lam[bin_t3]
        t3_dlam  = raw_t3.t3_dlam[bin_t3];  t3_flag  = raw_t3.t3_flag[bin_t3]
        t3_baseline    = raw_t3.t3_baseline[bin_t3]
        t3_maxbaseline = raw_t3.t3_maxbaseline[bin_t3]
        t3_sta_index   = raw_t3.t3_sta_index[:, bin_t3]
        nt3    = length(t3amp)
        nt3amp = count(i -> !isnan(t3amp[i]) && !isnan(t3amp_err[i]) && t3amp_err[i] > 0, eachindex(t3amp))
        nt3phi = count(i -> !isnan(t3phi[i]) && !isnan(t3phi_err[i]) && t3phi_err[i] > 0, eachindex(t3phi))
        off    = length(visamp) + nv2
        indx_t3_1 = collect(off .+ (1:nt3))
        indx_t3_2 = collect(off .+ (nt3+1:2nt3))
        indx_t3_3 = collect(off .+ (2nt3+1:3nt3))
        t3amp_corr_idx = isempty(raw_t3.t3amp_corr_idx) ? Int64[] : raw_t3.t3amp_corr_idx[bin_t3]
        t3phi_corr_idx = isempty(raw_t3.t3phi_corr_idx) ? Int64[] : raw_t3.t3phi_corr_idx[bin_t3]
    else
        t3amp = T[]; t3amp_err = T[]; t3phi = T[]; t3phi_err = T[]
        t3_mjd = T[]; t3_lam = T[]; t3_dlam = T[]; t3_flag = Bool[]
        t3_baseline = T[]; t3_maxbaseline = T[]; t3_sta_index = Matrix{Int64}(undef, 3, 0)
        nt3amp = 0; nt3phi = 0; indx_t3_1 = Int64[]; indx_t3_2 = Int64[]; indx_t3_3 = Int64[]
        t3amp_corr_idx = Int64[]; t3phi_corr_idx = Int64[]
    end
    t3amp_corr = raw_t3.t3amp_corr;  t3phi_corr = raw_t3.t3phi_corr

    # -- Flux --
    if use_flux
        flux      = raw_flux.flux[bin_flux];      flux_err  = raw_flux.flux_err[bin_flux]
        flux_mjd  = raw_flux.flux_mjd[bin_flux];  flux_lam  = raw_flux.flux_lam[bin_flux]
        flux_dlam = raw_flux.flux_dlam[bin_flux]; flux_flag = raw_flux.flux_flag[bin_flux]
        flux_sta_index = raw_flux.flux_sta_index[bin_flux]
        nflux     = length(flux)
        flux_calibrated = raw_flux.flux_calibrated
        flux_corr_idx = isempty(raw_flux.flux_corr_idx) ? Int64[] : raw_flux.flux_corr_idx[bin_flux]
    else
        flux = T[]; flux_err = T[]; flux_mjd = T[]; flux_lam = T[]
        flux_dlam = T[]; flux_flag = Bool[]; flux_sta_index = Int64[]; nflux = 0
        flux_calibrated = false
        flux_corr_idx = Int64[]
    end
    flux_corr = raw_flux.flux_corr

    # -- Assemble UV plane (2×nuv) --
    t3_uv = use_t3 ? hcat(
        hcat(raw_t3.t3_u1[bin_t3], raw_t3.t3_v1[bin_t3])',
        hcat(raw_t3.t3_u2[bin_t3], raw_t3.t3_v2[bin_t3])',
        hcat(raw_t3.t3_u3[bin_t3], raw_t3.t3_v3[bin_t3])') : Matrix{T}(undef, 2, 0)

    uv = hcat(vis_uv, v2_uv, t3_uv)
    uv_lam  = vcat(use_vis ? vis_lam  : T[], use_v2 ? v2_lam  : T[], use_t3 ? repeat(t3_lam,  3) : T[])
    uv_dlam = vcat(use_vis ? vis_dlam : T[], use_v2 ? v2_dlam : T[], use_t3 ? repeat(t3_dlam, 3) : T[])
    uv_mjd  = vcat(use_vis ? vis_mjd  : T[], use_v2 ? v2_mjd  : T[], use_t3 ? repeat(t3_mjd,  3) : T[])
    uv_baseline = vec(sqrt.(sum(uv.^2, dims=1)))
    nuv      = size(uv, 2)
    mean_mjd = isempty(uv_mjd) ? 0.0 : Float64(mean(uv_mjd))

    return BinData{T}(
        visamp, visamp_err, visphi, visphi_err,
        vis_mjd, vis_lam, vis_dlam, vis_flag, vis_baseline, vis_sta_index,
        indx_vis, nvisamp, nvisphi,
        v2, v2_err, v2_mjd, v2_lam, v2_dlam, v2_flag, v2_baseline, v2_sta_index,
        indx_v2, nv2,
        t3amp, t3amp_err, t3phi, t3phi_err,
        t3_mjd, t3_lam, t3_dlam, t3_flag, t3_baseline, t3_maxbaseline, t3_sta_index,
        indx_t3_1, indx_t3_2, indx_t3_3, nt3amp, nt3phi,
        flux, flux_err, flux_mjd, flux_lam, flux_dlam, flux_flag, flux_sta_index, nflux, flux_calibrated,
        uv, uv_lam, uv_dlam, uv_mjd, uv_baseline, nuv, mean_mjd,
        v2_corr,     v2_corr_idx,
        t3amp_corr,  t3amp_corr_idx,
        t3phi_corr,  t3phi_corr_idx,
        visamp_corr, visamp_corr_idx,
        visphi_corr, visphi_corr_idx,
        flux_corr,   flux_corr_idx,
        raw_vis.amptyp, raw_vis.phityp,
    )
end

# ---------------------------------------------------------------------------
# filter_bad_observables!: remove flagged/NaN/out-of-range data in-place,
# then prune the UV plane to only the points still referenced.
# ---------------------------------------------------------------------------
function filter_bad_observables!(bd::BinData{T};
        use_vis::Bool, use_v2::Bool, use_t3::Bool,
        force_full_vis::Bool, force_full_t3::Bool, special_filter_diffvis::Bool,
        cutoff_minv2::Real, cutoff_maxv2::Real,
        cutoff_mint3amp::Real, cutoff_maxt3amp::Real,
        filter_v2_snr_threshold::Real) where T

    good_uv_vis  = Int64[]; good_uv_v2   = Int64[]
    good_uv_t3_1 = Int64[]; good_uv_t3_2 = Int64[]; good_uv_t3_3 = Int64[]

    if use_vis
        va_ok = (.!isnan.(bd.visamp)) .& (.!isnan.(bd.visamp_err)) .& (bd.visamp_err .> 0)
        vp_ok = (.!isnan.(bd.visphi)) .& (.!isnan.(bd.visphi_err)) .& (bd.visphi_err .> 0)
        if special_filter_diffvis
            # Keep all VIS points; the cross-bin post-treatment will prune later
            vis_good = collect(1:length(bd.vis_flag))
        elseif force_full_vis
            vis_good = findall(.!bd.vis_flag .& (va_ok .& vp_ok))
        else
            vis_good = findall(.!bd.vis_flag .& (va_ok .| vp_ok))
        end
        good_uv_vis          = bd.indx_vis[vis_good]
        bd.visamp            = bd.visamp[vis_good];    bd.visamp_err    = bd.visamp_err[vis_good]
        bd.visphi            = bd.visphi[vis_good];    bd.visphi_err    = bd.visphi_err[vis_good]
        bd.nvisamp           = count(i -> !isnan(bd.visamp[i]) && !isnan(bd.visamp_err[i]) && bd.visamp_err[i] > 0, eachindex(bd.visamp))
        bd.nvisphi           = count(i -> !isnan(bd.visphi[i]) && !isnan(bd.visphi_err[i]) && bd.visphi_err[i] > 0, eachindex(bd.visphi))
        bd.vis_baseline      = bd.vis_baseline[vis_good]
        bd.vis_mjd           = bd.vis_mjd[vis_good];   bd.vis_lam       = bd.vis_lam[vis_good]
        bd.vis_dlam          = bd.vis_dlam[vis_good];  bd.vis_flag      = bd.vis_flag[vis_good]
        bd.vis_sta_index     = bd.vis_sta_index[:, vis_good]
        !isempty(bd.visamp_corr_idx) && (bd.visamp_corr_idx = bd.visamp_corr_idx[vis_good])
        !isempty(bd.visphi_corr_idx) && (bd.visphi_corr_idx = bd.visphi_corr_idx[vis_good])
    end

    if use_v2
        v2_good = findall(
            (.!bd.v2_flag) .& (bd.v2_err .> 0) .& (bd.v2_err .< 1) .&
            (bd.v2 .> cutoff_minv2) .& (bd.v2 .< cutoff_maxv2) .&
            .!isnan.(bd.v2) .& .!isnan.(bd.v2_err) .&
            (abs.(bd.v2 ./ bd.v2_err) .> filter_v2_snr_threshold))
        good_uv_v2           = bd.indx_v2[v2_good]
        bd.v2                = bd.v2[v2_good];         bd.v2_err        = bd.v2_err[v2_good]
        bd.v2_baseline       = bd.v2_baseline[v2_good]; bd.nv2          = length(bd.v2)
        bd.v2_mjd            = bd.v2_mjd[v2_good];    bd.v2_lam        = bd.v2_lam[v2_good]
        bd.v2_dlam           = bd.v2_dlam[v2_good];   bd.v2_flag       = bd.v2_flag[v2_good]
        bd.v2_sta_index      = bd.v2_sta_index[:, v2_good]
        !isempty(bd.v2_corr_idx) && (bd.v2_corr_idx = bd.v2_corr_idx[v2_good])
    end

    if use_t3
        ta_ok = (.!isnan.(bd.t3amp)) .& (.!isnan.(bd.t3amp_err)) .& (bd.t3amp_err .> 0)
        tp_ok = (.!isnan.(bd.t3phi)) .& (.!isnan.(bd.t3phi_err)) .& (bd.t3phi_err .> 0)
        t3_good = force_full_t3 ?
            findall(.!bd.t3_flag .& ta_ok .& tp_ok .& (bd.t3amp .> cutoff_mint3amp) .& (bd.t3amp .< cutoff_maxt3amp)) :
            findall(.!bd.t3_flag .& (ta_ok .| tp_ok))
        good_uv_t3_1         = bd.indx_t3_1[t3_good]
        good_uv_t3_2         = bd.indx_t3_2[t3_good]
        good_uv_t3_3         = bd.indx_t3_3[t3_good]
        bd.t3amp             = bd.t3amp[t3_good];      bd.t3amp_err     = bd.t3amp_err[t3_good]
        bd.t3phi             = bd.t3phi[t3_good];      bd.t3phi_err     = bd.t3phi_err[t3_good]
        bd.nt3amp            = count(i -> !isnan(bd.t3amp[i]) && !isnan(bd.t3amp_err[i]) && bd.t3amp_err[i] > 0, eachindex(bd.t3amp))
        bd.nt3phi            = count(i -> !isnan(bd.t3phi[i]) && !isnan(bd.t3phi_err[i]) && bd.t3phi_err[i] > 0, eachindex(bd.t3phi))
        bd.t3_baseline       = bd.t3_baseline[t3_good]; bd.t3_maxbaseline = bd.t3_maxbaseline[t3_good]
        bd.t3_mjd            = bd.t3_mjd[t3_good];    bd.t3_lam        = bd.t3_lam[t3_good]
        bd.t3_dlam           = bd.t3_dlam[t3_good];   bd.t3_flag       = bd.t3_flag[t3_good]
        bd.t3_sta_index      = bd.t3_sta_index[:, t3_good]
        !isempty(bd.t3amp_corr_idx) && (bd.t3amp_corr_idx = bd.t3amp_corr_idx[t3_good])
        !isempty(bd.t3phi_corr_idx) && (bd.t3phi_corr_idx = bd.t3phi_corr_idx[t3_good])
    end

    # Prune UV plane: keep only points still referenced by surviving observables
    uv_sel = falses(bd.nuv)
    isempty(good_uv_vis)  || (uv_sel[good_uv_vis]  .= true)
    isempty(good_uv_v2)   || (uv_sel[good_uv_v2]   .= true)
    isempty(good_uv_t3_1) || (uv_sel[good_uv_t3_1] .= true)
    isempty(good_uv_t3_2) || (uv_sel[good_uv_t3_2] .= true)
    isempty(good_uv_t3_3) || (uv_sel[good_uv_t3_3] .= true)

    iconv = cumsum(uv_sel); sel = findall(uv_sel)
    bd.uv          = bd.uv[:, sel];      bd.uv_lam      = bd.uv_lam[sel]
    bd.uv_dlam     = bd.uv_dlam[sel];    bd.uv_mjd      = bd.uv_mjd[sel]
    bd.uv_baseline = bd.uv_baseline[sel]; bd.nuv        = size(bd.uv, 2)

    isempty(good_uv_vis)  || (bd.indx_vis   = iconv[good_uv_vis])
    isempty(good_uv_v2)   || (bd.indx_v2    = iconv[good_uv_v2])
    isempty(good_uv_t3_1) || (bd.indx_t3_1  = iconv[good_uv_t3_1])
    isempty(good_uv_t3_2) || (bd.indx_t3_2  = iconv[good_uv_t3_2])
    isempty(good_uv_t3_3) || (bd.indx_t3_3  = iconv[good_uv_t3_3])
    bd.mean_mjd = isempty(bd.uv_mjd) ? 0.0 : Float64(mean(bd.uv_mjd))
    return bd
end

# ---------------------------------------------------------------------------
# remove_redundant_uv!(bd::BinData) — UV deduplication, BinData overload
# ---------------------------------------------------------------------------
function remove_redundant_uv!(bd::BinData{T}; uvtol::Float64=2e2) where T
    isempty(bd.uv) && return bd
    iconv, tokeep    = rm_redundance_kdtree(bd.uv, uvtol)
    bd.uv            = bd.uv[:, tokeep];  bd.uv_lam  = bd.uv_lam[tokeep]
    bd.uv_dlam       = bd.uv_dlam[tokeep]; bd.uv_mjd = bd.uv_mjd[tokeep]
    bd.uv_baseline   = bd.uv_baseline[tokeep]; bd.nuv = length(tokeep)
    !isempty(bd.indx_vis)   && (bd.indx_vis   = iconv[bd.indx_vis])
    !isempty(bd.indx_v2)    && (bd.indx_v2    = iconv[bd.indx_v2])
    !isempty(bd.indx_t3_1)  && (bd.indx_t3_1  = iconv[bd.indx_t3_1])
    !isempty(bd.indx_t3_2)  && (bd.indx_t3_2  = iconv[bd.indx_t3_2])
    !isempty(bd.indx_t3_3)  && (bd.indx_t3_3  = iconv[bd.indx_t3_3])
    return bd
end

# ---------------------------------------------------------------------------
# make_oidata: package a BinData + shared station info into OIdata{T}
# ---------------------------------------------------------------------------
function make_oidata(bd::BinData{T}, station_info, filename::String) where T
    return OIdata{T}(
        bd.visamp, bd.visamp_err, bd.visphi, bd.visphi_err,
        bd.vis_baseline, bd.vis_mjd, bd.vis_lam, bd.vis_dlam, bd.vis_flag,
        bd.v2, bd.v2_err, bd.v2_baseline, bd.v2_mjd, bd.mean_mjd,
        bd.v2_lam, bd.v2_dlam, bd.v2_flag,
        bd.t3amp, bd.t3amp_err, bd.t3phi, bd.t3phi_err,
        T[], T[],            # vonmises fields (computed separately)
        bd.t3_baseline, bd.t3_maxbaseline,
        bd.t3_mjd, bd.t3_lam, bd.t3_dlam, bd.t3_flag,
        bd.flux, bd.flux_err, bd.flux_mjd, bd.flux_lam, bd.flux_dlam, bd.flux_flag,
        bd.flux_sta_index, bd.flux_calibrated,
        bd.uv, bd.uv_lam, bd.uv_dlam, bd.uv_mjd, bd.uv_baseline,
        bd.nflux, bd.nvisamp, bd.nvisphi, bd.nv2, bd.nt3amp, bd.nt3phi, bd.nuv,
        bd.indx_vis, bd.indx_v2, bd.indx_t3_1, bd.indx_t3_2, bd.indx_t3_3,
        station_info.new_station_name, station_info.new_telescope_name, station_info.new_station_index,
        bd.vis_sta_index, bd.v2_sta_index, bd.t3_sta_index,
        bd.v2_corr,     bd.v2_corr_idx,
        bd.t3amp_corr,  bd.t3amp_corr_idx,
        bd.t3phi_corr,  bd.t3phi_corr_idx,
        bd.visamp_corr, bd.visamp_corr_idx,
        bd.visphi_corr, bd.visphi_corr_idx,
        bd.flux_corr,   bd.flux_corr_idx,
        bd.amptyp, bd.phityp,
        filename
    )
end

# ===========================================================================
# readoifits — main entry point
# ===========================================================================
"""
    readoifits(oifitsfile; kwargs...) -> Array{OIdata{T}, 2}

Read an OIFITS file and return a 2-D array of `OIdata{T}` indexed as
`[nwavbin, ntimebin]`. The simplest call returns a `1×1` array containing all
data in a single bin.

# Keyword arguments

## Target selection
- `targetname` — select a single target by name. Default: all targets combined.

## Spectral / temporal binning
- `spectralbin` — vector of `[λ_min, λ_max]` windows in metres. Default: single bin
  spanning all wavelengths.
- `temporalbin` — vector of `[mjd_min, mjd_max]` windows. Default: single bin
  spanning all epochs.
- `polychromatic` — if `true`, derive one spectral bin per instrument channel using
  midpoints between adjacent channel centres as boundaries. Overrides `spectralbin`.
- `merge_oi_wavelength` — if `true` and `polychromatic=true`, check whether all
  OI_WAVELENGTH tables share the same spectral channels (matching `eff_band` overlap)
  and, if so, deduplicate them to avoid redundant bins. Errors if tables are incompatible.
  Ignored when `spectralbin` is explicitly provided. Default: `false`.
- `splitting` — force multi-bin mode even when bin vectors equal `[[]]`.
- `get_specbin_file` — auto-derive spectral bins from the file. Default: `true`.
- `get_timebin_file` — auto-derive temporal bins from the file. Default: `true`.

## Observable selection
- `use_vis`, `use_v2`, `use_t3`, `use_flux` — load each observable type. Default: all `true`.
## Quality filtering
- `filter_bad_data` — apply quality cuts on load. Default: `true`.
- `force_full_vis` — require both visamp and visphi to be valid. Default: `false`.
- `force_full_t3` — require both t3amp and t3phi. Default: `false`.
- `cutoff_minv2`, `cutoff_maxv2` — V² range cut. Default: `(-1, 2.0)`.
- `cutoff_mint3amp`, `cutoff_maxt3amp` — T3amp range cut. Default: `(-1.0, 1.5)`.
- `filter_v2_snr_threshold` — minimum |V²/σ|. Default: `0.01`.
- `special_filter_diffvis` — differential visibility mode: keep only vis points common
  across all spectral bins.

## UV deduplication
- `redundance_remove` — merge UV points closer than `uvtol`. Default: `true`.
- `uvtol` — merge radius in cycles/rad (i.e. B/λ). Default: `200.0`.

## Numeric precision
- `T` — element type for all numeric arrays. Use `Float32` for ~50% memory reduction.
  Default: `Float64`.

## Output
- `warn` — print warnings about non-standard files. Default: `true`.
- `verbose` — print summary of loaded tables. Default: `false`.

# Example
```julia
data = readoifits("mystar.oifits")                         # all data, single bin
data = readoifits("mystar.oifits"; polychromatic=true)     # one bin per channel
data = readoifits("multi.oifits"; polychromatic=true, merge_oi_wavelength=true) # merge compatible tables
data = readoifits("mystar.oifits"; T=Float32)              # half memory
data = readoifits("multi.oifits"; targetname="Betelgeuse") # one target
```
"""
function readoifits(oifitsfile;
        targetname        = "",
        spectralbin       = [[]],
        temporalbin       = [[]],
        splitting         = false,
        polychromatic     = false,
        merge_oi_wavelength = false,
        get_specbin_file  = true,
        get_timebin_file  = true,
        redundance_remove = true,
        uvtol             = 2e2,
        filter_bad_data   = true,
        force_full_vis    = false,
        force_full_t3     = false,
        filter_v2_snr_threshold = 0.01,
        use_vis   = true, use_v2  = true, use_t3   = true,
        use_flux = true,
        cutoff_minv2    = -1,   cutoff_maxv2    = 2.0,
        cutoff_mint3amp = -1.0, cutoff_maxt3amp = 1.5,
        special_filter_diffvis = false,
        warn = true,
        verbose = true,
        T::Type{<:AbstractFloat} = Float64)

    if !isfile(oifitsfile)
        @error("readoifits could not locate the requested data file — please check path")
        return [[]]
    end

    ds = load_oifits(oifitsfile)

    # ---- Build OI_CORR lookup: corrname → symmetric SparseMatrixCSC ----------
    correlt = Dict{String, SparseMatrixCSC{T,Int64}}()
    try
        for ct in ds.correl
            correlt[ct.corrname] = _build_corr_sparse(ct; T)
        end
    catch
    end

    # ---- Validate and extract tables ----------------------------------------
    wavtables    = ds.instr;  isempty(wavtables)   && error("No OI_WAVELENGTH table in $oifitsfile")
    arraytables  = ds.array;  isempty(arraytables)  && error("No OI_ARRAY table in $oifitsfile")
    wavtableref  = [db.insname for db in wavtables]
    arraytableref = [db.arrname for db in arraytables]

    # ---- Target filtering ---------------------------------------------------
    tgt = ds.target
    all_target_ids   = [tgt[i].target_id for i in 1:length(tgt)]
    all_target_names = [tgt[i].target    for i in 1:length(tgt)]
    if !isempty(all_target_ids) && minimum(all_target_ids) == 0 && warn
        @warn("OI_TARGET does not follow OIFITSv2 standard — target indexing should start at 1, not 0.")
    end
    if targetname != ""
        targetid_filter = all_target_ids[findall(all_target_names .== targetname)]
        isempty(targetid_filter) && error("Target '$targetname' not found in $oifitsfile")
    else
        targetid_filter = unique(all_target_ids)
    end

    # ---- Observable table availability --------------------------------------
    fluxtables = ds.flux;  use_flux = use_flux && !isempty(fluxtables)
    vistables  = ds.vis;   use_vis  = use_vis  && !isempty(vistables)
    v2tables   = ds.vis2;  use_v2   = use_v2   && !isempty(v2tables)
    t3tables   = ds.t3;    use_t3   = use_t3   && !isempty(t3tables)

    # ---- Verbose: table summary ---------------------------------------------
    if verbose
        printstyled("File: $oifitsfile\n", color=:cyan)
        printstyled("  OI_ARRAY:      $(length(arraytables)) table(s)\n", color=:cyan)
        printstyled("  OI_WAVELENGTH: $(length(wavtables)) table(s)\n", color=:cyan)
        use_flux && printstyled("  OI_FLUX:       $(length(fluxtables)) table(s)\n", color=:cyan)
        use_vis  && printstyled("  OI_VIS:        $(length(vistables)) table(s)\n", color=:cyan)
        use_v2   && printstyled("  OI_VIS2:       $(length(v2tables)) table(s)\n", color=:cyan)
        use_t3   && printstyled("  OI_T3:         $(length(t3tables)) table(s)\n", color=:cyan)
        !isempty(correlt) && printstyled("  OI_CORR:       $(length(correlt)) matrix/matrices\n", color=:cyan)
    end

    # ---- Station indexing ---------------------------------------------------
    station_info = collect_station_info(arraytables, v2tables, t3tables, arraytableref; warn)
    ci = station_info.conversion_index
    so = station_info.station_index_offset

    # ---- Read and flatten all observable tables -----------------------------
    _ec = spzeros(T, 0, 0)   # empty corr matrix sentinel
    raw_flux = use_flux ? read_flux_tables(fluxtables, targetid_filter, wavtables, wavtableref, arraytableref, ci, so; correlt, T) :
               (; flux=T[], flux_err=T[], flux_mjd=T[], flux_lam=T[], flux_dlam=T[], flux_flag=Bool[], flux_sta_index=Int64[],
                  flux_calibrated=false, flux_corr=_ec, flux_corr_idx=Int64[])
    raw_vis  = use_vis  ? read_vis_tables(vistables,  targetid_filter, wavtables, wavtableref, arraytableref, ci, so; correlt, T) :
               (; visamp=T[], visamp_err=T[], visphi=T[], visphi_err=T[], vis_mjd=T[], vis_lam=T[], vis_dlam=T[], vis_flag=Bool[],
                  vis_uv=Matrix{T}(undef,0,2), vis_baseline=T[], vis_sta_index=Matrix{Int64}(undef,2,0),
                  visamp_corr=_ec, visamp_corr_idx=Int64[], visphi_corr=_ec, visphi_corr_idx=Int64[],
                  amptyp="", phityp="")
    raw_v2   = use_v2   ? read_v2_tables(v2tables,   targetid_filter, wavtables, wavtableref, arraytableref, ci, so; correlt, T) :
               (; v2=T[], v2_err=T[], v2_mjd=T[], v2_lam=T[], v2_dlam=T[], v2_flag=Bool[],
                  v2_uv=Matrix{T}(undef,0,2), v2_baseline=T[], v2_sta_index=Matrix{Int64}(undef,2,0),
                  v2_corr=_ec, v2_corr_idx=Int64[])
    raw_t3   = use_t3   ? read_t3_tables(t3tables,   targetid_filter, wavtables, wavtableref, arraytableref, ci, so; correlt, T) :
               (; t3amp=T[], t3amp_err=T[], t3phi=T[], t3phi_err=T[], t3_mjd=T[], t3_lam=T[], t3_dlam=T[], t3_flag=Bool[],
                  t3_u1=T[], t3_v1=T[], t3_u2=T[], t3_v2=T[], t3_u3=T[], t3_v3=T[], t3_baseline=T[], t3_maxbaseline=T[],
                  t3_sta_index=Matrix{Int64}(undef,3,0),
                  t3amp_corr=_ec, t3amp_corr_idx=Int64[], t3phi_corr=_ec, t3phi_corr_idx=Int64[])

    # ---- AMPTYP / PHITYP fallback via FITSIO headers -------------------------
    if use_vis && isempty(raw_vis.amptyp) && isempty(raw_vis.phityp)
        _at, _pt = _read_vis_header_keywords(oifitsfile)
        raw_vis = merge(raw_vis, (; amptyp=_at, phityp=_pt))
    end
    if use_vis && verbose && (!isempty(raw_vis.amptyp) || !isempty(raw_vis.phityp))
        printstyled("OI_VIS: AMPTYP=$(raw_vis.amptyp) PHITYP=$(raw_vis.phityp)\n", color=:cyan)
    end

    # ---- Determine spectral / temporal bins ---------------------------------
    splitting = splitting || polychromatic || temporalbin != [[]] || spectralbin != [[]]

    if temporalbin == [[]] && get_timebin_file
        mjd_all = vcat(
            use_vis  ? raw_vis.vis_mjd   : T[],
            use_v2   ? raw_v2.v2_mjd     : T[],
            use_t3   ? raw_t3.t3_mjd     : T[],
            use_flux ? raw_flux.flux_mjd : T[])
        temporalbin = [[minimum(mjd_all) - 0.001, maximum(mjd_all) + 0.001]]
    end

    if spectralbin == [[]] && get_specbin_file && !polychromatic
        lam_all  = vcat(use_vis ? raw_vis.vis_lam : T[], use_v2 ? raw_v2.v2_lam : T[], use_t3 ? raw_t3.t3_lam : T[], use_flux ? raw_flux.flux_lam : T[])
        dlam_all = vcat(use_vis ? raw_vis.vis_dlam : T[], use_v2 ? raw_v2.v2_dlam : T[], use_t3 ? raw_t3.t3_dlam : T[], use_flux ? raw_flux.flux_dlam : T[])
        spectralbin = [[minimum(lam_all) - 0.5*minimum(dlam_all[argmin(lam_all)]),
                        maximum(lam_all) + 0.5*maximum(dlam_all[argmax(lam_all)])]]
    end

    if polychromatic && get_specbin_file
        if length(wavtables) > 1 && merge_oi_wavelength
            eff_wave_all = merge_wav_tables(wavtables; T, verbose)
        else
            length(wavtables) > 1 &&
                @warn("Multiple OI_WAVELENGTH tables — set merge_oi_wavelength=true to auto-merge if channels match, or specify spectralbin manually.")
            eff_wave_all = sort(vcat([T.(db.eff_wave) for db in wavtables]...))
        end
        # Use midpoints between adjacent channel centers as bin boundaries.
        # This guarantees each stored lam value (= eff_wave[k]) falls in exactly one bin,
        # regardless of eff_band (which can be larger than the channel spacing).
        nch = length(eff_wave_all)
        if nch == 1
            hw = T(wavtables[1].eff_band[1]) / 2
            spectralbin = [[eff_wave_all[1] - hw, eff_wave_all[1] + hw]]
        else
            bounds = vcat(eff_wave_all[1] - (eff_wave_all[2] - eff_wave_all[1]) / 2,
                          (eff_wave_all[1:end-1] .+ eff_wave_all[2:end]) ./ 2,
                          eff_wave_all[end] + (eff_wave_all[end] - eff_wave_all[end-1]) / 2)
            spectralbin = [[bounds[k], bounds[k+1]] for k in 1:nch]
        end
    end

    nwavbin  = length(spectralbin)
    ntimebin = length(temporalbin)
    OIdataArr = Array{OIdata{T}}(undef, nwavbin, ntimebin)

    # ---- Main binning / filtering loop --------------------------------------
    for itimebin in 1:ntimebin, iwavbin in 1:nwavbin
        if splitting
            tlo, thi = temporalbin[itimebin][1], temporalbin[itimebin][2]
            wlo, whi = spectralbin[iwavbin][1],  spectralbin[iwavbin][2]
            bin_vis  = use_vis  ? (raw_vis.vis_mjd   .>= tlo .&& raw_vis.vis_mjd   .<= thi .&& raw_vis.vis_lam   .>= wlo .&& raw_vis.vis_lam   .<= whi) : Bool[]
            bin_v2   = use_v2   ? (raw_v2.v2_mjd     .>= tlo .&& raw_v2.v2_mjd     .<= thi .&& raw_v2.v2_lam     .>= wlo .&& raw_v2.v2_lam     .<= whi) : Bool[]
            bin_t3   = use_t3   ? (raw_t3.t3_mjd     .>= tlo .&& raw_t3.t3_mjd     .<= thi .&& raw_t3.t3_lam     .>= wlo .&& raw_t3.t3_lam     .<= whi) : Bool[]
            bin_flux = use_flux ? (raw_flux.flux_mjd .>= tlo .&& raw_flux.flux_mjd .<= thi .&& raw_flux.flux_lam .>= wlo .&& raw_flux.flux_lam .<= whi) : Bool[]
        else
            bin_vis  = use_vis  ? trues(length(raw_vis.visphi))   : Bool[]
            bin_v2   = use_v2   ? trues(length(raw_v2.v2))        : Bool[]
            bin_t3   = use_t3   ? trues(length(raw_t3.t3phi))     : Bool[]
            bin_flux = use_flux ? trues(length(raw_flux.flux))    : Bool[]
        end

        bd = slice_to_bin(raw_vis, raw_v2, raw_t3, raw_flux,
                          bin_vis, bin_v2, bin_t3, bin_flux;
                          use_vis, use_v2, use_t3, use_flux, T)

        filter_bad_data && filter_bad_observables!(bd;
            use_vis, use_v2, use_t3, force_full_vis, force_full_t3,
            special_filter_diffvis, cutoff_minv2, cutoff_maxv2,
            cutoff_mint3amp, cutoff_maxt3amp, filter_v2_snr_threshold)

        redundance_remove && remove_redundant_uv!(bd; uvtol)

        OIdataArr[iwavbin, itimebin] = make_oidata(bd, station_info, oifitsfile)
    end

    # ---- Post-treatment: cross-bin differential visibility filter -----------
    # Keep only VIS points that are valid (unflagged, finite err>0) across ALL spectral bins.
    if special_filter_diffvis
        for itimebin in 1:ntimebin
            # Build per-channel validity masks (all same length since per-bin filter kept everything)
            good_masks = [
                let d = OIdataArr[i,itimebin]
                    va_ok = (.!isnan.(d.visamp)) .& (.!isnan.(d.visamp_err)) .& (d.visamp_err .> 0)
                    vp_ok = (.!isnan.(d.visphi)) .& (.!isnan.(d.visphi_err)) .& (d.visphi_err .> 0)
                    .!d.vis_flag .& va_ok .& vp_ok
                end
                for i in 1:nwavbin]
            # Point must be good in ALL channels
            common_good = reduce(.&, good_masks)
            indx = findall(common_good)
            for i in 1:nwavbin
                d = OIdataArr[i,itimebin]
                d.nvisamp       = count(j -> !isnan(d.visamp[j]) && !isnan(d.visamp_err[j]) && d.visamp_err[j] > 0, indx)
                d.nvisphi       = count(j -> !isnan(d.visphi[j]) && !isnan(d.visphi_err[j]) && d.visphi_err[j] > 0, indx)
                d.visamp        = d.visamp[indx];        d.visamp_err    = d.visamp_err[indx]
                d.visphi        = d.visphi[indx];        d.visphi_err    = d.visphi_err[indx]
                d.vis_baseline  = d.vis_baseline[indx];  d.vis_mjd       = d.vis_mjd[indx]
                d.vis_lam       = d.vis_lam[indx];       d.vis_dlam      = d.vis_dlam[indx]
                d.vis_sta_index = d.vis_sta_index[:, indx]
                d.vis_flag      = d.vis_flag[indx];      d.indx_vis      = d.indx_vis[indx]
                !isempty(d.visamp_corr_idx) && (d.visamp_corr_idx = d.visamp_corr_idx[indx])
                !isempty(d.visphi_corr_idx) && (d.visphi_corr_idx = d.visphi_corr_idx[indx])
            end
        end
    end

    return OIdataArr
end

# ===========================================================================
# Convenience wrappers
# ===========================================================================

"""
    readoifits_multiepochs(oifitsfiles; polychromatic=false, kwargs...) -> Matrix{OIdata{T}}

Read a list of OIFITS files, one per epoch, and return a `Matrix{OIdata{T}}`
of size `(nwav, nepochs)`.

- Without `polychromatic`: each file yields a single spectral bin → `1 × nepochs`
- With `polychromatic=true`: each file is split into wavelength channels → `nwav × nepochs`

Mean MJD per epoch is available via `data[1,i].mean_mjd`.

Prints a summary line per file. Passes `filter_bad_data`, `force_full_t3`, and
`polychromatic` through to `readoifits`.
"""
function readoifits_multiepochs(oifitsfiles; filter_bad_data=true, force_full_t3=false,
                                polychromatic=false,
                                T::Type{<:AbstractFloat}=Float64)
    per_file = [readoifits(f; filter_bad_data, force_full_t3, polychromatic, T) for f in oifitsfiles]
    data = hcat(per_file...)   # each is nwav×1, result is nwav×nepochs
    for i in eachindex(oifitsfiles)
        println(oifitsfiles[i], "\t MJD: ", data[1,i].mean_mjd,
                "\t nV2 = ", data[1,i].nv2, "\t nT3amp = ", data[1,i].nt3amp,
                "\t nT3phi = ", data[1,i].nt3phi)
    end
    return data
end


# ===========================================================================
# FITS file utilities
# ===========================================================================

function readfits(fitsfile; normalize=false, vectorize=false)
    x = read((FITS(fitsfile))[1])
    normalize  && (x ./= sum(x))
    vectorize  && (x = vec(x))
    return x
end

function writefits(data, fitsfile; pixsize=-1)
    f = FITS(fitsfile, "w")
    if pixsize != -1
        header = FITSHeader(
            ["CDELT1","CDELT2","CRVAL1","CRVAL2","CRPIX1","CRPIX2"],
            [-(pixsize/1000.0)/206265.0, (pixsize/1000.0)/206265.0, 0.0, 0.0,
             size(data,1)/2, size(data,1)/2],
            ["Radians per Pixel","Radians per Pixel",
             "X-coordinate of reference pixel","Y-coordinate of reference pixel",
             "reference pixel in X","reference pixel in Y"])
        write(f, data, header=header)
    else
        write(f, data)
    end
    close(f)
end

function updatefits_aspro(fitsfile_in, fitsfile_out, res)
    f      = FITS(fitsfile_in)
    data   = read(f[1]);  header = read_header(f[1])
    header["CDELT1"] = -res/1000.0;  header["CDELT2"] =  res/1000.0
    header["CRVAL1"] = 0.0;          header["CRVAL2"] = 0.0
    header["CRPIX1"] = size(data,1)/2; header["CRPIX2"] = size(data,1)/2
    header["CUNIT1"] = "arcsec";     header["CUNIT2"] = "arcsec"
    set_comment!(header, "CDELT1", "Arcseconds per Pixel")
    set_comment!(header, "CDELT2", "Arcseconds per Pixel")
    set_comment!(header, "CRVAL1", "X-coordinate of reference pixel")
    set_comment!(header, "CRVAL2", "Y-coordinate of reference pixel")
    set_comment!(header, "CRPIX1", "reference pixel in X")
    set_comment!(header, "CRPIX2", "reference pixel in Y")
    fout = FITS(fitsfile_out, "w"); write(fout, data, header=header)
    close(f); close(fout)
end

# ===========================================================================
# Error inflation / oifits_prep
# ===========================================================================

function oifits_prep(data::OIdata{T};
        min_v2_err_add=0, min_v2_err_rel=0, v2_err_mult=1,
        min_t3amp_err_add=0, min_t3amp_err_rel=0, t3amp_err_mult=1,
        min_t3phi_err_add=0, t3phi_err_mult=1, quad=false) where T

    if quad
        data.v2_err = sqrt.((data.v2 .* min_v2_err_rel).^2 .+
                            (v2_err_mult .* data.v2_err).^2 .+ T(min_v2_err_add)^2)
    else
        temperr = v2_err_mult .* data.v2_err
        newerr  = abs.(data.v2 .* min_v2_err_rel) .+ min_v2_err_add
        temperr[newerr .> temperr] .= newerr[newerr .> temperr]
        data.v2_err = temperr
    end

    if quad
        data.t3amp_err = sqrt.((data.t3amp .* min_t3amp_err_rel).^2 .+
                               (t3amp_err_mult .* data.t3amp_err).^2 .+ T(min_t3amp_err_add)^2)
    else
        temperr = t3amp_err_mult .* data.t3amp_err
        newerr  = abs.(data.t3amp .* min_t3amp_err_rel) .+ min_t3amp_err_add
        temperr[newerr .> temperr] .= newerr[newerr .> temperr]
        data.t3amp_err = temperr
    end

    if quad
        data.t3phi_err = sqrt.((t3phi_err_mult .* data.t3phi_err).^2 .+ T(min_t3phi_err_add)^2)
    else
        temperr = t3phi_err_mult .* data.t3phi_err
        newerr  = T(min_t3phi_err_add) .* ones(T, length(data.t3phi_err))
        temperr[newerr .> temperr] .= newerr[newerr .> temperr]
        data.t3phi_err = temperr
    end
    return data
end

function oifits_prep(data::AbstractVector{<:OIdata}; kwargs...)
    oifits_prep.(data; kwargs...)
    return data
end

# ===========================================================================
# Utilities
# ===========================================================================

function list_oifits_targets(oifitsfile)
    if !isfile(oifitsfile)
        @error("Could not locate the requested data file — please check path")
        return [[]]
    end
    ds = load_oifits(oifitsfile)
    return unique([ds.target[i].target for i in 1:length(ds.target)])
end
