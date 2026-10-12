# The Tools tab's operations on OIFITS files and datasets: the header inspector, the validator
# and its repaired copy, and the splits.

"""
    oifits_headers(path) -> Vector{NamedTuple}

Every HDU of a FITS file, in order, as `(index, extname, kind, revn, insname, arrname, nrows,
dims, columns, cards)`:

- `kind` — `:primary`, `:image` or `:table`; `extname` is `""` when the HDU has none;
- `revn` — `OI_REVN` (0 when absent); `insname`, `arrname` — `INSNAME`, `ARRNAME` (`""`);
- `nrows` — rows of a table, 0 otherwise; `dims` — an image's axis lengths (`NAXISn`), empty for
  a table or a header-only HDU;
- `columns` — one `(name, format, unit, dims)` per table column, from `TTYPEn`, `TFORMn`,
  `TUNITn`, `TDIMn` (`""` when absent);
- `cards` — every header card in order, COMMENT and HISTORY included, as `(name, value,
  comment)`; a commentary card's text is its `value`, and its `comment` is `""`.

Works on any FITS file, images included, whether or not OIFITS.jl reads it. A ZIP archive holding
one FITS file (some pipelines deliver their `.fits` files zipped) shows that file.
"""
function oifits_headers(path::AbstractString)
    isfile(path) || throw(ArgumentError("$(path): no such file"))
    _is_zip(path) || return _headers(path)
    tmp, _ = _unzip_fits(path)
    try
        return _headers(tmp)
    finally
        rm(tmp; force = true)
    end
end

function _headers(path)
    AstroFITS.FitsFile(path) do ff
        map(enumerate(ff)) do (i, hdu)
            hdr = AstroFITS.FitsHeader(hdu)
            cards = [(name = String(c.name), value = _card_value(c),
                      comment = c.type == AstroFITS.FITS_COMMENT ? "" : String(c.comment))
                     for c in hdr if c.name != "END"]
            kw = Dict(c.name => c.value for c in cards)
            str(k) = (v = get(kw, k, ""); String(v isa AbstractString ? strip(v) : string(v)))
            istable = hdu isa AstroFITS.FitsTableHDU
            ncol = istable ? hdu.ncols : 0
            columns = [(name = str("TTYPE$k"), format = str("TFORM$k"), unit = str("TUNIT$k"),
                        dims = str("TDIM$k")) for k in 1:ncol]
            naxis = istable ? 0 : get(kw, "NAXIS", 0)
            dims = [Int(get(kw, "NAXIS$k", 0)) for k in 1:(naxis isa Integer ? naxis : 0)]
            revn = get(kw, "OI_REVN", 0)
            (index = i, extname = String(something(hdu.extname, "")),
             kind = i == 1 ? :primary : istable ? :table : :image,
             revn = revn isa Integer ? Int(revn) : 0, insname = str("INSNAME"),
             arrname = str("ARRNAME"), nrows = istable ? Int(hdu.nrows) : 0,
             dims, columns, cards)
        end
    end
end

# A card's value as a plain Julia value; commentary cards (COMMENT, HISTORY, blank) carry text.
function _card_value(c)
    c.type == AstroFITS.FITS_COMMENT && return String(c.comment)
    v = c.value()
    v isa AbstractString ? String(v) : v === nothing ? "" : v
end

"""
    validate_oifits(path) -> Vector{NamedTuple}

Check an OIFITS file and report what is wrong with it, as `(severity, table, message)` with
`severity` one of `:error` (the file breaks the standard or cannot be read), `:warning` (it reads,
but some data are suspect) and `:info`. An empty report means nothing was found.

The file is first read HDU by HDU (as [`oifits_headers`](@ref) does), then as `readoifits` opens
it ([`load_oifits`](@ref), table by table), and the data are checked:

- everything the opening had to work around, as warnings -- a table that does not read at the
  revision it declares (a missing column, a revision that does not exist), GRAVITY's `FLUX` read
  as OI_FLUX's `FLUXDATA`, one half of an amplitude/phase pair missing, rows whose TARGET_ID is not
  in OI_TARGET (dropped); a table that cannot be read at all, as an error. A ZIP archive is
  checked through the FITS file it holds, and said to be one;
- INSNAME / ARRNAME with no OI_WAVELENGTH / OI_ARRAY; two OI_WAVELENGTH tables with one INSNAME;
  OI_WAVELENGTH tables no data table uses; stations absent from their OI_ARRAY; TARGET_IDs absent
  from OI_TARGET;
- AMPTYP / PHITYP disagreeing between OI_VIS tables (`readoifits` keeps one label per dataset);
- per table: empty tables; unflagged points with NaN values or NaN / non-positive errors, or one
  half of a VIS / T3 pair carried with a non-finite error; the flagged fraction;
- what a repaired copy ([`repair_oifits`](@ref)) would not keep: HDUs other than the OIFITS
  tables, several targets.
"""
function validate_oifits(path::AbstractString)
    rep = NamedTuple{(:severity, :table, :message),Tuple{Symbol,String,String}}[]
    add(sev, tab, msg) = push!(rep, (severity = sev, table = String(tab), message = String(msg)))
    isfile(path) || (add(:error, "", "no such file"); return rep)
    hdus = try
        oifits_headers(path)
    catch e
        add(:error, "", "not a readable FITS file: " * sprint(showerror, e))
        return rep
    end
    label(h) = isempty(h.insname) ? "$(h.extname) (HDU $(h.index))" :
                                    "$(h.extname) $(h.insname) (HDU $(h.index))"

    # ── the tables, from their headers ──────────────────────────────────────────────
    for h in hdus
        h.index == 1 && continue
        startswith(h.extname, "OI_") ||
            add(:info, label(h), "not an OIFITS table; a repaired copy does not keep it")
        h.kind == :table && h.nrows == 0 && startswith(h.extname, "OI_") &&
            add(:warning, label(h), "empty table")
    end
    wav = [h for h in hdus if h.extname == "OI_WAVELENGTH"]
    for n in unique(h.insname for h in wav)
        count(h -> h.insname == n, wav) > 1 &&
            add(:error, "OI_WAVELENGTH", "several OI_WAVELENGTH tables with INSNAME \"$(n)\"")
    end

    # ── indices the standard numbers from 1 ──────────────────────────────────────────
    for (ext, col, what, n, revn) in _index_columns(path)
        bad = sort(unique(filter(<(1), n)))
        isempty(bad) && continue
        add(revn >= 2 ? :error : :warning, ext,
            "$(col) $(join(bad, ", ")): OIFITS v2 numbers $(what) from 1" *
            (revn >= 2 ? "" : " (this table is revision $(revn))"))
    end

    # ── as readoifits opens it, and what that had to work around ─────────────────────
    notes = String[]
    ds = try
        load_oifits(path; notes)
    catch e
        add(:error, "", "cannot be read: " * sprint(showerror, e))
        nothing
    end
    for n in notes
        if startswith(n, "a ZIP archive")
            add(:info, "", n)
        else
            tab, msg = split(n, ": "; limit = 2)
            add(:warning, tab, msg)
        end
    end
    ds === nothing && return rep

    datatabs = [(t, "OI_VIS") for t in ds.vis] ∪ [(t, "OI_VIS2") for t in ds.vis2] ∪
               [(t, "OI_T3") for t in ds.t3] ∪ [(t, "OI_FLUX") for t in ds.flux]
    tlabel(t, ext) = "$(ext) $(_safe_string(t, :insname))"
    insnames = Set(String(w.insname) for w in ds.instr)
    arrays = Dict(String(a.arrname) => Set(Int.(a.sta_index)) for a in ds.array)
    targets = Set(Int(ds.target[i].target_id) for i in 1:length(ds.target))
    used = Set{String}()
    for (t, ext) in datatabs
        ins = _safe_string(t, :insname); push!(used, ins)
        ins in insnames || add(:error, tlabel(t, ext), "INSNAME \"$(ins)\" has no OI_WAVELENGTH table")
        arr = _safe_string(t, :arrname)
        sta = _safe_field(t, :sta_index)
        if !isempty(arr) && !haskey(arrays, arr)
            add(:error, tlabel(t, ext), "ARRNAME \"$(arr)\" has no OI_ARRAY table")
        elseif !isempty(arr) && sta !== nothing
            miss = setdiff(Set(Int.(vec(sta))), arrays[arr])
            isempty(miss) || add(:warning, tlabel(t, ext),
                "stations $(join(sort(collect(miss)), ", ")) have no entry in OI_ARRAY \"$(arr)\" " *
                "(no telescope or station name)")
        end
        tid = setdiff(Set(Int.(t.target_id)), targets)
        isempty(tid) || add(:error, tlabel(t, ext), "TARGET_ID $(sort(collect(tid))) not in OI_TARGET")
        _check_values!(add, t, ext, tlabel(t, ext))
    end
    for w in ds.instr
        String(w.insname) in used ||
            add(:info, "OI_WAVELENGTH $(w.insname)", "no data table refers to it")
    end
    vis = collect(ds.vis)
    for (f, what) in ((:amptyp, "AMPTYP"), (:phityp, "PHITYP"))
        vals = unique(_canon_obstype(_safe_string(t, f)) for t in vis)
        length(vals) > 1 && add(:warning, "OI_VIS",
            "the OI_VIS tables disagree on $(what): " * join(repr.(vals), ", ") *
            "; readoifits keeps one label for the dataset (read them by instrument)")
    end
    length(targets) > 1 && add(:info, "OI_TARGET",
        "$(length(targets)) targets; a repaired copy holds one, so repair one target at a time")
    return rep
end

# Values and errors of one data table's unflagged points; the flagged fraction. A value with a
# NaN, infinite or non-positive error is suspect. A NaN value is how a table says "not
# measured": fine for one half of a VIS / T3 pair when the other half is there (GRAVITY's T3AMP),
# suspect when the point has no value at all and is not flagged.
function _check_values!(add, t, ext, lab)
    flag = Bool.(t.flag)
    pairs = ext == "OI_VIS"  ? (("VISAMP", t.visamp, t.visamperr), ("VISPHI", t.visphi, t.visphierr)) :
            ext == "OI_T3"   ? (("T3AMP", t.t3amp, t.t3amperr), ("T3PHI", t.t3phi, t.t3phierr)) :
            ext == "OI_VIS2" ? (("VIS2DATA", t.vis2data, t.vis2err),) :
                               (("FLUXDATA", t.fluxdata, t.fluxerr),)
    n = length(flag)
    n == 0 && return
    nf = count(flag)
    nf > 0 && add(:info, lab, "$(nf) of $(n) points flagged ($(round(100nf / n; digits = 1)) %)")
    ok = .!flag
    for (name, v, e) in pairs
        k = count(ok .& isfinite.(v) .& .!(isfinite.(e) .& (e .> 0)))
        k > 0 && add(:warning, lab, "$(k) unflagged $(name) values have a NaN, infinite or " *
                     "non-positive error")
    end
    novalue = reduce((a, b) -> a .& b, (.!isfinite.(v) for (_, v, _) in pairs))
    k = count(ok .& novalue)
    k > 0 && add(:warning, lab, "$(k) unflagged points carry no value")
end

"""
    repair_oifits(path, out; targetname = "", overwrite = false) -> out

Write a repaired copy of an OIFITS file: read by `readoifits` (every point kept, uv points not
merged) and written by [`write_oifits`](@ref), which rebuilds every table as OIFITS v2 — the
missing columns filled (e.g. FOV in OI_ARRAY), each instrument under its own name, OI_TARGET at
revision 2. Never over `path` itself.

What the copy does not keep: HDUs other than the OIFITS tables, columns OITOOLS does not read,
flagged points, OI_CORR, and whatever [`load_oifits`](@ref) had to leave out (rows whose
TARGET_ID is not in OI_TARGET). A ZIP archive's FITS file is repaired into a plain FITS file. A file with several targets is repaired one target at a time
(`targetname`). [`validate_oifits`](@ref) reports these before one repairs.
"""
function repair_oifits(path::AbstractString, out::AbstractString; targetname::AbstractString = "",
                       overwrite::Bool = false)
    _same_file(path, out) && throw(ArgumentError("the repaired copy cannot replace $(path)"))
    if isempty(targetname)
        names = list_oifits_targets(path)
        length(names) > 1 && throw(ArgumentError("$(path) holds $(length(names)) targets; " *
            "repair one at a time with targetname = one of " * join(repr.(names), ", ")))
    end
    d = readoifits(path; T = Float64, warn = false, verbose = false, targetname,
                   redundance_remove = false)
    return write_oifits(out, d; template = path, overwrite)
end

# Splitting a dataset into parts: by time, by instrument (from the points' table tags, so no
# re-reading), by target (re-read per target, since a dataset does not tag targets).

"""
    split_by_time(data; gap = 0.5, binwidth = nothing, edges = nothing) -> (parts, mjd)

Split an `OIdata` into epochs. The epochs come from the times of EVERY observable (V², T3, VIS,
flux), so a closure phase measured between two V² epochs is kept.

- `gap` (days, default): a new epoch wherever consecutive times are more than `gap` apart;
- `binwidth` (days): fixed windows from the first time, `[t₀ + kw, t₀ + (k+1)w)`;
- `edges`: given windows, a vector of `(lo, hi)` MJD pairs, inclusive (a point inside two
  windows goes to both).

`parts` holds one `OIdata` per non-empty epoch, made by [`filter_data`](@ref); `mjd` is the mean
time of each part's points. Points are assigned by their own Float64 MJDs (`v2_mjd`, …), so the
split is exact whatever `T` is. Within a part, the model time `uv_mjd` is stored in `T`: at
`Float32` it resolves ~5.6 min, which matters only for a `\$MJD` model over a short part.
"""
function split_by_time(data::OIdata; gap::Real = 0.5, binwidth = nothing, edges = nothing)
    t = sort!(vcat(data.v2_mjd, data.t3_mjd, data.vis_mjd, data.flux_mjd))
    isempty(t) && throw(ArgumentError("the data hold no point"))
    # windows as (lo, hi, closed_hi)
    windows = if edges !== nothing
        [(Float64(e[1]), Float64(e[2]), true) for e in edges]
    elseif binwidth !== nothing
        binwidth > 0 || throw(ArgumentError("binwidth must be positive"))
        n = floor(Int, (t[end] - t[1]) / binwidth) + 1
        [(t[1] + (k - 1) * binwidth, t[1] + k * binwidth, false) for k in 1:n]
    else
        cut = findall(>(gap), diff(t))
        [(t[a], t[b], true) for (a, b) in zip([1; cut .+ 1], [cut; length(t)])]
    end
    parts = OIdata[]; mjd = Float64[]
    for (lo, hi, closed) in windows
        inw(m) = lo <= m && (closed ? m <= hi : m < hi)
        out(v) = findall(!inw, v)
        n = count(inw, data.v2_mjd) + count(inw, data.t3_mjd) + count(inw, data.vis_mjd) +
            count(inw, data.flux_mjd)
        n == 0 && continue
        push!(parts, filter_data(data, [Int64[], out(data.vis_mjd), out(data.v2_mjd),
                                        out(data.t3_mjd), out(data.flux_mjd)]))
        push!(mjd, sum(m for v in (data.v2_mjd, data.t3_mjd, data.vis_mjd, data.flux_mjd)
                       for m in v if inw(m); init = 0.0) / n)
    end
    return [p for p in parts]::Vector{typeof(data)}, mjd
end

"""
    split_by_instrument(data) -> Vector{Pair{String, OIdata}}

Split an `OIdata` by instrument (the INSNAME of each point's table, `data.tables`), in the order
the tables were read, without re-reading the file. Works on a filtered, merged or binned dataset.
Points whose table is unknown are left out.
"""
function split_by_instrument(data::OIdata)
    ins(t) = t == 0 ? "" : data.tables[t].insname
    names = unique(filter(!isempty, [ins(t) for tab in (data.vis_tab, data.v2_tab, data.t3_tab,
                                                        data.flux_tab) for t in tab]))
    return [n => filter_data(data, [Int64[], findall(t -> ins(t) != n, data.vis_tab),
                                    findall(t -> ins(t) != n, data.v2_tab),
                                    findall(t -> ins(t) != n, data.t3_tab),
                                    findall(t -> ins(t) != n, data.flux_tab)]) for n in names]
end

"""
    split_by_target(path; kwargs...) -> Vector{Pair{String, Matrix{OIdata}}}

Read an OIFITS file once per target (`readoifits(path; targetname, kwargs...)`), in the order of
its OI_TARGET table.
"""
split_by_target(path::AbstractString; kwargs...) =
    [String(n) => readoifits(path; targetname = n, kwargs...) for n in list_oifits_targets(path)]

# OI_TARGET's TARGET_ID and OI_ARRAY's STA_INDEX columns as the file stores them, with each table's
# OI_REVN: `(table, column, what is numbered, values, revn)`. A ZIP archive's FITS file is read.
function _index_columns(path)
    tmp = _is_zip(path) ? first(_unzip_fits(path)) : nothing
    try
        AstroFITS.FitsFile(something(tmp, path)) do f
            out = Tuple{String,String,String,Vector{Int},Int}[]
            for h in f
                h isa AstroFITS.FitsTableHDU || continue
                ext = something(h.extname, "")
                col, what = ext == "OI_TARGET" ? ("TARGET_ID", "targets") :
                            ext == "OI_ARRAY"  ? ("STA_INDEX", "stations") : ("", "")
                (isempty(col) || !(col in h.column_names)) && continue
                arr = get(String, h, "ARRNAME", "")
                push!(out, (isempty(arr) ? ext : "$(ext) $(arr)", col, what,
                            Int.(vec(read(h, col))), get(Int, h, "OI_REVN", 0)))
            end
            out
        end
    catch
        Tuple{String,String,String,Vector{Int},Int}[]
    finally
        tmp === nothing || rm(tmp; force = true)
    end
end
