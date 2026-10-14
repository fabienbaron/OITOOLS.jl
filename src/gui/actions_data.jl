# Actions on data.
#
# Every function here follows the same contract: change the session, and append the exact
# OITOOLS call that would have done the same thing. The contract is enforced by the
# round-trip test, which exports the log, runs it, and compares against the session.

# readoifits' own defaults, so exported scripts show only what the user actually changed.
const _READOIFITS_DEFAULTS = Dict{Symbol,Any}(
    :targetname => "", :insname => "", :polychromatic => false, :merge_oi_wavelength => false,
    :filter_bad_data => true, :use_vis => true, :use_v2 => true, :use_t3 => true,
    :use_flux => true, :uvtol => 2e2, :warn => true, :verbose => true,
)

"""
    load_dataset!(session, path; kwargs...) -> DatasetEntry

Load an OIFITS file into the session and log the `readoifits` call that did it. `path` may
be a list of files, read as one dataset (`readoifits(files)`, which merges them point by point
through `merge_oidata`).

`kwargs` are passed straight through to `readoifits`, so anything it accepts works here
(`targetname`, `polychromatic`, `merge_oi_wavelength`, `filter_bad_data`, `T`, …).

Two behaviours worth knowing, both inherited from `readoifits`:

  * it returns an `Array{OIdata,2}` of `(nwavbin, ntimebin)`, kept in that shape here;
  * on a missing file it returns `[[]]` rather than throwing, which would surface as a
    baffling error much later — so this checks first and throws immediately.
"""
function load_dataset!(s::Session, path::Union{AbstractString,AbstractVector{<:AbstractString}};
                       kwargs...)
    files = path isa AbstractString ? [String(path)] : String.(collect(path))
    isempty(files) && throw(ArgumentError("no file to read"))
    for f in files
        isfile(f) || throw(ArgumentError("no such OIFITS file: $(repr(f))"))
    end
    kw = Dict{Symbol,Any}(kwargs)
    data = OITOOLS.readoifits(length(files) == 1 ? files[1] : files; kwargs...)
    data isa AbstractArray || throw(ArgumentError("readoifits did not return data for $(repr(path))"))

    name = _next_name(s.datasets, "data")
    entry = DatasetEntry(name, files[1], data, kw,
                         get(kw, :spectralbin, nothing), get(kw, :temporalbin, nothing), nothing,
                         "", length(files) == 1 ? String[] : files)
    push!(s.datasets, entry)

    # `Pair{Symbol,Any}[...]`, not a bare comprehension: over an EMPTY `kw` the comprehension
    # infers `Vector{Any}`, which `_kwargs(::Vector{<:Pair})` does not accept — so loading a
    # file with no read options at all threw a MethodError rather than logging the plain call.
    kwstr = _kwargs(Pair{Symbol,Any}[k => v for (k, v) in sort(collect(kw), by = first)];
                    defaults = _READOIFITS_DEFAULTS)
    call  = isempty(kwstr) ? "readoifits($(_literal(sources(entry))))" :
                             "readoifits($(_literal(sources(entry))); $kwstr)"
    log!(s, "$name = $call"; note = "load $(source_label(entry))", binding = name)
    return entry
end

"""
    filter_dataset!(session, name; read = Dict(), prep = Dict(), spec = "", kwargs...)
        -> DatasetEntry

Filter a loaded dataset into a new entry, leaving the original untouched, and log the calls
that do it.

- `kwargs` go to `set_data_filter`, bin by bin: every bin is filtered by its own discard lists
  and the dataset keeps its shape.
- `read`: options for `readoifits` (`insname`, `targetname`, `use_v2`, `ignore_flags`, …).
  What they select is decided when the file is read, so when any is given the file is read
  again with the source's own options plus these, and the filter applies to that.
- `prep`: error-bar floors and multipliers for `oifits_prep`, applied to the result. They
  change the data rather than remove any of it.
- `drop`: observables removed by `filter_data(…; drop)` (`:t3amp`, `:visphi`, …): one half of an
  amplitude/phase pair is marked "not measured", the other kept.
- `spec`: the Filter section's settings, kept on the entry so the section can show them again.

Filtering is destructive-by-copy, so keeping the source is what gives the GUI undo for free:
choosing the source again is the undo.
"""
function filter_dataset!(s::Session, name::AbstractString;
                         read::AbstractDict = Dict{Symbol,Any}(),
                         prep::AbstractDict = Dict{Symbol,Any}(),
                         drop::AbstractVector = Symbol[],
                         spec::AbstractString = "", kwargs...)
    src = dataset(s, name)
    kw  = Dict{Symbol,Any}(kwargs)
    ropts = merge(src.readopts, Dict{Symbol,Any}(read))
    base = isempty(read) ? src.data :
           OITOOLS.readoifits(sources(src); ropts..., warn = false, verbose = false)
    base isa AbstractArray || throw(ArgumentError("readoifits returned nothing for $(source_label(src))"))
    out = map(d -> OITOOLS.filter_data(d, OITOOLS.set_data_filter(d; kwargs...); drop), base)
    isempty(prep) || foreach(d -> OITOOLS.oifits_prep(d; prep...), out)

    newname = _next_name(s.datasets, "data")
    entry = DatasetEntry(newname, src.path, out, ropts, src.spectralbin, src.temporalbin,
                         src.name, String(spec), src.files)
    push!(s.datasets, entry)

    # Typed: over an empty dict a bare comprehension is a `Vector{Any}`, which `_kwargs` refuses.
    kwstr = _kwargs(Pair{Symbol,Any}[k => v for (k, v) in sort(collect(kw), by = first)])
    call  = isempty(kwstr) ? "set_data_filter(d)" : "set_data_filter(d; $kwstr)"
    from  = if isempty(read)
        name
    else
        rstr = _kwargs(Pair{Symbol,Any}[k => v for (k, v) in sort(collect(ropts), by = first)];
                       defaults = _READOIFITS_DEFAULTS)
        isempty(rstr) ? "readoifits($(_literal(sources(src))))" :
                        "readoifits($(_literal(sources(src))); $rstr)"
    end
    dstr = isempty(drop) ? "" : "; drop = $(_literal(Symbol.(drop)))"
    log!(s, "$(newname) = map(d -> filter_data(d, $call$dstr), $from)";
         note = "filter $(src.name)", binding = newname)
    if !isempty(prep)
        pstr = _kwargs(Pair{Symbol,Any}[k => v for (k, v) in sort(collect(prep), by = first)])
        log!(s, "foreach(d -> oifits_prep(d; $pstr), $newname)"; note = "error bars of $newname")
    end
    return entry
end

# ── the Filter section's settings ─────────────────────────────────────────────────────

const _ERR_ROWS = ("v2", "t3amp", "t3phi", "visamp", "visphi", "flux")
_is_phase(o) = o in ("t3phi", "visphi")

"""
    _nights(mjd) -> Vector{(; lo, hi, n)}

The observing nights among a set of epochs: runs of sorted MJDs with no gap longer than 12
hours, each with its first and last epoch and its number of epochs. A UT day is not a night: one
that crosses 0h UT is one night here, not two. Twelve hours separates nights at any site, while
nothing inside one night is that far apart.
"""
function _nights(mjd)
    t = sort(Float64.(vec(mjd)))
    out = NamedTuple{(:lo, :hi, :n), Tuple{Float64, Float64, Int}}[]
    isempty(t) && return out
    lo = prev = t[1]; n = 0
    for x in t
        if x - prev > 0.5
            push!(out, (; lo, hi = prev, n)); lo = x; n = 0
        end
        n += 1; prev = x
    end
    push!(out, (; lo, hi = prev, n))
    return out
end

"""
    _night_window(nights, i) -> (lo, hi)

The MJD window that keeps night `i` of `_nights`: halfway to the neighbouring nights (half a day
beyond the first and last), so an epoch at the edge of a night is never lost to rounding.
"""
function _night_window(ns, i)
    a = i == 1 ? ns[i].lo - 0.5 : (ns[i - 1].hi + ns[i].lo) / 2
    b = i == length(ns) ? ns[i].hi + 0.5 : (ns[i].hi + ns[i + 1].lo) / 2
    return (a, b)
end

"""
    filter_spec_kwargs(spec) -> (; read, filter, prep, drop) or a String

The Filter section's settings, as `key\tvalue` lines, turned into the four groups
`filter_dataset!` takes. A field left empty takes no part. A String is a refusal to show as it
is (a range whose minimum is above its maximum, a number that does not parse).

| keys | becomes |
|---|---|
| `wav_min`, `wav_max` (µm) | `wav_range` |
| `bl_min`, `bl_max` with `bl_unit` = `Mλ` or `m` | `baseline_range` or `baseline_length_range` |
| `nights` (nights kept, comma-separated: an MJD window `lo:hi`, or an MJD day) | `mjd_range`, one window per night |
| `insname`, `target`, `use_v2`, `use_t3`, `use_vis`, `use_flux`, `ignore_flags` | read options |
| `exclude_stations`, `exclude_baselines`, `exclude_triangles` (comma-separated) | the same keywords |
| `<obs>_add`, `<obs>_rel`, `<obs>_mult` for v2, t3amp, t3phi, visamp, visphi, flux | `oifits_prep` |
| `<obs>_cap`, `<obs>_caprel` | `max_<obs>_err`, `max_<obs>_err_rel` |
| `bad`, `snr`, `v2_min`, `v2_max`, `t3amp_min`, `t3amp_max` | the quality cuts |
| `remove_v2`, `remove_t3`, `remove_vis`, `remove_flux` (indices), `remove_uv` | the points picked by hand (`uv_bad` for uv points) |
| `drop` (comma-separated: `t3amp`, `t3phi`, `visamp`, `visphi`, …) | `filter_data`'s `drop` |
| `binning` (channels averaged into one; 1 = none), `binmode` (`block`, `pmoired`), `bincorr` | read options |
"""
function filter_spec_kwargs(spec::AbstractString)
    d = Dict{String,String}()
    for line in split(String(spec), '\n')
        f = split(line, '\t'; limit = 2)
        length(f) == 2 && !isempty(strip(f[2])) && (d[String(strip(f[1]))] = String(strip(f[2])))
    end
    num(k) = (v = tryparse(Float64, get(d, k, "")); v)
    for (k, v) in d
        (endswith(k, "_min") || endswith(k, "_max") || endswith(k, "_add") ||
         endswith(k, "_rel") || endswith(k, "_mult") || endswith(k, "_cap") ||
         endswith(k, "_caprel") || k == "snr") && tryparse(Float64, v) === nothing &&
            return "! $(k): $(repr(v)) is not a number"
    end
    list(k) = String[strip(x) for x in split(get(d, k, ""), ',') if !isempty(strip(x))]
    flag(k) = get(d, k, "") == "1"
    read = Dict{Symbol,Any}(); filt = Dict{Symbol,Any}(); prep = Dict{Symbol,Any}()

    function range!(key, lo, hi, scale, what, none_lo)
        a, b = num(lo), num(hi)
        a === nothing && b === nothing && return nothing
        a !== nothing && b !== nothing && a > b && return "! the $what range is empty: $a is above $b"
        filt[key] = Float64[a === nothing ? none_lo : a * scale, b === nothing ? 1e99 : b * scale]
        return nothing
    end
    r = range!(:wav_range, "wav_min", "wav_max", 1e-6, "wavelength", -1.0)
    r === nothing || return r
    r = get(d, "bl_unit", "Mλ") == "m" ?
        range!(:baseline_length_range, "bl_min", "bl_max", 1.0, "baseline", 0.0) :
        range!(:baseline_range, "bl_min", "bl_max", 1e6, "baseline", 0.0)
    r === nothing || return r
    if haskey(d, "nights")
        w = map(list("nights")) do x
            ab = tryparse.(Float64, split(x, ':'))
            length(ab) == 2 && !any(isnothing, ab) && return Float64[ab[1], ab[2]]
            k = tryparse(Int, x)
            k === nothing ? nothing : Float64[k, k + 1]
        end
        any(isnothing, w) && return "! nights: not a list of MJD windows or days"
        isempty(w) && return "! no night left: tick at least one"
        filt[:mjd_range] = sort(w)
    end

    haskey(d, "insname") && (read[:insname] = d["insname"])
    haskey(d, "target")  && (read[:targetname] = d["target"])
    for t in ("v2", "t3", "vis", "flux")
        get(d, "use_$t", "1") == "0" && (read[Symbol("use_$t")] = false)
    end
    flag("ignore_flags") && (read[:ignore_flags] = true)
    if haskey(d, "binning")
        n = tryparse(Int, d["binning"])
        (n === nothing || n < 1) && return "! average every: $(repr(d["binning"])) is not a channel count"
        if n > 1
            read[:binning] = n
            m = get(d, "binmode", "block")
            m in ("block", "pmoired") || return "! averaging mode: $(repr(m))"
            m == "block" || (read[:binmode] = Symbol(m))
            c = num("bincorr")
            if c !== nothing && m == "block"
                c >= 1 || return "! c: $(c) is below 1"
                c == 1 || (read[:bincorr] = c)
            end
        end
    end

    for (k, kk) in (("exclude_stations", :exclude_stations), ("exclude_baselines", :exclude_baselines),
                    ("exclude_triangles", :exclude_triangles))
        l = list(k); isempty(l) || (filt[kk] = l)
    end

    # Points picked by hand with the selection tools, as indices.
    for (k, kk) in (("remove_v2", :remove_v2), ("remove_t3", :remove_t3), ("remove_vis", :remove_vis),
                    ("remove_flux", :remove_flux), ("remove_uv", :uv_bad))
        haskey(d, k) || continue
        v = [tryparse(Int, x) for x in list(k)]
        any(isnothing, v) && return "! $(k): not a list of point numbers"
        isempty(v) || (filt[kk] = Int.(v))
    end

    for o in _ERR_ROWS
        a, m = num("$(o)_add"), num("$(o)_mult")
        a === nothing || (prep[Symbol("min_$(o)_err_add")] = a)
        m === nothing || (prep[Symbol("$(o)_err_mult")] = m)
        if !_is_phase(o)
            rr = num("$(o)_rel"); rr === nothing || (prep[Symbol("min_$(o)_err_rel")] = rr)
            cr = num("$(o)_caprel"); cr === nothing || (filt[Symbol("max_$(o)_err_rel")] = cr)
        end
        c = num("$(o)_cap"); c === nothing || (filt[Symbol("max_$(o)_err")] = c)
    end

    if flag("bad")
        filt[:filter_bad_data] = true
        for (k, kk) in (("snr", :filter_v2_snr_threshold), ("v2_min", :cutoff_minv2),
                        ("v2_max", :cutoff_maxv2), ("t3amp_min", :cutoff_mint3amp),
                        ("t3amp_max", :cutoff_maxt3amp))
            v = num(k); v === nothing || (filt[kk] = v)
        end
    end
    drop = Symbol[]
    if haskey(d, "drop")
        drop = Symbol.(filter(!isempty, strip.(split(d["drop"], ","))))
        bad = setdiff(drop, OITOOLS._DROPPABLE)
        isempty(bad) || return "! drop: unknown observable $(join(bad, ", "))"
    end
    return (; read, filter = filt, prep, drop)
end

"""
    dataset_summary(entry::DatasetEntry) -> Vector{Pair{String,String}}

What Explore's Dataset section shows, as `label => value` rows: target(s), instrument(s),
array and the stations the observables use, the MJD range and the nights it spans, the channels
and their wavelength range, the valid points per observable, AMPTYP / PHITYP, the binning, and
what a filtered dataset was filtered from.

Counted over every bin of the dataset, so a binned dataset reports what it holds, not its first
bin. The instruments are those the points come from (`OIdata.tables`); target and array names
are not kept in `OIdata`, so they are read from the file's own tables (a filtered or re-read
entry keeps its path).
"""
function dataset_summary(e::DatasetEntry)
    data = vec(e.data)
    d1 = first(data)
    rows = Pair{String,String}[]
    files = isempty(e.files) ? [e.path] : e.files
    length(files) > 1 && push!(rows, "files" => string(length(files), ": ",
                                                       _first_few(basename.(files))))
    hdrs = [h for h in (try OITOOLS.load_oifits(f) catch; nothing end for f in files)
            if h !== nothing]
    if !isempty(hdrs)
        tn = get(e.readopts, :targetname, "")
        targets = isempty(tn) ? unique(String[t.target for h in hdrs for t in h.target]) : [String(tn)]
        push!(rows, "target" => _first_few(targets))
        # the instruments the points come from (a dataset read with one `insname` holds one),
        # else the file's
        ins = unique(String[d.tables[t].insname for d in data
                            for t in vcat(d.v2_tab, d.t3_tab, d.vis_tab, d.flux_tab)
                            if 0 < t <= length(d.tables)])
        isempty(ins) && (ins = unique(String[i.insname for h in hdrs for i in h.instr]))
        push!(rows, "instrument" => _first_few(ins))
        push!(rows, "array" => _first_few(unique(String[a.arrname for h in hdrs for a in h.array])))
    end
    used = Int[]
    for d in data
        append!(used, vec(d.v2_sta_index), vec(d.t3_sta_index), vec(d.vis_sta_index))
    end
    used = unique(used)
    stations = [d1.sta_name[k] for k in eachindex(d1.sta_name) if d1.sta_index[k] in used]
    push!(rows, "stations" => string(length(stations), ": ", join(stations, " ")))
    mjd = reduce(vcat, (Float64.(d.uv_mjd) for d in data))
    if !isempty(mjd)
        lo, hi = extrema(mjd)
        nights = length(_nights(mjd))
        push!(rows, "MJD" => @sprintf("%.4f – %.4f  (%d night%s)", lo, hi, nights,
                                      nights == 1 ? "" : "s"))
    end
    lam = reduce(vcat, (Float64.(d.uv_lam) for d in data))
    if !isempty(lam)
        nch = length(unique(round.(lam; sigdigits = 10)))
        lo, hi = extrema(lam) .* 1e6
        push!(rows, "channels" => nch == 1 ? @sprintf("1 at %.4g µm", lo) :
                                   @sprintf("%d, %.4g – %.4g µm", nch, lo, hi))
    end
    tot(f) = sum(getfield(d, f) for d in data)
    counts = [("V²", :nv2), ("T3amp", :nt3amp), ("T3φ", :nt3phi), ("VISamp", :nvisamp),
              ("VISφ", :nvisphi), ("flux", :nflux)]
    push!(rows, "points" => join((string(l, " ", tot(f)) for (l, f) in counts if tot(f) > 0), " · "))
    # AMPTYP / PHITYP describe the OI_VIS tables, so only when there are some.
    if any(d -> length(d.visamp) > 0, data)
        push!(rows, "AMPTYP / PHITYP" => string(isempty(d1.amptyp) ? "—" : d1.amptyp, " / ",
                                                isempty(d1.phityp) ? "—" : d1.phityp))
    end
    size(e.data) == (1, 1) ||
        push!(rows, "bins" => string(size(e.data, 1), " wavelength × ", size(e.data, 2), " epoch"))
    e.derived_from === nothing || push!(rows, "filtered from" => e.derived_from)
    return rows
end

# A merged file can carry a dozen instrument tables; the panel names three and counts the rest.
_first_few(v::AbstractVector{<:AbstractString}; n::Int = 3) =
    length(v) <= n ? join(v, ", ") : join(v[1:n], ", ") * " + $(length(v) - n) more"
