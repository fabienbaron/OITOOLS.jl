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

Load an OIFITS file into the session and log the `readoifits` call that did it.

`kwargs` are passed straight through to `readoifits`, so anything it accepts works here
(`targetname`, `polychromatic`, `merge_oi_wavelength`, `filter_bad_data`, `T`, …).

Two behaviours worth knowing, both inherited from `readoifits`:

  * it returns an `Array{OIdata,2}` of `(nwavbin, ntimebin)`, kept in that shape here;
  * on a missing file it returns `[[]]` rather than throwing, which would surface as a
    baffling error much later — so this checks first and throws immediately.
"""
function load_dataset!(s::Session, path::AbstractString; kwargs...)
    isfile(path) || throw(ArgumentError("no such OIFITS file: $(repr(path))"))
    kw = Dict{Symbol,Any}(kwargs)
    data = OITOOLS.readoifits(path; kwargs...)
    data isa AbstractArray || throw(ArgumentError("readoifits did not return data for $(repr(path))"))

    name = _next_name(s.datasets, "data")
    entry = DatasetEntry(name, String(path), data, kw,
                         get(kw, :spectralbin, nothing), get(kw, :temporalbin, nothing), nothing)
    push!(s.datasets, entry)

    # `Pair{Symbol,Any}[...]`, not a bare comprehension: over an EMPTY `kw` the comprehension
    # infers `Vector{Any}`, which `_kwargs(::Vector{<:Pair})` does not accept — so loading a
    # file with no read options at all threw a MethodError rather than logging the plain call.
    kwstr = _kwargs(Pair{Symbol,Any}[k => v for (k, v) in sort(collect(kw), by = first)];
                    defaults = _READOIFITS_DEFAULTS)
    call  = isempty(kwstr) ? "readoifits($(_literal(path)))" :
                             "readoifits($(_literal(path)); $kwstr)"
    log!(s, "$name = $call"; note = "load $(basename(path))", binding = name)
    return entry
end

"""
    filter_dataset!(session, name; kwargs...) -> DatasetEntry

Apply `set_data_filter` + `filter_data` to a loaded dataset and add the result as a new
entry, leaving the original untouched.

Filtering in OITOOLS is destructive-by-copy, so keeping the source is what gives the GUI undo
for free: the filtered dataset is derived, and dropping it restores the original view.
Every bin is filtered, each by its own discard lists, and the dataset keeps its shape: a
binned dataset stays binned.
"""
function filter_dataset!(s::Session, name::AbstractString; kwargs...)
    src = dataset(s, name)
    kw  = Dict{Symbol,Any}(kwargs)
    out = map(d -> OITOOLS.filter_data(d, OITOOLS.set_data_filter(d; kwargs...)), src.data)

    newname = _next_name(s.datasets, "data")
    # The source's READ options: the filter is in the log, and a re-read of the file (rebinning)
    # must not be handed filter keywords as though readoifits took them.
    entry = DatasetEntry(newname, src.path, out, copy(src.readopts), src.spectralbin,
                         src.temporalbin, src.name)
    push!(s.datasets, entry)

    kwstr = _kwargs([k => v for (k, v) in sort(collect(kw), by = first)])
    call  = isempty(kwstr) ? "set_data_filter(d)" : "set_data_filter(d; $kwstr)"
    log!(s, "$(newname) = map(d -> filter_data(d, $call), $name)";
         note = "filter $(src.name)", binding = newname)
    return entry
end

"""
    dataset_summary(entry::DatasetEntry) -> Vector{Pair{String,String}}

What Explore's Dataset section shows, as `label => value` rows: target(s), instrument(s),
array and the stations the observables use, the MJD range and the nights it spans, the channels
and their wavelength range, the valid points per observable, AMPTYP / PHITYP, the binning, and
what a filtered dataset was filtered from.

Counted over every bin of the dataset, so a binned dataset reports what it holds, not its first
bin. Target, instrument and array names are not kept in `OIdata`, so they are read from the
file's own tables (a filtered or re-read entry keeps its path).
"""
function dataset_summary(e::DatasetEntry)
    data = vec(e.data)
    d1 = first(data)
    rows = Pair{String,String}[]
    hdr = try
        OITOOLS.load_oifits(e.path)
    catch
        nothing
    end
    if hdr !== nothing
        tn = get(e.readopts, :targetname, "")
        targets = isempty(tn) ? unique(String[t.target for t in hdr.target]) : [String(tn)]
        push!(rows, "target" => _first_few(targets))
        push!(rows, "instrument" => _first_few(unique(String[i.insname for i in hdr.instr])))
        push!(rows, "array" => _first_few(unique(String[a.arrname for a in hdr.array])))
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
        nights = length(unique(floor.(Int, mjd)))
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
    push!(rows, "AMPTYP / PHITYP" => string(isempty(d1.amptyp) ? "—" : d1.amptyp, " / ",
                                            isempty(d1.phityp) ? "—" : d1.phityp))
    size(e.data) == (1, 1) ||
        push!(rows, "bins" => string(size(e.data, 1), " wavelength × ", size(e.data, 2), " epoch"))
    e.derived_from === nothing || push!(rows, "filtered from" => e.derived_from)
    return rows
end

# A merged file can carry a dozen instrument tables; the panel names three and counts the rest.
_first_few(v::AbstractVector{<:AbstractString}; n::Int = 3) =
    length(v) <= n ? join(v, ", ") : join(v[1:n], ", ") * " + $(length(v) - n) more"
