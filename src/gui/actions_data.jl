# Actions on data.
#
# Every function here follows the same contract: change the session, and append the exact
# OITOOLS call that would have done the same thing. The contract is enforced by the
# round-trip test, which exports the log, runs it, and compares against the session.

# readoifits' own defaults, so exported scripts show only what the user actually changed.
const _READOIFITS_DEFAULTS = Dict{Symbol,Any}(
    :targetname => "", :polychromatic => false, :merge_oi_wavelength => false,
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

    kwstr = _kwargs([k => v for (k, v) in sort(collect(kw), by = first)];
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
Only the `[1,1]` bin is filtered — `filter_data` has no method for the 2-D array.
"""
function filter_dataset!(s::Session, name::AbstractString; kwargs...)
    src = dataset(s, name)
    kw  = Dict{Symbol,Any}(kwargs)
    d0  = src.data[1, 1]
    idx = OITOOLS.set_data_filter(d0; kwargs...)
    out = OITOOLS.filter_data(d0, idx)

    newname = _next_name(s.datasets, "data")
    entry = DatasetEntry(newname, src.path, reshape([out], 1, 1), kw, src.spectralbin,
                         src.temporalbin, src.name)
    push!(s.datasets, entry)

    kwstr = _kwargs([k => v for (k, v) in sort(collect(kw), by = first)])
    call  = isempty(kwstr) ? "set_data_filter($name[1,1])" :
                             "set_data_filter($name[1,1]; $kwstr)"
    log!(s, "$(newname) = reshape([filter_data($name[1,1], $call)], 1, 1)";
         note = "filter $(src.name)", binding = newname)
    return entry
end
