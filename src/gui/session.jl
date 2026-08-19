# The session: one object holding everything the four perspectives share.
#
# This is the whole reason OITOOLSGUI is one application rather than four. JMMC ships Aspro2,
# OIFITSExplorer, LITpro and OIImaging as separate programs, so moving a dataset between them
# means exporting and reimporting files. Here a dataset is loaded once and every perspective
# refers to the same entry.
#
# Provenance is recorded, not inferred: an ImageEntry knows which dataset, weights,
# regularizers and starting image produced it, so the Image perspective can say where a
# figure came from and the command log can replay it.

"""
    DatasetEntry

One loaded OIFITS file.

`data` is what `readoifits` returned — an `Array{OIdata,2}` of `(nwavbin, ntimebin)`, kept in
that shape deliberately. `readoifits` does not record the bin *definitions* anywhere, so
`spectralbin`/`temporalbin` keep what was asked for; without them a session cannot describe
its own binning.

`name` is the binding used in the exported script (`data1`, `data2`, …).
"""
mutable struct DatasetEntry
    name        :: String
    path        :: String
    data        :: Any                     # Array{OIdata,2} as returned by readoifits
    readopts    :: Dict{Symbol,Any}
    spectralbin :: Any
    temporalbin :: Any
    derived_from:: Union{Nothing,String}   # set when this came from filtering another entry
end

"""
    ModelEntry

A parametric model plus whatever has been fitted with it. `dict` and `free` are the OITOOLS
model dictionary and free-parameter list; `results` accumulates fit results so competing fits
can be compared side by side (§3.9 of the design).
"""
mutable struct ModelEntry
    name    :: String
    dict    :: Dict{String,Any}
    free    :: Vector{String}
    lb      :: Dict{String,Float64}
    ub      :: Dict{String,Float64}
    results :: Vector{Any}
end

"""
    ImageEntry

A reconstruction, with the provenance needed to reproduce it.

`image` is either a single array or, for the stochastic engines (SQUEEZE annealing/tempering,
OIVI), an ensemble — `posterior` then holds the per-pixel spread and the samples. Those
engines return a distribution, not a point estimate, and the display branches on it.
"""
mutable struct ImageEntry
    name      :: String
    image     :: Any
    pixsize   :: Float64
    engine    :: Symbol                     # :vmlmb, :bsmem, :bsdmm, :sparco, :squeeze, ...
    dataset   :: String                     # DatasetEntry.name it was reconstructed from
    weights   :: Vector{Float64}
    regs      :: Any
    posterior :: Any                        # nothing, or (; mean, sigma, samples)
end

"""
    Selection

Points brushed in one perspective, shared with the others. `uv` indexes into `OIdata.uv`,
which is what `set_data_filter(...; uv_bad = …)` accepts — so a lasso in the Explore
perspective feeds straight back into the data path.
"""
mutable struct Selection
    dataset :: String
    uv      :: Vector{Int}
    v2      :: Vector{Int}
    t3      :: Vector{Int}
end
Selection() = Selection("", Int[], Int[], Int[])

"""
    Session()

The application state. Perspectives are views onto this; switching tabs changes nothing here.
"""
mutable struct Session
    datasets  :: Vector{DatasetEntry}
    models    :: Vector{ModelEntry}
    images    :: Vector{ImageEntry}
    configs   :: Any                        # the NamedTuple from OITOOLS.list_configs()
    log       :: Vector{LogEntry}
    selection :: Selection
end

function Session()
    cfg = try
        OITOOLS.list_configs()
    catch err                               # a broken configs dir must not stop the app
        @warn "could not enumerate OITOOLS configs" err
        (; facilities = String[], combiners = String[], wavelengths = String[],
           targets = String[], unknown = String[], wave_combiner = Dict{String,String}())
    end
    Session(DatasetEntry[], ModelEntry[], ImageEntry[], cfg, LogEntry[], Selection())
end

function Base.show(io::IO, s::Session)
    print(io, "Session(", length(s.datasets), " dataset(s), ", length(s.models),
              " model(s), ", length(s.images), " image(s), ", length(s.log), " log step(s))")
end

"""
    dataset(session, name) -> DatasetEntry

Look an entry up by the name used in the exported script. Throws with the available names
rather than returning `nothing`, because every caller would have to check anyway.
"""
function dataset(s::Session, name::AbstractString)
    i = findfirst(d -> d.name == name, s.datasets)
    i === nothing && throw(KeyError("no dataset $(repr(name)); have " *
                                    string([d.name for d in s.datasets])))
    return s.datasets[i]
end

_next_name(existing, stem) = string(stem, length(existing) + 1)
