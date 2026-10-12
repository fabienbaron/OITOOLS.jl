# The Tools perspective: operations on files and datasets that are none of the other four —
# the header inspector, the validator and its repaired copy, channel averaging, and the splits
# (PLAN_Tools.md). A split or an average is a RE-READ of the dataset's file with one more
# `readoifits` keyword (`targetname`, `insname`, `temporalbin`, `binning`), so each part is an
# ordinary dataset: the Filter section works on it, and the exported script replays it.

# ── session actions ───────────────────────────────────────────────────────────────────

"""
    reread_dataset!(session, name, extra; label, binstats = nothing) -> DatasetEntry

Read dataset `name`'s file(s) again with its own read options plus `extra`, as a new dataset
labelled `label`, and log the `readoifits` call. A filtered dataset is refused: the re-read
would silently drop its cut. `binstats` (a Dict) receives the channel-averaging scatter check;
it is not a read option and is not kept.
"""
function reread_dataset!(s::Session, name::AbstractString, extra::AbstractDict;
                         label::AbstractString = "", binstats = nothing)
    src = dataset(s, name)
    src.derived_from === nothing ||
        throw(ArgumentError("$(src.name) is a filtered copy of $(src.derived_from); work on " *
                            "$(src.derived_from), then filter the result"))
    ropts = merge(src.readopts, Dict{Symbol,Any}(extra))
    data = OITOOLS.readoifits(sources(src); ropts..., warn = false, verbose = false, binstats)
    data isa AbstractArray || throw(ArgumentError("readoifits returned nothing for $(source_label(src))"))
    newname = _next_name(s.datasets, "data")
    entry = DatasetEntry(newname, src.path, data, ropts, get(ropts, :spectralbin, nothing),
                         get(ropts, :temporalbin, nothing), nothing, "", src.files,
                         isempty(src.label) ? String(label) : string(src.label, " · ", label))
    push!(s.datasets, entry)
    kwstr = _kwargs(Pair{Symbol,Any}[k => v for (k, v) in sort(collect(ropts), by = first)];
                    defaults = _READOIFITS_DEFAULTS)
    call = isempty(kwstr) ? "readoifits($(_literal(sources(src))))" :
                            "readoifits($(_literal(sources(src))); $kwstr)"
    log!(s, "$newname = $call"; note = "$(label) of $(src.name)", binding = newname)
    return entry
end

"""
    split_parts(entry, by; gap = 0.5, binwidth = 1.0) -> Vector{(label, extra)}

How dataset `entry` splits: one `(label, readoifits keywords)` per part. `by` is `:target`,
`:instrument`, `:gap` (a new epoch wherever the times are more than `gap` days apart) or
`:width` (fixed windows of `binwidth` days). Times come from every observable of every bin, and
each epoch's window is its own first and last time, so a re-read with that `temporalbin` keeps
exactly its points.
"""
function split_parts(e::DatasetEntry, by::Symbol; gap::Real = 0.5, binwidth::Real = 1.0)
    file = sources(e)
    file isa AbstractString ||
        throw(ArgumentError("$(e.name) is several files read as one; split each file instead"))
    if by === :target
        return [(String(n), Dict{Symbol,Any}(:targetname => String(n)))
                for n in OITOOLS.list_oifits_targets(file)]
    elseif by === :instrument
        names = String[]
        for d in e.data, t in d.tables
            isempty(t.insname) || t.insname in names || push!(names, t.insname)
        end
        return [(n, Dict{Symbol,Any}(:insname => n)) for n in names]
    elseif by === :gap || by === :width
        t = sort!(reduce(vcat, [vcat(d.v2_mjd, d.t3_mjd, d.vis_mjd, d.flux_mjd) for d in e.data]))
        isempty(t) && throw(ArgumentError("$(e.name) holds no point"))
        groups = if by === :gap
            cut = findall(>(gap), diff(t))
            [t[a:b] for (a, b) in zip([1; cut .+ 1], [cut; length(t)])]
        else
            binwidth > 0 || throw(ArgumentError("the window width must be positive"))
            k = floor.(Int, (t .- t[1]) ./ binwidth)
            [t[k .== j] for j in unique(k)]
        end
        day(m) = string(Date(1858, 11, 17) + Day(floor(Int, m)))
        labels = [day(g[1]) for g in groups]
        for (i, l) in enumerate(labels)      # two epochs on one UT day: tell them apart
            count(==(l), labels) > 1 && (labels[i] = string(l, " ", lpad(round(Int, 24 * (groups[i][1] % 1)), 2, '0'), "h"))
        end
        return [(labels[i], Dict{Symbol,Any}(:temporalbin => [[g[1], g[end]]]))
                for (i, g) in enumerate(groups)]
    end
    throw(ArgumentError("unknown split: $(by)"))
end

"""
    split_dataset!(session, name, by; gap, binwidth) -> Vector{DatasetEntry}

Split dataset `name` into one new dataset per part ([`split_parts`](@ref)), each a logged
re-read.
"""
split_dataset!(s::Session, name::AbstractString, by::Symbol; kw...) =
    [reread_dataset!(s, name, extra; label) for (label, extra) in split_parts(dataset(s, name), by; kw...)]

# ── shell calls ───────────────────────────────────────────────────────────────────────

"""The current dataset's file, for the Files section's "use the current dataset"; "" with none,
or when it is several files read as one."""
function shell_tools_source()
    e = current_dataset(_shell())
    (e === nothing || !isempty(e.files)) && return ""
    return e.path
end

"""
    shell_tools_headers(path) -> String

The file's HDUs, one per line: `index, kind, EXTNAME, INSNAME, OI_REVN, rows, columns-or-dims`
(tab-separated). An error string starting with `!` when it cannot be read.
"""
function shell_tools_headers(path::AbstractString)
    p = strip_file_url(String(path))
    console!(_shell(), "oifits_headers($(_literal(p)))"; kind = :cmd)
    h = OITOOLS.oifits_headers(p)
    return join([join((x.index, x.kind, x.extname, x.insname, x.revn, x.nrows,
                       x.kind === :table ? string(length(x.columns), " columns") :
                       isempty(x.dims) ? "header only" : join(x.dims, "×")), "\t") for x in h], "\n")
end

"""
    shell_tools_hdu(path, i) -> String

HDU `i` of the file as text: its columns (name, format, unit, dimensions) for a table, then
every header card in order, `NAME = value / comment`.
"""
function shell_tools_hdu(path::AbstractString, i::Real)     # QML sends numbers as Float64
    i = round(Int, i)
    h = OITOOLS.oifits_headers(strip_file_url(String(path)))
    1 <= i <= length(h) || return "! no HDU $(i)"
    x = h[i]
    io = IOBuffer()
    if x.kind === :table
        println(io, "Columns ($(x.nrows) rows):")
        for c in x.columns
            println(io, "  ", rpad(c.name, 16), " ", rpad(c.format, 8), " ",
                    rpad(isempty(c.unit) ? "" : "[$(c.unit)]", 10), isempty(c.dims) ? "" : "dims $(c.dims)")
        end
        println(io)
    end
    println(io, "Header ($(length(x.cards)) cards):")
    for c in x.cards
        if c.name in ("COMMENT", "HISTORY", "")
            println(io, rpad(c.name, 8), " ", c.value)
        else
            v = c.value isa AbstractString ? "'" * c.value * "'" : string(c.value)
            println(io, rpad(c.name, 8), " = ", v, isempty(c.comment) ? "" : " / " * c.comment)
        end
    end
    return String(take!(io))
end

"""
    shell_tools_validate(path) -> String

`validate_oifits` as lines `severity, table, message` (tab-separated); "" when nothing is found.
"""
function shell_tools_validate(path::AbstractString)
    sh = _shell()
    p = strip_file_url(String(path))
    console!(sh, "validate_oifits($(_literal(p)))"; kind = :cmd)
    r = OITOOLS.validate_oifits(p)
    ne = count(x -> x.severity === :error, r); nw = count(x -> x.severity === :warning, r)
    console!(sh, "$(ne) error$(ne == 1 ? "" : "s"), $(nw) warning$(nw == 1 ? "" : "s"), " *
                 "$(length(r) - ne - nw) note$(length(r) - ne - nw == 1 ? "" : "s")")
    return join([join((x.severity, x.table, x.message), "\t") for x in r], "\n")
end

"""The file's targets, one per line (the repaired copy is written one target at a time); ""
for a path that is not (yet) a readable OIFITS file -- the field asks on every keystroke."""
function shell_tools_targets(path::AbstractString)
    p = strip_file_url(String(path))
    (isfile(p) && !OITOOLS._is_zip(p)) || return ""
    try
        return join(String.(OITOOLS.list_oifits_targets(p)), "\n")
    catch
        return ""
    end
end

"""
    shell_tools_repair(path, out, target) -> String

Write a repaired copy of `path` to `out` (`repair_oifits`); `target` is "" for a one-target
file. Returns a status line, starting with `!` on failure.
"""
function shell_tools_repair(path::AbstractString, out::AbstractString, target::AbstractString = "")
    sh = _shell()
    p = strip_file_url(String(path)); o = strip_file_url(String(out))
    tk = isempty(target) ? "" : "; targetname = $(_literal(String(target)))"
    console!(sh, "repair_oifits($(_literal(p)), $(_literal(o))$(tk))"; kind = :cmd)
    try
        OITOOLS.repair_oifits(p, o; targetname = String(target), overwrite = true)
    catch err
        return console!(sh, "! could not repair $(basename(p)): " * _cause(err); kind = :err)
    end
    return console!(sh, "wrote the repaired copy " * basename(o))
end

"""
    shell_tools_bin(n, mode, c) -> String

Average every `n` channels of the current dataset into a new dataset (`readoifits(…; binning =
n, binmode, bincorr = c)`), `mode` `"block"` or `"pmoired"`. Averaging is not grouping: the
channels become fewer, each with the average of the ones it replaces. Returns a status line with
the scatter check (χ²/dof about a line in λ, per observable) for block averaging.
"""
function shell_tools_bin(n::Real, mode::AbstractString, c::Real = 1)
    n = round(Int, n)
    sh = _shell()
    e = current_dataset(sh)
    e === nothing && return "! no dataset loaded"
    n >= 2 || return "! average at least 2 channels"
    m = Symbol(mode)
    m in (:block, :pmoired) || return "! unknown averaging mode: " * mode
    extra = Dict{Symbol,Any}(:binning => Int(n))
    m === :block || (extra[:binmode] = m)
    c == 1 || (extra[:bincorr] = Float64(c))
    nlog = length(sh.session.log)
    stats = Dict{Symbol,Float64}()
    entry = try
        reread_dataset!(sh.session, e.name, extra; label = "binned ×$(n)", binstats = stats)
    catch err
        return console!(sh, "! could not average $(e.name): " * _cause(err); kind = :err)
    end
    for l in sh.session.log[nlog+1:end]
        console!(sh, l.code; kind = :cmd)
    end
    sh.current = length(sh.session.datasets)
    sh.status = try refresh_plot!(sh) catch err; "averaged into $(entry.name), but cannot plot it: " * _cause(err) end
    msg = "averaged ×$(n) into $(entry.name)" *
          (isempty(stats) ? "" : "; scatter χ²/dof " *
           join(["$(k) $(round(v; sigdigits = 3))" for (k, v) in sort(collect(stats))], ", "))
    return console!(sh, msg)
end

_split_by(by::AbstractString) =
    by == "gap" ? :gap : by == "width" ? :width : by == "target" ? :target :
    by == "instrument" ? :instrument : throw(ArgumentError("unknown split: " * by))

"""
    shell_tools_split_preview(by, param) -> String

The parts a split of the current dataset would make, one per line: `label, V² / T3 / VIS points`
(tab-separated), counted in the dataset as it is; a target's points are only known once read,
so a split by target shows "—". `by` is `gap`, `width` (both take `param` in days), `target` or
`instrument`.
"""
function shell_tools_split_preview(by::AbstractString, param::Real = 0.5)
    e = current_dataset(_shell())
    e === nothing && return "! no dataset loaded"
    b = _split_by(by)
    parts = split_parts(e, b; gap = param, binwidth = param)
    counts(x) = if haskey(x, :temporalbin)
        lo, hi = x[:temporalbin][1]
        inw(v) = count(m -> lo <= m <= hi, v)
        (sum(d -> inw(d.v2_mjd), e.data), sum(d -> inw(d.t3_mjd), e.data), sum(d -> inw(d.vis_mjd), e.data))
    elseif haskey(x, :insname)
        n = x[:insname]
        inst(d, tab) = count(t -> t != 0 && d.tables[t].insname == n, tab)
        (sum(d -> inst(d, d.v2_tab), e.data), sum(d -> inst(d, d.t3_tab), e.data), sum(d -> inst(d, d.vis_tab), e.data))
    else
        nothing
    end
    line(l, x) = (c = counts(x); c === nothing ? string(l, "\t—") :
                  string(l, "\t", c[1], " V², ", c[2], " T3, ", c[3], " VIS"))
    return join([line(l, x) for (l, x) in parts], "\n")
end

"""
    shell_tools_split(by, param, folder) -> String

Split the current dataset into new datasets; with `folder` not empty, also write each part to
`<folder>/<file>_<part>.oifits` (an existing file is not replaced, nor ever the source).
"""
function shell_tools_split(by::AbstractString, param::Real = 0.5, folder::AbstractString = "")
    sh = _shell()
    e = current_dataset(sh)
    e === nothing && return "! no dataset loaded"
    nlog = length(sh.session.log)
    entries = try
        split_dataset!(sh.session, e.name, _split_by(by); gap = param, binwidth = param)
    catch err
        return console!(sh, "! could not split $(e.name): " * _cause(err); kind = :err)
    end
    for l in sh.session.log[nlog+1:end]
        console!(sh, l.code; kind = :cmd)
    end
    written = 0
    dir = strip_file_url(String(folder))
    if !isempty(dir)
        stem = replace(basename(e.path), r"\.(oi)?fits$"i => "")
        for x in entries
            out = joinpath(dir, string(stem, "_", replace(x.label, r"[^A-Za-z0-9._+-]+" => "_"), ".oifits"))
            OITOOLS._same_file(out, e.path) && continue
            console!(sh, "write_oifits($(_literal(out)), $(x.name))"; kind = :cmd)
            try
                OITOOLS.write_oifits(out, x.data)      # the console's record only, as Export's
                written += 1
            catch err
                console!(sh, "! not written: " * _cause(err); kind = :err)
            end
        end
    end
    sh.current = length(sh.session.datasets)
    sh.status = try refresh_plot!(sh) catch err; "split, but cannot plot: " * _cause(err) end
    return console!(sh, "split $(e.name) into $(length(entries)) datasets" *
                        (isempty(dir) ? "" : ", $(written) written to $(dir)"))
end
