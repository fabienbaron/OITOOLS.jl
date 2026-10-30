# The Tools perspective: operations on files and datasets that are none of the other four —
# the header inspector, the validator, channel binning, the splits and merging, the cube viewer.
# OIFITS handling works from a FILE: binning or splitting reads it with one more `readoifits`
# keyword (`binning`, `spectralbin`, `temporalbin`, `insname`, `targetname`), so each result is
# an ordinary dataset — the Filter section works on it, and the exported script replays it.

# ── session actions ───────────────────────────────────────────────────────────────────

"""
    read_file_dataset!(session, path, extra; label, binstats = nothing) -> DatasetEntry

Read an OIFITS file with the read options `extra` into a new dataset labelled `label`, and log
the `readoifits` call. `binstats` (a Dict) receives the channel-averaging scatter check; it is not
a read option and is not kept.
"""
function read_file_dataset!(s::Session, path::AbstractString, extra::AbstractDict;
                            label::AbstractString = "", binstats = nothing)
    p = String(path)
    isfile(p) || throw(ArgumentError("no such file: $(repr(p))"))
    ropts = Dict{Symbol,Any}(extra)
    data = OITOOLS.readoifits(p; ropts..., warn = false, verbose = false, binstats)
    data isa AbstractArray || throw(ArgumentError("readoifits returned nothing for $(basename(p))"))
    newname = _next_name(s.datasets, "data")
    entry = DatasetEntry(; name = newname, path = p, data, readopts = ropts,
                         spectralbin = get(ropts, :spectralbin, nothing),
                         temporalbin = get(ropts, :temporalbin, nothing), label = String(label))
    push!(s.datasets, entry)
    kwstr = _kwargs(Pair{Symbol,Any}[k => v for (k, v) in sort(collect(ropts), by = first)];
                    defaults = _READOIFITS_DEFAULTS)
    call = isempty(kwstr) ? "readoifits($(_literal(p)))" : "readoifits($(_literal(p)); $kwstr)"
    log!(s, "$newname = $call"; note = "$(label) of $(basename(p))", binding = newname)
    return entry
end

# The file read with the default options, as split and preview work out their parts from it;
# kept while the file is unchanged.
const FILE_READ = Ref{Any}(nothing)
function _file_entry(path::AbstractString)
    p = String(path)
    isfile(p) || throw(ArgumentError("no such file: $(repr(p))"))
    c = FILE_READ[]
    (c !== nothing && c.path == p && c.mtime == mtime(p)) && return c.entry
    data = OITOOLS.readoifits(p; warn = false, verbose = false)
    entry = DatasetEntry(; name = "file", path = p, data)
    FILE_READ[] = (; path = p, mtime = mtime(p), entry)
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
    split_file!(session, path, by; gap, binwidth) -> Vector{DatasetEntry}

One new dataset per part of the file ([`split_parts`](@ref)), each a logged read.
"""
split_file!(s::Session, path::AbstractString, by::Symbol; kw...) =
    [read_file_dataset!(s, path, extra; label) for (label, extra) in split_parts(_file_entry(path), by; kw...)]

# ── shell calls ───────────────────────────────────────────────────────────────────────

"""The current dataset's file, for **Current**; "" with none, or when it is several files read as
one."""
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
                       (isempty(x.dims) || prod(x.dims) == 0) ? "header only" : join(x.dims, "×")),
                  "\t") for x in h], "\n")
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

# ── a table's contents ────────────────────────────────────────────────────────────────

# One cell of the Data view: its text, and whether its line is flagged (greyed).
struct TableCell
    text    :: String
    flagged :: Bool
end

"""
    OITableView <: AbstractMatrix{TableCell}

The Data view's model contents: the table `oifits_table` read, and `perm`, the lines shown (the
filter's survivors, in the sort's order). Cells are formatted when Qt asks for them, so a table
of hundreds of thousands of lines costs its columns once and nothing per line. `channel` picks
the channel an array column shows in the table as stored.
"""
mutable struct OITableView <: AbstractMatrix{TableCell}
    t       :: Any
    perm    :: Vector{Int}
    channel :: Int
    sortcol :: Int
    rev     :: Bool
    filter  :: String
    hide    :: Bool
end
OITableView() = OITableView(nothing, Int[], 1, 0, false, "", false)
Base.size(v::OITableView) = v.t === nothing ? (0, 1) : (length(v.perm), length(v.t.names))
Base.IndexStyle(::Type{OITableView}) = IndexCartesian()
function Base.getindex(v::OITableView, i::Int, j::Int)
    v.t === nothing && return TableCell("", false)
    r = v.perm[i]
    return TableCell(_cell_text(v.t.columns[j][r], v), _line_flagged(v, r))
end

# A line's flag: per channel, its own; in a table as stored, the shown channel's (a row with any
# channel flagged is not a flagged row).
function _line_flagged(v::OITableView, r)
    t = v.t
    (t.per_channel || t.nwave == 0) && return t.flagged[r]
    k = findfirst(==("FLAG"), t.names)
    k === nothing && return t.flagged[r]
    f = t.columns[k][r]
    return f isa AbstractVector ? Bool(f[clamp(v.channel, 1, length(f))]) : Bool(f)
end

_num_text(x::AbstractFloat) = isnan(x) ? "NaN" : @sprintf("%.6g", x)
_num_text(x::Bool) = x ? "T" : "F"
_num_text(x::Integer) = string(x)
_num_text(x) = strip(string(x))
function _cell_text(x, v)
    x isa AbstractVector || return String(_num_text(x))
    # an array column of the table as stored: the chosen channel, else first … last
    v.t.nwave > 0 && length(x) == v.t.nwave && return String(_num_text(x[clamp(v.channel, 1, length(x))]))
    length(x) <= 4 && return join(_num_text.(x), ", ")
    return string(_num_text(first(x)), " … ", _num_text(last(x)), " (", length(x), ")")
end

# The view, and the model QML shows it through; built by `gui()` before the window exists.
const OITABLE = Ref{Any}(nothing)
const OITABLE_MODEL = Ref{Any}(nothing)

"""
    build_oitable_model() -> JuliaItemModel

The Data view's model, over an empty [`OITableView`](@ref): roles `text` and `flagged`, the
column headers the names with their units below.
"""
function build_oitable_model()
    v = OITableView()
    m = QML.JuliaItemModel(v)
    QML.setheadergetter!(m, function (data, i, orientation, role)
        t = data.t
        (t === nothing || orientation != QML.Horizontal || !(1 <= i <= length(t.names))) && return string(i)
        isempty(t.units[i]) ? t.names[i] : t.names[i] * "\n" * t.units[i]
    end)
    OITABLE[] = v; OITABLE_MODEL[] = m
    return m
end

_table_refresh!() = OITABLE_MODEL[] === nothing || QML.force_model_update(OITABLE_MODEL[])

# The filter's terms: `key a–b` (λ, B, B/λ, MJD; an en dash, `..` or `-`) restricts a column to a
# range, any other word must appear in the stations (or, without them, in any text column).
function _filter_lines(t, text)
    lines = collect(1:t.nlines)
    toks = split(strip(text))
    keys = Dict("λ" => "λ", "lambda" => "λ", "wl" => "λ", "b" => "B", "b/λ" => "B/λ",
                "mjd" => "MJD")
    i = 1
    while i <= length(toks)
        tok = lowercase(toks[i])
        if haskey(keys, tok) && i < length(toks)
            m = match(r"^([-+]?[0-9.eE]+)\s*(?:–|\.\.|-)\s*([-+]?[0-9.eE]+)$", toks[i+1])
            k = findfirst(==(keys[tok]), t.names)
            if m !== nothing
                # a range on a column this table does not have (λ in a table as stored) applies
                # to nothing rather than matching no line
                if k !== nothing
                    lo, hi = parse(Float64, m[1]), parse(Float64, m[2])
                    c = t.columns[k]
                    filter!(r -> c[r] isa Real && lo <= c[r] <= hi, lines)
                end
                i += 2; continue
            end
        end
        word = toks[i]
        k = findfirst(==("stations"), t.names)
        cs = k === nothing ? [c for c in t.columns if eltype(c) <: AbstractString] : [t.columns[k]]
        filter!(r -> any(occursin(word, c[r]) for c in cs), lines)
        i += 1
    end
    return lines
end

function _table_apply!(v::OITableView)
    t = v.t
    t === nothing && (v.perm = Int[]; return v)
    lines = _filter_lines(t, v.filter)
    v.hide && filter!(r -> !_line_flagged(v, r), lines)
    if v.sortcol > 0
        c = t.columns[v.sortcol]
        key(r) = (x = c[r]; x isa AbstractVector ? (isempty(x) ? NaN : x[clamp(v.channel, 1, length(x))]) : x)
        if all(r -> key(r) isa Real, lines)
            # NaN last whichever the direction
            sort!(lines; by = r -> (isnan(Float64(key(r))), v.rev ? -Float64(key(r)) : Float64(key(r))))
        else
            sort!(lines; by = r -> string(key(r)), rev = v.rev)
        end
    end
    v.perm = lines
    return v
end

function _table_status(v)
    v.t === nothing && return ""
    nf = count(r -> _line_flagged(v, r), 1:v.t.nlines)
    ch = !v.t.per_channel && v.t.nwave > 1 ?
         @sprintf(" · channel %d at %.6g µm", v.channel, v.t.wl[clamp(v.channel, 1, v.t.nwave)]) : ""
    return string(length(v.perm), " of ", v.t.nlines, " lines", v.t.per_channel ? " (per channel)" : "",
                  nf > 0 ? ", $(nf) flagged" : "", ch)
end

"""
    shell_table_show(path, hdu, per_channel) -> String

Show table HDU `hdu` of `path` in the Data view ([`oifits_table`](@ref)): an OI data table per
channel unless `per_channel` is false. Returns `ok\tstatus\tnwave\tper_channel` or a message
beginning with `!`. The filter, the "hide flagged" choice and the channel carry over; the sort
is reset.
"""
function shell_table_show(path::AbstractString, hdu::Real, per_channel::Bool)
    v = OITABLE[]
    v === nothing && return "! no table view"
    p = strip_file_url(String(path))
    t = try
        OITOOLS.oifits_table(p, round(Int, hdu); per_channel)
    catch err
        return "! " * _cause(err)
    end
    v.t = t; v.sortcol = 0; v.rev = false
    _table_apply!(v); _table_refresh!()
    return join(("ok", _table_status(v), t.nwave, t.per_channel ? 1 : 0), "\t")
end

"""
    shell_table_sort(column) -> String

Sort the Data view by `column` (0-based, as QML numbers them); the same column again reverses
it. Returns the status line.
"""
function shell_table_sort(column::Real)
    v = OITABLE[]
    (v === nothing || v.t === nothing) && return ""
    c = round(Int, column) + 1
    v.rev = v.sortcol == c ? !v.rev : false
    v.sortcol = c
    _table_apply!(v); _table_refresh!()
    return _table_status(v)
end

"""
    shell_table_filter(text, hide_flagged) -> String

Filter the Data view: `text`'s words must appear in the stations, `λ 2.16–2.17`, `B 100–200`,
`B/λ a–b` or `mjd a–b` restrict a column to a range; `hide_flagged` drops flagged lines. Returns
the status line.
"""
function shell_table_filter(text::AbstractString, hide_flagged::Bool)
    v = OITABLE[]
    v === nothing && return ""
    v.filter = String(text); v.hide = hide_flagged
    _table_apply!(v); _table_refresh!()
    return _table_status(v)
end

"""
    shell_table_channel(k) -> String

The channel (1-based) an array column, and the flag, show in a table as stored. Returns the
status line.
"""
function shell_table_channel(k::Real)
    v = OITABLE[]
    (v === nothing || v.t === nothing) && return ""
    v.channel = clamp(round(Int, k), 1, max(v.t.nwave, 1))
    _table_apply!(v); _table_refresh!()
    return _table_status(v)
end

"""
    shell_table_widths() -> String

Each column's width in characters, `header\tcells` per line: the longer of its name and unit,
and the longest of its cells over up to 500 lines spread through the table (capped at 60), so a
column is as wide as what it holds and a table of any size costs the same.
"""
function shell_table_widths()
    v = OITABLE[]
    (v === nothing || v.t === nothing) && return ""
    t = v.t
    n = length(v.perm)
    sample = n <= 500 ? (1:n) : round.(Int, range(1, n; length = 500))
    out = String[]
    for j in eachindex(t.names)
        hdr = max(length(t.names[j]), length(t.units[j]))
        cells = isempty(sample) ? 0 : maximum(length(v[i, j].text) for i in sample)
        push!(out, string(hdr, "\t", min(cells, 60)))
    end
    return join(out, "\n")
end

"The comment of column `column` (0-based) of the Data view, for its tooltip."
function shell_table_comment(column::Real)
    v = OITABLE[]
    (v === nothing || v.t === nothing) && return ""
    c = round(Int, column) + 1
    1 <= c <= length(v.t.comments) || return ""
    return v.t.comments[c]
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

"""The file's targets, one per line (Export writes one target at a time); "" for a path that is
not (yet) a readable OIFITS file -- the field asks on every keystroke."""
function shell_tools_targets(path::AbstractString)
    p = strip_file_url(String(path))
    isfile(p) || return ""
    try
        return join(String.(OITOOLS.list_oifits_targets(p)), "\n")
    catch
        return ""
    end
end

"""
    shell_tools_export(path, out, target) -> String

Write the file `path` as clean OIFITS v2 to `out` (`repair_oifits`: every unflagged point kept
exactly, instruments under their names); `target` is "" for a one-target file.
"""
function shell_tools_export(path::AbstractString, out::AbstractString, target::AbstractString = "")
    sh = _shell()
    p = strip_file_url(String(path)); o = strip_file_url(String(out))
    tk = isempty(target) ? "" : "; targetname = $(_literal(String(target)))"
    console!(sh, "repair_oifits($(_literal(p)), $(_literal(o))$(tk))"; kind = :cmd)
    try
        OITOOLS.repair_oifits(p, o; targetname = String(target), overwrite = true)
    catch err
        return console!(sh, "! could not export $(basename(p)): " * _cause(err); kind = :err)
    end
    return console!(sh, "wrote " * basename(o))
end

# Select the newest dataset and draw it.
function _show_newest!(sh, what)
    sh.current = length(sh.session.datasets)
    sh.status = try refresh_plot!(sh) catch err; "$(what), but cannot plot it: " * _cause(err) end
end

# Run `f`, then echo to the console the session-log lines it added.
function _with_log(f, sh)
    nlog = length(sh.session.log)
    r = f()
    for l in sh.session.log[nlog+1:end]
        console!(sh, l.code; kind = :cmd)
    end
    return r
end

"""
    shell_tools_bin(path, n, mode, c) -> String

Average every `n` channels of the file into a new dataset (`readoifits(path; binning = n,
binmode, bincorr = c)`), `mode` `"block"` or `"pmoired"`. Averaging makes fewer channels, each
the average of the ones it replaces. Returns a status line with the scatter check (χ²/dof about a
line in λ, per observable) for block averaging.
"""
function shell_tools_bin(path::AbstractString, n::Real, mode::AbstractString, c::Real = 1)
    n = round(Int, n)
    sh = _shell()
    p = strip_file_url(String(path))
    isempty(p) && return "! choose an OIFITS file"
    n >= 2 || return "! average at least 2 channels"
    m = Symbol(mode)
    m in (:block, :pmoired) || return "! unknown averaging mode: " * mode
    extra = Dict{Symbol,Any}(:binning => n)
    m === :block || (extra[:binmode] = m)
    c == 1 || (extra[:bincorr] = Float64(c))
    stats = Dict{Symbol,Float64}()
    entry = try
        _with_log(sh) do
            read_file_dataset!(sh.session, p, extra; label = "binned ×$(n)", binstats = stats)
        end
    catch err
        return console!(sh, "! could not average $(basename(p)): " * _cause(err); kind = :err)
    end
    _show_newest!(sh, "averaged into $(entry.name)")
    msg = "averaged ×$(n) into $(entry.name)" *
          (isempty(stats) ? "" : "; scatter χ²/dof " *
           join(["$(k) $(round(v; sigdigits = 3))" for (k, v) in sort(collect(stats))], ", "))
    return console!(sh, msg)
end

"""
    group_windows(lams, n) -> Vector{Vector{Float64}}

`spectralbin` windows holding `n` consecutive channels each (the last may hold fewer), from the
distinct wavelengths `lams`: a window's edges are halfway to the neighbouring channels, so every
channel lies in exactly one window.
"""
function group_windows(lams::AbstractVector{<:Real}, n::Integer)
    ew = sort(Float64.(lams))
    ew = ew[[true; diff(ew) .> 1e-6 .* ew[2:end]]]
    m = length(ew)
    m == 0 && return Vector{Float64}[]
    half(i) = i < m ? (ew[i + 1] - ew[i]) / 2 : (m > 1 ? (ew[m] - ew[m - 1]) / 2 : ew[m] * 1e-6)
    lo(i) = i > 1 ? (ew[i - 1] + ew[i]) / 2 : ew[1] - half(1)
    return [[lo(a), ew[b] + half(b)] for (a, b) in ((a, min(a + n - 1, m)) for a in 1:n:m)]
end

"""
    shell_tools_group(path, n) -> String

Group every `n` consecutive channels of the file into a spectral bin, each channel kept
(`readoifits(path; spectralbin)`), as a new dataset: a polychromatic dataset for Imaging.
"""
function shell_tools_group(path::AbstractString, n::Real)
    n = round(Int, n)
    sh = _shell()
    p = strip_file_url(String(path))
    isempty(p) && return "! choose an OIFITS file"
    n >= 1 || return "! group at least 1 channel"
    w = try
        group_windows(reduce(vcat, [Float64.(d.uv_lam) for d in _file_entry(p).data]), n)
    catch err
        return console!(sh, "! could not group $(basename(p)): " * _cause(err); kind = :err)
    end
    isempty(w) && return "! the file holds no wavelength"
    entry = try
        _with_log(sh) do
            read_file_dataset!(sh.session, p, Dict{Symbol,Any}(:spectralbin => w); label = "grouped ×$(n)")
        end
    catch err
        return console!(sh, "! could not group $(basename(p)): " * _cause(err); kind = :err)
    end
    _show_newest!(sh, "grouped into $(entry.name)")
    return console!(sh, "grouped every $(n) channels of $(basename(p)) into $(entry.name): " *
                        "$(size(entry.data, 1)) spectral bins")
end

_split_by(by::AbstractString) =
    by == "gap" ? :gap : by == "width" ? :width : by == "target" ? :target :
    by == "instrument" ? :instrument : throw(ArgumentError("unknown split: " * by))

"""
    shell_tools_split_preview(path, by, param) -> String

The parts a split of the file would make, one per line: `label, V² / T3 / VIS points`
(tab-separated); a target's points are only known once read, so a split by target shows "—".
`by` is `gap`, `width` (both take `param` in days), `target` or `instrument`.
"""
function shell_tools_split_preview(path::AbstractString, by::AbstractString, param::Real = 0.5)
    p = strip_file_url(String(path))
    isempty(p) && return "! choose an OIFITS file"
    e = try _file_entry(p) catch err; return "! " * _cause(err) end
    parts = split_parts(e, _split_by(by); gap = param, binwidth = param)
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
    shell_tools_split(path, by, param, folder) -> String

Split the file into new datasets; with `folder` not empty, also write each part to
`<folder>/<file>_<part>.oifits` (an existing file is not replaced, nor ever the source).
"""
function shell_tools_split(path::AbstractString, by::AbstractString, param::Real = 0.5,
                           folder::AbstractString = "")
    sh = _shell()
    p = strip_file_url(String(path))
    isempty(p) && return "! choose an OIFITS file"
    entries = try
        _with_log(sh) do
            split_file!(sh.session, p, _split_by(by); gap = param, binwidth = param)
        end
    catch err
        return console!(sh, "! could not split $(basename(p)): " * _cause(err); kind = :err)
    end
    written = 0
    dir = strip_file_url(String(folder))
    if !isempty(dir)
        stem = replace(basename(p), r"\.(oi)?fits$"i => "")
        for x in entries
            out = joinpath(dir, string(stem, "_", replace(x.label, r"[^A-Za-z0-9._+-]+" => "_"), ".oifits"))
            OITOOLS._same_file(out, p) && continue
            tn = String(get(x.readopts, :targetname, ""))
            console!(sh, "write_oifits($(_literal(out)), $(x.name)" *
                         (isempty(tn) ? "" : "; targetname = $(_literal(tn))") * ")"; kind = :cmd)
            try
                OITOOLS.write_oifits(out, x.data; targetname = tn)      # the console's record only, as Export's
                written += 1
            catch err
                console!(sh, "! not written: " * _cause(err); kind = :err)
            end
        end
    end
    _show_newest!(sh, "split")
    return console!(sh, "split $(basename(p)) into $(length(entries)) datasets" *
                        (isempty(dir) ? "" : ", $(written) written to $(dir)"))
end

"""
    shell_tools_merge(paths) -> String

Read the files (one per line) as one dataset (`readoifits(files)`: stations matched by name,
every point kept), and select it.
"""
function shell_tools_merge(paths::AbstractString)
    sh = _shell()
    ps = unique([strip_file_url(String(strip(l))) for l in split(paths, '\n') if !isempty(strip(l))])
    length(ps) >= 2 || return console!(sh, "! choose at least two files to merge"; kind = :err)
    entry = try
        _with_log(sh) do
            load_dataset!(sh.session, ps; warn = false, verbose = false)
        end
    catch err
        return console!(sh, "! could not merge: " * _cause(err); kind = :err)
    end
    entry.label = "merged"
    _show_newest!(sh, "merged")
    d = vec(entry.data)
    return console!(sh, "merged $(length(ps)) files into $(entry.name): " *
                        "$(sum(x.nv2 for x in d)) V², $(sum(x.nt3phi for x in d)) T3φ, " *
                        "$(sum(x.nvisamp for x in d)) VIS")
end

# ── the cube viewer ───────────────────────────────────────────────────────────────────

# Its canvas, built by `gui()` with the others before the window exists, and the cube on show.
const CUBE_CANVAS = Ref{Any}(nothing)
const CUBE = Ref{Any}(nothing)

_cube_planes(m) = [(w, t) for t in axes(m.image, 4) for w in axes(m.image, 3)]

function _cube_label(m, w, t)
    parts = String[]
    w <= length(m.wavelengths) && isfinite(m.wavelengths[w]) &&
        push!(parts, @sprintf("λ %.5g µm", m.wavelengths[w] * 1e6))
    t <= length(m.mjds) && isfinite(m.mjds[t]) && push!(parts, @sprintf("MJD %.4f", m.mjds[t]))
    return isempty(parts) ? "plane $(w + (t - 1) * size(m.image, 3))" : join(parts, " · ")
end

"""
    _viewer_cube(path, sh) -> (; image, pixsize, wavelengths, mjds, path)

The cube to view: `read_model_image`'s reading where it accepts the file. Otherwise the planes as
stored, which a viewer can show although a model cannot use them: a third axis that does not
say what it is (no `CTYPE3`) is wavelength when `CUNIT3` is a length (planes at `CRVAL3 +
(k - CRPIX3)·CDELT3`), else unlabelled; the pixel scale from the WCS when there is one, 1 mas
(said in the console) when not.
"""
# Axis 3's wavelengths (metres) from a linear WCS whose CUNIT3 is a length, CTYPE3 or not; empty
# when the header does not give them.
function _axis3_wavelengths(kw, n)
    unit = Dict("m" => 1.0, "um" => 1e-6, "µm" => 1e-6, "micron" => 1e-6, "microns" => 1e-6,
                "nm" => 1e-9, "angstrom" => 1e-10)
    u = lowercase(strip(string(get(kw, "CUNIT3", ""))))
    (n > 1 && haskey(unit, u) && get(kw, "CRVAL3", nothing) isa Real &&
     get(kw, "CDELT3", nothing) isa Real) || return Float64[]
    c0, d, r = Float64(kw["CRVAL3"]), Float64(kw["CDELT3"]), Float64(get(kw, "CRPIX3", 1.0))
    return [(c0 + (k - r) * d) * unit[u] for k in 1:n]
end

function _viewer_cube(p::AbstractString, sh)
    try
        return OITOOLS.read_model_image(p)
    catch
    end
    raw = Float64.(OITOOLS.readfits(p))
    ndims(raw) in (2, 3, 4) || throw(ArgumentError("$(ndims(raw))-D data, not an image or a cube"))
    img = reshape(raw, size(raw, 1), size(raw, 2), size(raw, 3), size(raw, 4))
    kw = first(OITOOLS._fits_header(p))
    wl = _axis3_wavelengths(kw, size(img, 3))
    ps = OITOOLS.fits_pixsize(p)
    ps === nothing && console!(sh, "  no pixel scale in the header; shown at 1 mas per pixel")
    console!(sh, "  read as stored planes" * (isempty(wl) ? ", unlabelled" : ", wavelengths from CUNIT3/CRVAL3/CDELT3"))
    return (; image = img, pixsize = something(ps, 1.0), wavelengths = wl, mjds = Float64[],
              path = String(p))
end

"""
    shell_cube_open(path) -> String

Read a FITS image or cube (`read_model_image`: East first, its pixel scale, wavelengths and
epochs) for the cube viewer, and show its first plane. Returns one line per plane, its label.
"""
function shell_cube_open(path::AbstractString)
    sh = _shell()
    p = strip_file_url(String(path))
    console!(sh, "read_model_image($(_literal(p)))"; kind = :cmd)
    m = try
        _viewer_cube(p, sh)
    catch err
        return console!(sh, "! not an image or a cube: " * _cause(err); kind = :err)
    end
    CUBE[] = m
    ps = _cube_planes(m)
    console!(sh, "$(size(m.image, 1))×$(size(m.image, 2)) pixels of $(m.pixsize) mas, $(length(ps)) plane" *
                 (length(ps) == 1 ? "" : "s"))
    shell_cube_show(1)
    return join([_cube_label(m, w, t) for (w, t) in ps], "\n")
end

"""
    shell_cube_show(k) -> String

Show plane `k` (1-based, wavelengths first, then epochs) of the cube on the viewer's canvas.
"""
function shell_cube_show(k::Real)
    m = CUBE[]
    m === nothing && return "! no cube open"
    ps = _cube_planes(m)
    k = clamp(round(Int, k), 1, length(ps))
    w, t = ps[k]
    lab = _cube_label(m, w, t)
    c = CUBE_CANVAS[]
    c === nothing || show_image!(c, _stretched(m.image[:, :, w, t]), m.pixsize;
                                 label = _stretch_label(), title = basename(m.path) * " — " * lab)
    return lab
end

# The viewer's stretch: planes shown as max(I, 0)^p, p = 1 the plane as it is.
const CUBE_POWER = Ref(1.0)
_stretched(img) = CUBE_POWER[] == 1 ? img : max.(img, 0) .^ CUBE_POWER[]
_stretch_label() = CUBE_POWER[] == 1 ? "flux / pixel" : "(flux / pixel)^" * _g(CUBE_POWER[])

"""
    shell_cube_power(p) -> String

The viewer's stretch: planes shown as max(I, 0)^p (p > 0; 1 linear), for the single plane and
the planes side by side.
"""
function shell_cube_power(p::AbstractString)
    v = tryparse(Float64, strip(p))
    (v === nothing || !(v > 0)) && return "! the power must be a positive number"
    CUBE_POWER[] = v
    return "ok"
end

# ── channels side by side ─────────────────────────────────────────────────────────────

# The panels, built by `gui()` with the other canvases; what they show (channels, epoch).
const CUBE_TILES = Ref{Any}(nothing)
const CUBE_NTILES = 6
const CUBE_TILES_SHOWN = Ref{Any}(nothing)

"""
    build_cube_tiles(fig) -> NamedTuple

`CUBE_NTILES` image panels in a row over one shared colour bar, for the cube viewer's channels
side by side: `figure`, `axis` (the first panel), `axes`, and each panel's `xs`, `ys` and `data`
observables, the shared `crange`, `label`. Panels not in use are hidden and given no width.
"""
function build_cube_tiles(fig)
    n = CUBE_NTILES
    xs = [Makie.Observable(Float32[0, 1]) for _ in 1:n]
    ys = [Makie.Observable(Float32[0, 1]) for _ in 1:n]
    data = [Makie.Observable(zeros(Float32, 2, 2)) for _ in 1:n]
    crange = Makie.Observable((0.0f0, 1.0f0))
    label = Makie.Observable("(flux / pixel)^0.5")       # its glyphs met before the window exists
    axs = map(1:n) do i
        # the title seeded with a label's glyphs, as `build_canvas` seeds its colour bar's
        ax = Makie.Axis(fig[1, i]; aspect = Makie.DataAspect(), title = "λ 2.1661 µm · MJD 0",
                        xlabel = "α (mas)", ylabel = i == 1 ? "δ (mas)" : "")
        style_axis!(ax)
        Makie.heatmap!(ax, xs[i], ys[i], data[i]; colormap = LIVE_COLORMAP, colorrange = crange)
        ax
    end
    Makie.Colorbar(fig[2, 1:n]; colormap = LIVE_COLORMAP, limits = crange, label,
                   vertical = false, flipaxis = false)
    for i in 2:n
        Makie.colsize!(fig.layout, i, Makie.Fixed(0)); axs[i].blockscene.visible[] = false
    end
    return (; figure = fig, axis = axs[1], axes = axs, xs, ys, data, crange, label)
end

# The channels a list names: numbers with a decimal point or exponent are wavelengths in µm (the
# nearest channel), plain integers channel numbers; empty, four channels evenly spread.
function _tile_channels(m, spec::AbstractString)
    nλ = size(m.image, 3)
    toks = split(spec, r"[,;\s]+"; keepempty = false)
    isempty(toks) && return unique(round.(Int, range(1, nλ; length = min(4, nλ))))
    out = Int[]
    for s in toks
        v = tryparse(Float64, s)
        v === nothing && return "! not a channel or a wavelength: $(s)"
        if occursin(r"[.eE]", s)
            all(isfinite, m.wavelengths) ||
                return "! the cube's channels have no wavelengths; give channel numbers"
            push!(out, argmin(abs.(m.wavelengths .* 1e6 .- v)))
        else
            1 <= v <= nλ || return "! channel $(s): the cube has $(nλ)"
            push!(out, Int(v))
        end
    end
    length(out) > CUBE_NTILES && return "! at most $(CUBE_NTILES) channels side by side"
    return out
end

# Draw channels `ws` of epoch `t` on the panels `T`: one colour scale for all, the viewer's
# stretch, East to the left.
function _show_tiles!(T, m, ws, t)
    sc = live_plot_scale()
    imgs = [_stretched(m.image[:, :, w, t]) for w in ws]
    for i in 1:CUBE_NTILES
        ax = T.axes[i]
        if i > length(ws)
            Makie.colsize!(T.figure.layout, i, Makie.Fixed(0)); ax.blockscene.visible[] = false
            continue
        end
        Makie.colsize!(T.figure.layout, i, Makie.Auto())
        ax.blockscene.visible[] = true
        style_axis!(ax; scale = sc)
        # few ticks: six panels share the width
        ax.xticks = Makie.WilkinsonTicks(3); ax.yticks = Makie.WilkinsonTicks(3)
        img = imgs[i]
        hx, hy = size(img, 1) * m.pixsize / 2, size(img, 2) * m.pixsize / 2
        T.xs[i][] = Float32.(range(-hx, hx; length = size(img, 1)))
        T.ys[i][] = Float32.(range(-hy, hy; length = size(img, 2)))
        T.data[i][] = Float32.(img)
        ax.title = _cube_label(m, ws[i], t)
        Makie.limits!(ax, hx, -hx, -hy, hy)
    end
    mx = maximum(maximum, imgs)
    T.crange[] = (0.0f0, Float32(mx > 0 ? mx : 1))
    T.label[] = _stretch_label()
    return nothing
end

"""
    shell_cube_tiles(spec, t) -> String

Show the open cube's channels side by side, of epoch `t`: `spec` lists channel numbers or
wavelengths in µm (a decimal point), at most six; empty, four channels evenly spread. Returns
the channels' labels, tab-separated.
"""
function shell_cube_tiles(spec::AbstractString, t::Real)
    m = CUBE[]
    m === nothing && return "! no cube open"
    ws = _tile_channels(m, String(spec))
    ws isa String && return ws
    t = clamp(round(Int, t), 1, size(m.image, 4))
    T = CUBE_TILES[]
    T === nothing || _show_tiles!(T, m, ws, t)
    CUBE_TILES_SHOWN[] = (; ws, t)
    return join([_cube_label(m, w, t) for w in ws], "\t")
end

"""
    shell_cube_axes() -> String

The open cube's two channel axes, three lines: `nwave\tnepoch`, then each wavelength channel's
value in µm and each epoch's MJD, tab-separated (blank where the file does not say).
"""
function shell_cube_axes()
    m = CUBE[]
    m === nothing && return "! no cube open"
    nλ, nt = size(m.image, 3), size(m.image, 4)
    lab(v, k, f) = k <= length(v) && isfinite(v[k]) ? f(v[k]) : ""
    wl = [lab(m.wavelengths, k, x -> @sprintf("λ=%.5g µm", x * 1e6)) for k in 1:nλ]
    mj = [lab(m.mjds, k, x -> @sprintf("MJD %.4f", x)) for k in 1:nt]
    return join(["$(nλ)\t$(nt)", join(wl, "\t"), join(mj, "\t")], "\n")
end

"""
    shell_cube_show_at(w, t) -> String

Show wavelength channel `w` of epoch `t` (both 1-based, clamped) on the viewer's canvas.
"""
function shell_cube_show_at(w::Real, t::Real)
    m = CUBE[]
    m === nothing && return "! no cube open"
    nλ = size(m.image, 3)
    w = clamp(round(Int, w), 1, nλ); t = clamp(round(Int, t), 1, size(m.image, 4))
    return shell_cube_show((t - 1) * nλ + w)
end

# ── Set keywords / Save FITS ──────────────────────────────────────────────────────────

# What Set keywords applied, for Save FITS: the source file and `fits_set_keywords`' arguments.
const FITSKW = Ref{Any}(nothing)

const _PIX_MAS = Dict("mas" => 1.0, "µas" => 1e-3, "arcsec" => 1e3, "deg" => 3.6e6,
                      "rad" => 3.6e6 * 180 / π)
const _WAVE_M  = Dict("m" => 1.0, "µm" => 1e-6, "nm" => 1e-9, "Å" => 1e-10)
const _STEP_D  = Dict("d" => 1.0, "h" => 1 / 24, "min" => 1 / 1440, "s" => 1 / 86400)

_g(x) = @sprintf("%.10g", x)

"""
    shell_fits_keywords(path) -> String

What Set keywords starts from, as `key\tvalue` lines read from the file: `dims` (the array's
sizes), `pixsize` (mas; blank when the header has no scale), `eastleft` (1 when its WCS stores
East on the left), `axis3` (`wave` or `time`), `wstart`/`wstep` (µm) and `tstart`/`tstep` (MJD,
days), `refx`/`refy` (the reference pixel, CRPIX1/2), blank where the file does not say — axis 3's
wavelengths also from a CUNIT3/CRVAL3/CDELT3 without CTYPE3; then one `card\tNAME\t1\tvalue` (present) or `card\tNAME\t0\t` (missing) line per
WCS card of each axis (CTYPE, CUNIT, CRVAL, CDELT or a CD entry, CRPIX).
"""
function shell_fits_keywords(path::AbstractString)
    p = strip_file_url(String(path))
    dims = try
        OITOOLS.AstroFITS.FitsFile(f -> collect(f[1].data_size), p)
    catch err
        return "! cannot read $(basename(p)): " * _cause(err)
    end
    2 <= length(dims) <= 4 || return "! $(basename(p)): $(length(dims))-D data, not an image or a cube"
    h = first(OITOOLS._fits_header(p))
    ps = OITOOLS.fits_pixsize(p)
    cd1 = get(h, "CDELT1", 0.0)
    eastleft = startswith(uppercase(string(get(h, "CTYPE1", ""))), "RA") && cd1 isa Real && cd1 < 0
    axis3 = uppercase(strip(string(get(h, "CTYPE3", "")))) in ("TIME", "MJD") ? "time" : "wave"
    wl, mj = try
        m = OITOOLS.read_model_image(p)
        (m.wavelengths, m.mjds)
    catch
        (_axis3_wavelengths(h, length(dims) > 2 ? dims[3] : 1), Float64[])
    end
    first2(v, scale) = isempty(v) || !isfinite(v[1]) ? ("", "") :
        (_g(v[1] * scale), length(v) > 1 && isfinite(v[2]) ? _g((v[2] - v[1]) * scale) : "")
    ws, wd = first2(wl, 1e6)
    ts, td = first2(mj, 1.0)
    rp = OITOOLS.fits_refpix(p)
    # the WCS cards of each axis: present (with its value) or missing; a CD matrix entry
    # stands for CDELT on axes 1 and 2
    cards = String[]
    for k in eachindex(dims), c in ("CTYPE", "CUNIT", "CRVAL", "CDELT", "CRPIX")
        name = "$(c)$(k)"
        v = get(h, name, nothing)
        v === nothing && c == "CDELT" && k <= 2 && (v = get(h, "CD$(k)_$(k)", nothing))
        push!(cards, "card\t$(name)\t" * (v === nothing ? "0\t" : "1\t" * strip(string(v))))
    end
    return join(["dims\t" * join(dims, " "), "pixsize\t" * (ps === nothing ? "" : _g(ps)),
                 "eastleft\t" * (eastleft ? "1" : "0"), "axis3\t" * axis3,
                 "wstart\t" * ws, "wstep\t" * wd, "tstart\t" * ts, "tstep\t" * td,
                 "refx\t" * (rp === nothing ? "" : _g(rp[1])), "refy\t" * (rp === nothing ? "" : _g(rp[2])),
                 cards...], "\n")
end

"""
    shell_fits_apply(path, pixsize, pixunit, eastleft, axis3, wstart, wstep, wunit,
                     tstart, tstep, tunit, refx = "", refy = "") -> String

Set keywords: check the values against the file's planes, keep them for Save FITS, and show the
cube with them (nothing is written). Pixel size in `pixunit` (mas, µas, arcsec, deg, rad);
wavelengths `wstart + k·wstep` in `wunit` (m, µm, nm, Å); epochs `tstart + k·tstep`, the start an
MJD and the step in `tunit` (d, h, min, s); the reference pixel `refx`, `refy` in the file's
numbering (CRPIX1/2), both blank for the phase centre nx÷2+1. Returns the planes' labels, one per
line.
"""
function shell_fits_apply(path::AbstractString, pixsize::AbstractString, pixunit::AbstractString,
                          eastleft::Bool, axis3::AbstractString,
                          wstart::AbstractString, wstep::AbstractString, wunit::AbstractString,
                          tstart::AbstractString, tstep::AbstractString, tunit::AbstractString,
                          refx::AbstractString = "", refy::AbstractString = "")
    sh = _shell()
    p = strip_file_url(String(path))
    num(x) = isempty(strip(x)) ? nothing : tryparse(Float64, strip(x))
    ps = num(pixsize)
    (ps === nothing || !(ps > 0) || !haskey(_PIX_MAS, pixunit)) &&
        return "! the pixel size must be a positive number"
    ax = axis3 == "time" ? :time : :wave
    cube = try
        OITOOLS.fits_keyword_cube(p; east_left = eastleft, axis3 = ax)
    catch err
        return "! " * _cause(err)
    end
    nλ, nt = size(cube, 3), size(cube, 4)
    function axis(n, s0, ds, scale, what)
        n == 1 && num(s0) === nothing && return Float64[]
        a, d = num(s0), num(ds)
        a === nothing && return "! the $(what) start must be a number"
        n > 1 && (d === nothing || d == 0) && return "! the $(what) step must be a non-zero number"
        return [(a + (k - 1) * something(d, 0.0)) * scale for k in 1:n]
    end
    wl = axis(nλ, wstart, wstep, _WAVE_M[wunit], "wavelength")
    wl isa String && return wl
    t0 = num(tstart); dt = num(tstep)
    mj = nt == 1 && t0 === nothing ? Float64[] :
         t0 === nothing ? "! the epoch start (MJD) must be a number" :
         nt > 1 && (dt === nothing || dt == 0) ? "! the epoch step must be a non-zero number" :
         [t0 + (k - 1) * something(dt, 0.0) * _STEP_D[tunit] for k in 1:nt]
    mj isa String && return mj
    rx, ry = num(refx), num(refy)
    isempty(strip(refx)) && isempty(strip(refy)) || (rx !== nothing && ry !== nothing) ||
        return "! the reference pixel needs two numbers, or none for the centre"
    rp = rx === nothing ? nothing : (rx, ry)
    mas = ps * _PIX_MAS[pixunit]
    FITSKW[] = (; path = p, pixsize = mas, east_left = eastleft, axis3 = ax,
                  wavelengths = wl, mjds = mj, refpix = rp, nx = size(cube, 1))
    m = (; image = cube, pixsize = mas, wavelengths = wl, mjds = mj, path = p)
    CUBE[] = m
    console!(sh, "keywords set for $(basename(p)): $(_g(mas)) mas" *
                 (isempty(wl) ? "" : ", $(nλ) wavelength" * (nλ == 1 ? "" : "s")) *
                 (isempty(mj) ? "" : ", $(nt) epoch" * (nt == 1 ? "" : "s")) * " (not saved)")
    shell_cube_show(1)
    return join([_cube_label(m, w, t) for (w, t) in _cube_planes(m)], "\n")
end

"The file whose keywords Set keywords holds, or \"\"."
shell_fits_keyword_source() = FITSKW[] === nothing ? "" : FITSKW[].path

# A vector as the code that rebuilds it: an evenly spaced one as a range.
function _range_literal(v)
    length(v) == 1 && return "[$(_g(v[1]))]"
    return "collect(range($(_g(v[1])), step = $(_g(v[2] - v[1])), length = $(length(v))))"
end

"""
    shell_fits_save(out) -> String

Save FITS: write the file Set keywords was applied to, with those keywords, to `out`
(`fits_set_keywords`; `out` may be the file itself).
"""
function shell_fits_save(out::AbstractString)
    sh = _shell()
    k = FITSKW[]
    k === nothing && return "! set the keywords first"
    o = strip_file_url(String(out))
    args = "pixsize = $(_g(k.pixsize))" * (k.east_left ? ", east_left = true" : "") *
           (k.axis3 === :time ? ", axis3 = :time" : "") *
           (isempty(k.wavelengths) ? "" : ", wavelengths = " * _range_literal(k.wavelengths)) *
           (isempty(k.mjds) ? "" : ", mjds = " * _range_literal(k.mjds)) *
           ", refpix = " * (k.refpix === nothing ? "nothing" : "($(_g(k.refpix[1])), $(_g(k.refpix[2])))")
    console!(sh, "fits_set_keywords($(_literal(k.path)), $(_literal(o)); $(args))"; kind = :cmd)
    try
        OITOOLS.fits_set_keywords(k.path, o; pixsize = k.pixsize, east_left = k.east_left,
                                  axis3 = k.axis3, wavelengths = k.wavelengths, mjds = k.mjds,
                                  refpix = k.refpix)
    catch err
        return console!(sh, "! could not save $(basename(o)): " * _cause(err); kind = :err)
    end
    # an overwritten file is now East-first in the convention's axis order
    o == k.path && (FITSKW[] = merge(k, (; east_left = false, axis3 = :wave,
        refpix = k.refpix === nothing || !k.east_left ? k.refpix : (k.nx + 1 - k.refpix[1], k.refpix[2]))))
    return console!(sh, "wrote " * basename(o))
end

# ── fit each dataset ──────────────────────────────────────────────────────────────────

# The last "Fit each" table: one row per dataset.
const EPOCHS = Ref{Any}(nothing)

"""
    shell_fit_each(names, free_lines, constraint_lines, prior_lines,
                   v2, t3amp, t3phi, visamp, visphi, flux, diffvis,
                   optimiser, maxeval, uncertainty) -> String

Fit the panel's model to each dataset named (one per line), one after another on the worker, as a
single fit runs: the same free parameters, bounds, constraints, priors, observables and optimiser.
Every fit starts from the model as it stands in the table, so the results do not depend on the
order. Levenberg-Marquardt or an NLopt method; uncertainty `"jacobian"` or `"none"`. The table
(`shell_epochs_table`) is ready when the job finishes.
"""
function shell_fit_each(names::AbstractString, free_lines::AbstractString,
                        constraint_lines::AbstractString, prior_lines::AbstractString,
                        v2::Bool, t3amp::Bool, t3phi::Bool, visamp::Bool, visphi::Bool,
                        flux::Bool, diffvis::Bool, optimiser::AbstractString, maxeval::Real,
                        uncertainty::AbstractString = "jacobian")
    sh = _shell()
    sh.job === nothing || return "! a fit or a reconstruction is running"
    ns = [String(strip(l)) for l in split(names, '\n') if !isempty(strip(l))]
    isempty(ns) && return "! tick at least one dataset"
    es = try
        [dataset(sh.session, n) for n in ns]
    catch err
        return "! " * _cause(err)
    end
    md = copy(_model(sh).dict)
    isempty(md) && return "! no model"
    free, lb, ub = parse_free_lines(free_lines)
    isempty(free) && return "! no free parameters"
    cons = try
        parse_constraint_lines(constraint_lines)
    catch err
        return "! constraint: " * _cause(err)
    end
    priors  = parse_prior_lines(prior_lines)
    weights = fitting_weights(; v2, t3amp, t3phi, visamp, visphi, flux, diffvis)
    all(iszero, weights) && return "! no observables selected"
    opt = String(optimiser)
    opt in ("grid", "nested", "ensemble", "pigeons") &&
        return "! fit each dataset with Levenberg-Marquardt or an NLopt method, not $(opt == "grid" ? "a grid search" : "a sampler")"
    unc = String(uncertainty)
    unc in ("none", "jacobian") || return "! fit each dataset with a Jacobian uncertainty or none"
    for n in ns
        console!(sh, "> fit_model(model, $(free), $(n)[1, 1]; weights = $(weights), method = :$(opt))")
    end
    sh.fitlog = ""
    optsym = Symbol(opt)
    start_job!(sh, :epochs, function (stop)
        rows = Any[]
        for (n, e) in zip(ns, es)
            stop[] && break
            d = e.data[1, 1]
            println("── ", n, " (", source_label(e), ")")
            row = try
                r = opt == "lsqfit" ?
                    fit_model_lsqfit(md, free, d; lb, ub, weights, constraints = cons,
                                     maxIter = Int(maxeval), verb = false) :
                    fit_model(md, free, d; lb, ub, weights, priors, constraints = cons,
                              method = optsym, maxeval = Int(maxeval), verb = false)
                sig = unc == "jacobian" ? _jacobian_sigma(r) : nothing
                (; name = n, mjd = d.mean_mjd, chi2r = r.chi2r, ndof = r.ndof,
                   x = collect(r.x_opt), sigma = sig, err = "")
            catch err
                (; name = n, mjd = d.mean_mjd, chi2r = NaN, ndof = 0,
                   x = fill(NaN, length(free)), sigma = nothing, err = _cause(err))
            end
            println(isempty(row.err) ? @sprintf("   χ²r %.4f", row.chi2r) : "   failed: " * row.err)
            push!(rows, row)
        end
        return (; rows, free)
    end)
    return "fitting $(length(ns)) dataset" * (length(ns) == 1 ? "" : "s") * "…"
end

"""
    finish_epochs!(sh, res) -> String

Keep a finished "Fit each" table and report it in the console. **GUI thread only.**
"""
function finish_epochs!(sh::ShellState, res)
    if !res.ok
        msg = "! fitting each dataset failed: " * _cause(res.err)
        sh.fitlog = _collapse_cr(res.output) * "\n" * msg
        return console!(sh, msg)
    end
    sh.fitlog = _collapse_cr(res.output)
    EPOCHS[] = res.result
    for r in res.result.rows
        console!(sh, isempty(r.err) ?
                     @sprintf("  %-10s MJD %.4f   χ²r %.4f   ", r.name, r.mjd, r.chi2r) *
                     join([@sprintf("%s = %.6g", k, v) for (k, v) in zip(res.result.free, r.x)], "   ") :
                     "! $(r.name): " * r.err)
    end
    nok = count(r -> isempty(r.err), res.result.rows)
    return "fitted $(nok) of $(length(res.result.rows)) datasets"
end

"""
    shell_epochs_table() -> String

The last "Fit each" table, tab-separated with a header: `dataset, MJD, chi2r, ndof`, then each
free parameter and its σ (empty without an uncertainty).
"""
function shell_epochs_table()
    t = EPOCHS[]
    t === nothing && return ""
    head = vcat(["dataset", "MJD", "chi2r", "ndof"], [x for k in t.free for x in (k, "σ " * k)])
    lines = [join(head, "\t")]
    for r in t.rows
        cells = String[r.name, @sprintf("%.5f", r.mjd),
                       isfinite(r.chi2r) ? @sprintf("%.4f", r.chi2r) : "", string(r.ndof)]
        for j in eachindex(t.free)
            push!(cells, isfinite(r.x[j]) ? @sprintf("%.6g", r.x[j]) : "")
            push!(cells, r.sigma === nothing || !isfinite(r.sigma[j]) ? "" : @sprintf("%.3g", r.sigma[j]))
        end
        push!(lines, join(cells, "\t"))
    end
    return join(lines, "\n")
end

"""Write the last "Fit each" table to `path` (tab-separated)."""
function shell_save_epochs(path::AbstractString)
    t = shell_epochs_table()
    isempty(t) && return "! no table to save"
    p = strip_file_url(String(path))
    write(p, t * "\n")
    return console!(_shell(), "saved the table to " * basename(p))
end
