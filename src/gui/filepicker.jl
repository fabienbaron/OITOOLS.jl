# The file picker's data layer.
#
# Why this package has its own picker at all: `QtQuick.Dialogs.FileDialog` sets `visible =
# false` and emits `accepted`, but on some systems never unmaps its window — measured with
# `test/gui/filepicker_min.jl`, which is that dialog and nothing else. The dialog is a separate
# top-level window, and where Qt cannot get a GL surface for it (EGL failing over to zink, say)
# hiding it does not take it off the screen. Nothing in QML can close a window Qt will not
# close.
#
# A `Popup` is drawn on the parent window's own scene graph, so no second window is ever
# created and none can be left behind. That removes the failure mode by construction, and the
# price is that the toolkit no longer supplies the directory listing, the places, or the
# filtering. This file is that price.
#
# Everything returns a string, because that is what crosses into QML cleanly: rows separated by
# newlines, fields by tabs. Nothing here draws.

"Fields of one row of a listing, in order."
const PICKER_FIELDS = ("kind", "name", "size", "modified", "sortkey")

_human_size(n::Integer) =
    n < 1024          ? string(n, " B") :
    n < 1024^2        ? string(round(n / 1024;    digits = 1), " kB") :
    n < 1024^3        ? string(round(n / 1024^2;  digits = 1), " MB") :
                        string(round(n / 1024^3;  digits = 2), " GB")

function _modified(path)
    try
        t = Dates.unix2datetime(mtime(path))
        return Dates.format(t, "yyyy-mm-dd HH:MM")
    catch
        return ""
    end
end

"""
    picker_matches(name, patterns) -> Bool

Does this filename match any of the glob patterns (`"*.oifits"`, `"*"`)?

Only `*` is honoured, and only as a suffix or the whole pattern — which covers every filter
this application uses and avoids pulling in a glob library for the sake of one wildcard.
Matching is case-insensitive: `.FITS` is a FITS file.
"""
function picker_matches(name::AbstractString, patterns)
    isempty(patterns) && return true
    n = lowercase(name)
    for p in patterns
        p = lowercase(strip(String(p)))
        (isempty(p) || p == "*" || p == "*.*") && return true
        if startswith(p, "*")
            endswith(n, p[2:end]) && return true
        elseif n == p
            return true
        end
    end
    return false
end

"""
    picker_list(dir, patterns, show_hidden) -> String

One directory, as `kind\\tname\\tsize\\tmodified\\tsortkey` rows.

Directories come first and are never filtered out — a filter says which files you want to open,
not which folders you may walk through, and hiding folders that contain no matching file today
is how a picker becomes unusable. Directories are listed whatever the filter says.

`sortkey` is what QML sorts on without having to know the rule: directories sort before files,
then case-insensitively by name.
"""
function picker_list(dir::AbstractString, patterns::AbstractString = "*",
                     show_hidden::Bool = false)
    d = abspath(expanduser(String(dir)))
    isdir(d) || return ""
    pats = [strip(p) for p in split(patterns, ',') if !isempty(strip(p))]

    names = try
        readdir(d)
    catch err
        return "error\t" * sprint(showerror, err) * "\t\t\t0"
    end

    rows = String[]
    for n in sort(names; by = lowercase)
        (!show_hidden && startswith(n, ".")) && continue
        full = joinpath(d, n)
        isdirectory = try; isdir(full); catch; false; end
        if isdirectory
            push!(rows, join(("dir", n, "", _modified(full), "0" * lowercase(n)), "\t"))
        else
            picker_matches(n, pats) || continue
            sz = try; filesize(full); catch; 0; end
            push!(rows, join(("file", n, _human_size(sz), _modified(full),
                              "1" * lowercase(n)), "\t"))
        end
    end
    return join(rows, "\n")
end

"""
    picker_volumes() -> Vector{Tuple{String,String}}

The filesystem roots this platform expects to see listed, as `label => path` pairs.

Every place below `picker_places` is somewhere INSIDE one tree. On Unix that is enough, because
everything hangs off `/` and a user can walk there. On Windows it is not: `D:` is not reachable
from `C:` by walking up, so without this a file on another drive can only be opened by typing
its path — which is what a picker exists to avoid.

The drive scan is a bare `isdir` over the 26 letters. Windows offers `GetLogicalDrives`, but a
26-element loop needs no `ccall` and no dependency, and a drive that is mapped but not ready
(an empty optical drive, a disconnected network share) fails `isdir` — which is the answer we
want, since listing it would only produce an error on click.

Unix gets the mount points where removable media actually appear, and `/` itself, for the same
reason: they are not reachable by walking up from `\$HOME`.
"""
function picker_volumes()
    vols = Tuple{String,String}[]
    if Sys.iswindows()
        for c in 'A':'Z'
            d = string(c, ":\\")
            isdir(d) && push!(vols, (string(c, ":"), d))
        end
    elseif Sys.isapple()
        isdir("/Volumes") && for v in readdir("/Volumes"; join = true)
            isdir(v) && push!(vols, (basename(v), v))
        end
        push!(vols, ("Computer", "/"))
    else
        for m in ("/media/" * (get(ENV, "USER", "")), "/run/media/" * (get(ENV, "USER", "")),
                  "/media", "/mnt")
            isdir(m) && !isempty(readdir(m)) && push!(vols, (basename(m), m))
        end
        push!(vols, ("Filesystem", "/"))
    end
    return vols
end

"""
    picker_places() -> String

Shortcut folders, as `label\\tpath` rows.

The bundled data directory is included because it is where every example in this package
lives, and a picker that opens in the filesystem root taxes every single use.
"""
function picker_places()
    rows = Tuple{String,String}[]
    home = homedir()
    isdir(home) && push!(rows, ("Home", home))
    cwd = pwd()
    isdir(cwd) && cwd != home && push!(rows, ("Working directory", cwd))
    for (label, path) in picker_volumes()
        push!(rows, (label, path))
    end
    for (label, sub) in (("OITOOLS demo data", joinpath("demos", "data")),
                         # The 2004 contest set: small, monochromatic and well understood, and
                         # what most of the imaging examples and tests are written against.
                         ("Beauty Contest 2004", joinpath("demos", "data", "BC2004")),
                         # The only shipped files with differential phase and OI_FLUX, so the
                         # only ones on which the diffphi and flux views show anything.
                         ("Beauty Contest 2026", joinpath("demos", "data", "BC2026")),
                         # Starting points for the Model perspective, in the TOML format
                         # `read_model_file` takes -- a model to open and edit rather than one
                         # to build from an empty table.
                         ("OITOOLS models", joinpath("demos", "models")),
                         # PMOIRED's own models: Python dicts, which is what "Import PMOIRED…"
                         # reads. It already OPENS here; this is the way back after browsing
                         # somewhere else.
                         ("PMOIRED models", joinpath("demos", "data", "pmoired")))
        p = OITOOLS.resource(sub)
        (p !== nothing && isdir(p)) && push!(rows, (label, p))
    end
    return join((join(r, "\t") for r in rows), "\n")
end

"Absolute, normalised path. QML should never do path arithmetic itself."
picker_join(a, b) = abspath(joinpath(expanduser(String(a)), String(b)))

"The containing directory, or the same path when already at the root."
picker_parent(p) = (d = abspath(expanduser(String(p))); dirname(d) == d ? d : dirname(d))

"""
    picker_start(hint) -> String

Where to open. `hint` may be a directory, a file (its directory is used), or empty.

Falls back to the bundled data directory and then the working directory, so the picker always
opens somewhere with something in it.
"""
function picker_start(hint::AbstractString = "")
    h = String(hint)
    if !isempty(h)
        p = abspath(expanduser(replace(h, r"^file://" => "")))
        isdir(p) && return p
        isfile(p) && return dirname(p)
        isdir(dirname(p)) && return dirname(p)
    end
    forced = get(ENV, "OITOOLSGUI_DATA_DIR", "")
    isempty(forced) || (isdir(forced) && return abspath(forced))
    for sub in (joinpath("demos", "data"), joinpath("test", "gui", "data"))
        p = OITOOLS.resource(sub)
        p === nothing || return p
    end
    return pwd()
end

"""
    picker_examples(kind) -> String

Directory of the shipped examples for a kind of file, or `""` if there are none.

Opening an import dialog in the filesystem root asks the user to find something they have never
seen. `"pmoired"` points at models known to import: they were written by `dict_to_pmoired_file`
and read back, which matters because the importer is a transpiler and not a validator — a file
it cannot read fails at the point of use.
"""
function picker_examples(kind::AbstractString)
    k = lowercase(strip(String(kind)))
    sub = k == "pmoired" ? joinpath("demos", "data", "pmoired") :
          k == "model"   ? joinpath("demos", "models")          : ""
    isempty(sub) && return ""
    p = OITOOLS.resource(sub)
    return p === nothing ? "" : p
end

"`dir`, `file` or `none` — what QML needs to decide whether Open descends or accepts."
function picker_kind(p)
    q = abspath(expanduser(String(p)))
    isdir(q) && return "dir"
    isfile(q) && return "file"
    return "none"
end

"True when saving here would replace something, so the picker can ask first."
picker_would_overwrite(p) = isfile(abspath(expanduser(String(p)))) ? "1" : "0"
