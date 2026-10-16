# Where the shipped resources live.
#
# QML, the instrument configs, the demo data and the demo models are files on disk that the
# package reads at RUN time, not code that a sysimage can carry. Everything else in this
# package is compiled; these are not.
#
# The reason this file exists is that `pkgdir` is not an answer once the package is inside an
# application bundle. `PackageCompiler.create_app` bundles the sysimage, the artifacts and
# Julia's libraries, and writes only a stub `Project.toml` -- no package source is copied. So
# `pkgdir(OITOOLS)` then names a directory on the machine that BUILT the app, which does not
# exist on the machine running it, and every `joinpath(pkgdir(OITOOLS), ...)` silently misses.
#
# Measured, with each resource removed in turn from a copy of the repository:
#
#   src/configs   gone -> the GUI starts and quits cleanly, facilities 5 -> 0, no warning
#   src/gui/qml   gone -> fatal, "Failed to load QML file ..." out of loadqml
#   demos/data    gone -> silent, picker places 7 -> 3, examples ""
#
# The fatal one masks the silent ones, which is why all the call sites move together rather
# than one at a time in response to whatever the application complains about first.
#
# A BUNDLE MIRRORS THE REPOSITORY'S OWN RELATIVE PATHS -- `src/gui/qml`, `src/configs`,
# `demos/data`, `demos/models` under the resource root. That makes `resource` a root
# substitution and nothing else. A mapping from repository layout to bundle layout would be a
# second list to keep in step with the file picker's places, and it would drift the first time
# a place was added.

"Environment variable that names the resource root outright, overriding every other candidate."
const RESOURCE_DIR_VAR = "OITOOLS_RESOURCE_DIR"

"""
    resource_dir() -> String or nothing

The directory holding the shipped resources, or `nothing` if there is none.

The first existing one of:

 1. `\$OITOOLS_RESOURCE_DIR`;
 2. `<Sys.BINDIR>/../share/oitools`, in an application bundle;
 3. `pkgdir(OITOOLS)`, a development checkout.
"""
function resource_dir()
    for r in _resource_roots()
        isdir(r) && return r
    end
    return nothing
end

function _resource_roots()
    roots = String[]
    forced = get(ENV, RESOURCE_DIR_VAR, "")
    isempty(forced) || push!(roots, abspath(expanduser(forced)))
    push!(roots, normpath(joinpath(Sys.BINDIR, "..", "share", "oitools")))
    root = pkgdir(@__MODULE__)
    root === nothing || push!(roots, root)
    return roots
end

"""
    resource(parts...) -> String or nothing

One shipped file or directory, named by its path relative to the repository root, or `nothing`
when no root has it.

    resource("src", "gui", "qml", "Main.qml")
    resource("demos", "data")

`nothing` rather than a non-existent path, because every caller but one treats a missing
resource as "this feature has nothing to offer" and carries on -- a picker place is simply not
listed, an examples directory is empty. The exception is `Main.qml`, whose absence has no
sensible fallback and which therefore reports the path it looked for.
"""
function resource(parts::AbstractString...)
    sub = joinpath(parts...)
    for r in _resource_roots()
        p = joinpath(r, sub)
        ispath(p) && return p
    end
    return nothing
end
