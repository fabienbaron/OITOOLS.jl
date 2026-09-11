# Build the application bundle.
#
#     xvfb-run -a julia --project=bin app/build.jl          # headless
#     julia --project=bin app/build.jl                      # with a screen
#
# RUN IT WITH --project=bin, not --project=app. PackageCompiler is needed by this script and
# must NOT be a dependency of the application, or it would be bundled into it; `bin/` already
# carries it. What gets built is `app/`, which this script names explicitly.
#
# A DISPLAY IS NEEDED, because the precompile trace opens a GL context. Xvfb is enough.
#
# WHAT create_app DOES NOT DO, and this script therefore does:
#
#   * ship the package's data files. `bundle_project` writes a stub Project.toml and nothing
#     else -- no source, no QML, no configs, no demo data. They are staged into
#     share/oitools/, mirroring the repository's own relative paths, which is what
#     `OITOOLS.resource` expects to find under `<Sys.BINDIR>/../share/oitools`.
#   * seed Makie's font atlas. Without it the first launch tries to DOWNLOAD the atlas and,
#     failing that, re-renders every glyph.
#   * install a launcher. On Linux the windowing system has to be chosen before the process
#     starts; see app/launcher.sh.

using Pkg
using PackageCompiler

const ROOT   = normpath(joinpath(@__DIR__, ".."))
const APPSRC = joinpath(ROOT, "app")
const OUT    = get(ENV, "OITOOLS_APP_DIR", joinpath(ROOT, "build", "OITOOLS"))
const TRACE  = joinpath(ROOT, "bin", "trace.jl")
const STMTS  = joinpath(APPSRC, "precompile_statements.jl")

isfile(TRACE) || error("missing precompile trace: $TRACE")

# ── optional: link with LLD instead of BFD (Windows) ─────────────────────────
#
# OFF by default, and deliberately so. Set OITOOLS_FAST_LINK=1 to try it.
#
# The sysimage link is the longest step of a Windows build -- about 16 minutes measured,
# against a few on Linux. The reason is in the command PackageCompiler issues:
# `--export-all-symbols` over a 1.5 GB image, which makes BFD `ld` collect and sort an export
# table of millions of entries and emit an import library beside it. ELF has no such step,
# which is why the same link is quick on Linux. `--whole-archive` guarantees nothing can be
# dropped first, and BFD is single-threaded.
#
# LLD does the same job far faster and ships with Julia (`libexec/julia/lld`). Two details make
# it awkward rather than a flag:
#
#   * `-fuse-ld=lld` looks for a binary named `ld.lld`, and Julia's is a multi-flavour `lld`
#     that chooses its personality from argv[0]. So a copy under the right name is made in a
#     scratch directory and pointed at with `-B`.
#   * The only way into PackageCompiler is `JULIA_CC`, which `get_compiler_cmd` parses with
#     `Base.shell_split` -- and that treats a backslash as an escape, so a Windows path does
#     NOT survive it. Every path here is written with forward slashes, which Windows accepts
#     and shell_split leaves alone. This is what broke an earlier attempt.
#
# `JULIA_CC` applies to every link in the build, not only the slow one. If anything fails,
# unset OITOOLS_FAST_LINK and the default toolchain is used again.
if Sys.iswindows() && get(ENV, "OITOOLS_FAST_LINK", "0") != "0" && !haskey(ENV, "JULIA_CC")
    try
        lld = joinpath(Sys.BINDIR, "..", "libexec", "julia", "lld.exe")
        isfile(lld) || error("no lld at $lld")
        shim = mkpath(joinpath(tempdir(), "oitools-lld"))
        cp(lld, joinpath(shim, "ld.lld.exe"); force = true)
        gcc = first(PackageCompiler.get_compiler_cmd().exec)
        fwd(x) = replace(abspath(x), '\\' => '/')      # shell_split eats backslashes
        ENV["JULIA_CC"] = "$(fwd(gcc)) -fuse-ld=lld -B$(fwd(shim))"
        @info "OITOOLS_FAST_LINK: linking with LLD" JULIA_CC = ENV["JULIA_CC"]
    catch err
        @warn "could not set up LLD; using the default linker" err
        delete!(ENV, "JULIA_CC")
    end
end

# ── no console window on Windows ─────────────────────────────────────────────
#
# PackageCompiler links the C driver as a CONSOLE subsystem executable, so Windows opens a
# terminal and the GUI appears out of it. The subsystem is a two-byte field in the PE header,
# so it is patched on the finished executable rather than passed as a link flag.
#
# NOT via `JULIA_CC`: `get_compiler_cmd` parses that variable with `Base.shell_split`, which
# treats a backslash as an escape, so a Windows compiler path does not survive it -- the build
# fails trying to spawn `C:Usersbaron.julia…gcc.exe`. It would also apply the flag to the
# SYSIMAGE link, which is not the target.
#
# Every field is checked before the write and the current value must be 3 (console), so a
# layout this does not understand is left alone with a warning rather than corrupted.
"""
    set_windows_gui_subsystem!(exe) -> Bool

Flip a PE executable from the console subsystem to the windows one, in place.

`e_lfanew` at 0x3C gives the PE signature; the optional header follows the 4-byte signature and
the 20-byte COFF header, and `Subsystem` sits at offset 68 within it — the same place in PE32
and PE32+, since everything before it is fixed width in both.
"""
function set_windows_gui_subsystem!(exe::AbstractString)
    isfile(exe) || return false
    open(exe, "r+") do io
        read(io, 2) == b"MZ" || (@warn "not a PE file; console subsystem left as it is" exe; return false)
        seek(io, 0x3C); pe = Int(read(io, UInt32))
        seek(io, pe);   read(io, 4) == b"PE\0\0" || (@warn "no PE signature; left alone" exe; return false)
        opt = pe + 4 + 20
        seek(io, opt); magic = read(io, UInt16)
        magic in (0x10b, 0x20b) || (@warn "unknown optional header magic; left alone" exe magic; return false)
        seek(io, opt + 68); sub = read(io, UInt16)
        sub == 2 && return true                       # already a GUI binary
        sub == 3 || (@warn "unexpected subsystem; left alone" exe subsystem = sub; return false)
        seek(io, opt + 68); write(io, UInt16(2))
        return true
    end
end

@info "Building the application" out = OUT source = APPSRC
@info "This takes tens of minutes and several GB of scratch space."

t = @elapsed create_app(APPSRC, OUT;
                        force = true,
                        precompile_execution_file = TRACE,

                        # The trace above cannot call `gui()` -- its event loop never returns --
                        # so the window, the QML bridge and GLMakie's conversions for our types
                        # were compiled on the user's first launch. These statements are the
                        # union of three traces that DO reach them; see the file's own header
                        # for how to regenerate. Missing here is not an error, only a pause.
                        precompile_statements_file = isfile(STMTS) ? [STMTS] : String[],

                        # THESE TWO GO TOGETHER, and getting the pair wrong is fatal rather
                        # than merely wasteful.
                        #
                        # MKL_jll and IntelOpenMP_jll arrive through FFTW, whose provider is
                        # `fftw` (see FFT_FLAGS in oichi2.jl for why it stays that way), so
                        # FFTW loads neither. With `include_transitive_dependencies = true`
                        # they are compiled into the image anyway, and a JLL in the image runs
                        # its `__init__` whether or not anything calls it -- so omitting their
                        # lazy artifacts makes `IntelOpenMP_jll.__init__` call
                        # `find_artifact_dir` on a directory that was never bundled, and the
                        # application dies before `julia_main`. Measured in a sandbox with no
                        # depot; invisible on the build machine, where the artifact is still
                        # sitting in ~/.julia.
                        #
                        # `false` is what the option is for: the manual says it "only makes a
                        # difference if some packages do not load all their dependencies when
                        # themselves are loaded", and FFTW/MKL is exactly that case. The lazy
                        # artifacts can then be left out too -- 726 MB of MKL and IntelOpenMP
                        # that no line of this application executes.
                        include_transitive_dependencies = false,
                        include_lazy_artifacts = false,

                        # Only the stdlibs this project actually names. The risk the manual
                        # warns about is depending on one WITHOUT naming it -- `rand()` needs
                        # Random, `A * B` needs LinearAlgebra and Random both, because those
                        # stdlibs practise type piracy and merely loading them changes
                        # behaviour. The sandbox run is what checks it, since the failure would
                        # be a MethodError at run time rather than anything the build notices.
                        filter_stdlibs = true,

                        # -g0 stops DWARF being GENERATED. Julia's default is -g1 and
                        # PackageCompiler does not override it, which is where the image's
                        # ~300 MB of .debug_* comes from. It has to be suppressed here because
                        # it cannot be removed afterwards: `strip --strip-debug` produces a
                        # 1.54 GB image that then SIGSEGVs.
                        #
                        # Two heavier options are deliberately NOT used.
                        #
                        # --strip-metadata takes another 324 MB off the serialised heap, but it
                        # removes source locations from backtraces, and `julia_main`'s crash log
                        # is the only diagnostic a user can send back. 300 MB is not worth a
                        # report that says `none:11`.
                        #
                        # --strip-ir would take more still, and it is the one that breaks this
                        # package: derived model parameters and radial profiles are compiled at
                        # run time through RuntimeGeneratedFunctions and `eval`, and inlining
                        # into that new code needs the IR of what it calls.
                        sysimage_build_args = `-g0`,

                        # The default is a multiversioned target, which is what lets the binary
                        # run on a CPU other than this one. Do not narrow it: all the native
                        # code in the image is 0.21 GB of 1.84, so there is little to win here
                        # and portability to lose.
                        )

# ── the shipped resources ────────────────────────────────────────────────────
#
# Repository-relative paths, mirrored. `resource` is then a root substitution and nothing else,
# with no layout mapping to keep in step with the file picker's places.

const SHARE = joinpath(OUT, "share", "oitools")

"""
    treesize(dir) -> Int

Bytes a directory really occupies.

`filesize` FOLLOWS a symlink and reports its target, so summing it over `walkdir` counts every
linked file twice. A bundle is full of versioned `.so` links -- 805 of them here, 1.74 GB --
and totalling them that way overstated this bundle by more than half.
"""
treesize(dir) = sum(islink(p) ? 0 : filesize(p)
                    for (r, _, fs) in walkdir(dir) for p in (joinpath(r, f) for f in fs);
                    init = 0)

function stage(rel...)
    src = joinpath(ROOT, rel...)
    ispath(src) || (@warn "resource missing, not staged" src; return 0)
    dst = joinpath(SHARE, rel...)
    mkpath(dirname(dst))
    cp(src, dst; force = true)
    return treesize(dst)
end

# QMLMakie's own QML module ("Makie", supplying MakieArea) is registered by its `__init__`
# with `QML.add_import_path(joinpath(@__DIR__, "qml"))` -- and `@__DIR__` is expanded at
# PRECOMPILE time, so the path baked into the image is this machine's package directory. In a
# bundle it does not exist, Qt cannot find the module, and every .qml that imports it fails
# with `module "Makie" is not installed` before the window appears. 16 kB, and the application
# registers it from here at startup.
function stage_qmlmakie()
    pid = Base.PkgId(Base.UUID("08f9cac3-3b11-4f1c-9d88-d0e81c500f64"), "QMLMakie")
    src = Base.locate_package(pid)
    src === nothing && (@warn "QMLMakie not found; the Makie QML module is not staged"; return 0)
    from = joinpath(dirname(src), "qml")
    isdir(from) || (@warn "QMLMakie has no qml directory" from; return 0)
    to = joinpath(SHARE, "qml-modules")
    mkpath(dirname(to))
    cp(from, to; force = true)
    return treesize(to)
end

"""
    stage_x11_locale() -> bytes staged

Ship X11's compose tables, because the bundled libxkbcommon cannot find the host's.

`xkbcommon_jll` has its BinaryBuilder sandbox path compiled in as the locale directory --
`strings` on the shipped library shows `/workspace/destdir/share/X11/locale`, which exists on
no real machine. So on a WAYLAND session, where GLFW builds an xkbcommon compose table, every
launch prints

    xkbcommon: ERROR: [XKB-679] No Compose file for locale "en_US.UTF-8"
    Warning: GLFW.GLFWError(... "Wayland: Failed to create XKB compose table")

and dead keys and compose sequences stop working. Installing the host's libx11 does NOT help:
the data lands in /usr/share/X11/locale and the library is not looking there.

The library does read `XLOCALEDIR`, so `launcher.sh` points it at this copy. 1.9 MB, MIT/X11
licensed. Removing the library instead is not an option -- GLFW_jll, Qt6Base_jll,
Vulkan_Loader_jll and libdecor_jll all require it.

Same bug class as the QML import path QMLMakie bakes at precompile time: a build-time absolute
path that means nothing on the running machine.
"""
function stage_x11_locale()
    for src in ("/usr/share/X11/locale", "/usr/local/share/X11/locale")
        isfile(joinpath(src, "compose.dir")) || continue
        dst = joinpath(OUT, "share", "X11", "locale")
        mkpath(dirname(dst))
        cp(src, dst; force = true)
        return treesize(dst)
    end
    @warn "no X11 compose tables found to stage; dead keys will not work on a Wayland session " *
          "unless the user's machine has them at the standard path"
    return 0
end

staged = 0
staged += stage_qmlmakie()
staged += stage_x11_locale()
staged += stage("src", "configs")
staged += stage("src", "gui", "qml")
staged += stage("demos", "data")
staged += stage("demos", "models")
staged += stage("app", "assets")           # the icon, for the desktop entry and the installer
@info "Resources staged" dir = SHARE MB = round(staged / 1e6, digits = 1)

# ── Makie's font atlas ───────────────────────────────────────────────────────
#
# The trace above has already built it, into this machine's scratch space. Copy it beside the
# resources; the application seeds a per-user cache from there on first run, because the
# bundle's own depot is read-only once installed.

try
    Makie = Base.require(Base.PkgId(Base.UUID("ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"), "Makie"))
    cache = Base.invokelatest(Makie.get_cache_path)
    bins  = filter(f -> endswith(f, ".bin"), readdir(cache; join = true))
    if isempty(bins)
        @warn "no font atlas found; the first launch will render one" cache
    else
        dst = joinpath(SHARE, "makie-cache")
        mkpath(dst)
        for b in bins
            cp(b, joinpath(dst, basename(b)); force = true)
        end
        @info "Font atlas staged" files = length(bins) dir = dst
    end
catch err
    @warn "could not stage the font atlas; the first launch will render one" err
end

# ── the console window ───────────────────────────────────────────────────────

if Sys.iswindows()
    exe = joinpath(OUT, "bin", "OITOOLSApp.exe")
    if set_windows_gui_subsystem!(exe)
        @info "Windows: linked for the windows subsystem; no console window" exe
        @info "stdout and stderr therefore go nowhere — crash.log is the diagnostic" *
              " (see crash_log_path() in app/src/OITOOLSApp.jl)"
    end
end

# ── the launcher ─────────────────────────────────────────────────────────────

if Sys.islinux()
    cp(joinpath(APPSRC, "launcher.sh"), joinpath(OUT, "OITOOLS"); force = true)
    chmod(joinpath(OUT, "OITOOLS"), 0o755)
    # Menu entry, icon and the .oifits association. Run once after unpacking; it computes its
    # own paths, so the bundle can live anywhere and be moved by re-running it.
    cp(joinpath(APPSRC, "install-desktop.sh"), joinpath(OUT, "install-desktop.sh"); force = true)
    chmod(joinpath(OUT, "install-desktop.sh"), 0o755)
end

total = treesize(OUT)
@info "Done" seconds = round(t, digits = 1) total_MB = round(total / 1e6, digits = 1)
println("""

Start it with:

    $(joinpath(OUT, Sys.islinux() ? "OITOOLS" :
                    joinpath("bin", Sys.iswindows() ? "OITOOLSApp.exe" : "OITOOLSApp"))) [file.oifits]

$(Sys.islinux() ? """
On Linux go through the launcher rather than bin/OITOOLSApp: a bundle cannot choose its own
windowing system, because GLFW and Qt are both initialised from the sysimage before any code
in this package runs, and left alone on a Wayland session they land on different ones.
""" : """
There is no launcher on this platform and none is needed: GLFW defaults to the native
windowing system here, so the two halves agree without being told to.
""")""")
