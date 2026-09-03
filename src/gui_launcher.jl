# ===========================================================================
# oitoolsgui() — one call from a plain `using OITOOLS`.
#
# `bin/oitoolsgui.jl` exists because the graphics hints have to be set in a fixed order and
# each one has to happen BEFORE the package that reads it is loaded: Mesa and GLFW both read
# their configuration when the first OpenGL context is created, and Qt reads its platform when
# it starts. A script gets that ordering by writing the `using` lines in the right places.
#
# A function can get the same ordering, and for the same reason, as long as it is the thing
# that does the loading: `Base.require` at the right point in the body is exactly a `using` at
# the right point in a file. What it cannot do is undo a `using GLMakie` the caller already
# ran, so that case is detected and reported rather than silently producing a window on the
# wrong platform.
# ===========================================================================

const _GUI_PACKAGES = (("GLMakie",  "e9467ef8-e4e7-5192-8a1a-b1aee30e663a"),
                       ("QMLMakie", "08f9cac3-3b11-4f1c-9d88-d0e81c500f64"),
                       ("QML",      "2db162a6-7e43-52c3-8d84-290c1c42d82a"))

_pkgid(name, uuid) = Base.PkgId(Base.UUID(uuid), name)
_is_loaded(name, uuid) = haskey(Base.loaded_modules, _pkgid(name, uuid))

function _require_or_explain(name, uuid)
    try
        return Base.require(_pkgid(name, uuid))
    catch err
        error("""
            oitoolsgui() needs $name, and it is not available in the active environment.

            GLMakie, QMLMakie and QML are WEAK dependencies of OITOOLS, so `using OITOOLS`
            does not bring them in and they are not installed with it. Add them once:

                using Pkg; Pkg.add(["GLMakie", "QMLMakie", "QML"])

            Or use the launcher, whose environment already pins all three:

                julia --project=bin bin/oitoolsgui.jl [file.oifits]

            The original error was: $(sprint(showerror, err))
            """)
    end
end

"""
    oitoolsgui(files...; nautilus = true, verbose = true)

Open the GUI, loading everything it needs on the way.

Intended for the REPL: `using OITOOLS; oitoolsgui()`, optionally with OIFITS files to open in
the session. Equivalent to `bin/oitoolsgui.jl`, which stays the right entry point for a
launcher, a desktop icon or a sysimage.

GLMakie, QMLMakie and QML are weak dependencies, so they are not installed alongside OITOOLS —
`Pkg.add` them once into whatever environment you work in, or use `--project=bin`, which pins
known-good versions in its manifest.

!!! note "Call this before loading GLMakie yourself"
    The Mesa, GLFW and Qt platform hints are only read once, when the first OpenGL context is
    created and when Qt starts. This function sets them and then loads those packages, in that
    order. If GLMakie is already loaded the hints have been missed, and it says so rather than
    opening a window whose platform silently disagrees with Qt's.
"""
function oitoolsgui(files::AbstractString...; nautilus::Bool = true, verbose::Bool = true)
    if _is_loaded("GLMakie", "e9467ef8-e4e7-5192-8a1a-b1aee30e663a")
        @warn """
            GLMakie was already loaded, so the graphics hints could not be applied in time.
            On Wayland this means XWayland rather than native Wayland, and Qt may pick the
            other platform. Restart Julia and call oitoolsgui() before loading GLMakie.
            """
    end

    configure_graphics!(; verbose)                     # before the first GL context

    # GLFW_jll only dlopens libglfw, so this needs no display and creates no context.
    wl = try
        _require_or_explain("GLFW_jll", "0656b61e-2033-5cc2-a64a-77c0f6c09b89")
        # `invokelatest`, here and below: `Base.require` defines these methods AFTER this
        # function was compiled, so a direct call is a world-age error -- the method exists
        # and is unreachable from the world this frame is running in.
        Base.invokelatest(prefer_native_wayland!; verbose)   # before GLMakie
    catch err
        verbose && @debug "native Wayland not configured" err
        (; applied = false)
    end
    configure_qt_platform!(; match_x11 = !wl.applied, verbose)   # before Qt starts

    for (name, uuid) in _GUI_PACKAGES
        _require_or_explain(name, uuid)
    end
    # Cheap and pure Julia; without it "Nested sampling" is greyed out in the Model panel.
    nautilus && try
        _require_or_explain("Nautilus", "0c5b9d3e-7a41-4d2f-9e6c-8b1f4a2d5c73")
    catch err
        verbose && @debug "Nautilus not loaded; nested sampling stays disabled" err
    end

    ext = Base.get_extension(@__MODULE__, :OITOOLSGUIExt)
    ext === nothing && error("OITOOLSGUIExt did not load even though GLMakie, QMLMakie and " *
                             "QML are present — try `Base.retry_load_extensions()`")

    session = Base.invokelatest(ext.Session)
    for f in files
        isfile(f) ? Base.invokelatest(ext.load_dataset!, session, f;
                                      warn = false, verbose = false) :
                    @warn "not a file, skipping" f
    end
    return Base.invokelatest(gui, session)
end
