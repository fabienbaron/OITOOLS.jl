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

# Without these three there is no window at all.
const _GUI_PACKAGES = (("GLMakie",  "e9467ef8-e4e7-5192-8a1a-b1aee30e663a"),
                       ("QMLMakie", "08f9cac3-3b11-4f1c-9d88-d0e81c500f64"),
                       ("QML",      "2db162a6-7e43-52c3-8d84-290c1c42d82a"))

# Loaded when present, suggested once when not. The note is what the package does TODAY, which
# for two of them is less than their name suggests -- see the caveats in each line. The Python
# extensions (PythonPlot for matplotlib figures, PythonCall for UltraNest) are deliberately
# absent: they pull a Conda environment, which is a much larger thing to install than any of
# these, and the GUI works fully without them.
# The fourth field is how to install it. Nautilus is UNREGISTERED -- `Pkg.add("Nautilus")`
# cannot find it -- so it carries the same url `[sources]` uses in Project.toml, and a
# suggestion that told you to add it by name would simply fail.
const NAUTILUS_URL = "https://github.com/fabienbaron/Nautilus.jl.git"

const _GUI_OPTIONAL = (
    ("GLFW_jll",  "0656b61e-2033-5cc2-a64a-77c0f6c09b89",
     "native Wayland instead of XWayland (pure binary, no display needed to load)",
     :registered),
    ("Nautilus",  "0c5b9d3e-7a41-4d2f-9e6c-8b1f4a2d5c73",
     "enables \"Nested sampling\" in the Model panel, which is greyed out without it",
     :url),
    ("Pigeons",   "0eb8d820-af6a-4919-95ae-11206f830c31",
     "parallel tempering: reconstruct_squeeze_tempered, and the Image tab's \"Tempering\" " *
     "entry, which is greyed without it",
     :registered),
    # Unregistered AND without a remote, so neither `Pkg.add(name)` nor `Pkg.add(url=)` can
    # find it: it is developed from a local checkout. Loading it here is what ACTIVATES
    # OITOOLSVarInfExt -- an extension needs its trigger package loaded, not merely installed,
    # and the Image tab's "Variational inference" entry greys itself on exactly that check.
    ("VarInf",    "cd771066-68b7-4d72-889c-da7b55f139d5",
     "variational inference: the Image tab's \"Variational inference\" entry, which is " *
     "greyed without it",
     :path),
    ("PairPlots", "43a3c2be-4208-490b-832a-a21dcd55d7da",
     "plot_corner_makie for sampler posteriors, from the REPL. Not reachable from the GUI",
     :registered),
)

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
    _install_hint(name, how) -> String

The command that would install one optional package, in the form it actually needs.

Three kinds, because three kinds exist: registered (`Pkg.add` by name), a URL (Nautilus, which
is unregistered but has a remote), and a local path (VarInf, which has neither — it is
developed from a checkout). A suggestion that told you to `Pkg.add("VarInf")` would simply
fail, which is worse than saying nothing.
"""
function _install_hint(name, how)
    how === :registered && return "Pkg.add(\"$name\")"
    how === :url && return "Pkg.add(url = \"$(NAUTILUS_URL)\")"
    p = _sibling_checkout(name)
    return p === nothing ? "Pkg.develop(path = \"…/$name.jl\")   # a local checkout" :
                           "Pkg.develop(path = \"$p\")"
end

"""
    _sibling_checkout(name) -> String or nothing

A `name.jl` checkout beside this package, which is where an unregistered dependency developed
alongside OITOOLS normally sits. Returns `nothing` rather than guessing when it is not there.
"""
function _sibling_checkout(name)
    root = pkgdir(@__MODULE__)
    root === nothing && return nothing
    for c in (joinpath(dirname(root), name * ".jl"), joinpath(dirname(root), name))
        isfile(joinpath(c, "Project.toml")) && return c
    end
    return nothing
end

"""
    _offer_to_install(missing_pkgs; verbose) -> Bool

Ask whether to install the optional packages that are not there, and do it if told to.

**The default is no**, on Enter and on anything that is not a yes, because this is a question
asked while a window is opening and an unattended session must not stop on it. Only an
interactive terminal is asked at all: a script, a test run or a sysimage build gets the message
and carries on.

Returns whether anything was installed, so the caller knows to load the new packages.
"""
function _offer_to_install(missing_pkgs; verbose::Bool = true)
    lines = join(["    $name  --  $what" for (name, what, _) in missing_pkgs], "\n")
    cmds  = join(["    " * _install_hint(name, how) for (name, _, how) in missing_pkgs], "\n")
    interactive = isinteractive() && isa(stdin, Base.TTY)

    if !interactive
        verbose && @info """
            The GUI is running. These optional packages are not installed:

            $lines

            To add them:

                using Pkg
            $cmds

            Pass `optional = false` to skip this check entirely.
            """
        return false
    end

    printstyled("\nOptional packages that are not installed:\n"; bold = true)
    println(lines)
    print("\nInstall them now? [y/N] ")
    ans = try; lowercase(strip(readline())); catch; ""; end
    if !(ans in ("y", "yes"))
        verbose && @info "Skipped. To add them later:\n\n    using Pkg\n$cmds"
        return false
    end

    # Pkg is a dependency but is deliberately NOT loaded by `using OITOOLS` -- it is wanted
    # here, once, and nowhere else. `Base.require` is the same as a `using` at this point in
    # the file, and `invokelatest` is needed for the same reason as everywhere else in this
    # file: the methods arrive after this function was compiled.
    pkg = try
        Base.require(Base.PkgId(Base.UUID("44cfe95a-1eb2-52ea-b672-e2afdf69b78f"), "Pkg"))
    catch err
        @warn "could not load Pkg; install by hand" exception = err
        return false
    end

    ok = false
    for (name, _, how) in missing_pkgs
        try
            @info "installing $name"
            if how === :registered
                Base.invokelatest(pkg.add, name)
            elseif how === :url
                Base.invokelatest(pkg.add; url = NAUTILUS_URL)
            else
                p = _sibling_checkout(name)
                if p === nothing
                    @warn "$name has no registry entry and no checkout beside OITOOLS; " *
                          "install it by hand with Pkg.develop(path = \"…/$name.jl\")"
                else
                    Base.invokelatest(pkg.develop; path = p)
                end
            end
            ok = true
        catch err
            @warn "could not install $name" exception = err
        end
    end
    return ok
end

"""
    oitoolsgui(files...; optional = true, verbose = true)

Open the GUI, loading everything it needs on the way.

Intended for the REPL: `using OITOOLS; oitoolsgui()`, optionally with OIFITS files to open in
the session. Equivalent to `bin/oitoolsgui.jl`, which stays the right entry point for a
launcher, a desktop icon or a sysimage.

GLMakie, QMLMakie and QML are weak dependencies, so they are not installed alongside OITOOLS —
`Pkg.add` them once into whatever environment you work in, or use `--project=bin`, which pins
known-good versions in its manifest.

Four further packages are loaded when present and named in one message when not: `GLFW_jll`
(native Wayland), `Nautilus` (nested sampling in the Model panel), `Pigeons` (tempering) and
`PairPlots` (corner plots). None is needed for the window, and `optional = false` skips the
check. The Python extensions are deliberately not among them: matplotlib figures and UltraNest
need a Conda environment, which is a far larger install than any of these.

!!! note "Call this before loading GLMakie yourself"
    The Mesa, GLFW and Qt platform hints are only read once, when the first OpenGL context is
    created and when Qt starts. This function sets them and then loads those packages, in that
    order. If GLMakie is already loaded the hints have been missed, and it says so rather than
    opening a window whose platform silently disagrees with Qt's.
"""
function oitoolsgui(files::AbstractString...; optional::Bool = true, verbose::Bool = true)
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

    # The optional ones are never fatal: one message naming all of them, once, beats four
    # separate failures for things the GUI runs perfectly well without.
    if optional
        missing_pkgs = Tuple{String,String,Symbol}[]
        for (name, uuid, what, how) in _GUI_OPTIONAL
            try
                Base.require(_pkgid(name, uuid))
            catch
                push!(missing_pkgs, (name, what, how))
            end
        end
        if !isempty(missing_pkgs) && _offer_to_install(missing_pkgs; verbose)
            # Installing brings the trigger packages into the environment, but an extension
            # only activates once its trigger is LOADED. Both, then, and in that order.
            for (name, _, _) in missing_pkgs
                i = findfirst(t -> t[1] == name, _GUI_OPTIONAL)
                i === nothing && continue
                try; Base.require(_pkgid(name, _GUI_OPTIONAL[i][2])); catch; end
            end
            try; Base.retry_load_extensions(); catch; end
        end
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
