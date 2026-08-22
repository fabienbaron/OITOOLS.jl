# Installation

## Dependencies

| Package | Purpose |
|---------|---------|
| [OIFITS.jl](https://github.com/emmt/OIFITS.jl) | Read OIFITS v1/v2 files |
| NFFT | Non-uniform fast Fourier transform for imaging |
| OptimPackNextGen | VMLMB gradient optimizer for image reconstruction (Éric Thiébaut) |
| NearestNeighbors | KD-tree UV deduplication |
| SparseArrays | OI_CORR correlation matrix storage |
| NLopt | Local and global optimizers (Nelder–Mead, ISRES, …) |
| LsqFit | Levenberg–Marquardt least-squares fitting |
| UltraNest | Bayesian model selection via nested sampling |
| SpecialFunctions | Bessel functions for visibility calculations |
| PythonCall | Julia↔Python bridge (SIMBAD queries, UltraNest) |
| PythonPlot | **optional** — matplotlib plotting, see [Plotting](@ref plotting-install) |
| Crayons | Coloured terminal output |

## Step 1: Python packages

Nothing to do. OITOOLS declares its Python dependencies in `CondaPkg.toml`, and
[PythonCall](https://github.com/JuliaPy/PythonCall.jl) provisions them automatically into a
project-local environment the first time you `using OITOOLS`:

| Package | Used by |
|---|---|
| matplotlib-base | `fit_model_ultranest`'s corner plot, and the plotting extension below |
| astroquery | `magnitudes_from_simbad`, `ra_dec_from_simbad` |
| ultranest | `fit_model_ultranest` |
| scipy | imported at module level by `ultranest.plot` |

`matplotlib-base`, not `matplotlib`: on conda-forge the latter is a metapackage adding
PySide6, which pulls in Qt, libclang and libLLVM — about 510 MB unpacked for a GUI binding
nothing here imports. PythonPlot declares full `matplotlib` itself, so an interactive backend
arrives with it if you want one.

!!! tip "Reusing an existing Python"

    Set `JULIA_PYTHONCALL_EXE` to your interpreter to reuse an existing installation. Note
    that doing so opts out of `CondaPkg`, so you become responsible for keeping matplotlib,
    astroquery and ultranest installed and mutually consistent.

## Step 2: Julia packages

Add the EmmtRegistry (required for OptimPackNextGen, OIFITS, and related packages),
then install all dependencies:

```julia
using Pkg
pkg"registry add General"
pkg"registry add https://github.com/emmt/EmmtRegistry"
Pkg.add([
    "AstroTime", "CFITSIO", "Crayons", "Dates", "DelimitedFiles",
    "FFTW", "FITSIO", "Glob", "LaTeXStrings", "LinearAlgebra",
    "LsqFit", "Match", "NFFT", "NLopt", "NearestNeighbors",
    "Parameters", "Random", "SparseArrays",
    "SpecialFunctions", "Statistics", "UltraNest",
    "ArrayTools", "LazyAlgebra", "OptimPackNextGen", "OIFITS",
])
Pkg.add(url="https://github.com/fabienbaron/OITOOLS.jl.git")
```

Installation typically takes 2–10 minutes depending on your system.

Verify the installation:

```julia
using OITOOLS
```

## [Step 3: Plotting](@id plotting-install)

**The plotting functions live in a package extension and require PythonPlot.** `using OITOOLS`
alone gives you `uvplot`, `plot_v2`, `imdisp` and the rest as functions with no methods, and
calling one raises a `MethodError`.

```julia
using Pkg
Pkg.add("PythonPlot")

using OITOOLS, PythonPlot     # both: PythonPlot is what activates the extension
uvplot(data)
```

The extension (`OITOOLSPythonPlotExt`) loads automatically once PythonPlot is present; there
is nothing else to enable. It supplies every plotting entry point:

`uvplot`, `plot_v2`, `plot_t3phi`, `plot_t3amp`, `plot_visamp`, `plot_visphi`, `plot_diffphi`,
`plot_flux`, `plot_obs`, `plot_multi`, the `*_residuals` family, `plot_residuals`,
`plot_v2_multifile`, `plot_facility`, `imdisp`, `imdisp_multi`, `plot_ultranest_corner`,
`set_oiplot_defaults`, and the observation planners `gantt_onenight`, `obs_plan`,
`chara_plan`, `empty_night`.

Keeping them out of the core package means `using OITOOLS` starts no Python plotting stack:
it does not import matplotlib, and it maps no Qt into the process — which matters on a
headless machine, in CI, and beside any other Qt application.

To pick a non-interactive backend, for scripts or CI, set `MPLBACKEND` before loading
PythonPlot:

```julia
ENV["MPLBACKEND"] = "Agg"
using OITOOLS, PythonPlot
```

Observable metadata — `OBS_PLOT_SPECS`, `oiplot_colors`, `canonical_color` — stays in the core
package, so a front-end that draws with something other than matplotlib can share the axis
labels, groupings and palette without loading Python at all.


## [Step 4: The GUI (optional)](@id gui-install)

The GUI is a separate stack — GLMakie for drawing, Qt/QML for the window — and none of it is
loaded by `using OITOOLS`. It costs about 6.4 s of load time against the package's own 1.7 s,
so it is gated behind its own dependencies and no script, reduction or CI run pays for it.

From a clone of the repository:

```bash
julia --project=bin -e 'using Pkg; Pkg.instantiate()'   # first time only
julia --project=bin bin/oitoolsgui.jl                   # every time
```

An OIFITS file given on the command line is opened at startup:

```bash
julia --project=bin bin/oitoolsgui.jl demos/data/2004-data1.oifits
```

Three things about that command are deliberate.

**`--project=bin`, not `--project=.`.** GLMakie, QMLMakie and QML are *weak* dependencies of
OITOOLS, and weak dependencies are not loadable from the package's own environment.
`bin/Project.toml` declares them as ordinary ones. Under `--project=.` they either fail to load
or, worse, resolve from your default environment at whatever versions happen to be there.

**Run it from the repository root.** `bin/Project.toml` reaches the package with
`[sources] OITOOLS = {path = ".."}`, which is relative to `bin/`.

**`bin/Manifest.toml` is not in the repository.** The first `instantiate` resolves the GL and
Qt stack itself, within the bounds `Project.toml` declares (`GLMakie 0.13`, `Makie 0.24`,
`QML 0.13`, `QMLMakie 0.3`, `GLFW_jll 3.4`). Keep the resulting manifest if you want the same
versions back later.

### Starting it from your own session

```julia
using OITOOLS
configure_graphics!()            # before the first OpenGL context exists
configure_qt_platform!()         # before Qt starts

using GLFW_jll
prefer_native_wayland!()         # before `using GLMakie`

using GLMakie, QMLMakie, QML     # these activate the GUI extension
gui()
```

The order is not stylistic. Mesa reads its driver variables, and GLFW its platform hint, when
the first OpenGL context is created — so both calls have to happen above `using GLMakie`.
That is also why neither function lives in the GUI extension: loading GLMakie is what creates
that extension, so anything inside it would already be too late. `configure_graphics!` is in
the core package and needs nothing; `prefer_native_wayland!` needs only `libglfw` and lives in
a one-function extension that `GLFW_jll` alone activates.

All three are no-ops where they do not apply. `configure_graphics!` acts only on WSL, where
there is no `/dev/dri` render node and Mesa would otherwise fall through to software rendering.
`configure_qt_platform!` acts only in a Wayland session.

The Wayland pair is worth reading together, because the second decides what the third does.
Qt's Wayland backend does not release the windows it opens: a `QtQuick.Dialogs.FileDialog`
closes as far as QML is concerned — `visibleChanged -> false` arrives before `accepted`, and
`visible` stays false — while its window remains on screen. Measured with
`test/gui/filepicker_min.jl`, one mode per run: leaving the dialog to manage itself, forcing
`DontUseNativeDialog`, and destroying the QML object outright all leave the window up, and the
same run under `QT_QPA_PLATFORM=xcb` closes it properly. So `configure_qt_platform!` pins xcb
on a Wayland session.

That also settles the third call. `prefer_native_wayland!` exists because GLFW.jl hard-codes
X11, which would put GLMakie on XWayland while Qt ran native Wayland — two windowing systems in
one process. With Qt pinned to xcb, GLFW's X11 default already matches, so the function stands
down of its own accord. The cost of the pin is XWayland's: fractional HiDPI scaling is done by
the compositor rather than the application. On a 1:1 display it does not show.

Any of them can be switched off:

```bash
OITOOLSGUI_NO_GPU_SETUP=1      # skip the Mesa setup
OITOOLSGUI_NO_QT_PLATFORM=1    # stay on Qt's Wayland backend
QT_QPA_PLATFORM=wayland        # ...or say so explicitly; an existing setting is never touched
OITOOLSGUI_GLFW_X11=1          # keep GLMakie on XWayland
QSG_INFO=1                     # make Qt report which GL renderer it actually got
```

```@docs
gui
configure_graphics!
configure_qt_platform!
is_wsl
```

The GUI's own types (`Session`, `LiveCanvas`, `ShellState`) are defined inside the extension,
since a function can be declared in the core package and given methods later but a type cannot.
`gui()` builds its own `Session`, so most callers never need to name one; scripts that do reach
them through the extension module:

```julia
const GUI = Base.get_extension(OITOOLS, :OITOOLSGUIExt)
using .GUI
```
