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
| PythonCall | Julia↔Python bridge (UltraNest, matplotlib) |
| PythonPlot | **optional** — matplotlib plotting, see [Plotting](@ref plotting-install) |
| Downloads | *(standard library)* SIMBAD queries, over its TAP service |
| Crayons | Coloured terminal output |

## Step 1: Python packages

Nothing to do. OITOOLS declares its Python dependencies in `CondaPkg.toml`, and
[PythonCall](https://github.com/JuliaPy/PythonCall.jl) provisions them automatically into a
project-local environment the first time you `using OITOOLS`:

| Package | Used by |
|---|---|
| matplotlib-base | the corner plot of a nested-sampling fit, and the plotting extension below |
| ultranest | `fit_model_nested`'s `:ultranest` backend — optional; `Nautilus.jl` needs no Python |
| scipy | imported at module level by `ultranest.plot` |

`matplotlib-base`, not `matplotlib`: on conda-forge the latter is a metapackage adding
PySide6, which pulls in Qt, libclang and libLLVM — about 510 MB unpacked for a GUI binding
nothing here imports. PythonPlot declares full `matplotlib` itself, so an interactive backend
arrives with it if you want one.

!!! tip "Reusing an existing Python"

    Set `JULIA_PYTHONCALL_EXE` to your interpreter to reuse an existing installation. Note
    that doing so opts out of `CondaPkg`, so you become responsible for keeping matplotlib
    and ultranest installed and mutually consistent.

## What each optimiser needs

`using OITOOLS` loads no Python and no plotting stack, so several engines live behind a package
you have to load yourself. Each one below says exactly what.

An engine whose package is missing is not silently substituted: the function has no methods and
the call raises a `MethodError` naming it, or — for the nested samplers, which share one
entry point — an error listing what to load.

### Model fitting

| Engine | Call | Load |
|---|---|---|
| NLopt, gradient and derivative-free | `fit_model(...; method = :LD_LBFGS)` and any other NLopt symbol | `using OITOOLS` |
| Levenberg–Marquardt | `fit_model_lsqfit` | `using OITOOLS` |
| Grid search over two parameters | `chi2_map` | `using OITOOLS` |
| Bootstrap resampling | `bootstrap_fit` | `using OITOOLS` |
| Nested sampling — pure Julia | `fit_model_nested` | `using Nautilus` |
| Nested sampling — UltraNest | `fit_model_nested(...; backend = :ultranest)`, or `fit_model_ultranest` | `using PythonCall`, plus the Python `ultranest` package CondaPkg installs |

`fit_model` takes any NLopt method symbol — there is no whitelist — and `use_grad` follows the
`LD_` prefix. `bootstrap_fit` wraps an inner fitter (`:lsqfit` or `:nlopt`), so it needs nothing
beyond what that fitter needs.

Both nested samplers fill the same `fit_model_nested`; `nested_backend()` reports which is in
force and `set_nested_backend!` chooses when both are loaded. Their `logz` should agree within
`logzerr`, which is the cheapest check that either has converged.

### Imaging

| Engine | Call | Load |
|---|---|---|
| VMLMB + regularisers | `reconstruct` | `using OITOOLS` |
| Maximum entropy (BSMEM) | `reconstruct_bsmem` | `using OITOOLS` |
| ADMM (BSDMM) | `reconstruct_bsdmm` | `using OITOOLS` |
| SPARCO, chromatic hybrid | `reconstruct_sparco`, `reconstruct_hybrid` | `using OITOOLS` |
| Simulated annealing (SQUEEZE) | `reconstruct_squeeze` | `using OITOOLS` |
| Parallel tempering (SQUEEZE) | `reconstruct_squeeze_tempered` | `using Pigeons` |

### Drawing what they produce

| Figure | Call | Load |
|---|---|---|
| Corner plot from any posterior | `plot_corner_makie` | `using PairPlots` and a Makie backend |
| Corner plot from an UltraNest result | `plot_ultranest_corner` | `using PythonPlot` |
| Live SQUEEZE monitor | `reconstruct_squeeze(...; monitor = n)` | `using PythonPlot` |
| Every other figure | see [Plotting](@ref) | `using PythonPlot` for the matplotlib set, any Makie backend for the `_makie` set |

### The short version

```julia
using OITOOLS                                    # NLopt, LM, grid search, bootstrap, all imaging but tempering
using OITOOLS, Nautilus, PairPlots               # + nested sampling and its corner plot, no Python
using OITOOLS, PythonPlot                        # + matplotlib figures, UltraNest, the SQUEEZE monitor
using OITOOLS, Pigeons                           # + parallel-tempered SQUEEZE
```

## Faster GUI startup: a sysimage

Most of the delay between launching the GUI and seeing a window is Julia compiling Makie, Qt
and the reconstruction driver. A [PackageCompiler](https://github.com/JuliaLang/PackageCompiler.jl)
sysimage holds that work already done:

```bash
julia --project=bin bin/compile.jl                  # tens of minutes, roughly 1 GB
julia --project=bin --sysimage bin/oitoolsgui.so bin/oitoolsgui.jl file.oifits
```

`bin/trace.jl` is what gets traced, and it is deliberately the work a cold start pays for:
loading the stack, building every figure `gui()` builds, drawing into each one, one model fit,
one χ² map, one nested-sampling run and two iterations of `reconstruct`.

!!! note "Wayland and the sysimage"

    GLFW.jl initialises from its own `__init__`, which in a sysimage runs before any user code,
    so `JULIA_GLFW_PLATFORM` arrives too late and GLMakie stays on X11. The launcher notices and
    brings Qt to xcb to match, so the two halves still agree — but a sysimage run on Wayland
    uses XWayland where the plain launcher uses Wayland. Export the variable before Julia starts
    to get the native backend:

    ```bash
    JULIA_GLFW_PLATFORM=wayland julia --project=bin --sysimage bin/oitoolsgui.so bin/oitoolsgui.jl
    ```

!!! warning "A stale sysimage runs the old code"

    The image contains the version of OITOOLS it was built from. Change the package, the
    Manifest or Julia itself and you must rebuild, or you will be running code that is no
    longer on disk — and nothing will tell you. This is the one failure mode of the approach
    worth being careful about.

**No Python is in the image.** Neither PythonPlot nor PythonCall is passed to `create_sysimage`
and `bin/trace.jl` loads neither, so a session started on the sysimage carries no conda
environment: nested sampling comes from Nautilus.jl and plotting from Makie. The
matplotlib figures and `plot_ultranest_corner` are simply absent, and calling one raises a
`MethodError` naming it. `using PythonPlot` in that session brings them back — the sysimage
does not prevent it, it just does not carry it.

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
configure_graphics!()                                   # before the first OpenGL context

using GLFW_jll
wl = prefer_native_wayland!()                           # before `using GLMakie`
configure_qt_platform!(; match_x11 = !wl.applied)       # before Qt starts

using GLMakie, QMLMakie, QML     # these activate the GUI extension
gui()
```

**The session decides the windowing system.** On Wayland both Qt and GLMakie use Wayland; on
X11 both use X11. Neither of the last two calls does anything in the ordinary case — Qt already
follows the session on its own, and the only reason `prefer_native_wayland!` is needed at all is
that GLFW.jl hard-codes `PLATFORM_X11` on Linux (`glfw3.jl:499`), so GLMakie would land on
XWayland under a Wayland compositor. It sets `JULIA_GLFW_PLATFORM=wayland`, which GLFW.jl reads
at its own `Init`.

What must not happen is one half on each: two EGL display connections and two surface lifetimes
in one process. Only GLMakie's choice can fail, so Qt is the one that yields — that is what
`match_x11 = !wl.applied` says. GLMakie's choice fails in two cases: `OITOOLSGUI_GLFW_X11` is
set, or GLFW is already initialised, which is what a sysimage does.

The order is not stylistic. Mesa reads its driver variables when the first OpenGL context is
created, and GLFW.jl reads `JULIA_GLFW_PLATFORM` at its `Init`, so both calls have to happen
above `using GLMakie`. That is also why neither lives in the GUI extension: loading GLMakie is
what creates that extension, so anything inside it would already be too late.
`configure_graphics!` is in the core package and needs nothing; `prefer_native_wayland!` needs
only `libglfw` and lives in a one-function extension that `GLFW_jll` alone activates.

`configure_graphics!` acts only on WSL, where there is no `/dev/dri` render node and Mesa would
otherwise fall through to software rendering.

!!! note "Qt used to be pinned to XWayland"

    Qt's Wayland backend does not release the windows it opens, and a
    `QtQuick.Dialogs.FileDialog` stayed on screen after closing — measured in
    `test/gui/filepicker_min.jl`, where leaving the dialog to manage itself, forcing
    `DontUseNativeDialog` and destroying the QML object outright all failed, and only
    `QT_QPA_PLATFORM=xcb` worked. The GUI draws its own picker now and no QML imports
    `QtQuick.Dialogs`, so the pin is gone. It cost nothing on the GPU either way — measured on
    WSL, `GL_RENDERER` is the same D3D12/NVIDIA adapter on both backends — but XWayland gives up
    application-side fractional scaling.

Any of it can be switched off:

```bash
OITOOLSGUI_NO_GPU_SETUP=1      # skip the Mesa setup
OITOOLSGUI_GLFW_X11=1          # keep GLMakie on X11; Qt follows it to xcb
QT_QPA_PLATFORM=wayland        # say so explicitly; an existing setting is never touched
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
