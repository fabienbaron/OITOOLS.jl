# Graphical interface

One window with five perspectives over a single session — Exploring, Observing, Modeling,
Imaging and Tools — so a dataset moves from planning to exploration to fitting to imaging without
being exported and read back.

```julia
using OITOOLS
oitoolsgui()                  # or oitoolsgui("mydata.oifits")
```

Left of the Dataset box, **?** opens About (the OITOOLS and Julia versions, the licence, author
and acknowledgements) and **⚙** the appearance settings.

[`oitoolsgui`](@ref) is the whole launch: it loads GLMakie, QMLMakie and QML itself, in the
order the graphics hints require, then opens the window. Call it *before* loading GLMakie
yourself — Mesa, GLFW and Qt each read their configuration exactly once, when the first OpenGL
context is created and when Qt starts, so a `using GLMakie` that has already run cannot be
taken back. It warns in that case rather than opening a window whose platform silently
disagrees with Qt's.

Any optional package you do not have is named once, not treated as a failure: `GLFW_jll`
(native Wayland rather than XWayland), `Nautilus` (enables "Nested sampling" in the Model
panel), `Pigeons` and `PairPlots`. Pass `optional = false` to skip that check.

### Why these are not installed for you

GLMakie, QMLMakie and QML are **weak dependencies**. A weak dependency is declared in
`Project.toml` under `[weakdeps]` rather than `[deps]`: it is not installed with the package
and not loaded by `using OITOOLS`, but when *you* load it, the matching entry under
`[extensions]` fires and the code that needs it becomes available. That is what keeps
`using OITOOLS` free of Makie, of Qt and of Python — which matters on a headless machine, in a
script that only reads OIFITS files, and alongside any other application that maps its own Qt.

The cost is that you install them once yourself:

```julia
using Pkg; Pkg.add(["GLMakie", "QMLMakie", "QML"])
```

Weak dependencies also cannot be loaded from the package's own environment, so `--project=.`
either fails to find them or picks them up from your default environment at whatever versions
happen to be there. From a clone, `bin/` is an ordinary environment where all three are real
dependencies with a pinned manifest, and it is the right entry point for a desktop icon or a
sysimage:

```
julia --project=bin bin/oitoolsgui.jl [file.oifits]
```

The console at the bottom echoes the equivalent OITOOLS call for every action, so anything below
can be read off the window and pasted into a script.

## Choosing observables

Two levels, for two questions.

**What the dataset holds** — Exploring › Filter › **Tables** › *observables*: V², T3amp, T3φ,
VISamp, VISφ, flux. A box the data hold nothing of is disabled, and says so when hovered.
Applying the filter makes a new dataset without the unticked observables, which every other
perspective, an export and the script then see.

A T3 amplitude and its closure phase are one point, as are a VIS amplitude and its phase.
Unticking one half keeps the points for the other, with the unticked half "not measured": it is
not counted, not fitted, not drawn, and is written as NaN on export. Unticking both halves leaves
the table out when the file is read. The script line is
`filter_data(d, set_data_filter(d; …); drop = [:t3amp])`.

Exploring › Filter › **Spectral** › *average every N channels* averages channels when the file is
read (`binning`, *block* or *pmoired*, *c*), as Tools' spectral binning does.

**What a run uses** — the Observables boxes of Modeling (v2, t3amp, t3phi, visamp, visφ, diffvis,
flux) and of Imaging (v2, t3amp, t3phi, diffvis). They choose among what the dataset holds, for one
fit or one reconstruction, and leave the dataset as it is: fitting V² and closure phases only is
a tick, not a new dataset. The header χ²ᵣ is over the ticked observables, as the fit is.

A box is enabled when the dataset holds that observable, so a dataset filtered without T3amp
shows Modeling's and Imaging's t3amp boxes disabled.

## Tools

Operations on files and datasets that belong to none of the other perspectives, in two frames.

**OIFITS Handling** — every operation works from the file in the frame's field (**…** to choose
one, **Current** for the selected dataset's file); no dataset needs to be open. A file chosen
with **…** (or typed and Enter) shows its headers.

- **Headers**: every HDU of the file, as in FITS Image Handling below.
- **Validate**: errors (the file cannot be read, or breaks the standard where nothing can be read
  around — a revision-2 table numbering targets or stations from 0, for one), warnings (what was
  read around, stations with no entry in OI_ARRAY, values with unusable errors, OI_VIS tables
  disagreeing on AMPTYP / PHITYP, 0-based numbering in revision-1 tables) and notes (flagged
  fractions, HDUs and targets an export does not keep). [`validate_oifits`](@ref).
- **Spectral binning**: every *N* channels, *average* or *group*, into a new dataset. *average*
  makes one channel of N (`binning = N`): *block* takes 1/σ²-weighted means, circular means for
  phases, errors multiplied by √c for a correlation *c* between adjacent channels, and reports a
  scatter check (χ²/dof about a straight line in λ within each group); *pmoired* reproduces
  PMOIRED's `binning`. *group* puts N consecutive channels in each spectral bin and keeps every
  channel (`spectralbin`): a polychromatic dataset, as Imaging reconstructs one image per bin.
- **Temporal binning**: one dataset per epoch, a new epoch after a gap longer than the given
  days, or fixed windows of the given days (`temporalbin`). **Preview** lists the epochs with
  their point counts; **also write each epoch to** a folder writes `<file>_<epoch>.oifits` for
  each, never over the source or an existing file.
- **Split**: one dataset per instrument (`insname`) or per target (`targetname`).
- **Merge files…**: several OIFITS files (Ctrl- or Shift-click) read as one dataset, stations
  matched by name, every point kept.
- **Export OIFITS…**: the file written as clean OIFITS v2 — every point it does not flag, each
  instrument under its own name, one target at a time (a **target** box appears for a file with
  several). [`repair_oifits`](@ref). Exploring's **Export** writes the current dataset instead,
  as it stands: filtered, binned, split or merged.

Each new dataset is the file read with one more `readoifits` keyword: it appears in the Dataset
box, the Filter works on it, and the exported script rebuilds it.

**FITS Image Handling** works on the file in its own field (**…** lists `*.fits` files); a file
chosen with **…** (or typed and Enter) is shown at once, as **View** does.

- **Headers** of any FITS file (images, cubes and the HDUs OITOOLS does not read included): the
  HDUs above, the selected one's columns (format, unit, dimensions) and every header card below,
  in order; the split between the two is movable. A file inside a ZIP archive is read through the
  archive. [`oifits_headers`](@ref).
- **View** shows the file's image or cube one channel at a time, East left, North up. Beside it,
  **Channel** steps through the wavelength channels (◀ ▶ or the slider), with each channel's
  wavelength (λ= … µm) when the file gives it; a cube with epochs as well gets an **Epoch** slider
  below, with each epoch's MJD (a time series alone steps its epochs on Channel). **Save PNG**
  writes the plane on screen.
- **Set keywords…** says what the file's axes are: the pixel size (mas, µas, arcsec, deg, rad),
  whether the file stores East on the left, whether axis 3 is wavelength or time (axis 4 being
  the other), and each axis as a start and a step (wavelength in m, µm, nm or Å; time as an MJD
  and a step in d, h, min or s). It opens with what the file's header says, and lists the
  file's WCS cards axis by axis (CTYPE, CUNIT, CRVAL, CDELT, CRPIX), ✓ present (the value on
  hover) or ✗ missing; Apply shows the cube
  with those values. **Save FITS…** then writes the file, to a new name or over itself, with the
  keywords [`read_model_image`](@ref) reads and the file's other header cards
  ([`fits_set_keywords`](@ref)).

The frame on the left and the results on the right are separated by a movable split.

## Building a model

A model is a dictionary of `"component,parameter"` keys. The Model perspective edits that
dictionary directly, and every rule below is the parser's, not the panel's — the same model
written by hand in a script behaves identically.

Each parameter is in exactly one of three states, chosen by the mode selector on its row:

| mode | what the row holds | in the fit vector? |
|:--|:--|:--|
| **fixed** | a number | no |
| **free** | a number (the initial guess), plus lower and upper bounds and its index `i` | yes |
| **expr** | an expression over other parameters | no — a free parameter must be numeric |

The index matters: `x[i]` and `list_free_params[i]` are the same parameter throughout OITOOLS,
so the free rows *in order* are the vector every optimiser sees and every `x_opt` reports
against.

### Setting a parameter derived

The mode selector is the small dropdown in the **mode** column of the parameter's own row —
the same control that makes a parameter free. Choose **derived** on it and the expression
editor opens; type the expression and press OK. The row then shows the expression, and its
value becomes a read-only display of what the resolver computes, because editing a derived
value would mean nothing.

The conversions are not symmetric, and the editor says so as it makes them:

- **fixed → derived** and **free → derived** both open the editor. Going derived also drops the
  parameter from the fit vector, and the dialog states how many free parameters will remain,
  so the count never changes behind you.
- **derived → free** needs a number to start from. The resolver has already evaluated the
  expression, so the current value is offered as the seed rather than asking you to retype a
  number you can see.

A parameter cannot be free and derived at once: the resolver raises on a string in
`list_free_params`, so the selector makes the combination unreachable rather than reporting it
after a fit.

## Limb-darkened disks: choosing the law

The kind list offers one **limb-darkened disk**, not one entry per law, because every one of
them is a disc of some diameter — which law it wears is a question about the atmosphere, not
about the geometry. The law is chosen beside the kind when the component is created, and can be
changed afterwards from the **Limb darkening** group on the selected component.

| law | profile | coefficients |
|:--|:--|:--|
| linear | `1 - u(1-μ)` | `u` |
| quadratic | `1 - u(1-μ) - w(1-μ)²` | `u`, `w` |
| square root | `1 - u(1-μ) - w(1-√μ)` | `u`, `w` |
| power law | `μ^α` (Hestroffer) | `alpha` |
| four-parameter | `1 - Σₖ cₖ(1 - μ^{k/2})` (Claret) | `c1`…`c4` |

Changing the law **keeps the diameter and reseeds the coefficients**. That asymmetry is
deliberate: the diameter is the same physical quantity under every law, and carries over along
with its bounds and whether it is being fitted. The coefficients are not. `u` is the linear
coefficient in both the quadratic and the square-root law, but the two laws are different
profiles, so a value fitted under one is not a value under the other — carrying it across would
make it look like a measurement that survived the change.

The same caution applies to published tables. Claret tabulates linear, quadratic, square-root,
logarithmic and four-parameter coefficients separately, and a power-law α is a different
quantity again — typically well below the linear `u` for the same star. Fit the law you intend
to quote.

## Making the fluxes behave as flux ratios

Component fluxes (`f`) are fractions of the total, and a fit only means what you think it means
if they sum to 1. There are three ways to arrange that, and they are not equally good.

### Two components: derive one from the other

```julia
"star,f" => 0.7,                 # free
"disk,f" => "1 - $star,f"        # expr
```

Set `disk,f` to **expr** and type `1 - $star,f`. The sum is exactly 1 for any value the fitter
tries, and the redundant degree of freedom is gone — two free fluxes plus a normalisation is one
parameter more than the model has information for.

### Three or more, or anything chromatic: normalise by a global

```julia
"star,f0" => 3.0,                       # free, unnormalised weights
"disk,f0" => 1.0,
"F"       => "$star,f0 + $disk,f0",     # a global: a bare key, no component
"star,f"  => "$star,f0 / $F",
"disk,f"  => "$disk,f0 / $F",
```

The fluxes are ratios by construction, stay positive, and no component is singled out to absorb
the remainder. Weights of 3:1 give 0.75 and 0.25; 1:1 gives 0.5 and 0.5. This is the form the
SPARCO models in `demos/` use, and the one that extends to wavelength: put the spectral
dependence in each numerator and the same terms in `F`, and the ratios stay normalised at every
wavelength.

Globals live in their own section above the components, because the parser groups bare keys
under `__global__`.

### SPARCO: a star and an environment with their own spectral slopes

**SPARCO fluxes…** (Model Creation, beside **Fluxes sum to 1**) writes this form for two
components: the star's flux a power law of fixed index (−4, Rayleigh–Jeans), the environment's
of a free slope `env,slope`, and the global `fs`, the star's share at λ0 (1.6 µm by default):

```julia
"star,f" => "$fs * ($WL/1.6)^(-4) / ($fs * ($WL/1.6)^(-4) + (1 - $fs) * ($WL/1.6)^$env,slope)"
"env,f"  => "(1 - $fs) * ($WL/1.6)^$env,slope / ($fs * ($WL/1.6)^(-4) + (1 - $fs) * ($WL/1.6)^$env,slope)"
```

`fs` and the slope are freed with bounds [0, 1] and [−4, 4]; other components keep their own
fluxes. The environment can be any component — a Gaussian, a ring, or an image
(`demos/models/SPARCO_Star_Image.toml` pairs a star with `SPARCO_env.fits`, an inclined ring;
`SPARCO_Star_Disk.toml` with a Gaussian disc). Such a model is simulated as any model, in
Observing › Source › Model.

### A constraint, only if you want the sum to be soft

In the Constraints panel: **+ constraint**, then lhs `$star,f + $disk,f`, op `==`, rhs `1`,
tol `0.001`.

Prefer one of the two forms above for an exact identity like this one. A constraint is a *hard*
nonlinear constraint only under the NLopt optimisers; under Levenberg-Marquardt and nested
sampling it becomes a one-sided quadratic penalty, and a penalty stiff enough to win against a
strong χ² gradient will wreck the conditioning while one that is not simply loses — the fit then
settles where the two gradients balance, which is neither the constrained answer nor the
unconstrained one. Constraints earn their keep on relations a box cannot express, such as
`$ring,diamout > $ring,diamin`.

## Keywords: `$R`, `$MU`, `$WL`, `$MJD`

Four names are reserved and mean something without being parameters. They are **scoped**, and
the expression editor lists all four with the out-of-scope pair greyed out and the reason
given.

| keyword | meaning | valid |
|:--|:--|:--|
| `$R` | radius over the profile grid, mas | only inside a `profile` string |
| `$MU` | `sqrt(1 - (r/r_max)^2)`, for limb darkening | only inside a `profile` string |
| `$WL` | wavelength in **microns**, one value per uv point | only outside a profile |
| `$MJD` | Modified Julian Date, one value per uv point | only outside a profile |

Crossing the boundary is not reported when it is typed. It surfaces later as
`UndefVarError: WL not defined in OITOOLS`, which names neither the rule nor the expression
that broke it — and a plain misspelling produces the same shape of error. The editor checks
each `$` reference as you type and distinguishes the three cases, which is the only place they
are told apart.

### Wavelength and time

Any ordinary expression may use them:

```julia
"disk,ud" => "3.0 * ($WL / 1.6)^1.2"
```

One consequence is worth knowing: if **any** expression in the model mentions `$WL` or `$MJD`,
the resolver broadcasts *every* parameter into a per-uv-point vector. That is a real change of
regime, and the inspector says when it has happened. It also says which wavelength the values
in the table are shown at, since a broadcast parameter has no single value.

`$MJD` resolves from `data.uv_mjd`, which is `Float64` whatever the dataset's precision, so
epochs minutes apart stay distinct.

## Radial profiles

A component whose geometry key is `profile` is a freeform I(r), compiled to an AD-transparent
closure and Hankel-transformed. Add one with **+ component → radial profile**.

```julia
"disk,profile" => "exp(-($R / $scale)^2 / 2)",
"disk,scale"   => 1.5,      # discovered from the expression
"disk,udout"   => 8.0,      # sets the r grid
"disk,nr"      => 64,       # optional: fixes the number of quadrature nodes
"disk,f"       => 1.0,
```

Two things about profiles differ from every other component.

**Parameters are discovered, not declared.** Every name in the expression that is not `$R` or
`$MU` is a parameter, resolved against this component first and then against the globals. Typing
`$scale` creates the *need* for a parameter rather than referring to one, so the editor lists
what is missing and offers to create it. A typo therefore shows up as a spurious extra
parameter instead of a crash.

**The grid comes from a key you may not have meant.** Four keys can set `r_max`, tried in this
order:

| key | interpretation |
|:--|:--|
| `udout` | a **diameter** — halved |
| `diamout` | a **diameter** — halved |
| `diam` | a **diameter** — halved |
| `r_max` | already a **radius** — used as is |

If none is present the model will not build, and the error names all four. The panel prints the
grid it actually got — "r ∈ [0.0, 4.0] mas from `udout`/2, 100 points" — naming the key that
won, because a grid twice the size you expected is otherwise a silent puzzle.

**Templates.** `demos/models/YSO_*.toml` are the six models of Ibrahim et al. (2023, ApJ 947, 68)
for HD 190073 (V1295 Aql, `demos/data/V1295_Aql.oifits`), at the paper's epoch-A values: a star
with a Rayleigh–Jeans spectrum, a disc and an over-resolved halo with their own spectral slopes,
the disc's profile being

| template | disc |
|:--|:--|
| `YSO_Doughnut` | PMOIRED's doughnut, `1 - (2*($R - ($RMIN + $RMAX)/2)/($RMAX - $RMIN))^2` on an annulus |
| `YSO_DoubleSigmoid` | `1/(1 + exp(-($R - $rin)/$sigin)) * 1/(1 + exp(($R - $rout)/$sigout))` |
| `YSO_AlphaDoubleSigmoid` | the same with a floor `$ain` inside the inner edge |
| `YSO_SkewedRing` | the double sigmoid with a first azimuthal mode (`az amp1`, `az projang1`) |
| `YSO_SkewedRing2` | and a second (`az amp2`, `az projang2`) |
| `YSO_OffCenterStar` | the skewed ring, the star free to move (`star,x`, `star,y`) |

Beside the profile the panel draws **I(r) and V(B) together**. The brightness distribution is
what you write; the visibility is what the data constrains, and two profiles that look alike in
I(r) can be obviously different in V(B).

## What the parser understood

The inspector reports what `dict_to_model` made of the model, rather than what was typed. It
exists for one failure in particular: a key the parser does not recognise is **ignored**, so a
misspelt geometry key silently fits a different model. Those keys are listed in red.

Validation sits at the top: values outside their bounds, `lb >= ub`, and numeric flux fractions
that do not sum to 1 — which is the check the first section of this page is about.

## Azimuthal modes

A radial-profile component can carry azimuthal brightness modes, added in the panel below the
profile editor. Each is a **pair** of keys, and the parser errors if one is present without the
other, so the panel adds and removes them together:

```julia
"disk,az amp1"     => 0.15,
"disk,az projang1" => 40.0,
```

Both also appear as ordinary rows in the parameter table, so a mode can be freed and fitted
like any other parameter. New modes are seeded at zero amplitude, which leaves the component
exactly as it was — adding the structure makes the asymmetry available to a fit without moving
one that is already running.

!!! note "Angle convention"
    `az projang` follows PMOIRED's convention — `cos(n·(ψ + φ − π/2))`, matching `_Vazvar` in
    its `oimodels.py` term for term — so an angle means the same thing in both packages and a
    model moves between them unchanged.

## Orbits and radial velocities

**Geometry › orbit** puts a component on an orbit around the origin (`c,x = c,y = "orbit"`): its
position at each data point's date comes from the elements, added as rows of the parameter table
(`c,orb P`, `MJD0`, `e`, `omega`, `OMEGA`, `incl`, `a`; PMOIRED's keys and definitions — see the
API page, *Orbits and radial velocities*). A PMOIRED dict with `'2,x': 'orbit'` imports as it is;
its `'additional residuals'` function is dropped, and a global named like `P-42.38` becomes
`P_42_38`.

**RV table…** (Model Optimization, beside the observables) opens radial velocities, a CSV
`MJD,RV,sigma,star` (`a` primary, `b` secondary), and adds the velocity elements (`Ka`, `Kb`,
`gamma`) to the component on an orbit. **RVs** then fits them with the interferometry
(Levenberg-Marquardt, NLopt, and the bootstrap, which resamples them); the console gives χ² for
the interferometry and for the velocities. `demos/data/HD210763/rv_gallenne2023.csv` (PMOIRED's
EX6) and `demos/data/Spica/rv_spica.csv` are examples.

**Orbit**, under the fits table, draws the relative orbit on the sky (East left) with the
positions at the dataset's nights, and the radial velocities against phase with the model's
curves.

## Multi-start search

**Multi-start…** (Model Optimization, beside **Companion search…**) fits the model from many
starting points and lists the distinct minima they reach — PMOIRED's `gridFit`
(`multistart_fit`). One row per free parameter says where its starts come from: **start** (its
value in the table), **grid** (min, max, step: every node; several grids make every
combination), **uniform** (min, max) or **normal** (mean, σ), drawn once per start. Without a
grid, **Starts** sets how many; **Seed** fixes the draws. **Conditions**, one per line
(`sqrt($B,x^2 + $B,y^2) <= 12`), drop starts before fitting; **Δχ²ᵣ ≤** keeps only minima that
close to the best. A PMOIRED `expl` dict pasted under **PMOIRED expl** fills the rows (**Apply**).

Every start is a Levenberg-Marquardt fit with the table's bounds, the constraints, the ticked
observables and, when ticked, the radial velocities. The fits are sorted as PMOIRED sorts them:
those that failed, have no uncertainties or barely moved are set aside, the rest grouped into
minima; for a component on an orbit, minima whose MJD0 differ by whole periods are one. The
global minimum becomes the fit's row; the console lists the first minima.

The **χ² map** view then shows the search over two of its parameters (the boxes under the
plot): each start (`+`) joined to the minimum it reached, the minima coloured by χ²ᵣ, set-aside
starts as red crosses, the selected minimum ringed. The third box selects the minimum the fit's
row, **Adopt** and the correlations report.

## Fitting each dataset

**Fit each…**, beside **Fit**, runs the same fit — model, free parameters and bounds,
constraints, priors, observables and optimiser — on each ticked dataset, one after another, each
from the values in the table. Datasets made by Tools' temporal binning are the usual input: one
fit per night. Levenberg-Marquardt or an NLopt method; uncertainty from the Jacobian or none. The
result is the **Epochs** table, one row per dataset: its mean MJD, χ²ᵣ, ndof and each free
parameter with its σ, with **Copy** (tab-separated) and **Save…** (`.tsv`).

## Diagnostics: residuals, χ² map and SED

Under the fits table on the **Model Optimization** tab, one plot area shows whichever of three
diagnostics you ask for.

**Residuals** draws `(model − data)/σ` against baseline, one colour per observable, with the
±1 and ±3 lines marked. It is the residuals only — not the ten-panel data-and-model stack that
[`plot_residuals`](@ref) draws for a figure — because the question here is whether they are
centred on zero, whether one observable is systematically worse, and whether there is structure
against baseline that a reduced χ² of 3.2 does not show.

It is the model **as it stands in the table**, so it works before a fit as well as after one:
whether a starting point is anywhere near the data is worth knowing before spending an
optimiser on it. After a fit, **Adopt** writes the fitted values into the table and the
residuals redraw with them. While the panel is showing, it follows every edit.

**χ² map** shows the surface a grid search produced. It does not run one — running a grid
search *is* the Grid search optimiser, and a button that quietly ran another would be a
different fit under the same name.

**SED** draws total flux and per-component flux against wavelength, over the loaded dataset's
own wavelengths. It is available only for a model that references `$WL`: a model that does not
is the same model at every wavelength, and its SED is a set of flat lines.

!!! warning "`model_to_sed` takes metres"
    `wl_grid` is in metres, as `uv_lam` in the data and every `wl` argument; only the
    expressions see the wavelength in microns, as `$WL`. A grid in microns evaluates every
    `$WL` expression a million times too close in and returns fluxes that are wrong without
    raising anything.

Each panel has its own **Save PNG**, which writes the plot as it is framed on screen.
