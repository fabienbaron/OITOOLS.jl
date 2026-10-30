# Graphical interface

One window with five perspectives over one session: Exploring, Observing, Modeling, Imaging and
Tools.

```julia
using OITOOLS
oitoolsgui()                  # or oitoolsgui("mydata.oifits")
```

[`oitoolsgui`](@ref) loads GLMakie, QMLMakie and QML and opens the window. Call it before
`using GLMakie`: it sets graphics options that have no effect once GLMakie is loaded, and warns
if it is. Missing optional packages are listed once: `GLFW_jll` (native Wayland), `Nautilus`
(nested sampling), `Pigeons` (parallel tempering), `PairPlots` (corner plots). `optional = false`
skips the check.

GLMakie, QMLMakie and QML are weak dependencies of OITOOLS and are installed separately:

```julia
using Pkg; Pkg.add(["GLMakie", "QMLMakie", "QML"])
```

From a clone, `bin/` is an environment with these and the optional packages pinned:

```
julia --project=bin bin/oitoolsgui.jl [file.oifits]
```

**?** (left of the Dataset box) shows the versions, licence and acknowledgements; **⚙** the
settings: UI scale and font, the start-up tab, window and Qt scale, plot text, symbols and zoom,
the Observing time axis and elevation limits, and the work directory. **Save config** keeps them
for the next launch; **Reset to defaults** deletes the saved file. The console at the bottom
shows the OITOOLS call for each action.

## Choosing observables

- **Exploring › Filter › Tables › observables** (V², T3amp, T3φ, VISamp, VISφ, flux) removes
  observables from the dataset: applying the filter creates a new dataset without them. Boxes
  for observables the data do not hold are disabled. Removing one half of a T3 or VIS pair (an
  amplitude or its phase) keeps the other; the removed half is not counted, fitted or drawn,
  and is written as NaN on export (`filter_data(…; drop = [:t3amp])`).
- **Exploring › Filter › Spectral › average every N channels** averages channels on reading
  (`binning`; *block* or *pmoired*, *c*).
- The **Observables** boxes of Modeling (v2, t3amp, t3phi, visamp, visφ, diffvis, flux) and
  Imaging (v2, t3amp, t3phi, diffvis) choose which of the dataset's observables a fit or a
  reconstruction uses; the dataset is unchanged. The χ²ᵣ in the Modeling header uses the same
  selection.

## Targets (Observing)

**Enter** in a target's name resolves it on SIMBAD, as **SIMBAD** does, unless the row already
holds SIMBAD's answer for that name (its dot is green). **Save…** writes the list as CSV
(`name, ra, dec, V, J, H, K, simbad, cal, for`, coordinates in degrees, an empty magnitude for
one not known, `simbad` 1 for a row SIMBAD resolved, `cal` 1 for a calibrator, `for` the science
targets a calibrator is for, `;` between them); **Load…** reads it
back. A file with a header is read by column name, in any order, and its other columns are kept
and written back on Save (a target's model among them, below); one without is read in that order, and a file of names alone, one per
line, loads as unresolved targets. An emptied magnitude field is a magnitude not known. **Cal**
marks a target as a calibrator: the Gantt draws its bars in orange, and a field beside it takes
the names of the science targets it is for, on a line under it.
It opens on `demos/targets`, whose `bright_stars.csv` lists twelve bright stars across the year,
resolved on SIMBAD.

## Simulating targets (Observing)

**Model** on a target's row sets what that target is simulated from, picked from one menu: none,
the models held in the Modeling tab by name, a uniform disc of a given diameter, a FITS image
(with its pixel size, taken from the header when it has one, and the image's channels), or a
model file. The mark beside the button is ✓ for a target with a model, ○ for one without, and ✗
for a model that cannot be read (a model since removed, a missing file); the model is named on a
line under the target's name. Renaming a model in Modeling renames the reference. The CSV keeps
the model in the columns `sim_kind` (`session`, `modelfile`, `image`, `ud`), `sim_path`,
`sim_pixsize`, `sim_channels` and `sim_diam`.

**Simulate** observes every target marked ✓; the count is shown beside it. One target is
written to the output file. Several are written one file per target, `<output>_<name>.oifits`,
or, with **one file per target** unticked, all to the output file, its OI_TARGET listing every
target ([`join_targets`](@ref)). Each target takes its magnitude in the setup's band, or the
panel's **mag** when it has none; a magnitude typed in **mag** applies to every target. The k-th
target is seeded `seed + k − 1`. A calibrator is written with `CATEGORY = "CAL"` in `OI_TARGET`.
**Open in Explore** opens what was written, one dataset per target, each labelled with its
target's name.

**epochs from the plan**: a target the Gantt shows is observed at the times it is observable on
the chart, every **Step** minutes, instead of over the HA range.

**Use existing OIFITS as template** observes each target through the loaded dataset's uv points
and signal-to-noise instead, one file per target.

## POPs (Observing)

**Reference** is the delay reference: its cart stays put and the others follow it, so every POP
search, Gantt chart and delay chart checks the baselines to it. There is always one, W2 by
default; unticking it makes the most central telescope left the reference (the one whose
largest distance to the others, projected on the equator, is smallest).

**AutoPOPs** searches the POPs each time **Compute** runs, over the hours the target is above
the elevation limits and the horizon, in astronomical night, clear of the moon. The menu under
it chooses the search, in two sections, single target and multiple targets:

- **Single POPs**: one setting for the night (`best_pop`).
- **POP changes on a single line**: up to **changes ≤** POP changes, each taking **change**
  minutes (`pop_plan`).
- **POP changes & Telescope Drops on a single line**: the same, dropping up to **drop ≤**
  telescopes where no setting keeps them all in delay. The reference telescope is dropped only
  when no other drop gains the same time.
- **Three-row POP (time coverage)**: for the earliest, middle and latest
  third of the observable hours, the single setting with the most time in delay over it, one
  row each (a setting best in several thirds once). **Show details** and the delay carts show
  the earliest.
- **POP changes on a single line, all targets**: one sequence of POP settings for the whole
  list, up to **changes ≤** changes (`pop_plan_targets`). It maximises the sum, over the science
  targets, of the fraction of each one's observable hours kept, so a target up for an hour
  counts as much as one up all night. A science target with calibrators counts at a minute
  only when one of its own is in delay under the same POPs within **bracket ±** minutes;
  unticking **needs cal** scores the science targets alone. A calibrator never up or never in
  delay is left out (named in the console). The Gantt has a row per target, the calibrators in
  orange, and a POPs row on top with each setting and a line through every row at each change;
  clicking a target in the list picks the one **Show details** and the delay carts show.
  Without **AutoPOPs** every target is charted under the POPs set by hand.

Plans are ranked by the time with every selected telescope, then the time with one dropped
(the reference last), then with two, then by the fewest changes, and last by the fewest POPs
moved at them: among settings giving the same times, those keeping more POPs where they were.
**Three-row POP (time coverage)** picks its tied settings the same way.

On the Gantt a multi-POP plan is one bar per segment, lighter where a telescope is dropped. The
first POPs are written at the start, each change as **→** and the new POPs where it starts
(those that change in red), with **POP** in the bar; a telescope not used is a **×** in the
POPs. **DROP** marks the end of the segment before a telescope is lost, **−E1** the point it is
lost and **+E1** where it comes back.
**EL** marks the time above the soft elevation limit. **Show details** adds a POPs row and a
line through every row at each change.

**POP strategies** has the same strategy menu, a target menu (the target selected in the
list) and, for the strategies with changes, the minutes a POP change takes (the POPs panel's
**change**). A change of reference takes no time. Under them: **min stretch**, the shortest
time a POP setting is held, and a telescope dropped or added back (15 minutes; a shorter drop
takes its neighbour's telescopes where they are in delay); **night**, the part of the night
used, from and to in UT (`H` or `H:MM`, empty for the start or the end), which narrows the
astronomical night everywhere, the Gantt included; and for the list, **bracket ±** and
**needs cal**. For a single target's changes, **+ … min/POP** and **cart … m/s** set how long a
change takes: the longer of the POP change time plus the minutes per POP moved, and the delay
carts' travel to their new positions at that speed (empty, no time). **Precise** lists the plans
searched with them, about a second's search; **Search POPs**, **Compute** and AutoPOPs use the
POP change time alone. **Search POPs** lists that target's candidates: for a multi-POP strategy the ranked
plans (the plan, its changes, the telescopes dropped, and the hours with every telescope and in
all), for a single-setting one the best settings with their hour-angle and time ranges. The
times in the tables are on the clock chosen under **Time** when the search ran. Under the
tables a Gantt chart draws every candidate as a row, numbered as in the table. **Use** charts
the target with a candidate in place of AutoPOPs; **in use** gives it back. For the list, the
table holds the shared plans (the POP sequence, its changes, its score, and the hours of the
science targets and of all), and **Use** charts the list with one.

**Timeline**, under the tables, gives the night of the row clicked, or of the chart after
**Compute**, one row per step (`pop_plan_timeline`): each stretch observed with its POPs
(**×** for a telescope not used), its hour angles and duration; each POP change with the
telescopes whose POP moves; each telescope dropped or added; and the dark hours with nothing
observed, with the first constraint that excludes them: elevation-restricted, too close to the
Moon, or out of delay. Times follow the Gantt's clock. **Use** charts the target with the
selected candidate; **Copy** copies the timeline as text.

## Tools

Two frames; the results are shown on the right.

**OIFITS Handling** works on the file in its field: **…** to choose one, **Current** for the
selected dataset's file, or type a path and press Enter. No dataset needs to be open.

- **Headers**: choosing a file lists its HDUs. **Header** shows the selected HDU's cards;
  **Data** its contents ([`oifits_table`](@ref)):
  - OI_VIS2, OI_VIS, OI_T3 and OI_FLUX are shown **per channel** by default: one line per row
    and channel, preceded by λ, the station names, the baseline B (for T3 the longest of the
    three) and B/λ. Unticked, one line per row; array columns show the chosen **channel**.
  - Click a column header to sort (again for descending). The filter takes station names
    (`S1-E1`) and ranges (`λ 2.16–2.17`, `B 100–200`, `B/λ a–b`, `mjd a–b`). Flagged lines are
    grey; **hide flagged** removes them.
- **Validate** ([`validate_oifits`](@ref)): errors (unreadable file, standard violations such as
  0-based indices in a revision-2 table), warnings (what was read around, stations missing from
  OI_ARRAY, unusable errors, OI_VIS tables with different AMPTYP/PHITYP) and notes (flagged
  fractions, HDUs and targets an export drops).
- **Spectral binning**, into a new dataset: *average* N channels into one (`binning = N`;
  *block*: 1/σ²-weighted means, circular means for phases, errors × √c for an inter-channel
  correlation c, with a χ²/dof scatter check per group; *pmoired*: PMOIRED's `binning`), or
  *group* N channels per spectral bin, keeping all channels (`spectralbin`, a polychromatic
  dataset for Imaging).
- **Temporal binning** (`temporalbin`): one dataset per epoch, split at gaps longer than the
  given days or in fixed windows. **Preview** lists the epochs and their point counts; **also
  write each epoch to** writes `<file>_<epoch>.oifits` files to a folder, never overwriting.
- **Split**: one dataset per instrument (`insname`) or per target (`targetname`).
- **Merge files…**: several files (Ctrl- or Shift-click) read as one dataset, stations matched
  by name.
- **Export OIFITS…** ([`repair_oifits`](@ref)): the file as OIFITS v2, unflagged points only,
  one target at a time (a **target** box appears when there are several). Exploring's
  **Export** writes the current dataset as filtered.

New datasets appear in the Dataset box and are rebuilt by the exported script.

**FITS Image Handling** works on the file in its own field; choosing a file shows it.

- **Headers** of any FITS file, images included: the HDUs, the selected one's columns (format,
  unit, dimensions) and its cards. ZIP archives are read through. [`oifits_headers`](@ref).
- **View** shows the image or cube one plane at a time, East left: **Channel** steps through
  wavelengths (with λ when the file gives it), **Epoch** through epochs. **Side by side** shows up
  to six channels (numbers or λ in µm) on one colour scale; **power** sets a brightness stretch
  I^p. **Save PNG** writes the view.
- **Set keywords…** sets the pixel size (mas, µas, arcsec, deg, rad), whether East is on the
  left, the reference pixel, whether axis 3 is wavelength or time, and each axis's start and step.
  It opens with the file's values and lists its WCS cards (present ✓ / missing ✗). **Apply**
  shows the cube with those values; **Save FITS…** writes the file with the keywords
  [`read_model_image`](@ref) reads, keeping its other cards ([`fits_set_keywords`](@ref)).

## Building a model

The Modeling perspective edits the model dictionary (`"component,parameter"` keys); a model
written in a script behaves the same.

The session holds several models. The **Model** box at the top lists them and selects one; typing
a name there and pressing Enter renames the selected model (a name another model has is
refused). A dot beside the box marks unsaved edits. **+** adds an empty `model k`, **Copy**
duplicates the selected model (parameters, free list, bounds, constraints and priors), **−**
removes it (asking first when it has unsaved edits). **Open…** and **Import PMOIRED…** add the
model to the list and select it, or fill the selected model when it is empty and unedited; the
model is named after its file. **Save…** names a `model k` after the file and keeps a name given
by hand. Each model keeps its own constraints and priors. New, Open and Import keep the
bandwidth smearing in force when the model does not set it.

Each parameter is in one of three modes, set in the **mode** column of its row:

| mode | the row holds | fitted |
|:--|:--|:--|
| **fixed** | a number | no |
| **free** | a starting value, lower and upper bounds, and its index `i` in the fit vector | yes |
| **derived** | an expression over other parameters | no |

Free rows, in order, are the fit vector: `x[i]` is `list_free_params[i]`.

Choosing **derived** opens the expression editor; the row then shows the expression and its
computed value, read-only. Going from free to derived reduces the number of free parameters,
and the editor shows the new count. Going from derived to free starts from the computed value.

## Limb-darkened disks

The component kind **limb-darkened disk** takes a law, chosen when the component is created or
later under **Limb darkening**:

| law | profile | coefficients |
|:--|:--|:--|
| linear | `1 - u(1-μ)` | `u` |
| quadratic | `1 - u(1-μ) - w(1-μ)²` | `u`, `w` |
| square root | `1 - u(1-μ) - w(1-√μ)` | `u`, `w` |
| power law | `μ^α` (Hestroffer) | `alpha` |
| four-parameter | `1 - Σₖ cₖ(1 - μ^{k/2})` (Claret) | `c1`…`c4` |
| free-form | I(μ) through N nodes, `I(1) = 1` | `ld1`…`ld<N-1>` |

Changing the law keeps the diameter (value, bounds, free or fixed) and resets the coefficients
to the new law's defaults. Coefficients are not comparable between laws: `u` of the quadratic and
of the square-root law differ, and tabulated coefficients (e.g. Claret's) are per law.

**Free-form** starts from the profile the disc had: its nodes are that law's I(μ) at
μ = 0, 1/(N-1), …, so a fit resumes where the law's fit ended. Under the law menu:

- **nodes**: N, 5 by default, evenly spaced in μ. Changing it re-samples the current profile.
- **spline** or **straight lines** between the nodes (`ldspline` or `ldpiecewise`).
- **monotonic**: the free parameters become the ratios `ldq<k>` = I(μ_k)/I(μ_k+1) in [0, 1],
  and the nodes are their products, so I(μ) does not rise toward the limb at the nodes.

Visibilities measure few profile shapes: about one with data in the first lobe, two to three
into the second, four to five past the third null. More free nodes than that leave some
undetermined; **monotonic** or fewer nodes keep the fit meaningful.

**Free-form profile** (Model Optimization, shown when the model has a disc; needs VarInf) runs
`fit_ld_profile_vi` on the disc chosen there: I(μ) on a grid of nodes with a Gaussian-process
prior about the disc's current law, and a band from the posterior samples.

- **nodes**: the grid (24 by default).
- **GP length** (in μ) and **σ** (in log I): the Matérn-3/2 kernel of the prior.
- **geoVI** or **MGVI**, **samples**: antithetic pairs drawn for the band.
- **scale errors**: multiply the error bars by √χ²r when the profile's best fit has χ²r > 1.
- **bootstrap**: bagged VI, the run repeated on this many bootstrap replicates of the data
  (configurations drawn with replacement), their samples pooled; 0 for a single run.
- **Run** uses V² and T3amp of the ticked observables; closure phases are left out.

The result is the **Profile** view under the fits table: the band (16–84 %), the posterior
mean and the starting law (dashed), normalised to I(1) = 1 or ⟨I⟩ = 1. **Adopt** makes the disc
this profile (a straight-line free-form disc on the grid, nodes fixed); **Save profile…**
writes μ, r/R and both normalisations with their bands as CSV.

## Flux fractions

Component fluxes `f` are fractions of the total and should sum to 1.

**Two components**: derive one from the other.

```julia
"star,f" => 0.7,                 # free
"disk,f" => "1 - $star,f"        # derived
```

**Three or more, or chromatic fluxes**: normalise by a global (a key without a component).

```julia
"star,f0" => 3.0,                       # free weights
"disk,f0" => 1.0,
"F"       => "$star,f0 + $disk,f0",     # global
"star,f"  => "$star,f0 / $F",
"disk,f"  => "$disk,f0 / $F",
```

With chromatic weights in the numerators and in `F`, the fractions sum to 1 at every
wavelength. Globals are listed above the components.

**SPARCO fluxes…** (Model Creation) writes this form for a star with a fixed power-law index (−4,
Rayleigh–Jeans) and an environment with a free slope `env,slope`, `fs` being the star's fraction
at λ0 (1.6 µm by default):

```julia
"star,f" => "$fs * ($WL/1.6)^(-4) / ($fs * ($WL/1.6)^(-4) + (1 - $fs) * ($WL/1.6)^$env,slope)"
"env,f"  => "(1 - $fs) * ($WL/1.6)^$env,slope / ($fs * ($WL/1.6)^(-4) + (1 - $fs) * ($WL/1.6)^$env,slope)"
```

`fs` and the slope are freed with bounds [0, 1] and [−4, 4]. The environment can be any
component, an image included (`demos/models/SPARCO_Star_Image.toml`, `SPARCO_Star_Disk.toml`).

**A constraint** (**Fluxes sum to 1** adds `$star,f + $disk,f + … == 1` to the Constraints panel)
also works, but is enforced exactly only by the NLopt optimisers; Levenberg-Marquardt and the
samplers apply it as a quadratic penalty, which a steep χ² can overrule. Use constraints for
relations bounds cannot express, such as `$ring,diamout > $ring,diamin`.

## Keywords: `$R`, `$MU`, `$WL`, `$MJD`

| keyword | meaning | where |
|:--|:--|:--|
| `$R` | radius on the profile grid, mas | in a `profile` expression |
| `$MU` | `sqrt(1 - (r/r_max)^2)` | in a `profile` expression |
| `$WL` | wavelength, µm | outside profiles |
| `$MJD` | Modified Julian Date | outside profiles |

The expression editor checks each `$` reference while typing and says when a keyword is used
out of place or a name is unknown. Outside the editor, both errors appear as
`UndefVarError: WL not defined in OITOOLS`.

A parameter whose expression uses `$WL` or `$MJD`, or refers to one that does, has one value per
data point; the others have one value. A parameter of `$MJD` alone is computed once per epoch.
The parameter table shows such parameters at the dataset's first wavelength. `$MJD` comes from
`data.uv_mjd`, which is `Float64` in every dataset.

```julia
"disk,ud" => "3.0 * ($WL / 1.6)^1.2"
```

## Radial profiles

A component with a `profile` key is a free-form I(r), Hankel-transformed (**+ component → radial
profile**).

```julia
"disk,profile" => "exp(-($R / $scale)^2 / 2)",
"disk,scale"   => 1.5,      # a parameter of the expression
"disk,udout"   => 8.0,      # sets the r grid
"disk,nr"      => 64,       # optional: number of quadrature nodes
"disk,f"       => 1.0,
```

Names in the expression other than `$R` and `$MU` are parameters, looked up in the component and
then among the globals; the editor lists missing ones and offers to create them.

The grid's outer radius comes from the first of these keys present: `udout`, `diamout`, `diam`
(diameters, halved) or `r_max` (a radius). Without any of them the model does not build. The
panel shows the grid used and the key it came from.

`demos/models/YSO_*.toml` are the six models of Ibrahim et al. (2023, ApJ 947, 68) for HD 190073
(`demos/data/V1295_Aql.oifits`, epoch A): a star with a Rayleigh–Jeans spectrum, a disc and an
over-resolved halo, the disc being

| template | disc profile |
|:--|:--|
| `YSO_Doughnut` | PMOIRED's doughnut, `1 - (2*($R - ($RMIN + $RMAX)/2)/($RMAX - $RMIN))^2` on an annulus |
| `YSO_DoubleSigmoid` | `1/(1 + exp(-($R - $rin)/$sigin)) * 1/(1 + exp(($R - $rout)/$sigout))` |
| `YSO_AlphaDoubleSigmoid` | the same with a floor `$ain` inside the inner edge |
| `YSO_SkewedRing` | the double sigmoid with a first azimuthal mode (`az amp1`, `az projang1`) |
| `YSO_SkewedRing2` | and a second (`az amp2`, `az projang2`) |
| `YSO_OffCenterStar` | the skewed ring with a moving star (`star,x`, `star,y`) |

The panel plots I(r) and V(B) side by side.

### Azimuthal modes

A profile component can carry azimuthal modes, added below the profile editor as key pairs:

```julia
"disk,az amp1"     => 0.15,
"disk,az projang1" => 40.0,
```

Both keys are ordinary parameters. A new mode starts at zero amplitude. `az projang` follows
PMOIRED's convention, `cos(n·(ψ + φ − π/2))` (`_Vazvar` in its `oimodels.py`).

## What the parser understood

The inspector shows the model as `dict_to_model` reads it. Keys the parser does not recognise
are ignored by the fit and listed in red. Its validation lists values outside their bounds,
`lb >= ub`, and numeric flux fractions that do not sum to 1.

## Spectral lines

**Model Creation › Spectral lines**, for any component: **+ line** adds a line with a name, a
profile (Lorentzian, Gaussian, or PMOIRED's one-sided `truncexp` / `trunclor`), centre λ₀ (µm),
peak (negative for absorption) and width (nm), as parameters with PMOIRED's keys
(`c,line_<name>_wl0`, `_f`, `_<profile>`). **by velocity, rest** writes the centre as a rest
wavelength shifted by a global velocity `v_<component>_<name>` (km/s). **−** removes the line.
The SED view shows the lines. PMOIRED dicts with `line_*` keys import as they are.

## Orbits and radial velocities

**Geometry › orbit** puts a component on an orbit around the origin (`c,x = c,y = "orbit"`),
adding the elements as parameters (`c,orb P`, `MJD0`, `e`, `omega`, `OMEGA`, `incl`, `a`;
PMOIRED's keys and definitions, see the API page). PMOIRED dicts with `'x': 'orbit'` import as
they are; `'additional residuals'` is dropped, and a global such as `P-42.38` is renamed
`P_42_38`.

**RV table…** (Model Optimization) reads radial velocities from a CSV `MJD,RV,sigma,star` (`a`
primary, `b` secondary) and adds `Ka`, `Kb` and `gamma`. With **RVs** ticked, Levenberg-Marquardt,
NLopt and the bootstrap fit them with the interferometry; the console gives the χ² of each.
Examples: `demos/data/HD210763/rv_gallenne2023.csv`, `demos/data/Spica/rv_spica.csv`.

## Fitting

**Optimiser** (Model Optimization):

| entry | function | notes |
|:--|:--|:--|
| Levenberg-Marquardt | `fit_model_lsqfit` | |
| Nelder-Mead, L-BFGS, MMA, SLSQP, COBYLA | `fit_model` (NLopt) | MMA, SLSQP and COBYLA enforce constraints |
| Nested sampling | `fit_model_nested` | Nautilus or UltraNest; needs finite bounds |
| Ensemble slice sampler | `fit_model_ensemble` | settings: sweeps, seed |
| Parallel tempering | `fit_model_pigeons` | needs Pigeons; settings: rounds, chains, explorer (slice or MALA), reference (stabilised, variational or fixed), seed |
| Grid search | `chi2_map` | two parameters |

**Uncertainty**: *skip*, *Jacobian* (Levenberg-Marquardt's covariance), *bootstrap*
(`bootstrap_fit` from the fitted values, with nboot, mode, granularity and seed; not with
constraints, priors or samplers) or *posterior* (the spread of a sampler's posterior).

**Companion search…** searches for a point-source companion on a grid (`companion_search`);
**Multi-start…** fits from many starting points (below). **Similar fit on…** fits each ticked dataset
in turn (below).

Views under the fits table:

- **χ² map**: the grid search's map, or a companion or multi-start search.
- **Correlations**: the latest fit's parameter correlations.
- **Corner**: the latest fit's bootstrap replicates or posterior samples, the best fit marked
  (needs PairPlots).
- **SED**: total and per-component flux against the dataset's wavelengths, for models that use
  `$WL`.
- **Orbit**: the relative orbit on the sky with the dataset's epochs, and the radial velocities
  against phase.
- **Adopt** writes the selected fit's values into the table, when that fit was made on the
  selected model.

The fits table holds every fit of the session; with more than one model, each row names the
model it was made on.

Residuals are shown in Exploring (Compare › residuals). Each view has **Save PNG**.

!!! warning "`model_to_sed` takes metres"
    `wl_grid` is in metres, as `uv_lam` and every `wl` argument; only expressions see the
    wavelength in µm, as `$WL`.

## Multi-start search

**Multi-start…** (`multistart_fit`, PMOIRED's `gridFit`) fits the model from many starting
points and lists the distinct minima. One row per parameter, free ones first, then fixed ones
(greyed; a fixed parameter explored is held at each start's value, which gives a χ² profile).
Each row's starts are **start** (the table value), **grid** (min, max, step; several grids give
all combinations), **uniform** (min, max) or **normal** (mean, σ). Without a grid, **Starts** sets
their number; **Seed** fixes the draws. **Conditions** (one per line, e.g.
`sqrt($B,x^2 + $B,y^2) <= 12`) drop starts before fitting; **Δχ²ᵣ ≤** keeps minima within that
of the best. **PMOIRED expl** takes a PMOIRED `expl` dict (**Apply**).

Each start is a Levenberg-Marquardt fit with the table's bounds, the constraints, the ticked
observables and, if ticked, the radial velocities. As in PMOIRED, failed fits, fits without
uncertainties and fits that did not move are set aside and the rest grouped into minima; on an
orbit, minima whose MJD0 differ by whole periods are merged. The best minimum becomes the fit's
row.

The **χ² map** view shows the starts (`+`) joined to their minima, coloured by χ²ᵣ, set-aside
starts as red crosses. The boxes below the plot choose the two parameters and the minimum shown in
the fit's row.

## Fitting each dataset

**Similar fit on…** runs the current fit (model, bounds, constraints, priors, observables, optimiser)
on each ticked dataset, each from the table's values, typically the epochs made by Tools'
temporal binning. Levenberg-Marquardt or NLopt; Jacobian uncertainty or none. The **Epochs**
table lists each dataset's mean MJD, χ²ᵣ, ndof and fitted parameters with σ; **Copy**
(tab-separated) and **Save…** (`.tsv`).
