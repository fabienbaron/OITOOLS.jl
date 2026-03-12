# Simulating observations

OITOOLS can generate synthetic OIFITS datasets from a model or image, either by
copying the UV coverage of an existing file or by specifying observation times and
array geometry from scratch.

## From an existing OIFITS file

The simplest approach reuses the UV coverage and noise properties of a real dataset:

```julia
using OITOOLS
simulate_from_oifits("data/2004-data1.oifits", "data/sim.oifits";
                     image="data/2004true.fits", pixsize=0.101)
```

An analytic model can be used instead of an image:

```julia
model = create_model(create_component(type="ldlin", name="Star"))
dispatch_params([3.0, 0.15], model)
simulate_from_oifits("data/2004-data1.oifits", "data/sim.oifits";
                     model=model, pixsize=0.101)
```

See `example_simulate_observations_from_OIFITS.jl`.

## From observation times and array geometry

To simulate a full night of observations at a specific array:

```julia
using AstroTime
dates = collect(TAIEpoch(2018,8,13,3,0,0.0):15minutes:TAIEpoch(2018,8,13,8,30,0.0))
facility = read_facility_file("data/CHARA_new.txt")
target   = read_obs_file("data/default_obs.txt")
simulate_from_image("data/2004true.fits", 0.101, dates, facility, target, "data/sim.oifits")
```

See `example_simulate_observations_from_image.jl`.

## Observation planning

`example_chara_plan.jl` demonstrates how to check delay-line feasibility and
produce ASPRO-style Gantt charts for a given target and night:

```julia
gantt_onenight(facility, target, date)
```
