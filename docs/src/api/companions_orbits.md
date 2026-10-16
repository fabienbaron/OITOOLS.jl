# Companions and Orbits

## Companion search

| Function | Description |
|----------|-------------|
| `companion_search(model_dict, data; companion, free, ...)` | Grid search for a companion around any primary model (including an image): the companion fitted at every node, a χ² map, the best node and its significance |
| `companion_grid(model_dict, data; companion, step, rmin, rmax, mask)` | The grid `companion_search` uses, and its number of nodes |
| `companion_flux_map(search)` | The companion's fitted flux at every node |
| `detection_sigma(chi2r_null, chi2r_best, ndof)` | CANDID's significance convention: the null model's χ² against a χ² distribution, as Gaussian σ |

With only the companion's flux free (CANDID's `chi2Map`), the model is linear in it and each
node is a 1-D fit (BraDi, or Halley with a BraDi fallback). With other parameters free (the
primary's diameter, a companion diameter), each node runs VMLMB with bounds. The default grid
step is λ/(2B_max), out to the smearing radius Rλ/(2B_max); an image primary's bright pixels are
masked. Bandwidth smearing follows the model's `"smear"`.

```@docs
companion_search
companion_grid
CompanionSearch
companion_flux_map
detection_sigma
```

## Orbits and radial velocities

A component on an orbit takes PMOIRED's keys and PMOIRED's definition (`_orbit`):
`"c,x" = "c,y" = "orbit"`, with the elements `c,orb P` (days), `MJD0` (periastron, MJD), `e`,
`omega`, `OMEGA`, `incl` (degrees) and `a` (mas), or `M` (M☉) and `plx` (mas) in place of `a`.
The position is computed in Float64 at each data point's MJD, once per epoch. ω is the
secondary's, relative to the primary (orbitize!'s convention, 180° from Gallenne et al. 2023's).
Each element can be fixed, free or an expression.

Radial velocities ([`read_rv`](@ref), a CSV `MJD,RV,sigma,star` with star `a` or `b`) are fitted
with the interferometry by `fit_model_lsqfit(…; rv)`, `fit_model(…; rv)` and
`bootstrap_fit(…; rv)`, with the velocity elements `c,orb Ka`, `Kb`, `gamma`, or `M` and
`q` (= Mb/Ma) with `incl`. `chi2` and `ndof` count them; `chi2_rv` is their part.

```julia
m = Dict{String,Any}("1,ud" => 0.2, "2,f" => 0.37, "2,x" => "orbit", "2,y" => "orbit",
                     "2,orb P" => 42.38, "2,orb MJD0" => 59785.45, "2,orb e" => 0.623,
                     "2,orb omega" => 114.1, "2,orb OMEGA" => 77.4, "2,orb incl" => 71.0,
                     "2,orb a" => 3.75, "2,orb Ka" => 50.3, "2,orb Kb" => 58.8, "2,orb gamma" => 14.9)
rv = read_rv("demos/data/HD210763/rv_gallenne2023.csv")
r  = fit_model_lsqfit(m, ["2,orb e", "2,orb a", "2,orb Ka", "2,orb Kb", "2,orb gamma"], data; rv)
```

| Function | Description |
|----------|-------------|
| `orbit_x(t, P, MJD0, e, omega, OMEGA, incl, a)`, `orbit_y(…)` | The secondary's offset East and North (mas) at MJD `t` |
| `orbit_rv(t, star, P, MJD0, e, omega, K, gamma)` | A star's radial velocity (km/s) |
| `orbit_a(P, M, plx)`, `orbit_semiamplitudes(P, e, incl, M, q)` | Kepler's third law: the semi-major axis, and Ka, Kb |
| `kepler_E(M, e)` | The eccentric anomaly, with derivatives for dual numbers |
| `read_rv(path)`, `RVData` | Radial velocities from a CSV table |

```@docs
orbit_x
orbit_y
orbit_rv
orbit_a
orbit_semiamplitudes
kepler_E
read_rv
RVData
```
