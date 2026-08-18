#
# CHARA observation planning — demonstrates the functions in plan.jl
#
using OITOOLS, Dates
using PythonPlot

# ── Target ────────────────────────────────────────────────────────────────────
targetname = "AZ Cyg"
ra, dec = ra_dec_from_simbad(targetname)   # decimal degrees

# ── Facility & night ──────────────────────────────────────────────────────────
facility = read_facility_file("CHARA")
obsdate  = DateTime(2026, 6, 3)

# Telescope config: 1=use, 0=unavailable, 2=reference cart
config = [1, 1, 1, 1, 1, 2]   # S1 S2 E1 E2 W1 W2(ref)

# ── Night observability ───────────────────────────────────────────────────────
obs = night_observability(facility, ra, dec, obsdate;
                          alt_limit=35.0, alt_max=80.0, moon_min_sep=30.0)
printstyled("HA range (35° < el < 80°): ",
            obs.ha[obs.good_alt[1]], " to ", obs.ha[obs.good_alt[end]], "\n";
            color=:red)
println("Moon illumination: ", round(obs.moon_fli * 100, digits=1), "%")
println("Moon separation: ", round(minimum(obs.moon_sep), digits=1), "° – ",
        round(maximum(obs.moon_sep), digits=1), "°")

# ── Find best POP configurations ─────────────────────────────────────────────
results = best_pop(facility, dec, obs.ha, config; n_best=5, min_minutes=10)
print_pop_results(facility, config, results)

# ── Gantt plot for a specific POP assignment ──────────────────────────────────
pop = [1, 2, 4, 5, 2, 5]
obs_plan(targetname, facility, ra, dec, obsdate, pop, config; alt_limit=35.0)

# ── CHARA-plan style delay plot ───────────────────────────────────────────────
chara_plan(targetname, facility, ra, dec, obsdate, pop, config; alt_limit=35.0)
