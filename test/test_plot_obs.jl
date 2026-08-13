using OITOOLS, PythonPlot

datadir = joinpath(@__DIR__, "..", "demos", "data")

# --- Polychromatic data (many wavelengths, multiple MJDs) ---
println("=== Polychromatic data (OBJECT1_N) ===")
data = readoifits(joinpath(datadir, "BC2026", "OBJECT1_N.oifits"), warn=false)[1,1]
display(data)

println("\n--- plot_obs(data) ---")
plot_obs(data)

println("--- plot_obs(data, color=\"wav\") ---")
plot_obs(data, color="wav")

println("--- plot_obs(data, color=\"mjd\") ---")
plot_obs(data, color="mjd")

println("--- plot_obs(data, logplot=true) ---")
plot_obs(data, logplot=true)

println("--- plot_obs custom selection ---")
plot_obs(data; obs=["V2", "T3PHI", "VISAMP", "VISPHI"], color="wav")

# --- Monochromatic data (single wavelength, single MJD) ---
println("\n=== Monochromatic data (2004-data1) ===")
data2 = readoifits(joinpath(datadir, "2004-data1.oifits"), warn=false)[1,1]
display(data2)

println("--- plot_obs(data2, color=\"wav\") — should fall back to baseline ---")
plot_obs(data2, color="wav")

println("--- plot_obs(data2, color=\"mjd\") — should fall back to baseline ---")
plot_obs(data2, color="mjd")

println("--- plot_obs(data2, logplot=true) ---")
plot_obs(data2, logplot=true)

println("\nAll plot_obs tests passed.")
PythonPlot.close("all")
