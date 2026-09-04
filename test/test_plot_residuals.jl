using OITOOLS, PythonPlot

datadir = joinpath(@__DIR__, "..", "demos", "data")

# --- Monochromatic data ---
println("=== Monochromatic (2004-data1) ===")
data = readoifits(joinpath(datadir, "BC2004", "2004-data1.oifits"), warn=false)[1,1]
nx = 64; pixsize = 0.2
ft = setup_nfft(data, nx, pixsize)
x = gaussian2d(nx, nx, nx/6)

println("--- plot_residuals baseline ---")
plot_residuals(x, ft, data)

println("--- plot_residuals wav (should fall back) ---")
plot_residuals(x, ft, data; color="wav")

println("--- plot_residuals mjd (should fall back) ---")
plot_residuals(x, ft, data; color="mjd")

println("--- plot_residuals logplot ---")
plot_residuals(x, ft, data; logplot=true)

# --- Single-observable residual plots ---
println("\n--- plot_v2_residuals ---")
obs = image_to_obs(x, ft, data)
plot_v2_residuals(data, obs.v2)

println("--- plot_t3phi_residuals ---")
plot_t3phi_residuals(data, obs.t3phi)

println("--- plot_t3amp_residuals ---")
plot_t3amp_residuals(data, obs.t3amp)

# --- Polychromatic data ---
println("\n=== Polychromatic (OBJECT1_N) ===")
data2 = readoifits(joinpath(datadir, "BC2026", "OBJECT1_N.oifits"), warn=false)[1,1]
nx2 = 64; pixsize2 = 0.5
ft2 = setup_nfft(data2, nx2, pixsize2)
x2 = gaussian2d(nx2, nx2, nx2/6)

println("--- plot_residuals baseline ---")
plot_residuals(x2, ft2, data2)

println("--- plot_residuals wav ---")
plot_residuals(x2, ft2, data2; color="wav")

println("--- plot_residuals mjd ---")
plot_residuals(x2, ft2, data2; color="mjd")

println("--- plot_residuals logplot ---")
plot_residuals(x2, ft2, data2; logplot=true)

# --- Single-observable with color kwarg ---
println("\n--- plot_v2_residuals color=wav ---")
obs2 = image_to_obs(x2, ft2, data2)
plot_v2_residuals(data2, obs2.v2; color="wav")

println("--- plot_t3phi_residuals color=mjd ---")
plot_t3phi_residuals(data2, obs2.t3phi; color="mjd")

println("\nAll plot_residuals tests passed.")
PythonPlot.close("all")
