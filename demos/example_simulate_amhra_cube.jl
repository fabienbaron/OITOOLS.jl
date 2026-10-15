using OITOOLS
using Dates
# A Be star and its Keplerian disc in Brγ, from AMHRA's "Be model star + Keplerian disk"
# (Delaa et al. 2011), observed with CHARA: PMOIRED's example 4, first part, at CHARA.
#
# The cube is the one PMOIRED's example uses (PMOIRED_examples/DATA/BeModelsAMHRA): 256 × 256
# pixels of 2.658e-10 rad, 51 channels from 2.15962 µm in steps of 0.27 nm (R ≈ 8000).
# Its header has no CTYPE3, so its keywords are set once and the cube written anew.
cube_file = length(ARGS) > 0 ? ARGS[1] :
    joinpath(homedir(), "SOFTWARE", "PMOIRED_examples", "DATA", "BeModelsAMHRA",
             "output_Be_1675260915739.fits")
keyed = "./data/amhra_be_cube.fits"

# The pixel size in mas and one wavelength per plane; the file's reference pixel (CRPIX 128)
# is kept. The cube is read as stored, its first axis East — as PMOIRED reads it, which puts
# the disc at PA 120° (≡ -60°). AMHRA's own PA (60°) needs `east_left = true`.
fits_set_keywords(cube_file, keyed; pixsize = 2.65845e-10 * 180 / π * 3.6e6,
                  wavelengths = [2.1661e-6 + (k - 25) * 2.7e-10 for k in 1:51])
m = read_model_image(keyed)
println("cube: ", size(m.image), ", ", round(m.pixsize; digits = 5), " mas/pixel, λ ",
        round.(extrema(m.wavelengths) .* 1e6; digits = 5), " µm, offset ", m.offset, " mas")

# PMOIRED's target, (0h, -24°), at transit from CHARA (altitude 32°)
target = TargetConfig(target = "AMHRA_Be", raep0 = 0.0, decep0 = -24.0)
facility = read_facility_file("CHARA")
combiner = read_comb_file("MYSTIC")
day = DateTime(2024, 9, 22)
ts = collect(day:Minute(1):day + Hour(24))
_, ha = OITOOLS.hour_angle_calc(ts, facility.lon, target.raep0 / 15)
dates = [ts[argmin(abs.(mod.(ha .+ 12, 24) .- 12))]]

# The cube's own channels (`:image`), and PMOIRED's fixed relative errors: 0.01 % on the
# amplitudes and the flux, 0.01° on the phases.
out_file = "./data/amhra_be_chara.oifits"
simulate(facility, target, combiner, :image, dates, out_file; image = keyed, seed = 1,
         noise = FractionalNoise(v2 = 1e-4, visamp = 1e-4, t3amp = 1e-4, visphi = 0.01,
                                 t3phi = 0.01, flux = 1e-4))

data = readoifits(out_file; filter_bad_data = false, T = Float64)[1, 1]
println("V² $(data.nv2), VIS $(data.nvisamp), T3 $(data.nt3phi), flux $(data.nflux), ",
        length(unique(data.uv_lam)), " channels")
# With PythonPlot loaded: uvplot(data, color = "wav"), plot_t3phi(data), plot_visphi(data)
