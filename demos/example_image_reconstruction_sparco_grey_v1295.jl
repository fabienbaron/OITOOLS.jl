# SPARCO image reconstruction of V1295 Aql (HD 190073)
#
# Physical model:
#   - Star:        point source, flux fraction f_star, spectrum ∝ (λ/λ0)^(-4)
#   - Background:  over-resolved (V=0), flux fraction f_bg, spectrum ∝ (λ/λ0)^(-4)
#   - Environment: grey NFFT image, flux fraction 1 - f_star - f_bg,
#                  spectrum ∝ (λ/λ0)^(d_env)
#
# T3amp is disabled (weights=[1, 0, 1]) — this dataset has noisy T3amp.

using OITOOLS

# Load data and setup
oifitsfile = joinpath(@__DIR__, "data", "2019_v1295Aql.WL_SMOOTH.A.oifits")
data = readoifits(oifitsfile, filter_bad_data=true)[1,1]

pixsize = 0.125  # mas/pixel
nx = 64
ft = setup_nfft(data, nx, pixsize)

x_start = rand(nx, nx)
x_start ./= sum(x_start)

# Reconstruct
result = reconstruct_sparco(x_start, data, ft;
    lambda_ref = 1.6e-6,
    star_flux  = 0.48,
    bg_flux    = 0.0,
    d_env      = -0.1,
    weights    = [1.0, 0.0, 1.0],
    regularizers = [["tv", 1e2]],
    maxiter    = 200,
    rounds     = 3,
    verb       = true)

println("\nFinal chi2 = ", result.chi2)
println("Parameters: ", result.param_names)

imdisp(result.image, pixsize=pixsize, figtitle="V1295 Aql — SPARCO")
