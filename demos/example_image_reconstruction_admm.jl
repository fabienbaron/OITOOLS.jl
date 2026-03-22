using OITOOLS
using PyPlot

# Load data
oifitsfile = joinpath(@__DIR__, "data", "2004-data1.oifits")
data = readoifits(oifitsfile)
d = data[1,1]
nx = 64; pixsize = 0.2

ft = setup_ft(data, nx, pixsize)
nfft_plan = ft[1,1][1]

x_start = gaussian2d(nx, nx, nx / 6)
println("Starting chi2:")
image_to_chi2(x_start, ft, data; verb=true)

# ADMM with TV + centering
x_admm = reconstruct_admm(x_start, d, nfft_plan, ft, data, nx;
                            rho_v2=10000.0, rho_t3=10000.0,
                            rho_reg=10000.0, mu_reg=1e5, mu_cen=1e2,
                            reg_type=:tv, maxit=300)

println("\nFinal chi2:")
image_to_chi2(x_admm, ft, data; verb=true)
imdisp(x_admm, pixsize=pixsize, figtitle="ADMM reconstruction")
savefig("admm_result.png", dpi=150, bbox_inches="tight")
