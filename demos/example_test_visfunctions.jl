# Generate fake uv, compute visibilities for a model and invert to find the model Image
#include("../src/OITOOLS.jl"); using Main.OITOOLS
using OITOOLS
using PyPlot, NFFT
N = 201 # uv sampling
Bmax = 1000 #meters
λ = 1.6e-6 # meters

# UV GRID
x = collect(range(-Bmax, Bmax, length=N))/λ
uv = Array{Float64}(undef,2, N*N)
uv[1,:] = vec(repeat(x,1, N))
uv[2,:] = vec(repeat(x,1, N)')

# # UV POLAR GRID
# r = collect(range(0, 2*Bmax, length=N))/λ
# θ = collect(range(0, 360, length=360))/180*pi
# uv = Array{Float64}(undef,2, N*length(θ))
# uv[1,:] = vec(r'.*cos.(θ))
# uv[2,:] = vec(r'.*sin.(θ))


# Setup model
#V = Complex.(visibility_ud([10.0], uv)) # works
#V = Complex.(visibility_ldpow([10.0, 0.3], uv)) # bug
#V = Complex.(visibility_ldlin([10.0, 0.7], uv)) # does not seem limb darkened -> check expression
#V = Complex.(visibility_ldquad([10.0, 0.5, 0.5], uv))  #  works
#V = Complex.(visibility_ellipse_uniform([10.0,2.0,-45.0], uv))  # works
#V = Complex.(visibility_ellipse_quad([10.0,0.5, 0.5, 2.0,-45.0], uv)) # works
#V = Complex.(visibility_thin_ring([10.0, 45, 60],uv)) # works
#V = Complex.(visibility_Gaussian_ring([10.0, 90, 60, .1],uv))
V = Complex.(visibility_Gaussian_ring_az([10.0, 50, 60, .3, 20, 0, 0.5,0,0],uv))

# Image of the model
nx = 512       # number of pixels
pixsize = 0.05 # mas/pixel
fftplan  = plan_nfft(pixsize * (pi / 180.0) / 3600000.0*[-1;1].*uv, (nx,nx), 4, 2.0);
img = real.(nfft_adjoint(fftplan, V))
imdisp(img, pixscale=pixsize)
scatter(0,0, marker="*", color=:red)
