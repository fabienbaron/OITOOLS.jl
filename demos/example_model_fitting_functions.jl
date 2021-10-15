# Generate fake uv, compute visibilities for a model and invert to find the model Image
#include("../src/OITOOLS.jl"); using Main.OITOOLS
using OITOOLS
using PyPlot, FFTW, NFFT
N = 1024 # uv sampling
Bmax = 3000 #meters
λ = 1.6e-6 # meters  infrared H band

# UV GRID for irfft
# uv = Array{Float64}(undef,2, (div(N,2)+1)*N);
# x = collect(range(0, Bmax, length=div(N,2)+1))
# y = collect(range(-Bmax, Bmax, step=x[2]-x[1]))
# uv[1,:] = vec(repeat(x/λ,1, N))
# uv[2,:] = vec(repeat(y/λ,1, div(N,2)+1)')

uv = Array{Float64}(undef,2, N*N);
x = collect(range(-Bmax, Bmax, length=N))/λ
uv[1,:] = vec(repeat(x,1, N))
uv[2,:] = vec(repeat(x,1, N)')

function pad(mat)
    Nx,Ny=size(mat)
    padded = zeros(typeof(mat[1]),Nx*2,Ny*2)
    padded[div(Nx,2):div(Nx,2)+Nx-1,div(Ny,2):div(Ny,2)+Ny-1] = mat
    return padded
end


# # UV POLAR GRID
# r = collect(range(0, 2*Bmax, length=N))/λ
# θ = collect(range(0, 360, length=360))/180*pi
# uv = Array{Float64}(undef,2, N*length(θ))
# uv[1,:] = vec(r'.*cos.(θ))
# uv[2,:] = vec(r'.*sin.(θ))

# Setup Model imaging
nx = 256       # number of pixels
pixsize = 0.1 # mas/pixel
fftplan  = plan_nfft(pixsize * (pi / 180.0) / 3600000.0*[-1;1].*uv, (nx,nx), 4, 2.0);

# Setup model
V= Complex.(visibility_ud([10.0], uv)) # works
img = real.(nfft_adjoint(fftplan, V)); img = img.*(img .>0); imdisp(img, pixscale=pixsize); scatter(0,0, marker="*", color=:red)

#imdisp((abs.(fftshift(ifft(reshape(V,(N,N)))))).^.5)

V=  Complex.(visibility_annulus([1.0,2.0],uv))
img = real.(nfft_adjoint(fftplan, V)); img = img.*(img .>0); imdisp(img, pixscale=pixsize); scatter(0,0, marker="*", color=:red)

V = Complex.(visibility_ldpow([10.0, 0.3], uv)) # works
img = real.(nfft_adjoint(fftplan, V)); img = img.*(img .>0); imdisp(img, pixscale=pixsize); scatter(0,0, marker="*", color=:red)

V = Complex.(visibility_ldlin([10.0, 0.7], uv)) # works
img = real.(nfft_adjoint(fftplan, V)); img = img.*(img .>0); imdisp(img, pixscale=pixsize); scatter(0,0, marker="*", color=:red)

V = Complex.(visibility_ldsquareroot([10.0, 0.5, 0.3], uv))  #  works
img = real.(nfft_adjoint(fftplan, V)); img = img.*(img .>0); imdisp(img, pixscale=pixsize); scatter(0,0, marker="*", color=:red)

V = Complex.(visibility_ldquad([10.0, 0.5, 0.5], uv))  #  works
img = real.(nfft_adjoint(fftplan, V)); img = img.*(img .>0); imdisp(img, pixscale=pixsize); scatter(0,0, marker="*", color=:red)

V = Complex.(visibility_ellipse_uniform([10.0,2.0,-45.0], uv))  # works
img = real.(nfft_adjoint(fftplan, V)); img = img.*(img .>0); imdisp(img, pixscale=pixsize); scatter(0,0, marker="*", color=:red)

V = Complex.(visibility_ellipse_quad([10.0,0.5, 0.5, 2.0,-45.0], uv)) # works
img = real.(nfft_adjoint(fftplan, V)); img = img.*(img .>0); imdisp(img, pixscale=pixsize); scatter(0,0, marker="*", color=:red)

V = Complex.(visibility_thin_ring([10.0, 0, 0],uv)) # works
img = real.(nfft_adjoint(fftplan, V)); img = img.*(img .>0); imdisp(img, pixscale=pixsize); scatter(0,0, marker="*", color=:red)

V = Complex.(visibility_Gaussian_ring([10.0, 90, 60, .1],uv))
img = real.(nfft_adjoint(fftplan, V)); img = img.*(img .>0); imdisp(img, pixscale=pixsize); scatter(0,0, marker="*", color=:red)

V = Complex.(visibility_Lorentzian_ring([10.0, 90, 60, .1],uv))
img = real.(nfft_adjoint(fftplan, V)); img = img.*(img .>0); imdisp(img, pixscale=pixsize); scatter(0,0, marker="*", color=:red)

V = Complex.(visibility_Gaussian_ring_az([10.0, -90, 60, .15, 0.5, 0,0,0],uv))
img = real.(nfft_adjoint(fftplan, V)); img = img.*(img .>0); imdisp(img, pixscale=pixsize); scatter(0,0, marker="*", color=:red)

V = Complex.(visibility_Gaussian_ring_az([10.0, -90, 60, .15, 0, 0.5,0,0],uv))
img = real.(nfft_adjoint(fftplan, V)); img = img.*(img .>0); imdisp(img, pixscale=pixsize); scatter(0,0, marker="*", color=:red)

V = Complex.(visibility_Gaussian_ring_az([10.0, -90, 60, .15, 0, 0,0.5,0],uv))
img = real.(nfft_adjoint(fftplan, V)); img = img.*(img .>0); imdisp(img, pixscale=pixsize); scatter(0,0, marker="*", color=:red)

V = Complex.(visibility_Gaussian_ring_az([10.0, -90, 60, .15, 0, 0,0,0.5],uv))
img = real.(nfft_adjoint(fftplan, V)); img = img.*(img .>0); imdisp(img, pixscale=pixsize); scatter(0,0, marker="*", color=:red)

V = Complex.(visibility_GaussianLorentzian_ring_az([10.0, -90, 60, .15, 0, 0,0, 0, .4],uv))
img = real.(nfft_adjoint(fftplan, V)); img = img.*(img .>0); imdisp(img, pixscale=pixsize); scatter(0,0, marker="*", color=:red)


# Example with the new interface
c1 = create_component(type="ring", name="ring")
c2 = create_component(type="ud", name="star")
model = create_model(c1,c2);
dispatch_params([0.5, 10.0, -90, 60, .15, 0, 0,0, 0, .4],model)
model_to_image(model, nx=256, pixsize=0.1)


# model = create_model(create_component(type="ldpow", name="Model"));



function model_to_image2(nx, pixsize, model)
	s1=-1;
    s2=-1;
	N = 4*nx;
	uint = 1.0/(N/4*pixsize * pi/(180*3600000));
	# Compute visibility
    matin=zeros(ComplexF64,N,N)



    matin=zeros(ComplexF64,N,N)
    i=collect(0:Int(N/2))
    k=collect(0:Int(N/2))

    u = s2 * k * uint;
    v = s1 * i * uint;
    A = model_to_cvis(model,Array([u v]'),[1.6e-6]);
    matin[CartesianIndex.(i.+1,k.+1)] = A.*(-1).^(i + k)
    u = s2*k*uint;
    v = s1*(-i.-1)*uint;
    A = model_to_cvis(model,Array([u v]'),[1.6e-6]);
    matin[CartesianIndex.(k.+1, N.-i)] = A.*(-1).^((N-1).-i+k)

    u = s2*(-k.-1)*uint;
    v = s1*i*uint;
    A = model_to_cvis(model,Array([u v]'),[1.6e-6]);
    matin[CartesianIndex.(N.-k,i.+1)] = A.*(-1).^(i-k.+(N-1));

    u = s2*(-k.-1)*uint;
    v = s1*(-i.-1)*uint;
    A = model_to_cvis(model,Array([u v]'),[1.6e-6]);
    matin[CartesianIndex.(N.-k,N.-i)] = A.*(-1).^(-i-k);

	for i=0:Int(N/2)
	    for k=0:Int(N/2)
			v = s1 * i * uint;
			u = s2 * k * uint;
			A = model_to_cvis(model,Array([u v]'),[1.6e-6]);
			matin[k+1, i+1] =  A[1]  * (-1)^(i + k)
			v = s1*(-i-1)*uint;
			u = s2*k*uint;
            A = model_to_cvis(model,Array([u v]'),[1.6e-6]);
			matin[k+1, N-i-1+1] = A[1] * (-1)^((N-i-1) + k)
			v = s1*i*uint;
			u = s2*(-k-1)*uint;
            A = model_to_cvis(model,Array([u v]'),[1.6e-6]);
			matin[N-k-1+1,i+1] = A[1] * (-1)^(i+N-k-1);
			v = s1*(-i-1)*uint;
			u = s2*(-k-1)*uint;
            A = model_to_cvis(model,Array([u v]'),[1.6e-6]);
			matin[N-k-1+1,N-i-1+1] = A[1] * (-1)^(-i-k);
        end
    end

    #Do fft
    matout = fft(matin)
    #Downsample -- bad?
    out = zeros(Float64, nx,nx)
	for i=0:nx-1
		for k=0:nx-1
			out[i+1,k+1]=real(matout[4i+1, 4k+1] + matout[4i+1,4k+2] +
			matout[4i+1,4k+3] + matout[4i+1,4k+4] +
			matout[(4i+1)+1,4k+1] + matout[(4i+1)+1,4k+2]+
			matout[(4i+1)+1,4k+3] + matout[(4i+1)+1,4k+4]+
			matout[(4i+2)+1,4k+1] + matout[(4i+2)+1,4k+2]+
			matout[(4i+2)+1,4k+3] + matout[(4i+2)+1,4k+4]+
			matout[(4i+3)+1,4k+1] + matout[(4i+3)+1,4k+2]+
			matout[(4i+3)+1,4k+3] + matout[(4i+3)+1,4k+4])/16.;
        end
    end
    out = out .* (out .>0)
    out = out / sum(out)
return out
end
