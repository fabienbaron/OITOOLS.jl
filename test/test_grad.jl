using LinearAlgebra, Statistics, OITOOLS

function numgrad_1D(func;x=[], N=100, δ = 1e-6)
    if x==[]
        x = abs.(rand(Float64,N))
    else
        N = length(x)
    end
    numerical_g = zeros(length(x),length(func(x)))
    for i=1:N
        orig = x[i]
        x[i] = orig + 2*δ
        f2r = func(x)
        x[i] = orig + δ
        f1r = func(x)
        x[i] = orig - δ
        f1l = func(x)
        x[i] = orig - 2*δ
        f2l = func(x)
        numerical_g[i,:] .= (f2l-f2r+8*(f1r-f1l))
        x[i] = orig
    end
    numerical_g ./= (12*δ)
    return numerical_g;
end

oifitsfile = "./data/2004-data1.oifits"
pixsize = 0.4# size of a pixel in milliarcseconds
nx = 32 # width of image (number of pixels)
data = readoifits(oifitsfile)[1,1];
dft = setup_dft(data, nx, pixsize);
#initial image is a simple Gaussian
x = 10*vec(gaussian2d(nx,nx,nx/6));

# V2
f = xx->chi2_dft_f(xx, dft, data, weights = [1,0,0], verb=false, vonmises=false);
g = (xx,gg)->chi2_dft_fg(xx, gg, dft, data, weights = [1,0,0], verb=false, vonmises=false);
numerical_g = numgrad_1D(f,x=x,δ=1e-5);
analytic_g=similar(x);
g(x, analytic_g); analytic_g[:] = (analytic_g .- sum(vec(x).*analytic_g) / sum(x) ) / sum(x); # gradient correction to take into account the non-normalized image
#imdisp(reshape(numerical_g,nx,nx), figtitle="Numerical gradient")
#imdisp(reshape(analytic_g,nx,nx), figtitle="Analytic Gradient")
print("V2")
norm(numerical_g-analytic_g)
mean(numerical_g./analytic_g)

# T3AMP
f = xx->chi2_dft_f(xx, dft, data, weights = [0,1,0], verb=false, vonmises=false);
g = (xx,gg)->chi2_dft_fg(xx, gg, dft, data, weights = [0,1,0], verb=false, vonmises=false);
numerical_g = numgrad_1D(f,x=x,δ=1e-5);
analytic_g=similar(x);
g(x, analytic_g); analytic_g[:] = (analytic_g .- sum(vec(x).*analytic_g) / sum(x) ) / sum(x); # gradient correction to take into account the non-normalized image
#imdisp(reshape(numerical_g,nx,nx), figtitle="Numerical gradient")
#imdisp(reshape(analytic_g,nx,nx), figtitle="Analytic Gradient")
norm(numerical_g-analytic_g)
mean(numerical_g./analytic_g)

# T3PHI - Haniff expression
f = xx->chi2_dft_f(xx, dft, data, weights = [0,0,1], verb=false, vonmises=false);
g = (xx,gg)->chi2_dft_fg(xx, gg, dft, data, weights = [0,0,1], verb=false, vonmises=false);
numerical_g = numgrad_1D(f,x=x,δ=1e-5);
analytic_g=similar(x);
g(x, analytic_g); analytic_g[:] = (analytic_g .- sum(vec(x).*analytic_g) / sum(x) ) / sum(x); # gradient correction to take into account the non-normalized image
#imdisp(reshape(numerical_g,nx,nx), figtitle="Numerical gradient")
#imdisp(reshape(analytic_g,nx,nx), figtitle="Analytic Gradient")
norm(numerical_g-analytic_g)
mean(numerical_g./analytic_g)

# T3PHI - von Mises expression
data.t3phi_vonmises_err = gaussianwrapped_to_vonmises_fast.(data.t3phi_err/180*pi);
data.t3phi_vonmises_chi2_offset = log(2*pi) .+ 2*(logbesselI0.(data.t3phi_vonmises_err)-log.(data.t3phi_err/180*pi));
f = xx->chi2_dft_f(xx, dft, data, weights = [0,0,1], verb=false, vonmises=true);
g = (xx,gg)->chi2_dft_fg(xx, gg, dft, data, weights = [0,0,1], verb=false, vonmises=true);
numerical_g = numgrad_1D(f,x=x,δ=1e-5);
analytic_g=similar(x);
g(x, analytic_g); analytic_g[:] = (analytic_g .- sum(vec(x).*analytic_g) / sum(x) ) / sum(x); # gradient correction to take into account the non-normalized image
#imdisp(reshape(numerical_g,nx,nx), figtitle="Numerical gradient")
#imdisp(reshape(analytic_g,nx,nx), figtitle="Analytic Gradient")
norm(numerical_g-analytic_g)
mean(numerical_g./analytic_g)





#
# The following is leftover from development
#


# Conclusions:
# Confirmed good pairs:
#       chi2_visamp = norm(cvis_model)^2; ????
#          g_visamp = 2*real(transpose(dft)*conj(cvis_model))         1e-8

# angle(z)= atan2(Im(z),Re(z)) =  2 atan( Im(z)/(abs(z)+ Re(z))

# u = Im(z)/(abs(z)+ Re(z))

# d atan(u)/dz = 1/(1+u^2) . du/dz
# 1+u^2=2*abs(z)/(abs(z)+Re(z))   #OK
# 1/(1+u^2) =  (abs(z)+Re(z))/(2*abs(z)) # OK
#
# du/dz = d Im(z)/dz*1/(abs(z)+ Re(z))+Im(z)* d(1/(abs(z)+ Re(z)))/dz
# d(1/(abs(z)+ Re(z)))/dz =   - d(abs(z)+Re(z))/(abs(z)+ Re(z)))^2
#
# d atan(u)/dz = 1/(1+u^2) . du/dz = 1/(2*abs(z)) . [ d Im(z)/dz - Im(z)* d(abs(z)+Re(z))/dz*1/(abs(z)+ Re(z)))]
# = 1/(2*abs(z)) . [ -imag(H') - Im(z)/(abs(z)+ Re(z)) * d(abs(z)+Re(z))/dz)]


# von Mises
# chi2 = -2k*cos(t3phi-t3phi_data)

#f=y->angle.(H*y)
#f=y->imag.(H*y)   # derivative is -imag(H')
#f=y->real.(H*y)  # derivative is  real(H')
#f=y->abs.(H*y)
#  sqrt( Im^2 + Re^2)  =   d(Im^2+Re^2)/dz*1/(2*abs(z))
#                      =   (dIm Im + dRe Re)/(abs(z))
#                      =   (-imag(H') Im(z) + real(H') Re(z) )/abs(z)
#                      =   (((-imag(H'))'.*imag(z)+(real(H'))'.*real.(z))./abs.(z))'
#derivative            =    (-imag(H').*imag(z)'+real(H').*real.(z)')./abs.(z)'

# f=y->angle.(V)
# V=H*y/sum(y)
#diff                  = ( (imag(H).*real(V)-real(H).*imag.(V))./abs2.(V) )'
#                      =     imag(H.*(conj(V)./abs2.(V)))'
# then do diff_corrected =  (diff .- sum(y.*diff) / sum(y) ) / sum(y);

#
# function numgrad_1D(func;x=[], N=100, δ = 1e-6)
#     if x==[]
#         x = abs.(rand(Float64,N))
#     else
#         N = length(x)
#     end
#     numerical_g = zeros(length(x),length(func(x)))
#     for i=1:N
#         orig = x[i]
#         x[i] = orig + 2*δ
#         f2r = func(x)
#         x[i] = orig + δ
#         f1r = func(x)
#         x[i] = orig - δ
#         f1l = func(x)
#         x[i] = orig - 2*δ
#         f2l = func(x)
#         numerical_g[i,:] .= (f2l-f2r+8*(f1r-f1l))
#         x[i] = orig
#     end
#     numerical_g ./= (12*δ)
#     return numerical_g;
# end
#
# using LinearAlgebra, Statistics, OITOOLS, PyPlot
# oifitsfile = "./data/2004-data1.oifits"
# pixsize = 0.4 # size of a pixel in milliarcseconds
# nx = 32 # width of image (number of pixels)
# data = readoifits(oifitsfile)[1,1];
# dft = setup_dft(data, nx, pixsize);
# x = 10*vec(gaussian2d(nx,nx,nx/6));
# b=rand(size(dft,1))
#
# # nx = 32
# # nuv= 247
# # x=rand(nx*nx)
# # b = rand(nuv)
# # dft=rand(Complex{Float64}, nuv,nx*nx)
#
# # Argument
# V=dft*x/sum(x)
# f=y->angle.((dft*y)/sum(y))
# diffn = numgrad_1D(f,x=x, δ=1e-5);
# diff = imag(dft.*(conj(V)./abs2.(V)))';
# diff_corrected =  (diff .- sum(x.*diff) / sum(x) ) / sum(x);
# norm(diffn-diff_corrected)
# mean(diffn./diff_corrected)
# # Chi2 of Argument
# V=dft*x/sum(x)
# f=y->sum(abs2.(angle.((dft*y)/sum(y))-b))
# diffn = numgrad_1D(f,x=x, δ=1e-5)
# diff = 2*imag(transpose(dft)*(((angle.(V)-b))./abs2.(V).*conj(V)))
# diffg =  (diff .- sum(x.*diff) / sum(x) ) / sum(x);
# norm(diffn-diffg)
# mean(diffn./diffg)
#
# # Chi2 von Mises
# V=dft*x/sum(x)
# f=y->sum(-2*cos.(3*angle.((dft*y)/sum(y))-b))
# diffn = numgrad_1D(f,x=x, δ=1e-5)
# diff = 2*3*imag(transpose(dft)*(sin.(3*angle.(V)-b).*(conj(V)./abs2.(V))))
# diffg =  (diff .- sum(x.*diff) / sum(x) ) / sum(x);
# norm(diffn-diffg)
# mean(diffn./diffg)
#
# #Chi2 sparse
# V = image_to_vis_dft(x, dft); #or   V=dft*x/sum(x)
# indx1=data.indx_t3_1
# indx2=data.indx_t3_2
# indx3=data.indx_t3_3
# σ=data.t3phi_err
# b=data.t3phi
# #~, ~, t3phi_model = vis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
# f=y->norm(mod360(180/pi*angle.((dft*y)/sum(y))[indx1]+180/pi*angle.((dft*y)/sum(y))[indx2]+180/pi*angle.((dft*y)/sum(y))[indx3]-b)./σ)^2
# diffn = numgrad_1D(f,x=x, δ=1e-5)
# ~, ~, t3phi_model = vis_to_t3(V, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
#
# dt3=360/pi*(mod360(180/pi*angle.(V[indx1])+180/pi*angle.(V[indx2])+180/pi*angle.(V[indx3])-b)./σ.^2)
# dt3=360/pi*(mod360(t3phi_model-b)./σ.^2)
# diff = imag(transpose(dft[indx1,:])*(dt3./abs2.(V[indx1]).*conj(V[indx1]))
# + transpose(dft[indx2,:])*(dt3./abs2.(V[indx2]).*conj(V[indx2]))
# + transpose(dft[indx3,:])*(dt3./abs2.(V[indx3]).*conj(V[indx3])))
# diffg =  (diff .- sum(x.*diff) / sum(x) ) / sum(x);
# norm(diffn-diffg)
# mean(diffn./diffg)
