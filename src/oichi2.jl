using LinearAlgebra, SpecialFunctions, NFFT

function setup_dft(uv::Array{Float64,2}, nx, pixsize)
    scale_rad = pixsize * (pi / 180.0) / 3600000.0
    nuv = size(uv, 2)
    span = collect((2 * (1:nx) .- (nx + 1)) * scale_rad * pi)
    xvals = reshape(span, 1, nx, 1)
    yvals = reshape(span, 1, 1, nx)
    dft = reshape(cis.((uv[1, :] .* xvals) .- (uv[2, :] .* yvals)), nuv, nx^2)
    return dft
end

function setup_dft(data::OIdata, nx, pixsize)
    return setup_dft(data.uv, nx, pixsize)
end

function setup_nfft(data::OIdata, nx, pixsize)::Array{NFFT.NFFTPlan{2,0,Float64},1}
    scale_rad = pixsize * (pi / 180.0) / 3600000.0*[-1;1].*data.uv;
    fftplan_uv  = plan_nfft(scale_rad, (nx,nx), 4, 2.0);
    fftplan_vis  = plan_nfft(scale_rad[:, data.indx_vis], (nx,nx), 4, 2.0);
    fftplan_v2   = plan_nfft(scale_rad[:, data.indx_v2], (nx,nx), 4, 2.0);
    fftplan_t3_1 = plan_nfft(scale_rad[:, data.indx_t3_1], (nx,nx), 4, 2.0);
    fftplan_t3_2 = plan_nfft(scale_rad[:, data.indx_t3_2], (nx,nx), 4, 2.0);
    fftplan_t3_3 = plan_nfft(scale_rad[:, data.indx_t3_3], (nx,nx), 4, 2.0);
    return [fftplan_uv,fftplan_vis, fftplan_v2,fftplan_t3_1,fftplan_t3_2,fftplan_t3_3]
end


function setup_nfft(uv::Array{Float64,2}, indx_vis, indx_v2, indx_t3_1, indx_t3_2,indx_t3_3, nx, pixsize)::Array{NFFT.NFFTPlan{2,0,Float64},1}
    scale_rad = pixsize * (pi / 180.0) / 3600000.0*[-1;1].*uv;
    fftplan_uv  = plan_nfft(scale_rad, (nx,nx), 4, 2.0);
    fftplan_vis  = plan_nfft(scale_rad[:, indx_vis], (nx,nx), 4, 2.0);
    fftplan_v2   = plan_nfft(scale_rad[:, indx_v2], (nx,nx), 4, 2.0);
    fftplan_t3_1 = plan_nfft(scale_rad[:, indx_t3_1], (nx,nx), 4, 2.0);
    fftplan_t3_2 = plan_nfft(scale_rad[:, indx_t3_2], (nx,nx), 4, 2.0);
    fftplan_t3_3 = plan_nfft(scale_rad[:, indx_t3_3], (nx,nx), 4, 2.0);
    return [fftplan_uv,fftplan_vis,fftplan_v2,fftplan_t3_1,fftplan_t3_2,fftplan_t3_3]
end


function setup_nfft_t4(data, nx, pixsize)::Array{NFFT.NFFTPlan{2,0,Float64},1}
    scale_rad = pixsize * (pi / 180.0) / 3600000.0*[-1;1].*data.uv
    fftplan_uv  = plan_nfft(scale_rad, (nx,nx), 4, 2.0);
    fftplan_vis  = plan_nfft(scale_rad[:, data.indx_vis], (nx,nx), 4, 2.0);
    fftplan_v2   = plan_nfft(scale_rad[:, data.indx_v2], (nx,nx), 4, 2.0);
    fftplan_t3_1 = plan_nfft(scale_rad[:, data.indx_t3_1], (nx,nx), 4, 2.0);
    fftplan_t3_2 = plan_nfft(scale_rad[:, data.indx_t3_2], (nx,nx), 4, 2.0);
    fftplan_t3_3 = plan_nfft(scale_rad[:, data.indx_t3_3], (nx,nx), 4, 2.0);
    fftplan_t4_1 = plan_nfft(scale_rad[:, data.indx_t4_1], (nx,nx), 4, 2.0);
    fftplan_t4_2 = plan_nfft(scale_rad[:, data.indx_t4_2], (nx,nx), 4, 2.0);
    fftplan_t4_3 = plan_nfft(scale_rad[:, data.indx_t4_3], (nx,nx), 4, 2.0);
    fftplan_t4_4 = plan_nfft(scale_rad[:, data.indx_t4_4], (nx,nx), 4, 2.0);
    return [fftplan_uv,fftplan_vis,fftplan_v2,fftplan_t3_1,fftplan_t3_2,fftplan_t3_3, fftplan_t4_1,fftplan_t4_2,fftplan_t4_3, fftplan_t4_4]
end


function setup_nfft_multiepochs(data, nx, pixsize)::Array{Array{NFFT.NFFTPlan{2,0,Float64},1},1}
    nepochs = size(data,1);
    scale_rad = pixsize * (pi / 180.0) / 3600000.0;
    fftplan_multi = Array{Any}(undef, nepochs);
    for i=1:nepochs
        fftplan_multi[i]=setup_nfft(data[i], nx, pixsize);
    end
    return fftplan_multi
end


function setup_nfft_polychromatic(data, nx, pixsize)::Array{Array{NFFT.NFFTPlan{2,0,Float64},1},1}
    nwavs = size(data,1);
    scale_rad = pixsize * (pi / 180.0) / 3600000.0;
    fftplan_multi = Array{Any}(undef, nwavs);
    for i=1:nwavs
        fftplan_multi[i]=setup_nfft(data[i], nx, pixsize);
    end
    return fftplan_multi
end

#
# using FINUFFT
#
# function image_to_cvis_finufft(x::Array{Float64,2})
#   out2 = nufft2d2(pixsize * (pi / 180.0) / 3600000.0*data.uv[1,:], -pixsize * (pi / 180.0) / 3600000.0*data.uv[2,:],-1, 1e-6, Complex{Float64}.(x) / sum(x))
# end
#
# function image_to_cvis_finufft(x::Array{Float64,1})
# nx = Int64(sqrt(length(x)))
# out2 = nufft2d2(pixsize * (pi / 180.0) / 3600000.0*data.uv[1,:], -pixsize * (pi / 180.0) / 3600000.0*data.uv[2,:],-1, 1e-6, Complex{Float64}.(reshape(x,(nx,nx))) / sum(x))
# end


function mod360(x)
    mod.(mod.(x.+180.0,360.0).+360.0, 360.0) .- 180.0
end

function cvis_to_v2(cvis, indx)
    v2_model = abs2.(cvis[indx]);
end

function cvis_to_t3(cvis, indx1, indx2, indx3)
    t3 = cvis[indx1].*cvis[indx2].*cvis[indx3];
    t3amp = abs.(t3);
    t3phi = angle.(t3)*180.0/pi;
    return t3, t3amp, t3phi
end

function cvis_to_t4(cvis, indx1, indx2, indx3, indx4)
    t4 = cvis[indx1].*cvis[indx2]./(cvis[indx3]*conj(cvis[indx4]));
    t4amp = abs.(t4);
    t4phi = angle.(t4)*180.0/pi;
    return t4, t4amp, t4phi
end


function image_to_cvis_dft(x, dft)
    cvis_model = Array{Complex{Float64}}(undef,size(dft, 1));
    cvis_model = dft * vec(x) / sum(x);
end

function image_to_cvis_nfft(x::Array{Float64,1}, nfft_plan::NFFT.NFFTPlan)
    nx = Int64(sqrt(length(x)))
    cvis_model = nfft(nfft_plan, Complex{Float64}.(reshape(x,(nx,nx)))) / sum(x);
end

function image_to_cvis_nfft(x::Array{Float64,2}, nfft_plan::NFFT.NFFTPlan)
    cvis_model = nfft(nfft_plan, Complex{Float64}.(x)) / sum(x);
end

# Overload
function image_to_cvis_nfft(x::Array{Float64,1}, nfft_plan::Array{NFFT.NFFTPlan{2,0,Float64},1})
    flux = sum(x);
    nx = Int64(sqrt(length(x)))
    cvis_model = nfft(nfft_plan[1], Complex{Float64}.(reshape(x,(nx,nx)))) / flux;
end

function image_to_cvis_nfft(x::Array{Float64,2}, nfft_plan::Array{NFFT.NFFTPlan{2,0,Float64},1})
    flux = sum(x);
    cvis_model = nfft(nfft_plan[1], Complex{Float64}.(x)) / flux;
end


function image_to_cvis_finufft(x, data)
    flux = sum(x);
    if (ndims(x) == 1)
        nx = Int64(sqrt(length(x)))
        cvis_model = nfft(nfft_plan, Complex{Float64}.(reshape(x,(nx,nx)))) / flux;
    else
        cvis_model = nfft(nfft_plan, Complex{Float64}.(x)) / flux;
    end
end

#
# function chi2_vis_dft_fg(x, g, dft, data ) # criterion function plus its gradient w/r x
#   cvis_model = image_to_cvis_dft(x, dft);
#   # compute observables from all cvis
#   visamp_model = abs.(cvis_model);
#   visphi_model = angle.(cvis_model)*(180.0/pi);
#   chi2_visamp = norm((visamp_model - data.visamp)./data.visamp_err)^2;
#   chi2_visphi = norm(mod360(visphi_model - data.visphi)./data.visphi_err)^2;
#   # Original formulas
#   # g_visamp = 2.0*sum(((visamp_model-data.visamp)./data.visamp_err.^2).*real( conj(cvis_model./visamp_model).*dft),1);
#   # g_visphi = 360.0/pi*sum(((mod360(visphi_model-data.visphi)./data.visphi_err.^2)./abs2.(cvis_model)).*(-imag(cvis_model).*real(dft)+real(cvis_model).*imag(dft)),1);
#   # Improved formulas
#   g_visamp = 2.0*real(transpose(dft)*(conj(cvis_model./visamp_model).*(visamp_model-data.visamp)./data.visamp_err.^2))
#   g_visphi = 360.0/pi*imag(transpose(dft)*((mod360(visphi_model-data.visphi)./data.visphi_err.^2)./(visamp_model.^2).*conj(cvis_model)));
#   g[:] = g_visamp + g_visphi;
#   flux = sum(x);
#   g[:] = (g - sum(x.*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
# #  imdisp(x);
#   println("VISAMP: ", chi2_visamp/data.nvisamp, " VISPHI: ", chi2_visphi/data.nvisphi, " Flux: ", flux)
#   return chi2_visamp + chi2_visphi
# end

# function chi2_vis_nfft_fg(x, g, fftplan::FFTPLAN, data; mu = 1e7 ) # criterion function plus its gradient w/r x
#   # Total squared variation
#   # Total variation
#     ϵ=1e-10
#     nx = Int(sqrt(length(x)))
#     y = reshape(x, nx, nx);
#     yx = circshift(y,(0,-1));
#     yy = circshift(y,(-1,0));
#     tv_f = sum(sqrt.((y-yx).^2+(y-yy).^2+ϵ))
#     tv_g =  vec((2*y-yx-yy)./(sqrt.((y-yx).^2+(y-yy).^2+ϵ)));
#
#
#   cvis_model = image_to_cvis_nfft(x, fftplan.fftplan_uv);
#   # compute observables from all cvis
#   visamp_model = abs.(cvis_model);
#   visphi_model = angle.(cvis_model)*(180.0/pi);
#   chi2_visamp = norm((visamp_model - data.visamp)./data.visamp_err)^2;
#   chi2_visphi = norm(mod360(visphi_model - data.visphi)./data.visphi_err)^2;
#   g_visamp = 2.0*real(nfft_adjoint(fftplan.fftplan_uv,(cvis_model./visamp_model.*(visamp_model-data.visamp)./data.visamp_err.^2)));
#   g_visphi = 360.0/pi*-imag(nfft_adjoint(fftplan.fftplan_uv,cvis_model.*((mod360(visphi_model-data.visphi)./data.visphi_err.^2)./visamp_model.^2)));
#   g[:] = vec(g_visamp + g_visphi) +  mu * tv_g;
#   flux = sum(x);
#   g[:] = (g - sum(x.*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
#   println("VISAMP: ", chi2_visamp/data.nvisamp, " VISPHI: ", chi2_visphi/data.nvisphi, " Flux: ", flux)
#   return chi2_visamp + chi2_visphi +  mu * tv_f;
# end

function gaussian2d(n,m,sigma)
    g2d = [exp(-((X-(m/2)).^2+(Y-n/2).^2)/(2*sigma.^2)) for X=1:m, Y=1:n]
    return g2d
end

# function create_start_image(type, nx, pixsize, param)
#   if type="gaussian"
#   gaussian2d(nx,nx,sigma)
#   elseif type="disc"
#     x =

#   end
# return
# end


function reg_centering(x,g; verb = false) # takes a 1D array
    nx = Int(sqrt(length(x)))
    flux= sum(x)
    c = cdg(reshape(x,nx,nx))
    f = (c[1]-(nx+1)/2)^2+(c[2]-(nx+1)/2)^2
    xx = [mod(i-1,nx)+1 for i=1:nx*nx]
    yy = [div(i-1,nx)+1 for i=1:nx*nx]
    g[:] = 2*(c[1]-(nx+1)/2)*xx + 2*(c[2]-(nx+1)/2)*yy
    if verb == true
        print(" COG:", c, " REGC: ", f);
    end
    return f
end

function tvsq(x,tvsq_g; verb = false)
    # Total squared variation
    nx = Int(sqrt(length(x)))
    y = reshape(x, nx, nx);
    lx = circshift(y,(0,-1));
    ly = circshift(y,(-1,0));
    rx = circshift(y,(0,1));
    ry = circshift(y,(1,0));
    tvsq_f = norm(y-rx)^2+norm(y-ry)^2
    tvsq_g[:] = 2*vec(4*y-lx-ly-rx-ry)
    if verb == true
        print(" TVSQ:", tvsq_f);
    end
    return tvsq_f
end

function tv(x,tv_g; verb = false, ϵ=1e-8)
    # Total variation
    # TODO: - treat edges properly
    #       - check vs matrix implementation
    #       - check performance/precision vs FFT version
    nx = Int(sqrt(length(x)))
    y = reshape(x, nx, nx);
    xijplus1  = circshift(y, (0,-1))
    xijminus1 = circshift(y, (0,1))
    xiplus1j  = circshift(y, (-1,0))
    ximinus1j = circshift(y, (1,0))
    ximinus1jplus1 = circshift(y, (1,-1))
    xiplus1jminus1 = circshift(y, (-1,1))
    d1 = sqrt.((xiplus1j-y).^2  + (xijplus1-y).^2  .+ ϵ^2)
    d2 = sqrt.((y-ximinus1j).^2 + (ximinus1jplus1-ximinus1j).^2 .+ ϵ^2)
    d3 = sqrt.((xiplus1jminus1-xijminus1).^2 + (y-xijminus1).^2 .+ ϵ^2)
    tv_f = sum(d1.-ϵ)
    tv_g[:] = vec( (2*y - xiplus1j - xijplus1)./d1 + (y-ximinus1j)./d2 + (y-xijminus1)./d3 )

    if verb == true
        print(" TV:", tv_f);
    end
    return tv_f
end

function l1l2(x, g; verb = false, ϵ=1e-8, α = 2.0)
    # Isotropic L1-L2 / Fair loss function /  Mugnier-Brette expression
    nx = Int(sqrt(length(x)))
    y = reshape(x, nx, nx);
    xijplus1  = circshift(y, (0,-1))
    xijminus1 = circshift(y, (0,1))
    xiplus1j  = circshift(y, (-1,0))
    ximinus1j = circshift(y, (1,0))
    ximinus1jplus1 = circshift(y, (1,-1))
    xiplus1jminus1 = circshift(y, (-1,1))
    d1 = sqrt.((xiplus1j-y).^2  + (xijplus1-y).^2  .+ ϵ^2)
    d2 = sqrt.((y-ximinus1j).^2 + (ximinus1jplus1-ximinus1j).^2 .+ ϵ^2)
    d3 = sqrt.((xiplus1jminus1-xijminus1).^2 + (y-xijminus1).^2 .+ ϵ^2)
    f = α^2*sum(d1/α - log.(1.0 .+ d1/α ) .-ϵ)
    g[:] = α^2*vec( ( (2*y - xiplus1j - xijplus1)./d1 + (y-ximinus1j)./d2 + (y-xijminus1)./d3 )/α -(2*y - xiplus1j - xijplus1)./(d1.*(α .+ d1)) - (y-ximinus1j)./(d2.*(α .+ d2)) - (y-xijminus1)./(d3.*(α .+ d3)) )
    if verb == true
        print(" ℓ1ℓ2:", f);
    end
    return f
end

function l2sq(x,g; verb = false)
    f = sum(x.^2)
    g[:] =  2*x;
    if verb == true
        print(" ℓ2^2:", f);
    end
    return f
end


function l1hyp(x,g; verb = false,ϵ=1e-9)
    f = sum(sqrt.(x.^2 .+ϵ^2).-ϵ)
    g[:] =  x./sqrt.(x.^2 .+ϵ^2);
    if verb == true
        print(" ℓ1hyp:", f);
    end
    return f
end

function l1l2w(x,g; verb = false)
    f = sum(x-log.(1.0 .+ x))
    g[:] .=  1.0 .- 1.0./(1.0 .+x);
    if verb == true
        print(" ℓ1ℓ2w:", f);
    end
    return f
end

function entropy(x,g; verb = false, ϵ=1e-12)
    f = sum(x.*log.(abs.(x).+ϵ) - x)
    g[:] .=  log.(abs.(x).+ϵ);
    if verb == true
        print(" MAXENT:", f);
    end
    return f
end

function compactness(x,g; verb = false)
    nx = Int(sqrt(length(x)))
    y = repeat(collect(1:nx).-0.5*(nx-1),1,nx)
    yy = y.*y
    rr = vec(yy+yy')/(nx*nx)
    f = sum(x.*rr)
    g[:] .=  rr;
    if verb == true
        print(" compactness:", f);
    end
    return f
end

function numgrad_1D(func;x=[], N=400, δ = 1e-6)
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
        numerical_g[i,:] = (f2l-f2r+8*(f1r-f1l))
        x[i] = orig
    end
    numerical_g ./= (12*δ)
    return numerical_g;
end




function checkgrad_1D(func;x=[], N=400, δ = 1e-6)
    if x==[]
        x = abs.(rand(N))
    else
        N = length(x)
    end
    numerical_g = similar(x)
    analytic_g = similar(x)
    for i=1:N
        orig = x[i]
        x[i] = orig + 2*δ
        f2r = func(x,analytic_g)
        x[i] = orig + δ
        f1r = func(x,analytic_g)
        x[i] = orig - δ
        f1l = func(x,analytic_g)
        x[i] = orig - 2*δ
        f2l = func(x,analytic_g)
        numerical_g[i] = (f2l-f2r+8*(f1r-f1l))
        x[i] = orig
    end
    numerical_g ./= (12*δ)
    fref = func(x,analytic_g)
    numerator = norm(analytic_g-numerical_g)
    denominator = norm(analytic_g) + norm(numerical_g)
    difference = numerator / denominator
    print(difference, "\n")
    return (numerical_g, analytic_g);
end


function checkgrad_2D(func;x=[], N=36, M= 25, δ = 1e-6)
    if x==[]
        x = abs.(rand(N,M))
    else
        N,M = size(x)
    end
    numerical_g = similar(x)
    analytic_g = vec(similar(x))
    for i=1:N
        for j=1:M
            orig = x[i,j]
            x[i,j] = orig + δ
            f1r = func(x,analytic_g)
            x[i,j] = orig - δ
            f1l = func(x,analytic_g)
            x[i,j] = orig + 2*δ
            f2r = func(x,analytic_g)
            x[i,j] = orig - 2*δ
            f2l = func(x,analytic_g)
            x[i,j] = orig + 3*δ
            f3r = func(x,analytic_g)
            x[i,j] = orig - 3*δ
            f3l = func(x,analytic_g)
            numerical_g[i,j] = (f3r-f3l+9*(f2l-f2r)+45*(f1r-f1l))
            x[i,j] = orig
        end
    end
    numerical_g = vec(numerical_g)/(60*δ)
    fref = func(x,analytic_g)
    numerator = norm(analytic_g-numerical_g)
    denominator = norm(analytic_g) + norm(numerical_g)
    difference = numerator / denominator
    print(difference, "\n")
    return (numerical_g, analytic_g);
end

function trans_structnorm(x, g;verb=false, ϵ=1e-9)
    #this x is under the form (npix,nwavs)
    #but return the gradient as a 1D vector to use with Optimpack
    d = sqrt.(sum(x.^2, dims=2).+ϵ^2)
    f = sum(d.-ϵ);
    g[:] = vec(x./d);
    return f
end

function trans_tv(x, g; ζ=1e-13)
    # Transpectral or transtemporal regularization
    #this x is under the form (npix,nwavs)
    #but return the gradient as a 1D vector to use with Optimpack
    xlplus1_minus_xl = circshift(x, (0,-1))-x
    xl_minus_xminus1 = x-circshift(x, (0,1))
    d1 = sqrt.(xlplus1_minus_xl.^2 .+ ζ^2)
    d2 = sqrt.(xl_minus_xminus1.^2 .+ ζ^2)
    f = sum(d1.-ζ)
    g[:] = vec((-xlplus1_minus_xl./d1 + xl_minus_xminus1./d2))
    return f
end


function trans_tvsq(x, g;verb=false)
    # Transpectral or transtemporal regularization
    #this x is under the form (npix,nwavs)
    #but return the gradient as a 1D vector to use with Optimpack
    xlplus1_minus_xl = circshift(x, (0,-1))-x
    xl_minus_xminus1 = x-circshift(x, (0,1))
    d1 = xlplus1_minus_xl.^2
    d2 = xl_minus_xminus1.^2
    f = sum(d1)
    g[:] = 2*vec(-xlplus1_minus_xl + xl_minus_xminus1)
    return f
end

function trans_l1l2(x, g;verb=false, ζ=1e-13, δ=2.0)
    # Transpectral or transtemporal regularization using the L1L2 scheme
    #  x is under the form (npix,nwavs)
    #but return the gradient as a 1D vector to use with Optimpack
    xlplus1_minus_xl = circshift(x, (0,-1))-x
   # xlplus1_minus_xl[:,end] .= 0
    xl_minus_xminus1 = x-circshift(x, (0,1))
    # xl_minus_xminus1[:,1] .=0
    d1 = sqrt.(xlplus1_minus_xl.^2 .+ ζ^2)
    d2 = sqrt.(xl_minus_xminus1.^2 .+ ζ^2)
    f = sum(d1/δ - log.(1.0 .+d1/δ).-ζ)
    g[:] = vec((-xlplus1_minus_xl./d1 + xl_minus_xminus1./d2)/δ-((-xlplus1_minus_xl./d1)./(δ .+ d1) + (xl_minus_xminus1./d2)./(δ .+d2)))
    return f
end

function regularization(x, reg_g; printcolor = :black, regularizers=[], verb=true, goffset::Int64=1) # compound regularization
    reg_f = 0.0;
    for ireg =1:length(regularizers)
        temp_g = Array{Float64}(undef,length(x))
        if regularizers[ireg][1] == "centering"
            reg_f += regularizers[ireg][2]*reg_centering(x, temp_g, verb = verb)
        elseif regularizers[ireg][1] == "tv"
            reg_f += regularizers[ireg][2]*tv(x,temp_g, verb = verb)
        elseif regularizers[ireg][1] == "tvsq"
            reg_f += regularizers[ireg][2]*tvsq(x,temp_g, verb = verb)
        elseif regularizers[ireg][1] == "EPLL"
            reg_f += regularizers[ireg][2]*EPLL_fg(x,temp_g, regularizers[ireg][3])
        elseif regularizers[ireg][1] == "l1l2"
            reg_f += regularizers[ireg][2]*l1l2(x,temp_g, verb = verb, α = regularizers[ireg][3])
        elseif regularizers[ireg][1] == "l1l2w"
            reg_f += regularizers[ireg][2]*l1l2w(x,temp_g, verb = verb)
        elseif regularizers[ireg][1] == "l1hyp"
            reg_f += regularizers[ireg][2]*l1hyp(x,temp_g, verb = verb)
        elseif regularizers[ireg][1] == "l2sq"
            reg_f += regularizers[ireg][2]*l2sq(x,temp_g, verb = verb)
        elseif regularizers[ireg][1] == "compactness"
            reg_f += regularizers[ireg][2]*compactness(x,temp_g, verb = verb)
        elseif regularizers[ireg][1] == "entropy"
            reg_f += regularizers[ireg][2]*entropy(x,temp_g, verb = verb)
        end
        reg_g[goffset:end] += regularizers[ireg][2]*temp_g
    end
    if (verb==true)
        print("\n");
    end
    return reg_f
end


function chi2_dft_f(x::Array{Float64,1}, dft::Array{Complex{Float64},2}, data::OIdata; weights = [1.0,1.0,1.0], printcolor =:black, verb=true, vonmises=false)
    flux = sum(x);
    cvis_model = image_to_cvis_dft(x, dft);
    v2_model = cvis_to_v2(cvis_model, data.indx_v2);
    t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
    chi2_v2 = 0.0;
    chi2_t3amp = 0.0;
    chi2_t3phi = 0.0;
    if weights[1]>0
        chi2_v2 = norm((v2_model - data.v2)./data.v2_err)^2;
    end
    if weights[2]>0
        chi2_t3amp = norm((t3amp_model - data.t3amp)./data.t3amp_err)^2;
    end
    if weights[3]>0
        if vonmises == false
            chi2_t3phi = norm(mod360(t3phi_model - data.t3phi)./data.t3phi_err)^2;
        else
            chi2_t3phi =  sum(-2*data.t3phi_vonmises_err.*cos.((t3phi_model - data.t3phi)/180*pi) + data.t3phi_vonmises_chi2_offset)
        end
    end
    if verb==true
        printstyled("V2: $(chi2_v2/data.nv2) ", color=:red)
        printstyled("T3A: $(chi2_t3amp/data.nt3amp) ", color=:blue);
        printstyled("T3P: $(chi2_t3phi/data.nt3phi) ", color=:green);
        printstyled("Flux: $(flux) ", color=:black);
    end
    return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

function chi2_nfft_f(x::Array{Float64,1}, ftplan::Array{NFFT.NFFTPlan{2,0,Float64},1}, data::OIdata; weights = [1.0,1.0,1.0], cvis = [], printcolor =:black,  verb = false, vonmises=false)
    flux = sum(x);
    cvis_model = image_to_cvis_nfft(x, ftplan[1]);
    if length(cvis)>0
        #Note: cvis_model includes all the complex visibilities needed to compute V2, T3, etc.
        #      while the cvis variable is used to export visibility observables (e.g. diff vis or diff phi)
        cvis[:] = cvis_model[data.indx_vis]
    end
    v2_model = cvis_to_v2(cvis_model, data.indx_v2);
    t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
    chi2_v2 = 0.0;
    chi2_t3amp = 0.0;
    chi2_t3phi = 0.0;
    if weights[1]>0
        chi2_v2 = norm((v2_model - data.v2)./data.v2_err)^2;
    end

    if weights[2]>0
        chi2_t3amp = norm((t3amp_model - data.t3amp)./data.t3amp_err)^2;
    end
    if weights[3]>0
        if vonmises == false
            chi2_t3phi = norm(mod360(t3phi_model - data.t3phi)./data.t3phi_err)^2;
        else
            chi2_t3phi =  sum(-2*data.t3phi_vonmises_err.*cos.((t3phi_model - data.t3phi)/180*pi) + data.t3phi_vonmises_chi2_offset)
        end
    end
    if verb==true
        printstyled("V2: $(chi2_v2/data.nv2) ", color=:red)
        printstyled("T3A: $(chi2_t3amp/data.nt3amp) ", color=:blue);
        printstyled("T3P: $(chi2_t3phi/data.nt3phi) ", color=:green);
        printstyled("Flux: $(flux) ", color=:black);
    end
    return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

function chi2_dft_fg(x::Array{Float64,1}, g::Array{Float64,1}, dft::Array{Complex{Float64},2}, data::OIdata; weights = [1.0,1.0,1.0], printcolor =:black, verb=true, vonmises=false)
    flux = sum(x);
    cvis_model = image_to_cvis_dft(x, dft);
    v2_model = cvis_to_v2(cvis_model, data.indx_v2);
    t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
    chi2_v2 = 0.0;
    chi2_t3amp = 0.0;
    chi2_t3phi = 0.0;
    g_v2 = 0.0
    g_t3amp = 0.0
    g_t3phi = 0.0
    if (weights[1]>0)||(data.nv2==0)
        chi2_v2 = norm((v2_model - data.v2)./data.v2_err)^2;
        g_v2 = real(transpose(dft[data.indx_v2,:])*(4*((v2_model-data.v2)./data.v2_err.^2).*conj(cvis_model[data.indx_v2])));
    end

    if (weights[2]>0)||(data.nt3amp==0)
        chi2_t3amp = norm((t3amp_model - data.t3amp)./data.t3amp_err)^2;
        dT3 = 2.0*(t3amp_model-data.t3amp)./(data.t3amp_err.^2)
        g_t3amp = real(transpose(dft[data.indx_t3_1,:])*(dT3.*conj(cvis_model[data.indx_t3_1])./abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3])))+real(transpose(dft[data.indx_t3_2,:])*(dT3.*conj(cvis_model[data.indx_t3_2])./abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3]) ))+real(transpose(dft[data.indx_t3_3,:])*(dT3.*conj(cvis_model[data.indx_t3_3])./abs.(cvis_model[data.indx_t3_3]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]) ))
    end

    if (weights[3]>0)||(data.nt3phi==0)
        if vonmises == false
            chi2_t3phi = norm(mod360(t3phi_model - data.t3phi)./data.t3phi_err)^2;
            dT3 = 360.0/pi*mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2
            g_t3phi = imag(  transpose(dft[data.indx_t3_1,:])*(dT3./abs2.(cvis_model[data.indx_t3_1]).*conj(cvis_model[data.indx_t3_1]))+transpose(dft[data.indx_t3_2,:])*(dT3./abs2.(cvis_model[data.indx_t3_2]).*conj(cvis_model[data.indx_t3_2]))+transpose(dft[data.indx_t3_3,:])*(dT3./abs2.(cvis_model[data.indx_t3_3]).*conj(cvis_model[data.indx_t3_3])))
        else
            chi2_t3phi =  sum(-2*data.t3phi_vonmises_err.*cos.((t3phi_model - data.t3phi)/180*pi) + data.t3phi_vonmises_chi2_offset)
            dT3 = 2.0*data.t3phi_vonmises_err.*sin.((t3phi_model - data.t3phi)/180*pi)
            g_t3phi = imag(  transpose(dft[data.indx_t3_1,:])*(dT3./abs2.(cvis_model[data.indx_t3_1]).*conj(cvis_model[data.indx_t3_1]))+transpose(dft[data.indx_t3_2,:])*(dT3./abs2.(cvis_model[data.indx_t3_2]).*conj(cvis_model[data.indx_t3_2]))+transpose(dft[data.indx_t3_3,:])*(dT3./abs2.(cvis_model[data.indx_t3_3]).*conj(cvis_model[data.indx_t3_3])))
        end
    end
    g[:] = vec(weights[1]*g_v2 .+ weights[2]*g_t3amp .+ weights[3]*g_t3phi)
    # g[:] = (g - sum(vec(x).*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
    if verb==true
        printstyled("V2: $(chi2_v2/data.nv2) ", color=:red)
        printstyled("T3A: $(chi2_t3amp/data.nt3amp) ", color=:blue);
        printstyled("T3P: $(chi2_t3phi/data.nt3phi) ", color=:green);
        printstyled("Flux: $(flux) ", color=:black);
    end
    return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

function chi2_nfft_fg(x::Array{Float64,1}, g::Array{Float64,1}, ftplan::Array{NFFT.NFFTPlan{2,0,Float64},1}, data::OIdata; weights = [1.0,1.0,1.0], cvis = [], printcolor =:black,  verb = false, vonmises=false)
    flux = sum(x);
    cvis_model = image_to_cvis_nfft(x, ftplan[1]);
    if length(cvis)>0
        #Note: cvis_model includes all the complex visibilities needed to compute V2, T3, etc.
        #      while the cvis variable is used to export visibility observables (e.g. diff vis or diff phi)
        cvis[:] = cvis_model[data.indx_vis]
    end
    v2_model = cvis_to_v2(cvis_model, data.indx_v2);
    t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
    chi2_v2 = 0.0;
    chi2_t3amp = 0.0;
    chi2_t3phi = 0.0;
    g_v2 = 0.0
    g_t3amp = 0.0
    g_t3phi = 0.0
    if (weights[1]>0)||(data.nv2==0)
        chi2_v2 = norm((v2_model - data.v2)./data.v2_err)^2;
        g_v2 = real(nfft_adjoint(ftplan[3], (4*((v2_model-data.v2)./data.v2_err.^2).*cvis_model[data.indx_v2])));
    end

    if (weights[2]>0)||(data.nt3amp==0)
        chi2_t3amp = norm((t3amp_model - data.t3amp)./data.t3amp_err)^2;
        dT3 = 2.0*(t3amp_model-data.t3amp)./(data.t3amp_err.^2)
        g_t3amp = real(nfft_adjoint(ftplan[4], dT3.*cvis_model[data.indx_t3_1]./abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3]))) + real(nfft_adjoint(ftplan[5],dT3.*cvis_model[data.indx_t3_2]./abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3]))) + real(nfft_adjoint(ftplan[6],dT3.*cvis_model[data.indx_t3_3]./abs.(cvis_model[data.indx_t3_3]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2])))
    end

    if (weights[3]>0)||(data.nt3phi==0)
        if vonmises == false
            chi2_t3phi = norm(mod360(t3phi_model - data.t3phi)./data.t3phi_err)^2;
            dT3 = -360.0/pi*mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2
            g_t3phi = imag(nfft_adjoint(ftplan[4], dT3./abs2.(cvis_model[data.indx_t3_1]).*cvis_model[data.indx_t3_1])+nfft_adjoint(ftplan[5], dT3./abs2.(cvis_model[data.indx_t3_2]).*cvis_model[data.indx_t3_2])+nfft_adjoint(ftplan[6], dT3./abs2.(cvis_model[data.indx_t3_3]).*cvis_model[data.indx_t3_3]))
        else
            chi2_t3phi =  sum(-2*data.t3phi_vonmises_err.*cos.((t3phi_model - data.t3phi)/180*pi) + data.t3phi_vonmises_chi2_offset)
            dT3 = -2.0*data.t3phi_vonmises_err.*sin.((t3phi_model - data.t3phi)/180*pi)
            g_t3phi = imag( nfft_adjoint(ftplan[4], dT3./abs2.(cvis_model[data.indx_t3_1]).*cvis_model[data.indx_t3_1]) + nfft_adjoint(ftplan[5], dT3./abs2.(cvis_model[data.indx_t3_2]).*cvis_model[data.indx_t3_2]) + nfft_adjoint(ftplan[6], dT3./abs2.(cvis_model[data.indx_t3_3]).*cvis_model[data.indx_t3_3]))
        end
    end

    g[:] = vec(weights[1]*g_v2 .+ weights[2]*g_t3amp .+ weights[3]*g_t3phi)
    if verb==true
        printstyled("V2: $(chi2_v2/data.nv2) ", color=:red)
        printstyled("T3A: $(chi2_t3amp/data.nt3amp) ", color=:blue);
        printstyled("T3P: $(chi2_t3phi/data.nt3phi) ", color=:green);
        printstyled("Flux: $(flux) ", color=:black);
    end
    return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

function chi2_polychromatic_nfft_f(x::Array{Float64,1}, ft::Array{Array{NFFT.NFFTPlan{2,0,Float64},1},1}, data::Array{OIdata,1};weights = [1.0,1.0,1.0], printcolor= [], use_diffphases = false, verb = false)
    nwavs = length(ft);
    if printcolor == []
        printcolor=Array{Symbol}(undef,nwavs);
        printcolor[:] .= :black
    end
    npix = div(length(x),nwavs);
    cvis = fill((Complex{Float64}[]),nwavs);
    if use_diffphases == true
        for i=1:nwavs
            cvis[i] = Array{Complex{Float64},1}(undef, length(data[i].indx_vis))
        end
    end
    f = zeros(nwavs)
    for i=1:nwavs # weighted sum -- should probably do the computation in parallel
        if use_diffphases == true
            f[i] = chi2_nfft_f(x[1+(i-1)*npix:i*npix], ft[i], data[i], cvis=cvis[i], verb = verb, weights = weights);
        else
            f[i] = chi2_nfft_f(x[1+(i-1)*npix:i*npix], ft[i], data[i], verb = verb, weights = weights);
        end
        fr = f[i]/(data[i].nv2+data[i].nt3amp+data[i].nt3phi)
        if verb == true
            printstyled("Spectral channel $i chi2= $(f[i]) chi2r = $(fr)\n",color=printcolor[i]);
        end
    end
    chi2f = sum(f)
    ndof = sum([data[i].nv2+data[i].nt3amp+data[i].nt3phi for i=1:nwavs]);
    printstyled("All V2+T3AMP+T3PHI -- Chi2: $chi2f Chi2/dof: $(chi2f/ndof) \n", color=:red);
    # Differential phase
    if use_diffphases == true
        phi_model = angle.(hcat(cvis...))*180/pi
        diffphi_model = (nwavs*phi_model .- vec(sum(phi_model, dims=2))) /(nwavs-1)
        diffphi     = hcat([data[i].visphi for i=1:nwavs]...)
        diffphi_err = hcat([data[i].visphi_err for i=1:nwavs]...)
        chi2_diffphi = norm(mod360(diffphi_model-diffphi)./diffphi_err)^2
        chi2f += chi2_diffphi;
        print("Differential phase chi2r: $(chi2_diffphi/length(diffphi_model)) \n");
    end
    return chi2f;
end



function crit_dft_fg(x, g, dft, data; regularizers=[], weights = [1.0,1.0,1.0], verb = true) # regularization + negloglikelihood
    chi2_f = chi2_dft_fg(x,  g, dft, data, verb = verb, weights = weights ); # this resets g
    reg_f = regularization(x, g, regularizers=regularizers, verb = verb) ;# this adds to g
    flux = sum(x)
    g[:] = (g .- sum(vec(x).*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
    return chi2_f + reg_f;
end

function crit_nfft_fg(x::Array{Float64,1},g::Array{Float64,1}, fftplan::Array{NFFT.NFFTPlan{2,0,Float64},1}, data::OIdata; weights = [1.0,1.0,1.0], cvis = [], printcolor = :black, regularizers=[], verb = true)
    chi2_f = chi2_nfft_fg(x, g, fftplan, data, cvis = cvis, verb = verb, weights = weights);
    reg_f = regularization(x,g, regularizers=regularizers, printcolor = printcolor, verb = verb);
    flux = sum(x)
    g[:] = (g .- sum(vec(x).*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
    return chi2_f + reg_f;
end

function crit_nfft_f(x::Array{Float64,1}, fftplan::Array{NFFT.NFFTPlan{2,0,Float64},1}, data::OIdata; weights = [1.0,1.0,1.0], cvis = [], printcolor = :black, regularizers=[], verb = true)
    chi2_f = chi2_nfft_fg(x, g, fftplan, data, cvis = cvis, verb = verb, weights = weights );
    reg_f = regularization(x,g, regularizers=regularizers, printcolor = printcolor, verb = verb);
    return chi2_f + reg_f;
end

function crit_multitemporal_nfft_fg(x::Array{Float64,1}, g::Array{Float64,1}, ft::Array{Array{NFFT.NFFTPlan{2,0,Float64},1},1}, data::Array{OIdata,1};weights = [1.0,1.0,1.0], printcolor= [], epochs_weights=[],regularizers=[], verb = false)
    nepochs = length(ft);
    if epochs_weights == []
        epochs_weights=ones(Float64, nepochs);
    end
    if printcolor == []
        printcolor=Array{Symbol}(undef,nepochs);
        printcolor[:] .= :black
    end
    npix = div(length(x),nepochs);
    f = 0.0;
    for i=1:nepochs # weighted sum -- should probably do the computation in parallel
        tslice = 1+(i-1)*npix:i*npix; # temporal slice
        subg = Array{Float64}(undef, npix);
        printstyled("Epoch $i ",color=printcolor[i]);
        f += epochs_weights[i]*crit_nfft_fg(x[tslice], subg, ft[i], data[i], regularizers=regularizers[i], printcolor = printcolor[i], verb = verb, weights = weights);
        g[tslice] = epochs_weights[i]*subg
    end

    # cross temporal regularization
    if length(regularizers)>nepochs
        if (regularizers[nepochs+1][1][1] == "temporal_tvsq")  & (nepochs>1)
            y = reshape(x,(npix,nepochs))
            temporalf = sum( (y[:,2:end]-y[:,1:end-1]).^2 )
            tv_g = Array{Float64}(undef, npix,nepochs)
            if nepochs>2
                tv_g[:,1] = 2*(y[:,1] - y[:,2])
                tv_g[:,2:end-1] = 4*y[:,2:end-1]-2*(y[:,1:end-2]+y[:,3:end])
                tv_g[:,end] = 2*(y[:,end] - y[:,end-1])
            else
                tv_g[:,1] = 2*(y[:,1]-y[:,2]);
                tv_g[:,2] = 2*(y[:,2]-y[:,1]);
            end
            f+= regularizers[nepochs+1][1][2]*temporalf
            g[:] += regularizers[nepochs+1][1][2]*vec(tv_g);
            printstyled("Temporal regularization: $(regularizers[nepochs+1][1][2]*temporalf)\n", color=:yellow)
        end
    end
    return f;
end


function crit_polychromatic_nfft_fg(x::Array{Float64,1}, g::Array{Float64,1}, ft::Array{Array{NFFT.NFFTPlan{2,0,Float64},1},1}, data::Array{OIdata,1};weights = [1.0,1.0,1.0], printcolor= [], regularizers=[], use_diffphases = false, verb = false)
    nwavs = length(ft);
    if printcolor == []
        printcolor=Array{Symbol}(undef,nwavs);
        printcolor[:] .= :black
    end
    npix = div(length(x),nwavs);
    cvis = fill((Complex{Float64}[]),nwavs);
    if use_diffphases == true
        for i=1:nwavs
            cvis[i] = Array{Complex{Float64},1}(undef, length(data[i].indx_vis))
        end
    end

    f = 0.0;
    for i=1:nwavs # weighted sum -- should probably do the computation in parallel
        tslice = 1+(i-1)*npix:i*npix; # chromatic slice #TODO: use reshape instead ?
        subg = Array{Float64}(undef, npix);
        if verb == true
            printstyled("Spectral channel $i ",color=printcolor[i]);
        end
        f += crit_nfft_fg(x[tslice], subg, ft[i], data[i], regularizers=regularizers[i], cvis = cvis[i], printcolor = printcolor[i], verb = verb, weights = weights);
        g[tslice] = subg
    end
    ndof = sum([data[i].nv2+data[i].nt3amp+data[i].nt3phi for i=1:nwavs]);
    printstyled("Indpt images -- Crit: $f Crit/dof: $(f/ndof) \n", color=:red);
    # Differential phase
    #  if data.nvisphi > 0
    # Compute vis_ref
    if use_diffphases == true
        phi_model = angle.(hcat(cvis...))*180/pi
        diffphi_model = (nwavs*phi_model .- vec(sum(phi_model, dims=2))) /(nwavs-1)
        diffphi     = hcat([data[i].visphi for i=1:nwavs]...)
        diffphi_err = hcat([data[i].visphi_err for i=1:nwavs]...)
        chi2_diffphi = norm(mod360(diffphi_model-diffphi)./diffphi_err)^2
        f += chi2_diffphi;
        print("Differential phase chi2r: $(chi2_diffphi/length(diffphi_model)) \n");
        # Compute differential phase gradient
        d_diffphi = zeros(npix,nwavs)
        for i=1:nwavs
            d_diffphi[:,i] = -360.0/pi*vec(imag.(nfft_adjoint(ft[i][2], ((mod360(diffphi_model[:,i]-diffphi[:,i])./diffphi_err[:,i].^2)./abs2.(cvis[i])).*cvis[i])))
        end
        g_diffphi = vec(nwavs*d_diffphi .- vec(sum(d_diffphi, dims=2)))/(nwavs-1)
        #d_diffphi = [] # free local memory immediately after use
        g[:] += g_diffphi
    end

    # transspectral regularization
    if length(regularizers)>nwavs
        ntransreg = length(regularizers[nwavs+1]);
        y = reshape(x,(npix,nwavs))
        tg = Array{Float64}(undef, npix*nwavs)
        for i=1:ntransreg
            if (regularizers[nwavs+1][i][1] == "transspectral_tv")
                tf = trans_tv(y,tg)
                f    += regularizers[nwavs+1][i][2]*tf
                g[:] += regularizers[nwavs+1][i][2]*tg;
                printstyled("Trans-spectral TV: $(regularizers[nwavs+1][i][2]*tf)\n", color=:yellow)
            end
            if (regularizers[nwavs+1][i][1] == "transspectral_tvsq")
                tf = trans_tvsq(y, tg)
                f+= regularizers[nwavs+1][i][2]*tf
                g[:] += regularizers[nwavs+1][i][2]*tg
                printstyled("Trans-spectral squared TV: $(regularizers[nwavs+1][i][2]*tf)\n", color=:yellow)
            end
            if (regularizers[nwavs+1][i][1] == "transspectral_structnorm")
                tf = trans_structnorm(y,tg)
                f    += regularizers[nwavs+1][i][2]*tf
                g[:] += regularizers[nwavs+1][i][2]*tg;
                printstyled("Trans-spectral Structured Norm: $(regularizers[nwavs+1][i][2]*tf)\n", color=:yellow)
            end
            if (regularizers[nwavs+1][i][1] == "transspectral_l1l2")
                tf = trans_l1l2(y, tg, δ=regularizers[nwavs+1][i][3])
                f+= regularizers[nwavs+1][i][2]*tf
                g[:] += regularizers[nwavs+1][i][2]*tg
                printstyled("Trans-spectral l1l2 norm: $(regularizers[nwavs+1][i][2]*tf)\n", color=:yellow)
            end
        end
    end
    printstyled("Post trans -- Crit: $f Crit/dof: $(f/ndof) \n", color=:blue);
    return f;
end


function image_to_cvis_polychromatic_nfft(x::Array{Float64,1}, ft::Array{Array{NFFT.NFFTPlan{2,0,Float64},1},1})
    nwavs = length(ft);
    npix = div(length(x),nwavs);
    cvis = fill((Complex{Float64}[]),nwavs);
    for i=1:nwavs
        cvis[i] = image_to_cvis_nfft(x[1+(i-1)*npix:i*npix], ft[i]);
    end
    return cvis;
end

using OptimPackNextGen
function reconstruct(x_start::Array{Float64,1}, data::OIdata, ft; weights = [1.0,1.0,1.0], printcolor = :black, verb = false, maxiter = 100, regularizers =[], gtol=(0,1e-8))
    x_sol = []
    if typeof(ft) == Array{NFFT.NFFTPlan{2,0,Float64},1}
        crit = (x,g)->crit_nfft_fg(x, g, ft, data, regularizers=regularizers, verb = verb , weights = weights)
        x_sol = OptimPackNextGen.vmlmb(crit, x_start, verb=verb, lower=0, maxiter=maxiter, blmvm=false, gtol=gtol);
    else
        crit = (x,g)->crit_dft_fg(x, g, ft, data, regularizers=regularizers, verb = verb, weights = weights)
        x_sol = OptimPackNextGen.vmlmb(crit, x_start, verb=verb, lower=0, maxiter=maxiter, blmvm=false, gtol=gtol);
    end
    return x_sol
end


function reconstruct_multitemporal(x_start::Array{Float64,1}, data::Array{OIdata, 1}, ft; weights = [1.0,1.0,1.0], epochs_weights =[], printcolor= [], verb = true, maxiter = 100, regularizers =[], gtol=(0,1e-8))
    x_sol = []
    if typeof(ft) == Array{Array{NFFT.NFFTPlan{2,0,Float64},1},1}
        crit = (x,g)->crit_multitemporal_nfft_fg(x, g, ft, data, printcolor=printcolor, weights = weights, epochs_weights=epochs_weights, regularizers=regularizers, verb = verb)
        x_sol = OptimPackNextGen.vmlmb(crit, x_start, verb=verb, lower=0, maxiter=maxiter, blmvm=false, gtol=gtol);
    else
        error("Sorry, polytemporal DFT methods not implemented yet");
    end
    return x_sol
end



function reconstruct_polychromatic(x_start::Array{Float64,1}, data::Array{OIdata,1}, ft; weights = [1.0,1.0,1.0], printcolor= [], verb = true, use_diffphases = false, maxiter = 100, regularizers =[])
    x_sol = []
    if regularizers == []
        regularizers = fill([],length(data))
    end

    if typeof(ft) == Array{Array{NFFT.NFFTPlan{2,0,Float64},1},1}
        crit = (x,g)->crit_polychromatic_nfft_fg(x, g, ft, data, weights = weights, printcolor=printcolor, regularizers=regularizers, use_diffphases = use_diffphases, verb = verb)
        x_sol = OptimPackNextGen.vmlmb(crit, x_start, verb=true, lower=0, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
    else
        error("Sorry, polychromatic DFT methods not implemented yet");
    end
    return x_sol
end


function chi2_sparco_nfft_f(x::Array{Float64,1}, ftplan::Array{NFFT.NFFTPlan{2,0,Float64},1}, data::OIdata, params::Array{Float64,1}; verb = true, weights = [1.0,1.0,1.0] ) # criterion function for nfft
    # Image x is of length = N*N
    # Nparams are passed as an array

    # The chromatism is defined as follows
    #        f_star_0 * (lambda/ lambda_0)^-4 * V_star + (1-f_star_0 - f_bg_0 )*(lambda/lambda_0)^d_ind * V_env
    # V_tot = -------------------------------------------------------------------------------------------------
    #        (f_star_0 + f_bg_0) (lambda/ lambda_0)^-4 + (1 - f_star_λ0 - f_bg_0 )*(lambda/lambda_0)^d_ind
    # param[1] = fs0
    # param[2] = fbg0
    # param[3] = diameter of star = 2.776e-01
    # param[4] : fixed, d_ind environment power law
    # param[5] :  lambda_0 (fixed) = 1.65e-06
    # params=[0.8, 0.1, 0.5, -4.0, 1.65e-6]
    # Compute visibilty for model + image
    λ0 = params[5];
    λ = data.uv_lam
    α = (λ/λ0).^-4.0;
    β = (λ/λ0).^params[4];
    fluxstar = params[1]*α;
    fluxbg = params[2]*α;
    fluxenv = (1.0-params[1]-params[2])*β;
    Vstar = visibility_ud([params[3]], data.uv);
    Venv = image_to_cvis_nfft(x, ftplan[1]);
    cvis_model = (fluxstar.*Vstar + fluxenv.*Venv)./(fluxstar+fluxenv+fluxbg)
    v2_model = cvis_to_v2(cvis_model, data.indx_v2);
    t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);

    return chi2_v2 + chi2_t3amp + chi2_t3phi
end
# TBD: merge with previous function
function chi2_sparco_nfft_f(x::Array{Float64,1}, ftplan::Array{NFFT.NFFTPlan{2,0,Float64},1}, data::OIdata, nparams::Int64; verb = true ) # criterion function for nfft
    # x is of length = N*N+Nparams
    params=x[1:nparams]  # extract parameters
    λ0 = params[5];
    λ = data.uv_lam
    α = (λ/λ0).^-4.0;
    β = (λ/λ0).^params[4];
    fluxstar = params[1]*α;
    fluxbg = params[2]*α;
    fluxenv = (1.0-params[1]-params[2])*β;
    Vstar = visibility_ud([params[3]], data.uv);
    Venv = image_to_cvis_nfft(x[nparams+1:end], ftplan[1]);
    cvis_model = (fluxstar.*Vstar + fluxenv.*Venv)./(fluxstar+fluxenv+fluxbg)
    v2_model = cvis_to_v2(cvis_model, data.indx_v2);
    t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
    chi2_v2 = 0.0;
    chi2_t3amp = 0.0;
    chi2_t3phi = 0.0;
    if weights[1]>0
        chi2_v2 = norm((v2_model - data.v2)./data.v2_err)^2;
    end

    if weights[2]>0
        chi2_t3amp = norm((t3amp_model - data.t3amp)./data.t3amp_err)^2;
    end
    if weights[3]>0
            chi2_t3phi = norm(mod360(t3phi_model - data.t3phi)./data.t3phi_err)^2;
    end
    if verb==true
        printstyled("V2: $(chi2_v2/data.nv2) ", color=:red)
        printstyled("T3A: $(chi2_t3amp/data.nt3amp) ", color=:blue);
        printstyled("T3P: $(chi2_t3phi/data.nt3phi) ", color=:green);
    end
    return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

function chi2_sparco_nfft_fg(x::Array{Float64,1},  g::Array{Float64,1}, ftplan::Array{NFFT.NFFTPlan{2,0,Float64},1}, data::OIdata, nparams::Int64; verb = true, weights=[1.0,1.0,1.0] ) # criterion function for nfft
    params=x[1:nparams]  # extract parameters
    λ0 = params[5];
    λ = data.uv_lam
    α = (λ/λ0).^-4.0;
    β = (λ/λ0).^params[4];
    fluxstar = params[1]*α;
    fluxbg = params[2]*α;
    fluxenv = (1.0-params[1]-params[2])*β;
    Vstar = visibility_ud([params[3]], data.uv);
    Venv = image_to_cvis_nfft(x[nparams+1:end], ftplan[1]);

    cvis_model = (fluxstar.*Vstar + fluxenv.*Venv)./(fluxstar+fluxenv+fluxbg)
    v2_model = cvis_to_v2(cvis_model, data.indx_v2);
    t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);

    chi2_v2 = 0.0;
    chi2_t3amp = 0.0;
    chi2_t3phi = 0.0;
    if weights[1]>0
        chi2_v2 = norm((v2_model - data.v2)./data.v2_err)^2;
    end

    if weights[2]>0
        chi2_t3amp = norm((t3amp_model - data.t3amp)./data.t3amp_err)^2;
    end
    if weights[3]>0
            chi2_t3phi = norm(mod360(t3phi_model - data.t3phi)./data.t3phi_err)^2;
    end
    if verb==true
        printstyled("V2: $(chi2_v2/data.nv2) ", color=:red)
        printstyled("T3A: $(chi2_t3amp/data.nt3amp) ", color=:blue);
        printstyled("T3P: $(chi2_t3phi/data.nt3phi) ", color=:green);
    end

    # Derivative with respect to fs0 (param[1])
    u =  fluxstar.*Vstar + fluxenv.*Venv;
    v =  fluxstar + fluxenv + fluxbg;
    du = α.*Vstar - β.*Venv # du/dfs0
    dv = α - β# dv/dfs0
    dcvis_model_dfs0 = (du.*v-u.*dv)./(v.*v)

    # Derivative with respect to fg0 (param[2])
    du = -β.*Venv # du/dfg
    dv = α-β
    dcvis_model_dfg0 = (du.*v-u.*dv)./(v.*v)

    # Derivative with respect to diameter D (param[3])
    dVstar = dvisibility_ud([params[3]], data.uv);
    du = fluxstar.*dVstar
    dv = 0;
    dcvis_model_dD = (du.*v-u.*dv)./(v.*v)

    # Derivative with respect to spectral index
    du = log.(λ/λ0).*fluxenv.*Venv
    dv = log.(λ/λ0).*fluxenv
    dcvis_model_dindx = (du.*v-u.*dv)./(v.*v)

    # Fill up gradient of parameters
    g[1] = Δchi2(dcvis_model_dfs0, cvis_model, v2_model, t3_model, t3amp_model, t3phi_model, data);
    g[2] = Δchi2(dcvis_model_dfg0, cvis_model, v2_model, t3_model, t3amp_model, t3phi_model, data);
    g[3] = Δchi2(dcvis_model_dD, cvis_model, v2_model, t3_model, t3amp_model, t3phi_model, data);
    g[4] = Δchi2(dcvis_model_dindx, cvis_model, v2_model, t3_model, t3amp_model, t3phi_model, data);
    g[5] = 0.0;

    # Gradient with respect to pixel fluxes
    imratio = fluxenv./(fluxstar+fluxenv+fluxbg)

    g_v2 = real(nfft_adjoint(ftplan[3], (4*((v2_model-data.v2)./data.v2_err.^2).*cvis_model[data.indx_v2].*imratio[data.indx_v2])));
    g_t3amp = real(nfft_adjoint(ftplan[4], (2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*cvis_model[data.indx_t3_1].*imratio[data.indx_t3_1]./abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3]) )))
    + real(nfft_adjoint(ftplan[5], (2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*cvis_model[data.indx_t3_2].*imratio[data.indx_t3_2]./abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3]) )))
    + real(nfft_adjoint(ftplan[6], (2.0*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*cvis_model[data.indx_t3_3].*imratio[data.indx_t3_3]./abs.(cvis_model[data.indx_t3_3]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]) )))

    #dt3 = dcvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3] + cvis_model[data.indx_t3_1].*dcvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3] + cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2].*dcvis_model[data.indx_t3_3]
    #dt3phi = 360.0/pi*sum( (mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2).*imag(conj(t3).*dt3)./abs2.(t3))
    g_t3phi = -360.0/pi*imag(nfft_adjoint(ftplan[4], ((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*imratio[data.indx_t3_1].*conj(cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3]).*t3_model)
    + nfft_adjoint(ftplan[5], ((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*imratio[data.indx_t3_2].*conj(cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_3]).*t3_model)
    + nfft_adjoint(ftplan[6], ((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*imratio[data.indx_t3_3].*conj(cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2]).*t3_model))

    g[nparams+1:end] = vec(weights[1]*g_v2 .+ weights[2]*g_t3amp .+ weights[3]*g_t3phi)
    #Normalization done later
    #flux = sum(x[nparams+1:end])
    #g[nparams+1:end] = (g[nparams+1:end] .- sum(vec(x[nparams+1:end]).*g[nparams+1:end]) / flux ) / flux; # gradient correction to take into account the non-normalized image

    return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

function Δchi2(dcvis_model::Union{Array{Complex{Float64},1}, Array{Float64,1}}, cvis_model::Union{Array{Complex{Float64},1},Array{Float64,1}}, v2_model::Array{Float64,1}, t3_model::Array{Complex{Float64},1}, t3amp_model::Array{Float64,1}, t3phi_model::Array{Float64,1}, data::OIdata) # return gradient of chi2 when cvis and dvis_model are known
    dv2 = 4*sum( ( (v2_model-data.v2)./data.v2_err.^2 ).*real.(cvis_model[data.indx_v2].*dcvis_model[data.indx_v2]))
    dt3amp = 2.0*sum(((t3amp_model-data.t3amp)./data.t3amp_err.^2).* ( real(cvis_model[data.indx_t3_1].*conj(dcvis_model[data.indx_t3_1]))./abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3])
    + real(cvis_model[data.indx_t3_2].*conj(dcvis_model[data.indx_t3_2]))./abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3])
    + real(cvis_model[data.indx_t3_3].*conj(dcvis_model[data.indx_t3_3]))./abs.(cvis_model[data.indx_t3_3]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2])))
    t3 = cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3]
    dt3 = dcvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3] + cvis_model[data.indx_t3_1].*dcvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3] + cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2].*dcvis_model[data.indx_t3_3]
    dt3phi = 360.0/pi*sum( (mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2).*imag(conj(t3).*dt3)./abs2.(t3))
    # println("dv2:", dv2, " dt3amp:", dt3amp, " dt3phi:", dt3phi);
    return dv2+dt3amp+dt3phi
end

function crit_sparco_nfft_fg(x::Array{Float64,1},g::Array{Float64,1}, ftplan::Array{NFFT.NFFTPlan{2,0,Float64},1}, data::OIdata, nparams::Int64; weights = [1.0,1.0,1.0], cvis = [], printcolor = :black, regularizers=[], verb = true)
    chi2_f = chi2_sparco_nfft_fg(x,  g, ftplan, data, nparams, weights=weights)
    reg_f = regularization(x[nparams+1:end], g, regularizers=regularizers, printcolor = printcolor, verb = verb, goffset=nparams+1);
    # Gradient correction for the image (parameters are left untouched)
    flux = sum(x[nparams+1:end])
    g[nparams+1:end] = (g[nparams+1:end] .- sum(vec(x[nparams+1:end]).*g[nparams+1:end]) / flux ) / flux; # gradient correction to take into account the non-normalized image
    return chi2_f + reg_f;
end

function reconstruct_sparco_gray(x_start::Array{Float64,1}, params_start::Array{Float64,1}, data::OIdata, ft; printcolor = :black, verb = false, maxiter = 100, regularizers =[], weights=[1.0,1.0,1.0]) #grey environment
    x_sol = []
    crit = (x,g)->crit_sparco_nfft_fg(x, g, ft, data, length(params_start), regularizers =regularizers, verb = verb, weights=weights)
    sol = OptimPackNextGen.vmlmb(crit, [params_start;x_start], verb=verb, lower=0, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
    return (sol[1:length(params_start)], sol[length(params_start)+1:end])
end

# if Pkg.installed("Wavelets") !=nothing
#   using Wavelets
#   function W(mat; wavelet_bases=[WT.db1, WT.db2, WT.db3, WT.db4, WT.db5, WT.db6, WT.db7, WT.db8, WT.haar])
#         nx = Int(sqrt(length(mat)))
#       Wx = Array{Float64}(nx*nx,length(wavelet_bases));
#         for i=1:length(wavelet_bases)
#             Wx[:, i]=vec(dwt(reshape(mat,nx,nx), wavelet(wavelet_bases[i])));
#         end
#         return Wx;
#   end

#   function Wt(mat; wavelet_bases=[WT.db1, WT.db2, WT.db3, WT.db4, WT.db5, WT.db6, WT.db7, WT.db8, WT.haar])
#     nx = size(mat,2)
#     IWx = Array{Float64}(nx*nx,length(wavelet_bases));
#      for i=1:length(wavelet_bases)
#          IWx[:,i] = vec(idwt(reshape(mat[:,i],(nx,nx)), wavelet(wavelet_bases[i])));
#       end
#    return sum(IWx,2);#/length(wavelet_bases);
#   end


#   function regwav(x,wav_g; wavelet_bases=[WT.db1, WT.db2, WT.db3, WT.db4, WT.db5, WT.db6, WT.db7, WT.db8, WT.haar])
#     wav_f = norm(W(x,wavelet_bases))^2
#     wav_g[:] = 2*length(wavelet_bases)*x
#     return tv_f
#   end
#   end
