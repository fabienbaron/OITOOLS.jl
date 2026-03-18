using Statistics,LinearAlgebra, SparseArrays, SpecialFunctions, NFFT, Match, Printf



function Base.display(ft::Vector{<:AbstractArray{<:NFFTPlan}})
    nmulti = length(ft)
    println("OITOOLS NFFT transform for $nmulti channels");
    println("Access e.g. ft[1] if wanting to check an individual channel NFFT");
end

# setup_dft: generic in the floating-point type T of the uv coordinates.
# Returns a Matrix{Complex{T}} DFT kernel.
# Convention: V(u,v) = Σ I(x,y) exp(-2πi(u·x + v·y))  (astronomical standard)
function setup_dft(uv::Matrix{T}, nx, pixsize) where T<:AbstractFloat
    scale_rad = T(pixsize * (pi / 180.0) / 3600000.0)
    nuv = size(uv, 2)
    span = T.((2 * (1:nx) .- (nx + 1)) .* (scale_rad * pi))
    xvals = reshape(span, 1, nx, 1)
    yvals = reshape(span, 1, 1, nx)
    dft = Complex{T}.(reshape(cis.(.-((uv[1, :] .* xvals) .+ (uv[2, :] .* yvals))), nuv, nx^2))
    return dft
end

# Convenience: dispatch on OIdata{T} so the eltype flows through automatically.
function setup_dft(data::OIdata{T}, nx, pixsize) where T
    return setup_dft(data.uv, nx, pixsize)
end

# function setup_dft(data::Matrix{OIdata}, nx, pixsize)
#     if size(data) == (1,1)
#         return setup_dft(data[1,1].uv, nx, pixsize);
#     else 
#         error("Multidimensional DFT not implemented yet");
#     end
# end

# function setup_nfft(data::Matrix{OIdata}, nx, pixsize)
#     if size(data) == (1,1)
#         return setup_nfft(data[1,1], nx, pixsize);
#     else 
#         return setup_nfft_multi(data, nx, pixsize);
#     end
# end

# setup_nfft: dispatches on OIdata{T}; plan_nfft infers T from Matrix{T} uv coords,
# returning NFFTPlan{T,2,1}. Works for both Float32 and Float64.
# Convention: NFFT kernel is exp(-2πi k·x), so passing (u,v) directly gives
# V(u,v) = Σ I(x,y) exp(-2πi(u·x + v·y))  (astronomical standard)
function setup_nfft(data::OIdata{T}, nx, pixsize) where T
    scale_rad = T(pixsize * (pi / 180.0) / 3600000.0) .* data.uv
    fftplan_uv  = plan_nfft(scale_rad, (nx,nx), m=4, σ=2.0)
    fftplan_vis  = plan_nfft(scale_rad[:, data.indx_vis], (nx,nx), m=4, σ=2.0)
    fftplan_v2   = plan_nfft(scale_rad[:, data.indx_v2],  (nx,nx), m=4, σ=2.0)
    fftplan_t3_1 = plan_nfft(scale_rad[:, data.indx_t3_1],(nx,nx), m=4, σ=2.0)
    fftplan_t3_2 = plan_nfft(scale_rad[:, data.indx_t3_2],(nx,nx), m=4, σ=2.0)
    fftplan_t3_3 = plan_nfft(scale_rad[:, data.indx_t3_3],(nx,nx), m=4, σ=2.0)
    return [fftplan_uv, fftplan_vis, fftplan_v2, fftplan_t3_1, fftplan_t3_2, fftplan_t3_3]
end

# Raw-uv overload: caller must pass Matrix{T} explicitly typed.
function setup_nfft(uv::Matrix{T}, indx_vis, indx_v2, indx_t3_1, indx_t3_2, indx_t3_3, nx, pixsize) where T<:AbstractFloat
    scale_rad = T(pixsize * (pi / 180.0) / 3600000.0) .* uv
    fftplan_uv  = plan_nfft(scale_rad, (nx,nx), m=4, σ=2.0)
    fftplan_vis  = plan_nfft(scale_rad[:, indx_vis], (nx,nx), m=4, σ=2.0)
    fftplan_v2   = plan_nfft(scale_rad[:, indx_v2],  (nx,nx), m=4, σ=2.0)
    fftplan_t3_1 = plan_nfft(scale_rad[:, indx_t3_1],(nx,nx), m=4, σ=2.0)
    fftplan_t3_2 = plan_nfft(scale_rad[:, indx_t3_2],(nx,nx), m=4, σ=2.0)
    fftplan_t3_3 = plan_nfft(scale_rad[:, indx_t3_3],(nx,nx), m=4, σ=2.0)
    return [fftplan_uv, fftplan_vis, fftplan_v2, fftplan_t3_1, fftplan_t3_2, fftplan_t3_3]
end

# function setup_nfft_multi(data, nx, pixsize)
#     nwavs = size(data,1);
#     nepochs = size(data,2);
#     scale_rad = pixsize * (pi / 180.0) / 3600000.0;
#     fftplan_multi = Array{Array{NFFT.NFFTPlan{Float64, 2, 1}},2}(undef, nwavs, nepochs);
#     for i=1:nepochs
#         for j=1:nwavs
#             fftplan_multi[j,i]=setup_nfft(data[j,i], nx, pixsize);
#         end
#     end
#     return fftplan_multi
# end

function setup_nfft_multiepochs(data::AbstractVector{<:OIdata}, nx, pixsize)
    nepochs = length(data)
    fftplan_multi = Vector{Vector{NFFTPlan}}(undef, nepochs)
    for i in 1:nepochs
        fftplan_multi[i] = setup_nfft(data[i], nx, pixsize)
    end
    return fftplan_multi
end

function setup_nfft_polychromatic(data::AbstractVecOrMat{<:OIdata}, nx, pixsize)
    nwavs = size(data, 1)
    fftplan_multi = Vector{Vector{NFFTPlan}}(undef, nwavs)
    for i in 1:nwavs
        fftplan_multi[i] = setup_nfft(data[i], nx, pixsize)
    end
    return fftplan_multi
end

function setup_dft_polychromatic(data::AbstractVecOrMat{<:OIdata{T}}, nx, pixsize) where T
    nwavs = size(data, 1)
    fftplan_multi = Vector{Matrix{Complex{T}}}(undef, nwavs)
    for i in 1:nwavs
        fftplan_multi[i] = setup_dft(data[i], nx, pixsize)
    end
    return fftplan_multi
end


# function setup_nfft_polychromatic(uv::Array{Float64,3}, indx_vis, indx_v2, indx_t3_1, indx_t3_2,indx_t3_3, nx, pixsize)
#     nwavs = size(data,1);
#     scale_rad = pixsize * (pi / 180.0) / 3600000.0;
#     fftplan_multi = Array{Array{NFFT.NFFTPlan{Float64, 2, 1}}}(undef, nwavs);
#     for i=1:nwavs
#         fftplan_multi[i]=setup_nfft(uv[:,:,i], nx, pixsize);
#     end
#     return fftplan_multi
# end

function mod360(x)
    mod.(mod.(x.+180.0,360.0).+360.0, 360.0) .- 180.0
end

function vis_to_v2(cvis, indx)
    v2_model = abs2.(cvis[indx]);
end

function vis_to_t3(cvis, indx1, indx2, indx3)
    t3 = cvis[indx1].*cvis[indx2].*cvis[indx3];
    t3amp = abs.(t3);
    t3phi = angle.(t3)*180.0/pi;
    return t3, t3amp, t3phi
end

function observables(x, ft, data)
return image_to_obs(x, ft, data)
end

function image_to_obs(x, ft, data)
    cvis_model = image_to_vis(x, ft);
    v2_model = vis_to_v2(cvis_model, data.indx_v2);
    _, t3amp_model, t3phi_model = vis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
    return v2_model, t3amp_model, t3phi_model
end

# image_to_vis: generic overloads.
# T is inferred from the image array eltype; Complex{T} is used consistently.
function image_to_vis(x::AbstractMatrix{T}, dft::Matrix{Complex{T}}) where T<:AbstractFloat
    cvis_model = dft * vec(x / sum(x))
end

function image_to_vis(x::AbstractVector{T}, nfft_plan::NFFT.NFFTPlan{T}) where T<:AbstractFloat
    nx = Int64(sqrt(length(x)))
    cvis_model = nfft_plan * Complex{T}.(reshape(x / sum(x), (nx, nx)))
end

function image_to_vis(x::AbstractMatrix{T}, nfft_plan::NFFT.NFFTPlan{T}) where T<:AbstractFloat
    cvis_model = nfft_plan * Complex{T}.(x / sum(x))
end

function image_to_vis(x::AbstractVector{T}, nfft_plan::AbstractVector{<:NFFT.NFFTPlan{T}}) where T<:AbstractFloat
    nx = Int64(sqrt(length(x)))
    cvis_model = nfft_plan[1] * Complex{T}.(reshape(x / sum(x), (nx, nx)))
end

function image_to_vis(x::AbstractMatrix{T}, nfft_plan::AbstractVector{<:NFFT.NFFTPlan{T}}) where T<:AbstractFloat
    cvis_model = nfft_plan[1] * Complex{T}.(x / sum(x))
end

function image_to_v2(x, data, ft::AbstractVector{<:NFFT.NFFTPlan})
   return vis_to_v2(image_to_vis(x, ft), data.indx_v2)
end

function image_to_t3phi(x, data, ft::AbstractVector{<:NFFT.NFFTPlan})
    _, _, t3phi_model = vis_to_t3(image_to_vis(x, ft), data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3)
    return t3phi_model
 end
 
function image_to_t3amp(x, data, ft::AbstractVector{<:NFFT.NFFTPlan})
    _, t3amp_model,_ = vis_to_t3(image_to_vis(x, ft), data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3)
    return t3amp_model
 end
 
# function chi2_vis_dft_fg(x, g, dft, data ) # criterion function plus its gradient w/r x
#   cvis_model = image_to_vis_dft(x, dft);
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
#   cvis_model = image_to_vis(x, fftplan.fftplan_uv);
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
    return g2d/sum(g2d)
end

function disk(;npix::Int64=128, diameter::Union{Float64,Int64}=128, cent_x::Float64=64.5, cent_y::Float64=64.5, outside::Float64=0.0, centered::Bool = true)
        """
        Returns a 2D aperture of the desired diameter pixels, centered on (cent_x,cent_y) and on support npix X npix
        """
      if centered == true
        cent_x = (npix+1)/2
        cent_y = (npix+1)/2
      end
 
    x = (collect(1:npix) .- cent_x) / (diameter / 2.)
    y = (collect(1:npix) .- cent_y) / (diameter / 2.)
    xx = repeat(x,1,npix).^2
    rho = sqrt.(xx + xx')
    aperture = ones(size(rho))
    aperture[findall(rho.>1)] .= outside  # this is the aperture mask
   return aperture
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
    nx = size(x,1)
    flux= sum(x)
    c = cdg(x)
    f = (c[1]-(nx+1)/2)^2+(c[2]-(nx+1)/2)^2
    xx = [mod(i-1,nx)+1 for i=1:nx*nx]
    yy = [div(i-1,nx)+1 for i=1:nx*nx]
    g[:] = 2*(c[1]-(nx+1)/2)*xx + 2*(c[2]-(nx+1)/2)*yy
    if verb == true
        @printf(" COG: [%.2f, %.2f] REGC: %.2f", c[1], c[2], f);
    end
    return f
end

function visual_radial_params(angles, x, i)
    nx = size(x,1)
    xx = repeat( collect(1:nx) .- (nx-1)/2, 1, nx)
    yy = xx'
    ϕ = angles[1]/180*pi; #position angle
    inc = angles[2]/180*pi;
    rx = yy*cos(ϕ) + xx*sin(ϕ)
    ry = (-yy*sin(ϕ) + xx*cos(ϕ))*cos(inc)
    rr = vec(sqrt.(rx.^2 + ry.^2));
    indx = findall(rr.>i-1 .&& rr.<i+1);
    y=deepcopy(x); y[indx].=0;imdisp(y)
end

function setup_radial_reg(params, nx)
# Create radial operators
xx = repeat( collect(1:nx) .- (nx-1)/2, 1, nx); yy = xx';
ϕ = params[1]/180*pi; #position angle
inc = params[2]/180*pi;
rx = yy*cos(ϕ) + xx*sin(ϕ);
ry = (-yy*sin(ϕ) + xx*cos(ϕ))*cos(inc);
rr = vec(sqrt.(rx.^2 + ry.^2));
# Visual check
#indx= findall(rr.>i-1 .&& rr.<i+1) ; y=deepcopy(x); y[indx].=0;imdisp(y)
nrad = div(nx,2)-1
profile_mask = Array{Vector{Int64}}(undef, nrad)
profile_npix = zeros(Int64, nrad)
radH = Array{SparseMatrixCSC{Float64, Int64}}(undef, nrad)
radM = Array{SparseMatrixCSC{Float64, Int64}}(undef, nrad)
@Threads.threads for i=1:nrad
    profile_mask[i] = findall(rr.>i-1 .&& rr.<i+1);
    profile_npix[i] = length(profile_mask[i])
    P = sparse(1:profile_npix[i],profile_mask[i], ones(Float64,profile_npix[i]),profile_npix[i],nx*nx)
    O = ones(Float64,profile_npix[i])/profile_npix[i]
    radH[i] = (P .- O'*P)/sqrt(profile_npix[i]-1)
    radM[i] = O'*P
end
# Check: var(x[profile_mask[i]]) == norm(radH[i]*x)^2 to numerical precision
# Create big H so that sum([var(x[profile_mask[i]]) for i=1:nrad]) == norm(H*x)^2
H = vcat(radH...)
M = vcat(radM...)
# and big G
G = 2*reduce(+, radH[i]'*radH[i] for i=1:nrad)
return H, G, M
end

function radial_variance(x,g; H=[], G=[], verb = false)
# Takes in matrices to compute the total radial variance and its gradient
    nx = size(x,1)
    f = norm(H*vec(x))^2
    g[:] .=  reshape(G*vec(x),nx,nx);
    if verb == true
        @printf(" radialvar: %.3f", f);
    end
    return f
end

function tvsq(x,tvsq_g; verb = false)
    # Total squared variation
    nx = size(x,1)
    y = reshape(x, nx, nx);
    lx = circshift(y,(0,-1));
    ly = circshift(y,(-1,0));
    rx = circshift(y,(0,1));
    ry = circshift(y,(1,0));
    tvsq_f = norm(y-rx)^2+norm(y-ry)^2
    tvsq_g[:] = 2*vec(4*y-lx-ly-rx-ry)
    if verb == true
        @printf(" TVSQ: %.3f", tvsq_f);
    end
    return tvsq_f
end

function tv(x,tv_g; verb = false, ϵ=1e-8)
    # Total variation
    # TODO: - treat edges properly
    #       - check vs matrix implementation
    #       - check performance/precision vs FFT version
    nx = size(x,1)
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
        @printf(" TV: %.3f", tv_f);
    end
    return tv_f
end

function l1l2(x, g; verb = false, ϵ=1e-8, α = 1e-4)
    # Isotropic L1-L2 / Fair loss function /  Mugnier-Brette expression
    nx = size(x,1)
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
        @printf(" ℓ1ℓ2: %.3f", f);
    end
    return f
end

function l2sq(x,g; verb = false)
    f = sum(x.^2)
    g[:] =  2*x;
    if verb == true
        @printf(" ℓ2^2: %.3f", f);
    end
    return f
end


function l1hyp(x,g; verb = false,ϵ=1e-9)
    f = sum(sqrt.(x.^2 .+ϵ^2).-ϵ)
    g[:] =  x./sqrt.(x.^2 .+ϵ^2);
    if verb == true
        @printf(" ℓ1hyp: %.3f", f);
    end
    return f
end

function l1l2w(x,g; verb = false)
    f = sum(x-log.(1.0 .+ x))
    g[:] =  1.0 .- 1.0./(1.0 .+x);
    if verb == true
        @printf(" ℓ1ℓ2w: %.3f", f);
    end
    return f
end

function entropy(x,g; verb = false, ϵ=1e-12)
    f = sum(x.*log.(abs.(x).+ϵ) - x)
    g[:] =  log.(abs.(x).+ϵ);
    if verb == true
        @printf(" MAXENT: %.3f", f);
    end
    return f
end

function compactness(x,g; verb = false, w = 20.0) # w is the size in pixels of the soft-support 
    nx = size(x,1)
    yy = repeat(collect(1:nx).-0.5*(nx+1),1,nx).^2
    rr = (yy+yy')/(nx*nx)
    f = if isnothing(w)
        g[:] =  2*rr.*x
        sum(rr.*(x.^2))
    else
        soft_support = 1.0 ./ (1.0 .+ 2*rr/w^2)
        soft_support /= sum(soft_support)
        g[:] =  2*x./soft_support;
        sum((x.^2)./soft_support)
    end
    if verb == true
        @printf(" compactness: %.3f", f);
    end
    return f
end

function reg_support(x,g; prior=[], verb = false) # assumes prior is vec()
    mask = zeros(eltype(x), size(x))
    if prior !=[]
        mask = 1.0.-(prior.>0)
    end
    f = sum(mask.*(x.^2))
    g[:] =  2*mask.*x;
    if verb == true
        @printf(" support: %.3f", f);
    end
    return f
end

function trans_structnorm(x, g;verb=false, ϵ=1e-12)
    #this x is under the form (npix,nwavs)
    #but return the gradient as a 1D vector to use with Optimpack
    d = sqrt.(dropdims(sum(x.^2, dims=3), dims=3).+ϵ^2)
    f = sum(d.-ϵ);
    g[:] = x./d;
    return f
end

function trans_tv(x, g; ζ=1e-13)
    # Transpectral or transtemporal regularization
    #this x is under the form (npix,nwavs)
    #but return the gradient as a 1D vector to use with Optimpack
    xlplus1_minus_xl = circshift(x, (0,0,-1))-x
    xl_minus_xminus1 = x-circshift(x, (0,0,1))
    d1 = sqrt.(xlplus1_minus_xl.^2 .+ ζ^2)
    d2 = sqrt.(xl_minus_xminus1.^2 .+ ζ^2)
    f = sum(d1.-ζ)
    g[:] = -xlplus1_minus_xl./d1 + xl_minus_xminus1./d2
    return f
end


function trans_tvsq(x, g;verb=false)
    # Transpectral or transtemporal regularization
    #this x is under the form (npix,nwavs)
    #but return the gradient as a 1D vector to use with Optimpack
    xlplus1_minus_xl = circshift(x, (0,0,-1))-x
    xl_minus_xminus1 = x-circshift(x, (0,0,1))
    d1 = xlplus1_minus_xl.^2
    #d2 = xl_minus_xminus1.^2
    f = sum(d1)
    g[:] = 2*(-xlplus1_minus_xl + xl_minus_xminus1)
    return f
end

function trans_l1l2(x, g;verb=false, ζ=1e-13, δ=2.0)
    # Transpectral or transtemporal regularization using the L1L2 scheme
    #  x is under the form (npix,nwavs)
    #but return the gradient as a 1D vector to use with Optimpack
    xlplus1_minus_xl = circshift(x, (0,0,-1))-x
   # xlplus1_minus_xl[:,end] .= 0
    xl_minus_xminus1 = x-circshift(x, (0,0,1))
    # xl_minus_xminus1[:,1] .=0
    d1 = sqrt.(xlplus1_minus_xl.^2 .+ ζ^2)
    d2 = sqrt.(xl_minus_xminus1.^2 .+ ζ^2)
    f = sum(d1/δ - log.(1.0 .+d1/δ).-ζ)
    g[:] = vec((-xlplus1_minus_xl./d1 + xl_minus_xminus1./d2)/δ-((-xlplus1_minus_xl./d1)./(δ .+ d1) + (xl_minus_xminus1./d2)./(δ .+d2)))
    return f
end

function regularization(x, reg_g; printcolor = :black, regularizers=[], verb=true) # compound regularization
    reg_f = 0.0;
    if verb == true && !isempty(regularizers)
        print("\nReg:");
    end
    for ireg in regularizers
            temp_g = zeros(eltype(x), size(x))
            reg_f += @match ireg[1] begin 
                "centering"   => ireg[2]*reg_centering(x, temp_g; verb)
                "tv"          => ireg[2]*tv(x,temp_g; verb)
                "tvsq"        => ireg[2]*tvsq(x,temp_g; verb)
                "EPLL"        => ireg[2]*EPLL_fg(x,temp_g, ireg[3])
                "l1l2"        => ireg[2]*l1l2(x,temp_g; verb, α = ireg[3])
                "l1l2w"       => ireg[2]*l1l2w(x,temp_g; verb)
                "l1hyp"       => ireg[2]*l1hyp(x,temp_g; verb)
                "l2sq"        => ireg[2]*l2sq(x,temp_g; verb)
                "compactness" => ireg[2]*compactness(x,temp_g; verb, w = length(ireg) > 2 ? ireg[3] : nothing)
                "radialvar"   => ireg[2]*radial_variance(x,temp_g, H=ireg[3], G=ireg[4]; verb)
                "entropy"     => ireg[2]*entropy(x,temp_g; verb)
                "support"     => ireg[2]*reg_support(x, prior=ireg[3], temp_g; verb)
                _             => error("Unknown regularizer")
            end
            reg_g[:,:] += ireg[2]*temp_g
    end
    if (verb==true)
        print("\n");
    end
    return reg_f
end

# DFT version
function chi2_f(x::AbstractMatrix{<:AbstractFloat}, dft::AbstractMatrix{<:Complex}, data::OIdata; weights = [1.0,1.0,1.0],  cvis = [], printcolor =:black, verb=true, vonmises=false)
    flux = sum(x);
    cvis_model = image_to_vis(x, dft);
    if length(cvis)>0
        #Note: cvis_model includes all the complex visibilities needed to compute V2, T3, etc.
        #      while the cvis variable is used to export visibility observables (e.g. diff vis or diff phi)
        cvis[:] = cvis_model[data.indx_vis]
    end
    v2_model = vis_to_v2(cvis_model, data.indx_v2);
    t3_model, t3amp_model, t3phi_model = vis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
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
        printstyled(@sprintf("V2: %.2f ", chi2_v2/data.nv2), color=:red)
        printstyled(@sprintf("T3A: %.2f ", chi2_t3amp/data.nt3amp), color=:blue);
        printstyled(@sprintf("T3P: %.2f ", chi2_t3phi/data.nt3phi), color=:green);
        printstyled(@sprintf("Flux: %.4f ", flux), color=:black);
    end
    return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

# NFFT version
function chi2_f(x::AbstractMatrix{<:AbstractFloat}, ftplan::AbstractVector{<:NFFT.NFFTPlan}, data::OIdata; weights = [1.0,1.0,1.0], cvis = [], printcolor =:black,  verb = false, vonmises=false)
    flux = sum(x);
    cvis_model = image_to_vis(x, ftplan[1]);
    if length(cvis)>0
        #Note: cvis_model includes all the complex visibilities needed to compute V2, T3, etc.
        #      while the cvis variable is used to export visibility observables (e.g. diff vis or diff phi)
        cvis[:] = cvis_model[data.indx_vis]
    end
    v2_model = vis_to_v2(cvis_model, data.indx_v2);
    t3_model, t3amp_model, t3phi_model = vis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
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
        printstyled(@sprintf("V2: %.2f ", chi2_v2/data.nv2), color=:red)
        printstyled(@sprintf("T3A: %.2f ", chi2_t3amp/data.nt3amp), color=:blue);
        printstyled(@sprintf("T3P: %.2f ", chi2_t3phi/data.nt3phi), color=:green);
        printstyled(@sprintf("Flux: %.4f ", flux), color=:black);
    end
    return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

# DFT version
function chi2_fg(x::AbstractMatrix{<:AbstractFloat}, g::AbstractMatrix{<:AbstractFloat}, dft::AbstractMatrix{<:Complex}, data::OIdata; weights = [1.0,1.0,1.0],  cvis = [],  printcolor =:black, verb=true, vonmises=false)
    flux = sum(x);
    cvis_model = image_to_vis(x, dft);
    if length(cvis)>0
        #Note: cvis_model includes all the complex visibilities needed to compute V2, T3, etc.
        #      while the cvis variable is used to export visibility observables (e.g. diff vis or diff phi)
        cvis[:] = cvis_model[data.indx_vis]
    end
    v2_model = vis_to_v2(cvis_model, data.indx_v2);
    t3_model, t3amp_model, t3phi_model = vis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
    chi2_v2 = 0.0;
    chi2_t3amp = 0.0;
    chi2_t3phi = 0.0;
    g_v2 = 0.0
    g_t3amp = 0.0
    g_t3phi = 0.0
    if (weights[1]>0)&&(data.nv2>0)
        chi2_v2 = norm((v2_model - data.v2)./data.v2_err)^2;
        g_v2 = real(transpose(dft[data.indx_v2,:])*(4*((v2_model-data.v2)./data.v2_err.^2).*conj(cvis_model[data.indx_v2])));
    end

    if (weights[2]>0)&&(data.nt3amp>0)
        chi2_t3amp = norm((t3amp_model - data.t3amp)./data.t3amp_err)^2;
        dT3 = 2.0*(t3amp_model-data.t3amp)./(data.t3amp_err.^2)
        g_t3amp = real(transpose(dft[data.indx_t3_1,:])*(dT3.*conj(cvis_model[data.indx_t3_1])./abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3])))+real(transpose(dft[data.indx_t3_2,:])*(dT3.*conj(cvis_model[data.indx_t3_2])./abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3]) ))+real(transpose(dft[data.indx_t3_3,:])*(dT3.*conj(cvis_model[data.indx_t3_3])./abs.(cvis_model[data.indx_t3_3]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]) ))
    end

    if (weights[3]>0)&&(data.nt3phi>0)
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
    g[:] = weights[1]*g_v2 .+ weights[2]*g_t3amp .+ weights[3]*g_t3phi
    # g[:] = (g - sum(vec(x).*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
    if verb==true
        printstyled(@sprintf("V2: %.2f ", chi2_v2/data.nv2), color=:red)
        printstyled(@sprintf("T3A: %.2f ", chi2_t3amp/data.nt3amp), color=:blue);
        printstyled(@sprintf("T3P: %.2f ", chi2_t3phi/data.nt3phi), color=:green);
        printstyled(@sprintf("Flux: %.4f ", flux), color=:black);
    end
    return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

#NFFT version
function chi2_fg(x::AbstractMatrix{<:AbstractFloat}, g::AbstractMatrix{<:AbstractFloat}, ftplan::AbstractVector{<:NFFT.NFFTPlan}, data::OIdata; weights = [1.0,1.0,1.0], cvis = [], printcolor =:black,  verb = false, vonmises=false)
    flux = sum(x);
    cvis_model = image_to_vis(x, ftplan[1]);
    if length(cvis)>0
        #Note: cvis_model includes all the complex visibilities needed to compute V2, T3, etc.
        #      while the cvis variable is used to export visibility observables (e.g. diff vis or diff phi)
        cvis[:] = cvis_model[data.indx_vis]
    end
    v2_model = vis_to_v2(cvis_model, data.indx_v2);
    t3_model, t3amp_model, t3phi_model = vis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
    chi2_v2 = 0.0;
    chi2_t3amp = 0.0;
    chi2_t3phi = 0.0;
    g_v2 = 0.0
    g_t3amp = 0.0
    g_t3phi = 0.0
    if (weights[1]>0)&&(data.nv2>0)
        chi2_v2 = norm((v2_model - data.v2)./data.v2_err)^2;
        g_v2 = real(adjoint(ftplan[3])*(4*((v2_model-data.v2)./data.v2_err.^2).*cvis_model[data.indx_v2]));
    end

    if (weights[2]>0)&&(data.nt3amp>0)
        chi2_t3amp = norm((t3amp_model - data.t3amp)./data.t3amp_err)^2;
        dT3 = 2.0*(t3amp_model-data.t3amp)./(data.t3amp_err.^2)
        g_t3amp = real(adjoint(ftplan[4])*(dT3.*cvis_model[data.indx_t3_1]./abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3]))) + real(adjoint(ftplan[5])*(dT3.*cvis_model[data.indx_t3_2]./abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3]))) + real(adjoint(ftplan[6])*(dT3.*cvis_model[data.indx_t3_3]./abs.(cvis_model[data.indx_t3_3]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2])))
    end

    if (weights[3]>0)&&(data.nt3phi>0)
        if vonmises == false
            chi2_t3phi = norm(mod360(t3phi_model - data.t3phi)./data.t3phi_err)^2;
            dT3 = -360.0/pi*mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2
            g_t3phi = imag(adjoint(ftplan[4])*(dT3./abs2.(cvis_model[data.indx_t3_1]).*cvis_model[data.indx_t3_1])+adjoint(ftplan[5])*(dT3./abs2.(cvis_model[data.indx_t3_2]).*cvis_model[data.indx_t3_2])+ adjoint(ftplan[6])*(dT3./abs2.(cvis_model[data.indx_t3_3]).*cvis_model[data.indx_t3_3]))
        else
            chi2_t3phi =  sum(-2*data.t3phi_vonmises_err.*cos.((t3phi_model - data.t3phi)/180*pi) + data.t3phi_vonmises_chi2_offset)
            dT3 = -2.0*data.t3phi_vonmises_err.*sin.((t3phi_model - data.t3phi)/180*pi)
            g_t3phi = imag(adjoint(ftplan[4])*(dT3./abs2.(cvis_model[data.indx_t3_1]).*cvis_model[data.indx_t3_1]) + adjoint(ftplan[5])*(dT3./abs2.(cvis_model[data.indx_t3_2]).*cvis_model[data.indx_t3_2]) + adjoint(ftplan[6])*(dT3./abs2.(cvis_model[data.indx_t3_3]).*cvis_model[data.indx_t3_3]))
        end
    end

    g[:] = weights[1]*g_v2 .+ weights[2]*g_t3amp .+ weights[3]*g_t3phi
    if verb==true
        if (weights[1]>0)&&(data.nv2>0)
            printstyled(@sprintf("V2: %.2f ", chi2_v2/data.nv2), color=:red)
        else
            printstyled("V2: (N/A) ", color=:black) # disabled
        end
        printstyled(@sprintf("T3A: %.2f ", chi2_t3amp/data.nt3amp), color=:blue);
        printstyled(@sprintf("T3P: %.2f ", chi2_t3phi/data.nt3phi), color=:green);
        printstyled(@sprintf("Flux: %.4f ", flux), color=:black);
    end
    return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

# ---------------------------------------------------------------------------
# _polychromatic_vis_flux_chi2: shared helper for computing VIS and FLUX chi2
# contributions in polychromatic mode (used by both chi2_polychromatic_f and
# crit_polychromatic_fg).  Returns (chi2, ndof_extra) and the model quantities
# needed for gradient computation.
# ---------------------------------------------------------------------------
function _polychromatic_vis_flux_chi2(x, data, cvis, nwavs;
        use_diffphases, use_diffvisamp, use_abs_visamp, verb)
    chi2 = 0.0; ndof = 0
    diffphi_model = nothing; diffphi = nothing; diffphi_err = nothing
    cvis_ref_all = nothing
    # Differential observables (ASPRO2 convention)
    if use_diffphases || use_diffvisamp
        cvis_all = hcat(cvis...)
        cvis_sum = vec(sum(cvis_all, dims=2))
        if use_diffphases
            diffphi_model = zeros(size(cvis_all))
            cvis_ref_all  = zeros(ComplexF64, size(cvis_all))
            for i=1:nwavs
                cvis_ref_all[:,i] = (cvis_sum .- cvis_all[:,i]) ./ (nwavs - 1)
                diffphi_model[:,i] = angle.(cvis_all[:,i] ./ cvis_ref_all[:,i]) .* (180.0/pi)
            end
            diffphi     = hcat([data[i].visphi for i=1:nwavs]...)
            diffphi_err = hcat([data[i].visphi_err for i=1:nwavs]...)
            chi2_dp = norm(mod360(diffphi_model-diffphi)./diffphi_err)^2
            chi2 += chi2_dp;  ndof += sum(data[i].nvisphi for i in 1:nwavs)
            verb && printstyled(@sprintf("DiffPhi: %.2f ", chi2_dp/length(diffphi_model)), color=:magenta)
        end
        if use_diffvisamp
            diffvisamp_model = zeros(size(cvis_all))
            for i=1:nwavs
                cvis_ref = (cvis_sum .- cvis_all[:,i]) ./ (nwavs - 1)
                diffvisamp_model[:,i] = abs.(cvis_all[:,i] ./ cvis_ref)
            end
            mean_dva = vec(mean(diffvisamp_model, dims=2))
            diffvisamp_model ./= mean_dva
            diffvisamp     = hcat([data[i].visamp for i=1:nwavs]...)
            diffvisamp_err = hcat([data[i].visamp_err for i=1:nwavs]...)
            chi2_dva = norm((diffvisamp_model - diffvisamp)./diffvisamp_err)^2
            chi2 += chi2_dva;  ndof += sum(data[i].nvisamp for i in 1:nwavs)
            verb && printstyled(@sprintf("DiffVA: %.2f ", chi2_dva/length(diffvisamp_model)), color=:magenta)
        end
    end
    # Absolute VISAMP
    if use_abs_visamp && !use_diffvisamp
        abs_va_model = hcat([abs.(cvis[i]) for i=1:nwavs]...)
        abs_va       = hcat([data[i].visamp for i=1:nwavs]...)
        abs_va_err   = hcat([data[i].visamp_err for i=1:nwavs]...)
        chi2_ava = norm((abs_va_model - abs_va)./abs_va_err)^2
        chi2 += chi2_ava;  ndof += sum(data[i].nvisamp for i in 1:nwavs)
        verb && printstyled(@sprintf("VA: %.2f ", chi2_ava/length(abs_va_model)), color=:magenta)
    end
    # OI_FLUX
    C_flux = NaN; flux_residual = Float64[]; flux_chan_idx = Int[]
    use_flux = any(data[i].nflux > 0 for i in 1:nwavs)
    if use_flux
        fm = Float64[]; fd = Float64[]; fe = Float64[]
        for i in 1:nwavs
            d = data[i]; d.nflux > 0 || continue
            fi = sum(x[:,:,i])
            append!(fm, fill(fi, d.nflux)); append!(fd, d.flux)
            append!(fe, d.flux_err);        append!(flux_chan_idx, fill(i, d.nflux))
        end
        w = 1.0 ./ fe.^2
        C_flux = sum(fm .* fd .* w) / sum(fm.^2 .* w)
        flux_residual = (C_flux .* fm .- fd) ./ fe.^2
        chi2_fl = sum(((C_flux .* fm .- fd) ./ fe).^2)
        chi2 += chi2_fl;  ndof += length(fd)
        verb && printstyled(@sprintf("OI_FLUX: %.2f (C=%.4f) ", chi2_fl/length(fd), C_flux), color=:yellow)
    end
    return (; chi2, ndof, diffphi_model, diffphi, diffphi_err, cvis_ref_all,
              C_flux, flux_residual, flux_chan_idx, use_flux)
end

# ===========================================================================
# chi2_polychromatic_f — NFFT version
# ===========================================================================
function chi2_polychromatic_f(x::AbstractArray{<:AbstractFloat,3}, ft::AbstractVector{<:AbstractVector{<:NFFTPlan}}, data::AbstractVector{<:OIdata};weights = [1.0,1.0,1.0], printcolor= [], use_diffphases = false, verb = false)
    nwavs = length(ft);
    npix = size(x,1);
    if printcolor == []
        printcolor=[ :black for i=1:nwavs]
    end

    # Determine VIS observable types from data metadata
    amptyp = !isempty(data) ? data[1].amptyp : ""
    phityp = !isempty(data) ? data[1].phityp : ""
    # Auto-enable differential phases if PHITYP says so
    use_diffphases = use_diffphases || phityp == "differential"
    use_diffvisamp = amptyp == "differential"
    use_abs_visamp = amptyp == "absolute"

    cvis = fill((ComplexF64[]), nwavs);
    if use_diffphases || use_diffvisamp || use_abs_visamp
        for i=1:nwavs
            cvis[i] = Vector{ComplexF64}(undef, length(data[i].indx_vis))
        end
    end
    f = zeros(nwavs)
    for i=1:nwavs # weighted sum -- should probably do the computation in parallel
        if verb == true
            printstyled("Channel $i ",color=printcolor[i]);
        end

        if use_diffphases || use_diffvisamp || use_abs_visamp
            f[i] = chi2_f(x[:,:,i], ft[i], data[i], cvis=cvis[i], verb = verb, weights = weights);
        else
            f[i] = chi2_f(x[:,:,i], ft[i], data[i], verb = verb, weights = weights);
        end
        fr = f[i]/(weights[1]*data[i].nv2+weights[2]*data[i].nt3amp+weights[3]*data[i].nt3phi)
        if verb == true
            printstyled(@sprintf("Chi2r: %.2f Chi2: %.2f\n", fr, f[i]),color=printcolor[i]);
        end
    end
    chi2f = sum(f)
    ndof = Int(sum([weights[1]*data[i].nv2+weights[2]*data[i].nt3amp+weights[3]*data[i].nt3phi for i=1:nwavs]));
    if verb == true
    printstyled(@sprintf("All V2+T3A+T3P -- Chi2: %.2f Chi2r: %.2f\n", chi2f, chi2f/ndof), color=:red);
    end

    # VIS + FLUX chi2 via shared helper
    vis_flux = _polychromatic_vis_flux_chi2(x, data, cvis, nwavs;
        use_diffphases=use_diffphases, use_diffvisamp=use_diffvisamp,
        use_abs_visamp=use_abs_visamp, verb=verb)
    chi2f += vis_flux.chi2
    ndof  += vis_flux.ndof

    return chi2f/ndof;
end

# ===========================================================================
# chi2_polychromatic_f — DFT version
# ===========================================================================
function chi2_polychromatic_f(x::AbstractArray{<:AbstractFloat,3}, ft::AbstractVector{<:AbstractMatrix{<:Complex}}, data::AbstractVector{<:OIdata};weights = [1.0,1.0,1.0], printcolor= [], use_diffphases = false, verb = false)
    nwavs = length(ft);
    npix = size(x,1);
    if printcolor == []
        printcolor=[ :black for i=1:nwavs]
    end

    # Determine VIS observable types from data metadata
    amptyp = !isempty(data) ? data[1].amptyp : ""
    phityp = !isempty(data) ? data[1].phityp : ""
    use_diffphases = use_diffphases || phityp == "differential"
    use_diffvisamp = amptyp == "differential"
    use_abs_visamp = amptyp == "absolute"

    cvis = fill((ComplexF64[]), nwavs);
    if use_diffphases || use_diffvisamp || use_abs_visamp
        for i=1:nwavs
            cvis[i] = Vector{ComplexF64}(undef, length(data[i].indx_vis))
        end
    end
    f = zeros(nwavs)
    for i=1:nwavs
        if verb == true
            printstyled("Channel $i ",color=printcolor[i]);
        end

        if use_diffphases || use_diffvisamp || use_abs_visamp
            f[i] = chi2_f(x[:,:,i], ft[i], data[i], cvis=cvis[i], verb = verb, weights = weights);
        else
            f[i] = chi2_f(x[:,:,i], ft[i], data[i], verb = verb, weights = weights);
        end
        fr = f[i]/(weights[1]*data[i].nv2+weights[2]*data[i].nt3amp+weights[3]*data[i].nt3phi)
        if verb == true
            printstyled(@sprintf("Chi2r: %.2f Chi2: %.2f\n", fr, f[i]),color=printcolor[i]);
        end
    end
    chi2f = sum(f)
    ndof = Int(sum([weights[1]*data[i].nv2+weights[2]*data[i].nt3amp+weights[3]*data[i].nt3phi for i=1:nwavs]));
    if verb == true
    printstyled(@sprintf("All V2+T3A+T3P -- Chi2: %.2f Chi2r: %.2f\n", chi2f, chi2f/ndof), color=:red);
    end

    # VIS + FLUX chi2 via shared helper
    vis_flux = _polychromatic_vis_flux_chi2(x, data, cvis, nwavs;
        use_diffphases=use_diffphases, use_diffvisamp=use_diffvisamp,
        use_abs_visamp=use_abs_visamp, verb=verb)
    chi2f += vis_flux.chi2
    ndof  += vis_flux.ndof

    return chi2f/ndof;
end

function crit_fg(x,g::AbstractMatrix{<:AbstractFloat}, ft, data::OIdata; weights = [1.0,1.0,1.0], cvis = [], printcolor = :black, regularizers=[], verb = true)
    chi2 = chi2_fg(x, g, ft, data, cvis = cvis, verb = verb, weights = weights);
    reg = regularization(x, g, regularizers=regularizers, printcolor = printcolor, verb = verb);
    flux = sum(x)
    g[:] = (g .- sum(x.*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
    return chi2 + reg;
end

function crit_f(x::AbstractMatrix{<:AbstractFloat}, fftplan::AbstractVector{<:NFFT.NFFTPlan}, data::OIdata; weights = [1.0,1.0,1.0], cvis = [], printcolor = :black, regularizers=[], verb = true)
    chi2 = chi2_f(x, fftplan, data, cvis = cvis, verb = verb, weights = weights );
    g = zeros(eltype(x), size(x));
    reg = regularization(x, g,  regularizers=regularizers, printcolor = printcolor, verb = verb);
    return chi2 + reg;
end

function crit_multitemporal_fg(x::AbstractArray{<:AbstractFloat,3}, g::AbstractArray{<:AbstractFloat,3}, ft::AbstractVector{<:AbstractVector{<:NFFT.NFFTPlan}}, data::Array{OIdata,1};weights = [1.0,1.0,1.0], printcolor= [], epochs_weights=[],regularizers=[], verb = false)
    nepochs = length(ft);
    if epochs_weights == []
        epochs_weights=ones(eltype(x_start), nepochs);
    end
    if printcolor == []
        printcolor=Array{Symbol}(undef,nepochs);
        printcolor[:] .= :black
    end
    npix = div(length(x),nepochs);
    f = 0.0;
    for i=1:nepochs # weighted sum -- should probably do the computation in parallel
        subg = Array{eltype(x)}(undef, npix, npix);
        printstyled("Epoch $i ",color=printcolor[i]);
        f += epochs_weights[i]*crit_fg(x[:,:,i], subg, ft[i], data[i], regularizers=regularizers[i], printcolor = printcolor[i], verb = verb, weights = weights);
        g[tslice] = epochs_weights[i]*subg
    end

    # cross temporal regularization
    if length(regularizers)>nepochs
        if (regularizers[nepochs+1][1][1] == "temporal_tvsq")  & (nepochs>1)
            y = reshape(x,(npix,nepochs))
            temporalf = sum( (y[:,2:end]-y[:,1:end-1]).^2 )
            tv_g = Array{eltype(x)}(undef, npix, nepochs)
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

# ===========================================================================
# _polychromatic_vis_gradient! — shared VIS chi2 + gradient logic.
# Computes chi2 from differential phase, absolute VISAMP, and flux, and
# accumulates their gradients into g.  The `vis_adjoint_fn` closure handles
# the DFT/NFFT adjoint convention difference:
#   NFFT:  vis_adjoint_fn(rhs, i) = vec(adjoint(ft[i][2]) * rhs)
#   DFT:   vis_adjoint_fn(rhs, i) = vec(transpose(ft[i][data[i].indx_vis,:]) * conj(rhs))
# ===========================================================================
function _polychromatic_vis_gradient!(x, g, data, cvis, nwavs, npix, vis_adjoint_fn;
        use_diffphases, use_diffvisamp, use_abs_visamp, verb)
    f = 0.0

    if use_diffphases || use_diffvisamp
        cvis_all = hcat(cvis...)
        cvis_sum = vec(sum(cvis_all, dims=2))

        if use_diffphases
            diffphi_model = zeros(size(cvis_all))
            cvis_ref_all  = zeros(ComplexF64, size(cvis_all))
            for i=1:nwavs
                cvis_ref_all[:,i] = (cvis_sum .- cvis_all[:,i]) ./ (nwavs - 1)
                diffphi_model[:,i] = angle.(cvis_all[:,i] ./ cvis_ref_all[:,i]) .* (180.0/pi)
            end
            diffphi     = hcat([data[i].visphi for i=1:nwavs]...)
            diffphi_err = hcat([data[i].visphi_err for i=1:nwavs]...)
            chi2_diffphi = norm(mod360(diffphi_model-diffphi)./diffphi_err)^2
            f += chi2_diffphi
            verb && printstyled(@sprintf("DiffPhi: %.2f ", chi2_diffphi/length(diffphi_model)), color=:magenta)

            # Gradient: d(arg(V_i/V_ref_i))/d(pixel) via chain rule
            # Leading-order approximation: gradient through numerator V_i only
            # d(arg(V/Vref))/dpixel ≈ d(arg(V))/dpixel = Im((1/V) * dV/dpixel)
            # Following the T3phi pattern: rhs = scalar_weight .* V ./ |V|²
            for i=1:nwavs
                dT3 = -360.0/pi .* mod360(diffphi_model[:,i]-diffphi[:,i])./diffphi_err[:,i].^2
                rhs_vis = dT3 .* cvis_all[:,i] ./ abs2.(cvis_all[:,i])
                g_contrib = imag.(vis_adjoint_fn(rhs_vis, i))
                g[:,:,i] .+= reshape(g_contrib, npix, npix)
            end
        end
    end

    # Absolute VISAMP chi2 + gradient
    if use_abs_visamp && !use_diffvisamp
        for i=1:nwavs
            abs_visamp_model = abs.(cvis[i])
            abs_visamp       = data[i].visamp
            abs_visamp_err   = data[i].visamp_err
            residual = (abs_visamp_model .- abs_visamp) ./ abs_visamp_err.^2
            f += norm((abs_visamp_model .- abs_visamp) ./ abs_visamp_err)^2
            # d|V|/dpixel: Re(V/|V| * dV*/dpixel) for NFFT, Re(V*/|V| * dV/dpixel) for DFT
            # With our vis_adjoint convention, pass cvis (NFFT: adjoint conjugates; DFT: explicit conj in wrapper)
            rhs_vis = 2.0 .* residual .* cvis[i] ./ abs_visamp_model
            g[:,:,i] .+= reshape(real.(vis_adjoint_fn(rhs_vis, i)), npix, npix)
        end
    end

    # OI_FLUX chi2 + gradient
    use_flux = any(data[i].nflux > 0 for i in 1:nwavs)
    if use_flux
        flux_model_vec = Float64[]; flux_data_vec = Float64[]
        flux_err_vec = Float64[]; flux_chan_idx = Int[]
        for i in 1:nwavs
            d = data[i]; d.nflux > 0 || continue
            fi = sum(x[:,:,i])
            append!(flux_model_vec, fill(fi, d.nflux)); append!(flux_data_vec, d.flux)
            append!(flux_err_vec, d.flux_err); append!(flux_chan_idx, fill(i, d.nflux))
        end
        w = 1.0 ./ flux_err_vec.^2
        C_flux = sum(flux_model_vec .* flux_data_vec .* w) / sum(flux_model_vec.^2 .* w)
        residual_flux = (C_flux .* flux_model_vec .- flux_data_vec) ./ flux_err_vec.^2
        chi2_flux = sum(((C_flux .* flux_model_vec .- flux_data_vec) ./ flux_err_vec).^2)
        f += chi2_flux
        verb && printstyled(@sprintf("OI_FLUX: %.2f (C=%.4f) ", chi2_flux/length(flux_data_vec), C_flux), color=:yellow)
        for i in 1:nwavs
            idx = findall(flux_chan_idx .== i)
            isempty(idx) && continue
            g[:,:,i] .+= 2.0 * C_flux * sum(residual_flux[idx])
        end
    end
    return f
end

# ===========================================================================
# _polychromatic_transspectral_reg! — shared transspectral regularization
# ===========================================================================
function _polychromatic_transspectral_reg!(x, g, f, ndof, npix, nwavs, regularizers; verb=false)
    if length(regularizers) > nwavs
        ntransreg = length(regularizers[nwavs+1])
        tg = zeros(npix, npix, nwavs)
        for i=1:ntransreg
            if (regularizers[nwavs+1][i][1] == "transspectral_tv")
                tf = trans_tv(x,tg)
                f += regularizers[nwavs+1][i][2]*tf
                g[:,:,:] += regularizers[nwavs+1][i][2]*tg
                printstyled("Trans-spectral TV: $(regularizers[nwavs+1][i][2]*tf)\n", color=:yellow)
            end
            if (regularizers[nwavs+1][i][1] == "transspectral_tvsq")
                tf = trans_tvsq(x, tg)
                f += regularizers[nwavs+1][i][2]*tf
                g[:,:,:] += regularizers[nwavs+1][i][2]*tg
                printstyled("Trans-spectral squared TV: $(regularizers[nwavs+1][i][2]*tf)\n", color=:yellow)
            end
            if (regularizers[nwavs+1][i][1] == "transspectral_structnorm")
                tf = trans_structnorm(x,tg)
                f += regularizers[nwavs+1][i][2]*tf
                g[:,:,:] += regularizers[nwavs+1][i][2]*tg
                printstyled("Trans-spectral Structured Norm: $(regularizers[nwavs+1][i][2]*tf)\n", color=:yellow)
            end
            if (regularizers[nwavs+1][i][1] == "transspectral_l1l2")
                tf = trans_l1l2(x, tg, δ=regularizers[nwavs+1][i][3])
                f += regularizers[nwavs+1][i][2]*tf
                g[:,:,:] += regularizers[nwavs+1][i][2]*tg
                printstyled("Trans-spectral l1l2 norm: $(regularizers[nwavs+1][i][2]*tf)\n", color=:yellow)
            end
        end
        printstyled(@sprintf("Post trans -- Crit: %.2f Crit/dof: %.2f\n", f, f/ndof), color=:blue)
    end
    return f
end

# ===========================================================================
# crit_polychromatic_fg — NFFT version
# ===========================================================================
function crit_polychromatic_fg(x::AbstractArray{<:AbstractFloat,3}, g::AbstractArray{<:AbstractFloat,3}, ft::AbstractVector{<:AbstractVector{<:NFFTPlan}}, data::AbstractVector{<:OIdata};weights = [1.0,1.0,1.0], printcolor= [], regularizers=[], use_diffphases = false, verb = false)
    nwavs = length(ft); npix = size(x,1)
    printcolor == [] && (printcolor = [:black for i=1:nwavs])
    regularizers == [] && (regularizers = fill([], nwavs))

    amptyp = !isempty(data) ? data[1].amptyp : ""
    phityp = !isempty(data) ? data[1].phityp : ""
    use_diffphases = use_diffphases || phityp == "differential"
    use_diffvisamp = amptyp == "differential"
    use_abs_visamp = amptyp == "absolute"

    cvis = fill((ComplexF64[]), nwavs)
    if use_diffphases || use_diffvisamp || use_abs_visamp
        for i=1:nwavs; cvis[i] = Vector{ComplexF64}(undef, length(data[i].indx_vis)); end
    end

    f = 0.0
    for i=1:nwavs
        subg = zeros(eltype(x), npix, npix)
        verb && printstyled("Channel $i ",color=printcolor[i])
        f += crit_fg(x[:,:,i], subg, ft[i], data[i], regularizers=regularizers[i], cvis=cvis[i], printcolor=printcolor[i], verb=verb, weights=weights)
        g[:,:,i] = subg
    end
    ndof = Int(sum([weights[1]*data[i].nv2+weights[2]*data[i].nt3amp+weights[3]*data[i].nt3phi for i=1:nwavs]))

    # NFFT adjoint convention: adjoint(plan) already conjugates
    vis_adjoint_nfft = (rhs, i) -> vec(adjoint(ft[i][2]) * rhs)
    f += _polychromatic_vis_gradient!(x, g, data, cvis, nwavs, npix, vis_adjoint_nfft;
        use_diffphases=use_diffphases, use_diffvisamp=use_diffvisamp,
        use_abs_visamp=use_abs_visamp, verb=verb)

    f = _polychromatic_transspectral_reg!(x, g, f, ndof, npix, nwavs, regularizers; verb=verb)
    g[:] = g[:]/ndof
    return f/ndof
end

# ===========================================================================
# crit_polychromatic_fg — DFT version
# ===========================================================================
function crit_polychromatic_fg(x::AbstractArray{<:AbstractFloat,3}, g::AbstractArray{<:AbstractFloat,3}, ft::AbstractVector{<:AbstractMatrix{<:Complex}}, data::AbstractVector{<:OIdata};weights = [1.0,1.0,1.0], printcolor= [], regularizers=[], use_diffphases = false, verb = false)
    nwavs = length(ft); npix = size(x,1)
    printcolor == [] && (printcolor = [:black for i=1:nwavs])
    regularizers == [] && (regularizers = fill([], nwavs))

    amptyp = !isempty(data) ? data[1].amptyp : ""
    phityp = !isempty(data) ? data[1].phityp : ""
    use_diffphases = use_diffphases || phityp == "differential"
    use_diffvisamp = amptyp == "differential"
    use_abs_visamp = amptyp == "absolute"

    cvis = fill((ComplexF64[]), nwavs)
    if use_diffphases || use_diffvisamp || use_abs_visamp
        for i=1:nwavs; cvis[i] = Vector{ComplexF64}(undef, length(data[i].indx_vis)); end
    end

    f = 0.0
    for i=1:nwavs
        subg = zeros(eltype(x), npix, npix)
        verb && printstyled("Channel $i ",color=printcolor[i])
        f += crit_fg(x[:,:,i], subg, ft[i], data[i], regularizers=regularizers[i], cvis=cvis[i], printcolor=printcolor[i], verb=verb, weights=weights)
        g[:,:,i] = subg
    end
    ndof = Int(sum([weights[1]*data[i].nv2+weights[2]*data[i].nt3amp+weights[3]*data[i].nt3phi for i=1:nwavs]))

    # DFT adjoint convention: conj(transpose(F) * conj(rhs)) = Σ_j conj(F_jp) * rhs_j
    # This matches the NFFT adjoint operation exactly, so both Re() and Im() give identical results
    vis_adjoint_dft = (rhs, i) -> vec(conj.(transpose(ft[i][data[i].indx_vis, :]) * conj.(rhs)))
    f += _polychromatic_vis_gradient!(x, g, data, cvis, nwavs, npix, vis_adjoint_dft;
        use_diffphases=use_diffphases, use_diffvisamp=use_diffvisamp,
        use_abs_visamp=use_abs_visamp, verb=verb)

    f = _polychromatic_transspectral_reg!(x, g, f, ndof, npix, nwavs, regularizers; verb=verb)
    g[:] = g[:]/ndof
    return f/ndof
end

function image_to_vis(x::AbstractArray{<:AbstractFloat,3}, ft::Union{AbstractVector{<:AbstractVector{<:NFFTPlan}},AbstractVector{<:AbstractMatrix{<:Complex}}})
    nwavs = length(ft);
    npix = size(x,1);
    cvis = fill((ComplexF64[]), nwavs);
    Threads.@threads for i=1:nwavs
        cvis[i] = image_to_vis(x[:,:,i], ft[i]);
    end
    return cvis;
end

using OptimPackNextGen
function reconstruct(x_start::AbstractMatrix{<:AbstractFloat}, data::OIdata, ft; weights = [1.0,1.0,1.0], printcolor = :black, verb = false, maxiter = 100, regularizers =[], ftol= (0,1e-8), xtol=(0,1e-8), gtol=(0,1e-8))
    crit = (x,g)->crit_fg(x, g, ft, data, regularizers=regularizers, verb = verb , weights = weights)
    x_sol = OptimPackNextGen.vmlmb(crit, x_start, verb=verb, lower=0, maxiter=maxiter, blmvm=false, xtol = xtol, ftol = ftol, gtol=gtol);
    return x_sol
end

function reconstruct_multitemporal(x_start::AbstractArray{<:AbstractFloat,3}, data::AbstractVector{<:OIdata}, ft; weights = [1.0,1.0,1.0], epochs_weights =[], printcolor= [], verb = true, maxiter = 100, regularizers =[], ftol= (0,1e-8), xtol=(0,1e-8), gtol=(0,1e-8))
    x_sol = []
    if eltype(eltype(ft)) <: NFFTPlan
        crit = (x,g)->crit_multitemporal_fg(x, g, ft, data, printcolor=printcolor, weights = weights, epochs_weights=epochs_weights, regularizers=regularizers, verb = verb)
        x_sol = OptimPackNextGen.vmlmb(crit, x_start, verb=verb, lower=0, maxiter=maxiter, blmvm=false, xtol = xtol, ftol = ftol, gtol=gtol);
    else
        error("Sorry, polytemporal DFT methods not implemented yet");
    end
    return x_sol
end

function reconstruct_polychromatic(x_start::AbstractArray{<:AbstractFloat,3}, data::AbstractVector{<:OIdata}, ft; weights = [1.0,1.0,1.0], printcolor= [], verb = true, use_diffphases = false, maxiter = 100, regularizers =[], ftol= (0,1e-8), xtol=(0,1e-8), gtol=(0,1e-8))
    x_sol = []
    if regularizers == []
        regularizers = fill([],length(data))
    end
    crit = (x,g)->crit_polychromatic_fg(x, g, ft, data, weights = weights, printcolor=printcolor, regularizers=regularizers, use_diffphases = use_diffphases, verb = verb)
    x_sol = OptimPackNextGen.vmlmb(crit, x_start, verb=true, lower=0, maxiter=maxiter, blmvm=false, xtol = xtol, ftol = ftol, gtol=gtol);
    return x_sol
end


# ---------------------------------------------------------------------------
# SPARCO grey reconstruction
# ---------------------------------------------------------------------------
# The SPARCO model decomposes the total visibility as:
#
#        f_star * α * V_star + f_env * β * V_env
# V_tot = ----------------------------------------
#              f_star * α + f_bg * α + f_env * β
#
# where α = (λ/λ₀)^-4   (Rayleigh–Jeans stellar+background scaling)
#       β = (λ/λ₀)^(d-4) (environment power law, d = spectral index)
#       f_env = 1 - f_star - f_bg
#
# params[1] = f_star  (stellar flux fraction at λ₀)
# params[2] = f_bg    (incoherent background flux fraction at λ₀)
# params[3] = D       (stellar angular diameter, uniform disk)
# params[4] = d       (environment spectral index + 4)
# params[5] = λ₀      (reference wavelength, fixed)
# ---------------------------------------------------------------------------

# Shared helper: compute SPARCO composite visibility and intermediate quantities.
function _sparco_model(params, x_img, ftplan, data)
    λ0 = params[5]
    λ  = data.uv_lam
    α  = (λ/λ0).^-4.0
    β  = (λ/λ0).^(params[4]-4.0)
    fluxstar = params[1]*α
    fluxbg   = params[2]*α
    fluxenv  = (1.0-params[1]-params[2])*β
    Vstar = visibility_ud([params[3]], data.uv)
    Venv  = image_to_vis(x_img, ftplan[1])
    u = fluxstar.*Vstar .+ fluxenv.*Venv
    v = fluxstar .+ fluxenv .+ fluxbg
    cvis_model = u ./ v
    imratio = fluxenv ./ v
    return (; cvis_model, Vstar, Venv, fluxstar, fluxbg, fluxenv, α, β, u, v, imratio)
end

function chi2_sparco_f(x::AbstractMatrix{<:AbstractFloat}, params::AbstractVector{<:AbstractFloat},
        ftplan::AbstractVector{<:NFFT.NFFTPlan}, data::OIdata;
        verb=true, weights=[1.0,1.0,1.0], vonmises=false)
    m = _sparco_model(params, x, ftplan, data)
    v2_model = vis_to_v2(m.cvis_model, data.indx_v2)
    t3_model, t3amp_model, t3phi_model = vis_to_t3(m.cvis_model, data.indx_t3_1, data.indx_t3_2, data.indx_t3_3)
    chi2_v2 = 0.0; chi2_t3amp = 0.0; chi2_t3phi = 0.0
    if (weights[1]>0) && (data.nv2>0)
        chi2_v2 = norm((v2_model - data.v2)./data.v2_err)^2
    end
    if (weights[2]>0) && (data.nt3amp>0)
        chi2_t3amp = norm((t3amp_model - data.t3amp)./data.t3amp_err)^2
    end
    if (weights[3]>0) && (data.nt3phi>0)
        if !vonmises
            chi2_t3phi = norm(mod360(t3phi_model - data.t3phi)./data.t3phi_err)^2
        else
            chi2_t3phi = sum(-2*data.t3phi_vonmises_err.*cos.((t3phi_model - data.t3phi)/180*pi) .+ data.t3phi_vonmises_chi2_offset)
        end
    end
    if verb
        printstyled(@sprintf("V2: %.2f ", chi2_v2/data.nv2), color=:red)
        printstyled(@sprintf("T3A: %.2f ", chi2_t3amp/data.nt3amp), color=:blue)
        printstyled(@sprintf("T3P: %.2f ", chi2_t3phi/data.nt3phi), color=:green)
        printstyled(@sprintf("Flux: %.4f ", sum(x)), color=:black)
    end
    return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

using NLopt
function optimize_sparco_parameters(params_start, x::AbstractMatrix{<:AbstractFloat}, ft, data;
        weights=[1.0,1.0,1.0], lb=[0.0, 0.0, 0.0, -20.0], ub=[1.0, 1.0, 1.0, 20.0])
    nparams = length(params_start)-1
    f_params = (params, _)->chi2_sparco_f(x, [params;params_start[end]], ft, data; verb=false, weights)
    optimizer = Opt(:LN_NELDERMEAD, nparams)
    min_objective!(optimizer, f_params)
    lower_bounds!(optimizer, lb)
    upper_bounds!(optimizer, ub)
    minchi2, params_opt, ret = optimize(optimizer, params_start[1:nparams])
    return minchi2, [params_opt;params_start[end]], ret
end

# Gradient of chi2 w.r.t. a scalar parameter, given dcvis_model/dparam.
function _sparco_param_grad(dcvis_model, cvis_model, v2_model, t3_model, t3amp_model, t3phi_model,
        data; weights=[1.0,1.0,1.0], vonmises=false)
    dv2 = 0.0; dt3amp = 0.0; dt3phi = 0.0
    if (weights[1]>0) && (data.nv2>0)
        dv2 = 4*sum(((v2_model-data.v2)./data.v2_err.^2).*real.(cvis_model[data.indx_v2].*dcvis_model[data.indx_v2]))
    end
    if (weights[2]>0) && (data.nt3amp>0)
        dt3amp = 2.0*sum(((t3amp_model-data.t3amp)./data.t3amp_err.^2).*(
            real(cvis_model[data.indx_t3_1].*conj(dcvis_model[data.indx_t3_1]))./abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3])
          + real(cvis_model[data.indx_t3_2].*conj(dcvis_model[data.indx_t3_2]))./abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3])
          + real(cvis_model[data.indx_t3_3].*conj(dcvis_model[data.indx_t3_3]))./abs.(cvis_model[data.indx_t3_3]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2])))
    end
    if (weights[3]>0) && (data.nt3phi>0)
        t3 = cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3]
        dt3 = dcvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3] +
              cvis_model[data.indx_t3_1].*dcvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3] +
              cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2].*dcvis_model[data.indx_t3_3]
        if !vonmises
            dt3phi = 360.0/pi*sum((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2).*imag(conj(t3).*dt3)./abs2.(t3))
        else
            dt3phi = sum(2*data.t3phi_vonmises_err.*sin.((t3phi_model - data.t3phi)/180*pi).*(pi/180).*imag(conj(t3).*dt3)./abs2.(t3))
        end
    end
    return weights[1]*dv2 + weights[2]*dt3amp + weights[3]*dt3phi
end

function chi2_sparco_fg(x::AbstractVector{<:AbstractFloat}, g::AbstractVector{<:AbstractFloat},
        ftplan::AbstractVector{<:NFFT.NFFTPlan}, data::OIdata, nparams::Int64;
        verb=true, weights=[1.0,1.0,1.0], vonmises=false)
    params = x[1:nparams]
    nx = Int(sqrt(length(x)-nparams))
    x_img = reshape(x[nparams+1:end], nx, nx)

    m = _sparco_model(params, x_img, ftplan, data)
    v2_model = vis_to_v2(m.cvis_model, data.indx_v2)
    t3_model, t3amp_model, t3phi_model = vis_to_t3(m.cvis_model, data.indx_t3_1, data.indx_t3_2, data.indx_t3_3)

    # Chi2 computation
    chi2_v2 = 0.0; chi2_t3amp = 0.0; chi2_t3phi = 0.0
    g_v2 = 0.0; g_t3amp = 0.0; g_t3phi = 0.0
    if (weights[1]>0) && (data.nv2>0)
        chi2_v2 = norm((v2_model - data.v2)./data.v2_err)^2
        g_v2 = real(adjoint(ftplan[3])*(4*((v2_model-data.v2)./data.v2_err.^2).*m.cvis_model[data.indx_v2].*m.imratio[data.indx_v2]))
    end
    if (weights[2]>0) && (data.nt3amp>0)
        chi2_t3amp = norm((t3amp_model - data.t3amp)./data.t3amp_err)^2
        dT3 = 2.0*(t3amp_model-data.t3amp)./(data.t3amp_err.^2)
        g_t3amp = real(adjoint(ftplan[4])*(dT3.*m.cvis_model[data.indx_t3_1].*m.imratio[data.indx_t3_1]./abs.(m.cvis_model[data.indx_t3_1]).*abs.(m.cvis_model[data.indx_t3_2]).*abs.(m.cvis_model[data.indx_t3_3]))) +
                  real(adjoint(ftplan[5])*(dT3.*m.cvis_model[data.indx_t3_2].*m.imratio[data.indx_t3_2]./abs.(m.cvis_model[data.indx_t3_2]).*abs.(m.cvis_model[data.indx_t3_1]).*abs.(m.cvis_model[data.indx_t3_3]))) +
                  real(adjoint(ftplan[6])*(dT3.*m.cvis_model[data.indx_t3_3].*m.imratio[data.indx_t3_3]./abs.(m.cvis_model[data.indx_t3_3]).*abs.(m.cvis_model[data.indx_t3_1]).*abs.(m.cvis_model[data.indx_t3_2])))
    end
    if (weights[3]>0) && (data.nt3phi>0)
        if !vonmises
            chi2_t3phi = norm(mod360(t3phi_model - data.t3phi)./data.t3phi_err)^2
            dT3 = -360.0/pi*mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2
        else
            chi2_t3phi = sum(-2*data.t3phi_vonmises_err.*cos.((t3phi_model - data.t3phi)/180*pi) .+ data.t3phi_vonmises_chi2_offset)
            dT3 = -2.0*data.t3phi_vonmises_err.*sin.((t3phi_model - data.t3phi)/180*pi)
        end
        g_t3phi = imag(adjoint(ftplan[4])*(dT3./abs2.(m.cvis_model[data.indx_t3_1]).*m.cvis_model[data.indx_t3_1].*m.imratio[data.indx_t3_1]) +
                       adjoint(ftplan[5])*(dT3./abs2.(m.cvis_model[data.indx_t3_2]).*m.cvis_model[data.indx_t3_2].*m.imratio[data.indx_t3_2]) +
                       adjoint(ftplan[6])*(dT3./abs2.(m.cvis_model[data.indx_t3_3]).*m.cvis_model[data.indx_t3_3].*m.imratio[data.indx_t3_3]))
    end
    if verb
        printstyled(@sprintf("V2: %.2f ", chi2_v2/data.nv2), color=:red)
        printstyled(@sprintf("T3A: %.2f ", chi2_t3amp/data.nt3amp), color=:blue)
        printstyled(@sprintf("T3P: %.2f ", chi2_t3phi/data.nt3phi), color=:green)
        printstyled(@sprintf("Flux: %.4f ", sum(x_img)), color=:black)
    end

    # Parameter gradients via quotient rule on composite visibility
    du = m.α.*m.Vstar .- m.β.*m.Venv                             # du/dfs
    dv = m.α .- m.β                                               # dv/dfs
    dcvis_dfs = (du.*m.v .- m.u.*dv)./(m.v.*m.v)
    du_bg = .-m.β.*m.Venv                                         # du/dfg
    dv_bg = m.α .- m.β                                            # dv/dfg
    dcvis_dfg = (du_bg.*m.v .- m.u.*dv_bg)./(m.v.*m.v)
    dVstar = dvisibility_ud([params[3]], data.uv)
    dcvis_dD = (m.fluxstar.*dVstar)./m.v
    λ = data.uv_lam; logλr = log.(λ/params[5])
    du_ind = logλr.*m.fluxenv.*m.Venv
    dv_ind = logλr.*m.fluxenv
    dcvis_dind = (du_ind.*m.v .- m.u.*dv_ind)./(m.v.*m.v)

    pgkw = (; weights, vonmises)
    g[1] = _sparco_param_grad(dcvis_dfs,  m.cvis_model, v2_model, t3_model, t3amp_model, t3phi_model, data; pgkw...)
    g[2] = _sparco_param_grad(dcvis_dfg,  m.cvis_model, v2_model, t3_model, t3amp_model, t3phi_model, data; pgkw...)
    g[3] = _sparco_param_grad(dcvis_dD,   m.cvis_model, v2_model, t3_model, t3amp_model, t3phi_model, data; pgkw...)
    g[4] = _sparco_param_grad(dcvis_dind, m.cvis_model, v2_model, t3_model, t3amp_model, t3phi_model, data; pgkw...)
    g[5] = 0.0  # λ₀ is fixed

    # Image pixel gradient
    g[nparams+1:end] = vec(weights[1]*g_v2 .+ weights[2]*g_t3amp .+ weights[3]*g_t3phi)
    return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

function crit_sparco_fg(x::AbstractVector{<:AbstractFloat}, g::AbstractVector{<:AbstractFloat},
        ftplan::AbstractVector{<:NFFT.NFFTPlan}, data::OIdata, nparams::Int64;
        weights=[1.0,1.0,1.0], printcolor=:black, regularizers=[], verb=true, vonmises=false)
    chi2 = chi2_sparco_fg(x, g, ftplan, data, nparams, weights=weights, verb=verb, vonmises=vonmises)
    nx = Int(sqrt(length(x)-nparams))
    reg_g = zeros(eltype(x), nx, nx)
    reg_f = regularization(reshape(x[nparams+1:end], nx, nx), reg_g, regularizers=regularizers, printcolor=printcolor, verb=verb)
    g[nparams+1:end] += vec(reg_g)
    # Gradient correction for the image (parameters are left untouched)
    flux = sum(x[nparams+1:end])
    g[nparams+1:end] = (g[nparams+1:end] .- sum(vec(x[nparams+1:end]).*g[nparams+1:end]) / flux) / flux
    return chi2 + reg_f
end

function reconstruct_sparco_gray(x_start::AbstractMatrix{<:AbstractFloat}, params_start::AbstractVector{<:AbstractFloat},
        data::OIdata, ft; printcolor=:black, verb=false, maxiter=100,
        regularizers=[], weights=[1.0,1.0,1.0], vonmises=false,
        ftol=(0,1e-8), xtol=(0,1e-8), gtol=(0,1e-8))
    crit = (x,g)->crit_sparco_fg(x, g, ft, data, length(params_start),
        regularizers=regularizers, verb=verb, weights=weights, vonmises=vonmises)
    sol = OptimPackNextGen.vmlmb(crit, [params_start;vec(x_start)], verb=verb,
        lower=0, maxiter=maxiter, blmvm=false, xtol=xtol, ftol=ftol, gtol=gtol)
    nparams = length(params_start)
    return (sol[1:nparams], reshape(sol[nparams+1:end], size(x_start)))
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
