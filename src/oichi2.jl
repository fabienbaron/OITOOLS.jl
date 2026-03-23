using Statistics,LinearAlgebra, SparseArrays, SpecialFunctions, NFFT, Match, Printf



function _ft_cell_info(cell)
    if cell isa AbstractVector{<:NFFTPlan}
        plan = cell[1]  # fftplan_uv has all UV points
        nx = plan.N[1]
        nuv = size(plan.k, 2)
        return "NFFT", nx, nuv
    elseif cell isa AbstractMatrix{<:Complex}
        nuv = size(cell, 1)
        nx = round(Int, sqrt(size(cell, 2)))
        return "DFT", nx, nuv
    else
        return nothing, 0, 0
    end
end

function Base.display(ft::Vector{<:AbstractArray{<:NFFTPlan}})
    nmulti = length(ft)
    println("OITOOLS NFFT transform for $nmulti channels");
    for i in 1:nmulti
        _, nx, nuv = _ft_cell_info(ft[i])
        println("  Channel $i: $(nx)×$(nx) image, $nuv UV points")
    end
end

function Base.display(ft::AbstractMatrix)
    nwav, nepoch = size(ft)
    info = _ft_cell_info(ft[1])
    info[1] === nothing && (invoke(Base.display, Tuple{Any}, ft); return)
    mode, nx, nuv = info
    ncells = nwav * nepoch
    # Collect all UV counts
    nuv_all = [_ft_cell_info(ft[w, t])[3] for w in 1:nwav, t in 1:nepoch]
    println("OITOOLS $mode Fourier transform plans: $nwav wavelength(s) × $nepoch epoch(s)")
    println("  Image size: $(nx)×$(nx) pixels")
    if ncells == 1
        println("  UV points: $nuv")
    else
        println("  UV points: $(minimum(nuv_all))–$(maximum(nuv_all)) (total $(sum(nuv_all)))")
    end
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
    Base.depwarn("`setup_nfft_multiepochs` is deprecated, use `setup_ft(data, nx, pixsize)` instead.", :setup_nfft_multiepochs)
    nepochs = length(data)
    fftplan_multi = Vector{Vector{NFFTPlan}}(undef, nepochs)
    for i in 1:nepochs
        fftplan_multi[i] = setup_nfft(data[i], nx, pixsize)
    end
    return fftplan_multi
end

function setup_nfft_polychromatic(data::AbstractVecOrMat{<:OIdata}, nx, pixsize)
    Base.depwarn("`setup_nfft_polychromatic` is deprecated, use `setup_ft(data, nx, pixsize)` instead.", :setup_nfft_polychromatic)
    nwavs = size(data, 1)
    fftplan_multi = Vector{Vector{NFFTPlan}}(undef, nwavs)
    for i in 1:nwavs
        fftplan_multi[i] = setup_nfft(data[i], nx, pixsize)
    end
    return fftplan_multi
end

function setup_dft_polychromatic(data::AbstractVecOrMat{<:OIdata{T}}, nx, pixsize) where T
    Base.depwarn("`setup_dft_polychromatic` is deprecated, use `setup_ft(data, nx, pixsize; mode=\"dft\")` instead.", :setup_dft_polychromatic)
    nwavs = size(data, 1)
    fftplan_multi = Vector{Matrix{Complex{T}}}(undef, nwavs)
    for i in 1:nwavs
        fftplan_multi[i] = setup_dft(data[i], nx, pixsize)
    end
    return fftplan_multi
end


"""
    setup_ft(data, nx, pixsize; mode="nfft")

Set up Fourier transform plans for the given data matrix.

`data` is a `Matrix{OIdata}` of size `(nwav, nepoch)` (as returned by `readoifits`
or `readoifits_multiepochs`). Returns a matching matrix of FT plans.

# Keywords
- `mode="nfft"` — use `"nfft"` for Non-uniform FFT (fast) or `"dft"` for direct DFT (exact)
"""
function setup_ft(data::AbstractMatrix{<:OIdata}, nx, pixsize; mode::String="nfft")
    nwav, nepoch = size(data)
    if mode == "nfft"
        return [setup_nfft(data[w,t], nx, pixsize) for w in 1:nwav, t in 1:nepoch]
    elseif mode == "dft"
        return [setup_dft(data[w,t], nx, pixsize) for w in 1:nwav, t in 1:nepoch]
    else
        error("Unknown mode \"$mode\" — use \"nfft\" or \"dft\"")
    end
end

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
    cvis_model = image_to_vis(x, ft)

    # V²
    v2_model = data.nv2 > 0 ? vis_to_v2(cvis_model, data.indx_v2) : Float64[]

    # Triple products
    if data.nt3amp > 0 || data.nt3phi > 0
        _, t3amp_model, t3phi_model = vis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2, data.indx_t3_3)
    else
        t3amp_model = Float64[]
        t3phi_model = Float64[]
    end

    # Visibility amplitude & phase
    if data.nvisamp > 0 || data.nvisphi > 0
        Vvis = cvis_model[data.indx_vis]
        visamp_model = abs.(Vvis)
        visphi_model = angle.(Vvis) .* (180.0 / π)
    else
        visamp_model = Float64[]
        visphi_model = Float64[]
    end

    return (v2=v2_model, t3amp=t3amp_model, t3phi=t3phi_model,
            visamp=visamp_model, visphi=visphi_model)
end

"""
    image_to_residuals(x, ft, data)

Compute normalised residuals `(model - data) / error` for each observable type.
Returns a NamedTuple with fields `v2`, `t3amp`, `t3phi`, `visamp`, `visphi`.
Phase residuals (t3phi, visphi) are wrapped to [-180, 180] before dividing by error.
"""
function image_to_residuals(x, ft, data)
    obs = image_to_obs(x, ft, data)

    v2_res     = data.nv2 > 0 ? (obs.v2 .- data.v2) ./ data.v2_err : Float64[]
    t3amp_res  = data.nt3amp > 0 ? (obs.t3amp .- data.t3amp) ./ data.t3amp_err : Float64[]
    t3phi_res  = data.nt3phi > 0 ? mod360(obs.t3phi .- data.t3phi) ./ data.t3phi_err : Float64[]
    visamp_res = data.nvisamp > 0 ? (obs.visamp .- data.visamp) ./ data.visamp_err : Float64[]
    visphi_res = data.nvisphi > 0 ? mod360(obs.visphi .- data.visphi) ./ data.visphi_err : Float64[]

    return (v2=v2_res, t3amp=t3amp_res, t3phi=t3phi_res,
            visamp=visamp_res, visphi=visphi_res)
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
    # Transpectral TV: x is (npix, npix, nwavs)
    # Forward differences along dim 3, no wrap-around at boundaries
    dx = diff(x, dims=3)  # (npix, npix, nwavs-1)
    d = sqrt.(dx.^2 .+ ζ^2)
    f = sum(d .- ζ)
    # Gradient: chain rule scatters back to full (npix, npix, nwavs)
    q = dx ./ d  # (npix, npix, nwavs-1)
    g[:,:,1]       = -q[:,:,1]
    g[:,:,2:end-1] = q[:,:,1:end-1] .- q[:,:,2:end]
    g[:,:,end]     = q[:,:,end]
    return f
end

function trans_tvsq(x, g; verb=false)
    # Transpectral quadratic TV: x is (npix, npix, nwavs)
    dx = diff(x, dims=3)  # (npix, npix, nwavs-1)
    f = sum(dx.^2)
    g[:,:,1]       = -2 .* dx[:,:,1]
    g[:,:,2:end-1] = 2 .* (dx[:,:,1:end-1] .- dx[:,:,2:end])
    g[:,:,end]     = 2 .* dx[:,:,end]
    return f
end

function trans_l1l2(x, g; verb=false, ζ=1e-13, δ=2.0)
    # Transpectral L1-L2: x is (npix, npix, nwavs)
    dx = diff(x, dims=3)  # (npix, npix, nwavs-1)
    d = sqrt.(dx.^2 .+ ζ^2)
    # f = Σ (d/δ - log(1 + d/δ))
    f = sum(d ./ δ .- log.(1.0 .+ d ./ δ))
    # df/d(dx) = (dx/d) · (1/δ - 1/(δ+d)) = (dx/d) · d/(δ(δ+d)) = dx/(δ(δ+d))
    q = dx ./ (δ .* (δ .+ d))
    g[:,:,1]       = -q[:,:,1]
    g[:,:,2:end-1] = q[:,:,1:end-1] .- q[:,:,2:end]
    g[:,:,end]     = q[:,:,end]
    return f
end

# Spectral polynomial prior: penalizes deviation of each pixel's spectrum
# from its best-fit polynomial of degree `degree`.
# x is (npix, npix, nwavs), λ is a vector of nwavs wavelengths.
# f = Σ_pixel ||spectrum - polyfit||², g = 2(spectrum - polyfit) scattered back.
function trans_polychromatic_polyfit(x, g; degree::Int=1, λ::Vector{Float64}=Float64[])
    nx, ny, nwavs = size(x)
    if isempty(λ)
        λ = collect(range(0.0, 1.0, length=nwavs))
    end
    # Build Vandermonde matrix: V[i,j] = λ_i^(j-1), normalized
    λn = (λ .- mean(λ)) ./ max(std(λ), eps())
    V = hcat([λn.^k for k in 0:degree]...)  # nwavs × (degree+1)
    # Projection onto polynomial subspace: P = V (VᵀV)⁻¹ Vᵀ
    P = V * ((V' * V) \ V')  # nwavs × nwavs
    R = I - P  # residual projector
    # For each pixel, residual = R * spectrum, f = ||residual||²
    # Reshape x to (npix*npix, nwavs) for matrix ops
    X = reshape(x, nx*ny, nwavs)  # each row is a pixel's spectrum
    Xres = X * R'  # residuals: (npix*npix, nwavs)
    f = sum(Xres.^2)
    g[:] = reshape(2.0 .* Xres * R, nx, ny, nwavs)
    return f
end

# Spectral group spatial TV (vectorial/collaborative TV).
# Encourages edges to appear at the same spatial locations across all wavelengths.
# f = Σ_{i,j} √( Σ_λ [ (x_{i+1,j,λ} - x_{i,j,λ})² + (x_{i,j+1,λ} - x_{i,j,λ})² ] + ε² ) - ε
# x is (nx, ny, nwavs).
function trans_grouptv(x, g; ε=1e-12)
    nx, ny, nwavs = size(x)
    # Spatial forward differences with zero-padding (no wrap)
    dx = zeros(nx, ny, nwavs)  # ∂x/∂i
    dy = zeros(nx, ny, nwavs)  # ∂x/∂j
    dx[1:end-1,:,:] = x[2:end,:,:] .- x[1:end-1,:,:]
    dy[:,1:end-1,:] = x[:,2:end,:] .- x[:,1:end-1,:]
    # Sum of squared gradients across all wavelengths at each pixel
    S = dropdims(sum(dx.^2 .+ dy.^2, dims=3), dims=3) .+ ε^2  # (nx, ny)
    d = sqrt.(S)  # (nx, ny)
    f = sum(d .- ε)
    # Gradient: ∂f/∂(dx_{i,j,λ}) = dx_{i,j,λ} / d_{i,j}
    #           ∂f/∂(dy_{i,j,λ}) = dy_{i,j,λ} / d_{i,j}
    # Then scatter back through the forward-difference operator
    qdx = dx ./ d  # (nx, ny, nwavs), broadcasting d over dim 3
    qdy = dy ./ d
    # g = -div(qdx, qdy): adjoint of the forward-difference operator
    g .= 0.0
    # dx contribution: g[i,j] -= qdx[i,j]; g[i+1,j] += qdx[i,j]
    g[1:end-1,:,:] .-= qdx[1:end-1,:,:]
    g[2:end,:,:]   .+= qdx[1:end-1,:,:]
    # dy contribution: g[i,j] -= qdy[i,j]; g[i,j+1] += qdy[i,j]
    g[:,1:end-1,:] .-= qdy[:,1:end-1,:]
    g[:,2:end,:]   .+= qdy[:,1:end-1,:]
    return f
end

function regularization(x, reg_g; printcolor = :normal, regularizers=[], verb=true) # compound regularization
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
function chi2_f(x::AbstractMatrix{<:AbstractFloat}, dft::AbstractMatrix{<:Complex}, data::OIdata; weights = [1.0,1.0,1.0],  cvis = [], printcolor =:normal, verb=true, vonmises=false)
    Base.depwarn("`chi2_f` is deprecated, use `image_to_chi2` instead.", :chi2_f)
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
        printstyled(@sprintf("Flux: %.4f ", flux), color=:normal);
    end
    return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

# NFFT version
function chi2_f(x::AbstractMatrix{<:AbstractFloat}, ftplan::AbstractVector{<:NFFT.NFFTPlan}, data::OIdata; weights = [1.0,1.0,1.0], cvis = [], printcolor =:normal,  verb = false, vonmises=false)
    Base.depwarn("`chi2_f` is deprecated, use `image_to_chi2` instead.", :chi2_f)
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
        printstyled(@sprintf("Flux: %.4f ", flux), color=:normal);
    end
    return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

# DFT version
function chi2_fg(x::AbstractMatrix{<:AbstractFloat}, g::AbstractMatrix{<:AbstractFloat}, dft::AbstractMatrix{<:Complex}, data::OIdata; weights = [1.0,1.0,1.0],  cvis = [],  printcolor =:normal, verb=true, vonmises=false)
    Base.depwarn("`chi2_fg` is deprecated, use `image_to_chi2_fg` instead.", :chi2_fg)
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
        printstyled(@sprintf("Flux: %.4f ", flux), color=:normal);
    end
    return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

#NFFT version
function chi2_fg(x::AbstractMatrix{<:AbstractFloat}, g::AbstractMatrix{<:AbstractFloat}, ftplan::AbstractVector{<:NFFT.NFFTPlan}, data::OIdata; weights = [1.0,1.0,1.0], cvis = [], printcolor =:normal,  verb = false, vonmises=false)
    Base.depwarn("`chi2_fg` is deprecated, use `image_to_chi2_fg` instead.", :chi2_fg)
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
            printstyled("V2: (N/A) ", color=:normal) # disabled
        end
        printstyled(@sprintf("T3A: %.2f ", chi2_t3amp/data.nt3amp), color=:blue);
        printstyled(@sprintf("T3P: %.2f ", chi2_t3phi/data.nt3phi), color=:green);
        printstyled(@sprintf("Flux: %.4f ", flux), color=:normal);
    end
    return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

# ===========================================================================
# image_to_chi2 — unified chi2 (no regularization, no ndof normalization)
# ===========================================================================

"""
    image_to_chi2(x::AbstractArray{T,4}, ft::AbstractMatrix, data::AbstractMatrix{<:OIdata}; ...)

Compute the weighted chi-squared of a 4-D image cube `x[px, py, wav, epoch]`
against a `Matrix{OIdata}` of size `(nwav, nepoch)`.

Returns the **raw** chi-squared (not divided by ndof).  For `nwav > 1`,
cross-channel observables (differential phases, visibility amplitudes,
OI_FLUX) are included automatically based on data metadata.

See also: [`image_to_chi2_fg`](@ref), `image_to_obs`.
"""
function image_to_chi2(x::AbstractArray{<:AbstractFloat,4},
                       ft::AbstractMatrix,
                       data::AbstractMatrix{<:OIdata};
                       weights = [1.0, 1.0, 1.0],
                       use_diffphases = false,
                       verb = false,
                       vonmises = false)
    nwav, nepoch = size(data)
    npix = size(x, 1)
    chi2_total = 0.0

    for t in 1:nepoch
        # Per-channel chi2
        nwav_t = nwav
        amptyp = data[1,t].amptyp
        phityp = data[1,t].phityp
        udp = use_diffphases || phityp == "differential"
        udva = amptyp == "differential"
        uava = amptyp == "absolute"

        cvis = fill(ComplexF64[], nwav_t)
        need_cvis = udp || udva || uava
        if need_cvis
            for i in 1:nwav_t
                cvis[i] = Vector{ComplexF64}(undef, length(data[i,t].indx_vis))
            end
        end

        for w in 1:nwav_t
            verb && nepoch > 1 && printstyled("Epoch $t ", color=:normal)
            verb && nwav_t > 1 && printstyled("Channel $w ", color=:normal)
            if need_cvis
                chi2_total += chi2_f(x[:,:,w,t], ft[w,t], data[w,t];
                    cvis=cvis[w], verb=verb, weights=weights, vonmises=vonmises)
            else
                chi2_total += chi2_f(x[:,:,w,t], ft[w,t], data[w,t];
                    verb=verb, weights=weights, vonmises=vonmises)
            end
            # Per-channel absolute VISAMP (handled inside chi2_polychromatic_f but not in chi2_f)
            if uava && !udva && !isempty(cvis[w])
                va_model = abs.(cvis[w])
                chi2_va = norm((va_model - data[w,t].visamp) ./ data[w,t].visamp_err)^2
                chi2_total += chi2_va
                verb && printstyled(@sprintf("VA: %.2f ", chi2_va/data[w,t].nvisamp), color=:magenta)
            end
            verb && print("\n")
        end

        # Cross-channel terms for this epoch (nwav > 1 only)
        if nwav_t > 1
            data_epoch = @view data[:, t]
            vis_flux = _polychromatic_vis_flux_chi2(
                (@view x[:,:,:,t]), data_epoch, cvis, nwav_t;
                use_diffphases=udp, use_diffvisamp=udva,
                use_abs_visamp=false, verb=verb)
            chi2_total += vis_flux.chi2
        end
    end

    return chi2_total
end

"""
    image_to_chi2(x::AbstractMatrix, ft, data::OIdata; ...)

Mono convenience: compute chi-squared for a single 2-D image against a single `OIdata`.
"""
function image_to_chi2(x::AbstractMatrix{<:AbstractFloat}, ft, data::OIdata;
                       weights = [1.0, 1.0, 1.0], verb = false, vonmises = false)
    return chi2_f(x, ft, data; weights=weights, verb=verb, vonmises=vonmises)
end

# ===========================================================================
# image_to_chi2_fg — unified chi2 + gradient (with flux correction)
# ===========================================================================

"""
    image_to_chi2_fg(x::AbstractArray{T,4}, g::AbstractArray{T,4},
                     ft::AbstractMatrix, data::AbstractMatrix{<:OIdata}; ...)

Compute chi-squared and its gradient w.r.t. the raw (unnormalised) image pixels.

The gradient includes the flux-normalisation chain-rule correction:
for each `(w,t)` cell, `g[:,:,w,t]` is corrected so that
`g = (g_raw .- ⟨y, g_raw⟩) / S` where `y = x/S` and `S = sum(x)`.

Returns the **raw** chi-squared (not divided by ndof).
"""
function image_to_chi2_fg(x::AbstractArray{<:AbstractFloat,4},
                          g::AbstractArray{<:AbstractFloat,4},
                          ft::AbstractMatrix,
                          data::AbstractMatrix{<:OIdata};
                          weights = [1.0, 1.0, 1.0],
                          use_diffphases = false,
                          verb = false,
                          vonmises = false,
                          normalize = true)
    nwav, nepoch = size(data)
    npix = size(x, 1)
    g .= 0
    chi2_total = 0.0

    for t in 1:nepoch
        nwav_t = nwav
        amptyp = data[1,t].amptyp
        phityp = data[1,t].phityp
        udp = use_diffphases || phityp == "differential"
        udva = amptyp == "differential"
        uava = amptyp == "absolute"

        cvis = fill(ComplexF64[], nwav_t)
        need_cvis = udp || udva || uava
        if need_cvis
            for i in 1:nwav_t
                cvis[i] = Vector{ComplexF64}(undef, length(data[i,t].indx_vis))
            end
        end

        # Per-channel chi2 + gradient (raw adjoint, no flux correction)
        for w in 1:nwav_t
            subg = zeros(eltype(x), npix, npix)
            verb && nepoch > 1 && printstyled("Epoch $t ", color=:normal)
            verb && nwav_t > 1 && printstyled("Channel $w ", color=:normal)
            if need_cvis
                chi2_total += chi2_fg(x[:,:,w,t], subg, ft[w,t], data[w,t];
                    cvis=cvis[w], verb=verb, weights=weights, vonmises=vonmises)
            else
                chi2_total += chi2_fg(x[:,:,w,t], subg, ft[w,t], data[w,t];
                    verb=verb, weights=weights, vonmises=vonmises)
            end
            g[:,:,w,t] = subg
            verb && print("\n")
        end

        # Cross-channel gradient for this epoch (nwav > 1 only)
        if nwav_t > 1
            data_epoch = @view data[:, t]
            g_epoch = @view g[:,:,:,t]

            # Build vis_adjoint closure based on plan type
            ft_epoch = @view ft[:, t]
            if ft_epoch[1] isa AbstractVector{<:NFFTPlan}
                vis_adjoint = (rhs, i) -> vec(adjoint(ft_epoch[i][2]) * rhs)
            else
                vis_adjoint = (rhs, i) -> vec(conj.(transpose(ft_epoch[i][data_epoch[i].indx_vis, :]) * conj.(rhs)))
            end

            chi2_total += _polychromatic_vis_gradient!(
                (@view x[:,:,:,t]), g_epoch, data_epoch, cvis, nwav_t, npix, vis_adjoint;
                use_diffphases=udp, use_diffvisamp=udva,
                use_abs_visamp=uava, verb=verb)
        end

        # Apply per-cell flux correction
        if normalize
            for w in 1:nwav_t
                x_wt = @view x[:,:,w,t]
                g_wt = @view g[:,:,w,t]
                flux = sum(x_wt)
                g_wt .= (g_wt .- sum(x_wt .* g_wt) / flux) ./ flux
            end
        end
    end

    return chi2_total
end

"""
    image_to_chi2_fg(x::AbstractMatrix, g::AbstractMatrix, ft, data::OIdata; ...)

Mono convenience: compute chi-squared and flux-corrected gradient for a single
2-D image against a single `OIdata`.
"""
function image_to_chi2_fg(x::AbstractMatrix{<:AbstractFloat},
                          g::AbstractMatrix{<:AbstractFloat},
                          ft, data::OIdata;
                          weights = [1.0, 1.0, 1.0], verb = false, vonmises = false,
                          normalize = true)
    chi2 = chi2_fg(x, g, ft, data; weights=weights, verb=verb, vonmises=vonmises)
    if normalize
        flux = sum(x)
        g .= (g .- sum(x .* g) / flux) ./ flux
    end
    return chi2
end

# ===========================================================================
# Auto-detection convenience: accept Matrix{OIdata} + 2D image directly.
# For 1×1 data (mono), extracts the single element.
# For multi-channel data, wraps 2D image in 4D and calls unified method.
# ===========================================================================

function image_to_chi2(x::AbstractMatrix{<:AbstractFloat}, ft::AbstractMatrix,
                       data::AbstractMatrix{<:OIdata}; kwargs...)
    nwav, nepoch = size(data)
    if nwav == 1 && nepoch == 1
        return image_to_chi2(x, ft[1], data[1]; kwargs...)
    else
        nx = size(x, 1)
        return image_to_chi2(reshape(x, nx, nx, 1, 1),
            reshape([ft[1]], 1, 1), reshape([data[1]], 1, 1); kwargs...)
    end
end

function image_to_chi2_fg(x::AbstractMatrix{<:AbstractFloat}, g::AbstractMatrix{<:AbstractFloat},
                          ft::AbstractMatrix, data::AbstractMatrix{<:OIdata}; kwargs...)
    nwav, nepoch = size(data)
    if nwav == 1 && nepoch == 1
        return image_to_chi2_fg(x, g, ft[1], data[1]; kwargs...)
    else
        nx = size(x, 1)
        return image_to_chi2_fg(reshape(x, nx, nx, 1, 1), reshape(g, nx, nx, 1, 1),
            reshape([ft[1]], 1, 1), reshape([data[1]], 1, 1); kwargs...)
    end
end

function image_to_obs(x::AbstractMatrix, ft::AbstractMatrix, data::AbstractMatrix{<:OIdata})
    return image_to_obs(x, ft[1], data[1])
end

function image_to_residuals(x::AbstractMatrix, ft::AbstractMatrix, data::AbstractMatrix{<:OIdata})
    return image_to_residuals(x, ft[1], data[1])
end

function reconstruct(x_start::AbstractMatrix{<:AbstractFloat},
                     data::AbstractMatrix{<:OIdata}, ft::AbstractMatrix; kwargs...)
    npix = size(x_start, 1)
    x4 = reshape(x_start, npix, npix, 1, 1)
    x_sol = reconstruct(x4, data, ft; kwargs...)
    return x_sol[:,:,1,1]
end

function crit_fg(x::AbstractMatrix{<:AbstractFloat}, g::AbstractMatrix{<:AbstractFloat},
                 ft::AbstractMatrix, data::AbstractMatrix{<:OIdata}; kwargs...)
    nx = size(x, 1)
    g4 = reshape(g, nx, nx, 1, 1)
    f = crit_fg(reshape(x, nx, nx, 1, 1), g4, ft, data; kwargs...)
    return f
end

function crit_f(x::AbstractMatrix{<:AbstractFloat}, ft::AbstractMatrix,
                data::AbstractMatrix{<:OIdata}; kwargs...)
    nx = size(x, 1)
    return crit_f(reshape(x, nx, nx, 1, 1), ft, data; kwargs...)
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
    # These require the same baselines across all channels; we match by (sta1, sta2, mjd).
    if use_diffphases || use_diffvisamp
        # Build per-channel key→index maps: (sta1, sta2, round(mjd,6)) → position
        vis_keys = Vector{Vector{Tuple{Int,Int,Float64}}}(undef, nwavs)
        for i in 1:nwavs
            d = data[i]; nv = length(d.indx_vis)
            vis_keys[i] = [(d.vis_sta_index[1,j], d.vis_sta_index[2,j],
                            round(d.vis_mjd[j], digits=6)) for j in 1:nv]
        end
        # Common keys present in every channel (preserving order from channel 1)
        common_set = reduce(intersect, Set.(vis_keys))
        # Per-channel index masks for common baselines (ordered as in each channel)
        cidx = Vector{Vector{Int}}(undef, nwavs)
        for i in 1:nwavs
            cidx[i] = [j for (j, k) in enumerate(vis_keys[i]) if k in common_set]
        end
        ncommon = length(cidx[1])

        cvis_all = Matrix{ComplexF64}(undef, ncommon, nwavs)
        for i in 1:nwavs
            cvis_all[:, i] = cvis[i][cidx[i]]
        end
        cvis_sum = vec(sum(cvis_all, dims=2))
        if use_diffphases
            diffphi_model = zeros(ncommon, nwavs)
            cvis_ref_all  = zeros(ComplexF64, ncommon, nwavs)
            for i=1:nwavs
                cvis_ref_all[:,i] = (cvis_sum .- cvis_all[:,i]) ./ (nwavs - 1)
                diffphi_model[:,i] = angle.(cvis_all[:,i] ./ cvis_ref_all[:,i]) .* (180.0/pi)
            end
            diffphi     = hcat([data[i].visphi[cidx[i]] for i=1:nwavs]...)
            diffphi_err = hcat([data[i].visphi_err[cidx[i]] for i=1:nwavs]...)
            chi2_dp = norm(mod360(diffphi_model-diffphi)./diffphi_err)^2
            chi2 += chi2_dp;  ndof += length(diffphi_model)
            verb && printstyled(@sprintf("DiffPhi: %.2f ", chi2_dp/length(diffphi_model)), color=:magenta)
        end
        if use_diffvisamp
            diffvisamp_model = zeros(ncommon, nwavs)
            for i=1:nwavs
                cvis_ref = (cvis_sum .- cvis_all[:,i]) ./ (nwavs - 1)
                diffvisamp_model[:,i] = abs.(cvis_all[:,i] ./ cvis_ref)
            end
            mean_dva = vec(mean(diffvisamp_model, dims=2))
            diffvisamp_model ./= mean_dva
            diffvisamp     = hcat([data[i].visamp[cidx[i]] for i=1:nwavs]...)
            diffvisamp_err = hcat([data[i].visamp_err[cidx[i]] for i=1:nwavs]...)
            chi2_dva = norm((diffvisamp_model - diffvisamp)./diffvisamp_err)^2
            chi2 += chi2_dva;  ndof += length(diffvisamp_model)
            verb && printstyled(@sprintf("DiffVA: %.2f ", chi2_dva/length(diffvisamp_model)), color=:magenta)
        end
    end
    # Absolute VISAMP — channels may have different sizes, accumulate per-channel
    if use_abs_visamp && !use_diffvisamp
        chi2_ava = 0.0; ndof_ava = 0
        for i in 1:nwavs
            va_model = abs.(cvis[i])
            va_data  = data[i].visamp
            va_err   = data[i].visamp_err
            chi2_ava += norm((va_model - va_data) ./ va_err)^2
            ndof_ava += data[i].nvisamp
        end
        chi2 += chi2_ava;  ndof += ndof_ava
        verb && printstyled(@sprintf("VA: %.2f ", chi2_ava/ndof_ava), color=:magenta)
    end
    # OI_FLUX
    C_flux = 1.0; flux_residual = Float64[]; flux_chan_idx = Int[]
    use_flux = any(data[i].nflux > 0 for i in 1:nwavs)
    if use_flux
        fm = Float64[]; fd = Float64[]; fe = Float64[]
        for i in 1:nwavs
            d = data[i]; d.nflux > 0 || continue
            fi = sum(x[:,:,i])
            append!(fm, fill(fi, d.nflux)); append!(fd, d.flux)
            append!(fe, d.flux_err);        append!(flux_chan_idx, fill(i, d.nflux))
        end
        calibrated = data[1].flux_calibrated
        if calibrated
            flux_residual = (fm .- fd) ./ fe.^2
            chi2_fl = sum(((fm .- fd) ./ fe).^2)
            C_flux = 1.0
        else
            w = 1.0 ./ fe.^2
            C_flux = sum(fm .* fd .* w) / sum(fm.^2 .* w)
            flux_residual = (C_flux .* fm .- fd) ./ fe.^2
            chi2_fl = sum(((C_flux .* fm .- fd) ./ fe).^2)
        end
        chi2 += chi2_fl;  ndof += length(fd)
        if verb
            if calibrated
                printstyled(@sprintf("OI_FLUX: %.2f ", chi2_fl/length(fd)), color=:yellow)
            else
                printstyled(@sprintf("OI_FLUX: %.2f (C=%.4f) ", chi2_fl/length(fd), C_flux), color=:yellow)
            end
        end
    end
    return (; chi2, ndof, diffphi_model, diffphi, diffphi_err, cvis_ref_all,
              C_flux, flux_residual, flux_chan_idx, use_flux)
end

# ===========================================================================
# chi2_polychromatic_f — NFFT version
# ===========================================================================
function chi2_polychromatic_f(x::AbstractArray{<:AbstractFloat,3}, ft::AbstractVector{<:AbstractVector{<:NFFTPlan}}, data::AbstractVector{<:OIdata};weights = [1.0,1.0,1.0], printcolor= [], use_diffphases = false, verb = false)
    Base.depwarn("`chi2_polychromatic_f` is deprecated, use `image_to_chi2` with 4D images and Matrix{OIdata} instead.", :chi2_polychromatic_f)
    nwavs = length(ft);
    npix = size(x,1);
    if printcolor == []
        printcolor=[ :normal for i=1:nwavs]
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
    ndof = 0
    for i=1:nwavs
        if verb == true
            printstyled("Channel $i ",color=printcolor[i]);
        end

        if use_diffphases || use_diffvisamp || use_abs_visamp
            f[i] = chi2_f(x[:,:,i], ft[i], data[i], cvis=cvis[i], verb = verb, weights = weights);
        else
            f[i] = chi2_f(x[:,:,i], ft[i], data[i], verb = verb, weights = weights);
        end
        ndof_i = (weights[1]>0 ? data[i].nv2 : 0) + (weights[2]>0 ? data[i].nt3amp : 0) + (weights[3]>0 ? data[i].nt3phi : 0)

        # Per-channel absolute VISAMP
        if use_abs_visamp && !use_diffvisamp && !isempty(cvis[i])
            va_model = abs.(cvis[i])
            chi2_va = norm((va_model - data[i].visamp) ./ data[i].visamp_err)^2
            f[i] += chi2_va;  ndof_i += data[i].nvisamp
            verb && printstyled(@sprintf("VA: %.2f ", chi2_va/data[i].nvisamp), color=:magenta)
        end

        # Per-channel OI_FLUX (raw, no global C yet — printed after global fit below)
        ndof += ndof_i
        if verb == true
            printstyled(@sprintf("Chi2r: %.2f Chi2: %.2f\n", f[i]/ndof_i, f[i]),color=printcolor[i]);
        end
    end
    chi2f = sum(f)

    # Cross-channel differential observables + OI_FLUX
    vis_flux = _polychromatic_vis_flux_chi2(x, data, cvis, nwavs;
        use_diffphases=use_diffphases, use_diffvisamp=use_diffvisamp,
        use_abs_visamp=false, verb=verb)
    chi2f += vis_flux.chi2
    ndof  += vis_flux.ndof

    if verb == true
        printstyled(@sprintf("FULLCHI2: %.2f  Chi2r: %.2f  (ndof: %d)\n", chi2f, chi2f/ndof, ndof), color=:red)
    end

    return chi2f/ndof;
end

# ===========================================================================
# chi2_polychromatic_f — DFT version
# ===========================================================================
function chi2_polychromatic_f(x::AbstractArray{<:AbstractFloat,3}, ft::AbstractVector{<:AbstractMatrix{<:Complex}}, data::AbstractVector{<:OIdata};weights = [1.0,1.0,1.0], printcolor= [], use_diffphases = false, verb = false)
    Base.depwarn("`chi2_polychromatic_f` is deprecated, use `image_to_chi2` with 4D images and Matrix{OIdata} instead.", :chi2_polychromatic_f)
    nwavs = length(ft);
    npix = size(x,1);
    if printcolor == []
        printcolor=[ :normal for i=1:nwavs]
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
    ndof = 0
    for i=1:nwavs
        if verb == true
            printstyled("Channel $i ",color=printcolor[i]);
        end

        if use_diffphases || use_diffvisamp || use_abs_visamp
            f[i] = chi2_f(x[:,:,i], ft[i], data[i], cvis=cvis[i], verb = verb, weights = weights);
        else
            f[i] = chi2_f(x[:,:,i], ft[i], data[i], verb = verb, weights = weights);
        end
        ndof_i = (weights[1]>0 ? data[i].nv2 : 0) + (weights[2]>0 ? data[i].nt3amp : 0) + (weights[3]>0 ? data[i].nt3phi : 0)

        # Per-channel absolute VISAMP
        if use_abs_visamp && !use_diffvisamp && !isempty(cvis[i])
            va_model = abs.(cvis[i])
            chi2_va = norm((va_model - data[i].visamp) ./ data[i].visamp_err)^2
            f[i] += chi2_va;  ndof_i += data[i].nvisamp
            verb && printstyled(@sprintf("VA: %.2f ", chi2_va/data[i].nvisamp), color=:magenta)
        end

        ndof += ndof_i
        if verb == true
            printstyled(@sprintf("Chi2r: %.2f Chi2: %.2f\n", f[i]/ndof_i, f[i]),color=printcolor[i]);
        end
    end
    chi2f = sum(f)

    # Cross-channel differential observables + OI_FLUX
    vis_flux = _polychromatic_vis_flux_chi2(x, data, cvis, nwavs;
        use_diffphases=use_diffphases, use_diffvisamp=use_diffvisamp,
        use_abs_visamp=false, verb=verb)
    chi2f += vis_flux.chi2
    ndof  += vis_flux.ndof

    if verb == true
        printstyled(@sprintf("FULLCHI2: %.2f  Chi2r: %.2f  (ndof: %d)\n", chi2f, chi2f/ndof, ndof), color=:red)
    end

    return chi2f/ndof;
end

function crit_fg(x,g::AbstractMatrix{<:AbstractFloat}, ft, data::OIdata; weights = [1.0,1.0,1.0], cvis = [], printcolor = :normal, regularizers=[], verb = true)
    chi2 = chi2_fg(x, g, ft, data, cvis = cvis, verb = verb, weights = weights);
    reg = regularization(x, g, regularizers=regularizers, printcolor = printcolor, verb = verb);
    flux = sum(x)
    g[:] = (g .- sum(x.*g) / flux ) / flux; # gradient correction to take into account the non-normalized image
    return chi2 + reg;
end

function crit_f(x::AbstractMatrix{<:AbstractFloat}, fftplan::AbstractVector{<:NFFT.NFFTPlan}, data::OIdata; weights = [1.0,1.0,1.0], cvis = [], printcolor = :normal, regularizers=[], verb = true)
    chi2 = chi2_f(x, fftplan, data, cvis = cvis, verb = verb, weights = weights );
    g = zeros(eltype(x), size(x));
    reg = regularization(x, g,  regularizers=regularizers, printcolor = printcolor, verb = verb);
    return chi2 + reg;
end

# ===========================================================================
# Regularization helpers for unified crit_fg / crit_f
# ===========================================================================

# Expand regularizers to per-channel form.
# Input can be:
#   []                          → empty per-channel
#   [["tv",μ], ...]             → same regs for all channels
#   [[["tv",μ]], [["l1",μ]]]   → per-channel lists (length == nwav)
function _per_channel_regs(regularizers, nwav)
    isempty(regularizers) && return fill([], nwav)
    # If first element is a String, it's a single reg spec like ["tv", 1e-3]
    # i.e. regularizers is a flat list of reg tuples
    if regularizers[1] isa AbstractString || (regularizers[1] isa AbstractVector && regularizers[1][1] isa AbstractString)
        # Check if it's a single tuple ["tv", 1e-3] vs list of tuples [["tv", 1e-3], ...]
        if regularizers[1] isa AbstractString
            # Single reg spec like ["tv", 1e-3] — wrap in list
            return fill([regularizers], nwav)
        else
            # List of reg specs like [["tv", 1e-3], ["l1l2", 1e-3, 0.01]]
            return fill(regularizers, nwav)
        end
    elseif regularizers[1] isa AbstractVector
        # Per-channel lists
        length(regularizers) >= nwav || error("Per-channel regularizers must have length ≥ nwav ($nwav)")
        return regularizers
    else
        return fill(regularizers, nwav)
    end
end

# 2-D (spatial) regularization: applied to each (wav, epoch) cell independently
function _spatial_reg!(x3d, g3d, nwav, regularizers; verb=false)
    isempty(regularizers) && return 0.0
    T = eltype(x3d); npix = size(x3d, 1)
    cell_regs = _per_channel_regs(regularizers, nwav)
    f = 0.0
    for w in 1:nwav
        isempty(cell_regs[w]) && continue
        g_reg = zeros(T, npix, npix)
        f += regularization((@view x3d[:,:,w]), g_reg;
            regularizers=cell_regs[w], verb=verb)
        g3d[:,:,w] .+= g_reg
    end
    return f
end

# Transspectral regularization: cross-wavelength penalty per epoch
function _transspectral_reg!(x3d, g3d, nwav, transspectral_regularizers, ndof;
        verb=false, data=nothing)
    (nwav > 1 && !isempty(transspectral_regularizers)) || return 0.0
    npix = size(x3d, 1)
    regs = vcat(fill([], nwav), [transspectral_regularizers])
    # _polychromatic_transspectral_reg! takes accumulated f, returns f + Δf
    return _polychromatic_transspectral_reg!(x3d, g3d, 0.0, ndof, npix, nwav, regs;
        verb=verb, data=data)
end

# Temporal regularization: cross-epoch penalty
function _temporal_reg!(x4, g4, temporal_regularizers; verb=false)
    isempty(temporal_regularizers) && return 0.0
    nepoch = size(x4, 4)
    nepoch > 1 || return 0.0
    T = eltype(x4)
    npixall = size(x4,1) * size(x4,2) * size(x4,3)
    y = reshape(x4, npixall, nepoch)
    f = 0.0
    for treg in temporal_regularizers
        rname, μ = treg[1], treg[2]
        if rname == "temporal_tvsq"
            temporalf = sum((y[:,2:end] .- y[:,1:end-1]).^2)
            tv_g = zeros(T, npixall, nepoch)
            if nepoch > 2
                tv_g[:,1]       .= 2 .* (y[:,1] .- y[:,2])
                tv_g[:,2:end-1] .= 4 .* y[:,2:end-1] .- 2 .* (y[:,1:end-2] .+ y[:,3:end])
                tv_g[:,end]     .= 2 .* (y[:,end] .- y[:,end-1])
            else
                tv_g[:,1] .= 2 .* (y[:,1] .- y[:,2])
                tv_g[:,2] .= 2 .* (y[:,2] .- y[:,1])
            end
            f += μ * temporalf
            g4[:] .+= μ .* vec(tv_g)
            verb && printstyled(@sprintf("Temporal TV²: %.3f\n", μ * temporalf), color=:yellow)
        end
    end
    return f
end

# ===========================================================================
# Unified crit_fg and crit_f
# ===========================================================================

"""
    crit_fg(x::AbstractArray{T,4}, g::AbstractArray{T,4},
            ft::AbstractMatrix, data::AbstractMatrix{<:OIdata}; ...)

Unified criterion: `(χ² + regularization) / ndof`, with gradient.

The flux-normalisation chain-rule correction is applied to the combined
(χ² + regularization) gradient per cell, matching the original behaviour.

Returns the criterion value.  `g` is modified in place.

# Regularization structure

Three categories of regularizers, each a list of `["name", μ, ...]` tuples:

- **2-D / spatial** (`regularizers`): applied independently to each `(wavelength, epoch)`
  cell, e.g. `[["tv", 1e-3], ["l1l2", 1e-3, 0.01]]`
- **Transspectral** (`transspectral_regularizers`): cross-wavelength penalties applied
  per epoch when `nwav > 1`, e.g. `[["transspectral_tv", 0.1]]`
- **Temporal** (`temporal_regularizers`): cross-epoch penalties applied when
  `nepoch > 1`, e.g. `[["temporal_tvsq", 0.01]]`

# Keywords
- `weights = [1.0, 1.0, 1.0]`: relative weights for (V², T3amp, T3phi)
- `regularizers = []`: 2-D per-cell regularizers
- `transspectral_regularizers = []`: cross-wavelength regularizers
- `temporal_regularizers = []`: cross-epoch regularizers
- `epochs_weights = []`: per-epoch scaling (default: uniform)
- `use_diffphases = false`: force differential-phase fitting
- `verb = false`: verbose output
- `vonmises = false`: von Mises loss for closure phases
"""
function crit_fg(x4::AbstractArray{<:AbstractFloat,4},
                 g4::AbstractArray{<:AbstractFloat,4},
                 ft::AbstractMatrix,
                 data::AbstractMatrix{<:OIdata};
                 weights = [1.0, 1.0, 1.0],
                 regularizers = [],
                 transspectral_regularizers = [],
                 temporal_regularizers = [],
                 epochs_weights = [],
                 use_diffphases = false,
                 verb = false,
                 vonmises = false)
    nwav, nepoch = size(data)
    npix = size(x4, 1)
    T = eltype(x4)

    if isempty(epochs_weights)
        epochs_weights = ones(T, nepoch)
    end

    # Total degrees of freedom
    ndof = 0.0
    for t in 1:nepoch, w in 1:nwav
        d = data[w,t]
        ndof += weights[1]*d.nv2 + weights[2]*d.nt3amp + weights[3]*d.nt3phi
    end
    ndof = max(ndof, 1.0)

    f = 0.0
    g4 .= zero(T)

    for t in 1:nepoch
        g_epoch = zeros(T, npix, npix, nwav)

        verb && nepoch > 1 && printstyled("Epoch $t ", color=:cyan)

        # Chi2 gradient (raw adjoint, no flux correction yet)
        f_epoch = 0.0
        for w in 1:nwav
            subg = zeros(T, npix, npix)
            verb && nwav > 1 && printstyled("Channel $w ", color=:normal)
            f_epoch += chi2_fg((@view x4[:,:,w,t]), subg, ft[w,t], data[w,t];
                verb=verb, weights=weights, vonmises=vonmises)
            g_epoch[:,:,w] .= subg
            verb && print("\n")
        end

        # 2-D regularization (per cell)
        f_epoch += _spatial_reg!((@view x4[:,:,:,t]), g_epoch, nwav, regularizers; verb=verb)

        # Transspectral regularization (cross-wavelength)
        f_epoch += _transspectral_reg!((@view x4[:,:,:,t]), g_epoch, nwav,
            transspectral_regularizers, ndof; verb=verb, data=(@view data[:,t]))

        # Apply flux correction to combined (chi2 + reg) gradient per cell
        for w in 1:nwav
            x_wt = @view x4[:,:,w,t]
            g_wt = @view g_epoch[:,:,w]
            flux = sum(x_wt)
            g_wt .= (g_wt .- sum(x_wt .* g_wt) / flux) ./ flux
        end

        # Apply epoch weight and accumulate
        f += epochs_weights[t] * f_epoch
        g4[:,:,:,t] .+= epochs_weights[t] .* g_epoch
    end

    # Temporal regularization (cross-epoch)
    f += _temporal_reg!(x4, g4, temporal_regularizers; verb=verb)

    # ndof normalization
    f /= ndof
    g4 ./= ndof
    verb && printstyled(@sprintf("Crit/dof: %.4f\n", f), color=:blue)
    return f
end

"""
    crit_f(x::AbstractArray{T,4}, ft::AbstractMatrix,
           data::AbstractMatrix{<:OIdata}; ...)

Unified criterion: `(χ² + regularization) / ndof`, forward-only (no gradient).

Skips adjoint/gradient computation for speed.  Accepts the same keywords as
`crit_fg` (minus `g`).
"""
function crit_f(x4::AbstractArray{<:AbstractFloat,4},
                ft::AbstractMatrix,
                data::AbstractMatrix{<:OIdata};
                weights = [1.0, 1.0, 1.0],
                regularizers = [],
                transspectral_regularizers = [],
                temporal_regularizers = [],
                epochs_weights = [],
                use_diffphases = false,
                verb = false,
                vonmises = false)
    nwav, nepoch = size(data)
    npix = size(x4, 1)
    T = eltype(x4)

    if isempty(epochs_weights)
        epochs_weights = ones(T, nepoch)
    end

    ndof = 0.0
    for t in 1:nepoch, w in 1:nwav
        d = data[w,t]
        ndof += weights[1]*d.nv2 + weights[2]*d.nt3amp + weights[3]*d.nt3phi
    end
    ndof = max(ndof, 1.0)

    cell_regs = _per_channel_regs(regularizers, nwav)

    f = 0.0
    for t in 1:nepoch
        x_t = reshape((@view x4[:,:,:,t]), npix, npix, nwav, 1)
        ft_t = reshape((@view ft[:,t]), nwav, 1)
        data_t = reshape((@view data[:,t]), nwav, 1)

        verb && nepoch > 1 && printstyled("Epoch $t ", color=:cyan)

        # Chi2 forward-only (no gradient, no adjoint FFT)
        f_epoch = image_to_chi2(x_t, ft_t, data_t;
            weights=weights, use_diffphases=use_diffphases,
            verb=verb, vonmises=vonmises)

        # 2-D regularization (forward-only, gradient discarded)
        for w in 1:nwav
            if !isempty(cell_regs[w])
                g_dummy = zeros(T, npix, npix)
                f_epoch += regularization((@view x4[:,:,w,t]), g_dummy;
                    regularizers=cell_regs[w], verb=verb)
            end
        end

        # Transspectral regularization (forward-only, gradient discarded)
        if nwav > 1 && !isempty(transspectral_regularizers)
            g_dummy_3d = zeros(T, npix, npix, nwav)
            regs = vcat(fill([], nwav), [transspectral_regularizers])
            f_epoch = _polychromatic_transspectral_reg!(
                (@view x4[:,:,:,t]), g_dummy_3d,
                f_epoch, ndof, npix, nwav, regs;
                verb=verb, data=(@view data[:,t]))
        end

        f += epochs_weights[t] * f_epoch
    end

    # Temporal regularization (forward-only)
    if nepoch > 1 && !isempty(temporal_regularizers)
        npixall = npix * npix * nwav
        y = reshape(x4, npixall, nepoch)
        for treg in temporal_regularizers
            rname, μ = treg[1], treg[2]
            if rname == "temporal_tvsq"
                f += μ * sum((y[:,2:end] .- y[:,1:end-1]).^2)
                verb && printstyled(@sprintf("Temporal TV²: %.3f\n", μ * sum((y[:,2:end] .- y[:,1:end-1]).^2)), color=:yellow)
            end
        end
    end

    f /= ndof
    verb && printstyled(@sprintf("Crit/dof: %.4f\n", f), color=:blue)
    return f
end

const image_to_crit = crit_f

function crit_multitemporal_fg(x::AbstractArray{<:AbstractFloat,3}, g::AbstractArray{<:AbstractFloat,3}, ft::AbstractVector{<:AbstractVector{<:NFFT.NFFTPlan}}, data::Array{OIdata,1};weights = [1.0,1.0,1.0], printcolor= [], epochs_weights=[],regularizers=[], verb = false)
    Base.depwarn("`crit_multitemporal_fg` is deprecated, use `crit_fg` with 4D images and Matrix{OIdata} instead.", :crit_multitemporal_fg)
    nepochs = length(ft);
    if epochs_weights == []
        epochs_weights=ones(eltype(x_start), nepochs);
    end
    if printcolor == []
        printcolor=Array{Symbol}(undef,nepochs);
        printcolor[:] .= :normal
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
        # Match baselines across channels by (sta1, sta2, mjd) — same logic as chi2 helper
        vis_keys = Vector{Vector{Tuple{Int,Int,Float64}}}(undef, nwavs)
        for i in 1:nwavs
            d = data[i]; nv = length(d.indx_vis)
            vis_keys[i] = [(d.vis_sta_index[1,j], d.vis_sta_index[2,j],
                            round(d.vis_mjd[j], digits=6)) for j in 1:nv]
        end
        common_set = reduce(intersect, Set.(vis_keys))
        cidx = Vector{Vector{Int}}(undef, nwavs)
        for i in 1:nwavs
            cidx[i] = [j for (j, k) in enumerate(vis_keys[i]) if k in common_set]
        end
        ncommon = length(cidx[1])

        cvis_all = Matrix{ComplexF64}(undef, ncommon, nwavs)
        for i in 1:nwavs
            cvis_all[:, i] = cvis[i][cidx[i]]
        end
        cvis_sum = vec(sum(cvis_all, dims=2))

        if use_diffphases
            diffphi_model = zeros(ncommon, nwavs)
            cvis_ref_all  = zeros(ComplexF64, ncommon, nwavs)
            for i=1:nwavs
                cvis_ref_all[:,i] = (cvis_sum .- cvis_all[:,i]) ./ (nwavs - 1)
                diffphi_model[:,i] = angle.(cvis_all[:,i] ./ cvis_ref_all[:,i]) .* (180.0/pi)
            end
            diffphi     = hcat([data[i].visphi[cidx[i]] for i=1:nwavs]...)
            diffphi_err = hcat([data[i].visphi_err[cidx[i]] for i=1:nwavs]...)
            chi2_diffphi = norm(mod360(diffphi_model-diffphi)./diffphi_err)^2
            f += chi2_diffphi
            verb && printstyled(@sprintf("DiffPhi: %.2f ", chi2_diffphi/length(diffphi_model)), color=:magenta)

            # Gradient: d(arg(V_i/V_ref_i))/d(pixel) via chain rule
            # Leading-order approximation: gradient through numerator V_i only
            # d(arg(V/Vref))/dpixel ≈ d(arg(V))/dpixel = Im((1/V) * dV/dpixel)
            # Following the T3phi pattern: rhs = scalar_weight .* V ./ |V|²
            for i=1:nwavs
                dT3 = -360.0/pi .* mod360(diffphi_model[:,i]-diffphi[:,i])./diffphi_err[:,i].^2
                # Scatter common-baseline rhs back to full-length vis vector for adjoint
                rhs_full = zeros(ComplexF64, length(cvis[i]))
                rhs_full[cidx[i]] = dT3 .* cvis_all[:,i] ./ abs2.(cvis_all[:,i])
                g_contrib = imag.(vis_adjoint_fn(rhs_full, i))
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
        calibrated = data[1].flux_calibrated
        if calibrated
            C_flux = 1.0
            residual_flux = (flux_model_vec .- flux_data_vec) ./ flux_err_vec.^2
            chi2_flux = sum(((flux_model_vec .- flux_data_vec) ./ flux_err_vec).^2)
        else
            w = 1.0 ./ flux_err_vec.^2
            C_flux = sum(flux_model_vec .* flux_data_vec .* w) / sum(flux_model_vec.^2 .* w)
            residual_flux = (C_flux .* flux_model_vec .- flux_data_vec) ./ flux_err_vec.^2
            chi2_flux = sum(((C_flux .* flux_model_vec .- flux_data_vec) ./ flux_err_vec).^2)
        end
        f += chi2_flux
        if verb
            if calibrated
                printstyled(@sprintf("OI_FLUX: %.2f ", chi2_flux/length(flux_data_vec)), color=:yellow)
            else
                printstyled(@sprintf("OI_FLUX: %.2f (C=%.4f) ", chi2_flux/length(flux_data_vec), C_flux), color=:yellow)
            end
        end
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
function _polychromatic_transspectral_reg!(x, g, f, ndof, npix, nwavs, regularizers; verb=false, data=nothing)
    if length(regularizers) > nwavs
        ntransreg = length(regularizers[nwavs+1])
        tg = zeros(npix, npix, nwavs)
        for i=1:ntransreg
            rname = regularizers[nwavs+1][i][1]
            μ     = regularizers[nwavs+1][i][2]
            if rname == "transspectral_tv"
                tf = trans_tv(x, tg)
                f += μ * tf;  g[:,:,:] += μ * tg
                verb && printstyled(@sprintf(" ts_tv: %.3f", μ * tf), color=:yellow)
            elseif rname == "transspectral_tvsq"
                tf = trans_tvsq(x, tg)
                f += μ * tf;  g[:,:,:] += μ * tg
                verb && printstyled(@sprintf(" ts_tvsq: %.3f", μ * tf), color=:yellow)
            elseif rname == "transspectral_structnorm"
                tf = trans_structnorm(x, tg)
                f += μ * tf;  g[:,:,:] += μ * tg
                verb && printstyled(@sprintf(" ts_struct: %.3f", μ * tf), color=:yellow)
            elseif rname == "transspectral_l1l2"
                δ = regularizers[nwavs+1][i][3]
                tf = trans_l1l2(x, tg; δ=δ)
                f += μ * tf;  g[:,:,:] += μ * tg
                verb && printstyled(@sprintf(" ts_l1l2: %.3f", μ * tf), color=:yellow)
            elseif rname == "transspectral_grouptv"
                tf = trans_grouptv(x, tg)
                f += μ * tf;  g[:,:,:] += μ * tg
                verb && printstyled(@sprintf(" ts_gtv: %.3f", μ * tf), color=:yellow)
            elseif rname == "transspectral_poly"
                deg = length(regularizers[nwavs+1][i]) > 2 ? Int(regularizers[nwavs+1][i][3]) : 1
                λ = _get_wavelengths(data, nwavs)
                tf = trans_polychromatic_polyfit(x, tg; degree=deg, λ=λ)
                f += μ * tf;  g[:,:,:] += μ * tg
                verb && printstyled(@sprintf(" ts_poly(d=%d): %.3f", deg, μ * tf), color=:yellow)
            end
        end
        verb && printstyled(@sprintf("\nPost trans -- Crit: %.2f Crit/dof: %.2f\n", f, f/ndof), color=:blue)
    end
    return f
end

# Extract wavelength vector from data for polynomial regularizer
function _get_wavelengths(data, nwavs)
    if !isnothing(data) && length(data) >= nwavs
        return Float64[mean(data[i].uv_lam) for i in 1:nwavs]
    else
        return collect(range(0.0, 1.0, length=nwavs))
    end
end

# ===========================================================================
# crit_polychromatic_fg — NFFT version
# ===========================================================================
function crit_polychromatic_fg(x::AbstractArray{<:AbstractFloat,3}, g::AbstractArray{<:AbstractFloat,3}, ft::AbstractVector{<:AbstractVector{<:NFFTPlan}}, data::AbstractVector{<:OIdata};weights = [1.0,1.0,1.0], printcolor= [], regularizers=[], use_diffphases = false, verb = false)
    Base.depwarn("`crit_polychromatic_fg` is deprecated, use `crit_fg` with 4D images and Matrix{OIdata} instead.", :crit_polychromatic_fg)
    nwavs = length(ft); npix = size(x,1)
    printcolor == [] && (printcolor = [:normal for i=1:nwavs])
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

    f = _polychromatic_transspectral_reg!(x, g, f, ndof, npix, nwavs, regularizers; verb=verb, data=data)
    g[:] = g[:]/ndof
    return f/ndof
end

# ===========================================================================
# crit_polychromatic_fg — DFT version
# ===========================================================================
function crit_polychromatic_fg(x::AbstractArray{<:AbstractFloat,3}, g::AbstractArray{<:AbstractFloat,3}, ft::AbstractVector{<:AbstractMatrix{<:Complex}}, data::AbstractVector{<:OIdata};weights = [1.0,1.0,1.0], printcolor= [], regularizers=[], use_diffphases = false, verb = false)
    Base.depwarn("`crit_polychromatic_fg` is deprecated, use `crit_fg` with 4D images and Matrix{OIdata} instead.", :crit_polychromatic_fg)
    nwavs = length(ft); npix = size(x,1)
    printcolor == [] && (printcolor = [:normal for i=1:nwavs])
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

    f = _polychromatic_transspectral_reg!(x, g, f, ndof, npix, nwavs, regularizers; verb=verb, data=data)
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

"""
    reconstruct(x_start::AbstractArray{T,4}, data::AbstractMatrix{<:OIdata},
                ft::AbstractMatrix; ...)

Unified image reconstruction from optical interferometric data using VMLMB.

`x_start` is the initial image of shape `(nx, nx, nwav, nepoch)`.
`data` is a `Matrix{OIdata}` of size `(nwav, nepoch)` (from `readoifits` or
`readoifits_multiepochs`).  `ft` is the matching `Matrix` of FT plans from
`setup_ft`.

The criterion minimised is `(χ² + regularization) / ndof`.  The flux-normalisation
chain-rule correction is applied to the combined (χ² + regularization) gradient
per cell.

# Keywords
- `weights = [1.0, 1.0, 1.0]`: relative weights for (V², T3amp, T3phi)
- `regularizers = []`: per-cell regularizers — a list of `["name", μ, ...]` tuples
  applied identically to every `(wavelength, epoch)` cell
- `transspectral_regularizers = []`: cross-wavelength regularizers (per epoch),
  e.g. `[["transspectral_tv", 0.1]]`
- `temporal_regularizers = []`: cross-epoch regularizers,
  e.g. `[["temporal_tvsq", 0.01]]`
- `epochs_weights = []`: per-epoch scaling (default: uniform)
- `use_diffphases = false`: force differential-phase fitting
- `verb = false`: verbose output
- `maxiter = 100`: maximum VMLMB iterations
- `vonmises = false`: von Mises loss for closure phases
- `ftol = (0, 1e-8)`, `xtol = (0, 1e-8)`, `gtol = (0, 1e-8)`: VMLMB tolerances
"""
function reconstruct(x_start::AbstractArray{<:AbstractFloat,4},
                     data::AbstractMatrix{<:OIdata},
                     ft::AbstractMatrix;
                     weights = [1.0, 1.0, 1.0],
                     regularizers = [],
                     transspectral_regularizers = [],
                     temporal_regularizers = [],
                     epochs_weights = [],
                     use_diffphases = false,
                     verb = false,
                     maxiter = 100,
                     vonmises = false,
                     ftol = (0, 1e-8),
                     xtol = (0, 1e-8),
                     gtol = (0, 1e-8))
    _crit = (x4, g4) -> crit_fg(x4, g4, ft, data;
        weights=weights, regularizers=regularizers,
        transspectral_regularizers=transspectral_regularizers,
        temporal_regularizers=temporal_regularizers,
        epochs_weights=epochs_weights,
        use_diffphases=use_diffphases,
        verb=verb, vonmises=vonmises)
    x_sol = OptimPackNextGen.vmlmb(_crit, x_start, verb=verb, lower=0,
        maxiter=maxiter, blmvm=false, xtol=xtol, ftol=ftol, gtol=gtol)
    return x_sol
end

"""
    reconstruct(x_start::AbstractMatrix, data::OIdata, ft; ...)

Monochromatic convenience wrapper — delegates to the unified 4-D method.
See the 4-D method for the full keyword list.
"""
function reconstruct(x_start::AbstractMatrix{<:AbstractFloat}, data::OIdata, ft;
                     weights = [1.0, 1.0, 1.0], printcolor = :normal,
                     verb = false, maxiter = 100,
                     regularizers = [], vonmises = false,
                     ftol = (0, 1e-8), xtol = (0, 1e-8), gtol = (0, 1e-8))
    npix = size(x_start, 1)
    x4 = reshape(x_start, npix, npix, 1, 1)
    data_m = reshape([data], 1, 1)
    ft_m = reshape([ft], 1, 1)
    x_sol = reconstruct(x4, data_m, ft_m;
        weights=weights, regularizers=regularizers,
        verb=verb, maxiter=maxiter, vonmises=vonmises,
        ftol=ftol, xtol=xtol, gtol=gtol)
    return x_sol[:,:,1,1]
end

function reconstruct_multitemporal(x_start::AbstractArray{<:AbstractFloat,3}, data::AbstractVector{<:OIdata}, ft; weights = [1.0,1.0,1.0], epochs_weights =[], printcolor= [], verb = true, maxiter = 100, regularizers =[], ftol= (0,1e-8), xtol=(0,1e-8), gtol=(0,1e-8))
    Base.depwarn("`reconstruct_multitemporal` is deprecated, use `reconstruct` with 4D images and Matrix{OIdata} instead.", :reconstruct_multitemporal)
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
    Base.depwarn("`reconstruct_polychromatic` is deprecated, use `reconstruct` with 4D images and Matrix{OIdata} instead.", :reconstruct_polychromatic)
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
        printstyled(@sprintf("Flux: %.4f ", sum(x)), color=:normal)
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
        printstyled(@sprintf("Flux: %.4f ", sum(x_img)), color=:normal)
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
        weights=[1.0,1.0,1.0], printcolor=:normal, regularizers=[], verb=true, vonmises=false)
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

"""
    reconstruct_sparco_gray(x_start, params_start, data, ft; kwargs...)

SPARCO grey image reconstruction with a single chromatic point source.
Jointly optimizes `[params; image_pixels]` via VMLMB.
Returns `(params_final, x_final)`.
"""
function reconstruct_sparco_gray(x_start::AbstractMatrix{<:AbstractFloat}, params_start::AbstractVector{<:AbstractFloat},
        data::OIdata, ft; printcolor=:normal, verb=false, maxiter=100,
        regularizers=[], weights=[1.0,1.0,1.0], vonmises=false,
        ftol=(0,1e-8), xtol=(0,1e-8), gtol=(0,1e-8))
    crit = (x,g)->crit_sparco_fg(x, g, ft, data, length(params_start),
        regularizers=regularizers, verb=verb, weights=weights, vonmises=vonmises)
    sol = OptimPackNextGen.vmlmb(crit, [params_start;vec(x_start)], verb=verb,
        lower=0, maxiter=maxiter, blmvm=false, xtol=xtol, ftol=ftol, gtol=gtol)
    nparams = length(params_start)
    return (sol[1:nparams], reshape(sol[nparams+1:end], size(x_start)))
end

# ---------------------------------------------------------------------------
# Multi-source SPARCO reconstruction
# ---------------------------------------------------------------------------
# Generalizes SPARCO to N chromatic point sources + grey environment image.
#
#        Σ_k f_k·α_k·V_k  +  f_env·β·V_env
# V_tot = -------------------------------------------
#        Σ_k f_k·α_k  +  f_bg·α_bg  +  f_env·β
#
# where α_k  = (λ/λ₀)^(-dindex_k)     (per-source chromatic scaling)
#       α_bg = (λ/λ₀)^(-4)             (Rayleigh–Jeans background)
#       β    = (λ/λ₀)^(d_env - 4)      (environment power law)
#       f_env= 1 - Σ_k f_k - f_bg
#
# Per-source parameters (5 each):
#   flux_k, dindex_k, ra_k (mas), dec_k (mas), diam_k (mas)
#
# Global parameters (3):
#   f_bg, d_env, λ₀  (last is fixed)
#
# Total nparams = 5*nsources + 3
# ---------------------------------------------------------------------------

const _MAS2RAD_SPARCO = π / (180.0 * 3600.0 * 1000.0)

function _sparco_multi_model(params, nsources, x_img, ftplan, data)
    nparams = 5*nsources + 3
    λ0 = params[nparams]        # reference wavelength (fixed)
    f_bg  = params[nparams-2]
    d_env = params[nparams-1]
    λ  = data.uv_lam
    λr = λ ./ λ0

    # Background and environment chromatic scaling
    α_bg = λr.^(-4.0)
    β    = λr.^(d_env - 4.0)

    # Accumulate source contributions
    nuv = length(λ)
    sum_fα   = f_bg .* α_bg                      # denominator accumulator
    sum_fαV  = zeros(ComplexF64, nuv)             # numerator accumulator

    Vsources = Vector{Vector{ComplexF64}}(undef, nsources)
    αsources = Vector{Vector{Float64}}(undef, nsources)
    fluxsources = Vector{Vector{Float64}}(undef, nsources)

    u = data.uv[1,:]
    v = data.uv[2,:]
    for k in 1:nsources
        off = 5*(k-1)
        fk     = params[off+1]
        dindex = params[off+2]
        ra     = params[off+3]   # mas
        dec    = params[off+4]   # mas
        diam   = params[off+5]   # mas

        αk = λr.^(-dindex)
        Vud = visibility_ud([diam], data.uv)
        phase = @. -2π * _MAS2RAD_SPARCO * (u * ra + v * dec)
        Vk = Vud .* cis.(phase)

        αsources[k] = αk
        Vsources[k] = Vk
        fαk = fk .* αk
        fluxsources[k] = fαk
        sum_fα  .+= fαk
        sum_fαV .+= fαk .* Vk
    end

    # Environment contribution
    f_env_scalar = 1.0 - sum(params[5*(j-1)+1] for j in 1:nsources) - f_bg
    fluxenv = f_env_scalar .* β
    Venv = image_to_vis(x_img, ftplan[1])
    sum_fαV .+= fluxenv .* Venv
    sum_fα  .+= fluxenv

    cvis_model = sum_fαV ./ sum_fα
    imratio = fluxenv ./ sum_fα

    return (; cvis_model, Venv, Vsources, αsources, fluxsources, α_bg, β,
              fluxenv, f_env_scalar, sum_fαV, sum_fα, imratio, λr, u, v)
end

"""
    chi2_sparco_multi_f(x, params, nsources, ftplan, data; kwargs...)

Compute chi-squared for the multi-source SPARCO model (forward only, no gradient).
`params` layout: 5 per source `[flux, dindex, ra, dec, diam]` + 3 global `[f_bg, d_env, λ₀]`.
"""
function chi2_sparco_multi_f(x::AbstractMatrix{<:AbstractFloat}, params::AbstractVector{<:AbstractFloat},
        nsources::Int, ftplan::AbstractVector{<:NFFT.NFFTPlan}, data::OIdata;
        verb=true, weights=[1.0,1.0,1.0], vonmises=false)
    m = _sparco_multi_model(params, nsources, x, ftplan, data)
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
        printstyled(@sprintf("Flux: %.4f ", sum(x)), color=:normal)
    end
    return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

"""
    optimize_sparco_multi_parameters(params_start, nsources, x, ft, data; kwargs...)

Optimize multi-source SPARCO parameters with the image held fixed.
Uses NelderMead (gradient-free). Returns `(minchi2, params_opt, ret)`.
"""
function optimize_sparco_multi_parameters(params_start, nsources::Int,
        x::AbstractMatrix{<:AbstractFloat}, ft, data;
        weights=[1.0,1.0,1.0], lb=nothing, ub=nothing)
    # Default bounds: per source [flux 0-1, dindex -20..20, ra -500..500, dec -500..500, diam 0..50]
    #                 global [f_bg 0..1, d_env -20..20]   (λ₀ is fixed)
    nfree = 5*nsources + 2  # exclude λ₀
    if lb === nothing
        lb = Float64[]
        for _ in 1:nsources
            append!(lb, [0.0, -20.0, -500.0, -500.0, 0.0])
        end
        append!(lb, [0.0, -20.0])
    end
    if ub === nothing
        ub = Float64[]
        for _ in 1:nsources
            append!(ub, [1.0, 20.0, 500.0, 500.0, 50.0])
        end
        append!(ub, [1.0, 20.0])
    end
    λ0 = params_start[end]
    f_params = (p, _)->chi2_sparco_multi_f(x, [p; λ0], nsources, ft, data; verb=false, weights=weights)
    optimizer = Opt(:LN_NELDERMEAD, nfree)
    min_objective!(optimizer, f_params)
    lower_bounds!(optimizer, lb)
    upper_bounds!(optimizer, ub)
    minchi2, params_opt, ret = optimize(optimizer, params_start[1:nfree])
    return minchi2, [params_opt; λ0], ret
end

function chi2_sparco_multi_fg(x::AbstractVector{<:AbstractFloat}, g::AbstractVector{<:AbstractFloat},
        ftplan::AbstractVector{<:NFFT.NFFTPlan}, data::OIdata, nsources::Int, nparams::Int;
        verb=true, weights=[1.0,1.0,1.0], vonmises=false)
    params = x[1:nparams]
    nx = Int(sqrt(length(x)-nparams))
    x_img = reshape(x[nparams+1:end], nx, nx)

    m = _sparco_multi_model(params, nsources, x_img, ftplan, data)
    v2_model = vis_to_v2(m.cvis_model, data.indx_v2)
    t3_model, t3amp_model, t3phi_model = vis_to_t3(m.cvis_model, data.indx_t3_1, data.indx_t3_2, data.indx_t3_3)

    # Chi2 + image gradient (same structure as single-source SPARCO)
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
        printstyled(@sprintf("Flux: %.4f ", sum(x_img)), color=:normal)
    end

    # Parameter gradients via quotient rule
    pgkw = (; weights, vonmises)
    logλr = log.(m.λr)

    for k in 1:nsources
        off = 5*(k-1)
        fk     = params[off+1]
        dindex = params[off+2]
        ra     = params[off+3]
        dec    = params[off+4]
        diam   = params[off+5]
        αk     = m.αsources[k]
        Vk     = m.Vsources[k]
        fαk    = m.fluxsources[k]

        # ∂/∂fk: affects both source k and f_env (= 1 - Σf - fbg)
        du_dfk = αk .* Vk .- m.β .* m.Venv         # d(numerator)/dfk
        dv_dfk = αk .- m.β                           # d(denominator)/dfk
        dcvis_dfk = (du_dfk .* m.sum_fα .- m.sum_fαV .* dv_dfk) ./ (m.sum_fα.^2)
        g[off+1] = _sparco_param_grad(dcvis_dfk, m.cvis_model, v2_model, t3_model, t3amp_model, t3phi_model, data; pgkw...)

        # ∂/∂dindex_k
        du_dind = -logλr .* fαk .* Vk
        dv_dind = -logλr .* fαk
        dcvis_dind = (du_dind .* m.sum_fα .- m.sum_fαV .* dv_dind) ./ (m.sum_fα.^2)
        g[off+2] = _sparco_param_grad(dcvis_dind, m.cvis_model, v2_model, t3_model, t3amp_model, t3phi_model, data; pgkw...)

        # ∂/∂ra_k
        dphase_dra = @. -2π * _MAS2RAD_SPARCO * m.u
        dVk_dra = Vk .* (im .* dphase_dra)
        dcvis_dra = (fαk .* dVk_dra) ./ m.sum_fα
        g[off+3] = _sparco_param_grad(dcvis_dra, m.cvis_model, v2_model, t3_model, t3amp_model, t3phi_model, data; pgkw...)

        # ∂/∂dec_k
        dphase_ddec = @. -2π * _MAS2RAD_SPARCO * m.v
        dVk_ddec = Vk .* (im .* dphase_ddec)
        dcvis_ddec = (fαk .* dVk_ddec) ./ m.sum_fα
        g[off+4] = _sparco_param_grad(dcvis_ddec, m.cvis_model, v2_model, t3_model, t3amp_model, t3phi_model, data; pgkw...)

        # ∂/∂diam_k
        dVud = dvisibility_ud([diam], data.uv)
        phase = @. -2π * _MAS2RAD_SPARCO * (m.u * ra + m.v * dec)
        dVk_dD = dVud .* cis.(phase)
        dcvis_dD = (fαk .* dVk_dD) ./ m.sum_fα
        g[off+5] = _sparco_param_grad(dcvis_dD, m.cvis_model, v2_model, t3_model, t3amp_model, t3phi_model, data; pgkw...)
    end

    # ∂/∂f_bg
    goff = 5*nsources
    du_dfbg = .- m.β .* m.Venv                       # f_env decreases
    dv_dfbg = m.α_bg .- m.β
    dcvis_dfbg = (du_dfbg .* m.sum_fα .- m.sum_fαV .* dv_dfbg) ./ (m.sum_fα.^2)
    g[goff+1] = _sparco_param_grad(dcvis_dfbg, m.cvis_model, v2_model, t3_model, t3amp_model, t3phi_model, data; pgkw...)

    # ∂/∂d_env
    du_denv = logλr .* m.fluxenv .* m.Venv
    dv_denv = logλr .* m.fluxenv
    dcvis_denv = (du_denv .* m.sum_fα .- m.sum_fαV .* dv_denv) ./ (m.sum_fα.^2)
    g[goff+2] = _sparco_param_grad(dcvis_denv, m.cvis_model, v2_model, t3_model, t3amp_model, t3phi_model, data; pgkw...)

    # λ₀ is fixed
    g[goff+3] = 0.0

    # Image pixel gradient
    g[nparams+1:end] = vec(weights[1]*g_v2 .+ weights[2]*g_t3amp .+ weights[3]*g_t3phi)
    return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

function crit_sparco_multi_fg(x::AbstractVector{<:AbstractFloat}, g::AbstractVector{<:AbstractFloat},
        ftplan::AbstractVector{<:NFFT.NFFTPlan}, data::OIdata, nsources::Int, nparams::Int;
        weights=[1.0,1.0,1.0], printcolor=:normal, regularizers=[], verb=true, vonmises=false)
    chi2 = chi2_sparco_multi_fg(x, g, ftplan, data, nsources, nparams, weights=weights, verb=verb, vonmises=vonmises)
    nx = Int(sqrt(length(x)-nparams))
    reg_g = zeros(eltype(x), nx, nx)
    reg_f = regularization(reshape(x[nparams+1:end], nx, nx), reg_g, regularizers=regularizers, printcolor=printcolor, verb=verb)
    g[nparams+1:end] += vec(reg_g)
    # Flux-normalization correction for image gradient
    flux = sum(x[nparams+1:end])
    g[nparams+1:end] = (g[nparams+1:end] .- sum(vec(x[nparams+1:end]).*g[nparams+1:end]) / flux) / flux
    return chi2 + reg_f
end

"""
    reconstruct_sparco_multi(x_start, params_start, nsources, data, ft; kwargs...)

Multi-source SPARCO image reconstruction with N chromatic point sources.
Jointly optimizes `[params; image_pixels]` via VMLMB.
Returns `(params_final, x_final)`.
"""
function reconstruct_sparco_multi(x_start::AbstractMatrix{<:AbstractFloat},
        params_start::AbstractVector{<:AbstractFloat}, nsources::Int,
        data::OIdata, ft; printcolor=:normal, verb=false, maxiter=100,
        regularizers=[], weights=[1.0,1.0,1.0], vonmises=false,
        ftol=(0,1e-8), xtol=(0,1e-8), gtol=(0,1e-8))
    nparams = length(params_start)
    crit = (x,g)->crit_sparco_multi_fg(x, g, ft, data, nsources, nparams,
        regularizers=regularizers, verb=verb, weights=weights, vonmises=vonmises)
    sol = OptimPackNextGen.vmlmb(crit, [params_start;vec(x_start)], verb=verb,
        lower=0, maxiter=maxiter, blmvm=false, xtol=xtol, ftol=ftol, gtol=gtol)
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
