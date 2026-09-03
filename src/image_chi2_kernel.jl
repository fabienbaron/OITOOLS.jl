# ===========================================================================
# Image-reconstruction chi2 + gradient kernel.
#
# One adjoint transform, not one per observable: every observable accumulates into a single
# adjoint source `g_cvis`, which is transformed back to the image domain once. That is the same
# construction model fitting uses (`cvis_to_chi2_fg`), so the two halves of the package share
# one definition of the observables instead of two.
#
# The adjoint source accumulates at `W`, independent of the plan's precision `T`. This is not
# an optimisation: visibility-domain partial sums are large and cancel under the transform, and
# rounding them to Float32 first costs 20-60% of the gradient at nx = 512. `W = Float64`
# against a Float32 plan matches a full-Float64 reference to three or four significant figures.
# ===========================================================================

"""
    ImageChi2Cache{T,W}(nx, nuv)
    ImageChi2Cache(cell, data; work=Float64)

Preallocated scratch for repeated [`image_chi2_fg!`](@ref) calls against one Fourier cell.

`T` is the cell's precision; `W` the precision the adjoint source accumulates in. `W === T`
aliases the working visibilities onto the transform output, so a same-precision cache costs
nothing extra.

Reusing one cache across an optimisation keeps the criterion's allocations flat in `nx`: at
nx = 512 that is 276 KiB per evaluation against 29 MB for a fresh workspace each time.
"""
struct ImageChi2Cache{T<:AbstractFloat, W<:AbstractFloat}
    xn     :: Matrix{Complex{T}}   # normalised image, forward transform input
    V      :: Vector{Complex{T}}   # model visibilities at every uv point
    Vw     :: Vector{Complex{W}}   # V at working precision (aliases V when W === T)
    g_cvis :: Vector{Complex{W}}   # adjoint source, accumulated over all observables
    src    :: Vector{Complex{T}}   # g_cvis prepared for the adjoint, at plan precision
    gimg   :: Matrix{Complex{T}}   # adjoint output
    wts    :: Vector{W}            # weights padded to 7
end

function ImageChi2Cache{T,W}(nx::Integer, nuv::Integer) where {T<:AbstractFloat,W<:AbstractFloat}
    V  = zeros(Complex{T}, nuv)
    Vw = W === T ? V : zeros(Complex{W}, nuv)
    return ImageChi2Cache{T,W}(zeros(Complex{T}, nx, nx), V, Vw,
                               zeros(Complex{W}, nuv), zeros(Complex{T}, nuv),
                               zeros(Complex{T}, nx, nx), zeros(W, 7))
end

ImageChi2Cache(cell::NFFTCell{T}, data::OIdata; work::Type{W}=Float64) where {T,W} =
    ImageChi2Cache{T,W}(cell.uv.N[1], size(data.uv, 2))
ImageChi2Cache(cell::DFTCell{T}, data::OIdata; work::Type{W}=Float64) where {T,W} =
    ImageChi2Cache{T,W}(isqrt(size(cell, 2)), size(cell, 1))
ImageChi2Cache(cell::AbstractVector{<:NFFTPlan}, data::OIdata; work::Type{W}=Float64) where {W} =
    ImageChi2Cache{ft_eltype(cell),W}(cell[1].N[1], size(data.uv, 2))

_cache_fits(c::ImageChi2Cache, nx::Integer, nuv::Integer) =
    size(c.xn, 1) == nx && length(c.V) == nuv

# Forward and adjoint, dispatching on the cell type. The two conventions differ by a
# conjugation, which `cvis_to_chi2_fg`'s docstring states:
#   NFFT:  g_image = real(NFFT^H * conj(g_cvis))     -- adjoint(plan) conjugates
#   DFT:   g_image = real(transpose(DFT) * g_cvis)
_forward!(ws::ImageChi2Cache, cell::NFFTCell) = mul!(ws.V, cell.uv, ws.xn)
_forward!(ws::ImageChi2Cache, cell::AbstractVector{<:NFFTPlan}) = mul!(ws.V, cell[1], ws.xn)
_forward!(ws::ImageChi2Cache, cell::DFTCell) = mul!(ws.V, cell.kernel, vec(ws.xn))

_prep_src!(ws::ImageChi2Cache{T}, ::Union{NFFTCell,AbstractVector{<:NFFTPlan}}) where {T} =
    (@inbounds @simd for i in eachindex(ws.src, ws.g_cvis); ws.src[i] = Complex{T}(conj(ws.g_cvis[i])); end)
_prep_src!(ws::ImageChi2Cache{T}, ::DFTCell) where {T} =
    (@inbounds @simd for i in eachindex(ws.src, ws.g_cvis); ws.src[i] = Complex{T}(ws.g_cvis[i]); end)

_adjoint!(ws::ImageChi2Cache, cell::NFFTCell) = mul!(ws.gimg, adjoint(cell.uv), ws.src)
_adjoint!(ws::ImageChi2Cache, cell::AbstractVector{<:NFFTPlan}) = mul!(ws.gimg, adjoint(cell[1]), ws.src)
_adjoint!(ws::ImageChi2Cache, cell::DFTCell) = mul!(vec(ws.gimg), transpose(cell.kernel), ws.src)

"""
    image_chi2_fg!(g, x, cell, data, ws; weights, vonmises=false, cvis=[], verb=false) -> chi2

Chi-squared of image `x` against `data`, writing the gradient into `g`.

Fits V², T3amp, T3phi, visamp, visphi and OI_FLUX — slots 1 to 6 of `weights`. Slot 7,
differential phase, is cross-channel and is handled by `_polychromatic_vis_gradient!`.

`g` is the gradient with respect to the **normalised** image, so the caller's flux correction
applies unchanged. Pass `cvis` to receive the model visibilities at `data.indx_vis`.

!!! note "OI_FLUX cannot constrain a flux-normalised image"
    The flux gradient is uniform, and normalising projects uniform vectors out exactly. When
    the flux is uncalibrated the fitted scaling also makes chi2 independent of total flux
    outright. Weighting slot 6 therefore reports a chi2 but moves a normalised image by nothing.
"""
function image_chi2_fg!(g::AbstractMatrix{<:AbstractFloat},
                        x::AbstractMatrix{<:AbstractFloat},
                        cell, data::OIdata, ws::ImageChi2Cache{T,W};
                        weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                        vonmises::Bool = false,
                        cvis = [],
                        verb::Bool = false) where {T,W}

    flux = sum(x)

    @inbounds @simd for i in eachindex(ws.xn, x)
        ws.xn[i] = Complex{T}(x[i] / flux)
    end
    _forward!(ws, cell)
    if W !== T
        @inbounds @simd for i in eachindex(ws.Vw, ws.V)
            ws.Vw[i] = Complex{W}(ws.V[i])
        end
    end
    length(cvis) > 0 && (cvis[:] = ws.V[data.indx_vis])

    fill!(ws.wts, zero(W))
    @inbounds for i in 1:min(7, length(weights))
        ws.wts[i] = W(weights[i])
    end

    t = _chi2_terms(ws.Vw, data, ws.wts, vonmises)
    chi2 = _total_chi2(t, ws.wts)
    fill!(ws.g_cvis, zero(Complex{W}))
    _accumulate_g_cvis!(ws.g_cvis, ws.Vw, t, data, ws.wts)

    _prep_src!(ws, cell)
    _adjoint!(ws, cell)
    @inbounds @simd for i in eachindex(g, ws.gimg)
        g[i] = real(ws.gimg[i])
    end

    # OI_FLUX: the model is the image's own total flux, so its gradient is uniform.
    if ws.wts[6] > 0 && data.nflux > 0
        fl = _chi2_flux(flux, data)
        chi2 += ws.wts[6] * fl.chi2
        gflux = 2 * fl.C * sum(fl.residual) * ws.wts[6]
        @inbounds @simd for i in eachindex(g)
            g[i] += gflux
        end
    end

    verb && _print_chi2_terms(t, data, ws.wts, flux)
    return chi2
end

# Convenience method for callers with no cache to reuse: allocates one per call.
function image_chi2_fg!(g::AbstractMatrix{<:AbstractFloat},
                        x::AbstractMatrix{<:AbstractFloat},
                        cell, data::OIdata; work::Type=Float64, kwargs...)
    ws = ImageChi2Cache(cell, data; work=work)
    return image_chi2_fg!(g, x, cell, data, ws; kwargs...)
end

function _print_chi2_terms(t, data::OIdata, w::AbstractVector{<:Real}, flux::Real)
    w[1] > 0 && data.nv2 > 0     ? printstyled(@sprintf("V2: %.2f ", t.chi2_v2/data.nv2), color=:red) :
                                   printstyled("V2: (N/A) ", color=:normal)
    w[2] > 0 && data.nt3amp > 0  && printstyled(@sprintf("T3A: %.2f ", t.chi2_t3amp/data.nt3amp), color=:blue)
    w[3] > 0 && data.nt3phi > 0  && printstyled(@sprintf("T3P: %.2f ", t.chi2_t3phi/data.nt3phi), color=:green)
    w[4] > 0 && data.nvisamp > 0 && printstyled(@sprintf("VA: %.2f ", t.chi2_visamp/data.nvisamp), color=:cyan)
    w[5] > 0 && data.nvisphi > 0 && printstyled(@sprintf("VP: %.2f ", t.chi2_visphi/data.nvisphi), color=:magenta)
    printstyled(@sprintf("Flux: %.4f ", flux), color=:normal)
    return nothing
end

"""
    image_chi2_f(x, cell, data, ws; weights, vonmises=false, cvis=[], verb=false) -> chi2

Chi-squared of image `x` against `data`, without the gradient — one forward transform.

Fits the same slots as [`image_chi2_fg!`](@ref), so the objective a criterion reports and the
objective its gradient descends are one definition. Passing a weight vector longer than three
to a kernel that reads only three is how an image comes back ignoring an observable whose chi2
is being printed.
"""
function image_chi2_f(x::AbstractMatrix{<:AbstractFloat},
                      cell, data::OIdata, ws::ImageChi2Cache{T,W};
                      weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                      vonmises::Bool = false,
                      cvis = [],
                      verb::Bool = false) where {T,W}

    flux = sum(x)
    @inbounds @simd for i in eachindex(ws.xn, x)
        ws.xn[i] = Complex{T}(x[i] / flux)
    end
    _forward!(ws, cell)
    if W !== T
        @inbounds @simd for i in eachindex(ws.Vw, ws.V)
            ws.Vw[i] = Complex{W}(ws.V[i])
        end
    end
    length(cvis) > 0 && (cvis[:] = ws.V[data.indx_vis])

    fill!(ws.wts, zero(W))
    @inbounds for i in 1:min(7, length(weights))
        ws.wts[i] = W(weights[i])
    end

    t = _chi2_terms(ws.Vw, data, ws.wts, vonmises)
    chi2 = _total_chi2(t, ws.wts)
    if ws.wts[6] > 0 && data.nflux > 0
        chi2 += ws.wts[6] * _chi2_flux(flux, data).chi2
    end
    verb && _print_chi2_terms(t, data, ws.wts, flux)
    return chi2
end

function image_chi2_f(x::AbstractMatrix{<:AbstractFloat}, cell, data::OIdata;
                      work::Type=Float64, kwargs...)
    return image_chi2_f(x, cell, data, ImageChi2Cache(cell, data; work=work); kwargs...)
end
