# The Image perspective's data layer: one dataset in, one reconstruction out.
#
# Positivity first, deliberately. `reconstruct` optimises with VMLMB under `lower = 0`, so a
# non-negative image is a property of the optimiser rather than something a regulariser asks
# for — and that makes it the one setting always in force and never listed among the
# regularisers. Getting a run to work end to end with nothing else switched on is what proves
# the plumbing: geometry, Fourier plan, starting image, weights, and a criterion that actually
# descends. Regularisers are then a weight vector added to a working loop, not a debugging
# variable tangled up with it.
#
# Nothing here draws. As with `model.jl`, that is what lets it be tested against OITOOLS
# rather than against a screenshot.

"""
What a reconstruction needs before it can start.

`pixsize` is in mas, and `nx * pixsize` is the field of view. `auto_pixsize` samples the
longest baseline about three times over, which is the sensible default and rarely the final
answer — hence [`imaging_defaults`](@ref) rather than a hard-coded number.
"""
struct ImagingSetup
    nx        :: Int
    pixsize   :: Float64
    mode      :: Symbol        # :nfft or :dft
    startkind :: Symbol        # :dirac, :gaussian or :fits
    startfwhm :: Float64       # mas; 0 means "let gaussian_prior decide" (fov/5)
    startpath :: String
end

ImagingSetup(; nx = 64, pixsize = 0.25, mode = :nfft,
               startkind = :gaussian, startfwhm = 0.0, startpath = "") =
    ImagingSetup(nx, Float64(pixsize), Symbol(mode), Symbol(startkind),
                 Float64(startfwhm), String(startpath))

"""
    imaging_defaults(data; nx = 64) -> ImagingSetup

Geometry the data itself suggests: `auto_pixsize` for the pixel, `nx` left to the caller.

The field of view that falls out is `nx * pixsize`; whether it covers the source is a
judgement the panel shows rather than makes.
"""
function imaging_defaults(data; nx::Integer = 64)
    d = data isa AbstractArray ? data[1] : data
    ImagingSetup(; nx = Int(nx), pixsize = auto_pixsize(d))
end

"Field of view in mas."
fov(s::ImagingSetup) = s.nx * s.pixsize

"""
    start_image(setup, ft) -> Matrix

The image the optimiser starts from, in the precision of `ft`.

That last part is not a detail. The plans are `Float32` by default, and `image_to_vis`
dispatches on the element type matching the plan's — so a `Float64` starting image gives a
`MethodError` from deep inside the criterion rather than anything mentioning precision.
`reconstruct` casts internally, but every other call a panel makes (a χ² readout before the
run, say) would hit it, so the cast happens here, once.

A Dirac start is a single central pixel: the least committal image there is, and the usual
choice when nothing is known. A Gaussian start is smoother and converges faster on a
centrally-condensed source. Both carry unit total flux.
"""
function start_image(setup::ImagingSetup, ft)
    nx = setup.nx
    x = if setup.startkind === :dirac
        m = zeros(Float64, nx, nx)
        m[nx ÷ 2 + 1, nx ÷ 2 + 1] = 1.0
        m
    elseif setup.startkind === :gaussian
        fwhm = setup.startfwhm > 0 ? setup.startfwhm : fov(setup) / 5
        reshape(gaussian_prior(nx, setup.pixsize; fwhm_mas = fwhm), nx, nx)
    elseif setup.startkind === :fits
        isfile(setup.startpath) || error("No starting image at '$(setup.startpath)'")
        m = Float64.(readfits(setup.startpath))
        size(m) == (nx, nx) ||
            error("The starting image is $(size(m,1))×$(size(m,2)) but nx is $nx")
        s = sum(m)
        s > 0 || error("The starting image sums to $s; it must carry positive flux")
        m ./ s
    else
        error("Unknown starting image '$(setup.startkind)'; expected :dirac, :gaussian or :fits")
    end
    return to_ft_precision(x, ft)
end

"""
    imaging_weights(; v2, t3amp, t3phi) -> Vector{Float64}

The three-element weight vector `reconstruct` takes, from three tick boxes.

Three, not the seven model fitting uses. The difference is a standing source of confusion, and
building the vector here rather than in the panel is what keeps the panel from having to know
which engine wants which length.

`cvis` and `flux` have no entry at any weight: `reconstruct` cannot use them. That is a
property of the engine, not of the file, and the panel says so separately.
"""
imaging_weights(; v2::Bool = true, t3amp::Bool = true, t3phi::Bool = true) =
    Float64[v2 ? 1.0 : 0.0, t3amp ? 1.0 : 0.0, t3phi ? 1.0 : 0.0]

"""
What a finished run produced, and what it cost.

`chi2r_start` is kept beside `chi2r` because a reduced χ² on its own says very little: what
tells you the run did something is the two together. `flux` is the image sum, which is worth
reporting precisely because nothing constrains it — V² and closure phase are both invariant
under a global scaling, so a reconstruction from them alone is free to drift away from unit
flux, and a user who does not know that reads the number as a bug.
"""
struct ImagingResult
    image       :: Matrix{Float64}
    chi2r       :: Float64
    chi2r_start :: Float64
    flux        :: Float64
    maxiter     :: Int
    seconds     :: Float64
    setup       :: ImagingSetup
    weights     :: Vector{Float64}
end

function Base.show(io::IO, r::ImagingResult)
    print(io, "ImagingResult: ", r.setup.nx, "×", r.setup.nx, " at ",
          round(r.setup.pixsize; digits = 4), " mas  ",
          "χ²ᵣ ", round(r.chi2r_start; digits = 2), " → ", round(r.chi2r; digits = 2),
          "  flux ", round(r.flux; digits = 4),
          "  (", r.maxiter, " iter, ", round(r.seconds; digits = 2), " s)")
end

"""
    reconstruct_image(data, setup; weights, maxiter, regularizers, verb) -> ImagingResult

Run one reconstruction and report what it did.

`data` is the full `Array{OIdata,2}` that `readoifits` returns, not one element of it: the
Fourier plans are built per wavelength/time bin and `setup_ft` wants the same grid.

`regularizers` defaults to none, which leaves positivity as the only thing shaping the image.
That is a real reconstruction, not a placeholder — it is what an unregularised maximum
likelihood image is — and it is the configuration in which a broken criterion, plan or
weighting shows up as a χ² that does not move.
"""
function reconstruct_image(data::AbstractArray, setup::ImagingSetup;
                           weights      = imaging_weights(),
                           maxiter      ::Integer = 200,
                           regularizers = [],
                           verb         ::Bool = false)
    all(iszero, weights) &&
        error("Every observable is switched off; there is nothing to fit")

    ft = setup_ft(data, setup.nx, setup.pixsize; mode = String(setup.mode))
    x0 = start_image(setup, ft)

    chi2r_start = image_to_chi2(x0, ft, data; weights, verb = false)
    t = @elapsed x = reconstruct(x0, data, ft; weights, regularizers, maxiter, verb)
    chi2r = image_to_chi2(x, ft, data; weights, verb = false)

    img = Float64.(x)
    return ImagingResult(img, chi2r, chi2r_start, sum(img), Int(maxiter), t,
                         setup, Float64.(weights))
end
