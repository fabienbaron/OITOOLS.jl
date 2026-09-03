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

using Random

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
    startkind :: Symbol        # :dirac, :gaussian, :random or :fits
    startfwhm :: Float64       # mas; NEGATIVE means "choose it from the field of view"
    startpath :: String
    startseed :: Int           # for :random, so a run can be repeated
    engine    :: Symbol        # which reconstructor runs; see IMAGING_ENGINES
end

ImagingSetup(; nx = 64, pixsize = 0.25, mode = :nfft,
               startkind = :gaussian, startfwhm = AUTO_FWHM, startpath = "",
               startseed = 1, engine = :vmlmb) =
    ImagingSetup(nx, Float64(pixsize), Symbol(mode), Symbol(startkind),
                 Float64(startfwhm), String(startpath), Int(startseed), Symbol(engine))

"""
The reconstructors the Image perspective can run, and the call each one makes.

Every engine here is genuinely dispatched to. A name that appears in the panel but falls
through to `reconstruct` would be worse than an absent one: the run succeeds, the image looks
plausible, and the label on it is wrong.

`:tempering` and `:vi` are deliberately absent. Both need a package that is not a dependency --
Pigeons and OIVI -- so the panel blocks them with that reason rather than offering them here.
"""
const IMAGING_ENGINES = Dict{Symbol,String}(
    :vmlmb          => "reconstruct",
    :bsmem          => "reconstruct_bsmem",
    :bsdmm          => "reconstruct_bsdmm",
    :sparco         => "reconstruct_sparco",
    :squeeze        => "reconstruct_squeeze",
    :squeeze_sparco => "reconstruct_squeeze (model = SqueezeSparco)")

"""
The engines that reconstruct a WAVELENGTH CUBE rather than one grey image.

`reconstruct` takes a 4-D `(nx, nx, nwav, nepoch)` image with `transspectral_regularizers`;
`reconstruct_bsdmm` takes the same shape with `mu_group` coupling the channels. Every other
engine here takes one `OIdata` and one plan cell.

A set rather than a check on the function: whether an engine can do this is a property of the
engine, and `_require_mono` needs to answer it for engines it is not otherwise dispatching.
"""
const POLYCHROMATIC_ENGINES = Set([:vmlmb, :bsdmm])

"""
Sentinel for "work the starting width out from the field of view".

Negative, not zero. Zero is a width — a degenerate one, but a reader has no way to tell it from
a value someone meant, and a control showing `0` invites being left alone as though it said
something. A negative FWHM cannot be mistaken for a width anyone intended.
"""
const AUTO_FWHM = -1.0

"""
    imaging_defaults(data; nx = 64) -> ImagingSetup

Geometry the data itself suggests: `auto_pixsize` for the pixel, `nx` left to the caller.

The field of view that falls out is `nx * pixsize`; whether it covers the source is a
judgement the panel shows rather than makes.
"""
function imaging_defaults(data; nx::Integer = 64)
    # The FINEST pixel any bin needs, not bin 1's. `auto_pixsize` takes a single `OIdata` and
    # sizes the pixel from that bin's longest baseline; the shortest wavelength resolves the
    # most, so a cube reconstructed at bin 1's pixel would be under-sampled wherever the
    # baselines reach further. Taking the minimum costs nothing on a single bin.
    ds = data isa AbstractArray ? vec(data) : [data]
    ImagingSetup(; nx = Int(nx), pixsize = minimum(auto_pixsize(d) for d in ds))
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
        # `gaussian2d(nx, nx, nx/6)`, which is what the shipped demos start from — a Gaussian
        # of sigma nx/6 PIXELS, unit total flux. A positive `startfwhm` overrides it in mas.
        setup.startfwhm > 0 ?
            reshape(gaussian_prior(nx, setup.pixsize; fwhm_mas = setup.startfwhm), nx, nx) :
            Float64.(gaussian2d(nx, nx, nx / 6))
    elseif setup.startkind === :random
        # Noise, for asking whether the answer depends on where you started — the check a
        # smooth start cannot make.
        #
        # Uniform, not Gaussian. The optimiser is bounded below by zero, so the start has to be
        # non-negative; `rand` already is, whereas a Gaussian would need folding, and folding
        # it piles half the pixels up against the bound the optimiser is about to work against.
        #
        # Seeded, always. An unseeded start makes a reconstruction that cannot be repeated and
        # an exported script that does not reproduce its own figure.
        m = rand(Random.Xoshiro(setup.startseed), Float64, nx, nx)
        m ./ sum(m)
    elseif setup.startkind === :fits
        isfile(setup.startpath) || error("No starting image at '$(setup.startpath)'")
        m = Float64.(readfits(setup.startpath))
        size(m) == (nx, nx) ||
            error("The starting image is $(size(m,1))×$(size(m,2)) but nx is $nx")
        s = sum(m)
        s > 0 || error("The starting image sums to $s; it must carry positive flux")
        m ./ s
    else
        error("Unknown starting image '$(setup.startkind)'; " *
              "expected :dirac, :gaussian, :random or :fits")
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
    _cube_for_chi2(image, data) -> image in the shape `image_to_chi2` needs

A 2-D image against a single bin stays 2-D. Against several bins it is REPLICATED to
`(nx, nx, nwav, nepoch)`, which is what a grey image means when it is scored against
polychromatic data: the same brightness distribution in every channel. An image that is already
a cube is passed through.

Replicated, not reshaped to `(nx, nx, 1, 1)`. The 4-D criterion loops `w in 1:nwav` and indexes
`x[:, :, w, t]`, so a one-channel cube against several bins raises a `BoundsError` — the same
failure this is here to prevent, moved one function along.
"""
function _cube_for_chi2(image, data)
    ndims(image) == 2 || return image
    (data isa AbstractArray && length(data) > 1) || return image
    nwav, nepoch = size(data)
    return repeat(image, 1, 1, nwav, nepoch)
end

"""
    chi2_breakdown(image, ft, data, weights; chi2_of = nothing) -> Vector{NamedTuple}

Reduced χ² per observable: `(; name, chi2, n, chi2r, used)`.

Computed by asking the criterion itself, one observable at a time with a one-hot weight vector,
rather than by re-deriving the residuals here. The breakdown does exist inside `_chi2_f`, but
only as something it prints when `verb = true` — it returns the weighted sum alone — and a
panel that recomputed it by hand would be free to drift away from what is actually being
minimised.

Always from the FINAL image, never from whatever an engine printed on its way past: BSMEM's own
number is 505 where the image scores 884 on the same run, because it minimises a different
statistic. The image is what the panel displays, so the image is what the panel scores.

`chi2_of` replaces the criterion for a model the image alone does not describe. SPARCO's is the
image PLUS a parametric component, so `image_to_chi2` of its image is meaningless — measured at
χ²ᵣ ≈ 2.7e7 against the engine's own 1.2 — and `chi2_sparco_flat_f` is what scores it.

`used` records whether the observable was in the fit at all. An observable the user unticked
still has a χ² worth seeing: it says how well the reconstruction predicts data it never saw.
"""
function chi2_breakdown(image, ft, data, weights = [1.0, 1.0, 1.0]; chi2_of = nothing)
    ds = data isa AbstractArray ? vec(data) : [data]
    # `n` below counts points across EVERY bin, so the chi2 has to as well. A 2-D image against
    # a multi-bin `data` silently reduces to `ft[1], data[1]` inside `image_to_chi2`
    # (oichi2.jl), which would score one channel and divide it by the point count of all of
    # them -- a reduced chi2 wrong by roughly the number of bins, in the direction that
    # flatters the fit. Reshaping to the 4-D form picks the method that loops the cells.
    img = _cube_for_chi2(image, data)
    counts = (v2 = sum(d -> d.nv2, ds; init = 0),
              t3amp = sum(d -> d.nt3amp, ds; init = 0),
              t3phi = sum(d -> d.nt3phi, ds; init = 0))
    labels = ("V²", "T3amp", "T3φ")
    keys_  = (:v2, :t3amp, :t3phi)
    out = NamedTuple[]
    for i in 1:3
        n = counts[keys_[i]]
        n == 0 && continue
        onehot = [j == i ? 1.0 : 0.0 for j in 1:3]
        c2 = try
            chi2_of === nothing ?
                Float64(image_to_chi2(img, ft, data; weights = onehot, verb = false)) :
                Float64(chi2_of(onehot))
        catch err
            # Named, not swallowed. A bare `catch` here turned a precision or shape mismatch
            # into a NaN in the panel with nothing anywhere saying why -- and NaN in a χ²
            # column reads as a failed reconstruction rather than as a failed measurement of
            # one.
            @warn "could not compute the $(labels[i]) χ²" exception = (err, catch_backtrace())
            NaN
        end
        push!(out, (; name = labels[i], chi2 = c2, n = n, chi2r = c2 / n,
                      used = weights[i] > 0))
    end
    return out
end

"""
What a finished run produced, and what it cost.

`image_to_chi2` returns the RAW χ², not a reduced one — the sum over every observable, with no
division. Storing it under a name ending in `r` would be a quiet lie of the worst kind: on this
dataset the raw number is 315 where the reduced one is 0.69, and 315 reads as a fit that failed
completely. `chi2r` is a property here, computed against the points actually fitted.

`ndof` counts data points and is NOT reduced by any parameter count — an image has as many
parameters as pixels, so the usual correction is meaningless and pretending otherwise would be
worse than omitting it.

`chi2_start` is kept beside `chi2` because a χ² on its own says very little: what tells you the
run did something is the two together. `flux` is the image sum, which is worth
reporting precisely because nothing constrains it — V² and closure phase are both invariant
under a global scaling, so a reconstruction from them alone is free to drift away from unit
flux, and a user who does not know that reads the number as a bug.
"""
struct ImagingResult
    image       :: Matrix{Float64}
    chi2        :: Float64          # raw, as image_to_chi2 returns it
    chi2_start  :: Float64
    ndof        :: Int              # valid points in the observables actually fitted
    flux        :: Float64
    maxiter     :: Int
    seconds     :: Float64
    setup       :: ImagingSetup
    weights     :: Vector{Float64}
    breakdown   :: Vector{NamedTuple}   # reduced chi2 per observable
    # Whatever the engine returned beyond the image: SQUEEZE's diagnostics NamedTuple, SPARCO's
    # fitted parametric values. `nothing` for the engines that hand back an image and no more.
    extra       :: Any
    # The ensemble, for the engines that return a distribution rather than a point estimate:
    # `nothing`, or `(; mean, sigma, samples, source)`. Shaped to match `ImageEntry.posterior`
    # so a result can become a session entry unchanged.
    ensemble    :: Any
end

# The engine is on `setup`, so a result always says which reconstructor made it.
ImagingResult(image, chi2, chi2_start, ndof, flux, maxiter, seconds, setup, weights, breakdown) =
    ImagingResult(image, chi2, chi2_start, ndof, flux, maxiter, seconds, setup, weights,
                  breakdown, nothing, nothing)

"""
    result_ensemble(extra) -> nothing or (; mean, sigma, samples, source)

The ensemble an engine returned, or `nothing` for the engines that return one image.

Only the sampling engines produce one, and what they produce is worth naming precisely:
`reconstruct_squeeze` runs `nchains` independent chains and hands back
`diagnostics.images`, **one posterior mean per chain** (each already averaged over that chain's
post-burn-in samples). So the spread across them is CHAIN-TO-CHAIN disagreement, not a within-
chain posterior width, and the panel says so rather than calling it σ and leaving it at that.

The reconstruction itself returns the BEST chain, not this mean — `reconstruct_squeeze` ends
with `results[best].image`. The two are different images and the panel offers both.

`sigma` is `nothing` for a single chain, where there is no spread to report.
"""
function result_ensemble(extra)
    extra isa NamedTuple && haskey(extra, :images) || return nothing
    imgs = extra.images
    (imgs isa AbstractVector && !isempty(imgs)) || return nothing
    all(m -> m isa AbstractMatrix, imgs) || return nothing
    samples = [Float64.(m) for m in imgs]
    allequal(size.(samples)) || return nothing
    n = length(samples)
    mn = reduce(+, samples) ./ n
    sg = n > 1 ? sqrt.(reduce(+, ((s .- mn).^2 for s in samples)) ./ (n - 1)) : nothing
    return (; mean = mn, sigma = sg, samples,
              source = n == 1 ? "1 chain" : "$(n) chains")
end

"Reduced χ², against the points actually fitted. `NaN` when nothing was."
chi2r(r::ImagingResult) = r.ndof > 0 ? r.chi2 / r.ndof : NaN

"Reduced χ² of the starting image, on the same footing."
chi2r_start(r::ImagingResult) = r.ndof > 0 ? r.chi2_start / r.ndof : NaN

function Base.show(io::IO, r::ImagingResult)
    print(io, "ImagingResult: ", r.setup.nx, "×", r.setup.nx, " at ",
          round(r.setup.pixsize; digits = 4), " mas  ",
          "χ²ᵣ ", round(chi2r_start(r); digits = 3), " → ", round(chi2r(r); digits = 3),
          "  flux ", round(r.flux; digits = 4),
          "  (", r.maxiter, " iter, ", round(r.seconds; digits = 2), " s)")
end

"""
    ft_summary(ft, setup) -> String

One line describing a Fourier plan set: mode, geometry, channels and uv points.

`ft_info` prints and returns nothing, which suits the REPL and cannot fill a label. This is the
same facts on one line, and it is worth showing: the plan is the thing that decides what the
reconstruction can represent, and it is silently rebuilt whenever the geometry changes.
"""
function ft_summary(ft, setup::ImagingSetup)
    nwav, nepoch = size(ft)
    nuv = try
        c = ft[1]
        c isa NFFTCell ? size(c.uv.k, 2) : size(c, 1)
    catch
        0
    end
    chans = nwav * nepoch
    mode = uppercase(String(setup.mode))
    parts = ["$mode  $(setup.nx)×$(setup.nx)  $(round(setup.pixsize; digits = 4)) mas",
             "FOV $(round(fov(setup); digits = 2)) mas",
             chans == 1 ? "1 channel" : "$chans channels",
             "$nuv uv points"]
    return join(parts, "  ·  ")
end

"""
    ensure_ft!(cache, data, setup) -> (ft, summary, rebuilt)

The plan for this geometry, building it only when the geometry actually changed.

`cache` is whatever the last call returned, or `nothing`. Rebuilding an NFFT plan is the
expensive part of setting a reconstruction up — it is why the geometry controls can be moved
freely and only the final one costs anything.
"""
function ensure_ft!(cache, data, setup::ImagingSetup)
    key = (objectid(data), setup.nx, setup.pixsize, setup.mode)
    if cache !== nothing && cache.key == key
        return (cache.ft, cache.summary, false)
    end
    ft = setup_ft(data, setup.nx, setup.pixsize; mode = String(setup.mode))
    return (ft, ft_summary(ft, setup), true)
end

"""
    reconstruct_image(data, setup; weights, maxiter, regularizers, ft, verb) -> ImagingResult

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
                           ft           = nothing,
                           x_start      = nothing,
                           options      ::AbstractDict = Dict{String,String}(),
                           verb         ::Bool = false)
    regularizers = parse_regularizers(regularizers)
    all(iszero, weights) &&
        error("Every observable is switched off; there is nothing to fit")

    # A plan built by the panel when the geometry changed is reused; nothing is cached here,
    # since this function has no business owning state that outlives one call.
    ft = ft === nothing ?
         setup_ft(data, setup.nx, setup.pixsize; mode = String(setup.mode)) : ft
    # A supplied image continues from where a previous run stopped, which is what makes
    # "run from previous" more than a restart: VMLMB begins at that point rather than at a
    # generic start, so a run can be extended instead of repeated.
    x0 = if x_start === nothing
        start_image(setup, ft)
    else
        size(x_start) == (setup.nx, setup.nx) ||
            error("The previous image is $(size(x_start,1))×$(size(x_start,2)) but nx is " *
                  "$(setup.nx); change nx back or start fresh")
        to_ft_precision(Float64.(x_start), ft)
    end

    chi2_start = image_to_chi2(x0, ft, data; weights, verb = false)
    t = @elapsed (xraw, extra, own_chi2) = run_engine(setup.engine, x0, data, ft;
                                                     weights, regularizers, maxiter, verb, options)
    # The engines hand back different shapes and precisions -- a Float32 matrix from VMLMB, a
    # Float64 one from SQUEEZE and BSDMM, an (nx, nx, 1) slice from BSMEM. The shared χ² path
    # takes one shape, so normalise once here; comparing two engines on the same number is the
    # reason they share a panel at all.
    x = to_ft_precision(Float64.(xraw), ft)

    # `reconstruct_sparco` hands back the fitted `model` and `params` alongside the image, and
    # with those the SPARCO criterion scores the same thing the engine scored -- so the panel's
    # split adds up to the engine's own χ² instead of being dropped.
    sparco = extra isa NamedTuple && haskey(extra, :model) ? extra : nothing
    bd = if sparco === nothing
        chi2_breakdown(x, ft, data, Float64.(weights))
    else
        # `ft[1, 1]` unchanged: `NFFTCell <: AbstractVector{<:NFFTPlan}`, so the cell already IS
        # the plan vector the criterion wants -- the same value `reconstruct_sparco` passes it.
        ds = data isa AbstractArray ? first(data) : data
        chi2_breakdown(x, ft, data, Float64.(weights);
                       chi2_of = w -> chi2_sparco_flat_f(x, sparco.params, sparco.model,
                                                         ft[1, 1], ds;
                                                         w_name = "W", weights = w, verb = false))
    end
    # Only the observables that were actually in the criterion: dividing by points the fit
    # never saw would flatter it. The point counts do not depend on the image, so this is right
    # even for the engines whose χ² is not computed from the image.
    ndof = sum(b -> b.used ? b.n : 0, bd; init = 0)

    chi2 = own_chi2 === nothing ? Float64(image_to_chi2(x, ft, data; weights, verb = false)) :
                                  Float64(own_chi2)
    # What is left are the engines that score something the image alone does not predict and
    # give no way to split it -- SQUEEZE with a SPARCO component, whose parametric part lives
    # inside the sampler. They report their own χ² and no breakdown, rather than an image-only
    # split that would not sum to it.
    (own_chi2 === nothing || sparco !== nothing) || (bd = NamedTuple[])

    img = Float64.(x)
    return ImagingResult(img, chi2, Float64(chi2_start), ndof, sum(img),
                         Int(maxiter), t, setup, Float64.(weights), bd, extra,
                         result_ensemble(extra))
end

# ─────────────────────────────────────────────────────────────────────────────
# The engines
# ─────────────────────────────────────────────────────────────────────────────
#
# Each branch below makes ONE call to a real OITOOLS reconstructor and returns `(image, extra)`.
# The χ², the per-observable breakdown and the point count are then computed the same way for
# all of them, by `reconstruct_image` above -- so two engines can be compared on the same
# number, which is the whole point of having them in one panel.
#
# Options arrive as a `Dict{String,String}` from the panel rather than as a keyword each: the
# engines take wildly different settings (BSMEM a four-integer method vector, BSDMM three μ
# weights and two mode symbols, SQUEEZE move probabilities), and threading every one of them
# through a shared signature would put BSDMM's knobs in BSMEM's call.

_optstr(o, k, d::AbstractString = "") = (v = get(o, k, ""); isempty(v) ? d : v)
_optbool(o, k, d::Bool) = (v = get(o, k, ""); isempty(v) ? d : lowercase(v) in ("1", "true", "yes"))
_optint(o, k, d::Integer) = (v = get(o, k, ""); isempty(v) ? Int(d) : something(tryparse(Int, v), Int(d)))
_optreal(o, k, d::Real) = (v = get(o, k, ""); isempty(v) ? Float64(d) : something(tryparse(Float64, v), Float64(d)))

"""
The dataset's single bin, checked for a spread of wavelengths.

The test is on `uv_lam`, NOT on the number of wavelength BINS. OITOOLS keeps a chromatic
dataset as one `OIdata` whose points each carry their own λ, and that is the form SPARCO wants:
it weights each point by its wavelength, so a file read into one bin spanning 1.5-1.7 μm is
exactly right and `size(data, 1) == 1` says nothing about it. Splitting the same file into 35
bins instead makes `reconstruct_sparco` fail on mismatched per-channel point counts.
"""
function _require_chromatic(data, engine)
    d = _require_mono(data, engine)
    length(unique(d.uv_lam)) > 1 || error(
        "$(get(IMAGING_ENGINES, engine, engine)) separates its components by their spectral " *
        "indices, so it needs a spread of wavelengths; every point in this dataset is at " *
        "$(round(Float64(first(d.uv_lam)) * 1e6, digits = 4)) μm.")
    return d
end

"Every engine but VMLMB is monochromatic here; say so rather than failing inside the engine."
function _require_mono(data, engine)
    size(data) == (1, 1) || error(
        "$(get(IMAGING_ENGINES, engine, engine)) runs on one wavelength bin; this dataset has " *
        "$(size(data, 1))×$(size(data, 2)). Use an engine that reconstructs a cube — " *
        join(sort(string.(collect(POLYCHROMATIC_ENGINES))), " or ") *
        " — or read the file with a single bin.")
    return data[1, 1]
end

"""
    _check_bins(data, engine)

Refuse a multi-bin dataset for an engine that cannot reconstruct a cube.

Every engine needs this, INCLUDING the ones that can: until the cube is threaded all the way
through, `run_engine` still hands them a 2-D `x0`. VMLMB is the case that matters, because it
is the one engine `_require_mono` never guarded: `reconstruct(::Matrix, ::Matrix{OIdata}, …)`
reshapes the image to one channel but forwards every bin, and `crit_fg` then indexes
`x4[:,:,w,t]` for `w in 1:nwav` and throws a `BoundsError` from inside OITOOLS with nothing to
say which control caused it.
"""
function _check_bins(data, engine)
    size(data) == (1, 1) && return data
    error("$(get(IMAGING_ENGINES, engine, engine)) cannot yet reconstruct a wavelength cube; " *
          "this dataset has $(size(data, 1)) bins × $(size(data, 2)) epochs. " *
          "Read the file with a single bin.")
end

"""
    run_engine(engine, x0, data, ft; weights, regularizers, maxiter, verb, options)
      -> (image, extra, chi2)

Dispatch one reconstruction.

`extra` carries whatever the engine returned beyond the image — SQUEEZE's diagnostics, SPARCO's
fitted parameters — and `chi2` is the engine's own raw χ², or `nothing` to have it computed
from the image. Only the engines with a parametric component need to supply it: for those the
image is half the model, and χ² of the image alone is not a number about anything.
"""
function run_engine(engine::Symbol, x0, data, ft;
                    weights, regularizers, maxiter::Integer, verb::Bool,
                    options::AbstractDict = Dict{String,String}())
    o = options

    if engine === :vmlmb
        _check_bins(data, engine)
        return (reconstruct(x0, data, ft; weights, regularizers, maxiter, verb), nothing, nothing)

    elseif engine === :bsmem
        d = _require_mono(data, engine)
        # BSMEM ignores the spec list except for a `["mem", prior]` entry, and takes its
        # entropy mode as the four-integer `method` vector instead. Passing the panel's
        # regularisers here would silently do nothing, so only the prior is forwarded.
        regs = Any[]
        prior = _optstr(o, "prior")
        isempty(prior) || push!(regs, ["mem", prior_image(prior, size(x0, 1))])
        method = [_optint(o, "method1", 4), _optint(o, "method2", 1),
                  _optint(o, "method3", 1), _optint(o, "method4", 2)]
        # `ft[1, 1]`, not `ft`: BSMEM reads the geometry off `ft[1].N`, so it wants the CELL
        # for this bin, whose first element is the NFFT plan. Handing it the whole `OIft` makes
        # `ft[1]` the cell itself, which has no `N`.
        # `history` is BSMEM's own per-iteration record — entropy, χ², α, ω. The engines have no
        # callback API, but the ones that keep a history will fill a vector you pass in, and
        # that is structured data rather than scraped text.
        hist = NamedTuple[]
        img = reconstruct_bsmem(x0, d, ft[1, 1]; regularizers = regs, method,
                                maxiter = _optint(o, "maxiter", maxiter), verbose = verb,
                                history = hist)
        # (nx, nx, nwav) even for one channel.
        return (img[:, :, 1], (; history = hist), nothing)

    elseif engine === :bsdmm
        _check_bins(data, engine)
        # Not string specs at all: ADMM takes μ weights and two mode symbols.
        mu_reg, mu_cen = _optreal(o, "mu_reg", 0.0), _optreal(o, "mu_cen", 0.0)
        # Centering specifically, not just any block. ADMM splits the problem across proximal
        # blocks and a weight of zero creates none, so all-zero weights leave nothing to solve
        # -- but beyond that, the centering block is what pins the image against the translation
        # degeneracy of the bispectrum. Without it the image is free to wander the field, and
        # the run is unreproducible even when it converges.
        mu_cen > 0 || error(
            "BSDMM needs a non-zero mu_cen: the centering block is what holds the image " *
            "against the translation degeneracy of the bispectrum, and with every weight at " *
            "zero there are no ADMM blocks to solve at all.")
        img = reconstruct_bsdmm(x0, data, ft;
                                mu_reg, mu_cen,
                                reg_type = Symbol(_optstr(o, "reg_type", "tv")),
                                maxit    = _optint(o, "maxiter", maxiter),
                                verb, history = true)
        # With `history = true` ADMM returns `(image, trace)`; the trace carries the primal and
        # dual residuals and ρ, which are what actually say whether it is converging.
        img isa Tuple && return (img[1], (; history = img[2]), nothing)
        return (img, nothing, nothing)

    elseif engine === :sparco
        d = _require_chromatic(data, engine)
        r = reconstruct_sparco(x0, d, ft[1, 1];
                               lambda_ref = _optreal(o, "lambda0", 1.65e-6),
                               star_flux  = _optreal(o, "f_star", 0.5),
                               bg_flux    = _optreal(o, "f_bg", 0.0),
                               d_env      = _optreal(o, "env_indx", 0.0),
                               # `star_di`, not `ud`: this is the star's SPECTRAL INDEX in
                               # (lambda/lambda0)^(-star,di), legal over -20..20 and defaulting
                               # to 4. The annealing SPARCO's `ud` is a uniform-disc diameter in
                               # mas -- a different quantity, which is why they are different
                               # option keys rather than one field feeding both.
                               star_di    = _optreal(o, "star_di", 4.0),
                               weights, regularizers, maxiter = Int(maxiter),
                               rounds = _optint(o, "rounds", 3), verb)
        return (r.image, (; r.params, r.param_names, r.free_params, r.model), r.chi2)

    elseif engine === :squeeze || engine === :squeeze_sparco
        d = _require_mono(data, engine)
        # `SqueezeSparco` carries spectral indices for the star, the environment and the
        # background, so it separates them across the band exactly as `reconstruct_sparco`
        # does -- and is just as degenerate at a single wavelength.
        engine === :squeeze_sparco && _require_chromatic(data, engine)
        model = engine === :squeeze_sparco ?
            SqueezeSparco(; f_star   = _optreal(o, "f_star", 0.5),
                            ud       = _optreal(o, "ud", 0.0),
                            env_indx = _optreal(o, "env_indx", 0.0),
                            lambda0  = _optreal(o, "lambda0", 1.65e-6),
                            f_bg     = _optreal(o, "f_bg", 0.0),
                            bg_indx  = _optreal(o, "bg_indx", 0.0),
                            free     = _sparco_free(o)) : nothing
        prior = _optstr(o, "prior")
        img, diag = reconstruct_squeeze(x0, d, ft[1, 1];      # the bin's cell, as BSMEM takes
            weights, regularizers,
            nelements  = _optint(o, "nelements", 0),
            niter      = _optint(o, "niter", maxiter),
            nchains    = _optint(o, "nchains", 1),
            tmin       = _optreal(o, "tmin", 1.0),
            f_anywhere = _optreal(o, "f_anywhere", 0.05),
            f_copycat  = _optreal(o, "f_copycat", 0.1),
            # No centering under SPARCO, by default and for the same reason the panel greys
            # the controls: the parametric component sits at the centre by construction, so the
            # translation degeneracy the centering term breaks is already broken, and pulling
            # the image to the centre on top of that fights the star for the same position.
            cent_mult  = _optreal(o, "cent_mult", engine === :squeeze_sparco ? 0.0 : 1.0),
            auto_centering = _optbool(o, "auto_centering", engine !== :squeeze_sparco),
            prior_image = isempty(prior) ? nothing : prior_image(prior, size(x0, 1)),
            model,
            # The in-package monitor draws with matplotlib, and matplotlib called off the main
            # thread segfaults inside PythonCall -- measured. The GUI draws the result itself.
            monitor = 0, print_every = 0, verb,
            seed = _optint(o, "seed", 12345))
        # The plain sampler's image IS the whole model, so the shared χ² applies. With a SPARCO
        # component it is not, and the sampler's own reduced χ² is the honest number.
        own = engine === :squeeze_sparco ? _squeeze_chi2(diag) : nothing
        return (img, diag, own)
    end

    error("Unknown imaging engine $(repr(engine)); known ones are " *
          join(sort(string.(keys(IMAGING_ENGINES))), ", "))
end

"""
    prior_image(path, nx) -> Matrix{Float64}

A prior or mask image from a FITS file, checked against the reconstruction geometry.

BSMEM takes it as the MEM prior; SQUEEZE takes it as a per-pixel prior probability, where a
zero makes the pixel unreachable — the same file therefore doubles as a hard mask there. Both
need it on the image grid, so a size mismatch is an error rather than something to interpolate.
"""
function prior_image(path::AbstractString, nx::Integer)
    isfile(path) || error("No prior image at '$path'")
    m = Float64.(readfits(String(path)))
    size(m) == (nx, nx) ||
        error("The prior image is $(size(m,1))×$(size(m,2)) but nx is $nx")
    all(isfinite, m) || error("The prior image at '$path' has non-finite pixels")
    return m
end

"""
The sampler's own raw χ², from its diagnostics, or `nothing` when they do not carry one.

Per-chain quantities come back as vectors, so the best chain is the one to quote: the ensemble
mean includes chains that are still hot, and reporting it would understate a run that found the
answer on one chain.
"""
function _squeeze_chi2(diag)
    hasproperty(diag, :ndf) || return nothing
    r = if hasproperty(diag, :chi2r_last) && diag.chi2r_last isa AbstractVector &&
           !isempty(diag.chi2r_last)
        best = hasproperty(diag, :best_chain) ? diag.best_chain : argmin(diag.chi2r_last)
        Float64(diag.chi2r_last[clamp(Int(best), 1, length(diag.chi2r_last))])
    elseif hasproperty(diag, :chi2r_mean)
        m = diag.chi2r_mean
        Float64(m isa AbstractVector ? (isempty(m) ? NaN : minimum(m)) : m)
    else
        return nothing
    end
    return isfinite(r) ? r * Float64(diag.ndf) : nothing
end

"Which SPARCO parameters SQUEEZE samples, from the panel's tick boxes."
function _sparco_free(o)
    free = Symbol[]
    for (key, name) in (("free_f_star", :f_star), ("free_ud", :ud),
                        ("free_env_indx", :env_indx), ("free_f_bg", :f_bg),
                        ("free_bg_indx", :bg_indx))
        _optbool(o, key, key == "free_f_star") && push!(free, name)
    end
    # SqueezeSparco with nothing free is a fixed parametric component, which is legitimate but
    # is almost never what a tick-box panel meant to say.
    isempty(free) ? (:f_star,) : Tuple(free)
end


"""
    observable_availability(data) -> NamedTuple

Which observables the dataset actually contains, as six booleans.

Counted, not inferred from the presence of a table: `nv2` and friends are *valid-point*
counts, and a table whose points were all sanitised away (they get `err = Inf` to keep index
alignment) contains nothing to fit even though it exists. A panel that enabled a tick box on
the strength of the table would offer an observable with no data behind it.

`diffvis` additionally needs more than one wavelength — a differential phase is measured
against the other channels, so a monochromatic file cannot carry one whatever its OI_VIS says.

This answers only "does the file have it". Whether an *engine* can use it is a separate
question with a separate answer, and the two must not be merged: greying a box out for two
different reasons without distinguishing them is how a user concludes the file is broken.
"""
function observable_availability(data)
    ds = data isa AbstractArray ? vec(data) : [data]
    tot(f) = sum(d -> Int(getfield(d, f)), ds; init = 0)
    nwav = data isa AbstractArray ? size(data, 1) : 1
    return (; v2      = tot(:nv2)     > 0,
              t3amp   = tot(:nt3amp)  > 0,
              t3phi   = tot(:nt3phi)  > 0,
              cvis    = tot(:nvisamp) > 0 || tot(:nvisphi) > 0,
              flux    = tot(:nflux)   > 0,
              diffvis = nwav > 1 && tot(:nvisphi) > 0)
end

"The six flags as `name=0|1` pairs, which is the shape QML reads them back in."
observable_flags_string(data) =
    join(("$k=$(Int(v))" for (k, v) in pairs(observable_availability(data))), ",")


"""
    parse_regularizers(spec) -> Vector

Accept regularisers as `reconstruct` wants them, or as the panel has them.

`reconstruct` takes `["name", mu, extra...]` lists. A panel has strings, so a row arrives as
`"name,mu"` or `"name,mu,extra"` and is turned into that list here — one place that knows the
format, rather than the format being rebuilt wherever a run is started.

An empty spec is the positivity-only reconstruction, which is a real one: positivity is VMLMB's
`lower = 0`, not a regulariser, so it is in force either way.
"""
function parse_regularizers(spec)
    spec === nothing && return Any[]
    spec isa AbstractString || return collect(spec)
    out = Any[]
    for row in split(String(spec), ';')
        row = strip(row)
        isempty(row) && continue
        f = [strip(x) for x in split(row, ',')]
        isempty(f[1]) && continue
        entry = Any[String(f[1])]
        for v in f[2:end]
            isempty(v) && continue
            push!(entry, parse(Float64, v))
        end
        length(entry) == 1 &&
            error("Regulariser '$(f[1])' has no weight; give it as \"name,mu\"")
        push!(out, entry)
    end
    return out
end
