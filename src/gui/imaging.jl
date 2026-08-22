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
    chi2_breakdown(image, ft, data, weights) -> Vector{NamedTuple}

Reduced χ² per observable: `(; name, chi2, n, chi2r, used)`.

Computed by asking the criterion itself, one observable at a time with a one-hot weight vector,
rather than by re-deriving the residuals here. The breakdown does exist inside `_chi2_f`, but
only as something it prints when `verb = true` — it returns the weighted sum alone — and a
panel that recomputed it by hand would be free to drift away from what is actually being
minimised.

`used` records whether the observable was in the fit at all. An observable the user unticked
still has a χ² worth seeing: it says how well the reconstruction predicts data it never saw.
"""
function chi2_breakdown(image, ft, data, weights = [1.0, 1.0, 1.0])
    ds = data isa AbstractArray ? vec(data) : [data]
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
            Float64(image_to_chi2(image, ft, data; weights = onehot, verb = false))
        catch
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
end

# The engine is on `setup`, so a result always says which reconstructor made it.
ImagingResult(image, chi2, chi2_start, ndof, flux, maxiter, seconds, setup, weights, breakdown) =
    ImagingResult(image, chi2, chi2_start, ndof, flux, maxiter, seconds, setup, weights,
                  breakdown, nothing)

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

    bd  = chi2_breakdown(x, ft, data, Float64.(weights))
    # Only the observables that were actually in the criterion: dividing by points the fit
    # never saw would flatter it. The point counts do not depend on the image, so this is right
    # even for the engines whose χ² is not computed from the image.
    ndof = sum(b -> b.used ? b.n : 0, bd; init = 0)

    # SPARCO's model is the image PLUS a parametric component, so the image alone does not
    # predict the data and `image_to_chi2` of it is meaningless -- measured at χ²r ≈ 2.7e7
    # against the engine's own 1.2. Those engines report their own χ², and their per-observable
    # breakdown is dropped rather than shown as an image-only split that would not sum to it.
    chi2 = own_chi2 === nothing ? Float64(image_to_chi2(x, ft, data; weights, verb = false)) :
                                  Float64(own_chi2)
    own_chi2 === nothing || (bd = NamedTuple[])

    img = Float64.(x)
    return ImagingResult(img, chi2, Float64(chi2_start), ndof, sum(img),
                         Int(maxiter), t, setup, Float64.(weights), bd, extra)
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
        "$(size(data, 1))×$(size(data, 2)). Use VMLMB, or read the file with a single bin.")
    return data[1, 1]
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
        _require_mono(data, engine)
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
                               star_di    = _optreal(o, "ud", 4.0),
                               weights, regularizers, maxiter = Int(maxiter),
                               rounds = _optint(o, "rounds", 3), verb)
        return (r.image, (; r.params, r.param_names, r.free_params), r.chi2)

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
            cent_mult  = _optreal(o, "cent_mult", 1.0),
            auto_centering = _optbool(o, "auto_centering", true),
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
