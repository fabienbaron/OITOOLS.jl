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
end

ImagingSetup(; nx = 64, pixsize = 0.25, mode = :nfft,
               startkind = :gaussian, startfwhm = AUTO_FWHM, startpath = "",
               startseed = 1) =
    ImagingSetup(nx, Float64(pixsize), Symbol(mode), Symbol(startkind),
                 Float64(startfwhm), String(startpath), Int(startseed))

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
    t = @elapsed x = reconstruct(x0, data, ft; weights, regularizers, maxiter, verb)
    chi2 = image_to_chi2(x, ft, data; weights, verb = false)

    img = Float64.(x)
    bd  = chi2_breakdown(x, ft, data, Float64.(weights))
    # Only the observables that were actually in the criterion: dividing by points the fit
    # never saw would flatter it.
    ndof = sum(b -> b.used ? b.n : 0, bd; init = 0)
    return ImagingResult(img, Float64(chi2), Float64(chi2_start), ndof, sum(img),
                         Int(maxiter), t, setup, Float64.(weights), bd)
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
