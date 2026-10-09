# An image or a cube as a model component.
#
# The FITS convention is shared with ROTIR, whose `write_sky_cube` (src/sky_export.jl) is the
# reference writer. In short — the full statement is in `read_model_image`'s docstring:
#
#   axis 1 East, axis 2 North, phase centre at pixel nx÷2+1 (CRPIX = nx÷2+1, CRVAL = 0, the
#   primary's geometric centre), CDELT positive, CUNIT honoured; axis 3 WAVE in metres, axis 4
#   TIME in days (MJD, UTC); CHANNELS (EFF_WAVE, optional EFF_BAND) and EPOCHS (MJD) tables are
#   the authority over a linear WCS; NORMALIZ = PLANE | RELATIVE.

"""
    ModelImage

An image or cube read for use as a model component, in OITOOLS' layout.

| field | meaning |
|---|---|
| `image` | `(nx, ny, nλ, nt)`, dimension 1 East, dimension 2 North, `Float64` |
| `pixsize` | mas per pixel |
| `wavelengths` | length `nλ`, metres; `NaN` when the file does not say (a grey image) |
| `bandwidths` | length `nλ`, metres; `NaN` when the file does not say |
| `mjds` | length `nt`, MJD; `NaN` when the file does not say |
| `normalization` | `:plane` (each plane sums to 1) or `:relative` (the cube's own SED and light curve kept) |
| `offset` | `(east, north)` in mas of the file's reference pixel from the phase-centre pixel `nx÷2+1`; zero in the shared convention |
| `path` | the file it came from |
"""
struct ModelImage
    image         :: Array{Float64,4}
    pixsize       :: Float64
    wavelengths   :: Vector{Float64}
    bandwidths    :: Vector{Float64}
    mjds          :: Vector{Float64}
    normalization :: Symbol
    offset        :: NTuple{2,Float64}
    path          :: String
end

# Wavelength units a CUNIT3 may name, in metres.
const _WAVE_UNITS = Dict("m" => 1.0, "um" => 1e-6, "micron" => 1e-6, "nm" => 1e-9,
                         "angstrom" => 1e-10, "a" => 1e-10)

"""
    read_model_image(path; pixsize = nothing, normalization = nothing) -> ModelImage

Read a FITS image or cube for use as a model component, in the convention shared with ROTIR.

**Dimensions.** 2-D, 3-D and 4-D are accepted, a missing axis being length 1. A 3-D file is a
wavelength cube or a time series according to `CTYPE3` (`WAVE` or `TIME`); a 3-D file with
more than one plane and no `CTYPE3` is an error naming the keyword, because guessing would
fit the wrong planes to the data.

**Orientation.** The returned array always has dimension 1 East and dimension 2 North, which is
what the NFFT and every model component assume (`test/test_image_orientation.jl`):

- a file whose `CTYPE1` names RA is read by its WCS: `CDELT1 > 0` is already East-first;
  `CDELT1 < 0` — the usual astronomical image, East to the left — is flipped. `CDELT2 < 0` is
  flipped likewise;
- a file without `CTYPE1` is read as stored, East-first, whatever the sign of its `CDELT1`.

**Pixel scale** from `fits_pixsize` (CDELT or CD, with its `CUNIT`), or `pixsize` in mas,
which wins. Neither is an error: an image with no scale cannot be a model.

**Reference pixel.** The convention puts it at `nx÷2+1`, the NFFT's phase centre. A file that
puts it elsewhere is honoured: the difference becomes `offset`, applied as a shift when the
component is evaluated, so positions are always relative to the file's reference pixel.

**Wavelengths and epochs** come from the `CHANNELS` (`EFF_WAVE`, optional `EFF_BAND`) and
`EPOCHS` (`MJD`) tables when present, which are the authority. Otherwise from a linear WCS on
axes 3 and 4 (`CRVAL`, `CDELT`, `CRPIX`, `CUNIT`, `MJDREF`); `WAVEIRR`/`TIMEIRR = T` without
the table is an error.

**Normalisation** from `NORMALIZ`, or `normalization` (`:plane` or `:relative`), which wins:

- `:plane` (`PLANE`, and a file without the key): each plane divided by its own sum, so each
  sums to 1. The sum, never `sum(abs)` or the peak: a band-limited rendering rings slightly
  negative at a limb, and `sum(plane)` is its flux.
- `:relative` (`RELATIVE`): the cube's own SED and light curve kept, the whole cube divided by
  its brightest plane's sum so that plane sums to 1.

A plane summing to zero or less is an error.
"""
function read_model_image(path::AbstractString; pixsize = nothing, normalization = nothing)
    p = String(path)
    isfile(p) || throw(ArgumentError("no such FITS file: $(repr(p))"))
    f = FITS(p)
    raw, hdr, channels, epochs = try
        h1 = f[1]
        x  = read(h1)
        hd = read_header(h1)
        ch = _fits_table(f, "CHANNELS")
        ep = _fits_table(f, "EPOCHS")
        (x, hd, ch, ep)
    finally
        close(f)
    end

    nd = ndims(raw)
    2 <= nd <= 4 || throw(ArgumentError(
        "$(basename(p)): a model image has 2, 3 or 4 axes; this one has $nd"))
    key(k, default = nothing) = haskey(hdr, k) ? hdr[k] : default
    skey(k) = (v = key(k); v isa AbstractString ? uppercase(strip(v)) : "")

    # ── axes 3 and 4 ────────────────────────────────────────────────────────────────────
    img = Float64.(raw)
    if nd == 2
        img = reshape(img, size(img, 1), size(img, 2), 1, 1)
    elseif nd == 3
        n3 = size(img, 3)
        t3 = skey("CTYPE3")
        if t3 == "WAVE" || (isempty(t3) && n3 == 1)
            img = reshape(img, size(img, 1), size(img, 2), n3, 1)
        elseif t3 == "TIME" || t3 == "MJD"
            img = reshape(img, size(img, 1), size(img, 2), 1, n3)
        else
            throw(ArgumentError(
                "$(basename(p)): a 3-D file with $n3 planes needs CTYPE3 = 'WAVE' or 'TIME' to " *
                "say what its third axis is" * (isempty(t3) ? "; it has none" : "; it has '$t3'")))
        end
    end
    nλ, nt = size(img, 3), size(img, 4)

    # ── orientation ─────────────────────────────────────────────────────────────────────
    cd1 = key("CDELT1"); cd2 = key("CDELT2")
    if startswith(skey("CTYPE1"), "RA")
        cd1 isa Real && cd1 < 0 && (img = reverse(img; dims = 1))
        cd2 isa Real && cd2 < 0 && (img = reverse(img; dims = 2))
    end

    # ── pixel scale ─────────────────────────────────────────────────────────────────────
    ps = pixsize === nothing ? fits_pixsize(p) : Float64(pixsize)
    (ps === nothing || !(ps > 0)) && throw(ArgumentError(
        "$(basename(p)): no pixel scale in the header (CDELT/CD with CUNIT, or PIXSIZE); " *
        "give pixsize in mas"))

    # ── reference pixel ─────────────────────────────────────────────────────────────────
    nx, ny = size(img, 1), size(img, 2)
    cx, cy = nx ÷ 2 + 1, ny ÷ 2 + 1
    off = (0.0, 0.0)
    r1, r2 = key("CRPIX1"), key("CRPIX2")
    if r1 isa Real && r2 isa Real
        # The reference pixel in the East-first array. A flip maps index i to n+1-i, which
        # for an even n moves the phase-centre pixel by one: a standard astronomical image
        # with CRPIX at nx÷2+1 comes out with its reference one pixel off the centre, and
        # that pixel is carried as an offset rather than lost.
        rx = startswith(skey("CTYPE1"), "RA")  && cd1 isa Real && cd1 < 0 ? nx + 1 - r1 : r1
        ry = startswith(skey("CTYPE2"), "DEC") && cd2 isa Real && cd2 < 0 ? ny + 1 - r2 : r2
        # The file's reference pixel is the origin: a reference pixel EAST of the phase
        # centre means the image content sits WEST of the origin by the same amount.
        (rx != cx || ry != cy) && (off = (-(rx - cx) * ps, -(ry - cy) * ps))
    end

    # ── wavelengths ─────────────────────────────────────────────────────────────────────
    λ  = fill(NaN, nλ); dλ = fill(NaN, nλ)
    if channels !== nothing && haskey(channels, "EFF_WAVE")
        w = Float64.(channels["EFF_WAVE"])
        length(w) == nλ || throw(ArgumentError(
            "$(basename(p)): CHANNELS has $(length(w)) rows for $nλ planes"))
        λ .= w
        haskey(channels, "EFF_BAND") && (dλ .= Float64.(channels["EFF_BAND"]))
    elseif skey("CTYPE3") == "WAVE"
        key("WAVEIRR", false) == true && throw(ArgumentError(
            "$(basename(p)): WAVEIRR = T but no CHANNELS table gives the wavelengths"))
        unit = get(_WAVE_UNITS, lowercase(strip(string(key("CUNIT3", "m")))), nothing)
        unit === nothing && throw(ArgumentError(
            "$(basename(p)): CUNIT3 = '$(key("CUNIT3"))' is not a wavelength unit"))
        c0, d0, r0 = key("CRVAL3"), key("CDELT3", 0.0), key("CRPIX3", 1.0)
        c0 isa Real && (λ .= (c0 .+ ((1:nλ) .- r0) .* d0) .* unit)
    end

    # ── epochs ──────────────────────────────────────────────────────────────────────────
    mjd = fill(NaN, nt)
    if epochs !== nothing && haskey(epochs, "MJD")
        m = Float64.(epochs["MJD"])
        length(m) == nt || throw(ArgumentError(
            "$(basename(p)): EPOCHS has $(length(m)) rows for $nt epochs"))
        mjd .= m
    elseif skey(nd == 3 ? "CTYPE3" : "CTYPE4") in ("TIME", "MJD")
        ax = nd == 3 ? 3 : 4
        key("TIMEIRR", false) == true && throw(ArgumentError(
            "$(basename(p)): TIMEIRR = T but no EPOCHS table gives the MJDs"))
        u = lowercase(strip(string(key("CUNIT$ax", "d"))))
        scale = u in ("d", "day", "days") ? 1.0 : u in ("s", "sec") ? 1 / 86400 :
                throw(ArgumentError("$(basename(p)): CUNIT$ax = '$u' is not a time unit"))
        c0, d0, r0 = key("CRVAL$ax"), key("CDELT$ax", 0.0), key("CRPIX$ax", 1.0)
        c0 isa Real && (mjd .= Float64(key("MJDREF", 0.0)) .+ (c0 .+ ((1:nt) .- r0) .* d0) .* scale)
    end

    # ── normalisation ───────────────────────────────────────────────────────────────────
    norm = skey("NORMALIZ")
    sums = [sum(@view img[:, :, i, j]) for i in 1:nλ, j in 1:nt]
    any(<=(0), sums) && throw(ArgumentError(
        "$(basename(p)): a plane sums to $(minimum(sums)); a model image needs positive flux " *
        "in every plane"))
    normalization in (nothing, :plane, :relative) || throw(ArgumentError(
        "normalization must be :plane or :relative, not $(repr(normalization))"))
    nrm = normalization !== nothing ? normalization : norm == "RELATIVE" ? :relative : :plane
    if nrm === :plane
        norm == "PLANE" && maximum(abs.(sums .- 1)) > 1e-4 && @warn "$(basename(p)): " *
            "NORMALIZ = PLANE but planes sum to $(extrema(sums)); renormalised"
        for i in 1:nλ, j in 1:nt; img[:, :, i, j] ./= sums[i, j]; end
    else
        norm == "PLANE" && length(sums) > 1 && @warn "$(basename(p)): NORMALIZ = PLANE, so " *
            "every plane already sums to 1 and there is no spectrum or light curve to keep"
        img ./= maximum(sums)
    end

    return ModelImage(img, ps, λ, dλ, mjd, nrm, off, p)
end

# A binary table's columns by name, or `nothing` when the file has no such extension.
function _fits_table(f, extname)
    for h in f
        h isa FITSIO.TableHDU || continue
        name = try read_key(h, "EXTNAME")[1] catch; "" end
        uppercase(strip(String(name))) == extname || continue
        return Dict(c => read(h, c) for c in FITSIO.colnames(h))
    end
    return nothing
end

"""
    model_image_plane(mi, λ, mjd; wave_tol = nothing, mjd_tol = 0.5) -> (iλ, it)

The plane of `mi` that a data point at wavelength `λ` (metres) and epoch `mjd` is evaluated
with: the nearest in wavelength and the nearest in time, each within a tolerance, never an
interpolation between planes.

- **Wavelength.** A single plane, or planes with no wavelength, serve every λ — a grey image.
  Otherwise the nearest plane, which must lie within half its `EFF_BAND` when the file gives
  one, else within a quarter of the gap to its neighbouring plane. That accepts the drift of a
  channel between nights — MIRC-X's channel 1 is 1.515382 µm one night and 1.515206 µm five
  nights later, 1.2e-4 apart — and refuses a different channel: a cube holding channels 1 and 4
  serves neither 2 nor 3. Further away is an error naming the cube's channels, never a silent
  substitution of a neighbouring channel's image. `wave_tol` (relative) overrides both.
- **Epoch.** A single epoch, or epochs with no MJD, serve every time. Otherwise the nearest,
  within `mjd_tol` days (default 0.5: one night's points share an epoch).
"""
function model_image_plane(mi::ModelImage, λ::Real, mjd::Real; wave_tol = nothing,
                           mjd_tol::Real = 0.5)
    iλ = 1
    if length(mi.wavelengths) > 1 && all(isfinite, mi.wavelengths)
        iλ = argmin(abs.(mi.wavelengths .- λ))
        d  = abs(mi.wavelengths[iλ] - λ)
        jλ = iλ                    # not reassigned, so the generator below captures it unboxed
        gap = minimum(abs(mi.wavelengths[j] - mi.wavelengths[jλ])
                      for j in eachindex(mi.wavelengths) if j != jλ)
        tol = wave_tol !== nothing ? Float64(wave_tol) * λ :
              isfinite(mi.bandwidths[iλ]) ? mi.bandwidths[iλ] / 2 : gap / 4
        d <= tol || throw(ArgumentError(
            "λ = $(round(λ * 1e6; sigdigits = 6)) µm is not one of the image's channels " *
            "($(join(round.(mi.wavelengths .* 1e6; sigdigits = 6), ", ")) µm); the nearest is " *
            "$(round(d / λ; sigdigits = 2)) away in relative terms. Give wave_tol to accept it."))
    end
    it = 1
    if length(mi.mjds) > 1 && all(isfinite, mi.mjds)
        it = argmin(abs.(mi.mjds .- mjd))
        abs(mi.mjds[it] - mjd) <= mjd_tol || throw(ArgumentError(
            "MJD $(round(mjd; digits = 4)) is not within $mjd_tol d of any of the image's " *
            "epochs ($(join(round.(mi.mjds; digits = 4), ", ")))"))
    end
    return (iλ, it)
end

# ─────────────────────────────────────────────────────────────────────────────
# Visibilities of a model image, cached per uv set
# ─────────────────────────────────────────────────────────────────────────────

"""
    ImageVisCache()

The visibilities an image component has computed, keyed by the uv points, wavelengths and
epochs they were computed at. An image's visibilities depend on nothing a fit varies — flux
and position are applied afterwards — so a fit transforms the image once per dataset and
reuses the result at every evaluation. `nfft` counts the planes transformed so far.

Holds the last `_IMAGE_CACHE_SIZE` uv sets, oldest dropped first. Guarded by a lock, since a
`FlatModel` is shared between the threads of a batch likelihood.
"""
mutable struct ImageVisCache
    entries :: Vector{Tuple{Matrix{Float64},Vector{Float64},Vector{Float64},Vector{ComplexF64}}}
    nfft    :: Int
    warned  :: Bool
    lock    :: ReentrantLock
end
ImageVisCache() = ImageVisCache([], 0, false, ReentrantLock())

const _IMAGE_CACHE_SIZE = 16

# Whether a cached per-point vector is what `wl` or `mjd` would expand to.
_same_per_point(c, v::Nothing) = all(isnan, c)
_same_per_point(c, v::AbstractVector) = length(c) == length(v) && isequal(c, v)
_same_per_point(c, v) = all(x -> isequal(x, Float64(v)), c)

# A per-point vector from what `eval_model` was given for `wl` or `mjd`: a vector, a scalar,
# or nothing (NaN, as in the resolver).
_per_point(v, n) = v === nothing ? fill(NaN, n) :
                   v isa AbstractVector ? Float64.(v) : fill(Float64(v), n)

"""
    model_image_vis(mi, uv; wl = nothing, mjd = nothing, cache = nothing, name = "image")
        -> Vector{ComplexF64}

The complex visibilities of the model image `mi` at the uv points `uv` (2×N, cycles/rad), each
point evaluated with the plane `model_image_plane` picks for its wavelength `wl` (metres) and
epoch `mjd` — per point, or one value for all. A plane is used as stored: one summing to 1
gives V(0) = 1, a `:relative` plane gives its own flux. `mi.offset` is applied as a shift.

A cube needs `wl` (and an epoch series `mjd`) to choose planes, and a point outside the cube's
channels or epochs is an error from `model_image_plane`.

A uv point beyond the image's Nyquist frequency, 1/(2 pixsize), is outside what the pixel grid
can describe: its visibility is set to 0, and the first time it happens for a given `cache` a
warning names the component (`name`) and the baseline. Fitting data there needs a finer image.

With `cache`, a uv set seen before is returned without a transform.
"""
function model_image_vis(mi::ModelImage, uv::AbstractMatrix; wl = nothing, mjd = nothing,
                         cache::Union{Nothing,ImageVisCache} = nothing, name = "image")
    n = size(uv, 2)
    cache === nothing &&
        return _model_image_vis(mi, Matrix{Float64}(uv), _per_point(wl, n), _per_point(mjd, n),
                                nothing, name)
    lock(cache.lock) do
        # Compared in place, so a hit -- every evaluation of a fit after the first -- copies
        # nothing. Float32 points compare equal to the Float64 copies made of them.
        for e in cache.entries
            size(e[1]) == size(uv) && e[1] == uv && _same_per_point(e[2], wl) &&
                _same_per_point(e[3], mjd) && return e[4]
        end
        uvf = Matrix{Float64}(uv)
        λs  = _per_point(wl, n)
        ts  = _per_point(mjd, n)
        V = _model_image_vis(mi, uvf, λs, ts, cache, name)
        length(cache.entries) >= _IMAGE_CACHE_SIZE && popfirst!(cache.entries)
        push!(cache.entries, (uvf, λs, ts, V))
        return V
    end
end

function _model_image_vis(mi::ModelImage, uv::Matrix{Float64}, λs, ts, cache, name)
    nx, ny, nλ, nt = size(mi.image)
    n = size(uv, 2)
    chromatic = nλ > 1 && all(isfinite, mi.wavelengths)
    temporal  = nt > 1 && all(isfinite, mi.mjds)
    chromatic && any(!isfinite, λs) && throw(ArgumentError(
        "$name: the image has $nλ channels, so each uv point needs a wavelength to choose one"))
    temporal && any(!isfinite, ts) && throw(ArgumentError(
        "$name: the image has $nt epochs, so each uv point needs an MJD to choose one"))

    # Which plane each point uses, deciding once per distinct (λ, MJD).
    planes = Dict{Tuple{Float64,Float64},Tuple{Int,Int}}()
    groups = Dict{Tuple{Int,Int},Vector{Int}}()
    for i in 1:n
        key = (chromatic ? λs[i] : NaN, temporal ? ts[i] : NaN)
        p = get!(() -> model_image_plane(mi, key[1], key[2]), planes, key)
        push!(get!(groups, p, Int[]), i)
    end

    # The NFFT's nodes are in cycles per pixel and must lie in [-1/2, 1/2): exactly 1/2 is the
    # same frequency as -1/2, beyond it the pixel grid has nothing to say.
    k = uv .* (mi.pixsize / _MAS2RAD_PM)
    beyond = vec(any(abs.(k) .> 0.5; dims = 1))
    k[k .== 0.5] .= -0.5
    if any(beyond) && (cache === nothing || !cache.warned)
        bmax = maximum(sqrt.(uv[1, :] .^ 2 .+ uv[2, :] .^ 2)) / 1e6
        @warn "$name: uv points out to $(round(bmax; sigdigits = 4)) Mλ lie beyond the " *
              "image's Nyquist frequency ($(round(_MAS2RAD_PM / (2mi.pixsize) / 1e6; sigdigits = 4)) " *
              "Mλ at $(mi.pixsize) mas per pixel); their visibilities are set to 0"
        cache === nothing || (cache.warned = true)
    end

    V = zeros(ComplexF64, n)
    for ((iλ, it), idx) in groups
        keep = filter(i -> !beyond[i], idx)
        isempty(keep) && continue
        plan = plan_nfft(k[:, keep], (nx, ny), m = 4, σ = 2.0; fftflags = FFTW.ESTIMATE)
        V[keep] .= (plan * ComplexF64.(@view mi.image[:, :, iλ, it]))::Vector{ComplexF64}
        cache === nothing || (cache.nfft += 1)
    end

    ox, oy = mi.offset
    if ox != 0 || oy != 0
        @. V *= cis(-2π * (uv[1, :] * ox + uv[2, :] * oy) / _MAS2RAD_PM)
    end
    return V
end
