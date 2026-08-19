using NFFT
using DelimitedFiles
using FITSIO
using Random

import TOML

# ─── Configuration structs ────────────────────────────────────────────────────

Base.@kwdef mutable struct FacilityConfig
    name::String                = ""
    lat::Float64                = 0.0
    lon::Float64                = 0.0
    alt::Float64                = 0.0
    # Extra transmission factor applied on top of CombinerConfig.transmission. Combiner
    # configs transcribed from ASPRO already carry the array transmission end-to-end, so this
    # is 1.0 for them; it exists for arrays whose combiners are specified instrument-only.
    throughput::Float64         = 1.0
    wind_speed::Float64         = 10.0
    r0::Float64                 = 0.1     # Fried parameter at 500 nm, zenith [m]
    seeing::Float64             = 0.0     # arcsec at 500 nm; 0 => derived from r0
    t0::Float64                 = 0.0     # coherence time [s]; 0 => derived from r0/wind_speed
    ao::Union{AOConfig,Nothing} = nothing # adaptive optics, or nothing for an uncorrected array
    ntel::Int                   = 0
    tel_names::Vector{String}   = String[]
    sta_names::Vector{String}   = String[]
    tel_diams::Vector{Float64}  = Float64[]
    tel_gain::Vector{Float64}   = Float64[]
    sta_index::Vector{Int}      = Int[]
    sta_xyz::Matrix{Float64}    = zeros(Float64, 0, 3)
    delay_lengths::Vector{Float64} = Float64[]   # per-telescope delay line length (m), empty if unknown
end

Base.@kwdef mutable struct TargetConfig
    target_id::Int      = 1
    target::String      = "OBJECT"
    raep0::Float64      = 0.0
    decep0::Float64     = 0.0
    equinox::Float64    = 2000.0
    ra_err::Float64     = 0.0
    dec_err::Float64    = 0.0
    sysvel::Float64     = 0.0
    veltyp::String      = "LSR"
    veldef::String      = "OPTICAL"
    pmra::Float64       = 0.0
    pmdec::Float64      = 0.0
    pmra_err::Float64   = 0.0
    pmdec_err::Float64  = 0.0
    parallax::Float64   = 0.0
    para_err::Float64   = 0.0
    spectyp::String     = "G0"
end

"""
    CombinerConfig

Beam-combiner properties. The fields mirror ASPRO 2's `FocalInstrumentSetup` one for one, so
that a combiner config can be transcribed directly from ASPRO's interferometer XML
(`aspro-conf`, e.g. `CHARA.xml`) and reproduce its predictions.

- `transmission`: **end-to-end** transmission, array *and* instrument, excluding quantum
  efficiency, Strehl and atmosphere (ASPRO `transmission`). Because it already covers the
  array, `FacilityConfig.throughput` is 1.0 for any facility whose combiners are specified
  this way.
- `instrument_visibility`: instrumental contrast (ASPRO `instrumentVisibility`).
- `dit`: detector integration time in seconds — a **fixed** instrument property, shortened
  only when the detector would saturate (ASPRO `dit`).
- `total_int_time`: total on-source time per uv point, seconds (ASPRO
  `defaultTotalIntegrationTime`).
- `detector_saturation`: full-well in electrons per pixel (ASPRO `detectorSaturation`).
- `flux_frac_fringes` / `flux_frac_photometry`: how the light splits between the
  interferometric and photometric channels (ASPRO `fracFluxIn…`).
- `vis_cal_err`: systematic error on `|V|`, as a fraction. ASPRO stores this as
  `instrumentVisibilityBias` in **percent**, so `instrumentVisibilityBias = 1` is
  `vis_cal_err = 0.01`. The systematic on `V²` is twice this.
- `phase_cal_err`: systematic error on phases and closure phases, in degrees (ASPRO
  `instrumentPhaseBias`).
"""
Base.@kwdef mutable struct CombinerConfig
    name::String                    = ""
    transmission::Float64           = 0.01
    instrument_visibility::Float64  = 0.9
    dit::Float64                    = 0.01
    total_int_time::Float64         = 600.0
    detector_saturation::Float64    = 50000.0
    n_pix_fringe::Int               = 128
    n_pix_photometry::Int           = 4
    flux_frac_photometry::Float64   = 0.2
    flux_frac_fringes::Float64      = 0.8
    read_noise::Float64             = 10.0
    quantum_efficiency::Float64     = 1.0
    vis_cal_err::Float64            = 0.01
    phase_cal_err::Float64          = 1.0
    # "ao"          use the facility's AO model (the physical default)
    # "fixed_spica" reproduce ASPRO's hardcoded SPICA Strehl: 0.25 for seeing <= 0.7",
    #               0.15 for seeing <= 1.0", 0.10 otherwise, independent of magnitude and
    #               elevation. ASPRO applies this in NoiseService even though its own
    #               published Strehl plots use the AO model, so the two disagree; this exists
    #               to reproduce ASPRO's *noise*.
    strehl_model::String            = "ao"
end

Base.@kwdef mutable struct WaveConfig
    combiner::String            = ""
    mode::String                = ""
    λ::Vector{Float64}          = Float64[]
    δλ::Vector{Float64}         = Float64[]
end

# Backward-compatible type aliases
const facility_info = FacilityConfig
const obsv_info     = TargetConfig
const combiner_info = CombinerConfig
const wave_info     = WaveConfig

# ─── Config file resolution ───────────────────────────────────────────────────

const _CONFIGS_DIR = joinpath(@__DIR__, "configs")

function _resolve_config(path::AbstractString)
    isfile(path) && return path
    fallback = joinpath(_CONFIGS_DIR, basename(path))
    isfile(fallback) && return fallback
    # Try adding .toml extension
    for p in (path * ".toml", joinpath(_CONFIGS_DIR, basename(path) * ".toml"))
        isfile(p) && return p
    end
    return path  # let downstream code raise the error
end

# ─── TOML readers (primary) ──────────────────────────────────────────────────

function _read_facility_toml(path)
    d = TOML.parsefile(path)
    tels = get(d, "telescope", [])
    ntel = length(tels)
    xyz = zeros(ntel, 3)
    for (i, t) in enumerate(tels)
        xyz[i, :] .= Float64.(t["xyz"])
    end
    r0 = Float64(get(d, "r0", 0.1))
    wind_speed = Float64(get(d, "wind_speed", 10.0))
    # seeing and t0 are the physical inputs; fall back to the legacy r0/wind_speed pair so
    # configs written before the AO model keep working unchanged.
    seeing = Float64(get(d, "seeing", 0.0));  seeing <= 0 && (seeing = seeing_from_r0(r0))
    t0     = Float64(get(d, "t0", 0.0));      t0     <= 0 && (t0 = 0.31 * r0 / wind_speed)

    aod = get(d, "ao", nothing)
    ao = aod === nothing ? nothing : AOConfig(
        band          = String(get(aod, "band", "V")),
        n_subpupils   = Int(get(aod, "n_subpupils", 36)),
        n_actuators   = Int(get(aod, "n_actuators", 24)),
        dit           = Float64(get(aod, "dit", 2.0e-3)),
        ron           = Float64(get(aod, "ron", 5.0)),
        qe            = Float64(get(aod, "qe", 0.7)),
        transmission  = Float64(get(aod, "transmission", 0.0)),
        mag_offset    = Float64(get(aod, "mag_offset", 0.0)),
        strehl_max    = Float64(get(aod, "strehl_max", 0.0)),
    )

    FacilityConfig(
        name       = String(get(d, "name", "")),
        lat        = Float64(get(d, "lat", 0.0)),
        lon        = Float64(get(d, "lon", 0.0)),
        alt        = Float64(get(d, "alt", 0.0)),
        throughput = Float64(get(d, "throughput", 1.0)),
        wind_speed = wind_speed,
        r0         = r0,
        seeing     = seeing,
        t0         = t0,
        ao         = ao,
        ntel       = ntel,
        tel_names  = [String(get(t, "name", "T$i"))    for (i,t) in enumerate(tels)],
        sta_names  = [String(get(t, "station", get(t, "name", "S$i"))) for (i,t) in enumerate(tels)],
        tel_diams  = [Float64(get(t, "diameter", 1.0))  for t in tels],
        tel_gain   = [Float64(get(t, "gain", 1.0))      for t in tels],
        sta_index  = [Int(get(t, "index", i))            for (i,t) in enumerate(tels)],
        sta_xyz    = xyz,
        delay_lengths = [Float64(get(t, "delay_length", 0.0)) for t in tels],
    )
end

function _read_target_toml(path)
    d = TOML.parsefile(path)
    TargetConfig(
        target_id = Int(get(d, "target_id", 1)),
        target    = String(get(d, "target", "OBJECT")),
        raep0     = Float64(get(d, "raep0", 0.0)),
        decep0    = Float64(get(d, "decep0", 0.0)),
        equinox   = Float64(get(d, "equinox", 2000.0)),
        ra_err    = Float64(get(d, "ra_err", 0.0)),
        dec_err   = Float64(get(d, "dec_err", 0.0)),
        sysvel    = Float64(get(d, "sysvel", 0.0)),
        veltyp    = String(get(d, "veltyp", "LSR")),
        veldef    = String(get(d, "veldef", "OPTICAL")),
        pmra      = Float64(get(d, "pmra", 0.0)),
        pmdec     = Float64(get(d, "pmdec", 0.0)),
        pmra_err  = Float64(get(d, "pmra_err", 0.0)),
        pmdec_err = Float64(get(d, "pmdec_err", 0.0)),
        parallax  = Float64(get(d, "parallax", 0.0)),
        para_err  = Float64(get(d, "para_err", 0.0)),
        spectyp   = String(get(d, "spectyp", "G0")),
    )
end

function _read_combiner_toml(path)
    d = TOML.parsefile(path)
    # Accept the pre-ASPRO field names so combiner configs that have not been transcribed
    # yet keep working: `throughput_fringes` folded into a single end-to-end transmission,
    # `vis`, `incoh_int_time`, and an absolute `v2_cal_err`.
    transmission = if haskey(d, "transmission")
        Float64(d["transmission"])
    else
        Float64(get(d, "throughput_fringes", 0.1)) * Float64(get(d, "int_trans", 1.0))
    end
    vis_cal = haskey(d, "vis_cal_err") ? Float64(d["vis_cal_err"]) :
                                         0.5 * Float64(get(d, "v2_cal_err", 0.02))
    CombinerConfig(
        name                 = String(get(d, "name", "")),
        transmission         = transmission,
        instrument_visibility = Float64(get(d, "instrument_visibility", get(d, "vis", 0.9))),
        dit                  = Float64(get(d, "dit", 0.01)),
        total_int_time       = Float64(get(d, "total_int_time", get(d, "incoh_int_time", 600.0))),
        detector_saturation  = Float64(get(d, "detector_saturation", 50000.0)),
        n_pix_fringe         = Int(get(d, "n_pix_fringe", 128)),
        n_pix_photometry     = Int(get(d, "n_pix_photometry", 4)),
        flux_frac_photometry = Float64(get(d, "flux_frac_photometry", 0.2)),
        flux_frac_fringes    = Float64(get(d, "flux_frac_fringes", 0.8)),
        read_noise           = Float64(get(d, "read_noise", 10.0)),
        quantum_efficiency   = Float64(get(d, "quantum_efficiency", 1.0)),
        vis_cal_err          = vis_cal,
        phase_cal_err        = Float64(get(d, "phase_cal_err", 1.0)),
        strehl_model         = String(get(d, "strehl_model", "ao")),
    )
end

function _read_wave_toml(path)
    d = TOML.parsefile(path)
    WaveConfig(
        combiner = String(get(d, "combiner", "")),
        mode     = String(get(d, "mode", "")),
        λ        = Float64.(d["lambda"]),
        δλ       = Float64.(d["delta_lambda"]),
    )
end

# ─── Legacy .txt readers (backward compatibility) ────────────────────────────

function _read_facility_txt(path)
    f = readdlm(path)
    _get(key) = f[(LinearIndices(f .== key))[findall(f .== key)], 3]
    n = Int(_get("ntel")[1])
    _row(key, n) = f[(LinearIndices(f .== key))[findall(f .== key)], 3:n+2]
    _has(key) = !isempty(findall(f .== key))
    flat_xyz = Float64.(f[(LinearIndices(f .== "sta_xyz"))[findall(f .== "sta_xyz")], 3:n*3+2])
    xyz = reshape(vec(flat_xyz), 3, n)'  # (ntel, 3)
    # Legacy .txt files predate delay_lengths; without it the delay-line check in plan.jl
    # silently falls back to a hardcoded 45.7 m, so say so rather than failing quietly.
    _has("delay_lengths") || @warn "Facility file $path has no 'delay_lengths' row; delay-line checks will use a default length." maxlog=1
    FacilityConfig(
        name       = String(_get("name")[1]),
        lat        = Float64(_get("lat")[1]),
        lon        = Float64(_get("lon")[1]),
        alt        = Float64(_get("alt")[1]),
        throughput = Float64(_get("throughput")[1]),
        wind_speed = Float64(_get("wind_speed")[1]),
        r0         = Float64(_get("r0")[1]),
        ntel       = n,
        tel_names  = String.(vec(_row("tel_names", n))),
        sta_names  = String.(vec(_row("sta_names", n))),
        tel_diams  = Float64.(vec(_row("tel_diams", n))),
        tel_gain   = Float64.(vec(_row("tel_gain", n))),
        sta_index  = Int.(vec(_row("sta_index", n))),
        sta_xyz    = xyz,
        delay_lengths = _has("delay_lengths") ? Float64.(vec(_row("delay_lengths", n))) : Float64[],
    )
end

function _read_target_txt(path)
    f = readdlm(path)
    _get(key) = f[(LinearIndices(f .== key))[findall(f .== key)], 3]
    TargetConfig(
        target_id = Int(_get("target_id")[1]),
        target    = String(_get("target")[1]),
        raep0     = Float64(_get("raep0")[1]),
        decep0    = Float64(_get("decep0")[1]),
        equinox   = Float64(_get("equinox")[1]),
        ra_err    = Float64(_get("ra_err")[1]),
        dec_err   = Float64(_get("dec_err")[1]),
        sysvel    = Float64(_get("sysvel")[1]),
        veltyp    = String(_get("veltyp")[1]),
        veldef    = String(_get("veldef")[1]),
        pmra      = Float64(_get("pmra")[1]),
        pmdec     = Float64(_get("pmdec")[1]),
        pmra_err  = Float64(_get("pmra_err")[1]),
        pmdec_err = Float64(_get("pmdec_err")[1]),
        parallax  = Float64(_get("parallax")[1]),
        para_err  = Float64(_get("para_err")[1]),
        spectyp   = String(_get("spectyp")[1]),
    )
end

function _read_combiner_txt(path)
    f = readdlm(path)
    _get(key) = f[(LinearIndices(f .== key))[findall(f .== key)], 3]
    _has(key) = !isempty(findall(f .== key))
    _num(key, default) = _has(key) ? Float64(_get(key)[1]) : default
    # Legacy whitespace format predates the ASPRO-aligned schema; map what it does carry.
    transmission = _has("transmission") ? _num("transmission", 0.01) :
                   _num("throughput_fringes", 0.1) * _num("int_trans", 1.0)
    CombinerConfig(
        name                 = String(_get("name")[1]),
        transmission         = transmission,
        instrument_visibility = _has("instrument_visibility") ? _num("instrument_visibility", 0.9) : _num("vis", 0.9),
        dit                  = _num("dit", 0.01),
        total_int_time       = _has("total_int_time") ? _num("total_int_time", 600.0) : _num("incoh_int_time", 600.0),
        detector_saturation  = _num("detector_saturation", 50000.0),
        n_pix_fringe         = Int(_num("n_pix_fringe", 128)),
        n_pix_photometry     = Int(_num("n_pix_photometry", 4)),
        flux_frac_photometry = _num("flux_frac_photometry", 0.2),
        flux_frac_fringes    = _num("flux_frac_fringes", 0.8),
        read_noise           = _num("read_noise", 10.0),
        quantum_efficiency   = _num("quantum_efficiency", 1.0),
        vis_cal_err          = _has("vis_cal_err") ? _num("vis_cal_err", 0.01) : 0.5*_num("v2_cal_err", 0.02),
        phase_cal_err        = _num("phase_cal_err", 1.0),
    )
end

function _read_wave_txt(path)
    f = readdlm(path)
    _get(key) = String.(f[(LinearIndices(f .== key))[findall(f .== key)], 3])[1]
    WaveConfig(
        combiner = _get("combiner"),
        mode     = _get("mode"),
        λ        = Float64.(f[3:end, 1]),
        δλ       = Float64.(f[3:end, 2]),
    )
end

# ─── Public API: auto-detect format by extension ─────────────────────────────

"""
    list_configs([dir]) -> (; facilities, combiners, wavelengths, targets, unknown, wave_combiner)

Enumerate the configuration files that [`read_facility_file`](@ref), [`read_comb_file`](@ref),
[`read_wave_file`](@ref) and [`read_obs_file`](@ref) can load, classified by what they
actually contain. `dir` defaults to the configs shipped with OITOOLS.

Nothing about a filename says what kind of config it is, so the classification reads each
file: a **facility** has a site (`lat`), a **wavelength setup** has `lambda` and names its
`combiner`, a **target** has `raep0`, and anything else carrying `transmission`/`read_noise`
is a **combiner**. Files that parse but match nothing land in `unknown`, and so do files that
fail to parse — enumerating a directory should never throw.

`wave_combiner` maps each wavelength setup to the combiner it belongs to. That link exists
only inside the file, so it is the only way to offer, say, MIRC-X's spectral modes without
hard-coding them.

```jldoctest
julia> c = list_configs();

julia> c.facilities
5-element Vector{String}:
 "CHARA"
 "VLTI_AT_large"
 "VLTI_AT_medium"
 "VLTI_AT_small"
 "VLTI_UT"

julia> c.wave_combiner["MIRCX_LOWH"]
"MIRCX"
```
"""
function list_configs(dir::AbstractString = _CONFIGS_DIR)
    facilities  = String[]
    combiners   = String[]
    wavelengths = String[]
    targets     = String[]
    unknown     = String[]
    wave_combiner = Dict{String,String}()

    isdir(dir) || return (; facilities, combiners, wavelengths, targets, unknown, wave_combiner)

    for f in sort(readdir(dir))
        endswith(f, ".toml") || continue
        name = f[1:end-5]
        d = try
            TOML.parsefile(joinpath(dir, f))
        catch
            push!(unknown, name); continue
        end
        # Order matters: CHARA.toml carries a throughput/transmission key as well as a site,
        # so the facility test has to come before the combiner test.
        if haskey(d, "lat")
            push!(facilities, name)
        elseif haskey(d, "lambda") && haskey(d, "combiner")
            push!(wavelengths, name)
            wave_combiner[name] = string(d["combiner"])
        elseif haskey(d, "raep0") || haskey(d, "decep0")
            push!(targets, name)
        elseif any(haskey(d, k) for k in ("transmission", "read_noise", "int_trans"))
            push!(combiners, name)
        else
            push!(unknown, name)
        end
    end
    return (; facilities, combiners, wavelengths, targets, unknown, wave_combiner)
end


"""
    predict_errors(facility, combiner, wavelength; mag, visamp=1.0, elevation_deg=90.0,
                   mag_ao=nothing, channel=nothing) -> NamedTuple

Predicted per-channel uncertainties for one setup, **without running a simulation or writing
a file**.

`simulate` returns the result of writing an OIFITS and keeps nothing in memory, so the only
way to ask "what SNR would I get?" used to be to simulate to disk and read it back. This
runs the same noise model directly.

Returns, one entry per spectral channel of `wavelength`:

| field | meaning |
|---|---|
| `λ`, `δλ` | channel centre and width, m |
| `strehl` | coupling efficiency actually used |
| `dit` | integration time per frame, s (possibly shortened for saturation) |
| `nframes` | frames coadded |
| `nphot` | photons per telescope per DIT, fringe channel |
| `sigma_v2` | σ(V²), including the `vis_cal_err` systematic |
| `sigma_cp` | σ(closure phase), degrees, including `phase_cal_err` |
| `sigma_visamp` | σ of the visibility amplitude |
| `sigma_visphi` | σ(visibility phase), degrees |

`mag` follows the same rules as `simulate`: a scalar, a `Dict` of band ⇒ magnitude, or one
value per channel. `visamp` is the visibility amplitude at which to evaluate the errors —
errors depend on it, so the default of 1.0 gives the unresolved (best) case.

Mirrors `simulate`'s internals exactly; see `demos/validate_noise_model.jl`, which checks
this model against ASPRO 2's published curves.
"""
function predict_errors(facility::FacilityConfig, combiner::CombinerConfig,
                        wavelength::WaveConfig;
                        mag = 2.0, visamp::Real = 1.0, elevation_deg::Real = 90.0,
                        mag_ao::Union{Nothing,Real} = nothing)
    λ, δλ = wavelength.λ, wavelength.δλ
    nw    = length(λ)
    nw > 0 || throw(ArgumentError("wavelength config has no channels"))
    mags  = resolve_magnitudes(mag, λ)
    m_ao  = mag_ao === nothing ? _default_mag_ao(mags, λ, facility) : Float64(mag_ao)
    elev  = [Float64(elevation_deg)]

    N,  dit, nframes = photons_per_telescope(mags, λ, δλ, facility, combiner, elev, m_ao;
                                             channel = :fringes)
    Np, _,   _       = photons_per_telescope(mags, λ, δλ, facility, combiner, elev, m_ao;
                                             channel = :photometry)

    ntel  = facility.ntel
    v2st  = hcat([[i, j] for i in 1:ntel for j in (i+1):ntel]...)
    sq, vc, vk, vp = correlated_flux_coefficients(N, Np, combiner, ntel, v2st)

    strehl   = Vector{Float64}(undef, nw)
    nphot    = Vector{Float64}(undef, nw)
    sigma_v2 = Vector{Float64}(undef, nw)
    sigma_cp = Vector{Float64}(undef, nw)
    sigma_va = Vector{Float64}(undef, nw)
    sigma_vp = Vector{Float64}(undef, nw)

    for w in 1:nw
        σc = complex_vis_error(visamp, sq[1, 1, w], vc[1, 1, w], vk[1, 1, w], vp[1, 1, w],
                               nframes)
        v2 = visamp^2
        s2 = vis2_error_from_sigma(v2, σc)
        sigma_v2[w] = _with_systematics(s2, 2 * combiner.vis_cal_err * abs(v2))
        sigma_va[w] = _with_systematics(σc, combiner.vis_cal_err * abs(visamp))
        sigma_vp[w] = _with_systematics(min(180 / π * σc / max(abs(visamp), eps()), 180.0),
                                        combiner.phase_cal_err)
        # Closure phase: three baselines' statistical phase errors in quadrature, then the
        # calibration systematic added ONCE (not once per baseline). Same expression as
        # demos/validate_noise_model.jl, which is the form checked against ASPRO 2.
        σφ_stat     = 180 / π * sigma_v2[w] / 2
        sigma_cp[w] = _with_systematics(sqrt(3) * σφ_stat, combiner.phase_cal_err)
        strehl[w]   = combiner.strehl_model == "fixed_spica" ?
                      _spica_fixed_strehl(facility.seeing) :
                      coupling_efficiency(λ[w], _tel_diam(facility), facility.seeing,
                                          facility.t0, m_ao, facility.ao;
                                          elevation_deg = elevation_deg)
        nphot[w] = N[1, 1, w]
    end

    return (; λ, δλ, strehl, dit, nframes, nphot,
              sigma_v2, sigma_cp, sigma_visamp = sigma_va, sigma_visphi = sigma_vp)
end

# Mean telescope diameter, or 1 m when a config omits them.
_tel_diam(f::FacilityConfig) = isempty(f.tel_diams) ? 1.0 : sum(f.tel_diams) / length(f.tel_diams)

# simulate derives the AO guide magnitude from `mag` and the AO band when mag_ao is absent.
function _default_mag_ao(mags::AbstractVector, λ::AbstractVector, facility::FacilityConfig)
    facility.ao === nothing && return mags[cld(length(mags), 2)]
    b = band_by_name(facility.ao.band)
    b === nothing && return mags[cld(length(mags), 2)]
    return resolve_magnitudes_at(mags, λ, b.lambda)
end

# Value of the per-channel magnitude vector nearest a given wavelength.
function resolve_magnitudes_at(mags::AbstractVector, λ::AbstractVector, λ0::Real)
    _, i = findmin(abs.(λ .- λ0))
    return mags[i]
end

"""
    read_facility_file(path) -> FacilityConfig

Read an interferometric array: site coordinates, telescopes, station positions, delay-line
lengths, seeing and AO parameters.

`path` may be a bare name (`"CHARA"`), a name with extension, or a full path. Bare names are
looked up in OITOOLS' own `src/configs` directory, so the shipped configs work without a
path; see [`list_configs`](@ref) for what is available. A `.toml` file is read with the
current schema, anything else with the legacy whitespace-delimited reader.
"""
function read_facility_file(path)
    p = _resolve_config(path)
    endswith(p, ".toml") ? _read_facility_toml(p) : _read_facility_txt(p)
end

"""
    read_obs_file(path) -> TargetConfig

Read a target definition (name, RA/Dec in degrees, proper motion, spectral type) — the
OI_TARGET header `simulate` writes.

`path` may be a bare name (`"CHARA"`), a name with extension, or a full path. Bare names are
looked up in OITOOLS' own `src/configs` directory, so the shipped configs work without a
path; see [`list_configs`](@ref) for what is available. A `.toml` file is read with the
current schema, anything else with the legacy whitespace-delimited reader.
"""
function read_obs_file(path)
    p = _resolve_config(path)
    endswith(p, ".toml") ? _read_target_toml(p) : _read_target_txt(p)
end

"""
    read_comb_file(path) -> CombinerConfig

Read a beam combiner: transmission, detector properties, integration time and the
calibration error floors used by the noise model.

`path` may be a bare name (`"CHARA"`), a name with extension, or a full path. Bare names are
looked up in OITOOLS' own `src/configs` directory, so the shipped configs work without a
path; see [`list_configs`](@ref) for what is available. A `.toml` file is read with the
current schema, anything else with the legacy whitespace-delimited reader.
"""
function read_comb_file(path)
    p = _resolve_config(path)
    endswith(p, ".toml") ? _read_combiner_toml(p) : _read_combiner_txt(p)
end

"""
    read_wave_file(path) -> WaveConfig

Read a spectral setup: per-channel wavelengths and bandwidths, plus the name of the
combiner it belongs to.

`path` may be a bare name (`"CHARA"`), a name with extension, or a full path. Bare names are
looked up in OITOOLS' own `src/configs` directory, so the shipped configs work without a
path; see [`list_configs`](@ref) for what is available. A `.toml` file is read with the
current schema, anything else with the legacy whitespace-delimited reader.
"""
function read_wave_file(path)
    p = _resolve_config(path)
    endswith(p, ".toml") ? _read_wave_toml(p) : _read_wave_txt(p)
end


# zero_point_flux and the photometric band table live in atmosphere.jl, included first.

"""
    resolve_magnitudes(mag, λ) -> Vector{Float64}

Expand a magnitude specification onto the `length(λ)` spectral channels.

`mag` may be

- a `Real`: the same magnitude in every channel, i.e. the source is assumed to have Vega
  colours. Fine within one band, wrong across a mode like MIRC-X `Low_J` that spans J and H.
- an `AbstractDict` mapping band names to magnitudes, e.g. `Dict("J"=>2.1,"H"=>1.8)`, or
  what [`magnitudes_from_simbad`](@ref) returns. Interpolated linearly in log(λ) between the
  band centres and held flat outside them.
- an `AbstractVector` of length `length(λ)`: one magnitude per channel.
"""
function resolve_magnitudes(mag, λ)
    nw = length(λ)
    if mag isa Real
        return fill(Float64(mag), nw)
    elseif mag isa AbstractVector
        length(mag) == nw || throw(ArgumentError(
            "mag has $(length(mag)) entries but there are $nw spectral channels"))
        return Float64.(collect(mag))
    elseif mag isa AbstractDict
        pts = Tuple{Float64,Float64}[]
        for (k, v) in mag
            b = band_by_name(String(k))
            b === nothing && continue
            push!(pts, (b.lambda, Float64(v)))
        end
        isempty(pts) && throw(ArgumentError(
            "none of the magnitude keys $(collect(keys(mag))) match a known photometric band"))
        sort!(pts, by=first)
        out = Vector{Float64}(undef, nw)
        for w in 1:nw
            l = λ[w]
            if l <= pts[1][1]
                out[w] = pts[1][2]
            elseif l >= pts[end][1]
                out[w] = pts[end][2]
            else
                k = findlast(p -> p[1] <= l, pts)
                l1, m1 = pts[k]; l2, m2 = pts[k+1]
                t = (log(l) - log(l1)) / (log(l2) - log(l1))
                out[w] = m1 + t * (m2 - m1)
            end
        end
        return out
    end
    throw(ArgumentError("mag must be a Real, a Dict of band=>magnitude, or a Vector per channel"))
end

"""
    photons_per_telescope(mags, λ, δλ, facility, combiner, elevation_deg, mag_ao;
                          channel=:fringes, spec_weight=nothing)

Photons detected per telescope, per spectral channel, per detector frame, as an
`(ntel, nhours, nwavs)` array, together with the frame time and frame count.

The budget is stated once, here, and nowhere else:

    N = F0(λ) · 10^(-m(λ)/2.5) · A_tel · δλ · DIT
        · T_atm(λ)             atmospheric transmission (see `atm_transmission`)
        · facility.throughput  telescopes and beam train only
        · combiner.transmission   end-to-end array + instrument (ASPRO convention)
        · flux_frac            fraction split into this channel
        · QE
        · S(λ, elevation, m_ao)   Strehl, or the seeing-limited coupling with no AO

Note in particular that `facility.throughput` no longer doubles as a coupling term: fibre
injection is the Strehl now, so a facility that has AO must carry an `[ao]` block or it will
be modelled as an uncorrected aperture.
"""
function photons_per_telescope(mags, λ, δλ, facility, combiner, elevation_deg, mag_ao;
                               channel::Symbol=:fringes, spec_weight=nothing)
    ntel   = facility.ntel
    nwavs  = length(λ)
    nhours = length(elevation_deg)

    flux_frac = channel === :fringes ? combiner.flux_frac_fringes :
                                       combiner.flux_frac_photometry

    QE      = combiner.quantum_efficiency
    t_total = combiner.total_int_time
    seeing  = facility.seeing
    t0      = facility.t0

    w_spec = spec_weight === nothing ? ones(nwavs) : spec_weight

    # Photons per second per telescope, before the flux split, for saturation accounting.
    rate = zeros(ntel, nhours, nwavs)
    for w in 1:nwavs
        F0     = zero_point_flux(λ[w])                     # ph/s/m^2/um
        δλ_μm  = δλ[w] * 1e6
        T_atm  = atm_transmission(λ[w])
        base   = F0 * 10^(-mags[w] / 2.5) * δλ_μm *
                 T_atm * facility.throughput * combiner.transmission * QE * w_spec[w]
        for h in 1:nhours, t in 1:ntel
            D = facility.tel_diams[t]
            S = combiner.strehl_model == "fixed_spica" ? _spica_fixed_strehl(seeing) :
                coupling_efficiency(λ[w], D, seeing, t0, mag_ao, facility.ao;
                                    elevation_deg=elevation_deg[h])
            rate[t, h, w] = base * (π/4 * D^2) * S
        end
    end

    # DIT is a fixed instrument property, shortened only to avoid saturating the detector.
    # The worst case is the peak pixel of the interferometric channel, which sees the flux of
    # every telescope spread over nbPixInterferometry pixels.
    peak_per_sec = combiner.flux_frac_fringes * ntel *
                   maximum(rate) / max(combiner.n_pix_fringe, 1)
    max_dit = peak_per_sec > 0 ? combiner.detector_saturation / peak_per_sec : Inf
    DIT = min(combiner.dit, max_dit, t_total)
    if DIT < combiner.dit
        @info "DIT shortened from $(combiner.dit) s to $(round(DIT, sigdigits=3)) s to avoid detector saturation" maxlog=1
    end
    N_frames = max(1.0, t_total / DIT)

    return rate .* (DIT * flux_frac), DIT, N_frames
end

# ASPRO's hardcoded SPICA Strehl (NoiseService.initParameters), pending CHARA AO
# characterisation in the visible. Independent of magnitude and elevation.
_spica_fixed_strehl(seeing) = seeing <= 0.7 ? 0.25 : (seeing <= 1.0 ? 0.15 : 0.10)

# Combine a statistical error with a systematic floor. One rule, used for every observable,
# so that V2, VISAMP, VISPHI, T3AMP and T3PHI stay on the same footing.
_with_systematics(σ_stat, σ_sys) = sqrt(σ_stat^2 + σ_sys^2)

"""
    correlated_flux_coefficients(N_interf, N_phot, combiner, ntel, v2_stations)

Coefficients of ASPRO's squared-correlated-flux error model, each an
`(nv2, nhours, nwavs)` array:

- `sq_coef`:   squared correlated flux at `V² = 1`, i.e. `(sqrt(N_i·N_j)·vinst)²`
- `var_coef`:  coefficient multiplying the squared correlated flux in its variance
- `var_const`: the constant part of that variance
- `v2phot`:    variance the photometric channels contribute to `V²`

The variance of the squared correlated flux is `var = sq_coef·V²·var_coef + var_const`. The
constant part is what gives `σ(V²)` a floor as `V² → 0`: on a fully resolved baseline the
error stops shrinking rather than vanishing.

The photon-noise term counts **all** telescopes feeding the interferometric channel, not just
the two forming the baseline — an all-in-one combiner overlaps every beam, so every beam
contributes photon noise to every baseline.
"""
function correlated_flux_coefficients(N_interf, N_phot, combiner, ntel, v2_stations)
    nv2 = size(v2_stations, 2)
    _, nhours, nwavs = size(N_interf)
    σ_ron  = combiner.read_noise
    n_pix  = combiner.n_pix_fringe
    n_pixp = combiner.n_pix_photometry
    vinst  = combiner.instrument_visibility

    sq_coef   = zeros(nv2, nhours, nwavs)
    var_coef  = zeros(nv2, nhours, nwavs)
    var_const = zeros(nv2, nhours, nwavs)
    v2phot    = zeros(nv2, nhours, nwavs)

    ron2 = n_pix * σ_ron^2
    for w in 1:nwavs, h in 1:nhours
        N_tot = zero(eltype(N_interf))
        for t in 1:ntel
            N_tot += N_interf[t, h, w]
        end
        vc = 2.0 * (N_tot + ron2) + 4.0
        vk = N_tot^2 + N_tot * (1.0 + 2.0 * ron2) + (3.0 + n_pix) * n_pix * σ_ron^4
        for b in 1:nv2
            i, j = v2_stations[1, b], v2_stations[2, b]
            sq_coef[b, h, w]   = N_interf[i, h, w] * N_interf[j, h, w] * vinst^2
            var_coef[b, h, w]  = vc
            var_const[b, h, w] = vk

            if combiner.flux_frac_photometry > 0
                Np = 0.5 * (N_phot[i, h, w] + N_phot[j, h, w])
                if Np > 0
                    σp = sqrt(Np + n_pixp * σ_ron^2)
                    v2phot[b, h, w] = 2.0 * (σp / Np)^2
                end
            end
        end
    end
    return sq_coef, var_coef, var_const, v2phot
end

"""
    complex_vis_error(visamp, sq_coef, var_coef, var_const, v2phot, N_frames)

Per-component error `σ_c` of the complex visibility, i.e. the standard deviation of the real
and of the imaginary part of a circular complex Gaussian perturbation.

This, and not `σ(V²)`, is the primitive quantity of the noise model. Every observable is
derived from one perturbed complex visibility, so every error bar has to be derived from the
same `σ_c` or the written uncertainties will not describe the scatter actually present.

`σ_c` is fixed by matching the two limits of ASPRO's correlated-flux variance. For a circular
complex Gaussian, `Var(|V+n|²) = 4V²σ_c² + 4σ_c⁴`, so

- as `V → 0`, `σ(V²) → 2σ_c²`, which must equal the additive floor
  `sqrt(var_const)/(sq_coef·sqrt(N_frames))`
- for large `V`, `σ(V²) → 2Vσ_c`, which must equal `V·sqrt(var_coef/(sq_coef·N_frames))`

giving the two terms below, plus a photometric-normalisation term that scales with `|V|`
because a photometric error is multiplicative on `V²`.

Deriving `σ_c` the other way round — as `σ(V²)/(2|V|)` — diverges on a resolved baseline and
injects noise far larger than the error bars claim.
"""
function complex_vis_error(visamp, sq_coef, var_coef, var_const, v2phot, N_frames)
    sq_coef <= 0 && return 10.0
    σ2  = 0.25 * var_coef / (sq_coef * N_frames)                    # multiplicative / photon
    σ2 += 0.5 * sqrt(max(var_const, 0.0)) / (sq_coef * sqrt(N_frames))  # additive floor
    σ2 += 0.25 * v2phot * visamp^2 / N_frames                       # photometric channels
    return min(sqrt(max(σ2, 0.0)), 10.0)
end

"""
    vis2_error_from_sigma(v2, σ_c)

`σ(V²)` implied by a circular complex Gaussian of per-component width `σ_c`:
`sqrt(4V²σ_c² + 4σ_c⁴) = 2σ_c·sqrt(V² + σ_c²)`.

Using this rather than an independent expression is what makes the written `VIS2ERR` match
the scatter that the sampler actually produces.
"""
vis2_error_from_sigma(v2, σ_c) = 2 * σ_c * sqrt(max(abs(v2), 0.0) + σ_c^2)


#CODES FOR SIMULATING OIFITS BASED ON INPUT IMAGE AND EITHER INPUT OIFITS OR HOUR ANGLES

"""
    vis_to_t3_closure(cvis, indx1, indx2, indx3)

Bispectrum and closure phase for a triangle whose three legs are indexed as `(i,j)`, `(j,k)`
and `(i,k)` — the convention `get_t3_baselines` produces.

The closure phase is `φ_ij + φ_jk + φ_ki`, and the third leg stored is `(i,k) = -(k,i)`, so
the third visibility has to be **conjugated**: `V_ij · V_jk · conj(V_ik)`.

This differs from `vis_to_t3` in oichi2.jl, which takes no conjugate because `readoifits`
gives it indices pointing at the already-negated third leg `u3 = -(u1+u2)`. Both are correct
for their own index convention; using the wrong one flips the sign of the closure phase while
leaving V2 and T3AMP untouched, which is exactly why it can go unnoticed.
"""
function vis_to_t3_closure(cvis::AbstractVector{Complex{T}}, indx1, indx2, indx3) where T<:AbstractFloat
    t3 = cvis[indx1] .* cvis[indx2] .* conj.(cvis[indx3])
    return t3, abs.(t3), angle.(t3) .* T(180/pi)
end


function get_v2_baselines(N,station_xyz,tel_names)
    # determine V2 Baselines and make necessary arrays
    nv2 = Int64(N*(N-1)/2);
    v2_baselines = Array{Float64}(undef,3,nv2);
    v2_stations  = Array{Int64}(undef,2,nv2);
    #v2_stations_nonredun=Array{Int64}(undef,2,nv2);
    v2_indx      = Array{Int64}(undef,nv2);
    baseline_name = Array{String}(undef,nv2);
    ind = 1
    for i=1:N-1
      for j=i+1:N
            v2_baselines[:,ind] .= station_xyz[j,:]-station_xyz[i,:];
            v2_stations[:,ind] = [i,j];
            v2_indx[ind] = ind;
            ind += 1
      end
    end
    for i=1:nv2
        baseline_name[i]=string(tel_names[v2_stations[1,i]],"-",tel_names[v2_stations[2,i]])
    end
    return nv2,v2_baselines,v2_stations,v2_indx,baseline_name
end

function v2mapt3(i,j,v2_stations)
    # really dumb way of doing this
    # there is probably an analytic formula for this
    listi = findall(v2_stations[1,:].==i);
    return listi[findfirst(v2_stations[2,listi].==j)]
end


function get_t3_baselines(N,station_xyz,v2_stations)
    # determine T3 Baselines and make neccessary
    nt3 = binomial(N,3);# V2 Baselines
    t3_baselines=Array{Float64}(undef,3,3,nt3);
    t3_stations  = Array{Int64}(undef,3, nt3);
    t3_indx_1  = Array{Int64}(undef,nt3);
    t3_indx_2  = Array{Int64}(undef,nt3);
    t3_indx_3  = Array{Int64}(undef,nt3);
    ind = 1;
    for i=1:N
      for j=i+1:N-1
        for k=j+1:N
            t3_baselines[:, 1, ind] .= station_xyz[j,:]-station_xyz[i,:];
            t3_baselines[:, 2, ind] .= station_xyz[k,:]-station_xyz[j,:];
            t3_baselines[:, 3, ind] .= station_xyz[i,:]-station_xyz[k,:];
            t3_stations[:,ind]=[i,j,k];
            t3_indx_1[ind] = v2mapt3(i,j,v2_stations);
            t3_indx_2[ind] = v2mapt3(j,k,v2_stations);
            t3_indx_3[ind] = v2mapt3(i,k,v2_stations);
            ind += 1
        end
      end
    end
    return nt3,t3_baselines,t3_stations,t3_indx_1,t3_indx_2,t3_indx_3
end

function get_uv(l, h, λ, δ::Float64, baselines)
    # Use following expression only if there are missing baselines somewhere
    # baselines = hcat(v2_baselines, t3_baselines[:, 1, :], t3_baselines[:, 2, :], t3_baselines[:, 3, :])
    # Expression to use for pure simulation where the full complement of v2 and t3 are created
    nhours = length(h);
    nuv = size(baselines, 2);
    # Now compute the UV coordinates, again according to the APIS++ standards.
    u_M = -sin.(l)*sin.(h) .*baselines[1,:]  .+ cos.(h) .* baselines[2,:]+cos.(l)*sin.(h).*baselines[3,:];
    v_M = (sin.(l)*sin(δ)*cos.(h).+cos.(l)*cos(δ)).* baselines[1,:]  + sin(δ)*sin.(h) .*baselines[2,:]   +(-cos.(l)*cos.(h)*sin(δ).+sin.(l)*cos(δ)) .* baselines[3,:];
    w_M = (-sin.(l)*cos(δ)*cos.(h).+cos.(l)*sin(δ)).* baselines[1,:] - cos(δ)*sin.(h) .*baselines[2,:]   + (cos.(l)*cos.(h)*cos(δ).+sin.(l)sin(δ))  .* baselines[3,:];
    # proj baselines to (uv wav)
    u = reshape((1 ./λ)'.*vec(u_M), (nuv,nhours,length(λ)));
    v = reshape((1 ./λ)'.*vec(v_M), (nuv,nhours,length(λ)));
    w = reshape((1 ./λ)'.*vec(w_M), (nuv,nhours,length(λ)));
    uv = vcat(vec(u)',vec(v)')
    return nuv,uv,u_M,v_M,w_M
end

function get_uv(l, h, λ, δ::Matrix{Float64}, baselines)
    # Use following expression only if there are missing baselines somewhere
    # baselines = hcat(v2_baselines, t3_baselines[:, 1, :], t3_baselines[:, 2, :], t3_baselines[:, 3, :])
    # Expression to use for pure simulation where the full complement of v2 and t3 are created
    nhours = length(h);
    nuv = size(baselines, 2);
    # Now compute the UV coordinates, again according to the APIS++ standards.
    # If δ(=DECLINATION) is a vector, this is because DEC changes with HA, i.e. it's not a star but e.g. a satellite
    u_M = -sin(l)*sin.(h) .*baselines[1,:]  .+ cos.(h) .* baselines[2,:]+cos(l)*sin.(h).*baselines[3,:];
    v_M = ( sin(l)*sin.(δ).*cos.(h).+cos(l).*cos.(δ)).* baselines[1,:]  + sin.(δ).*sin.(h) .*baselines[2,:]   +(-cos(l)*cos.(h).*sin.(δ).+sin(l).*cos.(δ)) .* baselines[3,:];
    w_M = (-sin(l)*cos.(δ).*cos.(h).+cos(l).*sin.(δ)).* baselines[1,:] - cos.(δ).*sin.(h) .*baselines[2,:]   + (cos(l)*cos.(h).*cos.(δ).+sin(l).*sin.(δ)) .* baselines[3,:];
    # proj baselines to (uv wav)
    u = reshape((1 ./λ)'.*vec(u_M), (nuv,nhours,length(λ)));
    v = reshape((1 ./λ)'.*vec(v_M), (nuv,nhours,length(λ)));
    w = reshape((1 ./λ)'.*vec(w_M), (nuv,nhours,length(λ)));
    uv = vcat(vec(u)',vec(v)')
    return nuv,uv,u_M,v_M,w_M
end







function get_uv_indxes(nhours,nuv,nv2,nt3,v2_indx,t3_indx_1,t3_indx_2,t3_indx_3,nw,uv)
    v2_indx_M = repeat(v2_indx,1,nhours)+repeat(Int64.(collect(range(0,stop=nuv*(nhours-1),length=nhours)))',nv2)
    t3_indx_1_M = repeat(t3_indx_1,1,nhours)+repeat(Int64.(collect(range(0,stop=nuv*(nhours-1),length=nhours)))',nt3)
    t3_indx_2_M = repeat(t3_indx_2,1,nhours)+repeat(Int64.(collect(range(0,stop=nuv*(nhours-1),length=nhours)))',nt3)
    t3_indx_3_M = repeat(t3_indx_3,1,nhours)+repeat(Int64.(collect(range(0,stop=nuv*(nhours-1),length=nhours)))',nt3)
    #indexes separated by wavelength channel
    v2_indx_w = vec(repeat(vec(v2_indx_M),1,nw)+repeat(Int64.(collect(range(0,stop=nuv*nhours*(nw-1),length=nw)))',nv2*nhours))
    t3_indx_1_w = vec(repeat(vec(t3_indx_1_M),1,nw)+repeat(Int64.(collect(range(0,stop=nuv*nhours*(nw-1),length=nw)))',nt3*nhours))
    t3_indx_2_w = vec(repeat(vec(t3_indx_2_M),1,nw)+repeat(Int64.(collect(range(0,stop=nuv*nhours*(nw-1),length=nw)))',nt3*nhours))
    t3_indx_3_w = vec(repeat(vec(t3_indx_3_M),1,nw)+repeat(Int64.(collect(range(0,stop=nuv*nhours*(nw-1),length=nw)))',nt3*nhours))
    t3_baseline = (sqrt.(uv[1,t3_indx_1_w].^2 + uv[2,t3_indx_1_w].^2).*sqrt.(uv[1,t3_indx_2_w].^2 + uv[2,t3_indx_2_w].^2).*sqrt.(uv[1,t3_indx_3_w].^2 + uv[2,t3_indx_3_w].^2)).^(1 ./ 3.);
    return v2_indx_M,t3_indx_1_M,t3_indx_2_M,t3_indx_3_M,v2_indx_w,t3_indx_1_w,t3_indx_2_w,t3_indx_3_w
 end

function simulate(facility,target,combiner,wavelength,dates,out_file; image::Union{String, Array{Float64,1}, Array{Float64,2}, Array{Float64, 3}, Array{Float64,4}}="", pixsize::Float64=0.1, flat_model::Union{FlatModel,Nothing}=nothing, flat_params::Vector{Float64}=Float64[], dft=false,
                  noise::Bool=true, nonoise::Union{Nothing,Bool}=nothing, debias::Bool=true,
                  mag=2.0, mag_ao::Union{Nothing,Real}=nothing,
                  n_samples::Int=100, rng::Union{Nothing,AbstractRNG}=nothing, seed::Union{Nothing,Integer}=nothing,
                  observability=nothing)
    #simulate an observation using input hour angles, info about array and combiner, and input image
    ntel=facility.ntel;

    # `nonoise=true` is the old spelling of `noise=false`.
    if nonoise !== nothing
        Base.depwarn("simulate(...; nonoise=$(nonoise)) is deprecated, use noise=$(!nonoise) instead", :simulate)
        noise = !nonoise
    end
    _rng = rng !== nothing ? rng : (seed !== nothing ? MersenneTwister(seed) : Random.default_rng())

    # Observability filtering is strictly opt-in: by default `simulate` is a pure uv-coverage
    # simulator and every epoch you pass in is used, observable or not.
    if observability !== nothing
        sel = observable_epochs(facility, target, collect(dates); observability...)
        r = sel.report
        if r.n_out == 0
            error("observability filter removed all $(r.n_in) epochs " *
                  "($(r.n_dropped_elevation) below/above elevation limits, $(r.n_dropped_delay) outside delay lines)")
        end
        if r.n_out < r.n_in
            @info "observability: keeping $(r.n_out)/$(r.n_in) epochs " *
                  "(dropped $(r.n_dropped_elevation) on elevation, $(r.n_dropped_delay) on delay lines)"
        end
        dates = sel.dates
    end

    # RA and DEC are in (decimal) degrees, per the OIFITS standard for OI_TARGET.
    dec = target.decep0/180*pi
    ra = target.raep0/180*pi

    # hour_angle_calc wants RA in hours, hence the /15.
    lst, hour_angles = hour_angle_calc(dates, facility.lon, target.raep0/15);
    nhours = length(hour_angles);
    # Altitude at each epoch drives the airmass-dependent Strehl and photon budget.
    elevation_deg, _ = alt_az(target.decep0, facility.lat, hour_angles)
    h_rad = hour_angles' .* pi / 12;
    l = facility.lat/180*pi;
    λ = wavelength.λ;
    δλ = wavelength.δλ;
    nwavs = length(λ)

    station_xyz = facility.sta_xyz  # (ntel, 3) matrix

    # Find physical baselines and triangles combinations
    nv2,v2_baselines,v2_stations,v2_indx,baseline_name         = get_v2_baselines(ntel,station_xyz,facility.tel_names);
    nt3,t3_baselines,t3_stations,t3_indx_1,t3_indx_2,t3_indx_3 = get_t3_baselines(ntel,station_xyz,v2_stations);

    # Compute uv coverage: nuv is the number of uv points, uv is (u,v) in Mλ, (u_M, v_M, w_M) are in meters
    nuv, uv, u_M, v_M, w_M = get_uv(l,h_rad,λ,dec,v2_baselines)

    # Setup indexing for OIFITS:  _M -> in meters, _w -> scaled by λ
    v2_indx_M,t3_indx_1_M,t3_indx_2_M,t3_indx_3_M,v2_indx_w,t3_indx_1_w,t3_indx_2_w,t3_indx_3_w = get_uv_indxes(nhours,nuv,nv2,nt3,v2_indx,t3_indx_1,t3_indx_2,t3_indx_3,nwavs,uv)

    # Compute complex visibilities from either a truth image or a model
    cvis_model = ComplexF64[]
    # Per-channel relative flux of the source, normalised to a mean of 1. Multiplies the
    # photon count per channel and fills OI_FLUX. Stays at 1 unless an image cube supplies a
    # spectrum: normalising each channel of a cube (which is required to get V right) would
    # otherwise throw the source SED away entirely.
    spec_weight = ones(nwavs)
    # Determine if image or model

    if (((typeof(image) == String) && (image !="")) || ((typeof(image) != String) && (image !="")) )
        # TODO: Update polychromatic and dynamic imaging and non square images
        x = image isa AbstractString ? readfits(image) : Float64.(image)
        nx = size(x,1); 
        if ndims(x) == 2 # Monochromatic image
            println("Simulating monochromatic observations...")
            x = x/sum(x);
            ft = setup_nfft(uv, v2_indx_w, v2_indx_w, t3_indx_1_w, t3_indx_2_w, t3_indx_3_w, nx, pixsize);
            cvis_model = image_to_vis(x, ft);
        elseif ndims(x) == 3 # at the moment, we assume polychromatic
            println("Simulating polychromatic observations...")
            # Capture the spectrum before normalising each channel to unit total flux.
            chan_flux = vec(sum(x, dims=(1,2)))
            size(x,3) == nwavs || throw(DimensionMismatch(
                "image cube has $(size(x,3)) planes but the wavelength table has $nwavs channels"))
            mean_flux = sum(chan_flux) / nwavs
            spec_weight = mean_flux > 0 ? chan_flux ./ mean_flux : ones(nwavs)
            x = x ./ sum(x, dims=(1,2))
            # uv_lambda = reshape(uv, 2*nuv*nhours,nw) 
            # k=8; range=(k-1)*nuv*nhours+1:k*nuv*nhours; norm(uv_lambda[:,:,k]-uv[:,range]); 
            # We're basically redoing setup_nfft_polychromatic here
            #      ft = Array{Array{NFFT.NFFTPlan{Float64, 2, 1}}}(undef, nwavs);
            cvis = fill((Complex{Float64}[]),nwavs);
            if dft==true
                println("Using DFT")
                #ft = Array{Array{ComplexF64,2}}(undef, nwavs);
                for i=1:nwavs
                    ft = setup_dft(uv[:,(i-1)*nuv*nhours+1:i*nuv*nhours], nx, pixsize);
                    cvis[i] = image_to_vis(x[:,:,i], ft);
                end
            else
                println("Using NFFT")
                for i=1:nwavs
                    #ft[i]=setup_nfft(uv[:,(i-1)*nuv*nhours+1:i*nuv*nhours], vec(v2_indx_M), vec(v2_indx_M), vec(t3_indx_1_M), vec(t3_indx_2_M), vec(t3_indx_3_M), nx, pixsize);
                    ft = setup_nfft(uv[:,(i-1)*nuv*nhours+1:i*nuv*nhours], vec(v2_indx_M), vec(v2_indx_M), vec(t3_indx_1_M), vec(t3_indx_2_M), vec(t3_indx_3_M), nx, pixsize);
                    cvis[i] = image_to_vis(x[:,:,i], ft);
                end
            end
            cvis_model = vec(hcat(cvis...))
        end
    elseif flat_model !== nothing
        # Flat-dict parametric model
        # Build per-UV-point wavelength and MJD vectors matching the uv ordering
        # uv is (nuv, nhours, nwavs) flattened: baseline fastest, wavelength slowest
        wl_per_uv  = repeat(λ, inner=nuv*nhours)                          # (nuv*nhours*nwavs,)
        mjd_per_uv = repeat(repeat(datetime_to_mjd.(dates), inner=nuv), outer=nwavs)
        cvis_model = eval_model(flat_model, flat_params, uv; wl=wl_per_uv, mjd=mjd_per_uv)
    else
        @warn("No image nor model definition in call to simulate()")
        @warn("Will generate zero visibilities")
        cvis_model = zeros(ComplexF64, size(uv,2))
    end

    # ── Photon budget ────────────────────────────────────────────────────────
    mags   = resolve_magnitudes(mag, λ)
    # The AO wavefront sensor guides on the science target unless told otherwise, and it
    # senses in its own band, not the science band.
    _ao_mag = if mag_ao !== nothing
        Float64(mag_ao)
    elseif mag isa AbstractDict && facility.ao !== nothing && haskey(mag, facility.ao.band)
        Float64(mag[facility.ao.band])
    elseif facility.ao !== nothing
        b = band_by_name(facility.ao.band)
        resolve_magnitudes(mag, [b === nothing ? 0.55e-6 : b.lambda])[1]
    else
        mags[1]
    end

    N_interf, DIT, N_frames = photons_per_telescope(mags, λ, δλ, facility, combiner,
                                                    elevation_deg, _ao_mag;
                                                    channel=:fringes, spec_weight=spec_weight)
    N_phot, _, _            = photons_per_telescope(mags, λ, δλ, facility, combiner,
                                                    elevation_deg, _ao_mag;
                                                    channel=:photometry, spec_weight=spec_weight)

    # Error model coefficients per (baseline, epoch, channel). Flatten to the ordering
    # simulate uses everywhere: baseline fastest, then epoch, then wavelength.
    sq_coef3, var_coef3, var_const3, v2phot3 =
        correlated_flux_coefficients(N_interf, N_phot, combiner, ntel, v2_stations)
    sq_coef   = vec(sq_coef3)
    var_coef  = vec(var_coef3)
    var_const = vec(var_const3)
    v2phot    = vec(v2phot3)

    # ── Error bars from the noiseless model ──────────────────────────────────
    v2_true     = vis_to_v2(cvis_model, v2_indx_w)
    visamp_true = abs.(cvis_model[v2_indx_w])

    # The complex-visibility error is the primitive: the noise is drawn from it and every
    # error bar below is derived from it, so the written uncertainties describe the scatter.
    σ_cvis = complex_vis_error.(visamp_true, sq_coef, var_coef, var_const, v2phot, N_frames)
    # Scatter into cvis index space. v2_indx_w happens to be the identity here (nuv == nv2),
    # but do not rely on that: the triangle indices below index cvis directly.
    σ_cvis_full = zeros(length(cvis_model))
    σ_cvis_full[v2_indx_w] .= σ_cvis

    v2_stat = vis2_error_from_sigma.(v2_true, σ_cvis)
    v2_model_err = _with_systematics.(v2_stat, 2 * combiner.vis_cal_err .* abs.(v2_true))
    v2_model_err = min.(v2_model_err, 100.0)   # cap runaway errors, as ASPRO does

    visamp_model_err = _with_systematics.(σ_cvis, combiner.vis_cal_err .* visamp_true)
    # Phase error of a complex visibility with amplitude |V| and per-component error sigma.
    visphi_stat_deg  = (180.0/π) .* σ_cvis ./ max.(visamp_true, 1e-6)
    dphi_model_err   = _with_systematics.(min.(visphi_stat_deg, 180.0), combiner.phase_cal_err)

    # Closure quantities: add the three baselines' relative errors in quadrature.
    rel1 = σ_cvis_full[t3_indx_1_w] ./ max.(abs.(cvis_model[t3_indx_1_w]), 1e-6)
    rel2 = σ_cvis_full[t3_indx_2_w] ./ max.(abs.(cvis_model[t3_indx_2_w]), 1e-6)
    rel3 = σ_cvis_full[t3_indx_3_w] ./ max.(abs.(cvis_model[t3_indx_3_w]), 1e-6)
    rel_t3 = sqrt.(rel1.^2 .+ rel2.^2 .+ rel3.^2)

    _, t3amp_true, _ = vis_to_t3_closure(cvis_model, t3_indx_1_w, t3_indx_2_w, t3_indx_3_w)
    t3amp_model_err = _with_systematics.(abs.(t3amp_true) .* rel_t3,
                                         3 * combiner.vis_cal_err .* abs.(t3amp_true))
    t3phi_model_err = _with_systematics.(min.((180.0/π) .* rel_t3, 180.0), combiner.phase_cal_err)

    # OI_FLUX: per-telescope spectrophotometry, normalised so the mean channel is 1.
    ron2_p       = combiner.n_pix_photometry * combiner.read_noise^2
    snr_tel_3d   = sqrt.(N_frames) .* N_phot ./ sqrt.(max.(N_phot .+ ron2_p, 1.0))
    nobs_flux    = ntel * nhours
    flux_true_3d = repeat(reshape(spec_weight, 1, 1, nwavs), ntel, nhours, 1)
    flux_true    = vec(flux_true_3d)
    flux_model_err = vec(flux_true_3d ./ max.(snr_tel_3d, 1e-3))

    # ── Draw the noise ONCE, on the complex visibility ───────────────────────
    # Every observable is then derived from the same perturbed cvis, so VISAMP^2 == VIS2 and
    # T3PHI is exactly the sum of the three baseline phases. Drawing each observable
    # independently, as this code used to, breaks both relations.
    cvis_noisy = cvis_model
    flux_model = copy(flux_true)
    if noise
        cvis_noisy = cvis_model .+ σ_cvis_full .* (randn(_rng, length(cvis_model)) .+
                                            im .* randn(_rng, length(cvis_model)))
        flux_model .+= flux_model_err .* randn(_rng, length(flux_model))
    end

    v2_model = vis_to_v2(cvis_noisy, v2_indx_w)
    # A squared visibility built from a noisy complex visibility is biased high by 2*sigma_c^2
    # (E|V+n|^2 = V^2 + 2 sigma_c^2). Real pipelines subtract this, so by default so do we:
    # otherwise every simulated V2 sits ~0.5 sigma above the truth and biases anything fitted
    # to it. Pass debias=false to keep the raw, biased estimator.
    if noise && debias
        v2_model = v2_model .- 2 .* σ_cvis.^2
    end
    _, t3amp_model, t3phi_model = vis_to_t3_closure(cvis_noisy, t3_indx_1_w, t3_indx_2_w, t3_indx_3_w)
    # Debias the amplitude with the same correction so that VISAMP^2 == VIS2 still holds
    # exactly: the two must stay the same measurement expressed two ways.
    visamp_model = (noise && debias) ? sqrt.(max.(v2_model, 0.0)) : abs.(cvis_noisy[v2_indx_w])
    visphi_abs   = angle.(cvis_noisy[v2_indx_w]) .* (180.0 / π)

    # Differential phase: subtract the mean phase over the other channels of the same
    # baseline and epoch. Excluding the channel itself is what makes this a differential
    # measurement rather than a slightly-shrunk absolute one.
    phi_3d = reshape(visphi_abs, nv2, nhours, nwavs)
    dphi_model = if nwavs > 1
        tot = sum(phi_3d, dims=3)
        vec(phi_3d .- (tot .- phi_3d) ./ (nwavs - 1))
    else
        zeros(length(visphi_abs))
    end

    # Wrap phases into (-180, 180].
    _wrap180(x) = x - 360.0 * floor((x + 180.0) / 360.0)
    t3phi_model = _wrap180.(t3phi_model)
    dphi_model  = _wrap180.(dphi_model)

    # Replace the analytic T3 error bars with sampled ones.
    #
    # The analytic closure errors are a small-error expansion: they hold while sigma/|V| << 1
    # and fail badly outside it. Measured against the realised scatter on a resolved disc:
    #
    #     median SNR(V2)      6.3     1.0     0.09    0.02
    #     T3AMP z-std analytic  1.01    1.90   25.9    152
    #     T3AMP z-std sampled   1.08    1.09    1.08    1.08
    #
    # i.e. two orders of magnitude wrong exactly in the faint regime worth simulating, for no
    # measurable extra time. Hence sampling is the default; pass n_samples=0 for the analytic
    # form (cheaper, and adequate when every baseline is well detected).
    #
    # The spread is taken about the noiseless truth rather than about the sample mean, so the
    # error bar covers T3AMP's bias as well as its scatter. That is what makes a chi2 against
    # the true model come out at 1.
    if n_samples > 0
        _, t3amp_tr, t3phi_tr = vis_to_t3_closure(cvis_model, t3_indx_1_w, t3_indx_2_w, t3_indx_3_w)
        t3amp_s2 = zeros(length(t3amp_model_err))
        t3phi_s2 = zeros(length(t3phi_model_err))
        for _ in 1:n_samples
            c = cvis_model .+ σ_cvis_full .* (randn(_rng, length(cvis_model)) .+
                                       im .* randn(_rng, length(cvis_model)))
            _, a, p = vis_to_t3_closure(c, t3_indx_1_w, t3_indx_2_w, t3_indx_3_w)
            t3amp_s2 .+= (a .- t3amp_tr).^2
            t3phi_s2 .+= _wrap180.(p .- t3phi_tr).^2
        end
        t3amp_model_err = sqrt.(t3amp_s2 ./ n_samples)
        t3phi_model_err = sqrt.(t3phi_s2 ./ n_samples)
        t3amp_model_err = _with_systematics.(t3amp_model_err, 3 * combiner.vis_cal_err .* abs.(t3amp_true))
        t3phi_model_err = _with_systematics.(min.(t3phi_model_err, 180.0), combiner.phase_cal_err)
    end

    # Flag points whose error bars make them meaningless rather than writing garbage.
    v2_flag_flat     = (.!isfinite.(v2_model))     .| (v2_model_err     .> 10.0)
    t3amp_flag_flat  = (.!isfinite.(t3amp_model))  .| (t3amp_model_err  .> 10.0)
    t3phi_flag_flat  = (.!isfinite.(t3phi_model))  .| (t3phi_model_err  .> 180.0)
    t3_flag_flat     = t3amp_flag_flat .| t3phi_flag_flat
    vis_flag_flat    = (.!isfinite.(visamp_model)) .| (visamp_model_err .> 10.0)
    flux_flag_flat   = (.!isfinite.(flux_model))   .| (flux_model_err   .> 10.0)

    # --- Build OIFITS data structures and save ---
    nobs_v2  = nv2 * nhours
    nobs_t3  = nt3 * nhours
    mjd_vals = datetime_to_mjd.(dates)
    date_obs = string(Dates.Date(dates[1]))
    tid      = target.target_id

    # Reshape flat vectors to (nwavs, nobs) matrices for OIFITS tables.
    # simulate ordering is (element, hour, wavelength-slowest), reshape via
    # (n1, nhours, nwavs) → (nwavs, n1*nhours) by permuting dims.
    function to_oifits_matrix(flat, n1, nh, nw)
        return reshape(permutedims(reshape(flat, n1, nh, nw), (3, 1, 2)), nw, n1*nh)
    end

    v2_matrix       = to_oifits_matrix(v2_model,        nv2, nhours, nwavs)
    v2_err_matrix   = to_oifits_matrix(v2_model_err,    nv2, nhours, nwavs)
    t3amp_matrix    = to_oifits_matrix(t3amp_model,     nt3, nhours, nwavs)
    t3amperr_matrix = to_oifits_matrix(t3amp_model_err, nt3, nhours, nwavs)
    t3phi_matrix    = to_oifits_matrix(t3phi_model,     nt3, nhours, nwavs)
    t3phierr_matrix = to_oifits_matrix(t3phi_model_err, nt3, nhours, nwavs)
    visamp_matrix   = to_oifits_matrix(visamp_model,     nv2, nhours, nwavs)
    visamperr_matrix= to_oifits_matrix(visamp_model_err, nv2, nhours, nwavs)
    dphi_matrix     = to_oifits_matrix(dphi_model,       nv2, nhours, nwavs)
    dphierr_matrix  = to_oifits_matrix(dphi_model_err,   nv2, nhours, nwavs)
    flux_matrix     = to_oifits_matrix(flux_model,       ntel, nhours, nwavs)
    fluxerr_matrix  = to_oifits_matrix(flux_model_err,   ntel, nhours, nwavs)
    v2_flag_matrix   = Matrix{Bool}(to_oifits_matrix(v2_flag_flat,   nv2,  nhours, nwavs))
    t3_flag_matrix   = Matrix{Bool}(to_oifits_matrix(t3_flag_flat,   nt3,  nhours, nwavs))
    vis_flag_matrix  = Matrix{Bool}(to_oifits_matrix(vis_flag_flat,  nv2,  nhours, nwavs))
    flux_flag_matrix = Matrix{Bool}(to_oifits_matrix(flux_flag_flat, ntel, nhours, nwavs))

    # Per-observation vectors (meter-baseline coords, MJD, station indices)
    ucoord_vis2 = vec(u_M[v2_indx_M])
    vcoord_vis2 = vec(v_M[v2_indx_M])
    u1coord     = vec(u_M[t3_indx_1_M])
    v1coord     = vec(v_M[t3_indx_1_M])
    u2coord     = vec(u_M[t3_indx_2_M])
    v2coord     = vec(v_M[t3_indx_2_M])
    mjd_vis2    = repeat(mjd_vals, inner=nv2)
    mjd_t3      = repeat(mjd_vals, inner=nt3)
    mjd_flux    = repeat(mjd_vals, inner=ntel)
    sta_idx_arr = facility.sta_index

    # OI_TARGET
    tgt = OIFITS.OI_TARGET([OIFITS.OITargetEntry(
        tid, target.target,
        target.raep0,  target.decep0,
        target.equinox,
        target.ra_err, target.dec_err,
        target.sysvel, target.veltyp,
        target.veldef,
        target.pmra,   target.pmdec,
        target.pmra_err, target.pmdec_err,
        target.parallax, target.para_err,
        target.spectyp,  "SCI")]; revn=1)

    # OI_ARRAY
    arr = OIFITS.OI_ARRAY(undef)
    arr.revn    = 1
    arr.arrname = facility.name
    arr.frame   = "GEOCENTRIC"
    arr.arrayx  = 0.0; arr.arrayy = 0.0; arr.arrayz = 0.0
    arr.tel_name  = facility.tel_names
    arr.sta_name  = facility.sta_names
    arr.sta_index = sta_idx_arr
    arr.diameter  = facility.tel_diams
    arr.staxyz    = collect(facility.sta_xyz')  # (ntel,3)' → (3,ntel)

    # OI_WAVELENGTH
    ins = OIFITS.OI_WAVELENGTH(undef)
    ins.revn     = 1
    ins.insname  = string(wavelength.combiner, "_", wavelength.mode)
    ins.eff_wave = Float64.(λ)
    ins.eff_band = Float64.(δλ)

    # OI_VIS2
    db_vis2 = OIFITS.OI_VIS2(undef)
    db_vis2.revn      = 2
    db_vis2.date_obs  = date_obs
    db_vis2.arrname   = arr.arrname
    db_vis2.insname   = ins.insname
    db_vis2.target_id = fill(tid, nobs_v2)
    db_vis2.time      = zeros(nobs_v2)
    db_vis2.mjd       = mjd_vis2
    db_vis2.int_time  = ones(nobs_v2)
    db_vis2.vis2data  = v2_matrix
    db_vis2.vis2err   = v2_err_matrix
    db_vis2.ucoord    = ucoord_vis2
    db_vis2.vcoord    = vcoord_vis2
    db_vis2.sta_index = repeat(v2_stations, 1, nhours)
    db_vis2.flag      = v2_flag_matrix

    # OI_VIS (differential visibility amplitude and phase)
    db_vis = OIFITS.OI_VIS(undef)
    db_vis.revn      = 2
    db_vis.date_obs  = date_obs
    db_vis.arrname   = arr.arrname
    db_vis.insname   = ins.insname
    db_vis.amptyp    = "ABSOLUTE"
    db_vis.phityp    = "DIFFERENTIAL"
    db_vis.target_id = fill(tid, nobs_v2)
    db_vis.time      = zeros(nobs_v2)
    db_vis.mjd       = mjd_vis2
    db_vis.int_time  = ones(nobs_v2)
    db_vis.visamp    = visamp_matrix
    db_vis.visamperr = visamperr_matrix
    db_vis.visphi    = dphi_matrix
    db_vis.visphierr = dphierr_matrix
    db_vis.ucoord    = ucoord_vis2
    db_vis.vcoord    = vcoord_vis2
    db_vis.sta_index = repeat(v2_stations, 1, nhours)
    db_vis.flag      = vis_flag_matrix

    # OI_T3
    db_t3 = OIFITS.OI_T3(undef)
    db_t3.revn      = 2
    db_t3.date_obs  = date_obs
    db_t3.arrname   = arr.arrname
    db_t3.insname   = ins.insname
    db_t3.target_id = fill(tid, nobs_t3)
    db_t3.time      = zeros(nobs_t3)
    db_t3.mjd       = mjd_t3
    db_t3.int_time  = ones(nobs_t3)
    db_t3.t3amp     = t3amp_matrix
    db_t3.t3amperr  = t3amperr_matrix
    db_t3.t3phi     = t3phi_matrix
    db_t3.t3phierr  = t3phierr_matrix
    db_t3.u1coord   = u1coord
    db_t3.v1coord   = v1coord
    db_t3.u2coord   = u2coord
    db_t3.v2coord   = v2coord
    db_t3.sta_index = repeat(t3_stations, 1, nhours)
    db_t3.flag      = t3_flag_matrix

    # OI_FLUX (per-telescope spectrophotometry)
    db_flux = OIFITS.OI_FLUX(undef)
    db_flux.revn      = 1
    db_flux.date_obs  = date_obs
    db_flux.arrname   = arr.arrname
    db_flux.insname   = ins.insname
    db_flux.calstat   = "C"
    db_flux.fov       = 0.0
    db_flux.fovtype   = "FWHM"
    db_flux.target_id = fill(tid, nobs_flux)
    db_flux.mjd       = mjd_flux
    db_flux.int_time  = ones(nobs_flux)
    db_flux.fluxdata  = flux_matrix
    db_flux.fluxerr   = fluxerr_matrix
    db_flux.sta_index = repeat(sta_idx_arr, nhours)
    db_flux.flag      = flux_flag_matrix

    # Assemble and write
    ds = OIFITS.OIDataSet(tgt, arr, ins, db_vis2, db_vis, db_t3, db_flux)
    OIFITS.write(out_file, ds; overwrite=true)
end

function simulate_from_oifits(in_oifits, out_file; mode="copy_errors", errors=[],
        image::Union{String, Array{Float64,1}, Array{Float64,2}, Array{Float64,3}, Array{Float64,4}}="",
        pixsize::Float64=0.1, flat_model::Union{FlatModel,Nothing}=nothing, flat_params::Vector{Float64}=Float64[])

    # Read input OIFITS data for UV coordinates and errors
    data = readoifits(in_oifits, filter_bad_data=false)[1]

    # Compute complex visibilities from truth image or model
    if (((typeof(image) == String) && (image != "")) || ((typeof(image) != String) && (image != "")))
        x = image isa AbstractString ? readfits(image) : Float64.(image)
        nx = size(x, 1)
        x = vec(x) / sum(x)
        ft = setup_nfft(data, nx, pixsize)
        cvis_model = image_to_vis(to_ft_precision(x, ft), ft)
    elseif flat_model !== nothing
        cvis_model = eval_model(flat_model, flat_params, data.uv)
    else
        error("Bad image or model definition in call to simulate()")
    end

    # Compute model observables
    v2_model = vis_to_v2(cvis_model, data.indx_v2)
    _, t3amp_model, t3phi_model = vis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2, data.indx_t3_3)

    # Compute errors based on mode
    if mode == "copy_errors"
        v2_model_err = copy(data.v2_err)
        t3amp_model_err = copy(data.t3amp_err)
        t3phi_model_err = copy(data.t3phi_err)
    elseif mode == "copy_snr"
        v2_model_err = abs.(v2_model ./ data.v2 .* data.v2_err)
        gooddata = readoifits(in_oifits, filter_bad_data=true)[1]
        t3amp_model_err = abs.(t3amp_model ./ (data.t3amp ./ data.t3amp_err))
        t3phi_model_err = max.(abs.(t3phi_model ./ (data.t3phi ./ data.t3phi_err)), minimum(gooddata.t3phi_err))
    elseif mode == "noise_model"
        v2_model_err = errors.v2_multit * v2_model .+ errors.v2_addit
        t3amp_model_err = errors.t3amp_multit * t3amp_model .+ errors.t3amp_addit
        t3phi_model_err = zeros(length(t3phi_model)) .+ errors.t3phi_addit
    end

    # Fix NaN errors
    isbad = findall(isnan.(t3amp_model_err))
    t3amp_model_err[isbad] .= abs.(t3amp_model[isbad]) / 10
    isbad = findall(isnan.(t3phi_model_err))
    t3phi_model_err[isbad] .= abs.(t3phi_model[isbad]) / 10 .+ 1.0

    # Add noise
    v2_model .+= v2_model_err .* randn(length(v2_model))
    t3amp_model .+= t3amp_model_err .* randn(length(t3amp_model))
    t3phi_model .+= t3phi_model_err .* randn(length(t3phi_model))

    # Read input OIFITS as OIDataSet, modify data fields, write out
    ds = OIFITS.OIDataSet(in_oifits)

    # Replace V2 data+errors in all OI_VIS2 blocks
    v2_offset = 0
    for db in ds.vis2
        nrows = length(db.ucoord)
        nwavs = size(db.vis2data, 1)
        n = nwavs * nrows
        db.vis2data = reshape(v2_model[v2_offset+1:v2_offset+n], nwavs, nrows)
        db.vis2err  = reshape(v2_model_err[v2_offset+1:v2_offset+n], nwavs, nrows)
        v2_offset += n
    end

    # Replace T3 data+errors in all OI_T3 blocks
    t3_offset = 0
    for db in ds.t3
        nrows = length(db.u1coord)
        nwavs = size(db.t3amp, 1)
        n = nwavs * nrows
        db.t3amp    = reshape(t3amp_model[t3_offset+1:t3_offset+n], nwavs, nrows)
        db.t3amperr = reshape(t3amp_model_err[t3_offset+1:t3_offset+n], nwavs, nrows)
        db.t3phi    = reshape(t3phi_model[t3_offset+1:t3_offset+n], nwavs, nrows)
        db.t3phierr = reshape(t3phi_model_err[t3_offset+1:t3_offset+n], nwavs, nrows)
        t3_offset += n
    end

    OIFITS.write(out_file, ds; overwrite=true)
end
