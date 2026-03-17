using NFFT
using DelimitedFiles
using FITSIO
using CFITSIO
import TOML

# ─── Configuration structs ────────────────────────────────────────────────────

Base.@kwdef mutable struct FacilityConfig
    name::String                = ""
    lat::Float64                = 0.0
    lon::Float64                = 0.0
    alt::Float64                = 0.0
    throughput::Float64         = 1.0
    wind_speed::Float64         = 10.0
    r0::Float64                 = 0.1
    ntel::Int                   = 0
    tel_names::Vector{String}   = String[]
    sta_names::Vector{String}   = String[]
    tel_diams::Vector{Float64}  = Float64[]
    tel_gain::Vector{Float64}   = Float64[]
    sta_index::Vector{Int}      = Int[]
    sta_xyz::Matrix{Float64}    = zeros(Float64, 0, 3)
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

Base.@kwdef mutable struct CombinerConfig
    name::String                    = ""
    int_trans::Float64              = 1.0
    vis::Float64                    = 0.9
    n_pix_fringe::Int               = 128
    n_pix_photometry::Int           = 4
    flux_frac_photometry::Float64   = 0.5
    flux_frac_fringes::Float64      = 0.5
    throughput_photometry::Float64  = 0.1
    throughput_fringes::Float64     = 0.1
    n_splits::Int                   = 4
    read_noise::Float64             = 10.0
    quantum_efficiency::Float64     = 0.7
    v2_cal_err::Float64             = 0.005
    phase_cal_err::Float64          = 1.0
    incoh_int_time::Float64         = 100.0
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
    FacilityConfig(
        name       = String(get(d, "name", "")),
        lat        = Float64(get(d, "lat", 0.0)),
        lon        = Float64(get(d, "lon", 0.0)),
        alt        = Float64(get(d, "alt", 0.0)),
        throughput = Float64(get(d, "throughput", 1.0)),
        wind_speed = Float64(get(d, "wind_speed", 10.0)),
        r0         = Float64(get(d, "r0", 0.1)),
        ntel       = ntel,
        tel_names  = [String(get(t, "name", "T$i"))    for (i,t) in enumerate(tels)],
        sta_names  = [String(get(t, "station", get(t, "name", "S$i"))) for (i,t) in enumerate(tels)],
        tel_diams  = [Float64(get(t, "diameter", 1.0))  for t in tels],
        tel_gain   = [Float64(get(t, "gain", 1.0))      for t in tels],
        sta_index  = [Int(get(t, "index", i))            for (i,t) in enumerate(tels)],
        sta_xyz    = xyz,
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
    CombinerConfig(
        name                = String(get(d, "name", "")),
        int_trans           = Float64(get(d, "int_trans", 1.0)),
        vis                 = Float64(get(d, "vis", 0.9)),
        n_pix_fringe        = Int(get(d, "n_pix_fringe", 128)),
        n_pix_photometry    = Int(get(d, "n_pix_photometry", 4)),
        flux_frac_photometry = Float64(get(d, "flux_frac_photometry", 0.5)),
        flux_frac_fringes   = Float64(get(d, "flux_frac_fringes", 0.5)),
        throughput_photometry = Float64(get(d, "throughput_photometry", 0.1)),
        throughput_fringes  = Float64(get(d, "throughput_fringes", 0.1)),
        n_splits            = Int(get(d, "n_splits", 4)),
        read_noise          = Float64(get(d, "read_noise", 10.0)),
        quantum_efficiency  = Float64(get(d, "quantum_efficiency", 0.7)),
        v2_cal_err          = Float64(get(d, "v2_cal_err", 0.005)),
        phase_cal_err       = Float64(get(d, "phase_cal_err", 1.0)),
        incoh_int_time      = Float64(get(d, "incoh_int_time", 100.0)),
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
    flat_xyz = Float64.(f[(LinearIndices(f .== "sta_xyz"))[findall(f .== "sta_xyz")], 3:n*3+2])
    xyz = reshape(vec(flat_xyz), 3, n)'  # (ntel, 3)
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
    CombinerConfig(
        name                = String(_get("name")[1]),
        int_trans           = Float64(_get("int_trans")[1]),
        vis                 = Float64(_get("vis")[1]),
        n_pix_fringe        = Int(_get("n_pix_fringe")[1]),
        n_pix_photometry    = Int(_get("n_pix_photometry")[1]),
        flux_frac_photometry = Float64(_get("flux_frac_photometry")[1]),
        flux_frac_fringes   = Float64(_get("flux_frac_fringes")[1]),
        throughput_photometry = Float64(_get("throughput_photometry")[1]),
        throughput_fringes  = Float64(_get("throughput_fringes")[1]),
        n_splits            = Int(_get("n_splits")[1]),
        read_noise          = Float64(_get("read_noise")[1]),
        quantum_efficiency  = Float64(_get("quantum_efficiency")[1]),
        v2_cal_err          = Float64(_get("v2_cal_err")[1]),
        phase_cal_err       = Float64(_get("phase_cal_err")[1]),
        incoh_int_time      = Float64(_get("incoh_int_time")[1]),
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

function read_facility_file(path)
    p = _resolve_config(path)
    endswith(p, ".toml") ? _read_facility_toml(p) : _read_facility_txt(p)
end

function read_obs_file(path)
    p = _resolve_config(path)
    endswith(p, ".toml") ? _read_target_toml(p) : _read_target_txt(p)
end

function read_comb_file(path)
    p = _resolve_config(path)
    endswith(p, ".toml") ? _read_combiner_toml(p) : _read_combiner_txt(p)
end

function read_wave_file(path)
    p = _resolve_config(path)
    endswith(p, ".toml") ? _read_wave_toml(p) : _read_wave_txt(p)
end


# Zero-point flux for a mag=0 star in photons/s/m²/μm
# Log-linear interpolation from standard photometric bands (Cohen et al. 2003, Bessell et al. 1998)
const _F0_λ_μm = [0.55, 0.70, 0.90, 1.25, 1.65, 2.20, 3.50, 4.80, 10.0]
const _F0_phot = [3.6e10, 1.8e10, 8.4e9, 3.1e9, 1.15e9, 4.3e8, 8.1e7, 2.1e7, 9.7e5]

function zero_point_flux(λ_meters)
    λ_μm = clamp(λ_meters * 1e6, _F0_λ_μm[1], _F0_λ_μm[end])
    logλ = log.(_F0_λ_μm)
    logF = log.(_F0_phot)
    for i in 1:length(_F0_λ_μm)-1
        if λ_μm <= _F0_λ_μm[i+1]
            t = (log(λ_μm) - logλ[i]) / (logλ[i+1] - logλ[i])
            return exp(logF[i] + t * (logF[i+1] - logF[i]))
        end
    end
    return _F0_phot[end]
end

function compute_telescope_snr(mag, λ, δλ, facility, combiner, ntel, nwavs)
    # Single-telescope photometric SNR per telescope per wavelength
    # Returns (ntel, nwavs) matrix
    r0_ref   = facility.r0
    v_wind   = facility.wind_speed
    T_atm    = facility.throughput
    T_phot   = combiner.throughput_photometry
    QE       = combiner.quantum_efficiency
    n_pix_ph = combiner.n_pix_photometry
    σ_ron    = combiner.read_noise
    t_incoh  = combiner.incoh_int_time
    T_total  = T_atm * T_phot * QE
    snr = zeros(ntel, nwavs)
    for w in 1:nwavs
        F0 = zero_point_flux(λ[w])
        r0_λ = r0_ref * (λ[w] / 0.5e-6)^(6/5)
        τ0 = 0.31 * r0_λ / v_wind
        DIT = min(τ0, t_incoh)
        N_frames = max(1.0, t_incoh / DIT)
        δλ_μm = δλ[w] * 1e6
        for t in 1:ntel
            D_t = facility.tel_diams[t]
            η_t = min(1.0, (r0_λ / D_t)^2)
            Nphot = F0 * 10^(-mag / 2.5) * (π / 4 * D_t^2) * δλ_μm * DIT * T_total * η_t
            N_noise = Nphot + n_pix_ph * σ_ron^2
            snr[t, w] = sqrt(max(Nphot, 0.0)) / sqrt(max(N_noise, 1.0)) * sqrt(N_frames)
        end
    end
    return snr
end

function compute_baseline_snr(mag, λ, δλ, facility, combiner, v2_stations, nv2, nwavs)
    # Photometric SNR (for unresolved source V=1) per baseline per wavelength
    # Returns (nv2, nwavs) matrix
    r0_ref   = facility.r0
    v_wind   = facility.wind_speed
    T_atm    = facility.throughput
    T_inst   = combiner.throughput_fringes
    QE       = combiner.quantum_efficiency
    n_pix    = combiner.n_pix_fringe
    n_splits = combiner.n_splits
    σ_ron    = combiner.read_noise
    t_incoh  = combiner.incoh_int_time
    T_total  = T_atm * T_inst * QE
    snr = zeros(nv2, nwavs)
    for w in 1:nwavs
        F0 = zero_point_flux(λ[w])
        r0_λ = r0_ref * (λ[w] / 0.5e-6)^(6/5)
        τ0 = 0.31 * r0_λ / v_wind
        DIT = min(τ0, t_incoh)
        N_frames = max(1.0, t_incoh / DIT)
        δλ_μm = δλ[w] * 1e6
        for b in 1:nv2
            i, j = v2_stations[1, b], v2_stations[2, b]
            D_i = facility.tel_diams[i]
            D_j = facility.tel_diams[j]
            η_i = min(1.0, (r0_λ / D_i)^2)
            η_j = min(1.0, (r0_λ / D_j)^2)
            Nphot_i = F0 * 10^(-mag / 2.5) * (π / 4 * D_i^2) * δλ_μm * DIT * T_total * η_i
            Nphot_j = F0 * 10^(-mag / 2.5) * (π / 4 * D_j^2) * δλ_μm * DIT * T_total * η_j
            N_noise = (Nphot_i + Nphot_j) / n_splits + n_pix * σ_ron^2
            snr[b, w] = sqrt(max(Nphot_i * Nphot_j, 0.0)) / sqrt(max(N_noise, 1.0)) * sqrt(N_frames)
        end
    end
    return snr
end


#CODES FOR SIMULATING OIFITS BASED ON INPUT IMAGE AND EITHER INPUT OIFITS OR HOUR ANGLES

#Functions used in main Functions
# vis_to_t3_conj is removed — use vis_to_t3 from oichi2.jl instead.
# The conjugate was a workaround for a sign convention mismatch that is now fixed.


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

function simulate(facility,target,combiner,wavelength,dates,out_file; image::Union{String, Array{Float64,1}, Array{Float64,2}, Array{Float64, 3}, Array{Float64,4}}="", pixsize::Float64=0.1, flat_model::Union{FlatModel,Nothing}=nothing, flat_params::Vector{Float64}=Float64[], dft=false,nonoise=false, mag::Float64=2.0)
    #simulate an observation using input hour angles, info about array and combiner, and input image
    ntel=facility.ntel;

    # Override RA and DEC if not fixed (e.g. Sat) -- input RA and DEC should be in (decimal) degrees
    dec = target.decep0/180*pi
    ra = target.raep0/180*pi

    lst, hour_angles = hour_angle_calc(dates,facility.lon, ra*180/pi);
    nhours = length(hour_angles);
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
    # Determine if image or model

    if (((typeof(image) == String) && (image !="")) || ((typeof(image) != String) && (image !="")) )
        # TODO: Update polychromatic and dynamic imaging and non square images
        x = readfits(image)
        nx = size(x,1); 
        if ndims(x) == 2 # Monochromatic image
            println("Simulating monochromatic observations...")
            x = x/sum(x);
            ft = setup_nfft(uv, v2_indx_w, v2_indx_w, t3_indx_1_w, t3_indx_2_w, t3_indx_3_w, nx, pixsize);
            cvis_model = image_to_vis(x, ft);
        elseif ndims(x) == 3 # at the moment, we assume polychromatic
            println("Simulating polychromatic observations...")
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
        mjd_per_uv = repeat(repeat(Float64.(value.(modified_julian.(dates))), inner=nuv), outer=nwavs)
        cvis_model = eval_model(flat_model, flat_params, uv; wl=wl_per_uv, mjd=mjd_per_uv)
    else
        @warn("No image nor model definition in call to simulate()")
        @warn("Will generate zero visibilities")
        cvis_model = zeros(ComplexF64, size(uv,2))
    end

    # Compute true values of observables
    v2_model = vis_to_v2(cvis_model, v2_indx_w);
    t3_model, t3amp_model, t3phi_model = vis_to_t3(cvis_model, t3_indx_1_w, t3_indx_2_w, t3_indx_3_w);

    # Compute uncertainties — physical noise model from magnitude, instrument, and atmosphere
    snr0 = compute_baseline_snr(mag, λ, δλ, facility, combiner, v2_stations, nv2, nwavs)
    # Expand (nv2, nwavs) → flat vector matching v2_model ordering (baseline fastest, hour, wavelength)
    snr0_v2 = vec(repeat(max.(snr0, 1.0), nhours, 1))

    # V² error: σ(V²) ≈ 2|V|/SNR_0 + calibration systematic
    v2_model_err = 2.0 .* sqrt.(abs.(v2_model)) ./ snr0_v2 .+ combiner.v2_cal_err .* abs.(v2_model)

    # T3 errors: get SNR for each triangle baseline
    snr0_t3_1 = max.(snr0[t3_indx_1, :], 1.0)  # (nt3, nwavs)
    snr0_t3_2 = max.(snr0[t3_indx_2, :], 1.0)
    snr0_t3_3 = max.(snr0[t3_indx_3, :], 1.0)

    # |V| per triangle baseline from v2_model
    v2_per_bw = reshape(v2_model, nv2, nhours, nwavs)[:, 1, :]  # (nv2, nwavs)
    V_t3_1 = sqrt.(clamp.(v2_per_bw[t3_indx_1, :], 1e-20, Inf))
    V_t3_2 = sqrt.(clamp.(v2_per_bw[t3_indx_2, :], 1e-20, Inf))
    V_t3_3 = sqrt.(clamp.(v2_per_bw[t3_indx_3, :], 1e-20, Inf))

    # Inverse SNR² sum for the three baselines
    inv_snr2_sum = 1.0 ./ (V_t3_1 .* snr0_t3_1).^2 .+ 1.0 ./ (V_t3_2 .* snr0_t3_2).^2 .+ 1.0 ./ (V_t3_3 .* snr0_t3_3).^2

    # Closure phase error (degrees) with calibration floor
    sigma_t3phi_bw = sqrt.((180.0/π)^2 .* inv_snr2_sum .+ combiner.phase_cal_err^2)
    t3phi_model_err = vec(repeat(sigma_t3phi_bw, nhours, 1))

    # T3 amplitude error
    t3amp_model_err = abs.(t3amp_model) .* vec(repeat(sqrt.(inv_snr2_sum), nhours, 1))

    # OI_VIS observables: visibility amplitude and differential phase
    visamp_model = abs.(cvis_model[v2_indx_w])
    visphi_abs   = angle.(cvis_model[v2_indx_w]) .* (180.0 / π)  # absolute phase [deg]

    # Differential phase: subtract mean phase across wavelengths per baseline/hour
    phi_3d = reshape(visphi_abs, nv2, nhours, nwavs)
    phi_ref = sum(phi_3d, dims=3) ./ nwavs
    dphi_model = vec(phi_3d .- phi_ref)

    # Visibility amplitude error: σ(|V|) ≈ 1/SNR_0 + calibration systematic
    visamp_model_err = 1.0 ./ snr0_v2 .+ combiner.v2_cal_err .* visamp_model

    # Differential phase error [deg]: σ(φ) ≈ (180/π)/SNR_0 + calibration floor
    dphi_model_err = (180.0 / π) ./ snr0_v2 .+ combiner.phase_cal_err

    # OI_FLUX: per-telescope photometric flux (normalised to 1)
    snr_tel = compute_telescope_snr(mag, λ, δλ, facility, combiner, ntel, nwavs)
    nobs_flux = ntel * nhours
    flux_model     = ones(nobs_flux * nwavs)  # normalised total flux
    flux_model_err = vec(repeat(1.0 ./ max.(snr_tel, 1.0), nhours, 1))

    # Add noise
    if nonoise==true
        v2_model_err[:] .= 1.0
        t3amp_model_err[:] .= 1.0
        t3phi_model_err[:] .= 1.0
        visamp_model_err[:] .= 1.0
        dphi_model_err[:] .= 1.0
        flux_model_err[:] .= 1.0
    else
        v2_model       += v2_model_err     .* randn(length(v2_model))
        t3amp_model    += t3amp_model_err  .* randn(length(t3amp_model))
        t3phi_model    += t3phi_model_err  .* randn(length(t3phi_model))
        visamp_model   += visamp_model_err .* randn(length(visamp_model))
        dphi_model     += dphi_model_err   .* randn(length(dphi_model))
        flux_model     += flux_model_err   .* randn(length(flux_model))
    end

    # --- Build OIFITS data structures and save ---
    nobs_v2  = nv2 * nhours
    nobs_t3  = nt3 * nhours
    mjd_vals = Float64.(value.(modified_julian.(dates)))
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
    db_vis2.flag      = fill(false, nwavs, nobs_v2)

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
    db_vis.flag      = fill(false, nwavs, nobs_v2)

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
    db_t3.flag      = fill(false, nwavs, nobs_t3)

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
    db_flux.flag      = fill(false, nwavs, nobs_flux)

    # Assemble and write
    ds = OIFITS.OIDataSet(tgt, arr, ins, db_vis2, db_vis, db_t3, db_flux)
    OIFITS.write(out_file, ds; overwrite=true)
end

function simulate_from_oifits(in_oifits, out_file; mode="copy_errors", errors=[],  image::Union{String, Array{Float64,1}, Array{Float64,2}, Array{Float64, 3}, Array{Float64,4}}="", pixsize::Float64=0.1, flat_model::Union{FlatModel,Nothing}=nothing, flat_params::Vector{Float64}=Float64[])
 
    outfilename= string("!", out_file)
    data = (readoifits(in_oifits,filter_bad_data=false))[1,1];
    oifits = FITS(in_oifits);
    #setup simulation
    nuv = data.nuv
    # Compute complex visibilities from either a truth image or a model
    cvis_model = ComplexF64[]
    # Determine if image or model
    if (((typeof(image) == String) && (image !="")) || ((typeof(image) != String) && (image !="")) )
        # Truth image
        # TODO: Polychromatic and dynamic imaging
        x = readfits(image)
        nx = size(x,1);
        x = vec(x)/sum(x);
        ft = setup_nfft(data, nx, pixsize);
        cvis_model = image_to_vis(x, ft);
    elseif flat_model !== nothing
        # Flat-dict parametric model
        cvis_model = eval_model(flat_model, flat_params, data.uv)
    else
        error("Bad image or model definition in call to simulate()")
    end

    f = fits_create_file(outfilename);
    #setup initial table
    copy_oi_header(f,oifits[1]);
    ioiarraytables=[]
    iotargettables=[]
    swavetables=[]
    iwavetables=[]
    svistables=[]
    ivistables=[]
    st3tables=[]
    it3tables=[]

    #setup new table arrays

    for itable=2:length(oifits)
        if (read_key(oifits[itable],"EXTNAME")[1] == "OI_ARRAY")
            push!(ioiarraytables,itable)
        end
        if (read_key(oifits[itable],"EXTNAME")[1] == "OI_TARGET")
            push!(iotargettables,itable)
        end
        if (read_key(oifits[itable],"EXTNAME")[1] == "OI_WAVELENGTH")
            push!(swavetables,read_key(oifits[itable],"INSNAME")[1])
            push!(iwavetables,itable)
        end
        if (read_key(oifits[itable],"EXTNAME")[1] == "OI_VIS2")
            push!(svistables,read_key(oifits[itable],"INSNAME")[1])
            push!(ivistables,itable)
        end
        if (read_key(oifits[itable],"EXTNAME")[1] == "OI_T3")
            push!(st3tables,read_key(oifits[itable],"INSNAME")[1])
            push!(it3tables,itable)
        end
    end
    #TO GET OTHER PIECES OF INFO
    #get OI_ARRAY details (currently only handles one Array)
    for itable=1:length(ioiarraytables)
        telnames=read(oifits[ioiarraytables[itable]],"TEL_NAME")
                 read(oifits[ioiarraytables[1]],"TEL_NAME")
        sta_names=read(oifits[ioiarraytables[itable]],"STA_NAME")
        sta_index=read(oifits[ioiarraytables[itable]],"STA_INDEX")
        tel_diams= read(oifits[ioiarraytables[itable]],"DIAMETER")
        staxyz=read(oifits[ioiarraytables[itable]],"STAXYZ");
        station_xyz=transpose(staxyz);
        N=length(read(oifits[ioiarraytables[itable]],"TEL_NAME"))
        if read_key(oifits[ioiarraytables[itable]],"OI_REVN")[1] == 2
            print,"HERE"
            fov=read(oifits[ioiarraytables[itable]],"FOV")
            fovtype=read(oifits[ioiarraytables[itable]],"FOVTYPE")
            oi_array=[telnames,sta_names,sta_index,tel_diams,staxyz,fov,fovtype];
            copy_oi_array(f,oifits[ioiarraytables[itable]],oi_array);
        else
            oi_array=[telnames,sta_names,sta_index,tel_diams,staxyz];
            copy_oi_array(f,oifits[ioiarraytables[itable]],oi_array);
        end
    end

    #get target array details

    target_id=read(oifits["OI_TARGET"],"TARGET_ID");
    target=read(oifits["OI_TARGET"],"TARGET");
    raep0=read(oifits["OI_TARGET"],"RAEP0");
    decep0=read(oifits["OI_TARGET"],"DECEP0");
    equinox=read(oifits["OI_TARGET"],"EQUINOX");
    ra_err=read(oifits["OI_TARGET"],"RA_ERR");
    dec_err=read(oifits["OI_TARGET"],"DEC_ERR");
    sysvel=-read(oifits["OI_TARGET"],"SYSVEL");
    veltyp=read(oifits["OI_TARGET"],"VELTYP");
    veldef=read(oifits["OI_TARGET"],"VELDEF");
    pmra=read(oifits["OI_TARGET"],"PMRA");
    pmdec=read(oifits["OI_TARGET"],"PMDEC");
    pmra_err=read(oifits["OI_TARGET"],"PMRA_ERR");
    pmdec_err=read(oifits["OI_TARGET"],"PMDEC_ERR");
    parallax=read(oifits["OI_TARGET"],"PARALLAX"); #COULD PULL FROM GAIA
    para_err=read(oifits["OI_TARGET"],"PARA_ERR");
    spectyp=read(oifits["OI_TARGET"],"SPECTYP");
    target_array=[target_id,target,raep0,decep0,equinox,ra_err,dec_err,sysvel,veltyp,veldef,pmra,pmdec,pmra_err,pmdec_err,parallax,para_err,spectyp];
    copy_oi_target(f,oifits[iotargettables[length(iotargettables)]],target_array);
    #GET WAVE TABLES
    for itable=1:length(iwavetables)
        eff_wave=read(oifits[iwavetables[itable]],"EFF_WAVE")
        eff_band=read(oifits[iwavetables[itable]],"EFF_BAND")
        wave_array=[eff_wave,eff_band];
        copy_oi_wavelength(f,oifits[iwavetables[itable]],wave_array);
    end

    #GET V2 INFO
    v2_model = vis_to_v2(cvis_model, data.indx_v2 )# based on uv points
    v2_model_err = Float64[]
    if mode == "copy_errors"
        v2_model_err = data.v2_err
    elseif mode == "copy_snr"
        v2_model_err = abs.(v2_model./data.v2.*data.v2_err)
    elseif mode == "noise_model"
        v2_model_err = errors.v2_multit*v2_model .+ errors.v2_addit;
    end

    v2_model += v2_model_err.*randn(length(v2_model))

    #Now to fill the tables_
    visindxstart=1 
    for itable=1:length(ivistables); #for each table
        v2_model_table=[]
        v2_model_err_table=[]
        uvis=read(oifits[ivistables[itable]],"UCOORD");
        vvis=read(oifits[ivistables[itable]],"VCOORD")
        wavetablenum=iwavetables[findall(swavetables.==svistables[itable])]
        eff_wave= read(oifits[wavetablenum[1]],"EFF_WAVE")
        eff_band = read(oifits[wavetablenum[1]],"EFF_BAND")
        nwavs=length(eff_wave);
        for hour=1:length(uvis)
            for wave=1:nwavs
                v2_model_table=push!(v2_model_table,v2_model[visindxstart+wave-1]);
                v2_model_err_table=push!(v2_model_err_table,v2_model_err[visindxstart+wave-1]);
            end
            visindxstart=visindxstart+nwavs
        end
        target_id_vis2=read(oifits[ivistables[itable]],"TARGET_ID");
        time_vis2=read(oifits[ivistables[itable]],"TIME");
        mjd_vis2=read(oifits[ivistables[itable]],"MJD");
        int_time_vis2=read(oifits[ivistables[itable]],"INT_TIME");
        v2_model_stations=read(oifits[ivistables[itable]],"STA_INDEX");
        nrowvis=Int(length(v2_model_table)/nwavs)
        v2_model_reshape=reshape(v2_model_table,(nwavs,nrowvis));
        v2_model_err_reshape=reshape(v2_model_err_table,(nwavs,nrowvis));
        flag_vis2=read(oifits[ivistables[itable]],"FLAG")
        vis2_array=[target_id_vis2,time_vis2,mjd_vis2,int_time_vis2,v2_model_reshape,v2_model_err_reshape,uvis,vvis,v2_model_stations,flag_vis2'];
        copy_oi_vis2(f,oifits[ivistables[itable]],vis2_array);
    end

    #GET T3 NOW
    t3_model, t3amp_model, t3phi_model = vis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2, data.indx_t3_3);
    if mode == "copy_errors"
        t3amp_model_err = data.t3amp_err
        t3phi_model_err = data.t3phi_err
        isbad = findall(isnan.(t3amp_model_err))
        t3amp_model_err[isbad] .= abs.(t3amp_model[isbad])/10
        isbad = findall(isnan.(t3phi_model_err))
        t3phi_model_err[isbad] .= abs.(t3phi_model[isbad])/10 .+ 1.0 
    elseif mode == "copy_snr"
        gooddata = (readoifits(in_oifits,filter_bad_data=true))[1,1];
        t3amp_model_err = abs.(t3amp_model./(data.t3amp./data.t3amp_err))
        t3phi_model_err = max.(abs.(t3phi_model./(data.t3phi./data.t3phi_err)), minimum(gooddata.t3phi_err))
        isbad = findall(isnan.(t3amp_model_err))
        t3amp_model_err[isbad] .= abs.(t3amp_model[isbad])/10
        isbad = findall(isnan.(t3phi_model_err))
        t3phi_model_err[isbad] .= abs.(t3phi_model[isbad])/10 .+ 1.0 
    elseif mode == "noise_model"
        t3amp_model_err = errors.t3amp_multit*t3amp_model .+ errors.t3amp_addit;   
        t3phi_model_err = zeros(length(t3phi_model)) .+ errors.t3phi_addit;
    end
    t3amp_model += t3amp_model_err.*randn(length(t3amp_model)) 
    t3phi_model += t3phi_model_err.*randn(length(t3phi_model))
    t3indxstart=1
    for itable=1:length(it3tables); #for each table
        t3amp_model_table=[]
        t3amp_model_err_table=[]
        t3phi_model_table=[]
        t3phi_model_err_table=[]
        u1=read(oifits[it3tables[itable]],"U1COORD");
        v1=read(oifits[it3tables[itable]],"V1COORD")
        u2=read(oifits[it3tables[itable]],"U2COORD");
        v2=read(oifits[it3tables[itable]],"V2COORD")
        wavetablenum=iwavetables[findall(swavetables.==st3tables[itable])]
        eff_wave= read(oifits[wavetablenum[1]],"EFF_WAVE")
        eff_band = read(oifits[wavetablenum[1]],"EFF_BAND")
        nwavs=length(eff_wave);
        for hour=1:length(u1)
            for wave=1:nwavs
                t3amp_model_table=push!(t3amp_model_table,t3amp_model[t3indxstart+wave-1]);
                t3amp_model_err_table=push!(t3amp_model_err_table,t3amp_model_err[t3indxstart+wave-1]);
                t3phi_model_table=push!(t3phi_model_table,t3phi_model[t3indxstart+wave-1]);
                t3phi_model_err_table=push!(t3phi_model_err_table,t3phi_model_err[t3indxstart+wave-1]);
            end
            t3indxstart=t3indxstart+nwavs
        end
        target_id_t3=read(oifits[it3tables[itable]],"TARGET_ID");
        time_t3=read(oifits[it3tables[itable]],"TIME");
        mjd_t3=read(oifits[it3tables[itable]],"MJD");
        int_time_t3=read(oifits[it3tables[itable]],"INT_TIME");
        t3_model_stations=read(oifits[it3tables[itable]],"STA_INDEX");
        nrowt3=Int(length(t3amp_model_table)/nwavs)
        t3amp_model_reshape=reshape(t3amp_model_table,(nwavs,nrowt3));
        t3amp_model_err_reshape=reshape(t3amp_model_err_table,(nwavs,nrowt3));
        t3phi_model_reshape=reshape(t3phi_model_table,(nwavs,nrowt3));
        t3phi_model_err_reshape=reshape(t3phi_model_err_table,(nwavs,nrowt3));
        flag_t3=read(oifits[it3tables[itable]],"FLAG")
        t3_array=[target_id_t3,time_t3,mjd_t3,int_time_t3,t3amp_model_reshape,t3amp_model_err_reshape,t3phi_model_reshape,t3phi_model_err_reshape,u1,v1,u2,v2,t3_model_stations,flag_t3'];
        copy_oi_t3(f,oifits[it3tables[itable]],t3_array);
    end
    fits_close_file(f)
end
