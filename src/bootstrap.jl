# bootstrap.jl
#
# Uncertainty estimation by resampling, for parametric model fitting.
#
#   resample_blocks(data, blk)  Nonparametric block bootstrap.  Keeps the data
#                               values untouched and instead resamples which
#                               observations are used, in blocks of
#                               (MJD, telescope configuration) — the same
#                               atomic unit as PMOIRED.  Answers: "how much
#                               does the answer depend on which observations I
#                               happen to have?"  Captures correlated
#                               calibration errors and mis-stated error bars,
#                               at the price of needing many independent
#                               blocks.
#
# `bootstrap_fit` drives it over many replicates and summarises the resulting
# parameter distribution (median, 16/84 percentiles, covariance).
#
#   perturb_data(data)          Adds Gaussian noise drawn from the error bars.
#                               A simulation utility, NOT an uncertainty
#                               estimator: refitting such replicates only
#                               reproduces the analytic covariance, and is blind
#                               to exactly the error sources that matter (see
#                               demos/bootstrap_validation).

using Random, Statistics, Printf


# ─────────────────────────────────────────────────────────────────────────────
# Parametric Monte Carlo
# ─────────────────────────────────────────────────────────────────────────────

"""
    perturb_data(data::OIdata; rng=Random.default_rng()) -> OIdata

Return a copy of `data` with independent Gaussian noise added to every
observable, each draw scaled by that point's error bar.  uv coordinates,
flags and the number of data points are unchanged.

Useful for building noisy realisations of a dataset — simulation, end-to-end
tests, sensitivity checks.

!!! warning "Not an uncertainty estimator"
    Refitting `perturb_data` replicates is a parametric Monte Carlo, not a
    bootstrap: it assumes the quoted error bars are correct and uncorrelated,
    so it can only ever reproduce the analytic covariance of the fit.  In the
    validation study under `demos/bootstrap_validation` it was never better
    than `fit_model_lsqfit`'s covariance and, when calibration systematics were
    present, it understated the true scatter by a factor 11 to 25.  Use
    [`bootstrap_fit`](@ref) for uncertainties.
"""
function perturb_data(data::OIdata; rng::AbstractRNG=Random.default_rng())
    d = deepcopy(data)
    d.v2      .+= d.v2_err      .* randn(rng, length(d.v2))
    d.t3amp   .+= d.t3amp_err   .* randn(rng, length(d.t3amp))
    d.t3phi   .+= d.t3phi_err   .* randn(rng, length(d.t3phi))
    d.visamp  .+= d.visamp_err  .* randn(rng, length(d.visamp))
    d.visphi  .+= d.visphi_err  .* randn(rng, length(d.visphi))
    if length(d.flux) > 0
        d.flux .+= d.flux_err   .* randn(rng, length(d.flux))
    end
    return d
end

"""
    resample_data(data::OIdata)

Deprecated alias for [`perturb_data`](@ref).

Despite the name, this function never resampled anything: it adds Gaussian
noise drawn from the error bars.  For a genuine bootstrap (resampling which
observations are used) see [`bootstrap_fit`](@ref) and
[`resample_blocks`](@ref).
"""
function resample_data(data::OIdata; kwargs...)
    @warn("`resample_data` adds Gaussian noise from the error bars — it is a " *
          "parametric Monte Carlo, not a bootstrap. Use `perturb_data` for that, " *
          "or `bootstrap_fit` for a genuine block bootstrap.", maxlog=1)
    return perturb_data(data; kwargs...)
end


# ─────────────────────────────────────────────────────────────────────────────
# Data blocks
# ─────────────────────────────────────────────────────────────────────────────

"""
    DataBlocks

Partition of an `OIdata` into resampling units.  Each block holds the row
indices it owns in each observable table (`OI_VIS`, `OI_VIS2`, `OI_T3`,
`OI_FLUX`).  Built by [`data_blocks`](@ref).
"""
struct DataBlocks
    keys        ::Vector{String}
    idx_vis     ::Vector{Vector{Int}}
    idx_v2      ::Vector{Vector{Int}}
    idx_t3      ::Vector{Vector{Int}}
    idx_flux    ::Vector{Vector{Int}}
    granularity ::Symbol
end

Base.length(b::DataBlocks) = length(b.keys)

function Base.show(io::IO, b::DataBlocks)
    nvis  = sum(length, b.idx_vis;  init=0)
    nv2   = sum(length, b.idx_v2;   init=0)
    nt3   = sum(length, b.idx_t3;   init=0)
    nflux = sum(length, b.idx_flux; init=0)
    npts  = nvis + nv2 + nt3 + nflux
    @printf(io, "DataBlocks: %d blocks (%s), %d points [vis=%d v2=%d t3=%d flux=%d]",
            length(b), string(b.granularity), npts, nvis, nv2, nt3, nflux)
    if length(b) > 0
        @printf(io, "\n  %.1f points per block on average", npts / length(b))
    end
end

_mjd_str(mjd, digits) = string(round(Float64(mjd), digits=digits))

"""
    data_blocks(data::OIdata; granularity=:config, mjd_digits=5) -> DataBlocks

Partition `data` into resampling blocks.

# Granularity
- `:config` (default) — one block per (MJD, telescope configuration), where
  the configuration is the baseline for V²/VIS, the triangle for T3 and the
  telescope for FLUX.  All wavelength channels of a block stay together, and
  observables sharing a configuration at the same MJD (e.g. `|V|` and `V²` on
  the same baseline) stay together too.  This is PMOIRED's "spectral vector".
- `:epoch` — one block per MJD: an entire observation, all baselines and all
  observables, is kept or dropped as a unit.  Coarser, more conservative;
  needs many epochs to be usable.
- `:point` — one block per data point.  The classical i.i.d. bootstrap; it
  destroys the spectral and per-configuration correlation structure and will
  therefore *underestimate* uncertainties on real data.  Provided for
  comparison studies.

`mjd_digits` sets the rounding used to group MJDs into a common epoch
(default 5 decimals ≈ 0.86 s, matching PMOIRED).

!!! note "MJD precision"
    Blocking needs the MJDs at full precision: near MJD 55000 the representable
    `Float32` values are 3.9e-3 d = 5.6 minutes apart, so exposures closer than
    that would collapse onto one value and their blocks would merge.  `readoifits`
    therefore stores `v2_mjd`, `t3_mjd`, `vis_mjd` and `flux_mjd` as `Float64`
    whatever `T` is — the block structure is identical at `T=Float32` and
    `T=Float64`.  (`uv_mjd` follows `T`: it is the array handed to the model
    evaluator as `\$MJD`, where the precision is irrelevant.)
"""
function data_blocks(data::OIdata; granularity::Symbol=:config, mjd_digits::Int=5)
    granularity in (:config, :epoch, :point) ||
        error("data_blocks: granularity must be :config, :epoch or :point (got :$granularity)")

    keys     = String[]
    lookup   = Dict{String,Int}()
    idx_vis  = Vector{Int}[]
    idx_v2   = Vector{Int}[]
    idx_t3   = Vector{Int}[]
    idx_flux = Vector{Int}[]

    function block_id!(key::String)
        get!(lookup, key) do
            push!(keys, key)
            push!(idx_vis,  Int[]); push!(idx_v2,   Int[])
            push!(idx_t3,   Int[]); push!(idx_flux, Int[])
            length(keys)
        end
    end

    function key_of(mjd, sta::AbstractVector{<:Integer}, tag::String, i::Int)
        granularity === :point && return "$tag#$i"
        granularity === :epoch && return _mjd_str(mjd, mjd_digits)
        return _mjd_str(mjd, mjd_digits) * "|" * join(sort(collect(sta)), "-")
    end

    # ── OI_VIS ───────────────────────────────────────────────────────────────
    nvis = max(length(data.visamp), length(data.visphi))
    for i in 1:nvis
        sta = size(data.vis_sta_index, 2) >= i ? data.vis_sta_index[:, i] : Int[]
        mjd = length(data.vis_mjd) >= i ? data.vis_mjd[i] : 0.0
        push!(idx_vis[block_id!(key_of(mjd, sta, "vis", i))], i)
    end

    # ── OI_VIS2 ──────────────────────────────────────────────────────────────
    for i in eachindex(data.v2)
        sta = size(data.v2_sta_index, 2) >= i ? data.v2_sta_index[:, i] : Int[]
        mjd = length(data.v2_mjd) >= i ? data.v2_mjd[i] : 0.0
        push!(idx_v2[block_id!(key_of(mjd, sta, "v2", i))], i)
    end

    # ── OI_T3 ────────────────────────────────────────────────────────────────
    nt3 = max(length(data.t3amp), length(data.t3phi))
    for i in 1:nt3
        sta = size(data.t3_sta_index, 2) >= i ? data.t3_sta_index[:, i] : Int[]
        mjd = length(data.t3_mjd) >= i ? data.t3_mjd[i] : 0.0
        push!(idx_t3[block_id!(key_of(mjd, sta, "t3", i))], i)
    end

    # ── OI_FLUX ──────────────────────────────────────────────────────────────
    for i in eachindex(data.flux)
        sta = length(data.flux_sta_index) >= i ? [data.flux_sta_index[i]] : Int[]
        mjd = length(data.flux_mjd) >= i ? data.flux_mjd[i] : 0.0
        push!(idx_flux[block_id!(key_of(mjd, sta, "flux", i))], i)
    end

    return DataBlocks(keys, idx_vis, idx_v2, idx_t3, idx_flux, granularity)
end


# ─────────────────────────────────────────────────────────────────────────────
# Block resampling
# ─────────────────────────────────────────────────────────────────────────────

_reidx(v::AbstractVector, idx) = isempty(v) ? v : v[idx]
_reidx(m::AbstractMatrix, idx) = size(m, 2) == 0 ? m : m[:, idx]

function _expand(lists::Vector{Vector{Int}}, counts::Vector{Int})
    out = Int[]
    for b in eachindex(counts), _ in 1:counts[b]
        append!(out, lists[b])
    end
    sort!(out)
    return out
end

"""
    block_counts(nblocks, mode; rng) -> Vector{Int}

Multiplicity drawn for each block under the given resampling `mode`.  Exposed
mainly for testing and for the calibration study of the different schemes.

# Modes
- `:replacement` (default) — draw `nblocks` blocks with replacement.  This is
  the textbook nonparametric bootstrap; multiplicities are ≈ Poisson(1), i.e.
  mean 1 and variance 1.
- `:halfsample` — keep a random half of the blocks, each once.  Balanced
  repeated replication: correctly calibrated (the finite-population correction
  cancels the factor 2 from halving the data), and cheaper since each replicate
  fit uses half the data.  In the validation study under
  `demos/bootstrap_validation` this was the best-calibrated scheme whenever
  blocks were numerous (ratio 0.94-1.06); with very few blocks it turns
  conservative (1.24-1.56), since each replicate then fits very little data.
- `:pmoired` — PMOIRED's scheme: build two independent copies of the dataset,
  drop a random half of the blocks in each, and fit the union.  Multiplicities
  are 0/1/2 with probabilities 1/4, 1/2, 1/4, so their variance is 1/2 rather
  than 1.

!!! warning "`:pmoired` is biased low by √2 — for comparison only"
    Because its block multiplicities have variance 1/2 instead of 1, this
    scheme returns uncertainties that are too small by roughly a factor √2.
    Measured against simulated truth (`demos/bootstrap_validation`) it gives
    0.62–0.79 of the true parameter scatter in every regime tested, with a 1σ
    coverage of ~0.50 instead of 0.683.  Use it to reproduce or cross-check a
    PMOIRED result, never to quote an error bar: use `:replacement` (default)
    or `:halfsample` for that.

`resample_blocks` additionally accepts `mode=:weights`, the multiplier
(Bayesian) bootstrap, which draws continuous weights rather than integer
multiplicities — see [`block_weights`](@ref).
"""
function block_counts(nblocks::Int, mode::Symbol; rng::AbstractRNG=Random.default_rng())
    nblocks > 0 || error("block_counts: no blocks to resample")
    counts = zeros(Int, nblocks)
    if mode === :replacement
        for _ in 1:nblocks
            counts[rand(rng, 1:nblocks)] += 1
        end
    elseif mode === :halfsample
        p = randperm(rng, nblocks)
        counts[p[1:(nblocks ÷ 2)]] .= 1
    elseif mode === :pmoired
        @warn("mode=:pmoired reproduces PMOIRED's two-half-samples scheme, whose " *
              "block multiplicities have variance 1/2 instead of 1: the resulting " *
              "uncertainties are too small by about √2 (measured 0.62-0.79 of the " *
              "true scatter). Use it for comparison only — quote error bars from " *
              ":replacement or :halfsample.", maxlog=1)
        for _ in 1:2
            p = randperm(rng, nblocks)
            counts[p[(nblocks ÷ 2 + 1):end]] .+= 1
        end
    else
        error("block_counts: unknown mode :$mode (use :replacement, :halfsample or :pmoired)")
    end
    return counts
end

"""
    block_weights(nblocks; rng) -> Vector{Float64}

Random weights for the multiplier ("Bayesian") bootstrap: `nblocks` draws from
Dirichlet(1, …, 1) rescaled to mean 1, i.e. mean 1 and variance 1 − 1/nblocks,
the same first two moments as the multiplicity of a draw with replacement.

Unlike `block_counts`, these are continuous: no block is ever dropped entirely
and none is duplicated, so a replicate costs exactly one dataset instead of
one-and-a-bit.  See [`apply_block_weights`](@ref).
"""
function block_weights(nblocks::Int; rng::AbstractRNG=Random.default_rng())
    nblocks > 0 || error("block_weights: no blocks to resample")
    w = randexp(rng, nblocks)
    return w .* (nblocks / sum(w))
end

"""
    apply_block_weights(data, blocks, w) -> OIdata

Multiplier bootstrap replicate: instead of duplicating and dropping blocks, give
block `i` the weight `w[i]` in the χ² by dividing its error bars by `sqrt(w[i])`.

For a weighted least-squares fit this is the standard multiplier (or Bayesian)
bootstrap: with `mean(w) = 1` and `var(w) = 1` it has the same asymptotic
behaviour as resampling blocks with replacement, at the cost of one dataset per
replicate and with no discrete jumps as blocks enter and leave.

The reduced χ² of such a replicate is not meaningful; the parameter scatter is.
"""
function apply_block_weights(data::OIdata, blocks::DataBlocks, w::AbstractVector{<:Real})
    length(w) == length(blocks) ||
        error("apply_block_weights: w has length $(length(w)), expected $(length(blocks))")
    all(x -> x > 0, w) || error("apply_block_weights: weights must be strictly positive")

    d = deepcopy(data)
    for b in eachindex(w)
        s = sqrt(w[b])
        for i in blocks.idx_vis[b]
            d.visamp_err[i] /= s
            d.visphi_err[i] /= s
        end
        for i in blocks.idx_v2[b]
            d.v2_err[i] /= s
        end
        for i in blocks.idx_t3[b]
            d.t3amp_err[i] /= s
            d.t3phi_err[i] /= s
            if !isempty(d.t3phi_vonmises_err)
                d.t3phi_vonmises_err[i] *= w[b]
                d.t3phi_vonmises_chi2_offset[i] *= w[b]
            end
        end
        for i in blocks.idx_flux[b]
            d.flux_err[i] /= s
        end
    end
    return d
end

"""
    resample_blocks(data, blocks; mode=:replacement, rng=Random.default_rng()) -> OIdata
    resample_blocks(data; granularity=:config, mode=:replacement, ...) -> OIdata

Bootstrap replicate of `data`: the observable *values* and their error bars are
left untouched, and the blocks of `blocks` are drawn according to `mode` (see
[`block_counts`](@ref)).  A block drawn twice contributes its data twice; a
block not drawn is absent.

The uv table is carried over unchanged — duplicated data points simply share
the same uv row — so the model evaluation cost per replicate is unchanged.
"""
function resample_blocks(data::OIdata, blocks::DataBlocks;
                         mode::Symbol=:replacement,
                         rng::AbstractRNG=Random.default_rng())
    if mode === :weights
        return apply_block_weights(data, blocks, block_weights(length(blocks); rng=rng))
    end
    counts = block_counts(length(blocks), mode; rng=rng)
    return apply_block_counts(data, blocks, counts)
end

function resample_blocks(data::OIdata; granularity::Symbol=:config, mjd_digits::Int=5,
                         mode::Symbol=:replacement,
                         rng::AbstractRNG=Random.default_rng())
    blocks = data_blocks(data; granularity=granularity, mjd_digits=mjd_digits)
    return resample_blocks(data, blocks; mode=mode, rng=rng)
end

"""
    apply_block_counts(data, blocks, counts) -> OIdata

Build the `OIdata` in which each block of `blocks` appears `counts[i]` times.
"""
function apply_block_counts(data::OIdata, blocks::DataBlocks, counts::Vector{Int})
    length(counts) == length(blocks) ||
        error("apply_block_counts: counts has length $(length(counts)), expected $(length(blocks))")

    iv  = _expand(blocks.idx_vis,  counts)
    i2  = _expand(blocks.idx_v2,   counts)
    i3  = _expand(blocks.idx_t3,   counts)
    ifl = _expand(blocks.idx_flux, counts)

    d = deepcopy(data)

    # ── OI_VIS ───────────────────────────────────────────────────────────────
    d.visamp        = _reidx(data.visamp,          iv)
    d.visamp_err    = _reidx(data.visamp_err,      iv)
    d.visphi        = _reidx(data.visphi,          iv)
    d.visphi_err    = _reidx(data.visphi_err,      iv)
    d.vis_baseline  = _reidx(data.vis_baseline,    iv)
    d.vis_mjd       = _reidx(data.vis_mjd,         iv)
    d.vis_lam       = _reidx(data.vis_lam,         iv)
    d.vis_dlam      = _reidx(data.vis_dlam,        iv)
    d.vis_flag      = _reidx(data.vis_flag,        iv)
    d.indx_vis      = _reidx(data.indx_vis,        iv)
    d.vis_sta_index = _reidx(data.vis_sta_index,   iv)
    d.visamp_corr_idx = _reidx(data.visamp_corr_idx, iv)
    d.visphi_corr_idx = _reidx(data.visphi_corr_idx, iv)

    # ── OI_VIS2 ──────────────────────────────────────────────────────────────
    d.v2            = _reidx(data.v2,              i2)
    d.v2_err        = _reidx(data.v2_err,          i2)
    d.v2_baseline   = _reidx(data.v2_baseline,     i2)
    d.v2_mjd        = _reidx(data.v2_mjd,          i2)
    d.v2_lam        = _reidx(data.v2_lam,          i2)
    d.v2_dlam       = _reidx(data.v2_dlam,         i2)
    d.v2_flag       = _reidx(data.v2_flag,         i2)
    d.indx_v2       = _reidx(data.indx_v2,         i2)
    d.v2_sta_index  = _reidx(data.v2_sta_index,    i2)
    d.v2_corr_idx   = _reidx(data.v2_corr_idx,     i2)

    # ── OI_T3 ────────────────────────────────────────────────────────────────
    d.t3amp                     = _reidx(data.t3amp,                      i3)
    d.t3amp_err                 = _reidx(data.t3amp_err,                  i3)
    d.t3phi                     = _reidx(data.t3phi,                      i3)
    d.t3phi_err                 = _reidx(data.t3phi_err,                  i3)
    d.t3phi_vonmises_err        = _reidx(data.t3phi_vonmises_err,         i3)
    d.t3phi_vonmises_chi2_offset= _reidx(data.t3phi_vonmises_chi2_offset, i3)
    d.t3_baseline               = _reidx(data.t3_baseline,                i3)
    d.t3_maxbaseline            = _reidx(data.t3_maxbaseline,             i3)
    d.t3_mjd                    = _reidx(data.t3_mjd,                     i3)
    d.t3_lam                    = _reidx(data.t3_lam,                     i3)
    d.t3_dlam                   = _reidx(data.t3_dlam,                    i3)
    d.t3_flag                   = _reidx(data.t3_flag,                    i3)
    d.indx_t3_1                 = _reidx(data.indx_t3_1,                  i3)
    d.indx_t3_2                 = _reidx(data.indx_t3_2,                  i3)
    d.indx_t3_3                 = _reidx(data.indx_t3_3,                  i3)
    d.t3_sta_index              = _reidx(data.t3_sta_index,               i3)
    d.t3amp_corr_idx            = _reidx(data.t3amp_corr_idx,             i3)
    d.t3phi_corr_idx            = _reidx(data.t3phi_corr_idx,             i3)

    # ── OI_FLUX ──────────────────────────────────────────────────────────────
    d.flux           = _reidx(data.flux,           ifl)
    d.flux_err       = _reidx(data.flux_err,       ifl)
    d.flux_mjd       = _reidx(data.flux_mjd,       ifl)
    d.flux_lam       = _reidx(data.flux_lam,       ifl)
    d.flux_dlam      = _reidx(data.flux_dlam,      ifl)
    d.flux_flag      = _reidx(data.flux_flag,      ifl)
    d.flux_sta_index = _reidx(data.flux_sta_index, ifl)
    d.flux_corr_idx  = _reidx(data.flux_corr_idx,  ifl)

    # ── Counts (same validity rules as readoifits) ───────────────────────────
    _nvalid(x, e) = count(i -> !isnan(x[i]) && !isnan(e[i]) && e[i] > 0, eachindex(x))
    d.nv2     = length(d.v2)
    d.nflux   = length(d.flux)
    d.nvisamp = isempty(d.visamp) ? 0 : _nvalid(d.visamp, d.visamp_err)
    d.nvisphi = isempty(d.visphi) ? 0 : _nvalid(d.visphi, d.visphi_err)
    d.nt3amp  = isempty(d.t3amp)  ? 0 : _nvalid(d.t3amp,  d.t3amp_err)
    d.nt3phi  = isempty(d.t3phi)  ? 0 : _nvalid(d.t3phi,  d.t3phi_err)

    return d
end


# ─────────────────────────────────────────────────────────────────────────────
# Bootstrap result
# ─────────────────────────────────────────────────────────────────────────────

"""
    BootstrapResult

Outcome of [`bootstrap_fit`](@ref).

# Fields
- `samples` — `nboot × npar` matrix of per-replicate best-fit parameters
- `list_free_params` — parameter names, matching the columns of `samples`
- `x_opt` — best fit to the full dataset
- `median`, `sigma`, `sigma_minus`, `sigma_plus` — median and 16/84 percentile
  half-widths of the (masked) bootstrap distribution
- `covar`, `correlation` — parameter covariance and correlation matrices
- `chi2r` — reduced χ² of each replicate
- `mask` — replicates kept after sigma/χ² clipping
- `nblocks`, `mode`, `granularity` — resampling setup
- `nfailed` — replicates whose fit threw an error
"""
struct BootstrapResult
    samples          ::Matrix{Float64}
    list_free_params ::Vector{String}
    x_opt            ::Vector{Float64}
    median           ::Vector{Float64}
    sigma            ::Vector{Float64}
    sigma_minus      ::Vector{Float64}
    sigma_plus       ::Vector{Float64}
    covar            ::Matrix{Float64}
    correlation      ::Matrix{Float64}
    chi2r            ::Vector{Float64}
    mask             ::Vector{Bool}
    nblocks          ::Int
    mode             ::Symbol
    granularity      ::Symbol
    nfailed          ::Int
end

function Base.show(io::IO, r::BootstrapResult)
    @printf(io, "BootstrapResult: %d/%d replicates used  (mode=:%s, blocks=:%s, nblocks=%d)\n",
            count(r.mask), size(r.samples, 1), string(r.mode), string(r.granularity), r.nblocks)
    r.nfailed > 0 && @printf(io, "  %d replicate fit(s) failed\n", r.nfailed)
    for (i, p) in enumerate(r.list_free_params)
        @printf(io, "  %-20s = %10.5f  +%.5f -%.5f   (fit: %.5f)\n",
                p, r.median[i], r.sigma_plus[i], r.sigma_minus[i], r.x_opt[i])
    end
end


# ─────────────────────────────────────────────────────────────────────────────
# Bootstrap driver
# ─────────────────────────────────────────────────────────────────────────────

"""
    bootstrap_fit(model_dict, list_free_params, data; kwargs...) -> BootstrapResult

Estimate parameter uncertainties by nonparametric block bootstrap: refit the
model to `nboot` replicates in which the blocks of data (by default one per
MJD and telescope configuration) are resampled, and summarise the scatter of
the best-fit parameters.

Unlike the analytic covariance of `fit_model_lsqfit`, this does not assume the
quoted error bars are correct or uncorrelated: correlated calibration errors
and mis-stated error bars show up as extra scatter between replicates.  It
does assume the blocks are numerous and exchangeable — with a single snapshot
(a handful of blocks) the estimate is unreliable, and `perturb_data` based
Monte Carlo or the analytic errors are the better tool.

# Arguments
- `model_dict::Dict{String}` — flat parameter dictionary
- `list_free_params::Vector{String}` — names of the free parameters
- `data::OIdata` — interferometric data

# Keywords
- `nboot` — number of replicates (default 200; use ≥1000 for stable percentiles)
- `mode` — `:replacement` (default), `:halfsample` or `:weights`; see
  [`block_counts`](@ref) and [`block_weights`](@ref).  `:pmoired` reproduces
  PMOIRED's scheme and is biased low by about √2 — comparison only
- `granularity` — `:config` (default), `:epoch` or `:point`; see [`data_blocks`](@ref)
- `fitter` — `:lsqfit` (default, Levenberg-Marquardt with analytic Jacobian) or `:nlopt`
- `lb`, `ub`, `weights`, `vonmises`, `nB_workspace` — passed through to the fitter
- `method`, `maxeval` — NLopt settings (`fitter=:nlopt`)
- `maxIter` — LsqFit iterations (`fitter=:lsqfit`)
- `start_scatter` — perturb each replicate's starting point by this many σ of
  the full-data fit (default 0; PMOIRED uses the equivalent of 1). Non-zero
  values mix optimiser scatter into the uncertainty and are meant for
  comparison studies only.
- `sigma_clipping` — reject replicates beyond this many σ from the median
  (default `nothing`; PMOIRED's plots default to 4.5)
- `chi2r_max` — reject replicates whose reduced χ² exceeds this value
- `seed` — RNG seed; replicate `i` uses `Xoshiro(seed + i)` so results are
  reproducible regardless of thread scheduling
- `threaded` — run replicates on all available Julia threads (default `true`)
- `verb` — progress reporting (default `true`)

# Example
```julia
data  = readoifits("data/AlphaCenA.oifits")[1,1]
model = Dict{String,Any}("star,ldlin" => 8.0, "star,u" => 0.3, "star,f" => 1.0)
boot  = bootstrap_fit(model, ["star,ldlin", "star,u"], data;
                      lb=Dict("star,ldlin"=>5.0, "star,u"=>0.0),
                      ub=Dict("star,ldlin"=>12.0, "star,u"=>1.0),
                      weights=[1.0, 0.0, 0.0], nboot=500)
```
"""
function bootstrap_fit(model_dict       ::Dict{String},
                       list_free_params ::Vector{String},
                       data             ::OIdata;
    nboot          ::Int    = 200,
    mode           ::Symbol = :replacement,
    granularity    ::Symbol = :config,
    mjd_digits     ::Int    = 5,
    fitter         ::Symbol = :lsqfit,
    lb                      = Dict{String,Float64}(),
    ub                      = Dict{String,Float64}(),
    weights                 = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    vonmises       ::Bool   = false,
    nB_workspace            = nothing,
    method         ::Symbol = :LD_LBFGS,
    maxeval        ::Int    = 2000,
    maxIter        ::Int    = 200,
    start_scatter  ::Real   = 0.0,
    sigma_clipping          = nothing,
    chi2r_max               = nothing,
    seed                    = nothing,
    threaded       ::Bool   = true,
    verb           ::Bool   = true,
)
    fitter in (:lsqfit, :nlopt) ||
        error("bootstrap_fit: fitter must be :lsqfit or :nlopt (got :$fitter)")
    npar = length(list_free_params)
    nB   = something(nB_workspace, size(data.uv, 2))

    # ── 1. Fit the full dataset ──────────────────────────────────────────────
    verb && println("Fitting the full dataset...")
    x_opt, x_sigma = if fitter == :lsqfit
        f = fit_model_lsqfit(model_dict, list_free_params, data;
                             lb=lb, ub=ub, weights=weights, vonmises=vonmises,
                             nB_workspace=nB, maxIter=maxIter)
        f.x_opt, f.stderror
    else
        f = fit_model(model_dict, list_free_params, data;
                      lb=lb, ub=ub, weights=weights, vonmises=vonmises,
                      nB_workspace=nB, method=method, maxeval=maxeval)
        f.x_opt, fill(NaN, npar)
    end

    if start_scatter > 0 && any(isnan, x_sigma)
        # need per-parameter scales for the starting-point perturbation
        f = fit_model_lsqfit(model_dict, list_free_params, data;
                             lb=lb, ub=ub, weights=weights, vonmises=vonmises,
                             nB_workspace=nB, maxIter=maxIter)
        x_sigma = f.stderror
    end

    # ── 2. Build the blocks ──────────────────────────────────────────────────
    blocks = data_blocks(data; granularity=granularity, mjd_digits=mjd_digits)
    nblocks = length(blocks)
    verb && println(blocks)
    if nblocks < 10 && granularity !== :point
        @warn "Only $nblocks resampling blocks: bootstrap uncertainties will be " *
              "poorly determined. Consider granularity=:config, or the analytic " *
              "errors from fit_model_lsqfit / a perturb_data Monte Carlo."
    end

    # ── 3. Replicates ────────────────────────────────────────────────────────
    seed0    = seed === nothing ? rand(UInt32) : UInt64(seed)
    samples  = fill(NaN, nboot, npar)
    chi2r    = fill(NaN, nboot)
    nworkers = threaded ? min(Threads.nthreads(), nboot) : 1

    # One compiled model per worker: FlatModel carries mutable Hankel
    # workspaces, so it cannot be shared between threads.  Compiled serially
    # to keep the runtime-generated-function cache single-threaded.
    models = [parse_model(model_dict, list_free_params; nB_workspace=nB) for _ in 1:nworkers]

    done = Threads.Atomic{Int}(0)
    step = max(nboot ÷ 10, 1)

    function run_chunk(w::Int, chunk)
        fm = models[w]
        for i in chunk
            rng = Xoshiro(seed0 + i)
            try
                rd = resample_blocks(data, blocks; mode=mode, rng=rng)
                x0 = start_scatter > 0 ?
                     x_opt .+ start_scatter .* x_sigma .* randn(rng, npar) : x_opt
                if fitter == :lsqfit
                    r = fit_model_lsqfit(fm, x0, rd; lb=lb, ub=ub, weights=weights,
                                         vonmises=vonmises, maxIter=maxIter)
                    samples[i, :] = r.x_opt
                    chi2r[i]      = r.chi2r
                else
                    r = fit_model(fm, x0, rd; lb=lb, ub=ub, weights=weights,
                                  vonmises=vonmises, method=method, maxeval=maxeval)
                    samples[i, :] = r.x_opt
                    chi2r[i]      = r.chi2r
                end
            catch e
                e isa InterruptException && rethrow(e)
                # leave the row as NaN — counted as a failure below
            end
            if verb
                n = Threads.atomic_add!(done, 1) + 1
                (n % step == 0 || n == nboot) &&
                    @printf("\rBootstrap: %d/%d replicates", n, nboot)
            end
        end
    end

    verb && println("Running $nboot bootstrap fits on $nworkers thread(s)...")
    if nworkers == 1
        run_chunk(1, 1:nboot)
    else
        chunks = [w:nworkers:nboot for w in 1:nworkers]
        tasks  = [Threads.@spawn run_chunk(w, chunks[w]) for w in 1:nworkers]
        foreach(wait, tasks)
    end
    verb && println()

    # ── 4. Clipping ──────────────────────────────────────────────────────────
    mask = [all(isfinite, @view samples[i, :]) for i in 1:nboot]
    nfailed = count(!, mask)
    if chi2r_max !== nothing
        mask .&= (isfinite.(chi2r) .& (chi2r .<= chi2r_max))
    end
    if sigma_clipping !== nothing
        for j in 1:npar
            keep = copy(mask)
            for _ in 1:3
                any(keep) || break
                x   = @view samples[:, j]
                med = median(x[keep])
                s   = 0.5 * (quantile(x[keep], 0.84) - quantile(x[keep], 0.16))
                s > 0 || break
                keep = mask .& (abs.(x .- med) .<= sigma_clipping * s)
            end
            mask .&= keep
        end
    end
    count(mask) >= 2 ||
        error("bootstrap_fit: fewer than 2 usable replicates ($(count(mask))/$nboot); " *
              "check the fit setup, the bounds, or relax the clipping")

    # ── 5. Summary statistics ────────────────────────────────────────────────
    X = samples[mask, :]
    med  = [median(X[:, j])              for j in 1:npar]
    p16  = [quantile(X[:, j], 0.16)      for j in 1:npar]
    p84  = [quantile(X[:, j], 0.84)      for j in 1:npar]
    sm   = med .- p16
    sp   = p84 .- med
    sig  = 0.5 .* (sm .+ sp)
    C    = npar > 1 ? cov(X) : reshape([var(X[:, 1])], 1, 1)
    R    = npar > 1 ? cor(X) : ones(1, 1)

    if verb
        println()
        for (j, p) in enumerate(list_free_params)
            @printf("  %-20s = %10.5f  +%.5f -%.5f\n", p, med[j], sp[j], sm[j])
        end
    end

    return BootstrapResult(samples, list_free_params, x_opt, med, sig, sm, sp,
                           C, R, chi2r, mask, nblocks, mode, granularity, nfailed)
end
