# ─────────────────────────────────────────────────────────────────────────────
# cvis_to_chi2_noalloc — append this to src/chi2_flat.jl.
#
# WHY THIS EXISTS ALONGSIDE `cvis_to_chi2_f`
#
# `cvis_to_chi2_f` already gives chi2 without the gradient, and is the right
# choice almost everywhere — it also handles the flux term and does proper
# `_chi2_eltype` precision handling.  But it routes through `_chi2_terms`, which
# allocates on nearly every line:
#
#     v2_model    = abs2.(V[data.indx_v2])        # gather + broadcast temporary
#     V1          = V[data.indx_t3_1]             # gather
#     t3amp_model = abs.(t3)                      # temporary
#     dphi        = _mod360(t3phi_model .- ...)   # two more
#
# That is ~10 heap allocations per call, several of them nuv-length arrays.  The
# SQUEEZE sampler calls this ONCE PER PROPOSAL — millions of times per run — so
# the allocation and GC pressure dominate the arithmetic.  Measured on the
# original OITOOLS: 0 bytes and ~9x faster than `cvis_to_chi2_fg`.
#
# So: use `cvis_to_chi2_f` for ordinary work; use this in inner sampler loops.
# The two must agree to round-off, which `test_squeeze.jl` checks.
#
# The chi2 accumulator is Float64 regardless of the data precision `T`: it is a
# sum of ~ndata positive terms, so Float64 costs one scalar and removes
# summation error from the Metropolis acceptance test.
# ─────────────────────────────────────────────────────────────────────────────

@inline _mod360s(x::Float64) = mod(mod(x + 180.0, 360.0) + 360.0, 360.0) - 180.0

const _RAD2DEG = 180.0 / π
const _DEG2RAD = π / 180.0

"""
    _weights_tuple(weights) -> NTuple{7,Float64}

Pad/truncate a weight vector to a stack-allocated 7-tuple, so the sampler's inner
loop never re-allocates (`_pad_weights` returns a heap `Vector`).
"""
function _weights_tuple(weights::AbstractVector{<:Real})
    n = length(weights)
    return ntuple(i -> i <= n ? Float64(weights[i]) : 0.0, Val(7))
end

function _cvis_to_chi2(V::AbstractVector{<:Complex}, data::OIdata,
                       w::NTuple{7,Float64}, vonmises::Bool)
    chi2 = 0.0

    # ── V2 ────────────────────────────────────────────────────────────────────
    if w[1] > 0 && data.nv2 > 0
        s = 0.0
        idx = data.indx_v2
        @inbounds for i in eachindex(idx)
            d = Float64(abs2(V[idx[i]])) - Float64(data.v2[i])
            e = Float64(data.v2_err[i])
            s += d * d / (e * e)
        end
        chi2 += w[1] * s
    end

    # ── T3 (one pass computes the triple product for both amp and phi) ────────
    do_t3amp = w[2] > 0 && data.nt3amp > 0
    do_t3phi = w[3] > 0 && data.nt3phi > 0
    if do_t3amp || do_t3phi
        s_amp = 0.0
        s_phi = 0.0
        i1 = data.indx_t3_1; i2 = data.indx_t3_2; i3 = data.indx_t3_3
        @inbounds for i in eachindex(i1)
            t3 = V[i1[i]] * V[i2[i]] * V[i3[i]]
            if do_t3amp
                d = Float64(abs(t3)) - Float64(data.t3amp[i])
                e = Float64(data.t3amp_err[i])
                s_amp += d * d / (e * e)
            end
            if do_t3phi
                phi = Float64(angle(t3)) * _RAD2DEG
                if !vonmises
                    dphi = _mod360s(phi - Float64(data.t3phi[i]))
                    e = Float64(data.t3phi_err[i])
                    s_phi += dphi * dphi / (e * e)
                else
                    dphi_rad = (phi - Float64(data.t3phi[i])) * _DEG2RAD
                    s_phi += -2.0 * Float64(data.t3phi_vonmises_err[i]) * cos(dphi_rad) +
                             Float64(data.t3phi_vonmises_chi2_offset[i])
                end
            end
        end
        chi2 += w[2] * s_amp + w[3] * s_phi
    end

    # ── Absolute VISAMP / VISPHI ──────────────────────────────────────────────
    do_visamp = w[4] > 0 && data.nvisamp > 0
    do_visphi = w[5] > 0 && data.nvisphi > 0
    if do_visamp || do_visphi
        s_amp = 0.0
        s_phi = 0.0
        idx = data.indx_vis
        @inbounds for i in eachindex(idx)
            Vv = V[idx[i]]
            if do_visamp
                d = Float64(abs(Vv)) - Float64(data.visamp[i])
                e = Float64(data.visamp_err[i])
                s_amp += d * d / (e * e)
            end
            if do_visphi
                dphi = _mod360s(Float64(angle(Vv)) * _RAD2DEG - Float64(data.visphi[i]))
                e = Float64(data.visphi_err[i])
                s_phi += dphi * dphi / (e * e)
            end
        end
        chi2 += w[4] * s_amp + w[5] * s_phi
    end

    return chi2
end

"""
    cvis_to_chi2_noalloc(V, data; weights=[1,1,1,0,0,0,0], vonmises=false)

Weighted chi2 of model complex visibilities `V` against `data`, computed **without
allocating**.

Returns the same value as [`cvis_to_chi2_f`](@ref) (and as
`cvis_to_chi2_fg(V, data; weights, vonmises)[1]`) to round-off, but fuses every
observable into a single pass with no gathered temporaries.  Intended for inner
sampler loops that evaluate chi2 millions of times — see `reconstruct_squeeze`.
Prefer [`cvis_to_chi2_f`](@ref) for ordinary use: it is clearer and additionally
handles the flux term (`weights[6]`), which this does not.
"""
function cvis_to_chi2_noalloc(V::AbstractVector{<:Complex}, data::OIdata;
                              weights::AbstractVector{<:Real} = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
                              vonmises::Bool = false)
    return _cvis_to_chi2(V, data, _weights_tuple(weights), vonmises)
end
