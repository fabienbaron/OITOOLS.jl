# Spectral vectors: the points of one dataset that a single row of an OIFITS table holds — one
# baseline (triangle, telescope) at one time, across the wavelengths. The measurement model
# references the model's phases and amplitudes per vector, as the pipeline referenced the data
# (`PLAN_C1.md`).

"""
    SpectralGroups

The spectral vectors of a dataset, each a vector of indices into the observable's own arrays
(`visphi`, `v2`, `t3phi`, `flux`), ordered by wavelength: `vis`, `v2`, `t3`, `flux`.
"""
struct SpectralGroups
    vis  :: Vector{Vector{Int}}
    v2   :: Vector{Vector{Int}}
    t3   :: Vector{Vector{Int}}
    flux :: Vector{Vector{Int}}
end

# Indices grouped by `key(i)` in first-seen order, each group ordered by wavelength `lam`.
function _groups(n, key, lam)
    order = Dict{Any,Int}(); out = Vector{Int}[]
    for i in 1:n
        k = key(i)
        g = get!(order, k) do
            push!(out, Int[]); length(out)
        end
        push!(out[g], i)
    end
    for g in out
        sort!(g; by = i -> lam[i])
    end
    return out
end

"""
    spectral_groups(data::OIdata) -> SpectralGroups

The dataset's spectral vectors: VIS and V² points sharing their table, their two stations and
their MJD; T3 points their three stations; flux points their station (0 for a calibrated
spectrum) — each vector ordered by wavelength.
"""
function spectral_groups(d::OIdata)
    vis  = _groups(length(d.vis_mjd), i -> (d.vis_tab[i], d.vis_sta_index[1, i], d.vis_sta_index[2, i], d.vis_mjd[i]), d.vis_lam)
    v2   = _groups(length(d.v2_mjd), i -> (d.v2_tab[i], d.v2_sta_index[1, i], d.v2_sta_index[2, i], d.v2_mjd[i]), d.v2_lam)
    t3   = _groups(length(d.t3_mjd), i -> (d.t3_tab[i], d.t3_sta_index[1, i], d.t3_sta_index[2, i],
                                           d.t3_sta_index[3, i], d.t3_mjd[i]), d.t3_lam)
    flux = _groups(length(d.flux_mjd), i -> (d.flux_tab[i], d.flux_sta_index[i], d.flux_mjd[i]), d.flux_lam)
    return SpectralGroups(vis, v2, t3, flux)
end
