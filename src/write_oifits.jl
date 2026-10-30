# Writing an OIdata back out as OIFITS v2.
#
# Built table by table rather than by editing the file the data came from: OIFITS.jl cannot
# write back several readable files (a missing CATEGORY in OI_TARGET, FOV in OI_ARRAY, a
# CFITSIO 402 on others), while a dataset assembled from scratch -- as `simulate` assembles
# its own -- writes. What OIdata does not keep (the target's coordinates, the stations'
# positions) is copied from the original file when it can still be read.

"""
    write_oifits(path, data; template = :auto, targetname = "", overwrite = false) -> path

Write an `OIdata`, or an array of them (the bins `readoifits` returns), to an OIFITS v2 file
that `readoifits` reads back to the same points.

- One OI_WAVELENGTH table per instrument the points were read from (`data.tables`), named
  by its INSNAME and holding the distinct wavelengths of that instrument's points, every bin
  included; points whose table is unknown go to an instrument named `"OITOOLS"`.
- Per instrument, one OI_VIS2, OI_T3 and OI_FLUX table, and one OI_VIS table per AMPTYP /
  PHITYP pair. One row per epoch, stations and baseline (or triangle), spanning the
  instrument's channels; a channel without a point is written NaN and flagged, so it is
  dropped again on reading. OI_FLUX likewise, per epoch and station.
- OI_TARGET is the original file's target named `targetname`, or its first target when
  `targetname` is empty, and OI_ARRAY the data's own stations, their diameters and positions
  copied from the original file where the names match; both from `template` (by default the
  file the data was read from), or minimal when it cannot be read.

An existing file is refused unless `overwrite = true`.
"""
function write_oifits(path::AbstractString, data; template = :auto, targetname::AbstractString = "",
                      overwrite::Bool = false)
    bins = data isa OIdata ? [data] : [d for d in vec(collect(data))]
    isempty(bins) && throw(ArgumentError("no data to write"))
    d1 = first(bins)
    (!overwrite && isfile(path)) &&
        throw(ArgumentError("$(path) exists; pass overwrite = true to replace it"))
    tpl = template === :auto ? (isfile(d1.filename) ? d1.filename : nothing) : template
    hdr = tpl === nothing ? nothing : (try load_oifits(tpl) catch; nothing end)

    # ── instruments and their wavelength grids ──────────────────────────────────────
    # A point's instrument is its table's INSNAME; instruments are written in the order their
    # tables were read.
    insof(d, t) = (t == 0 || isempty(d.tables[t].insname)) ? "OITOOLS" : d.tables[t].insname
    instruments = String[]
    for d in bins, tab in (d.vis_tab, d.v2_tab, d.t3_tab, d.flux_tab), t in unique(tab)
        n = insof(d, t); n in instruments || push!(instruments, n)
    end
    lamkey(l) = round(Float64(l); sigdigits = 7)
    band = Dict(n => Dict{Float64,Tuple{Float64,Float64}}() for n in instruments)  # key => (λ, Δλ) as read
    for d in bins
        for (l, dl, tab) in ((d.v2_lam, d.v2_dlam, d.v2_tab), (d.t3_lam, d.t3_dlam, d.t3_tab),
                             (d.vis_lam, d.vis_dlam, d.vis_tab), (d.flux_lam, d.flux_dlam, d.flux_tab))
            for i in eachindex(l)
                get!(band[insof(d, tab[i])], lamkey(l[i]), (Float64(l[i]), Float64(dl[i])))
            end
        end
    end
    isempty(instruments) && throw(ArgumentError("the data hold no V², T3, VIS or flux point"))
    lamkeys = Dict(n => sort(collect(keys(band[n]))) for n in instruments)
    chan    = Dict(n => Dict(w => k for (k, w) in enumerate(lamkeys[n])) for n in instruments)

    # ── stations ─────────────────────────────────────────────────────────────────────
    sta_index = Int.(d1.sta_index)
    sta_name  = String.(d1.sta_name)
    tel_name  = length(d1.tel_name) == length(sta_name) ? String.(d1.tel_name) : copy(sta_name)
    diam = zeros(length(sta_index)); xyz = zeros(3, length(sta_index))
    arrname = "OITOOLS"
    if hdr !== nothing && !isempty(hdr.array)
        arrname = String(hdr.array[1].arrname)
        for a in hdr.array, (k, n) in enumerate(sta_name)
            j = findfirst(==(n), String.(a.sta_name))
            j === nothing && continue
            diam[k] = Float64(a.diameter[j])
            xyz[:, k] .= Float64.(a.staxyz[:, j])
        end
    end

    # ── rows: one per epoch, stations and baseline, across the channels ─────────────
    m4(x) = round(Float64(x); digits = 4)          # metres, to group a baseline's channels
    uvm(d, k, l) = (Float64(d.uv[1, k]) * Float64(l), Float64(d.uv[2, k]) * Float64(l))

    # A row per key; a key whose channel is already taken (the file repeated a baseline at one
    # epoch) opens a further row, so no point overwrites another.
    function gather(items)        # items: (key, channel, values tuple, row info)
        rows = Dict{Any,Vector{Int}}(); used = Set{Tuple{Int,Int}}()
        info = Any[]; where = Vector{Int}(undef, length(items))
        for (n, it) in enumerate(items)
            cands = get!(rows, it[1], Int[])
            r = findfirst(r -> (r, it[2]) ∉ used, cands)
            if r === nothing
                push!(info, it[4]); push!(cands, length(info)); r = length(info)
            else
                r = cands[r]
            end
            push!(used, (r, it[2])); where[n] = r
        end
        return where, info
    end
    nanmat(nw, nr) = fill(NaN, nw, nr)
    trues_(nw, nr) = fill(true, nw, nr)
    # One half of an amplitude/phase pair that the file did not carry is held by `readoifits` as
    # 0 ± Inf and not counted (GRAVITY's T3AMP and VISPHI, for instance). It is written back as
    # the NaN it was, not as a value with an infinite error, which a reader would count.
    half(v, e) = isfinite(e) ? (v, e) : (NaN, NaN)

    # Items per table to write: V2, T3 and flux per instrument, VIS per instrument, AMPTYP and
    # PHITYP (from the point's table, else the dataset's).
    v2items = Dict{String,Vector{Any}}(); t3items = Dict{String,Vector{Any}}()
    fluxitems = Dict{String,Vector{Any}}(); visitems = Dict{Tuple{String,String,String},Vector{Any}}()
    additem!(dict, key, it) = push!(get!(dict, key, Any[]), it)
    vistyp(d, t) = t == 0 ? (d.amptyp, d.phityp) : (d.tables[t].amptyp, d.tables[t].phityp)
    for d in bins
        for i in eachindex(d.v2)
            n = insof(d, d.v2_tab[i])
            k = d.indx_v2[i]; u, v = uvm(d, k, d.v2_lam[i])
            s = Int.(d.v2_sta_index[:, i])
            additem!(v2items, n, ((d.v2_mjd[i], s[1], s[2], m4(u), m4(v)), chan[n][lamkey(d.v2_lam[i])],
                              (d.v2[i], d.v2_err[i]), (d.v2_mjd[i], u, v, s)))
        end
        for i in eachindex(d.visamp)
            n = insof(d, d.vis_tab[i])
            k = d.indx_vis[i]; u, v = uvm(d, k, d.vis_lam[i])
            s = Int.(d.vis_sta_index[:, i])
            additem!(visitems, (n, vistyp(d, d.vis_tab[i])...),
                 ((d.vis_mjd[i], s[1], s[2], m4(u), m4(v)), chan[n][lamkey(d.vis_lam[i])],
                  (half(d.visamp[i], d.visamp_err[i])..., half(d.visphi[i], d.visphi_err[i])...),
                  (d.vis_mjd[i], u, v, s)))
        end
        for i in eachindex(d.t3phi)
            n = insof(d, d.t3_tab[i])
            l = d.t3_lam[i]
            u1, v1 = uvm(d, d.indx_t3_1[i], l); u2, v2 = uvm(d, d.indx_t3_2[i], l)
            s = Int.(d.t3_sta_index[:, i])
            additem!(t3items, n, ((d.t3_mjd[i], s..., m4(u1), m4(v1), m4(u2), m4(v2)), chan[n][lamkey(l)],
                              (half(d.t3amp[i], d.t3amp_err[i])..., half(d.t3phi[i], d.t3phi_err[i])...),
                              (d.t3_mjd[i], u1, v1, u2, v2, s)))
        end
        for i in eachindex(d.flux)
            n = insof(d, d.flux_tab[i])
            s = Int(d.flux_sta_index[i])
            additem!(fluxitems, n, ((d.flux_mjd[i], s), chan[n][lamkey(d.flux_lam[i])],
                                (d.flux[i], d.flux_err[i]), (d.flux_mjd[i], s)))
        end
    end

    mjd_all = Float64[]
    for it in (v2items, visitems, t3items, fluxitems), (_, xs) in it, x in xs
        push!(mjd_all, Float64(x[4][1]))
    end
    date_obs = string(Date(1858, 11, 17) + Day(floor(Int, minimum(mjd_all))))

    # ── the tables ───────────────────────────────────────────────────────────────────
    t = hdr === nothing || isempty(hdr.target) ? nothing :
        hdr.target[something(findfirst(e -> String(e.target) == targetname, hdr.target.list), 1)]
    tgt = OIFITS.OI_TARGET([t === nothing ?
        OIFITS.OITargetEntry(1, "unknown", 0.0, 0.0, 2000.0, 0.0, 0.0, 0.0, "UNKNOWN",
                             "OPTICAL", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "UNKNOWN", "SCI") :
        OIFITS.OITargetEntry(1, String(t.target), t.raep0, t.decep0, t.equinox, t.ra_err,
                             t.dec_err, t.sysvel, String(t.veltyp), String(t.veldef), t.pmra,
                             t.pmdec, t.pmra_err, t.pmdec_err, t.parallax, t.para_err,
                             String(t.spectyp), "SCI")]; revn = 2)

    arr = OIFITS.OI_ARRAY(undef)
    arr.revn = 2; arr.arrname = arrname; arr.frame = "GEOCENTRIC"
    arr.arrayx = 0.0; arr.arrayy = 0.0; arr.arrayz = 0.0
    arr.tel_name = tel_name; arr.sta_name = sta_name; arr.sta_index = sta_index
    arr.diameter = diam; arr.staxyz = xyz
    arr.fov = zeros(length(sta_index)); arr.fovtype = fill("RADIUS", length(sta_index))

    blocks = Any[tgt, arr]
    for n in instruments
        ins = OIFITS.OI_WAVELENGTH(undef)
        ins.revn = 2; ins.insname = n
        ins.eff_wave = [band[n][w][1] for w in lamkeys[n]]
        ins.eff_band = [band[n][w][2] for w in lamkeys[n]]
        push!(blocks, ins)
    end
    nwof(n) = length(lamkeys[n])

    for n in instruments
        haskey(v2items, n) || continue
        items = v2items[n]; nw = nwof(n)
        rows, info = gather(items); nr = length(info)
        val = nanmat(nw, nr); err = nanmat(nw, nr); flg = trues_(nw, nr)
        for ((_, c, x, _), r) in zip(items, rows)
            val[c, r] = x[1]; err[c, r] = x[2]; flg[c, r] = false
        end
        db = OIFITS.OI_VIS2(undef)
        db.revn = 2; db.date_obs = date_obs; db.arrname = arrname; db.insname = n
        db.target_id = ones(Int, nr); db.time = zeros(nr); db.int_time = ones(nr)
        db.mjd = [Float64(x[1]) for x in info]
        db.ucoord = [x[2] for x in info]; db.vcoord = [x[3] for x in info]
        db.sta_index = reduce(hcat, [x[4] for x in info])
        db.vis2data = val; db.vis2err = err; db.flag = flg
        push!(blocks, db)
    end
    for key in sort(collect(keys(visitems)); by = k -> (findfirst(==(k[1]), instruments), k[2], k[3]))
        n, amptyp, phityp = key
        items = visitems[key]; nw = nwof(n)
        rows, info = gather(items); nr = length(info)
        amp = nanmat(nw, nr); ae = nanmat(nw, nr); phi = nanmat(nw, nr); pe = nanmat(nw, nr)
        flg = trues_(nw, nr)
        for ((_, c, x, _), r) in zip(items, rows)
            amp[c, r] = x[1]; ae[c, r] = x[2]; phi[c, r] = x[3]; pe[c, r] = x[4]; flg[c, r] = false
        end
        db = OIFITS.OI_VIS(undef)
        db.revn = 2; db.date_obs = date_obs; db.arrname = arrname; db.insname = n
        db.amptyp = amptyp; db.phityp = phityp
        db.target_id = ones(Int, nr); db.time = zeros(nr); db.int_time = ones(nr)
        db.mjd = [Float64(x[1]) for x in info]
        db.ucoord = [x[2] for x in info]; db.vcoord = [x[3] for x in info]
        db.sta_index = reduce(hcat, [x[4] for x in info])
        db.visamp = amp; db.visamperr = ae; db.visphi = phi; db.visphierr = pe; db.flag = flg
        push!(blocks, db)
    end
    for n in instruments
        haskey(t3items, n) || continue
        items = t3items[n]; nw = nwof(n)
        rows, info = gather(items); nr = length(info)
        amp = nanmat(nw, nr); ae = nanmat(nw, nr); phi = nanmat(nw, nr); pe = nanmat(nw, nr)
        flg = trues_(nw, nr)
        for ((_, c, x, _), r) in zip(items, rows)
            amp[c, r] = x[1]; ae[c, r] = x[2]; phi[c, r] = x[3]; pe[c, r] = x[4]; flg[c, r] = false
        end
        db = OIFITS.OI_T3(undef)
        db.revn = 2; db.date_obs = date_obs; db.arrname = arrname; db.insname = n
        db.target_id = ones(Int, nr); db.time = zeros(nr); db.int_time = ones(nr)
        db.mjd = [Float64(x[1]) for x in info]
        db.u1coord = [x[2] for x in info]; db.v1coord = [x[3] for x in info]
        db.u2coord = [x[4] for x in info]; db.v2coord = [x[5] for x in info]
        db.sta_index = reduce(hcat, [x[6] for x in info])
        db.t3amp = amp; db.t3amperr = ae; db.t3phi = phi; db.t3phierr = pe; db.flag = flg
        push!(blocks, db)
    end
    for n in instruments
        haskey(fluxitems, n) || continue
        items = fluxitems[n]; nw = nwof(n)
        rows, info = gather(items); nr = length(info)
        val = nanmat(nw, nr); err = nanmat(nw, nr); flg = trues_(nw, nr)
        for ((_, c, x, _), r) in zip(items, rows)
            val[c, r] = x[1]; err[c, r] = x[2]; flg[c, r] = false
        end
        db = OIFITS.OI_FLUX(undef)
        db.revn = 1; db.date_obs = date_obs; db.arrname = arrname; db.insname = n
        db.calstat = d1.flux_calibrated ? "C" : "U"; db.fov = 0.0; db.fovtype = "RADIUS"
        db.target_id = ones(Int, nr); db.int_time = ones(nr)
        db.mjd = [Float64(x[1]) for x in info]
        db.sta_index = [x[2] for x in info]
        db.fluxdata = val; db.fluxerr = err; db.flag = flg
        push!(blocks, db)
    end

    OIFITS.write(path, OIFITS.OIDataSet(blocks...); overwrite = true)
    return path
end
