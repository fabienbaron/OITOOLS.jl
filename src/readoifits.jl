using FITSIO, OIFITS, Statistics, SparseArrays, Crayons

# TO DO : * use the splat operator to simplify things
#         * add differential and complex visibilities
#         * why are we negating u coordinates ?
#         * check outer=[1,8]
#         * sort indexing of uv points when t3amp and t3phi differ
#         * define booleans "use_v2, use_t3amp" and apply to all file
#         * check if t3_uv_mjd and similar are necessary

using NearestNeighbors
function rm_redundance_kdtree(uv,uvtol)
  #@inbounds uv[:,findall(uv[1,:].<0)] *= -1; # change any (-u,v) -> (u,-v)
  # Note: good idea but we would need to complex conjugate the data point in this case
  indx_redundance = collect(1:size(uv,2));
  kdtree = KDTree(uv);
  for value in indx_redundance
    redundance = inrange(kdtree,uv[:,value],uvtol);
    @inbounds indx_redundance[redundance] .= minimum(redundance);
  end
  tokeep = unique(indx_redundance) # we only need the tokeep points in the uv plane'
  indx_red_conv = indexin(indx_redundance, tokeep)
  return indx_red_conv, tokeep  ;
end

function tablemerge(tabtomerge)
    return vcat([vec(tabtomerge[i]) for i=1:length(tabtomerge)]...);
end


mutable struct OIdata
    # Complex visibilities
    visamp::Array{Float64,1}
    visamp_err::Array{Float64,1}
    visphi::Array{Float64,1}
    visphi_err::Array{Float64,1}
    vis_baseline::Array{Float64,1}
    vis_mjd::Array{Float64,1}
    vis_lam::Array{Float64,1}
    vis_dlam::Array{Float64,1}
    vis_flag::Array{Bool,1}
    # V2
    v2::Array{Float64,1}
    v2_err::Array{Float64,1}
    v2_baseline::Array{Float64,1}
    v2_mjd::Array{Float64,1}
    mean_mjd::Float64
    v2_lam::Array{Float64,1}
    v2_dlam::Array{Float64,1}
    v2_flag::Array{Bool,1}
    # T3
    t3amp::Array{Float64,1}
    t3amp_err::Array{Float64,1}
    t3phi::Array{Float64,1}
    t3phi_err::Array{Float64,1}
    t3phi_vonmises_err::Array{Float64,1}
    t3phi_vonmises_chi2_offset::Array{Float64,1}
    t3_baseline::Array{Float64,1}
    t3_maxbaseline::Array{Float64,1}
    t3_mjd::Array{Float64,1}
    t3_lam::Array{Float64,1}
    t3_dlam::Array{Float64,1}
    t3_flag::Array{Bool,1}
    #T4
#    t4amp::Array{Float64,1}
#    t4amp_err::Array{Float64,1}
#    t4phi::Array{Float64,1}
#    t4phi_err::Array{Float64,1}
#    t4_baseline::Array{Float64,1}
#    t4_maxbaseline::Array{Float64,1}
#    t4_mjd::Array{Float64,1}
#    t4_lam::Array{Float64,1}
#    t4_dlam::Array{Float64,1}
#    t4_flag::Array{Bool,1}
    #OIFlux
    flux::Array{Float64,1}
    flux_err::Array{Float64,1}
    flux_mjd::Array{Float64,1}
    flux_lam::Array{Float64,1}
    flux_dlam::Array{Float64,1}
    flux_flag::Array{Bool,1}
    flux_sta_index::Array{Int64,1}
    #UV coverage
    uv::Array{Float64,2}
    uv_lam::Array{Float64,1}
    uv_dlam::Array{Float64,1}
    uv_mjd::Array{Float64,1}
    uv_baseline::Array{Float64,1}
    # Data product sizes
    nflux::Int64
    nvisamp::Int64
    nvisphi::Int64
    nv2::Int64
    nt3amp::Int64
    nt3phi::Int64
    # nt4amp::Int64
    # nt4phi::Int64
    nuv::Int64
    # Indexing logic
    indx_vis::Array{Int64,1}
    indx_v2::Array{Int64,1}
    indx_t3_1::Array{Int64,1}
    indx_t3_2::Array{Int64,1}
    indx_t3_3::Array{Int64,1}
    # indx_t4_1::Array{Int64,1}
    # indx_t4_2::Array{Int64,1}
    # indx_t4_3::Array{Int64,1}
    # indx_t4_4::Array{Int64,1}
    sta_name::Array{String,1}
    tel_name::Array{String,1}
    sta_index::Array{Int64,1}
    vis_sta_index::Array{Int64,2}
    v2_sta_index::Array{Int64,2}
    t3_sta_index::Array{Int64,2}
    # t4_sta_index::Array{Int64,2}
    filename::String
end


function Base.display(data::OIdata)
println("Mean MJD: $(data.mean_mjd)")
println("Wavelength range: $(minimum(data.uv_lam)) - $(maximum(data.uv_lam))")
println("nflux: $(data.nflux) | nuv: $(data.nuv) | nvisamp: $(data.nvisamp) | nvisphi: $(data.nvisphi) | nv2: $(data.nv2) | nt3amp: $(data.nt3amp) | nt3phi: $(data.nt3phi)")
end

function Base.display(data::Array{OIdata})
println("Original data file: $(data[1].filename)")
println("Number of wavelength bins: $(size(data,1))")
println("Number of time/epoch bins: $(size(data,2))")
printcolors=[196, 166, 208, 220, 148 , 112, 28, 20, 92, 165]
if  size(data,1) == (1,1)
    display(data[1,1])
elseif (size(data,1) < 11) && (size(data,2) == 1)
    for i=1:size(data,1) 
    print(Crayon(foreground = printcolors[i], bold = true), "Wavelength: $i/$(size(data,1))\n")
    display(data[i,1])
    end
end
end


function set_data_filter(data::OIdata; wav_range::Union{Array{Float64,1}, Array{Array{Float64,1}}} = [-1.0, 1e99], mjd_range::Union{Array{Float64,1}, Array{Array{Float64,1}}} = [-1.0, 1e99], baseline_range::Array{Float64}=[0,1e99],
    filter_bad_data = false, filter_vis = true, filter_v2 = true, filter_t3amp = true, filter_t3phi = true,
    cutoff_minv2 = -1, cutoff_maxv2 = 2.0, cutoff_mint3amp = -1.0, cutoff_maxt3amp = 1.5, special_filter_diffvis=false, force_full_vis = false, force_full_t3 = false, filter_v2_snr_threshold=0.01, uv_bad=Int64[], filter_visphi=false, filter_visamp=false)

    if !isempty(wav_range)
        if typeof(wav_range) == Array{Float64,1}
            wav_range = [wav_range]
        end
    end

    if !isempty(mjd_range)
        if typeof(mjd_range) == Array{Float64,1}
            mjd_range = [mjd_range]
        end
    end

    # Select points (to keep) in the uv plane
    # one can directly deselect/disable uv points by setting the list of bad uv point in uv_bad
    # One can disable the selection of e.g. V2 observables by setting filter_v2 = false
    # Can select wavelengths in a range
    # Can select MJDs in a range
    use_visphi = (data.nvisphi > 0) && (filter_visphi == true)
    use_visamp = (data.nvisamp > 0) && (filter_visamp == true)
    use_vis = use_visphi || use_visamp
    use_v2 = (data.nv2 > 0) && (filter_v2 == true)
    use_t3amp = (data.nt3amp > 0) && (filter_t3amp == true)
    use_t3phi = (data.nt3phi > 0) && (filter_t3phi == true)
    use_t3 = use_t3phi || use_t3amp
#    use_t4amp = (data.nt4amp > 0) && (filter_t4amp == true)
#    use_t4phi = (data.nt4phi > 0) && (filter_t4phi == true)
#    use_t4 = use_t4phi || use_t4amp

    # Bad points are easier to track/concatenate/union than good points are to intersect
    vis_bad = Int64[]
    v2_bad = Int64[]
    t3_bad = Int64[]
#    t4_bad = Int64[]

    if filter_bad_data == true
        if use_vis
            visamp_good =  (.!isnan.(data.visamp )) .& (.!isnan.(data.visamp_err )) .& (data.visamp_err.>0.0)
            visphi_good =  (.!isnan.(data.visphi ) ).& (.!isnan.(data.visphi_err )) .& (data.visphi_err.>0.0)
            vis_good = []
            if force_full_vis == false
                vis_good = findall(.!data.vis_flag .& (visamp_good .| visphi_good) )
            else
                vis_good = findall(.!data.vis_flag .& (visamp_good .& visphi_good) )
            end

            if special_filter_diffvis==true # if we use the diff vis filter, we need to do this after all the data has been read in
                # hence here we won't actually filter anything
                vis_good = findall( data.vis_flag.!=2)
            end
            vis_bad = setdiff(collect(1:length(data.vis_flag)), vis_good)
        end
        if use_v2
            # Filter OBVIOUSLY bad V2 data
            v2_good = findall(  (data.v2_flag.==false) .& (data.v2_err.>0.0).& (data.v2_err.<1.0) .& (data.v2.>cutoff_minv2) .& (data.v2.<cutoff_maxv2) .& .!isnan.(data.v2) .& .!isnan.(data.v2_err) .& (abs.(data.v2./data.v2_err).>filter_v2_snr_threshold))
            v2_bad  = setdiff(collect(1:length(data.v2_flag)), v2_good)
        end

        if use_t3
            t3amp_good =  (.!isnan.(data.t3amp )) .& (.!isnan.(data.t3amp_err )) .& (data.t3amp_err.>0.0)
            t3phi_good =  (.!isnan.(data.t3phi )) .& (.!isnan.(data.t3phi_err )) .& (data.t3phi_err.>0.0)

            #if force_full_t3 is set to "true", then we require both t3amp and t3phi to be defined
            t3_good = Int64[]
            if force_full_t3 == false
                t3_good = findall(.!data.t3_flag .& (t3amp_good .| t3phi_good))
            else
                t3_good = findall(.!data.t3_flag .& (t3amp_good .& t3phi_good) .& (data.t3amp.>cutoff_mint3amp) .& (data.t3amp.<cutoff_maxt3amp) )
            end
            t3_bad  = setdiff(collect(1:length(data.t3_flag)), t3_good)
        end

        # if use_t4 == true
        #     t4amp_good =  (.!isnan.(data.t4amp )) .& (.!isnan.(data.t4amp_err )) .& (data.t4amp_err.>0.0)
        #     t4phi_good =  (.!isnan.(data.t4phi )) .& (.!isnan.(data.t4phi_err )) .& (data.t4phi_err.>0.0)
        #     # t4_good = []
        #     force_full_t4 = true;
        #     if force_full_t4 == false
        #         t4_good = findall(.!data.t4_flag .& (t4amp_good .| t4phi_good) )
        #     else
        #         t4_good = findall(.!data.t4_flag .& (t4amp_good .& t4phi_good) )
        #     end
        #     t4_bad  = setdiff(collect(1:length(data.t4_flag)), t4_good)
        # end
    end

    # Filtering the uv plane (this will filter the observables too)
    # We define the uv indexes we keep as uv_good
    if baseline_range != [0,1e99]
        uv_bad = union(uv_bad, findall(data.uv_baseline.<baseline_range[1]), findall(data.uv_baseline.>baseline_range[2]))
    end
    uv_good = setdiff(collect(1:data.nuv), uv_bad) # select all as good except in the uv_bad list
    if !isempty(mjd_range)
        uv_good = intersect(uv_good, vcat([findall( mjd_range[i][1]  .<= data.uv_mjd .<=  mjd_range[i][2]) for i=1:length(mjd_range)]...))
    end
    if !isempty(wav_range)
        uv_good = intersect(uv_good, vcat([findall( wav_range[i][1]  .<= data.uv_lam .<=  wav_range[i][2]) for i=1:length(wav_range)]...))
    end


    uv_bad = setdiff(collect(1:data.nuv), uv_good)

    if (data.nvisamp>0 || data.nvisphi >0)
        vis_bad = union(vis_bad, findall([data.indx_vis[i] ∉ uv_good for i=1:length(data.indx_vis)]))
    end
    if data.nv2>0
        v2_bad = union(v2_bad, findall([data.indx_v2[i] ∉ uv_good for i=1:length(data.indx_v2)]))
    end
    if (data.nt3amp>0 || data.nt3phi>0)
        t3_bad = union(t3_bad, findall([((data.indx_t3_1[i] ∉ uv_good) || (data.indx_t3_2[i] ∉ uv_good) ||  (data.indx_t3_3[i] ∉ uv_good)) for i=1:length(data.indx_t3_1)]))
    end
    # if (data.nt4amp>0 || data.nt4phi>0)
    #     t4_bad = union(t4_bad, findall([((data.indx_t4_1[i] ∉ uv_good) || (data.indx_t4_2[i] ∉ uv_good) ||  (data.indx_t4_3[i] ∉ uv_good) || (data.indx_t4_4[i] ∉ uv_good)) for i=1:length(data.indx_t3_1)]))
    # end
    return [uv_bad, vis_bad, v2_bad, t3_bad] #[uv_bad, vis_bad, v2_bad, t3_bad, t4_bad]
end


function filter_data(data_in::OIdata, indexes_to_discard = Int64[])
    data = deepcopy(data_in)
    good_uv_vis = Int64[]
    good_uv_v2 = Int64[]
    good_uv_t3_1 = Int64[]
    good_uv_t3_2 = Int64[]
    good_uv_t3_3 = Int64[]
    # good_uv_t4_1 = Int64[]
    # good_uv_t4_2 = Int64[]
    # good_uv_t4_3 = Int64[]
    # good_uv_t4_4 = Int64[]

    if (data.nvisamp>0 || data.nvisphi >0)
        vis_good = setdiff(collect(1:length(data.indx_vis)), indexes_to_discard[2])
        good_uv_vis        = data.indx_vis[vis_good];
        data.visamp        = data.visamp[vis_good];
        data.visamp_err    = data.visamp_err[vis_good];
        data.visphi        = data.visphi[vis_good];
        data.visphi_err    = data.visphi_err[vis_good];
        data.vis_baseline  = data.vis_baseline[vis_good];
        data.vis_mjd       = data.vis_mjd[vis_good];
        data.vis_lam       = data.vis_lam[vis_good];
        data.vis_dlam      = data.vis_dlam[vis_good];
        data.vis_flag      = data.vis_flag[vis_good];
        data.vis_sta_index = data.vis_sta_index[:,vis_good];
        data.nvisamp       = length(data.visamp);
        data.nvisphi       = length(data.visphi);
    end


    if (data.nv2>0)
        v2_good = setdiff(collect(1:length(data.indx_v2)), indexes_to_discard[3])
        good_uv_v2         = data.indx_v2[v2_good]
        data.v2            = data.v2[v2_good]
        data.v2_err        = data.v2_err[v2_good]
        data.v2_baseline   = data.v2_baseline[v2_good]
        data.v2_mjd        = data.v2_mjd[v2_good]
        data.v2_lam        = data.v2_lam[v2_good]
        data.v2_dlam       = data.v2_dlam[v2_good]
        data.v2_flag       = data.v2_flag[v2_good]
        data.v2_sta_index  = data.v2_sta_index[:,v2_good]
        data.nv2           = length(data.v2)
    end

    if (data.nt3amp>0 || data.nt3phi >0)
        t3_good = setdiff(collect(1:length(data.indx_t3_1)), indexes_to_discard[4])
        good_uv_t3_1 = data.indx_t3_1[t3_good]
        good_uv_t3_2 = data.indx_t3_2[t3_good]
        good_uv_t3_3 = data.indx_t3_3[t3_good]
        data.t3amp = data.t3amp[t3_good]
        data.t3amp_err = data.t3amp_err[t3_good]
        data.t3phi = data.t3phi[t3_good]
        data.t3phi_err = data.t3phi_err[t3_good]
        data.t3_baseline  = data.t3_baseline[t3_good]
        data.t3_maxbaseline  = data.t3_maxbaseline[t3_good]
        data.t3_mjd  = data.t3_mjd[t3_good]
        data.t3_lam  = data.t3_lam[t3_good]
        data.t3_dlam = data.t3_dlam[t3_good]
        data.t3_flag = data.t3_flag[t3_good]
        data.t3_sta_index = data.t3_sta_index[:,t3_good]
        data.nt3amp = length(data.t3amp)
        data.nt3phi = length(data.t3phi)
    end

    # if (data.nt4amp>0 || data.nt4phi >0)
    #     t3_good = setdiff(collect(1:length(data.indx_t4_1)), indexes_to_discard[5])
    #     good_uv_t4_1 = indx_t4_1[t4_good]
    #     good_uv_t4_2 = indx_t4_2[t4_good]
    #     good_uv_t4_3 = indx_t4_3[t4_good]
    #     good_uv_t4_4 = indx_t4_4[t4_good]
    #     data.t4amp         = data.t4amp[t4_good]
    #     data.t4amp_err     = data.t4amp_err[t4_good]
    #     data.t4phi         = data.t4phi[t4_good]
    #     data.t4phi_err     = data.t4phi_err[t4_good]
    #     data.t4_baseline   = data.t4_baseline[t4_good]
    #     data.t4_maxbaseline= data.t4_maxbaseline[t4_good]
    #     data.t4_mjd        = data.t4_mjd[t4_good]
    #     data.t4_lam        = data.t4_lam[t4_good]
    #     data.t4_dlam       = data.t4_dlam[t4_good]
    #     data.t4_flag       = data.t4_flag[t4_good]
    #     data.t4_sta_index  = data.t4_sta_index[:,t4_good]
    #     data.nt4amp        = length(data.t4amp)
    #     data.nt4phi        = length(data.t4phi)
    # end

    # uv points filtering
    uv_select  = Array{Bool}(undef, data.nuv)
    uv_select[:]  .= false;
    uv_select[good_uv_vis] .= true
    uv_select[good_uv_v2] .= true
    uv_select[good_uv_t3_1] .= true
    uv_select[good_uv_t3_2] .= true
    uv_select[good_uv_t3_3] .= true
    # uv_select[good_uv_t4_1] .= true
    # uv_select[good_uv_t4_2] .= true
    # uv_select[good_uv_t4_3] .= true
    # uv_select[good_uv_t4_4] .= true
    uv_select[indexes_to_discard[1]] .= false # disable wanted bad uv

    indx_conv = Array{Int64}(undef, length(uv_select))
    acc = 0;
    for i=1:length(uv_select)
        if uv_select[i]
            acc+=1;
        end
        indx_conv[i]=acc;
    end

    indx_uv_sel = findall(uv_select.==true)
    data.uv = data.uv[:,indx_uv_sel]
    data.uv_lam = data.uv_lam[indx_uv_sel]
    data.uv_dlam = data.uv_dlam[indx_uv_sel]
    data.uv_mjd = data.uv_mjd[indx_uv_sel]
    data.uv_baseline = data.uv_baseline[indx_uv_sel]
    data.nuv = size(data.uv,2)
    data.indx_vis  = indx_conv[good_uv_vis]
    data.indx_v2   = indx_conv[good_uv_v2]
    data.indx_t3_1 = indx_conv[good_uv_t3_1]
    data.indx_t3_2 = indx_conv[good_uv_t3_2]
    data.indx_t3_3 = indx_conv[good_uv_t3_3]
    # data.indx_t4_1 = indx_conv[good_uv_t4_1]
    # data.indx_t4_2 = indx_conv[good_uv_t4_2]
    # data.indx_t4_3 = indx_conv[good_uv_t4_3]
    # data.indx_t4_4 = indx_conv[good_uv_t4_4]
    return data
end



function remove_redundant_uv!(data::OIdata; uvtol=2e2)
# removes redundant uv points and reassigns observable indexes accordingly
     indx_red_conv, to_keep = rm_redundance_kdtree(data.uv,uvtol);
     data.uv           =  data.uv[:,to_keep]
     data.uv_lam       =  data.uv_lam[to_keep]
     data.uv_dlam      =  data.uv_dlam[to_keep]
     data.uv_mjd       =  data.uv_mjd[to_keep]
     data.uv_baseline  = data.uv_baseline[to_keep];
     data.nuv          = length(to_keep)
     if (data.nvisphi>0 || data.nvisamp > 0)
         data.indx_vis = indx_red_conv[data.indx_vis];
     end
     if data.nv2>0
         data.indx_v2 = indx_red_conv[data.indx_v2];
     end
     if (data.nt3amp>0 || data.nt3phi>0)
         data.indx_t3_1 = indx_red_conv[data.indx_t3_1];
         data.indx_t3_2 = indx_red_conv[data.indx_t3_2];
         data.indx_t3_3 = indx_red_conv[data.indx_t3_3];
     end
     # if (data.nt4amp>0 || data.nt4phi>0)
     #     data.indx_t4_1 = indx_red_conv[data.indx_t4_1];
     #     data.indx_t4_2 = indx_red_conv[data.indx_t4_2];
     #     data.indx_t4_3 = indx_red_conv[data.indx_t4_3];
     #     data.indx_t4_4 = indx_red_conv[data.indx_t4_4];
     # end
     return data
end


function readoifits(oifitsfile; targetname ="", spectralbin=[[]], temporalbin=[[]], splitting = false,  polychromatic = false, get_specbin_file=true, get_timebin_file=true,redundance_remove=true,uvtol=2e2, filter_bad_data= true, force_full_vis = false,force_full_t3 = false, filter_v2_snr_threshold=0.01, use_vis = true, use_v2 = true, use_t3 = true, use_t4 = true, use_flux = true, cutoff_minv2 = -1, cutoff_maxv2 = 2.0, cutoff_mint3amp = -1.0, cutoff_maxt3amp = 1.5, special_filter_diffvis=false, verbose=true)

    #TODO: rethink indexing by station -- there should be a station pair for each uv point, then v2/t3/etc. stations are indexed with index_v2, etc.

    #  targetname =""; spectralbin=[[]]; temporalbin=[[]]; splitting = false;  polychromatic = false; get_specbin_file=true; get_timebin_file=true;redundance_remove=false;uvtol=1.e3; filter_bad_data= true; force_full_vis = false;force_full_t3 = false; filter_v2_snr_threshold=0.01 ;use_vis = true; use_v2 = true; use_t3 = true; use_t4 = true; cutoff_minv2 = -1; cutoff_maxv2 = 2.0; cutoff_mint3amp = -1.0; cutoff_maxt3amp = 1.5; use_flux=false;
    if !isfile(oifitsfile)
        @error("readoifits could not locate the requested data file -- please check path\n")
        return [[]];
    end

    tables = OIFITS.load(oifitsfile);


    wavtables = OIFITS.select(tables,"OI_WAVELENGTH");
    if length(wavtables) ==0
        error("No OI_WAVELENGTH table in $oifitsfile")
    end
    wavtableref = [wavtables[i].insname for i=1:length(wavtables)];
    # In case there are multiple targets per file (NPOI, some old MIRC), select only the wanted data
    targetid_filter = [];
    targettables = OIFITS.select(tables, "OI_TARGET");
    if minimum(vcat([targettables[i].target_id for i=1:length(targettables)]...)) == 0
        for i=1:length(targettables)
            if targettables[i].revn!=1 # only OIFITSv2 enforces indexing
                if verbose == true
                    @warn("OI_TARGET table $i does not follow the OIFITSv2 standard - target indexing should start at 1, not 0.")
                end
            end
        end
    end

    if(targetname !="")
        if length(targettables)>1
            error("The OIFITS file has several target tables -- please specify which target you want");
        else
            targettables = targettables[1];
        end
        targetid_filter = targettables.target_id[findall(targettables.target .==targetname)]; # match target ids to target name
    else # we select everything = all target tables and all targets in these
        targetid_filter = unique(vcat([targettables[i].target_id for i=1:length(targettables)]...));
    end


    fluxtables = []
    flux_ntables = 0
    if use_flux == true
        fluxtables = OIFITS.select(tables,"OI_FLUX");
        flux_ntables = length(fluxtables);
        if flux_ntables == 0
            use_flux = false;
        end
    end

    vistables = []
    vis_ntables = 0
    if use_vis == true
        vistables = OIFITS.select(tables,"OI_VIS");
        vis_ntables = length(vistables);
        if vis_ntables == 0
            use_vis = false;
        end
    end

    v2tables = []
    v2_ntables = 0
    if use_v2 == true
        v2tables = OIFITS.select(tables,"OI_VIS2");
        v2_ntables = length(v2tables);
        if v2_ntables == 0
            use_v2 = false;
        end
    end

    t3tables = []
    t3_ntables = 0
    if use_t3 == true
        t3tables = OIFITS.select(tables,"OI_T3");
        t3_ntables = length(t3tables);
        if t3_ntables == 0
            use_t3 = false;
        end
    end

    # t4tables = []
    # t4_ntables = 0
    # if use_t4 == true
    #     t4tables = OIFITS.select(tables,"OI_T4");
    #     t4_ntables = length(t4tables);
    #     if t4_ntables == 0
    #         use_t4 = false;
    #     end
    # end


    # OI_ARRAY
    arraytables=OIFITS.select(tables,"OI_ARRAY")
    if length(arraytables) ==0
        error("No OI_ARRAY table in $oifitsfile")
    end
    arraytableref = [arraytables[i].arrname for i=1:length(arraytables)];
    array_ntables=length(arraytables)

    #get info from array tables_  #TBD -> update with knowledge from above
    vis_sta_index=Array{Array{Int64,2}}(undef, array_ntables);
    v2_sta_index=Array{Array{Int64,2}}(undef, array_ntables);
    t3_sta_index=Array{Array{Int64,3}}(undef, array_ntables);
    t4_sta_index=Array{Array{Int64,4}}(undef, array_ntables);

    station_name=Array{Array{String,1}}(undef, array_ntables)
    telescope_name=Array{Array{String,1}}(undef, array_ntables);
    station_index=Array{Array{Int64,1}}(undef, array_ntables);

    for itable = 1:array_ntables
        station_name[itable] = arraytables[itable].sta_name; # station_names
        telescope_name[itable] = arraytables[itable].tel_name; # tel_names
        station_index[itable] = arraytables[itable].sta_index; # station_indexes for matchin names to indexes in v2 and t3
    end

    station_index_offset = 0
    if minimum(vcat(station_index...)) == 0  #determine if compliant with OIFITS format (min index = 1,not 0)
        revn = maximum([arraytables[i].revn for i=1:length(arraytables)])
        if verbose==true
            if revn==2
                @warn("This file does not follow the OIFITSv2 standard - station indexing should start at 1, not 0.")
            else
                @warn("OIFITSv1 detected. To be compliant with OIFITSv2 standard, OITOOLS will internally reindex stations from 1.")
            end
        end
        station_index_offset = 1
    end

    # START OF INDEXING LOGIC (this should catch a lot of errors due to non-compliance with OIFITS v2)
    # At CHARA, station names = telescope names (for the moment)
    # At VLTI, stations and telescopes are different
    # This code can handle simple situations (pure CHARA or VLTI data) as well as CHARA + VLTI combined data
    # as long as a give station always keep the same indexes in merged data sets


    #
    # Test if there are unknown station indexes
    #

    # TODO: currently using only V2 & T3 tables, add for other data products
    unknown_station_names = String[]
    unknown_tel_names = String[]
    unknown_station_indexes = Int64[]

    for itable = 1:v2_ntables
        iarray = findall(v2tables[itable].arrname .== arraytableref)
        if length(iarray)>0 # Will fail if no corresponding OI_ARRAY table
            iarray = iarray[1]
            corresp_station_indexes = arraytables[iarray].sta_index
            for jj in unique(v2tables[itable].sta_index) # Will fail if non-existent indexes in V2 tables
                if !(jj in corresp_station_indexes)
                    @warn("V2 table $itable refers to station index $jj, non existent in OI_ARRAY=$(v2tables[itable].arrname); available indexes are $(corresp_station_indexes)")
                    push!(unknown_station_names, string("UNKN",jj))
                    push!(unknown_tel_names, string("UNKN",jj))
                    push!(unknown_station_indexes, jj);
                end
            end
        else
            @warn("V2 table $itable is missing its corresponding OI_ARRAY $(v2tables[itable].arrname)")
        end
    end

    for itable = 1:t3_ntables
        iarray = findall(t3tables[itable].arrname .== arraytableref)
        if length(iarray)>0 # Will fail if no corresponding OI_ARRAY table
            iarray = iarray[1]
            corresp_station_indexes = arraytables[iarray].sta_index
            for jj in unique(t3tables[itable].sta_index) # Will fail if non-existent indexes in T3 tables
                if !(jj in corresp_station_indexes)
                    @warn("T3 table $itable refers to station index $jj, non existent in OI_ARRAY=$(t3tables[itable].arrname); available indexes are $(corresp_station_indexes)")
                    push!(unknown_station_names, string("UNKN",jj))
                    push!(unknown_tel_names, string("UNKN",jj))
                    push!(unknown_station_indexes, jj);
                end
            end
        else
            @warn("T3 table $itable is missing its corresponding OI_ARRAY $(t3tables[itable].arrname)")
        end
    end

    if unknown_station_names != []
        # Keep only non-redundant and sort them
        non_redundant_unknown_stations = indexin(unique(unknown_station_indexes), unknown_station_indexes)
        unknown_station_indexes = unknown_station_indexes[non_redundant_unknown_stations]
        unknown_station_names = unknown_station_names[non_redundant_unknown_stations]
        unknown_tel_names = unknown_tel_names[non_redundant_unknown_stations]
        sorted_non_redundant_stations = sortperm(unknown_station_indexes)
        unknown_station_indexes = unknown_station_indexes[sorted_non_redundant_stations]
        unknown_station_names = unknown_station_names[sorted_non_redundant_stations]
        unknown_tel_names = unknown_tel_names[sorted_non_redundant_stations]
        @warn("Unknown stations, creating station names $unknown_station_names with corresponding indexes = $unknown_station_indexes \n");
    end

    # To simplify things, we merge known and unknown (missing) stations
    station_names_all = vec(vcat(station_name..., unknown_station_names))
    telescope_names_all = vec(vcat(telescope_name..., unknown_tel_names))
    station_indexes_all = vec(vcat(station_index..., unknown_station_indexes))

    # Check if index use is consistent
    list_stations = unique(station_names_all)
    nstations =  length(list_stations)
    #println("Counting ", nstations, " unique stations");

    new_station_name = Array{String}(undef, nstations)
    new_telescope_name = Array{String}(undef, nstations)
    new_station_index = zeros(Int64,nstations)

    for istation=1:length(list_stations)
        name = list_stations[istation]
        loc = findall(station_names_all .== name)
        tel = unique(telescope_names_all[loc])
        if length(tel)>1
            @warn("OIFITS file indicate that several telescopes can be located the same station -- Plots may be misleading")
        end
        tel=tel[1]
        # Renumber all stations
        new_station_name[istation]   = name
        new_telescope_name[istation] = tel
        new_station_index[istation] = istation
    end
    # we will need to convert the old indexes into the new ones
    conversion_index = spzeros(Int64, array_ntables, maximum(station_indexes_all)+station_index_offset)

    # Existing indexes in OI_ARRAY
    for itable = 1:array_ntables
        for istation = 1:length(station_name[itable])
            name = station_name[itable][istation];
            oldindx = station_index[itable][istation]
            newindx = findall(new_station_name .== name)[1]
            conversion_index[itable,station_index_offset+oldindx] = newindx;
        end
        # TODO: check logic
        if unknown_station_names != []
            nmax = sum(conversion_index[itable,:] .!=0) ; # length(station_name) ?
            conversion_index[itable, unknown_station_indexes] = nmax+1:size(conversion_index,2);
        end
    end

    # END OF STATION INDEXING LOGIC

    # Quick OI-ARRAY check
    # TODO: update so that it's not redundant with previous checks
    #       update to use visphi and flux
    all_oitables_names = unique(vcat((arraytables[i].arrname for i=1:length(arraytables))...))
    used_oiarray_tables = unique(vcat([v2tables[itable].arrname for itable = 1:v2_ntables], [t3tables[itable].arrname for itable = 1:t3_ntables]))
    if length(used_oiarray_tables)>length(all_oitables_names)
        missing_oiarray_tables =  used_oiarray_tables[.![used_oiarray_tables[i] in all_oitables_names for i=1:length(used_oiarray_tables)]]
        @warn("Missing at least $(length(used_oiarray_tables)-length(all_oitables_names)) OI-ARRAY tables in this file - won't be able to import stations properly although uv coverage will be fine.")
        @warn("Missing tables are: $(missing_oiarray_tables)");
    end

    if use_flux
        flux_ntables = length(fluxtables)
        flux = Array{Array{Float64,2}}(undef,flux_ntables)
        flux_err = Array{Array{Float64,2}}(undef,flux_ntables)
        flux_mjd = Array{Array{Float64,2}}(undef,flux_ntables)
        flux_lam = Array{Array{Float64,2}}(undef,flux_ntables)
        flux_dlam = Array{Array{Float64,2}}(undef,flux_ntables)
        flux_flag = Array{Array{Float64,2}}(undef,flux_ntables)
        flux_sta_index=Array{Array{Int64,2}}(undef,flux_ntables);
        for itable = 1:flux_ntables
            flux[itable] = fluxtables[itable].fluxdata
            flux_err[itable] = fluxtables[itable].fluxerr
            flux_flag[itable] = fluxtables[itable].flag
            flux_mjd[itable] = repeat(fluxtables[itable].mjd', size(flux[itable],1));
            iarray = findall(fluxtables[itable].arrname .== arraytableref)
            if length(iarray)>0
                flux_sta_index[itable]=conversion_index[iarray[1], station_index_offset.+repeat(fluxtables[itable].sta_index', size(flux[itable],1))];
            else
                flux_sta_index[itable]=1000 .+station_index_offset.+repeat(fluxtables[itable].sta_index', size(flux[itable],1));
            end
            iwavtable = findall(fluxtables[itable].insname .== wavtableref);
            if (length(iwavtable) != 1)
                error("Wave table confusion or Missing table\n");
            end
            flux_lam[itable] = repeat(wavtables[iwavtable[1]].eff_wave, outer=[1,size(flux[itable],2)]); # spectral channels
            flux_dlam[itable] = repeat(wavtables[iwavtable[1]].eff_band, outer=[1,size(flux[itable],2)]); # width of spectral channels
        end
    end

    if use_vis
    visamp_old = Array{Array{Float64,2}}(undef,vis_ntables);
    visamp_err_old = Array{Array{Float64,2}}(undef,vis_ntables);
    visphi_old = Array{Array{Float64,2}}(undef,vis_ntables);
    visphi_err_old = Array{Array{Float64,2}}(undef,vis_ntables);
    vis_ucoord_old = Array{Array{Float64,1}}(undef,vis_ntables);
    vis_vcoord_old = Array{Array{Float64,1}}(undef,vis_ntables);
    vis_mjd_old = Array{Array{Float64,2}}(undef,vis_ntables);
    vis_lam_old = Array{Array{Float64,2}}(undef,vis_ntables);
    vis_dlam_old = Array{Array{Float64,2}}(undef,vis_ntables);
    vis_flag_old = Array{Array{Bool,2}}(undef,vis_ntables);
    vis_u_old = Array{Array{Float64,1}}(undef,vis_ntables);
    vis_v_old = Array{Array{Float64,1}}(undef,vis_ntables);
    vis_uv_old = Array{Array{Float64,2}}(undef,vis_ntables);
    vis_baseline_old = Array{Array{Float64,1}}(undef,vis_ntables);
    vis_sta_index_old=Array{Array{Int64,2}}(undef, vis_ntables);
    for itable = 1:vis_ntables
        vis_targetid_filter = findall(sum([vistables[itable].target_id .==targetid_filter[i] for i=1:length(targetid_filter)],dims=1)[1].>0);
        visamp_old[itable] = vistables[itable].visamp[:,vis_targetid_filter];
        visamp_err_old[itable] = vistables[itable].visamperr[:,vis_targetid_filter];
        visphi_old[itable] = vistables[itable].visphi[:,vis_targetid_filter];
        visphi_err_old[itable] = vistables[itable].visphierr[:,vis_targetid_filter];
        vis_ucoord_old[itable] = vistables[itable].ucoord[vis_targetid_filter];
        vis_vcoord_old[itable] = vistables[itable].vcoord[vis_targetid_filter];
        vis_mjd_old[itable] = repeat(vistables[itable].mjd[vis_targetid_filter]', outer=[size(visamp_old[itable],1),1]); # Modified Julian Date (JD - 2400000.5)
        iarray = findall(vistables[itable].arrname .== arraytableref)
        if length(iarray)>0
            vis_sta_index_old[itable]= conversion_index[iarray[1], station_index_offset.+repeat(vistables[itable].sta_index[:,vis_targetid_filter],outer=[size(visamp_old[itable],1),1])];
        else
            vis_sta_index_old[itable]= 1000 .+station_index_offset.+repeat(vistables[itable].sta_index[:,vis_targetid_filter],outer=[size(visamp_old[itable],1),1]);
        end
        whichwav = findall(vistables[itable].insname.==wavtableref);
        vis_lam_old[itable] = repeat(wavtables[whichwav[1]].eff_wave, outer=[1,size(visamp_old[itable],2)]); # spectral channels
        vis_dlam_old[itable] = repeat(wavtables[whichwav[1]].eff_band, outer=[1,size(visamp_old[itable],2)]); # width of spectral channels
        vis_flag_old[itable] = vistables[itable].flag[:,vis_targetid_filter]; # flag for vis table
        nvis_lam_old = length(vis_lam_old[itable][:,1]);
        vis_u_old[itable] = Float64[];
        vis_v_old[itable] = Float64[];
        for u = 1:length(vis_ucoord_old[itable])
            vis_u_old[itable] = vcat(vis_u_old[itable],vis_ucoord_old[itable][u]./vis_lam_old[itable][:,1]);
            vis_v_old[itable] = vcat(vis_v_old[itable],vis_vcoord_old[itable][u]./vis_lam_old[itable][:,1]);
        end
        vis_uv_old[itable] = hcat(vec(vis_u_old[itable]),vec(vis_v_old[itable]));
        vis_baseline_old[itable] = vec(sqrt.(vis_u_old[itable].^2 + vis_v_old[itable].^2));
    end
    end

    if use_v2
    # get V2 data from tables
    v2_old = Array{Array{Float64,2}}(undef, v2_ntables);
    v2_err_old = Array{Array{Float64,2}}(undef,v2_ntables);
    v2_ucoord_old = Array{Array{Float64,1}}(undef,v2_ntables);
    v2_vcoord_old = Array{Array{Float64,1}}(undef,v2_ntables);
    v2_mjd_old = Array{Array{Float64,2}}(undef,v2_ntables);
    v2_lam_old = Array{Array{Float64,2}}(undef,v2_ntables);
    v2_dlam_old = Array{Array{Float64,2}}(undef,v2_ntables);
    v2_flag_old = Array{Array{Bool,2}}(undef,v2_ntables);
    v2_u_old = Array{Array{Float64,1}}(undef,v2_ntables);
    v2_v_old = Array{Array{Float64,1}}(undef,v2_ntables);
    v2_uv_old = Array{Array{Float64,2}}(undef,v2_ntables);
    v2_baseline_old = Array{Array{Float64,1}}(undef,v2_ntables);
    v2_sta_index_old=Array{Array{Int64,2}}(undef, v2_ntables);
    for itable = 1:v2_ntables
        v2_targetid_filter = findall(sum([v2tables[itable].target_id .==targetid_filter[i] for i=1:length(targetid_filter)],dims=1)[1].>0);
        v2_old[itable] = v2tables[itable].vis2data[:,v2_targetid_filter]; # Visibility squared
        v2_err_old[itable] = v2tables[itable].vis2err[:,v2_targetid_filter]; # error in Visibility squared
        v2_ucoord_old[itable] = v2tables[itable].ucoord[v2_targetid_filter]; # u coordinate in uv plane
        v2_vcoord_old[itable] = v2tables[itable].vcoord[v2_targetid_filter]; #  v coordinate in uv plane
        v2_mjd_old[itable] = repeat(v2tables[itable].mjd[v2_targetid_filter]', outer=[size(v2_old[itable],1),1]); # Modified Julian Date (JD - 2400000.5)
        iarray = findall(v2tables[itable].arrname .== arraytableref)
        if length(iarray)>0
            v2_sta_index_old[itable]=conversion_index[iarray[1], station_index_offset.+repeat(v2tables[itable].sta_index[:,v2_targetid_filter],outer=[size(v2_old[itable],1),1])];
        else # TO DO (or just hope people use correct OIFITS !)
            v2_sta_index_old[itable]=1000 .+station_index_offset.+repeat(v2tables[itable]sta_index[:,v2_targetid_filter],outer=[size(v2_old[itable],1),1]);
        end
        whichwav = findall(v2tables[itable].insname .== wavtableref);
        if (length(whichwav) != 1)
            error("Wave table confusion -- Missing table ?\n");
        end
        v2_lam_old[itable] = repeat(wavtables[whichwav[1]].eff_wave,  outer=[1,size(v2_old[itable],2)]); # spectral channels
        v2_dlam_old[itable] = repeat(wavtables[whichwav[1]].eff_band, outer=[1,size(v2_old[itable],2)]); # width of spectral channels
        v2_flag_old[itable] = v2tables[itable].flag[:,v2_targetid_filter]; # flag for v2 table
        nv2_lam_old = length(v2_lam_old[itable][:,1]);
        v2_u_old[itable] = Float64[];
        v2_v_old[itable] = Float64[];
        for u = 1:length(v2_ucoord_old[itable])
            v2_u_old[itable] = vcat(v2_u_old[itable],v2_ucoord_old[itable][u]./v2_lam_old[itable][:,1]);
            v2_v_old[itable] = vcat(v2_v_old[itable],v2_vcoord_old[itable][u]./v2_lam_old[itable][:,1]);
        end
        v2_uv_old[itable] = hcat(vec(v2_u_old[itable]),vec(v2_v_old[itable]));
        v2_baseline_old[itable] = vec(sqrt.(v2_u_old[itable].^2 + v2_v_old[itable].^2));
    end
    end

    if use_t3
    # same with T3, VIS
    # Get T3 data from tables
    t3amp_old = Array{Array{Float64,2}}(undef,t3_ntables);
    t3amp_err_old = Array{Array{Float64,2}}(undef,t3_ntables);
    t3phi_old = Array{Array{Float64,2}}(undef,t3_ntables);
    t3phi_err_old = Array{Array{Float64,2}}(undef,t3_ntables);
    t3_u1coord_old = Array{Array{Float64,1}}(undef,t3_ntables);
    t3_v1coord_old = Array{Array{Float64,1}}(undef,t3_ntables);
    t3_u2coord_old = Array{Array{Float64,1}}(undef,t3_ntables);
    t3_v2coord_old = Array{Array{Float64,1}}(undef,t3_ntables);
    t3_u3coord_old = Array{Array{Float64,1}}(undef,t3_ntables);
    t3_v3coord_old = Array{Array{Float64,1}}(undef,t3_ntables);
    t3_mjd_old = Array{Array{Float64,2}}(undef,t3_ntables);
    t3_lam_old = Array{Array{Float64,2}}(undef,t3_ntables);
    t3_dlam_old = Array{Array{Float64,2}}(undef,t3_ntables);
    t3_flag_old = Array{Array{Bool,2}}(undef,t3_ntables);
    t3_u1_old = Array{Array{Float64,1}}(undef,t3_ntables);
    t3_v1_old = Array{Array{Float64,1}}(undef,t3_ntables);
    t3_u2_old = Array{Array{Float64,1}}(undef,t3_ntables);
    t3_v2_old = Array{Array{Float64,1}}(undef,t3_ntables);
    t3_u3_old = Array{Array{Float64,1}}(undef,t3_ntables);
    t3_v3_old = Array{Array{Float64,1}}(undef,t3_ntables);
    t3_baseline_old = Array{Array{Float64,1}}(undef,t3_ntables);
    t3_maxbaseline_old = Array{Array{Float64,1}}(undef,t3_ntables);
    t3_sta_index_old=Array{Array{Int64,2}}(undef, t3_ntables);
    for itable = 1:t3_ntables
        t3_targetid_filter = findall(sum([t3tables[itable].target_id .==targetid_filter[i] for i=1:length(targetid_filter)],dims=1)[1].>0);
        t3amp_old[itable] = t3tables[itable].t3amp[:,t3_targetid_filter];
        t3amp_err_old[itable] = t3tables[itable].t3amperr[:,t3_targetid_filter];
        t3phi_old[itable] = t3tables[itable].t3phi[:,t3_targetid_filter];
        t3phi_err_old[itable] = t3tables[itable].t3phierr[:,t3_targetid_filter];
        t3_u1coord_old[itable] = t3tables[itable].u1coord[t3_targetid_filter];
        t3_v1coord_old[itable] = t3tables[itable].v1coord[t3_targetid_filter];
        t3_u2coord_old[itable] = t3tables[itable].u2coord[t3_targetid_filter];
        t3_v2coord_old[itable] = t3tables[itable].v2coord[t3_targetid_filter];
        t3_u3coord_old[itable] = -(t3_u1coord_old[itable] + t3_u2coord_old[itable]); # the minus takes care of complex conjugate
        t3_v3coord_old[itable] = -(t3_v1coord_old[itable] + t3_v2coord_old[itable]);
        t3_mjd_old[itable] = repeat(t3tables[itable].mjd[t3_targetid_filter]', outer=[size(t3amp_old[itable],1),1]); # Modified Julian Date (JD - 2400000.5)
        iarray = findall(t3tables[itable].arrname .== arraytableref)
        if length(iarray)>0
            t3_sta_index_old[itable]= conversion_index[iarray[1], station_index_offset.+repeat(t3tables[itable].sta_index[:,t3_targetid_filter],outer=[size(t3amp_old[itable],1),1])];
        else
            t3_sta_index_old[itable]= 1000 .+station_index_offset.+repeat(t3tables[itable].sta_index[:,t3_targetid_filter],outer=[size(t3amp_old[itable],1),1]);
        end
        whichwav = findall(t3tables[itable].insname .==wavtableref);
        t3_lam_old[itable] = repeat(wavtables[whichwav[1]].eff_wave, outer=[1,size(t3amp_old[itable],2)]); # spectral channels
        t3_dlam_old[itable] = repeat(wavtables[whichwav[1]].eff_band, outer=[1,size(t3amp_old[itable],2)]); # width of spectral channels
        t3_flag_old[itable] = t3tables[itable].flag[:,t3_targetid_filter]; # flag for t3 table
        nt3_lam_old = length(t3_lam_old[itable][:,1]);

        t3_u1_old[itable] = Float64[];
        t3_v1_old[itable] = Float64[];
        t3_u2_old[itable] = Float64[];
        t3_v2_old[itable] = Float64[];
        t3_u3_old[itable] = Float64[];
        t3_v3_old[itable] = Float64[];
        for u = 1:length(t3_u1coord_old[itable])
            t3_u1_old[itable] = vcat(t3_u1_old[itable],t3_u1coord_old[itable][u]./t3_lam_old[itable][:,1]);
            t3_v1_old[itable] = vcat(t3_v1_old[itable],t3_v1coord_old[itable][u]./t3_lam_old[itable][:,1]);
            t3_u2_old[itable] = vcat(t3_u2_old[itable],t3_u2coord_old[itable][u]./t3_lam_old[itable][:,1]);
            t3_v2_old[itable] = vcat(t3_v2_old[itable],t3_v2coord_old[itable][u]./t3_lam_old[itable][:,1]);
            t3_u3_old[itable] = vcat(t3_u3_old[itable],t3_u3coord_old[itable][u]./t3_lam_old[itable][:,1]);
            t3_v3_old[itable] = vcat(t3_v3_old[itable],t3_v3coord_old[itable][u]./t3_lam_old[itable][:,1]);
        end

        t3_baseline_old[itable] = vec((sqrt.(t3_u1_old[itable].^2 + t3_v1_old[itable].^2).*sqrt.(t3_u2_old[itable].^2 + t3_v2_old[itable].^2).*
        sqrt.(t3_u3_old[itable].^2 + t3_v3_old[itable].^2)).^(1.0/3.0));
        t3_maxbaseline_old[itable] = vec(max.(sqrt.(t3_u1_old[itable].^2 + t3_v1_old[itable].^2), sqrt.(t3_u2_old[itable].^2 + t3_v2_old[itable].^2), sqrt.(t3_u3_old[itable].^2 + t3_v3_old[itable].^2)));
    end
    end

    # if use_t4 == true
    #     # same with T4
    #     # Get T4 data from tables
    #     t4amp_old = Array{Array{Float64,2}}(undef,t4_ntables);
    #     t4amp_err_old = Array{Array{Float64,2}}(undef,t4_ntables);
    #     t4phi_old = Array{Array{Float64,2}}(undef,t4_ntables);
    #     t4phi_err_old = Array{Array{Float64,2}}(undef,t4_ntables);
    #     t4_u1coord_old = Array{Array{Float64,1}}(undef,t4_ntables);
    #     t4_v1coord_old = Array{Array{Float64,1}}(undef,t4_ntables);
    #     t4_u2coord_old = Array{Array{Float64,1}}(undef,t4_ntables);
    #     t4_v2coord_old = Array{Array{Float64,1}}(undef,t4_ntables);
    #     t4_u3coord_old = Array{Array{Float64,1}}(undef,t4_ntables);
    #     t4_v3coord_old = Array{Array{Float64,1}}(undef,t4_ntables);
    #     t4_u4coord_old = Array{Array{Float64,1}}(undef,t4_ntables);
    #     t4_v4coord_old = Array{Array{Float64,1}}(undef,t4_ntables);
    #     t4_mjd_old = Array{Array{Float64,2}}(undef,t4_ntables);
    #     t4_lam_old = Array{Array{Float64,2}}(undef,t4_ntables);
    #     t4_dlam_old = Array{Array{Float64,2}}(undef,t4_ntables);
    #     t4_flag_old = Array{Array{Bool,2}}(undef,t4_ntables);
    #     t4_u1_old = Array{Array{Float64,1}}(undef,t4_ntables);
    #     t4_v1_old = Array{Array{Float64,1}}(undef,t4_ntables);
    #     t4_u2_old = Array{Array{Float64,1}}(undef,t4_ntables);
    #     t4_v2_old = Array{Array{Float64,1}}(undef,t4_ntables);
    #     t4_u3_old = Array{Array{Float64,1}}(undef,t4_ntables);
    #     t4_v3_old = Array{Array{Float64,1}}(undef,t4_ntables);
    #     t4_u4_old = Array{Array{Float64,1}}(undef,t4_ntables);
    #     t4_v4_old = Array{Array{Float64,1}}(undef,t4_ntables);
    #     t4_baseline_old = Array{Array{Float64,1}}(undef,t4_ntables);
    #     t4_maxbaseline_old = Array{Array{Float64,1}}(undef,t4_ntables);
    #     t4_sta_index_old=Array{Array{Int64,2}}(undef, t4_ntables);
    #     for itable = 1:t4_ntables
    #         t4_targetid_filter = findall(sum([t4tables[itable].target_id .==targetid_filter[i] for i=1:length(targetid_filter)],dims=1)[1].>0);
    #         t4amp_old[itable] = t4tables[itable].t4amp[:,t4_targetid_filter];
    #         t4amp_err_old[itable] = t4tables[itable].t4amperr[:,t4_targetid_filter];
    #         t4phi_old[itable] = t4tables[itable].t4phi[:,t4_targetid_filter];
    #         t4phi_err_old[itable] = t4tables[itable].t4phierr[:,t4_targetid_filter];
    #         t4_u1coord_old[itable] = t4tables[itable].u1coord[t4_targetid_filter];
    #         t4_v1coord_old[itable] = t4tables[itable].v1coord[t4_targetid_filter];
    #         t4_u2coord_old[itable] = t4tables[itable].u2coord[t4_targetid_filter];
    #         t4_v2coord_old[itable] = t4tables[itable].v2coord[t4_targetid_filter];
    #         t4_u3coord_old[itable] = t4tables[itable].u3coord[t4_targetid_filter];
    #         t4_v3coord_old[itable] = t4tables[itable].v3coord[t4_targetid_filter];
    #         t4_u4coord_old[itable] = -(t4_u1coord_old[itable] + t4_u2coord_old[itable] + t4_u3coord_old[itable]); # the minus takes care of complex conjugate
    #         t4_v4coord_old[itable] = -(t4_v1coord_old[itable] + t4_v2coord_old[itable] + t4_v3coord_old[itable]);
    #         t4_mjd_old[itable] = repeat(t4tables[itable].mjd[t4_targetid_filter]', outer=[size(t4amp_old[itable],1),1]); # Modified Julian Date (JD - 2400000.5)
    #         iarray = findall(t4tables[itable].arrname .== arraytableref)
    #         if length(iarray)>0
    #             t4_sta_index_old[itable]=conversion_index[iarray[1], station_index_offset.+repeat(t4tables[itable].sta_index[:,t4_targetid_filter],outer=[size(t4amp_old[itable],1),1])];
    #         else
    #             t4_sta_index_old[itable]=1000 .+station_index_offset.+repeat(t4tables[itable].sta_index[:,t4_targetid_filter],outer=[size(t4amp_old[itable],1),1]);
    #         end
    #         whichwav = findall(t4tables[itable].insname .==wavtableref);
    #         t4_lam_old[itable] = repeat(wavtables[whichwav[1]].eff_wave, outer=[1,size(t4amp_old[itable],2)]); # spectral channels
    #         t4_dlam_old[itable] = repeat(wavtables[whichwav[1]].eff_band,outer=[1,size(t4amp_old[itable],2)]); # width of spectral channels
    #         t4_flag_old[itable] = t4tables[itable].flag[:,t4_targetid_filter]; # flag for t4 table
    #         nt4_lam_old = length(t4_lam_old[itable][:,1]);
    #
    #         t4_u1_old[itable] = Float64[];
    #         t4_v1_old[itable] = Float64[];
    #         t4_u2_old[itable] = Float64[];
    #         t4_v2_old[itable] = Float64[];
    #         t4_u3_old[itable] = Float64[];
    #         t4_v3_old[itable] = Float64[];
    #         t4_u4_old[itable] = Float64[];
    #         t4_v4_old[itable] = Float64[];
    #         for u = 1:length(t4_u1coord_old[itable])
    #             t4_u1_old[itable] = vcat(t4_u1_old[itable],t4_u1coord_old[itable][u]./t4_lam_old[itable][:,1]);
    #             t4_v1_old[itable] = vcat(t4_v1_old[itable],t4_v1coord_old[itable][u]./t4_lam_old[itable][:,1]);
    #             t4_u2_old[itable] = vcat(t4_u2_old[itable],t4_u2coord_old[itable][u]./t4_lam_old[itable][:,1]);
    #             t4_v2_old[itable] = vcat(t4_v2_old[itable],t4_v2coord_old[itable][u]./t4_lam_old[itable][:,1]);
    #             t4_u3_old[itable] = vcat(t4_u3_old[itable],t4_u3coord_old[itable][u]./t4_lam_old[itable][:,1]);
    #             t4_v3_old[itable] = vcat(t4_v3_old[itable],t4_v3coord_old[itable][u]./t4_lam_old[itable][:,1]);
    #             t4_u4_old[itable] = vcat(t4_u4_old[itable],t4_u4coord_old[itable][u]./t4_lam_old[itable][:,1]);
    #             t4_v4_old[itable] = vcat(t4_v4_old[itable],t4_v4coord_old[itable][u]./t4_lam_old[itable][:,1]);
    #         end
    #
    #         t4_baseline_old[itable] = vec((sqrt.(t4_u1_old[itable].^2 + t4_v1_old[itable].^2).*sqrt.(t4_u2_old[itable].^2 + t4_v2_old[itable].^2).*
    #         sqrt.(t4_u3_old[itable].^2 + t4_v3_old[itable].^2).*sqrt.(t4_u4_old[itable].^2 + t4_v4_old[itable].^2)).^(1.0/4.0));
    #         t4_maxbaseline_old[itable] = vec(max.(sqrt.(t4_u1_old[itable].^2 + t4_v1_old[itable].^2), sqrt.(t4_u2_old[itable].^2 + t4_v2_old[itable].^2), sqrt.(t4_u3_old[itable].^2 + t4_v3_old[itable].^2), sqrt.(t4_u4_old[itable].^2 + t4_v4_old[itable].^2)));
    #     end
    # end

    # combine data from all tables into a single array for each observable
    if use_flux
        flux_all = tablemerge(flux)
        flux_err_all = tablemerge(flux_err)
        flux_mjd_all = tablemerge(flux_mjd)
        flux_lam_all = tablemerge(flux_lam)
        flux_dlam_all = tablemerge(flux_dlam)
        flux_flag_all = tablemerge(flux_flag)
        flux_sta_index_all = tablemerge(flux_sta_index)
    end

    if use_vis
    visamp_all = tablemerge(visamp_old);
    visamp_err_all =  tablemerge(visamp_err_old);
    visphi_all = tablemerge(visphi_old);
    visphi_err_all =  tablemerge(visphi_err_old);
    vis_mjd_all = tablemerge(vis_mjd_old);
    vis_lam_all = tablemerge(vis_lam_old);
    vis_dlam_all = tablemerge(vis_dlam_old);
    vis_flag_all = tablemerge(vis_flag_old);
    vis_uv_all = vcat(vis_uv_old...)
    vis_baseline_all = tablemerge(vis_baseline_old);
    vis_sta_index_all= hcat([ reshape(vis_sta_index_old[i], 2, div(length(vis_sta_index_old[i]), 2)) for i=1:length(vis_sta_index_old) ]...)
    end

    if use_v2
    v2_all = tablemerge(v2_old);
    v2_err_all =  tablemerge(v2_err_old);
    v2_mjd_all = tablemerge(v2_mjd_old);
    v2_lam_all = tablemerge(v2_lam_old);
    v2_dlam_all = tablemerge(v2_dlam_old);
    v2_flag_all = tablemerge(v2_flag_old);
    v2_uv_all = vcat(v2_uv_old...)
    v2_baseline_all = tablemerge(v2_baseline_old);
    v2_sta_index_all= hcat([ reshape(v2_sta_index_old[i], 2, div(length(v2_sta_index_old[i]), 2)) for i=1:length(v2_sta_index_old) ]...)
    end

    if use_t3
    t3amp_all = tablemerge(t3amp_old);
    t3amp_err_all = tablemerge(t3amp_err_old);
    t3phi_all = tablemerge(t3phi_old);
    t3phi_err_all = tablemerge(t3phi_err_old);
    t3_mjd_all = tablemerge(t3_mjd_old);
    t3_lam_all = tablemerge(t3_lam_old);
    t3_dlam_all = tablemerge(t3_dlam_old);
    t3_flag_all = tablemerge(t3_flag_old);
    t3_u1_all = vcat(t3_u1_old...);
    t3_v1_all = vcat(t3_v1_old...);
    t3_u2_all = vcat(t3_u2_old...);
    t3_v2_all = vcat(t3_v2_old...);
    t3_u3_all = vcat(t3_u3_old...);
    t3_v3_all = vcat(t3_v3_old...);
    t3_baseline_all = tablemerge(t3_baseline_old);
    t3_maxbaseline_all = tablemerge(t3_maxbaseline_old);
    t3_sta_index_all= hcat([ reshape(t3_sta_index_old[i], 3, div(length(t3_sta_index_old[i]), 3)) for i=1:length(t3_sta_index_old) ]...)
    t3_uv_all = cat(hcat(t3_u1_all, t3_v1_all), hcat(t3_u2_all, t3_v2_all),hcat(t3_u3_all, t3_v3_all), dims=3);
    end

    # if use_t4
    #     t4amp_all = tablemerge(t4amp_old);
    #     t4amp_err_all = tablemerge(t4amp_err_old);
    #     t4phi_all = tablemerge(t4phi_old);
    #     t4phi_err_all = tablemerge(t4phi_err_old);
    #     t4_mjd_all = tablemerge(t4_mjd_old);
    #     t4_lam_all = tablemerge(t4_lam_old);
    #     t4_dlam_all = tablemerge(t4_dlam_old);
    #     t4_flag_all = tablemerge(t4_flag_old);
    #     t4_u1_all = vcat(t4_u1_old...);
    #     t4_v1_all = vcat(t4_v1_old...);
    #     t4_u2_all = vcat(t4_u2_old...);
    #     t4_v2_all = vcat(t4_v2_old...);
    #     t4_u3_all = vcat(t4_u3_old...);
    #     t4_v3_all = vcat(t4_v3_old...);
    #     t4_u4_all = vcat(t4_u4_old...);
    #     t4_v4_all = vcat(t4_v4_old...);
    #     t4_baseline_all = tablemerge(t4_baseline_old);
    #     t4_maxbaseline_all = tablemerge(t4_maxbaseline_old);
    #     t4_sta_index_all= hcat([ reshape(t4_sta_index_old[i], 4, div(length(t4_sta_index_old[i]), 4)) for i=1:length(t4_sta_index_old) ]...)
    #     t4_uv_all = cat(hcat(t4_u1_all, t4_v1_all), hcat(t4_u2_all, t4_v2_all),hcat(t4_u3_all, t4_v3_all), hcat(t4_u4_all, t4_v4_all), dims=3);
    # end

    #
    # Data splitting logic
    #

    if (polychromatic == true)||(temporalbin != [[]])||(spectralbin != [[]])
        splitting = true
    end

    # calculate default timebin if user picks timebin = [[]]
    if ((temporalbin == [[]]) && (get_timebin_file == true))
        temporalbin = [[]]
        mjd_all = [];
        if use_vis
            mjd_all = vcat(mjd_all,vis_mjd_all)
        end
        if use_v2
            mjd_all = vcat(mjd_all,v2_mjd_all)
        end
        if use_t3
            mjd_all = vcat(mjd_all, t3_mjd_all)
        end
        # if use_t4
        #     mjd_all = vcat(mjd_all, t4_mjd_all)
        # end
        if use_flux
            mjd_all = vcat(mjd_all, flux_mjd_all)
        end
        temporalbin[1] = [minimum(mjd_all)-0.001,maximum(mjd_all)+0.001]; # start & end mjd
    end

    # get spectralbin if get_spectralbin_from_file == true
    if ((spectralbin == [[]]) && (get_specbin_file == true) && (polychromatic == false))
        dlam_all = []
        lam_all = []
        if use_vis
            lam_all = vcat(lam_all, vis_lam_all)
            dlam_all = vcat(dlam_all, vis_dlam_all)
        end
        if use_v2
            lam_all = vcat(lam_all, v2_lam_all)
            dlam_all = vcat(dlam_all, v2_dlam_all)
        end
        if use_t3
            lam_all = vcat(lam_all, t3_lam_all)
            dlam_all = vcat(dlam_all, t3_dlam_all)
        end
        # if use_t4
        #     lam_all = vcat(lam_all, t4_lam_all)
        #     dlam_all = vcat(dlam_all, t4_dlam_all)
        # end
        if use_flux
            lam_all = vcat(lam_all, flux_lam_all)
            dlam_all = vcat(dlam_all, flux_dlam_all)
        end
        spectralbin[1] = vcat(spectralbin[1], minimum(lam_all)-minimum(dlam_all[argmin(lam_all)])*0.5, maximum(lam_all)+maximum(dlam_all[argmax(lam_all)])*0.5);
    end


    if ((polychromatic == true) && (get_specbin_file == true))
        if length(wavtables)>1
            @warn("There are multiple OI_WAVELENGTH tables in this file. Please specify spectralbin to select spectral channels. Wavelength selection is currently based on straight merging of tables");
        end
        #            Double check no band overlap
        wavarray = vcat([sort(hcat(wavtables[i].eff_wave-wavtables[i].eff_band/2, wavtables[i].eff_wave+wavtables[i].eff_band/2),dims=1) for i=1:length(wavtables)]...)
        overlap = sum(wavarray[2:end,2].<wavarray[1:end-1,1])
        if overlap>0
             @warn("Wavebands are overlapping - polychromatic imaging channels selection will be based on narrower bands.");
             wavarray = vcat([hcat(wavtables[i].eff_wave-wavtables[i].eff_band/2.2, wavtables[i].eff_wave+wavtables[i].eff_band/2.2) for i=1:length(wavtables)]...)
        end
        spectralbin = [wavarray[i,:] for i=1:size(wavarray,1)];
    end

    # number of spectral bins
    nwavbin = length(spectralbin);
    # number of temporal bins
    ntimebin = length(temporalbin);

    OIdataArr = Array{OIdata}(undef, nwavbin,ntimebin);

    # Define new arrays so that they get binned properly

    #TODO: replace some undef with zeros
    mean_mjd = Array{Float64}(undef, nwavbin,ntimebin);
    uv = Array{Array{Float64,2}}(undef, nwavbin,ntimebin);
    uv_lam = Array{Array{Float64,1}}(undef, nwavbin,ntimebin);
    uv_dlam = Array{Array{Float64,1}}(undef, nwavbin,ntimebin);
    uv_mjd = Array{Array{Float64,1}}(undef, nwavbin,ntimebin);
    uv_baseline = Array{Array{Float64,1}}(undef, nwavbin,ntimebin);
    #  full_sta_index = fill((vcat(Int64[]',Int64[]')),nwavbin,ntimebin);
    nuv = Array{Int64}(undef,nwavbin,ntimebin);
    nvisamp = zeros(Int64,nwavbin,ntimebin);
    nvisphi = zeros(Int64,nwavbin,ntimebin);
    vis_sta_index=fill((vcat(Int64[]',Int64[]')),nwavbin,ntimebin);
    visamp = fill((Float64[]),nwavbin,ntimebin);
    visphi = fill((Float64[]),nwavbin,ntimebin);
    visamp_err = fill((Float64[]),nwavbin,ntimebin);
    visphi_err = fill((Float64[]),nwavbin,ntimebin);
    vis_mjd = fill((Float64[]),nwavbin,ntimebin);
    vis_lam = fill((Float64[]),nwavbin,ntimebin);
    vis_dlam = fill((Float64[]),nwavbin,ntimebin);
    vis_flag = fill((Bool[]),nwavbin,ntimebin);
    vis_uv = fill((vcat(Float64[]',Float64[]')),nwavbin,ntimebin);
    vis_baseline = fill((Float64[]),nwavbin,ntimebin);
    indx_vis = fill(Int64[],nwavbin,ntimebin);

    nv2 = zeros(Int64,nwavbin,ntimebin);
    v2_sta_index=fill((vcat(Int64[]',Int64[]')),nwavbin,ntimebin);
    v2 = fill((Float64[]),nwavbin,ntimebin);
    v2_err = fill((Float64[]),nwavbin,ntimebin);
    v2_mjd = fill((Float64[]),nwavbin,ntimebin);
    v2_lam = fill((Float64[]),nwavbin,ntimebin);
    v2_dlam = fill((Float64[]),nwavbin,ntimebin);
    v2_flag = fill((Bool[]),nwavbin,ntimebin);
    v2_uv = fill((vcat(Float64[]',Float64[]')),nwavbin,ntimebin);
    v2_baseline = fill((Float64[]),nwavbin,ntimebin);
    indx_v2 = fill(Int64[],nwavbin,ntimebin);

    nt3amp = zeros(Int64,nwavbin,ntimebin);
    nt3phi = zeros(Int64,nwavbin,ntimebin);
    t3amp = fill((Float64[]),nwavbin,ntimebin);
    t3amp_err = fill((Float64[]),nwavbin,ntimebin);
    t3phi = fill((Float64[]),nwavbin,ntimebin);
    t3phi_err = fill((Float64[]),nwavbin,ntimebin);
    t3_mjd = fill((Float64[]),nwavbin,ntimebin);
    t3_lam = fill((Float64[]),nwavbin,ntimebin);
    t3_dlam = fill((Float64[]),nwavbin,ntimebin);
    t3_flag = fill((Bool[]),nwavbin,ntimebin);
    t3_uv = fill((vcat(Float64[]',Float64[]')),nwavbin,ntimebin);
    t3_baseline = fill((Float64[]),nwavbin,ntimebin);
    t3_maxbaseline = fill((Float64[]),nwavbin,ntimebin);
    indx_t3_1 = fill(Int64[],nwavbin,ntimebin);
    indx_t3_2 = fill(Int64[],nwavbin,ntimebin);
    indx_t3_3 = fill(Int64[],nwavbin,ntimebin);

    t3_sta_index=fill((vcat(Int64[]',Int64[]',Int64[]')),nwavbin,ntimebin);

    # nt4amp = zeros(Int64, nwavbin,ntimebin);
    # nt4phi = zeros(Int64, nwavbin,ntimebin);
    # t4amp = fill((Float64[]),nwavbin,ntimebin);
    # t4amp_err = fill((Float64[]),nwavbin,ntimebin);
    # t4phi = fill((Float64[]),nwavbin,ntimebin);
    # t4phi_err = fill((Float64[]),nwavbin,ntimebin);
    # t4_mjd = fill((Float64[]),nwavbin,ntimebin);
    # t4_lam = fill((Float64[]),nwavbin,ntimebin);
    # t4_dlam = fill((Float64[]),nwavbin,ntimebin);
    # t4_flag = fill((Bool[]),nwavbin,ntimebin);
    # t4_uv = fill((vcat(Float64[]',Float64[]')),nwavbin,ntimebin);
    # t4_baseline = fill((Float64[]),nwavbin,ntimebin);
    # t4_maxbaseline = fill((Float64[]),nwavbin,ntimebin);
    # indx_t4_1 = fill(Int64[],nwavbin,ntimebin);
    # indx_t4_2 = fill(Int64[],nwavbin,ntimebin);
    # indx_t4_3 = fill(Int64[],nwavbin,ntimebin);
    # indx_t4_4 = fill(Int64[],nwavbin,ntimebin);
    # t4_sta_index=fill((vcat(Int64[]',Int64[]',Int64[]',Int64[]')),nwavbin,ntimebin);

    nflux = zeros(Int64,nwavbin,ntimebin);
    flux = fill((Float64[]),nwavbin,ntimebin);
    flux_err = fill((Float64[]),nwavbin,ntimebin);
    flux_mjd = fill((Float64[]),nwavbin,ntimebin);
    flux_lam = fill((Float64[]),nwavbin,ntimebin);
    flux_dlam = fill((Float64[]),nwavbin,ntimebin);
    flux_flag =  fill((Bool[]),nwavbin,ntimebin);
    flux_sta_index= fill(Int64[],nwavbin,ntimebin);

    # New iteration for splitting data
    for itimebin = 1:ntimebin
        # combine data to one bin
        for iwavbin = 1:nwavbin

            # OPERATION 1
            # Binning according to spectral & temporal channel
            #
            bin_vis=Bool[];bin_v2=Bool[];bin_t3=Bool[]; bin_t4=Bool[]
            if splitting == true
                if use_vis
                bin_vis = (vis_mjd_all.<=temporalbin[itimebin][2]).&(vis_mjd_all.>=temporalbin[itimebin][1]).&(vis_lam_all.<=spectralbin[iwavbin][2]).&(vis_lam_all.>=spectralbin[iwavbin][1]);
                end
                if use_v2
                bin_v2 = (v2_mjd_all.<=temporalbin[itimebin][2]).&(v2_mjd_all.>=temporalbin[itimebin][1]).&(v2_lam_all.<=spectralbin[iwavbin][2]).&(v2_lam_all.>=spectralbin[iwavbin][1]);
                end
                if use_t3
                    bin_t3 = (t3_mjd_all.<=temporalbin[itimebin][2]).&(t3_mjd_all.>=temporalbin[itimebin][1]).&(t3_lam_all.<=spectralbin[iwavbin][2]).&(t3_lam_all.>=spectralbin[iwavbin][1]);
                end
                # if use_t4
                #     bin_t4 = (t4_mjd_all.<=temporalbin[itimebin][2]).&(t4_mjd_all.>=temporalbin[itimebin][1]).&(t4_lam_all.<=spectralbin[iwavbin][2]).&(t4_lam_all.>=spectralbin[iwavbin][1]);
                # end
                if use_flux
                bin_flux = (flux_mjd_all.<=temporalbin[itimebin][2]).&(flux_mjd_all.>=temporalbin[itimebin][1]).&(flux_lam_all.<=spectralbin[iwavbin][2]).&(flux_lam_all.>=spectralbin[iwavbin][1]);
                end
            else # select all
                if use_vis
                    bin_vis = Bool.(ones(length(visphi_all)))
                end
                if use_v2
                    bin_v2 = Bool.(ones(length(v2_all)))
                end
                if use_t3
                    bin_t3 = Bool.(ones(length(t3phi_all)))
                end
                # if use_t4
                #     bin_t4 = Bool.(ones(length(t4amp_all)))
                # end
                if use_flux
                    bin_flux = Bool.(ones(length(flux_all)))
                end
            end
            if use_vis
            visamp[iwavbin,itimebin] = visamp_all[bin_vis];
            visamp_err[iwavbin,itimebin] = visamp_err_all[bin_vis];
            visphi[iwavbin,itimebin] = visphi_all[bin_vis];
            visphi_err[iwavbin,itimebin] = visphi_err_all[bin_vis];
            vis_mjd[iwavbin,itimebin] = vis_mjd_all[bin_vis];
            vis_lam[iwavbin,itimebin] = vis_lam_all[bin_vis];
            vis_dlam[iwavbin,itimebin] = vis_dlam_all[bin_vis];
            vis_flag[iwavbin,itimebin] = vis_flag_all[bin_vis];
            vis_uv[iwavbin,itimebin] = hcat(vis_uv_all[bin_vis,1],vis_uv_all[bin_vis,2])';
            vis_sta_index[iwavbin,itimebin]= vis_sta_index_all[:,bin_vis];
            vis_baseline[iwavbin,itimebin] = vis_baseline_all[bin_vis];
            nvisamp[iwavbin,itimebin] = length(visamp[iwavbin,itimebin]);
            nvisphi[iwavbin,itimebin] = length(visphi[iwavbin,itimebin]);
            indx_vis[iwavbin,itimebin] = collect(1:nvisamp[iwavbin,itimebin]);
            end

            if use_flux
                flux[iwavbin,itimebin] = flux_all[bin_flux];
                flux_err[iwavbin,itimebin] = flux_err_all[bin_flux];
                flux_mjd[iwavbin,itimebin] = flux_mjd_all[bin_flux];
                flux_lam[iwavbin,itimebin] = flux_lam_all[bin_flux];
                flux_dlam[iwavbin,itimebin] = flux_dlam_all[bin_flux];
                flux_flag[iwavbin,itimebin] = flux_flag_all[bin_flux];
                flux_sta_index[iwavbin,itimebin] = flux_sta_index_all[bin_flux];
                nflux[iwavbin,itimebin] = length(flux_all[bin_flux])
            end

            if use_v2
            v2[iwavbin,itimebin] = v2_all[bin_v2];
            v2_err[iwavbin,itimebin] = v2_err_all[bin_v2];
            v2_mjd[iwavbin,itimebin] = v2_mjd_all[bin_v2];
            v2_lam[iwavbin,itimebin] = v2_lam_all[bin_v2];
            v2_dlam[iwavbin,itimebin] = v2_dlam_all[bin_v2];
            v2_flag[iwavbin,itimebin] = v2_flag_all[bin_v2];
            v2_uv[iwavbin,itimebin] = hcat(v2_uv_all[bin_v2,1],v2_uv_all[bin_v2,2])';
            v2_baseline[iwavbin,itimebin] = v2_baseline_all[bin_v2];
            v2_sta_index[iwavbin,itimebin]= v2_sta_index_all[:,bin_v2];
            nv2[iwavbin,itimebin] = length(v2[iwavbin,itimebin]);
            indx_v2[iwavbin,itimebin] = collect(nvisamp[iwavbin,itimebin].+(1:nv2[iwavbin,itimebin]));
            end

            if use_t3
            t3amp[iwavbin,itimebin] = t3amp_all[bin_t3];
            t3amp_err[iwavbin,itimebin] = t3amp_err_all[bin_t3];
            t3phi[iwavbin,itimebin] = t3phi_all[bin_t3];
            t3phi_err[iwavbin,itimebin] = t3phi_err_all[bin_t3];
            t3_mjd[iwavbin,itimebin] = t3_mjd_all[bin_t3];
            t3_lam[iwavbin,itimebin] = t3_lam_all[bin_t3];
            t3_dlam[iwavbin,itimebin] = t3_dlam_all[bin_t3];
            t3_flag[iwavbin,itimebin] = t3_flag_all[bin_t3];
            t3_uv[iwavbin,itimebin] = hcat(vec(t3_uv_all[bin_t3,1,:]),vec(t3_uv_all[bin_t3,2,:]))';
            t3_baseline[iwavbin,itimebin] = t3_baseline_all[bin_t3];
            t3_maxbaseline[iwavbin,itimebin] = t3_maxbaseline_all[bin_t3];
            t3_sta_index[iwavbin,itimebin]= t3_sta_index_all[:,bin_t3];
            nt3amp[iwavbin,itimebin] = length(t3amp[iwavbin,itimebin]);
            nt3phi[iwavbin,itimebin] = length(t3phi[iwavbin,itimebin]);
            indx_t3_1[iwavbin,itimebin] = collect(nvisamp[iwavbin,itimebin]+nv2[iwavbin,itimebin].+(1:nt3amp[iwavbin,itimebin]));
            indx_t3_2[iwavbin,itimebin] = collect(nvisamp[iwavbin,itimebin]+nv2[iwavbin,itimebin].+(nt3amp[iwavbin,itimebin]+1:2*nt3amp[iwavbin,itimebin]));
            indx_t3_3[iwavbin,itimebin] = collect(nvisamp[iwavbin,itimebin]+nv2[iwavbin,itimebin].+(2*nt3amp[iwavbin,itimebin]+1:3*nt3amp[iwavbin,itimebin]));
            end

            # if use_t4 == true
            #     t4amp[iwavbin,itimebin] = t4amp_all[bin_t4];
            #     t4amp_err[iwavbin,itimebin] = t4amp_err_all[bin_t4];
            #     t4phi[iwavbin,itimebin] = t4phi_all[bin_t4];
            #     t4phi_err[iwavbin,itimebin] = t4phi_err_all[bin_t4];
            #     t4_mjd[iwavbin,itimebin] = t4_mjd_all[bin_t4];
            #     t4_lam[iwavbin,itimebin] = t4_lam_all[bin_t4];
            #     t4_dlam[iwavbin,itimebin] = t4_dlam_all[bin_t4];
            #     t4_flag[iwavbin,itimebin] = t4_flag_all[bin_t4];
            #     t4_uv[iwavbin,itimebin] = hcat(vec(t4_uv_all[bin_t4,1,:]),vec(t4_uv_all[bin_t4,2,:]))';
            #     t4_baseline[iwavbin,itimebin] = t4_baseline_all[bin_t4];
            #     t4_maxbaseline[iwavbin,itimebin] = t4_maxbaseline_all[bin_t4];
            #     t4_sta_index[iwavbin,itimebin]=t4_sta_index_all[:,bin_t4];
            #     nt4amp[iwavbin,itimebin] = length(t4amp[iwavbin,itimebin]);
            #     nt4phi[iwavbin,itimebin] = length(t4phi[iwavbin,itimebin]);
            #     indx_t4_1[iwavbin,itimebin] = collect(nvisamp[iwavbin,itimebin]+nv2[iwavbin,itimebin]+nt3amp[iwavbin,itimebin].+(0                          +1:nt4amp[iwavbin,itimebin]));
            #     indx_t4_2[iwavbin,itimebin] = collect(nvisamp[iwavbin,itimebin]+nv2[iwavbin,itimebin]+nt3amp[iwavbin,itimebin].+(  nt4amp[iwavbin,itimebin]+1:2*nt4amp[iwavbin,itimebin]));
            #     indx_t4_3[iwavbin,itimebin] = collect(nvisamp[iwavbin,itimebin]+nv2[iwavbin,itimebin]+nt3amp[iwavbin,itimebin].+(2*nt4amp[iwavbin,itimebin]+1:3*nt4amp[iwavbin,itimebin]));
            #     indx_t4_4[iwavbin,itimebin] = collect(nvisamp[iwavbin,itimebin]+nv2[iwavbin,itimebin]+nt3amp[iwavbin,itimebin].+(3*nt4amp[iwavbin,itimebin]+1:4*nt4amp[iwavbin,itimebin]));
            # end
            #
            # if use_t4 == true
            #     uv[iwavbin,itimebin] = hcat(vis_uv[iwavbin,itimebin], v2_uv[iwavbin,itimebin],t3_uv[iwavbin,itimebin],t4_uv[iwavbin,itimebin]);
            #     uv_lam[iwavbin,itimebin]  = vcat(vis_lam[iwavbin,itimebin],v2_lam[iwavbin,itimebin], repeat(t3_lam[iwavbin,itimebin],3), repeat(t4_lam[iwavbin,itimebin],4));
            #     uv_dlam[iwavbin,itimebin] = vcat(vis_dlam[iwavbin,itimebin],v2_dlam[iwavbin,itimebin],repeat(t3_dlam[iwavbin,itimebin],3),repeat(t4_dlam[iwavbin,itimebin],4));
            #     uv_mjd[iwavbin,itimebin]  = vcat(vis_mjd[iwavbin,itimebin],v2_mjd[iwavbin,itimebin],repeat(t3_mjd[iwavbin,itimebin],3),repeat(t4_mjd[iwavbin,itimebin],4));
            # else
                uv[iwavbin,itimebin] = hcat(vis_uv[iwavbin,itimebin], v2_uv[iwavbin,itimebin],t3_uv[iwavbin,itimebin]);
                uv_lam[iwavbin,itimebin]  = vcat(vis_lam[iwavbin,itimebin], v2_lam[iwavbin,itimebin], repeat(t3_lam[iwavbin,itimebin],3));
                uv_dlam[iwavbin,itimebin] = vcat(vis_dlam[iwavbin,itimebin], v2_dlam[iwavbin,itimebin],repeat(t3_dlam[iwavbin,itimebin],3));
                uv_mjd[iwavbin,itimebin]  = vcat(vis_mjd[iwavbin,itimebin], v2_mjd[iwavbin,itimebin],repeat(t3_mjd[iwavbin,itimebin],3));
            #end
            uv_baseline[iwavbin,itimebin]  = vec(sqrt.(sum(uv[iwavbin,itimebin].^2,dims=1)));
            nuv[iwavbin,itimebin] = size(uv[iwavbin,itimebin],2);
            mean_mjd[iwavbin,itimebin] = mean(uv_mjd[iwavbin,itimebin]); # TODO obviously will fail if any NaN

            # OPERATION 2
            # Filtering bad points
            #
            if (filter_bad_data==true) # TODO: move out and make its own function

                if use_vis

                visamp_good =  (.!isnan.(visamp[iwavbin,itimebin] )) .& (.!isnan.(visamp_err[iwavbin,itimebin] )) .& (visamp_err[iwavbin,itimebin].>0.0)
                visphi_good =  (.!isnan.(visphi[iwavbin,itimebin] ) ).& (.!isnan.(visphi_err[iwavbin,itimebin] )) .& (visphi_err[iwavbin,itimebin].>0.0)

                vis_good = []
                if force_full_vis == false
                    vis_good = findall(.!vis_flag[iwavbin,itimebin] .& (visamp_good .| visphi_good) )
                else
                    vis_good = findall(.!vis_flag[iwavbin,itimebin] .& (visamp_good .& visphi_good) )
                end

                if special_filter_diffvis==true # if we use the diff vis filter, we need to do this after all the data has been read in
                                                # hence at this location we won't be actually filter anything yet
                    vis_good = findall( vis_flag[iwavbin,itimebin].!=2)
                end
                good_uv_vis = indx_vis[iwavbin,itimebin][vis_good];
                visamp[iwavbin,itimebin]        = visamp[iwavbin,itimebin][vis_good];
                visamp_err[iwavbin,itimebin]    = visamp_err[iwavbin,itimebin][vis_good];
                visphi[iwavbin,itimebin]        = visphi[iwavbin,itimebin][vis_good];
                visphi_err[iwavbin,itimebin]    = visphi_err[iwavbin,itimebin][vis_good];
                nvisamp[iwavbin,itimebin]       = length(visamp[iwavbin,itimebin]);
                nvisphi[iwavbin,itimebin]       = length(visphi[iwavbin,itimebin]);
                vis_baseline[iwavbin,itimebin]  = vis_baseline[iwavbin,itimebin][vis_good];
                vis_mjd[iwavbin,itimebin]       = vis_mjd[iwavbin,itimebin][vis_good];
                vis_lam[iwavbin,itimebin]       = vis_lam[iwavbin,itimebin][vis_good];
                vis_dlam[iwavbin,itimebin]      = vis_dlam[iwavbin,itimebin][vis_good];
                vis_flag[iwavbin,itimebin]      = vis_flag[iwavbin,itimebin][vis_good];
                vis_sta_index[iwavbin,itimebin] = vis_sta_index[iwavbin,itimebin][:,vis_good];
                end

                if use_v2
                # Filter OBVIOUSLY bad V2 data
                v2_good = findall(  (v2_flag[iwavbin,itimebin].==false) .& (v2_err[iwavbin,itimebin].>0.0)
                .& (v2_err[iwavbin,itimebin].<1.0) .& (v2[iwavbin,itimebin].>cutoff_minv2)
                .& (v2[iwavbin,itimebin].<cutoff_maxv2)
                .& .!isnan.(v2[iwavbin,itimebin]) .& .!isnan.(v2_err[iwavbin,itimebin])
                .& (abs.(v2[iwavbin,itimebin]./v2_err[iwavbin,itimebin]).>filter_v2_snr_threshold))

                good_uv_v2 = indx_v2[iwavbin,itimebin][v2_good]

                v2[iwavbin,itimebin]            = v2[iwavbin,itimebin][v2_good]
                v2_err[iwavbin,itimebin]        = v2_err[iwavbin,itimebin][v2_good]
                v2_baseline[iwavbin,itimebin]   = v2_baseline[iwavbin,itimebin][v2_good]
                nv2[iwavbin,itimebin]           = length(v2[iwavbin,itimebin])
                v2_mjd[iwavbin,itimebin]        = v2_mjd[iwavbin,itimebin][v2_good]
                v2_lam[iwavbin,itimebin]        = v2_lam[iwavbin,itimebin][v2_good]
                v2_dlam[iwavbin,itimebin]       = v2_dlam[iwavbin,itimebin][v2_good]
                v2_flag[iwavbin,itimebin]       = v2_flag[iwavbin,itimebin][v2_good]
                v2_sta_index[iwavbin,itimebin]  = v2_sta_index[iwavbin,itimebin][:,v2_good]
                end



                if use_t3
                t3amp_good =  (.!isnan.(t3amp[iwavbin,itimebin] )) .& (.!isnan.(t3amp_err[iwavbin,itimebin] )) .& (t3amp_err[iwavbin,itimebin].>0.0)
                t3phi_good =  (.!isnan.(t3phi[iwavbin,itimebin] )) .& (.!isnan.(t3phi_err[iwavbin,itimebin] )) .& (t3phi_err[iwavbin,itimebin].>0.0)

                #if force_full_t3 is set to "true", then we require both t3amp and t3phi to be defined
                t3_good = []
                if force_full_t3 == false
                    t3_good = findall(.!t3_flag[iwavbin,itimebin] .& (t3amp_good .| t3phi_good))
                else
                    t3_good = findall(.!t3_flag[iwavbin,itimebin] .& (t3amp_good .& t3phi_good) .& (t3amp[iwavbin,itimebin].>cutoff_mint3amp) .& (t3amp[iwavbin,itimebin].<cutoff_maxt3amp) )
                end
                good_uv_t3_1 = indx_t3_1[iwavbin,itimebin][t3_good]
                good_uv_t3_2 = indx_t3_2[iwavbin,itimebin][t3_good]
                good_uv_t3_3 = indx_t3_3[iwavbin,itimebin][t3_good]
                t3amp[iwavbin,itimebin] = t3amp[iwavbin,itimebin][t3_good]
                t3amp_err[iwavbin,itimebin] = t3amp_err[iwavbin,itimebin][t3_good]
                t3phi[iwavbin,itimebin] = t3phi[iwavbin,itimebin][t3_good]
                t3phi_err[iwavbin,itimebin] = t3phi_err[iwavbin,itimebin][t3_good]
                nt3amp[iwavbin,itimebin] = length(t3amp[iwavbin,itimebin])
                nt3phi[iwavbin,itimebin] = length(t3phi[iwavbin,itimebin])
                t3_baseline[iwavbin,itimebin]  = t3_baseline[iwavbin,itimebin][t3_good]
                t3_maxbaseline[iwavbin,itimebin]  = t3_maxbaseline[iwavbin,itimebin][t3_good]
                t3_mjd[iwavbin,itimebin]  = t3_mjd[iwavbin,itimebin][t3_good]
                t3_lam[iwavbin,itimebin]  = t3_lam[iwavbin,itimebin][t3_good]
                t3_dlam[iwavbin,itimebin] = t3_dlam[iwavbin,itimebin][t3_good]
                t3_flag[iwavbin,itimebin] = t3_flag[iwavbin,itimebin][t3_good]
                t3_sta_index[iwavbin,itimebin] = t3_sta_index[iwavbin,itimebin][:,t3_good]
                end

                # if use_t4 == true
                # t4amp_good =  (.!isnan.(t4amp[iwavbin,itimebin] )) .& (.!isnan.(t4amp_err[iwavbin,itimebin] )) .& (t4amp_err[iwavbin,itimebin].>0.0)
                # t4phi_good =  (.!isnan.(t4phi[iwavbin,itimebin] )) .& (.!isnan.(t4phi_err[iwavbin,itimebin] )) .& (t4phi_err[iwavbin,itimebin].>0.0)
                # # t4_good = []
                # force_full_t4 = true;
                # if force_full_t4 == false
                #     t4_good = findall(.!t4_flag[iwavbin,itimebin] .& (t4amp_good .| t4phi_good) )
                # else
                #     t4_good = findall(.!t4_flag[iwavbin,itimebin] .& (t4amp_good .& t4phi_good) )
                # end
                # good_uv_t4_1 = indx_t4_1[iwavbin,itimebin][t4_good]
                # good_uv_t4_2 = indx_t4_2[iwavbin,itimebin][t4_good]
                # good_uv_t4_3 = indx_t4_3[iwavbin,itimebin][t4_good]
                # good_uv_t4_4 = indx_t4_4[iwavbin,itimebin][t4_good]
                #
                # t4amp[iwavbin,itimebin]         = t4amp[iwavbin,itimebin][t4_good]
                # t4amp_err[iwavbin,itimebin]     = t4amp_err[iwavbin,itimebin][t4_good]
                # t4phi[iwavbin,itimebin]         = t4phi[iwavbin,itimebin][t4_good]
                # t4phi_err[iwavbin,itimebin]     = t4phi_err[iwavbin,itimebin][t4_good]
                # nt4amp[iwavbin,itimebin]        = length(t4amp[iwavbin,itimebin])
                # nt4phi[iwavbin,itimebin]        = length(t4phi[iwavbin,itimebin])
                # t4_baseline[iwavbin,itimebin]   = t4_baseline[iwavbin,itimebin][t4_good]
                # t4_maxbaseline[iwavbin,itimebin]= t4_maxbaseline[iwavbin,itimebin][t4_good]
                # t4_mjd[iwavbin,itimebin]        = t4_mjd[iwavbin,itimebin][t4_good]
                # t4_lam[iwavbin,itimebin]        = t4_lam[iwavbin,itimebin][t4_good]
                # t4_dlam[iwavbin,itimebin]       = t4_dlam[iwavbin,itimebin][t4_good]
                # t4_flag[iwavbin,itimebin]       = t4_flag[iwavbin,itimebin][t4_good]
                # t4_sta_index[iwavbin,itimebin]  = t4_sta_index[iwavbin,itimebin][:,t4_good]
                # end

                # uv points filtering
                uv_select  = Array{Bool}(undef, size(uv[iwavbin,itimebin],2))
                uv_select[:]  .= false;
                if use_vis
                uv_select[good_uv_vis] .= true
                end
                if use_v2
                uv_select[good_uv_v2] .= true
                end
                if use_t3==true
                uv_select[good_uv_t3_1] .= true
                uv_select[good_uv_t3_2] .= true
                uv_select[good_uv_t3_3] .= true
                end
                # if use_t4 == true
                # uv_select[good_uv_t4_1] .= true
                # uv_select[good_uv_t4_2] .= true
                # uv_select[good_uv_t4_3] .= true
                # uv_select[good_uv_t4_4] .= true
                # end
                #indx_conv = [sum(uv_select[1:i]) for i=1:length(uv_select)] # Performance pitfall
                indx_conv = Array{Int64}(undef, length(uv_select))
                acc = 0;
                for i=1:length(uv_select)
                    if uv_select[i]
                        acc+=1;
                    end
                    indx_conv[i]=acc;
                end

                indx_uv_sel = findall(uv_select.==true)
                uv[iwavbin,itimebin] = uv[iwavbin,itimebin][:,indx_uv_sel]
                uv_lam[iwavbin,itimebin] = uv_lam[iwavbin,itimebin][indx_uv_sel]
                uv_dlam[iwavbin,itimebin] = uv_dlam[iwavbin,itimebin][indx_uv_sel]
                uv_mjd[iwavbin,itimebin] = uv_mjd[iwavbin,itimebin][indx_uv_sel]
                uv_baseline[iwavbin,itimebin] = uv_baseline[iwavbin,itimebin][indx_uv_sel]
                nuv[iwavbin,itimebin] = size(uv[iwavbin,itimebin],2)
                if use_vis
                indx_vis[iwavbin,itimebin]  = indx_conv[good_uv_vis]
                end
                if use_v2
                indx_v2[iwavbin,itimebin]   = indx_conv[good_uv_v2]
                end
                if use_t3 == true
                indx_t3_1[iwavbin,itimebin] = indx_conv[good_uv_t3_1]
                indx_t3_2[iwavbin,itimebin] = indx_conv[good_uv_t3_2]
                indx_t3_3[iwavbin,itimebin] = indx_conv[good_uv_t3_3]
                end
                # if use_t4 == true
                # indx_t4_1[iwavbin,itimebin] = indx_conv[good_uv_t4_1]
                # indx_t4_2[iwavbin,itimebin] = indx_conv[good_uv_t4_2]
                # indx_t4_3[iwavbin,itimebin] = indx_conv[good_uv_t4_3]
                # indx_t4_4[iwavbin,itimebin] = indx_conv[good_uv_t4_4]
                # end
            end


            if (redundance_remove == true) # Remove duplicate uv points accross different data products (V2, T3, etc.)
                #TODO: find a better heuristic for uvtol, currently uvtol=1000 works for VLTI & CHARA near-infrared
                indx_red_conv, tokeep = rm_redundance_kdtree(uv[iwavbin,itimebin],uvtol);
                uv[iwavbin,itimebin] = uv[iwavbin,itimebin][:,tokeep]
                uv_lam[iwavbin,itimebin]  =  uv_lam[iwavbin,itimebin][tokeep]
                uv_dlam[iwavbin,itimebin] =  uv_dlam[iwavbin,itimebin][tokeep]
                uv_mjd[iwavbin,itimebin] =  uv_mjd[iwavbin,itimebin][tokeep]
                uv_baseline[iwavbin,itimebin]  = uv_baseline[iwavbin,itimebin][tokeep];
                nuv[iwavbin,itimebin] = size(uv[iwavbin,itimebin],2);
                if use_vis
                    indx_vis[iwavbin,itimebin] = indx_red_conv[indx_vis[iwavbin,itimebin]];
                end
                if use_v2
                indx_v2[iwavbin,itimebin] = indx_red_conv[indx_v2[iwavbin,itimebin]];
                end

                if use_t3
                    indx_t3_1[iwavbin,itimebin] = indx_red_conv[indx_t3_1[iwavbin,itimebin]];
                    indx_t3_2[iwavbin,itimebin] = indx_red_conv[indx_t3_2[iwavbin,itimebin]];
                    indx_t3_3[iwavbin,itimebin] = indx_red_conv[indx_t3_3[iwavbin,itimebin]];
                end

                # if use_t4
                #     indx_t4_1[iwavbin,itimebin] = indx_red_conv[indx_t4_1[iwavbin,itimebin]];
                #     indx_t4_2[iwavbin,itimebin] = indx_red_conv[indx_t4_2[iwavbin,itimebin]];
                #     indx_t4_3[iwavbin,itimebin] = indx_red_conv[indx_t4_3[iwavbin,itimebin]];
                #     indx_t4_4[iwavbin,itimebin] = indx_red_conv[indx_t4_4[iwavbin,itimebin]];
                # end
            end

#            OIdataArr[iwavbin,itimebin] = OIdata( visamp[iwavbin,itimebin], visamp_err[iwavbin,itimebin], visphi[iwavbin,itimebin], visphi_err[iwavbin,itimebin], vis_baseline[iwavbin,itimebin], vis_mjd[iwavbin,itimebin], vis_lam[iwavbin,itimebin], vis_dlam[iwavbin,itimebin], vis_flag[iwavbin,itimebin], v2[iwavbin,itimebin], v2_err[iwavbin,itimebin], v2_baseline[iwavbin,itimebin], v2_mjd[iwavbin,itimebin],
#            mean_mjd[iwavbin,itimebin], v2_lam[iwavbin,itimebin], v2_dlam[iwavbin,itimebin], v2_flag[iwavbin,itimebin], t3amp[iwavbin,itimebin], t3amp_err[iwavbin,itimebin], t3phi[iwavbin,itimebin], t3phi_err[iwavbin,itimebin], [], [], t3_baseline[iwavbin,itimebin],t3_maxbaseline[iwavbin,itimebin], t3_mjd[iwavbin,itimebin], t3_lam[iwavbin,itimebin], t3_dlam[iwavbin,itimebin], t3_flag[iwavbin,itimebin], t4amp[iwavbin,itimebin], t4amp_err[iwavbin,itimebin], t4phi[iwavbin,itimebin], t4phi_err[iwavbin,itimebin], t4_baseline[iwavbin,itimebin],t4_maxbaseline[iwavbin,itimebin],t4_mjd[iwavbin,itimebin], t4_lam[iwavbin,itimebin], t4_dlam[iwavbin,itimebin], t4_flag[iwavbin,itimebin],flux[iwavbin,itimebin], flux_err[iwavbin,itimebin], flux_mjd[iwavbin,itimebin], flux_lam[iwavbin,itimebin], flux_dlam[iwavbin,itimebin], flux_flag[iwavbin,itimebin], flux_sta_index[iwavbin,itimebin],
#            uv[iwavbin,itimebin], uv_lam[iwavbin,itimebin], uv_dlam[iwavbin,itimebin],uv_mjd[iwavbin,itimebin], uv_baseline[iwavbin,itimebin], nflux[iwavbin,itimebin], nvisamp[iwavbin,itimebin], nvisphi[iwavbin,itimebin], nv2[iwavbin,itimebin], nt3amp[iwavbin,itimebin], nt3phi[iwavbin,itimebin], nt4amp[iwavbin,itimebin], nt4phi[iwavbin,itimebin], nuv[iwavbin,itimebin], indx_vis[iwavbin,itimebin], indx_v2[iwavbin,itimebin],
#            indx_t3_1[iwavbin,itimebin], indx_t3_2[iwavbin,itimebin], indx_t3_3[iwavbin,itimebin],indx_t4_1[iwavbin,itimebin], indx_t4_2[iwavbin,itimebin], indx_t4_3[iwavbin,itimebin],indx_t4_4[iwavbin,itimebin],new_station_name,new_telescope_name,new_station_index,vis_sta_index[iwavbin,itimebin],v2_sta_index[iwavbin,itimebin],t3_sta_index[iwavbin,itimebin], t4_sta_index[iwavbin,itimebin],oifitsfile);
            OIdataArr[iwavbin,itimebin] = OIdata( visamp[iwavbin,itimebin], visamp_err[iwavbin,itimebin], visphi[iwavbin,itimebin], visphi_err[iwavbin,itimebin], vis_baseline[iwavbin,itimebin], vis_mjd[iwavbin,itimebin], vis_lam[iwavbin,itimebin], vis_dlam[iwavbin,itimebin], vis_flag[iwavbin,itimebin], v2[iwavbin,itimebin], v2_err[iwavbin,itimebin], v2_baseline[iwavbin,itimebin], v2_mjd[iwavbin,itimebin],
            mean_mjd[iwavbin,itimebin], v2_lam[iwavbin,itimebin], v2_dlam[iwavbin,itimebin], v2_flag[iwavbin,itimebin], t3amp[iwavbin,itimebin], t3amp_err[iwavbin,itimebin], t3phi[iwavbin,itimebin], t3phi_err[iwavbin,itimebin], [], [], t3_baseline[iwavbin,itimebin],t3_maxbaseline[iwavbin,itimebin], t3_mjd[iwavbin,itimebin], t3_lam[iwavbin,itimebin], t3_dlam[iwavbin,itimebin], t3_flag[iwavbin,itimebin],flux[iwavbin,itimebin], flux_err[iwavbin,itimebin], flux_mjd[iwavbin,itimebin], flux_lam[iwavbin,itimebin], flux_dlam[iwavbin,itimebin], flux_flag[iwavbin,itimebin], flux_sta_index[iwavbin,itimebin], uv[iwavbin,itimebin], uv_lam[iwavbin,itimebin], uv_dlam[iwavbin,itimebin],uv_mjd[iwavbin,itimebin], uv_baseline[iwavbin,itimebin], nflux[iwavbin,itimebin], nvisamp[iwavbin,itimebin], nvisphi[iwavbin,itimebin], nv2[iwavbin,itimebin], nt3amp[iwavbin,itimebin], nt3phi[iwavbin,itimebin], nuv[iwavbin,itimebin], indx_vis[iwavbin,itimebin], indx_v2[iwavbin,itimebin],
            indx_t3_1[iwavbin,itimebin], indx_t3_2[iwavbin,itimebin], indx_t3_3[iwavbin,itimebin],new_station_name,new_telescope_name,new_station_index,vis_sta_index[iwavbin,itimebin],v2_sta_index[iwavbin,itimebin],t3_sta_index[iwavbin,itimebin],oifitsfile);

        end
    end

    # Post-treatment of data
    #
    if special_filter_diffvis==true # remove diffphases when data in some spectral channels have been flagged
    for itimebin = 1:ntimebin
        indx=findall(vec((prod(hcat([(1 .-OIdataArr[i, itimebin].vis_flag) for i=1:nwavbin]...),dims=2))).==1) # diffphi baselines for which no spectral channels are flagged
        for i=1:nwavbin
            OIdataArr[i,itimebin].nvisphi = length(indx)
            OIdataArr[i,itimebin].visphi = OIdataArr[i,itimebin].visphi[indx]
            OIdataArr[i,itimebin].visphi_err = OIdataArr[i,itimebin].visphi_err[indx]
            OIdataArr[i,itimebin].vis_baseline = OIdataArr[i,itimebin].vis_baseline[indx]
            OIdataArr[i,itimebin].vis_mjd = OIdataArr[i,itimebin].vis_mjd[indx]
            OIdataArr[i,itimebin].vis_lam = OIdataArr[i,itimebin].vis_lam[indx]
            OIdataArr[i,itimebin].vis_dlam = OIdataArr[i,itimebin].vis_dlam[indx]
            OIdataArr[i,itimebin].vis_sta_index = OIdataArr[i,itimebin].vis_sta_index[:,indx]
            OIdataArr[i,itimebin].vis_flag= OIdataArr[i,itimebin].vis_flag[indx]
            OIdataArr[i,itimebin].indx_vis = OIdataArr[i,itimebin].indx_vis[indx]
        end
    end
    end

    return OIdataArr;
end

function readoifits_multiepochs(oifitsfiles; filter_bad_data=true, force_full_t3 = false) # read multiple files, each containing a single epochs
    nepochs = length(oifitsfiles);
    tepochs = Array{Float64}(undef, nepochs);
    data = Array{OIdata}(undef, nepochs);
    for i=1:nepochs
        data[i] = readoifits(oifitsfiles[i], filter_bad_data=filter_bad_data, force_full_t3 =force_full_t3 )[1,1];
        tepochs[i] = data[i].mean_mjd;
        println(oifitsfiles[i], "\t MJD: ", tepochs[i], "\t nV2 = ", data[i].nv2, "\t nT3amp = ", data[i].nt3amp, "\t nT3phi = ", data[i].nt3phi);
    end
    return nepochs, tepochs, data
end

function readoifits_multicolors(oifitsfiles; filter_bad_data=false,  force_full_t3 = false) # read multiple files, each containing a single wavelength
    nwavs = length(oifitsfiles);
    data = Array{OIdata}(undef, nwavs);
    for i=1:nwavs
        data[i] = readoifits(oifitsfiles[i], filter_bad_data=filter_bad_data, force_full_t3 =force_full_t3 )[1,1];
        println(oifitsfiles[i], "\t nV2 = ", data[i].nv2, "\t nT3amp = ", data[i].nt3amp, "\t nT3phi = ", data[i].nt3phi);
    end
    return data
end

# period in days
function time_split(mjd,period;mjd_start=mjd[1])
    timebins = (maximum(mjd) - mjd_start)/(period);
    itimebin = Int(ceil(timebins));
    temporalbin = [[],[]];
    temporalbin[1] = [mjd_start,mjd_start+period];
    temporalbin[2] = [mjd_start+period,mjd_start+2*period];
    for i = 3:itimebin
        temporalbin[2] = vcat(temporalbin[2],mjd_start+(i-1)*period);
        temporalbin[2] = vcat(temporalbin[2],mjd_start+(i)*period);
    end
    return temporalbin
end

function readfits(fitsfile; normalize = false, vectorize=false)
    x = (read((FITS(fitsfile))[1]))
    if normalize == true
        x ./= sum(x)
    end
    if vectorize == true
        x = vec(x)
    end
    return x;
end

function writefits(data, fitsfile;pixsize=-1)
    """
    pixsize should be input as mas
    """
    f = FITS(fitsfile, "w");
    if pixsize!=-1
        header = FITSHeader(["CDELT1","CDELT2","CRVAL1","CRVAL2","CRPIX1","CRPIX2"],[-(pixsize/1000.0)/(206265.0),(pixsize/1000.0)/(206265.0),0.0,0.0,(size(data)[1]/2),(size(data)[1]/2)],["Radians per Pixel","Radians per Pixel","X-coordinate of reference pixel","Y-coordinate of reference pixel","reference pixel in X","reference pixel in Y"])
        write(f, data,header=header);
    else
        write(f, data);
    end
    close(f);
end

function updatefits_aspro(fitsfile_in,fitsfile_out,res)
    """
    res should be input as mas
    """
    f = FITS(fitsfile_in);
    data = read(f[1])
    header = read_header(f[1])
    header["CDELT1"] = -res/1000.0
    header["CDELT2"] = res/1000.0
    header["CRVAL1"] = 0.0
    header["CRVAL2"] = 0.0
    header["CRPIX1"] = (size(data)[1]/2)
    header["CRPIX2"] = (size(data)[1]/2)
    header["CUNIT1"] = "arcsec"
    header["CUNIT2"] = "arcsec"
    set_comment!(header,"CDELT1","Arcseconds per Pixel")
    set_comment!(header,"CDELT2","Arcseconds per Pixel")
    set_comment!(header,"CRVAL1","X-coordinate of reference pixel")
    set_comment!(header,"CRVAL2","Y-coordinate of reference pixel")
    set_comment!(header,"CRPIX1","reference pixel in X")
    set_comment!(header,"CRPIX2","reference pixel in Y")
    fout = FITS(fitsfile_out,"w")
    write(fout,data,header=header);
    close(f)
    close(fout)
end


function oifits_prep(data::OIdata; min_v2_err_add = 0.0, min_v2_err_rel = 0.0 , v2_err_mult = 1.0, min_t3amp_err_add = 0.0,  min_t3amp_err_rel = 0.0, t3amp_err_mult = 1.0, min_t3phi_err_add = 0.0, t3phi_err_mult = 1.0, quad = false)
# e.g. MIRC from Monnier et al. https://arxiv.org/pdf/1211.6055.pdf
# min_v2_err_add = 2e-4, min_v2_err_rel = 0.066, min_t3amp_err_add = 1e-5, min_t3amp_err_rel = 0.1, min_t3phi_err_add = 1.0

# Prep V2
if quad == false
    temperr  = v2_err_mult*data.v2_err
    newerr   = abs.(data.v2*min_v2_err_rel) .+ min_v2_err_add
    newerrin = findall(newerr.>temperr)
    temperr[newerrin] = newerr[newerrin]
else
   temperr=sqrt.( (data.v2*min_v2_err_rel).^2.  + (v2_err_mult*data.v2_err).^2 .+ min_v2_err_add^2 )
end

data.v2_err = temperr

# Prep t3amp
if quad == false
    temperr  = t3amp_err_mult*data.t3amp_err
    newerr   = abs.(data.t3amp*min_t3amp_err_rel) .+ min_t3amp_err_add
    newerrin = findall(newerr.>temperr)
    temperr[newerrin] = newerr[newerrin]
else
   temperr=sqrt.( (data.t3amp*min_t3amp_err_rel).^2.  + (t3amp_err_mult*data.t3amp_err).^2 .+ min_t3amp_err_add^2 )
end

data.t3amp_err = temperr

# Prep t3phi -- need to be in degrees
if quad == true
    temperr = t3phi_err_mult*data.t3phi_err
    newerr  = min_t3phi_err_add*ones(length(data.t3phi_err))
    newerrin = findall(newerr.>temperr)
    temperr[newerrin] = newerr[newerrin]
else
    temperr = sqrt.( (t3phi_err_mult*data.t3phi_err).^2 .+ min_t3phi_err_add^2)
end
data.t3phi_err = temperr

return data
end

function oifits_prep(data::Array{OIdata,1};kwargs...) # TODO: prep diff visibilities here
    for i=1:length(data)
        oifits_prep(data[i], kwargs...)
    end
return data
end

function list_oifits_targets(oifitsfile)
    if !isfile(oifitsfile)
        @error("Could not locate the requested data file -- please check path\n")
        return [[]];
    end
    tables = OIFITS.load(oifitsfile);
    targettables = OIFITS.select(tables, "OI_TARGET");
    list_target = unique(vcat([targettables[i].target for i=1:length(targettables)]...));
    return list_target
end
