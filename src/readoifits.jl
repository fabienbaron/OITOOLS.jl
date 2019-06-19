using FITSIO

using OIFITS
#include("../../OIFITS.jl/src/OIFITS.jl"); using Main.OIFITS; #modified for T4

using Statistics

# TO DO : * use the splat operator to simplify things
#         * add differential and complex visibilities
#         * why are we negating u coordinates ?
#         * check outer=[1,8]
#         * sort indexing of uv points when t3amp and t3phi differ
#         * define booleans "use_v2, use_t3amp" and apply to all file
#         * check if t3_uv_mjd and similar are necessary

include("remove_redundance.jl")

function tablemerge(tabtomerge)
return vcat([vec(tabtomerge[i]) for i=1:length(tabtomerge)]...);
end


mutable struct OIdata
  visamp::Array{Any,1}
  visamp_err::Array{Any,1}
  visphi::Array{Any,1}
  visphi_err::Array{Any,1}
  v2::Array{Float64,1}
  v2_err::Array{Float64,1}
  v2_baseline::Array{Float64,1}
  v2_mjd::Array{Float64,1}
  mean_mjd::Float64
  v2_lam::Array{Float64,1}
  v2_dlam::Array{Float64,1}
  v2_flag::Array{Bool,1}
  t3amp::Array{Float64,1}
  t3amp_err::Array{Float64,1}
  t3phi::Array{Float64,1}
  t3phi_err::Array{Float64,1}
  t3_baseline::Array{Float64,1}
  t3_maxbaseline::Array{Float64,1}
  t3_mjd::Array{Float64,1}
  t3_lam::Array{Float64,1}
  t3_dlam::Array{Float64,1}
  t3_flag::Array{Bool,1}
  t4amp::Array{Float64,1}
  t4amp_err::Array{Float64,1}
  t4phi::Array{Float64,1}
  t4phi_err::Array{Float64,1}
  t4_baseline::Array{Float64,1}
  t4_maxbaseline::Array{Float64,1}
  t4_mjd::Array{Float64,1}
  t4_lam::Array{Float64,1}
  t4_dlam::Array{Float64,1}
  t4_flag::Array{Bool,1}
  uv::Array{Float64,2}
  uv_lam::Array{Float64,1}
  uv_dlam::Array{Float64,1}
  uv_mjd::Array{Float64,1}
  uv_baseline::Array{Float64,1}
  nvisamp::Int64
  nvisphi::Int64
  nv2::Int64
  nt3amp::Int64
  nt3phi::Int64
  nt4amp::Int64
  nt4phi::Int64
  nuv::Int64
  indx_v2::Array{Int64,1}
  indx_t3_1::Array{Int64,1}
  indx_t3_2::Array{Int64,1}
  indx_t3_3::Array{Int64,1}
  indx_t4_1::Array{Int64,1}
  indx_t4_2::Array{Int64,1}
  indx_t4_3::Array{Int64,1}
  indx_t4_4::Array{Int64,1}
  sta_name::Array{String,1}
  tel_name::Array{String,1}
  sta_index::Array{Int64,1}
  v2_sta_index::Array{Int64,2}
  t3_sta_index::Array{Int64,2}
  t4_sta_index::Array{Int64,2}
  filename::String
end

function readoifits(oifitsfile; targetname ="", spectralbin=[[]], temporalbin=[[]], splitting = false,  polychromatic = false, get_specbin_file=true, get_timebin_file=true,redundance_chk=false,uvtol=1.e3, filter_bad_data= true, force_full_t3 = false, filter_v2_snr_threshold=0.01, use_v2 = true, use_t3 = true, use_t4 = true)

  if !isfile(oifitsfile)
    @warn("Could not find file")
    return [];
  end

  tables = OIFITS.load(oifitsfile);
  wavtable = OIFITS.select(tables,"OI_WAVELENGTH");
  wavtableref = [wavtable[i][:insname] for i=1:length(wavtable)];
  # In case there are multiple targets per file (NPOI, some old MIRC), select only the wanted data
  targetid_filter = [];
  targettables = OIFITS.select(tables, "OI_TARGET");
  if(targetname !="")
    if length(targettables)>1
      error("The OIFITS file has several target tables -- import not implemented yet for this case");
    else
      targettables = targettables[1];
    end
    targetid_filter = targettables[:target_id][findall(targettables[:target].==targetname)]; # match target ids to target name
  else # we select everything = all target tables and all targets in these

    targetid_filter = unique(vcat([targettables[i][:target_id] for i=1:length(targettables)]...));
  end

  arraytables=OIFITS.select(tables,"OI_ARRAY")
  array_ntables=length(arraytables)

#get info from array tables_
  station_name=Array{Array{String,1}}(undef, array_ntables)
  telescope_name=Array{Array{String,1}}(undef, array_ntables);
  station_index=Array{Array{Int64,1}}(undef, array_ntables);
  v2_sta_index=Array{Array{Int64,2}}(undef, array_ntables);
  t3_sta_index=Array{Array{Int64,3}}(undef, array_ntables);
  t4_sta_index=Array{Array{Int64,4}}(undef, array_ntables);

  for itable = 1:array_ntables
      station_name[itable] = arraytables[itable][:sta_name]; # station_names
      telescope_name[itable] = arraytables[itable][:tel_name]; # tel_names
      station_index[itable] = arraytables[itable][:sta_index]; # station_indexes for matchin names to indexes in v2 and t3
  end

  min_sta_index = minimum(vcat(station_index...)); #determine if compliant with OIFITS format (min index = 1,not 0)
  if min_sta_index == 0
    @warn("This file does not follow the oifits standard - station indexing should start at 1, not zero")
   for itable = 1:array_ntables
     station_index[itable] .+= 1;
   end
  end


  v2table = OIFITS.select(tables,"OI_VIS2");
  v2_ntables = length(v2table);

  if v2_ntables == 0
    use_v2 = false;
  end

  t3table = OIFITS.select(tables,"OI_T3");
  t3_ntables = length(t3table);

  if t3_ntables == 0
    use_v2 = false;
  end

  t4table = OIFITS.select(tables,"OI_T4");
  t4_ntables = length(t4table);

  if t4_ntables == 0
    use_t4 = false;
  end

  # Quick OI-ARRAY check
  all_oitables_names = unique(vcat((arraytables[i][:arrname] for i=1:length(arraytables))...))
  used_oiarray_tables = unique(vcat([v2table[itable][:arrname] for itable = 1:v2_ntables], [t3table[itable][:arrname] for itable = 1:t3_ntables]))
  if length(used_oiarray_tables)>length(all_oitables_names)
    missing_oiarray_tables =  used_tables[.![used_tables[i] in oitables_names for i=1:length(used_tables)]]
    @warn("Missing at least $(length(used_oiarray_tables)-length(all_oitables_names)) OI-ARRAY tables in this file - won't be able to import stations properly although uv coverage will be fine.")
    @warn("Missing tables are: $(missing_oiarray_tables)");
  end

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
    v2_targetid_filter = findall(sum([v2table[itable][:target_id].==targetid_filter[i] for i=1:length(targetid_filter)],dims=1)[1].>0);
    v2_old[itable] = v2table[itable][:vis2data][:,v2_targetid_filter]; # Visibility squared
    v2_err_old[itable] = v2table[itable][:vis2err][:,v2_targetid_filter]; # error in Visibility squared
    v2_ucoord_old[itable] = -v2table[itable][:ucoord][v2_targetid_filter]; # u coordinate in uv plane
    v2_vcoord_old[itable] = v2table[itable][:vcoord][v2_targetid_filter]; #  v coordinate in uv plane
    v2_mjd_old[itable] = repeat(v2table[itable][:mjd][v2_targetid_filter]', outer=[size(v2_old[itable],1),1]); # Modified Julian Date (JD - 2400000.5)
    v2_sta_index_old[itable]=repeat(v2table[itable][:sta_index][:,v2_targetid_filter],outer=[size(v2_old[itable],1),1]);
    whichwav = findall(v2table[itable][:insname].==wavtableref);
    if (length(whichwav) != 1)
      error("Wave table confusion -- Missing table ?\n");
    end
    v2_lam_old[itable] = repeat(wavtable[whichwav[1]][:eff_wave],  outer=[1,size(v2_old[itable],2)]); # spectral channels
    v2_dlam_old[itable] = repeat(wavtable[whichwav[1]][:eff_band], outer=[1,size(v2_old[itable],2)]); # width of spectral channels
    v2_flag_old[itable] = v2table[itable][:flag][:,v2_targetid_filter]; # flag for v2 table
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
    t3_targetid_filter = findall(sum([t3table[itable][:target_id].==targetid_filter[i] for i=1:length(targetid_filter)],dims=1)[1].>0);
    t3amp_old[itable] = t3table[itable][:t3amp][:,t3_targetid_filter];
    t3amp_err_old[itable] = t3table[itable][:t3amperr][:,t3_targetid_filter];
    t3phi_old[itable] = t3table[itable][:t3phi][:,t3_targetid_filter];
    t3phi_err_old[itable] = t3table[itable][:t3phierr][:,t3_targetid_filter];
    t3_u1coord_old[itable] = -t3table[itable][:u1coord][t3_targetid_filter];
    t3_v1coord_old[itable] = t3table[itable][:v1coord][t3_targetid_filter];
    t3_u2coord_old[itable] = -t3table[itable][:u2coord][t3_targetid_filter];
    t3_v2coord_old[itable] = t3table[itable][:v2coord][t3_targetid_filter];
    t3_u3coord_old[itable] = -(t3_u1coord_old[itable] + t3_u2coord_old[itable]); # the minus takes care of complex conjugate
    t3_v3coord_old[itable] = -(t3_v1coord_old[itable] + t3_v2coord_old[itable]);
    t3_mjd_old[itable] = repeat(t3table[itable][:mjd][t3_targetid_filter]', outer=[size(t3amp_old[itable],1),1]); # Modified Julian Date (JD - 2400000.5)
    t3_sta_index_old[itable]= repeat(t3table[itable][:sta_index][:,t3_targetid_filter],outer=[size(t3amp_old[itable],1),1]);
    whichwav = findall(t3table[itable][:insname].==wavtableref);
    t3_lam_old[itable] = repeat(wavtable[whichwav[1]][:eff_wave], outer=[1,size(t3amp_old[itable],2)]); # spectral channels
    t3_dlam_old[itable] = repeat(wavtable[whichwav[1]][:eff_band], outer=[1,size(t3amp_old[itable],2)]); # width of spectral channels
    t3_flag_old[itable] = t3table[itable][:flag][:,t3_targetid_filter]; # flag for t3 table
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

  if use_t4 == true
  # same with T4
  # Get T4 data from tables
  t4amp_old = Array{Array{Float64,2}}(undef,t4_ntables);
  t4amp_err_old = Array{Array{Float64,2}}(undef,t4_ntables);
  t4phi_old = Array{Array{Float64,2}}(undef,t4_ntables);
  t4phi_err_old = Array{Array{Float64,2}}(undef,t4_ntables);
  t4_u1coord_old = Array{Array{Float64,1}}(undef,t4_ntables);
  t4_v1coord_old = Array{Array{Float64,1}}(undef,t4_ntables);
  t4_u2coord_old = Array{Array{Float64,1}}(undef,t4_ntables);
  t4_v2coord_old = Array{Array{Float64,1}}(undef,t4_ntables);
  t4_u3coord_old = Array{Array{Float64,1}}(undef,t4_ntables);
  t4_v3coord_old = Array{Array{Float64,1}}(undef,t4_ntables);
  t4_u4coord_old = Array{Array{Float64,1}}(undef,t4_ntables);
  t4_v4coord_old = Array{Array{Float64,1}}(undef,t4_ntables);
  t4_mjd_old = Array{Array{Float64,2}}(undef,t4_ntables);
  t4_lam_old = Array{Array{Float64,2}}(undef,t4_ntables);
  t4_dlam_old = Array{Array{Float64,2}}(undef,t4_ntables);
  t4_flag_old = Array{Array{Bool,2}}(undef,t4_ntables);
  t4_u1_old = Array{Array{Float64,1}}(undef,t4_ntables);
  t4_v1_old = Array{Array{Float64,1}}(undef,t4_ntables);
  t4_u2_old = Array{Array{Float64,1}}(undef,t4_ntables);
  t4_v2_old = Array{Array{Float64,1}}(undef,t4_ntables);
  t4_u3_old = Array{Array{Float64,1}}(undef,t4_ntables);
  t4_v3_old = Array{Array{Float64,1}}(undef,t4_ntables);
  t4_u4_old = Array{Array{Float64,1}}(undef,t4_ntables);
  t4_v4_old = Array{Array{Float64,1}}(undef,t4_ntables);
  t4_baseline_old = Array{Array{Float64,1}}(undef,t4_ntables);
  t4_maxbaseline_old = Array{Array{Float64,1}}(undef,t4_ntables);
  t4_sta_index_old=Array{Array{Int64,2}}(undef, t4_ntables);
  for itable = 1:t4_ntables
    t4_targetid_filter = findall(sum([t4table[itable][:target_id].==targetid_filter[i] for i=1:length(targetid_filter)],dims=1)[1].>0);
    t4amp_old[itable] = t4table[itable][:t4amp][:,t4_targetid_filter];
    t4amp_err_old[itable] = t4table[itable][:t4amperr][:,t4_targetid_filter];
    t4phi_old[itable] = t4table[itable][:t4phi][:,t4_targetid_filter];
    t4phi_err_old[itable] = t4table[itable][:t4phierr][:,t4_targetid_filter];
    t4_u1coord_old[itable] = -t4table[itable][:u1coord][t4_targetid_filter];
    t4_v1coord_old[itable] = t4table[itable][:v1coord][t4_targetid_filter];
    t4_u2coord_old[itable] = -t4table[itable][:u2coord][t4_targetid_filter];
    t4_v2coord_old[itable] = t4table[itable][:v2coord][t4_targetid_filter];
    t4_u3coord_old[itable] = -t4table[itable][:u3coord][t4_targetid_filter];
    t4_v3coord_old[itable] = t4table[itable][:v3coord][t4_targetid_filter];
    t4_u4coord_old[itable] = -(t4_u1coord_old[itable] + t4_u2coord_old[itable] + t4_u3coord_old[itable]); # the minus takes care of complex conjugate
    t4_v4coord_old[itable] = -(t4_v1coord_old[itable] + t4_v2coord_old[itable] + t4_v3coord_old[itable]);
    t4_mjd_old[itable] = repeat(t4table[itable][:mjd][t4_targetid_filter]', outer=[size(t4amp_old[itable],1),1]); # Modified Julian Date (JD - 2400000.5)
    t4_sta_index_old[itable]=repeat(t4table[itable][:sta_index][:,t4_targetid_filter],outer=[size(t4amp_old[itable],1),1]);
    whichwav = findall(t4table[itable][:insname].==wavtableref);
    t4_lam_old[itable] = repeat(wavtable[whichwav[1]][:eff_wave], outer=[1,size(t4amp_old[itable],2)]); # spectral channels
    t4_dlam_old[itable] = repeat(wavtable[whichwav[1]][:eff_band],outer=[1,size(t4amp_old[itable],2)]); # width of spectral channels
    t4_flag_old[itable] = t4table[itable][:flag][:,t4_targetid_filter]; # flag for t4 table
    nt4_lam_old = length(t4_lam_old[itable][:,1]);

    t4_u1_old[itable] = Float64[];
    t4_v1_old[itable] = Float64[];
    t4_u2_old[itable] = Float64[];
    t4_v2_old[itable] = Float64[];
    t4_u3_old[itable] = Float64[];
    t4_v3_old[itable] = Float64[];
    t4_u4_old[itable] = Float64[];
    t4_v4_old[itable] = Float64[];
    for u = 1:length(t4_u1coord_old[itable])
      t4_u1_old[itable] = vcat(t4_u1_old[itable],t4_u1coord_old[itable][u]./t4_lam_old[itable][:,1]);
      t4_v1_old[itable] = vcat(t4_v1_old[itable],t4_v1coord_old[itable][u]./t4_lam_old[itable][:,1]);
      t4_u2_old[itable] = vcat(t4_u2_old[itable],t4_u2coord_old[itable][u]./t4_lam_old[itable][:,1]);
      t4_v2_old[itable] = vcat(t4_v2_old[itable],t4_v2coord_old[itable][u]./t4_lam_old[itable][:,1]);
      t4_u3_old[itable] = vcat(t4_u3_old[itable],t4_u3coord_old[itable][u]./t4_lam_old[itable][:,1]);
      t4_v3_old[itable] = vcat(t4_v3_old[itable],t4_v3coord_old[itable][u]./t4_lam_old[itable][:,1]);
      t4_u4_old[itable] = vcat(t4_u4_old[itable],t4_u4coord_old[itable][u]./t4_lam_old[itable][:,1]);
      t4_v4_old[itable] = vcat(t4_v4_old[itable],t4_v4coord_old[itable][u]./t4_lam_old[itable][:,1]);
    end

    t4_baseline_old[itable] = vec((sqrt.(t4_u1_old[itable].^2 + t4_v1_old[itable].^2).*sqrt.(t4_u2_old[itable].^2 + t4_v2_old[itable].^2).*
    sqrt.(t4_u3_old[itable].^2 + t4_v3_old[itable].^2).*sqrt.(t4_u4_old[itable].^2 + t4_v4_old[itable].^2)).^(1.0/4.0));
    t4_maxbaseline_old[itable] = vec(max.(sqrt.(t4_u1_old[itable].^2 + t4_v1_old[itable].^2), sqrt.(t4_u2_old[itable].^2 + t4_v2_old[itable].^2), sqrt.(t4_u3_old[itable].^2 + t4_v3_old[itable].^2), sqrt.(t4_u4_old[itable].^2 + t4_v4_old[itable].^2)));
  end
  end

  # combine data from various tables into one
  # This linearizes the tables

  v2_all = tablemerge(v2_old);
  v2_err_all =  tablemerge(v2_err_old);
  v2_mjd_all = tablemerge(v2_mjd_old);
  v2_lam_all = tablemerge(v2_lam_old);
  v2_dlam_all = tablemerge(v2_dlam_old);
  v2_flag_all = tablemerge(v2_flag_old);
  v2_uv_all = vcat(v2_uv_old...)
  v2_baseline_all = tablemerge(v2_baseline_old);
  v2_sta_index_all= hcat([ reshape(v2_sta_index_old[i], 2, div(length(v2_sta_index_old[i]), 2)) for i=1:length(v2_sta_index_old) ]...)

  # Fix to v2_sta_index_all since it's a frequent mistake to have 0 there
  # OIFITS standard says 1 should be the minimum
  if min_sta_index == 0
   v2_sta_index_all .+= 1;
  end

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


  # to do: handle unusual cases such as no v2
  if min_sta_index == 0
   t3_sta_index_all .+= 1;
  end
  t3_uv_all = cat(hcat(t3_u1_all, t3_v1_all), hcat(t3_u2_all, t3_v2_all),hcat(t3_u3_all, t3_v3_all), dims=3);

if use_t4 == true
  t4amp_all = tablemerge(t4amp_old);
  t4amp_err_all = tablemerge(t4amp_err_old);
  t4phi_all = tablemerge(t4phi_old);
  t4phi_err_all = tablemerge(t4phi_err_old);
  t4_mjd_all = tablemerge(t4_mjd_old);
  t4_lam_all = tablemerge(t4_lam_old);
  t4_dlam_all = tablemerge(t4_dlam_old);
  t4_flag_all = tablemerge(t4_flag_old);
  t4_u1_all = vcat(t4_u1_old...);
  t4_v1_all = vcat(t4_v1_old...);
  t4_u2_all = vcat(t4_u2_old...);
  t4_v2_all = vcat(t4_v2_old...);
  t4_u3_all = vcat(t4_u3_old...);
  t4_v3_all = vcat(t4_v3_old...);
  t4_u4_all = vcat(t4_u4_old...);
  t4_v4_all = vcat(t4_v4_old...);
  t4_baseline_all = tablemerge(t4_baseline_old);
  t4_maxbaseline_all = tablemerge(t4_maxbaseline_old);
  t4_sta_index_all= hcat([ reshape(t4_sta_index_old[i], 4, div(length(t4_sta_index_old[i]), 4)) for i=1:length(t4_sta_index_old) ]...)
  if min_sta_index == 0
     t4_sta_index_all .+= 1;
  end
  t4_uv_all = cat(hcat(t4_u1_all, t4_v1_all), hcat(t4_u2_all, t4_v2_all),hcat(t4_u3_all, t4_v3_all), hcat(t4_u4_all, t4_v4_all), dims=3);
end

#
# Data splitting logic
#

if (temporalbin != [[]])||(spectralbin != [[]]|| (polychromatic == true))
  splitting = true
end

  # calculate default timebin if user picks timebin = [[]]
  if ((temporalbin == [[]]) && (get_timebin_file == true))
    temporalbin = [[]]
    temporalbin[1] = [minimum(v2_mjd_all)-0.001,maximum(v2_mjd_all)+0.001]; # start & end mjd
  end

  # get spectralbin if get_spectralbin_from_file == true

  if ((spectralbin == [[]]) && (get_specbin_file == true) && (polychromatic == false))
    spectralbin[1] = vcat(spectralbin[1], minimum(v2_lam_all)-minimum(v2_dlam_all[argmin(v2_lam_all)])*0.5, maximum(v2_lam_all)+maximum(v2_dlam_all[argmax(v2_lam_all)])*0.5);
  end

  if ((polychromatic == true) && (get_specbin_file == true))
    if length(wavtable)>1
      @warn("Sorry, multiple OI_WAVELENGTH table, I do not know how to pick spectral channels. I'll make an attempt though...");
      #nwav = unique(vcat([length(wavtable[i][:eff_wave]) for i=1:length(wavtable)]...)) # are we always using the same mode ?
      #allwav = sort(vcat([wavtable[i][:eff_wave] for i=1:length(wavtable)]...))
    else
      wavarray = hcat(wavtable[1][:eff_wave]-wavtable[1][:eff_band]/2, wavtable[1][:eff_wave]+wavtable[1][:eff_band]/2);
      spectralbin = [wavarray[i,:] for i=1:size(wavarray,1)];
    end
  end

  # number of spectral bins
  nwavbin = length(spectralbin);
  # number of temporal bins
  ntimebin = length(temporalbin);

  OIdataArr = Array{OIdata}(undef, nwavbin,ntimebin);

  # Define new arrays so that they get binned properly

  #TODO: replace some undef with zeros
  mean_mjd = Array{Float64}(undef, nwavbin,ntimebin);
  full_uv = Array{Array{Float64,2}}(undef, nwavbin,ntimebin);
  full_uv_lam = Array{Array{Float64,1}}(undef, nwavbin,ntimebin);
  full_uv_dlam = Array{Array{Float64,1}}(undef, nwavbin,ntimebin);
  full_uv_mjd = Array{Array{Float64,1}}(undef, nwavbin,ntimebin);
  full_uv_baseline = Array{Array{Float64,1}}(undef, nwavbin,ntimebin);
  full_sta_index = fill((vcat(Int64[]',Int64[]')),nwavbin,ntimebin);
  nuv = Array{Int64}(undef,nwavbin,ntimebin);

  nv2 = Array{Int64}(undef,nwavbin,ntimebin);
  v2_sta_index_new=fill((vcat(Int64[]',Int64[]')),nwavbin,ntimebin);
  v2_new = fill((Float64[]),nwavbin,ntimebin);
  v2_err_new = fill((Float64[]),nwavbin,ntimebin);
  v2_mjd_new = fill((Float64[]),nwavbin,ntimebin);
  v2_lam_new = fill((Float64[]),nwavbin,ntimebin);
  v2_dlam_new = fill((Float64[]),nwavbin,ntimebin);
  v2_flag_new = fill((Bool[]),nwavbin,ntimebin);
  v2_uv_new = fill((vcat(Float64[]',Float64[]')),nwavbin,ntimebin);
  v2_baseline_new = fill((Float64[]),nwavbin,ntimebin);
  indx_v2 = Array{Array{Int64,1}}(undef,nwavbin,ntimebin);

  nt3amp = Array{Int64}(undef,nwavbin,ntimebin);
  nt3phi = Array{Int64}(undef,nwavbin,ntimebin);
  t3amp_new = fill((Float64[]),nwavbin,ntimebin);
  t3amp_err_new = fill((Float64[]),nwavbin,ntimebin);
  t3phi_new = fill((Float64[]),nwavbin,ntimebin);
  t3phi_err_new = fill((Float64[]),nwavbin,ntimebin);
  t3_mjd_new = fill((Float64[]),nwavbin,ntimebin);
  t3_lam_new = fill((Float64[]),nwavbin,ntimebin);
  t3_dlam_new = fill((Float64[]),nwavbin,ntimebin);
  t3_flag_new = fill((Bool[]),nwavbin,ntimebin);
  t3_uv_new = fill((vcat(Float64[]',Float64[]')),nwavbin,ntimebin);
  t3_baseline_new = fill((Float64[]),nwavbin,ntimebin);
  t3_maxbaseline_new = fill((Float64[]),nwavbin,ntimebin);
  indx_t3_1 = Array{Array{Int64,1}}(undef,nwavbin,ntimebin);
  indx_t3_2 = Array{Array{Int64,1}}(undef,nwavbin,ntimebin);
  indx_t3_3 = Array{Array{Int64,1}}(undef,nwavbin,ntimebin);
  t3_sta_index_new=fill((vcat(Int64[]',Int64[]',Int64[]')),nwavbin,ntimebin);

  nt4amp = zeros(Int64, nwavbin,ntimebin);
  nt4phi = zeros(Int64, nwavbin,ntimebin);
  t4amp_new = fill((Float64[]),nwavbin,ntimebin);
  t4amp_err_new = fill((Float64[]),nwavbin,ntimebin);
  t4phi_new = fill((Float64[]),nwavbin,ntimebin);
  t4phi_err_new = fill((Float64[]),nwavbin,ntimebin);
  t4_mjd_new = fill((Float64[]),nwavbin,ntimebin);
  t4_lam_new = fill((Float64[]),nwavbin,ntimebin);
  t4_dlam_new = fill((Float64[]),nwavbin,ntimebin);
  t4_flag_new = fill((Bool[]),nwavbin,ntimebin);
  t4_uv_new = fill((vcat(Float64[]',Float64[]')),nwavbin,ntimebin);
  t4_baseline_new = fill((Float64[]),nwavbin,ntimebin);
  t4_maxbaseline_new = fill((Float64[]),nwavbin,ntimebin);
  indx_t4_1 = fill(Int64[],nwavbin,ntimebin);
  indx_t4_2 = fill(Int64[],nwavbin,ntimebin);
  indx_t4_3 = fill(Int64[],nwavbin,ntimebin);
  indx_t4_4 = fill(Int64[], nwavbin,ntimebin);
  t4_sta_index_new=fill((vcat(Int64[]',Int64[]',Int64[]',Int64[]')),nwavbin,ntimebin);

  # New iteration for splitting data
  for itimebin = 1:ntimebin
    # combine data to one bin
    for iwavbin = 1:nwavbin
      # Binning filter
      bin_v2=Bool[];bin_t3=Bool[]; bin_t4=Bool[]
      if splitting == true
        bin_v2 = (v2_mjd_all.<=temporalbin[itimebin][2]).&(v2_mjd_all.>=temporalbin[itimebin][1]).&(v2_lam_all.<=spectralbin[iwavbin][2]).&(v2_lam_all.>=spectralbin[iwavbin][1]);
        bin_t3 = (t3_mjd_all.<=temporalbin[itimebin][2]).&(t3_mjd_all.>=temporalbin[itimebin][1]).&(t3_lam_all.<=spectralbin[iwavbin][2]).&(t3_lam_all.>=spectralbin[iwavbin][1]);
        if use_t4 == true
          bin_t4 = (t4_mjd_all.<=temporalbin[itimebin][2]).&(t4_mjd_all.>=temporalbin[itimebin][1]).&(t4_lam_all.<=spectralbin[iwavbin][2]).&(t4_lam_all.>=spectralbin[iwavbin][1]);
        end
      else
        bin_v2 = Bool.(ones(length(v2_all)))
        bin_t3 = Bool.(ones(length(t3amp_all)))
        if use_t4 == true
          bin_t4 = Bool.(ones(length(t4amp_all)))
        end
      end

      v2_new[iwavbin,itimebin] = v2_all[bin_v2];
      v2_err_new[iwavbin,itimebin] = v2_err_all[bin_v2];
      v2_mjd_new[iwavbin,itimebin] = v2_mjd_all[bin_v2];
      v2_lam_new[iwavbin,itimebin] = v2_lam_all[bin_v2];
      v2_dlam_new[iwavbin,itimebin] = v2_dlam_all[bin_v2];
      v2_flag_new[iwavbin,itimebin] = v2_flag_all[bin_v2];
      v2_uv_new[iwavbin,itimebin] = hcat(v2_uv_all[bin_v2,1],v2_uv_all[bin_v2,2])';
      v2_baseline_new[iwavbin,itimebin] = v2_baseline_all[bin_v2];
      v2_sta_index_new[iwavbin,itimebin]= v2_sta_index_all[:,bin_v2];
      nv2[iwavbin,itimebin] = length(v2_new[iwavbin,itimebin]);
      indx_v2[iwavbin,itimebin] = collect(1:nv2[iwavbin,itimebin]);


      t3amp_new[iwavbin,itimebin] = t3amp_all[bin_t3];
      t3amp_err_new[iwavbin,itimebin] = t3amp_err_all[bin_t3];
      t3phi_new[iwavbin,itimebin] = t3phi_all[bin_t3];
      t3phi_err_new[iwavbin,itimebin] = t3phi_err_all[bin_t3];
      t3_mjd_new[iwavbin,itimebin] = t3_mjd_all[bin_t3];
      t3_lam_new[iwavbin,itimebin] = t3_lam_all[bin_t3];
      t3_dlam_new[iwavbin,itimebin] = t3_dlam_all[bin_t3];
      t3_flag_new[iwavbin,itimebin] = t3_flag_all[bin_t3];
      t3_uv_new[iwavbin,itimebin] = hcat(vec(t3_uv_all[bin_t3,1,:]),vec(t3_uv_all[bin_t3,2,:]))';
      t3_baseline_new[iwavbin,itimebin] = t3_baseline_all[bin_t3];
      t3_maxbaseline_new[iwavbin,itimebin] = t3_maxbaseline_all[bin_t3];
      t3_sta_index_new[iwavbin,itimebin]= t3_sta_index_all[:,bin_t3];
      nt3amp[iwavbin,itimebin] = length(t3amp_new[iwavbin,itimebin]);
      nt3phi[iwavbin,itimebin] = length(t3phi_new[iwavbin,itimebin]);
      indx_t3_1[iwavbin,itimebin] = collect(nv2[iwavbin,itimebin].+(1:nt3amp[iwavbin,itimebin]));
      indx_t3_2[iwavbin,itimebin] = collect(nv2[iwavbin,itimebin].+(nt3amp[iwavbin,itimebin]+1:2*nt3amp[iwavbin,itimebin]));
      indx_t3_3[iwavbin,itimebin] = collect(nv2[iwavbin,itimebin].+(2*nt3amp[iwavbin,itimebin]+1:3*nt3amp[iwavbin,itimebin]));


      if use_t4 == true
      t4amp_new[iwavbin,itimebin] = t4amp_all[bin_t4];
      t4amp_err_new[iwavbin,itimebin] = t4amp_err_all[bin_t4];
      t4phi_new[iwavbin,itimebin] = t4phi_all[bin_t4];
      t4phi_err_new[iwavbin,itimebin] = t4phi_err_all[bin_t4];
      t4_mjd_new[iwavbin,itimebin] = t4_mjd_all[bin_t4];
      t4_lam_new[iwavbin,itimebin] = t4_lam_all[bin_t4];
      t4_dlam_new[iwavbin,itimebin] = t4_dlam_all[bin_t4];
      t4_flag_new[iwavbin,itimebin] = t4_flag_all[bin_t4];
      t4_uv_new[iwavbin,itimebin] = hcat(vec(t4_uv_all[bin_t4,1,:]),vec(t4_uv_all[bin_t4,2,:]))';
      t4_baseline_new[iwavbin,itimebin] = t4_baseline_all[bin_t4];
      t4_maxbaseline_new[iwavbin,itimebin] = t4_maxbaseline_all[bin_t4];
      t4_sta_index_new[iwavbin,itimebin]=t4_sta_index_all[:,bin_t4];
      nt4amp[iwavbin,itimebin] = length(t4amp_new[iwavbin,itimebin]);
      nt4phi[iwavbin,itimebin] = length(t4phi_new[iwavbin,itimebin]);
      indx_t4_1[iwavbin,itimebin] = collect(nv2[iwavbin,itimebin]+nt3amp[iwavbin,itimebin].+(0                          +1:nt4amp[iwavbin,itimebin]));
      indx_t4_2[iwavbin,itimebin] = collect(nv2[iwavbin,itimebin]+nt3amp[iwavbin,itimebin].+(  nt4amp[iwavbin,itimebin]+1:2*nt4amp[iwavbin,itimebin]));
      indx_t4_3[iwavbin,itimebin] = collect(nv2[iwavbin,itimebin]+nt3amp[iwavbin,itimebin].+(2*nt4amp[iwavbin,itimebin]+1:3*nt4amp[iwavbin,itimebin]));
      indx_t4_4[iwavbin,itimebin] = collect(nv2[iwavbin,itimebin]+nt3amp[iwavbin,itimebin].+(3*nt4amp[iwavbin,itimebin]+1:4*nt4amp[iwavbin,itimebin]));
      end

      if use_t4 == true
        full_uv[iwavbin,itimebin] = hcat(v2_uv_new[iwavbin,itimebin],t3_uv_new[iwavbin,itimebin],t4_uv_new[iwavbin,itimebin]);
        full_uv_lam[iwavbin,itimebin]  = vcat(v2_lam_new[iwavbin,itimebin], repeat(t3_lam_new[iwavbin,itimebin],3), repeat(t4_lam_new[iwavbin,itimebin],4));
        full_uv_dlam[iwavbin,itimebin] = vcat(v2_dlam_new[iwavbin,itimebin],repeat(t3_dlam_new[iwavbin,itimebin],3),repeat(t4_dlam_new[iwavbin,itimebin],4));
        full_uv_mjd[iwavbin,itimebin]  = vcat(v2_mjd_new[iwavbin,itimebin],repeat(t3_mjd_new[iwavbin,itimebin],3),repeat(t4_mjd_new[iwavbin,itimebin],4));
      else
        full_uv[iwavbin,itimebin] = hcat(v2_uv_new[iwavbin,itimebin],t3_uv_new[iwavbin,itimebin]);
        full_uv_lam[iwavbin,itimebin]  = vcat(v2_lam_new[iwavbin,itimebin], repeat(t3_lam_new[iwavbin,itimebin],3));
        full_uv_dlam[iwavbin,itimebin] = vcat(v2_dlam_new[iwavbin,itimebin],repeat(t3_dlam_new[iwavbin,itimebin],3));
        full_uv_mjd[iwavbin,itimebin]  = vcat(v2_mjd_new[iwavbin,itimebin],repeat(t3_mjd_new[iwavbin,itimebin],3));
      end
      full_uv_baseline[iwavbin,itimebin]  = vec(sqrt.(sum(full_uv[iwavbin,itimebin].^2,dims=1)));
      nuv[iwavbin,itimebin] = size(full_uv[iwavbin,itimebin],2);
      mean_mjd[iwavbin,itimebin] = mean(full_uv_mjd[iwavbin,itimebin]); # TODO obviously will fail if any NaN

      if (redundance_chk == true) # temp fix?
        full_uv[iwavbin,itimebin], indx_redun = rm_redundance_kdtree(full_uv[iwavbin,itimebin],uvtol);
        #TODO: fix full_uv_lam, dlam, and mjd, full_uv_baseline
        nuv[iwavbin,itimebin] = size(full_uv[iwavbin,itimebin],2);
        indx_v2[iwavbin,itimebin] = indx_redun[indx_v2[iwavbin,itimebin]];
        indx_t3_1[iwavbin,itimebin] = indx_redun[indx_t3_1[iwavbin,itimebin]];
        indx_t3_2[iwavbin,itimebin] = indx_redun[indx_t3_2[iwavbin,itimebin]];
        indx_t3_3[iwavbin,itimebin] = indx_redun[indx_t3_3[iwavbin,itimebin]];

        if use_t4 == true
          indx_t4_1[iwavbin,itimebin] = indx_redun[indx_t4_1[iwavbin,itimebin]];
          indx_t4_2[iwavbin,itimebin] = indx_redun[indx_t4_2[iwavbin,itimebin]];
          indx_t4_3[iwavbin,itimebin] = indx_redun[indx_t4_3[iwavbin,itimebin]];
          indx_t4_4[iwavbin,itimebin] = indx_redun[indx_t4_4[iwavbin,itimebin]];
        end
      end

      OIdataArr[iwavbin,itimebin] = OIdata([],[], [], [], v2_new[iwavbin,itimebin], v2_err_new[iwavbin,itimebin], v2_baseline_new[iwavbin,itimebin], v2_mjd_new[iwavbin,itimebin],
      mean_mjd[iwavbin,itimebin], v2_lam_new[iwavbin,itimebin], v2_dlam_new[iwavbin,itimebin], v2_flag_new[iwavbin,itimebin], t3amp_new[iwavbin,itimebin],
      t3amp_err_new[iwavbin,itimebin], t3phi_new[iwavbin,itimebin], t3phi_err_new[iwavbin,itimebin], t3_baseline_new[iwavbin,itimebin],t3_maxbaseline_new[iwavbin,itimebin],
      t3_mjd_new[iwavbin,itimebin], t3_lam_new[iwavbin,itimebin], t3_dlam_new[iwavbin,itimebin], t3_flag_new[iwavbin,itimebin], t4amp_new[iwavbin,itimebin], t4amp_err_new[iwavbin,itimebin], t4phi_new[iwavbin,itimebin], t4phi_err_new[iwavbin,itimebin], t4_baseline_new[iwavbin,itimebin],t4_maxbaseline_new[iwavbin,itimebin],
      t4_mjd_new[iwavbin,itimebin], t4_lam_new[iwavbin,itimebin], t4_dlam_new[iwavbin,itimebin], t4_flag_new[iwavbin,itimebin],
      full_uv[iwavbin,itimebin], full_uv_lam[iwavbin,itimebin], full_uv_dlam[iwavbin,itimebin],full_uv_mjd[iwavbin,itimebin], full_uv_baseline[iwavbin,itimebin], 0, 0, nv2[iwavbin,itimebin], nt3amp[iwavbin,itimebin], nt3phi[iwavbin,itimebin], nt4amp[iwavbin,itimebin], nt4phi[iwavbin,itimebin], nuv[iwavbin,itimebin], indx_v2[iwavbin,itimebin],
      indx_t3_1[iwavbin,itimebin], indx_t3_2[iwavbin,itimebin], indx_t3_3[iwavbin,itimebin],indx_t4_1[iwavbin,itimebin], indx_t4_2[iwavbin,itimebin], indx_t4_3[iwavbin,itimebin],indx_t4_4[iwavbin,itimebin],station_name[1],telescope_name[1],station_index[:][1],v2_sta_index_new[iwavbin,itimebin],t3_sta_index_new[iwavbin,itimebin], t4_sta_index_new[iwavbin,itimebin],oifitsfile);

      if (filter_bad_data==true) # TODO: move out and make its own function
        # Filter OBVIOUSLY bad V2 data
        v2_good = findall(  (OIdataArr[iwavbin,itimebin].v2_flag.==false) .& (OIdataArr[iwavbin,itimebin].v2_err.>0.0)
        .& (OIdataArr[iwavbin,itimebin].v2_err.<1.0) .& (OIdataArr[iwavbin,itimebin].v2.>-0.2)
        .& (OIdataArr[iwavbin,itimebin].v2.<1.2)
        .& .!isnan.(OIdataArr[iwavbin,itimebin].v2) .& .!isnan.(OIdataArr[iwavbin,itimebin].v2_err)
        .& (abs.(OIdataArr[iwavbin,itimebin].v2./OIdataArr[iwavbin,itimebin].v2_err).>filter_v2_snr_threshold))

        good_uv_v2 = OIdataArr[iwavbin,itimebin].indx_v2[v2_good]

        OIdataArr[iwavbin,itimebin].v2 = OIdataArr[iwavbin,itimebin].v2[v2_good]
        OIdataArr[iwavbin,itimebin].v2_err = OIdataArr[iwavbin,itimebin].v2_err[v2_good]
        OIdataArr[iwavbin,itimebin].v2_baseline = OIdataArr[iwavbin,itimebin].v2_baseline[v2_good]
        OIdataArr[iwavbin,itimebin].nv2 = length(OIdataArr[iwavbin,itimebin].v2)
        OIdataArr[iwavbin,itimebin].v2_mjd  = OIdataArr[iwavbin,itimebin].v2_mjd[v2_good]
        OIdataArr[iwavbin,itimebin].v2_lam  = OIdataArr[iwavbin,itimebin].v2_lam[v2_good]
        OIdataArr[iwavbin,itimebin].v2_dlam = OIdataArr[iwavbin,itimebin].v2_dlam[v2_good]
        OIdataArr[iwavbin,itimebin].v2_flag = OIdataArr[iwavbin,itimebin].v2_flag[v2_good]
        OIdataArr[iwavbin,itimebin].v2_sta_index = OIdataArr[iwavbin,itimebin].v2_sta_index[:,v2_good]

        t3amp_good =  (.!isnan.(OIdataArr[iwavbin,itimebin].t3amp )) .& (.!isnan.(OIdataArr[iwavbin,itimebin].t3amp_err )) .& (OIdataArr[iwavbin,itimebin].t3amp_err.>0.0)
        t3phi_good =  (.!isnan.(OIdataArr[iwavbin,itimebin].t3phi ) ).& (.!isnan.(OIdataArr[iwavbin,itimebin].t3phi_err )) .& (OIdataArr[iwavbin,itimebin].t3phi_err.>0.0)

        t4amp_good =  (.!isnan.(OIdataArr[iwavbin,itimebin].t4amp )) .& (.!isnan.(OIdataArr[iwavbin,itimebin].t4amp_err )) .& (OIdataArr[iwavbin,itimebin].t4amp_err.>0.0)
        t4phi_good =  (.!isnan.(OIdataArr[iwavbin,itimebin].t4phi ) ).& (.!isnan.(OIdataArr[iwavbin,itimebin].t4phi_err )) .& (OIdataArr[iwavbin,itimebin].t4phi_err.>0.0)


        #if force_full_t3 is set to "true", then we require both t3amp and t3phi to be defined
        t3_good = []
        if force_full_t3 == false
          t3_good = findall(.!OIdataArr[iwavbin,itimebin].t3_flag .& (t3amp_good .| t3phi_good) )
        else
          t3_good = findall(.!OIdataArr[iwavbin,itimebin].t3_flag .& (t3amp_good .& t3phi_good) )
        end
        good_uv_t3_1 = OIdataArr[iwavbin,itimebin].indx_t3_1[t3_good]
        good_uv_t3_2 = OIdataArr[iwavbin,itimebin].indx_t3_2[t3_good]
        good_uv_t3_3 = OIdataArr[iwavbin,itimebin].indx_t3_3[t3_good]
        OIdataArr[iwavbin,itimebin].t3amp = OIdataArr[iwavbin,itimebin].t3amp[t3_good]
        OIdataArr[iwavbin,itimebin].t3amp_err = OIdataArr[iwavbin,itimebin].t3amp_err[t3_good]
        OIdataArr[iwavbin,itimebin].t3phi = OIdataArr[iwavbin,itimebin].t3phi[t3_good]
        OIdataArr[iwavbin,itimebin].t3phi_err = OIdataArr[iwavbin,itimebin].t3phi_err[t3_good]
        OIdataArr[iwavbin,itimebin].nt3amp = length(OIdataArr[iwavbin,itimebin].t3amp)
        OIdataArr[iwavbin,itimebin].nt3phi = length(OIdataArr[iwavbin,itimebin].t3phi)
        OIdataArr[iwavbin,itimebin].t3_baseline  = OIdataArr[iwavbin,itimebin].t3_baseline[t3_good]
        OIdataArr[iwavbin,itimebin].t3_maxbaseline  = OIdataArr[iwavbin,itimebin].t3_maxbaseline[t3_good]
        OIdataArr[iwavbin,itimebin].t3_mjd  = OIdataArr[iwavbin,itimebin].t3_mjd[t3_good]
        OIdataArr[iwavbin,itimebin].t3_lam  = OIdataArr[iwavbin,itimebin].t3_lam[t3_good]
        OIdataArr[iwavbin,itimebin].t3_dlam = OIdataArr[iwavbin,itimebin].t3_dlam[t3_good]
        OIdataArr[iwavbin,itimebin].t3_flag = OIdataArr[iwavbin,itimebin].t3_flag[t3_good]
        OIdataArr[iwavbin,itimebin].t3_sta_index = OIdataArr[iwavbin,itimebin].t3_sta_index[:,t3_good]

        # t4_good = []
        force_full_t4 = true;
        if force_full_t4 == false
          t4_good = findall(.!OIdataArr[iwavbin,itimebin].t4_flag .& (t4amp_good .| t4phi_good) )
        else
          t4_good = findall(.!OIdataArr[iwavbin,itimebin].t4_flag .& (t4amp_good .& t4phi_good) )
        end
        good_uv_t4_1 = OIdataArr[iwavbin,itimebin].indx_t4_1[t4_good]
        good_uv_t4_2 = OIdataArr[iwavbin,itimebin].indx_t4_2[t4_good]
        good_uv_t4_3 = OIdataArr[iwavbin,itimebin].indx_t4_3[t4_good]
        good_uv_t4_4 = OIdataArr[iwavbin,itimebin].indx_t4_4[t4_good]

        OIdataArr[iwavbin,itimebin].t4amp = OIdataArr[iwavbin,itimebin].t4amp[t4_good]
        OIdataArr[iwavbin,itimebin].t4amp_err = OIdataArr[iwavbin,itimebin].t4amp_err[t4_good]
        OIdataArr[iwavbin,itimebin].t4phi = OIdataArr[iwavbin,itimebin].t4phi[t4_good]
        OIdataArr[iwavbin,itimebin].t4phi_err = OIdataArr[iwavbin,itimebin].t4phi_err[t4_good]
        OIdataArr[iwavbin,itimebin].nt4amp = length(OIdataArr[iwavbin,itimebin].t4amp)
        OIdataArr[iwavbin,itimebin].nt4phi = length(OIdataArr[iwavbin,itimebin].t4phi)
        OIdataArr[iwavbin,itimebin].t4_baseline  = OIdataArr[iwavbin,itimebin].t4_baseline[t4_good]
        OIdataArr[iwavbin,itimebin].t4_maxbaseline  = OIdataArr[iwavbin,itimebin].t4_maxbaseline[t4_good]
        OIdataArr[iwavbin,itimebin].t4_mjd  = OIdataArr[iwavbin,itimebin].t4_mjd[t4_good]
        OIdataArr[iwavbin,itimebin].t4_lam  = OIdataArr[iwavbin,itimebin].t4_lam[t4_good]
        OIdataArr[iwavbin,itimebin].t4_dlam = OIdataArr[iwavbin,itimebin].t4_dlam[t4_good]
        OIdataArr[iwavbin,itimebin].t4_flag = OIdataArr[iwavbin,itimebin].t4_flag[t4_good]
        OIdataArr[iwavbin,itimebin].t4_sta_index = OIdataArr[iwavbin,itimebin].t4_sta_index[:,t4_good]


        # uv points filtering
        uv_select  = Array{Bool}(undef, size(OIdataArr[iwavbin,itimebin].uv,2))
        uv_select[:]  .= false;
        uv_select[good_uv_v2] .= true
        uv_select[good_uv_t3_1] .= true
        uv_select[good_uv_t3_2] .= true
        uv_select[good_uv_t3_3] .= true

        uv_select[good_uv_t4_1] .= true
        uv_select[good_uv_t4_2] .= true
        uv_select[good_uv_t4_3] .= true
        uv_select[good_uv_t4_4] .= true

        indx_conv = [sum(uv_select[1:i]) for i=1:length(uv_select)]
        OIdataArr[iwavbin,itimebin].uv = OIdataArr[iwavbin,itimebin].uv[:,findall(uv_select.==true)]
        OIdataArr[iwavbin,itimebin].uv_lam = OIdataArr[iwavbin,itimebin].uv_lam[findall(uv_select.==true)]
        OIdataArr[iwavbin,itimebin].uv_dlam = OIdataArr[iwavbin,itimebin].uv_dlam[findall(uv_select.==true)]
        OIdataArr[iwavbin,itimebin].uv_mjd = OIdataArr[iwavbin,itimebin].uv_mjd[findall(uv_select.==true)]
        OIdataArr[iwavbin,itimebin].uv_baseline = OIdataArr[iwavbin,itimebin].uv_baseline[findall(uv_select.==true)]
        OIdataArr[iwavbin,itimebin].nuv = size(OIdataArr[iwavbin,itimebin].uv,2)
        OIdataArr[iwavbin,itimebin].indx_v2 =   indx_conv[good_uv_v2]
        OIdataArr[iwavbin,itimebin].indx_t3_1 = indx_conv[good_uv_t3_1]
        OIdataArr[iwavbin,itimebin].indx_t3_2 = indx_conv[good_uv_t3_2]
        OIdataArr[iwavbin,itimebin].indx_t3_3 = indx_conv[good_uv_t3_3]

        OIdataArr[iwavbin,itimebin].indx_t4_1 = indx_conv[good_uv_t4_1]
        OIdataArr[iwavbin,itimebin].indx_t4_2 = indx_conv[good_uv_t4_2]
        OIdataArr[iwavbin,itimebin].indx_t4_3 = indx_conv[good_uv_t4_3]
        OIdataArr[iwavbin,itimebin].indx_t4_4 = indx_conv[good_uv_t4_4]
      end

    end
  end

  return OIdataArr;
end

function readoifits_multiepochs(oifitsfiles; filter_bad_data=false,  force_full_t3 = false) # read multiple files, each containing a single epochs
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

function writefits(data, fitsfile)
f = FITS(fitsfile, "w");
write(f, data);
close(f);
end
