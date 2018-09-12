include("../OITOOLS.jl/oitools.jl")
using NFFT
#Functions used in main Functions


function cvis_to_t3_conj(cvis, indx1, indx2, indx3)
    #get t3 from caculated visibilities
    #because of ordering of v_ij j>i we need to use conjugate (we use conj(t3_13), rather than t3_31)
 # t3 = cvis[indx1].*cvis[indx2].*conj(cvis[indx3]);
  t3 = cvis[indx1].*cvis[indx2].*conj(cvis[indx3]);
  t3amp = abs.(t3);
  t3phi = angle.(t3)*180./pi;
  return t3, t3amp, t3phi
end

function v2mapt3(i,j)
    # really dumb way of doing this
    # there is probably an analytic formula for this
    listi = find(v2_stations[1,:].==i);
    return listi[findfirst(v2_stations[2,listi].==j)]
end

function get_v2_baselines(N,station_xyz)
    # determine V2 Baselines and make necessary arrays
    nv2 = Int64(N*(N-1)/2);
    v2_baselines = Array{Float64}(3,nv2);
    v2_stations  = Array{Int64}(2,nv2);
    v2_stations_nonredun=Array{Int64}(2,nv2);
    v2_indx      = Array{Int64}(nv2);
    baseline_list = Array{String}(nv2);
    indx = 1
    for i=1:N+function
      for j=i+1:N
        #if(j!=i)
            v2_baselines[:,indx] = station_xyz[j]-station_xyz[i];
            v2_stations[:,indx] = [i,j];
            v2_indx[indx] = indx;
            indx += 1
    #    end
      end
    end
    return nv2,v2_baselines,v2_stations,v2_stations_nonredun,v2_indx,baseline_list,indx

end

function get_t3_baselines(N,station_xyz)
    # determine T3 Baselines and make neccessary https://mail.google.com/mail/u/0/#inboxarrays
    nt3 = binomial(N,3);# V2 Baselines
    t3_baselines=Array{Float64}(3,3,nt3);
    t3_stations  = Array{Int64}(3, nt3);
    t3_indx_1  = Array{Int64}(nt3);
    t3_indx_2  = Array{Int64}(nt3);
    t3_indx_3  = Array{Int64}(nt3);
    indx = 1;
    for i=1:N
      for j=i+1:N-1
        for k=j+1:N
            t3_baselines[:, 1, indx] = station_xyz[j]-station_xyz[i];
            t3_baselines[:, 2, indx] = station_xyz[k]-station_xyz[j];
            t3_baselines[:, 3, indx] = station_xyz[i]-station_xyz[k];
            t3_stations[:,indx]=[i,j,k];
            t3_indx_1[indx] = v2mapt3(i,j);
            t3_indx_2[indx] = v2mapt3(j,k);
            t3_indx_3[indx] = v2mapt3(i,k);
            indx += 1
        end
      end
    end
    return nt3,t3_baselines,t3_stations,t3_indx_1,t3_indx_2,t3_indx_3,indx
end

function get_uv(l,h,λ,δ,v2_baselines)
    # Merge all uv
    #Use following expression only if there are missing baselines somewhere
    #vis_baselines = hcat(v2_baselines, t3_baselines[:, 1, :], t3_baselines[:, 2, :], t3_baselines[:, 3, :])
    #Expression to use for pure simulation where the full complement of v2 and t3 are created
    vis_baselines = deepcopy(v2_baselines)
    nuv = size(vis_baselines, 2);
    # Baselines to proj Baselines (uv m)
    # Now compute the UV coordinates, again according to the APIS++ standards.
    u_M = -sin(l)*sin.(h) .* vis_baselines[1,:] .+ cos.(h) .* vis_baselines[2,:]+cos(l)*sin.(h).*vis_baselines[3,:];
    v_M = (sin(l)*sin(δ)*cos.(h).+cos(l)*cos(δ)) .* vis_baselines[1,:] + sin(δ)*sin.(h) .* vis_baselines[2,:]+(-cos(l)*cos.(h)*sin(δ).+sin(l)*cos(δ)) .* vis_baselines[3,:];
    w_M =  (-sin(l)*cos(δ)* cos.(h).+cos(l)*sin(δ)).* vis_baselines[1,:] - cos(δ)* sin.(h) .* vis_baselines[2,:]  .+ (cos(l)*cos.(h)*cos(δ).+sin(l)sin(δ)) .* vis_baselines[3,:];
    u_M = -u_M
    v_M = -v_M
    # proj baselines to (uv wav)
    u = reshape((1./λ)'.*vec(u_M), (nuv,nhours,length(λ)));
    v = reshape((1./λ)'.*vec(v_M), (nuv,nhours,length(λ)));
    w = reshape((1./λ)'.*vec(w_M), (nuv,nhours,length(λ)));
    uv = vcat(vec(u)',vec(v)')
    return nuv,uv,u_M,v_M,w_M
end

function get_uv_indxes(nhours,nuv,nv2,nt3,v2_indx,t3_indx_1,t3_indx_2,t3_indx_3)
    v2_indx_M = repmat(v2_indx,1,nhours)+repmat(Int64.(collect(linspace(0,nuv*(nhours-1),nhours)))',nv2)
    t3_indx_1_M = repmat(t3_indx_1,1,nhours)+repmat(Int64.(collect(linspace(0,nuv*(nhours-1),nhours)))',nt3)
    t3_indx_2_M = repmat(t3_indx_2,1,nhours)+repmat(Int64.(collect(linspace(0,nuv*(nhours-1),nhours)))',nt3)
    t3_indx_3_M = repmat(t3_indx_3,1,nhours)+repmat(Int64.(collect(linspace(0,nuv*(nhours-1),nhours)))',nt3)
    #indexes separated by wavelength channel
    v2_indx_w = vec(repmat(vec(v2_indx_M),1,nw)+repmat(Int64.(collect(linspace(0,nuv*nhours*(nw-1),nw)))',nv2*nhours))
    t3_indx_1_w = vec(repmat(vec(t3_indx_1_M),1,nw)+repmat(Int64.(collect(linspace(0,nuv*nhours*(nw-1),nw)))',nt3*nhours))
    t3_indx_2_w = vec(repmat(vec(t3_indx_2_M),1,nw)+repmat(Int64.(collect(linspace(0,nuv*nhours*(nw-1),nw)))',nt3*nhours))
    t3_indx_3_w = vec(repmat(vec(t3_indx_3_M),1,nw)+repmat(Int64.(collect(linspace(0,nuv*nhours*(nw-1),nw)))',nt3*nhours))
    t3_baseline = (sqrt.(uv[1,t3_indx_1_w].^2 + uv[2,t3_indx_1_w].^2).*
        sqrt.(uv[1,t3_indx_2_w].^2 + uv[2,t3_indx_2_w].^2).*
        sqrt.(uv[1,t3_indx_3_w].^2 + uv[2,t3_indx_3_w].^2)).^(1./3.);
    return v2_indx_M,t3_indx_1_M,t3_indx_2_M,t3_indx_3_M,v2_indx_w,t3_indx_1_w,t3_indx_2_w,t3_indx_3_w
 end

 function get_simulated_image(image_file,pixsize,uv,v2_indx_w,t3_indx_1_w, t3_indx_2_w, t3_indx_3_w)
     #simulate observation using input image file and pixsize
     x = (read((FITS(image_file))[1])); x=x[:,end:-1:1]; nx = (size(x))[1]; x=vec(x)/sum(x);
     #nfft_plan = setup_nfft(-uv, nx, pixsize)
     #cvis_model = image_to_cvis_nfft(x, nfft_plan)
     dft = setup_dft(-uv, nx, pixsize);
     cvis_model = image_to_cvis_dft(x, dft);
     v2_model = cvis_to_v2(cvis_model, v2_indx_w);
     t3_model, t3amp_model, t3phi_model = cvis_to_t3_conj(cvis_model, t3_indx_1_w, t3_indx_2_w, t3_indx_3_w);
     #Add noise
     v2_model_err = 1.5*v2_model+1e-5
     v2_model += v2_model_err.*randn(length(v2_model))
     t3amp_model_err =1.5*t3amp_model+1e-6
     t3amp_model += t3amp_model_err.*randn(length(t3amp_model))
     t3phi_model_err = zeros(length(t3phi_model))+2. # degree  -- there is another way of setting this with Haniff formula
    t3phi_model += t3phi_model_err.*randn(length(t3phi_model))
    return v2_model,v2_model_err,t3amp_model,t3amp_model_err,t3phi_model,t3phi_model_err
 end

function prep_arrays(obsv_info)
    oi_array=[tel_names,sta_names,sta_index,tel_diams,staxyz];
    target_array=[target_id,target,raep0,decep0,equinox,ra_err,dec_err,sysvel,veltyp,veldef,pmra,pmdec,pmra_err,pmdec_err,parallax,para_err,spectyp];
    wave_array=[eff_wave,eff_band];
    vis2_array=[target_id_vis2,time_vis2,mjd_vis2,int_time_vis2,v2_model,v2_model_err,vec(ucoord_vis2),vec(vcoord_vis2),v2_model_stations,flag_vis2];
    t3_array=[target_id_t3,time_t3,mjd_t3,int_time_t3,t3amp_model,t3amp_model_err,t3phi_model,t3phi_model_err,vec(u1coord),vec(v1coord),vec(u2coord),vec(v2coord),t3_model_stations,flag_t3];
    return oi_array,target_array,wave_array,vis2_array,t3oi_array
 end

function read_array_file(array_file)
    info=readdlm(array_file)
    array_name=info[find(info.=="name"),3]
    lat=info[find(info.=="lat"),3]
    lon=info[find(info.=="lon"),3]
    alt=info[find(info.=="alt"),3]
    coord=info[find(info.=="coord"),3]
    throughput=info[find(info.=="throughput"),3]
    wind_speed=info[find(info.=="wind_speed"),3]
    r0=info[find(info.=="r0"),3]
    ntel=info[find(info.=="ntel")+function,3]
    tel_diams=info[find(info.=="tel_diams"),3:ntel[1]+2]
    tel_gain=info[find(info.=="tel_gain"),3:ntel[1]+2]
    tel_names=info[find(info.=="tel_names"),3:ntel[1]+2]
    sta_names=info[find(info.=="sta_names"),3:ntel[1]+2]
    sta_index=info[find(info.=="sta_index"),3:ntel[1]+2]
    sta_xyz=info[find(info.=="sta_xyz"),3:ntel[1]*3+2]
    station_xyz=Array{Float64}(ntel[1],3)
    for i=1:ntel[1]
        station_xyz[i,1:3]=sta_xyz[(i*3-2):i*3]
    end
    return array_name,lat,lon,alt,coord,throughput,wind_speed,r0,ntel,tel_diams,tel_gain,tel_names,sta_names,sta_index,station_xyz
 end

function read_obs_file(obsv_file)
    info=readdlm(obsv_file)
    target_id=info[find(info.=="target_id"),3]
    target=info[find(info.=="target"),3]
    raep0=info[find(info.=="raep0"),3]
    decep0=info[find(info.=="decep0"),3]
    equinox=info[find(info.=="equinox"),3]
    ra_err=info[find(info.=="ra_err"),3]
    dec_err = info[find(info.=="dec_err"),3]
    sysvel = info[find(info.=="sysvel"),3]
    veltyp =  info[find(info.=="veltyp"),3]
    veldef = info[find(info.=="veldef"),3]
    pmra = info[find(info.=="pmra"),3]
    pmdec =  info[find(info.=="pmdec"),3]
    pmra_err =  info[find(info.=="pmra_err"),3]
    pmdec_err =  info[find(info.=="pmdec_err"),3]
    parallax =  info[find(info.=="parallax"),3]
    para_err =   info[find(info.=="para_err"),3]
    spectyp =  info[find(info.=="spectyp"),3]
    return target_id,target,raep0,decep0,equinox,ra_err,dec_err,sysvel,veltyp,veldef,pmra,pmdec,pmra_err,pmdec_err,parallax,para_err,spectyp
 end

function read_comb_file(comb_file)
    info=readdlm(comb_file)
    name =info[find(info.=="name"),3]
    int_trans = info[find(info.=="int_trans"),3]
    vis = info[find(info.=="vis"),3]
    n_pix_fringe =info[find(info.=="n_pix_fringe"),3]
    n_pix_photometry = info[find(info.=="n_pix_photometry"),3]
    flux_frac_photometry = info[find(info.=="flux_frac_photometry"),3]
    flux_frac_fringes = info[find(info.=="flux_frac_fringes"),3]
    throughput_photometry = info[find(info.=="throughput_photometry"),3]
    throughput_fringes = info[find(info.=="throughput_fringes"),3]
    n_splits = info[find(info.=="n_splits"),3]
    read_noise = info[find(info.=="read_noise"),3]
    quantum_efficiency = info[find(info.=="quantum_efficiency"),3]
    v2_cal_err = info[find(info.=="v2_cal_err"),3]
    phase_cal_err = info[find(info.=="phase_cal_err"),3]
    incoh_int_time = info[find(info.=="incoh_int_time"),3]
    return info,name,int_trans,vis,n_pix_fringe,n_pix_photometry,flux_frac_photometry,flux_frac_fringes,throughput_photometry,throughput_fringes,n_splits,read_noise,quantum_efficiency,v2_cal_err,phase_cal_err,incoh_int_time
 end

function read_wave_file(wave_file)
    info=readdlm(wave_file)
    combiner =info[find(info.=="combiner"),3]
    mode = info[find(info.=="mode"),3]
    λ=info[3:length(info[1]),1]
    δλ=info[3:length(info[1]),2]
    return combiner,mode,λ,δλ
 end

#user input
δ=+18.594234/180*pi #ce tau 10 mas 32w 0.37 pix

#include("chara_config.jl")
#include("npoi_config.jl")
hour_angles = linspace(-6,6,20);

function simulate_ha(array_config_file,obsv_info_file,comb_file,wave_file,hour_angles,δ,image_file,pixsize,outfilename)
    #get from config file: N, station_xyz,λ
     array_name,lat,lon,alt,coord,throughput,wind_speed,r0,ntel,tel_diams,tel_gain,tel_names,sta_names,sta_index,station_xyz=read_array_file(array_config_file)
     target_id,target,raep0,decep0,equinox,ra_err,dec_err,sysvel,veltyp,veldef,pmra,pmdec,pmra_err,pmdec_err,parallax,para_err,spectyp=read_obs_file(obsv_info_file)
     info,name,int_trans,vis,n_pix_fringe,n_pix_photometry,flux_frac_photometry,flux_frac_fringes,throughput_photometry,throughput_fringes,n_splits,read_noise,quantum_efficiency,v2v2_cal_err,phase_cal_err,incoh_int_time=read_comb_file(comb_file)
     combiner,mode,λ,δλ=read_wave_file(wave_file)
    N=ntel[1]
    nhours = length(hour_angles);
    h = hour_angles' * pi / 12;
    l=lat/180*pi;
    nw=length(λ)
    N=ntel[1]
    staxyz=Array{Float64}(3,N);
    for i=1:N
            staxyz[:,i]=station_xyz[i,:];
    end
    nv2,v2_baselines,v2_stations,v2_stations_nonredun,v2_indx,baseline_list,indx=get_v2_baselines(N,station_xyz);
    nt3,t3_baselines,t3_stations,t3_indx_1,t3_indx_2,t3_indx_3,indx=get_t3_baselines(N,station_xyz);
    nuv,uv,u_M,v_M,w_M=get_uv(l,h,λ,δ,v2_baselines)
    v2_indx_M,t3_indx_1_M,t3_indx_2_M,t3_indx_3_M,v2_indx_w,t3_indx_1_w,t3_indx_2_w,t3_indx_3_w=get_uv_indxes(nhours,nuv,nv2,nt3,v2_indx,t3_indx_1,t3_indx_2,t3_indx_3)
    v2_model,v2_model_err,t3amp_model,t3amp_model_err,t3phi_model,t3phi_model_err=get_simulated_image(image_file,pixsize,uv,v2_indx_w,t3_indx_1_w, t3_indx_2_w, t3_indx_3_w)

    #setup arrays for OIFITS format
    sta_names=tel_names
    sta_index=Int64.(collect(linspace(1,N,N)))

    #input telescope data


    target_id_vis2=ones(nv2*nhours)
    time_vis2=ones(nv2*nhours) #change
    mjd_vis2=ones(nv2*nhours)*58297.  #change
    int_time_vis2=ones(nv2*nhours) #exchange
    flag_vis2=fill(false,nv2*nhours,1)
    #need to get vis2,vis2err,u,v,sta_index from DATA
    target_id_t3=ones(nt3*nhours)
    time_t3=ones(nt3*nhours) #change
    mjd_t3=ones(nt3*nhours)*58297.  #change
    int_time_t3=ones(nt3*nhours) #exchange;
    flag_t3=fill(false,(nt3*nhours),1);

    eff_wave= λ
    eff_band = δλ

    ucoord_vis2 = u_M[v2_indx_M]
    vcoord_vis2 = v_M[v2_indx_M]
    u1coord = u_M[t3_indx_1_M]
    v1coord = v_M[t3_indx_1_M]
    u2coord = u_M[t3_indx_2_M]
    v2coord = v_M[t3_indx_2_M]

    v2_model=reshape(reshape(v2_model,(nhours,nv2,nw)),(nhours*nv2,nw));
    v2_model_err=reshape(reshape(v2_model_err,(nhours,nv2,nw)),(nhours*nv2,nw));

    t3amp_model=reshape(reshape(t3amp_model,(nhours,nt3,nw)),(nhours*nt3,nw));
    t3amp_model_err=reshape(reshape(t3amp_model_err,(nhours,nt3,nw)),(nhours*nt3,nw));

    t3phi_model=reshape(reshape(t3phi_model,(nhours,nt3,nw)),(nhours*nt3,nw));
    t3phi_model_err=reshape(reshape(t3phi_model_err,(nhours,nt3,nw)),(nhours*nt3,nw));

    v2_model_stations=repmat(v2_stations,1,nhours);
    t3_model_stations=repmat(t3_stations,1,nhours);
    v2_model=transpose(v2_model);
    v2_model_err=transpose(v2_model_err);
    t3amp_model=transpose(t3amp_model);
    t3amp_model_err=transpose(t3amp_model_err);
    t3phi_model=transpose(t3phi_model);
    t3phi_model_err=transpose(t3phi_model_err);

    oi_array=[tel_names,sta_names,sta_index,tel_diams,staxyz];
    target_array=[convert(Int16,target_id[1]),convert(String,target[1]),convert(Float64,raep0[1]),convert(Float64,decep0[1]),convert(Int16,equinox[1]),convert(Int16,ra_err[1]),convert(Int16,dec_err[1]),convert(Int16,sysvel[1]),convert(String,veltyp[1]),convert(String,veldef[1]),convert(Int16,pmra[1]),convert(Int16,pmdec[1]),convert(Int16,pmra_err[1]),convert(Int16,pmdec_err[1]),convert(Float64,parallax[1]),convert(Float64,para_err[1]),convert(String,spectyp[1])];
    wave_array=[eff_wave,eff_band];
    vis2_array=[target_id_vis2,time_vis2,mjd_vis2,int_time_vis2,v2_model,v2_model_err,vec(ucoord_vis2),vec(vcoord_vis2),v2_model_stations,flag_vis2];
    t3_array=[target_id_t3,time_t3,mjd_t3,int_time_t3,t3amp_model,t3amp_model_err,t3phi_model,t3phi_model_err,vec(u1coord),vec(v1coord),vec(u2coord),vec(v2coord),t3_model_stations,flag_t3];


    f = fits_create_file(outfilename);
    write_oi_header(f,1);
    write_oi_array(f,oi_array);
    write_oi_target(f,target_array);
    write_oi_wavelength(f,wave_array);
    write_oi_vis2(f,vis2_array);
    write_oi_t3(f,t3_array);
    fits_close_file(f);

end
