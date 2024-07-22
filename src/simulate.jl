using NFFT
using DelimitedFiles
using FITSIO
using CFITSIO

mutable struct facility_info
    facility_name::Array{Any,1}
    lat::Array{Any,1}
    lon::Array{Any,1}
    atl::Array{Any,1}
    coord::Array{Any,1}
    throughput::Array{Any,1}
    wind_speed::Array{Any,1}
    r0::Array{Any,1}
    ntel::Array{Any,1}
    tel_diams::Array{Any,2}
    tel_gain::Array{Any,2}
    tel_names::Array{Any,2}
    sta_names::Array{Any,2}
    sta_index::Array{Any,2}
    sta_xyz::Array{Float64,2}
end

mutable struct obsv_info
    target_id::Array{Any,1}
    target::Array{Any,1}
    raep0::Array{Any,1}
    decep0::Array{Any,1}
    equinox::Array{Any,1}
    ra_err::Array{Any,1}
    dec_err::Array{Any,1}
    sysvel::Array{Any,1}
    veltyp::Array{Any,1}
    veldef::Array{Any,1}
    pmra::Array{Any,1}
    pmdec::Array{Any,1}
    pmra_err::Array{Any,1}
    pmdec_err::Array{Any,1}
    parallax::Array{Any,1}
    para_err::Array{Any,1}
    spectyp::Array{Any,1}
end

mutable struct combiner_info
    comb_name::Array{Any,1}
    int_trans::Array{Any,1}
    vis::Array{Any,1}
    n_pix::Array{Any,1}
    n_pix_photometry::Array{Any,1}
    flux_frac_photometry::Array{Any,1}
    flux_frac_fringes::Array{Any,1}
    throughput_photometry::Array{Any,1}
    throughput_fringes::Array{Any,1}
    n_splits::Array{Any,1}
    read_noise::Array{Any,1}
    quantum_efficiency::Array{Any,1}
    v2_cal_err::Array{Any,1}
    phase_cal_err::Array{Any,1}
    incoh_int_time::Array{Any,1}
end

mutable struct wave_info
    combiner_name::String
    combiner_mode::String
    λ::Array{Float64,1}
    δλ::Array{Float64,1}
end

mutable struct error_info
    v2_multit::Float64
    v2_addit::Float64
    t3amp_multit::Float64
    t3amp_addit::Float64
    t3phi_multit::Float64
    t3phi_addit::Float64
end

function read_facility_file(facility_file)
    facility=readdlm(facility_file)
    facility_name=facility[(LinearIndices(facility.=="name"))[findall(facility.=="name")],3]
    latitude=facility[(LinearIndices(facility.=="lat"))[findall(facility.=="lat")],3]
    longitude=facility[(LinearIndices(facility.=="lon"))[findall(facility.=="lon")],3]
    altitude=facility[(LinearIndices(facility.=="alt"))[findall(facility.=="alt")],3]
    coordinates=facility[(LinearIndices(facility.=="coord"))[findall(facility.=="coord")],3]
    facility_throughput=facility[(LinearIndices(facility.=="throughput"))[findall(facility.=="throughput")],3]
    windspeed=facility[(LinearIndices(facility.=="wind_speed"))[findall(facility.=="wind_speed")],3]
    r_0=facility[(LinearIndices(facility.=="r0"))[findall(facility.=="r0")],3]
    n_tel=facility[(LinearIndices(facility.=="ntel"))[findall(facility.=="ntel")],3]
    teldiams=facility[(LinearIndices(facility.=="tel_diams"))[findall(facility.=="tel_diams")],3:n_tel[1]+2]
    telgain=facility[(LinearIndices(facility.=="tel_gain"))[findall(facility.=="tel_gain")],3:n_tel[1]+2]
    telnames=facility[(LinearIndices(facility.=="tel_names"))[findall(facility.=="tel_names")],3:n_tel[1]+2]
    stanames=facility[(LinearIndices(facility.=="sta_names"))[findall(facility.=="sta_names")],3:n_tel[1]+2]
    staindex=facility[(LinearIndices(facility.=="sta_index"))[findall(facility.=="sta_index")],3:n_tel[1]+2]
    staxyz=facility[(LinearIndices(facility.=="sta_xyz"))[findall(facility.=="sta_xyz")],3:n_tel[1]*3+2]
    facility_out=facility_info(facility_name,latitude,longitude,altitude,coordinates,facility_throughput,windspeed,r_0,n_tel,teldiams,telgain,telnames,stanames,staindex,staxyz)
    return facility_out
 end



function read_obs_file(obsv_file)
    obs=readdlm(obsv_file)
    targetid=obs[(LinearIndices(obs.=="target_id"))[findall(obs.=="target_id")],3]
    target_name=obs[(LinearIndices(obs.=="target"))[findall(obs.=="target")],3]
    raep_0=obs[(LinearIndices(obs.=="raep0"))[findall(obs.=="raep0")],3]
    decep_0=obs[(LinearIndices(obs.=="decep0"))[findall(obs.=="decep0")],3]
    equi=obs[(LinearIndices(obs.=="equinox"))[findall(obs.=="equinox")],3]
    raerr=obs[(LinearIndices(obs.=="ra_err"))[findall(obs.=="ra_err")],3]
    decerr = obs[(LinearIndices(obs.=="dec_err"))[findall(obs.=="dec_err")],3]
    sys_vel = obs[(LinearIndices(obs.=="sysvel"))[findall(obs.=="sysvel")],3]
    vel_typ =  obs[(LinearIndices(obs.=="veltyp"))[findall(obs.=="veltyp")],3]
    vel_def = obs[(LinearIndices(obs.=="veldef"))[findall(obs.=="veldef")],3]
    pm_ra = obs[(LinearIndices(obs.=="pmra"))[findall(obs.=="pmra")],3]
    pm_dec =  obs[(LinearIndices(obs.=="pmdec"))[findall(obs.=="pmdec")],3]
    pm_ra_err =  obs[(LinearIndices(obs.=="pmra_err"))[findall(obs.=="pmra_err")],3]
    pm_dec_err =  obs[(LinearIndices(obs.=="pmdec_err"))[findall(obs.=="pmdec_err")],3]
    para =  obs[(LinearIndices(obs.=="parallax"))[findall(obs.=="parallax")],3]
    par_err =   obs[(LinearIndices(obs.=="para_err"))[findall(obs.=="para_err")],3]
    spec_typ =  obs[(LinearIndices(obs.=="spectyp"))[findall(obs.=="spectyp")],3]
    obs_out=obsv_info(targetid,target_name,raep_0,decep_0,equi,raerr,decerr,sys_vel,vel_typ,vel_def,pm_ra,pm_dec,pm_ra_err,pm_dec_err,para,par_err,spec_typ)
    return obs_out
 end

function read_comb_file(comb_file)
    comb=readdlm(comb_file)
    name =comb[(LinearIndices(comb.=="name"))[findall(comb.=="name")],3]
    int_trans = comb[(LinearIndices(comb.=="int_trans"))[findall(comb.=="int_trans")],3]
    vis = comb[(LinearIndices(comb.=="vis"))[findall(comb.=="vis")],3]
    n_pix_fringe =comb[(LinearIndices(comb.=="n_pix_fringe"))[findall(comb.=="n_pix_fringe")],3]
    n_pix_photometry = comb[(LinearIndices(comb.=="n_pix_photometry"))[findall(comb.=="n_pix_photometry")],3]
    flux_frac_photometry = comb[(LinearIndices(comb.=="flux_frac_photometry"))[findall(comb.=="flux_frac_photometry")],3]
    flux_frac_fringes = comb[(LinearIndices(comb.=="flux_frac_fringes"))[findall(comb.=="flux_frac_fringes")],3]
    throughput_photometry = comb[(LinearIndices(comb.=="throughput_photometry"))[findall(comb.=="throughput_photometry")],3]
    throughput_fringes = comb[(LinearIndices(comb.=="throughput_fringes"))[findall(comb.=="throughput_fringes")],3]
    n_splits = comb[(LinearIndices(comb.=="n_splits"))[findall(comb.=="n_splits")],3]
    read_noise = comb[(LinearIndices(comb.=="read_noise"))[findall(comb.=="read_noise")],3]
    quantum_efficiency = comb[(LinearIndices(comb.=="quantum_efficiency"))[findall(comb.=="quantum_efficiency")],3]
    v2_cal_err = comb[(LinearIndices(comb.=="v2_cal_err"))[findall(comb.=="v2_cal_err")],3]
    phase_cal_err = comb[(LinearIndices(comb.=="phase_cal_err"))[findall(comb.=="phase_cal_err")],3]
    incoh_int_time = comb[(LinearIndices(comb.=="incoh_int_time"))[findall(comb.=="incoh_int_time")],3]
    return combiner_info(name,int_trans,vis,n_pix_fringe,n_pix_photometry,flux_frac_photometry,flux_frac_fringes,throughput_photometry,throughput_fringes,n_splits,read_noise,quantum_efficiency,v2_cal_err,phase_cal_err,incoh_int_time)
 end

function read_wave_file(wave_file)
    wave       = readdlm(wave_file)
    combiner   = String.(wave[(LinearIndices(wave.=="combiner"))[findall(wave.=="combiner")],3])[1]
    mode       = String.(wave[(LinearIndices(wave.=="mode"))[findall(wave.=="mode")],3])[1]
    λ          = Float64.(wave[3:size(wave)[1],1])
    δλ         = Float64.(wave[3:size(wave)[1],2])
    return wave_info(combiner,mode,λ,δλ)
 end

function define_errors(v2mult,v2add,t3ampmult,t3ampadd,t3phimult,t3phiadd)
 return error_info(v2mult,v2add,t3ampmult,t3ampadd,t3phimult,t3phiadd)
end


#CODES FOR SIMULATING OIFITS BASED ON INPUT IMAGE AND EITHER INPUT OIFITS OR HOUR ANGLES

#Functions used in main Functions
function vis_to_t3_conj(cvis, indx1, indx2, indx3)
    #get t3 from caculated visibilities
    #because of ordering of v_ij j>i we need to use conjugate (we use conj(t3_13), rather than t3_31)
 # t3 = cvis[indx1].*cvis[indx2].*conj(cvis[indx3]);
  t3 = cvis[indx1].*cvis[indx2].*conj(cvis[indx3]);
  t3amp = abs.(t3);
  t3phi = angle.(t3)*180 ./pi;
  return t3, t3amp, t3phi
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
    for i=1:N+1
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

function get_uv(l, h, λ, δ, baselines)
    #Use following expression only if there are missing baselines somewhere
    # baselines = hcat(v2_baselines, t3_baselines[:, 1, :], t3_baselines[:, 2, :], t3_baselines[:, 3, :])
    #Expression to use for pure simulation where the full complement of v2 and t3 are created
    nhours = length(h);
    nuv = size(baselines, 2);
    # Now compute the UV coordinates, again according to the APIS++ standards.
    u_M = -sin.(l)*sin.(h) .*baselines[1,:]  .+ cos.(h) .* baselines[2,:]+cos.(l)*sin.(h).*baselines[3,:];
    v_M = (sin.(l)*sin(δ)*cos.(h).+cos.(l)*cos(δ)).* baselines[1,:]  + sin(δ)*sin.(h) .*baselines[2,:]   +(-cos.(l)*cos.(h)*sin(δ).+sin.(l)*cos(δ)) .* baselines[3,:];
    w_M = (-sin.(l)*cos(δ)*cos.(h).+cos.(l)*sin(δ)).* baselines[1,:] - cos(δ)*sin.(h) .*baselines[2,:]   + (cos.(l)*cos.(h)*cos(δ).+sin.(l)sin(δ))  .* baselines[3,:];
    # proj baselines to (uv wav)
    u = reshape((1 ./λ)'.*vec(-u_M), (nuv,nhours,length(λ))); #TODO: we have -u_M for the moment, need to fix nfft first
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

function simulate(facility,target,combiner,wave,dates,errors,out_file; image::Union{String, Array{Float64,1}, Array{Float64,2}, Array{Float64, 3}, Array{Float64,4}}="", pixsize::Float64=0.1, model::OImodel=create_model())
    outfilename= string("!", out_file)
    #simulate an observation using input hour angles, info about array and combiner, and input image
    ntel=facility.ntel[1];
    lst, hour_angles = hour_angle_calc(dates,facility.lon[1],target.raep0[1]);
    nhours = length(hour_angles);
    h_rad = hour_angles' .* pi / 12;
    δ = target.decep0[1]/180*pi;
    l = facility.lat[1]/180*pi;
    λ = wave.λ;
    δλ = wave.δλ;
    nw = length(λ)

    station_xyz=zeros(Float64,ntel,3) #✓
    for i=1:ntel
        station_xyz[i,1:3]=facility.sta_xyz[(i*3-2):i*3]#✓
    end

    # Find physical baselines and triangles combinations
    nv2,v2_baselines,v2_stations,v2_indx,baseline_name         = get_v2_baselines(ntel,station_xyz,facility.tel_names);
    nt3,t3_baselines,t3_stations,t3_indx_1,t3_indx_2,t3_indx_3 = get_t3_baselines(ntel,station_xyz,v2_stations);

    # Compute uv coverage: nuv is the number of uv points, uv is (u,v) in Mλ, (u_M, v_M, w_M) are in meters
    nuv, uv, u_M, v_M, w_M = get_uv(l,h_rad,λ,δ,v2_baselines)

    # Setup indexing for OIFITS:  _M -> in meters, _w -> scaled by λ
    v2_indx_M,t3_indx_1_M,t3_indx_2_M,t3_indx_3_M,v2_indx_w,t3_indx_1_w,t3_indx_2_w,t3_indx_3_w = get_uv_indxes(nhours,nuv,nv2,nt3,v2_indx,t3_indx_1,t3_indx_2,t3_indx_3,nw,uv)

    # Compute complex visibilities from either a truth image or a model
    cvis_model = ComplexF64[]
    # Determine if image or model

    if (((typeof(image) == String) && (image !="")) || ((typeof(image) != String) && (image !="")) )
        # Truth image
        # TODO: Polychromatic and dynamic imaging
        x = readfits(image)
        nx = size(x,1);
        x = vec(x)/sum(x);
        ft = setup_nfft(uv, v2_indx_w, v2_indx_w, t3_indx_1_w, t3_indx_2_w, t3_indx_3_w, nx, pixsize);
        cvis_model = image_to_vis(x, ft);
    elseif model.components !=[]
        # Model
        cvis_model = model_to_vis(model, uv, λ);
    else
        error("Bad image or model definition in call to simulate()");
    end

    # Compute true values of observables
    v2_model = vis_to_v2(cvis_model, v2_indx_w);
    t3_model, t3amp_model, t3phi_model = vis_to_t3(cvis_model, t3_indx_1_w, t3_indx_2_w, t3_indx_3_w);

    # Compute uncertainties
    v2_model_err = errors.v2_multit*v2_model .+ errors.v2_addit;
    t3amp_model_err = errors.t3amp_multit*t3amp_model .+ errors.t3amp_addit;
    t3phi_model_err = zeros(length(t3phi_model)) .+ errors.t3phi_addit; # degree  -- there is another way of setting this with Haniff's formula

    # Add errors
    v2_model    += v2_model_err.*randn(length(v2_model));
    t3amp_model += t3amp_model_err.*randn(length(t3amp_model));
    t3phi_model += t3phi_model_err.*randn(length(t3phi_model));


    #setup arrays for OIFITS format
    sta_names=facility.tel_names
    sta_index=Int64.(collect(range(1,step=1,length=ntel)))

    #input telescope data
    target_id_vis2   = ones(nv2*nhours).*target.target_id[1]
    time_vis2        = zeros(nv2*nhours)     # OIFITS v2 requires zeros
    mjd_vis2         = repeat(value.(modified_julian.(dates)), nv2) # TOCHECK (could be transposed)
    int_time_vis2    = ones(Float64,nv2*nhours)      # TODO
    flag_vis2        = fill(false,nv2*nhours,1)
    #need to get vis2,vis2err,u,v,sta_index from DATA
    target_id_t3     = ones(nt3*nhours).*target.target_id[1]
    time_t3          = zeros(nt3*nhours) #change
    mjd_t3           = repeat(value.(modified_julian.(dates)), nt3)  #change
    int_time_t3      = ones(Float64, nt3*nhours)      #  TODO;
    flag_t3          = fill(false,(nt3*nhours),1);

    ucoord_vis2 = u_M[v2_indx_M]
    vcoord_vis2 = v_M[v2_indx_M]
    u1coord = u_M[t3_indx_1_M]
    v1coord = v_M[t3_indx_1_M]
    u2coord = u_M[t3_indx_2_M]
    v2coord = v_M[t3_indx_2_M]

    v2_model          = reshape(reshape(v2_model,(nhours,nv2,nw)),(nhours*nv2,nw))';
    v2_model_err      = reshape(reshape(v2_model_err,(nhours,nv2,nw)),(nhours*nv2,nw))';
    t3amp_model       = reshape(reshape(t3amp_model,(nhours,nt3,nw)),(nhours*nt3,nw))';
    t3amp_model_err   = reshape(reshape(t3amp_model_err,(nhours,nt3,nw)),(nhours*nt3,nw))';
    t3phi_model       = reshape(reshape(t3phi_model,(nhours,nt3,nw)),(nhours*nt3,nw))';
    t3phi_model_err   = reshape(reshape(t3phi_model_err,(nhours,nt3,nw)),(nhours*nt3,nw))';
    v2_model_stations = repeat(v2_stations,1,nhours);
    t3_model_stations = repeat(t3_stations,1,nhours);

    oiarray     = [facility.tel_names,sta_names,facility.sta_index,facility.tel_diams,station_xyz'];
    oitarget    = [target.target_id[1],target.target[1],target.raep0[1],target.decep0[1],target.equinox[1],target.ra_err[1],target.dec_err[1],target.sysvel[1],target.veltyp[1],target.veldef[1],target.pmra[1],target.pmdec[1],target.pmra_err[1],target.pmdec_err[1],target.parallax[1],target.para_err[1],target.spectyp[1]];
    oiwavelength= [λ,δλ];
    oivis2      = [target_id_vis2,time_vis2,mjd_vis2,int_time_vis2,v2_model,v2_model_err,vec(ucoord_vis2),vec(vcoord_vis2),v2_model_stations,flag_vis2];
    oit3        = [target_id_t3,time_t3,mjd_t3,int_time_t3,t3amp_model,t3amp_model_err,t3phi_model,t3phi_model_err,vec(u1coord),vec(v1coord),vec(u2coord),vec(v2coord),t3_model_stations,flag_t3];

    # Write everything
    f = fits_create_file(outfilename);
    write_oi_header(f,1);
    write_oi_array(f,oiarray);
    write_oi_target(f,oitarget)
    write_oi_wavelength(f,oiwavelength);
    write_oi_vis2(f,oivis2);
    write_oi_t3(f,oit3);
    fits_close_file(f);
end

function simulate_from_oifits(in_oifits, out_file; mode="copy_errors", errors=[],  image::Union{String, Array{Float64,1}, Array{Float64,2}, Array{Float64, 3}, Array{Float64,4}}="", pixsize::Float64=0.1, model::OImodel=create_model())
 
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
    elseif model.components !=[]
        # Model
        cvis_model = model_to_vis(model, data.uv, data.uv_lam);
    else
        error("Bad image or model definition in call to simulate()");
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
        nw=length(eff_wave);
        for hour=1:length(uvis)
            for wave=1:nw
                v2_model_table=push!(v2_model_table,v2_model[visindxstart+wave-1]);
                v2_model_err_table=push!(v2_model_err_table,v2_model_err[visindxstart+wave-1]);
            end
            visindxstart=visindxstart+nw
        end
        target_id_vis2=read(oifits[ivistables[itable]],"TARGET_ID");
        time_vis2=read(oifits[ivistables[itable]],"TIME");
        mjd_vis2=read(oifits[ivistables[itable]],"MJD");
        int_time_vis2=read(oifits[ivistables[itable]],"INT_TIME");
        v2_model_stations=read(oifits[ivistables[itable]],"STA_INDEX");
        nrowvis=Int(length(v2_model_table)/nw)
        v2_model_reshape=reshape(v2_model_table,(nw,nrowvis));
        v2_model_err_reshape=reshape(v2_model_err_table,(nw,nrowvis));
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
        nw=length(eff_wave);
        for hour=1:length(u1)
            for wave=1:nw
                t3amp_model_table=push!(t3amp_model_table,t3amp_model[t3indxstart+wave-1]);
                t3amp_model_err_table=push!(t3amp_model_err_table,t3amp_model_err[t3indxstart+wave-1]);
                t3phi_model_table=push!(t3phi_model_table,t3phi_model[t3indxstart+wave-1]);
                t3phi_model_err_table=push!(t3phi_model_err_table,t3phi_model_err[t3indxstart+wave-1]);
            end
            t3indxstart=t3indxstart+nw
        end
        target_id_t3=read(oifits[it3tables[itable]],"TARGET_ID");
        time_t3=read(oifits[it3tables[itable]],"TIME");
        mjd_t3=read(oifits[it3tables[itable]],"MJD");
        int_time_t3=read(oifits[it3tables[itable]],"INT_TIME");
        t3_model_stations=read(oifits[it3tables[itable]],"STA_INDEX");
        nrowt3=Int(length(t3amp_model_table)/nw)
        t3amp_model_reshape=reshape(t3amp_model_table,(nw,nrowt3));
        t3amp_model_err_reshape=reshape(t3amp_model_err_table,(nw,nrowt3));
        t3phi_model_reshape=reshape(t3phi_model_table,(nw,nrowt3));
        t3phi_model_err_reshape=reshape(t3phi_model_err_table,(nw,nrowt3));
        flag_t3=read(oifits[it3tables[itable]],"FLAG")
        t3_array=[target_id_t3,time_t3,mjd_t3,int_time_t3,t3amp_model_reshape,t3amp_model_err_reshape,t3phi_model_reshape,t3phi_model_err_reshape,u1,v1,u2,v2,t3_model_stations,flag_t3'];
        copy_oi_t3(f,oifits[it3tables[itable]],t3_array);
    end
    fits_close_file(f)
end
