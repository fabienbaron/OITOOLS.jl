using NFFT
using DelimitedFiles
using FITSIO
using FITSIO.Libcfitsio

struct facility_info
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

struct obsv_info
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

struct combiner_info
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

struct wave_info
    combiner_name::Array{Any,1}
    combiner_mode::Array{Any,1}
    lam::Array{Any,1}
    del_lam::Array{Any,1}
end

struct error_struct
    v2_multit::Float64
    v2_addit::Float64
    t3amp_multit::Float64
    t3amp_addit::Float64
    t3phi_multit::Float64
    t3phi_addit::Float64
end

function read_facility_file(facility_file,facility_struct)
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
    facility_out=facility_struct(facility_name,latitude,longitude,altitude,coordinates,facility_throughput,windspeed,r_0,n_tel,teldiams,telgain,telnames,stanames,staindex,staxyz)
    return facility_out
 end



function read_obs_file(obsv_file,obs_structure)
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
    obs_out=obs_structure(targetid,target_name,raep_0,decep_0,equi,raerr,decerr,sys_vel,vel_typ,vel_def,pm_ra,pm_dec,pm_ra_err,pm_dec_err,para,par_err,spec_typ)
    return obs_out
 end

function read_comb_file(comb_file,comb_structure)
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
    comb_out=comb_structure(name,int_trans,vis,n_pix_fringe,n_pix_photometry,flux_frac_photometry,flux_frac_fringes,throughput_photometry,throughput_fringes,n_splits,read_noise,quantum_efficiency,v2_cal_err,phase_cal_err,incoh_int_time)
    return comb_out
 end

function read_wave_file(wave_file,wave_structure)
    wave=readdlm(wave_file)
    combiner =wave[(LinearIndices(wave.=="combiner"))[findall(wave.=="combiner")],3]
    mode = wave[(LinearIndices(wave.=="mode"))[findall(wave.=="mode")],3]
    lambda=wave[3:size(wave)[1],1]
    del_lambda=wave[3:size(wave)[1],2]
    wave_out=wave_structure(combiner,mode,lambda,del_lambda)
    return wave_out
 end

function define_errors(errorstruct,v2mult,v2add,t3ampmult,t3ampadd,t3phimult,t3phiadd)
        error_out=errorstruct(v2mult,v2add,t3ampmult,t3ampadd,t3phimult,t3phiadd)
        return error_out
end


function writefits_aspro(data, fitsfile;res=0.05)
    """
    res should be input as mas
    """
    f = FITS(fitsfile, "w");
    header = FITSHeader(["CDELT1","CDELT2","CRVAL1","CRVAL2","CRPIX1","CRPIX2"],[-(res/1000.0)/(206265.0),(res/1000.0)/(206265.0),0.0,0.0,(size(data)[1]/2),(size(data)[1]/2)],["Radians per Pixel","Radians per Pixel","X-coordinate of reference pixel","Y-coordinate of reference pixel","reference pixel in X","reference pixel in Y"])
    write(f, data,header=header);
    close(f);
end

function updatefits_aspro(fitsfile_in,fitsfile_out;res=0.05)
    """
    res should be input as mas
    """
    f = FITS(fitsfile_in);
    data = read(f[1])
    header = read_header(f[1])
    header["CDELT1"] = -(res/1000.0)/(206265.0)
    header["CDELT2"] = (res/1000.0)/(206265.0)
    header["CRVAL1"] = 0.0
    header["CRVAL2"] = 0.0
    header["CRPIX1"] = (size(data)[1]/2)
    header["CRPIX2"] = (size(data)[1]/2)
    set_comment!(header,"CDELT1","Radians per Pixel")
    set_comment!(header,"CDELT2","Radians per Pixel")
    set_comment!(header,"CRVAL1","X-coordinate of reference pixel")
    set_comment!(header,"CRVAL2","Y-coordinate of reference pixel")
    set_comment!(header,"CRPIX1","reference pixel in X")
    set_comment!(header,"CRPIX2","reference pixel in Y")
    fout = FITS(fitsfile_out,"w")
    write(fout,data,header=header);
    close(f)
    close(fout)
end



function hour_angle_calc(dates::Union{Array{Any,2},Array{Float64,2}},longitude::Float64, ra::Float64;dst="no",ldir="W",timezone="UTC")

"""
This function calculates and returns the hour angle for the desired object given a RA, time, and longitude
of observations.  The program assumes UTC but changing the timezone argument to "local" will adjust the
input times accordingly.

Arguments:
Manual Inputs:
* date [2018  3 5 21 13 56.7; 2018 3 5 21 14 12.7] day,month,year,hour,min,sec: Correction to UTC is done within code if necessary/dictated.
* longitude: degree (in decimal form) value of the longitudinal location of the telescope.
* ra: [h,m,s] array of the right ascension of the desired object.

Optional Inputs:
* dst: ["yes","no"] whether daylight savings time needs to be accounted for.
* ldir: ["W","E"] defines the positive direction for change in longitude (should be kept to W for best applicability).
* timezone: ["UTC","local"] states whether the time inputs are UTC or not.  If not the code will adjust the times as needed with the given longitude.

Accuracy:
* With preliminary testing the LST returned is accurate to within a few minutes when compared with other calculators (MORE TESTING AND QUANTITATIVE ERROR NEEDS TO BE ESTABLISHED).
"""


if ldir == "W"
    alpha = -1.
elseif ldir == "E"
    alpha = 1.
end

years=Int.(dates[:,1])
months=Int.(dates[:,2])
days = Int.(dates[:,3])
hours=Int.(dates[:,4])
minutes=Int.(dates[:,5])
seconds=Float64.(dates[:,6])

#rah=ra[1]
#ram=ra[2]
#ras=ra[3]

#raht = rah+(ram/60.)+(ras/3600.) #Converts the RA array into a hour decimal
raht=ra
h_ad = alpha*longitude/15 #Measures the hours offset due to longitude


if timezone != "UTC"
    if dst=="no"
        hours .-= h_ad
    elseif dst=="yes"
        h_ad += 1
        hours .-= h_ad
    end
    wrap_hours_over = findall(hours.>24)
    hours[wrap_hours_over] .-= 24
    days[wrap_hours_over] .+= 1

    wrap_hours_under = findall(hours.<0)
    hours[wrap_hours_under] +=24
    days[wrap_hours_under] -= 1
end

#Below: the code calculates first the Julian Date given the input time and then determines the GMST based
#on this JD.  This is then converted to LST and finally to Local Hour Angle (LHA or HA).  The final result
#is in terms of hours for both LST and HA.
jdn = floor.((1461*(years .+4800 .+(months.-14)/12))/4+(367*(months .-2-12*((months.-14)/12)))/12-(3*((years .+4900 +(months.-14)/12)/100))/4 .+days.-32075)
jdn = jdn + ((hours.-12)/24)+(minutes/1440)+(seconds/86400)
jd0 = jdn .-2451545.0
#gmst = 6.697374558 + (0.06570982441908*(jd0-(hour/24.))) + (1.00273790935*hour) + 0.000026*(jd0/36525)^2
t = jd0/36525.0
H = 24110.54841 .+ 8640184.812866*t + 0.093104*(t.^2) - (6.2*10.0^(-6))*(t.^3)
w = 1.00273790935 .+ (5.9*10.0^(-11)*t)
tt = hours*3600 + minutes*60 + seconds
gmst = H + w.*tt #seconds
gmst = gmst/3600 #hours

gmst_over = findall(gmst.>24)
gmst[gmst_over] -= (24*floor.(gmst[gmst_over]/24))
lst = gmst .+ h_ad

lst_under = findall(lst.<0)
lst_over = findall(lst.>24)
lst[lst_under] .+= 24
lst[lst_over] -= (24*floor.(lst[lst_over]/24))

hour_angle = lst .-raht
return lst,hour_angle
end





function opd_limits(dec,ha,base;latitude=34.224722)
    """
    dec should be input as [deg,am,as]
    ha should be in hours and may be input as array
    base is the difference in baselines (delX i, delY j, delZ k) where i,j,k
        point respecitvely towards the East, North, and Meridian(?).
    res should be input as mas
    latitude is latitude of chara in degrees
    """
    ha = ha .* 15.0
    dec = dec[1]+(dec[2]/60.0)+(dec[3]/3600.0)
    dtr = pi/180.0 #degrees to radians
    opd = zeros(length(ha))
    for i in range(1,length(ha))
        alt = asin.(sin(dec*dtr)*sin(latitude*dtr)+(cos(dec*dtr)*cos(latitude*dtr)*cos.(ha[i]*dtr)))/dtr #returns alt in degree
        az = atan((-1.0*cos(dec*dtr)*sin.(ha[i]*dtr))/(sin(dec*dtr)*cos(latitude*dtr)-(cos(dec*dtr)*cos.(ha[i]*dtr)*sin(latitude*dtr))))/dtr
        s = [cos.(alt*dtr)*cos.(az*dtr),cos.(alt*dtr)*sin.(az*dtr),sin.(alt*dtr)]
        opd[i] = dot(base,s) #in meters (?)
    end
    return opd
end



#CODES FOR SIMULATING OIFITS BASED ON INPUT IMAGE AND EITHER INPUT OIFITS OR HOUR ANGLES

#Functions used in main Functions
function cvis_to_t3_conj(cvis, indx1, indx2, indx3)
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
    v2_stations_nonredun=Array{Int64}(undef,2,nv2);
    v2_indx      = Array{Int64}(undef,nv2);
    baseline_list = Array{String}(undef,nv2);
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
        baseline_list[i]=string(tel_names[v2_stations[1,i]],"-",tel_names[v2_stations[2,i]])
    end
    return nv2,v2_baselines,v2_stations,v2_stations_nonredun,v2_indx,baseline_list,ind
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
    return nt3,t3_baselines,t3_stations,t3_indx_1,t3_indx_2,t3_indx_3,ind
end

function get_uv(l,h,λ,δ,v2_baselines,nhours)
    #Use following expression only if there are missing baselines somewhere
    #vis_baselines = hcat(v2_baselines, t3_baselines[:, 1, :], t3_baselines[:, 2, :], t3_baselines[:, 3, :])
    #Expression to use for pure simulation where the full complement of v2 and t3 are created
    vis_baselines = deepcopy(v2_baselines)
    nuv = size(vis_baselines, 2);
    # Now compute the UV coordinates, again according to the APIS++ standards.
    u_M = -sin.(l)*sin.(h) .*vis_baselines[1,:]  .+ cos.(h) .* vis_baselines[2,:]+cos.(l)*sin.(h).*vis_baselines[3,:];

    v_M = (sin.(l)*sin(δ)*cos.(h).+cos.(l)*cos(δ)) .* vis_baselines[1,:] + sin(δ)*sin.(h) .* vis_baselines[2,:]+(-cos.(l)*cos.(h)*sin(δ).+sin.(l)*cos(δ)) .* vis_baselines[3,:];

    w_M =  (sin.(l)*cos(δ)* cos.(h).+cos.(l)*sin(δ)).* vis_baselines[1,:] - cos(δ)* sin.(h) .* vis_baselines[2,:]  .+ (cos.(l)*cos.(h)*cos(δ).+sin.(l)sin(δ)) .* vis_baselines[3,:];
    # proj baselines to (uv wav)
    u_M=-u_M
    v_M=v_M
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
    t3_baseline = (sqrt.(uv[1,t3_indx_1_w].^2 + uv[2,t3_indx_1_w].^2).*
        sqrt.(uv[1,t3_indx_2_w].^2 + uv[2,t3_indx_2_w].^2).*
        sqrt.(uv[1,t3_indx_3_w].^2 + uv[2,t3_indx_3_w].^2)).^(1 ./ 3.);
    return v2_indx_M,t3_indx_1_M,t3_indx_2_M,t3_indx_3_M,v2_indx_w,t3_indx_1_w,t3_indx_2_w,t3_indx_3_w
 end

 function compute_observables(image,pixsize,uv,v2_indx_w,t3_indx_1_w, t3_indx_2_w, t3_indx_3_w,error_struc)
     nx = size(image,1);
     x = vec(image)/sum(image);
     dft = setup_dft(uv, nx, pixsize);
     cvis_model = image_to_cvis_dft(x, dft);
     v2_model = cvis_to_v2(cvis_model, v2_indx_w);
     t3_model, t3amp_model, t3phi_model = cvis_to_t3_conj(cvis_model, t3_indx_1_w, t3_indx_2_w, t3_indx_3_w);
     #Add noise
     v2_model_err = error_struc.v2_multit*v2_model .+ error_struc.v2_addit
     t3amp_model_err = error_struc.t3amp_multit*t3amp_model .+ error_struc.t3amp_addit
     t3phi_model_err = zeros(length(t3phi_model)) .+ error_struc.t3phi_addit # degree  -- there is another way of setting this with Haniff formula
     v2_model    += v2_model_err.*randn(length(v2_model))
     t3amp_model += t3amp_model_err.*randn(length(t3amp_model))
     t3phi_model += t3phi_model_err.*randn(length(t3phi_model))
    return v2_model,v2_model_err,t3amp_model,t3amp_model_err,t3phi_model,t3phi_model_err
 end

function prep_arrays(_info)
    oi_array=[tel_names,sta_names,sta_index,tel_diams,staxyz];
    target_array=[target_id,target,raep0,decep0,equinox,ra_err,dec_err,sysvel,veltyp,veldef,pmra,pmdec,pmra_err,pmdec_err,parallax,para_err,spectyp];
    wave_array=[eff_wave,eff_band];
    vis2_array=[target_id_vis2,time_vis2,mjd_vis2,int_time_vis2,v2_model,v2_model_err,vec(ucoord_vis2),vec(vcoord_vis2),v2_model_stations,flag_vis2];
    t3_array=[target_id_t3,time_t3,mjd_t3,int_time_t3,t3amp_model,t3amp_model_err,t3phi_model,t3phi_model_err,vec(u1coord),vec(v1coord),vec(u2coord),vec(v2coord),t3_model_stations,flag_t3];
    return oi_array,target_array,wave_array,vis2_array,t3oi_array
 end


#user input
#δ=+18.594234/180*pi #ce tau 10 mas 32w 0.37 pix

#include("chara_config.jl")
#include("npoi_config.jl")
#hour_angles = range(-6,6,20);

function simulate_ha(facility,observatory,combiner,wave_info_out,hour_angles,image_file,pixsize,errors,outfilename)
    #simulate an observation using input hour angles, info about array and combiner, and input image


    ntel=facility.ntel[1] #✓
    nhours = length(hour_angles); #✓
    h = hour_angles' .* pi / 12; #✓
    δ=observatory.decep0[1]/180*pi #✓
    l=facility.lat[1]/180*pi; #✓
    λ=wave_info_out.lam#✓
    δλ=wave_info_out.del_lam #✓
    nw=length(λ)

    station_xyz=Array{Float64}(undef,ntel,3) #✓
    staxyz=Array{Float64}(undef,3,ntel); #✓
    for i=1:ntel[1] #✓
        station_xyz[i,1:3]=facility.sta_xyz[(i*3-2):i*3]#✓
    end


    for i=1:ntel #✓
            staxyz[:,i]=station_xyz[i,:];#✓
    end
#✓

    nv2,v2_baselines,v2_stations,v2_stations_nonredun,v2_indx,baseline_list,ind=get_v2_baselines(ntel,station_xyz,facility.tel_names);
    nt3,t3_baselines,t3_stations,t3_indx_1,t3_indx_2,t3_indx_3,ind=get_t3_baselines(ntel,station_xyz,v2_stations);
    nuv,uv,u_M,v_M,w_M=get_uv(l,h,λ,δ,v2_baselines,nhours)
    v2_indx_M,t3_indx_1_M,t3_indx_2_M,t3_indx_3_M,v2_indx_w,t3_indx_1_w,t3_indx_2_w,t3_indx_3_w=get_uv_indxes(nhours,nuv,nv2,nt3,v2_indx,t3_indx_1,t3_indx_2,t3_indx_3,nw,uv)
#✓

#simulate observation using input image file and pixsize
    x = readfits(image_file)
    #flip image for an unknown reason
    x=x[end:-1:1,:]
    # TODO: "wavelength to image" vector, "epoch to image" vector
    # TODO: test the size of x, handle temporal and polychromatic loops
    v2_model,v2_model_err,t3amp_model,t3amp_model_err,t3phi_model,t3phi_model_err=compute_observables(x,pixsize,uv,v2_indx_w,t3_indx_1_w, t3_indx_2_w, t3_indx_3_w,errors)

    #setup arrays for OIFITS format
    sta_names=facility.tel_names
    sta_index=Int64.(collect(range(1,step=1,length=ntel)))

    #input telescope data
    target_id_vis2=ones(nv2*nhours).*observatory.target_id[1]
    time_vis2=ones(nv2*nhours) #change
    mjd_vis2=ones(nv2*nhours)*58297.  #change
    int_time_vis2=ones(nv2*nhours) #exchange
    flag_vis2=fill(false,nv2*nhours,1)
    #need to get vis2,vis2err,u,v,sta_index from DATA
    target_id_t3=ones(nt3*nhours).*observatory.target_id[1]
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

    v2_model_stations=repeat(v2_stations,1,nhours);
    t3_model_stations=repeat(t3_stations,1,nhours);

    v2_model=transpose(v2_model);
    v2_model_err=transpose(v2_model_err);
    t3amp_model=transpose(t3amp_model);
    t3amp_model_err=transpose(t3amp_model_err);
    t3phi_model=transpose(t3phi_model);
    t3phi_model_err=transpose(t3phi_model_err);

    oi_array=[facility.tel_names,sta_names,sta_index,facility.tel_diams,staxyz];
    target_array=[observatory.target_id[1],observatory.target[1],observatory.raep0[1],observatory.decep0[1],observatory.equinox[1],observatory.ra_err[1],observatory.dec_err[1],observatory.sysvel[1],observatory.veltyp[1],observatory.veldef[1],observatory.pmra[1],observatory.pmdec[1],observatory.pmra_err[1],observatory.pmdec_err[1],observatory.parallax[1],observatory.para_err[1],observatory.spectyp[1]];
    wave_array=[λ,δλ];
    vis2_array=[target_id_vis2,time_vis2,mjd_vis2,int_time_vis2,v2_model,v2_model_err,vec(ucoord_vis2),vec(vcoord_vis2),v2_model_stations,flag_vis2];
    t3_array=[target_id_t3,time_t3,mjd_t3,int_time_t3,t3amp_model,t3amp_model_err,t3phi_model,t3phi_model_err,vec(u1coord),vec(v1coord),vec(u2coord),vec(v2coord),t3_model_stations,flag_t3];


    f = fits_create_file(outfilename);
    write_oi_header(f,1);
    write_oi_array(f,[facility.tel_names,sta_names,facility.sta_index,facility.tel_diams,staxyz]);
    #write_oi_target(f,[convert(Int16,target_id[1]),convert(String,target[1]),convert(Float64,raep0[1]),convert(Float64,decep0[1]),convert(Int16,equinox[1]),convert(Int16,ra_err[1]),convert(Int16,dec_err[1]),convert(Int16,sysvel[1]),convert(String,veltyp[1]),convert(String,veldef[1]),convert(Int16,pmra[1]),convert(Int16,pmdec[1]),convert(Int16,pmra_err[1]),convert(Int16,pmdec_err[1]),convert(Float64,parallax[1]),convert(Float64,para_err[1]),convert(String,spectyp[1])]);
    write_oi_target(f,target_array)
    write_oi_wavelength(f,wave_array);
    write_oi_vis2(f,[target_id_vis2,time_vis2,mjd_vis2,int_time_vis2,v2_model,v2_model_err,vec(ucoord_vis2),vec(vcoord_vis2),v2_model_stations,flag_vis2]);
    write_oi_t3(f,[target_id_t3,time_t3,mjd_t3,int_time_t3,t3amp_model,t3amp_model_err,t3phi_model,t3phi_model_err,vec(u1coord),vec(v1coord),vec(u2coord),vec(v2coord),t3_model_stations,flag_t3]);
    fits_close_file(f);

end

function simulate_obs(oifitsin,outfilename,fitsfiles,pixsize;dft=false,nfft=false)
    #simulate observation from input oifits and input image.
    if typeof(fitsfiles)== String
        num_files = 1.0
    else
        num_files=length(fitsfiles)
    end


    data = (readoifits(oifitsin))[1,1]; # data can be split by wavelength, time, etc.
    oifits=FITS(oifitsin);

    #setup simulation
    nuv = data.nuv

    #x = (read((FITS(fitsfiles))[1])); x=x[:,end:-1:1]; nx = (size(x))[1]; x=vec(x)/sum(x); #read in images
    x = (read((FITS(fitsfiles))[1])); nx = (size(x))[1]; x=vec(x)/sum(x); #read in images
    #nfft_plan = setup_nfft(-uv, nx, pixsize)
    #cvis_model = image_to_cvis_nfft(x, nfft_plan)
    cvis_model=[]
    if dft==true
        ft = setup_dft(data, nx, pixsize);
        cvis_model = image_to_cvis_dft(x, ft);
    end

    if nfft==true
        ft = setup_nfft(data, nx, pixsize);
        cvis_model = image_to_cvis_nfft(x, ft);
    end
    #setup new file
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


    #get OI_ARRY details (currently only handles one Array)
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
    #copy_oi_target(f,oifits[ioiarraytables[length(ioiarraytables)]+1],target_array);
    #copy_oi_target(f,oifits[ioiarraytables[length(ioiarraytables)]-1],target_array); comment out to deal with EXTNAME issue check this
    copy_oi_target(f,oifits[iotargettables[length(iotargettables)]],target_array);
    #GET WAVE TABLES
    for itable=1:length(iwavetables)
        eff_wave=read(oifits[iwavetables[itable]],"EFF_WAVE")
        eff_band=read(oifits[iwavetables[itable]],"EFF_BAND")
        wave_array=[eff_wave,eff_band];
        copy_oi_wavelength(f,oifits[iwavetables[itable]],wave_array);
    end

    #GET V2 INFO

    v2_model = cvis_to_v2(cvis_model, data.indx_v2 )# based on uv points
    v2_model_err=data.v2_err
    # v2_model_err = v2_model_true./abs.(data.v2./data.v2_err)

    v2_model += v2_model_err.*randn(length(v2_model))
    #Now to fill the tables_
    #uvis_lam=[];
    #vvis_lam=[];
    #v2_stations=[];
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
    t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2, data.indx_t3_3);

    #t3amp_model_err =0.007*t3amp_model+1e-6
    t3amp_model_err=data.t3amp_err
    #t3amp_model_err=abs.(t3amp_model./(data.t3amp./data.t3amp_err))
    t3amp_model += t3amp_model_err.*randn(length(t3amp_model))

    #t3phi_model_err = zeros(length(t3phi_model))+2. # degree  -- there is another way of setting this with Haniff formula
    t3phi_model_err=data.t3phi_err
    #t3phi_model_err=abs.(t3phi_model./(data.t3phi./data.t3phi_err))
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
