
#NEED TO ADD ABILTY TO HAVE DIFFERENT INSNAME IN OUTPUT IF INPUT HAS DIFFERENT INSNAME...basically just copy and pass header.


function simulate_obs(oifitsin,outfilename,fitsfiles,pixsize;dft=false,nfft=true)

    if typeof(fitsfiles)== String
        num_files = 1.0
    else
        num_files=length(fitsfiles)
    end


    data = (readoifits(oifitsin))[1,1]; # data can be split by wavelength, time, etc.
    oifits=FITS(oifitsin);

    #setup simulation
    nuv = data.nuv

    x = (read((FITS(fitsfiles))[1])); x=x[:,end:-1:1]; nx = (size(x))[1]; x=vec(x)/sum(x); #read in images

    #nfft_plan = setup_nfft(-uv, nx, pixsize)
    #cvis_model = image_to_cvis_nfft(x, nfft_plan)
    if dft==true
        ft = setup_dft(data, nx, pixsize);
    end

    if nfft==true
        ft = setup_nfft(data, nx, pixsize);
    end
    cvis_model = image_to_cvis_dft(x, ft);
    #setup new file
    f = fits_create_file(outfilename);

    #setup initial table
    copy_oi_header(f,oifits[1]);

    ioiarraytables=[]
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
        sta_names=read(oifits[ioiarraytables[itable]],"STA_NAME")
        sta_index=read(oifits[ioiarraytables[itable]],"STA_INDEX")
        tel_diams= read(oifits[ioiarraytables[itable]],"DIAMETER")
        staxyz=read(oifits[ioiarraytables[itable]],"STAXYZ");
        station_xyz=transpose(staxyz);
        N=length(read(oifits[ioiarraytables[itable]],"TEL_NAME"))
        oi_array=[telnames,sta_names,sta_index,tel_diams,staxyz];
        copy_oi_array(f,oifits[ioiarraytables[itable]],oi_array);
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
    copy_oi_target(f,oifits[ioiarraytables[length(ioiarraytables)]-1],target_array);
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
