using OITOOLS
using FITSIO

function make_disk(npix,scale,d_mean,majorminor,pa)
    minor=2*d_mean/(1.0+majorminor)
    major=minor*majorminor
    x=repeat(collect(-(npix-1)*scale/2:scale:(npix-1)*scale/2),1,npix)'
    y=x'
    θ=pa*π/180
    x=(x.*cos.(θ).-y.*sin.(θ))'
    y=(x.*sin(θ).+y.*cos(θ))'
    x=x.*d_mean/minor
    y=(y.*d_mean/major)
    r=sqrt.(x.^2 .+ y.^2)
    r[LinearIndices((r.<d_mean/2))[findall((r.<d_mean/2))]].=1
    r[LinearIndices((r.!=1))[findall((r.!=1))]].=0
    return r
end

function place_spot(disk,scale,x_d_spot,y_d_spot,dist_from_edge_mas,pa_loc,pa_spot,ratio)
    #make spot_params
    x_rad=x_d_spot/(2.0)
    y_rad=y_d_spot/(2.0)
    d_mean=(x_d_spot+y_d_spot)/2.0
    bounds=[x_d_spot,y_d_spot]
    npix=Int(round(maximum(bounds)/scale))
    major=maximum(bounds)
    minor=minimum(bounds)
    x=repeat(collect(-(npix-1)*scale/2:scale:(npix-1)*scale/2),1,npix)'
    y=x'
    θ=pa_spot*π/180
    x=(x.*cos.(θ).-y.*sin.(θ))'
    y=(x.*sin(θ).+y.*cos(θ))'
    x=x.*d_mean/minor
    y=(y.*d_mean/major)
    r=sqrt.(x.^2 .+ y.^2)
    r[LinearIndices((r.>d_mean/2))[findall((r.>d_mean/2))]].=0
    r[(div(size(x)[1],2)+1),(div(size(x)[1],2)+1)]=1
    #r[LinearIndices((r.<d_mean/2))[findall((r.<d_mean/2))]].=1
    #r[LinearIndices((r.!=1))[findall((r.!=1))]].=0
    r[LinearIndices((r.!=0))[findall((r.!=0))]].=1
    spot=r*ratio
    #place spot
    delta_y = (dist_from_edge_mas-d_mean/2)*sin((pa_loc+90)*π/180.)
    delta_x = -(dist_from_edge_mas-d_mean/2)*cos((pa_loc+90)*π/180.)
    location_x=Int(round(size(disk)[1]/2+delta_x/scale))
    location_y=Int(round(size(disk)[1]/2+delta_y/scale))
    start_spot_x=Int(round(location_x-npix/2,RoundDown))
    start_spot_y=Int(round(location_y-npix/2,RoundDown))
    end_spot_x=Int(round(location_x+npix/2,RoundDown))
    end_spot_y=Int(round(location_y+npix/2,RoundDown))
    #apply only where there is not a spot and it is in the star.
    disk[start_spot_x:end_spot_x-1,start_spot_y:end_spot_y-1].=disk[start_spot_x:end_spot_x-1,start_spot_y:end_spot_y-1].-spot.*(disk[start_spot_x:end_spot_x-1,start_spot_y:end_spot_y-1].==1)
    return disk
end

function make_spot_move(star_params,spot_params,ra,longitude,error_array,config_files,simulation_params,hours,minutes)
    scale=Float64(star_params[2])
    npix=Int(star_params[1])
    facility_out=read_facility_file(config_files[1],facility_info)
    observatory_out=read_obs_file(config_files[2],obsv_info)
    combiner_out=read_comb_file(config_files[3],combiner_info)
    wave_out=read_wave_file(config_files[4],wave_info)
    errors=define_errors(error_struct,error_array[1],error_array[2],error_array[3],error_array[4],error_array[5],error_array[6])

    #make simulated observations
    speed=simulation_params[1] #mas/day
    span=simulation_params[2] #days
    months=span/30
    daymonth=simulation_params[3]
    steps=months*daymonth
    start_month=simulation_params[4]
    start_day=simulation_params[5]
    year=simulation_params[6]
    hours=simulation_params[7]
    minutes=simulation_params[8]
    nspots=size(spot_params)[1]
    filelist=Array{String}(undef,Int(daymonth*months),2)
    images=Array{Float64}(undef,npix,npix,Int(daymonth*months))
    #disk=make_disk(star_params[1],star_params[2],star_params[3],star_params[4],star_params[5])
    for i=1:months
        for j=1:daymonth
            entry=1
            dates=Array{Float64,2}(undef,Int(size(hours)[1]*size(minutes)[1]),6)
            for k=1:size(hours)[1]
                for l=1:size(minutes)[1]
                    dates[Int(entry),:]=[year start_month+(i-1) start_day+(j-1) hours[Int(i),Int(k)] minutes[Int(l)] 0 ]
                    #println(year," ",start_month+(i-1)," ", start_day+(j-1), " ", hours[Int(i),Int(k)], " " minutes[Int(l)], " ", 0)
                    entry+=1
                end
            end
            number=Int((j+(i-1)*daymonth))
            println(number)
            filename="./data/starspot"*string(number)
            out_file="!"*filename*".oifits"
            lsts,hour_angles=hour_angle_calc(dates,longitude,ra)
            disk=make_disk(npix,scale,star_params[3],star_params[4],star_params[5])
            for m=1:nspots
                disk=place_spot(disk,scale,spot_params[m,1],spot_params[m,2],spot_params[m,3]*speed*(j+(i-1)*daymonth),spot_params[m,4],spot_params[m,5],spot_params[m,6])
            end
            image_file=filename*".fits"
            f=FITS(image_file,"w")
            write(f,disk)
            close(f)
            imdisp(disk)
            filelist[Int(j+(i-1)*daymonth),1]=filename*".oifits"
            filelist[Int(j+(i-1)*daymonth),2]=filename*".fits"
            simulate_ha(facility_out,observatory_out,combiner_out,wave_out,hour_angles,image_file,scale,errors,out_file)
            images[:,:,Int(j+(i-1)*daymonth)]=disk
        end
    end
    return filelist,images
end
