# gather common display tasks
using PyPlot,PyCall
using LaTeXStrings
PyDict(pyimport("matplotlib")."rcParams")["font.family"]=["serif"]
PyDict(pyimport("matplotlib")."rcParams")["xtick.major.size"]=[6]
PyDict(pyimport("matplotlib")."rcParams")["ytick.major.size"]=[6]
PyDict(pyimport("matplotlib")."rcParams")["xtick.minor.size"]=[6]
PyDict(pyimport("matplotlib")."rcParams")["ytick.minor.size"]=[6]
PyDict(pyimport("matplotlib")."rcParams")["xtick.major.width"]=[1]
PyDict(pyimport("matplotlib")."rcParams")["ytick.major.width"]=[1]
PyDict(pyimport("matplotlib")."rcParams")["xtick.minor.width"]=[1]
PyDict(pyimport("matplotlib")."rcParams")["ytick.minor.width"]=[1]
PyDict(pyimport("matplotlib")."rcParams")["lines.markeredgewidth"]=[1]
PyDict(pyimport("matplotlib")."rcParams")["legend.numpoints"]=[1]
PyDict(pyimport("matplotlib")."rcParams")["legend.handletextpad"]=[0.3]

global colors=["black", "gold","chartreuse","blue","red", "pink","lightgray","darkorange","darkgreen","aqua",
"fuchsia","saddlebrown","dimgray","darkslateblue","violet","indigo","blue","dodgerblue",
"sienna","olive","purple","darkorchid","tomato","darkturquoise","steelblue","seagreen","darkgoldenrod","darkseagreen"]

global markers=["o","s","v","P","*","x","^","D","p",1,"<","H","X","4",4,"_","1",6,"8","d",9]


#xclicks=Array{Float64,1}(undef,1)
#yclicks=Array{Float64,1}(undef,1)
left_click=1
double_click=true
right_click=3
middle_click=2
delete_x=[]
delete_y=[]
delete_err=[]
filename=""
#allow removing of numpoints
function onclickv2(event)
    clicktype=event.button
    xdat=deepcopy(data.v2_baseline)./1e6
    ydat=deepcopy(data.v2)
    errdat=deepcopy(data.v2_err)
    if clicktype == left_click
        if event.dblclick == double_click
            ax=gca()
            ymin,ymax=ax.get_ylim()
            xmin,xmax=ax.get_xlim()
            normfactor=abs(xmax-xmin)/abs(ymax-ymin)
            xclick=event.xdata
            yclick=event.ydata
            diffx=abs.(xdat.-xclick)
            diffy=abs.(ydat.-yclick)
            diffr=sqrt.((xdat.-xclick).^2+((ydat.-yclick).*normfactor).^2)
            indx_remove=indmin(diffr)
            closestx=xdat[indx_remove]
            closesty=ydat[indx_remove]
            closesterr=errdat[indx_remove]
            println(closestx," ",xclick," ",closesty," ",yclick)
            errorbar(closestx,closesty,yerr=closesterr,fmt="o", markersize=3,color="Red")
            push!(delete_x,closestx)
            push!(delete_y,closesty)
            push!(delete_err,closesterr)
            println(delete_x)
        end
    elseif clicktype == right_click
        println("Cancelling!")
        clf()
        errorbar(xdat,ydat,yerr=errdat,fmt="o", markersize=3,color="Black")
    elseif clicktype == middle_click
        filter!(a->a∉delete_x,xdat)
        filter!(a->a∉delete_y,ydat)
        filter!(a->a∉delete_err,errdat)
        clf()
        errorbar(xdat,ydat,yerr=errdat,fmt="o", markersize=3,color="Black")
        # save

    end
end

function onclickidentify(event)
    clicktype=event.button
    if clicktype == left_click    #left click to id
            xclick=event.xdata
            yclick=event.ydata
            ax=gca()
            ymin,ymax=ax.get_ylim()
            xmin,xmax=ax.get_xlim()
            normfactor=abs(xmax-xmin)/abs(ymax-ymin)
            xdat=v2base./1e6
            ydat=v2value
            errdat=v2err
            #diffx=abs.(xdat.-xclicks[length(xclicks)])
            #diffy=abs.(ydat.-yclicks[length(xclicks)])
            diffr=sqrt.((xdat.-xclick).^2+((ydat.-yclick).*normfactor).^2)
            indx_id=argmin(diffr)
            closestx=xdat[indx_id]
            closesty=ydat[indx_id]
            closesterr=errdat[indx_id]
            clickbaseval=clickbase[:,indx_id].+1
            clicksec,clickmin,clickhour,clickday,clickmonth,clickyear=mjd_to_utc(clickmjd[indx_id])
            printstyled("Frequency: ",closestx," Squared Vis Amp (V2): ",closesty," V2 Error: ",closesterr,"\n",color=:blue)
            printstyled("Baseline: ",clickname[clickbaseval[1]],"-",clickname[clickbaseval[2]],"\n",color=:red)
            if length(clickfile)!=1
                printstyled("Filename: ",clickfile[indx_id],"\n",color=:green)
            else
                printstyled("Filename: ",clickfile[1],"\n",color=:green)
            end
            printstyled("MJD: ",clickmjd[indx_id],"\n",color=:black)
            printstyled("Year: ", clickyear," Month: ",clickmonth," Day: ",clickday,"\n",color=:blue)
            printstyled("Hour: ",clickhour, " Min: ",clickmin, " Sec: ", clicksec,"\n",color=:red )
    #elseif clicktype == right_click
    #    fig.canvas.mpl_disconnect(cid)
    end
end


# Overloaded uvplot functions
function uvplot(uv::Array{Float64,2};filename="")
    u = -uv[1,:]/1e6
    v = uv[2,:]/1e6
    fig = figure("UV plot",figsize=(8,8),facecolor="White")
    clf();
    ax = gca()
    markeredgewidth=0.1
    ax.locator_params(axis ="y", nbins=20)
    ax.locator_params(axis ="x", nbins=20)
    axis("equal")
    scatter(u, v,alpha=1.0, s = 12.0,color="Black")
    scatter(-u, -v,alpha=1.0, s = 12.0, color="Black")
    title("UV coverage")
    xlabel(L"U (M$\lambda$)")
    ylabel(L"V (M$\lambda$)")
    ax.grid(true,which="both",color="LightGrey");
    tight_layout();
    if filename !=""
        savefig(filename)
    end
end

function uvplot(data::OIdata;fancy=false,filename="")
    uvplot(data.uv,data.nv2,data.tel_name,data.v2_sta_index,data.v2_lam,fancy=fancy,filename=filename)
end

function uvplot(uv::Array{Float64,2},nv2::Int64,tel_name::Array{String,1},v2_sta_index::Array{Int64,2},v2_lam::Array{Float64,1};fancy=false,filename="")
    u = -uv[1,:]/1e6
    v = uv[2,:]/1e6
    fig = figure("UV plot",figsize=(8,8),facecolor="White")
    clf();
    ax = gca()
    markeredgewidth=0.1
    ax.locator_params(axis ="y", nbins=20)
    ax.locator_params(axis ="x", nbins=20)
    axis("equal")
    if fancy == true
        baseline_list=get_baseline_list_v2(nv2,tel_name,v2_sta_index)
        baseline=sort(unique(baseline_list))
        for i=1:length(baseline)
            loc=findall(baseline_list->baseline_list==baseline[i],baseline_list)
            #scatter(uv[1,loc[1]:loc[length(loc)]].*v2_lam[loc[1]:loc[length(loc)]], uv[2,loc[1]:loc[length(loc)]].*v2_lam[loc[1]:loc[length(loc)]],alpha=0.5, color=colors[i],label=baseline)
            #scatter(-uv[1,loc[1]:loc[length(loc)]].*v2_lam[loc[1]:loc[length(loc)]], -uv[2,loc[1]:loc[length(loc)]].*v2_lam[loc[1]:loc[length(loc)]],alpha=0.5, color=colors[i])
            scatter( u[loc],  v[loc], alpha=1.0, s=12.0, color=colors[i],label=baseline[i])
            scatter(-u[loc], -v[loc], alpha=1.0, s=12.0, color=colors[i])
        end
        ax.legend(fontsize=8, fancybox=true, shadow=true, ncol=3,loc="best")
    else
        scatter(u, v,alpha=1.0, s = 12.0,color="Black")
        scatter(-u, -v,alpha=1.0, s = 12.0, color="Black")
    end
    title("UV coverage")
    xlabel(L"U (M$\lambda$)")
    ylabel(L"V (M$\lambda$)")
    ax.grid(true,which="both",color="LightGrey");
    tight_layout();
    if filename !=""
        savefig(filename)
    end
end


function v2plot_modelvsdata(baseline_v2::Array{Float64,1},v2_data::Array{Float64,1},v2_data_err::Array{Float64,1}, v2_model::Array{Float64,1}; logplot = false) #plots V2 data vs v2 model
    fig = figure("V2 plot - Model vs Data",figsize=(8,8),facecolor="White")
    clf();
    subplot(211)
    ax = gca();
    if logplot==true
        ax.set_yscale("log")
    end
    errorbar(baseline_v2/1e6,v2_data,yerr=v2_data_err,fmt="o", markersize=2,color="Black")
    plot(baseline_v2/1e6, v2_model, color="Red", linestyle="none", marker="o", markersize=3)
    title("Squared Visbility Amplitudes - Model vs data plot")
    #xlabel(L"Baseline (M$\lambda$)")
    ylabel("Squared Visibility Amplitudes")
    ax.grid(true);
    subplot(212)
    plot(baseline_v2/1e6, (v2_model - v2_data)./v2_data_err,color="Black", linestyle="none", marker="o", markersize=3)
    xlabel(L"Baseline (M$\lambda$)")
    ylabel("Residuals (number of sigma)")
    ax = gca();
    ax.grid(true,which="both",color="LightGrey")
    tight_layout()
end



# This draws a continuous line based on the analytic function
function v2plot_modelvsfunc(data::OIdata, visfunc, params; drawpoints = false, yrange=[], drawfunc = true, logplot = false) #plots V2 data vs v2 model
# Compute model points (discrete)
baseline_v2 = data.v2_baseline;
v2_data = data.v2;
v2_data_err = data.v2_err;
cvis_model = visfunc(params,sqrt.(data.uv[1,:].^2+data.uv[2,:].^2));
v2_model = cvis_to_v2(cvis_model, data.indx_v2); # model points
# Compute model curve (continous)
r = sqrt.(data.uv[1,data.indx_v2].^2+data.uv[2,data.indx_v2].^2)
xrange = range(minimum(r),maximum(r),step=(maximum(r)-minimum(r))/1000);
cvis_func = visfunc(params,xrange);
v2_func = abs2.(cvis_func);


fig = figure("V2 plot - Model vs Data",figsize=(8,8),facecolor="White")
clf();
subplot(211)
ax = gca();
if logplot==true
ax.set_yscale("log")
end

if yrange !=[]
ylim((yrange[1], yrange[2]))
end

errorbar(baseline_v2/1e6,v2_data,yerr=v2_data_err,fmt="o", markersize=2,color="Black")
if drawpoints == true
plot(baseline_v2/1e6, v2_model, color="Red", linestyle="none", marker="o", markersize=3)
end

if drawfunc == true
plot(xrange/1e6, v2_func, color="Red", linestyle="-", markersize=3)
end

    title("Squared Visbility Amplitudes - Model vs data plot")
    #xlabel(L"Baseline (M$\lambda$)")
    ylabel("Squared Visibility Amplitudes")
    ax.grid(true,which="both",color="LightGrey")
    subplot(212)
    plot(baseline_v2/1e6, (v2_model - v2_data)./v2_data_err,color="Black", linestyle="none", marker="o", markersize=3)
    xlabel(L"Baseline (M$\lambda$)")
    ylabel("Residuals (number of sigma)")
    ax = gca();
    ax.grid(true,which="both",color="LightGrey")
    tight_layout()
end

function v2plot(data::OIdata;logplot=false,remove=false,idpoint=false,clean=false,fancy=false,markopt=false,ledgendcount=1)
        if idpoint==true
            global v2base=data.v2_baseline
            global v2value=data.v2
            global v2err=data.v2_err
            global clickbase=data.v2_sta_index
            global clickname=data.tel_name
            global clickmjd=data.v2_mjd
            global clickfile=[]
            push!(clickfile,data.filename)
            v2plot(data.v2_baseline,data.v2,data.v2_err,data.nv2,data.tel_name,data.v2_sta_index,logplot=logplot,remove=remove,clean=clean,idpoint=idpoint,fancy=fancy,markopt=markopt)
        else
            v2plot(data.v2_baseline,data.v2,data.v2_err,data.nv2,data.tel_name,data.v2_sta_index,logplot=logplot,remove=remove,clean=clean,idpoint=idpoint,fancy=fancy,markopt=markopt)
        end
end

function v2plot_multifile(data::Array{OIdata,1}; logplot = false, remove = false,idpoint=false,clean=false,filename="")
    global v2base=[]
    global v2value=[]
    global v2err=[]
    global clickbase=Array{Int64,2}
    global clickname=[]
    global clickmjd=[]
    global clickfile=[]
    axiscount=0
    testaxis=0
    fig = figure("V2 data",figsize=(10,5),facecolor="White");
    if clean == true
        clf();
    end
    ax = gca();
    if logplot==true
        ax.set_yscale("log")
    end
    if remove == true # No support for removing points across multimple nights just yet, although this could be VERY useful....
        fig.canvas.mpl_connect("button_press_event",onclickv2)
    end
    for i=1:length(data) #plot each night
        nv2=data[i].nv2
        tel_name=data[i].tel_name
        v2_sta_index=data[i].v2_sta_index
        baseline_v2=data[i].v2_baseline
        v2_data=data[i].v2
        v2_data_err=data[i].v2_err
        if idpoint == true
            v2base=vcat(v2base,baseline_v2)
            v2value=vcat(v2value,v2_data)
            v2err=vcat(v2err,v2_data_err)
            if i==1
                clickbase=cat(data[i].v2_sta_index,dims=2)
            else
                clickbase=cat(clickbase,data[i].v2_sta_index,dims=2)
            end
            clickname=vcat(clickname,data[i].tel_name)
            clickmjd=vcat(clickmjd,data[i].v2_mjd)
            basearray=Array{String,1}(undef,data[i].nv2)
            basearray[1:data[i].nv2].=data[i].filename
            if i==1
                clickfile=cat(basearray,dims=1)
            else
                clickfile=cat(clickfile,basearray,dims=1)
            end
        end
        baseline_list=get_baseline_list_v2(nv2,tel_name,v2_sta_index)
        for j=1:length(unique(baseline_list))
            baseline=unique(baseline_list)[j]
            loc=findall(baseline_list->baseline_list==baseline,baseline_list)
            errorbar(baseline_v2[loc]/1e6,v2_data[loc],yerr=v2_data_err[loc],fmt="o",marker=markers[j], markeredgecolor=colors[i],color=colors[i],markersize=3,ecolor="Gainsboro",elinewidth=1.0,label=baseline)
            if axiscount==0
                if (length(unique(baseline_list)))==15
                    ax.legend(bbox_to_anchor=[0.925,1.0],loc=2,borderaxespad=0)
                    testaxis=1
                end
            end
        end
        if testaxis ==1
            axiscount=1
        end
    end
    title("Squared Visibility Amplitude Data")
    xlabel(L"Baseline (M$\lambda$)")
    ylabel("Squared Visibility Amplitudes")
    ax.grid(true,which="both",color="LightGrey")
    tight_layout()
    if idpoint==true
        cid=fig.canvas.mpl_connect("button_press_event",onclickidentify)
    end
    if filename !=""
        savefig(filename)
    end
end



function v2plot(baseline_v2::Array{Float64,1},v2_data::Array{Float64,1},v2_data_err::Array{Float64,1}, nv2::Int64,tel_name::Array{String,1},v2_sta_index::Array{Int64,2}; logplot = false, remove = false,idpoint=false,clean=true,color="Black",fancy=false,markopt=false,ledgendcount=1) # plots v2 data only
    fig = figure("V2 data",figsize=(10,5),facecolor="White");
    if clean == true
        clf();
    end
    ax = gca();
    if logplot==true
        ax.set_yscale("log")
    end
    if remove == true
        fig.canvas.mpl_connect("button_press_event",onclickv2)
    end
    if fancy == true
        baseline_list=get_baseline_list_v2(nv2,tel_name,v2_sta_index)
        baseline=sort(unique(baseline_list))
        for i=1:length(baseline)
            loc=findall(baseline_list->baseline_list==baseline[i],baseline_list)
            if markopt == false
                errorbar(baseline_v2[loc]/1e6,v2_data[loc],yerr=v2_data_err[loc],fmt="o", markeredgecolor=colors[i],markersize=3,ecolor="Gainsboro",color=colors[i],elinewidth=1.0,label=baseline[i])
            else
                errorbar(baseline_v2[loc]/1e6,v2_data[loc],yerr=v2_data_err[loc],fmt="o",marker=markers[i], markeredgecolor=color,markersize=3,ecolor="Gainsboro",color=color,elinewidth=1.0,label=baseline[i])
            end
        end
        ax.legend(fontsize=8, fancybox=true, shadow=true, ncol=3,loc="best")
    else
        errorbar(baseline_v2/1e6,v2_data,yerr=v2_data_err,fmt="o", markersize=3,color=color,ecolor="Gainsboro",elinewidth=1.0)
    end
    title("Squared Visibility Amplitude Data")
    xlabel(L"Baseline (M$\lambda$)")
    ylabel("Squared Visibility Amplitudes")
    ax.grid(true,which="both",color="LightGrey")
    tight_layout()
    if idpoint==true
        cid=fig.canvas.mpl_connect("button_press_event",onclickidentify)
    end
end


function t3phiplot(data::OIdata;fancy=false,filename="")
    t3phiplot(data.t3_maxbaseline,data.t3phi,data.t3phi_err,data.nt3phi,data.tel_name,data.t3_sta_index;fancy=fancy,filename=filename);
end

function t3phiplot(baseline_t3,t3phi_data,t3phi_data_err,nt3phi,tel_name,t3_sta_index;color="Black",fancy=false,filename="") # plots v2 data only
  fig = figure("Closure phase data",figsize=(10,5),facecolor="White");
  clf();
  ax=gca();
  if fancy == true
      baseline_list=get_baseline_list_t3phi(nt3phi,tel_name,t3_sta_index)
baseline = unique(baseline_list)
      for i=1:length(baseline)
          loc=findall(baseline_list->baseline_list==baseline[i],baseline_list)
          errorbar(baseline_t3[loc]/1e6,t3phi_data[loc],yerr=t3phi_data_err[loc],fmt="o",markeredgecolor=colors[i],color=colors[i], markersize=3,ecolor="Gainsboro",elinewidth=1.0,label=baseline[i])
      end
      ax.legend(fontsize=8, fancybox=true, shadow=true, ncol=3,loc="best")
  else
      errorbar(baseline_t3/1e6,t3phi_data,yerr=t3phi_data_err,fmt="o", markersize=3,color=color, ecolor="Gainsboro",elinewidth=1.0)
  end
  title("Closure phase data")
  xlabel(L"Maximum Baseline (M$\lambda$)")
  ylabel("Closure phase (degrees)")
  ax.grid(true,which="both",color="LightGrey")
  tight_layout()
  if filename !=""
      savefig(filename)
  end
end


function imdisp(image; colormap = "gist_heat", pixscale = -1.0, tickinterval = 10, use_colorbar = false, beamsize = -1, beamlocation = [.9, .9])
 fig = figure("Image",figsize=(6,6),facecolor="White")
 nx=ny=-1;
 pixmode = false;
 if pixscale == -1
     pixmode = true;
     pixscale = 1
 end
 img = []
 if ndims(image) ==1
   ny=nx=Int64(sqrt(length(image)))
   img = imshow(rotl90(reshape(image,nx,nx))/maximum(image), ColorMap(colormap), interpolation="none", extent=[-0.5*nx*pixscale,0.5*nx*pixscale,-0.5*ny*pixscale,0.5*ny*pixscale]); # uses Monnier's orientation
 else
   nx,ny = size(image);
   img = imshow(rotl90(image)/maximum(image), ColorMap(colormap), interpolation="none", extent=[-0.5*nx*pixscale,0.5*nx*pixscale,-0.5*ny*pixscale,0.5*ny*pixscale]); # uses Monnier's orientation
 end
 if pixmode == false
 xlabel("RA (mas)")
 ylabel("DEC (mas)")
end

 ax = gca()
 ax.set_aspect("equal")
 mx = matplotlib.ticker.MultipleLocator(tickinterval) # Define interval of minor ticks
 ax.xaxis.set_minor_locator(mx) # Set interval of minor ticks
 ax.yaxis.set_minor_locator(mx) # Set interval of minor ticks
 ax.xaxis.set_tick_params(which="major",length=10,width=2)
 ax.xaxis.set_tick_params(which="minor",length=5,width=1)
 ax.yaxis.set_tick_params(which="major",length=10,width=2)
 ax.yaxis.set_tick_params(which="minor",length=5,width=1)


 if use_colorbar == true
     divider = pyimport("mpl_toolkits.axes_grid1").make_axes_locatable(ax)
     cax = divider.append_axes("right", size="5%", pad=0.2)
     colorbar(img, cax=cax, ax=ax)
 end

  if beamsize > 0
   c = matplotlib.patches.Circle((0.5*nx*pixscale*beamlocation[1],-0.5*ny*pixscale*beamlocation[2]),beamsize,fc="white",ec="white",linewidth=.5)
   ax.add_artist(c)
  end
 tight_layout()
end

#TODO: work for rectangular
function imdisp_polychromatic(image_vector::Union{Array{Float64,1}, Array{Float64,2},Array{Float64,3}}; imtitle="Polychromatic image", nwavs = 1, colormap = "gist_heat", pixscale = -1.0, tickinterval = 10, use_colorbar = false, beamsize = -1, beamlocation = [.9, .9])
if typeof(image_vector)==Array{Float64,2}
    nwavs = size(image_vector,2)
elseif typeof(image_vector)==Array{Float64,3}
    nwavs = size(image_vector,3)
end

fig = figure(imtitle,figsize=(nwavs*10,4),facecolor="White")
clf();
images_all =reshape(image_vector, (div(length(vec(image_vector)),nwavs), nwavs))
  for i=1:nwavs
    plotnum = 100+nwavs*10+i
    subplot(plotnum)
    title("Wave $i")
    image = images_all[:,i]
    nx=ny=-1;
    pixmode = false;
    if pixscale == -1
      pixmode = true;
      pixscale = 1
    end
    ny=nx=Int64.(sqrt(length(image)))
    img = imshow(rotl90(reshape(image,nx,nx))/maximum(image), ColorMap(colormap), interpolation="none", extent=[-0.5*nx*pixscale,0.5*nx*pixscale,-0.5*ny*pixscale,0.5*ny*pixscale]); # uses Monnier's orientation
    if pixmode == false
        xlabel("RA (mas)")
        ylabel("DEC (mas)")
    end

    ax = gca()
    ax.set_aspect("equal")
    mx = matplotlib.ticker.MultipleLocator(tickinterval) # Define interval of minor ticks
    ax.xaxis.set_minor_locator(mx) # Set interval of minor ticks
    ax.yaxis.set_minor_locator(mx) # Set interval of minor ticks
    ax.xaxis.set_tick_params(which="major",length=10,width=2)
    ax.xaxis.set_tick_params(which="minor",length=5,width=1)
    ax.yaxis.set_tick_params(which="major",length=10,width=2)
    ax.yaxis.set_tick_params(which="minor",length=5,width=1)

    if use_colorbar == true
     divider = pyimport("mpl_toolkits.axes_grid1").make_axes_locatable(ax)
     cax = divider.append_axes("right", size="5%", pad=0.2)
     colorbar(img, cax=cax, ax=ax)
    end

   if beamsize > 0
    c = matplotlib.patches.Circle((0.5*nx*pixscale*beamlocation[1],-0.5*ny*pixscale*beamlocation[2]),beamsize,fc="white",ec="white",linewidth=.5)
    ax.add_artist(c)
   end
  tight_layout()
  end
end

# TODO: rework for julia 1.0+
function imdisp_temporal(image_vector, nepochs; colormap = "gist_heat", name="Time-variable images",pixscale = -1.0, tickinterval = 10, use_colorbar = false, beamsize = -1, beamlocation = [.9, .9])
  fig = figure(name,figsize=(nepochs*10,6+round(nepochs/3)),facecolor="White")
  images_all =reshape(image_vector, (div(length(image_vector),nepochs), nepochs))
  cols=6
  rows=div(nepochs,cols)+1
  for i=1:nepochs
    #plotnum = 100*(div(nepochs,9)+1)
    fig.add_subplot(rows,cols,i)
    #subplot()
    title("Epoch $i")
    image = images_all[:,i]
    nx=ny=-1;
    pixmode = false;
    if pixscale == -1
      pixmode = true;
      pixscale = 1
  end
  img = []
  if ndims(image) ==1
    ny=nx=Int64(sqrt(length(image)))
    img = imshow(rotl90(reshape(image,nx,nx)), ColorMap(colormap), interpolation="none", extent=[-0.5*nx*pixscale,0.5*nx*pixscale,-0.5*ny*pixscale,0.5*ny*pixscale]); # uses Monnier's orientation
  else
    nx,ny = size(image);
    img = imshow(rotl90(image), ColorMap(colormap), interpolation="none", extent=[-0.5*nx*pixscale,0.5*nx*pixscale,-0.5*ny*pixscale,0.5*ny*pixscale]); # uses Monnier's orientation
  end
  if pixmode == false
  xlabel("RA (mas)")
  ylabel("DEC (mas)")
 end

  ax = gca()
  ax.set_aspect("equal")
  mx = matplotlib.ticker.MultipleLocator(tickinterval) # Define interval of minor ticks
  ax.xaxis.set_minor_locator(mx) # Set interval of minor ticks
  ax.yaxis.set_minor_locator(mx) # Set interval of minor ticks
  ax.xaxis.set_tick_params(which="major",length=10,width=2)
  ax.xaxis.set_tick_params(which="minor",length=5,width=1)
  ax.yaxis.set_tick_params(which="major",length=10,width=2)
  ax.yaxis.set_tick_params(which="minor",length=5,width=1)

   if use_colorbar == true
       divider = pyimport("mpl_toolkits.axes_grid1").make_axes_locatable(ax)
       cax = divider.append_axes("right", size="5%", pad=0.2)
       colorbar(img, cax=cax, ax=ax)
   end


   if beamsize > 0
    c = matplotlib.patches.Circle((0.5*nx*pixscale*beamlocation[1],-0.5*ny*pixscale*beamlocation[2]),beamsize,fc="white",ec="white",linewidth=.5)
    ax.add_artist(c)
   end
  tight_layout()
  end
end

function get_baseline_list_v2(nv2,tel_names,v2_stations)
    baseline_list=Array{String}(undef,nv2)
    for i=1:nv2
        baseline_list[i]=string(tel_names[v2_stations[1,i]],"-",tel_names[v2_stations[2,i]])
    end
    return baseline_list
end

function get_baseline_list_t3phi(nt3phi,tel_names,t3_stations)
    baseline_list=Array{String}(undef,nt3phi)
    for i=1:nt3phi
        baseline_list[i]=string(tel_names[t3_stations[1,i]],"-",tel_names[t3_stations[2,i]],"-",tel_names[t3_stations[3,i]])
    end
    return baseline_list
end
