#
# TO DO: plot_t3phivsdata, plot_t3phivsmodel, everything should be able to take as input a 1D array of OIdata
#

# gather common display tasks
using PyPlot,PyCall, LaTeXStrings, Statistics

function set_oiplot_defaults()
    PyDict(pyimport("matplotlib")."rcParams")["font.family"]=["serif"]
    PyDict(pyimport("matplotlib")."rcParams")["font.size"]=[14]
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
    #PyDict(pyimport("matplotlib")."rcParams")["agg.path.chunksize"]=[10000]
end

const global oiplot_colors=["black", "gold","chartreuse","blue","red", "pink","lightgray","darkorange","darkgreen","aqua",
"fuchsia","saddlebrown","dimgray","darkslateblue","violet","indigo","blue","dodgerblue",
"sienna","olive","purple","darkorchid","tomato","darkturquoise","steelblue","seagreen","darkgoldenrod","darkseagreen","salmon","slategray","lime","coral","maroon","mistyrose","sandybrown","tan","olivedrab"]

const global oiplot_markers=["o","s","v","P","*","x","^","D","p",1,"<","H","X","4",4,"_","1",6,"8","d",9]


# #xclicks=Array{Float64,1}(undef,1)
# #yclicks=Array{Float64,1}(undef,1)
# left_click=1
# double_click=true
# right_click=3
# middle_click=2
# delete_x=[]
# delete_y=[]
# delete_err=[]
# filename=""
# #allow removing of numpoints
# function onclickv2(event)
#     clicktype=event.button
#     xdat=deepcopy(data.v2_baseline)./1e6
#     ydat=deepcopy(data.v2)
#     errdat=deepcopy(data.v2_err)
#     if clicktype == left_click
#         if event.dblclick == double_click
#             ax=gca()
#             ymin,ymax=ax.get_ylim()
#             xmin,xmax=ax.get_xlim()
#             normfactor=abs(xmax-xmin)/abs(ymax-ymin)
#             xclick=event.xdata
#             yclick=event.ydata
#             diffx=abs.(xdat.-xclick)
#             diffy=abs.(ydat.-yclick)
#             diffr=sqrt.((xdat.-xclick).^2+((ydat.-yclick).*normfactor).^2)
#             indx_remove=indmin(diffr)
#             closestx=xdat[indx_remove]
#             closesty=ydat[indx_remove]
#             closesterr=errdat[indx_remove]
#             println(closestx," ",xclick," ",closesty," ",yclick)
#             errorbar(closestx,closesty,yerr=closesterr,fmt="o", markersize=3,color="Red")
#             push!(delete_x,closestx)
#             push!(delete_y,closesty)
#             push!(delete_err,closesterr)
#             println(delete_x)
#         end
#     elseif clicktype == right_click
#         println("Cancelling!")
#         clf()
#         errorbar(xdat,ydat,yerr=errdat,fmt="o", markersize=3,color="Black")
#     elseif clicktype == middle_click
#         filter!(a->a∉delete_x,xdat)
#         filter!(a->a∉delete_y,ydat)
#         filter!(a->a∉delete_err,errdat)
#         clf()
#         errorbar(xdat,ydat,yerr=errdat,fmt="o", markersize=3,color="Black")
#         # save
#     end
# end

# function onclickidentify(event)
#     clicktype=event.button
#     if clicktype == left_click    #left click to id
#         xclick=event.xdata
#         yclick=event.ydata
#         if !isnothing(xclick) & !isnothing(yclick)
#             ax=gca()
#             ymin,ymax=ax.get_ylim()
#             xmin,xmax=ax.get_xlim()
#             normfactor=abs(xmax-xmin)/abs(ymax-ymin)
#             xdat = v2base./1e6
#             ydat = v2value
#             errdat = v2err
#             indx_id = argmin(sqrt.((xdat.-xclick).^2+((ydat.-yclick).*normfactor).^2))
#             clickbaseval=clickbase[:,indx_id]
#             printstyled("----------------------------\n",color=:black);
#             if length(clickfile)!=1
#                 printstyled("Filename: ",clickfile[indx_id],"\n",color=:orange)
#             else
#                 printstyled("Filename: ",clickfile[1],"\n",color=:orange)
#             end
#             printstyled("Radial frequency: ",xdat[indx_id]," V2: ",ydat[indx_id]," V2_err: ",errdat[indx_id],"\n",color=:blue);
#             printstyled("λ: ", clicklam[indx_id], " δλ: ", clickdlam[indx_id], "\n",color=:green);
#             printstyled("Baseline: ",clickname[clickbaseval[1]],"-",clickname[clickbaseval[2]],"\n",color=:red);
#             printstyled("MJD: ",clickmjd[indx_id], " UT: ", Dates.format(mjd_to_utdate(clickmjd[indx_id]),"Y-m-d H:M:S"), "\n",color=:red )
#             #elseif clicktype == right_click
#             #    fig.canvas.mpl_disconnect(cid
#         end
#     end
# end


# Overloaded uvplot functions
function uvplot(uv::Array{Float64,2})
    u = uv[1,:]/1e6
    v = uv[2,:]/1e6
    fig = figure("UV plot",figsize=(8,8),facecolor="White")
    set_oiplot_defaults()
    clf();
    ax = gca()
    markeredgewidth=0.1
    ax.locator_params(axis ="y", nbins=20)
    ax.locator_params(axis ="x", nbins=20)
    ax.set_aspect("equal")
    scatter(u, v,alpha=1.0, s = 12.0,color="Black")
    scatter(-u, -v,alpha=1.0, s = 12.0, color="Black")
    title("UV coverage")
    xlabel(L"U (M$\lambda$)")
    ylabel(L"V (M$\lambda$)")
    ax.grid(true,which="both",color="Grey", linestyle=":");
    tight_layout();
end

function uvplot(data::Union{OIdata,Array{OIdata,1}, Array{OIdata,2}};color::String="baseline",filename="", figsize=(10,10), minuv= -1e99, maxuv= 1e99, square = true, legend_below = true, figtitle = "", windowtitle="", cmap="Spectral_r", flipx = false)
    set_oiplot_defaults()
    if typeof(data)==OIdata
        data = [data]
    end
    if typeof(data)==Array{OIdata,2}
        data = vec(data)
    end
    nuv = sum(data[i].nuv for i=1:length(data))
    mean_mjd = mean(data[i].mean_mjd for i=1:length(data))
    if windowtitle==""
        windowtitle = string("Mean MJD: $(round(mean_mjd*100)/100), nuv: $(nuv)")
    end
    fig = figure(windowtitle,figsize=figsize,facecolor="White")
    set_oiplot_defaults()
    clf();
    ax = gca()
    ax.locator_params(axis ="y", nbins=10)
    ax.locator_params(axis ="x", nbins=10)
    ax.set_aspect(1.0)
    if (color == "baseline" || color =="base" || color =="bases" || color =="baselines") # we need to identify corresponding baselines #TBD --> could be offloaded to readoifits
        baseline_list_v2 = [get_baseline_names(data[n].sta_name,data[n].v2_sta_index) for n=1:length(data)];
        baseline_list_t3 = [get_triple_baselines_names(data[n].sta_name,data[n].t3_sta_index) for n=1:length(data)];
        baseline=sort(unique(vcat(vcat(baseline_list_v2...),vec(hcat(baseline_list_t3...)))))
        if length(baseline)>length(oiplot_colors)
            @warn("I ran out of colors!")
        end
        indx_v2 = [data[n].indx_v2 for n=1:length(data)]
        indx_t3 = [hcat(data[n].indx_t3_1,data[n].indx_t3_2, data[n].indx_t3_3)' for n=1:length(data)]
        for i=1:length(baseline)
            loc =  [vcat(indx_v2[n][findall(baseline_list_v2[n] .== baseline[i])], indx_t3[n][findall(baseline_list_t3[n] .== baseline[i])]) for n=1:length(data)]
            u = vcat([data[n].uv[1,loc[n]] for n=1:length(data)]...)/1e6;
            v = vcat([data[n].uv[2,loc[n]] for n=1:length(data)]...)/1e6;
            scatter( u,  v, alpha=1.0, s=12.0, color=oiplot_colors[i],label=baseline[i]) #TBD: handle case where length(baseline)>length(oiplot_colors)
            scatter(-u, -v, alpha=1.0, s=12.0, color=oiplot_colors[i])
        end
               
        if legend_below == false
            ax.legend(fontsize=10, fancybox=true, shadow=true, ncol=3, loc="upper right")
        else
            ax.legend(fontsize=10, fancybox=true, shadow=true, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.10));
            tight_layout();
        end
    elseif (color == "wavelength" || color == "wav" || color =="wavs" || color =="wavelengths")
        u = vcat([data[n].uv[1,:]/1e6 for n=1:length(data)]...)
        v = vcat([data[n].uv[2,:]/1e6 for n=1:length(data)]...)
        wavcol = vcat([data[n].uv_lam*1e6 for n=1:length(data)]...)
        wavs=sort(unique(wavcol))
        scatter([u;-u], [v;-v], alpha=1.0, s = 12.0, c=[wavcol;wavcol], cmap=cmap)
        cbar = colorbar(ax=ax, aspect=50, orientation="horizontal", label="Wavelength (μm)", pad=0.1, fraction=0.02)
        if length(wavs)==2
            cbar_range = round.(collect(range(ceil(minimum(wavcol)*100)/100, floor(maximum(wavcol)*100)/100,length=2))*100)/100
            cbar.set_ticks(cbar_range)
            cbar.set_ticklabels(cbar_range)
        elseif length(wavs)==1
            cbar_range = wavs
            cbar.set_ticks(cbar_range)
            cbar.set_ticklabels(cbar_range)
        else # >= 9
            nwavs = min(length(wavs), 10) # max 10 divisions
            minrange = ceil(minimum(wavcol)*100)/100
            maxrange = floor(maximum(wavcol)*100)/100
            step = round((maxrange-minrange)/nwavs*100)/100
            cbar_range = collect(range(minrange, maxrange, step=step))
            cbar.set_ticks(cbar_range)
        end
    elseif (color == "mjd" || color == "time")
        u = vcat([data[n].uv[1,:]/1e6 for n=1:length(data)]...)
        v = vcat([data[n].uv[2,:]/1e6 for n=1:length(data)]...)
        mjdcol = vcat([data[n].uv_mjd for n=1:length(data)]...)
        scatter([u;-u], [v;-v], alpha=1.0, s = 12.0, c=[mjdcol;mjdcol], cmap=cmap)
        cbar = colorbar(ax=ax, aspect=50, orientation="horizontal", label="MJD", pad=0.1, fraction=0.02)
        mjds=unique(mjdcol)
        if length(mjds)<5
            cbar.set_ticks(round.(sort(mjds)*100)/100)
            cbar.set_ticklabels(round.(sort(mjds)*100)/100)
        else
            cbar_range = round.(collect(range(ceil(minimum(mjdcol)*100)/100, floor(maximum(mjdcol)*100)/100,length=5))*100)/100
            cbar.set_ticks(cbar_range)
            cbar.set_ticklabels(cbar_range)
        end
    else
        u = vcat([data[n].uv[1,:]/1e6 for n=1:length(data)]...)
        v = vcat([data[n].uv[2,:]/1e6 for n=1:length(data)]...)
        scatter(u, v,alpha=1.0, s = 12.0,color="Black")
        scatter(-u, -v,alpha=1.0, s = 12.0, color="Black")
    end
    if square==true
        if minuv<-1e98
            minuv = minimum([ax.get_xlim()[1],ax.get_ylim()[1]] )
        end
        if maxuv>1e98
            maxuv = maximum([ax.get_xlim()[2],ax.get_ylim()[2]])
        end     
        ax.set_xlim((minuv,maxuv))
        ax.set_ylim((minuv,maxuv))
    end
    title(figtitle)
    xlabel(L"U (M$\lambda$)")
    ylabel(L"V (M$\lambda$)")
    ax.grid(true,which="both",color="Grey",linestyle=":");
    ax.set_aspect("equal")
    if flipx == true
        ax.invert_xaxis()
    end
    if filename !=""
        savefig(filename)
    end
    tight_layout();
    show(block=false)
end

# This draws a continuous line based on the analytic function
function plot_v2_vs_func(data::OIdata, model::OImodel, params; drawpoints = false, yrange=[], drawfunc = true, logplot = false) #plots V2 data vs v2 model
    set_oiplot_defaults()
    # Compute model points (discrete)
    baseline_v2 = data.v2_baseline;
    v2_data = data.v2;
    v2_data_err = data.v2_err;
    update_model(model)
    dispatch_params(params, model);
    cvis_model = model_to_vis(model, data)
    v2_model = cvis_to_v2(cvis_model, data.indx_v2); # model points
    # Compute model curve (pseudo-continous) by setting up cyclindrical uv plane
    r = sqrt.(data.uv[1,data.indx_v2].^2+data.uv[2,data.indx_v2].^2)
    r_range = collect(range(minimum(r),maximum(r),length=100));
    r_proj = collect(range(minimum(r),maximum(r),length=100));
    θ_range = collect(range(-pi,pi,length=100));
    uv = repeat(hcat(vec(r_range'.*cos.(θ_range)), vec(r_range'.*sin.(θ_range)))',1,100)
    λ = repeat(collect(range(minimum(data.uv_lam),maximum(data.uv_lam),step=100)), size(uv,2));
    cvis_func = model_to_vis(model,uv, λ)
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
        plot(r_proj/1e6, v2_func, color="Red", linestyle="-", markersize=3)
    end

    title("Squared Visbility Amplitudes - Model vs data plot")
    #xlabel(L"Baseline (M$\lambda$)")
    ylabel("Squared Visibility Amplitudes")
    ax.grid(true,which="both",color="Grey",linestyle=":")
    subplot(212)
    plot(baseline_v2/1e6, (v2_model - v2_data)./v2_data_err,color="Black", linestyle="none", marker="o", markersize=3)
    xlabel(L"Baseline (M$\lambda$)")
    ylabel("Residuals (number of sigma)")
    ax = gca();
    ax.grid(true,which="both",color="Grey",linestyle=":")
    tight_layout()
    show(block=false)
end

function plot_v2(data::Union{OIdata,Array{OIdata,1},Array{OIdata,2}};figsize=(12, 6), logplot = false, remove = false,idpoint=false,clean=true,color::String="baseline",markopt=false, legend_below=false, figtitle="")
    set_oiplot_defaults()

    if idpoint==true # interactive plot, click to identify point
        global v2base=data.v2_baseline
        global v2value=data.v2
        global v2err=data.v2_err
        global clickbase=data.v2_sta_index
        global clickname=data.sta_name
        global clickmjd=data.v2_mjd
        global clicklam=data.v2_lam
        global clickdlam=data.v2_dlam
        global clickfile=[]
        push!(clickfile,data.filename)
    end
    if typeof(data)==OIdata
        data = [data]
    end
    if typeof(data)==Array{OIdata,2}
        data = vec(data)
    end
    fig = figure(string(figtitle, "V2 data"),figsize=figsize,facecolor="White");
    if clean == true # do not overplot on existing window by default
        clf();
    end
    ax = gca();
    if logplot==true
        ax.set_yscale("log")
    end
    #if remove == true
    #    fig.canvas.mpl_connect("button_press_event",onclickv2)
    #end
    if (color == "baseline" || color =="base"|| color =="bases" || color =="baselines") # we need to identify corresponding baselines #TBD --> could be offloaded to readoifits
        baseline_list_v2 = [get_baseline_names(data[n].sta_name,data[n].v2_sta_index) for n=1:length(data)];
        baseline=sort(unique(vcat(baseline_list_v2...)))
        for i=1:length(baseline)
            loc = [findall(baseline_list_v2[n] .== baseline[i]) for n=1:length(data)]
            baseline_v2 = vcat([data[n].v2_baseline[loc[n]] for n=1:length(data)]...)/1e6
            v2 = vcat([data[n].v2[loc[n]] for n=1:length(data)]...)
            v2_err = vcat([data[n].v2_err[loc[n]] for n=1:length(data)]...)
            if markopt == false
                errorbar(baseline_v2,v2,yerr=v2_err,fmt="o", markeredgecolor=oiplot_colors[i],markersize=3,ecolor="Gainsboro",color=oiplot_colors[i],elinewidth=1.0,label=baseline[i])
            else
                errorbar(baseline_v2,v2,yerr=v2_err,fmt="o",marker=oiplot_markers[i], markeredgecolor="Black",markersize=3,ecolor="Gainsboro",color="Black",elinewidth=1.0,label=baseline[i])
            end
        end
        if legend_below == false
            ax.legend(fontsize=8, fancybox=true, shadow=true, ncol=4,loc="best")
        else
            ax.legend(fontsize=8, fancybox=true, shadow=true, ncol=8,loc="upper center", bbox_to_anchor=(0.5, -0.15))
        end
    elseif (color == "wavelength" || color == "wav" || color =="wavs" || color =="wavelengths")
        wavcol = vcat([data[n].uv_lam[data[n].indx_v2]*1e6 for n=1:length(data)]...)
        baseline_v2 = vcat([data[n].v2_baseline for n=1:length(data)]...)/1e6
        v2 = vcat([data[n].v2 for n=1:length(data)]...);
        v2_err = vcat([data[n].v2_err for n=1:length(data)]...);
        sc = scatter(baseline_v2,v2, c=wavcol, cmap="Spectral_r", alpha=1.0, s=6.0, zorder=100)
        el = errorbar(baseline_v2,v2,yerr=v2_err,fmt="none", marker="none",ecolor="Gainsboro", elinewidth=1.0, zorder=0)
        cbar = colorbar(sc, aspect=50, orientation="horizontal", label="Wavelength (μm)", pad=0.18, fraction=0.05)
        cbar_range = floor.(collect(range(minimum(wavcol), maximum(wavcol), length=11))*100)/100
        cbar.set_ticks(cbar_range)
    else
        baseline_v2 = vcat([data[n].v2_baseline for n=1:length(data)]...)/1e6
        v2 = vcat([data[n].v2 for n=1:length(data)]...);
        v2_err = vcat([data[n].v2_err for n=1:length(data)]...);
        errorbar(baseline_v2,v2,yerr=v2_err,fmt="o", markersize=3,color="Black",ecolor="Gainsboro",elinewidth=1.0);
    end
    title("Squared Visibility Amplitude Data")
    xlabel(L"Baseline (M$\lambda$)")
    ylabel("Squared Visibility Amplitudes")
    ax.grid(true,which="both",color="Grey", linestyle=":")
    tight_layout()
    #if idpoint==true
    #    cid=fig.canvas.mpl_connect("button_press_event",onclickidentify)
    #end
    show(block=false)
end


function plot_t3phi(data::Union{OIdata,Array{OIdata,1},Array{OIdata,2}}; figsize=(12,6), color::String="baseline",markopt=false, legend_below=false, t3base="max")
    set_oiplot_defaults()

    if typeof(data)==OIdata
        data = [data]
    end
    if typeof(data)==Array{OIdata,2}
        data = vec(data)
    end
    fig = figure("Closure phase data",figsize=figsize,facecolor="White");
    clf();
    ax=gca();
    if (color == "baseline" || color =="base"|| color =="bases" || color =="baselines")
        baseline_list_t3 = [get_triplet_names(data[n].sta_name,data[n].t3_sta_index) for n=1:length(data)];
        baseline=sort(unique(vcat(baseline_list_t3...)))
        #indx_t3 = [hcat(data[n].indx_t3_1,data[n].indx_t3_2, data[n].indx_t3_3)' for n=1:length(data)]
        for i=1:length(baseline)
            loc =  [findall(baseline_list_t3[n] .== baseline[i]) for n=1:length(data)]
            if t3base=="max"
                baseline_t3 = vcat([data[n].t3_maxbaseline[loc[n]] for n=1:length(data)]...)/1e6
            elseif t3base=="geom"
                baseline_t3 = vcat([data[n].t3_baseline[loc[n]] for n=1:length(data)]...)/1e6
            end
            t3phi = vcat([data[n].t3phi[loc[n]] for n=1:length(data)]...)
            t3phi_err = vcat([data[n].t3phi_err[loc[n]] for n=1:length(data)]...)
            errorbar(baseline_t3,t3phi,yerr=t3phi_err,fmt="o",markeredgecolor=oiplot_colors[i],color=oiplot_colors[i], markersize=3,ecolor="Gainsboro",elinewidth=1.0,label=baseline[i])
        end
        if legend_below == false
            ax.legend(fontsize=8, fancybox=true, shadow=true, ncol=4,loc="best")
        else
            ax.legend(fontsize=8, fancybox=true, shadow=true, ncol=5,loc="upper center", bbox_to_anchor=(0.5, -0.15))
        end
    elseif  (color == "wavelength" || color == "wav" || color =="wavs" || color =="wavelengths")
        wavcol = vcat([data[n].uv_lam[data[n].indx_t3_1]*1e6 for n=1:length(data)]...)
        if t3base=="max"
            baseline_t3 = vcat([data[n].t3_maxbaseline for n=1:length(data)]...)/1e6
        elseif t3base=="geom"
            baseline_t3 = vcat([data[n].t3_baseline for n=1:length(data)]...)/1e6
        end
        t3phi = vcat([data[n].t3phi for n=1:length(data)]...);
        t3phi_err = vcat([data[n].t3phi_err for n=1:length(data)]...);
        sc = scatter(baseline_t3, t3phi, c=wavcol, cmap="Spectral_r", alpha=1.0, s=6.0, zorder=100)
        el = errorbar(baseline_t3, t3phi,yerr=t3phi_err,fmt="none", marker="none",ecolor="Gainsboro", elinewidth=1.0, zorder=0)
        cbar = colorbar(sc, aspect=50, orientation="horizontal", label="Wavelength (μm)", pad=0.18, fraction=0.05)
        cbar_range = floor.(collect(range(minimum(wavcol), maximum(wavcol), length=11))*100)/100
        cbar.set_ticks(cbar_range)
    else
        if t3base=="max"
            baseline_t3 = vcat([data[n].t3_maxbaseline for n=1:length(data)]...)/1e6
        elseif t3base=="geom"
            baseline_t3 = vcat([data[n].t3_baseline for n=1:length(data)]...)/1e6
        end
        t3phi = vcat([data[n].t3phi for n=1:length(data)]...);
        t3phi_err = vcat([data[n].t3phi_err for n=1:length(data)]...);
        errorbar(baseline_t3,t3phi,yerr=t3phi_err,fmt="o", markersize=3,color="Black", ecolor="Gainsboro",elinewidth=1.0)
    end
    title("Closure phase data")
    if t3base=="max"
        xlabel(L"Maximum Baseline (M$\lambda$)")
    elseif t3base=="geom"
        xlabel(L"Geometric Mean Baseline (M$\lambda$)")
    end
    ylabel("Closure phase (degrees)")
    ax.grid(true,which="both",color="Grey",linestyle=":")
    tight_layout()
    show(block=false)
end

function plot_v2_residuals(x, data::OIdata, ft::Array{NFFT.NFFTPlan{Float64, 2, 1}, 1}; logplot = false, y_range=[], res_range=[])
    plot_v2_residuals(data, image_to_v2(x, data, ft), logplot = logplot, y_range=y_range, res_range=res_range);
end

function plot_t3phi_residuals(x, data::OIdata, ft::Array{NFFT.NFFTPlan{Float64, 2, 1}, 1}; logplot = false, y_range=[], res_range=[])
    plot_t3phi_residuals(data, image_to_t3phi(x, data, ft), logplot = logplot, y_range=y_range, res_range=res_range);
end

function plot_t3amp_residuals(x, data::OIdata, ft::Array{NFFT.NFFTPlan{Float64, 2, 1}, 1}; logplot = false, y_range=[], res_range=[])
    plot_t3amp_residuals(data, image_to_t3amp(x, data, ft), logplot = logplot, y_range=y_range, res_range=res_range);
end


function plot_v2_residuals(data::OIdata, v2_model::Array{Float64,1}; figsize=(12,12), logplot = false, y_range=[], res_range=[]) #plots V2 data vs v2 model
    set_oiplot_defaults()

    fig = figure("V2 plot - Model vs Data",figsize=figsize,facecolor="White");
    fig.subplots(nrows=2, ncols=1, sharex=true)
    ax = subplot(2, 1, 1)
    if logplot==true
        ax.set_yscale("log")
    end
    errorbar(data.v2_baseline/1e6,data.v2,yerr=data.v2_err,fmt="o", markersize=1, color="Grey", ecolor="Grey")
    plot(data.v2_baseline/1e6, v2_model, color="Red", linestyle="none", marker="o", markersize=2)
    title("Squared Visibility Amplitudes - Model vs data plot")
    if y_range != []
        ylim(y_range)
    end
    ylabel("Squared Visibility Amplitudes")
    ax.grid(true,which="both",color="Grey",linestyle=":");
    ax = subplot(2, 1, 2)
    ax.axhline(color="Grey")
    plot(data.v2_baseline/1e6, (v2_model - data.v2)./data.v2_err, color="Grey", linestyle="none", marker="o", markersize=2)
    xlabel(L"Baseline (M$\lambda$)")
    ylabel("Residuals (number of sigma)")
    if res_range !=[]
        ylim(res_range)
    end
    ax.grid(true,which="both",color="Grey",linestyle=":")
    tight_layout()
    show(block=false)
end

function plot_t3phi_residuals(data::OIdata, t3phi_model::Array{Float64,1}; figsize=(12,12), logplot = false, y_range=[],  res_range=[]) #plots V2 data vs v2 model
    set_oiplot_defaults()

    fig = figure("Closure phase plot - Model vs Data",figsize=figsize,facecolor="White");
    fig.subplots(nrows=2, ncols=1, sharex=true)
    ax = subplot(2, 1, 1)
    if logplot==true
        ax.set_yscale("log")
    end
    errorbar(data.t3_maxbaseline/1e6,data.t3phi,yerr=data.t3phi_err,fmt="o", markersize=1, color="Grey", ecolor="Grey")
    plot(data.t3_maxbaseline/1e6, t3phi_model, color="Red", linestyle="none", marker="o", markersize=2)
    title("Closure phases - Model vs data plot")
    if y_range != []
        ylim(y_range)
    end
    ylabel("Closure phases (degrees)")
    ax.grid(true,which="both",color="Grey",linestyle=":");
    ax = subplot(2, 1, 2)
    ax.axhline(color="Grey")
    plot(data.t3_maxbaseline/1e6, mod360(t3phi_model - data.t3phi)./data.t3phi_err,color="Grey", linestyle="none", marker="o", markersize=1)
    if res_range !=[]
        ylim(res_range)
    end
    xlabel(L"Max Baseline (M$\lambda$)")
    ylabel("Residuals (number of sigma)")
    ax.grid(true,which="both",color="Grey",linestyle=":")
    tight_layout()
    show(block=false)
end

function plot_t3amp_residuals(data::OIdata, t3amp_model::Array{Float64,1}; figsize=(12,12), logplot = false, y_range=[],  res_range=[]) #plots V2 data vs v2 model
    set_oiplot_defaults()

    fig = figure("Triple amplitudes plot - Model vs Data",figsize=figsize,facecolor="White");
    fig.subplots(nrows=2, ncols=1, sharex=true)
    ax = subplot(2, 1, 1)
    if logplot==true
        ax.set_yscale("log")
    end
    errorbar(data.t3_maxbaseline/1e6,data.t3amp,yerr=data.t3amp_err,fmt="o", markersize=1, color="Grey", ecolor="Grey")
    plot(data.t3_maxbaseline/1e6, t3amp_model, color="Red", linestyle="none", marker="o", markersize=2)
    title("Triple amplitudes - Model vs data plot")
    if y_range != []
        ylim(y_range)
    end
    ylabel("Triple amplitudes")
    ax.grid(true,which="both",color="Grey",linestyle=":");
    ax = subplot(2, 1, 2)
    ax.axhline(color="Grey")
    plot(data.t3_maxbaseline/1e6, (t3amp_model - data.t3amp)./data.t3amp_err,color="Grey", linestyle="none", marker="o", markersize=2)
    if res_range !=[]
        ylim(res_range)
    end
    xlabel(L"Max Baseline (M$\lambda$)")
    ylabel("Residuals (number of sigma)")
    ax = gca();
    ax.grid(true,which="both",color="Grey",linestyle=":")
    tight_layout()
    show(block=false)
end


function plot_t3amp(data::Union{OIdata,Array{OIdata,1},Array{OIdata,2}}; figsize=(12,6), color::String="baseline",markopt=false, legend_below=false, t3base="max")
    set_oiplot_defaults()

    if typeof(data)==OIdata
        data = [data]
    end
    if typeof(data)==Array{OIdata,2}
        data = vec(data)
    end
    fig = figure("Triple amplitude data",figsize=(12,6),facecolor="White");
    clf();
    ax=gca();
    baseline_t3 = []

    if (color == "baseline" || color =="base"|| color =="bases" || color =="baselines")
        baseline_list_t3 = [get_triplet_names(data[n].sta_name,data[n].t3_sta_index) for n=1:length(data)];
        baseline=sort(unique(vcat(baseline_list_t3...)))
        #indx_t3 = [hcat(data[n].indx_t3_1,data[n].indx_t3_2, data[n].indx_t3_3)' for n=1:length(data)]
        for i=1:length(baseline)
            loc =  [findall(baseline_list_t3[n] .== baseline[i]) for n=1:length(data)]
            if t3base=="max"
                baseline_t3 = vcat([data[n].t3_maxbaseline[loc[n]] for n=1:length(data)]...)/1e6
            elseif t3base=="geom"
                baseline_t3 = vcat([data[n].t3_baseline[loc[n]] for n=1:length(data)]...)/1e6
            end
            t3amp = vcat([data[n].t3amp[loc[n]] for n=1:length(data)]...)
            t3amp_err = vcat([data[n].t3amp_err[loc[n]] for n=1:length(data)]...)
            errorbar(baseline_t3,t3amp,yerr=t3amp_err,fmt="o",markeredgecolor=oiplot_colors[i],color=oiplot_colors[i], markersize=3,ecolor="Gainsboro",elinewidth=1.0,label=baseline[i])
        end
        if legend_below == false
            ax.legend(fontsize=8, fancybox=true, shadow=true, ncol=4,loc="best")
        else
            ax.legend(fontsize=8, fancybox=true, shadow=true, ncol=8,loc="upper center", bbox_to_anchor=(0.5, -0.15))
        end
    elseif  (color == "wavelength" || color == "wav" || color =="wavs" || color =="wavelengths")
        wavcol = vcat([data[n].uv_lam[data[n].indx_t3_1]*1e6 for n=1:length(data)]...)
        
        if t3base=="max"
            baseline_t3 = vcat([data[n].t3_maxbaseline for n=1:length(data)]...)/1e6
        elseif t3base=="geom"
            baseline_t3 = vcat([data[n].t3_baseline for n=1:length(data)]...)/1e6
        end
        t3amp = vcat([data[n].t3amp for n=1:length(data)]...);
        t3amp_err = vcat([data[n].t3amp_err for n=1:length(data)]...);
        sc = scatter(baseline_t3, t3amp, c=wavcol, cmap="Spectral_r", alpha=1.0, s=6.0, zorder=100)
        el = errorbar(baseline_t3, t3amp,yerr=t3amp_err,fmt="none", marker="none",ecolor="Gainsboro", elinewidth=1.0, zorder=0)
        cbar = colorbar(sc, aspect=50, orientation="horizontal", label="Wavelength (μm)", pad=0.18, fraction=0.05)
        cbar_range = floor.(collect(range(minimum(wavcol), maximum(wavcol), length=11))*100)/100
        cbar.set_ticks(cbar_range)
    else
        if t3base=="max"
            baseline_t3 = vcat([data[n].t3_maxbaseline for n=1:length(data)]...)/1e6
        elseif t3base=="geom"
            baseline_t3 = vcat([data[n].t3_baseline for n=1:length(data)]...)/1e6
        end
        t3amp = vcat([data[n].t3amp for n=1:length(data)]...);
        t3amp_err = vcat([data[n].t3amp_err for n=1:length(data)]...);
        errorbar(baseline_t3,t3amp,yerr=t3amp_err,fmt="o", markersize=3,color="Black", ecolor="Gainsboro",elinewidth=1.0)
    end
    title("Triple amplitude data")
    if t3base=="max"
        xlabel(L"Maximum Baseline (M$\lambda$)")
    elseif t3base=="geom"
        xlabel(L"Geometric Mean Baseline (M$\lambda$)")
    end
    ylabel("Triple amplitude")
    ax.grid(true,which="both",color="Grey",linestyle=":")
    tight_layout()
    show(block=false)
end

function plot_flux(data::Union{OIdata,Array{OIdata,1}}; color="Black")
    list_stations = sort(unique(vcat([data[n].flux_sta_index for n=1:length(data)]...)))
    for i=1:length(list_stations)
        for n=1:length(data)
            indx = findall(data[n].flux_sta_index .== list_stations[i])
            wavcol = data[n].flux_lam[indx]*1e6
            mjd = data[n].flux_mjd[indx]
            flux = data[n].flux[indx]
            flux_err = data[n].flux_err[indx]
            el = errorbar(wavcol,flux,yerr=flux_err,fmt="none", marker="none",ecolor="Gainsboro", elinewidth=1.0, zorder=0)
        end
    end
end

function plot_visphi(data::Union{OIdata,Array{OIdata,1}}; color::String="baseline",markopt=false, legend_below=false)
    if typeof(data)==OIdata
        data = [data]
    end
    if typeof(data)==Array{OIdata,2}
        data = vec(data)
    end
    fig = figure("Visibility phase data",figsize=(10,5),facecolor="White");
    clf();
    ax=gca();
    if (color == "baseline" || color =="base")
        baseline_list_vis = [get_baseline_names(data[n].sta_name,data[n].vis_sta_index) for n=1:length(data)];
        baseline=sort(unique(vcat(baseline_list_vis...)))
        for i=1:length(baseline)
            loc =  [findall(baseline_list_vis[n] .== baseline[i]) for n=1:length(data)]
            baseline_vis = vcat([data[n].vis_baseline[loc[n]] for n=1:length(data)]...)/1e6
            visphi = vcat([data[n].visphi[loc[n]] for n=1:length(data)]...)
            visphi_err = vcat([data[n].visphi_err[loc[n]] for n=1:length(data)]...)
            errorbar(baseline_vis,visphi,yerr=visphi_err,fmt="o",markeredgecolor=oiplot_colors[i],color=oiplot_colors[i], markersize=3,ecolor="Gainsboro",elinewidth=1.0,label=baseline[i])
        end
        if legend_below == false
            ax.legend(fontsize=8, fancybox=true, shadow=true, ncol=4,loc="best")
        else
            ax.legend(fontsize=8, fancybox=true, shadow=true, ncol=8,loc="upper center", bbox_to_anchor=(0.5, -0.15))
        end
    elseif (color == "wavelength" || color == "wav")
        wavcol = vcat([data[n].uv_lam[data[n].indx_vis]*1e6 for n=1:length(data)]...)
        baseline_vis = vcat([data[n].vis_baseline for n=1:length(data)]...)/1e6
        visphi = vcat([data[n].visphi for n=1:length(data)]...);
        visphi_err = vcat([data[n].visphi_err for n=1:length(data)]...);
        sc = scatter(baseline_vis, visphi, c=wavcol, cmap="Spectral_r", alpha=1.0, s=6.0, zorder=100)
        el = errorbar(baseline_vis, visphi,yerr=visphi_err,fmt="none", marker="none",ecolor="Gainsboro", elinewidth=1.0, zorder=0)
        cbar = colorbar(sc, aspect=50, orientation="horizontal", label="Wavelength (μm)", pad=0.18, fraction=0.05)
        cbar_range = floor.(collect(range(minimum(wavcol), maximum(wavcol), length=11))*100)/100
        cbar.set_ticks(cbar_range)
    else
        baseline_vis = vcat([data[n].vis_baseline for n=1:length(data)]...)/1e6
        visphi = vcat([data[n].visphi for n=1:length(data)]...);
        visphi_err = vcat([data[n].visphi_err for n=1:length(data)]...);
        errorbar(baseline_vis,visphi,yerr=visphi_err,fmt="o", markersize=3,color="Black", ecolor="Gainsboro",elinewidth=1.0)
    end
    title("Visibility phase phase data")
    xlabel(L"Maximum Baseline (M$\lambda$)")
    ylabel("Visibility phase (degrees)")
    ax.grid(true,which="both",color="Grey",linestyle=":")
    tight_layout()
    if filename !=""
        savefig(filename)
    end
    show(block=false)
end



function plot_diffphi(data::Array{OIdata,1}; color="Black",markopt=false, legend_below=false, filename="")
    #
    # Note: this is a special kind of plot, which doesn't follow the classic plotting recipe
    #
    baseline_list_vis = [get_baseline_names(data[n].sta_name,data[n].vis_sta_index) for n=1:length(data)];
    baseline=sort(unique(vcat(baseline_list_vis...)))
    # Creating one subplot per baseline
    fig ,ax=  plt.subplots(num="Differential phase data",nrows=length(baseline), sharex=true,figsize=(10,5),facecolor="White")
    suptitle("Differential phase data")
    subplots_adjust(hspace=0.0)
    mx=matplotlib[:ticker][:MultipleLocator](20)
    for i=1:length(baseline)
        title(baseline[i], x=0.9, y=0.75)
        loc =  [findall(baseline_list_vis[n] .== baseline[i]) for n=1:length(data)]
        baseline_vis = vcat([data[n].vis_baseline[loc[n]] for n=1:length(data)]...)/1e6
        wavcol = vcat([data[n].uv_lam[data[n].indx_vis[loc[n]]]*1e9 for n=1:length(data)]...)
        visphi = vcat([data[n].visphi[loc[n]] for n=1:length(data)]...)
        visphi_err = vcat([data[n].visphi_err[loc[n]] for n=1:length(data)]...)
        plt.axes(ax[i])
        ax[i][:xaxis][:set_minor_locator](mx)
        errorbar(wavcol,visphi,yerr=visphi_err,fmt="o",markersize=0.5,ecolor="Gainsboro",elinewidth=.5)
        if i==length(baseline)
            xlabel("λ (nm)")
        end
        ylabel("Δϕ (°)")
        grid(true,which="both",color="Grey",linestyle=":")
    end
    ax[length(baseline)][:tick_params](axis="x", which="major", length=10.0)
    ax[length(baseline)][:tick_params](axis="x", which="minor", length=5.0)
    tight_layout();
    if filename !=""
        savefig(filename)
    end
    show(block=false)
end


function plot_v2_multifile(data::Array{OIdata,1}; logplot = false, remove = false,idpoint=false,clean=false,filename="",legend_below=true)
    global v2base=[]
    global v2value=[]
    global v2err=[]
    global clickbase=Array{Int64,2}
    global clickname=[]
    global clickmjd=[]
    global clicklam=[]
    global clickdlam=[]
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
    #if remove == true # No support for removing points across multimple nights just yet, although this could be VERY useful....
    #    fig.canvas.mpl_connect("button_press_event",onclickv2)
    #end
    for i=1:length(data) #plot each night
        nv2=data[i].nv2
        sta_name=data[i].sta_name
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
            clickname=vcat(clickname,data[i].sta_name)
            clickmjd=vcat(clickmjd,data[i].v2_mjd)
            clicklam=vcat(clicklam,data[i].v2_lam)
            clickdlam=vcat(clickdlam,data[i].v2_dlam)
            basearray=Array{String,1}(undef,data[i].nv2)
            basearray[1:data[i].nv2].=data[i].filename
            if i==1
                clickfile=cat(basearray,dims=1)
            else
                clickfile=cat(clickfile,basearray,dims=1)
            end
        end
        baseline_list=get_baseline_names(sta_name,v2_sta_index)
        for j=1:length(unique(baseline_list))
            baseline=unique(baseline_list)[j]
            loc=findall(baseline_list->baseline_list==baseline,baseline_list)
            errorbar(baseline_v2[loc]/1e6,v2_data[loc],yerr=v2_data_err[loc],fmt="o",marker=oiplot_markers[j], markeredgecolor=oiplot_colors[i],color=oiplot_colors[i],markersize=3,ecolor="Gainsboro",elinewidth=1.0,label=baseline)
            if axiscount==0
                if (length(unique(baseline_list)))==15
                    ax.legend(fontsize=8, fancybox=true, shadow=true, ncol=8,loc="upper center", bbox_to_anchor=(0.5, -0.10));
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
    ax.grid(true,which="both",color="Grey",linestyle=":")
    tight_layout()
    
    #if idpoint==true
    #    cid=fig.canvas.mpl_connect("button_press_event",onclickidentify)
    #end
    if filename !=""
        savefig(filename)
    end
    show(block=false)
end



function plot_v2_and_t3phi_wav(data::Union{OIdata,Array{OIdata,1},Array{OIdata,2}};figsize=(10,10), logplot = false, clean=true,t3base="max",markopt=false, legend_below=false, figtitle="")
    set_oiplot_defaults()
    if typeof(data)==OIdata
            data = [data]
        end
        if typeof(data)==Array{OIdata,2}
            data = vec(data)
        end
        fig = figure(string(figtitle, "V2 & T3phi data"),figsize=figsize,facecolor="White");
        if clean == true # do not overplot on existing window by default
            clf();
        end
    #
    # V2
    #
        axes = fig.subplots(2,1)
        subplot(2,1,1)   
        ax = gca();
        if logplot==true
                ax.set_yscale("log")
            end
        
            xlabel(L"Baseline (M$\lambda$)")
            ylabel("Squared Visibility Amplitudes")
            grid(true,which="both",color="Grey", linestyle=":")
            tight_layout() 
            wavcol = vcat([data[n].uv_lam[data[n].indx_v2]*1e6 for n=1:length(data)]...)
            baseline_v2 = vcat([data[n].v2_baseline for n=1:length(data)]...)/1e6
            v2 = vcat([data[n].v2 for n=1:length(data)]...);
            v2_err = vcat([data[n].v2_err for n=1:length(data)]...);
            sc = scatter(baseline_v2,v2, c=wavcol, cmap="Spectral_r", alpha=1.0, s=6.0, zorder=100)
            el = errorbar(baseline_v2,v2,yerr=v2_err,fmt="none", marker="none",ecolor="Gainsboro", elinewidth=1.0, zorder=0)
            
    #
    # T3phi
    #
        subplot(2,1,2)   
            wavcol = vcat([data[n].uv_lam[data[n].indx_t3_1]*1e6 for n=1:length(data)]...)
            if t3base=="max"
                baseline_t3 = vcat([data[n].t3_maxbaseline for n=1:length(data)]...)/1e6
            elseif t3base=="geom"
                baseline_t3 = vcat([data[n].t3_baseline for n=1:length(data)]...)/1e6
            end
            t3phi = vcat([data[n].t3phi for n=1:length(data)]...);
            t3phi_err = vcat([data[n].t3phi_err for n=1:length(data)]...);
            sc = scatter(baseline_t3, t3phi, c=wavcol, cmap="Spectral_r", alpha=1.0, s=6.0, zorder=100)
            el = errorbar(baseline_t3, t3phi,yerr=t3phi_err,fmt="none", marker="none",ecolor="Gainsboro", elinewidth=1.0, zorder=0)
           
        if t3base=="max"
            xlabel(L"Maximum Baseline (M$\lambda$)")
        elseif t3base=="geom"
            xlabel(L"Geometric Mean Baseline (M$\lambda$)")
        end
        ylabel("Closure phase (degrees)")
       # ax = gca();
        grid(true,which="both",color="Grey",linestyle=":")
        tight_layout()
        cbar = colorbar(sc, aspect=50, orientation="horizontal", label="Wavelength (μm)", pad=0.18, fraction=0.05)
        minrange = ceil(minimum(wavcol)*100)/100
        maxrange = floor(maximum(wavcol)*100)/100
        step = round((maxrange-minrange)*10)/100
        cbar_range = collect(range(minrange, maxrange, step=step))
        cbar.set_ticks(cbar_range)
        show(block=false)
end
    



function imdisp(image; figtitle="OITOOLS image", colormap = "gist_heat", pixsize = -1.0, tickinterval = 0.5, use_colorbar = false, beamsize = -1, beamlocation = [])
    fig = figure(figtitle,figsize=(6,6),facecolor="White")
    clf();
    nx=ny=-1;
    pixmode = false;
    if pixsize == -1
        pixmode = true;
        pixsize = 1
    end
    scaling_factor = maximum(image);
    if abs.(scaling_factor) <  1e-20
        scaling_factor = 1.0;
        @warn("Maximum of image < tol");
    end

    img = []
    if ndims(image) ==1
        ny=nx=Int64(sqrt(length(image)))
        img = imshow(rotl90(reshape(image,nx,nx))/scaling_factor, ColorMap(colormap), interpolation="none", extent=[0.5*nx*pixsize,-0.5*nx*pixsize, -0.5*ny*pixsize,0.5*ny*pixsize]); # uses Monnier's orientation
    else
        nx,ny = size(image);
        img = imshow(rotl90(image)/scaling_factor, ColorMap(colormap), interpolation="none", extent=[0.5*nx*pixsize,-0.5*nx*pixsize,-0.5*ny*pixsize,0.5*ny*pixsize]); # uses Monnier's orientation
    end
    if pixmode == false
        xlabel("x ← E (mas)")
        ylabel("y → N (mas)")
    end

    ax = gca()
    ax.set_aspect("equal")
    if (size(image,1)*pixsize)>1000
        tickinterval = 50.0
    elseif (size(image,1)*pixsize)>100
        tickinterval = 5.0
        # ax.invert_xaxis();
    end
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
        if beamlocation == []
            beamlocation = [.8, .8]
        end
        c = matplotlib.patches.Circle((0.5*nx*pixsize*beamlocation[1],-0.5*ny*pixsize*beamlocation[2]),0.5*beamsize,fc="white",ec="white",linewidth=0)
        ax.add_artist(c)
    end
    tight_layout()
    show(block=false)
end

#TODO: work for rectangular
function imdisp_polychromatic(image_vector::Union{Array{Float64,1}, Array{Float64,2},Array{Float64,3}}; wavs = [], figtitle="Polychromatic image", nwavs = 1, colormap = "gist_heat", pixsize = -1.0, tickinterval = 10, use_colorbar = false, beamsize = -1, beamlocation = [.9, .9])
    if typeof(image_vector)==Array{Float64,2}
        nwavs = size(image_vector,2)
    elseif typeof(image_vector)==Array{Float64,3}
        nwavs = size(image_vector,3)
    end
    nside = ceil(Int64,sqrt(nwavs))

    fig = figure(figtitle,figsize=(10,10),facecolor="White")
    clf();
    images_all =reshape(image_vector, (div(length(vec(image_vector)),nwavs), nwavs))
    for i=1:nwavs
        fig.add_subplot(nside,nside,i)
        if wavs==[]
            title("λ: $i")
        else
            title("λ: $(wavs[i])")
        end

        image = images_all[:,i]
        nx=ny=-1;
        pixmode = false;
        if pixsize == -1
            pixmode = true;
            pixsize = 1
        end
        ny=nx=Int64.(sqrt(length(image)))
        img = imshow(rotl90(reshape(image,nx,nx))/maximum(image), ColorMap(colormap), interpolation="none", extent=[0.5*nx*pixsize,-0.5*nx*pixsize,-0.5*ny*pixsize,0.5*ny*pixsize]); # uses Monnier's orientation
        # if pixmode == false
        #     xlabel("RA (mas)")
        #     ylabel("DEC (mas)")
        # end
        #
        # ax = gca()
        # ax.set_aspect("equal")
        # mx = matplotlib.ticker.MultipleLocator(tickinterval) # Define interval of minor ticks
        # ax.xaxis.set_minor_locator(mx) # Set interval of minor ticks
        # ax.yaxis.set_minor_locator(mx) # Set interval of minor ticks
        # ax.xaxis.set_tick_params(which="major",length=5,width=2)
        # ax.xaxis.set_tick_params(which="minor",length=2,width=1)
        # ax.yaxis.set_tick_params(which="major",length=5,width=2)
        # ax.yaxis.set_tick_params(which="minor",length=2,width=1)
        #
        # if use_colorbar == true
        #     divider = pyimport("mpl_toolkits.axes_grid1").make_axes_locatable(ax)
        #     cax = divider.append_axes("right", size="5%", pad=0.2)
        #     colorbar(img, cax=cax, ax=ax)
        # end
        #
        # if beamsize > 0
        #     c = matplotlib.patches.Circle((0.5*nx*pixsize*beamlocation[1],-0.5*ny*pixsize*beamlocation[2]),beamsize,fc="white",ec="white",linewidth=.5)
        #     ax.add_artist(c)
        # end
        tight_layout()
    end
    show(block=false)
end

# TODO: rework for julia 1.0+
function imdisp_temporal(image_vector, nepochs; colormap = "gist_heat", name="Time-variable images",pixsize = -1.0, tickinterval = 10, use_colorbar = false, beamsize = -1, beamlocation = [.9, .9])
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
        if pixsize == -1
            pixmode = true;
            pixsize = 1
        end
        img = []
        if ndims(image) ==1
            ny=nx=Int64(sqrt(length(image)))
            img = imshow(rotl90(reshape(image,nx,nx)), ColorMap(colormap), interpolation="none", extent=[-0.5*nx*pixsize,0.5*nx*pixsize,-0.5*ny*pixsize,0.5*ny*pixsize]); # uses Monnier's orientation
        else
            nx,ny = size(image);
            img = imshow(rotl90(image), ColorMap(colormap), interpolation="none", extent=[-0.5*nx*pixsize,0.5*nx*pixsize,-0.5*ny*pixsize,0.5*ny*pixsize]); # uses Monnier's orientation
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
            c = matplotlib.patches.Circle((0.5*nx*pixsize*beamlocation[1],-0.5*ny*pixsize*beamlocation[2]),beamsize,fc="white",ec="white",linewidth=.5)
            ax.add_artist(c)
        end
        tight_layout()
    end
    show(block=false)
end


function get_baseline_names(sta_names,sta_indx)
    nbaselines = size(sta_indx,2)
    baseline_names=Array{String}(undef,nbaselines)
    for i=1:nbaselines
        baseline_names[i]=string(sta_names[sta_indx[1,i]],"-",sta_names[sta_indx[2,i]])
    end
    return baseline_names
end

function get_triplet_names(sta_names, sta_indx)
    nt3=size(sta_indx,2)
    triplet_names=Array{String}(undef,nt3)
    for i=1:nt3
        triplet_names[i]=string(sta_names[sta_indx[1,i]],"-",sta_names[sta_indx[2,i]],"-",sta_names[sta_indx[3,i]])
    end
    return triplet_names
end

function get_triple_baselines_names(sta_names, sta_indx)
    nt3=size(sta_indx,2)
    baseline_names=Array{String}(undef, 3, nt3)
    for i=1:nt3
        baseline_names[1, i]=string(sta_names[sta_indx[1,i]],"-",sta_names[sta_indx[2,i]])
        baseline_names[2, i]=string(sta_names[sta_indx[2,i]],"-",sta_names[sta_indx[3,i]])
        baseline_names[3, i]=string(sta_names[sta_indx[1,i]],"-",sta_names[sta_indx[3,i]])
    end
    return baseline_names
end


function plot_facility(facility)
    coords=reshape(facility.sta_xyz, 3, facility.ntel[1])
    scatter(coords[1,:], coords[2,:])
    axis("equal")
    for i=1:facility.ntel[1]
        text(coords[1,i]+1, coords[2,i]+1, facility.sta_names[i])
    end
end
    