# gather common display tasks
using PyPlot,PyCall

#@pyimport mpl_toolkits.axes_grid1 as axgrid

# double check by plotting uv coverage

function uvplot(uv)
u = uv[1,:]
v = uv[2,:]
fig = figure("UV plot",figsize=(10,10),facecolor="White")
ax = axes()
scatter(u, v,alpha=0.5, color="Black")
scatter(-u, -v,alpha=0.5, color="Black")
title("UV coverage")
xlabel("U")
ylabel("V")
grid("on")
tight_layout();
PyPlot.draw();PyPlot.pause(0.5); # this is used to see plots when running code in batch mode
end


function v2plot_modelvsdata(baseline_v2,v2_data,v2_data_err, v2_model; logplot = false) #plots V2 data vs v2 model
fig = figure("V2 plot - Model vs Data",figsize=(10,10),facecolor="White")
subplot(211)
ax = gca()
if logplot==true
ax[:set_yscale]("log")
end
errorbar(baseline_v2,v2_data,yerr=v2_data_err,fmt="o", markersize=2,color="Black")
plot(baseline_v2, v2_model, color="Red", linestyle="none", marker="o", markersize=2)
title("V2 - Model vs data plot")
xlabel("Baseline")
ylabel("V2")
grid("on")
subplot(212)
plot(baseline_v2, (v2_model - v2_data)./v2_data_err,color="Black", linestyle="none", marker="o", markersize=2)
xlabel("Baseline")
ylabel("Residuals (number of sigma)")
grid("on");
tight_layout()
PyPlot.show();PyPlot.pause(0.05);  # this is used to see plots when running code in batch mode
end


function v2plot(baseline_v2,v2_data,v2_data_err) # plots v2 data only
fig = figure("V2 data",figsize=(15,8),facecolor="White")
  errorbar(baseline_v2,v2_data,yerr=v2_data_err,fmt="o", color="Red")
  title("V2 data")
  xlabel("Baseline")
  ylabel("V2")
  grid("on")
  tight_layout()
  PyPlot.show();PyPlot.pause(0.05);  # this is used to see plots when running code in batch mode
end


function t3phiplot(baseline_t3,t3phi_data,t3phi_data_err) # plots v2 data only
  fig = figure("Closure phase data",figsize=(10,10),facecolor="White")
  errorbar(baseline_t3,t3phi_data,yerr=t3phi_data_err,fmt="o", color="Red")
  title("Closure phase data")
  xlabel("Baseline")
  ylabel("Closure phase (degrees)")
  tight_layout()
  PyPlot.show();PyPlot.pause(0.5);  # this is used to see plots when running code in batch mode
end


function imdisp(image, pixellation = -1)
 fig = figure("Image",figsize=(10,10),facecolor="White")
# if pixellation < 0 -> no pixellation entered -> do not draw in milliarcseconds
 nx=Int64(sqrt(length(image)))
 #ax = gca()
 imshow(rotl90(reshape(image,nx,nx)), ColorMap("hot"), interpolation="none"); # uses Monnier's orientation
 #divider = axgrid.make_axes_locatable(ax)
 #cax = divider[:append_axes]("right", size="5%", pad=0.05)
 #colorbar(image, cax=cax)
 tight_layout()
 PyPlot.draw();PyPlot.pause(0.05);
end
