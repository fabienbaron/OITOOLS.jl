# gather common display tasks
using PyPlot,PyCall

#@pyimport mpl_toolkits.axes_grid1 as axgrid

# double check by plotting uv coverage

function uvplot(u,v)
fig = figure("UV plot",figsize=(10,10),facecolor="White")
ax = axes()
scatter(u, v,alpha=0.5)
scatter(-u, -v,alpha=0.5)
title("UV coverage")
xlabel("U")
ylabel("V")
grid("on")
PyPlot.draw();PyPlot.pause(0.5); # this is used to see plots when running code in batch mode
end


function v2plot_modelvsdata(baseline_v2,v2_data,v2_data_err, v2_model) #plots V2 data vs v2 model
fig = figure("V2 plot - Model vs Data",figsize=(10,10),facecolor="White")
subplot(211)
errorbar(baseline_v2,v2_data,yerr=v2_data_err,fmt="o", color="Red")
scatter(baseline_v2, v2_model, color="Blue")
title("V2 - Model vs data plot")
xlabel("Baseline")
ylabel("V2")
subplot(212)
scatter(baseline_v2, (v2_model - v2_data)./v2_data_err)
xlabel("Baseline")
ylabel("Residuals (number of sigma)")
PyPlot.show();PyPlot.pause(1);  # this is used to see plots when running code in batch mode
end


function v2plot(baseline_v2,v2_data,v2_data_err) # plots v2 data only
  fig = figure("V2 data",figsize=(10,10),facecolor="White")
  errorbar(baseline_v2,v2_data,yerr=v2_data_err,fmt="o", color="Red")
  title("V2 data")
  xlabel("Baseline")
  ylabel("V2")
  PyPlot.show();PyPlot.pause(0.5);  # this is used to see plots when running code in batch mode
end


function t3phiplot(baseline_t3,t3phi_data,t3phi_data_err) # plots v2 data only
  fig = figure("Closure phase data",figsize=(10,10),facecolor="White")
  errorbar(baseline_t3,t3phi_data,yerr=t3phi_data_err,fmt="o", color="Red")
  title("Closure phase data")
  xlabel("Baseline")
  ylabel("Closure phase (degrees)")
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
 PyPlot.draw();PyPlot.pause(0.05);
end
