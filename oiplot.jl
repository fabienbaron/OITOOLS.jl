# gather common display tasks
using PyPlot,PyCall
PyDict(pyimport("matplotlib")["rcParams"])["font.family"]=["serif"]
#PyDict(pyimport("matplotlib")["rcParams"])["mathtext.fontset"]=["custom"]
PyDict(pyimport("matplotlib")["rcParams"])["xtick.major.size"]=[6]
PyDict(pyimport("matplotlib")["rcParams"])["ytick.major.size"]=[6]
PyDict(pyimport("matplotlib")["rcParams"])["xtick.minor.size"]=[6]
PyDict(pyimport("matplotlib")["rcParams"])["ytick.minor.size"]=[6]
PyDict(pyimport("matplotlib")["rcParams"])["xtick.major.width"]=[1]
PyDict(pyimport("matplotlib")["rcParams"])["ytick.major.width"]=[1]
PyDict(pyimport("matplotlib")["rcParams"])["xtick.minor.width"]=[1]
PyDict(pyimport("matplotlib")["rcParams"])["ytick.minor.width"]=[1]
PyDict(pyimport("matplotlib")["rcParams"])["lines.markeredgewidth"]=[1]
PyDict(pyimport("matplotlib")["rcParams"])["legend.numpoints"]=[1]
#PyDict(pyimport("matplotlib")["rcParams"])["legend.frameon"]=["False"]
PyDict(pyimport("matplotlib")["rcParams"])["legend.handletextpad"]=[0.3]
#@pyimport mpl_toolkits.axes_grid1 as axgrid

# double check by plotting uv coverage

function uvplot(uv)
u = uv[1,:]/1e6
v = uv[2,:]/1e6
fig = figure("UV plot",figsize=(10,10),facecolor="White")
clf();
ax = axes()
minorticks_on
markeredgewidth=0.1
ax[:locator_params](axis ="y", nbins=20)
ax[:locator_params](axis ="x", nbins=20)
scatter(u, v,alpha=0.5, color="Black")
scatter(-u, -v,alpha=0.5, color="Black")
title("UV coverage")
xlabel(L"U (M$\lambda$)")
ylabel(L"V (M$\lambda$)")
grid("on")
tight_layout();
PyPlot.draw();PyPlot.pause(0.5); # this is used to see plots when running code in batch mode
end


function v2plot_modelvsdata(baseline_v2,v2_data,v2_data_err, v2_model; logplot = false) #plots V2 data vs v2 model
fig = figure("V2 plot - Model vs Data",figsize=(10,10),facecolor="White")
clf();
subplot(211)
ax = gca();
if logplot==true
ax[:set_yscale]("log")
end
errorbar(baseline_v2/1e6,v2_data,yerr=v2_data_err,fmt="o", markersize=2,color="Black")
plot(baseline_v2/1e6, v2_model, color="Red", linestyle="none", marker="o", markersize=3)
title("Squared Visbility Amplitudes - Model vs data plot")
#xlabel(L"Baseline (M$\lambda$)")
ylabel("Squared Visibility Amplitudes")
grid("on")
subplot(212)
plot(baseline_v2/1e6, (v2_model - v2_data)./v2_data_err,color="Black", linestyle="none", marker="o", markersize=3)
xlabel(L"Baseline (M$\lambda$)")
ylabel("Residuals (number of sigma)")
grid("on");
tight_layout()
PyPlot.show();PyPlot.pause(0.05);  # this is used to see plots when running code in batch mode
end


function v2plot(baseline_v2,v2_data,v2_data_err; logplot = false) # plots v2 data only
fig = figure("V2 data",figsize=(10,5),facecolor="White");
clf();
ax = gca();
if logplot==true
ax[:set_yscale]("log")
end
errorbar(baseline_v2/1e6,v2_data,yerr=v2_data_err,fmt="o", markersize=3,color="Black")
title("Squared Visibility Amplitude Data")
xlabel(L"Baseline (M$\lambda$)")
ylabel("Squared Visibility Amplitudes")
grid("on")
tight_layout()
PyPlot.show();PyPlot.pause(0.05);  # this is used to see plots when running code in batch mode
end

function t3phiplot(baseline_t3,t3phi_data,t3phi_data_err) # plots v2 data only
  fig = figure("Closure phase data",figsize=(10,5),facecolor="White");
  clf();
  errorbar(baseline_t3/1e6,t3phi_data,yerr=t3phi_data_err,fmt="o", markersize=3,color="Black")
  title("Closure phase data")
  xlabel(L"Maximum Baseline (M$\lambda$)")
  ylabel("Closure phase (degrees)")
  grid("on")
  tight_layout()
  PyPlot.show();PyPlot.pause(0.5);  # this is used to see plots when running code in batch mode
end


function imdisp(image, pixellation = -1)
 fig = figure("Image",figsize=(5,5),facecolor="White")
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
