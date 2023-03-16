using OITOOLS, PyPlot
#
#  Simple diameter vs linear darkening coefficient grid search
#  This shows the correlation between the two parameters
#

oifitsfile = "./data/AlphaCenA.oifits";
data = (readoifits(oifitsfile))[1,1]; # data can be split by wavelength, time, etc.
diameters=range(1.0,10.0, step=0.1) #mas
ld=range(0.0,1.0, step=0.01) # limb darkening
chi2 = zeros(length(diameters), length(ld))
model = create_model(create_component(type="ldlin", name="Model"));
for i= 1:length(diameters)
    for j=1:length(ld)
     chi2[i,j] = model_to_chi2(data, model, [diameters[i], ld[j]], weights=[1.0,0,0])
 end
end
res = findmin(chi2)
print("Best chi2: ", res[1], " Diameter:", diameters[res[2][1]], " LD:", ld[res[2][2]])
# Plot the chi2 map
imshow(chi2.^.05, ColorMap("gist_heat_r"), interpolation="none")
xlabel("Limb Darkening Coefficient")
ylabel("Diameter (mas)")
x_ticks = collect(range(1, length(ld), step=Int(10)))
y_ticks = collect(range(1, length(diameters), step=Int(10)))
xticks(x_ticks, ld[x_ticks])
yticks(y_ticks, diameters[y_ticks])
ax=gca()
ax[:tick_params](axis="x", which="major", length=8.0)
ax[:tick_params](axis="y", which="major", length=8.0)

cbar = colorbar(ax=ax, aspect=50, orientation="vertical", label=L"Reduced $\chi^2$ values")
