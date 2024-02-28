using OITOOLS
using PyPlot
data=readoifits("./data/MYSTIC_L2.2023Oct25.rho_Cas.MIRCX_IDL.lbd_all.AVG10m.oifits")
plot_v2_and_t3phi_wav(data, logplot=true, figsize=(14,14))
savefig("plot.svg")
#plot_v2(data, color="wav", logplot=true)
#plot_t3phi(data, color="wav")
uvplot(data, color="baseline")
uvplot(data, color="mjd")
uvplot(data, color="wav")

