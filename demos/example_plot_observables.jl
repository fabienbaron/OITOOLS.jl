using OITOOLS
#
# EXAMPLE 1: read image and data, and compare the observables
#
# read the data file
oifitsfile = "./data/rho_Cas_example.oifits"
data = readoifits(oifitsfile)
# Display the data
uvplot(data, color="bases");
uvplot(data, color="wav");
uvplot(data, color="mjd");

plot_v2(data,color="bases")
plot_v2(data,color="wav")
plot_t3phi(data);
plot_t3phi(data, t3base="geom");

#plot_v2_and_t3phi_wav(data, logplot=true, figsize=(14,14))
