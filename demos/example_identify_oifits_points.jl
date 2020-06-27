#
# Colored plots by baselines, left click identifies the point
#
using OITOOLS
oifitsin = "./data/AZCYG2011_FINAL2018.oifits"
data1 = readoifits(oifitsin,filter_bad_data=true)[1,1];

oifitsin ="./data/AZCYG2014NEW.oifits"
data2 = readoifits(oifitsin,filter_bad_data=true)[1,1];

merged_data =[data1,data2];

#plot one v2 with baselines colored logscale --left click identifies the point
v2plot(data1,logplot=true, idpoint = true, legend_below=true)

#plot two files by baseline, with different colors for each and different symbols for baselines, and id points on click
#this also indicates from which file the point comes from
v2plot_multifile(merged_data, logplot = true, idpoint=true)
