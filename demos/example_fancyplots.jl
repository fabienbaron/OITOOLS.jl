oifitsin = "./data/AZCYG2011_FINAL2018.oifits"
data1 = readoifits(oifitsin,filter_bad_data=true)[1,1];

oifitsin ="./data/AZCYG2014NEW.oifits"
data2 = readoifits(oifitsin,filter_bad_data=true)[1,1];

data =[data1,data2]

#plot one file uv by baseline colorful
uvplot(data1,fancy=true,filename="")
#plot one v2 with baselines colored logscale an
v2plot(data1;logplot=true,fancy=true)
#plot two files by baseline with different colors for each and different symbols for baselines and id points on click
v2plot_multifile(data, logplot = true, idpoint=true,filename="")

#plot t3phi with baselines colored
t3phiplot(data1;fancy=true,filename="")
