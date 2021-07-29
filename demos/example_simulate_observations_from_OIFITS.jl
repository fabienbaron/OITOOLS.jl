using OITOOLS
#
#EXAMPLE 9
#Simulate an observation using an input oifits file as base
# The uv coverage and SNR of data points will be copied

image_file="./data/2004true.fits"
pixsize=0.101

in_oifits  = "./data/2004-data1.oifits"
out_oifits  = "./data/2004-simulated.oifits"

simulate_from_oifits(in_oifits,out_oifits,image=image_file,pixsize=pixsize)

#Compare simulated data to input
data1 = readoifits(in_oifits)
data2 = readoifits(out_oifits)
